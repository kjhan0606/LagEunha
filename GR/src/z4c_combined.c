/**
 * @file z4c_combined.c
 * @brief Combined Z4c DG-FVM Implementation
 * 
 * All source files combined for easy compilation.
 * Compile with: gcc -O3 -o z4c_dg z4c_combined.c -lm
 * 
 * References:
 * [1] Baumgarte & Shapiro, "Numerical Relativity", Cambridge (2010)
 * [2] Dumbser et al., J. Comput. Phys. 348, 70-117 (2017)
 * [3] Hilditch et al., Phys. Rev. D 88, 084057 (2013)
 * [4] Bernuzzi & Hilditch, Phys. Rev. D 81, 084003 (2010)
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * CONFIGURATION AND CONSTANTS
 * ============================================================================*/
#define DIM 3
#define N_GAUSS 5
#define MAX_FACES 24
#define MAX_NEIGHBORS 24
#define EPS 1.0e-14

#define N_Z4C_VARS 22
#define N_AUX_VARS 33
#define N_TOTAL_VARS (N_Z4C_VARS + N_AUX_VARS)

/* Variable indices - Ref: [3] Hilditch et al. */
#define IDX_PHI 0
#define IDX_GT_XX 1
#define IDX_GT_XY 2
#define IDX_GT_XZ 3
#define IDX_GT_YY 4
#define IDX_GT_YZ 5
#define IDX_GT_ZZ 6
#define IDX_K 7
#define IDX_AT_XX 8
#define IDX_AT_XY 9
#define IDX_AT_XZ 10
#define IDX_AT_YY 11
#define IDX_AT_YZ 12
#define IDX_AT_ZZ 13
#define IDX_GTFUNC_X 14
#define IDX_GTFUNC_Y 15
#define IDX_GTFUNC_Z 16
#define IDX_THETA 17
#define IDX_ALPHA 18
#define IDX_BETA_X 19
#define IDX_BETA_Y 20
#define IDX_BETA_Z 21
#define IDX_AUX_A_X 22
#define IDX_AUX_A_Y 23
#define IDX_AUX_A_Z 24
#define IDX_AUX_P_X 52
#define IDX_AUX_P_Y 53
#define IDX_AUX_P_Z 54

/* ============================================================================
 * DATA STRUCTURES
 * ============================================================================*/
typedef struct { double x, y, z; } Vec3;

typedef struct {
    double kappa1, kappa2, kappa3, eta;
} Z4cDampingParams;

typedef enum { GAUGE_1_PLUS_LOG, GAUGE_HARMONIC } LapseGauge;
typedef enum { SHIFT_ZERO, SHIFT_GAMMA_DRIVER } ShiftGauge;

typedef struct {
    LapseGauge lapse_type;
    ShiftGauge shift_type;
    double mu_L, mu_S, eta_beta;
} GaugeParams;

typedef struct {
    int neighbor_id;
    Vec3 centroid, normal;
    double area;
} VoronoiFace;

typedef struct {
    int id;
    Vec3 centroid;
    double volume;
    int n_faces;
    VoronoiFace faces[MAX_FACES];
    int n_neighbors;
    int neighbors[MAX_NEIGHBORS];
    double U[N_TOTAL_VARS];
} VoronoiCell;

typedef struct {
    int n_cells;
    VoronoiCell *cells;
    Vec3 bbox_min, bbox_max;
    double min_cell_size;
} VoronoiMesh;

typedef struct {
    int n_points;
    double points[N_GAUSS], weights[N_GAUSS];
} GaussQuadrature;

typedef struct {
    Vec3 domain_min, domain_max;
    int nx, ny, nz;
    double cfl, t_final, dt_output;
    Z4cDampingParams damping;
    GaugeParams gauge;
    char output_dir[256];
} SimConfig;

typedef struct {
    SimConfig config;
    VoronoiMesh mesh;
    GaussQuadrature gauss;
    double time, dt;
    int step;
    double *rhs;
    double hamiltonian_L2, z4_constraint_L2;
} Simulation;

/* ============================================================================
 * VECTOR UTILITIES
 * ============================================================================*/
double vec3_dot(const Vec3 *a, const Vec3 *b) { return a->x*b->x + a->y*b->y + a->z*b->z; }
double vec3_norm(const Vec3 *a) { return sqrt(vec3_dot(a, a)); }

/* ============================================================================
 * SYMMETRIC TENSOR UTILITIES - Ref: [1] Appendix A
 * ============================================================================*/
int sym_idx(int i, int j) {
    if (i > j) { int tmp = i; i = j; j = tmp; }
    static const int idx_map[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
    return idx_map[i][j];
}

void sym_to_full(double full[3][3], const double sym[6]) {
    full[0][0] = sym[0]; full[0][1] = sym[1]; full[0][2] = sym[2];
    full[1][0] = sym[1]; full[1][1] = sym[3]; full[1][2] = sym[4];
    full[2][0] = sym[2]; full[2][1] = sym[4]; full[2][2] = sym[5];
}

double det_sym(const double sym[6]) {
    double full[3][3]; sym_to_full(full, sym);
    return full[0][0]*(full[1][1]*full[2][2]-full[1][2]*full[2][1])
         - full[0][1]*(full[1][0]*full[2][2]-full[1][2]*full[2][0])
         + full[0][2]*(full[1][0]*full[2][1]-full[1][1]*full[2][0]);
}

void inv_sym(double inv[6], const double sym[6]) {
    double full[3][3], inv_full[3][3];
    sym_to_full(full, sym);
    double det = det_sym(sym);
    double inv_det = 1.0 / (det + EPS);
    
    inv_full[0][0] = (full[1][1]*full[2][2] - full[1][2]*full[2][1]) * inv_det;
    inv_full[0][1] = (full[0][2]*full[2][1] - full[0][1]*full[2][2]) * inv_det;
    inv_full[0][2] = (full[0][1]*full[1][2] - full[0][2]*full[1][1]) * inv_det;
    inv_full[1][1] = (full[0][0]*full[2][2] - full[0][2]*full[2][0]) * inv_det;
    inv_full[1][2] = (full[0][2]*full[1][0] - full[0][0]*full[1][2]) * inv_det;
    inv_full[2][2] = (full[0][0]*full[1][1] - full[0][1]*full[1][0]) * inv_det;
    
    inv[0] = inv_full[0][0]; inv[1] = inv_full[0][1]; inv[2] = inv_full[0][2];
    inv[3] = inv_full[1][1]; inv[4] = inv_full[1][2]; inv[5] = inv_full[2][2];
}

/* ============================================================================
 * GAUSS-LEGENDRE QUADRATURE - Ref: Hesthaven & Warburton Appendix C
 * ============================================================================*/
double legendre_P(int n, double x) {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    double P_nm2 = 1.0, P_nm1 = x, P_n = 0.0;
    for (int k = 2; k <= n; k++) {
        P_n = ((2.0*k-1.0)*x*P_nm1 - (k-1.0)*P_nm2) / (double)k;
        P_nm2 = P_nm1; P_nm1 = P_n;
    }
    return P_n;
}

double legendre_dP(int n, double x) {
    if (n == 0) return 0.0;
    if (n == 1) return 1.0;
    double P_n = legendre_P(n, x), P_nm1 = legendre_P(n-1, x);
    if (fabs(x*x - 1.0) < EPS) return (x > 0 ? 1 : -1) * 0.5 * n * (n + 1.0);
    return n * (P_nm1 - x*P_n) / (1.0 - x*x);
}

void gauss_legendre_init(GaussQuadrature *gauss, int n) {
    gauss->n_points = n;
    for (int k = 0; k < n; k++) {
        double x = cos(M_PI * (4.0*(k+1) - 1.0) / (4.0*n + 2.0));
        for (int iter = 0; iter < 100; iter++) {
            double P = legendre_P(n, x), dP = legendre_dP(n, x);
            double dx = -P / (dP + EPS);
            x += dx;
            if (fabs(dx) < 1e-15) break;
        }
        gauss->points[k] = x;
        double dP = legendre_dP(n, x);
        gauss->weights[k] = 2.0 / ((1.0 - x*x) * dP * dP);
    }
}

/* ============================================================================
 * VORONOI MESH (CARTESIAN APPROXIMATION)
 * ============================================================================*/
int voronoi_mesh_from_cartesian(VoronoiMesh *mesh, Vec3 domain_min, Vec3 domain_max, int nx, int ny, int nz) {
    mesh->bbox_min = domain_min;
    mesh->bbox_max = domain_max;
    mesh->n_cells = nx * ny * nz;
    
    double dx = (domain_max.x - domain_min.x) / nx;
    double dy = (domain_max.y - domain_min.y) / ny;
    double dz = (domain_max.z - domain_min.z) / nz;
    mesh->min_cell_size = fmin(dx, fmin(dy, dz));
    
    mesh->cells = (VoronoiCell *)calloc(mesh->n_cells, sizeof(VoronoiCell));
    if (!mesh->cells) return -1;
    
    int cell_id = 0;
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                VoronoiCell *cell = &mesh->cells[cell_id];
                cell->id = cell_id;
                cell->centroid.x = domain_min.x + (i + 0.5) * dx;
                cell->centroid.y = domain_min.y + (j + 0.5) * dy;
                cell->centroid.z = domain_min.z + (k + 0.5) * dz;
                cell->volume = dx * dy * dz;
                cell->n_faces = 6;
                
                /* Face 0: x- */
                cell->faces[0].centroid = (Vec3){cell->centroid.x - 0.5*dx, cell->centroid.y, cell->centroid.z};
                cell->faces[0].normal = (Vec3){-1.0, 0.0, 0.0};
                cell->faces[0].area = dy * dz;
                cell->faces[0].neighbor_id = (i > 0) ? cell_id - 1 : -1;
                
                /* Face 1: x+ */
                cell->faces[1].centroid = (Vec3){cell->centroid.x + 0.5*dx, cell->centroid.y, cell->centroid.z};
                cell->faces[1].normal = (Vec3){1.0, 0.0, 0.0};
                cell->faces[1].area = dy * dz;
                cell->faces[1].neighbor_id = (i < nx-1) ? cell_id + 1 : -1;
                
                /* Face 2: y- */
                cell->faces[2].centroid = (Vec3){cell->centroid.x, cell->centroid.y - 0.5*dy, cell->centroid.z};
                cell->faces[2].normal = (Vec3){0.0, -1.0, 0.0};
                cell->faces[2].area = dx * dz;
                cell->faces[2].neighbor_id = (j > 0) ? cell_id - nx : -1;
                
                /* Face 3: y+ */
                cell->faces[3].centroid = (Vec3){cell->centroid.x, cell->centroid.y + 0.5*dy, cell->centroid.z};
                cell->faces[3].normal = (Vec3){0.0, 1.0, 0.0};
                cell->faces[3].area = dx * dz;
                cell->faces[3].neighbor_id = (j < ny-1) ? cell_id + nx : -1;
                
                /* Face 4: z- */
                cell->faces[4].centroid = (Vec3){cell->centroid.x, cell->centroid.y, cell->centroid.z - 0.5*dz};
                cell->faces[4].normal = (Vec3){0.0, 0.0, -1.0};
                cell->faces[4].area = dx * dy;
                cell->faces[4].neighbor_id = (k > 0) ? cell_id - nx*ny : -1;
                
                /* Face 5: z+ */
                cell->faces[5].centroid = (Vec3){cell->centroid.x, cell->centroid.y, cell->centroid.z + 0.5*dz};
                cell->faces[5].normal = (Vec3){0.0, 0.0, 1.0};
                cell->faces[5].area = dx * dy;
                cell->faces[5].neighbor_id = (k < nz-1) ? cell_id + nx*ny : -1;
                
                cell_id++;
            }
        }
    }
    printf("Mesh created: %d cells (%dx%dx%d)\n", mesh->n_cells, nx, ny, nz);
    return 0;
}

void voronoi_mesh_free(VoronoiMesh *mesh) { free(mesh->cells); }

/* ============================================================================
 * Z4c SOURCE TERMS - Ref: [3] Hilditch et al., Eqs. (13)-(23)
 * ============================================================================*/
void z4c_compute_source(double *source, const double *U, const Z4cDampingParams *damp, const GaugeParams *gauge) {
    for (int i = 0; i < N_TOTAL_VARS; i++) source[i] = 0.0;
    
    double phi = U[IDX_PHI];
    double exp_m4phi = exp(-4.0 * phi);
    
    double gt[6] = {U[IDX_GT_XX], U[IDX_GT_XY], U[IDX_GT_XZ], U[IDX_GT_YY], U[IDX_GT_YZ], U[IDX_GT_ZZ]};
    double gt_inv[6]; inv_sym(gt_inv, gt);
    double gt_inv_full[3][3]; sym_to_full(gt_inv_full, gt_inv);
    
    double K = U[IDX_K];
    double At[6] = {U[IDX_AT_XX], U[IDX_AT_XY], U[IDX_AT_XZ], U[IDX_AT_YY], U[IDX_AT_YZ], U[IDX_AT_ZZ]};
    double At_full[3][3]; sym_to_full(At_full, At);
    
    double Theta = U[IDX_THETA];
    double alpha = U[IDX_ALPHA];
    double beta[3] = {U[IDX_BETA_X], U[IDX_BETA_Y], U[IDX_BETA_Z]};
    
    /* Compute tilde_A^{ij} and A_ij A^{ij} */
    double At_sq = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double At_up_ij = 0.0;
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    At_up_ij += gt_inv_full[i][k] * gt_inv_full[j][l] * At_full[k][l];
            At_sq += At_full[i][j] * At_up_ij;
        }
    }
    
    double kappa1 = damp->kappa1, kappa2 = damp->kappa2;
    
    /* Eq. (13): d_t phi = -1/6 alpha K + 1/6 div_beta */
    source[IDX_PHI] = -alpha * K / 6.0;
    
    /* Eq. (14): d_t tilde_gamma_ij = -2 alpha tilde_A_ij + ... */
    for (int s = 0; s < 6; s++)
        source[IDX_GT_XX + s] = -2.0 * alpha * At[s];
    
    /* Eq. (15): d_t K = alpha (A_ij A^ij + K^2/3) + alpha kappa1 (1-kappa2) Theta */
    source[IDX_K] = alpha * (At_sq + K*K/3.0) + alpha * kappa1 * (1.0 - kappa2) * Theta;
    
    /* Eq. (16): d_t tilde_A_ij = alpha K tilde_A_ij - 2 alpha tilde_A_il tilde_A^l_j + ... */
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j <= i; j++) {
            int s = sym_idx(i, j);
            double src = alpha * K * At[s];
            for (int l = 0; l < 3; l++)
                for (int m = 0; m < 3; m++)
                    src -= 2.0 * alpha * At[sym_idx(i,l)] * gt_inv_full[l][m] * At[sym_idx(m,j)];
            source[IDX_AT_XX + s] = src;
        }
    }
    
    /* Eq. (18): d_t Theta = -alpha kappa1 (2 + kappa2) Theta + ... */
    source[IDX_THETA] = -alpha * kappa1 * (2.0 + kappa2) * Theta;
    
    /* Gauge: 1+log lapse - Ref: [1] Eq. (4.60) */
    if (gauge->lapse_type == GAUGE_1_PLUS_LOG)
        source[IDX_ALPHA] = -2.0 * alpha * K;
    
    /* Gamma-driver shift - Ref: [1] Eq. (4.82) */
    if (gauge->shift_type == SHIFT_GAMMA_DRIVER) {
        double Gt[3] = {U[IDX_GTFUNC_X], U[IDX_GTFUNC_Y], U[IDX_GTFUNC_Z]};
        for (int i = 0; i < 3; i++)
            source[IDX_BETA_X + i] = 0.75 * alpha * alpha * Gt[i] - gauge->eta_beta * beta[i];
    }
}

/* ============================================================================
 * Z4c FLUX - Ref: [2] Dumbser et al., Section 2.2
 * ============================================================================*/
void z4c_compute_flux(double *flux, const double *U, int d) {
    for (int i = 0; i < N_TOTAL_VARS; i++) flux[i] = 0.0;
    
    double beta_d = U[IDX_BETA_X + d];
    
    /* Advection flux: beta^d * U */
    flux[IDX_PHI] = beta_d * U[IDX_PHI];
    for (int s = 0; s < 6; s++) flux[IDX_GT_XX + s] = beta_d * U[IDX_GT_XX + s];
    flux[IDX_K] = beta_d * U[IDX_K];
    for (int s = 0; s < 6; s++) flux[IDX_AT_XX + s] = beta_d * U[IDX_AT_XX + s];
    for (int i = 0; i < 3; i++) flux[IDX_GTFUNC_X + i] = beta_d * U[IDX_GTFUNC_X + i];
    flux[IDX_THETA] = beta_d * U[IDX_THETA];
    flux[IDX_ALPHA] = beta_d * U[IDX_ALPHA];
    for (int i = 0; i < 3; i++) flux[IDX_BETA_X + i] = beta_d * U[IDX_BETA_X + i];
}

/* ============================================================================
 * CHARACTERISTIC SPEEDS - Ref: [2] Section 2.3
 * ============================================================================*/
void z4c_characteristic_speeds(double *speeds, const double *U, int d) {
    double alpha = U[IDX_ALPHA];
    double beta_d = U[IDX_BETA_X + d];
    double phi = U[IDX_PHI];
    
    double gt[6] = {U[IDX_GT_XX], U[IDX_GT_XY], U[IDX_GT_XZ], U[IDX_GT_YY], U[IDX_GT_YZ], U[IDX_GT_ZZ]};
    double gt_inv[6]; inv_sym(gt_inv, gt);
    
    double gt_inv_dd = (d == 0) ? gt_inv[0] : (d == 1) ? gt_inv[3] : gt_inv[5];
    double v_light = alpha * sqrt(fabs(exp(-4.0*phi) * gt_inv_dd) + EPS);
    
    speeds[0] = beta_d - v_light;
    speeds[1] = beta_d + v_light;
    if (speeds[0] > 0.0) speeds[0] = 0.0;
    if (speeds[1] < 0.0) speeds[1] = 0.0;
}

/* ============================================================================
 * LOCAL LAX-FRIEDRICHS FLUX - Ref: [2] Eq. (35)
 * ============================================================================*/
void riemann_llf(double *flux_num, const double *U_L, const double *U_R, const Vec3 *n, const SimConfig *config) {
    double flux_L[N_TOTAL_VARS], flux_R[N_TOTAL_VARS];
    double fx[N_TOTAL_VARS], fy[N_TOTAL_VARS], fz[N_TOTAL_VARS];
    
    z4c_compute_flux(fx, U_L, 0); z4c_compute_flux(fy, U_L, 1); z4c_compute_flux(fz, U_L, 2);
    for (int i = 0; i < N_TOTAL_VARS; i++) flux_L[i] = fx[i]*n->x + fy[i]*n->y + fz[i]*n->z;
    
    z4c_compute_flux(fx, U_R, 0); z4c_compute_flux(fy, U_R, 1); z4c_compute_flux(fz, U_R, 2);
    for (int i = 0; i < N_TOTAL_VARS; i++) flux_R[i] = fx[i]*n->x + fy[i]*n->y + fz[i]*n->z;
    
    double lambda_max = 0.0;
    for (int d = 0; d < 3; d++) {
        double speeds_L[2], speeds_R[2];
        z4c_characteristic_speeds(speeds_L, U_L, d);
        z4c_characteristic_speeds(speeds_R, U_R, d);
        double n_d = (d == 0) ? n->x : (d == 1) ? n->y : n->z;
        double s_max = fmax(fmax(fabs(speeds_L[0]), fabs(speeds_L[1])),
                          fmax(fabs(speeds_R[0]), fabs(speeds_R[1])));
        lambda_max += s_max * fabs(n_d);
    }
    if (lambda_max < EPS) lambda_max = 1.0;
    
    for (int i = 0; i < N_TOTAL_VARS; i++)
        flux_num[i] = 0.5*(flux_L[i] + flux_R[i]) - 0.5*lambda_max*(U_R[i] - U_L[i]);
}

/* ============================================================================
 * BOUNDARY CONDITIONS - Ref: [16] Ruiz et al.
 * ============================================================================*/
void boundary_outflow(double *U_ghost, const double *U_int, const Vec3 *n) {
    for (int i = 0; i < N_TOTAL_VARS; i++) U_ghost[i] = U_int[i];
}

void enforce_constraints(Simulation *sim) {
    for (int c = 0; c < sim->mesh.n_cells; c++) {
        double *U = sim->mesh.cells[c].U;
        
        /* Enforce det(tilde_gamma) = 1 - Ref: [3] Eq. (3.5) */
        double gt[6] = {U[IDX_GT_XX], U[IDX_GT_XY], U[IDX_GT_XZ], U[IDX_GT_YY], U[IDX_GT_YZ], U[IDX_GT_ZZ]};
        double det = det_sym(gt);
        if (fabs(det - 1.0) > 1e-10) {
            double scale = pow(det, -1.0/3.0);
            for (int s = 0; s < 6; s++) U[IDX_GT_XX + s] *= scale;
            U[IDX_PHI] += log(det) / 12.0;
        }
        
        /* Enforce tr(tilde_A) = 0 - Ref: [3] Eq. (3.6) */
        double gt_inv[6]; inv_sym(gt_inv, gt);
        double At[6] = {U[IDX_AT_XX], U[IDX_AT_XY], U[IDX_AT_XZ], U[IDX_AT_YY], U[IDX_AT_YZ], U[IDX_AT_ZZ]};
        double trace = gt_inv[0]*At[0] + gt_inv[3]*At[3] + gt_inv[5]*At[5]
                     + 2.0*(gt_inv[1]*At[1] + gt_inv[2]*At[2] + gt_inv[4]*At[4]);
        if (fabs(trace) > 1e-14)
            for (int s = 0; s < 6; s++) U[IDX_AT_XX + s] -= (trace/3.0) * gt[s];
        
        /* Ensure alpha > 0 */
        if (U[IDX_ALPHA] < 1e-6) U[IDX_ALPHA] = 1e-6;
    }
}

/* ============================================================================
 * COMPUTE RHS - Ref: [2] Dumbser et al., Section 3
 * ============================================================================*/
void z4c_compute_rhs(Simulation *sim) {
    VoronoiMesh *mesh = &sim->mesh;
    
    for (int c = 0; c < mesh->n_cells; c++) {
        VoronoiCell *cell = &mesh->cells[c];
        double *rhs = &sim->rhs[c * N_TOTAL_VARS];
        
        /* 1. Source terms */
        z4c_compute_source(rhs, cell->U, &sim->config.damping, &sim->config.gauge);
        for (int i = 0; i < N_TOTAL_VARS; i++) rhs[i] *= cell->volume;
        
        /* 2. Surface flux integral */
        for (int f = 0; f < cell->n_faces; f++) {
            VoronoiFace *face = &cell->faces[f];
            double U_ext[N_TOTAL_VARS];
            
            if (face->neighbor_id >= 0)
                memcpy(U_ext, mesh->cells[face->neighbor_id].U, N_TOTAL_VARS * sizeof(double));
            else
                boundary_outflow(U_ext, cell->U, &face->normal);
            
            double flux_num[N_TOTAL_VARS];
            riemann_llf(flux_num, cell->U, U_ext, &face->normal, &sim->config);
            
            for (int i = 0; i < N_TOTAL_VARS; i++)
                rhs[i] -= flux_num[i] * face->area;
        }
        
        /* 3. Apply inverse mass matrix (1/volume for FVM) */
        for (int i = 0; i < N_TOTAL_VARS; i++)
            rhs[i] /= cell->volume;
    }
}

/* ============================================================================
 * TIME STEPPING - Ref: [12] Shu & Osher SSP-RK3
 * ============================================================================*/
double compute_dt(Simulation *sim) {
    double lambda_max = 0.0;
    for (int c = 0; c < sim->mesh.n_cells; c++) {
        for (int d = 0; d < 3; d++) {
            double speeds[2];
            z4c_characteristic_speeds(speeds, sim->mesh.cells[c].U, d);
            double s = fmax(fabs(speeds[0]), fabs(speeds[1]));
            if (s > lambda_max) lambda_max = s;
        }
    }
    if (lambda_max < EPS) lambda_max = 1.0;
    return sim->config.cfl * sim->mesh.min_cell_size / lambda_max;
}

void time_step_rk3_ssp(Simulation *sim) {
    int n = sim->mesh.n_cells * N_TOTAL_VARS;
    double *U_n = (double *)malloc(n * sizeof(double));
    
    /* Save U^n */
    for (int c = 0; c < sim->mesh.n_cells; c++)
        memcpy(&U_n[c*N_TOTAL_VARS], sim->mesh.cells[c].U, N_TOTAL_VARS * sizeof(double));
    
    /* Stage 1: U^{(1)} = U^n + dt * L(U^n) */
    z4c_compute_rhs(sim);
    for (int c = 0; c < sim->mesh.n_cells; c++)
        for (int i = 0; i < N_TOTAL_VARS; i++)
            sim->mesh.cells[c].U[i] += sim->dt * sim->rhs[c*N_TOTAL_VARS + i];
    enforce_constraints(sim);
    
    /* Stage 2: U^{(2)} = 3/4 U^n + 1/4 U^{(1)} + 1/4 dt L(U^{(1)}) */
    z4c_compute_rhs(sim);
    for (int c = 0; c < sim->mesh.n_cells; c++)
        for (int i = 0; i < N_TOTAL_VARS; i++)
            sim->mesh.cells[c].U[i] = 0.75*U_n[c*N_TOTAL_VARS+i] + 0.25*sim->mesh.cells[c].U[i]
                                    + 0.25*sim->dt*sim->rhs[c*N_TOTAL_VARS+i];
    enforce_constraints(sim);
    
    /* Stage 3: U^{n+1} = 1/3 U^n + 2/3 U^{(2)} + 2/3 dt L(U^{(2)}) */
    z4c_compute_rhs(sim);
    for (int c = 0; c < sim->mesh.n_cells; c++)
        for (int i = 0; i < N_TOTAL_VARS; i++)
            sim->mesh.cells[c].U[i] = (1.0/3.0)*U_n[c*N_TOTAL_VARS+i] + (2.0/3.0)*sim->mesh.cells[c].U[i]
                                    + (2.0/3.0)*sim->dt*sim->rhs[c*N_TOTAL_VARS+i];
    enforce_constraints(sim);
    
    free(U_n);
    sim->time += sim->dt;
    sim->step++;
}

/* ============================================================================
 * INITIAL DATA - Ref: [14] Alcubierre et al. for gauge wave
 * ============================================================================*/
void initial_data_flat(Simulation *sim) {
    for (int c = 0; c < sim->mesh.n_cells; c++) {
        double *U = sim->mesh.cells[c].U;
        for (int i = 0; i < N_TOTAL_VARS; i++) U[i] = 0.0;
        U[IDX_PHI] = 0.0;
        U[IDX_GT_XX] = U[IDX_GT_YY] = U[IDX_GT_ZZ] = 1.0;
        U[IDX_ALPHA] = 1.0;
    }
    printf("Initialized: Flat spacetime\n");
}

void initial_data_gauge_wave(Simulation *sim, double A) {
    double d = 1.0;
    for (int c = 0; c < sim->mesh.n_cells; c++) {
        double *U = sim->mesh.cells[c].U;
        double x = sim->mesh.cells[c].centroid.x;
        
        for (int i = 0; i < N_TOTAL_VARS; i++) U[i] = 0.0;
        
        /* H(x) = 1 - A*sin(2*pi*x/d) - Ref: [14] Eq. (4.2) */
        double H = 1.0 - A * sin(2.0 * M_PI * x / d);
        double dH = -A * (2.0 * M_PI / d) * cos(2.0 * M_PI * x / d);
        if (H <= 0.0) H = 0.01;
        
        U[IDX_PHI] = log(H) / 12.0;
        U[IDX_GT_XX] = pow(H, 2.0/3.0);
        U[IDX_GT_YY] = U[IDX_GT_ZZ] = pow(H, -1.0/3.0);
        
        double alpha = sqrt(H);
        double K_xx = -dH / (2.0 * alpha);
        double K = K_xx / H;
        U[IDX_K] = K;
        
        double exp_m4phi = 1.0 / pow(H, 1.0/3.0);
        U[IDX_AT_XX] = exp_m4phi * (2.0/3.0) * K_xx;
        U[IDX_AT_YY] = U[IDX_AT_ZZ] = -exp_m4phi * (1.0/3.0) * K;
        
        U[IDX_ALPHA] = alpha;
        U[IDX_AUX_A_X] = 0.5 * dH / H;
        U[IDX_AUX_P_X] = dH / (12.0 * H);
    }
    printf("Initialized: Gauge wave (A=%.4f)\n", A);
}

/* ============================================================================
 * CONSTRAINT MONITORING - Ref: [1] Eqs. (2.125)-(2.126)
 * ============================================================================*/
double compute_hamiltonian(const double *U) {
    double K = U[IDX_K];
    double At[6] = {U[IDX_AT_XX], U[IDX_AT_XY], U[IDX_AT_XZ], U[IDX_AT_YY], U[IDX_AT_YZ], U[IDX_AT_ZZ]};
    double gt[6] = {U[IDX_GT_XX], U[IDX_GT_XY], U[IDX_GT_XZ], U[IDX_GT_YY], U[IDX_GT_YZ], U[IDX_GT_ZZ]};
    double gt_inv[6]; inv_sym(gt_inv, gt);
    double gt_inv_full[3][3], At_full[3][3];
    sym_to_full(gt_inv_full, gt_inv);
    sym_to_full(At_full, At);
    
    double At_sq = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double At_up = 0.0;
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    At_up += gt_inv_full[i][k] * gt_inv_full[j][l] * At_full[k][l];
            At_sq += At_full[i][j] * At_up;
        }
    return (2.0/3.0)*K*K - At_sq;
}

void compute_constraints(Simulation *sim) {
    double H_sum = 0.0, Theta_sum = 0.0, vol = 0.0;
    for (int c = 0; c < sim->mesh.n_cells; c++) {
        double *U = sim->mesh.cells[c].U;
        double v = sim->mesh.cells[c].volume;
        double H = compute_hamiltonian(U);
        H_sum += H*H*v;
        Theta_sum += U[IDX_THETA]*U[IDX_THETA]*v;
        vol += v;
    }
    sim->hamiltonian_L2 = sqrt(H_sum / vol);
    sim->z4_constraint_L2 = sqrt(Theta_sum / vol);
}

/* ============================================================================
 * OUTPUT - VTK format
 * ============================================================================*/
void output_vtk(Simulation *sim, const char *filename) {
    char fname[512]; snprintf(fname, sizeof(fname), "%s.vtk", filename);
    FILE *fp = fopen(fname, "w");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", fname); return; }
    
    fprintf(fp, "# vtk DataFile Version 3.0\nZ4c t=%.6e\nASCII\nDATASET UNSTRUCTURED_GRID\n", sim->time);
    fprintf(fp, "POINTS %d double\n", sim->mesh.n_cells);
    for (int c = 0; c < sim->mesh.n_cells; c++)
        fprintf(fp, "%.10e %.10e %.10e\n", sim->mesh.cells[c].centroid.x,
                sim->mesh.cells[c].centroid.y, sim->mesh.cells[c].centroid.z);
    
    fprintf(fp, "\nCELLS %d %d\n", sim->mesh.n_cells, 2*sim->mesh.n_cells);
    for (int c = 0; c < sim->mesh.n_cells; c++) fprintf(fp, "1 %d\n", c);
    fprintf(fp, "\nCELL_TYPES %d\n", sim->mesh.n_cells);
    for (int c = 0; c < sim->mesh.n_cells; c++) fprintf(fp, "1\n");
    
    fprintf(fp, "\nCELL_DATA %d\n", sim->mesh.n_cells);
    fprintf(fp, "SCALARS alpha double 1\nLOOKUP_TABLE default\n");
    for (int c = 0; c < sim->mesh.n_cells; c++) fprintf(fp, "%.10e\n", sim->mesh.cells[c].U[IDX_ALPHA]);
    fprintf(fp, "\nSCALARS phi double 1\nLOOKUP_TABLE default\n");
    for (int c = 0; c < sim->mesh.n_cells; c++) fprintf(fp, "%.10e\n", sim->mesh.cells[c].U[IDX_PHI]);
    fprintf(fp, "\nSCALARS K double 1\nLOOKUP_TABLE default\n");
    for (int c = 0; c < sim->mesh.n_cells; c++) fprintf(fp, "%.10e\n", sim->mesh.cells[c].U[IDX_K]);
    fprintf(fp, "\nSCALARS Theta double 1\nLOOKUP_TABLE default\n");
    for (int c = 0; c < sim->mesh.n_cells; c++) fprintf(fp, "%.10e\n", sim->mesh.cells[c].U[IDX_THETA]);
    
    fclose(fp);
    printf("Wrote: %s\n", fname);
}

/* ============================================================================
 * SIMULATION DRIVER
 * ============================================================================*/
void config_defaults(SimConfig *config) {
    config->domain_min = (Vec3){-0.5, -0.5, -0.5};
    config->domain_max = (Vec3){0.5, 0.5, 0.5};
    config->nx = config->ny = config->nz = 10;
    config->cfl = 0.4;
    config->t_final = 1.0;
    config->dt_output = 0.1;
    config->damping = (Z4cDampingParams){0.02, 0.0, 0.5, 2.0};
    config->gauge.lapse_type = GAUGE_1_PLUS_LOG;
    config->gauge.shift_type = SHIFT_GAMMA_DRIVER;
    config->gauge.eta_beta = 2.0;
    strcpy(config->output_dir, "./output");
}

int simulation_init(Simulation *sim, SimConfig *config) {
    sim->config = *config;
    if (voronoi_mesh_from_cartesian(&sim->mesh, config->domain_min, config->domain_max,
                                    config->nx, config->ny, config->nz) != 0) return -1;
    gauss_legendre_init(&sim->gauss, N_GAUSS);
    sim->rhs = (double *)calloc(sim->mesh.n_cells * N_TOTAL_VARS, sizeof(double));
    sim->time = 0.0; sim->step = 0;
    return 0;
}

void simulation_free(Simulation *sim) {
    voronoi_mesh_free(&sim->mesh);
    free(sim->rhs);
}

int simulation_run(Simulation *sim) {
    char cmd[512]; snprintf(cmd, sizeof(cmd), "mkdir -p %s", sim->config.output_dir);
    system(cmd);
    
    char fname[512];
    snprintf(fname, sizeof(fname), "%s/z4c_%05d", sim->config.output_dir, 0);
    output_vtk(sim, fname);
    
    compute_constraints(sim);
    printf("t=%.4e: ||H||=%.4e, ||Theta||=%.4e\n", sim->time, sim->hamiltonian_L2, sim->z4_constraint_L2);
    
    double next_output = sim->config.dt_output;
    int out_count = 1;
    
    while (sim->time < sim->config.t_final) {
        sim->dt = compute_dt(sim);
        if (sim->time + sim->dt > sim->config.t_final) sim->dt = sim->config.t_final - sim->time;
        if (sim->time + sim->dt > next_output) sim->dt = next_output - sim->time;
        
        time_step_rk3_ssp(sim);
        
        if (sim->step % 100 == 0)
            printf("Step %d: t=%.4e, dt=%.4e\n", sim->step, sim->time, sim->dt);
        
        if (sim->time >= next_output - 1e-14) {
            compute_constraints(sim);
            printf("Output %d: t=%.4e, ||H||=%.4e, ||Theta||=%.4e\n",
                   out_count, sim->time, sim->hamiltonian_L2, sim->z4_constraint_L2);
            snprintf(fname, sizeof(fname), "%s/z4c_%05d", sim->config.output_dir, out_count);
            output_vtk(sim, fname);
            next_output += sim->config.dt_output;
            out_count++;
        }
        
        /* Check for NaN */
        for (int c = 0; c < sim->mesh.n_cells; c++)
            for (int v = 0; v < N_TOTAL_VARS; v++)
                if (isnan(sim->mesh.cells[c].U[v])) {
                    fprintf(stderr, "NaN at cell %d var %d\n", c, v);
                    return -1;
                }
    }
    
    printf("\nSimulation complete: t=%.4e, steps=%d\n", sim->time, sim->step);
    return 0;
}

/* ============================================================================
 * MAIN
 * ============================================================================*/
int main(int argc, char *argv[]) {
    printf("╔════════════════════════════════════════════════════════════╗\n");
    printf("║     Z4c DG-FVM Numerical Relativity Code                   ║\n");
    printf("║     Refs: [1] Baumgarte & Shapiro (2010)                   ║\n");
    printf("║           [2] Dumbser et al., JCP 348 (2017)               ║\n");
    printf("║           [3] Hilditch et al., PRD 88, 084057 (2013)       ║\n");
    printf("╚════════════════════════════════════════════════════════════╝\n\n");
    
    SimConfig config;
    config_defaults(&config);
    
    char test[64] = "gauge_wave";
    double amplitude = 0.01;
    
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-t") == 0 && i+1 < argc) strcpy(test, argv[++i]);
        else if (strcmp(argv[i], "-n") == 0 && i+1 < argc) config.nx = config.ny = config.nz = atoi(argv[++i]);
        else if (strcmp(argv[i], "-T") == 0 && i+1 < argc) config.t_final = atof(argv[++i]);
        else if (strcmp(argv[i], "-A") == 0 && i+1 < argc) amplitude = atof(argv[++i]);
        else if (strcmp(argv[i], "-o") == 0 && i+1 < argc) strcpy(config.output_dir, argv[++i]);
        else if (strcmp(argv[i], "-h") == 0) {
            printf("Usage: %s [-t test] [-n cells] [-T tfinal] [-A amp] [-o dir]\n", argv[0]);
            printf("  -t: flat, gauge_wave (default)\n");
            return 0;
        }
    }
    
    printf("Test: %s, Grid: %dx%dx%d, T_final: %.2f\n\n", test, config.nx, config.ny, config.nz, config.t_final);
    
    Simulation sim;
    if (simulation_init(&sim, &config) != 0) return 1;
    
    if (strcmp(test, "flat") == 0) initial_data_flat(&sim);
    else initial_data_gauge_wave(&sim, amplitude);
    
    clock_t start = clock();
    int result = simulation_run(&sim);
    double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
    
    printf("\nWall time: %.2f seconds\n", elapsed);
    printf("Result: %s\n", result == 0 ? "SUCCESS" : "FAILED");
    
    simulation_free(&sim);
    return result;
}
