/**
 * @file z4c_dg.h
 * @brief Discontinuous Galerkin Finite Volume Method for Z4c Numerical Relativity
 * 
 * Key References:
 * [1] Baumgarte & Shapiro, "Numerical Relativity", Cambridge (2010)
 * [2] Dumbser et al., J. Comput. Phys. 348, 70-117 (2017) - ADER-DG for GR
 * [3] Hilditch et al., Phys. Rev. D 88, 084057 (2013) - Z4c formulation
 * [4] Bernuzzi & Hilditch, Phys. Rev. D 81, 084003 (2010) - Constraint damping
 */

#ifndef Z4C_DG_H
#define Z4C_DG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

/* ============ CONFIGURATION ============ */
#define DIM 3
#define MAX_POLY_DEGREE 4
#define N_GAUSS 5
#define MAX_FACES 24
#define MAX_NEIGHBORS 24
#define MAX_DOF_PER_CELL 125
#define EPS 1.0e-14
#define C_LIGHT 1.0

/* ============ Z4c VARIABLES ============ */
/* Primary: phi, tilde_gamma_ij(6), K, tilde_A_ij(6), tilde_Gamma^i(3), Theta, alpha, beta^i(3) = 22 */
/* Auxiliary: A_i(3), B^i_j(9), D_kij(18), P_i(3) = 33 */
#define N_Z4C_VARS 22
#define N_AUX_VARS 33
#define N_TOTAL_VARS (N_Z4C_VARS + N_AUX_VARS)

/* Variable indices - Reference: [3] Hilditch et al. */
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

/* Auxiliary indices - Reference: [2] Dumbser et al., Eqs. (15)-(18) */
#define IDX_AUX_A_X 22
#define IDX_AUX_A_Y 23
#define IDX_AUX_A_Z 24
#define IDX_AUX_B_XX 25
#define IDX_AUX_B_XY 26
#define IDX_AUX_B_XZ 27
#define IDX_AUX_B_YX 28
#define IDX_AUX_B_YY 29
#define IDX_AUX_B_YZ 30
#define IDX_AUX_B_ZX 31
#define IDX_AUX_B_ZY 32
#define IDX_AUX_B_ZZ 33
#define IDX_AUX_D_BASE 34
#define IDX_AUX_P_X 52
#define IDX_AUX_P_Y 53
#define IDX_AUX_P_Z 54

/* ============ STRUCTURES ============ */

typedef struct { double x, y, z; } Vec3;

/* Constraint damping - Reference: [4] Bernuzzi & Hilditch, Eq. (2.7) */
typedef struct {
    double kappa1, kappa2, kappa3, eta;
} Z4cDampingParams;

/* Gauge conditions - Reference: [1] Chapter 4.2 */
typedef enum { GAUGE_GEODESIC, GAUGE_HARMONIC, GAUGE_1_PLUS_LOG, GAUGE_SHOCK_AVOIDING } LapseGauge;
typedef enum { SHIFT_ZERO, SHIFT_GAMMA_DRIVER, SHIFT_HARMONIC } ShiftGauge;

typedef struct {
    LapseGauge lapse_type;
    ShiftGauge shift_type;
    double mu_L, mu_S, eta_beta;
} GaugeParams;

typedef struct {
    int neighbor_id;
    Vec3 centroid, normal;
    double area;
    int n_vertices;
    Vec3 *vertices;
} VoronoiFace;

typedef struct {
    int id;
    Vec3 generator, centroid;
    double volume;
    int n_faces;
    VoronoiFace faces[MAX_FACES];
    int n_neighbors;
    int neighbors[MAX_NEIGHBORS];
    int n_dof;
    double *U, *U_modal;
} VoronoiCell;

typedef struct {
    int n_cells, n_boundary_cells;
    VoronoiCell *cells;
    Vec3 bbox_min, bbox_max;
    double min_cell_size, max_cell_size;
} VoronoiMesh;

typedef struct {
    int n_points;
    double points[N_GAUSS], weights[N_GAUSS];
} GaussQuadrature;

typedef struct {
    int poly_degree, n_basis_1d, n_basis_3d;
    double **phi_at_gauss, **dphi_at_gauss;
    double *mass_matrix, *inv_mass_matrix;
    double *stiffness_x, *stiffness_y, *stiffness_z;
} DGBasis;

typedef enum { RIEMANN_LLF, RIEMANN_UPWIND, RIEMANN_CENTRAL } RiemannSolverType;
typedef enum { TIME_EULER, TIME_RK2, TIME_RK3_SSP, TIME_RK4 } TimeIntegrator;

typedef struct {
    double domain_size[DIM];
    Vec3 domain_min, domain_max;
    int n_cells_target, poly_degree;
    TimeIntegrator time_integrator;
    double cfl, t_final, dt_output;
    Z4cDampingParams damping;
    GaugeParams gauge;
    RiemannSolverType riemann_solver;
    bool use_limiting;
    double limiter_threshold;
    char output_dir[256];
    int output_every_n_steps;
} SimConfig;

typedef struct {
    SimConfig config;
    VoronoiMesh mesh;
    DGBasis basis;
    GaussQuadrature gauss;
    double time;
    int step;
    double dt;
    double *rhs, *flux_buffer, *U_temp;
    double hamiltonian_L2, momentum_L2[DIM], z4_constraint_L2;
} Simulation;

/* ============ FUNCTION DECLARATIONS ============ */

/* Initialization */
int simulation_init(Simulation *sim, SimConfig *config);
void simulation_free(Simulation *sim);
void config_set_defaults(SimConfig *config);
int simulation_run(Simulation *sim);

/* Mesh */
int voronoi_mesh_create(VoronoiMesh *mesh, Vec3 domain_min, Vec3 domain_max, int n_cells_target);
void voronoi_mesh_free(VoronoiMesh *mesh);
int voronoi_mesh_from_cartesian(VoronoiMesh *mesh, Vec3 domain_min, Vec3 domain_max, int nx, int ny, int nz);

/* DG Basis */
int dg_basis_init(DGBasis *basis, int poly_degree);
void dg_basis_free(DGBasis *basis);
double legendre_P(int n, double x);
double legendre_dP(int n, double x);
void gauss_legendre_init(GaussQuadrature *gauss, int n_points);

/* Z4c Evolution - Reference: [3] Hilditch et al., Eqs. (13)-(23) */
void z4c_compute_rhs(Simulation *sim);
void z4c_compute_source(double *source, const double *U, const Z4cDampingParams *damp, const GaugeParams *gauge);
void z4c_compute_flux(double *flux, const double *U, int direction);
void z4c_characteristic_speeds(double *speeds, const double *U, int direction);

/* Riemann Solvers - Reference: [2] Section 3 */
void riemann_llf(double *flux_num, const double *U_left, const double *U_right, const Vec3 *normal, const SimConfig *config);
void compute_numerical_flux(double *flux_num, const double *U_left, const double *U_right, const Vec3 *normal, const SimConfig *config);
void compute_face_flux_integral(double *flux_integral, const VoronoiCell *cell, const VoronoiMesh *mesh, const DGBasis *basis, const GaussQuadrature *gauss, const SimConfig *config);

/* Time Integration - Reference: [12] Shu & Osher SSP-RK methods */
void time_step(Simulation *sim);
void time_step_euler(Simulation *sim);
void time_step_rk3_ssp(Simulation *sim);
void time_step_rk4(Simulation *sim);
double compute_dt(Simulation *sim);

/* Initial Data - Reference: [14] Alcubierre et al. gauge wave, [1] Section 2.4 Schwarzschild */
void initial_data_flat_spacetime(Simulation *sim);
void initial_data_gauge_wave(Simulation *sim, double amplitude);
void initial_data_schwarzschild(Simulation *sim, double mass);
void initial_data_robust_stability(Simulation *sim, double amplitude);
void initial_data_linear_wave(Simulation *sim, double amplitude);

/* Boundary Conditions - Reference: [16] Ruiz et al. */
void apply_boundary_conditions(Simulation *sim);
void boundary_outflow(double *U_ghost, const double *U_interior, const Vec3 *normal);
void boundary_reflecting(double *U_ghost, const double *U_interior, const Vec3 *normal);
void enforce_conformal_metric_constraint(Simulation *sim);
void enforce_traceless_constraint(Simulation *sim);
void enforce_lapse_positivity(Simulation *sim);

/* Constraints - Reference: [1] Eqs. (2.125)-(2.126) */
void compute_constraints(Simulation *sim);
double compute_hamiltonian_constraint(const double *U);
void compute_momentum_constraint(double *momentum, const double *U);
void print_constraint_summary(const Simulation *sim);

/* Output */
void output_vtk(Simulation *sim, const char *filename);
void output_1d_slice(Simulation *sim, const char *filename, double y_slice, double z_slice);
void output_time_series(Simulation *sim, const char *filename);
void output_checkpoint(Simulation *sim, const char *filename);
int load_checkpoint(Simulation *sim, const char *filename);

/* Utilities */
double vec3_dot(const Vec3 *a, const Vec3 *b);
Vec3 vec3_cross(const Vec3 *a, const Vec3 *b);
double vec3_norm(const Vec3 *a);
Vec3 vec3_normalize(const Vec3 *a);
Vec3 vec3_add(const Vec3 *a, const Vec3 *b);
Vec3 vec3_sub(const Vec3 *a, const Vec3 *b);
Vec3 vec3_scale(const Vec3 *a, double s);
int sym_idx(int i, int j);
void sym_to_full(double full[DIM][DIM], const double sym[6]);
double det_sym(const double sym[6]);
void inv_sym(double inv[6], const double sym[6]);

/* Tests */
int test_flat_spacetime(int n_cells, int n_steps);
int test_gauge_wave(double amplitude, int n_cells, double t_final);

#endif /* Z4C_DG_H */
