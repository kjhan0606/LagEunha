/* ================================================================
 *  exam_gpu.cu
 *
 *  CUDA kernels and GPU memory management for Voronoi force loop.
 *  Phase 2 implementation: av_mode=0 (Monaghan AV) force kernel.
 *  Phases 3-4 will add av_mode=1 (NS stress) and av_mode=5 (HLLC).
 *
 *  Compiled with nvcc.
 * ================================================================ */

#ifdef USE_CUDA

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <cuda_runtime.h>

/* CUB for device-wide reduction */
#include <cub/cub.cuh>
#include "gfs_pair.h"
#include "gfs_velocity_limiter.h"
#include "gfs_thermal64.h"

extern "C" {
#include "exam_gpu.h"
}

static __device__ unsigned long long d_gfs_pair_hits;
static __device__ unsigned long long d_gfs_candidate_fallback_cells;
static __device__ unsigned long long d_gfs_invalid_raw_volume_cells;
static __device__ unsigned long long d_gfs_invalid_velocity_fit_cells;
static __device__ unsigned long long d_gfs_tess_geometry_errors;
static __device__ unsigned long long d_gfs_stencil_expanded_cells;
static __device__ unsigned long long d_gfs_stencil_halo_clipped_cells;
static long long h_gfs_pair_hits = 0;
static __device__ unsigned long long d_varpress_geom_max_key;
static double h_varpress_geom_max_rate=0;
static long long h_varpress_geom_max_id=-1;
extern "C" int gpu_varpress_w0_enabled(void) {
#ifdef GFS_VARPRESS_W0
    return 1;
#else
    return 0;
#endif
}
extern "C" double gpu_varpress_geom_max_rate(void) { return h_varpress_geom_max_rate; }
extern "C" long long gpu_varpress_geom_max_id(void) { return h_varpress_geom_max_id; }
static int h_geometry_trial=0,h_geometry_failed=0;
extern "C" void gpu_geometry_trial_begin(void) { h_geometry_trial=1; h_geometry_failed=0; }
extern "C" void gpu_geometry_trial_end(void) { h_geometry_trial=0; h_geometry_failed=0; }
extern "C" int gpu_geometry_trial_active(void) { return h_geometry_trial; }
extern "C" int gpu_geometry_trial_failed(void) { return h_geometry_failed; }
extern "C" void gpu_geometry_trial_mark_failed(void) { h_geometry_failed=1; }
/* Geometry from the latest host tessellation launch, for optional dumps. */
static int h_diag_tess_mx = 0, h_diag_tess_my = 0;
static double h_diag_tess_cellsize = 0.0;
static double h_diag_tess_xmin = 0.0, h_diag_tess_ymin = 0.0;
static double h_diag_global_w2_max = 0.0;

/* ================================================================
 *  Error checking macro
 * ================================================================ */
#define CUDA_CHECK(call) do { \
    cudaError_t err = (call); \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(1); \
    } \
} while(0)

/* ================================================================
 *  Device inline: get2dUpqradRk4 (Laguerre face velocity)
 *
 *  Returns (v_face - v_i): the face velocity relative to particle i.
 *  Accounts for Laguerre weight evolution (w2, w2old) and
 *  anti-aliasing correction.
 * ================================================================ */
__device__ __forceinline__
void dev_get2dUpqradRk4(
    /* particle i */
    double xi, double yi, double vxi, double vyi,
    float w2i, float w2oldi, float csi,
    /* particle j */
    double xj, double yj, double vxj, double vyj,
    float w2j, float w2oldj, float csj,
    /* params */
    double dtold, int w2lin,
    /* output: v_face - v_i */
    double *uradx, double *urady)
{
    double upqx = vxj - vxi;
    double upqy = vyj - vyi;
    double erx = xj - xi;
    double ery = yj - yi;
    double dpq2 = erx * erx + ery * ery;
    double dpq = sqrt(dpq2);
    double dpq_inv = 1.0 / dpq;
    erx *= dpq_inv;
    ery *= dpq_inv;

    double er_dot_upq = erx * upqx + ery * upqy;

    double wp2 = (double)w2i,  wq2 = (double)w2j;
    double wpold2 = (double)w2oldi, wqold2 = (double)w2oldj;

    double fact1 = 0.5 * (1.0 + (wp2 - wq2) / dpq2);

    double fact2;
    if (w2lin) {
        /* GFS_W2_LINEAR: w2old carries Wdot (mirrors gfs_upqrad_rk4). */
        fact2 = 0.5 * (wpold2 - wqold2) * dpq_inv;
    } else {
    double dwpdt = (sqrt(wp2) - sqrt(wpold2)) / dtold;
    double dwqdt = (sqrt(wq2) - sqrt(wqold2)) / dtold;

    /* Clamp dw/dt by sound speed */
    double vpw, vqw;
    if (dwpdt > 0) vpw = fmin((double)csi, dwpdt);
    else           vpw = fmax(-(double)csi, dwpdt);
    if (dwqdt > 0) vqw = fmin((double)csj, dwqdt);
    else           vqw = fmax(-(double)csj, dwqdt);

    fact2 = (sqrt(wp2) * vpw - sqrt(wq2) * vqw) * dpq_inv;
    }
    fact2 -= (wp2 - wq2) / dpq2 * er_dot_upq;

    *uradx = fact1 * upqx + fact2 * erx;
    *urady = fact1 * upqy + fact2 * ery;
}

/* ================================================================
 *  Device inline: HLLC solver
 * ================================================================ */
/* Fix C v2: extreme P-ratio guard with pmin floor discriminator (mirrors CPU) */
#ifndef HLLC_VACUUM_PRATIO
#define HLLC_VACUUM_PRATIO 100.0
#endif
#ifndef HLLC_VACUUM_PMIN_GUARD
#define HLLC_VACUUM_PMIN_GUARD 0.01
#endif

/* Acceleration CFL: dt ≤ ACC_CFL_FRAC · sqrt(h_i / |a_i|).
 * Mirrors CPU exam.c — caps free-fall time across cell under current force. */
#ifndef ACC_CFL_FRAC
#define ACC_CFL_FRAC 0.25
#endif

/* dK-aware CFL: dt ≤ DKCFL_FRAC · K / |dK| (mirrors CPU exam.c).
 * Mode 2 mitigation — preserves conservation, no dK clipping. */
#ifndef DKCFL_FRAC
#define DKCFL_FRAC 0.1
#endif
#ifndef DKCFL_KFLOOR
#define DKCFL_KFLOOR 1.0e-10
#endif

__device__ __forceinline__
void dev_hll_star_state(
    double rhoL, double pL, double vnL, double cL,
    double rhoR, double pR, double vnR, double cR,
    double Gamma,
    double *pstar, double *vnstar)
{
    /* HLL average of the conserved state, then the EOS pressure of that average.
     * Same state as hll_star_state in exam.c. */
    double SL = fmin(vnL - cL, vnR - cR);
    double SR = fmax(vnL + cL, vnR + cR);
    /* No upwind selection here. This is the HLL average state between the
     * two waves, used as the face pressure and velocity of a moving
     * (Lagrangian) face, not an Eulerian flux at x/t = 0. The old
     * "SL >= 0 -> left state, SR <= 0 -> right state" test is made in
     * the frame of the caller, which is the lab frame for the extreme
     * face: in an orbiting disk (|v| >> c) it returned the hot cell's own
     * P and v_n, so a cold neighbour expanded against the hot pressure and
     * went to ie < 0 (Kepler floor runaway, see code_review_grokbot.md).
     * The average state is Galilean covariant (P* invariant, v* shifts
     * with the frame), so dropping the test gives the frame-independent
     * answer, and it is identical wherever SL < 0 < SR already held. */
    double gm1 = Gamma - 1.0;
    if (gm1 < 1.0e-8) gm1 = 1.0e-8;
    double rhoLs = fmax(rhoL, 1.0e-30);
    double rhoRs = fmax(rhoR, 1.0e-30);
    double EL = pL / (gm1 * rhoLs) + 0.5 * vnL * vnL;
    double ER = pR / (gm1 * rhoRs) + 0.5 * vnR * vnR;
    double EdenL = rhoL * EL, EdenR = rhoR * ER;
    double denom = SR - SL;
    if (fabs(denom) < 1.0e-30) {
        *pstar = 0.5 * (pL + pR);
        *vnstar = 0.5 * (vnL + vnR);
        return;
    }
    double FrL = rhoL * vnL, FrR = rhoR * vnR;
    double FmL = rhoL * vnL * vnL + pL, FmR = rhoR * vnR * vnR + pR;
    double FeL = vnL * (EdenL + pL), FeR = vnR * (EdenR + pR);
    double rhoS = (SR * rhoR - SL * rhoL - (FrR - FrL)) / denom;
    double momS = (SR * (rhoR * vnR) - SL * (rhoL * vnL) - (FmR - FmL)) / denom;
    double ES   = (SR * EdenR - SL * EdenL - (FeR - FeL)) / denom;
    if (!(rhoS > 1.0e-30) || isnan(rhoS) || isinf(rhoS)
            || isnan(momS) || isinf(momS) || isnan(ES) || isinf(ES)) {
        *pstar = 0.5 * (pL + pR);
        *vnstar = 0.5 * (vnL + vnR);
        return;
    }
    double vnS = momS / rhoS;
    double eint = ES / rhoS - 0.5 * vnS * vnS;
    double PS = gm1 * rhoS * eint;
    if (!(PS > 0) || isnan(PS) || isinf(PS) || isnan(vnS) || isinf(vnS)) {
        *pstar = 0.5 * (pL + pR);
        *vnstar = 0.5 * (vnL + vnR);
        return;
    }
    *pstar = PS;
    *vnstar = vnS;
}

__device__ __forceinline__
void dev_hllc_face_2d(
    double rhoL, double pL, double vnL, double cL,
    double rhoR, double pR, double vnR, double cR,
    double Gamma, int phase1,
    double *pstar, double *vnstar)
{
    /* True HLLC: ports GIZMO get_wavespeeds_and_pressure_star (reimann.h:527).
     * 3-stage: Gaburov/Davis → Roe fallback → PVRS-Rusanov last resort. */
    const double tiny = 1.0e-30;
    double S_L, S_R, S_M, P_M;
    double cmax = cL > cR ? cL : cR;

    /* Phase 1: HLL average only on a large pressure ratio. A uniformly
     * cold face stays on HLLC. Runs before Fix C. */
    if (phase1) {
        double pmin = pL < pR ? pL : pR;
        double pmax = pL > pR ? pL : pR;
        if (pmax > 100.0 * fmax(pmin, 1.0e-30)) {
            dev_hll_star_state(rhoL, pL, vnL, cL, rhoR, pR, vnR, cR,
                               Gamma, pstar, vnstar);
            return;
        }
    }

    /* Fix C v2: extreme P-ratio guard (HLLE fallback), only when both sides above floor */
    {
        double pmin = pL < pR ? pL : pR;
        double pmax = pL > pR ? pL : pR;
        if (pmin > HLLC_VACUUM_PMIN_GUARD && pmax > HLLC_VACUUM_PRATIO * pmin) {
            *pstar  = pmin;
            *vnstar = 0.5*(vnL + vnR);
            return;
        }
    }

    if ((vnR - vnL) > cmax) {
        *pstar = tiny; *vnstar = 0.5*(vnL + vnR);
        return;
    }

    /* Primary: Gaburov/Davis */
    {
        double vmin = vnL < vnR ? vnL : vnR;
        double vmax = vnL > vnR ? vnL : vnR;
        S_L = vmin - cmax;
        S_R = vmax + cmax;
        double rho_wt_L = rhoL * (S_L - vnL);
        double rho_wt_R = rhoR * (S_R - vnR);
        double denom = rho_wt_L - rho_wt_R;
        if (fabs(denom) > tiny) {
            S_M = ((pR - pL) + rho_wt_L*vnL - rho_wt_R*vnR) / denom;
            P_M = (pL*rho_wt_R - pR*rho_wt_L + rho_wt_L*rho_wt_R*(vnR - vnL))
                  / (rho_wt_R - rho_wt_L);
            if (P_M > 0 && !isnan(P_M) && !isinf(P_M)) {
                *pstar = P_M; *vnstar = S_M; return;
            }
        }
    }

    /* Roe fallback */
    {
        double sqL = sqrt(rhoL), sqR = sqrt(rhoR);
        double sq_inv = 1.0 / (sqL + sqR + tiny);
        double vn_roe = (sqL*vnL + sqR*vnR) * sq_inv;
        double h_L = Gamma*pL/((Gamma - 1.0)*rhoL) + 0.5*vnL*vnL;
        double h_R = Gamma*pR/((Gamma - 1.0)*rhoR) + 0.5*vnR*vnR;
        double h_roe = (sqL*h_L + sqR*h_R) * sq_inv;
        double c2_roe = (Gamma - 1.0) * (h_roe - 0.5*vn_roe*vn_roe);
        if (c2_roe < tiny) c2_roe = tiny;
        double c_roe = sqrt(c2_roe);

        double srA = vnR + cR, srB = vn_roe + c_roe;
        double slA = vnL - cL, slB = vn_roe - c_roe;
        S_R = srA > srB ? srA : srB;
        S_L = slA < slB ? slA : slB;

        double rho_wt_R =  rhoR * (S_R - vnR);
        double rho_wt_L = -rhoL * (S_L - vnL);
        double denom = rho_wt_R + rho_wt_L;
        if (fabs(denom) > tiny) {
            S_M = (rho_wt_R*vnR + rho_wt_L*vnL + (pL - pR)) / denom;
            P_M = rhoL*(vnL - S_L)*(vnL - S_M) + pL;
            if (P_M > 0 && !isnan(P_M) && !isinf(P_M)) {
                *pstar = P_M; *vnstar = S_M; return;
            }
        }
    }

    /* Last resort: PVRS + Rusanov */
    {
        P_M = 0.5*((pL + pR) + (vnL - vnR)*0.25*(rhoL + rhoR)*(cL + cR));
        S_M = 0.5*(vnR + vnL) + 2.0*(pL - pR)/((rhoL + rhoR)*(cL + cR) + tiny);
        double a1 = fabs(vnL - cL), a2 = fabs(vnR - cR);
        double a3 = fabs(vnL + cL), a4 = fabs(vnR + cR);
        double S_plus = a1;
        if (a2 > S_plus) S_plus = a2;
        if (a3 > S_plus) S_plus = a3;
        if (a4 > S_plus) S_plus = a4;
        if (S_M < -S_plus) S_M = -S_plus;
        if (S_M >  S_plus) S_M =  S_plus;
        if (P_M < tiny) P_M = tiny;
        *pstar = P_M; *vnstar = S_M;
    }
}

__device__ __forceinline__
void dev_hllc_face_2d_rest_frame(
    double rhoL, double pL, double vnL_lab, double cL,
    double rhoR, double pR, double vnR_lab, double cR,
    double wn, double Gamma, int phase1,
    double *pstar, double *vnstar_lab)
{
    double vnL = vnL_lab - wn;
    double vnR = vnR_lab - wn;
    double pst, vnst;
    dev_hllc_face_2d(rhoL, pL, vnL, cL, rhoR, pR, vnR, cR, Gamma, phase1, &pst, &vnst);
    *pstar      = pst;
    *vnstar_lab = vnst + wn;
}

/* ================================================================
 *  Main force kernel: 1 thread per real particle
 *
 *  Each thread loops over its faces [face_offset[i], face_offset[i+1]).
 *  Computes: ax, ay, die, dt_min, vsig_max.
 *
 *  Supports all av_modes (0, 1, 2, 3, 5) via if-branching.
 *  Since av_mode is uniform across all threads, no warp divergence.
 * ================================================================ */
__global__ __launch_bounds__(256, 2)
void getAccVoro2DBlend_kernel(
    /* Face CSR */
    const int *__restrict__ face_offset,
    const double *__restrict__ c1x, const double *__restrict__ c1y,
    const double *__restrict__ c2x, const double *__restrict__ c2y,
    const int *__restrict__ neighbor_idx,
    const int *__restrict__ kp_idx, const int *__restrict__ km_idx,
    const int *__restrict__ is_ghost,
    /* Particle SoA — read */
    const double *__restrict__ px, const double *__restrict__ py,
    const double *__restrict__ pvx, const double *__restrict__ pvy,
    const double *__restrict__ pmass,
    const float *__restrict__ pden, const float *__restrict__ ppressure,
    const float *__restrict__ pcsound, const float *__restrict__ pvolume,
    const float *__restrict__ pie_in, const float *__restrict__ pw2,
    const float *__restrict__ pw2old,
    /* Stress fields (may be NULL if av_mode==0 && nu_phys<=0) */
    const double *__restrict__ pgUxx, const double *__restrict__ pgUxy,
    const double *__restrict__ pgUyx, const double *__restrict__ pgUyy,
    const double *__restrict__ pdPdx, const double *__restrict__ pdPdy,
    const double *__restrict__ pdRhodx, const double *__restrict__ pdRhody,
    const double *__restrict__ ptauxx, const double *__restrict__ ptauxy,
    const double *__restrict__ ptauyy,
    const double *__restrict__ palpha_cd, const double *__restrict__ pdivv,
    /* Hyperviscosity Pass 2 input: per-particle velocity Laplacian */
    const double *__restrict__ plap_vx, const double *__restrict__ plap_vy,
    /* Output arrays */
    double *__restrict__ ax_out, double *__restrict__ ay_out,
    float *__restrict__ die_out, float *__restrict__ dK_out,
    float *__restrict__ dt_out,
    double *__restrict__ vsig_max_out,
    double *__restrict__ vsmoothx_out, double *__restrict__ vsmoothy_out,
    /* Optional per-face audit output, NULL in ordinary runs. */
    double *__restrict__ audit_work,
    double *__restrict__ audit_pressure,
    double *__restrict__ audit_udotds,
    /* Physics parameters */
    int n_particles,
    int av_mode, int use_muscl,
    double OoA, double Courant, double Gamma,
    double alphavis, double betavis, double etavis, double epsvis,
    double nu_phys, double prandtl,
    double cd_amax, double blend_theta, double dtold,
    double hyperv_alpha, double hyperv_force_cap,
    int phase1, int laguerre_rot, int face_charge, int geom_fv, int w2lin,
    double var_geom_cfl, int var_geom_audit,
    int thermal64, const double *__restrict__ pU64,
    const double *__restrict__ pV64, double *__restrict__ out_Udot64)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_particles) return;

    /* Load owner particle to registers */
    double ibp_x  = px[i],  ibp_y  = py[i];
    double ibp_vx = pvx[i], ibp_vy = pvy[i];
    double ibp_mass = pmass[i];
    double ibp_den      = (double)pden[i];
    double ibp_pressure = (double)ppressure[i];
    double ibp_csound   = (double)pcsound[i];
    double ibp_volume   = (double)pvolume[i];
    if (thermal64) {
        ibp_volume = pV64[i];
        ibp_den = ibp_mass / ibp_volume;
        ibp_pressure = (Gamma - 1.0) * pU64[i] / ibp_volume;
        ibp_csound = sqrt(Gamma * ibp_pressure / ibp_den);
    }
    float  ibp_w2     = pw2[i];
    float  ibp_w2old  = pw2old[i];

    double fx = 0, fy = 0, die = 0;
    /* Dissipative-only contributions (viscous heating + heat conduction).
     * Mirrors CPU dK_diss (exam.c:3829..4459); converted to dK after the loop. */
    double dK_diss = 0;
    float  my_dt = 1.0e10f;
    double my_vsig_max = 0;
    /* XSPH: face-area-weighted velocity-difference accumulator
     * (real-neighbor faces only; ghost mirrors excluded). */
    double xsph_vx_acc = 0, xsph_vy_acc = 0, xsph_w_acc = 0;
    /* Hyperviscosity Pass 2: A_face/d_ij × (lap_v_j - lap_v_i) accumulator */
    double biharm_vx_acc = 0, biharm_vy_acc = 0;
    double ibp_lap_vx = (hyperv_alpha > 0) ? plap_vx[i] : 0.0;
    double ibp_lap_vy = (hyperv_alpha > 0) ? plap_vy[i] : 0.0;
    double dte = 0;
    double thermal_rate64 = 0;
#ifdef GFS_VARPRESS_W0
    double var_sweep_abs=0;
#endif

    int f_begin = face_offset[i];
    int f_end   = face_offset[i + 1];

    for (int f = f_begin; f < f_end; f++) {
        if (audit_work) {
            audit_work[f] = 0.0;
            audit_pressure[f] = 0.0;
            audit_udotds[f] = 0.0;
        }
        int j = neighbor_idx[f];
        if (j < 0) continue;

        int jbp_is_ghost = is_ghost[f];
#ifdef GFS_VARPRESS_W0
        /* This experimental variational branch is defined only for a
         * smooth, unweighted, interior material cell. Fail visibly if an
         * unsupported face is presented instead of silently evolving it. */
        if (jbp_is_ghost || ibp_w2 != 0.0f || pw2[j] != 0.0f) {
            ax_out[i] = NAN; ay_out[i] = NAN; die_out[i] = NAN;
            dt_out[i] = NAN; return;
        }
#endif

        /* Load neighbor particle */
        double jbp_x  = px[j],  jbp_y  = py[j];
        double jbp_vx = pvx[j], jbp_vy = pvy[j];
        double jbp_den      = (double)pden[j];
        double jbp_pressure = (double)ppressure[j];
        double jbp_csound   = (double)pcsound[j];
        if (thermal64) {
            jbp_den = pmass[j] / pV64[j];
            jbp_pressure = (Gamma - 1.0) * pU64[j] / pV64[j];
            jbp_csound = sqrt(Gamma * jbp_pressure / jbp_den);
        }
        float  jbp_w2     = pw2[j];
        float  jbp_w2old  = pw2old[j];

        /* Face geometry from corners (local/relative coords) */
        double line_x = c2x[f] - c1x[f];
        double line_y = c2y[f] - c1y[f];

        /* dS = outward-pointing area vector = (line_y, -line_x) */
        /* voro2D_norm rotates line 90° CCW: (y,-x), normalized, then
           scaled by facearea. Since the cross product check always
           gives c.z > 0 for these line vectors, dS = (line_y, -line_x). */
        double dSx = line_y;
        double dSy = -line_x;
        double facearea = sqrt(dSx * dSx + dSy * dSy);

        /* XSPH accumulation: real (non-ghost) neighbors only.
         * Wall mirror velocities are reflected and would pull the
         * particle off-axis. */
        if (!jbp_is_ghost) {
            xsph_vx_acc += facearea * (jbp_vx - ibp_vx);
            xsph_vy_acc += facearea * (jbp_vy - ibp_vy);
            xsph_w_acc  += facearea;
        }

        /* Hyperviscosity Pass 2: accumulate A_face/d_ij × (lap_v_j - lap_v_i)
         * for real neighbors only. */
        if (hyperv_alpha > 0 && !jbp_is_ghost) {
            double drx_h = jbp_x - ibp_x;
            double dry_h = jbp_y - ibp_y;
            double dramp_h = sqrt(drx_h * drx_h + dry_h * dry_h);
            if (dramp_h > 0) {
                double invd_bh = facearea / dramp_h;
                biharm_vx_acc += invd_bh * (plap_vx[j] - ibp_lap_vx);
                biharm_vy_acc += invd_bh * (plap_vy[j] - ibp_lap_vy);
            }
        }

        /* Inter-particle direction */
        double drx = jbp_x - ibp_x;
        double dry = jbp_y - ibp_y;
        double dramp = sqrt(drx * drx + dry * dry);
        double dramp_inv = 1.0 / (dramp + 1.0e-30);
        double erx = drx * dramp_inv;
        double ery = dry * dramp_inv;

        /* Face velocity (always needed for energy accumulation) */
        double uradx, urady;
        dev_get2dUpqradRk4(
            ibp_x, ibp_y, ibp_vx, ibp_vy, ibp_w2, ibp_w2old, (float)ibp_csound,
            jbp_x, jbp_y, jbp_vx, jbp_vy, jbp_w2, jbp_w2old, (float)jbp_csound,
            dtold, w2lin, &uradx, &urady);
        /* Centroid correction for a pure Voronoi face (w = 0). */
        if (ibp_w2 == 0.0f && jbp_w2 == 0.0f && dramp > 0.0) {
            double cx = 0.5 * (c1x[f] + c2x[f]);
            double cy = 0.5 * (c1y[f] + c2y[f]);
            double dxc = cx - 0.5 * drx;
            double dyc = cy - 0.5 * dry;
            double delta = -((jbp_vx - ibp_vx) * dxc + (jbp_vy - ibp_vy) * dyc) / dramp;
            uradx += delta * erx;
            urady += delta * ery;
        } else if (laguerre_rot && dramp > 0.0) {
            /* GFS_LAGUERRE_ROTATION: Laguerre face, reference point at the
             * anchor fact1 (x_j - x_i) that dev_get2dUpqradRk4 moves with. */
            double fact1 = 0.5 * (1.0 + ((double)ibp_w2 - (double)jbp_w2) / (drx * drx + dry * dry));
            double cx = 0.5 * (c1x[f] + c2x[f]);
            double cy = 0.5 * (c1y[f] + c2y[f]);
            double dxc = cx - fact1 * drx;
            double dyc = cy - fact1 * dry;
            double delta = -((jbp_vx - ibp_vx) * dxc + (jbp_vy - ibp_vy) * dyc) / dramp;
            uradx += delta * erx;
            urady += delta * ery;
        }

        double pi_total;
        double tau_dot_dS_x = 0, tau_dot_dS_y = 0;
        /* Fix A+B: dissipative AV pressure extracted from pi_total.
         * Sum of (Monaghan AV + CD10 + HLLC-jump-excess) >= 0 when active.
         * Used downstream to compute pair-symmetric eps_AV >= 0 heating. */
        double p_av_for_heat = 0;
        /* PdV speed. Ordinary faces keep the geometric Voronoi speed.
         * Extreme faces (and the non-MUSCL cell HLL) use the Riemann normal speed. */
        int riemann_vstar = 0;
        double riemann_vn = 0, riemann_nx = 0, riemann_ny = 0;
        int phase1_cell = 0;
        double phase1_vn = 0, phase1_nx = 0, phase1_ny = 0;

        /* ======== av_mode dispatch ======== */
        if (av_mode == 0 && nu_phys <= 0) {
            /* --- Monaghan AV only (no NS stress) --- */
            /* M(n,m) pressure */
            int kp = kp_idx[f];
            int km = km_idx[f];
            double p_kp = (kp >= 0) ? (double)ppressure[kp] : ibp_pressure;
            double p_km = (km >= 0) ? (double)ppressure[km] : jbp_pressure;
            double w = OoA / 3.0;
            pi_total = (0.5 - w) * (ibp_pressure + jbp_pressure)
                     + w * (p_kp + p_km);

            /* Monaghan AV */
            double uijx = jbp_vx - ibp_vx;
            double uijy = jbp_vy - ibp_vy;
            double rvel = erx * uijx + ery * uijy;
            if (rvel < 0) {
                double wcomp = sqrt((double)ibp_w2) + sqrt((double)jbp_w2);
                double scaleFactor = (wcomp > etavis ? wcomp : etavis);
                double drampScale = dramp / scaleFactor;
                double mu = rvel /* Voronoi-AV-reform: no length regularization */;
                double meanden = 0.5 * (ibp_den + jbp_den);
                double meanCsound = 0.5 * (ibp_csound + jbp_csound);
                pi_total += (-alphavis * meanCsound * mu
                             + betavis * mu * mu) * meanden;
            }

        } else if (jbp_is_ghost) {
            /* --- Ghost face: M(n,m) + Monaghan AV, no NS/HLLC --- */
            int kp = kp_idx[f];
            int km = km_idx[f];
            double p_kp = (kp >= 0) ? (double)ppressure[kp] : ibp_pressure;
            double p_km = (km >= 0) ? (double)ppressure[km] : jbp_pressure;
            double w = OoA / 3.0;
            pi_total = (0.5 - w) * (ibp_pressure + jbp_pressure)
                     + w * (p_kp + p_km);

            double uijx = jbp_vx - ibp_vx;
            double uijy = jbp_vy - ibp_vy;
            double rvel = erx * uijx + ery * uijy;
            if (rvel < 0) {
                double wcomp = sqrt((double)ibp_w2) + sqrt((double)jbp_w2);
                double scaleFactor = (wcomp > etavis ? wcomp : etavis);
                double drampScale = dramp / scaleFactor;
                double mu = rvel /* Voronoi-AV-reform: no length regularization */;
                double meanden = 0.5 * (ibp_den + jbp_den);
                double meanCsound = 0.5 * (ibp_csound + jbp_csound);
                pi_total += (-alphavis * meanCsound * mu
                             + betavis * mu * mu) * meanden;
            }

        } else if (av_mode == 5) {
            /* --- Arepo-Laguerre: face rest frame HLLC + MUSCL --- */
            double ds_mag_inv = 1.0 / (facearea + 1.0e-30);
            double nx_hat = dSx * ds_mag_inv;
            double ny_hat = dSy * ds_mag_inv;

            /* Reuse pre-computed face velocity */
            double wx = ibp_vx + uradx;
            double wy = ibp_vy + urady;
            double wn = wx * nx_hat + wy * ny_hat;

            double vnL_lab = ibp_vx * nx_hat + ibp_vy * ny_hat;
            double vnR_lab = jbp_vx * nx_hat + jbp_vy * ny_hat;
            double pL = ibp_pressure;
            double pR = jbp_pressure;
            double rhoL = ibp_den;
            double rhoR = jbp_den;

            if (use_muscl) {
                double xf_rel_i = 0.5 * (c1x[f] + c2x[f]);
                double yf_rel_i = 0.5 * (c1y[f] + c2y[f]);
                double dx_ij = jbp_x - ibp_x;
                double dy_ij = jbp_y - ibp_y;
                double dx_iF = xf_rel_i;
                double dy_iF = yf_rel_i;
                double dx_jF = xf_rel_i - dx_ij;
                double dy_jF = yf_rel_i - dy_ij;

                double drho_i = pdRhodx[i] * dx_iF + pdRhody[i] * dy_iF;
                double drho_j = pdRhodx[j] * dx_jF + pdRhody[j] * dy_jF;
                double dp_i = pdPdx[i] * dx_iF + pdPdy[i] * dy_iF;
                double dp_j = pdPdx[j] * dx_jF + pdPdy[j] * dy_jF;
                double dvxL = pgUxx[i] * dx_iF + pgUxy[i] * dy_iF;
                double dvyL = pgUyx[i] * dx_iF + pgUyy[i] * dy_iF;
                double dvnL0 = nx_hat * dvxL + ny_hat * dvyL;
                double dvxR = pgUxx[j] * dx_jF + pgUxy[j] * dy_jF;
                double dvyR = pgUyx[j] * dx_jF + pgUyy[j] * dy_jF;
                double dvnR0 = nx_hat * dvxR + ny_hat * dvyR;

                rhoL += drho_i;
                rhoR += drho_j;
                pL += dp_i;
                pR += dp_j;
                vnL_lab += dvnL0;
                vnR_lab += dvnR0;
            }

            if (pL < 1.0e-10) pL = 1.0e-10;
            if (pR < 1.0e-10) pR = 1.0e-10;
            if (rhoL < 1.0e-10) rhoL = 1.0e-10;
            if (rhoR < 1.0e-10) rhoR = 1.0e-10;

            double pst, vnst_lab;
            dev_hllc_face_2d_rest_frame(rhoL, pL, vnL_lab, ibp_csound,
                                        rhoR, pR, vnR_lab, jbp_csound,
                                        wn, Gamma, phase1, &pst, &vnst_lab);
            pi_total = pst;
            /* Same HLL state from both ends: cell-centered values along i→j.
             * MUSCL left/right states are not swaps of each other, so the
             * shock speed was applied on one end of the face only. */
            if (phase1 && dramp > 0) {
                double pmin_c = fmin(ibp_pressure, jbp_pressure);
                double pmax_c = fmax(ibp_pressure, jbp_pressure);
                if (pmax_c > 100.0 * fmax(pmin_c, 1.0e-30)) {
                    double ps, vns;
                    dev_hll_star_state(
                        ibp_den, ibp_pressure, ibp_vx * erx + ibp_vy * ery, ibp_csound,
                        jbp_den, jbp_pressure, jbp_vx * erx + jbp_vy * ery, jbp_csound,
                        Gamma, &ps, &vns);
                    pi_total = ps;
                    riemann_vstar = 1;
                    riemann_vn = vns;
                    riemann_nx = erx;
                    riemann_ny = ery;
                }
            }

            /* Optional Monaghan AV for contact noise */
            if (alphavis > 0) {
                double uijx = jbp_vx - ibp_vx;
                double uijy = jbp_vy - ibp_vy;
                double rvel = erx * uijx + ery * uijy;
                if (rvel < 0) {
                    double wcomp = sqrt((double)ibp_w2) + sqrt((double)jbp_w2);
                    double scaleFactor = (wcomp > etavis ? wcomp : etavis);
                    double drampScale = dramp / scaleFactor;
                    double mu = rvel /* Voronoi-AV-reform: no length regularization */;
                    double meanden = 0.5 * (ibp_den + jbp_den);
                    double meanCsound = 0.5 * (ibp_csound + jbp_csound);
                    pi_total += (-alphavis * meanCsound * mu
                                 + betavis * mu * mu) * meanden;
                }
            }

        } else if (av_mode == 3) {
            /* --- DEPRECATED: pure HLLC lab frame --- */
            double ds_mag_inv = 1.0 / (facearea + 1.0e-30);
            double nx_hat = dSx * ds_mag_inv;
            double ny_hat = dSy * ds_mag_inv;
            double vnL = ibp_vx * nx_hat + ibp_vy * ny_hat;
            double vnR = jbp_vx * nx_hat + jbp_vy * ny_hat;
            double pL = ibp_pressure, pR = jbp_pressure;

            if (use_muscl) {
                double xf_rel_i = 0.5 * (c1x[f] + c2x[f]);
                double yf_rel_i = 0.5 * (c1y[f] + c2y[f]);
                double dx_ij = jbp_x - ibp_x;
                double dy_ij = jbp_y - ibp_y;
                double dx_iF = xf_rel_i, dy_iF = yf_rel_i;
                double dx_jF = xf_rel_i - dx_ij, dy_jF = yf_rel_i - dy_ij;

                double dp_i = pdPdx[i] * dx_iF + pdPdy[i] * dy_iF;
                double dp_j = pdPdx[j] * dx_jF + pdPdy[j] * dy_jF;
                pL += dp_i; pR += dp_j;

                double dvnL0 = (pgUxx[i] * nx_hat + pgUxy[i] * ny_hat) * dx_iF
                             + (pgUyx[i] * nx_hat + pgUyy[i] * ny_hat) * dy_iF;
                double dvnR0 = (pgUxx[j] * nx_hat + pgUxy[j] * ny_hat) * dx_jF
                             + (pgUyx[j] * nx_hat + pgUyy[j] * ny_hat) * dy_jF;
                vnL += dvnL0; vnR += dvnR0;

                double pmin = fmin(ibp_pressure, jbp_pressure);
                double pmax = fmax(ibp_pressure, jbp_pressure);
                pL = fmax(pmin, fmin(pmax, pL));
                pR = fmax(pmin, fmin(pmax, pR));
                double vnmin = fmin(vnL - dvnL0, vnR - dvnR0);
                double vnmax = fmax(vnL - dvnL0, vnR - dvnR0);
                vnL = fmax(vnmin, fmin(vnmax, vnL));
                vnR = fmax(vnmin, fmin(vnmax, vnR));
            }
            if (pL < 1.0e-10) pL = 1.0e-10;
            if (pR < 1.0e-10) pR = 1.0e-10;

            double pst, vnst;
            dev_hllc_face_2d(ibp_den, pL, vnL, ibp_csound,
                             jbp_den, pR, vnR, jbp_csound,
                             Gamma, phase1, &pst, &vnst);
            pi_total = pst;

            if (cd_amax > 0) {
                double dvn_hllc = vnR - vnL;
                if (dvn_hllc < 0) {
                    double alpha_face = 0.5 * (palpha_cd[i] + palpha_cd[j]);
                    double rho_mean = 0.5 * (ibp_den + jbp_den);
                    double vsig_cd = ibp_csound + jbp_csound - fmin(0.0, dvn_hllc);
                    pi_total += 0.5 * alpha_face * vsig_cd * rho_mean * (-dvn_hllc);
                }
            }
            if (alphavis > 0) {
                double uijx = jbp_vx - ibp_vx;
                double uijy = jbp_vy - ibp_vy;
                double rvel = erx * uijx + ery * uijy;
                if (rvel < 0) {
                    double wcomp = sqrt((double)ibp_w2) + sqrt((double)jbp_w2);
                    double scaleFactor = (wcomp > etavis ? wcomp : etavis);
                    double drampScale = dramp / scaleFactor;
                    double mu = rvel /* Voronoi-AV-reform: no length regularization */;
                    double meanden = 0.5 * (ibp_den + jbp_den);
                    double meanCsound = 0.5 * (ibp_csound + jbp_csound);
                    pi_total += (-alphavis * meanCsound * mu
                                 + betavis * mu * mu) * meanden;
                }
            }

        } else {
            /* --- av_mode == 0 (with nu_phys>0), 1, or 2: NS stress path --- */
            /* M(n,m) pressure */
            int kp = kp_idx[f];
            int km = km_idx[f];
            double p_kp = (kp >= 0) ? (double)ppressure[kp] : ibp_pressure;
            double p_km = (km >= 0) ? (double)ppressure[km] : jbp_pressure;
            double w = OoA / 3.0;
            double p_mnm = (0.5 - w) * (ibp_pressure + jbp_pressure)
                         + w * (p_kp + p_km);

            /* NS stress traction τ·dS */
            if (av_mode == 1 || (av_mode == 0 && nu_phys > 0)) {
                double txx_face, txy_face, tyy_face;
                if (OoA > 0 && kp >= 0 && km >= 0) {
                    /* M(n,m) interpolation of stress */
                    double t_kp, t_km;
                    t_kp = ptauxx[kp]; t_km = ptauxx[km];
                    txx_face = (0.5 - w) * (ptauxx[i] + ptauxx[j]) + w * (t_kp + t_km);
                    t_kp = ptauxy[kp]; t_km = ptauxy[km];
                    txy_face = (0.5 - w) * (ptauxy[i] + ptauxy[j]) + w * (t_kp + t_km);
                    t_kp = ptauyy[kp]; t_km = ptauyy[km];
                    tyy_face = (0.5 - w) * (ptauyy[i] + ptauyy[j]) + w * (t_kp + t_km);
                } else {
                    txx_face = 0.5 * (ptauxx[i] + ptauxx[j]);
                    txy_face = 0.5 * (ptauxy[i] + ptauxy[j]);
                    tyy_face = 0.5 * (ptauyy[i] + ptauyy[j]);
                }
                /* τ_code = -τ_NS, negate for correct sign */
                tau_dot_dS_x = -(txx_face * dSx + txy_face * dSy);
                tau_dot_dS_y = -(txy_face * dSx + tyy_face * dSy);
            }

            if (av_mode == 0 || av_mode == 1) {
                /* M(n,m) + NS stress + optional Monaghan AV (Balsara-modulated) */
                pi_total = p_mnm;
                if (alphavis > 0) {
                    double uijx = jbp_vx - ibp_vx;
                    double uijy = jbp_vy - ibp_vy;
                    double rvel = erx * uijx + ery * uijy;
                    if (rvel < 0) {
                        double wcomp = sqrt((double)ibp_w2) + sqrt((double)jbp_w2);
                        double scaleFactor = (wcomp > etavis ? wcomp : etavis);
                        double drampScale = dramp / scaleFactor;
                        double mu = rvel /* Voronoi-AV-reform: no length regularization */;
                        double meanden = 0.5 * (ibp_den + jbp_den);
                        double meanCsound = 0.5 * (ibp_csound + jbp_csound);
                        double divv_i  = pdivv[i];
                        double curlv_i = pgUyx[i] - pgUxy[i];
                        double divv_j  = pdivv[j];
                        double curlv_j = pgUyx[j] - pgUxy[j];
                        double xi_b    = 1.0e-4 * meanCsound / dramp;
                        double f_i = fabs(divv_i) / (fabs(divv_i) + fabs(curlv_i) + xi_b);
                        double f_j = fabs(divv_j) / (fabs(divv_j) + fabs(curlv_j) + xi_b);
                        double f_balsara = 0.5 * (f_i + f_j);
                        double mon_av_pi = f_balsara * (-alphavis * meanCsound * mu
                                                         + betavis * mu * mu) * meanden;
                        pi_total += mon_av_pi;
                        if (mon_av_pi > 0) p_av_for_heat += mon_av_pi;
                    }
                }
                /* CD10 viscous pressure (Cullen-Dehnen 2010) — strong-shock
                 * detector that survives Voronoi cell collapse. */
                if (av_mode == 1 && cd_amax > 0) {
                    double uijx_cd = jbp_vx - ibp_vx;
                    double uijy_cd = jbp_vy - ibp_vy;
                    double rvel_cd = erx * uijx_cd + ery * uijy_cd;
                    if (rvel_cd < 0) {
                        double alpha_face = 0.5 * (palpha_cd[i] + palpha_cd[j]);
                        double rho_mean   = 0.5 * (ibp_den + jbp_den);
                        double vsig_cd    = ibp_csound + jbp_csound - rvel_cd;
                        double cd10_pi = 0.5 * alpha_face * vsig_cd * rho_mean * (-rvel_cd);
                        pi_total += cd10_pi;
                        p_av_for_heat += cd10_pi;
                    }
                }
            } else {
                /* av_mode == 2: two-tier blend */
                double alpha_max_pq = fmax(palpha_cd[i], palpha_cd[j]);
                double theta_pq = 0;
                if (pdivv[i] < -1.0e-10 || pdivv[j] < -1.0e-10)
                    theta_pq = blend_theta;
                double f_pq = fmin(1.0, alpha_max_pq / (cd_amax + 1.0e-30) + theta_pq);

                double p_hllc = p_mnm;
                double Pi_cd10 = 0;

                if (f_pq > 1.0e-6) {
                    double ds_mag_inv = 1.0 / (facearea + 1.0e-30);
                    double nx_hat = dSx * ds_mag_inv;
                    double ny_hat = dSy * ds_mag_inv;
                    double vnL = ibp_vx * nx_hat + ibp_vy * ny_hat;
                    double vnR = jbp_vx * nx_hat + jbp_vy * ny_hat;

                    /* Contact detection */
                    double rho_mean_det = 0.5 * (ibp_den + jbp_den);
                    double p_mean_det = 0.5 * (ibp_pressure + jbp_pressure);
                    double drho_rel = fabs(ibp_den - jbp_den) / (rho_mean_det + 1.0e-30);
                    double dp_rel = fabs(ibp_pressure - jbp_pressure) / (p_mean_det + 1.0e-30);
                    int is_pure_contact = (drho_rel > 0.1) && (dp_rel < 0.02);

                    double pL = ibp_pressure, pR = jbp_pressure;
                    double rhoL = ibp_den, rhoR = jbp_den;
                    double cL = ibp_csound, cR = jbp_csound;

                    if (use_muscl) {
                        double xf_rel_i = 0.5 * (c1x[f] + c2x[f]);
                        double yf_rel_i = 0.5 * (c1y[f] + c2y[f]);
                        double dx_ij = jbp_x - ibp_x, dy_ij = jbp_y - ibp_y;
                        double dx_iF = xf_rel_i, dy_iF = yf_rel_i;
                        double dx_jF = xf_rel_i - dx_ij, dy_jF = yf_rel_i - dy_ij;

                        double dp_i = pdPdx[i] * dx_iF + pdPdy[i] * dy_iF;
                        double dp_j = pdPdx[j] * dx_jF + pdPdy[j] * dy_jF;
                        pL += dp_i; pR += dp_j;

                        double dvnL0 = 0, dvnR0 = 0;
                        if (!is_pure_contact) {
                            dvnL0 = (pgUxx[i]*nx_hat + pgUxy[i]*ny_hat)*dx_iF
                                   + (pgUyx[i]*nx_hat + pgUyy[i]*ny_hat)*dy_iF;
                            dvnR0 = (pgUxx[j]*nx_hat + pgUxy[j]*ny_hat)*dx_jF
                                   + (pgUyx[j]*nx_hat + pgUyy[j]*ny_hat)*dy_jF;
                            vnL += dvnL0; vnR += dvnR0;
                        }

                        double pmin = fmin(ibp_pressure, jbp_pressure);
                        double pmax = fmax(ibp_pressure, jbp_pressure);
                        pL = fmax(pmin, fmin(pmax, pL));
                        pR = fmax(pmin, fmin(pmax, pR));
                        double vnmin = fmin(vnL - dvnL0, vnR - dvnR0);
                        double vnmax = fmax(vnL - dvnL0, vnR - dvnR0);
                        vnL = fmax(vnmin, fmin(vnmax, vnL));
                        vnR = fmax(vnmin, fmin(vnmax, vnR));
                    }

                    pL = use_muscl ? pL : ibp_pressure;
                    pR = use_muscl ? pR : jbp_pressure;

                    double wx = 0.5 * (ibp_vx + jbp_vx);
                    double wy = 0.5 * (ibp_vy + jbp_vy);
                    double wn = wx * nx_hat + wy * ny_hat;

                    double pst, vnst_lab;
                    dev_hllc_face_2d_rest_frame(rhoL, pL, vnL, cL,
                                                rhoR, pR, vnR, cR,
                                                wn, Gamma, phase1, &pst, &vnst_lab);
                    p_hllc = pst;

                    if (f_pq <= 0.5) {
                        double dvn = vnR - vnL;
                        if (dvn < 0) {
                            double alpha_face = 0.5 * (palpha_cd[i] + palpha_cd[j]);
                            double rho_mean = 0.5 * (rhoL + rhoR);
                            double vsig_cd = cL + cR - fmin(0.0, dvn);
                            Pi_cd10 = 0.5 * alpha_face * vsig_cd * rho_mean * (-dvn);
                        }
                    }
                }

                pi_total = (1.0 - f_pq) * p_mnm + f_pq * (p_hllc + Pi_cd10);
                tau_dot_dS_x *= (1.0 - f_pq);
                tau_dot_dS_y *= (1.0 - f_pq);
                /* HLLC star excess pressure over centered M(n,m) is the
                 * dissipative part of the HLLC face flux. Together with
                 * Pi_cd10 (always non-negative when active) this is the
                 * total dissipative AV pressure for the av_mode=2 blend. */
                double Pi_hllc_extra = (p_hllc > p_mnm ? (p_hllc - p_mnm) : 0);
                p_av_for_heat = f_pq * (Pi_hllc_extra + Pi_cd10);
            }
        }

        /* ======== Accumulate forces and energy rates ======== */
        /* uradx, urady already computed before av_mode dispatch */

        /* Cell-centered HLL when MUSCL is off. Ghost faces stay on the wall path.
         * With MUSCL, av_mode 5 already built P* from the limited face states. */
        if (phase1 && !use_muscl && !jbp_is_ghost && facearea > 0) {
            double invA = 1.0 / facearea;
            double nx = dSx * invA;
            double ny = dSy * invA;
            double ps, vns;
            p_av_for_heat = 0;
            tau_dot_dS_x = 0;
            tau_dot_dS_y = 0;
            dev_hll_star_state(
                ibp_den, ibp_pressure, ibp_vx * nx + ibp_vy * ny, ibp_csound,
                jbp_den, jbp_pressure, jbp_vx * nx + jbp_vy * ny, jbp_csound,
                Gamma, &ps, &vns);
            pi_total = ps;
            phase1_cell = 1;
            phase1_vn = vns;
            phase1_nx = nx;
            phase1_ny = ny;
        }

        /* Local pair restoring pressure. Stays on at zero relative velocity
         * and is part of the face pressure, so the work remains in dE. */
        if (!jbp_is_ghost) {
            double vclose = -((jbp_vx - ibp_vx) * erx + (jbp_vy - ibp_vy) * ery);
            if (vclose < 0.0) vclose = 0.0;
            double qpair = (jbp_vx - ibp_vx) * erx + (jbp_vy - ibp_vy) * ery;
            double p_pair = gfs_pair_pressure(
                dramp, ibp_volume, (double)pvolume[j],
                ibp_den, jbp_den, ibp_csound, jbp_csound, vclose);
            p_pair = gfs_pair_work_limit(
                p_pair, facearea, dtold,
                ibp_mass, (double)pmass[j], qpair,
                (double)pie_in[i], (double)pie_in[j], 6.0);
            /* Charge cap (gfs_pair.h): s_i from the face velocity used for
             * the work below, s_j = qA - s_i. Same as exam.c. */
            if (p_pair > 0.0) {
                double uax, uay;
                if (phase1_cell && !geom_fv) {
                    uax = phase1_vn * phase1_nx; uay = phase1_vn * phase1_ny;
                } else if (riemann_vstar && !geom_fv) {
                    uax = riemann_vn * riemann_nx; uay = riemann_vn * riemann_ny;
                } else {
                    uax = ibp_vx + uradx; uay = ibp_vy + urady;
                }
                double s_i = (uax - ibp_vx) * dSx + (uay - ibp_vy) * dSy;
                double qA  = (jbp_vx - ibp_vx) * dSx + (jbp_vy - ibp_vy) * dSy;
                p_pair = gfs_pair_charge_limit(p_pair, s_i, qA - s_i, dtold,
                    (double)pie_in[i], (double)pie_in[j], 6.0);
            }
            if (p_pair > 0.0) {
                pi_total += p_pair;
                atomicAdd(&d_gfs_pair_hits, 1ull);
            }
        }

        /* GFS_FACE_CHARGE_LIMIT (default off): same cap on the whole face
         * pressure, see exam.c. */
        if (face_charge && !jbp_is_ghost && pi_total > 0.0) {
            double uax, uay;
            if (phase1_cell && !geom_fv) {
                uax = phase1_vn * phase1_nx; uay = phase1_vn * phase1_ny;
            } else if (riemann_vstar && !geom_fv) {
                uax = riemann_vn * riemann_nx; uay = riemann_vn * riemann_ny;
            } else {
                uax = ibp_vx + uradx; uay = ibp_vy + urady;
            }
            double s_i = (uax - ibp_vx) * dSx + (uay - ibp_vy) * dSy;
            double qA  = (jbp_vx - ibp_vx) * dSx + (jbp_vy - ibp_vy) * dSy;
            double pc = gfs_pair_charge_limit(pi_total, s_i, qA - s_i, dtold,
                (double)pie_in[i], (double)pie_in[j], 6.0);
            if (pc < pi_total) {
                double f = pc / pi_total;
                p_av_for_heat *= f;
                pi_total = pc;
            }
        }

        /* Fix A+B: compression gate on AV.
         *   pi_total = p_rev + p_av_for_heat (decomposition)
         *     p_rev          := p_mnm (face-symmetric, reversible)
         *     p_av_for_heat  := Monaghan AV + CD10 + HLLC-jump-excess (dissipative)
         *   FACE GATE: if u_n_compress > 0 (approaching cells):
         *     - F_AV acts; eps_AV_face = p_av * u_n * |dS| (>=0)
         *     - die_AV += 1/2 eps_AV  (pair-symmetric: j gets the other half)
         *   else: strip p_av_for_heat from pi_total (force AND heat both off). */
        double u_n_compress;
        {
            double dvx_av = jbp_vx - ibp_vx;
            double dvy_av = jbp_vy - ibp_vy;
            u_n_compress = -(dvx_av * erx + dvy_av * ery);
        }
        double eps_av_face = 0;
        if (p_av_for_heat > 0 && u_n_compress > 0) {
            eps_av_face = p_av_for_heat * u_n_compress * facearea;
        } else if (p_av_for_heat != 0) {
            pi_total -= p_av_for_heat;
            p_av_for_heat = 0;
        }

        /* Reversible PdV from face-symmetric reference pressure.
         *   p_rev_face = pi_total - p_av_for_heat
         * In av_mode 3/5/ghost: p_av_for_heat=0 -> p_rev_face = pi_total
         *   (preserves OLD behavior for those paths). */
        {
            double p_rev_face = pi_total - p_av_for_heat;
            die += -p_rev_face * (uradx * dSx + urady * dSy);
        }
        /* Pair-symmetric AV heating (1/2 here, 1/2 in j's loop) — >=0 by construction */
        if (eps_av_face > 0) {
            double half_eps = 0.5 * eps_av_face;
            die     += half_eps;
            dK_diss += half_eps;
        }
        /* NS viscous heating: (τ·dS) · (v_face - v_i) — dissipative */
        {
            double die_visc_face = tau_dot_dS_x * uradx + tau_dot_dS_y * urady;
            die     += die_visc_face;
            dK_diss += die_visc_face;
        }

        /* Heat conduction — dissipative.
         * dramp is floored at heff (= cell radius) — same floor used by CFL —
         * to avoid divide-by-near-zero when Voronoi tessellation produces
         * degenerate face geometry (Sedov2D rings IC at Nx=512 hit dramp~1e-26
         * and detonated K to ~1e22 in one step). */
        if (!phase1 && nu_phys > 0 && prandtl > 0 && !jbp_is_ghost
            && ibp_den > 0 && jbp_den > 0
            && ibp_volume > 0 && facearea > 0) {
            double chi = nu_phys / prandtl;
            double Ti = ibp_pressure / ibp_den;
            double Tj = jbp_pressure / jbp_den;
            double rho_face = 0.5 * (ibp_den + jbp_den);
            /* dramp floor at heff (= cell radius) — same floor used by CFL. */
            double dramp_cond = fmax(dramp, 0.25 * sqrt(ibp_volume));
            double die_cond_face = chi * rho_face * (Tj - Ti) / dramp_cond * facearea;
            /* OOB diagnostic: log anomalous conduction faces. Sedov2D Nx>=256
             * explodes K to ~1e23 in step 1 — pinpoint which i,j pair is bad. */
            if (fabs(die_cond_face) > 1.0e6) {
                printf("[GPU_COND_DIAG] i=%d j=%d Pi=%g Pj=%g rhoi=%g rhoj=%g "
                       "dramp=%g dramp_cond=%g facearea=%g die_cond=%g jghost=%d\n",
                       i, j, ibp_pressure, jbp_pressure, ibp_den, jbp_den,
                       dramp, dramp_cond, facearea, die_cond_face, jbp_is_ghost);
            }
            die     += die_cond_face;
            dK_diss += die_cond_face;
        }

#ifdef GFS_VARPRESS_W0
        /* For an unweighted Voronoi pair, w_n = v_i·(n/2+a)
         * +v_j·(n/2-a), a=(x_face-(x_i+x_j)/2)/d. The discrete
         * pressure force B^T P is the adjoint of geometric Vdot=Bv.
         * Pair reciprocity is exact only for identical opposite faces; the
         * retained native 5x5 geometry is checked separately. */
        const double face_cx = 0.5 * (c1x[f] + c2x[f]);
        const double face_cy = 0.5 * (c1y[f] + c2y[f]);
        const double off_x = face_cx - 0.5 * drx;
        const double off_y = face_cy - 0.5 * dry;
        const double b_ix = 0.5*dSx + facearea*off_x*dramp_inv;
        const double b_iy = 0.5*dSy + facearea*off_y*dramp_inv;
        const double dp = ibp_pressure - jbp_pressure;
        const double var_fx = dp*b_ix;
        const double var_fy = dp*b_iy;
        const double var_vdot = uradx*dSx + urady*dSy;
        if (thermal64) thermal_rate64 -= ibp_pressure * var_vdot;
        if (var_geom_cfl>0 || var_geom_audit) var_sweep_abs+=fabs(var_vdot);
#endif

        /* Phase 1 stores dE/dt. Ordinary faces move at the geometric speed
         * v_i + u_rad. Extreme faces move at the lab-frame Riemann normal speed. */
        if (phase1) {
            double uax, uay;
            if (phase1_cell && !geom_fv) {
                uax = phase1_vn * phase1_nx;
                uay = phase1_vn * phase1_ny;
            } else if (riemann_vstar && !geom_fv) {
                uax = riemann_vn * riemann_nx;
                uay = riemann_vn * riemann_ny;
            } else {
                uax = ibp_vx + uradx;
                uay = ibp_vy + urady;
            }
#ifdef GFS_VARPRESS_W0
            const double face_work = ibp_vx*var_fx + ibp_vy*var_fy
                                   - ibp_pressure*var_vdot;
            dte += face_work;
            if (audit_work) {
                audit_work[f] = face_work;
                audit_pressure[f] = ibp_pressure;
                audit_udotds[f] = var_vdot;
            }
#else
            double udotds = uax * dSx + uay * dSy;
            double face_work = -pi_total * udotds
                 + tau_dot_dS_x * uax + tau_dot_dS_y * uay;
            dte += face_work;
            if (audit_work) {
                audit_work[f] = face_work;
                audit_pressure[f] = pi_total;
                audit_udotds[f] = udotds;
            }
#endif
        }

        /* Pressure force consistent with the selected work law. */
#ifdef GFS_VARPRESS_W0
        fx += var_fx;
        fy += var_fy;
#else
        fx += -pi_total * dSx + tau_dot_dS_x;
        fy += -pi_total * dSy + tau_dot_dS_y;
#endif

        /* Signal velocity for CFL */
        double dvx = jbp_vx - ibp_vx;
        double dvy = jbp_vy - ibp_vy;
        double VdotR = dvx * erx + dvy * ery;
        double vsig = jbp_csound + ibp_csound - fmin(0.0, VdotR);

        /* Real faces use the true separation. Flooring dramp up to the
         * cell radius lets a flattened pair step across each other.
         * Ghost mirrors keep the floor. */
        double heff = 0.25 * sqrt(ibp_volume);
        double dramp_cfl = jbp_is_ghost ? fmax(dramp, heff) : dramp;
        float dt_face = (float)(2.0 * Courant * dramp_cfl / vsig);
        if (dt_face < my_dt) my_dt = dt_face;

        /* Track vsig_max for CD10 */
        if (av_mode >= 1) {
            double vsig_cd = jbp_is_ghost ?
                (ibp_csound + jbp_csound) : vsig;
            if (vsig_cd > my_vsig_max) my_vsig_max = vsig_cd;
        }
    } /* end face loop */

#ifdef GFS_VARPRESS_W0
    /* Bound swept volume before face cancellation. This captures rapid
     * bisector rotation as well as compression/expansion, independently
     * of the normal acoustic signal speed. No force or heat is changed. */
    if (var_sweep_abs>0 && ibp_volume>0) {
        double rate=var_sweep_abs/ibp_volume;
        if (var_geom_cfl>0) my_dt=fmin((double)my_dt,var_geom_cfl/rate);
        if (var_geom_audit) {
            unsigned long long key=((unsigned long long)__float_as_uint((float)rate)<<32)
                                  |(unsigned int)i;
            atomicMax(&d_varpress_geom_max_key,key);
        }
    }
#endif

    /* Viscous CFL */
    if (av_mode >= 1 || nu_phys > 0) {
        double h_i = sqrt(ibp_volume);
        double nu_cd_i = (av_mode >= 1) ?
            palpha_cd[i] * h_i * ibp_csound : 0;
        double chi = (nu_phys > 0 && prandtl > 0) ? nu_phys / prandtl : 0;
        double nu_eff = fmax(fmax(nu_phys, nu_cd_i), chi);
        if (nu_eff > 0) {
            float dt_visc = (float)(0.5 * h_i * h_i / nu_eff);
            if (dt_visc < my_dt) my_dt = dt_visc;
        }
    }

    /* Hyperviscosity force: dv/dt -= ν₄ × biharm_v_i / V_i
     * ν₄ = α₄ × cs × Δx³ with Δx = sqrt(V_i) (cell-local scale).
     * Selectively damps high-k modes (rate ∝ ν₄ k⁴).
     * Optional cap: |a_HV| ≤ cap·cs²/Δx (CFL accel) — prevents lap_v runaway
     * from feeding back into a destabilizing force. */
    if (hyperv_alpha > 0 && ibp_volume > 0) {
        double invVol_hv = 1.0 / ibp_volume;
        double dx_eff = sqrt(ibp_volume);
        double nu4 = hyperv_alpha * ibp_csound * dx_eff * dx_eff * dx_eff;
        double a_hv_x = -nu4 * biharm_vx_acc * invVol_hv;
        double a_hv_y = -nu4 * biharm_vy_acc * invVol_hv;
        if (hyperv_force_cap > 0) {
            double a_hv_mag = sqrt(a_hv_x*a_hv_x + a_hv_y*a_hv_y);
            double a_hv_max = hyperv_force_cap * ibp_csound * ibp_csound / dx_eff;
            if (a_hv_mag > a_hv_max && a_hv_mag > 0) {
                double s = a_hv_max / a_hv_mag;
                a_hv_x *= s; a_hv_y *= s;
            }
        }
        fx += ibp_mass * a_hv_x;
        fy += ibp_mass * a_hv_y;
    }

    /* Acceleration CFL: dt ≤ ACC_CFL_FRAC · sqrt(h_i / |a_i|).
     * Mode 2 (Sedov shock-front): |a| ~ 1e+10 → free-fall time across cell
     * collapses. Constrain dt so cell displacement under a stays small
     * fraction of h. Preserves conservation (only shrinks dt). */
    {
        double ax_i = fx / ibp_mass;
        double ay_i = fy / ibp_mass;
        double amag2 = ax_i*ax_i + ay_i*ay_i;
        if (amag2 > 0 && ibp_volume > 0) {
            double amag = sqrt(amag2);
            double h_i = sqrt(ibp_volume);
            double dt_acc = ACC_CFL_FRAC * sqrt(h_i / amag);
            if (dt_acc < (double)my_dt) my_dt = (float)dt_acc;
        }
    }

    /* Fix A: per-cell entropy gate — dK_diss must be >= 0.
     * The face-pair NS viscous term Σ_faces τ·dS·(v_f-v_i) is positive
     * in the well-resolved limit (= q_NS·V >= 0), but individual face
     * terms can be negative, and discretization noise/boundary effects
     * can drive the cell sum slightly negative.  Clamping enforces the
     * 2nd law per cell. Phase 1 keeps dK = 0 and skips this CFL. */
    if (!phase1 && dK_diss < 0) dK_diss = 0;
    /* Convert dissipative die to entropy rate dK/dt (mirrors exam.c:4818-4830):
     *   dK/dt = (γ-1) · q_diss / ρ^γ                    (per unit volume)
     *         = (γ-1) · dK_diss / (V_cell · ρ^γ)
     *         = (γ-1) · dK_diss / (m · ρ^(γ-1))         (m = ρ·V)
     * Always written; CPU integrator uses it iff entropy_mode==1. */
    double dK_val = 0.0;
    if (!phase1 && isfinite(dK_diss) && ibp_mass > 0 && ibp_den > 0) {
        /* T8 (T4 component): floor ρ in dK conversion to prevent low-ρ
         * amplifier blowing up. Pairs with T7b entropy diffusion. */
        double rho_eff = fmax(ibp_den, 0.01);
        double rho_pow = pow(rho_eff, Gamma - 1.0);
        double cand = (Gamma - 1.0) * dK_diss / (ibp_mass * rho_pow);
        if (isfinite(cand)) dK_val = cand;
    }

    /* dK-aware CFL: dt ≤ DKCFL_FRAC · K / |dK|.
     * K derived locally from P/ρ^γ (avoids extra kernel input). */
    {
        double dK_abs = fabs(dK_val);
        if (!phase1 && dK_abs > 0 && ibp_den > 0) {
            double rho_pow_K = pow(fmax(ibp_den, 1.0e-30), Gamma);
            double K_eff = ibp_pressure / rho_pow_K;
            if (K_eff < DKCFL_KFLOOR) K_eff = DKCFL_KFLOOR;
            double dt_dK = DKCFL_FRAC * K_eff / dK_abs;
            if (dt_dK < (double)my_dt) my_dt = (float)dt_dK;
        }
    }

    /* Write outputs */
    ax_out[i]       = fx / ibp_mass;
    ay_out[i]       = fy / ibp_mass;
    /* Phase 1: die is dE/dt. The CPU kick rebuilds ie from E − kinetic. */
    die_out[i]      = (float)(phase1 ? dte : die);
    if (thermal64) out_Udot64[i] = thermal_rate64;
    dK_out[i]       = phase1 ? 0.f : (float)dK_val;
    dt_out[i]       = my_dt;
    vsig_max_out[i] = my_vsig_max;

    /* XSPH face-averaged velocity offset (zero if no real neighbors) */
    if (xsph_w_acc > 0) {
        vsmoothx_out[i] = xsph_vx_acc / xsph_w_acc;
        vsmoothy_out[i] = xsph_vy_acc / xsph_w_acc;
    } else {
        vsmoothx_out[i] = 0;
        vsmoothy_out[i] = 0;
    }
}

/* ================================================================
 *  GPU Memory Management
 * ================================================================ */

/* Helper: allocate pinned host + device arrays for SoA */
static void alloc_particle_soa(ParticleSoA *h, ParticleSoA *d, int n)
{
    size_t sd = n * sizeof(double);
    size_t sf = n * sizeof(float);
    size_t sl = n * sizeof(long long);

#define ALLOC_PAIR_D(field) \
    CUDA_CHECK(cudaMallocHost(&h->field, sd)); \
    CUDA_CHECK(cudaMalloc(&d->field, sd))
#define ALLOC_PAIR_F(field) \
    CUDA_CHECK(cudaMallocHost(&h->field, sf)); \
    CUDA_CHECK(cudaMalloc(&d->field, sf))

    ALLOC_PAIR_D(x);       ALLOC_PAIR_D(y);
    ALLOC_PAIR_D(vx);      ALLOC_PAIR_D(vy);
    ALLOC_PAIR_D(mass);
    /* indx for ghost detection */
    CUDA_CHECK(cudaMallocHost(&h->indx, sl));
    CUDA_CHECK(cudaMalloc(&d->indx, sl));
    ALLOC_PAIR_F(den);      ALLOC_PAIR_F(pressure);
    ALLOC_PAIR_F(csound);   ALLOC_PAIR_F(volume);
    ALLOC_PAIR_F(ie);       ALLOC_PAIR_F(w2);      ALLOC_PAIR_F(w2old);
    ALLOC_PAIR_F(w2ceil);   ALLOC_PAIR_F(avgNeighP);
    ALLOC_PAIR_F(K);
    ALLOC_PAIR_D(gUxx);     ALLOC_PAIR_D(gUxy);
    ALLOC_PAIR_D(gUyx);     ALLOC_PAIR_D(gUyy);
    ALLOC_PAIR_D(dPdx);     ALLOC_PAIR_D(dPdy);
    ALLOC_PAIR_D(dRhodx);   ALLOC_PAIR_D(dRhody);
    ALLOC_PAIR_D(tauxx);    ALLOC_PAIR_D(tauxy);   ALLOC_PAIR_D(tauyy);
    ALLOC_PAIR_D(alpha_cd); ALLOC_PAIR_D(divv);
    /* LagMFM fields */
    ALLOC_PAIR_D(E_inv_xx); ALLOC_PAIR_D(E_inv_xy);
    ALLOC_PAIR_D(E_inv_yx); ALLOC_PAIR_D(E_inv_yy);
    ALLOC_PAIR_D(h_mfm);

    ALLOC_PAIR_D(ax_out);   ALLOC_PAIR_D(ay_out);
    ALLOC_PAIR_F(die_out);  ALLOC_PAIR_F(dK_out);  ALLOC_PAIR_F(dt_out);
    ALLOC_PAIR_D(vsig_max_out);
    ALLOC_PAIR_D(vsmoothx_out); ALLOC_PAIR_D(vsmoothy_out);
    ALLOC_PAIR_D(lap_vx);       ALLOC_PAIR_D(lap_vy);

#undef ALLOC_PAIR_D
#undef ALLOC_PAIR_F
}

static void free_particle_soa(ParticleSoA *h, ParticleSoA *d)
{
#define FREE_PAIR_H(field) if(h->field) { cudaFreeHost(h->field); h->field = NULL; }
#define FREE_PAIR_D(field) if(d->field) { cudaFree(d->field); d->field = NULL; }
#define FREE_PAIR(field) FREE_PAIR_H(field); FREE_PAIR_D(field)

    FREE_PAIR(x);       FREE_PAIR(y);
    FREE_PAIR(vx);      FREE_PAIR(vy);
    FREE_PAIR(mass);
    FREE_PAIR(indx);
    FREE_PAIR(den);      FREE_PAIR(pressure);
    FREE_PAIR(csound);   FREE_PAIR(volume);
    FREE_PAIR(ie);       FREE_PAIR(w2);      FREE_PAIR(w2old);
    FREE_PAIR(w2ceil);   FREE_PAIR(avgNeighP);
    FREE_PAIR(K);
    FREE_PAIR(gUxx);     FREE_PAIR(gUxy);
    FREE_PAIR(gUyx);     FREE_PAIR(gUyy);
    FREE_PAIR(dPdx);     FREE_PAIR(dPdy);
    FREE_PAIR(dRhodx);   FREE_PAIR(dRhody);
    FREE_PAIR(tauxx);    FREE_PAIR(tauxy);   FREE_PAIR(tauyy);
    FREE_PAIR(alpha_cd); FREE_PAIR(divv);
    FREE_PAIR(E_inv_xx); FREE_PAIR(E_inv_xy);
    FREE_PAIR(E_inv_yx); FREE_PAIR(E_inv_yy);
    FREE_PAIR(h_mfm);

    FREE_PAIR(ax_out);   FREE_PAIR(ay_out);
    FREE_PAIR(die_out);  FREE_PAIR(dK_out);  FREE_PAIR(dt_out);
    FREE_PAIR(vsig_max_out);
    FREE_PAIR(vsmoothx_out); FREE_PAIR(vsmoothy_out);
    FREE_PAIR(lap_vx);       FREE_PAIR(lap_vy);

#undef FREE_PAIR
#undef FREE_PAIR_H
#undef FREE_PAIR_D
}

static void alloc_face_csr(FaceCSR *h, FaceCSR *d, int max_p, int max_f)
{
    size_t si = max_f * sizeof(int);
    size_t sd = max_f * sizeof(double);
    size_t so = (max_p + 1) * sizeof(int);

    CUDA_CHECK(cudaMallocHost(&h->face_offset, so));
    CUDA_CHECK(cudaMalloc(&d->face_offset, so));

#define ALLOC_F_D(field) CUDA_CHECK(cudaMallocHost(&h->field, sd)); CUDA_CHECK(cudaMalloc(&d->field, sd))
#define ALLOC_F_I(field) CUDA_CHECK(cudaMallocHost(&h->field, si)); CUDA_CHECK(cudaMalloc(&d->field, si))

    ALLOC_F_D(c1x);  ALLOC_F_D(c1y);
    ALLOC_F_D(c2x);  ALLOC_F_D(c2y);
    ALLOC_F_I(neighbor_idx);
    ALLOC_F_I(kp_idx);  ALLOC_F_I(km_idx);
    ALLOC_F_I(is_ghost);

#undef ALLOC_F_D
#undef ALLOC_F_I
}

static void free_face_csr(FaceCSR *h, FaceCSR *d)
{
#define FREE_CSR(field) \
    if(h->field) { cudaFreeHost(h->field); h->field = NULL; } \
    if(d->field) { cudaFree(d->field); d->field = NULL; }

    FREE_CSR(face_offset);
    FREE_CSR(c1x);  FREE_CSR(c1y);
    FREE_CSR(c2x);  FREE_CSR(c2y);
    FREE_CSR(neighbor_idx);
    FREE_CSR(kp_idx);  FREE_CSR(km_idx);
    FREE_CSR(is_ghost);

#undef FREE_CSR
}

extern "C"
void gpu_init(GPUContext *ctx, int max_particles, int max_faces, int mpi_rank)
{
    memset(ctx, 0, sizeof(GPUContext));
    ctx->max_particles = max_particles;
    ctx->max_faces     = max_faces;

    /* --- Select GPU device: rank % n_devices --- */
    {
        int n_devices = 0;
        CUDA_CHECK(cudaGetDeviceCount(&n_devices));
        if (n_devices == 0) {
            fprintf(stderr, "[GPU] FATAL: no CUDA devices found\n");
            exit(1);
        }
        int device_id = mpi_rank % n_devices;
        CUDA_CHECK(cudaSetDevice(device_id));

        cudaDeviceProp prop;
        CUDA_CHECK(cudaGetDeviceProperties(&prop, device_id));
        fprintf(stderr, "[GPU] rank %d -> device %d/%d: %s (%.0f MB)\n",
                mpi_rank, device_id, n_devices, prop.name,
                prop.totalGlobalMem / 1048576.0);
    }

    alloc_particle_soa(&ctx->h_parts, &ctx->d_parts, max_particles);
    alloc_face_csr(&ctx->h_faces, &ctx->d_faces, max_particles, max_faces);
    CUDA_CHECK(cudaMalloc(&ctx->d_tess_den_input, max_particles * sizeof(float)));

    /* --- Create CUDA stream --- */
    {
        cudaStream_t s;
        CUDA_CHECK(cudaStreamCreate(&s));
        ctx->stream = (void *)s;
    }

    /* --- Pre-allocate CUB scratch for DeviceReduce::Min --- */
    {
        CUDA_CHECK(cudaMalloc((float **)&ctx->d_cub_min_out, sizeof(float)));
        CUDA_CHECK(cudaMallocHost((float **)&ctx->h_cub_min_out, sizeof(float)));

        void *d_temp = NULL;
        size_t temp_bytes = 0;
        cub::DeviceReduce::Min(d_temp, temp_bytes,
                               ctx->d_parts.dt_out,
                               (float *)ctx->d_cub_min_out,
                               max_particles,
                               (cudaStream_t)ctx->stream);
        ctx->cub_temp_bytes = temp_bytes;
        CUDA_CHECK(cudaMalloc(&ctx->d_cub_temp, temp_bytes));
    }

    ctx->initialized = 1;
    fprintf(stderr, "[GPU] Initialized: max_particles=%d, max_faces=%d\n",
            max_particles, max_faces);
}

extern "C"
void gpu_free(GPUContext *ctx)
{
    if (!ctx->initialized) return;
    free_particle_soa(&ctx->h_parts, &ctx->d_parts);
    free_face_csr(&ctx->h_faces, &ctx->d_faces);
    gpu_free_tess_buffers(ctx);
    gpu_free_tree_buffers(ctx);
    if (ctx->d_tess_den_input) {
        cudaFree(ctx->d_tess_den_input);
        ctx->d_tess_den_input = NULL;
    }
    if (ctx->d_cub_temp) { cudaFree(ctx->d_cub_temp); ctx->d_cub_temp = NULL; }
    if (ctx->d_cub_min_out) { cudaFree(ctx->d_cub_min_out); ctx->d_cub_min_out = NULL; }
    if (ctx->h_cub_min_out) { cudaFreeHost(ctx->h_cub_min_out); ctx->h_cub_min_out = NULL; }
    if (ctx->stream) { cudaStreamDestroy((cudaStream_t)ctx->stream); ctx->stream = NULL; }
    ctx->cub_temp_bytes = 0;
    ctx->initialized = 0;
}

/* ================================================================
 *  Upload / Download
 * ================================================================ */

extern "C"
void gpu_upload_particles(GPUContext *ctx, int n, int has_stress)
{
    cudaStream_t s = (cudaStream_t)ctx->stream;
    size_t sd = n * sizeof(double);
    size_t sf = n * sizeof(float);
    size_t sl = n * sizeof(long long);

#define UP_D(field) CUDA_CHECK(cudaMemcpyAsync(ctx->d_parts.field, ctx->h_parts.field, sd, cudaMemcpyHostToDevice, s))
#define UP_F(field) CUDA_CHECK(cudaMemcpyAsync(ctx->d_parts.field, ctx->h_parts.field, sf, cudaMemcpyHostToDevice, s))

    UP_D(x);  UP_D(y);  UP_D(vx);  UP_D(vy);  UP_D(mass);
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_parts.indx, ctx->h_parts.indx, sl, cudaMemcpyHostToDevice, s));
    UP_F(den);  UP_F(pressure);  UP_F(csound);  UP_F(volume);
    UP_F(ie);  UP_F(w2);  UP_F(w2old);

    if (has_stress) {
        UP_D(gUxx);  UP_D(gUxy);  UP_D(gUyx);  UP_D(gUyy);
        UP_D(dPdx);  UP_D(dPdy);  UP_D(dRhodx);  UP_D(dRhody);
        UP_D(tauxx);  UP_D(tauxy);  UP_D(tauyy);
        UP_D(alpha_cd);  UP_D(divv);
        UP_D(lap_vx);  UP_D(lap_vy);
        UP_F(K);  /* entropy variable; mirrors sbp->stress.K */
        if (gfs_thermal64_requested()) {
            UP_D(E_inv_xx); UP_D(E_inv_xy); UP_D(E_inv_yx); UP_D(E_inv_yy);
        }
    }

#undef UP_D
#undef UP_F
}

extern "C"
void gpu_upload_faces(GPUContext *ctx, int n_particles, int n_faces)
{
    cudaStream_t s = (cudaStream_t)ctx->stream;
    size_t si = n_faces * sizeof(int);
    size_t sd = n_faces * sizeof(double);
    size_t so = (n_particles + 1) * sizeof(int);

    CUDA_CHECK(cudaMemcpyAsync(ctx->d_faces.face_offset, ctx->h_faces.face_offset, so, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_faces.c1x, ctx->h_faces.c1x, sd, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_faces.c1y, ctx->h_faces.c1y, sd, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_faces.c2x, ctx->h_faces.c2x, sd, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_faces.c2y, ctx->h_faces.c2y, sd, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_faces.neighbor_idx, ctx->h_faces.neighbor_idx, si, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_faces.kp_idx, ctx->h_faces.kp_idx, si, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_faces.km_idx, ctx->h_faces.km_idx, si, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_faces.is_ghost, ctx->h_faces.is_ghost, si, cudaMemcpyHostToDevice, s));
}

extern "C"
void gpu_download_face_csr(GPUContext *ctx, int n_particles, int n_faces,
    int *h_face_offset, double *h_c1x, double *h_c1y, double *h_c2x, double *h_c2y,
    int *h_neighbor_idx, int *h_kp_idx, int *h_km_idx, int *h_is_ghost)
{
    cudaStream_t s = (cudaStream_t)ctx->stream;
    size_t si = n_faces * sizeof(int);
    size_t sd = n_faces * sizeof(double);
    size_t so = (n_particles + 1) * sizeof(int);
    CUDA_CHECK(cudaMemcpyAsync(h_face_offset, ctx->d_faces.face_offset, so, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(h_c1x, ctx->d_faces.c1x, sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(h_c1y, ctx->d_faces.c1y, sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(h_c2x, ctx->d_faces.c2x, sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(h_c2y, ctx->d_faces.c2y, sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(h_neighbor_idx, ctx->d_faces.neighbor_idx, si, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(h_kp_idx, ctx->d_faces.kp_idx, si, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(h_km_idx, ctx->d_faces.km_idx, si, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(h_is_ghost, ctx->d_faces.is_ghost, si, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaStreamSynchronize(s));
}

extern "C"
void gpu_download_results(GPUContext *ctx, int n)
{
    cudaStream_t s = (cudaStream_t)ctx->stream;
    size_t sd = n * sizeof(double);
    size_t sf = n * sizeof(float);

    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.ax_out, ctx->d_parts.ax_out, sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.ay_out, ctx->d_parts.ay_out, sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.die_out, ctx->d_parts.die_out, sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dK_out, ctx->d_parts.dK_out, sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dt_out, ctx->d_parts.dt_out, sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.vsig_max_out, ctx->d_parts.vsig_max_out, sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.vsmoothx_out, ctx->d_parts.vsmoothx_out, sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.vsmoothy_out, ctx->d_parts.vsmoothy_out, sd, cudaMemcpyDeviceToHost, s));
    if (gfs_thermal64_requested())
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.E_inv_yx, ctx->d_parts.E_inv_yx, sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaStreamSynchronize(s));
}

#ifdef GFS_VARPRESS_W0
static __device__ unsigned long long d_local_periodic_sync_errors;
__global__ void refresh_local_periodic_fields_kernel(ParticleSoA p,
    int nbp,int n_total,double boxsize)
{
    int h=nbp+blockIdx.x*blockDim.x+threadIdx.x;
    if(h>=n_total) return;
    long long id=p.indx[h];
    if(id<0 || id>=nbp || p.indx[(int)id]!=id) {
        atomicAdd(&d_local_periodic_sync_errors,1ull); return;
    }
    int j=(int)id;
    double dx=p.x[h]-p.x[j],dy=p.y[h]-p.y[j];
    double tol=128.0*DBL_EPSILON*boxsize;
    if(fabs(dx-boxsize*nearbyint(dx/boxsize))>tol ||
       fabs(dy-boxsize*nearbyint(dy/boxsize))>tol ||
       p.vx[h]!=p.vx[j] || p.vy[h]!=p.vy[j] || p.mass[h]!=p.mass[j] ||
       p.ie[h]!=p.ie[j] || p.w2[h]!=0.0f || p.w2[j]!=0.0f) {
        atomicAdd(&d_local_periodic_sync_errors,1ull); return;
    }
#define GFS_ALIAS(field) p.field[h]=p.field[j]
    GFS_ALIAS(den); GFS_ALIAS(pressure); GFS_ALIAS(csound); GFS_ALIAS(volume);
    GFS_ALIAS(w2ceil); GFS_ALIAS(avgNeighP); GFS_ALIAS(K);
    GFS_ALIAS(gUxx); GFS_ALIAS(gUxy); GFS_ALIAS(gUyx); GFS_ALIAS(gUyy);
    GFS_ALIAS(dPdx); GFS_ALIAS(dPdy); GFS_ALIAS(dRhodx); GFS_ALIAS(dRhody);
    GFS_ALIAS(tauxx); GFS_ALIAS(tauxy); GFS_ALIAS(tauyy); GFS_ALIAS(divv);
    GFS_ALIAS(lap_vx); GFS_ALIAS(lap_vy); GFS_ALIAS(alpha_cd);
#undef GFS_ALIAS
}
#endif

extern "C" int gpu_refresh_local_periodic_fields(GPUContext *ctx,
    int nbp,int n_total,double boxsize)
{
#ifdef GFS_VARPRESS_W0
    if(n_total<nbp || !(boxsize>0)) return 0;
    if(n_total==nbp) return 1;
    cudaStream_t stream=(cudaStream_t)ctx->stream;
    unsigned long long errors=0;
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_local_periodic_sync_errors,&errors,sizeof(errors),0,cudaMemcpyHostToDevice,stream));
    refresh_local_periodic_fields_kernel<<<(n_total-nbp+255)/256,256,0,stream>>>(ctx->d_parts,nbp,n_total,boxsize);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyFromSymbolAsync(&errors,d_local_periodic_sync_errors,sizeof(errors),0,cudaMemcpyDeviceToHost,stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));
    /* A failed alias can leave partial device copies. The caller then does
     * the full host refresh/upload before any force evaluation. */
    return errors==0;
#else
    (void)ctx;(void)nbp;(void)n_total;(void)boxsize;
    return 0;
#endif
}

/* Host-only structural snapshot used with GFS_GPU_FACE_DUMP_CALL.
 * It reads the polygons and linked-cell data produced by the tessellation
 * kernel. It does not feed any values back into the solver. */
static void dump_kepler_tess_snapshot(GPUContext *ctx, int n_particles,
                                      int n_total, int mpi_rank, int force_call,
                                      int mx, int my, double cellsize,
                                      double xmin, double ymin)
{
    const double cx0 = getenv("EUNHA_KEPLER_CX") ? atof(getenv("EUNHA_KEPLER_CX")) : 12.0;
    const double cy0 = getenv("EUNHA_KEPLER_CY") ? atof(getenv("EUNHA_KEPLER_CY")) : 12.0;
    const double r_cut = 2.5;
    PolygonOut *poly = (PolygonOut *)malloc((size_t)n_particles * sizeof(PolygonOut));
    float *w2 = (float *)malloc((size_t)n_total * sizeof(float));
    if (!poly || !w2) {
        fprintf(stderr, "[KEPLER_TESS_DUMP] rank=%d allocation failed\n", mpi_rank);
        free(poly); free(w2);
        exit(1);
    }
    CUDA_CHECK(cudaMemcpy(poly, ctx->d_polygons,
                          (size_t)n_particles * sizeof(PolygonOut), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(w2, ctx->d_parts.w2,
                          (size_t)n_total * sizeof(float), cudaMemcpyDeviceToHost));

    char path[160];
    snprintf(path, sizeof(path), "kepler_tess_diag_r%d_call%d.tsv", mpi_rank, force_call);
    FILE *out = fopen(path, "w");
    if (!out) {
        fprintf(stderr, "[KEPLER_TESS_DUMP] rank=%d could not open %s\n", mpi_rank, path);
        free(poly); free(w2);
        exit(1);
    }
    const int max_candidates = ctx->h_cells.n_entries;
    fprintf(out, "# rank=%d call=%d nbp=%d n_total=%d mx=%d my=%d cellsize=%.17g xmin=%.17g ymin=%.17g n_cells=%d n_entries=%d center=%.17g,%.17g r_cut=%.17g max_neigh=%d\n",
            mpi_rank, force_call, n_particles, n_total, mx, my,
            cellsize, xmin, ymin,
            ctx->h_cells.n_cells, ctx->h_cells.n_entries, cx0, cy0, r_cut, max_candidates);
    fprintf(out, "# CELL rank owner_soa owner_id x y binx biny nc raw_area maxR2 floor_active n_raw n_stored n_truncated min_stencil_drad2 n_outside_le_R2 min_out_id min_out_soa min_out_drad2 min_out_x min_out_y min_out_binx min_out_biny\n");
    fprintf(out, "# VERT rank owner_soa owner_id corner xrel yrel edge_neigh_soa edge_neigh_id xworld yworld\n");
    fprintf(out, "# CAND rank owner_soa owner_id sorted_order soa id dx dy w2 drad2 binx biny\n");

    if (mx <= 0 || my <= 0 || !(cellsize > 0.0) ||
        ctx->h_cells.n_cells != mx * my) {
        fprintf(out, "[KEPLER_TESS_DUMP] rank=%d inconsistent cell grid: mx=%d my=%d n_cells=%d\n",
                mpi_rank, mx, my, ctx->h_cells.n_cells);
        fclose(out); free(poly); free(w2);
        return;
    }

    typedef struct {
        int soa, bx, by;
        double dx, dy, w2, drad2;
    } CandidateDump;
    CandidateDump *cand = (CandidateDump *)malloc((size_t)max_candidates * sizeof(*cand));
    if (!cand) {
        fprintf(stderr, "[KEPLER_TESS_DUMP] rank=%d candidate allocation failed (%d entries)\n",
                mpi_rank, max_candidates);
        fclose(out); free(poly); free(w2);
        exit(1);
    }
    const double volume_floor = 1.0e-3 * cellsize * cellsize;

    for (int i = 0; i < n_particles; i++) {
        const double xi = ctx->h_parts.x[i], yi = ctx->h_parts.y[i];
        const double rx = xi - cx0, ry = yi - cy0;
        if (rx * rx + ry * ry >= r_cut * r_cut) continue;

        int ix = (int)((xi - xmin) / cellsize);
        int iy = (int)((yi - ymin) / cellsize);
        if (ix < 0) ix = 0; if (ix >= mx) ix = mx - 1;
        if (iy < 0) iy = 0; if (iy >= my) iy = my - 1;

        const long long idi = ctx->h_parts.indx[i];
        const PolygonOut *p = &poly[i];
        const int nc = (p->n_corners < 0) ? 0 :
                       (p->n_corners > GPU_MAX_CORNERS ? GPU_MAX_CORNERS : p->n_corners);
        double raw_area2 = 0.0, maxR2 = 0.0;
        for (int e = 0; e < nc; e++) {
            int en = (e + 1) % nc;
            raw_area2 += p->cx[e] * p->cy[en] - p->cy[e] * p->cx[en];
            double r2 = p->cx[e] * p->cx[e] + p->cy[e] * p->cy[e];
            if (r2 > maxR2) maxR2 = r2;
        }
        double raw_area = 0.5 * fabs(raw_area2);

        int n_stored = 0;
        long long n_raw = 0;
        double min_stencil_drad2 = DBL_MAX;
        for (int dyc = -2; dyc <= 2; dyc++) {
            int cy = iy + dyc;
            if (cy < 0 || cy >= my) continue;
            for (int dxc = -2; dxc <= 2; dxc++) {
                int cx = ix + dxc;
                if (cx < 0 || cx >= mx) continue;
                int cell_id = cx + mx * cy;
                int begin = ctx->h_cells.cell_offset[cell_id];
                int end = ctx->h_cells.cell_offset[cell_id + 1];
                for (int k = begin; k < end; k++) {
                    int j = ctx->h_cells.particle_index[k];
                    if (j < 0 || j >= n_total || ctx->h_parts.indx[j] == idi) continue;
                    double dx = ctx->h_parts.x[j] - xi;
                    double dy = ctx->h_parts.y[j] - yi;
                    double dist2 = dx * dx + dy * dy;
                    double w2j = (double)w2[j];
                    double wdiff = ((double)w2[i] - w2j) / dist2;
                    double wfrac = 0.5 + 0.5 * wdiff;
                    if (wfrac < 0.05) wfrac = 0.05; else if (wfrac > 0.95) wfrac = 0.95;
                    /* Sort/filter by the actual clipped bisector offset. */
                    double drad2 = wfrac * wfrac * dist2;
                    n_raw++;
                    if (n_stored < max_candidates) {
                        cand[n_stored].soa = j; cand[n_stored].bx = cx; cand[n_stored].by = cy;
                        cand[n_stored].dx = dx; cand[n_stored].dy = dy;
                        cand[n_stored].w2 = w2j; cand[n_stored].drad2 = drad2;
                        if (drad2 < min_stencil_drad2) min_stencil_drad2 = drad2;
                        n_stored++;
                    }
                }
            }
        }
        /* Match the GPU's stable insertion sort on the stored candidate prefix. */
        for (int a = 1; a < n_stored; a++) {
            CandidateDump key = cand[a];
            int b = a - 1;
            while (b >= 0 && cand[b].drad2 > key.drad2) {
                cand[b + 1] = cand[b]; b--;
            }
            cand[b + 1] = key;
        }

        long long n_outside_le_R2 = 0;
        double min_out_drad2 = DBL_MAX, min_out_dx = 0, min_out_dy = 0;
        long long min_out_id = -1;
        int min_out_soa = -1, min_out_bx = -1, min_out_by = -1;
        for (int j = 0; j < n_total; j++) {
            if (ctx->h_parts.indx[j] == idi) continue;
            int jx = (int)((ctx->h_parts.x[j] - xmin) / cellsize);
            int jy = (int)((ctx->h_parts.y[j] - ymin) / cellsize);
            if (jx < 0) jx = 0; if (jx >= mx) jx = mx - 1;
            if (jy < 0) jy = 0; if (jy >= my) jy = my - 1;
            if (abs(jx - ix) <= 2 && abs(jy - iy) <= 2) continue;
            double dx = ctx->h_parts.x[j] - xi, dy = ctx->h_parts.y[j] - yi;
            double dist2 = dx * dx + dy * dy;
            double w2j = (double)w2[j];
            double wdiff = ((double)w2[i] - w2j) / dist2;
            double wfrac = 0.5 + 0.5 * wdiff;
            if (wfrac < 0.05) wfrac = 0.05; else if (wfrac > 0.95) wfrac = 0.95;
            double drad2 = wfrac * wfrac * dist2;
            if (drad2 <= maxR2) n_outside_le_R2++;
            if (drad2 < min_out_drad2) {
                min_out_drad2 = drad2; min_out_dx = dx; min_out_dy = dy;
                min_out_id = ctx->h_parts.indx[j]; min_out_soa = j; min_out_bx = jx; min_out_by = jy;
            }
        }

        fprintf(out, "CELL\t%d\t%d\t%lld\t%.17g\t%.17g\t%d\t%d\t%d\t%.17g\t%.17g\t%d\t%lld\t%d\t%lld\t%.17g\t%lld\t%lld\t%d\t%.17g\t%.17g\t%.17g\t%d\t%d\n",
                mpi_rank, i, idi, xi, yi, ix, iy, nc, raw_area, maxR2,
                raw_area < volume_floor, n_raw, n_stored, n_raw - n_stored,
                min_stencil_drad2, n_outside_le_R2, min_out_id, min_out_soa,
                min_out_drad2, xi + min_out_dx, yi + min_out_dy,
                min_out_bx, min_out_by);

        for (int e = 0; e < nc; e++) {
            int j = p->edge_neigh[e];
            long long idj = (j >= 0 && j < n_total) ? ctx->h_parts.indx[j] : -1;
            fprintf(out, "VERT\t%d\t%d\t%lld\t%d\t%.17g\t%.17g\t%d\t%lld\t%.17g\t%.17g\n",
                    mpi_rank, i, idi, e, p->cx[e], p->cy[e], j, idj,
                    xi + p->cx[e], yi + p->cy[e]);
        }
        for (int a = 0; a < n_stored; a++) {
            int j = cand[a].soa;
            fprintf(out, "CAND\t%d\t%d\t%lld\t%d\t%d\t%lld\t%.17g\t%.17g\t%.17g\t%.17g\t%d\t%d\n",
                    mpi_rank, i, idi, a, j, ctx->h_parts.indx[j],
                    cand[a].dx, cand[a].dy, cand[a].w2, cand[a].drad2,
                    cand[a].bx, cand[a].by);
        }
    }

    /* Full-domain audit of generators omitted by the fixed 5x5 stencil.
     * The detailed dump above is intentionally limited to the Kepler core;
     * this compact audit covers every owned cell and asks whether an omitted
     * bisector would actually cut its GPU polygon.  For a weighted face, a
     * conservative search radius follows from the clamped wfrac in [0.05,.95]
     * and the global maximum w2. */
    {
        char audit_path[160];
        snprintf(audit_path, sizeof(audit_path),
                 "kepler_stencil_audit_r%d_call%d.tsv", mpi_rank, force_call);
        FILE *audit = fopen(audit_path, "w");
        if (!audit) {
            fprintf(stderr, "[KEPLER_STENCIL_AUDIT] rank=%d could not open %s\n",
                    mpi_rank, audit_path);
            free(cand); free(poly); free(w2);
            exit(1);
        }
        fprintf(audit, "# rank owner_soa owner_id x y area maxR2 search_radius n_outside_possible n_outside_cut n_outside_any_cut max_violation worst_neighbor_id worst_neighbor_x worst_neighbor_y halo_search_clipped\n");
        const double global_w2_max = h_diag_global_w2_max;

        long long audit_possible_total = 0;
        long long audit_cut_total = 0;
        long long audit_any_cut_total = 0;
        int audit_cells_with_cut = 0;
        int audit_halo_clipped_cells = 0;
        double audit_max_violation = 0.0;

        for (int i = 0; i < n_particles; i++) {
            const double xi = ctx->h_parts.x[i], yi = ctx->h_parts.y[i];
            int ix = (int)((xi - xmin) / cellsize);
            int iy = (int)((yi - ymin) / cellsize);
            if (ix < 0) ix = 0; if (ix >= mx) ix = mx - 1;
            if (iy < 0) iy = 0; if (iy >= my) iy = my - 1;

            const PolygonOut *p = &poly[i];
            const int nc = (p->n_corners < 0) ? 0 :
                           (p->n_corners > GPU_MAX_CORNERS ? GPU_MAX_CORNERS : p->n_corners);
            double maxR2 = 0.0, area2 = 0.0;
            for (int e = 0; e < nc; e++) {
                const int en = (e + 1) % nc;
                area2 += p->cx[e] * p->cy[en] - p->cy[e] * p->cx[en];
                const double r2 = p->cx[e] * p->cx[e] + p->cy[e] * p->cy[e];
                if (r2 > maxR2) maxR2 = r2;
            }
            const double radius = sqrt(maxR2);
            const double deltaW = fmax(0.0, global_w2_max - (double)w2[i]);
            const double d_clamped = (deltaW > 0.0)
                ? fmin(sqrt(deltaW / 0.9), 20.0 * radius) : 0.0;
            const double d_unclamped = radius + sqrt(maxR2 + deltaW);
            double search_radius = fmax(2.0 * radius, fmax(d_clamped, d_unclamped));
            if (!(search_radius >= 0.0) || !isfinite(search_radius))
                search_radius = hypot((double)mx * cellsize, (double)my * cellsize);
            int pad = (int)ceil(search_radius / cellsize) + 1;
            if (pad < 3) pad = 3;
            if (pad > (mx > my ? mx : my)) pad = mx > my ? mx : my;
            const double grid_xmax = xmin + (double)mx * cellsize;
            const double grid_ymax = ymin + (double)my * cellsize;
            const int halo_search_clipped =
                (xi - search_radius < xmin || xi + search_radius > grid_xmax ||
                 yi - search_radius < ymin || yi + search_radius > grid_ymax);
            if (halo_search_clipped) audit_halo_clipped_cells++;
            int xlo = ix - pad, xhi = ix + pad;
            int ylo = iy - pad, yhi = iy + pad;
            if (xlo < 0) xlo = 0; if (xhi >= mx) xhi = mx - 1;
            if (ylo < 0) ylo = 0; if (yhi >= my) yhi = my - 1;

            long long n_possible = 0, n_cut = 0, n_any_cut = 0;
            double max_violation = 0.0;
            long long worst_id = -1;
            double worst_x = 0.0, worst_y = 0.0;
            const long long owner_id = ctx->h_parts.indx[i];
            for (int cy = ylo; cy <= yhi; cy++) {
                for (int cx = xlo; cx <= xhi; cx++) {
                    if (abs(cx - ix) <= 2 && abs(cy - iy) <= 2) continue;
                    const int cell_id = cx + mx * cy;
                    const int begin = ctx->h_cells.cell_offset[cell_id];
                    const int end = ctx->h_cells.cell_offset[cell_id + 1];
                    for (int k = begin; k < end; k++) {
                        const int j = ctx->h_cells.particle_index[k];
                        if (j < 0 || j >= n_total || ctx->h_parts.indx[j] == owner_id) continue;
                        const double dx = ctx->h_parts.x[j] - xi;
                        const double dy = ctx->h_parts.y[j] - yi;
                        const double dist2 = dx * dx + dy * dy;
                        if (!(dist2 > 0.0)) continue;
                        double wfrac = 0.5 + 0.5 * ((double)w2[i] - (double)w2[j]) / dist2;
                        if (wfrac < 0.05) wfrac = 0.05;
                        else if (wfrac > 0.95) wfrac = 0.95;
                        const double drad2 = wfrac * wfrac * dist2;
                        if (drad2 > maxR2) continue;
                        n_possible++;

                        const double threshold = wfrac * dist2;
                        double violation = -DBL_MAX;
                        for (int e = 0; e < nc; e++) {
                            const double v = p->cx[e] * dx + p->cy[e] * dy - threshold;
                            if (v > violation) violation = v;
                        }
                        if (violation > 0.0) n_any_cut++;
                        const double tol = 64.0 * DBL_EPSILON * fmax(1.0, dist2);
                        if (violation > tol) {
                            n_cut++;
                            if (violation > max_violation) {
                                max_violation = violation;
                                worst_id = ctx->h_parts.indx[j];
                                worst_x = ctx->h_parts.x[j];
                                worst_y = ctx->h_parts.y[j];
                            }
                        }
                    }
                }
            }
            audit_possible_total += n_possible;
            audit_cut_total += n_cut;
            audit_any_cut_total += n_any_cut;
            if (n_cut > 0) audit_cells_with_cut++;
            if (max_violation > audit_max_violation) audit_max_violation = max_violation;
            fprintf(audit, "%d\t%d\t%lld\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%lld\t%lld\t%lld\t%.17g\t%lld\t%.17g\t%.17g\t%d\n",
                    mpi_rank, i, owner_id, xi, yi, 0.5 * fabs(area2), maxR2,
                    search_radius, n_possible, n_cut, n_any_cut, max_violation,
                    worst_id, worst_x, worst_y, halo_search_clipped);
        }
        fclose(audit);
        fprintf(stderr, "[GPU_STENCIL_AUDIT] rank=%d file=%s cells=%d possible_5x5=%lld cut_5x5=%lld any_cut_5x5=%lld cells_with_cut_5x5=%d max_violation=%.6e halo_search_clipped=%d\n",
                mpi_rank, audit_path, n_particles, audit_possible_total,
                audit_cut_total, audit_any_cut_total, audit_cells_with_cut,
                audit_max_violation, audit_halo_clipped_cells);
    }

    fclose(out);
    fprintf(stderr, "[KEPLER_TESS_DUMP] rank=%d call=%d nbp=%d file=%s\n",
            mpi_rank, force_call, n_particles, path);
    free(cand); free(poly); free(w2);
}

/* ================================================================
 *  Kernel launch + CFL reduction
 * ================================================================ */

extern "C"
long long gpu_take_pair_hits(void)
{
    return h_gfs_pair_hits;
}

/* Bounded, read-only force-input capture. Actual face fluxes are still written
 * by the production force kernel. No branch or force arithmetic is changed. */
static int gfs_face_dump_selected(long long id, const long long *ids, int n)
{
    if (!n) return 1;
    for (int k=0; k<n; k++) if (id==ids[k]) return 1;
    return 0;
}

static void gfs_dump_force_inputs(GPUContext *ctx, int n_particles, int n_total,
    int rank, int call, const GPUPhysicsParams *params,
    const long long *ids, int nids)
{
    if (!nids) return;
    if (params->av_mode!=5 || !params->use_muscl || !params->phase1) {
        fprintf(stderr,"[FACE_INPUT] selected input capture requires phase1 MUSCL av_mode5\n"); exit(1);
    }
    unsigned char *needed=(unsigned char *)calloc((size_t)n_total,1);
    if (!needed) { fprintf(stderr,"[FACE_INPUT] allocation failed\n"); exit(1); }
    for (int i=0; i<n_particles; i++) {
        if (!gfs_face_dump_selected(ctx->h_parts.indx[i],ids,nids)) continue;
        needed[i]=1;
        for (int f=ctx->h_faces.face_offset[i]; f<ctx->h_faces.face_offset[i+1]; f++) {
            int j=ctx->h_faces.neighbor_idx[f];
            if (j>=0 && j<n_total) needed[j]=1;
        }
    }
    /* CUDA arrays contain the exact EOS/reconstruction inputs consumed by the
     * force kernel. Host particle fields from the prior RK stage are not used. */
#define GFS_CAPTURE_FIELD(field,type) CUDA_CHECK(cudaMemcpy(ctx->h_parts.field,ctx->d_parts.field,(size_t)n_total*sizeof(type),cudaMemcpyDeviceToHost))
    GFS_CAPTURE_FIELD(vx,double); GFS_CAPTURE_FIELD(vy,double);
    GFS_CAPTURE_FIELD(mass,double);
    GFS_CAPTURE_FIELD(den,float); GFS_CAPTURE_FIELD(pressure,float);
    GFS_CAPTURE_FIELD(csound,float); GFS_CAPTURE_FIELD(volume,float);
    GFS_CAPTURE_FIELD(ie,float); GFS_CAPTURE_FIELD(w2,float); GFS_CAPTURE_FIELD(w2old,float);
    GFS_CAPTURE_FIELD(gUxx,double); GFS_CAPTURE_FIELD(gUxy,double);
    GFS_CAPTURE_FIELD(gUyx,double); GFS_CAPTURE_FIELD(gUyy,double);
    GFS_CAPTURE_FIELD(dPdx,double); GFS_CAPTURE_FIELD(dPdy,double);
    GFS_CAPTURE_FIELD(dRhodx,double); GFS_CAPTURE_FIELD(dRhody,double);
#undef GFS_CAPTURE_FIELD
    char path[160]; snprintf(path,sizeof(path),"face_input_r%d_call%d.tsv",rank,call);
    FILE *out=fopen(path,"wx");
    if (!out) { fprintf(stderr,"[FACE_INPUT] cannot create %s\n",path); exit(1); }
    fprintf(out,"index\tid\tx\ty\tvx\tvy\tmass\tV\trho\tP\tcs\tU\tW\tWold\tBxx\tBxy\tByx\tByy\tdPdx\tdPdy\tdRhodx\tdRhody\n");
    for (int i=0; i<n_total; i++) if (needed[i]) {
        const ParticleSoA *h=&ctx->h_parts;
        fprintf(out,"%d\t%lld\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
            i,h->indx[i],h->x[i],h->y[i],h->vx[i],h->vy[i],h->mass[i],
            (double)h->volume[i],(double)h->den[i],(double)h->pressure[i],
            (double)h->csound[i],(double)h->ie[i],(double)h->w2[i],(double)h->w2old[i],
            h->gUxx[i],h->gUxy[i],h->gUyx[i],h->gUyy[i],h->dPdx[i],h->dPdy[i],h->dRhodx[i],h->dRhody[i]);
    }
    if (ferror(out) || fclose(out)) { fprintf(stderr,"[FACE_INPUT] write failed\n"); exit(1); }
    free(needed);
    snprintf(path,sizeof(path),"face_input_meta_r%d_call%d.json",rank,call);
    out=fopen(path,"wx");
    if (!out) { fprintf(stderr,"[FACE_INPUT] cannot create %s\n",path); exit(1); }
    fprintf(out,"{\"rank\":%d,\"call\":%d,\"Gamma\":%.17g,\"dtold\":%.17g,\"av_mode\":%d,\"use_muscl\":%d,\"phase1\":%d,\"geom_fv\":%d,\"face_charge\":%d,\"w2lin\":%d,\"alphavis\":%.17g,\"betavis\":%.17g,\"nu_phys\":%.17g,\"cd_amax\":%.17g,\"hyperv_alpha\":%.17g}\n",
        rank,call,params->Gamma,params->dtold,params->av_mode,params->use_muscl,
        params->phase1,params->geom_fv,params->face_charge,params->w2lin,
        params->alphavis,params->betavis,params->nu_phys,params->cd_amax,params->hyperv_alpha);
    if (ferror(out) || fclose(out)) { fprintf(stderr,"[FACE_INPUT] meta write failed\n"); exit(1); }
}

#include "gfs_ale_shadow.inc"

extern "C" void gfs_rk4_geometry_abort_snapshot(void);

extern "C"
double gpu_launch_force_kernel(GPUContext *ctx, int n_particles,
                               int n_total, int mpi_rank,
                               const GPUPhysicsParams *params)
{
    static int force_call = 0;
    static int dump_call = -1;
    static int dump_first=0,dump_last=0,dump_nids=0;
    static long long dump_ids[64];
    force_call++;
    double var_geom_cfl=0;
    int var_geom_audit=0;
#ifdef GFS_VARPRESS_W0
    const char *geom_cfl_env=getenv("GFS_VARPRESS_GEOM_CFL");
    if (geom_cfl_env && geom_cfl_env[0]) {
        char *end=NULL; var_geom_cfl=strtod(geom_cfl_env,&end);
        if (!end || *end || !isfinite(var_geom_cfl) || var_geom_cfl<0 || var_geom_cfl>1) {
            fprintf(stderr,"[VARPRESS_GEOM_CFL] requires a coefficient in [0,1]\n"); exit(1);
        }
    }
    const char *geom_audit_env=getenv("GFS_VARPRESS_GEOM_AUDIT");
    var_geom_audit=geom_audit_env && atoi(geom_audit_env)>0;
    if (force_call==1)
        fprintf(stderr,"[VARPRESS_GEOM_CFL] coefficient=%.17g audit=%d\n",var_geom_cfl,var_geom_audit);
#endif
#ifdef GFS_VARPRESS_W0
    if (!params->phase1 || params->av_mode != 5 || !params->use_muscl ||
        params->alphavis != 0.0 || params->betavis != 0.0 ||
        params->nu_phys != 0.0 || params->cd_amax != 0.0 ||
        params->hyperv_alpha != 0.0 || params->face_charge ||
        params->w2lin || params->laguerre_rot) {
        fprintf(stderr, "[VARPRESS_W0] requires phase1/MUSCL/av5, W0, Re=inf, no AV/weight rates\n");
        exit(1);
    }
    if (force_call == 1)
        fprintf(stderr, "[VARPRESS_W0] experimental smooth fixed-mass variational pressure active\n");
#endif
    if (dump_call < 0) {
        const char *s = getenv("GFS_GPU_FACE_DUMP_CALL");
        dump_call = (s && s[0]) ? atoi(s) : 0;
        const char *cursor=getenv("GFS_GPU_FACE_DUMP_IDS");
        while (cursor && *cursor) {
            while (*cursor==',' || *cursor==' ' || *cursor=='\t') cursor++;
            if (!*cursor) break;
            if (*cursor<'0' || *cursor>'9' || dump_nids>=64) {
                fprintf(stderr,"[FACE_INPUT] invalid/over-capacity ID list\n"); exit(1);
            }
            errno=0; char *end=NULL; unsigned long long id=strtoull(cursor,&end,10);
            if (errno || id>(unsigned long long)MAX_INDEX_GPU ||
                (*end && *end!=',' && *end!=' ' && *end!='\t')) {
                fprintf(stderr,"[FACE_INPUT] malformed ID\n"); exit(1);
            }
            dump_ids[dump_nids++]=(long long)id; cursor=end;
        }
        const char *first=getenv("GFS_GPU_FACE_DUMP_FIRST");
        const char *last=getenv("GFS_GPU_FACE_DUMP_LAST");
        if (first || last) {
            char *e1=NULL,*e2=NULL;
            long a=first?strtol(first,&e1,10):0, b=last?strtol(last,&e2,10):0;
            if (!first || !last || !first[0] || !last[0] || *e1 || *e2 ||
                a<1 || b<a || b-a>=16 || b>1000000000L || !dump_nids) {
                fprintf(stderr,"[FACE_INPUT] range needs positive FIRST/LAST, <=16 calls, and explicit IDs\n"); exit(1);
            }
            dump_first=(int)a; dump_last=(int)b;
        }
    }
    int n_faces = ctx->d_faces.n_faces_total;
    int do_dump = n_faces>0 && ((dump_call>0 && force_call==dump_call) ||
        (dump_first>0 && force_call>=dump_first && force_call<=dump_last));
    const char *shadow_env=getenv("GFS_ALE_SHADOW_CALL");
    int do_shadow=shadow_env && atoi(shadow_env)==force_call;
    double *d_audit_work = NULL, *d_audit_pressure = NULL, *d_audit_udotds = NULL;
    if (do_dump || do_shadow) {
        CUDA_CHECK(cudaMalloc(&d_audit_work, n_faces * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_audit_pressure, n_faces * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_audit_udotds, n_faces * sizeof(double)));
    }
    int blockSize = 256;
    int gridSize = (n_particles + blockSize - 1) / blockSize;
    int thermal64 = gfs_thermal64_requested();

    cudaStream_t s = (cudaStream_t)ctx->stream;
    {
        unsigned long long z = 0;
        CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_pair_hits, &z, sizeof(z), 0,
                                           cudaMemcpyHostToDevice, s));
        if (var_geom_audit)
            CUDA_CHECK(cudaMemcpyToSymbolAsync(d_varpress_geom_max_key,&z,sizeof(z),0,cudaMemcpyHostToDevice,s));
    }

    getAccVoro2DBlend_kernel<<<gridSize, blockSize, 0, s>>>(
        /* Face CSR */
        ctx->d_faces.face_offset,
        ctx->d_faces.c1x, ctx->d_faces.c1y,
        ctx->d_faces.c2x, ctx->d_faces.c2y,
        ctx->d_faces.neighbor_idx,
        ctx->d_faces.kp_idx, ctx->d_faces.km_idx,
        ctx->d_faces.is_ghost,
        /* Particle SoA */
        ctx->d_parts.x, ctx->d_parts.y,
        ctx->d_parts.vx, ctx->d_parts.vy,
        ctx->d_parts.mass,
        ctx->d_parts.den, ctx->d_parts.pressure,
        ctx->d_parts.csound, ctx->d_parts.volume,
        ctx->d_parts.ie, ctx->d_parts.w2, ctx->d_parts.w2old,
        /* Stress fields */
        ctx->d_parts.gUxx, ctx->d_parts.gUxy,
        ctx->d_parts.gUyx, ctx->d_parts.gUyy,
        ctx->d_parts.dPdx, ctx->d_parts.dPdy,
        ctx->d_parts.dRhodx, ctx->d_parts.dRhody,
        ctx->d_parts.tauxx, ctx->d_parts.tauxy, ctx->d_parts.tauyy,
        ctx->d_parts.alpha_cd, ctx->d_parts.divv,
        /* Hyperviscosity Pass 2 input */
        ctx->d_parts.lap_vx, ctx->d_parts.lap_vy,
        /* Output */
        ctx->d_parts.ax_out, ctx->d_parts.ay_out,
        ctx->d_parts.die_out, ctx->d_parts.dK_out,
        ctx->d_parts.dt_out,
        ctx->d_parts.vsig_max_out,
        ctx->d_parts.vsmoothx_out, ctx->d_parts.vsmoothy_out,
        d_audit_work, d_audit_pressure, d_audit_udotds,
        /* Physics */
        n_particles,
        params->av_mode, params->use_muscl,
        params->OoA, params->Courant, params->Gamma,
        params->alphavis, params->betavis,
        params->etavis, params->epsvis,
        params->nu_phys, params->prandtl,
        params->cd_amax, params->blend_theta, params->dtold,
        params->hyperv_alpha, params->hyperv_force_cap,
        params->phase1, params->laguerre_rot, params->face_charge, params->geom_fv,
        params->w2lin,var_geom_cfl,var_geom_audit,
        thermal64,ctx->d_parts.E_inv_xx,ctx->d_parts.E_inv_xy,ctx->d_parts.E_inv_yx);

    CUDA_CHECK(cudaGetLastError());

    /* CFL reduction using pre-allocated CUB scratch */
    cub::DeviceReduce::Min(ctx->d_cub_temp, ctx->cub_temp_bytes,
                           ctx->d_parts.dt_out,
                           (float *)ctx->d_cub_min_out,
                           n_particles, s);
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_cub_min_out, ctx->d_cub_min_out,
                                sizeof(float), cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaStreamSynchronize(s));
    {
        unsigned long long h = 0;
        CUDA_CHECK(cudaMemcpyFromSymbol(&h, d_gfs_pair_hits, sizeof(h)));
        h_gfs_pair_hits = (long long)h;
    }

    h_varpress_geom_max_rate=0; h_varpress_geom_max_id=-1;
    if (var_geom_audit) {
        unsigned long long key=0;
        CUDA_CHECK(cudaMemcpyFromSymbol(&key,d_varpress_geom_max_key,sizeof(key)));
        unsigned int bits=(unsigned int)(key>>32), index=(unsigned int)key;
        float rate=0; memcpy(&rate,&bits,sizeof(rate));
        h_varpress_geom_max_rate=(double)rate;
        if (index<(unsigned int)n_particles) h_varpress_geom_max_id=ctx->h_parts.indx[index];
    }

    if(do_shadow) gfs_ale_shadow_capture(ctx,n_particles,n_total,mpi_rank,force_call,params,d_audit_pressure,d_audit_work);

    if (do_dump) {
        size_t bytes = (size_t)n_faces * sizeof(double);
        double *h_work = (double *)malloc(bytes);
        double *h_pressure = (double *)malloc(bytes);
        double *h_udotds = (double *)malloc(bytes);
        CUDA_CHECK(cudaMemcpy(h_work, d_audit_work, bytes, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(h_pressure, d_audit_pressure, bytes, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(h_udotds, d_audit_udotds, bytes, cudaMemcpyDeviceToHost));
        gpu_download_face_csr(ctx, n_particles, n_faces,
            ctx->h_faces.face_offset,
            ctx->h_faces.c1x, ctx->h_faces.c1y,
            ctx->h_faces.c2x, ctx->h_faces.c2y,
            ctx->h_faces.neighbor_idx,
            ctx->h_faces.kp_idx, ctx->h_faces.km_idx,
            ctx->h_faces.is_ghost);
        CUDA_CHECK(cudaMemcpy(ctx->h_parts.indx, ctx->d_parts.indx,
            (size_t)n_total * sizeof(long long), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(ctx->h_parts.x, ctx->d_parts.x,
            (size_t)n_total * sizeof(double), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(ctx->h_parts.y, ctx->d_parts.y,
            (size_t)n_total * sizeof(double), cudaMemcpyDeviceToHost));
        gfs_dump_force_inputs(ctx,n_particles,n_total,mpi_rank,force_call,params,dump_ids,dump_nids);
        /* Face flux/geometry diagnostics do not require the much larger
         * all-candidate structural dump on a locally refined mesh. */
        const char *tess_dump = getenv("GFS_GPU_TESS_DUMP");
        if (!(tess_dump && atoi(tess_dump)==0))
            dump_kepler_tess_snapshot(ctx, n_particles, n_total, mpi_rank, force_call,
                                      h_diag_tess_mx, h_diag_tess_my,
                                      h_diag_tess_cellsize,
                                      h_diag_tess_xmin, h_diag_tess_ymin);

        char path[128];
        snprintf(path, sizeof(path), "face_dump_r%d_call%d.tsv", mpi_rank, force_call);
        FILE *out = fopen(path, "w");
        if (!out) {
            fprintf(stderr, "[FACE_DUMP] rank=%d could not open %s\n", mpi_rank, path);
            exit(1);
        }
        fprintf(out, "rank\tcall\ti\tj\tid_i\tid_j\twall\tpstar\tudotds\twork\tdSx\tdSy\tx1\ty1\tx2\ty2\n");
        long long nrows = 0;
        for (int i = 0; i < n_particles; i++) {
            if (!gfs_face_dump_selected(ctx->h_parts.indx[i],dump_ids,dump_nids)) continue;
            for (int f = ctx->h_faces.face_offset[i];
                 f < ctx->h_faces.face_offset[i + 1]; f++) {
                int j = ctx->h_faces.neighbor_idx[f];
                if (j < 0 || j >= n_total) continue;
                double x0 = ctx->h_parts.x[i], y0 = ctx->h_parts.y[i];
                double c1x = ctx->h_faces.c1x[f], c1y = ctx->h_faces.c1y[f];
                double c2x = ctx->h_faces.c2x[f], c2y = ctx->h_faces.c2y[f];
                fprintf(out, "%d\t%d\t%d\t%d\t%lld\t%lld\t%d\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
                    mpi_rank, force_call, i, j,
                    ctx->h_parts.indx[i], ctx->h_parts.indx[j],
                    ctx->h_faces.is_ghost[f],
                    h_pressure[f], h_udotds[f], h_work[f],
                    c2y-c1y, -(c2x-c1x),
                    x0+c1x, y0+c1y, x0+c2x, y0+c2y);
                nrows++;
            }
        }
        fclose(out);
        fprintf(stderr, "[FACE_DUMP] rank=%d call=%d rows=%lld file=%s\n",
            mpi_rank, force_call, nrows, path);
        free(h_work); free(h_pressure); free(h_udotds);
    }
    if(do_dump || do_shadow){
        CUDA_CHECK(cudaFree(d_audit_work));
        CUDA_CHECK(cudaFree(d_audit_pressure));
        CUDA_CHECK(cudaFree(d_audit_udotds));
    }

    return (double)(*(float *)ctx->h_cub_min_out);
}

/* ================================================================

/* ================================================================
 *  GPU Tessellation: Voronoi/Laguerre on GPU
 *
 *  Three kernels:
 *    A. voronoi_tessellate_kernel — Sutherland-Hodgman polygon clipping
 *    B. pressure_stress_kernel   — EOS + NS stress from gradients
 *    C. scatter_faces_kernel     — PolygonOut → FaceCSR on device
 *
 *  Orchestrated by gpu_tessellate_and_build_faces().
 * ================================================================ */

/* ---- Tessellation buffer management ---- */

extern "C"
void gpu_alloc_tess_buffers(GPUContext *ctx, int max_particles,
                            int max_cells, int max_cell_entries)
{
    cudaStream_t s = (cudaStream_t)ctx->stream;

    /* Cell CSR */
    ctx->max_cells = max_cells;
    ctx->max_cell_entries = max_cell_entries;
    size_t so = (max_cells + 1) * sizeof(int);
    size_t se = max_cell_entries * sizeof(int);
    CUDA_CHECK(cudaMallocHost(&ctx->h_cells.cell_offset, so));
    CUDA_CHECK(cudaMalloc(&ctx->d_cells.cell_offset, so));
    CUDA_CHECK(cudaMallocHost(&ctx->h_cells.particle_index, se));
    CUDA_CHECK(cudaMalloc(&ctx->d_cells.particle_index, se));

    /* PolygonOut per particle */
    CUDA_CHECK(cudaMalloc(&ctx->d_polygons,
               max_particles * sizeof(PolygonOut)));

    /* Face count and offset arrays */
    CUDA_CHECK(cudaMalloc(&ctx->d_face_count,  max_particles * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&ctx->d_face_offset, (max_particles + 1) * sizeof(int)));

    /* CUB scratch for ExclusiveScan */
    {
        void *d_temp = NULL;
        size_t temp_bytes = 0;
        cub::DeviceScan::ExclusiveSum(d_temp, temp_bytes,
                                       ctx->d_face_count,
                                       ctx->d_face_offset,
                                       max_particles, s);
        ctx->cub_scan_temp_bytes = temp_bytes;
        CUDA_CHECK(cudaMalloc(&ctx->d_cub_scan_temp, temp_bytes));
    }

    ctx->tess_faces_on_device = 0;

    fprintf(stderr, "[GPU] Tessellation buffers: max_cells=%d, max_entries=%d, "
            "PolygonOut=%.1f MB\n",
            max_cells, max_cell_entries,
            max_particles * sizeof(PolygonOut) / 1048576.0);
}

extern "C"
void gpu_free_tess_buffers(GPUContext *ctx)
{
    if (ctx->h_cells.cell_offset)    { cudaFreeHost(ctx->h_cells.cell_offset);    ctx->h_cells.cell_offset = NULL; }
    if (ctx->d_cells.cell_offset)    { cudaFree(ctx->d_cells.cell_offset);        ctx->d_cells.cell_offset = NULL; }
    if (ctx->h_cells.particle_index) { cudaFreeHost(ctx->h_cells.particle_index); ctx->h_cells.particle_index = NULL; }
    if (ctx->d_cells.particle_index) { cudaFree(ctx->d_cells.particle_index);     ctx->d_cells.particle_index = NULL; }
    if (ctx->d_polygons)    { cudaFree(ctx->d_polygons);    ctx->d_polygons = NULL; }
    if (ctx->d_face_count)  { cudaFree(ctx->d_face_count);  ctx->d_face_count = NULL; }
    if (ctx->d_face_offset) { cudaFree(ctx->d_face_offset); ctx->d_face_offset = NULL; }
    if (ctx->d_cub_scan_temp) { cudaFree(ctx->d_cub_scan_temp); ctx->d_cub_scan_temp = NULL; }
    ctx->cub_scan_temp_bytes = 0;
    ctx->max_cells = 0;
    ctx->max_cell_entries = 0;
}

/* ================================================================
 *  Option C: OST k-d tree NN search — buffer management
 * ================================================================ */
extern "C"
void gpu_alloc_tree_buffers(GPUContext *ctx, int max_tree_nodes)
{
    if (max_tree_nodes <= ctx->max_tree_nodes) return;  /* already large enough */

    /* Free previous allocation if any (idempotent realloc) */
    if (ctx->d_tree_type) { cudaFree(ctx->d_tree_type);     ctx->d_tree_type = NULL; }
    if (ctx->d_tree_sibling)  { cudaFree(ctx->d_tree_sibling);  ctx->d_tree_sibling = NULL; }
    if (ctx->d_tree_daughter) { cudaFree(ctx->d_tree_daughter); ctx->d_tree_daughter = NULL; }
    if (ctx->d_tree_ptl_idx)  { cudaFree(ctx->d_tree_ptl_idx);  ctx->d_tree_ptl_idx = NULL; }
    if (ctx->d_tree_cx)       { cudaFree(ctx->d_tree_cx);       ctx->d_tree_cx = NULL; }
    if (ctx->d_tree_cy)       { cudaFree(ctx->d_tree_cy);       ctx->d_tree_cy = NULL; }
    if (ctx->d_tree_size)     { cudaFree(ctx->d_tree_size);     ctx->d_tree_size = NULL; }
    if (ctx->h_tree_type)     { cudaFreeHost(ctx->h_tree_type);     ctx->h_tree_type = NULL; }
    if (ctx->h_tree_sibling)  { cudaFreeHost(ctx->h_tree_sibling);  ctx->h_tree_sibling = NULL; }
    if (ctx->h_tree_daughter) { cudaFreeHost(ctx->h_tree_daughter); ctx->h_tree_daughter = NULL; }
    if (ctx->h_tree_ptl_idx)  { cudaFreeHost(ctx->h_tree_ptl_idx);  ctx->h_tree_ptl_idx = NULL; }
    if (ctx->h_tree_cx)       { cudaFreeHost(ctx->h_tree_cx);       ctx->h_tree_cx = NULL; }
    if (ctx->h_tree_cy)       { cudaFreeHost(ctx->h_tree_cy);       ctx->h_tree_cy = NULL; }
    if (ctx->h_tree_size)     { cudaFreeHost(ctx->h_tree_size);     ctx->h_tree_size = NULL; }

    size_t si = max_tree_nodes * sizeof(int);
    size_t sf = max_tree_nodes * sizeof(float);

    CUDA_CHECK(cudaMalloc(&ctx->d_tree_type,     si));
    CUDA_CHECK(cudaMalloc(&ctx->d_tree_sibling,  si));
    CUDA_CHECK(cudaMalloc(&ctx->d_tree_daughter, si));
    CUDA_CHECK(cudaMalloc(&ctx->d_tree_ptl_idx,  si));
    CUDA_CHECK(cudaMalloc(&ctx->d_tree_cx,       sf));
    CUDA_CHECK(cudaMalloc(&ctx->d_tree_cy,       sf));
    CUDA_CHECK(cudaMalloc(&ctx->d_tree_size,     sf));

    CUDA_CHECK(cudaMallocHost(&ctx->h_tree_type,     si));
    CUDA_CHECK(cudaMallocHost(&ctx->h_tree_sibling,  si));
    CUDA_CHECK(cudaMallocHost(&ctx->h_tree_daughter, si));
    CUDA_CHECK(cudaMallocHost(&ctx->h_tree_ptl_idx,  si));
    CUDA_CHECK(cudaMallocHost(&ctx->h_tree_cx,       sf));
    CUDA_CHECK(cudaMallocHost(&ctx->h_tree_cy,       sf));
    CUDA_CHECK(cudaMallocHost(&ctx->h_tree_size,     sf));

    ctx->max_tree_nodes = max_tree_nodes;
    fprintf(stderr, "[GPU] Tree buffers allocated: %d nodes (%.1f MB total)\n",
            max_tree_nodes, max_tree_nodes * (4*sizeof(int)+3*sizeof(float)) / 1048576.0);
}

extern "C"
void gpu_free_tree_buffers(GPUContext *ctx)
{
    if (ctx->d_tree_type)     { cudaFree(ctx->d_tree_type);     ctx->d_tree_type = NULL; }
    if (ctx->d_tree_sibling)  { cudaFree(ctx->d_tree_sibling);  ctx->d_tree_sibling = NULL; }
    if (ctx->d_tree_daughter) { cudaFree(ctx->d_tree_daughter); ctx->d_tree_daughter = NULL; }
    if (ctx->d_tree_ptl_idx)  { cudaFree(ctx->d_tree_ptl_idx);  ctx->d_tree_ptl_idx = NULL; }
    if (ctx->d_tree_cx)       { cudaFree(ctx->d_tree_cx);       ctx->d_tree_cx = NULL; }
    if (ctx->d_tree_cy)       { cudaFree(ctx->d_tree_cy);       ctx->d_tree_cy = NULL; }
    if (ctx->d_tree_size)     { cudaFree(ctx->d_tree_size);     ctx->d_tree_size = NULL; }
    if (ctx->h_tree_type)     { cudaFreeHost(ctx->h_tree_type);     ctx->h_tree_type = NULL; }
    if (ctx->h_tree_sibling)  { cudaFreeHost(ctx->h_tree_sibling);  ctx->h_tree_sibling = NULL; }
    if (ctx->h_tree_daughter) { cudaFreeHost(ctx->h_tree_daughter); ctx->h_tree_daughter = NULL; }
    if (ctx->h_tree_ptl_idx)  { cudaFreeHost(ctx->h_tree_ptl_idx);  ctx->h_tree_ptl_idx = NULL; }
    if (ctx->h_tree_cx)       { cudaFreeHost(ctx->h_tree_cx);       ctx->h_tree_cx = NULL; }
    if (ctx->h_tree_cy)       { cudaFreeHost(ctx->h_tree_cy);       ctx->h_tree_cy = NULL; }
    if (ctx->h_tree_size)     { cudaFreeHost(ctx->h_tree_size);     ctx->h_tree_size = NULL; }
    ctx->max_tree_nodes = 0;
    ctx->n_tree_nodes = 0;
}

extern "C"
void gpu_upload_tree(GPUContext *ctx, int n_nodes)
{
    cudaStream_t s = (cudaStream_t)ctx->stream;
    size_t si = n_nodes * sizeof(int);
    size_t sf = n_nodes * sizeof(float);

    CUDA_CHECK(cudaMemcpyAsync(ctx->d_tree_type,     ctx->h_tree_type,     si, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_tree_sibling,  ctx->h_tree_sibling,  si, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_tree_daughter, ctx->h_tree_daughter, si, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_tree_ptl_idx,  ctx->h_tree_ptl_idx,  si, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_tree_cx,       ctx->h_tree_cx,       sf, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_tree_cy,       ctx->h_tree_cy,       sf, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_tree_size,     ctx->h_tree_size,     sf, cudaMemcpyHostToDevice, s));

    ctx->n_tree_nodes = n_nodes;
}

extern "C"
void gpu_upload_cells(GPUContext *ctx, int n_cells, int n_entries)
{
    cudaStream_t s = (cudaStream_t)ctx->stream;
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_cells.cell_offset, ctx->h_cells.cell_offset,
               (n_cells + 1) * sizeof(int), cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_cells.particle_index, ctx->h_cells.particle_index,
               n_entries * sizeof(int), cudaMemcpyHostToDevice, s));
    ctx->d_cells.n_cells = n_cells;
    ctx->d_cells.n_entries = n_entries;
}

/* ================================================================
 *  Kernel A: voronoi_tessellate_kernel
 *
 *  1 thread per real particle. Builds Voronoi/Laguerre polygon via
 *  Sutherland-Hodgman clipping. Computes volume, density, Green-Gauss
 *  gradients, avgNeighP, w2ceil. Writes PolygonOut for face scatter.
 * ================================================================ */

/* Legacy inline capacity, retained only to count dense seed stencils.
 * Candidate clipping now streams directly from CellCSR; no candidate is
 * truncated and no per-thread arrays of this size are allocated. */
#define GPU_MAX_NEIGH 512
#ifdef GFS_FAST_SEED3
#define GPU_SEED_PAD 1
#else
#define GPU_SEED_PAD 2
#endif

/* Remove only exactly coincident consecutive vertices. The outgoing
 * support belongs to the last emitted copy of that point. Keeping both
 * copies can exhaust the polygon capacity during dense grid clipping,
 * even when the actual convex cell has only a few nonzero edges. */
__device__ __forceinline__
int gpu_raw_append_vertex(double x, double y, int support,
    double *out_x, double *out_y, int *out_support, int *count)
{
    int n=*count;
    if(n>0 && out_x[n-1]==x && out_y[n-1]==y){
        out_support[n-1]=support;
        return 1;
    }
    if(n>=GPU_MAX_CORNERS) return 0;
    out_x[n]=x;out_y[n]=y;out_support[n]=support;*count=n+1;
    return 1;
}

__device__ __forceinline__
int gpu_clip_voronoi_halfplane(
    double nx_n, double ny_n, double wfrac, int neighbor_soa,
    double *&poly_x, double *&poly_y, int *&poly_neigh, int nc,
    double *&tmp_x, double *&tmp_y, int *&tmp_neigh, int exact_intersections)
{
    const double dist2 = nx_n * nx_n + ny_n * ny_n;
    const double threshold = wfrac * dist2;
    int out_count = 0;

    if(exact_intersections){
        int valid=1;
        for(int e=0;e<nc && valid;e++){
            int en=(e+1)%nc;
            double ax=poly_x[e],ay=poly_y[e],bx=poly_x[en],by=poly_y[en];
            double da=ax*nx_n+ay*ny_n,db=bx*nx_n+by*ny_n;
            int ain=da<=threshold,bin=db<=threshold;
            if(ain)
                valid=gpu_raw_append_vertex(ax,ay,poly_neigh[e],tmp_x,tmp_y,tmp_neigh,&out_count);
            if(ain!=bin && valid){
                double denom=db-da;
                if(denom==0 || !isfinite(denom)){
                    atomicAdd(&d_gfs_tess_geometry_errors,1ull);
                    return -1;
                }
                double t=fmax(0.0,fmin(1.0,(threshold-da)/denom));
                double x=ax+t*(bx-ax),y=ay+t*(by-ay);
                valid=gpu_raw_append_vertex(x,y,ain?neighbor_soa:poly_neigh[e],
                    tmp_x,tmp_y,tmp_neigh,&out_count);
            }
        }
        if(!valid){
            unsigned long long errors=atomicAdd(&d_gfs_tess_geometry_errors,1ull);
            if(errors<8)
                printf("[GPU_CLIP_CAPACITY] nc=%d capacity=%d neighbor_soa=%d normal=(%.17g,%.17g) threshold=%.17g\n",nc,GPU_MAX_CORNERS,neighbor_soa,nx_n,ny_n,threshold);
            return -1;
        }
        if(out_count>1 && tmp_x[0]==tmp_x[out_count-1] && tmp_y[0]==tmp_y[out_count-1])
            out_count--;
        /* The output polygon is already complete in the scratch buffers.
         * Swap the views instead of copying three arrays back after every
         * accepted halfplane. The clipping order and arithmetic stay exact. */
        double *swap_x=poly_x; poly_x=tmp_x; tmp_x=swap_x;
        double *swap_y=poly_y; poly_y=tmp_y; tmp_y=swap_y;
        int *swap_neigh=poly_neigh; poly_neigh=tmp_neigh; tmp_neigh=swap_neigh;
        return out_count;
    }

    for (int e = 0; e < nc; e++) {
        const int e_next = (e + 1) % nc;
        const double ax = poly_x[e], ay = poly_y[e];
        const double bx = poly_x[e_next], by = poly_y[e_next];
        const double da = ax * nx_n + ay * ny_n;
        const double db = bx * nx_n + by * ny_n;
        const int a_in = (da <= threshold);
        const int b_in = (db <= threshold);

        if (a_in && b_in) {
            if (out_count < GPU_MAX_CORNERS) {
                tmp_x[out_count] = ax;
                tmp_y[out_count] = ay;
                tmp_neigh[out_count] = poly_neigh[e];
                out_count++;
            } else atomicAdd(&d_gfs_tess_geometry_errors,1ull);
        } else if (a_in && !b_in) {
            if (out_count < GPU_MAX_CORNERS) {
                tmp_x[out_count] = ax;
                tmp_y[out_count] = ay;
                tmp_neigh[out_count] = poly_neigh[e];
                out_count++;
            } else atomicAdd(&d_gfs_tess_geometry_errors,1ull);
            const double denom = db - da;
            if (exact_intersections && (denom == 0.0 || !isfinite(denom)))
                atomicAdd(&d_gfs_tess_geometry_errors, 1ull);
            double t = exact_intersections ? ((denom != 0.0) ? (threshold-da)/denom : 0.0)
                       : ((fabs(denom) > 1e-30) ? (threshold-da)/denom : wfrac);
            t = exact_intersections ? fmax(0.0,fmin(1.0,t)) : fmax(1e-9,fmin(1.0-1e-9,t));
            if (out_count < GPU_MAX_CORNERS) {
                tmp_x[out_count] = ax + t * (bx - ax);
                tmp_y[out_count] = ay + t * (by - ay);
                tmp_neigh[out_count] = neighbor_soa;
                out_count++;
            } else atomicAdd(&d_gfs_tess_geometry_errors,1ull);
        } else if (!a_in && b_in) {
            const double denom = db - da;
            if (exact_intersections && (denom == 0.0 || !isfinite(denom)))
                atomicAdd(&d_gfs_tess_geometry_errors, 1ull);
            double t = exact_intersections ? ((denom != 0.0) ? (threshold-da)/denom : 0.0)
                       : ((fabs(denom) > 1e-30) ? (threshold-da)/denom : wfrac);
            t = exact_intersections ? fmax(0.0,fmin(1.0,t)) : fmax(1e-9,fmin(1.0-1e-9,t));
            if (out_count < GPU_MAX_CORNERS) {
                tmp_x[out_count] = ax + t * (bx - ax);
                tmp_y[out_count] = ay + t * (by - ay);
                tmp_neigh[out_count] = poly_neigh[e];
                out_count++;
            } else atomicAdd(&d_gfs_tess_geometry_errors,1ull);
        }
    }

    for (int e = 0; e < out_count; e++) {
        poly_x[e] = tmp_x[e];
        poly_y[e] = tmp_y[e];
        poly_neigh[e] = tmp_neigh[e];
    }
    return out_count;
}

/* Exact ordinary-Voronoi fast path: equal zero power weights place the
 * bisector halfway between the generators, so avoid a double division.
 * Nonzero or unequal weights retain the original expression and clamp. */
__device__ __forceinline__
double gpu_laguerre_wfrac_minmax(double w2i, float w2j, double dist2)
{
    if (w2i == 0.0 && w2j == 0.0f) return 0.5;
    return fmax(0.05, fmin(0.95,
        0.5 + 0.5 * (w2i - (double)w2j) / dist2));
}

__device__ __forceinline__
double gpu_laguerre_wfrac_ifclamp(double w2i, float w2j, double dist2)
{
    if (w2i == 0.0 && w2j == 0.0f) return 0.5;
    double wfrac = 0.5 + 0.5 * (w2i - (double)w2j) / dist2;
    if (wfrac < 0.05) wfrac = 0.05;
    else if (wfrac > 0.95) wfrac = 0.95;
    return wfrac;
}

__global__ __launch_bounds__(256, 2)
void voronoi_tessellate_kernel(
    /* Particle SoA */
    const double *__restrict__ px, const double *__restrict__ py,
    const double *__restrict__ pvx, const double *__restrict__ pvy,
    const double *__restrict__ pmass,
    const long long *__restrict__ pindx,
    const float *__restrict__ pie, const float *__restrict__ pw2,
    const float *__restrict__ ppressure_in,
    const float *__restrict__ pden_in,
    /* Cell CSR */
    const int *__restrict__ cell_offset,
    const int *__restrict__ cell_particle_index,
    /* Output: tessellation results in SoA */
    float *__restrict__ out_den, float *__restrict__ out_pressure,
    float *__restrict__ out_volume, float *__restrict__ out_w2ceil,
    float *__restrict__ out_avgNeighP,
    /* Output: gradients */
    double *__restrict__ out_gUxx, double *__restrict__ out_gUxy,
    double *__restrict__ out_gUyx, double *__restrict__ out_gUyy,
    double *__restrict__ out_dPdx, double *__restrict__ out_dPdy,
    double *__restrict__ out_dRhodx, double *__restrict__ out_dRhody,
    /* Output: hyperviscosity Pass 1 (velocity Laplacian) */
    double *__restrict__ out_lap_vx, double *__restrict__ out_lap_vy,
    /* Output: polygon for face scatter */
    PolygonOut *__restrict__ polygons,
    int *__restrict__ face_count_out,
    /* Grid parameters */
    int nbp, int n_total, int mx, int my,
    double cellsize, double xmin, double ymin,
    double boxsize, double global_w2_max, int all_weights_zero,
    double OoA, int av_mode,
    int gradient_method, int raw_volume, int limit_velocity, int merge_small_faces,
    int thermal64, double *__restrict__ out_V64)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nbp) return;

    double xi = px[i], yi = py[i];
    double w2i = all_weights_zero ? 0.0 : (double)pw2[i];
    long long indx_i = pindx[i];
    double massi = pmass[i];

    /* Determine cell (ix, iy) from position */
    int ix = (int)((xi - xmin) / cellsize);
    int iy = (int)((yi - ymin) / cellsize);
    if (ix < 0) ix = 0;  if (ix >= mx) ix = mx - 1;
    if (iy < 0) iy = 0;  if (iy >= my) iy = my - 1;

    /* ---- Initialize box polygon ---- */
    double max_halfbox = 0.5 * boxsize;
    double halfbox = raw_volume >= 2 ? fmin(max_halfbox,cellsize) : max_halfbox;
    /* 4 box corners: (-hb,-hb), (+hb,-hb), (+hb,+hb), (-hb,+hb) */
    double poly_x_storage[GPU_MAX_CORNERS], poly_y_storage[GPU_MAX_CORNERS];
    int    poly_neigh_storage[GPU_MAX_CORNERS];
    int    nc = 4, stencil_pad = GPU_SEED_PAD, seed_candidates = 0, box_retry = 0;
    double maxdist2;
    double tmp_x_storage[GPU_MAX_CORNERS], tmp_y_storage[GPU_MAX_CORNERS];
    int    tmp_neigh_storage[GPU_MAX_CORNERS];
    double *poly_x=poly_x_storage,*poly_y=poly_y_storage;
    int *poly_neigh=poly_neigh_storage;
    double *tmp_x=tmp_x_storage,*tmp_y=tmp_y_storage;
    int *tmp_neigh=tmp_neigh_storage;

    /* Start locally to avoid cancellation when a cell of width~h is clipped
     * out of a nominal box thousands of times wider. A local box is accepted
     * only when no artificial box support remains. If it remains, double
     * the box and repeat the complete search. Convexity then certifies that
     * the accepted polygon is not truncated by the initial box. */
    for (;;) {
    nc=4;
    /* Same init as Voro2D_InitializeCorner: corner i = (((i+1)%4)/2, i/2)*box - halfbox */
    poly_x[0] = -halfbox; poly_y[0] = -halfbox; poly_neigh[0] = -1;
    poly_x[1] =  halfbox; poly_y[1] = -halfbox; poly_neigh[1] = -2;
    poly_x[2] =  halfbox; poly_y[2] =  halfbox; poly_neigh[2] = -3;
    poly_x[3] = -halfbox; poly_y[3] =  halfbox; poly_neigh[3] = -4;

    maxdist2 = 2.0 * halfbox * halfbox;

    /* ---- Sutherland-Hodgman clipping ---- */

    if(raw_volume>=2){
        /* Dense local refinement can give a distant row dozens of active
         * intermediate supports before the nearest generators are visited.
         * Tighten the starting superset with up to eight real halfplanes
         * first. This changes only processing order: every seed/adaptive
         * candidate is still examined below, and the corner cap still
         * aborts rather than truncating a genuinely complex polygon. */
        int nearest[8]; double offset2[8];
        for(int s=0;s<8;s++){nearest[s]=-1;offset2[s]=DBL_MAX;}
        for(int yc=-GPU_SEED_PAD;yc<=GPU_SEED_PAD;yc++){
            int cy=iy+yc;if(cy<0 || cy>=my)continue;
            for(int xc=-GPU_SEED_PAD;xc<=GPU_SEED_PAD;xc++){
                int cx=ix+xc;if(cx<0 || cx>=mx)continue;
                int cell=cx+mx*cy;
                for(int k=cell_offset[cell];k<cell_offset[cell+1];k++){
                    int j=cell_particle_index[k];if(pindx[j]==indx_i)continue;
                    double dx=px[j]-xi,dy=py[j]-yi,d2=dx*dx+dy*dy;
                    if(!(d2>0))continue;
                    double wf=all_weights_zero ? 0.5
                        : gpu_laguerre_wfrac_minmax(w2i,pw2[j],d2);
                    double dr2=wf*wf*d2;
                    int s=fabs(dx)>=fabs(dy)
                        ? (dx>=0?(dy>=0?0:7):(dy>=0?3:4))
                        : (dy>=0?(dx>=0?1:2):(dx>=0?6:5));
                    if(dr2<offset2[s]){offset2[s]=dr2;nearest[s]=j;}
                }
            }
        }
        for(int s=0;s<8 && nc>=3;s++){
            int j=nearest[s];if(j<0 || offset2[s]>maxdist2)continue;
            double dx=px[j]-xi,dy=py[j]-yi,d2=dx*dx+dy*dy;
            double wf=all_weights_zero ? 0.5
                : gpu_laguerre_wfrac_minmax(w2i,pw2[j],d2);
            nc=gpu_clip_voronoi_halfplane(dx,dy,wf,j,poly_x,poly_y,poly_neigh,nc,
                                         tmp_x,tmp_y,tmp_neigh,1);
            maxdist2=0;
            for(int e=0;e<nc;e++)maxdist2=fmax(maxdist2,poly_x[e]*poly_x[e]+poly_y[e]*poly_y[e]);
        }
    }

    /* Stream the selected seed in linked-cell order. The old
     * stored prefix and overflow pass together used this same sequence.
     * Keeping only the current half-plane removes 18 KiB of per-thread
     * candidate arrays, while retaining the polygon-radius rejection and
     * the complete adaptive search below. */
    seed_candidates = 0;
    for (int dy_c = -GPU_SEED_PAD; dy_c <= GPU_SEED_PAD && nc >= 3; dy_c++) {
        int cy = iy + dy_c;
        if (cy < 0 || cy >= my) continue;
        for (int dx_c = -GPU_SEED_PAD; dx_c <= GPU_SEED_PAD && nc >= 3; dx_c++) {
            int cx = ix + dx_c;
            if (cx < 0 || cx >= mx) continue;
            int cell_id = cx + mx * cy;
            int c_begin = cell_offset[cell_id];
            int c_end = cell_offset[cell_id + 1];
            for (int k = c_begin; k < c_end && nc >= 3; k++) {
                int j = cell_particle_index[k];
                if (pindx[j] == indx_i) continue;
                seed_candidates++;

                double dx = px[j] - xi;
                double dy_val = py[j] - yi;
                double dist2 = dx * dx + dy_val * dy_val;
                double wfrac = all_weights_zero ? 0.5
                    : gpu_laguerre_wfrac_ifclamp(w2i,pw2[j],dist2);
                double drad2 = wfrac * wfrac * dist2;
                if (drad2 > maxdist2) continue;

                nc = gpu_clip_voronoi_halfplane(dx, dy_val, wfrac, j,
                                                poly_x, poly_y, poly_neigh, nc,
                                                tmp_x, tmp_y, tmp_neigh, raw_volume >= 2);
                if (nc < 3) break;

                double new_max = 0.0;
                for (int e = 0; e < nc; e++) {
                    double cd = poly_x[e] * poly_x[e] + poly_y[e] * poly_y[e];
                    if (cd > new_max) new_max = cd;
                }
                maxdist2 = new_max;
            }
        }
    }

    /* Expand beyond the seed stencil when the current polygon can still
     * be cut by a more distant generator.  The seed polygon is a superset of
     * the complete cell, so its circumradius gives a safe search bound.  For
     * a clamped Laguerre bisector, the global weight bound gives the exact
     * conservative Euclidean radius needed to include every possible cut. */
    stencil_pad = GPU_SEED_PAD;
    if (nc >= 3) {
        const double radius = sqrt(fmax(0.0, maxdist2));
        const double deltaW = fmax(0.0, global_w2_max - w2i);
        const double d_clamped = (deltaW > 0.0)
            ? fmin(sqrt(deltaW / 0.9), 20.0 * radius) : 0.0;
        const double d_unclamped = radius + sqrt(maxdist2 + deltaW);
        double search_radius = fmax(2.0 * radius,
                                    fmax(d_clamped, d_unclamped));
        if (!(search_radius >= 0.0) || !isfinite(search_radius))
            search_radius = hypot((double)mx * cellsize,
                                  (double)my * cellsize);
#ifdef GFS_FAST_SEED3
        /* The 3x3 clipped polygon is a superset. All omitted bins lie at
         * least this far from the owner, even when their contents are
         * unknown. A strict radius gap certifies completeness. */
        const double bin_x0 = xmin + (double)ix * cellsize;
        const double bin_y0 = ymin + (double)iy * cellsize;
        const double min_outside_distance = fmin(
            fmin(xi - (bin_x0 - cellsize), (bin_x0 + 2.0*cellsize) - xi),
            fmin(yi - (bin_y0 - cellsize), (bin_y0 + 2.0*cellsize) - yi));
        if (search_radius < min_outside_distance) {
            stencil_pad = 1;
        } else {
#endif
        stencil_pad = (int)ceil(search_radius / cellsize) + 1;
        if (stencil_pad < GPU_SEED_PAD) stencil_pad = GPU_SEED_PAD;
        if (stencil_pad > (mx > my ? mx : my))
            stencil_pad = mx > my ? mx : my;
#ifdef GFS_FAST_SEED3
        }
#endif

        if (stencil_pad > GPU_SEED_PAD) {
            int ylo = iy - stencil_pad, yhi = iy + stencil_pad;
            int xlo = ix - stencil_pad, xhi = ix + stencil_pad;
            if (ylo < 0) ylo = 0; if (yhi >= my) yhi = my - 1;
            if (xlo < 0) xlo = 0; if (xhi >= mx) xhi = mx - 1;
            for (int cy = ylo; cy <= yhi && nc >= 3; cy++) {
                for (int cx = xlo; cx <= xhi && nc >= 3; cx++) {
                    if (abs(cx - ix) <= GPU_SEED_PAD && abs(cy - iy) <= GPU_SEED_PAD) continue;
                    const int cell_id = cx + mx * cy;
                    const int c_begin = cell_offset[cell_id];
                    const int c_end = cell_offset[cell_id + 1];
                    for (int k = c_begin; k < c_end && nc >= 3; k++) {
                        const int j = cell_particle_index[k];
                        if (pindx[j] == indx_i) continue;
                        const double dx = px[j] - xi;
                        const double dy_val = py[j] - yi;
                        const double dist2 = dx * dx + dy_val * dy_val;
                        if (!(dist2 > 0.0)) continue;
                        double wfrac = all_weights_zero ? 0.5
                            : gpu_laguerre_wfrac_ifclamp(w2i,pw2[j],dist2);
                        if (wfrac * wfrac * dist2 > maxdist2) continue;

                        nc = gpu_clip_voronoi_halfplane(dx, dy_val, wfrac, j,
                                                        poly_x, poly_y, poly_neigh, nc,
                                                        tmp_x, tmp_y, tmp_neigh, raw_volume >= 2);
                        if (nc < 3) break;
                        double new_max = 0.0;
                        for (int e = 0; e < nc; e++) {
                            const double cd = poly_x[e] * poly_x[e]
                                            + poly_y[e] * poly_y[e];
                            if (cd > new_max) new_max = cd;
                        }
                        maxdist2 = new_max;
                    }
                }
            }
        }
    }

    if (raw_volume < 2 || nc < 3) break;
    int touches_box=0;
    for (int e=0;e<nc;e++)
        if (poly_neigh[e]<0 || fabs(poly_x[e])>=halfbox*(1-32*DBL_EPSILON) ||
            fabs(poly_y[e])>=halfbox*(1-32*DBL_EPSILON)) touches_box=1;
    if (!touches_box) break;
    if (halfbox>=max_halfbox || ++box_retry>32) {
        atomicAdd(&d_gfs_tess_geometry_errors,1ull);
        face_count_out[i]=0;
        return;
    }
    halfbox=fmin(max_halfbox,2.0*halfbox);
    } /* adaptive initial box */
    if (seed_candidates > GPU_MAX_NEIGH)
        atomicAdd(&d_gfs_candidate_fallback_cells, 1ull);
    if (stencil_pad > GPU_SEED_PAD)
        atomicAdd(&d_gfs_stencil_expanded_cells, 1ull);

    if(raw_volume>=2){
        /* A convex power-cell has one connected edge per support plane.
         * Repeated clipping can leave consecutive collinear segments with
         * the same outgoing support. They must be one Riemann face: using
         * different segment midpoints at the two cell ends gives different
         * nonlinear pressures even when summed area vectors agree. Remove
         * the redundant vertex by support identity, with no length cutoff. */
        int changed=1;
        while(changed && nc>=3){
            changed=0;
            for(int e=0;e<nc;e++){
                int prev=(e+nc-1)%nc;
                if(poly_neigh[e]>=0 && poly_neigh[e]==poly_neigh[prev]){
                    for(int k=e;k+1<nc;k++){
                        poly_x[k]=poly_x[k+1];poly_y[k]=poly_y[k+1];poly_neigh[k]=poly_neigh[k+1];
                    }
                    nc--;changed=1;break;
                }
            }
        }
    }

    /* ---- Merge near-coincident vertices (degenerate faces at grid corners) ----
     * On a regular/near-regular grid, 4+ Voronoi cells meet at (near) a
     * single point. Sutherland-Hodgman clipping produces two nearly-
     * identical vertices, creating a tiny face.  These cause CFL dt
     * collapse.  Use RELATIVE threshold: edge < 1e-4 * max_edge. */
    if (!raw_volume && merge_small_faces) {
        double max_edge2 = 0;
        for (int e = 0; e < nc; e++) {
            int en = (e + 1) % nc;
            double dx_m = poly_x[en] - poly_x[e];
            double dy_m = poly_y[en] - poly_y[e];
            double e2 = dx_m * dx_m + dy_m * dy_m;
            if (e2 > max_edge2) max_edge2 = e2;
        }
        double merge_eps2 = max_edge2 * 1e-6;  /* (1e-3 * max_edge)^2 */
        if (merge_eps2 < 1e-30) merge_eps2 = 1e-30;

        int new_nc = 0;
        for (int e = 0; e < nc; e++) {
            int en = (e + 1) % nc;
            double dx_m = poly_x[en] - poly_x[e];
            double dy_m = poly_y[en] - poly_y[e];
            if (dx_m * dx_m + dy_m * dy_m > merge_eps2) {
                tmp_x[new_nc] = poly_x[e];
                tmp_y[new_nc] = poly_y[e];
                tmp_neigh[new_nc] = poly_neigh[e];
                new_nc++;
            }
        }
        if (new_nc >= 3 && new_nc < nc) {
            nc = new_nc;
            for (int e = 0; e < nc; e++) {
                poly_x[e] = tmp_x[e];
                poly_y[e] = tmp_y[e];
                poly_neigh[e] = tmp_neigh[e];
            }
        }
    }

    /* ---- Compute volume via shoelace formula ---- */
    double area = 0;
    for (int e = 0; e < nc; e++) {
        int en = (e + 1) % nc;
        area += poly_x[e] * poly_y[en] - poly_y[e] * poly_x[en];
    }
    double volume = 0.5 * fabs(area);
    if (raw_volume >= 2) {
        /* Each positive edge must remain on its recorded support plane.
         * This catches mislabeled/backward/ill-resolved micro-faces rather
         * than treating positive shoelace area as a topology certificate. */
        for (int e=0;e<nc;e++) {
            int jn=poly_neigh[e], en=(e+1)%nc;
            if (jn < 0) continue;
            double ex=poly_x[en]-poly_x[e], ey=poly_y[en]-poly_y[e];
            if (!(ex*ex+ey*ey > 0.0)) continue;
            double dx=px[jn]-xi, dy=py[jn]-yi, d2=dx*dx+dy*dy;
            double wf=all_weights_zero ? 0.5
                : gpu_laguerre_wfrac_minmax(w2i,pw2[jn],d2);
            double rhs=wf*d2;
            double scale=fabs(rhs)+fabs(dx)*(fabs(poly_x[e])+fabs(poly_x[en]))+
                         fabs(dy)*(fabs(poly_y[e])+fabs(poly_y[en]));
            double tol=1024.0*DBL_EPSILON*fmax(scale,DBL_MIN);
#ifdef GFS_VARPRESS_W0
            /* Clipping/intersection roundoff is measured in the cell's
             * Cartesian coordinates. Include that construction scale even
             * when a nearly axis-aligned support has a tiny dot product. */
            double vertex_scale=fmax(hypot(poly_x[e],poly_y[e]),
                                     hypot(poly_x[en],poly_y[en]));
            double construction_scale=hypot(dx,dy)*fmax(vertex_scale,sqrt(volume));
            tol=fmax(tol,1024.0*DBL_EPSILON*construction_scale);
#endif
            double err=fmax(fabs(dx*poly_x[e]+dy*poly_y[e]-rhs),
                            fabs(dx*poly_x[en]+dy*poly_y[en]-rhs));
            /* A zero-measure Voronoi edge at a cocircular grid corner can
             * reverse by one coordinate ulp. Bound subtraction/product
             * roundoff explicitly; retain every vertex and its closure.
             * A resolved reverse edge remains an error. */
            double orientation=ey*dx-ex*dy;
            double orientation_error=8.0*DBL_EPSILON*(
                fabs(dx)*(fabs(poly_y[e])+fabs(poly_y[en]))+
                fabs(dy)*(fabs(poly_x[e])+fabs(poly_x[en]))+
                fabs(ey*dx)+fabs(ex*dy));
#ifdef GFS_VARPRESS_W0
            /* A face shorter than local coordinate roundoff has no
             * resolvable orientation. Keep its vertices and support check;
             * only exempt the sign of this zero-measure edge. */
            double local_coord=fmax(fmax(fabs(poly_x[e]),fabs(poly_y[e])),
                                    fmax(fabs(poly_x[en]),fabs(poly_y[en])));
            double unresolved_length=32.0*DBL_EPSILON*fmax(local_coord,sqrt(volume));
            bool unresolved_orientation=hypot(ex,ey)<=unresolved_length;
#else
            bool unresolved_orientation=false;
#endif
            if (!isfinite(err) || err>tol ||
                (orientation < -orientation_error && !unresolved_orientation)) {
                unsigned long long count=atomicAdd(&d_gfs_tess_geometry_errors,1ull);
                if (count<8)
                    printf("[GPU_SUPPORT_ERROR] i=%d id=%lld edge=%d neighbor=%lld err=%.17g tol=%.17g orientation=%.17g edge=(%.17g,%.17g) n=(%.17g,%.17g) a=(%.17g,%.17g) b=(%.17g,%.17g) area=%.17g\n",
                        i,(long long)pindx[i],e,(long long)pindx[jn],err,tol,ey*dx-ex*dy,ex,ey,dx,dy,
                        poly_x[e],poly_y[e],poly_x[en],poly_y[en],area);
                face_count_out[i]=0;
                return;
            }
        }
    }
    if (raw_volume && (!(area > 0.0) || !isfinite(volume) || nc < 3)) {
        atomicAdd(&d_gfs_invalid_raw_volume_cells, 1ull);
        face_count_out[i] = 0;
        return;
    }
    /* Volume floor at 0.1% of nominal cell area (= cellsize²).  Sutherland-
     * Hodgman clipping at corners of rings IC (Nx≥256) sometimes produces
     * polygon area ~1e-30 → density = mass/V → 1.5e+25 garbage → step-1
     * K explosion via conduction.  See bug_gpu_sedov_zerovol.md. */
    double volume_floor = 1.0e-3 * cellsize * cellsize;
    if (!raw_volume && volume < volume_floor) {
        printf("[GPU_VOL_DIAG] i=%d nc=%d V_raw=%g V_floor=%g cellsize=%g\n",
               i, nc, volume, volume_floor, cellsize);
        volume = volume_floor;
    }
    double density = massi / volume;

    /* ---- w2ceil and avgNeighboringPressure ---- */
    double w2ceil_val = 1e20;
    double avgNP_sum = 0, avgNP_len = 0;
    for (int e = 0; e < nc; e++) {
        int en = (e + 1) % nc;
        int j_neigh = poly_neigh[e];
        if (j_neigh >= 0) {
            /* w2ceil = min over real faces of (dist2_j + w2_j) */
            double jdx = px[j_neigh] - xi;
            double jdy = py[j_neigh] - yi;
            double jdist2 = jdx * jdx + jdy * jdy;
            double val = jdist2 + (all_weights_zero ? 0.0 : (double)pw2[j_neigh]);
            if (val < w2ceil_val) w2ceil_val = val;

            /* avgNeighP: weighted by face edge length */
            double edx = poly_x[en] - poly_x[e];
            double edy = poly_y[en] - poly_y[e];
            double ds = sqrt(edx * edx + edy * edy);
            avgNP_sum += ds * (double)ppressure_in[j_neigh];
            avgNP_len += ds;
        }
    }
    double avgNP = (avgNP_len > 0) ? avgNP_sum / avgNP_len : (double)ppressure_in[i];

    /* ---- Gradients: Green-Gauss (gradient_method==0) or Pakmor 2016 LSF (==1) ---- */
    int use_lsf = (gradient_method == 1);
    /* Green-Gauss accumulators */
    double sum_gUxx = 0, sum_gUxy = 0, sum_gUyx = 0, sum_gUyy = 0;
    double sum_dPdx = 0, sum_dPdy = 0;
    double sum_dRhodx = 0, sum_dRhody = 0;
    /* LSF normal-equation matrix M = Σ w (Δr Δr^T) and per-Q RHS */
    double Mxx = 0, Mxy = 0, Myy = 0;
    double bxVx = 0, byVx = 0, bxVy = 0, byVy = 0;
    double bxP  = 0, byP  = 0, bxR  = 0, byR  = 0;
    /* Hyperviscosity Pass 1 accumulators (real neighbors only) */
    double sum_lap_vx = 0, sum_lap_vy = 0;

    for (int e = 0; e < nc; e++) {
        int en = (e + 1) % nc;
        int j_neigh = poly_neigh[e];
        if (j_neigh < 0) continue;  /* box boundary edge, skip */

        /* dS = (c2y - c1y, -(c2x - c1x)) */
        double dSx = poly_y[en] - poly_y[e];
        double dSy = -(poly_x[en] - poly_x[e]);

        double dx_ij = px[j_neigh] - xi;
        double dy_ij = py[j_neigh] - yi;
        double d2_ij = dx_ij * dx_ij + dy_ij * dy_ij;

        /* Hyperviscosity Pass 1: accumulate (v_j - v_i) weighted by A_face/d_ij;
         * exclude ghost mirrors so wall noise doesn't leak into biharmonic. */
        if (pindx[j_neigh] != MAX_INDEX_GPU && d2_ij > 0) {
            double A_face_lap = sqrt(dSx * dSx + dSy * dSy);
            double invd_lap = A_face_lap / sqrt(d2_ij);
            sum_lap_vx += invd_lap * (pvx[j_neigh] - pvx[i]);
            sum_lap_vy += invd_lap * (pvy[j_neigh] - pvy[i]);
        }

        if (use_lsf) {
            /* Pakmor 2016 LSF: weight = A_face / |Δr|^2 */
            double A_face = sqrt(dSx * dSx + dSy * dSy);
            double w = (d2_ij > 0) ? A_face / d2_ij : 0.0;
            double dvx = pvx[j_neigh] - pvx[i];
            double dvy = pvy[j_neigh] - pvy[i];
            double dP  = (double)ppressure_in[j_neigh] - (double)ppressure_in[i];
            double dR  = (double)pden_in[j_neigh] - (double)pden_in[i];
            Mxx += w * dx_ij * dx_ij;
            Mxy += w * dx_ij * dy_ij;
            Myy += w * dy_ij * dy_ij;
            bxVx += w * dx_ij * dvx;  byVx += w * dy_ij * dvx;
            bxVy += w * dx_ij * dvy;  byVy += w * dy_ij * dvy;
            bxP  += w * dx_ij * dP;   byP  += w * dy_ij * dP;
            bxR  += w * dx_ij * dR;   byR  += w * dy_ij * dR;
        } else {
            double wfrac_g = (all_weights_zero || !(d2_ij > 0)) ? 0.5
                : 0.5 + 0.5 * (w2i - (double)pw2[j_neigh]) / d2_ij;

            int is_ghost = (pindx[j_neigh] == MAX_INDEX_GPU);

            double vx_face, vy_face, P_face, Rho_face;

            if (is_ghost) {
                /* Wall ghost: simple weighted average */
                double vxi = pvx[i], vyi = pvy[i];
                double vxj = pvx[j_neigh], vyj = pvy[j_neigh];
                vx_face  = vxi + wfrac_g * (vxj - vxi);
                vy_face  = vyi + wfrac_g * (vyj - vyi);
                P_face   = (double)ppressure_in[i] + wfrac_g * ((double)ppressure_in[j_neigh] - (double)ppressure_in[i]);
                Rho_face = (double)pden_in[i] + wfrac_g * ((double)pden_in[j_neigh] - (double)pden_in[i]);
            } else {
                int kp_neigh = poly_neigh[en];
                int km_neigh = poly_neigh[(e - 1 + nc) % nc];

                if (OoA > 0 && kp_neigh >= 0 && km_neigh >= 0) {
                    double OoA3 = OoA / 3.0;
                    double vxi = pvx[i], vyi = pvy[i];
                    double vxj = pvx[j_neigh], vyj = pvy[j_neigh];
                    double pi_val = (double)ppressure_in[i], pj_val = (double)ppressure_in[j_neigh];
                    double rhoi = (double)pden_in[i], rhoj = (double)pden_in[j_neigh];

                    vx_face  = vxi + wfrac_g * (vxj - vxi)
                             + OoA3 * (pvx[kp_neigh] + pvx[km_neigh] - vxi - vxj);
                    vy_face  = vyi + wfrac_g * (vyj - vyi)
                             + OoA3 * (pvy[kp_neigh] + pvy[km_neigh] - vyi - vyj);
                    P_face   = pi_val + wfrac_g * (pj_val - pi_val)
                             + OoA3 * ((double)ppressure_in[kp_neigh] + (double)ppressure_in[km_neigh] - pi_val - pj_val);
                    Rho_face = rhoi + wfrac_g * (rhoj - rhoi)
                             + OoA3 * ((double)pden_in[kp_neigh] + (double)pden_in[km_neigh] - rhoi - rhoj);
                } else {
                    double vxi = pvx[i], vyi = pvy[i];
                    double vxj = pvx[j_neigh], vyj = pvy[j_neigh];
                    vx_face  = vxi + wfrac_g * (vxj - vxi);
                    vy_face  = vyi + wfrac_g * (vyj - vyi);
                    P_face   = (double)ppressure_in[i] + wfrac_g * ((double)ppressure_in[j_neigh] - (double)ppressure_in[i]);
                    Rho_face = (double)pden_in[i] + wfrac_g * ((double)pden_in[j_neigh] - (double)pden_in[i]);
                }
            }

            sum_gUxx += vx_face * dSx;
            sum_gUxy += vx_face * dSy;
            sum_gUyx += vy_face * dSx;
            sum_gUyy += vy_face * dSy;
            sum_dPdx += P_face * dSx;
            sum_dPdy += P_face * dSy;
            sum_dRhodx += Rho_face * dSx;
            sum_dRhody += Rho_face * dSy;
        }
    }

    double gUxx, gUxy, gUyx, gUyy;
    double dPdx_val, dPdy_val, dRhodx_val, dRhody_val;
    if (use_lsf) {
        double det = Mxx * Myy - Mxy * Mxy;
        if (fabs(det) > 1e-30) {
            double invDet = 1.0 / det;
            gUxx = ( Myy*bxVx - Mxy*byVx) * invDet;
            gUxy = (-Mxy*bxVx + Mxx*byVx) * invDet;
            gUyx = ( Myy*bxVy - Mxy*byVy) * invDet;
            gUyy = (-Mxy*bxVy + Mxx*byVy) * invDet;
            dPdx_val   = ( Myy*bxP - Mxy*byP) * invDet;
            dPdy_val   = (-Mxy*bxP + Mxx*byP) * invDet;
            dRhodx_val = ( Myy*bxR - Mxy*byR) * invDet;
            dRhody_val = (-Mxy*bxR + Mxx*byR) * invDet;
        } else {
            gUxx = gUxy = gUyx = gUyy = 0;
            dPdx_val = dPdy_val = dRhodx_val = dRhody_val = 0;
        }
    } else {
        double invVol = 1.0 / volume;
        gUxx = sum_gUxx * invVol; gUxy = sum_gUxy * invVol;
        gUyx = sum_gUyx * invVol; gUyy = sum_gUyy * invVol;
        dPdx_val = sum_dPdx * invVol; dPdy_val = sum_dPdy * invVol;
        dRhodx_val = sum_dRhodx * invVol; dRhody_val = sum_dRhody * invVol;
    }
    double divv_val = gUxx + gUyy;

    /* ---- TVD slope limiter (av_mode==5) ---- */
    if (av_mode == 5) {
        double rho_i = density, vx_i = pvx[i], vy_i = pvy[i];
        double P_i = (double)ppressure_in[i];  /* use input pressure for limiter */
        double rho_max = rho_i, rho_min = rho_i;
        double vx_max = vx_i, vx_min = vx_i;
        double vy_max = vy_i, vy_min = vy_i;
        double P_max = P_i, P_min = P_i;

        /* Pass A: min/max over real neighbors */
        for (int e = 0; e < nc; e++) {
            int jn = poly_neigh[e];
            if (jn < 0) continue;
            if (pindx[jn] == MAX_INDEX_GPU) continue;
            double rj = (double)pden_in[jn];
            double vxj = pvx[jn], vyj = pvy[jn];
            double pj = (double)ppressure_in[jn];
            if (rj > rho_max) rho_max = rj;  if (rj < rho_min) rho_min = rj;
            if (vxj > vx_max) vx_max = vxj;  if (vxj < vx_min) vx_min = vxj;
            if (vyj > vy_max) vy_max = vyj;  if (vyj < vy_min) vy_min = vyj;
            if (pj > P_max) P_max = pj;      if (pj < P_min) P_min = pj;
        }

        /* Pass B: tighten alpha */
        double alpha_rho = 1.0, alpha_vx = 1.0, alpha_vy = 1.0, alpha_P = 1.0;
        for (int e = 0; e < nc; e++) {
            if (poly_neigh[e] < 0) continue;
            int en = (e + 1) % nc;
            double dxf = 0.5 * (poly_x[e] + poly_x[en]);
            double dyf = 0.5 * (poly_y[e] + poly_y[en]);

            double drhoF = dRhodx_val * dxf + dRhody_val * dyf;
            double dvxF  = gUxx * dxf + gUxy * dyf;
            double dvyF  = gUyx * dxf + gUyy * dyf;
            double dPF   = dPdx_val * dxf + dPdy_val * dyf;

            #define TVD_LIMIT_GPU(alpha, delta, Vi, Vmax, Vmin) do { \
                if (delta > 1e-30) { \
                    double r = (Vmax - Vi) / delta; \
                    if (r < 0) r = 0; \
                    if (r < alpha) alpha = r; \
                } else if (delta < -1e-30) { \
                    double r = (Vmin - Vi) / delta; \
                    if (r < 0) r = 0; \
                    if (r < alpha) alpha = r; \
                } \
            } while(0)

            TVD_LIMIT_GPU(alpha_rho, drhoF, rho_i, rho_max, rho_min);
            TVD_LIMIT_GPU(alpha_vx,  dvxF,  vx_i,  vx_max,  vx_min);
            TVD_LIMIT_GPU(alpha_vy,  dvyF,  vy_i,  vy_max,  vy_min);
            TVD_LIMIT_GPU(alpha_P,   dPF,   P_i,   P_max,   P_min);
            #undef TVD_LIMIT_GPU
        }

        /* Apply limiter */
        dRhodx_val *= alpha_rho;  dRhody_val *= alpha_rho;
        if (limit_velocity == 2) {
            /* Use every real face neighbor. tmp_x/y are no longer needed
             * for clipping; keep geometry intact and reuse this scratch. */
            double sample_x[GPU_MAX_CORNERS], sample_y[GPU_MAX_CORNERS];
            int nsample = 0;
            double input_scale = hypot(pvx[i],pvy[i]);
            for (int e = 0; e < nc; e++) {
                int jn = poly_neigh[e];
                if (jn < 0 || pindx[jn] == MAX_INDEX_GPU) continue;
                sample_x[nsample] = px[jn] - xi;
                sample_y[nsample] = py[jn] - yi;
                tmp_x[nsample] = pvx[jn] - pvx[i];
                tmp_y[nsample] = pvy[jn] - pvy[i];
                input_scale = fmax(input_scale,hypot(pvx[jn],pvy[jn]));
                nsample++;
            }
            double alpha_v = gfs_velocity_midpoint_alpha(nsample,
                sample_x, sample_y, tmp_x, tmp_y, gUxx, gUxy, gUyx, gUyy,
                gfs_velocity_matrix_condition(Mxx,Mxy,Myy),input_scale);
            if (!(alpha_v >= 0.0)) {
                atomicAdd(&d_gfs_invalid_velocity_fit_cells,1ull);
                face_count_out[i] = 0;
                return;
            }
            gUxx *= alpha_v; gUxy *= alpha_v;
            gUyx *= alpha_v; gUyy *= alpha_v;
        } else if (limit_velocity) {
            gUxx *= alpha_vx;     gUxy *= alpha_vx;
            gUyx *= alpha_vy;     gUyy *= alpha_vy;
        }
        dPdx_val *= alpha_P;      dPdy_val *= alpha_P;
        divv_val = gUxx + gUyy;  /* recompute after limiting */
    }

    /* ---- Write outputs ---- */
    out_volume[i]   = (float)volume;
    if (thermal64) out_V64[i] = volume;
    out_den[i]      = (float)density;
    out_w2ceil[i]   = (float)w2ceil_val;
    out_avgNeighP[i]= (float)avgNP;

    out_gUxx[i] = gUxx;  out_gUxy[i] = gUxy;
    out_gUyx[i] = gUyx;  out_gUyy[i] = gUyy;
    out_dPdx[i] = dPdx_val;  out_dPdy[i] = dPdy_val;
    out_dRhodx[i] = dRhodx_val;  out_dRhody[i] = dRhody_val;

    /* Hyperviscosity Pass 1 finalize: lap_v = (1/V) Σ A_face/d_ij × Δv */
    {
        double invVol_lap = (volume > 0) ? 1.0 / volume : 0.0;
        out_lap_vx[i] = sum_lap_vx * invVol_lap;
        out_lap_vy[i] = sum_lap_vy * invVol_lap;
    }

    /* Write PolygonOut for face scatter */
    PolygonOut *po = &polygons[i];
    po->n_corners = nc;
    po->stencil_pad = stencil_pad;
    int n_real_faces = 0;
    for (int e = 0; e < nc; e++) {
        po->cx[e] = poly_x[e];
        po->cy[e] = poly_y[e];
        po->edge_neigh[e] = poly_neigh[e];
        int en = (e+1)%nc;
        double ex = poly_x[en]-poly_x[e], ey = poly_y[en]-poly_y[e];
        if (poly_neigh[e] >= 0 && (raw_volume < 2 || ex*ex+ey*ey > 0.0)) n_real_faces++;
    }
    po->n_real_faces = n_real_faces;
    if (nc >= 3) {
        double final_maxdist2 = 0.0;
        for (int e = 0; e < nc; e++) {
            const double cd = poly_x[e] * poly_x[e] + poly_y[e] * poly_y[e];
            if (cd > final_maxdist2) final_maxdist2 = cd;
        }
        const double final_radius = sqrt(final_maxdist2);
        const double deltaW = fmax(0.0, global_w2_max - w2i);
        const double d_clamped = (deltaW > 0.0)
            ? fmin(sqrt(deltaW / 0.9), 20.0 * final_radius) : 0.0;
        const double d_unclamped = final_radius + sqrt(final_maxdist2 + deltaW);
        const double final_search_radius = fmax(2.0 * final_radius,
                                                 fmax(d_clamped, d_unclamped));
        const double grid_xmax = xmin + (double)mx * cellsize;
        const double grid_ymax = ymin + (double)my * cellsize;
        po->halo_clipped =
            (xi - final_search_radius < xmin || xi + final_search_radius > grid_xmax ||
             yi - final_search_radius < ymin || yi + final_search_radius > grid_ymax);
        if (po->halo_clipped)
            atomicAdd(&d_gfs_stencil_halo_clipped_cells, 1ull);
    } else {
        po->halo_clipped = 0;
    }
    face_count_out[i] = n_real_faces;
}

/* ================================================================
 *  Kernel B: pressure_stress_kernel
 *
 *  Computes P, cs from EOS, NS stress from gradients.
 *  1 thread per real particle, trivial arithmetic.
 * ================================================================ */
__global__
void pressure_stress_kernel(
    const float *__restrict__ pie, const float *__restrict__ pvolume,
    const float *__restrict__ pden,
    float *__restrict__ pK,         /* in/out: K is clamped to K_floor when entropy_mode==1 */
    const double *__restrict__ pgUxx, const double *__restrict__ pgUxy,
    const double *__restrict__ pgUyx, const double *__restrict__ pgUyy,
    const double *__restrict__ palpha_cd,
    float *__restrict__ out_pressure, float *__restrict__ out_csound,
    double *__restrict__ out_tauxx, double *__restrict__ out_tauxy,
    double *__restrict__ out_tauyy,
    double *__restrict__ out_divv,
    int nbp, int av_mode, double Gamma, double nu_phys, double cd_amax,
    int entropy_mode, double K_floor,
    int thermal64, const double *__restrict__ pU64,
    const double *__restrict__ pV64, const double *__restrict__ pmass)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nbp) return;

    double ie_val = (double)pie[i];
    double vol    = (double)pvolume[i];
    double den    = (double)pden[i];
    if (thermal64) {
        vol = pV64[i];
        ie_val = pU64[i];
        den = pmass[i] / vol;
    }

    double P;
    if (entropy_mode == 1) {
        /* K-EOS: P = K · ρ^γ. Mirrors CPU updateDenW2Pressure2DBlend
         * (exam.c:2686-2692): clamp K to floor, then derive P. */
        double K_val = (double)pK[i];
        if (K_val < K_floor) { K_val = K_floor; pK[i] = (float)K_val; }
        P = K_val * pow(den, Gamma);
    } else {
        P = ie_val / vol * (Gamma - 1.0);
    }
    if (P <= 0) P = 1e-6;
    double cs = sqrt(Gamma * P / den);

    out_pressure[i] = (float)P;
    out_csound[i]   = (float)cs;

    double gxx = pgUxx[i], gxy = pgUxy[i];
    double gyx = pgUyx[i], gyy = pgUyy[i];
    double divv_val = gxx + gyy;
    out_divv[i] = divv_val;

    if (av_mode == 0 && nu_phys > 0) {
        out_tauxx[i] = -nu_phys * den * (2.0 * gxx - (2.0/3.0) * divv_val);
        out_tauxy[i] = -nu_phys * den * (gxy + gyx);
        out_tauyy[i] = -nu_phys * den * (2.0 * gyy - (2.0/3.0) * divv_val);
    } else if (av_mode == 1) {
        /* NS-stress uses ONLY physical viscosity. CD10 enters as viscous *pressure*
         * Π_CD10 = ½ α vsig ρ̄ (-Δv) in the face loop (line ~612), not as NS-stress
         * traction. The old nu_cd = α_cd·sqrt(V)·c_s injected an explicit 1/Nx factor
         * that, combined with the 1/m ~ Nx² in dK conversion, drove K-runaway at
         * high Nx. This matches av_mode>=2 design intent (CD10 as viscous pressure). */
        double nu = (nu_phys > 0) ? nu_phys : 0.0;
        out_tauxx[i] = -nu * den * (2.0 * gxx - (2.0/3.0) * divv_val);
        out_tauxy[i] = -nu * den * (gxy + gyx);
        out_tauyy[i] = -nu * den * (2.0 * gyy - (2.0/3.0) * divv_val);
    } else {
        out_tauxx[i] = 0;
        out_tauxy[i] = 0;
        out_tauyy[i] = 0;
    }
}

/* ================================================================
 *  Kernel C: scatter_faces_kernel
 *
 *  Reads PolygonOut per particle, writes face data to FaceCSR
 *  using prefix-summed offsets.
 * ================================================================ */
__global__
void scatter_faces_kernel(
    const PolygonOut *__restrict__ polygons,
    const int *__restrict__ face_offset,
    const long long *__restrict__ pindx,
    /* FaceCSR output arrays */
    double *__restrict__ f_c1x, double *__restrict__ f_c1y,
    double *__restrict__ f_c2x, double *__restrict__ f_c2y,
    int *__restrict__ f_neighbor_idx,
    int *__restrict__ f_kp_idx, int *__restrict__ f_km_idx,
    int *__restrict__ f_is_ghost,
    int nbp, int omit_zero_edges)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nbp) return;

    const PolygonOut *po = &polygons[i];
    int nc = po->n_corners;
    int base = face_offset[i];
    int local_f = 0;

    for (int e = 0; e < nc; e++) {
        int j_neigh = po->edge_neigh[e];
        if (j_neigh < 0) continue;  /* skip box boundary edges */

        int en = (e + 1) % nc;
        double ex=po->cx[en]-po->cx[e], ey=po->cy[en]-po->cy[e];
        if (omit_zero_edges && !(ex*ex+ey*ey > 0.0)) continue;
        int dst = base + local_f;

        f_c1x[dst] = po->cx[e];
        f_c1y[dst] = po->cy[e];
        f_c2x[dst] = po->cx[en];
        f_c2y[dst] = po->cy[en];
        f_neighbor_idx[dst] = j_neigh;

        /* kp = edge_neigh of next edge, km = edge_neigh of prev edge */
        int kp = po->edge_neigh[en];
        int km = po->edge_neigh[(e - 1 + nc) % nc];
        f_kp_idx[dst] = (kp >= 0) ? kp : -1;
        f_km_idx[dst] = (km >= 0) ? km : -1;

        f_is_ghost[dst] = (pindx[j_neigh] == MAX_INDEX_GPU) ? 1 : 0;
        local_f++;
    }
}

/* ================================================================
 *  Orchestrator: gpu_tessellate_and_build_faces
 *
 *  Runs kernels A → B → prefix sum → C, producing FaceCSR on device.
 * ================================================================ */
extern "C"
void gpu_tessellate_and_build_faces(GPUContext *ctx, int nbp,
    int mx, int my, double cellsize, double xmin, double ymin,
    double boxsize, double global_w2_max, int all_weights_zero,
    double OoA, int av_mode, double Gamma,
    double nu_phys, double prandtl, double cd_amax,
    int gradient_method,
    int entropy_mode, double K_floor)
{
    h_diag_tess_mx = mx; h_diag_tess_my = my;
    h_diag_tess_cellsize = cellsize;
    h_diag_tess_xmin = xmin; h_diag_tess_ymin = ymin;
    h_diag_global_w2_max = global_w2_max;
    cudaStream_t s = (cudaStream_t)ctx->stream;
    int blockSize = 256;
    int gridSize = (nbp + blockSize - 1) / blockSize;
    int n_total = ctx->d_parts.n_total;
    int thermal64 = gfs_thermal64_requested();
    static unsigned long long tess_profile_call = 0;
    tess_profile_call++;
    const char *tess_profile_env = getenv("GFS_GPU_TESS_PROFILE");
    int tess_profile = tess_profile_env && atoi(tess_profile_env) > 0
        && (tess_profile_call <= 8 || tess_profile_call % 1024 == 0);
    cudaEvent_t tess_events[6];
    if (tess_profile) {
        for (int k = 0; k < 6; k++) CUDA_CHECK(cudaEventCreate(&tess_events[k]));
    }

    /* Causal controls for the explicit affine patch test only. Production
     * retains both limiters. Rho/P positivity limiting is never bypassed. */
    const char *probe_v = getenv("GFS_AFFINE_UNLIMITED_VELOCITY");
    const char *probe_m = getenv("GFS_AFFINE_UNMERGED_FACES");
    const char *velocity_env = getenv("GFS_VELOCITY_LIMITER");
    int limit_velocity = (velocity_env && atoi(velocity_env) > 0) ? 2 : 1;
    if (probe_v && atoi(probe_v) > 0) limit_velocity = 0;
    if (limit_velocity == 2 && gradient_method != 1) {
        fprintf(stderr,"[VELOCITY_LIMITER] requires LSF gradient_method=1\n"); exit(1);
    }
    const char *raw_env = getenv("GFS_GPU_RAW_FACES");
    int raw_geometry = (raw_env && atoi(raw_env) > 0) ? 2 : 0;
    int merge_small_faces = !(probe_m && atoi(probe_m) > 0);
    if (raw_geometry) merge_small_faces = 0;
    if ((probe_v && atoi(probe_v)>0) || (probe_m && atoi(probe_m)>0)) {
        const char *ic = getenv("EUNHA_IC");
        if (!ic || strcmp(ic, "affine_probe") != 0) {
            fprintf(stderr, "[AFFINE_PROBE] diagnostic controls require EUNHA_IC=affine_probe\n");
            exit(1);
        }
        static int banner = 0;
        if (!banner) {
            fprintf(stderr, "[AFFINE_PROBE] velocity_limiter=%d merge_small_faces=%d; rho/P limiting retained\n",
                    limit_velocity, merge_small_faces);
            banner = 1;
        }
    }

    unsigned long long zero_fallbacks = 0;
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_invalid_raw_volume_cells,&zero_fallbacks,sizeof(zero_fallbacks),0,cudaMemcpyHostToDevice,s));
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_invalid_velocity_fit_cells,&zero_fallbacks,sizeof(zero_fallbacks),0,cudaMemcpyHostToDevice,s));
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_tess_geometry_errors,&zero_fallbacks,sizeof(zero_fallbacks),0,cudaMemcpyHostToDevice,s));
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_candidate_fallback_cells, &zero_fallbacks, sizeof(zero_fallbacks), 0,
                                       cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_stencil_expanded_cells, &zero_fallbacks, sizeof(zero_fallbacks), 0,
                                       cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_stencil_halo_clipped_cells, &zero_fallbacks, sizeof(zero_fallbacks), 0,
                                       cudaMemcpyHostToDevice, s));

    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[0], s));
    /* Kernel A reads neighbor densities while publishing each cell's new
     * density. Passing d_parts.den as both input and output made those reads
     * depend on block execution order, and violated the restrict contract.
     * Freeze exactly the uploaded stage input (including MPI padding); the
     * copy and kernel share a stream, so every read sees one time level. */
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_tess_den_input, ctx->d_parts.den,
                               n_total * sizeof(float),
                               cudaMemcpyDeviceToDevice, s));

    /* Kernel A: Voronoi tessellation + gradients */
    voronoi_tessellate_kernel<<<gridSize, blockSize, 0, s>>>(
        ctx->d_parts.x, ctx->d_parts.y,
        ctx->d_parts.vx, ctx->d_parts.vy,
        ctx->d_parts.mass, ctx->d_parts.indx,
        ctx->d_parts.ie, ctx->d_parts.w2,
        ctx->d_parts.pressure, ctx->d_tess_den_input,
        ctx->d_cells.cell_offset, ctx->d_cells.particle_index,
        /* outputs */
        ctx->d_parts.den, ctx->d_parts.pressure,
        ctx->d_parts.volume, ctx->d_parts.w2ceil, ctx->d_parts.avgNeighP,
        ctx->d_parts.gUxx, ctx->d_parts.gUxy,
        ctx->d_parts.gUyx, ctx->d_parts.gUyy,
        ctx->d_parts.dPdx, ctx->d_parts.dPdy,
        ctx->d_parts.dRhodx, ctx->d_parts.dRhody,
        ctx->d_parts.lap_vx, ctx->d_parts.lap_vy,
        ctx->d_polygons, ctx->d_face_count,
        nbp, n_total, mx, my, cellsize, xmin, ymin, boxsize, global_w2_max,
        all_weights_zero, OoA, av_mode,
        gradient_method, raw_geometry, limit_velocity, merge_small_faces,
        thermal64, ctx->d_parts.E_inv_xy);
    CUDA_CHECK(cudaGetLastError());
    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[1], s));

    /* Kernel B: Pressure, sound speed, NS stress */
    pressure_stress_kernel<<<gridSize, blockSize, 0, s>>>(
        ctx->d_parts.ie, ctx->d_parts.volume, ctx->d_parts.den,
        ctx->d_parts.K,
        ctx->d_parts.gUxx, ctx->d_parts.gUxy,
        ctx->d_parts.gUyx, ctx->d_parts.gUyy,
        ctx->d_parts.alpha_cd,
        ctx->d_parts.pressure, ctx->d_parts.csound,
        ctx->d_parts.tauxx, ctx->d_parts.tauxy, ctx->d_parts.tauyy,
        ctx->d_parts.divv,
        nbp, av_mode, Gamma, nu_phys, cd_amax,
        entropy_mode, K_floor, thermal64,
        ctx->d_parts.E_inv_xx, ctx->d_parts.E_inv_xy, ctx->d_parts.mass);
    CUDA_CHECK(cudaGetLastError());
    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[2], s));

    /* CUB prefix sum: face_count → face_offset */
    cub::DeviceScan::ExclusiveSum(ctx->d_cub_scan_temp, ctx->cub_scan_temp_bytes,
                                   ctx->d_face_count, ctx->d_face_offset,
                                   nbp, s);
    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[3], s));

    /* Copy total face count (last element) to set face_offset[nbp] */
    /* We need face_offset[nbp] = face_offset[nbp-1] + face_count[nbp-1].
       CUB ExclusiveSum of N elements produces N outputs. We need N+1.
       Compute total via: copy face_offset[nbp-1] + face_count[nbp-1] → face_offset[nbp].
       Simpler: use InclusiveSum of face_count → gives cumulative sums,
       then shift. But ExclusiveSum is more standard. We'll just compute the
       total on host after a small copy. */
    {
        /* Alternative: run a tiny kernel to set face_offset[nbp] */
        /* For simplicity, sync and read last values */
        int h_last_offset, h_last_count;
        CUDA_CHECK(cudaMemcpyAsync(&h_last_offset, ctx->d_face_offset + nbp - 1,
                   sizeof(int), cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(&h_last_count, ctx->d_face_count + nbp - 1,
                   sizeof(int), cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaStreamSynchronize(s));
        unsigned long long h_candidate_fallbacks = 0;
        CUDA_CHECK(cudaMemcpyFromSymbol(&h_candidate_fallbacks, d_gfs_candidate_fallback_cells,
                                       sizeof(h_candidate_fallbacks)));
        if (h_candidate_fallbacks > 0)
            fprintf(stderr, "[GPU_CANDIDATE_STREAM] cells_over_legacy_cap=%llu legacy_cap=%d all_candidates_accounted_for=1\n",
                    h_candidate_fallbacks, GPU_MAX_NEIGH);
        unsigned long long h_stencil_expanded = 0, h_halo_clipped = 0;
        CUDA_CHECK(cudaMemcpyFromSymbol(&h_stencil_expanded, d_gfs_stencil_expanded_cells,
                                       sizeof(h_stencil_expanded)));
        CUDA_CHECK(cudaMemcpyFromSymbol(&h_halo_clipped, d_gfs_stencil_halo_clipped_cells,
                                       sizeof(h_halo_clipped)));
        if (h_stencil_expanded > 0 || h_halo_clipped > 0)
            fprintf(stderr, "[GPU_STENCIL_EXPANSION] expanded_cells=%llu halo_clipped_cells=%llu\n",
                    h_stencil_expanded, h_halo_clipped);
        unsigned long long invalid_volume=0,invalid_fit=0,geometry_errors=0;
        CUDA_CHECK(cudaMemcpyFromSymbol(&invalid_volume,d_gfs_invalid_raw_volume_cells,sizeof(invalid_volume)));
        CUDA_CHECK(cudaMemcpyFromSymbol(&invalid_fit,d_gfs_invalid_velocity_fit_cells,sizeof(invalid_fit)));
        CUDA_CHECK(cudaMemcpyFromSymbol(&geometry_errors,d_gfs_tess_geometry_errors,sizeof(geometry_errors)));
        if (invalid_volume || invalid_fit || geometry_errors || (raw_geometry && h_halo_clipped)) {
            fprintf(stderr,"[GPU_TESS_INVALID] volume=%llu velocity_fit=%llu geometry=%llu halo=%llu\n",
                    invalid_volume,invalid_fit,geometry_errors,h_halo_clipped);
            if (h_geometry_trial) {
                h_geometry_failed=1; ctx->tess_faces_on_device=0; return;
            }
#ifdef GFS_VARPRESS_W0
            gfs_rk4_geometry_abort_snapshot();
#endif
            exit(1);
        }
        int total_faces = h_last_offset + h_last_count;

        /* Set face_offset[nbp] on device */
        CUDA_CHECK(cudaMemcpyAsync(ctx->d_face_offset + nbp, &total_faces,
                   sizeof(int), cudaMemcpyHostToDevice, s));

        /* Update FaceCSR metadata */
        ctx->d_faces.n_particles = nbp;
        ctx->d_faces.n_faces_total = total_faces;

        /* Check capacity (should not happen with proper allocation) */
        if (total_faces > ctx->max_faces) {
            fprintf(stderr, "[GPU TESS] total_faces=%d > max_faces=%d; refusing truncated faces\n",total_faces,ctx->max_faces);
            if (h_geometry_trial) {
                h_geometry_failed=1; ctx->tess_faces_on_device=0; return;
            }
#ifdef GFS_VARPRESS_W0
            gfs_rk4_geometry_abort_snapshot();
#endif
            exit(1);
        }
    }

    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[4], s));
    /* Kernel C: Scatter faces → FaceCSR */
    scatter_faces_kernel<<<gridSize, blockSize, 0, s>>>(
        ctx->d_polygons, ctx->d_face_offset, ctx->d_parts.indx,
        ctx->d_faces.c1x, ctx->d_faces.c1y,
        ctx->d_faces.c2x, ctx->d_faces.c2y,
        ctx->d_faces.neighbor_idx,
        ctx->d_faces.kp_idx, ctx->d_faces.km_idx,
        ctx->d_faces.is_ghost,
        nbp, raw_geometry >= 2);
    CUDA_CHECK(cudaGetLastError());

    /* Copy face_offset to FaceCSR (device already has it, just set pointer) */
    /* The force kernel uses d_faces.face_offset. Copy d_face_offset → d_faces.face_offset */
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_faces.face_offset, ctx->d_face_offset,
               (nbp + 1) * sizeof(int), cudaMemcpyDeviceToDevice, s));

    if (tess_profile) CUDA_CHECK(cudaEventRecord(tess_events[5], s));
    CUDA_CHECK(cudaStreamSynchronize(s));
    if (tess_profile) {
        float ms[5];
        for (int k = 0; k < 5; k++)
            CUDA_CHECK(cudaEventElapsedTime(&ms[k], tess_events[k], tess_events[k+1]));
        fprintf(stderr, "[GPU tess device profile] call=%llu density_copy_A=%.3f B=%.3f scan=%.3f host_roundtrip=%.3f scatter=%.3f ms\n",
            tess_profile_call, ms[0], ms[1], ms[2], ms[3], ms[4]);
        for (int k = 0; k < 6; k++) CUDA_CHECK(cudaEventDestroy(tess_events[k]));
    }
    ctx->tess_faces_on_device = 1;
}

/* ================================================================
 *  gpu_download_tess_results: device SoA → pinned host SoA
 *
 *  Downloads the tessellation outputs (den, P, cs, volume, w2ceil,
 *  avgNeighP, gradients, stress) so CPU can write them back to AoS.
 * ================================================================ */
/* End-of-step geometry uses the same complete candidate search, but keeps
 * raw positive area: CPU exam2dUpdateVol neither merges small faces nor applies
 * the force-stage volume floor. No pressure/EOS or face-CSR update is performed.
 * Return failure to the MPI wrapper before any host particle is changed. */
extern "C"
int gpu_tessellate_raw_volume(GPUContext *ctx, int nbp,
    int mx, int my, double cellsize, double xmin, double ymin,
    double boxsize, double global_w2_max, double OoA, int av_mode,
    int gradient_method)
{
    ctx->tess_faces_on_device = 0;
    if (nbp <= 0) return 0;
    cudaStream_t s = (cudaStream_t)ctx->stream;
    int blockSize = 256, gridSize = (nbp + blockSize - 1) / blockSize;
    int n_total = ctx->d_parts.n_total;
    unsigned long long zero = 0;
    const char *raw_env=getenv("GFS_GPU_RAW_FACES");
    int raw_geometry=(raw_env && atoi(raw_env)>0) ? 2 : 1;
    int thermal64=gfs_thermal64_requested();
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_tess_geometry_errors,&zero,sizeof(zero),0,cudaMemcpyHostToDevice,s));
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_invalid_raw_volume_cells, &zero,
        sizeof(zero), 0, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_stencil_halo_clipped_cells, &zero,
        sizeof(zero), 0, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_stencil_expanded_cells, &zero,
        sizeof(zero), 0, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_candidate_fallback_cells, &zero,
        sizeof(zero), 0, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_tess_den_input, ctx->d_parts.den,
        n_total * sizeof(float), cudaMemcpyDeviceToDevice, s));
    voronoi_tessellate_kernel<<<gridSize, blockSize, 0, s>>>(
        ctx->d_parts.x, ctx->d_parts.y,
        ctx->d_parts.vx, ctx->d_parts.vy,
        ctx->d_parts.mass, ctx->d_parts.indx,
        ctx->d_parts.ie, ctx->d_parts.w2,
        ctx->d_parts.pressure, ctx->d_tess_den_input,
        ctx->d_cells.cell_offset, ctx->d_cells.particle_index,
        /* outputs */
        ctx->d_parts.den, ctx->d_parts.pressure,
        ctx->d_parts.volume, ctx->d_parts.w2ceil, ctx->d_parts.avgNeighP,
        ctx->d_parts.gUxx, ctx->d_parts.gUxy,
        ctx->d_parts.gUyx, ctx->d_parts.gUyy,
        ctx->d_parts.dPdx, ctx->d_parts.dPdy,
        ctx->d_parts.dRhodx, ctx->d_parts.dRhody,
        ctx->d_parts.lap_vx, ctx->d_parts.lap_vy,
        ctx->d_polygons, ctx->d_face_count,
        nbp, n_total, mx, my, cellsize, xmin, ymin, boxsize, global_w2_max,
        0, OoA, av_mode,
        gradient_method, raw_geometry, 1, 0,
        thermal64, ctx->d_parts.E_inv_xy);
    CUDA_CHECK(cudaGetLastError());
    unsigned long long invalid = 0, clipped = 0, geometry_errors=0;
    CUDA_CHECK(cudaMemcpyFromSymbolAsync(&geometry_errors,d_gfs_tess_geometry_errors,sizeof(geometry_errors),0,cudaMemcpyDeviceToHost,s));
    CUDA_CHECK(cudaMemcpyFromSymbolAsync(&invalid, d_gfs_invalid_raw_volume_cells,
        sizeof(invalid), 0, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyFromSymbolAsync(&clipped, d_gfs_stencil_halo_clipped_cells,
        sizeof(clipped), 0, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.volume, ctx->d_parts.volume,
        nbp * sizeof(float), cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.avgNeighP, ctx->d_parts.avgNeighP,
        nbp * sizeof(float), cudaMemcpyDeviceToHost, s));
    if (thermal64)
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.E_inv_xy, ctx->d_parts.E_inv_xy,
            nbp * sizeof(double), cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaStreamSynchronize(s));
    if (invalid || clipped || geometry_errors)
        fprintf(stderr, "[GPU_FINAL_VOLUME] invalid=%llu halo_clipped=%llu geometry=%llu\n",
                invalid, clipped,geometry_errors);
    return (invalid || clipped || geometry_errors) ? 1 : 0;
}

extern "C"
void gpu_download_tess_results(GPUContext *ctx, int nbp, int has_stress)
{
    cudaStream_t s = (cudaStream_t)ctx->stream;
    size_t sf = nbp * sizeof(float);
    size_t sd = nbp * sizeof(double);
    int thermal64 = gfs_thermal64_requested();

    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.den,       ctx->d_parts.den,       sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.pressure,  ctx->d_parts.pressure,  sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.csound,    ctx->d_parts.csound,    sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.volume,    ctx->d_parts.volume,    sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.w2ceil,    ctx->d_parts.w2ceil,    sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.avgNeighP, ctx->d_parts.avgNeighP, sf, cudaMemcpyDeviceToHost, s));

    if (has_stress) {
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.gUxx,    ctx->d_parts.gUxx,    sd, cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.gUxy,    ctx->d_parts.gUxy,    sd, cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.gUyx,    ctx->d_parts.gUyx,    sd, cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.gUyy,    ctx->d_parts.gUyy,    sd, cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dPdx,    ctx->d_parts.dPdx,    sd, cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dPdy,    ctx->d_parts.dPdy,    sd, cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dRhodx,  ctx->d_parts.dRhodx,  sd, cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dRhody,  ctx->d_parts.dRhody,  sd, cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.tauxx,   ctx->d_parts.tauxx,   sd, cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.tauxy,   ctx->d_parts.tauxy,   sd, cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.tauyy,   ctx->d_parts.tauyy,   sd, cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.divv,    ctx->d_parts.divv,    sd, cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.lap_vx,  ctx->d_parts.lap_vx,  sd, cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.lap_vy,  ctx->d_parts.lap_vy,  sd, cudaMemcpyDeviceToHost, s));
        /* K may have been clamped to K_floor by pressure_stress_kernel */
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.K,       ctx->d_parts.K,       sf, cudaMemcpyDeviceToHost, s));
        if (gfs_thermal64_requested())
            CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.E_inv_xy, ctx->d_parts.E_inv_xy, sd, cudaMemcpyDeviceToHost, s));
    }

    CUDA_CHECK(cudaStreamSynchronize(s));
}

/* ================================================================
 *  LagMFM (av_mode=4) upload/download
 * ================================================================ */

extern "C"
void gpu_upload_lagmfm_fields(GPUContext *ctx, int n)
{
    cudaStream_t s = (cudaStream_t)ctx->stream;
    size_t sd = n * sizeof(double);

    CUDA_CHECK(cudaMemcpyAsync(ctx->d_parts.E_inv_xx, ctx->h_parts.E_inv_xx, sd, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_parts.E_inv_xy, ctx->h_parts.E_inv_xy, sd, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_parts.E_inv_yx, ctx->h_parts.E_inv_yx, sd, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_parts.E_inv_yy, ctx->h_parts.E_inv_yy, sd, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_parts.h_mfm,    ctx->h_parts.h_mfm,    sd, cudaMemcpyHostToDevice, s));
}

extern "C"
void gpu_download_lagmfm_density(GPUContext *ctx, int n)
{
    cudaStream_t s = (cudaStream_t)ctx->stream;
    size_t sd = n * sizeof(double);
    size_t sf = n * sizeof(float);

    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.den,     ctx->d_parts.den,     sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.volume,  ctx->d_parts.volume,  sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.pressure,ctx->d_parts.pressure,sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.csound,  ctx->d_parts.csound,  sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.E_inv_xx,ctx->d_parts.E_inv_xx,sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.E_inv_xy,ctx->d_parts.E_inv_xy,sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.E_inv_yx,ctx->d_parts.E_inv_yx,sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.E_inv_yy,ctx->d_parts.E_inv_yy,sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.h_mfm,   ctx->d_parts.h_mfm,   sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.gUxx,    ctx->d_parts.gUxx,    sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.gUxy,    ctx->d_parts.gUxy,    sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.gUyx,    ctx->d_parts.gUyx,    sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.gUyy,    ctx->d_parts.gUyy,    sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dPdx,    ctx->d_parts.dPdx,    sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dPdy,    ctx->d_parts.dPdy,    sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dRhodx,  ctx->d_parts.dRhodx,  sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dRhody,  ctx->d_parts.dRhody,  sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.divv,    ctx->d_parts.divv,    sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.tauxx,   ctx->d_parts.tauxx,   sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.tauxy,   ctx->d_parts.tauxy,   sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.tauyy,   ctx->d_parts.tauyy,   sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.w2ceil,  ctx->d_parts.w2ceil,  sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaStreamSynchronize(s));
}

/* ================================================================
 *  LagMFM (av_mode=4) Wendland C2 kernel (2D), device inline
 *
 *  W(r,h) = (7/(4πh²)) · (1-q/2)⁴ · (1+2q),  q=r/h ∈ [0,2]
 * ================================================================ */
__device__ __forceinline__
double dev_wendland_c2_2d(double r, double h)
{
    if (h <= 0) return 0;
    double q = r / h;
    if (q >= 2.0) return 0;
    double u = 1.0 - 0.5 * q;
    double u2 = u * u;
    double u4 = u2 * u2;
    double C = 7.0 / (4.0 * M_PI * h * h);
    return C * u4 * (1.0 + 2.0 * q);
}

/* 2D cubic spline (M4), Hopkins (2015) fiducial; σ = 10/(7π h²). */
__device__ __forceinline__
double dev_cspline_2d(double r, double h)
{
    if (h <= 0) return 0;
    double q = r / h;
    if (q >= 2.0) return 0;
    double sigma = 10.0 / (7.0 * M_PI * h * h);
    if (q < 1.0) return sigma * (1.0 - 1.5*q*q + 0.75*q*q*q);
    double tq = 2.0 - q;
    return sigma * 0.25 * tq*tq*tq;
}

__device__ __forceinline__
double dev_dWdr_cspline_2d(double r, double h)
{
    if (h <= 0) return 0;
    double q = r / h;
    if (q >= 2.0) return 0;
    double sigma = 10.0 / (7.0 * M_PI * h * h);
    if (q < 1.0) return (sigma/h) * (-3.0*q + 2.25*q*q);
    double tq = 2.0 - q;
    return -(sigma/h) * 0.75 * tq*tq;
}

/* ================================================================
 *  GPU nearest-neighbor kernel: replaces det2d_dpqRK4 for all av_modes.
 *  For each real particle, finds min distance² to any neighbor in the
 *  cell linked list → w2ceil.  Optionally clips w2 ≤ w2ceil.
 *
 *  Uses CellCSR (3×3 cell neighborhood) instead of k-d tree.
 * ================================================================ */
__global__ __launch_bounds__(256, 4)
void nearest_neighbor_kernel(
    const double *__restrict__ px, const double *__restrict__ py,
    const long long *__restrict__ pindx,
    const int *__restrict__ cell_offset,
    const int *__restrict__ particle_index,
    float *__restrict__ w2ceil_out,
    float *__restrict__ w2_out,
    int n_particles, int n_total,
    double cellsize, int mx, int my,
    double xmin, double ymin,
    float kappa)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_particles) return;
    if (pindx[i] == MAX_INDEX_GPU) return;

    double xi = px[i], yi = py[i];
    double invcs = 1.0 / cellsize;
    int cix = (int)((xi - xmin) * invcs);
    int ciy = (int)((yi - ymin) * invcs);
    if (cix < 0) cix = 0; if (cix >= mx) cix = mx - 1;
    if (ciy < 0) ciy = 0; if (ciy >= my) ciy = my - 1;

    double min_dist2 = 1e30;

    for (int jy = ciy - 1; jy <= ciy + 1; jy++) {
        if (jy < 0 || jy >= my) continue;
        for (int jx = cix - 1; jx <= cix + 1; jx++) {
            if (jx < 0 || jx >= mx) continue;
            int cell_id = jx + mx * jy;
            int k_start = cell_offset[cell_id];
            int k_end   = cell_offset[cell_id + 1];
            for (int k = k_start; k < k_end; k++) {
                int j = particle_index[k];
                if (j == i || j < 0 || j >= n_total) continue;
                double dx = px[j] - xi;
                double dy = py[j] - yi;
                double d2 = dx * dx + dy * dy;
                if (d2 > 0 && d2 < min_dist2) min_dist2 = d2;
            }
        }
    }

    w2ceil_out[i] = (float)min_dist2;
    if (kappa == 0.0f)
        w2_out[i] = 0.0f;
    else if (kappa > 0.0f)
        w2_out[i] = fminf(w2_out[i], (float)min_dist2);
}

/* ================================================================
 *  LagMFM density kernel: adaptive h, kernel density, E_inv, gradients,
 *  EOS, NS stress.  1 thread per real particle.
 *
 *  Direct port of updateDenW2Pressure2D_LagMFM (exam.c:2098-2337).
 *  Neighbor search via CellCSR.
 * ================================================================ */
__global__ __launch_bounds__(256, 2)
void lagmfm_density_kernel(
    /* Particle SoA — input */
    const double *__restrict__ px, const double *__restrict__ py,
    const double *__restrict__ pvx, const double *__restrict__ pvy,
    const double *__restrict__ pmass,
    const long long *__restrict__ pindx,
    const float *__restrict__ ppressure_in,
    const float *__restrict__ pvolume_in,
    const float *__restrict__ pie_in,
    const double *__restrict__ palpha_cd,
    /* Cell CSR */
    const int *__restrict__ cell_offset,
    const int *__restrict__ particle_index,
    /* Output arrays */
    float *__restrict__ den_out,
    float *__restrict__ volume_out,
    float *__restrict__ pressure_out,
    float *__restrict__ csound_out,
    double *__restrict__ E_inv_xx_out, double *__restrict__ E_inv_xy_out,
    double *__restrict__ E_inv_yx_out, double *__restrict__ E_inv_yy_out,
    double *__restrict__ h_mfm_out,
    double *__restrict__ gUxx_out, double *__restrict__ gUxy_out,
    double *__restrict__ gUyx_out, double *__restrict__ gUyy_out,
    double *__restrict__ dPdx_out, double *__restrict__ dPdy_out,
    double *__restrict__ dRhodx_out, double *__restrict__ dRhody_out,
    double *__restrict__ divv_out,
    double *__restrict__ tauxx_out, double *__restrict__ tauxy_out,
    double *__restrict__ tauyy_out,
    float *__restrict__ w2ceil_out,
    /* Parameters */
    int n_particles, int n_total,
    double eta, int h_iter_max, double h_tol,
    double cellsize, int mx, int my,
    double xmin, double ymin,
    double Gamma, double nu_phys, int pure_gizmo)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_particles) return;

    /* Skip wall ghost particles */
    if (pindx[i] == MAX_INDEX_GPU) return;

    double xi = px[i], yi = py[i];
    double ivx = pvx[i], ivy = pvy[i];
    double imass = pmass[i];
    double iP = (double)ppressure_in[i];
    double iV = (double)pvolume_in[i];

    double invcs = 1.0 / cellsize;

    /* Initial h from previous volume */
    double Vprev = iV > 0 ? iV : cellsize * cellsize;
    double h_i = eta * sqrt(Vprev);

    /* Cell index of particle i */
    int cix = (int)((xi - xmin) * invcs);
    int ciy = (int)((yi - ymin) * invcs);
    if (cix < 0) cix = 0; if (cix >= mx) cix = mx - 1;
    if (ciy < 0) ciy = 0; if (ciy >= my) ciy = my - 1;

    double n_i = 0;
    double rho_ker = 0;
    double Exx = 0, Exy = 0, Eyy = 0;
    double Sxx_vx = 0, Sxy_vx = 0;
    double Sxx_vy = 0, Sxy_vy = 0;
    double Sxx_P = 0, Sxy_P = 0;
    /* Hopkins (2015) B31 neighbor extrema for limiter */
    double P_min = iP, P_max = iP;
    double vx_min = ivx, vx_max = ivx;
    double vy_min = ivy, vy_max = ivy;
    double rho_min_local = 0, rho_max_local = 0; /* filled below */
    /* Initialize rho bounds after we have rho_i (deferred: use mass*n_prev
     * estimate by reading any neighbor).  Safer: initialize to +inf/-inf
     * and overwrite once we include self or at first neighbor. */
    rho_min_local =  1e300;
    rho_max_local = -1e300;

    /* Newton iteration with bisection fallback (Hopkins 2015). */
    double h_lo = 0, h_hi = 0;
    double dh_prev = 0;
    int use_bisect = 0;

    for (int h_it = 0; h_it < h_iter_max; h_it++) {
        n_i = 0;
        rho_ker = 0;
        Exx = Exy = Eyy = 0;
        Sxx_vx = Sxy_vx = 0;
        Sxx_vy = Sxy_vy = 0;
        Sxx_P = Sxy_P = 0;
        P_min = P_max = iP;
        vx_min = vx_max = ivx;
        vy_min = vy_max = ivy;
        rho_min_local =  1e300;
        rho_max_local = -1e300;

        double two_h = 2.0 * h_i;
        double two_h2 = two_h * two_h;
        int nbr = (int)ceil(two_h * invcs);
        if (nbr < 1) nbr = 1;
        if (nbr > 10) nbr = 10;

        /* Neighbor loop via cell grid */
        for (int jy = ciy - nbr; jy <= ciy + nbr; jy++) {
            if (jy < 0 || jy >= my) continue;
            for (int jx = cix - nbr; jx <= cix + nbr; jx++) {
                if (jx < 0 || jx >= mx) continue;
                int cell_id = jx + mx * jy;
                int k_start = cell_offset[cell_id];
                int k_end   = cell_offset[cell_id + 1];
                for (int k = k_start; k < k_end; k++) {
                    int j = particle_index[k];
                    if (j < 0 || j >= n_total) continue;
                    double dx = px[j] - xi;
                    double dy = py[j] - yi;
                    double r2 = dx * dx + dy * dy;
                    if (r2 < two_h2) {
                        double r = sqrt(r2);
                        double Wij = dev_wendland_c2_2d(r, h_i);
                        n_i += Wij;
                        rho_ker += Wij * pmass[j];
                        Exx += Wij * dx * dx;
                        Exy += Wij * dx * dy;
                        Eyy += Wij * dy * dy;
                        double dvx = pvx[j] - ivx;
                        double dvy = pvy[j] - ivy;
                        double dP  = (double)ppressure_in[j] - iP;
                        Sxx_vx += Wij * dvx * dx;
                        Sxy_vx += Wij * dvx * dy;
                        Sxx_vy += Wij * dvy * dx;
                        Sxy_vy += Wij * dvy * dy;
                        Sxx_P  += Wij * dP  * dx;
                        Sxy_P  += Wij * dP  * dy;
                        /* Neighbor extrema for B31 limiter.
                         * ρ_j ≈ m_j / V_j_prev (lagged, consistent with pden). */
                        double Pj = (double)ppressure_in[j];
                        if (Pj < P_min) P_min = Pj; if (Pj > P_max) P_max = Pj;
                        if (pvx[j] < vx_min) vx_min = pvx[j];
                        if (pvx[j] > vx_max) vx_max = pvx[j];
                        if (pvy[j] < vy_min) vy_min = pvy[j];
                        if (pvy[j] > vy_max) vy_max = pvy[j];
                        double rho_j_est = pmass[j] / (double)pvolume_in[j];
                        if (rho_j_est < rho_min_local) rho_min_local = rho_j_est;
                        if (rho_j_est > rho_max_local) rho_max_local = rho_j_est;
                    }
                }
            }
        }

        if (n_i <= 0) break;  /* isolated */

        double V_new = 1.0 / n_i;
        double h_target = eta * sqrt(V_new);
        double R = h_i - h_target;
        if (R > 0) h_hi = h_i; else h_lo = h_i;

        double h_new;
        if (use_bisect && h_lo > 0 && h_hi > 0) {
            h_new = 0.5*(h_lo + h_hi);
        } else {
            h_new = h_target;
            double dh = h_new - h_i;
            if (h_it > 0 && dh*dh_prev < 0 &&
                fabs(dh) > 0.5*fabs(dh_prev) &&
                h_lo > 0 && h_hi > 0) {
                use_bisect = 1;
                h_new = 0.5*(h_lo + h_hi);
            }
            dh_prev = dh;
        }
        double adh = fabs(h_new - h_i);
        h_i = h_new;
        if (adh < h_tol * h_i) break;
    }

    /* Isolated particle — retain previous state */
    if (n_i <= 0) return;

    /* Finalize density / volume */
    double Vi = 1.0 / n_i;

    /* Invert E to get matrix-weighted gradient operator */
    double det = Exx * Eyy - Exy * Exy;
    double adet = fabs(det);
    double nrm2 = Exx * Exx + 2.0 * Exy * Exy + Eyy * Eyy;
    double tol = 1.0e-14 * nrm2 + 1.0e-30;

    double Einv_xx = 0, Einv_xy = 0, Einv_yx = 0, Einv_yy = 0;
    if (adet >= tol) {
        double invdet = 1.0 / det;
        Einv_xx =  Eyy * invdet;
        Einv_xy = -Exy * invdet;
        Einv_yx = -Exy * invdet;
        Einv_yy =  Exx * invdet;
    }

    /* Gradients: ∇Q|_i^α = Σ_β (E_inv)^{αβ} · S_β */
    double gux = Einv_xx * Sxx_vx + Einv_xy * Sxy_vx;
    double guy = Einv_yx * Sxx_vx + Einv_yy * Sxy_vx;
    double gvx = Einv_xx * Sxx_vy + Einv_xy * Sxy_vy;
    double gvy = Einv_yx * Sxx_vy + Einv_yy * Sxy_vy;
    double dpx = Einv_xx * Sxx_P  + Einv_xy * Sxy_P;
    double dpy = Einv_yx * Sxx_P  + Einv_yy * Sxy_P;

    /* Density and EOS (mass-weighted: ρ_i = Σ_j m_j W_ij) */
    double rho_val = rho_ker;
    double ie_val = (double)pie_in[i];
    double P_val = ie_val / Vi * (Gamma - 1.0);
    if (P_val <= 0) {
        P_val = 1e-6;
        ie_val = P_val * Vi / (Gamma - 1.0);
    }
    double cs_val = sqrt(Gamma * P_val / rho_val);

    /* Post-Newton ρ gradient pass (needs final n_i → ρ_i).
     * Re-walk neighbors with converged h_i and accumulate ρ-moments. */
    double Sxx_Rho_f = 0, Sxy_Rho_f = 0;
    double rho_min_f = rho_val, rho_max_f = rho_val;
    {
        double two_h = 2.0 * h_i;
        double two_h2 = two_h * two_h;
        int nbr = (int)ceil(two_h * invcs);
        if (nbr < 1) nbr = 1;
        if (nbr > 10) nbr = 10;
        for (int jy = ciy - nbr; jy <= ciy + nbr; jy++) {
            if (jy < 0 || jy >= my) continue;
            for (int jx = cix - nbr; jx <= cix + nbr; jx++) {
                if (jx < 0 || jx >= mx) continue;
                int cell_id = jx + mx * jy;
                int ks = cell_offset[cell_id];
                int ke = cell_offset[cell_id + 1];
                for (int k = ks; k < ke; k++) {
                    int j = particle_index[k];
                    if (j < 0 || j >= n_total) continue;
                    double dx = px[j] - xi;
                    double dy = py[j] - yi;
                    double r2 = dx * dx + dy * dy;
                    if (r2 < two_h2) {
                        double rr = sqrt(r2);
                        double Wij = dev_wendland_c2_2d(rr, h_i);
                        double rho_j = pmass[j] / (double)pvolume_in[j];
                        double dRho = rho_j - rho_val;
                        Sxx_Rho_f += Wij * dRho * dx;
                        Sxy_Rho_f += Wij * dRho * dy;
                        if (rho_j < rho_min_f) rho_min_f = rho_j;
                        if (rho_j > rho_max_f) rho_max_f = rho_j;
                    }
                }
            }
        }
    }
    double dRhodx_val = Einv_xx * Sxx_Rho_f + Einv_xy * Sxy_Rho_f;
    double dRhody_val = Einv_yx * Sxx_Rho_f + Einv_yy * Sxy_Rho_f;

    /* Hopkins (2015) B31 gradient-scaling limiter with ADAPTIVE β_i.
     * β_i = max(β_min, β_max · min(1, N_crit/N_cond)), N_cond = λ_max/λ_min(E).
     * In degenerate neighborhoods (rolled-up vortex edges) β drops to β_min=1
     * → tight limiter → suppresses grid-noise sub-mode amplification.
     * Matches Hopkins Appendix B (β_min=1, β_max=2, N_crit=10). */
    double Ncond_i;
    {
        double tr = Exx + Eyy;
        double disc = (Exx - Eyy) * (Exx - Eyy) + 4.0 * Exy * Exy;
        double sqrtdisc = sqrt(disc > 0 ? disc : 0);
        double lam_max = 0.5 * (tr + sqrtdisc);
        double lam_min = 0.5 * (tr - sqrtdisc);
        if (det <= 0 || lam_min <= 0) Ncond_i = 1.0e30;
        else                          Ncond_i = lam_max / lam_min;
    }
    {
        double recon_scale = h_i;
        double N_crit = 10.0;
        double beta_lim = fmax(1.0, 2.0 * fmin(1.0, N_crit / Ncond_i));
        double gmag_P = sqrt(dpx*dpx + dpy*dpy) * recon_scale;
        if (gmag_P > 0) {
            double a_up = beta_lim * (P_max - P_val) / gmag_P;
            double a_dn = beta_lim * (P_val - P_min) / gmag_P;
            double alpha = fmin(1.0, fmin(a_up, a_dn));
            if (alpha < 0) alpha = 0;
            dpx *= alpha; dpy *= alpha;
        }
        double gmag_vx = sqrt(gux*gux + guy*guy) * recon_scale;
        if (gmag_vx > 0) {
            double a_up = beta_lim * (vx_max - ivx) / gmag_vx;
            double a_dn = beta_lim * (ivx - vx_min) / gmag_vx;
            double alpha = fmin(1.0, fmin(a_up, a_dn));
            if (alpha < 0) alpha = 0;
            gux *= alpha; guy *= alpha;
        }
        double gmag_vy = sqrt(gvx*gvx + gvy*gvy) * recon_scale;
        if (gmag_vy > 0) {
            double a_up = beta_lim * (vy_max - ivy) / gmag_vy;
            double a_dn = beta_lim * (ivy - vy_min) / gmag_vy;
            double alpha = fmin(1.0, fmin(a_up, a_dn));
            if (alpha < 0) alpha = 0;
            gvx *= alpha; gvy *= alpha;
        }
        double gmag_rho = sqrt(dRhodx_val*dRhodx_val + dRhody_val*dRhody_val) * recon_scale;
        if (gmag_rho > 0) {
            double a_up = beta_lim * (rho_max_f - rho_val) / gmag_rho;
            double a_dn = beta_lim * (rho_val - rho_min_f) / gmag_rho;
            double alpha = fmin(1.0, fmin(a_up, a_dn));
            if (alpha < 0) alpha = 0;
            dRhodx_val *= alpha; dRhody_val *= alpha;
        }
    }
    double divv_val = gux + gvy;

    /* NS stress: zeroed in pure-GIZMO mode (inviscid Euler).
     * Otherwise (3D-slab trace convention, matches CPU). */
    double tau_xx, tau_xy, tau_yy;
    if (pure_gizmo) {
        tau_xx = 0; tau_xy = 0; tau_yy = 0;
        (void)nu_phys; (void)palpha_cd;
    } else {
        /* NS-stress: physical viscosity only. CD10 NS-stress path removed (caused
         * Nx² K-runaway via h·c_s factor). CD10 still acts as viscous pressure. */
        (void)palpha_cd;
        double nu = (nu_phys > 0) ? nu_phys : 0.0;
        tau_xx = -nu * rho_val * (2.0 * gux - (2.0 / 3.0) * divv_val);
        tau_xy = -nu * rho_val * (guy + gvx);
        tau_yy = -nu * rho_val * (2.0 * gvy - (2.0 / 3.0) * divv_val);
    }

    /* Write outputs */
    den_out[i]      = (float)rho_val;
    volume_out[i]   = (float)Vi;
    pressure_out[i] = (float)P_val;
    csound_out[i]   = (float)cs_val;
    E_inv_xx_out[i] = Einv_xx;
    E_inv_xy_out[i] = Einv_xy;
    E_inv_yx_out[i] = Einv_yx;
    E_inv_yy_out[i] = Einv_yy;
    h_mfm_out[i]    = h_i;
    gUxx_out[i]     = gux;
    gUxy_out[i]     = guy;
    gUyx_out[i]     = gvx;
    gUyy_out[i]     = gvy;
    dPdx_out[i]     = dpx;
    dPdy_out[i]     = dpy;
    dRhodx_out[i]   = dRhodx_val;
    dRhody_out[i]   = dRhody_val;
    divv_out[i]     = divv_val;
    tauxx_out[i]    = tau_xx;
    tauxy_out[i]    = tau_xy;
    tauyy_out[i]    = tau_yy;

    /* Nearest-neighbor distance² for w2ceil (replaces det2d_dpqRK4 on GPU).
     * Reuse the final Newton-converged neighbor loop to find min distance. */
    {
        double min_dist2 = 1e30;
        double two_h = 2.0 * h_i;
        double two_h2 = two_h * two_h;
        int nbr = (int)ceil(two_h * invcs);
        if (nbr < 1) nbr = 1;
        if (nbr > 10) nbr = 10;

        for (int jy = ciy - nbr; jy <= ciy + nbr; jy++) {
            if (jy < 0 || jy >= my) continue;
            for (int jx = cix - nbr; jx <= cix + nbr; jx++) {
                if (jx < 0 || jx >= mx) continue;
                int cell_id = jx + mx * jy;
                int k_start = cell_offset[cell_id];
                int k_end   = cell_offset[cell_id + 1];
                for (int k = k_start; k < k_end; k++) {
                    int j = particle_index[k];
                    if (j == i || j < 0 || j >= n_total) continue;
                    double dxx = px[j] - xi;
                    double dyy = py[j] - yi;
                    double d2  = dxx * dxx + dyy * dyy;
                    if (d2 > 0 && d2 < min_dist2) min_dist2 = d2;
                }
            }
        }
        w2ceil_out[i] = (float)min_dist2;
    }
}

/* ================================================================
 *  LagMFM force kernel: MFM effective faces, MUSCL+clamp, HLLC
 *  face-rest-frame, Monaghan AV, NS traction, CFL.
 *
 *  Direct port of getAccVoro2D_LagMFM (exam.c:3129-3525).
 *  1 thread per real particle, neighbor search via CellCSR.
 * ================================================================ */
__global__ __launch_bounds__(256, 2)
void lagmfm_force_kernel(
    /* Particle SoA — input */
    const double *__restrict__ px, const double *__restrict__ py,
    const double *__restrict__ pvx, const double *__restrict__ pvy,
    const double *__restrict__ pmass,
    const long long *__restrict__ pindx,
    const float *__restrict__ pden, const float *__restrict__ ppressure,
    const float *__restrict__ pcsound, const float *__restrict__ pvolume,
    const float *__restrict__ pw2,
    const double *__restrict__ pE_inv_xx, const double *__restrict__ pE_inv_xy,
    const double *__restrict__ pE_inv_yx, const double *__restrict__ pE_inv_yy,
    const double *__restrict__ ph_mfm,
    const double *__restrict__ pgUxx, const double *__restrict__ pgUxy,
    const double *__restrict__ pgUyx, const double *__restrict__ pgUyy,
    const double *__restrict__ pdPdx, const double *__restrict__ pdPdy,
    const double *__restrict__ pdRhodx, const double *__restrict__ pdRhody,
    const double *__restrict__ ptauxx, const double *__restrict__ ptauxy,
    const double *__restrict__ ptauyy,
    const double *__restrict__ palpha_cd,
    /* Cell CSR */
    const int *__restrict__ cell_offset,
    const int *__restrict__ particle_index,
    /* Output arrays */
    double *__restrict__ ax_out, double *__restrict__ ay_out,
    float *__restrict__ die_out, float *__restrict__ dt_out,
    double *__restrict__ vsig_max_out,
    /* Parameters */
    int n_particles, int n_total,
    double cellsize, int mx, int my,
    double xmin, double ymin,
    double Courant, double Gamma,
    double alphavis, double betavis, double etavis, double epsvis,
    double nu_phys, double eta, int pure_gizmo)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_particles) return;

    if (pindx[i] == MAX_INDEX_GPU) return;

    /* Load particle i */
    double xi = px[i], yi = py[i];
    double ivx = pvx[i], ivy = pvy[i];
    double imass = pmass[i];
    double iP   = (double)ppressure[i];
    double irho = (double)pden[i];
    double ics  = (double)pcsound[i];
    double iV   = (double)pvolume[i];
    double iw2  = (double)pw2[i];
    double Ei_xx = pE_inv_xx[i], Ei_xy = pE_inv_xy[i];
    double Ei_yx = pE_inv_yx[i], Ei_yy = pE_inv_yy[i];
    double h_i = ph_mfm[i];
    if (h_i <= 0) h_i = eta * sqrt(iV > 0 ? iV : cellsize * cellsize);

    double i_gUxx = pgUxx[i], i_gUxy = pgUxy[i];
    double i_gUyx = pgUyx[i], i_gUyy = pgUyy[i];
    double i_dPdx = pdPdx[i], i_dPdy = pdPdy[i];
    double i_dRhodx = pdRhodx[i], i_dRhody = pdRhody[i];
    double i_tauxx = ptauxx[i], i_tauxy = ptauxy[i], i_tauyy = ptauyy[i];

    double invcs = 1.0 / cellsize;

    int cix = (int)((xi - xmin) * invcs);
    int ciy = (int)((yi - ymin) * invcs);
    if (cix < 0) cix = 0; if (cix >= mx) cix = mx - 1;
    if (ciy < 0) ciy = 0; if (ciy >= my) ciy = my - 1;

    double fx = 0, fy = 0, die = 0;
    float  my_dt = 1.0e10f;
    double my_vsig_max = 0;

    double two_h_i = 2.0 * h_i;
    int nbr = (int)ceil(two_h_i * invcs) + 1;
    if (nbr < 1) nbr = 1;
    if (nbr > 10) nbr = 10;

    for (int jy = ciy - nbr; jy <= ciy + nbr; jy++) {
        if (jy < 0 || jy >= my) continue;
        for (int jx = cix - nbr; jx <= cix + nbr; jx++) {
            if (jx < 0 || jx >= mx) continue;
            int cell_id = jx + mx * jy;
            int k_start = cell_offset[cell_id];
            int k_end   = cell_offset[cell_id + 1];
            for (int k = k_start; k < k_end; k++) {
                int j = particle_index[k];
                if (j < 0 || j >= n_total || j == i) continue;
                int jwall = (pindx[j] == MAX_INDEX_GPU);

                double dx = px[j] - xi;
                double dy = py[j] - yi;
                double r2 = dx * dx + dy * dy;

                double h_j = ph_mfm[j];
                if (h_j <= 0) h_j = h_i;
                double h_max = fmax(h_i, h_j);
                double two_h_max = 2.0 * h_max;
                if (r2 > two_h_max * two_h_max) continue;
                double r = sqrt(r2);

                /* Effective face area vector.  GIZMO volume-disparity test:
                 * |V_i-V_j|/min(V)/NUMDIMS > 1.25 (2D: 2.5) → symmetric weight. */
                double W_ij = dev_wendland_c2_2d(r, h_i);
                double W_ji = dev_wendland_c2_2d(r, h_j);
                double Vj  = (double)pvolume[j];
                double Vi2, Vj2;
                {
                    double Vmin = (iV < Vj) ? iV : Vj;
                    double Vabs = iV - Vj; if (Vabs < 0) Vabs = -Vabs;
                    if (Vmin > 0 && (Vabs/Vmin) > 2.5) {
                        double denom = iV*W_ij + Vj*W_ji;
                        double wcomm = (denom > 0) ? (iV*Vj*(W_ij+W_ji)/denom) : iV;
                        Vi2 = Vj2 = wcomm;
                    } else {
                        Vi2 = iV;
                        Vj2 = Vj;
                    }
                }

                double Ax_i = Vi2 * W_ij * (Ei_xx * dx + Ei_xy * dy);
                double Ay_i = Vi2 * W_ij * (Ei_yx * dx + Ei_yy * dy);
                double Ax_j, Ay_j;
                if (jwall) {
                    Ax_j = 0; Ay_j = 0;
                } else {
                    double Ej_xx = pE_inv_xx[j], Ej_xy = pE_inv_xy[j];
                    double Ej_yx = pE_inv_yx[j], Ej_yy = pE_inv_yy[j];
                    Ax_j = Vj2 * W_ji * (Ej_xx * dx + Ej_xy * dy);
                    Ay_j = Vj2 * W_ji * (Ej_yx * dx + Ej_yy * dy);
                }
                double Ax = Ax_i + Ax_j;
                double Ay = Ay_i + Ay_j;
                double Amag2 = Ax * Ax + Ay * Ay;
                if (Amag2 < 1e-30) continue;
                double Amag = sqrt(Amag2);
                double nx_hat = Ax / Amag;
                double ny_hat = Ay / Amag;

                double jvx = pvx[j], jvy = pvy[j];
                double jP  = (double)ppressure[j];
                double jrho = (double)pden[j];
                double jcs  = (double)pcsound[j];

                /* GIZMO default equal-weight midpoint (hydro_core_meshless.h:86-102
                 * with s_star_ij=0): reconstruct from particle center to 0.5·r. */
                double hdx = 0.5 * dx, hdy = 0.5 * dy;
                double pL   = iP   + i_dPdx * hdx   + i_dPdy * hdy;
                double pR   = jP   - pdPdx[j] * hdx - pdPdy[j] * hdy;
                double vxL  = ivx  + i_gUxx * hdx   + i_gUxy * hdy;
                double vxR  = jvx  - pgUxx[j] * hdx - pgUxy[j] * hdy;
                double vyL  = ivy  + i_gUyx * hdx   + i_gUyy * hdy;
                double vyR  = jvy  - pgUyx[j] * hdx - pgUyy[j] * hdy;
                double rhoL = irho + i_dRhodx * hdx + i_dRhody * hdy;
                double rhoR = jrho - pdRhodx[j]*hdx - pdRhody[j]*hdy;

                /* GIZMO asymmetric pair-wise limiter (reimann.h:188-237).
                 * ψ_1 = 0.5 (fac_minmax), ψ_2 = 0.375 (fac_meddev).
                 * L = from i, R = from j (our convention; role-swap for Q_i vs Q_j ordering). */
                if (!jwall) {
                    #define GIZMO_PAIRLIM(QL, QR, Qi, Qj) do {                  \
                        if ((Qi) != (Qj)) {                                     \
                            double _Qmed = 0.5 * ((Qi) + (Qj));                 \
                            double _Qmin, _Qmax;                                 \
                            if ((Qi) < (Qj)) { _Qmin=(Qi); _Qmax=(Qj); }         \
                            else             { _Qmin=(Qj); _Qmax=(Qi); }         \
                            double _dQ = _Qmax - _Qmin;                          \
                            double _f1 = 0.5   * _dQ;                            \
                            double _f2 = 0.375 * _dQ;                            \
                            double _Qmax_eff = _Qmax + _f1;                      \
                            double _Qmin_eff = _Qmin - _f1;                      \
                            double _Qmed_max = _Qmed + _f2;                      \
                            double _Qmed_min = _Qmed - _f2;                      \
                            if (_Qmed_max > _Qmax_eff) _Qmed_max = _Qmax_eff;    \
                            if (_Qmed_min < _Qmin_eff) _Qmed_min = _Qmin_eff;    \
                            if ((Qi) < (Qj)) {                                   \
                                if ((QL) < _Qmin_eff) (QL) = _Qmin_eff;          \
                                if ((QL) > _Qmed_max) (QL) = _Qmed_max;          \
                                if ((QR) > _Qmax_eff) (QR) = _Qmax_eff;          \
                                if ((QR) < _Qmed_min) (QR) = _Qmed_min;          \
                            } else {                                             \
                                if ((QL) > _Qmax_eff) (QL) = _Qmax_eff;          \
                                if ((QL) < _Qmed_min) (QL) = _Qmed_min;          \
                                if ((QR) < _Qmin_eff) (QR) = _Qmin_eff;          \
                                if ((QR) > _Qmed_max) (QR) = _Qmed_max;          \
                            }                                                    \
                        }                                                        \
                    } while(0)
                    GIZMO_PAIRLIM(pL,   pR,   iP,   jP);
                    GIZMO_PAIRLIM(vxL,  vxR,  ivx,  jvx);
                    GIZMO_PAIRLIM(vyL,  vyR,  ivy,  jvy);
                    GIZMO_PAIRLIM(rhoL, rhoR, irho, jrho);
                    #undef GIZMO_PAIRLIM
                }

                /* Wall: skip reconstruction on both sides */
                if (jwall) {
                    pL   = iP;   vxL = ivx; vyL = ivy; rhoL = irho;
                    pR   = jP;   vxR = jvx; vyR = jvy; rhoR = jrho;
                }
                /* Safety floor for reconstructed ρ, P */
                if (rhoL <= 0) rhoL = irho;
                if (rhoR <= 0) rhoR = jrho;
                if (pL   <= 0) pL   = iP;
                if (pR   <= 0) pR   = jP;

                double vnL_lab = vxL * nx_hat + vyL * ny_hat;
                double vnR_lab = vxR * nx_hat + vyR * ny_hat;
                double cL = sqrt(Gamma * pL / rhoL);
                double cR = sqrt(Gamma * pR / rhoR);

                /* GIZMO default equal-weight face velocity. */
                double wx = 0.5 * (ivx + jvx);
                double wy = 0.5 * (ivy + jvy);
                double wn = wx * nx_hat + wy * ny_hat;

                double pstar, vnstar_lab;
                if (jwall) {
                    pstar = 0.5 * (pL + pR);
                    vnstar_lab = 0;
                } else {
                    dev_hllc_face_2d_rest_frame(
                        rhoL, pL, vnL_lab, cL,
                        rhoR, pR, vnR_lab, cR,
                        wn, Gamma, 0, &pstar, &vnstar_lab);
                }

                /* Non-GIZMO: Monaghan AV add-on.  Pure mode skips. */
                if (!pure_gizmo && alphavis > 0 && !jwall) {
                    double uij_x = jvx - ivx;
                    double uij_y = jvy - ivy;
                    double er_x = dx / r, er_y = dy / r;
                    double rvel = er_x * uij_x + er_y * uij_y;
                    if (rvel < 0) {
                        double wcomp = sqrt(iw2) + sqrt((double)pw2[j]);
                        double scaleFactor = (wcomp > etavis ? wcomp : etavis);
                        double drampScale = r / scaleFactor;
                        double mu = rvel /* Voronoi-AV-reform: no length regularization */;
                        double meanden = 0.5 * (irho + jrho);
                        double meanCs  = 0.5 * (ics + jcs);
                        pstar += (-alphavis * meanCs * mu + betavis * mu * mu) * meanden;
                    }
                }

                /* Non-GIZMO: NS viscous traction. */
                double tau_A_x = 0, tau_A_y = 0;
                if (!pure_gizmo && !jwall) {
                    double txx_f = 0.5 * (i_tauxx + ptauxx[j]);
                    double txy_f = 0.5 * (i_tauxy + ptauxy[j]);
                    double tyy_f = 0.5 * (i_tauyy + ptauyy[j]);
                    tau_A_x = -(txx_f * Ax + txy_f * Ay);
                    tau_A_y = -(txy_f * Ax + tyy_f * Ay);
                }

                /* Force */
                fx += -pstar * Ax + tau_A_x;
                fy += -pstar * Ay + tau_A_y;

                /* Full face velocity: Hopkins (2015) §2.5 passive tangential
                 * advection (upwind across contact wave). */
                double tx_hat = -ny_hat;
                double ty_hat =  nx_hat;
                double vtL = vxL * tx_hat + vyL * ty_hat;
                double vtR = vxR * tx_hat + vyR * ty_hat;
                double vt_up = (vnstar_lab - wn) >= 0 ? vtL : vtR;
                double vfx = vnstar_lab * nx_hat + vt_up * tx_hat;
                double vfy = vnstar_lab * ny_hat + vt_up * ty_hat;

                if (jwall) {
                    double vti = ivx * tx_hat + ivy * ty_hat;
                    vfx = vti * tx_hat;
                    vfy = vti * ty_hat;
                }

                /* Energy */
                die += -pstar * ((vfx - ivx) * Ax + (vfy - ivy) * Ay);
                die += tau_A_x * (vfx - ivx) + tau_A_y * (vfy - ivy);

                /* CFL (Hopkins 2015, eq. 9): dt = C · h_i / vsig */
                double dvx_sig = jvx - ivx;
                double dvy_sig = jvy - ivy;
                double er_x = dx / r, er_y = dy / r;
                double VdotR = dvx_sig * er_x + dvy_sig * er_y;
                double vsig = cL + cR - fmin(0.0, VdotR);
                float dt_face = (float)(Courant * h_i / vsig);
                if (dt_face < my_dt) my_dt = dt_face;

                if (!pure_gizmo && !jwall && vsig > my_vsig_max)
                    my_vsig_max = vsig;
            }
        }
    }

    /* Viscous CFL (NS + CD10; non-GIZMO).  Pure mode skips. */
    if (!pure_gizmo) {
        double h_eff = sqrt(iV);
        double nu_cd = palpha_cd[i] * h_eff * ics;
        double nu_eff = fmax(nu_phys, nu_cd);
        if (nu_eff > 0) {
            float dt_visc = (float)(0.5 * h_eff * h_eff / nu_eff);
            if (dt_visc < my_dt) my_dt = dt_visc;
        }
    }

    /* Write outputs */
    ax_out[i]       = fx / imass;
    ay_out[i]       = fy / imass;
    die_out[i]      = (float)die;
    dt_out[i]       = my_dt;
    vsig_max_out[i] = my_vsig_max;
}

/* ================================================================
 *  LagMFM kernel launch wrappers
 * ================================================================ */

extern "C"
void gpu_download_w2ceil(GPUContext *ctx, int n)
{
    cudaStream_t s = (cudaStream_t)ctx->stream;
    size_t sf = n * sizeof(float);
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.w2ceil, ctx->d_parts.w2ceil, sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.w2,     ctx->d_parts.w2,     sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaStreamSynchronize(s));
}

/* ================================================================
 *  Option C: OST tree-walk nearest-neighbor kernel.
 *
 *  Mirrors CPU find_GNearest (gnnost.c:320) using a flat index-based
 *  representation of the OST tree.  Stack-less while-loop with
 *  bounding-sphere pruning (dist - nodesize > mindist → skip subtree).
 *
 *  Robust to arbitrary cell expansion: prunes by tree topology, not by
 *  fixed-radius cell stencil.  Solves the Sedov rarefaction failure mode
 *  where sqrt(Vcell) > 3*cellsize breaks CellCSR 3×3 search.
 *
 *  1 thread / real particle.  Each thread independently walks the tree.
 * ================================================================ */
__global__ __launch_bounds__(256, 4)
void nearest_neighbor_tree_kernel(
    /* Query particles */
    const double *__restrict__ px, const double *__restrict__ py,
    const long long *__restrict__ pindx,
    /* Tree SoA */
    const int   *__restrict__ tree_type,
    const int   *__restrict__ tree_sibling,
    const int   *__restrict__ tree_daughter,
    const int   *__restrict__ tree_ptl_idx,
    const float *__restrict__ tree_cx,
    const float *__restrict__ tree_cy,
    const float *__restrict__ tree_size,
    int root_idx, int n_tree_nodes,
    /* Output */
    float *__restrict__ w2ceil_out,
    float *__restrict__ w2_out,
    int n_particles,
    float kappa)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_particles) return;
    if (pindx[i] == MAX_INDEX_GPU) return;

    float qx = (float)px[i];
    float qy = (float)py[i];
    /* Self skip: when ptl[] is built by det2d_dpqRK4_GPU, real particle i
     * occupies ptl[i] with indx=i (matches SoA index). */
    int self_ptl_idx = i;

    float mindist2 = 1.0e23f;
    float mindist  = 1.0e23f;

    int cur = root_idx;
    /* Safety bound on iterations to guard against any cycles (should not
     * happen by construction).  Use 4*n_tree_nodes as a generous cap. */
    int iter = 0;
    int iter_max = 4 * n_tree_nodes;

    while (cur >= 0 && cur < n_tree_nodes && iter < iter_max) {
        ++iter;
        int type = tree_type[cur];
        if (type == 0) {
            /* TREE node: test bounding sphere */
            float cx = tree_cx[cur];
            float cy = tree_cy[cur];
            float dx = qx - cx;
            float dy = qy - cy;
            float dist = sqrtf(dx*dx + dy*dy);
            if (dist - tree_size[cur] > mindist) {
                cur = tree_sibling[cur];     /* prune subtree */
            } else {
                cur = tree_daughter[cur];    /* descend */
            }
        } else {
            /* PTL leaf */
            int pidx = tree_ptl_idx[cur];
            if (pidx != self_ptl_idx) {
                float cx = tree_cx[cur];
                float cy = tree_cy[cur];
                float dx = qx - cx;
                float dy = qy - cy;
                float d2 = dx*dx + dy*dy;
                if (d2 < mindist2 && d2 > 0.0f) {
                    mindist2 = d2;
                    mindist  = sqrtf(d2);
                }
            }
            cur = tree_sibling[cur];
        }
    }

    w2ceil_out[i] = mindist2;
    if (kappa == 0.0f)
        w2_out[i] = 0.0f;
    else if (kappa > 0.0f)
        w2_out[i] = fminf(w2_out[i], mindist2);
}

extern "C"
void gpu_launch_nearest_neighbor_tree(GPUContext *ctx, int n_particles,
                                       int n_total, float kappa)
{
    int blockSize = 256;
    int gridSize = (n_particles + blockSize - 1) / blockSize;
    cudaStream_t s = (cudaStream_t)ctx->stream;

    nearest_neighbor_tree_kernel<<<gridSize, blockSize, 0, s>>>(
        ctx->d_parts.x, ctx->d_parts.y,
        ctx->d_parts.indx,
        ctx->d_tree_type,
        ctx->d_tree_sibling,
        ctx->d_tree_daughter,
        ctx->d_tree_ptl_idx,
        ctx->d_tree_cx, ctx->d_tree_cy,
        ctx->d_tree_size,
        ctx->tree_root_idx, ctx->n_tree_nodes,
        ctx->d_parts.w2ceil, ctx->d_parts.w2,
        n_particles,
        kappa);

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(s));
}

extern "C"
void gpu_launch_nearest_neighbor(GPUContext *ctx, int n_particles, int n_total,
                                  const GPUPhysicsParams *params, float kappa)
{
    int blockSize = 256;
    int gridSize = (n_particles + blockSize - 1) / blockSize;
    cudaStream_t s = (cudaStream_t)ctx->stream;

    nearest_neighbor_kernel<<<gridSize, blockSize, 0, s>>>(
        ctx->d_parts.x, ctx->d_parts.y,
        ctx->d_parts.indx,
        ctx->d_cells.cell_offset, ctx->d_cells.particle_index,
        ctx->d_parts.w2ceil, ctx->d_parts.w2,
        n_particles, n_total,
        params->cellsize, params->mx, params->my,
        params->xmin, params->ymin,
        kappa);

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(s));
}

extern "C"
void gpu_launch_lagmfm_density_kernel(GPUContext *ctx, int n_particles,
                                       int n_total,
                                       const GPUPhysicsParams *params)
{
    int blockSize = 256;
    int gridSize = (n_particles + blockSize - 1) / blockSize;
    cudaStream_t s = (cudaStream_t)ctx->stream;

    lagmfm_density_kernel<<<gridSize, blockSize, 0, s>>>(
        ctx->d_parts.x, ctx->d_parts.y,
        ctx->d_parts.vx, ctx->d_parts.vy,
        ctx->d_parts.mass, ctx->d_parts.indx,
        ctx->d_parts.pressure, ctx->d_parts.volume,
        ctx->d_parts.ie, ctx->d_parts.alpha_cd,
        /* Cell CSR */
        ctx->d_cells.cell_offset, ctx->d_cells.particle_index,
        /* Output */
        ctx->d_parts.den, ctx->d_parts.volume,
        ctx->d_parts.pressure, ctx->d_parts.csound,
        ctx->d_parts.E_inv_xx, ctx->d_parts.E_inv_xy,
        ctx->d_parts.E_inv_yx, ctx->d_parts.E_inv_yy,
        ctx->d_parts.h_mfm,
        ctx->d_parts.gUxx, ctx->d_parts.gUxy,
        ctx->d_parts.gUyx, ctx->d_parts.gUyy,
        ctx->d_parts.dPdx, ctx->d_parts.dPdy,
        ctx->d_parts.dRhodx, ctx->d_parts.dRhody,
        ctx->d_parts.divv,
        ctx->d_parts.tauxx, ctx->d_parts.tauxy, ctx->d_parts.tauyy,
        ctx->d_parts.w2ceil,
        /* Parameters */
        n_particles, n_total,
        params->lagmfm_eta, params->lagmfm_h_iter, params->lagmfm_h_tol,
        params->cellsize, params->mx, params->my,
        params->xmin, params->ymin,
        params->Gamma, params->nu_phys, params->lagmfm_pure_gizmo);

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(s));
}

extern "C"
double gpu_launch_lagmfm_force_kernel(GPUContext *ctx, int n_particles,
                                       int n_total,
                                       const GPUPhysicsParams *params)
{
    int blockSize = 256;
    int gridSize = (n_particles + blockSize - 1) / blockSize;
    cudaStream_t s = (cudaStream_t)ctx->stream;

    lagmfm_force_kernel<<<gridSize, blockSize, 0, s>>>(
        ctx->d_parts.x, ctx->d_parts.y,
        ctx->d_parts.vx, ctx->d_parts.vy,
        ctx->d_parts.mass, ctx->d_parts.indx,
        ctx->d_parts.den, ctx->d_parts.pressure,
        ctx->d_parts.csound, ctx->d_parts.volume,
        ctx->d_parts.w2,
        ctx->d_parts.E_inv_xx, ctx->d_parts.E_inv_xy,
        ctx->d_parts.E_inv_yx, ctx->d_parts.E_inv_yy,
        ctx->d_parts.h_mfm,
        ctx->d_parts.gUxx, ctx->d_parts.gUxy,
        ctx->d_parts.gUyx, ctx->d_parts.gUyy,
        ctx->d_parts.dPdx, ctx->d_parts.dPdy,
        ctx->d_parts.dRhodx, ctx->d_parts.dRhody,
        ctx->d_parts.tauxx, ctx->d_parts.tauxy, ctx->d_parts.tauyy,
        ctx->d_parts.alpha_cd,
        /* Cell CSR */
        ctx->d_cells.cell_offset, ctx->d_cells.particle_index,
        /* Output */
        ctx->d_parts.ax_out, ctx->d_parts.ay_out,
        ctx->d_parts.die_out, ctx->d_parts.dt_out,
        ctx->d_parts.vsig_max_out,
        /* Parameters */
        n_particles, n_total,
        params->cellsize, params->mx, params->my,
        params->xmin, params->ymin,
        params->Courant, params->Gamma,
        params->alphavis, params->betavis,
        params->etavis, params->epsvis,
        params->nu_phys, params->lagmfm_eta, params->lagmfm_pure_gizmo);

    CUDA_CHECK(cudaGetLastError());

    /* CFL reduction */
    cub::DeviceReduce::Min(ctx->d_cub_temp, ctx->cub_temp_bytes,
                           ctx->d_parts.dt_out,
                           (float *)ctx->d_cub_min_out,
                           n_particles, s);
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_cub_min_out, ctx->d_cub_min_out,
                                sizeof(float), cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaStreamSynchronize(s));

    return (double)(*(float *)ctx->h_cub_min_out);
}

#include "gfs_moment_gpu.inc"
#include "gfs_dg_gpu.inc"
#include "gfs_enriched_gpu.inc"
#include "gfs_enriched_coupling_gpu.inc"

#endif /* USE_CUDA */
