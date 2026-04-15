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
#include <math.h>
#include <float.h>
#include <cuda_runtime.h>

/* CUB for device-wide reduction */
#include <cub/cub.cuh>

extern "C" {
#include "exam_gpu.h"
}

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
    double dtold,
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

    double dwpdt = (sqrt(wp2) - sqrt(wpold2)) / dtold;
    double dwqdt = (sqrt(wq2) - sqrt(wqold2)) / dtold;

    /* Clamp dw/dt by sound speed */
    double vpw, vqw;
    if (dwpdt > 0) vpw = fmin((double)csi, dwpdt);
    else           vpw = fmax(-(double)csi, dwpdt);
    if (dwqdt > 0) vqw = fmin((double)csj, dwqdt);
    else           vqw = fmax(-(double)csj, dwqdt);

    double fact2 = (sqrt(wp2) * vpw - sqrt(wq2) * vqw) * dpq_inv;
    fact2 -= (wp2 - wq2) / dpq2 * er_dot_upq;

    *uradx = fact1 * upqx + fact2 * erx;
    *urady = fact1 * upqy + fact2 * ery;
}

/* ================================================================
 *  Device inline: HLLC solver
 * ================================================================ */
__device__ __forceinline__
void dev_hllc_face_2d(
    double rhoL, double pL, double vnL, double cL,
    double rhoR, double pR, double vnR, double cR,
    double Gamma,
    double *pstar, double *vnstar)
{
    double ZL = rhoL * cL, ZR = rhoR * cR;
    double GP1 = Gamma + 1.0;

    double p_pvrs = (ZR * pL + ZL * pR + ZL * ZR * (vnL - vnR)) / (ZL + ZR);
    if (p_pvrs < 0) p_pvrs = 0;

    double qL = 1.0, qR = 1.0;
    if (p_pvrs > pL)
        qL = sqrt(1.0 + GP1 / (2.0 * Gamma) * (p_pvrs / pL - 1.0));
    if (p_pvrs > pR)
        qR = sqrt(1.0 + GP1 / (2.0 * Gamma) * (p_pvrs / pR - 1.0));

    double WL = rhoL * cL * qL;
    double WR = rhoR * cR * qR;
    double Ws = WL + WR;

    *pstar  = (WR * pL + WL * pR + WL * WR * (vnL - vnR)) / Ws;
    *vnstar = (WL * vnL + WR * vnR + pL - pR) / Ws;
}

__device__ __forceinline__
void dev_hllc_face_2d_rest_frame(
    double rhoL, double pL, double vnL_lab, double cL,
    double rhoR, double pR, double vnR_lab, double cR,
    double wn, double Gamma,
    double *pstar, double *vnstar_lab)
{
    double vnL = vnL_lab - wn;
    double vnR = vnR_lab - wn;
    double pst, vnst;
    dev_hllc_face_2d(rhoL, pL, vnL, cL, rhoR, pR, vnR, cR, Gamma, &pst, &vnst);
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
    /* Output arrays */
    double *__restrict__ ax_out, double *__restrict__ ay_out,
    float *__restrict__ die_out, float *__restrict__ dt_out,
    double *__restrict__ vsig_max_out,
    /* Physics parameters */
    int n_particles,
    int av_mode, int use_muscl,
    double OoA, double Courant, double Gamma,
    double alphavis, double betavis, double etavis, double epsvis,
    double nu_phys, double prandtl,
    double cd_amax, double blend_theta, double dtold)
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
    float  ibp_w2     = pw2[i];
    float  ibp_w2old  = pw2old[i];

    double fx = 0, fy = 0, die = 0;
    float  my_dt = 1.0e10f;
    double my_vsig_max = 0;

    int f_begin = face_offset[i];
    int f_end   = face_offset[i + 1];

    for (int f = f_begin; f < f_end; f++) {
        int j = neighbor_idx[f];
        if (j < 0) continue;

        int jbp_is_ghost = is_ghost[f];

        /* Load neighbor particle */
        double jbp_x  = px[j],  jbp_y  = py[j];
        double jbp_vx = pvx[j], jbp_vy = pvy[j];
        double jbp_den      = (double)pden[j];
        double jbp_pressure = (double)ppressure[j];
        double jbp_csound   = (double)pcsound[j];
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
            dtold, &uradx, &urady);

        double pi_total;
        double tau_dot_dS_x = 0, tau_dot_dS_y = 0;

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
                double scaleFactor = (wcomp == 0 ? etavis : wcomp);
                double drampScale = dramp / scaleFactor;
                double mu = rvel / (drampScale + epsvis / drampScale);
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
                double scaleFactor = (wcomp == 0 ? etavis : wcomp);
                double drampScale = dramp / scaleFactor;
                double mu = rvel / (drampScale + epsvis / drampScale);
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
                                        wn, Gamma, &pst, &vnst_lab);
            pi_total = pst;

            /* Optional Monaghan AV for contact noise */
            if (alphavis > 0) {
                double uijx = jbp_vx - ibp_vx;
                double uijy = jbp_vy - ibp_vy;
                double rvel = erx * uijx + ery * uijy;
                if (rvel < 0) {
                    double wcomp = sqrt((double)ibp_w2) + sqrt((double)jbp_w2);
                    double scaleFactor = (wcomp == 0 ? etavis : wcomp);
                    double drampScale = dramp / scaleFactor;
                    double mu = rvel / (drampScale + epsvis / drampScale);
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
                             Gamma, &pst, &vnst);
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
                    double scaleFactor = (wcomp == 0 ? etavis : wcomp);
                    double drampScale = dramp / scaleFactor;
                    double mu = rvel / (drampScale + epsvis / drampScale);
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
                /* M(n,m) + NS stress + optional Monaghan AV */
                pi_total = p_mnm;
                if (alphavis > 0) {
                    double uijx = jbp_vx - ibp_vx;
                    double uijy = jbp_vy - ibp_vy;
                    double rvel = erx * uijx + ery * uijy;
                    if (rvel < 0) {
                        double wcomp = sqrt((double)ibp_w2) + sqrt((double)jbp_w2);
                        double scaleFactor = (wcomp == 0 ? etavis : wcomp);
                        double drampScale = dramp / scaleFactor;
                        double mu = rvel / (drampScale + epsvis / drampScale);
                        double meanden = 0.5 * (ibp_den + jbp_den);
                        double meanCsound = 0.5 * (ibp_csound + jbp_csound);
                        pi_total += (-alphavis * meanCsound * mu
                                     + betavis * mu * mu) * meanden;
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
                                                wn, Gamma, &pst, &vnst_lab);
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
            }
        }

        /* ======== Accumulate forces and energy rates ======== */
        /* uradx, urady already computed before av_mode dispatch */

        /* Pressure work: -p * (v_face - v_i) · dS */
        die += -pi_total * (uradx * dSx + urady * dSy);
        /* Viscous heating: (τ·dS) · (v_face - v_i) */
        die += tau_dot_dS_x * uradx + tau_dot_dS_y * urady;

        /* Heat conduction */
        if (nu_phys > 0 && prandtl > 0 && !jbp_is_ghost) {
            double chi = nu_phys / prandtl;
            double Ti = ibp_pressure / ibp_den;
            double Tj = jbp_pressure / jbp_den;
            double rho_face = 0.5 * (ibp_den + jbp_den);
            die += chi * rho_face * (Tj - Ti) / dramp * facearea;
        }

        /* Force: -p·dS + τ·dS */
        fx += -pi_total * dSx + tau_dot_dS_x;
        fy += -pi_total * dSy + tau_dot_dS_y;

        /* Signal velocity for CFL */
        double dvx = jbp_vx - ibp_vx;
        double dvy = jbp_vy - ibp_vy;
        double VdotR = dvx * erx + dvy * ery;
        double vsig = jbp_csound + ibp_csound - fmin(0.0, VdotR);

        /* CFL with floor */
        double heff = 0.25 * sqrt(ibp_volume);
        double dramp_cfl = fmax(dramp, heff);
        float dt_face = (float)(2.0 * Courant * dramp_cfl / vsig);
        if (dt_face < my_dt) my_dt = dt_face;

        /* Track vsig_max for CD10 */
        if (av_mode >= 1) {
            double vsig_cd = jbp_is_ghost ?
                (ibp_csound + jbp_csound) : vsig;
            if (vsig_cd > my_vsig_max) my_vsig_max = vsig_cd;
        }
    } /* end face loop */

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

    /* Write outputs */
    ax_out[i]       = fx / ibp_mass;
    ay_out[i]       = fy / ibp_mass;
    die_out[i]      = (float)die;
    dt_out[i]       = my_dt;
    vsig_max_out[i] = my_vsig_max;
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
    ALLOC_PAIR_F(die_out);  ALLOC_PAIR_F(dt_out);
    ALLOC_PAIR_D(vsig_max_out);

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
    FREE_PAIR(die_out);  FREE_PAIR(dt_out);
    FREE_PAIR(vsig_max_out);

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
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dt_out, ctx->d_parts.dt_out, sf, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.vsig_max_out, ctx->d_parts.vsig_max_out, sd, cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaStreamSynchronize(s));
}

/* ================================================================
 *  Kernel launch + CFL reduction
 * ================================================================ */

extern "C"
double gpu_launch_force_kernel(GPUContext *ctx, int n_particles,
                               const GPUPhysicsParams *params)
{
    int blockSize = 256;
    int gridSize = (n_particles + blockSize - 1) / blockSize;

    cudaStream_t s = (cudaStream_t)ctx->stream;

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
        /* Output */
        ctx->d_parts.ax_out, ctx->d_parts.ay_out,
        ctx->d_parts.die_out, ctx->d_parts.dt_out,
        ctx->d_parts.vsig_max_out,
        /* Physics */
        n_particles,
        params->av_mode, params->use_muscl,
        params->OoA, params->Courant, params->Gamma,
        params->alphavis, params->betavis,
        params->etavis, params->epsvis,
        params->nu_phys, params->prandtl,
        params->cd_amax, params->blend_theta, params->dtold);

    CUDA_CHECK(cudaGetLastError());

    /* CFL reduction using pre-allocated CUB scratch */
    cub::DeviceReduce::Min(ctx->d_cub_temp, ctx->cub_temp_bytes,
                           ctx->d_parts.dt_out,
                           (float *)ctx->d_cub_min_out,
                           n_particles, s);
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_cub_min_out, ctx->d_cub_min_out,
                                sizeof(float), cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaStreamSynchronize(s));

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

/* Max neighbors gathered from 3×3 stencil */
#define GPU_MAX_NEIGH 256

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
    /* Output: polygon for face scatter */
    PolygonOut *__restrict__ polygons,
    int *__restrict__ face_count_out,
    /* Grid parameters */
    int nbp, int n_total, int mx, int my,
    double cellsize, double xmin, double ymin,
    double boxsize, double OoA, int av_mode)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nbp) return;

    double xi = px[i], yi = py[i];
    double w2i = (double)pw2[i];
    long long indx_i = pindx[i];
    double massi = pmass[i];

    /* Determine cell (ix, iy) from position */
    int ix = (int)((xi - xmin) / cellsize);
    int iy = (int)((yi - ymin) / cellsize);
    if (ix < 0) ix = 0;  if (ix >= mx) ix = mx - 1;
    if (iy < 0) iy = 0;  if (iy >= my) iy = my - 1;

    /* ---- Gather neighbors from 3×3 stencil ---- */
    double n_dx[GPU_MAX_NEIGH], n_dy[GPU_MAX_NEIGH];
    double n_w2[GPU_MAX_NEIGH], n_drad2[GPU_MAX_NEIGH];
    int    n_soa_idx[GPU_MAX_NEIGH];
    int    n_count = 0;

    for (int dy_c = -1; dy_c <= 1; dy_c++) {
        int cy = iy + dy_c;
        if (cy < 0 || cy >= my) continue;
        for (int dx_c = -1; dx_c <= 1; dx_c++) {
            int cx = ix + dx_c;
            if (cx < 0 || cx >= mx) continue;
            int cell_id = cx + mx * cy;
            int c_begin = cell_offset[cell_id];
            int c_end   = cell_offset[cell_id + 1];
            for (int k = c_begin; k < c_end && n_count < GPU_MAX_NEIGH; k++) {
                int j = cell_particle_index[k];
                if (pindx[j] == indx_i) continue;  /* skip self */
                double dx = px[j] - xi;
                double dy_val = py[j] - yi;
                double dist2 = dx * dx + dy_val * dy_val;
                double w2j = (double)pw2[j];
                double wdiff = (w2i - w2j) / dist2;
                double drad2 = 0.25 * dist2 * (1.0 + wdiff) * (1.0 + wdiff);
                n_dx[n_count] = dx;
                n_dy[n_count] = dy_val;
                n_w2[n_count] = w2j;
                n_drad2[n_count] = drad2;
                n_soa_idx[n_count] = j;
                n_count++;
            }
        }
    }

    /* ---- Insertion sort by drad2 ---- */
    for (int a = 1; a < n_count; a++) {
        double key_d = n_drad2[a];
        double key_dx = n_dx[a], key_dy = n_dy[a], key_w2 = n_w2[a];
        int key_idx = n_soa_idx[a];
        int b = a - 1;
        while (b >= 0 && n_drad2[b] > key_d) {
            n_drad2[b+1] = n_drad2[b];
            n_dx[b+1] = n_dx[b]; n_dy[b+1] = n_dy[b];
            n_w2[b+1] = n_w2[b]; n_soa_idx[b+1] = n_soa_idx[b];
            b--;
        }
        n_drad2[b+1] = key_d;
        n_dx[b+1] = key_dx; n_dy[b+1] = key_dy;
        n_w2[b+1] = key_w2; n_soa_idx[b+1] = key_idx;
    }

    /* ---- Initialize box polygon ---- */
    double halfbox = 0.5 * boxsize;
    /* 4 box corners: (-hb,-hb), (+hb,-hb), (+hb,+hb), (-hb,+hb) */
    double poly_x[GPU_MAX_CORNERS], poly_y[GPU_MAX_CORNERS];
    int    poly_neigh[GPU_MAX_CORNERS];
    int    nc = 4;
    /* Same init as Voro2D_InitializeCorner: corner i = (((i+1)%4)/2, i/2)*box - halfbox */
    poly_x[0] = -halfbox; poly_y[0] = -halfbox; poly_neigh[0] = -1;
    poly_x[1] =  halfbox; poly_y[1] = -halfbox; poly_neigh[1] = -2;
    poly_x[2] =  halfbox; poly_y[2] =  halfbox; poly_neigh[2] = -3;
    poly_x[3] = -halfbox; poly_y[3] =  halfbox; poly_neigh[3] = -4;

    double maxdist2 = halfbox * sqrt(3.0);

    /* ---- Sutherland-Hodgman clipping ---- */
    double tmp_x[GPU_MAX_CORNERS], tmp_y[GPU_MAX_CORNERS];
    int    tmp_neigh[GPU_MAX_CORNERS];

    for (int ni = 0; ni < n_count; ni++) {
        if (n_drad2[ni] > maxdist2) break;

        double nx_n = n_dx[ni], ny_n = n_dy[ni];
        double dist2 = nx_n * nx_n + ny_n * ny_n;
        double wfrac = 0.5 + 0.5 * (w2i - n_w2[ni]) / dist2;
        double threshold = wfrac * dist2;  /* wfrac * dot(neigh, neigh) */
        int j_soa = n_soa_idx[ni];

        int out_count = 0;
        for (int e = 0; e < nc; e++) {
            int e_next = (e + 1) % nc;
            double ax = poly_x[e],    ay = poly_y[e];
            double bx = poly_x[e_next], by = poly_y[e_next];
            double da = ax * nx_n + ay * ny_n;  /* dot(A, neigh) */
            double db = bx * nx_n + by * ny_n;  /* dot(B, neigh) */
            int a_in = (da <= threshold);
            int b_in = (db <= threshold);

            if (a_in && b_in) {
                /* Both inside: output A with original edge_neigh */
                if (out_count < GPU_MAX_CORNERS) {
                    tmp_x[out_count] = ax;
                    tmp_y[out_count] = ay;
                    tmp_neigh[out_count] = poly_neigh[e];
                    out_count++;
                }
            } else if (a_in && !b_in) {
                /* A in, B out: output A, then intersection with neighbor ni */
                if (out_count < GPU_MAX_CORNERS) {
                    tmp_x[out_count] = ax;
                    tmp_y[out_count] = ay;
                    tmp_neigh[out_count] = poly_neigh[e];
                    out_count++;
                }
                /* Intersection: t = (threshold - da) / (db - da) */
                double denom = db - da;
                double t = (fabs(denom) > 1e-30) ? (threshold - da) / denom : wfrac;
                t = fmax(1e-9, fmin(1.0 - 1e-9, t));
                if (out_count < GPU_MAX_CORNERS) {
                    tmp_x[out_count] = ax + t * (bx - ax);
                    tmp_y[out_count] = ay + t * (by - ay);
                    tmp_neigh[out_count] = j_soa;  /* new face starts */
                    out_count++;
                }
            } else if (!a_in && b_in) {
                /* A out, B in: output intersection with original edge_neigh */
                double denom = db - da;
                double t = (fabs(denom) > 1e-30) ? (threshold - da) / denom : wfrac;
                t = fmax(1e-9, fmin(1.0 - 1e-9, t));
                if (out_count < GPU_MAX_CORNERS) {
                    tmp_x[out_count] = ax + t * (bx - ax);
                    tmp_y[out_count] = ay + t * (by - ay);
                    tmp_neigh[out_count] = poly_neigh[e];
                    out_count++;
                }
            }
            /* both out: skip */
        }

        nc = out_count;
        for (int e = 0; e < nc; e++) {
            poly_x[e] = tmp_x[e];
            poly_y[e] = tmp_y[e];
            poly_neigh[e] = tmp_neigh[e];
        }

        if (nc < 3) break;  /* degenerate polygon */

        /* Update maxdist2 from farthest surviving corner */
        double new_max = 0;
        for (int e = 0; e < nc; e++) {
            double cd = poly_x[e] * poly_x[e] + poly_y[e] * poly_y[e];
            if (cd > new_max) new_max = cd;
        }
        maxdist2 = new_max;
    }

    /* ---- Merge near-coincident vertices (degenerate faces at grid corners) ----
     * On a regular/near-regular grid, 4+ Voronoi cells meet at (near) a
     * single point. Sutherland-Hodgman clipping produces two nearly-
     * identical vertices, creating a tiny face.  These cause CFL dt
     * collapse.  Use RELATIVE threshold: edge < 1e-4 * max_edge. */
    {
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
    if (volume < 1e-30) volume = 1e-30;
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
            double val = jdist2 + (double)pw2[j_neigh];
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

    /* ---- Green-Gauss gradients ---- */
    double sum_gUxx = 0, sum_gUxy = 0, sum_gUyx = 0, sum_gUyy = 0;
    double sum_dPdx = 0, sum_dPdy = 0;
    double sum_dRhodx = 0, sum_dRhody = 0;

    for (int e = 0; e < nc; e++) {
        int en = (e + 1) % nc;
        int j_neigh = poly_neigh[e];
        if (j_neigh < 0) continue;  /* box boundary edge, skip */

        /* dS = (c2y - c1y, -(c2x - c1x)) */
        double dSx = poly_y[en] - poly_y[e];
        double dSy = -(poly_x[en] - poly_x[e]);

        /* Laguerre-aware face interpolation weight */
        double dx_ij = px[j_neigh] - xi;
        double dy_ij = py[j_neigh] - yi;
        double d2_ij = dx_ij * dx_ij + dy_ij * dy_ij;
        double wfrac_g = (d2_ij > 0) ?
            0.5 + 0.5 * (w2i - (double)pw2[j_neigh]) / d2_ij : 0.5;

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
            /* Check M(n,m) stencil neighbors: kp = edge_neigh[(e+1)%nc], km = edge_neigh[(e-1+nc)%nc] */
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

    double invVol = 1.0 / volume;
    double gUxx = sum_gUxx * invVol, gUxy = sum_gUxy * invVol;
    double gUyx = sum_gUyx * invVol, gUyy = sum_gUyy * invVol;
    double divv_val = gUxx + gUyy;
    double dPdx_val = sum_dPdx * invVol, dPdy_val = sum_dPdy * invVol;
    double dRhodx_val = sum_dRhodx * invVol, dRhody_val = sum_dRhody * invVol;

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
        gUxx *= alpha_vx;         gUxy *= alpha_vx;
        gUyx *= alpha_vy;         gUyy *= alpha_vy;
        dPdx_val *= alpha_P;      dPdy_val *= alpha_P;
        divv_val = gUxx + gUyy;  /* recompute after limiting */
    }

    /* ---- Write outputs ---- */
    out_volume[i]   = (float)volume;
    out_den[i]      = (float)density;
    out_w2ceil[i]   = (float)w2ceil_val;
    out_avgNeighP[i]= (float)avgNP;

    out_gUxx[i] = gUxx;  out_gUxy[i] = gUxy;
    out_gUyx[i] = gUyx;  out_gUyy[i] = gUyy;
    out_dPdx[i] = dPdx_val;  out_dPdy[i] = dPdy_val;
    out_dRhodx[i] = dRhodx_val;  out_dRhody[i] = dRhody_val;

    /* Write PolygonOut for face scatter */
    PolygonOut *po = &polygons[i];
    po->n_corners = nc;
    int n_real_faces = 0;
    for (int e = 0; e < nc; e++) {
        po->cx[e] = poly_x[e];
        po->cy[e] = poly_y[e];
        po->edge_neigh[e] = poly_neigh[e];
        if (poly_neigh[e] >= 0) n_real_faces++;
    }
    po->n_real_faces = n_real_faces;
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
    const double *__restrict__ pgUxx, const double *__restrict__ pgUxy,
    const double *__restrict__ pgUyx, const double *__restrict__ pgUyy,
    const double *__restrict__ palpha_cd,
    float *__restrict__ out_pressure, float *__restrict__ out_csound,
    double *__restrict__ out_tauxx, double *__restrict__ out_tauxy,
    double *__restrict__ out_tauyy,
    double *__restrict__ out_divv,
    int nbp, int av_mode, double Gamma, double nu_phys, double cd_amax)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nbp) return;

    double ie_val = (double)pie[i];
    double vol    = (double)pvolume[i];
    double den    = (double)pden[i];

    double P = ie_val / vol * (Gamma - 1.0);
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
        double h = sqrt(vol);
        double nu_cd = palpha_cd[i] * h * cs;
        double nu = (nu_phys > 0) ? fmax(nu_phys, nu_cd) : nu_cd;
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
    int nbp)
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
    double boxsize, double OoA, int av_mode, double Gamma,
    double nu_phys, double prandtl, double cd_amax)
{
    cudaStream_t s = (cudaStream_t)ctx->stream;
    int blockSize = 256;
    int gridSize = (nbp + blockSize - 1) / blockSize;
    int n_total = ctx->d_parts.n_total;

    /* Kernel A: Voronoi tessellation + gradients */
    voronoi_tessellate_kernel<<<gridSize, blockSize, 0, s>>>(
        ctx->d_parts.x, ctx->d_parts.y,
        ctx->d_parts.vx, ctx->d_parts.vy,
        ctx->d_parts.mass, ctx->d_parts.indx,
        ctx->d_parts.ie, ctx->d_parts.w2,
        ctx->d_parts.pressure, ctx->d_parts.den,
        ctx->d_cells.cell_offset, ctx->d_cells.particle_index,
        /* outputs */
        ctx->d_parts.den, ctx->d_parts.pressure,
        ctx->d_parts.volume, ctx->d_parts.w2ceil, ctx->d_parts.avgNeighP,
        ctx->d_parts.gUxx, ctx->d_parts.gUxy,
        ctx->d_parts.gUyx, ctx->d_parts.gUyy,
        ctx->d_parts.dPdx, ctx->d_parts.dPdy,
        ctx->d_parts.dRhodx, ctx->d_parts.dRhody,
        ctx->d_polygons, ctx->d_face_count,
        nbp, n_total, mx, my, cellsize, xmin, ymin, boxsize, OoA, av_mode);
    CUDA_CHECK(cudaGetLastError());

    /* Kernel B: Pressure, sound speed, NS stress */
    pressure_stress_kernel<<<gridSize, blockSize, 0, s>>>(
        ctx->d_parts.ie, ctx->d_parts.volume, ctx->d_parts.den,
        ctx->d_parts.gUxx, ctx->d_parts.gUxy,
        ctx->d_parts.gUyx, ctx->d_parts.gUyy,
        ctx->d_parts.alpha_cd,
        ctx->d_parts.pressure, ctx->d_parts.csound,
        ctx->d_parts.tauxx, ctx->d_parts.tauxy, ctx->d_parts.tauyy,
        ctx->d_parts.divv,
        nbp, av_mode, Gamma, nu_phys, cd_amax);
    CUDA_CHECK(cudaGetLastError());

    /* CUB prefix sum: face_count → face_offset */
    cub::DeviceScan::ExclusiveSum(ctx->d_cub_scan_temp, ctx->cub_scan_temp_bytes,
                                   ctx->d_face_count, ctx->d_face_offset,
                                   nbp, s);

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
        int total_faces = h_last_offset + h_last_count;

        /* Set face_offset[nbp] on device */
        CUDA_CHECK(cudaMemcpyAsync(ctx->d_face_offset + nbp, &total_faces,
                   sizeof(int), cudaMemcpyHostToDevice, s));

        /* Update FaceCSR metadata */
        ctx->d_faces.n_particles = nbp;
        ctx->d_faces.n_faces_total = total_faces;

        /* Check capacity (should not happen with proper allocation) */
        if (total_faces > ctx->max_faces) {
            fprintf(stderr, "[GPU TESS] WARNING: total_faces=%d > max_faces=%d, "
                    "faces may be truncated!\n", total_faces, ctx->max_faces);
        }
    }

    /* Kernel C: Scatter faces → FaceCSR */
    scatter_faces_kernel<<<gridSize, blockSize, 0, s>>>(
        ctx->d_polygons, ctx->d_face_offset, ctx->d_parts.indx,
        ctx->d_faces.c1x, ctx->d_faces.c1y,
        ctx->d_faces.c2x, ctx->d_faces.c2y,
        ctx->d_faces.neighbor_idx,
        ctx->d_faces.kp_idx, ctx->d_faces.km_idx,
        ctx->d_faces.is_ghost,
        nbp);
    CUDA_CHECK(cudaGetLastError());

    /* Copy face_offset to FaceCSR (device already has it, just set pointer) */
    /* The force kernel uses d_faces.face_offset. Copy d_face_offset → d_faces.face_offset */
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_faces.face_offset, ctx->d_face_offset,
               (nbp + 1) * sizeof(int), cudaMemcpyDeviceToDevice, s));

    CUDA_CHECK(cudaStreamSynchronize(s));
    ctx->tess_faces_on_device = 1;
}

/* ================================================================
 *  gpu_download_tess_results: device SoA → pinned host SoA
 *
 *  Downloads the tessellation outputs (den, P, cs, volume, w2ceil,
 *  avgNeighP, gradients, stress) so CPU can write them back to AoS.
 * ================================================================ */
extern "C"
void gpu_download_tess_results(GPUContext *ctx, int nbp, int has_stress)
{
    cudaStream_t s = (cudaStream_t)ctx->stream;
    size_t sf = nbp * sizeof(float);
    size_t sd = nbp * sizeof(double);

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
    if (kappa >= 0.0f) {
        w2_out[i] = fminf(w2_out[i], (float)min_dist2);
    }
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
    double *__restrict__ divv_out,
    double *__restrict__ tauxx_out, double *__restrict__ tauxy_out,
    double *__restrict__ tauyy_out,
    float *__restrict__ w2ceil_out,
    /* Parameters */
    int n_particles, int n_total,
    double eta, int h_iter_max, double h_tol,
    double cellsize, int mx, int my,
    double xmin, double ymin,
    double Gamma, double nu_phys)
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
    double Exx = 0, Exy = 0, Eyy = 0;
    double Sxx_vx = 0, Sxy_vx = 0;
    double Sxx_vy = 0, Sxy_vy = 0;
    double Sxx_P = 0, Sxy_P = 0;

    /* Newton iteration for adaptive h */
    for (int h_it = 0; h_it < h_iter_max; h_it++) {
        n_i = 0;
        Exx = Exy = Eyy = 0;
        Sxx_vx = Sxy_vx = 0;
        Sxx_vy = Sxy_vy = 0;
        Sxx_P = Sxy_P = 0;

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
                    }
                }
            }
        }

        if (n_i <= 0) break;  /* isolated */

        double V_new = 1.0 / n_i;
        double h_new = eta * sqrt(V_new);
        double dh = fabs(h_new - h_i);
        if (dh < h_tol * h_i) {
            h_i = h_new;
            break;
        }
        h_i = h_new;
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
    double divv_val = gux + gvy;

    /* EOS */
    double ie_val = (double)pie_in[i];
    double P_val = ie_val / Vi * (Gamma - 1.0);
    if (P_val <= 0) {
        P_val = 1e-6;
        ie_val = P_val * Vi / (Gamma - 1.0);
    }
    double rho_val = imass * n_i;
    double cs_val = sqrt(Gamma * P_val / rho_val);

    /* NS stress (3D-slab trace convention, matches CPU) */
    double h_eff = sqrt(Vi);
    double alpha_cd_val = palpha_cd[i];
    double nu_cd = alpha_cd_val * h_eff * cs_val;
    double nu = (nu_phys > 0) ? fmax(nu_phys, nu_cd) : nu_cd;
    double tau_xx = -nu * rho_val * (2.0 * gux - (2.0 / 3.0) * divv_val);
    double tau_xy = -nu * rho_val * (guy + gvx);
    double tau_yy = -nu * rho_val * (2.0 * gvy - (2.0 / 3.0) * divv_val);

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
    double nu_phys, double eta)
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

                /* Effective face area vector */
                double W_ij = dev_wendland_c2_2d(r, h_i);
                double W_ji = dev_wendland_c2_2d(r, h_j);
                double Vi2 = iV;       /* Hopkins (2015): V_i, not V_i^2 */
                double Vj  = (double)pvolume[j];
                double Vj2 = Vj;       /* Hopkins (2015): V_j, not V_j^2 */

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

                /* MUSCL reconstruction */
                double hdx = 0.5 * dx, hdy = 0.5 * dy;
                double pL  = iP  + i_dPdx * hdx + i_dPdy * hdy;
                double pR  = jP  - pdPdx[j] * hdx - pdPdy[j] * hdy;
                double vxL = ivx + i_gUxx * hdx + i_gUxy * hdy;
                double vxR = jvx - pgUxx[j] * hdx - pgUxy[j] * hdy;
                double vyL = ivy + i_gUyx * hdx + i_gUyy * hdy;
                double vyR = jvy - pgUyx[j] * hdx - pgUyy[j] * hdy;

                /* Pair-wise clamp */
                #define DEV_CLAMP(Qrec,Qi,Qj) do { \
                    double qmn = fmin(Qi,Qj), qmx = fmax(Qi,Qj); \
                    if (Qrec < qmn) Qrec = qmn; \
                    if (Qrec > qmx) Qrec = qmx; \
                } while(0)
                DEV_CLAMP(pL,  iP,  jP);
                DEV_CLAMP(pR,  iP,  jP);
                DEV_CLAMP(vxL, ivx, jvx);
                DEV_CLAMP(vxR, ivx, jvx);
                DEV_CLAMP(vyL, ivy, jvy);
                DEV_CLAMP(vyR, ivy, jvy);
                #undef DEV_CLAMP

                /* Wall: skip reconstruction */
                if (jwall) {
                    pL  = iP;  vxL = ivx; vyL = ivy;
                    pR  = jP;  vxR = jvx; vyR = jvy;
                }

                double vnL_lab = vxL * nx_hat + vyL * ny_hat;
                double vnR_lab = vxR * nx_hat + vyR * ny_hat;
                double rhoL = irho, rhoR = jrho;
                double cL = ics, cR = jcs;

                /* HLLC face rest frame */
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
                        wn, Gamma, &pstar, &vnstar_lab);
                }

                /* Optional Monaghan AV */
                if (alphavis > 0 && !jwall) {
                    double uij_x = jvx - ivx;
                    double uij_y = jvy - ivy;
                    double er_x = dx / r, er_y = dy / r;
                    double rvel = er_x * uij_x + er_y * uij_y;
                    if (rvel < 0) {
                        double wcomp = sqrt(iw2) + sqrt((double)pw2[j]);
                        double scaleFactor = (wcomp == 0 ? etavis : wcomp);
                        double drampScale = r / scaleFactor;
                        double mu = rvel / (drampScale + epsvis / drampScale);
                        double meanden = 0.5 * (irho + jrho);
                        double meanCs  = 0.5 * (ics + jcs);
                        pstar += (-alphavis * meanCs * mu + betavis * mu * mu) * meanden;
                    }
                }

                /* NS viscous traction */
                double tau_A_x = 0, tau_A_y = 0;
                if (!jwall) {
                    double txx_f = 0.5 * (i_tauxx + ptauxx[j]);
                    double txy_f = 0.5 * (i_tauxy + ptauxy[j]);
                    double tyy_f = 0.5 * (i_tauyy + ptauyy[j]);
                    tau_A_x = -(txx_f * Ax + txy_f * Ay);
                    tau_A_y = -(txy_f * Ax + tyy_f * Ay);
                }

                /* Force */
                fx += -pstar * Ax + tau_A_x;
                fy += -pstar * Ay + tau_A_y;

                /* Full face velocity */
                double tx_hat = -ny_hat;
                double ty_hat =  nx_hat;
                double vt_avg = 0.5 * ((ivx + jvx) * tx_hat + (ivy + jvy) * ty_hat);
                double vfx = vnstar_lab * nx_hat + vt_avg * tx_hat;
                double vfy = vnstar_lab * ny_hat + vt_avg * ty_hat;

                if (jwall) {
                    /* Free-slip: tangential = v_i, normal = 0 */
                    double vti = ivx * tx_hat + ivy * ty_hat;
                    vfx = vti * tx_hat;
                    vfy = vti * ty_hat;
                }

                /* Energy */
                die += -pstar * ((vfx - ivx) * Ax + (vfy - ivy) * Ay);
                die += tau_A_x * (vfx - ivx) + tau_A_y * (vfy - ivy);

                /* Signal velocity for CFL */
                double dvx_sig = jvx - ivx;
                double dvy_sig = jvy - ivy;
                double er_x = dx / r, er_y = dy / r;
                double VdotR = dvx_sig * er_x + dvy_sig * er_y;
                double vsig = cL + cR - fmin(0.0, VdotR);
                double heff = sqrt(iV);
                double dramp_cfl = fmax(r, heff);
                float dt_face = (float)(2.0 * Courant * dramp_cfl / vsig);
                if (dt_face < my_dt) my_dt = dt_face;

                if (!jwall && vsig > my_vsig_max)
                    my_vsig_max = vsig;
            }
        }
    }

    /* Viscous CFL */
    {
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
        ctx->d_parts.divv,
        ctx->d_parts.tauxx, ctx->d_parts.tauxy, ctx->d_parts.tauyy,
        ctx->d_parts.w2ceil,
        /* Parameters */
        n_particles, n_total,
        params->lagmfm_eta, params->lagmfm_h_iter, params->lagmfm_h_tol,
        params->cellsize, params->mx, params->my,
        params->xmin, params->ymin,
        params->Gamma, params->nu_phys);

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
        params->nu_phys, params->lagmfm_eta);

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

#endif /* USE_CUDA */
