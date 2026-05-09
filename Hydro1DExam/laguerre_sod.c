/*
 * laguerre_sod.c  -  1D Hydro Test Suite
 * Laguerre vs Voronoi vs SPH  |  RK4 integration
 * Problems: Sod, Blast, Shu-Osher, Noh, Lax, DblFan, Collision, Contact
 *
 * gcc -O2 -o laguerre_sod laguerre_sod.c -lm
 * ./laguerre_sod [N]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

/* ===== Fixed constants ===== */
#define NMAX 5000
#define CFL     0.3
#define SPH_A   1.0
#define SPH_B   2.0
#define SPH_ETA 2.0

/* ===== Problem / BC enums ===== */
enum { PROB_SOD=0, PROB_BLAST=1, PROB_SHUOSHER=2, PROB_NOH=3,
       PROB_LAX=4, PROB_DBLFAN=5, PROB_COLLISION=6, PROB_CONTACT=7 };
enum { BC_REFLECT=0, BC_TRANSMIT=1 };

/* ===== Runtime problem globals ===== */
static int    problem_id;
static int    bc_left, bc_right;
static int    n_regions;
static int    has_exact_sol;
static double GAM, GM1, GP1;
static double x0_dom, x1_dom, xmid1, xmid2, tend;
static double rhoL_ic, vL_ic, pL_ic;
static double rhoM_ic, vM_ic, pM_ic;   /* Blast middle region */
static double rhoR_ic, vR_ic, pR_ic;
static double ps_exact, vs_exact;       /* exact Riemann star state */
static double sph_hmin, sph_hmax;       /* SPH smoothing length limits */

/* ===== Mode enum ===== */
enum { LAGUERRE=0, VORONOI=1, SPH_MODE=2 };

/* ===== AV parameters (runtime-configurable for sweep) ===== */
static double av_alpha = 3.0;    /* Monaghan alpha */
static double av_beta  = 4.0;    /* Monaghan beta  */
static double av_eps   = 0.01;   /* Monaghan softening */
static double face_dp  = 0.0;    /* dp/Z face velocity correction (B: vf) */
static double face_dv  = 0.0;    /* dv*Z face pressure correction (B: pf) */
static double alpha_u  = 0.0;    /* artificial thermal conductivity coeff */
static double w_coeff  = 0.4;    /* Laguerre w coefficient */

/* Laguerre weight mode */
enum { WM_PRESSURE=0, WM_ENTHALPY=1, WM_SPECVOL=2, WM_CSOUND=3 };
static int    weight_mode = WM_PRESSURE;
static double h0_ref = 1.0;             /* enthalpy reference */
static double v0_ref = 1.0;             /* specific volume reference */
static double c0_ref = 1.0;             /* sound speed reference */
static double face_shift_max = 1.0;     /* face shift limiter (fraction of dx) */
static double w2_rate_max   = 0.0;     /* max fractional change of w^2 per step (0=off) */
static double p_exp        = 0.0;     /* pressure exponent for WM_PRESSURE (0=use (γ-1)/γ) */
static double w_floor_frac = 0.0;     /* w floor as fraction of meand (0=off) */
static int    use_adiab_w  = 0;       /* adiabatic w correction: w²*=(Vvor/Vlag)^(γ-1) */
static double shock_w_damp = 0.0;    /* shock-switch: w²/=(1+κ*S) where S=compress*h/c (0=off) */
static double shock_p_ratio = 0.0;  /* pressure-ratio threshold for shock_w_damp (0=always active) */
static double w_relax_tau  = 0.0;    /* w relaxation timescale in units of h/c_ref (0=off) */
static double w_grad_max   = 0.0;    /* max |w_i-w_j|/dx gradient limiter (0=off) */
static double shock_th = 0.0;    /* w2 shock detector threshold (D: 0=off) */
static double av_al_mn = 0.0;    /* AV alpha minimum (smooth region) */
static double av_bt_mn = 0.0;    /* AV beta minimum (smooth region) */
static int    use_muscl= 0;      /* MUSCL reconstruction (C: 0=off) */
static int    use_riemann=0;     /* 1=acoustic, 2=exact, 3=HLLC */
static int    av_mode = 0;       /* 0=Monaghan, 1=Price signal, 2=CD10 */
static double nu_phys  = 0.0;   /* fixed kinematic viscosity for NS stress (0=off) */
static double nu_coeff = 0.0;   /* adaptive: nu += nu_coeff*h*cs (0=off) */

/* Cullen-Dehnen (2010) switch parameters */
static double cd_ell  = 0.1;     /* decay length parameter */
static double cd_amax = 2.0;     /* alpha max */
static double cd_amin = 0.0;     /* alpha min */
static double alpha_cd[NMAX];    /* per-cell time-dependent alpha */
static double divv_old[NMAX];    /* div v from previous timestep */
static double vsig_max[NMAX];    /* max signal velocity per cell */

/* global y-range tracking */
static double gmin[4]={ 1e30, 1e30, 1e30, 1e30};
static double gmax[4]={-1e30,-1e30,-1e30,-1e30};
static int cmp_dbl(const void *a, const void *b){
    double da=*(const double*)a, db=*(const double*)b;
    return (da>db)-(da<db);
}

/* ===== Reference solution arrays ===== */
static int    nref;
static double ref_x[NMAX], ref_rho[NMAX], ref_v[NMAX], ref_p[NMAX];

/* ===== Particle ===== */
typedef struct {
    double x, v, e, m;
    double rho, p, c, w, vol;
    double av, ae;
} Pt;

static int    N, mode;
static double meand;
static double p0_ref=0.2;
static Pt     P[NMAX];
static double face[NMAX+1], pf[NMAX+1], vf[NMAX+1], qf[NMAX+1];
static double hh[NMAX];
static double s0x[NMAX], s0v[NMAX], s0e[NMAX];
static double kx[4][NMAX], kv[4][NMAX], ke[4][NMAX];
static double ow2[NMAX];    /* w^2 at start of previous timestep */
static double prev_dt=0;    /* previous timestep size */
static int    use_AA=0;     /* 1=finite-diff AA, 2=analytical AA */
static double aa_max_frac = 1.0;  /* AA limiter as fraction of Cs */
static double divv[NMAX];   /* velocity divergence from previous RK4 stage */

/* ===== EOS ===== */
static inline double eos_p(double r,double e){return GM1*r*e;}
static inline double eos_e(double r,double p){return p/(GM1*r);}
static inline double eos_c(double r,double p){
    return(r>0&&p>0)?sqrt(GAM*p/r):1e-10;
}

static inline double w2Measure(double pressure){
    double t=w_coeff*meand*pow(pressure/p0_ref, GM1/GAM);
    return t*t;
}
static inline double w2Measure_gen(int i){
    double val;
    switch(weight_mode){
    case WM_PRESSURE: {
        double alpha = (p_exp > 0) ? p_exp : GM1/GAM;
        val = pow(fmax(P[i].p,1e-14)/p0_ref, alpha);
        break; }
    case WM_ENTHALPY: {
        double h = P[i].e + P[i].p/fmax(P[i].rho,1e-14);
        val = h / h0_ref;
        break; }
    case WM_SPECVOL:
        val = (1.0/fmax(P[i].rho,1e-14)) / v0_ref;
        break;
    case WM_CSOUND:
        val = P[i].c / c0_ref;
        break;
    default:
        val = pow(fmax(P[i].p,1e-14)/p0_ref, GM1/GAM);
    }
    double t = w_coeff * meand * val;
    if(w_floor_frac > 0){
        double wf = w_floor_frac * meand;
        if(t < wf) t = wf;
    }
    return t*t;
}

/* ============================================================
 *  SPH cubic spline kernel (1D)
 * ============================================================ */
static double Wk(double r, double h){
    double q=fabs(r)/h;
    if(q>=2) return 0;
    double s=2.0/(3.0*h);
    if(q<=1) return s*(1.0-1.5*q*q+0.75*q*q*q);
    double t=2.0-q; return s*0.25*t*t*t;
}
static double dWk(double xij, double h){
    double r=fabs(xij), q=r/h;
    if(q>=2.0||r<1e-30) return 0;
    double sg=(xij>0)?1.0:-1.0, s=2.0/(3.0*h*h);
    if(q<=1) return s*(-3.0*q+2.25*q*q)*sg;
    double t=2.0-q; return s*(-0.75*t*t)*sg;
}

/* ============================================================
 *  Cell geometry
 * ============================================================ */
static void cells_vor(void){
    for(int i=0;i<N-1;i++)
        face[i+1]=0.5*(P[i].x+P[i+1].x);
    /* Boundary faces: fixed for reflecting, extrapolated for transmissive */
    if(bc_left == BC_TRANSMIT && N >= 2)
        face[0] = P[0].x - 0.5*(P[1].x - P[0].x);
    else
        face[0] = x0_dom;
    if(bc_right == BC_TRANSMIT && N >= 2)
        face[N] = P[N-1].x + 0.5*(P[N-1].x - P[N-2].x);
    else
        face[N] = x1_dom;
    for(int i=0;i<N;i++){
        P[i].vol=face[i+1]-face[i];
        if(P[i].vol<1e-15) P[i].vol=1e-15;
        P[i].rho=P[i].m/P[i].vol;
        P[i].p=eos_p(P[i].rho,P[i].e);
        P[i].c=eos_c(P[i].rho,P[i].p);
        P[i].w=0;
    }
}

static void cells_lag(void){
    /* Boundary faces: fixed for reflecting, extrapolated for transmissive */
    if(bc_left == BC_TRANSMIT && N >= 2)
        face[0] = P[0].x - 0.5*(P[1].x - P[0].x);
    else
        face[0] = x0_dom;
    if(bc_right == BC_TRANSMIT && N >= 2)
        face[N] = P[N-1].x + 0.5*(P[N-1].x - P[N-2].x);
    else
        face[N] = x1_dom;

    /* Voronoi volumes for adiabatic w correction */
    double vvol[NMAX];
    if(use_adiab_w && w_coeff > 0){
        double vf0 = x0_dom;
        for(int i=0;i<N;i++){
            double vf1 = (i<N-1) ? 0.5*(P[i].x+P[i+1].x) : x1_dom;
            vvol[i] = vf1 - vf0;
            if(vvol[i] < 1e-15) vvol[i] = 1e-15;
            vf0 = vf1;
        }
    }

    for(int i=0;i<N;i++){
        double dleft  = (i>0)   ? (P[i].x-P[i-1].x) : 2*(P[i].x-x0_dom);
        double dright = (i<N-1) ? (P[i+1].x-P[i].x) : 2*(x1_dom-P[i].x);
        double ceil2  = fmin(dleft, dright); ceil2 *= ceil2;
        double w2     = w2Measure_gen(i);

        /* adiabatic correction: w² *= (V_vor / V_lag_prev)^(γ-1) */
        if(use_adiab_w && w_coeff > 0 && P[i].vol > 1e-15){
            double ratio = vvol[i] / P[i].vol;
            if(ratio > 0.1 && ratio < 10.0)
                w2 *= pow(ratio, GM1);
        }

        /* Method 1: Shock-switch — damp w in compression regions
         * S = |compress| * h / c_s  (dimensionless shock strength)
         * w² /= (1 + κ * S²)  →  0 at strong shocks
         * Optional: only activate when local pressure ratio > shock_p_ratio */
        if(shock_w_damp > 0){
            int do_damp = 1;
            if(shock_p_ratio > 0){
                /* check max pressure ratio with neighbors */
                double pi = fmax(P[i].p, 1e-30);
                double prmax = 1.0;
                if(i > 0){
                    double pn = fmax(P[i-1].p, 1e-30);
                    double r = (pi > pn) ? pi/pn : pn/pi;
                    if(r > prmax) prmax = r;
                }
                if(i < N-1){
                    double pn = fmax(P[i+1].p, 1e-30);
                    double r = (pi > pn) ? pi/pn : pn/pi;
                    if(r > prmax) prmax = r;
                }
                do_damp = (prmax > shock_p_ratio);
            }
            if(do_damp){
                double compress = fmax(0.0, -divv[i]);
                double ci = fmax(P[i].c, 1e-10);
                double S = compress * meand / ci;
                w2 /= (1.0 + shock_w_damp * S * S);
            }
        }

        /* Method 2: w relaxation — exponential smoothing toward target
         * w²_actual = ow2 + α*(w²_target - ow2),  α = min(1, dt/τ)
         * Starts from w=0 (ow2=0 at first step) → gradual buildup */
        if(w_relax_tau > 0 && prev_dt > 0){
            double ci = fmax(P[i].c, 1e-10);
            double tau = w_relax_tau * meand / ci;
            double alpha_r = fmin(1.0, prev_dt / tau);
            double w2_old = (ow2[i] > 0) ? ow2[i] : 0.0;
            w2 = w2_old + alpha_r * (w2 - w2_old);
        }

        w2 = fmin(w2, ceil2);
        if(w2_rate_max > 0 && prev_dt > 0 && ow2[i] > 0){
            double w2_hi = ow2[i] * (1.0 + w2_rate_max);
            double w2_lo = ow2[i] * (1.0 - w2_rate_max);
            if(w2_lo < 0) w2_lo = 0;
            if(w2 > w2_hi) w2 = w2_hi;
            if(w2 < w2_lo) w2 = w2_lo;
        }
        P[i].w = sqrt(w2);
    }

    /* Method 3: w gradient limiter — reduce w when neighbor difference too large
     * Iterative: limit |w_i - w_{i+1}| / dx ≤ w_grad_max */
    if(w_grad_max > 0){
        for(int pass=0; pass<3; pass++){
            for(int i=0;i<N-1;i++){
                double dx = P[i+1].x - P[i].x;
                if(dx < 1e-15) dx = 1e-15;
                double wi = P[i].w, wj = P[i+1].w;
                double dw = fabs(wi - wj);
                double maxdw = w_grad_max * dx;
                if(dw > maxdw){
                    double mean_w = 0.5*(wi + wj);
                    double scale = maxdw / dw;
                    P[i].w   = mean_w + scale * (wi - mean_w);
                    P[i+1].w = mean_w + scale * (wj - mean_w);
                    if(P[i].w   < 0) P[i].w   = 0;
                    if(P[i+1].w < 0) P[i+1].w = 0;
                }
            }
        }
    }
    for(int i=0;i<N-1;i++){
        double dx = P[i+1].x - P[i].x;
        if(dx < 1e-15) dx = 1e-15;
        double w2i = P[i].w*P[i].w, w2j = P[i+1].w*P[i+1].w;
        double shift = (w2i-w2j)/(2*dx);
        double max_s = face_shift_max * dx;
        if(shift >  max_s) shift =  max_s;
        if(shift < -max_s) shift = -max_s;
        face[i+1] = 0.5*(P[i].x+P[i+1].x) + shift;
        if(face[i+1] < P[i].x  +1e-12) face[i+1] = P[i].x  +1e-12;
        if(face[i+1] > P[i+1].x-1e-12) face[i+1] = P[i+1].x-1e-12;
    }
    for(int i=0;i<N;i++){
        P[i].vol=face[i+1]-face[i];
        if(P[i].vol<1e-15) P[i].vol=1e-15;
        P[i].rho=P[i].m/P[i].vol;
        P[i].p=eos_p(P[i].rho,P[i].e);
        P[i].c=eos_c(P[i].rho,P[i].p);
    }
}

/* forward declarations for face Riemann solvers (defined after exact_riemann) */
static void exact_riemann_face(double rhoL, double pL, double vL, double cL,
                               double rhoR, double pR, double vR, double cR,
                               double *pstar, double *vstar);
static void hllc_face(double rhoL, double pL, double vL, double cL,
                      double rhoR, double pR, double vR, double cR,
                      double *pstar, double *vstar);

/* ============================================================
 *  Face states: Monaghan AV + dp/Z correction + thermal cond.
 * ============================================================ */
static void face_st(void){
    for(int i=1;i<N;i++){
        double denij = 0.5*(P[i-1].rho + P[i].rho);
        double Cs    = 0.5*(P[i-1].c   + P[i].c);
        double dx    = P[i].x - P[i-1].x;
        if(dx < 1e-30) dx = 1e-30;

        /* --- C: MUSCL reconstruction or piecewise constant --- */
        double pL, pR, vL, vR;
        if(use_muscl && mode!=SPH_MODE){
            /* minmod slope limiter */
            double sp_i=0, sp_im=0, sv_i=0, sv_im=0;
            if(i>=2){
                double dp_l = P[i-1].p - P[i-2].p;
                double dp_r = P[i].p   - P[i-1].p;
                if(dp_l*dp_r > 0) sp_im = (fabs(dp_l)<fabs(dp_r))? dp_l : dp_r;
                double dv_l = P[i-1].v - P[i-2].v;
                double dv_r = P[i].v   - P[i-1].v;
                if(dv_l*dv_r > 0) sv_im = (fabs(dv_l)<fabs(dv_r))? dv_l : dv_r;
            }
            if(i<=N-2){
                double dp_l = P[i].p   - P[i-1].p;
                double dp_r = P[i+1].p - P[i].p;
                if(dp_l*dp_r > 0) sp_i = (fabs(dp_l)<fabs(dp_r))? dp_l : dp_r;
                double dv_l = P[i].v   - P[i-1].v;
                double dv_r = P[i+1].v - P[i].v;
                if(dv_l*dv_r > 0) sv_i = (fabs(dv_l)<fabs(dv_r))? dv_l : dv_r;
            }
            pL = P[i-1].p + 0.5*sp_im;
            pR = P[i].p   - 0.5*sp_i;
            vL = P[i-1].v + 0.5*sv_im;
            vR = P[i].v   - 0.5*sv_i;
        } else {
            pL = P[i-1].p;  pR = P[i].p;
            vL = P[i-1].v;  vR = P[i].v;
        }

        /* --- Face states --- */
        if(use_riemann >= 2 && mode!=SPH_MODE){
            double rL = P[i-1].rho, rR = P[i].rho;
            double csL = P[i-1].c, csR = P[i].c;
            if(use_riemann == 2)
                exact_riemann_face(rL, pL, vL, csL, rR, pR, vR, csR,
                                   &pf[i], &vf[i]);
            else
                hllc_face(rL, pL, vL, csL, rR, pR, vR, csR,
                          &pf[i], &vf[i]);
        } else if(use_riemann == 1 && mode!=SPH_MODE){
            /* Full acoustic Riemann with separate Z_L, Z_R */
            double ZL = fmax(P[i-1].rho * P[i-1].c, 1e-30);
            double ZR = fmax(P[i].rho   * P[i].c,   1e-30);
            double Zsum = ZL + ZR;
            pf[i] = (ZR*pL + ZL*pR + ZL*ZR*(vL - vR)) / Zsum;
            vf[i] = (ZL*vL + ZR*vR + pL - pR) / Zsum;
        } else {
            pf[i] = 0.5*(pL + pR);
            vf[i] = 0.5*(vL + vR);
            /* partial Riemann corrections */
            if((face_dp > 0 || face_dv > 0) && mode!=SPH_MODE){
                double Z = fmax(denij*Cs, 1e-30);
                if(face_dp > 0) vf[i] += face_dp * (pL - pR) / Z;
                if(face_dv > 0) pf[i] += face_dv * (vL - vR) * Z;
            }
        }

        /* AA correction: df/dt for Laguerre face */
        if(use_AA && mode==LAGUERRE){
            double w2p = P[i-1].w * P[i-1].w;
            double w2q = P[i].w   * P[i].w;
            double vpq = P[i].v - P[i-1].v;
            double dw2p_dt=0, dw2q_dt=0;

            if(use_AA==1 && prev_dt > 0){
                dw2p_dt = (w2p - ow2[i-1]) / prev_dt;
                dw2q_dt = (w2q - ow2[i])   / prev_dt;
            } else if(use_AA==2){
                /* analytical dw^2/dt = -beta * w^2 * div(v)
                 * PRESSURE:  beta = 2*(gamma-1)
                 * ENTHALPY:  beta = 2*(gamma-1)
                 * CSOUND:    beta = (gamma-1)
                 * SPECVOL:   beta = -2              */
                double beta;
                switch(weight_mode){
                case WM_CSOUND:  beta = GM1;       break;
                case WM_SPECVOL: beta = -2.0;      break;
                default:         beta = 2.0 * GM1;  break; /* PRESSURE, ENTHALPY */
                }
                dw2p_dt = -beta * w2p * divv[i-1];
                dw2q_dt = -beta * w2q * divv[i];
            }

            double AA = (dw2p_dt - dw2q_dt) / (2*dx)
                      - (w2p - w2q) * vpq / (2*dx*dx);
            if(aa_max_frac > 0){
                double AA_max = aa_max_frac * Cs;
                if(AA >  AA_max) AA =  AA_max;
                if(AA < -AA_max) AA = -AA_max;
            }
            vf[i] += AA;
        }

        /* Artificial viscosity (compression only) */
        if(av_mode == 1){
            double dv_star = vL - vR;
            if(dv_star > 0){
                double v_sig = P[i-1].c + P[i].c + av_beta * dv_star;
                double av_price = 0.5 * av_alpha * v_sig * denij * dv_star;
                if(mode==LAGUERRE){
                    double wp = P[i-1].w, wq = P[i].w, wsum = wp+wq;
                    if(wsum > 0.5*fabs(dx))
                        pf[i] += sqrt(fabs(dx)/wsum) * av_price;
                    else
                        pf[i] += av_price;
                } else {
                    pf[i] += av_price;
                }
            }
        } else if(av_mode == 2){
            double dv_star = vL - vR;
            double v_sig = P[i-1].c + P[i].c + av_beta * fmax(dv_star, 0.0);
            if(v_sig > vsig_max[i-1]) vsig_max[i-1] = v_sig;
            if(v_sig > vsig_max[i])   vsig_max[i]   = v_sig;
            if(dv_star > 0){
                double alpha_face = 0.5 * (alpha_cd[i-1] + alpha_cd[i]);
                double av_cd = 0.5 * alpha_face * v_sig * denij * dv_star;
                if(mode==LAGUERRE){
                    double wp = P[i-1].w, wq = P[i].w, wsum = wp+wq;
                    if(wsum > 0.5*fabs(dx))
                        pf[i] += sqrt(fabs(dx)/wsum) * av_cd;
                    else
                        pf[i] += av_cd;
                } else {
                    pf[i] += av_cd;
                }
            }
        } else {
            /* Monaghan AV */
            double dv_half = 0.5*(P[i].v - P[i-1].v);
            double rvel    = dv_half * dx;
            double rscale  = 0.25*(P[i-1].vol + P[i].vol);
            double muij    = rvel * rscale / (dx*dx + av_eps*av_eps*rscale*rscale);

            if(muij < 0){
                double al = av_alpha, bt = av_beta;
                if(mode==LAGUERRE && shock_th > 0){
                    double w2L = P[i-1].w * P[i-1].w;
                    double w2R = P[i].w   * P[i].w;
                    double w2avg = 0.5*(w2L + w2R);
                    double sind = (w2avg > 1e-30) ?
                        fabs(w2L - w2R) / w2avg : 0.0;
                    double fac = fmin(sind / shock_th, 1.0);
                    al = av_al_mn + (av_alpha - av_al_mn) * fac;
                    bt = av_bt_mn + (av_beta  - av_bt_mn) * fac;
                }
                double av_base = (-al*Cs*muij + bt*muij*muij)*denij;
                if(mode==LAGUERRE){
                    double wp = P[i-1].w, wq = P[i].w, wsum = wp+wq;
                    if(wsum > 0.5*fabs(dx))
                        pf[i] += sqrt(fabs(dx)/wsum) * av_base;
                    else
                        pf[i] += av_base;
                } else {
                    pf[i] += av_base;
                }
            }
        }

        /* NS viscous stress: tau_xx = (4/3)*nu*rho*dv/dx
         * p_eff = p_face - tau_xx  (tension in expansion, extra pressure in compression) */
        {
            double nu = nu_phys;
            if(nu_coeff > 0)
                nu += nu_coeff * 0.5*(P[i-1].vol + P[i].vol) * Cs;
            if(nu > 0){
                double dvdx = (P[i].v - P[i-1].v) / dx;
                pf[i] -= (4.0/3.0) * nu * denij * dvdx;
            }
        }

        /* Floor: allow negative pf when NS stress is active (physical tension) */
        if(pf[i] < 0 && nu_phys <= 0 && nu_coeff <= 0) pf[i] = 0;
        qf[i] = 0;
    }
    /* Boundary conditions */
    if(bc_left == BC_REFLECT){
        pf[0] = P[0].p   - P[0].rho  *P[0].c  *P[0].v;   vf[0]=0;
    } else {
        pf[0] = P[0].p;   vf[0] = P[0].v;
    }
    if(bc_right == BC_REFLECT){
        pf[N] = P[N-1].p + P[N-1].rho*P[N-1].c*P[N-1].v;  vf[N]=0;
    } else {
        pf[N] = P[N-1].p;  vf[N] = P[N-1].v;
    }
    if(pf[0]<0) pf[0]=0;
    if(pf[N]<0) pf[N]=0;
    qf[0] = 0;  qf[N] = 0;
}

/* ============================================================
 *  RHS: mesh modes
 * ============================================================ */
static void rhs_mesh(void){
    if(mode==VORONOI) cells_vor(); else cells_lag();
    face_st();
    for(int i=0;i<N;i++){
        P[i].av = (pf[i] - pf[i+1]) / P[i].m;
        P[i].ae = (pf[i]*(vf[i]-P[i].v) - pf[i+1]*(vf[i+1]-P[i].v)) / P[i].m;
        P[i].ae += (qf[i] - qf[i+1]) / P[i].m;
        divv[i] = (vf[i+1] - vf[i]) / fmax(P[i].vol, 1e-15);
    }
}

/* ============================================================
 *  SPH
 * ============================================================ */
static void sph_density(void){
    for(int i=0;i<N;i++){
        hh[i]=SPH_ETA*P[i].m/fmax(P[i].rho,1e-10);
        if(hh[i]<sph_hmin) hh[i]=sph_hmin;
        if(hh[i]>sph_hmax) hh[i]=sph_hmax;
    }
    for(int i=0;i<N;i++){
        double rho=0;
        for(int j=0;j<N;j++){
            double hij=0.5*(hh[i]+hh[j]);
            rho+=P[j].m*Wk(P[i].x - P[j].x, hij);
            if(bc_left == BC_REFLECT)
                rho+=P[j].m*Wk(P[i].x - (2*x0_dom - P[j].x), hij);
            if(bc_right == BC_REFLECT)
                rho+=P[j].m*Wk(P[i].x - (2*x1_dom - P[j].x), hij);
        }
        P[i].rho=fmax(rho,1e-10);
        P[i].p=eos_p(P[i].rho,P[i].e);
        P[i].c=eos_c(P[i].rho,P[i].p);
        P[i].vol=P[i].m/P[i].rho;
    }
    for(int i=0;i<N;i++){
        hh[i]=SPH_ETA*P[i].m/P[i].rho;
        if(hh[i]<sph_hmin) hh[i]=sph_hmin;
        if(hh[i]>sph_hmax) hh[i]=sph_hmax;
    }
}

static void sph_pair(int i,double xj,double vj,double mj,
                     double rhoj,double pj,double cj,double hj)
{
    double xij=P[i].x-xj, vij=P[i].v-vj;
    double hij=0.5*(hh[i]+hj);
    double dw=dWk(xij,hij);
    if(dw==0) return;
    double Pi2=P[i].p/(P[i].rho*P[i].rho);
    double Pj2=pj/(rhoj*rhoj);
    double Piav=0;
    if(vij*xij<0){
        double mu=hij*vij*xij/(xij*xij+0.01*hij*hij);
        double cb=0.5*(P[i].c+cj);
        double rb=0.5*(P[i].rho+rhoj);
        Piav=(-SPH_A*cb*mu+SPH_B*mu*mu)/rb;
    }
    double T=Pi2+Pj2+Piav;
    P[i].av -= mj*T*dw;
    P[i].ae += 0.5*mj*T*vij*dw;
}

static void rhs_sph(void){
    sph_density();
    for(int i=0;i<N;i++){ P[i].av=0; P[i].ae=0; }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            if(i!=j)
                sph_pair(i, P[j].x,P[j].v,P[j].m,
                         P[j].rho,P[j].p,P[j].c,hh[j]);
            if(bc_left == BC_REFLECT)
                sph_pair(i, 2*x0_dom-P[j].x,-(P[j].v),P[j].m,
                         P[j].rho,P[j].p,P[j].c,hh[j]);
            if(bc_right == BC_REFLECT)
                sph_pair(i, 2*x1_dom-P[j].x,-(P[j].v),P[j].m,
                         P[j].rho,P[j].p,P[j].c,hh[j]);
        }
    }
}

/* ===== Dispatcher ===== */
static void compute_rhs(void){
    if(mode==SPH_MODE) rhs_sph(); else rhs_mesh();
}
static void update_derived(void){
    if(mode==VORONOI) cells_vor();
    else if(mode==LAGUERRE) cells_lag();
    else sph_density();
}

/* ============================================================
 *  Cullen-Dehnen alpha update (called once per RK4 step)
 * ============================================================ */
static void update_alpha_cd(double dt){
    for(int i=0; i<N; i++){
        double h = P[i].vol;
        double A = fmax(0.0, -(divv[i] - divv_old[i]) / dt);
        double vs = vsig_max[i];
        double alpha_loc = cd_amax * h*h * A / (vs*vs + h*h * A + 1e-30);
        double tau = h / (2.0 * cd_ell * vs + 1e-30);
        if(alpha_loc > alpha_cd[i])
            alpha_cd[i] = alpha_loc;
        else
            alpha_cd[i] = alpha_loc + (alpha_cd[i] - alpha_loc) * exp(-dt/tau);
        if(alpha_cd[i] < cd_amin) alpha_cd[i] = cd_amin;
        if(alpha_cd[i] > cd_amax) alpha_cd[i] = cd_amax;
        divv_old[i] = divv[i];
    }
}

/* ============================================================
 *  CFL timestep
 * ============================================================ */
static double get_dt(void){
    double dt=1e30;
    for(int i=0;i<N;i++){
        double s=P[i].c+fabs(P[i].v);
        if(s<1e-30) continue;
        double dx=(mode==SPH_MODE)?hh[i]:P[i].vol;
        double d=CFL*dx/s;
        if(d<dt) dt=d;
    }
    return dt;
}

/* ============================================================
 *  RK4 step
 * ============================================================ */
static void rk4_step(double dt){
    if(mode==LAGUERRE && w_coeff > 0){
        for(int i=0;i<N;i++) ow2[i] = P[i].w * P[i].w;
        prev_dt = dt;
    }
    for(int i=0;i<N;i++){
        s0x[i]=P[i].x; s0v[i]=P[i].v; s0e[i]=P[i].e;
    }
    for(int st=0;st<4;st++){
        if(st>0){
            double f=(st<3)?0.5:1.0;
            int k=st-1;
            for(int i=0;i<N;i++){
                P[i].x=s0x[i]+f*dt*kx[k][i];
                P[i].v=s0v[i]+f*dt*kv[k][i];
                P[i].e=s0e[i]+f*dt*ke[k][i];
                if(P[i].e<1e-14) P[i].e=1e-14;
            }
        }
        if(av_mode==2)
            for(int i=0;i<N;i++) vsig_max[i]=0;
        compute_rhs();
        for(int i=0;i<N;i++){
            kx[st][i]=P[i].v;
            kv[st][i]=P[i].av;
            ke[st][i]=P[i].ae;
        }
    }
    for(int i=0;i<N;i++){
        P[i].x=s0x[i]+dt/6.0*(kx[0][i]+2*kx[1][i]+2*kx[2][i]+kx[3][i]);
        P[i].v=s0v[i]+dt/6.0*(kv[0][i]+2*kv[1][i]+2*kv[2][i]+kv[3][i]);
        P[i].e=s0e[i]+dt/6.0*(ke[0][i]+2*ke[1][i]+2*ke[2][i]+ke[3][i]);
        if(P[i].e<1e-14) P[i].e=1e-14;
    }
    update_derived();
    if(av_mode==2)
        update_alpha_cd(dt);
}

/* ============================================================
 *  Problem setup
 * ============================================================ */
static void setup_problem(int pid){
    problem_id = pid;
    switch(pid){
    case PROB_SOD:
        GAM=1.4; GM1=0.4; GP1=2.4;
        x0_dom=0; x1_dom=1; xmid1=0.5;
        n_regions=2;
        rhoL_ic=1; vL_ic=0; pL_ic=1;
        rhoR_ic=0.125; vR_ic=0; pR_ic=0.1;
        bc_left=BC_REFLECT; bc_right=BC_REFLECT;
        tend=0.2;
        has_exact_sol=1;
        break;
    case PROB_BLAST:
        GAM=1.4; GM1=0.4; GP1=2.4;
        x0_dom=0; x1_dom=1; xmid1=0.1; xmid2=0.9;
        n_regions=3;
        rhoL_ic=1; vL_ic=0; pL_ic=1000;
        rhoM_ic=1; vM_ic=0; pM_ic=0.01;
        rhoR_ic=1; vR_ic=0; pR_ic=100;
        bc_left=BC_REFLECT; bc_right=BC_REFLECT;
        tend=0.038;
        has_exact_sol=0;
        break;
    case PROB_SHUOSHER:
        GAM=1.4; GM1=0.4; GP1=2.4;
        x0_dom=-5; x1_dom=5; xmid1=-4;
        n_regions=2;
        rhoL_ic=3.857143; vL_ic=2.629369; pL_ic=10.33333;
        rhoR_ic=1; vR_ic=0; pR_ic=1;
        bc_left=BC_TRANSMIT; bc_right=BC_TRANSMIT;
        tend=1.8;
        has_exact_sol=0;
        break;
    case PROB_NOH:
        GAM=5.0/3.0; GM1=2.0/3.0; GP1=8.0/3.0;
        x0_dom=0; x1_dom=1; xmid1=1.0; /* single region */
        n_regions=1;
        rhoL_ic=1; vL_ic=-1; pL_ic=1e-6;
        rhoR_ic=1; vR_ic=-1; pR_ic=1e-6;
        bc_left=BC_REFLECT; bc_right=BC_TRANSMIT;
        tend=0.6;
        has_exact_sol=1;
        break;
    case PROB_LAX:
        GAM=1.4; GM1=0.4; GP1=2.4;
        x0_dom=0; x1_dom=1; xmid1=0.5;
        n_regions=2;
        rhoL_ic=0.445; vL_ic=0.698; pL_ic=3.528;
        rhoR_ic=0.5;   vR_ic=0;     pR_ic=0.571;
        bc_left=BC_TRANSMIT; bc_right=BC_TRANSMIT;
        tend=0.13;
        has_exact_sol=1;
        break;
    case PROB_DBLFAN:
        GAM=1.4; GM1=0.4; GP1=2.4;
        x0_dom=0; x1_dom=1; xmid1=0.5;
        n_regions=2;
        rhoL_ic=1; vL_ic=-2; pL_ic=0.4;
        rhoR_ic=1; vR_ic=2;  pR_ic=0.4;
        bc_left=BC_TRANSMIT; bc_right=BC_TRANSMIT;
        tend=0.15;
        has_exact_sol=1;
        break;
    case PROB_COLLISION:
        GAM=1.4; GM1=0.4; GP1=2.4;
        x0_dom=0; x1_dom=1; xmid1=0.4;
        n_regions=2;
        rhoL_ic=5.99924; vL_ic=19.5975;  pL_ic=460.894;
        rhoR_ic=5.99242; vR_ic=-6.19633; pR_ic=46.095;
        bc_left=BC_TRANSMIT; bc_right=BC_TRANSMIT;
        tend=0.035;
        has_exact_sol=1;
        break;
    case PROB_CONTACT:
        GAM=1.4; GM1=0.4; GP1=2.4;
        x0_dom=0; x1_dom=1; xmid1=0.5;
        n_regions=2;
        rhoL_ic=1.0;   vL_ic=0; pL_ic=1.0;
        rhoR_ic=0.125; vR_ic=0; pR_ic=1.0;
        bc_left=BC_TRANSMIT; bc_right=BC_TRANSMIT;
        tend=0.25;
        has_exact_sol=1;
        break;
    }
    sph_hmin = 0.003 * (x1_dom - x0_dom);
    sph_hmax = 0.1   * (x1_dom - x0_dom);

    /* adaptive reference values for weight scaling */
    double pmax_ic = fmax(pL_ic, pR_ic);
    if(pid == PROB_BLAST) pmax_ic = fmax(pmax_ic, pM_ic);
    double pmin_ic = fmin(pL_ic, pR_ic);
    if(pid == PROB_BLAST) pmin_ic = fmin(pmin_ic, pM_ic);
    p0_ref = sqrt(pmax_ic * fmax(pmin_ic, 1e-14));  /* geometric mean */

    double hL = eos_e(rhoL_ic, pL_ic) + pL_ic/fmax(rhoL_ic,1e-14);
    double hR = eos_e(rhoR_ic, pR_ic) + pR_ic/fmax(rhoR_ic,1e-14);
    h0_ref = 0.5*(hL + hR);
    if(h0_ref < 1e-10) h0_ref = 1e-10;

    v0_ref = 0.5*(1.0/fmax(rhoL_ic,1e-14) + 1.0/fmax(rhoR_ic,1e-14));

    double cL_ic = eos_c(rhoL_ic, pL_ic);
    double cR_ic = eos_c(rhoR_ic, pR_ic);
    c0_ref = sqrt(cL_ic * cR_ic);  /* geometric mean */
    if(c0_ref < 1e-10) c0_ref = 1e-10;
    if(v0_ref < 1e-10) v0_ref = 1e-10;
}

/* ============================================================
 *  Initialization
 * ============================================================ */
static int ic_type = 1;  /* default: equal-mass */
static int NL_actual, NR_actual;

static void init_two_region(int n){
    N = n;
    if(ic_type == 1){
        double total_mass = rhoL_ic*(xmid1-x0_dom) + rhoR_ic*(x1_dom-xmid1);
        double mp = total_mass / N;
        double dxL = mp / rhoL_ic;
        double dxR = (rhoR_ic > 1e-30) ? mp / rhoR_ic : 0;
        int NL = (int)(0.5 + (xmid1-x0_dom)/dxL);
        if(NL > N-1 && rhoR_ic > 1e-30) NL = N-1;
        if(NL > N) NL = N;
        if(NL < 1) NL = 1;
        int NR = N - NL;
        NL_actual = NL; NR_actual = NR;
        meand = (x1_dom-x0_dom)/N;

        for(int i=0; i<NL; i++){
            P[i].x = x0_dom + (i+0.5)*dxL;
            P[i].v = vL_ic; P[i].m = mp;
            P[i].rho = rhoL_ic; P[i].p = pL_ic;
            P[i].e = eos_e(rhoL_ic, pL_ic);
            P[i].c = eos_c(rhoL_ic, pL_ic);
            P[i].vol = dxL; P[i].w = 0;
            P[i].av = P[i].ae = 0;
        }
        for(int i=0; i<NR; i++){
            int k = NL + i;
            P[k].x = xmid1 + (i+0.5)*dxR;
            P[k].v = vR_ic; P[k].m = mp;
            P[k].rho = rhoR_ic; P[k].p = pR_ic;
            P[k].e = eos_e(rhoR_ic, pR_ic);
            P[k].c = eos_c(rhoR_ic, pR_ic);
            P[k].vol = dxR; P[k].w = 0;
            P[k].av = P[k].ae = 0;
        }
    } else {
        double dx = (x1_dom-x0_dom)/N;
        meand = dx;
        NL_actual = N/2; NR_actual = N - NL_actual;
        for(int i=0; i<N; i++){
            P[i].x = x0_dom+(i+0.5)*dx;
            P[i].v = 0; P[i].av = P[i].ae = 0;
            if(P[i].x < xmid1){
                P[i].rho=rhoL_ic; P[i].p=pL_ic; P[i].m=rhoL_ic*dx;
                P[i].v=vL_ic;
            } else {
                P[i].rho=rhoR_ic; P[i].p=pR_ic; P[i].m=rhoR_ic*dx;
                P[i].v=vR_ic;
            }
            P[i].e = eos_e(P[i].rho, P[i].p);
            P[i].c = eos_c(P[i].rho, P[i].p);
            P[i].vol = dx; P[i].w = 0;
        }
    }
    memset(divv, 0, sizeof(double)*N);
    memset(divv_old, 0, sizeof(double)*N);
    memset(vsig_max, 0, sizeof(double)*N);
    for(int i=0;i<N;i++) alpha_cd[i] = cd_amin;
    update_derived();
}

static void init_three_region(int n){
    N = n;
    double total_mass = rhoL_ic*(xmid1-x0_dom) + rhoM_ic*(xmid2-xmid1)
                      + rhoR_ic*(x1_dom-xmid2);
    double mp = total_mass / N;
    double dxL = mp/rhoL_ic, dxM = mp/rhoM_ic, dxR = mp/rhoR_ic;
    int NL = (int)(0.5 + (xmid1-x0_dom)/dxL);
    int NM = (int)(0.5 + (xmid2-xmid1)/dxM);
    if(NL < 0) NL = 0;
    if(NM < 0) NM = 0;
    int NR = N - NL - NM;
    if(NR < 0){ NM += NR; NR = 0; }
    NL_actual = NL; NR_actual = NR;
    meand = (x1_dom - x0_dom) / N;

    int k = 0;
    for(int i=0; i<NL; i++){
        P[k].x = x0_dom + (i+0.5)*dxL;
        P[k].v = vL_ic; P[k].m = mp;
        P[k].rho = rhoL_ic; P[k].p = pL_ic;
        P[k].e = eos_e(rhoL_ic, pL_ic);
        P[k].c = eos_c(rhoL_ic, pL_ic);
        P[k].vol = dxL; P[k].w = 0;
        P[k].av = P[k].ae = 0;
        k++;
    }
    for(int i=0; i<NM; i++){
        P[k].x = xmid1 + (i+0.5)*dxM;
        P[k].v = vM_ic; P[k].m = mp;
        P[k].rho = rhoM_ic; P[k].p = pM_ic;
        P[k].e = eos_e(rhoM_ic, pM_ic);
        P[k].c = eos_c(rhoM_ic, pM_ic);
        P[k].vol = dxM; P[k].w = 0;
        P[k].av = P[k].ae = 0;
        k++;
    }
    for(int i=0; i<NR; i++){
        P[k].x = xmid2 + (i+0.5)*dxR;
        P[k].v = vR_ic; P[k].m = mp;
        P[k].rho = rhoR_ic; P[k].p = pR_ic;
        P[k].e = eos_e(rhoR_ic, pR_ic);
        P[k].c = eos_c(rhoR_ic, pR_ic);
        P[k].vol = dxR; P[k].w = 0;
        P[k].av = P[k].ae = 0;
        k++;
    }
    memset(divv, 0, sizeof(double)*N);
    memset(divv_old, 0, sizeof(double)*N);
    memset(vsig_max, 0, sizeof(double)*N);
    for(int i=0; i<N; i++) alpha_cd[i] = cd_amin;
    update_derived();
}

static void init_shuosher(int n){
    N = n;
    double mass_L = rhoL_ic * (xmid1 - x0_dom);
    /* ∫_{-4}^{5} (1+0.2sin(5πx)) dx = 9 + 0.2/(5π)[cos(-20π)-cos(25π)]
       cos(-20π)=1, cos(25π)=-1 → 9 + 0.2/(5π)*2 = 9 + 0.08/π */
    double mass_R = (x1_dom - xmid1)
                  + 0.2/(5.0*M_PI)*(cos(5.0*M_PI*xmid1) - cos(5.0*M_PI*x1_dom));
    double total_mass = mass_L + mass_R;
    double mp = total_mass / N;

    int NL = (int)(0.5 + mass_L / mp);
    if(NL < 1) NL = 1;
    if(NL > N-1) NL = N-1;
    int NR = N - NL;
    NL_actual = NL; NR_actual = NR;
    meand = (x1_dom - x0_dom) / N;

    /* Left region: uniform density */
    double dxL = (xmid1 - x0_dom) / NL;
    for(int i=0; i<NL; i++){
        P[i].x = x0_dom + (i+0.5)*dxL;
        P[i].v = vL_ic; P[i].m = mp;
        P[i].rho = rhoL_ic; P[i].p = pL_ic;
        P[i].e = eos_e(rhoL_ic, pL_ic);
        P[i].c = eos_c(rhoL_ic, pL_ic);
        P[i].vol = dxL; P[i].w = 0;
        P[i].av = P[i].ae = 0;
    }

    /* Right region: sinusoidal density, Newton iteration for equal-mass */
    /* F(x) = ∫_{xmid1}^{x} (1+0.2sin(5πs)) ds
            = (x-xmid1) + 0.2/(5π)[cos(5π·xmid1) - cos(5πx)] */
    double cos0 = cos(5.0*M_PI*xmid1);
    double x_prev = xmid1;
    for(int j=0; j<NR; j++){
        double target = (j+0.5) * mp;
        double x = x_prev + mp;  /* initial guess */
        for(int it=0; it<100; it++){
            double Fx = (x - xmid1) + 0.2/(5.0*M_PI)*(cos0 - cos(5.0*M_PI*x));
            double dFx = 1.0 + 0.2*sin(5.0*M_PI*x);
            if(dFx < 1e-30) dFx = 1e-30;
            double dx = -(Fx - target) / dFx;
            x += dx;
            if(fabs(dx) < 1e-14) break;
        }
        int k = NL + j;
        double rho_local = 1.0 + 0.2*sin(5.0*M_PI*x);
        if(rho_local < 0.1) rho_local = 0.1;
        P[k].x = x;
        P[k].v = vR_ic; P[k].m = mp;
        P[k].rho = rho_local; P[k].p = pR_ic;
        P[k].e = eos_e(rho_local, pR_ic);
        P[k].c = eos_c(rho_local, pR_ic);
        P[k].vol = mp / rho_local; P[k].w = 0;
        P[k].av = P[k].ae = 0;
        x_prev = x;
    }

    memset(divv, 0, sizeof(double)*N);
    memset(divv_old, 0, sizeof(double)*N);
    memset(vsig_max, 0, sizeof(double)*N);
    for(int i=0; i<N; i++) alpha_cd[i] = cd_amin;
    update_derived();
}

static void init(int n){
    if(problem_id == PROB_BLAST)
        init_three_region(n);
    else if(problem_id == PROB_SHUOSHER)
        init_shuosher(n);
    else
        init_two_region(n);
}

/* ============================================================
 *  Exact Riemann solver (generalized for runtime IC vars)
 * ============================================================ */
static double rfunc(double p,double rK,double pK,double cK){
    if(p<=pK) return(2*cK/GM1)*(pow(p/pK,GM1/(2*GAM))-1);
    double A=2.0/(GP1*rK), B=GM1/GP1*pK;
    return(p-pK)*sqrt(A/(p+B));
}
static double rfunc_d(double p,double rK,double pK,double cK){
    if(p<=pK) return(1.0/(rK*cK))*pow(p/pK,-GP1/(2*GAM));
    double A=2.0/(GP1*rK), B=GM1/GP1*pK;
    double s=sqrt(A/(p+B));
    return s*(1-0.5*(p-pK)/(p+B));
}
static void exact_riemann(double *ps,double *vs){
    double cL=eos_c(rhoL_ic,pL_ic), cR=eos_c(rhoR_ic,pR_ic);
    double p=0.5*(pL_ic+pR_ic);
    for(int it=0;it<200;it++){
        double f=rfunc(p,rhoL_ic,pL_ic,cL)+rfunc(p,rhoR_ic,pR_ic,cR)+(vR_ic-vL_ic);
        double fp=rfunc_d(p,rhoL_ic,pL_ic,cL)+rfunc_d(p,rhoR_ic,pR_ic,cR);
        if(fabs(fp)<1e-30) break;
        double dp=-f/fp; p+=dp;
        if(p<1e-10)p=1e-10;
        if(fabs(dp)<1e-12*(1+p)) break;
    }
    *ps=p;
    *vs=0.5*(vL_ic+vR_ic)+0.5*(rfunc(p,rhoR_ic,pR_ic,cR)-rfunc(p,rhoL_ic,pL_ic,cL));
}

/* ============================================================
 *  Face-level exact Riemann solver
 * ============================================================ */
static void exact_riemann_face(double rhoL, double pL, double vL, double cL,
                               double rhoR, double pR, double vR, double cR,
                               double *pstar, double *vstar)
{
    double z = GM1 / (2.0*GAM);
    double pguess = pow((cL + cR - 0.5*GM1*(vR - vL)) /
                        (cL/pow(pL,z) + cR/pow(pR,z)), 1.0/z);
    if(pguess < 1e-10) pguess = 1e-10;

    double p = pguess;
    int converged = 0;
    for(int it = 0; it < 50; it++){
        double fL = rfunc(p, rhoL, pL, cL);
        double fR = rfunc(p, rhoR, pR, cR);
        double f  = fL + fR + (vR - vL);
        double fp = rfunc_d(p, rhoL, pL, cL) + rfunc_d(p, rhoR, pR, cR);
        if(fabs(fp) < 1e-30) break;
        double dp = -f / fp;
        p += dp;
        if(p < 1e-10) p = 1e-10;
        if(fabs(dp) < 1e-8*(1.0 + p)){ converged = 1; break; }
    }

    if(!converged){
        double ZL = rhoL * cL, ZR = rhoR * cR;
        double Zs = ZL + ZR;
        *pstar = (ZR*pL + ZL*pR + ZL*ZR*(vL - vR)) / Zs;
        *vstar = (ZL*vL + ZR*vR + pL - pR) / Zs;
        return;
    }

    *pstar = p;
    *vstar = 0.5*(vL + vR) + 0.5*(rfunc(p, rhoR, pR, cR) - rfunc(p, rhoL, pL, cL));
}

/* ============================================================
 *  HLLC face solver (TRUE HLLC with Gaburov/Davis wave speeds +
 *  Roe-averaged fallback + PVRS-Rusanov last resort, ports GIZMO
 *  reimann.h:527-593.  Previous impl here was TSRS disguised as HLLC.)
 * ============================================================ */
static void hllc_face(double rhoL, double pL, double vL, double cL,
                      double rhoR, double pR, double vR, double cR,
                      double *pstar, double *vstar)
{
    const double tiny = 1.0e-30;
    double S_L, S_R, S_M, P_M;
    double cmax = cL > cR ? cL : cR;

    if((vR - vL) > cmax){
        *pstar = tiny; *vstar = 0.5*(vL + vR);
        return;
    }

    /* Primary: Gaburov/Davis */
    {
        double vmin = vL < vR ? vL : vR;
        double vmax = vL > vR ? vL : vR;
        S_L = vmin - cmax;
        S_R = vmax + cmax;
        double rho_wt_L = rhoL*(S_L - vL);
        double rho_wt_R = rhoR*(S_R - vR);
        double denom = rho_wt_L - rho_wt_R;
        if(fabs(denom) > tiny){
            S_M = ((pR - pL) + rho_wt_L*vL - rho_wt_R*vR) / denom;
            P_M = (pL*rho_wt_R - pR*rho_wt_L + rho_wt_L*rho_wt_R*(vR - vL))
                  / (rho_wt_R - rho_wt_L);
            if(P_M > 0 && !isnan(P_M) && !isinf(P_M)){
                *pstar = P_M; *vstar = S_M; return;
            }
        }
    }

    /* Roe fallback */
    {
        double sqL = sqrt(rhoL), sqR = sqrt(rhoR);
        double sq_inv = 1.0/(sqL + sqR + tiny);
        double v_roe = (sqL*vL + sqR*vR)*sq_inv;
        double h_L = GAM*pL/((GAM - 1.0)*rhoL) + 0.5*vL*vL;
        double h_R = GAM*pR/((GAM - 1.0)*rhoR) + 0.5*vR*vR;
        double h_roe = (sqL*h_L + sqR*h_R)*sq_inv;
        double c2_roe = (GAM - 1.0)*(h_roe - 0.5*v_roe*v_roe);
        if(c2_roe < tiny) c2_roe = tiny;
        double c_roe = sqrt(c2_roe);

        double srA = vR + cR, srB = v_roe + c_roe;
        double slA = vL - cL, slB = v_roe - c_roe;
        S_R = srA > srB ? srA : srB;
        S_L = slA < slB ? slA : slB;

        double rho_wt_R =  rhoR*(S_R - vR);
        double rho_wt_L = -rhoL*(S_L - vL);
        double denom = rho_wt_R + rho_wt_L;
        if(fabs(denom) > tiny){
            S_M = (rho_wt_R*vR + rho_wt_L*vL + (pL - pR))/denom;
            P_M = rhoL*(vL - S_L)*(vL - S_M) + pL;
            if(P_M > 0 && !isnan(P_M) && !isinf(P_M)){
                *pstar = P_M; *vstar = S_M; return;
            }
        }
    }

    /* Last resort */
    {
        P_M = 0.5*((pL + pR) + (vL - vR)*0.25*(rhoL + rhoR)*(cL + cR));
        S_M = 0.5*(vR + vL) + 2.0*(pL - pR)/((rhoL + rhoR)*(cL + cR) + tiny);
        double a1 = fabs(vL - cL), a2 = fabs(vR - cR);
        double a3 = fabs(vL + cL), a4 = fabs(vR + cR);
        double S_plus = a1;
        if(a2 > S_plus) S_plus = a2;
        if(a3 > S_plus) S_plus = a3;
        if(a4 > S_plus) S_plus = a4;
        if(S_M < -S_plus) S_M = -S_plus;
        if(S_M >  S_plus) S_M =  S_plus;
        if(P_M < tiny) P_M = tiny;
        *pstar = P_M; *vstar = S_M;
    }
}

/* ============================================================
 *  Exact sample (generalized for runtime IC vars)
 * ============================================================ */
static void exact_sample(double x,double t,
                         double *rho,double *vel,double *pre)
{
    double ps=ps_exact, vs=vs_exact;
    double cL=eos_c(rhoL_ic,pL_ic), cR=eos_c(rhoR_ic,pR_ic);
    double S=(x-xmid1)/t;
    /* Left wave */
    if(ps<=pL_ic){
        double csL=cL*pow(ps/pL_ic,GM1/(2*GAM));
        if(S<=vL_ic-cL)      { *rho=rhoL_ic; *vel=vL_ic; *pre=pL_ic; return; }
        if(S<=vs-csL){
            double cf=(2.0/GP1)*(cL+0.5*GM1*(vL_ic-S));
            *vel=(2.0/GP1)*(cL+S)+(GM1/GP1)*vL_ic;
            *rho=rhoL_ic*pow(cf/cL,2.0/GM1);
            *pre=pL_ic  *pow(cf/cL,2.0*GAM/GM1);
            return;
        }
    } else {
        double SL=vL_ic-cL*sqrt(GP1/(2*GAM)*ps/pL_ic+GM1/(2*GAM));
        if(S<=SL) { *rho=rhoL_ic; *vel=vL_ic; *pre=pL_ic; return; }
    }
    if(S<=vs){
        *vel=vs; *pre=ps;
        if(ps<=pL_ic) *rho=rhoL_ic*pow(ps/pL_ic,1.0/GAM);
        else { double r=ps/pL_ic; *rho=rhoL_ic*(r+GM1/GP1)/(GM1/GP1*r+1); }
        return;
    }
    /* Right wave */
    if(ps<=pR_ic){
        double csR=cR*pow(ps/pR_ic,GM1/(2*GAM));
        if(S<=vs+csR) { *rho=rhoR_ic*pow(ps/pR_ic,1.0/GAM); *vel=vs; *pre=ps; return; }
        if(S<=vR_ic+cR){
            double cf=(2.0/GP1)*(cR+0.5*GM1*(S-vR_ic));
            *vel=(2.0/GP1)*(-cR+S)+(GM1/GP1)*vR_ic;
            *rho=rhoR_ic*pow(cf/cR,2.0/GM1);
            *pre=pR_ic  *pow(cf/cR,2.0*GAM/GM1);
            return;
        }
        *rho=rhoR_ic; *vel=vR_ic; *pre=pR_ic;
    } else {
        double SR=vR_ic+cR*sqrt(GP1/(2*GAM)*ps/pR_ic+GM1/(2*GAM));
        if(S<=SR){
            double r=ps/pR_ic;
            *rho=rhoR_ic*(r+GM1/GP1)/(GM1/GP1*r+1);
            *vel=vs; *pre=ps;
        } else { *rho=rhoR_ic; *vel=vR_ic; *pre=pR_ic; }
    }
}

/* ============================================================
 *  Noh exact solution
 * ============================================================ */
static void exact_noh(double x, double t, double *rho, double *vel, double *pre){
    double xs = t/3.0;
    if(x < xs){
        *rho = 4.0;
        *vel = 0.0;
        *pre = 4.0/3.0;
    } else {
        *rho = 1.0;
        *vel = -1.0;
        *pre = 0.0;
    }
}

/* ============================================================
 *  Reference solution I/O
 * ============================================================ */
static void read_reference(const char *filename){
    FILE *fp = fopen(filename, "r");
    if(!fp){ fprintf(stderr,"ERROR: cannot open %s\n",filename); nref=0; return; }
    char buf[256];
    nref = 0;
    while(fgets(buf, sizeof(buf), fp)){
        if(buf[0]=='#') continue;
        double x, rho, v, p;
        if(sscanf(buf, "%lf %lf %lf %lf", &x, &rho, &v, &p)==4){
            ref_x[nref] = x; ref_rho[nref] = rho;
            ref_v[nref] = v; ref_p[nref] = p;
            nref++;
            if(nref >= NMAX) break;
        }
    }
    fclose(fp);
}

static double interp_ref(double x, const double *ref_y){
    if(nref < 2) return 0;
    if(x <= ref_x[0]) return ref_y[0];
    if(x >= ref_x[nref-1]) return ref_y[nref-1];
    int lo=0, hi=nref-1;
    while(hi - lo > 1){
        int mid = (lo+hi)/2;
        if(ref_x[mid] <= x) lo = mid; else hi = mid;
    }
    double t = (x - ref_x[lo]) / (ref_x[hi] - ref_x[lo] + 1e-30);
    return ref_y[lo] + t*(ref_y[hi] - ref_y[lo]);
}

static void eval_solution(double x, double t, double *rho, double *vel, double *pre){
    if(problem_id == PROB_SOD || problem_id == PROB_LAX ||
       problem_id == PROB_DBLFAN || problem_id == PROB_COLLISION ||
       problem_id == PROB_CONTACT){
        exact_sample(x, t, rho, vel, pre);
    } else if(problem_id == PROB_NOH){
        exact_noh(x, t, rho, vel, pre);
    } else {
        *rho = interp_ref(x, ref_rho);
        *vel = interp_ref(x, ref_v);
        *pre = interp_ref(x, ref_p);
    }
}

static void write_exact_file(const char *filename){
    FILE *fp = fopen(filename, "w");
    int np = 1001;
    for(int i=0; i<np; i++){
        double x = x0_dom + (x1_dom - x0_dom)*i/(np-1.0);
        double re, ve, pe;
        eval_solution(x, tend, &re, &ve, &pe);
        double ee = (GM1*re > 1e-30) ? pe/(GM1*re) : 0;
        fprintf(fp,"%.10e %.10e %.10e %.10e %.10e\n",x,re,ve,pe,ee);
    }
    fclose(fp);
}

/* ============================================================
 *  Run one simulation, return metrics
 * ============================================================ */
typedef struct {
    double l1_rho, l1_vel, l1_pre;
    double p_osc;
    double eng_err;
    int    nstep;
    int    failed;
} Metrics;

static Metrics run_sim(int smode, int np, const char *outfile, int verbose)
{
    Metrics met = {0};
    mode = smode;
    init(np);

    /* initial total energy (computed numerically) */
    double eng0=0;
    for(int i=0;i<N;i++)
        eng0 += P[i].m*(0.5*P[i].v*P[i].v+P[i].e);

    clock_t t0 = clock();
    double t=0;
    int max_steps = (np > 1000) ? 200000 : 50000;
    int diag_blast = (problem_id==PROB_BLAST && verbose);
    while(t<tend-1e-14){
        double dt=get_dt();
        if(t+dt>tend) dt=tend-t;
        rk4_step(dt);
        t+=dt; met.nstep++;
        if(diag_blast && (met.nstep%2000==0 || met.nstep<=5)){
            double vmin=1e30, vmax=-1e30;
            int imin=-1;
            for(int i=0;i<N;i++){
                if(P[i].vol<vmin){ vmin=P[i].vol; imin=i; }
                if(P[i].vol>vmax) vmax=P[i].vol;
            }
            fprintf(stderr,"    [diag] step=%d t=%.4e dt=%.4e vmin=%.4e(i=%d,x=%.4f,p=%.2e) vmax=%.4e\n",
                    met.nstep, t, dt, vmin, imin, P[imin].x, P[imin].p, vmax);
        }
        if(met.nstep >= max_steps){
            fprintf(stderr,"  WARNING: max steps (%d) at t=%.4e/%.4e\n",
                    max_steps, t, tend);
            if(diag_blast){
                double vmin=1e30; int imin=-1;
                for(int i=0;i<N;i++) if(P[i].vol<vmin){ vmin=P[i].vol; imin=i; }
                fprintf(stderr,"    [diag] FAIL vmin=%.4e i=%d x=%.4f rho=%.4e p=%.4e c=%.4e w=%.4e\n",
                        vmin, imin, P[imin].x, P[imin].rho, P[imin].p, P[imin].c, P[imin].w);
            }
            met.failed = 1;
            break;
        }
    }
    double wall_ms = 1000.0*(clock()-t0)/CLOCKS_PER_SEC;

    /* conservation */
    double eng=0;
    for(int i=0;i<N;i++)
        eng += P[i].m*(0.5*P[i].v*P[i].v+P[i].e);
    met.eng_err = (eng0>1e-30) ? fabs(eng-eng0)/eng0 : fabs(eng-eng0);

    /* L1 error */
    for(int i=0;i<N;i++){
        /* For reference-based problems, skip particles outside reference range */
        if(!has_exact_sol && nref >= 2 &&
           (P[i].x < ref_x[0] || P[i].x > ref_x[nref-1]))
            continue;
        double re,ve,pe;
        eval_solution(P[i].x,tend,&re,&ve,&pe);
        double dv=P[i].vol;
        met.l1_rho += fabs(P[i].rho-re)*dv;
        met.l1_vel += fabs(P[i].v  -ve)*dv;
        met.l1_pre += fabs(P[i].p  -pe)*dv;
    }

    /* pressure oscillation (Sod only) */
    if(problem_id == PROB_SOD){
        double sum=0, sum2=0; int cnt=0;
        for(int i=0;i<N;i++){
            if(P[i].x > 0.69 && P[i].x < 0.82){
                sum  += P[i].p;
                sum2 += P[i].p*P[i].p;
                cnt++;
            }
        }
        if(cnt>1){
            double mean = sum/cnt;
            met.p_osc = sqrt(sum2/cnt - mean*mean);
        }
    }

    if(verbose){
        const char *lbl[]={"Laguerre","Voronoi","SPH"};
        fprintf(stderr,"  %s: %d steps  L1_rho=%.4e  p_osc=%.4e  eng=%.2e  %.1fms\n",
                lbl[smode], met.nstep, met.l1_rho, met.p_osc, met.eng_err, wall_ms);
    }

    /* write output (skip if failed) */
    if(!met.failed && outfile){
        FILE *fp=fopen(outfile,"w");
        fprintf(fp,"# x  rho  v  p  e  w  vol\n");
        for(int i=0;i<N;i++){
            double w_out=(smode==SPH_MODE)?hh[i]:P[i].w;
            fprintf(fp,"%.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
                    P[i].x,P[i].rho,P[i].v,P[i].p,P[i].e,w_out,P[i].vol);
        }
        fclose(fp);
    }
    return met;
}

/* ============================================================
 *  Grid plot: 4 rows (variables) x nmeth columns (methods)
 * ============================================================ */
static void write_plot_grid(const char *gpfile, const char *pngfile,
                            const char *ref_file, const char *title,
                            int nmeth, const char **labels,
                            const char **files, const char **colors,
                            double xlo, double xhi,
                            double ylo4[4], double yhi4[4])
{
    FILE *fp = fopen(gpfile, "w");
    const char *var[] = {"Density", "Velocity", "Pressure", "Internal Energy"};
    int ec[] = {2, 3, 4, 5};

    fprintf(fp,
        "set terminal pngcairo size %d,1200 enhanced font 'Arial,11'\n"
        "set output '%s'\n"
        "set multiplot layout 4,%d title '%s' font ',14'\n"
        "set key top right font ',9'\n\n",
        450*nmeth, pngfile, nmeth, title);

    for(int r = 0; r < 4; r++){
        for(int c = 0; c < nmeth; c++){
            if(r == 0)
                fprintf(fp, "set title '%s' font ',12'\n", labels[c]);
            else
                fprintf(fp, "unset title\n");
            if(c == 0)
                fprintf(fp, "set ylabel '%s'\n", var[r]);
            else
                fprintf(fp, "unset ylabel\n");
            fprintf(fp, "set xlabel 'x'\nset xrange [%g:%g]\n", xlo, xhi);
            fprintf(fp, "set yrange [%.6g:%.6g]\n", ylo4[r], yhi4[r]);
            fprintf(fp,
                "plot '%s' u 1:%d w l lw 2 "
                "lc rgb '#000000' t 'Ref', \\\n"
                "     '%s' u 1:%d w p pt 2 ps 1.0 lw 2 "
                "lc rgb '%s' t '%s'\n\n",
                ref_file, ec[r], files[c], ec[r], colors[c], labels[c]);
        }
    }
    fprintf(fp, "unset multiplot\n");
    fclose(fp);
}

/* ============================================================
 *  Method configuration helpers
 * ============================================================ */
static void set_method_hllc_cd10(void){
    w_coeff=0; use_riemann=3; use_muscl=1; av_mode=2;
    cd_amax=1.0; cd_ell=0.05; cd_amin=0; av_beta=1.0;
    shock_th=0; face_dp=0; face_dv=0; prev_dt=0;
    use_AA=0; alpha_u=0; nu_phys=0; nu_coeff=0;
}
static void set_method_hllc_price(void){
    w_coeff=0; use_riemann=3; use_muscl=1; av_mode=1;
    av_alpha=0.25; av_beta=1.0;
    shock_th=0; face_dp=0; face_dv=0; prev_dt=0;
    use_AA=0; alpha_u=0; nu_phys=0; nu_coeff=0;
}
static void set_method_price(void){
    w_coeff=0; use_riemann=0; use_muscl=1; av_mode=1;
    av_alpha=1.5; av_beta=2.0;
    shock_th=0; face_dp=0; face_dv=0; prev_dt=0;
    use_AA=0; alpha_u=0; nu_phys=0; nu_coeff=0;
}
static void set_method_voronoi(void){
    w_coeff=0; use_riemann=0; use_muscl=0; av_mode=0;
    av_alpha=3; av_beta=4; face_dp=0.5; face_dv=0;
    shock_th=0; prev_dt=0;
    use_AA=0; alpha_u=0; nu_phys=0; nu_coeff=0;
}
static void set_method_sph(void){
    w_coeff=0; use_riemann=0; use_muscl=0; av_mode=0;
    av_alpha=1; av_beta=2; face_dp=0; face_dv=0;
    shock_th=0; prev_dt=0;
    use_AA=0; alpha_u=0; nu_phys=0; nu_coeff=0;
}
/* M(n,m) + NS stress + dp/Z + Monaghan AV — matches multi-D recipe */
static void set_method_ns_av(void){
    w_coeff=0; use_riemann=0; use_muscl=0; av_mode=0;
    av_alpha=0.5; av_beta=4.0; face_dp=0.25; face_dv=0;
    nu_coeff=0.10; nu_phys=0;
    shock_th=0; prev_dt=0;
    use_AA=0; alpha_u=0;
}

/* ============================================================
 *  main
 * ============================================================ */
int main(int argc, char **argv)
{
    /* ---- Setup Sod problem for Phase 1 sweeps ---- */
    setup_problem(PROB_SOD);
    exact_riemann(&ps_exact,&vs_exact);

    int np_sweep = 200;
    ic_type = 1; use_AA = 0; alpha_u = 0; use_muscl = 0;

    /* ---- Phase 1a: Full Acoustic Riemann sweep ---- */
    fprintf(stderr,"\n===== Full Acoustic Riemann (N=%d) =====\n", np_sweep);
    fprintf(stderr,"%5s %4s %4s | %10s %10s %10s\n",
            "w","a_av","b_av","L1_rho","p_osc","eng");
    fprintf(stderr,"------|----|----|------------|------------|----------\n");

    double best_score=1e30;
    double best_wc=0.2, best_al=0, best_bt=0;

    double wc_v[] = {0.0, 0.1, 0.2, 0.3, 0.4};
    double al_v[] = {0.0, 0.5, 1.0, 2.0};
    double bt_v[] = {0.0, 0.5, 1.0, 2.0};

    for(int iw=0; iw<5; iw++)
    for(int ia=0; ia<4; ia++){
        double al = al_v[ia], bt = bt_v[ia];
        w_coeff = wc_v[iw];
        use_AA = (w_coeff > 0) ? 1 : 0;
        use_riemann = 1; shock_th = 0; av_mode = 0;
        av_alpha = al; av_beta = bt;
        face_dp = 0; face_dv = 0; prev_dt = 0;
        Metrics m = run_sim(LAGUERRE, np_sweep, "/dev/null", 0);
        double sc = m.l1_rho + m.p_osc;
        fprintf(stderr," %4.1f %4.1f %4.1f | %10.4e %10.4e %10.4e%s\n",
                wc_v[iw], al, bt, m.l1_rho, m.p_osc, m.eng_err,
                (sc<best_score)?" *":"");
        if(sc < best_score){
            best_score = sc;
            best_wc = wc_v[iw]; best_al = al; best_bt = bt;
        }
    }

    /* Also test Voronoi with full Riemann */
    fprintf(stderr,"\n--- Voronoi with full Riemann ---\n");
    for(int ia=0; ia<4; ia++){
        w_coeff = 0;
        use_riemann = 1; shock_th = 0; av_mode = 0;
        av_alpha = al_v[ia]; av_beta = bt_v[ia];
        face_dp = 0; face_dv = 0; prev_dt = 0;
        Metrics m = run_sim(VORONOI, np_sweep, "/dev/null", 0);
        fprintf(stderr,"  Vor a=%3.1f b=%3.1f | %10.4e %10.4e %10.4e\n",
                al_v[ia], bt_v[ia], m.l1_rho, m.p_osc, m.eng_err);
    }

    fprintf(stderr,"\nBest acoustic: w=%.1f av(%.1f,%.1f) score=%.4e\n",
            best_wc, best_al, best_bt, best_score);

    /* ---- Phase 1b: Exact & HLLC Riemann sweep ---- */
    double wc_v2[] = {0.0, 0.1, 0.2, 0.3, 0.4};
    double al_v2[] = {0.0, 0.5, 1.0};

    double best_exact_score=1e30, best_exact_wc=0.2, best_exact_al=0, best_exact_bt=0;
    double best_hllc_score=1e30,  best_hllc_wc=0.2,  best_hllc_al=0,  best_hllc_bt=0;

    for(int ri = 2; ri <= 3; ri++){
        const char *rlabel = (ri==2) ? "Exact" : "HLLC";
        double *b_score = (ri==2) ? &best_exact_score : &best_hllc_score;
        double *b_wc    = (ri==2) ? &best_exact_wc    : &best_hllc_wc;
        double *b_al    = (ri==2) ? &best_exact_al    : &best_hllc_al;
        double *b_bt    = (ri==2) ? &best_exact_bt    : &best_hllc_bt;

        fprintf(stderr,"\n===== %s Riemann (N=%d) =====\n", rlabel, np_sweep);
        fprintf(stderr,"%5s %4s %4s | %10s %10s %10s\n",
                "w","a_av","b_av","L1_rho","p_osc","eng");
        fprintf(stderr,"------|----|----|------------|------------|----------\n");

        for(int iw=0; iw<5; iw++)
        for(int ia=0; ia<3; ia++){
            double al = al_v2[ia], bt = al_v2[ia];
            w_coeff = wc_v2[iw];
            use_AA = (w_coeff > 0) ? 1 : 0;
            use_riemann = ri; shock_th = 0; av_mode = 0; use_muscl = 0;
            av_alpha = al; av_beta = bt;
            face_dp = 0; face_dv = 0; prev_dt = 0;
            Metrics m = run_sim(LAGUERRE, np_sweep, "/dev/null", 0);
            double sc = m.l1_rho + m.p_osc;
            fprintf(stderr," %4.1f %4.1f %4.1f | %10.4e %10.4e %10.4e%s\n",
                    wc_v2[iw], al, bt, m.l1_rho, m.p_osc, m.eng_err,
                    (sc < *b_score)?" *":"");
            if(sc < *b_score){
                *b_score = sc;
                *b_wc = wc_v2[iw]; *b_al = al; *b_bt = bt;
            }
        }
        fprintf(stderr,"\nBest %s: w=%.1f av(%.1f,%.1f) score=%.4e\n",
                rlabel, *b_wc, *b_al, *b_bt, *b_score);
    }

    /* ---- Phase 1c: Price signal AV sweep ---- */
    double best_price_score=1e30, best_price_wc=0, best_price_al=1, best_price_bt=2;
    double best_price_l1=1e30, best_price_posc=1e30;
    {
        double pw_v[]  = {0.0, 0.2, 0.4};
        double pal_v[] = {0.5, 1.0, 1.5, 2.0};
        double pbt_v[] = {1.0, 2.0, 3.0};

        fprintf(stderr,"\n===== Price Signal AV (no Riemann, N=%d) =====\n", np_sweep);
        fprintf(stderr,"%5s %5s %5s | %10s %10s %10s\n",
                "w","alpha","beta","L1_rho","p_osc","eng");
        fprintf(stderr,"------|-----|-----|------------|------------|----------\n");

        for(int iw=0; iw<3; iw++)
        for(int ia=0; ia<4; ia++)
        for(int ib=0; ib<3; ib++){
            w_coeff = pw_v[iw];
            use_AA = (w_coeff > 0) ? 1 : 0;
            av_alpha = pal_v[ia]; av_beta = pbt_v[ib];
            use_riemann = 0; use_muscl = 1; av_mode = 1;
            shock_th = 0; face_dp = 0; face_dv = 0; prev_dt = 0;
            Metrics m = run_sim(LAGUERRE, np_sweep, "/dev/null", 0);
            double sc = m.l1_rho + m.p_osc;
            fprintf(stderr," %4.1f %5.1f %5.1f | %10.4e %10.4e %10.4e%s\n",
                    pw_v[iw], pal_v[ia], pbt_v[ib],
                    m.l1_rho, m.p_osc, m.eng_err,
                    (sc < best_price_score)?" *":"");
            if(sc < best_price_score){
                best_price_score = sc;
                best_price_wc = pw_v[iw];
                best_price_al = pal_v[ia];
                best_price_bt = pbt_v[ib];
                best_price_l1 = m.l1_rho;
                best_price_posc = m.p_osc;
            }
        }
        fprintf(stderr,"\nBest Price AV: w=%.1f alpha=%.1f beta=%.1f "
                "L1=%.4e p_osc=%.4e score=%.4e\n",
                best_price_wc, best_price_al, best_price_bt,
                best_price_l1, best_price_posc, best_price_score);
    }

    /* ---- Phase 1d: Price signal AV + HLLC sweep ---- */
    double best_ph_score=1e30, best_ph_wc=0, best_ph_al=0.5, best_ph_bt=1;
    double best_ph_l1=1e30, best_ph_posc=1e30;
    {
        double phw_v[]  = {0.0, 0.2};
        double phal_v[] = {0.0, 0.25, 0.5, 1.0};
        double phbt_v[] = {1.0, 2.0};

        fprintf(stderr,"\n===== Price Signal AV + HLLC (N=%d) =====\n", np_sweep);
        fprintf(stderr,"%5s %5s %5s | %10s %10s %10s\n",
                "w","alpha","beta","L1_rho","p_osc","eng");
        fprintf(stderr,"------|-----|-----|------------|------------|----------\n");

        for(int iw=0; iw<2; iw++)
        for(int ia=0; ia<4; ia++)
        for(int ib=0; ib<2; ib++){
            w_coeff = phw_v[iw];
            use_AA = (w_coeff > 0) ? 1 : 0;
            av_alpha = phal_v[ia]; av_beta = phbt_v[ib];
            use_riemann = 3; use_muscl = 1; av_mode = 1;
            shock_th = 0; face_dp = 0; face_dv = 0; prev_dt = 0;
            Metrics m = run_sim(LAGUERRE, np_sweep, "/dev/null", 0);
            double sc = m.l1_rho + m.p_osc;
            fprintf(stderr," %4.1f %5.2f %5.1f | %10.4e %10.4e %10.4e%s\n",
                    phw_v[iw], phal_v[ia], phbt_v[ib],
                    m.l1_rho, m.p_osc, m.eng_err,
                    (sc < best_ph_score)?" *":"");
            if(sc < best_ph_score){
                best_ph_score = sc;
                best_ph_wc = phw_v[iw];
                best_ph_al = phal_v[ia];
                best_ph_bt = phbt_v[ib];
                best_ph_l1 = m.l1_rho;
                best_ph_posc = m.p_osc;
            }
        }
        fprintf(stderr,"\nBest HLLC+Price: w=%.1f alpha=%.2f beta=%.1f "
                "L1=%.4e p_osc=%.4e score=%.4e\n",
                best_ph_wc, best_ph_al, best_ph_bt,
                best_ph_l1, best_ph_posc, best_ph_score);
    }

    /* ---- Phase 1e: CD10 sweep (no Riemann) ---- */
    double best_cd_score=1e30, best_cd_wc=0, best_cd_bt=1, best_cd_amax=2, best_cd_ell=0.1;
    double best_cd_l1=1e30, best_cd_posc=1e30;
    {
        double cdw_v[]   = {0.0, 0.2};
        double cdax_v[]  = {1.0, 2.0};
        double cdel_v[]  = {0.05, 0.1, 0.2};
        double cdbt_v[]  = {1.0, 2.0};

        fprintf(stderr,"\n===== CD10 (no Riemann, N=%d) =====\n", np_sweep);
        fprintf(stderr,"%5s %5s %5s %5s | %10s %10s %10s\n",
                "w","amax","ell","beta","L1_rho","p_osc","eng");
        fprintf(stderr,"------|-----|-----|-----|------------|------------|----------\n");

        for(int iw=0; iw<2; iw++)
        for(int iax=0; iax<2; iax++)
        for(int iel=0; iel<3; iel++)
        for(int ib=0; ib<2; ib++){
            w_coeff = cdw_v[iw];
            use_AA = (w_coeff > 0) ? 1 : 0;
            cd_amax = cdax_v[iax]; cd_ell = cdel_v[iel];
            av_beta = cdbt_v[ib]; cd_amin = 0;
            use_riemann = 0; use_muscl = 1; av_mode = 2;
            shock_th = 0; face_dp = 0; face_dv = 0; prev_dt = 0;
            Metrics m = run_sim(LAGUERRE, np_sweep, "/dev/null", 0);
            double sc = m.l1_rho + m.p_osc;
            fprintf(stderr," %4.1f %5.1f %5.2f %5.1f | %10.4e %10.4e %10.4e%s\n",
                    cdw_v[iw], cdax_v[iax], cdel_v[iel], cdbt_v[ib],
                    m.l1_rho, m.p_osc, m.eng_err,
                    (sc < best_cd_score)?" *":"");
            if(sc < best_cd_score){
                best_cd_score = sc;
                best_cd_wc = cdw_v[iw]; best_cd_amax = cdax_v[iax];
                best_cd_ell = cdel_v[iel]; best_cd_bt = cdbt_v[ib];
                best_cd_l1 = m.l1_rho; best_cd_posc = m.p_osc;
            }
        }
        fprintf(stderr,"\nBest CD10: w=%.1f amax=%.1f ell=%.2f beta=%.1f "
                "L1=%.4e p_osc=%.4e score=%.4e\n",
                best_cd_wc, best_cd_amax, best_cd_ell, best_cd_bt,
                best_cd_l1, best_cd_posc, best_cd_score);
    }

    /* ---- Phase 1f: CD10 + HLLC sweep ---- */
    double best_cdh_score=1e30, best_cdh_amax=1, best_cdh_ell=0.1, best_cdh_bt=1;
    double best_cdh_l1=1e30, best_cdh_posc=1e30;
    {
        double cdhax_v[] = {0.5, 1.0, 2.0};
        double cdhel_v[] = {0.05, 0.1, 0.2};
        double cdhbt_v[] = {1.0, 2.0};

        fprintf(stderr,"\n===== CD10 + HLLC (N=%d) =====\n", np_sweep);
        fprintf(stderr,"%5s %5s %5s | %10s %10s %10s\n",
                "amax","ell","beta","L1_rho","p_osc","eng");
        fprintf(stderr,"------|-----|-----|------------|------------|----------\n");

        for(int iax=0; iax<3; iax++)
        for(int iel=0; iel<3; iel++)
        for(int ib=0; ib<2; ib++){
            w_coeff = 0.0;
            cd_amax = cdhax_v[iax]; cd_ell = cdhel_v[iel];
            av_beta = cdhbt_v[ib]; cd_amin = 0;
            use_riemann = 3; use_muscl = 1; av_mode = 2;
            shock_th = 0; face_dp = 0; face_dv = 0; prev_dt = 0;
            Metrics m = run_sim(LAGUERRE, np_sweep, "/dev/null", 0);
            double sc = m.l1_rho + m.p_osc;
            fprintf(stderr," %5.1f %5.2f %5.1f | %10.4e %10.4e %10.4e%s\n",
                    cdhax_v[iax], cdhel_v[iel], cdhbt_v[ib],
                    m.l1_rho, m.p_osc, m.eng_err,
                    (sc < best_cdh_score)?" *":"");
            if(sc < best_cdh_score){
                best_cdh_score = sc;
                best_cdh_amax = cdhax_v[iax]; best_cdh_ell = cdhel_v[iel];
                best_cdh_bt = cdhbt_v[ib];
                best_cdh_l1 = m.l1_rho; best_cdh_posc = m.p_osc;
            }
        }
        fprintf(stderr,"\nBest HLLC+CD10: amax=%.1f ell=%.2f beta=%.1f "
                "L1=%.4e p_osc=%.4e score=%.4e\n",
                best_cdh_amax, best_cdh_ell, best_cdh_bt,
                best_cdh_l1, best_cdh_posc, best_cdh_score);
    }

    /* reset av_mode for Monaghan-based runs */
    av_mode = 0; use_muscl = 0; nu_phys = 0; nu_coeff = 0;

    /* ---- Summary ---- */
    fprintf(stderr,"\n===== Comparison =====\n");
    fprintf(stderr,"  Acoustic best:  w=%.1f av(%.1f,%.1f) score=%.4e\n",
            best_wc, best_al, best_bt, best_score);
    fprintf(stderr,"  Exact best:     w=%.1f av(%.1f,%.1f) score=%.4e\n",
            best_exact_wc, best_exact_al, best_exact_bt, best_exact_score);
    fprintf(stderr,"  HLLC best:      w=%.1f av(%.1f,%.1f) score=%.4e\n",
            best_hllc_wc, best_hllc_al, best_hllc_bt, best_hllc_score);
    fprintf(stderr,"  Price AV best:  w=%.1f a=%.1f b=%.1f score=%.4e\n",
            best_price_wc, best_price_al, best_price_bt, best_price_score);
    fprintf(stderr,"  HLLC+Price:     w=%.1f a=%.2f b=%.1f score=%.4e\n",
            best_ph_wc, best_ph_al, best_ph_bt, best_ph_score);
    fprintf(stderr,"  CD10 best:      w=%.1f amax=%.1f ell=%.2f b=%.1f score=%.4e\n",
            best_cd_wc, best_cd_amax, best_cd_ell, best_cd_bt, best_cd_score);
    fprintf(stderr,"  HLLC+CD10:      amax=%.1f ell=%.2f b=%.1f score=%.4e\n",
            best_cdh_amax, best_cdh_ell, best_cdh_bt, best_cdh_score);

    /* ============================================================
     *  Phase W: Laguerre Weight Experiments
     *  weight_mode × w_coeff × face_shift_max × 3 methods × 4 problems
     * ============================================================ */
    {
        const char *wm_name[] = {"PRESS", "ENTH", "SPECV", "CSND"};
        double wc_exp[]  = {0.0, 0.1, 0.2, 0.3, 0.5};
        double fs_exp[]  = {0.2, 0.4, 1.0};  /* face shift limit */
        int n_wc = 5, n_fs = 3, n_wm = 4;

        /* 3 Laguerre methods to test */
        typedef struct { const char *name; void (*setup)(void); } WMethod;
        WMethod wmeth[] = {
            {"HLLC+CD10",  set_method_hllc_cd10},
            {"HLLC+Price", set_method_hllc_price},
            {"Price AV",   set_method_price}
        };
        int n_wmeth = 3;

        fprintf(stderr,"\n===== Phase W: Laguerre Weight Experiments (N=200) =====\n");
        fprintf(stderr,"%-10s %-5s %5s %5s | %-10s %-10s %-10s %-10s\n",
                "Method","WMode","w_c","fs_lim","Sod","Blast","ShuOsher","Noh");
        fprintf(stderr,"----------|-----|-----|-----|----------|----------|----------|----------\n");

        /* w=0 baselines first */
        for(int im=0; im<n_wmeth; im++){
            double l1_base[4];
            for(int pid=0; pid<4; pid++){
                setup_problem(pid);
                if(!has_exact_sol){
                    char rf[128];
                    sprintf(rf, "%s/ref.dat", (const char*[]){"sod","blast","shuosher","noh"}[pid]);
                    read_reference(rf);
                } else if(pid==PROB_SOD){
                    exact_riemann(&ps_exact,&vs_exact);
                }
                wmeth[im].setup();
                w_coeff = 0;
                use_AA = 0;
                weight_mode = WM_PRESSURE;
                face_shift_max = 1.0;
                Metrics m = run_sim(LAGUERRE, 200, "/dev/null", 0);
                l1_base[pid] = m.failed ? -1 : m.l1_rho;
            }
            fprintf(stderr,"%-10s %-5s %5.1f %5.1f | %10.4e %10.4e %10.4e %10.4e  <<BASE>>\n",
                    wmeth[im].name, "w=0", 0.0, 1.0,
                    l1_base[0], l1_base[1], l1_base[2], l1_base[3]);
        }

        fprintf(stderr,"----------|-----|-----|-----|----------|----------|----------|----------\n");

        /* sweep all combinations */
        for(int im=0; im<n_wmeth; im++)
        for(int iwm=0; iwm<n_wm; iwm++)
        for(int iwc=1; iwc<n_wc; iwc++)  /* skip w=0 (already done) */
        for(int ifs=0; ifs<n_fs; ifs++){
            double l1[4];
            int any_failed = 0;
            for(int pid=0; pid<4; pid++){
                setup_problem(pid);
                if(!has_exact_sol){
                    char rf[128];
                    sprintf(rf, "%s/ref.dat", (const char*[]){"sod","blast","shuosher","noh"}[pid]);
                    read_reference(rf);
                } else if(pid==PROB_SOD){
                    exact_riemann(&ps_exact,&vs_exact);
                }
                wmeth[im].setup();
                w_coeff = wc_exp[iwc];
                use_AA = (w_coeff > 0) ? 1 : 0;
                weight_mode = iwm;
                face_shift_max = fs_exp[ifs];
                Metrics m = run_sim(LAGUERRE, 200, "/dev/null", 0);
                l1[pid] = m.failed ? -1 : m.l1_rho;
                if(m.failed) any_failed = 1;
            }
            fprintf(stderr,"%-10s %-5s %5.1f %5.1f | %10.4e %10.4e %10.4e %10.4e\n",
                    wmeth[im].name, wm_name[iwm],
                    wc_exp[iwc], fs_exp[ifs],
                    l1[0], l1[1], l1[2], l1[3]);
        }
        /* ---- Multi-problem: AA + gentle stabilization ---- */
        fprintf(stderr,"\n--- Multi-Problem: AA + gentle stabilization (HLLC+CD10) ---\n");
        {
            const char *pname2[] = {"sod","blast","shuosher","noh"};

            struct { double wc, fs, sd, pr, gm, rm; int mode; const char *lbl; const char *fname; } cfgs[] = {
                {0.0,  1.0,  0, 0, 0, 0, LAGUERRE, "Base (w=0)",      "base"},
                {0.10, 0.2,  0, 0, 0, 0, LAGUERRE, "AA only w=0.10",  "aa_only"},
                {0.08, 0.2, 15, 0, 0, 0, LAGUERRE, "AA+Ks15 w=0.08",  "aa_ks15"},
                {0.0,  1.0,  0, 0, 0, 0, SPH_MODE, "SPH",             "sph"},
            };
            int n_cfgs = sizeof(cfgs)/sizeof(cfgs[0]);

            fprintf(stderr,"%-16s %5s %4s %5s %5s %5s | %-10s %-10s %-10s %-10s\n",
                    "Label","w","fs","s_dmp","gmax","rmax","Sod","Blast","ShuOsher","Noh");
            fprintf(stderr,"----------------|-----|----|----- |----- |-----|----------|----------|----------|----------\n");

            for(int ic=0; ic<n_cfgs; ic++){
                double l1[4];
                for(int pid=0; pid<4; pid++){
                    setup_problem(pid);
                    if(!has_exact_sol){
                        char rf[128];
                        sprintf(rf, "%s/ref.dat", pname2[pid]);
                        read_reference(rf);
                    } else if(pid==PROB_SOD){
                        exact_riemann(&ps_exact,&vs_exact);
                    }
                    if(cfgs[ic].mode == SPH_MODE){
                        set_method_sph();
                    } else {
                        set_method_hllc_cd10();
                        w_coeff = cfgs[ic].wc;
                        weight_mode = WM_SPECVOL;
                        face_shift_max = cfgs[ic].fs;
                        use_AA = (w_coeff > 0) ? 1 : 0;
                        aa_max_frac = 1.0;
                        shock_w_damp = cfgs[ic].sd;
                        shock_p_ratio = cfgs[ic].pr;
                        w_grad_max = cfgs[ic].gm;
                        w2_rate_max = cfgs[ic].rm;
                        w_relax_tau = 0; use_adiab_w = 0;
                        p_exp = 0; w_floor_frac = 0;
                    }
                    char outf[256];
                    sprintf(outf, "%s/%s.dat", pname2[pid], cfgs[ic].fname);
                    mkdir(pname2[pid], 0755);
                    Metrics m = run_sim(cfgs[ic].mode, 200, outf, 0);
                    l1[pid] = m.failed ? -1.0 : m.l1_rho;
                }
                fprintf(stderr,"%-16s %5.2f %4.1f %5.1f %5.1f %5.1f | %10.4e %10.4e %10.4e %10.4e\n",
                        cfgs[ic].lbl, cfgs[ic].wc, cfgs[ic].fs,
                        cfgs[ic].sd, cfgs[ic].gm, cfgs[ic].rm,
                        l1[0], l1[1], l1[2], l1[3]);
            }

            shock_w_damp = 0; shock_p_ratio = 0; w_relax_tau = 0; w_grad_max = 0;
            use_adiab_w = 0; aa_max_frac = 1.0; use_AA = 0;
            w2_rate_max = 0; w_floor_frac = 0;
        }
        fprintf(stderr,"--- End Multi-Problem Test ---\n\n");

        fprintf(stderr,"===== End Phase W =====\n\n");
    }

    /* ============================================================
     *  Multi-problem test suite (5 methods × 8 problems)
     * ============================================================ */
    const char *pname[]  = {"sod", "blast", "shuosher", "noh",
                            "lax", "dblfan", "collision", "contact"};
    const char *ptitle[] = {"Sod Shock Tube", "Interacting Blast Waves",
                            "Shu-Osher Shock-Entropy", "Noh Strong Shock",
                            "Lax Shock Tube", "Double Rarefaction (Toro #2)",
                            "Two-Shock Collision (Toro #4)", "Stationary Contact"};
    const char *mname[]  = {"hllc_cd10", "hllc_price", "price", "voronoi", "sph", "ns_av"};
    const char *mlabel[] = {"HLLC+CD10", "HLLC+Price", "Price AV",
                            "Vor M(3,4)+dpZ", "SPH", "NS+AV"};
    const char *mcolor[] = {"#a65628", "#377eb8", "#e41a1c", "#4daf4a", "#984ea3", "#ff7f00"};
    int         mmode[]  = {LAGUERRE, LAGUERRE, LAGUERRE, VORONOI, SPH_MODE, LAGUERRE};

    typedef void (*msetup_fn)(void);
    msetup_fn msetup[] = {set_method_hllc_cd10, set_method_hllc_price,
                          set_method_price, set_method_voronoi, set_method_sph,
                          set_method_ns_av};

    int np_test = 200;
    int nprob = 8, nmeth = 6;

    for(int pid = 0; pid < nprob; pid++){
        setup_problem(pid);

        /* create subdirectory */
        mkdir(pname[pid], 0755);
        fprintf(stderr,"\n====== %s (N=%d) -> %s/ ======\n",
                ptitle[pid], np_test, pname[pid]);

        /* exact solution or high-res reference */
        char ref_file[128];
        sprintf(ref_file, "%s/ref.dat", pname[pid]);

        if(has_exact_sol){
            if(pid != PROB_NOH)
                exact_riemann(&ps_exact, &vs_exact);
        } else {
            /* generate N=4000 reference with HLLC+CD10 */
            set_method_hllc_cd10();
            fprintf(stderr,"  Generating N=4000 reference...\n");
            run_sim(LAGUERRE, 4000, ref_file, 1);
            read_reference(ref_file);
        }

        /* reset y-range */
        for(int k=0;k<4;k++){ gmin[k]=1e30; gmax[k]=-1e30; }

        /* run 5 methods, track which succeed */
        char files[6][128];
        int mfailed[6] = {0};
        for(int mi=0; mi<nmeth; mi++){
            msetup[mi]();
            sprintf(files[mi], "%s/%s.dat", pname[pid], mname[mi]);
            Metrics m = run_sim(mmode[mi], np_test, files[mi], 1);
            mfailed[mi] = m.failed;
        }

        /* write exact/reference for plotting */
        if(has_exact_sol)
            write_exact_file(ref_file);
        /* else: already written by N=4000 run */

        /* update y-range from reference file */
        double evals[NMAX]; int ne=0;
        {
            FILE *rf = fopen(ref_file, "r");
            char buf[256];
            while(rf && fgets(buf, sizeof(buf), rf)){
                if(buf[0]=='#') continue;
                double x, rho, vel, pre, eng;
                if(sscanf(buf, "%lf %lf %lf %lf %lf", &x, &rho, &vel, &pre, &eng)>=5){
                    double vv[]={rho,vel,pre,eng};
                    for(int k=0;k<4;k++){
                        if(vv[k]<gmin[k]) gmin[k]=vv[k];
                        if(vv[k]>gmax[k]) gmax[k]=vv[k];
                    }
                    if(ne<NMAX) evals[ne++]=eng;
                }
            }
            if(rf) fclose(rf);
        }
        /* For internal energy: detect bimodal gap (e.g. Blast wall artifacts)
         * Clip only if: gap > 30% of range AND < 15% of values above gap */
        if(ne>10){
            qsort(evals, ne, sizeof(double), cmp_dbl);
            double erange = evals[ne-1]-evals[0];
            if(erange > 1e-10){
                double max_gap=0; int gap_idx=-1;
                for(int i=1;i<ne;i++){
                    double g=evals[i]-evals[i-1];
                    if(g>max_gap){ max_gap=g; gap_idx=i; }
                }
                double frac_above = (double)(ne-gap_idx)/ne;
                if(max_gap > 0.3*erange && frac_above < 0.15){
                    double below=evals[gap_idx-1];
                    double ipad=0.12*(below-evals[0]);
                    if(ipad<0.1) ipad=0.1;
                    gmin[3]=evals[0];
                    gmax[3]=below+ipad;
                }
            }
        }

        /* expand y-range with method data, capped at 50% beyond ref range */
        double gmin_base[4], gmax_base[4];
        for(int k=0;k<4;k++){ gmin_base[k]=gmin[k]; gmax_base[k]=gmax[k]; }
        for(int mi=0; mi<nmeth; mi++){
            if(mfailed[mi]) continue;
            FILE *mf = fopen(files[mi], "r");
            if(!mf) continue;
            char buf[256];
            while(fgets(buf, sizeof(buf), mf)){
                if(buf[0]=='#') continue;
                double x, rho, vel, pre, eng;
                if(sscanf(buf, "%lf %lf %lf %lf %lf", &x,&rho,&vel,&pre,&eng)>=5){
                    double vv[]={rho,vel,pre,eng};
                    for(int k=0;k<4;k++){
                        if(vv[k]<gmin[k]) gmin[k]=vv[k];
                        if(vv[k]>gmax[k]) gmax[k]=vv[k];
                    }
                }
            }
            fclose(mf);
        }
        /* cap expansion: don't let method outliers expand more than 50% of base range */
        for(int k=0;k<4;k++){
            double base_range = gmax_base[k]-gmin_base[k];
            if(base_range<1e-10) base_range=0.1;
            double lo_limit = gmin_base[k] - 0.5*base_range;
            double hi_limit = gmax_base[k] + 0.5*base_range;
            if(gmin[k] < lo_limit) gmin[k] = lo_limit;
            if(gmax[k] > hi_limit) gmax[k] = hi_limit;
        }

        double ylo[4], yhi[4];
        for(int k=0;k<4;k++){
            double pad=0.08*(gmax[k]-gmin[k]);
            if(pad<1e-10) pad=0.1;
            ylo[k]=gmin[k]-pad;
            yhi[k]=gmax[k]+pad;
        }
        /* clamp non-negative quantities */
        if(ylo[0]<0) ylo[0]=0; /* density */
        if(ylo[2]<0) ylo[2]=0; /* pressure */
        if(ylo[3]<0) ylo[3]=0; /* internal energy */

        /* generate gnuplot script (exclude failed methods) */
        char gpf[128], pngf[128], title[256];
        sprintf(gpf, "%s/plot_%d.gp", pname[pid], np_test);
        sprintf(pngf, "%s/result_%d.png", pname[pid], np_test);
        sprintf(title, "%s (N=%d)", ptitle[pid], np_test);

        const char *fptrs[6], *lptrs[6], *cptrs[6];
        int nplot = 0;
        for(int mi=0; mi<nmeth; mi++){
            if(!mfailed[mi]){
                fptrs[nplot] = files[mi];
                lptrs[nplot] = mlabel[mi];
                cptrs[nplot] = mcolor[mi];
                nplot++;
            }
        }

        write_plot_grid(gpf, pngf, ref_file, title,
                        nplot, lptrs, fptrs, cptrs,
                        x0_dom, x1_dom, ylo, yhi);

        fprintf(stderr,"Plot: gnuplot %s -> %s\n", gpf, pngf);
    }

    /* ============================================================
     *  Phase NS: NS stress + dp/Z + Monaghan AV sweep
     *  In 1D, simple average has no upwinding (unlike multi-D M(n,m)
     *  which has transverse terms). dp/Z correction provides the
     *  1D equivalent of multi-D M(n,m) implicit upwinding.
     * ============================================================ */
    {
        double nc_v[]  = {0.0, 0.02, 0.05, 0.1};
        double dp_v[]  = {0.0, 0.25, 0.5};
        double al_ns[] = {0.5, 1.0, 2.0, 3.0};
        double bt_ns[] = {0.0, 2.0, 4.0};
        int n_nc = 4, n_dp = 3, n_al = 4, n_bt = 3;
        const char *pname_ns[] = {"sod","blast","shuosher","noh",
                                  "lax","dblfan","collision","contact"};

        fprintf(stderr,"\n===== Phase NS: NS + dp/Z + Monaghan AV (N=200) =====\n");
        fprintf(stderr,"%6s %4s %4s %4s | %10s %10s %10s %10s %10s %10s %10s %10s\n",
                "nu_c","dpZ","a_av","b_av","Sod","Blast","ShuOsher","Noh",
                "Lax","DblFan","Collision","Contact");
        fprintf(stderr,"------|----|----|----|-");
        for(int k=0;k<8;k++) fprintf(stderr,"----------|-");
        fprintf(stderr,"\n");

        double best_ns_score = 1e30;
        double best_nc = 0, best_dpz = 0.5, best_aln = 3, best_btn = 4;

        for(int inc=0; inc<n_nc; inc++)
        for(int idp=0; idp<n_dp; idp++)
        for(int ia=0; ia<n_al; ia++)
        for(int ib=0; ib<n_bt; ib++){
            double l1[8];
            int any_fail = 0;
            for(int pid=0; pid<8; pid++){
                setup_problem(pid);
                if(!has_exact_sol){
                    char rf[128];
                    sprintf(rf, "%s/ref.dat", pname_ns[pid]);
                    read_reference(rf);
                } else if(pid != PROB_NOH){
                    exact_riemann(&ps_exact,&vs_exact);
                }
                w_coeff=0; use_riemann=0; use_muscl=0; av_mode=0;
                face_dv=0; shock_th=0; prev_dt=0;
                use_AA=0; alpha_u=0;
                nu_coeff = nc_v[inc]; nu_phys = 0;
                face_dp = dp_v[idp];
                av_alpha = al_ns[ia]; av_beta = bt_ns[ib];
                Metrics m = run_sim(LAGUERRE, 200, "/dev/null", 0);
                l1[pid] = m.failed ? -1.0 : m.l1_rho;
                if(m.failed) any_fail = 1;
            }
            double score = 0;
            for(int k=0;k<8;k++){
                if(l1[k] < 0) { score = 1e30; break; }
                score += l1[k];
            }
            fprintf(stderr," %5.3f %4.2f %4.1f %4.1f |",
                    nc_v[inc],dp_v[idp],al_ns[ia],bt_ns[ib]);
            for(int k=0;k<8;k++) fprintf(stderr," %10.4e",l1[k]);
            fprintf(stderr,"%s\n",(score<best_ns_score)?" *":"");
            if(score < best_ns_score){
                best_ns_score = score;
                best_nc = nc_v[inc]; best_dpz = dp_v[idp];
                best_aln = al_ns[ia]; best_btn = bt_ns[ib];
            }
        }
        fprintf(stderr,"\nBest NS+AV: nu_coeff=%.3f dp/Z=%.2f alpha=%.1f beta=%.1f score=%.4e\n",
                best_nc, best_dpz, best_aln, best_btn, best_ns_score);
        nu_phys = 0; nu_coeff = 0;
    }

    return 0;
}
