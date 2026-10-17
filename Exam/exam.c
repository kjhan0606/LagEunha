#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>

#ifdef _OPENMP
#include<omp.h>
#endif

#include "eunha.h"

#include "voro.h"
#include "voro_eunha.h"
#include "nnost.h"
#include "gnnost.h"
#include "exam.h"
#include "exam2d.h"
#include "gfs_pair.h"

/* Dual energy for the RK4 blend path. GFS_DUAL_ENERGY = eta > 0 turns it
 * on (default off). The entropy K = P/rho^gamma is stress.K. It is set from
 * P and rho on the first RK4 call when it is not positive, and SEDOV_PHASE1
 * keeps it adiabatic (dK = 0). A cell whose ie falls below eta * ie_K,
 * ie_K = K rho^gamma V/(gamma-1), takes ie_K instead of the P = 1e-6 floor.
 * A shock only raises ie above ie_K, so shocked gas is not touched. */
static double gfs_dual_energy_eta(void){
	static int cached = 0;
	static double eta = 0;
	if(!cached){
		const char *s = getenv("GFS_DUAL_ENERGY");
		eta = (s && s[0]) ? atof(s) : 0.0;
		if(!(eta > 0)) eta = 0;
		cached = 1;
	}
	return eta;
}

/* GIZMO-style entropy switch and half-limit for the RK4 blend path.
 * Both off by default. They act once per step, on the host, after the
 * final RK4 combination and the new Voronoi volumes, where GIZMO's kick
 * applies them (kicks.c, ENERGY_ENTROPY_SWITCH_IS_ACTIVE and the
 * "dEnt < 0.5 u" line). So they work on GPU runs too.
 *
 * GFS_ENTROPY_SWITCH=1   gravity-dominated cold cells, where
 *                        coef * m |g| h > max(0.5 ie_n, ie_E),
 *                        take the adiabatic ie = ie_n (V_n/V)^(gamma-1)
 *                        (K = P/rho^gamma held over the step) instead of
 *                        the energy-equation value ie_E. g is gravity only
 *                        (kepler_accel + GAS_ACC), h = sqrt(V).
 * GFS_ES_COEF            coef, default 0.01 (GIZMO).
 * GFS_HALF_LIMIT         0/1. ie < 0.5 ie_n becomes 0.5 ie_n. No floor
 *                        refill. Default: same as GFS_ENTROPY_SWITCH.
 *
 * ie_n and V_n are held in stress.E_inv_xx / E_inv_xy, which only the
 * LagMFM path (av_mode 4) uses, so they migrate with the particle. The
 * switch is disabled for av_mode 4 and entropy_mode 1. */
#define GFS_ES_IE0(p) ((p).stress.E_inv_xx)
#define GFS_ES_V0(p)  ((p).stress.E_inv_xy)
static int gfs_es_on = -1, gfs_half_on = -1;
static double gfs_es_coef = 0.01;
static void gfs_es_config(void){
	if(gfs_es_on >= 0) return;
	const char *s = getenv("GFS_ENTROPY_SWITCH");
	gfs_es_on = (s && s[0] && atoi(s) > 0) ? 1 : 0;
	s = getenv("GFS_ES_COEF");
	if(s && s[0] && atof(s) > 0) gfs_es_coef = atof(s);
	s = getenv("GFS_HALF_LIMIT");
	gfs_half_on = (s && s[0]) ? (atoi(s) > 0 ? 1 : 0) : gfs_es_on;
}

/* GFS_LAGUERRE_ROTATION=1 extends the face-rotation correction of the face
 * normal speed (voronoi_face_rotation) to Laguerre faces (w2_i != w2_j or
 * either nonzero). Default 0: the correction stays w = 0 only, as before.
 * A Laguerre face is still perpendicular to e = (x_j - x_i)/d, so
 *   u_n(c_f) = u_n(a) - u_ij·(c_f - a)/d,
 * where a = x_i + fact1 (x_j - x_i), fact1 = (1 + (w2_i - w2_j)/d^2)/2, is
 * the anchor point whose velocity get2dUpqradRk4 returns. For w = 0,
 * fact1 = 1/2 and a is the midpoint. The GPU force path reads the same
 * variable (exam_gpu_extract.c, GPUPhysicsParams.laguerre_rot). */
/* GFS_FACE_CHARGE_LIMIT=1 (default 0): apply gfs_pair_charge_limit to the
 * whole face pressure, not only the pair part. A cell that the face sweeps
 * outward (s_i = (u_face - v_i)·dS > 0) then pays at most ie/6 per step on
 * that face. Near the Kepler centre the face-centroid rotation term makes a
 * cold cell expand against a hot neighbour's pressure (P* ~ P_j ~ 100 P_i),
 * which drains it to the stage floor. The cap is symmetric in i and j
 * (same P on both ends), so momentum and energy stay pairwise conserved;
 * it weakens the force between the two cells while it is active. Policy,
 * so off by default. GPU: GPUPhysicsParams.face_charge. */
/* GFS_GEOM_FACE_VEL=1 (default 0): the work and energy terms use the
 * geometric face velocity (midpoint velocity + rotation term, the velocity
 * the Voronoi face actually moves with) on every face. HLL/HLLC still
 * supply p*; their v_n* is no longer used as the face velocity on the
 * phase1 extreme and riemann_vstar faces. GPU: GPUPhysicsParams.geom_fv. */
static int gfs_geom_fv_on(void){
	static int on = -1;
	if(on < 0){ const char *e = getenv("GFS_GEOM_FACE_VEL"); on = (e && atoi(e) > 0) ? 1 : 0; }
	return on;
}
static int gfs_face_charge_on(void){
	static int on = -1;
	if(on < 0){
		const char *s = getenv("GFS_FACE_CHARGE_LIMIT");
		on = (s && s[0] && atoi(s) > 0) ? 1 : 0;
	}
	return on;
}

static int gfs_laguerre_rotation_on(void){
	static int on = -1;
	if(on < 0){
		const char *s = getenv("GFS_LAGUERRE_ROTATION");
		on = (s && s[0] && atoi(s) > 0) ? 1 : 0;
	}
	return on;
}

/* GFS_W2_LINEAR=1 (default 0; CPU RK4 blend path, GAS kappa > 0 only).
 * Books the face motion caused by the weight change in P dV.
 * Default path: W is recomputed from each RK4 stage state, the face speed
 * uses (w_stage - w_n)/dt_old clamped by c_s (about half of the weight
 * motion over a step), and W jumps at the end of the step after the
 * end-of-step volume update. None of the jump enters P dV.
 * With the flag, each particle gets one rate per step,
 *   Wdot = (W_target - W_n)/dt_old,  |Wdot| <= 2 c_s max(w_n, w_target),
 * with W_target from getw2forHydroParticle + applyW2Controls on the state
 * at the start of the step (a one-step relaxation toward the target).
 * The target is capped at w2ceil/2 (margin for the nearest-neighbour
 * ceiling). Stage k uses W = W_n + c_k Dt Wdot (c = 0, 1/2, 1/2, 1), the step ends
 * at W_n + Dt Wdot before the volume update, and the face speed carries
 * exactly (Wdot_i - Wdot_j)/(2d) (w2old = w2 in the stages, so the lagged
 * term in get2dUpqradRk4 is zero). Wdot is held in w2hydro while active.
 * GPU: GPUPhysicsParams.w2lin, and the uploaded w2old array carries Wdot
 * (exam_gpu_extract.c fillParticleSoA, dev_get2dUpqradRk4). */
static int gfs_w2lin_active = 0;
int gfs_w2lin_is_active(void){ return gfs_w2lin_active; }
static int gfs_w2lin_on(void){
	static int on = -1;
	if(on < 0){
		const char *s = getenv("GFS_W2_LINEAR");
		on = (s && s[0] && atoi(s) > 0) ? 1 : 0;
	}
	return on;
}
/* GFS_DYN_W2=1 (default 0; scheme 3, Laguerre + dynamic weights,
 * code_review.md 2(c); CPU RK4 blend path, av_mode 5, kappa = 0, 1 rank).
 * Each RK4 stage solves the weighted least squares
 *   min sum_f a_f (Wdot_p - Wdot_q - b_f)^2,
 *   b_f = 2 d (v*_f - n.ubar - rot_f) + (W_p - W_q) ddot/d,
 * (graph Laplacian L_a Wdot = B^T A b, Jacobi-PCG, zero mean) with v*_f the
 * HLL face speed (cell states along e; GFS_DYN_W2_MID=1: mid-plane speed on
 * ordinary faces). a_f = A_f * GFS_DYN_W2_RATIO on extreme faces
 * (P_max > GFS_DYN_W2_PX * P_min, default 100), A_f elsewhere.
 * W is a state variable: stage W = W_n + c Dt Wdot_prev, final
 * W_n + Dt (k1 + 2k2 + 2k3 + k4)/6, shifted so min W = 0 (only differences
 * enter the geometry). The face speed carries (Wdot_p - Wdot_q)/(2d) and the
 * rotation term, and the work uses this realised speed (as GFS_GEOM_FACE_VEL).
 * Cell validity |W_p - W_q| < (1 - GFS_DYN_W2_DELTA) d^2 (default 0.2):
 * Dt limit from the K1 rates, ((1-delta) d^2 - |dW|)/(d|dW|/dt - (1-delta) d(d^2)/dt),
 * logged as bound=1 / nbind; past (1 - 2 delta) d^2 a face switches to a
 * barrier target (keep dW/d^2, relax over dt_old) with weight
 * a_f * GFS_DYN_W2_BARRIER (default 1000), logged as nbar.
 * Other knobs: GFS_DYN_W2_TOL (CG rel. tol, 1e-10), GFS_DYN_W2_MAXIT (2000),
 * GFS_DYN_W2_OMP=1 (OpenMP face collection).
 * [DYNW2] per step: ||r||/||b|| on extreme and ordinary faces (sums over
 * the 4 stages), rms (realised - v*) on extreme faces from the solve
 * (rms_dv_ext) and from the force loop (rms_dv_ext_force), rms_dv_ord,
 * CG iterations, max |dW|/d^2, invalid faces, Dt limit.
 * Single rank only (the solve is global; not active when NID > 1).
 * No GPU path: CPU only, inactive when the GPU is enabled. */
/* GFS_TORQUE_AUDIT=1 (default 0; diagnostic only, CPU face loop, results
 * unchanged). Torque about (EUNHA_KEPLER_CX, EUNHA_KEPLER_CY) of the face
 * pressure force -p dS on each particle, split into the Riemann part on
 * ordinary faces (HLLC), the Riemann part on extreme faces (HLL average
 * state) and the pair pressure, summed in radius bins of the particle.
 * Faces left without a neighbour (edges of the bounding cell polygon) are
 * counted with their length. Rank 0 prints [TORQ] each step: cumulative
 * integral of the stage-averaged torque over dt, per bin and term. */
#define TQ_NB 9
static const double tq_edges[TQ_NB-1] = {0.5, 0.55, 0.6, 0.75, 1.0, 1.5, 1.9, 2.1};
static double sfl_bin[TQ_NB];
static double gfs_tnow;
/* GFS_FLOOR_LOG=1: one line per stage pressure-floor reset (first 200 at t >= GFS_FLOOR_LOG_T1) */
static void gfs_sfl_log(treevorork4particletype *b, double ie_before){
	static int on = -1, n = 0;
	if(on < 0){ const char *s = getenv("GFS_FLOOR_LOG"); on = (s && s[0] == '1') ? 1 : 0; }
	static double t1 = -1;
	if(t1 < 0){ const char *s = getenv("GFS_FLOOR_LOG_T1"); t1 = (s && s[0]) ? atof(s) : 0; }
	if(!on || n >= 200 || gfs_tnow < t1) return;
	double cx = 0, cy = 0;
	{ const char *a = getenv("EUNHA_KEPLER_CX"), *c = getenv("EUNHA_KEPLER_CY"); if(a) cx = atof(a); if(c) cy = atof(c); }
	double dx = (double)b->x - cx, dy = (double)b->y - cy, R = sqrt(dx*dx + dy*dy);
	double vr = (dx*(double)b->vx + dy*(double)b->vy)/(R > 0 ? R : 1), vp = (dx*(double)b->vy - dy*(double)b->vx)/(R > 0 ? R : 1);
#ifdef _OPENMP
#pragma omp critical (sfl_log)
#endif
	{
		n++;
		fprintf(stderr, "[SFLLOG] t=%.6f id=%ld x=%.5f y=%.5f R=%.5f vr=%.4e vphi=%.4e den=%.4e vol=%.4e mass=%.4e ie_before=%.4e\n",
				gfs_tnow, (long)PINDX(b), (double)b->x, (double)b->y, R, vr, vp, (double)b->den, (double)b->volume, (double)b->mass, ie_before);
	}
}   /* stage pressure-floor resets per bin (GFS_TORQUE_AUDIT) */
/* M1: face torque split by (compression/expansion of the reconstructed normal
 * jump) x (sign of r_i x e_ij): k = 2*compress + (arm > 0) */
static double tqs_acc[TQ_NB][4], tqs_cum[TQ_NB][4], tqs_n[TQ_NB][4], tqs_pp[TQ_NB][4], tqs_du[TQ_NB][4];
/* EUNHA_LIM_POINT=anchor: the limiter checks the reconstruction at the
 * point of the face on the generator line (phi d e, the midpoint for W = 0)
 * instead of the face centroid; a linear field is never limited there
 * (Zier & Springel 2022 §4.3; Hopkins 2015 App. A). Fluxes keep the
 * centroid. Default: centroid. */
static int gfs_lim_anchor_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("EUNHA_LIM_POINT"); on = (s && strcmp(s, "anchor") == 0) ? 1 : 0; }
	return on;
}
static int gfs_gauss2_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("EUNHA_FACE_GAUSS2"); on = (s && atoi(s) > 0) ? 1 : 0; }
	return on;
}
/* EUNHA_LIM_PAIR=1 (with EUNHA_REC_TCOMMON=1): no cell-wide monotonicity
 * factor. Each side's anchor value q_i + g_i.a is clamped per face to the
 * pair interval [min(q_i,q_j), max(q_i,q_j)] (Hopkins 2015 App. B
 * structure), so a free-surface cell whose truncated stencil makes it the
 * extremum no longer loses its whole gradient (the 0/delta degeneracy,
 * Astra 14th reply); a linear field is never clamped. The common
 * tangential increment from the anchor to the centroid is kept, and each
 * side's state at the evaluation point is bounded by the values of the
 * four generators of the two Delaunay triangles sharing the face (i, j,
 * k+, k-): linear interpolation on the face stays inside, rho and P stay
 * positive, and an ill-conditioned LSF gradient (close pair dominating
 * the A/d^2 weights) cannot push both sides to a common unphysical state.
 * Centroid positivity of the cell gradients of rho and P is kept. */
static int gfs_lim_pair_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("EUNHA_LIM_PAIR"); on = (s && atoi(s) > 0) ? 1 : 0; }
	return on;
}
/* EUNHA_REC_TCOMMON=1 (av_mode 5 MUSCL): each side extrapolates to the
 * anchor a = phi e with its own limited gradient and from a to the face
 * centroid with the pair mean gradient. The face jump is then the anchor
 * jump checked by the limiter; the term (g_j - g_i).(c - a), unbounded for
 * close pairs with long faces (256^2 Kepler t=17.08), is removed. Linear
 * fields with exact gradients reconstruct as before. Astra 13th reply. */
static int gfs_rec_tcommon_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("EUNHA_REC_TCOMMON"); on = (s && atoi(s) > 0) ? 1 : 0; }
	return on;
}
/* EUNHA_LSF_ORDER=2: quadratic least-squares fit (gradient + Hessian) over
 * all candidates within r_s = min(2.5 sqrt(V_i), searched radius), so a
 * one-sided stencil carries no O(h) curvature bias in the gradient. Falls
 * back to the linear face LSF when fewer than 8 points or ill-conditioned. */
static int gfs_lsf_order(void){
	static int o = -1;
	if(o < 0){ const char *s = getenv("EUNHA_LSF_ORDER"); o = (s && atoi(s) == 2) ? 2 : 1; }
	return o;
}
/* EUNHA_LIM_SMOOTH=1 (needs EUNHA_LSF_ORDER=2): in cells where the quadratic
 * fit explains the data (weighted rms residual <= 0.1 rms of the
 * differences), the limiter range [Qmin, Qmax] is widened by the curvature
 * over the stencil, max_j |e_j^T H e_j|/2, which is O(h^2) and zero for a
 * linear field; smooth extrema are then not clipped. */
static int gfs_lim_smooth_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("EUNHA_LIM_SMOOTH"); on = (s && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static long gfs_q2_n = 0, gfs_q2_fail = 0, gfs_q2_smooth[4] = {0, 0, 0, 0};
/* GFS_ORACLE_EDGE=1 (diagnostic, Kepler only): cells with a vacuum arc get
 * the exact Keplerian velocity gradient and no velocity limiting. */
/* modes: 1 = exact gradient + no velocity limiting, 2 = exact gradient
 * only (limiter as usual), 3 = no velocity limiting only (LSF gradient) */
static int gfs_oracle_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_ORACLE_EDGE"); on = (s && s[0]) ? atoi(s) : 0; if(on < 0) on = 0; }
	return on;
}
static double tq_acc[TQ_NB][3], tq_cum[TQ_NB][3], tq_cap_n[TQ_NB], tq_cap_len[TQ_NB];
static long tq_ncalls = 0;
static double lq_n[TQ_NB], lq_asum[TQ_NB], lq_low[TQ_NB];            /* limiter: cells, sum min(a_vx,a_vy), count <0.5 */
static double jq_n[TQ_NB], jq_cell[TQ_NB], jq_rec[TQ_NB], jq_keep[TQ_NB]; /* faces: |du_n| cell, reconstructed, count rec>0.5 cell */
static int gfs_tq_on(void){
	static int on = -1;
	if(on < 0){ const char *e = getenv("GFS_TORQUE_AUDIT"); on = (e && atoi(e) > 0) ? 1 : 0; }
	return on;
}
static void tq_center(double *cx, double *cy){
	const char *a = getenv("EUNHA_KEPLER_CX"), *b = getenv("EUNHA_KEPLER_CY");
	*cx = (a && a[0]) ? atof(a) : 0.0; *cy = (b && b[0]) ? atof(b) : 0.0;
}
static int tq_bin(double R){ int k = 0; while(k < TQ_NB-1 && R >= tq_edges[k]) k++; return k; }
static int gfs_dynw2_active = 0;
/* GFS_SCHEME_H=1 (with GFS_DYN_W2=1; default 0): hybrid scheme (H) of
 * Docs/gfs_scheme3_formulation.md §3, §19. Geometric (Laguerre) volume is
 * the thermodynamic volume; the face force carries the VPH tangential term
 * -A (P_j - P_i) t/d; the cell total-energy rate is v_i.f - P_i A w_f
 * + theta_i D_f with D_f = A[(P_i-P_j) delta + (P* - Ptilde)(u_i - u_j)],
 * so total energy is conserved for any Wdot. Wdot follows the Godunov
 * volume rate with relaxation to the Voronoi volume (tau_V, GFS_H_TAUV)
 * and error feedback (tau_C, GFS_H_TAUC). Validity |dW| < d^2 is measured,
 * not enforced ([SCHH]). CPU, one rank (as GFS_DYN_W2). */
/* GFS_EDGE_DUMP=1 (default 0, diagnostic): per-face records for particles
 * with r < GFS_EDGE_DUMP_R (centre EUNHA_KEPLER_CX/CY) during the time
 * windows [T1, T1+W) and [T2, T2+W) (GFS_EDGE_DUMP_T1/T2/W), appended to
 * edge_faces.txt. Time is accumulated by the RK4 driver (gfs_tnow). */
static double gfs_tnow = 0; static long gfs_edge_call = 0;
/* EUNHA_VORO_VACBC=1 (with EUNHA_VORO_RHO; default 0): the exposed arcs are a
 * physical vacuum boundary for density and pressure (rho_b = P_b = 0).
 * The LSF gradient of rho and P gets the exact arc integrals
 *   M_b = int r r^T/|r|^2 ds,  b_b = int r (0 - Q_i)/|r|^2 ds  (|r| = R),
 * and the limiter's lower bound for rho and P includes 0 on arc cells.
 * Velocity gets no boundary value (a free surface is not a wall).
 * Docs/consult_vacuum_bc_design_reply_astra.md (A). */
static int gfs_vacbc_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("EUNHA_VORO_VACBC"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static int gfs_twopass_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("EUNHA_TWOPASS_GRAD"); on = (s && s[0]) ? (atoi(s) > 0) : 1; }
	return on;
}
static long gfs_cellcheck_n = 0, gfs_cellcheck_fail = 0;
/* GFS_TOPO_AUDIT=1: face-neighbour set of every particle at RK stage 1;
 * stages 2-4 count particles whose set differs (topology change inside a
 * step, where the RK4 right-hand side is discontinuous). Keyed by PINDX. */
static int gfs_rk_stage = 0;
static long *gfs_topo_ids = NULL; static unsigned char *gfs_topo_cnt = NULL; static long gfs_topo_cap = 0;
static long gfs_topo_chg = 0, gfs_topo_n = 0;
static int gfs_topo_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_TOPO_AUDIT"); on = (s && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static int gfs_cmp_long(const void *a, const void *b){ long x = *(const long*)a, y = *(const long*)b; return (x > y) - (x < y); }   /* cells that may be cut by generators beyond the searched stencil */
static long gfs_recfloor_n = 0;   /* HLLC input floor (1e-10) activations, av_mode 5 face loop */
static long gfs_rf_negP = 0, gfs_rf_tinyP = 0, gfs_rf_negR = 0, gfs_rf_tinyR = 0, gfs_rf_arc = 0, gfs_rf_faces = 0;
static double gfs_rf_APs_floor = 0, gfs_rf_APs_all = 0;
/* magnitude of negative reconstructions relative to the cell value: bins
 * |q_rec|/q_cell in [0,1e-12), [1e-12,1e-6), [1e-6,1e-2), [1e-2,1), >=1 */
static long gfs_rf_mag[2][5];
static int gfs_lim_dump_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_LIM_DUMP"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static int gfs_edge_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_EDGE_DUMP"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
/* GFS_TCF=1: thermodynamically compatible face flux (Docs/plan_thermo_compatible_flux.md,
 * Docs/tcf_derivation.md; W=0 only). Reversible force from cell pressures,
 * dissipation Q = Z_eff minmod(du_rec, du_cell), antisymmetric energy exchange. */
static int gfs_tcf_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_TCF"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static int gfs_tcf_q0(void){   /* GFS_TCF_Q0=1: Q = 0 (reversible G2 test) */
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_TCF_Q0"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static FILE *tcf_fdump = NULL;   /* <file>.faces: per-face record */
static FILE *tcf_dfp = NULL;
static int tcf_dcall0 = -1, tcf_dcall1 = -1;
/* GFS_TCF_DUMP=<file>: right-hand side (and <file>.faces) of force calls
 * GFS_TCF_DUMP_CALL..GFS_TCF_DUMP_LAST (default 1..1), one file per call
 * (suffix .<call> when more than one), then exit. Called at the serial entry. */
static void tcf_dump_switch(int call){
	const char *b = getenv("GFS_TCF_DUMP");
	if(!b || !b[0]) return;
	if(tcf_dcall0 < 0){ const char *e_ = getenv("GFS_TCF_DUMP_CALL"); tcf_dcall0 = (e_ && e_[0]) ? atoi(e_) : 1;
		e_ = getenv("GFS_TCF_DUMP_LAST"); tcf_dcall1 = (e_ && e_[0]) ? atoi(e_) : tcf_dcall0; }
	if(tcf_dfp){ fclose(tcf_dfp); tcf_dfp = NULL; }
	if(tcf_fdump){ fclose(tcf_fdump); tcf_fdump = NULL; }
	if(call > tcf_dcall1){ fprintf(stderr, "[TCF] dump written\n"); exit(0); }
	if(call < tcf_dcall0) return;
	char fn[1024];
	if(tcf_dcall1 > tcf_dcall0) snprintf(fn, sizeof(fn), "%s.%d", b, call); else snprintf(fn, sizeof(fn), "%s", b);
	tcf_dfp = fopen(fn, "w");
	if(tcf_dfp) fprintf(tcf_dfp, "# id x y vx vy P V m fx fy dE TdS dVgeom sumAn_x sumAn_y sumAe_x sumAe_y maxbis/d max|e x n| minA/d nface gUxx gUxy gUyx gUyy\n");
	strncat(fn, ".faces", sizeof(fn) - strlen(fn) - 1);
	tcf_fdump = fopen(fn, "w");
	if(tcf_fdump) fprintf(tcf_fdump, "# id_i id_j x_i y_i A d ex ey du_rec du_cell Z c_i uL_rec uR_rec Q\n");
}
static FILE *tcf_dump_fp(void){ return tcf_dfp; }
/* GFS_TCF_SACC=<file>: per-cell RK4 integral of d ln K/dt = T dS/U from the
 * face dissipation (G2 entropy-error test), rewritten after every step. */
static double *tcf_rate = NULL, *tcf_sacc = NULL; static size_t tcf_cap = 0;
static double gfs_dt_cap = 0, gfs_step_dtmin = 1e300;   /* step rejection (exam2d_vph_rk4_int_blend) */
static int gfs_attempt_failed = 0, gfs_in_retry = 0;   /* set by any stage floor / non-positive state */
static double gfs_dtb[4] = {1e300, 1e300, 1e300, 1e300}; static long gfs_dtb_cnt[4] = {0, 0, 0, 0};   /* face CFL, face geometry, dt3, HG */
static const char *tcf_sacc_file(void){ const char *s = getenv("GFS_TCF_SACC"); return (s && s[0]) ? s : NULL; }
static void tcf_sacc_ensure(treevorostressrk4particletype *sbp, int np){
	size_t mx = 0, k; int i;
	for(i = 0; i < np; i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(sbp+i)); if(id != (size_t)MAX_INDEX && id + 1 > mx) mx = id + 1; }
	if(mx <= tcf_cap) return;
	tcf_rate = (double*)realloc(tcf_rate, mx*sizeof(double)); tcf_sacc = (double*)realloc(tcf_sacc, mx*sizeof(double));
	for(k = tcf_cap; k < mx; k++) tcf_rate[k] = tcf_sacc[k] = 0;
	tcf_cap = mx;
}
static void tcf_sacc_stage(treevorostressrk4particletype *sbp, int np, double wdt){
	int i;
	for(i = 0; i < np; i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(sbp+i)); if(id < tcf_cap) tcf_sacc[id] += wdt*tcf_rate[id]; }
}
static void tcf_sacc_write(treevorostressrk4particletype *sbp, int np){
	FILE *fp = fopen(tcf_sacc_file(), "w"); int i;
	if(!fp) return;
	fprintf(fp, "# id int_0^t TdS/U dt\n");
	for(i = 0; i < np; i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(sbp+i)); if(id < tcf_cap) fprintf(fp, "%ld %.17g\n", (long)id, tcf_sacc[id]); }
	fclose(fp);
}
static int gfs_tcf_notan(void){   /* GFS_TCF_NOTAN=1: diagnostic, tangential force removed (Astra 26) */
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_TCF_NOTAN"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static double tcf_tanwork = 0;
/* GFS_TCF_ZH=1: impedance floor Z >= Z_H = (rho_h/2)[a + sqrt(a^2 + c^2)], a = (gamma+1)|m|/8,
 * rho_h = 2 rho_i rho_j/(rho_i + rho_j), c = max(c_i, c_j): symmetric two-shock Hugoniot form
 * (Docs/subzonal_derivation.md 3.2, section 6); m = minmod(du_rec, du_cell). */
static int gfs_tcf_zh(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_TCF_ZH"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
/* GFS_TCF_QLAW=2: Q = Z m^2/du_cell, m = minmod(du_rec, du_cell) (0 when du_cell = 0).
 * Heat A Q du_cell = A Z m^2 >= 0 is second order in the reconstruction jump when
 * |du_rec| < |du_cell| (mean shear across the face, Fable consult 26) and equals
 * the minmod law when |du_rec| >= |du_cell|. Default 1: Q = Z m. */
static int gfs_tcf_qlaw(void){
	static int q = -1;
	if(q < 0){ const char *s = getenv("GFS_TCF_QLAW"); q = (s && s[0]) ? atoi(s) : 1; if(q != 2) q = 1; }
	return q;
}   /* sum_i v_i . F_t,i of the removed tangential force (cell adiabatic residual) */
/* GFS_SUBZONE=1: subzonal (centroid-referenced) potential Phi on top of the cell-volume
 * forces (Docs/subzonal_derivation.md 3.2, signed; reference Docs/subzonal_checks/sub2d.py).
 * W=0, no disc, one MPI rank; cells with cap or wall edges get Phi = 0 (counted).
 * Cell pressure P_i0 = (gamma-1) ie/(V g), g = sum_f (B_f(c)/V)(h_f(c)/h_f(x))^(gamma-1). */
static int gfs_subzone_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_SUBZONE"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static double *sub_g = NULL, *sub_fx = NULL, *sub_fy = NULL, *sub_dphi = NULL, *sub_phi = NULL;
static size_t sub_cap = 0;
static long sub_ncap = 0, sub_ncell = 0, sub_nbad = 0;
static double sub_vf_last = 0;   /* sum_i v_i.F_Phi,i at the last force call (same velocities as dPhi) */
static void sub_ensure(SimParameters *simpar){
	size_t ps = TVORORK4_DDINFO(simpar)[0].n_size, mx = 0, k; char *b = (char*)VORORK4_TBP(simpar); int i;
	for(i = 0; i < VORO_NP(simpar); i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(b + i*ps)); if(id != (size_t)MAX_INDEX && id + 1 > mx) mx = id + 1; }
	if(mx <= sub_cap) return;
	sub_g = (double*)realloc(sub_g, mx*sizeof(double)); sub_fx = (double*)realloc(sub_fx, mx*sizeof(double));
	sub_fy = (double*)realloc(sub_fy, mx*sizeof(double)); sub_dphi = (double*)realloc(sub_dphi, mx*sizeof(double));
	sub_phi = (double*)realloc(sub_phi, mx*sizeof(double));
	for(k = sub_cap; k < mx; k++){ sub_g[k] = 1; sub_fx[k] = sub_fy[k] = sub_dphi[k] = sub_phi[k] = 0; }
	sub_cap = mx;
}
static double sub_gfac(treevorork4particletype *p){
	size_t id = (size_t)PINDX(p);
	return (gfs_subzone_on() && id < sub_cap && sub_g[id] > 0) ? sub_g[id] : 1.0;
}
#define SUB_MAXC 256
typedef struct {
	int n;
	double ax[SUB_MAXC], ay[SUB_MAXC], nx[SUB_MAXC], ny[SUB_MAXC];
	treevorork4particletype *up[SUB_MAXC];   /* neighbour across edge k = (corner k, corner k+1) */
} sub_poly;
/* corners relative to the generator and edge normals: toward the neighbour on a pair
 * face; the outward geometric normal on a cap edge (box side, moves rigidly with the
 * generator; S->up = NULL). 0 for an arc edge (disc), a wall ghost, or too many corners. */
/* GFS_HG_ARC=1 (Astra audit 15): HG in disc-clipped cells. Arc edges (VORO_ARC_EDGE) become NULL
 * neighbours in sub_poly; hg_cell then uses the exact clipped area/centroid (circular-arc moments),
 * a generic centroid Jacobian (Voronoi vertices + bisector-circle points), skips arc edges in the face
 * modes, and adds one cell-owned mode per connected arc component with xi = log(B(c)/W). Split only. */
static int gfs_hg_arc(void){ static int on = -1; if(on < 0){ const char *s_ = getenv("GFS_HG_ARC"); on = (s_ && s_[0] && atoi(s_) > 0) ? 1 : 0; } return on; }
static double hg_arc_R = 0, hg_arc_chk = 0; static long hg_arc_nmode = 0, hg_arc_ncell = 0, hg_arc_fail = 0;
static int sub_build(Voro2D_Corner *vc, Voro2D_point *nw, treevorork4particletype *ibp, sub_poly *S){
	Voro2D_Corner *t = vc; int n = 0;
	do {
		if(n >= SUB_MAXC || (t->upperrelated == VORO_ARC_EDGE && !gfs_hg_arc())) return 0;
		S->ax[n] = (double)t->x; S->ay[n] = (double)t->y;
		if(t->upperrelated >= 0){
			treevorork4particletype *nb = (treevorork4particletype*)nw[t->upperrelated].bp;
			if(!nb || PINDX(nb) == MAX_INDEX) return 0;
			double dx = (double)nb->x - (double)ibp->x, dy = (double)nb->y - (double)ibp->y, d = sqrt(dx*dx + dy*dy);
			if(!(d > 0)) return 0;
			S->nx[n] = dx/d; S->ny[n] = dy/d; S->up[n] = nb;
		} else {
			double ex = (double)t->upperlink->x - (double)t->x, ey = (double)t->upperlink->y - (double)t->y, L = sqrt(ex*ex + ey*ey);
			if(!(L > 0)) return 0;
			S->nx[n] = ey/L; S->ny[n] = -ex/L;   /* outward for a CCW polygon */
			S->up[n] = NULL;
		}
		n++;
		t = t->upperlink;
	} while(t != vc);
	S->n = n;
	return n >= 3;
}
/* g, Phi = U0 (g - 1), and (if gv) the gradients of Phi w.r.t. the generator (gP) and corners (gv) */
static int sub_eval(const sub_poly *S, double m, double P0, double gam,
		double *g_out, double *phi_out, double *gP, double *gvx, double *gvy){
	int n = S->n, k;
	double V = 0, Jx = 0, Jy = 0;
	for(k = 0; k < n; k++){
		int k1 = (k+1) % n; double cr = S->ax[k]*S->ay[k1] - S->ay[k]*S->ax[k1];
		V += 0.5*cr; Jx += cr*(S->ax[k] + S->ax[k1])/6.0; Jy += cr*(S->ay[k] + S->ay[k1])/6.0;
	}
	if(!(V > 0)) return 0;
	double cx = Jx/V, cy = Jy/V, g = 0;
	double r[SUB_MAXC], B[SUB_MAXC];
	for(k = 0; k < n; k++){
		int k1 = (k+1) % n;
		double hx = S->ax[k]*S->nx[k] + S->ay[k]*S->ny[k];
		double hc = (S->ax[k] - cx)*S->nx[k] + (S->ay[k] - cy)*S->ny[k];
		if(!(hx > 0)) return 0;
		r[k] = hc/hx;
		B[k] = 0.5*((S->ax[k] - cx)*(S->ay[k1] - cy) - (S->ay[k] - cy)*(S->ax[k1] - cx));
		g += (B[k]/V)*pow(r[k] > 0 ? r[k] : 0, gam - 1);
	}
	*g_out = g;
	if(phi_out) *phi_out = P0*V/(gam - 1)*(g - 1);
	if(!gvx) return 1;
	double rho = m/V, gcx = 0, gcy = 0, HM = 0;
	for(k = 0; k < n; k++){ gvx[k] = gvy[k] = 0; }
	gP[0] = gP[1] = 0;
	for(k = 0; k < n; k++){
		int k1 = (k+1) % n;
		double rk = r[k] > 0 ? r[k] : 0;
		double Pf = P0*pow(rk, gam), rf = rho*rk, H = (rf > 0) ? gam/(gam - 1)*Pf/rf : 0, M = m*B[k]/V;
		HM += H*M;
		double ax = S->ax[k], ay = S->ay[k], bx = S->ax[k1], by = S->ay[k1];
		/* -P_f dW_f, W_f = 1/2 a x b (generator at the origin): R(w) = (w_y, -w_x), L(u) = (-u_y, u_x) */
		gvx[k]  += -Pf*0.5*by;   gvy[k]  += -Pf*0.5*(-bx);
		gvx[k1] += -Pf*0.5*(-ay); gvy[k1] += -Pf*0.5*ax;
		gP[0] += -Pf*0.5*(ay - by); gP[1] += -Pf*0.5*(-(ax - bx));
		/* (m/V) H_f dB_f */
		double sF = m/V*H;
		gvx[k]  += sF*0.5*(by - cy);    gvy[k]  += sF*0.5*(-(bx - cx));
		gvx[k1] += sF*0.5*(-(ay - cy)); gvy[k1] += sF*0.5*(ax - cx);
		gcx += sF*0.5*(ay - by); gcy += sF*0.5*(-(ax - bx));
	}
	for(k = 0; k < n; k++){
		int kn = (k+1) % n, kp = (k+n-1) % n;
		double ax = S->ax[k], ay = S->ay[k], nx_ = S->ax[kn], ny_ = S->ay[kn], px = S->ax[kp], py = S->ay[kp];
		double dVx = 0.5*(ny_ - py), dVy = 0.5*(-(nx_ - px));
		/* -(sum H M / V) dV, and +P0 dV removes the single-state part (Phi = U_sub - U0) */
		gvx[k] += (-HM/V + P0)*dVx; gvy[k] += (-HM/V + P0)*dVy;
		/* centroid chain: dc = (dJ - c dV)/V */
		double ck = ax*ny_ - ay*nx_, cp = px*ay - py*ax;
		double J00 = ((ax + nx_)*ny_ + ck + (px + ax)*(-py) + cp)/6.0, J01 = ((ax + nx_)*(-nx_) + (px + ax)*px)/6.0;
		double J10 = ((ay + ny_)*ny_ + (py + ay)*(-py))/6.0,             J11 = ((ay + ny_)*(-nx_) + ck + (py + ay)*px + cp)/6.0;
		J00 = (J00 - cx*dVx)/V; J01 = (J01 - cx*dVy)/V; J10 = (J10 - cy*dVx)/V; J11 = (J11 - cy*dVy)/V;
		gvx[k] += J00*gcx + J10*gcy; gvy[k] += J01*gcx + J11*gcy;
	}
	return 1;
}
/* GFS_SHARED_GEOM=1 (design_astra.md step 1): every pair face is used with the segment
 * common to both cells (the two cells cut the shared bisector with their own boxes or
 * discs); the rest of each cell's edge is vacuum. Endpoints are recorded per cell in the
 * stage geometry pass (updateDenW2Pressure2DBlend pass 0, same cells as the force loop). */
static int gfs_shared_geom_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_SHARED_GEOM"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
#define SHF_K 64
static int *shf_n = NULL; static long *shf_j = NULL; static double *shf_g = NULL;   /* g: dx dy ax ay bx by (rel. to cell) */
static size_t shf_cap = 0;
static long shf_nlook = 0, shf_nmiss = 0, shf_ncut = 0;
static double shf_maxcut = 0;   /* largest cut fraction of a face with A > 1e-6 d (shared-geometry mismatch) */
static void shf_ensure(SimParameters *simpar){
	size_t ps = TVORORK4_DDINFO(simpar)[0].n_size, mx = 0, k; char *b = (char*)VORORK4_TBP(simpar); int i;
	for(i = 0; i < VORO_NP(simpar); i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(b + i*ps)); if(id != (size_t)MAX_INDEX && id + 1 > mx) mx = id + 1; }
	if(mx <= shf_cap) return;
	shf_n = (int*)realloc(shf_n, mx*sizeof(int)); shf_j = (long*)realloc(shf_j, mx*SHF_K*sizeof(long));
	shf_g = (double*)realloc(shf_g, mx*SHF_K*6*sizeof(double));
	for(k = shf_cap; k < mx; k++) shf_n[k] = 0;
	shf_cap = mx;
}
static void shf_record(Voro2D_Corner *vc, Voro2D_point *nw, treevorork4particletype *ibp){
	size_t id = (size_t)PINDX(ibp); Voro2D_Corner *t = vc; int n = 0;
	if(id >= shf_cap) return;
	do {
		if(t->upperrelated >= 0 && n < SHF_K){
			treevorork4particletype *nb = (treevorork4particletype*)nw[t->upperrelated].bp;
			if(nb && PINDX(nb) != MAX_INDEX){
				double *g = shf_g + (id*SHF_K + n)*6;
				shf_j[id*SHF_K + n] = (long)PINDX(nb);
				g[0] = (double)nb->x - (double)ibp->x; g[1] = (double)nb->y - (double)ibp->y;
				g[2] = (double)t->x; g[3] = (double)t->y; g[4] = (double)t->upperlink->x; g[5] = (double)t->upperlink->y;
				n++;
			}
		}
		t = t->upperlink;
	} while(t != vc);
	shf_n[id] = n;
}
/* common segment of face (i -> j): own endpoints a, b and the separation d (all relative
 * to x_i). Returns the common length and its centroid (relative to x_i). */
static double shf_common(long idi, long idj, double dx, double dy, double ax, double ay, double bx, double by,
		double *cx, double *cy){
	double d = sqrt(dx*dx + dy*dy), tx = -dy/d, ty = dx/d, mx = 0.5*dx, my = 0.5*dy;
	double s1 = (ax - mx)*tx + (ay - my)*ty, s2 = (bx - mx)*tx + (by - my)*ty;
	double lo = fmin(s1, s2), hi = fmax(s1, s2);
	*cx = mx + tx*0.5*(lo + hi); *cy = my + ty*0.5*(lo + hi);
	if((size_t)idj >= shf_cap) return hi - lo;
	int k, n = shf_n[idj], found = 0;
	for(k = 0; k < n; k++){
		if(shf_j[(size_t)idj*SHF_K + k] != idi) continue;
		double *g = shf_g + ((size_t)idj*SHF_K + k)*6;
		if(fabs(g[0] + dx) + fabs(g[1] + dy) > 1e-9*d) continue;   /* another periodic image */
		/* j's endpoints relative to x_i */
		double s3 = (g[2] + dx - mx)*tx + (g[3] + dy - my)*ty, s4 = (g[4] + dx - mx)*tx + (g[5] + dy - my)*ty;
		double l2 = fmax(lo, fmin(s3, s4)), h2 = fmin(hi, fmax(s3, s4));
		if(hi - lo > 1e-6*d){
			double cutf = ((hi - lo) - fmax(h2 - l2, 0))/(hi - lo);
			if(cutf > shf_maxcut){
#ifdef _OPENMP
#pragma omp critical (shf_maxcut)
#endif
				{ if(cutf > shf_maxcut) shf_maxcut = cutf; }
			}
		}
		if(l2 != lo || h2 != hi){
#ifdef _OPENMP
#pragma omp atomic
#endif
			shf_ncut++;
		}
		lo = l2; hi = h2; found = 1; break;
	}
#ifdef _OPENMP
#pragma omp atomic
#endif
	shf_nlook++;
	if(!found){
#ifdef _OPENMP
#pragma omp atomic
#endif
		shf_nmiss++;
		return 0;   /* the face exists on one side only: vacuum for both */
	}
	if(hi <= lo) return 0;
	*cx = mx + tx*0.5*(lo + hi); *cy = my + ty*0.5*(lo + hi);
	return hi - lo;
}
/* GFS_TCF_FACE_Z=1 (design step 5): Q impedance from the face subzonal states
 * rho_if = rho_i h_f(c_i)/h_f(x_i) = rho_i (1 - 2 delta_i.e/d), delta_i = c_i - x_i,
 * P_if = P_i (rho_if/rho_i)^gamma; Z = max(Z_Davis, Z_H), Z_Davis = mu (c_M + max(m,0)),
 * Z_H = mu [a + sqrt(a^2 + c_M^2)], a = (gamma+1)|m|/8, mu = rho_L rho_R/(rho_L + rho_R). */
static int gfs_tcf_face_z(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_TCF_FACE_Z"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
/* GFS_PERIODIC_COMPLETE=1 (design step 2): a periodic run must not have active cap edges */
static int gfs_periodic_complete(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_PERIODIC_COMPLETE"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static double *cen_dx = NULL, *cen_dy = NULL; static size_t cen_cap = 0;
static long pc_ncap = 0; static double pc_vsum = 0;
static void cen_ensure(SimParameters *simpar){
	size_t ps = TVORORK4_DDINFO(simpar)[0].n_size, mx = 0, k; char *b = (char*)VORORK4_TBP(simpar); int i;
	for(i = 0; i < VORO_NP(simpar); i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(b + i*ps)); if(id != (size_t)MAX_INDEX && id + 1 > mx) mx = id + 1; }
	if(mx <= cen_cap) return;
	cen_dx = (double*)realloc(cen_dx, mx*sizeof(double)); cen_dy = (double*)realloc(cen_dy, mx*sizeof(double));
	for(k = cen_cap; k < mx; k++) cen_dx[k] = cen_dy[k] = 0;
	cen_cap = mx;
}
/* GFS_TCF_QBP=1 (Astra audit 12, needs EUNHA_LSF_ORDER=2): Q input m = minmod(delta_b', du_cell),
 * delta_b' = du_cell + d/2 e^T (G_i + G_j) e with the UNLIMITED quadratic-fit velocity gradients
 * (cached before the limiter). Rotation invariant; zero for linear and quadratic velocity fields. */
static __thread double q2r_cd2[4096], q2r_cph[4096][9], q2r_cdq[4096][4], q2r_Qc[4096][9]; static __thread int q2r_cord[4096];
/* GFS_Q3=1 (Astra audit 14, item 1): cubic fit (9 rotation-covariant columns) for the unlimited velocity
 * gradient used by QBP, in every cell (stateless, so restart stays exact); quadratic fallback if rank < 9. */
static int gfs_q3(void){ static int on = -1; if(on < 0){ const char *s_ = getenv("GFS_Q3"); on = (s_ && s_[0] && atoi(s_) > 0) ? 1 : 0; } return on; }
static long q3_n = 0, q3_fb = 0;
static int gfs_qbp_chi(void){ static int on = -1; if(on < 0){ const char *s_ = getenv("GFS_QBP_CHI"); on = (s_ && s_[0] && atoi(s_) > 0) ? 1 : 0; } return on; }
static int gfs_q2_rank(void){ static int on = -1; if(on < 0){ const char *s_ = getenv("GFS_Q2_RANK"); on = (s_ && s_[0] && atoi(s_) > 0) ? 1 : 0; } return on; }
static int gfs_tcf_qbp(void){ static int on = -1; if(on < 0){ const char *s_ = getenv("GFS_TCF_QBP"); on = (s_ && s_[0] && atoi(s_) > 0) ? 1 : 0; } return on; }
static double *qbp_g = NULL; static char *qbp_ok = NULL; static size_t qbp_cap = 0; static long qbp_nfail = 0, qbp_nface = 0;
/* GFS_QBP_EDGE_INNER=1 (Astra audit 20 (b)): a cell with a vacuum arc has a one-sided velocity fit; on a face
 * between an arc cell and a cell without arc, delta and chi use the arc-free cell's gradient only (same rule on both sides). */
/* =1: the arc-free gradient as is (rejected: first-order midpoint error in curved shear);
 * =2: the arc-free cell's quadratic fit Taylor-extrapolated to the arc cell, G' = G + H.(x_arc - x_in), exact for quadratic v */
static int gfs_qbp_edge_inner(void){ static int on = -1; if(on < 0){ const char *s_ = getenv("GFS_QBP_EDGE_INNER"); on = (s_ && s_[0]) ? atoi(s_) : 0; if(on < 0) on = 0; } return on; }
static double *qbp_h = NULL;   /* velocity Hessians (vx: Hxx Hxy Hyy, vy: Hxx Hxy Hyy); zero for the linear fallback */
static double *qbp_t = NULL;   /* =3: velocity third derivatives (vx: Txxx Txxy Txyy Tyyy, vy: same) from the GFS_Q3 cubic fit; zero otherwise */
static char *qbp_arc = NULL; static long qbp_narc = 0, qbp_nmixed = 0;
/* qbp_ok: 0 invalid, 1 unlimited linear LSF gradient (fallback, Astra audit 13 (b)), 2 unlimited quadratic fit */
static long qbp_src[3] = {0, 0, 0}, qbp_q2fail_npt = 0, qbp_q2fail_piv = 0;
static void qbp_ensure(SimParameters *simpar){
	size_t ps = TVORORK4_DDINFO(simpar)[0].n_size, mx = 0, k; char *b = (char*)VORORK4_TBP(simpar); int i;
	for(i = 0; i < VORO_NP(simpar); i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(b + i*ps)); if(id != (size_t)MAX_INDEX && id + 1 > mx) mx = id + 1; }
	if(mx > qbp_cap){ qbp_g = (double*)realloc(qbp_g, 4*mx*sizeof(double)); qbp_ok = (char*)realloc(qbp_ok, mx); qbp_arc = (char*)realloc(qbp_arc, mx); qbp_h = (double*)realloc(qbp_h, 6*mx*sizeof(double)); qbp_t = (double*)realloc(qbp_t, 8*mx*sizeof(double)); qbp_cap = mx; }
	for(k = 0; k < qbp_cap; k++){ qbp_ok[k] = 0; qbp_arc[k] = 0; }
	for(k = 0; k < 6*qbp_cap; k++) qbp_h[k] = 0;
	for(k = 0; k < 8*qbp_cap; k++) qbp_t[k] = 0;
}
/* polygon centroid offset c - x (corners relative to the generator) and cap-edge count */
static int cen_poly(Voro2D_Corner *vc, double *dx, double *dy){
	Voro2D_Corner *t = vc; double V = 0, Jx = 0, Jy = 0; int ncap = 0;
	do {
		Voro2D_Corner *u = t->upperlink;
		double cr = (double)t->x*(double)u->y - (double)t->y*(double)u->x;
		V += 0.5*cr; Jx += cr*((double)t->x + (double)u->x)/6.0; Jy += cr*((double)t->y + (double)u->y)/6.0;
		if(t->upperrelated < 0) ncap++;
		t = u;
	} while(t != vc);
	*dx = (V > 0) ? Jx/V : 0; *dy = (V > 0) ? Jy/V : 0;
	return ncap;
}
/* GFS_HG_AFFINE=1 (design_astra.md steps 3-4 with the audit corrections): affine-filtered
 * hourglass dissipation, no potential. Per pair face xi = ln(1-u), u = 2 delta.Delta/d^2
 * (delta = c_i - x_i, Delta = x_j - x_i); g_b = grad_b xi through the centroid (circumcentre
 * Jacobians) plus the explicit (normal-rotation) terms; 2-hop weighted LSQ affine filter
 * l_b = g_b - w_b M S^-1 ytilde_b; r = sum l_b.v_b; C = rho_f h W [b|r| + sqrt(b^2 r^2 + c_f^2)];
 * F_b = -C r l_b, heat C r^2 in cell i. Reference: Docs/subzonal_checks/hg2d_ref.py.
 * Needs GFS_SHARED_GEOM (2-hop data) and one MPI rank; cells with cap/arc edges are skipped. */
static int gfs_lim_proj(void){ static int on = -1; if(on < 0){ const char *s_ = getenv("GFS_LIM_PROJ"); on = (s_ && s_[0] && atoi(s_) > 0) ? 1 : 0; } return on; }
static int gfs_hg_affine(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_HG_AFFINE"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static double *hg_fx = NULL, *hg_fy = NULL, *hg_H = NULL; static size_t hg_cap = 0;
static treevorork4particletype **hg_ptr = NULL;   /* id -> local particle */
static long hg_ncell = 0, hg_nskip = 0, hg_npinv = 0; static double hg_vf_last = 0, hg_H_last = 0;
static double *hg_lam = NULL;   /* row-sum bound of M^-1/2 (D) M^-1/2 per particle (RK4 damping limit) */
static double hg_dtlim = 1e300;
static double hg_worst[12] = {0};   /* exact single-face eigenvalue k sum|z|^2, d/sqrtV, 1-u, |g M^-1/2|, |l M^-1/2|, cond, |r|h/c, ns, id, C/acou */
/* face operators of the current force evaluation (Astra audit 2: spectrum diagnosis) */
static long *hgo_id = NULL; static double *hgo_l = NULL; static long hgo_n = 0, hgo_cap = 0;
static long *hgf_start = NULL, *hgf_len = NULL; static double *hgf_k = NULL; static long hgf_n = 0, hgf_cap = 0;
static double hg_spec_lam = 0, hg_spec_res = 0, hg_spec_bound = 0;
static int gfs_hg_spec(void){ static int on = -1; if(on < 0){ const char *s = getenv("GFS_HG_SPEC"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; } return on; }
static void hgo_reset(long ncell_est){
	if(!gfs_hg_spec()) return;
	long need_f = ncell_est*16 + 64, need_o = need_f*48;
	if(need_f > hgf_cap){ hgf_start = (long*)realloc(hgf_start, (need_f + 1)*sizeof(long)); hgf_len = (long*)realloc(hgf_len, need_f*sizeof(long)); hgf_k = (double*)realloc(hgf_k, need_f*sizeof(double)); hgf_cap = need_f; }
	if(need_o > hgo_cap){ hgo_id = (long*)realloc(hgo_id, need_o*sizeof(long)); hgo_l = (double*)realloc(hgo_l, 2*need_o*sizeof(double)); hgo_cap = need_o; }
	hgo_n = 0; hgf_n = 0;
}
   /* largest per-face Lambda contribution: lam, d/sqrtV, 1-u, max|l| sqrtV, C/(rho c h W), |r| h/c, ns, id */
static long fz_nclamp = 0;
static void bp_reset(void); static int gfs_bprime(void); static int hgm_building;
static void hg_ensure(SimParameters *simpar){
	size_t ps = TVORORK4_DDINFO(simpar)[0].n_size, mx = 0, k; char *b = (char*)VORORK4_TBP(simpar); int i;
	for(i = 0; i < VORO_NP(simpar); i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(b + i*ps)); if(id != (size_t)MAX_INDEX && id + 1 > mx) mx = id + 1; }
	if(mx > hg_cap){
		hg_fx = (double*)realloc(hg_fx, mx*sizeof(double)); hg_fy = (double*)realloc(hg_fy, mx*sizeof(double));
		hg_H = (double*)realloc(hg_H, mx*sizeof(double));
		hg_ptr = (treevorork4particletype**)realloc(hg_ptr, mx*sizeof(void*));
		hg_lam = (double*)realloc(hg_lam, mx*sizeof(double));
		hg_cap = mx;
	}
	for(k = 0; k < hg_cap; k++){ hg_fx[k] = hg_fy[k] = hg_H[k] = hg_lam[k] = 0; hg_ptr[k] = NULL; }
	for(i = 0; i < VORO_NP(simpar); i++){
		treevorork4particletype *q = (treevorork4particletype*)(b + i*ps); size_t id = (size_t)PINDX(q);
		if(id < hg_cap) hg_ptr[id] = q;
	}
	if(hgm_building && gfs_bprime()) bp_reset();
}
/* circumcentre Jacobians of the corner (i, p1, p2) (p relative to x_i, v the corner): J_i, J_1, J_2 (row-major 2x2) */
static int hg_vjac(double p1x, double p1y, double p2x, double p2y, double vx, double vy, double *Ji, double *J1, double *J2){
	double a00 = 2*p1x, a01 = 2*p1y, a10 = 2*p2x, a11 = 2*p2y, det = a00*a11 - a01*a10;
	if(!(fabs(det) > 0)) return 0;
	double i00 = a11/det, i01 = -a01/det, i10 = -a10/det, i11 = a00/det;
	Ji[0] = (i00 + i01)*2*vx; Ji[1] = (i00 + i01)*2*vy; Ji[2] = (i10 + i11)*2*vx; Ji[3] = (i10 + i11)*2*vy;
	double qx = 2*(p1x - vx), qy = 2*(p1y - vy), rx = 2*(p2x - vx), ry = 2*(p2y - vy);
	J1[0] = i00*qx; J1[1] = i00*qy; J1[2] = i10*qx; J1[3] = i10*qy;
	J2[0] = i01*rx; J2[1] = i01*ry; J2[2] = i11*rx; J2[3] = i11*ry;
	return 1;
}
/* GFS_HG_SPLIT=1 (Astra audit 3): HG dissipation integrated as a separate step,
 * HG(dt/2) -> RK4 without HG (dt) -> HG(dt/2); each mode (owner cell i, face f) is advanced
 * by the exact scalar flow dr/dt = -alpha C(|r|, U_i) r with the kinetic-energy loss put
 * into U_i; modes applied forward then backward with tau/2 each, in a fixed order. */
static int gfs_hg_split(void){ static int on = -1; if(on < 0){ const char *s_ = getenv("GFS_HG_SPLIT"); on = (s_ && s_[0] && atoi(s_) > 0) ? 1 : 0; } return on; }
static int hgm_building;
/* GFS_BULK_AFFINE=1 (Astra design 5): bulk dissipation on the compressive part of the
 * 2-hop QR affine divergence d_i = tr G_i = L_i v, eta_i = m h [c + (gamma+1)/4 h |a|],
 * h = sqrt(V); integrated as a one-sided split mode (only d < 0). Requires GFS_HG_SPLIT. */
static int gfs_bulk_affine(void){ static int on = -1; if(on < 0){ const char *s_ = getenv("GFS_BULK_AFFINE"); on = (s_ && s_[0] && atoi(s_) > 0) ? 1 : 0; } return on; }
static long bulk_nmode = 0, bulk_nskip = 0; static double bulk_heat = 0;
static double bulk_test_err = 0;
/* GFS_BPRIME=1 (Astra design 6, sign corrected): per shared face the non-affine normal-velocity residual
 * delta_f = (v_i - v_j).e + d e^T (G_i + G_j)/2 e (G from the same 2-hop QR), C_f = A_f Z_f with face-state
 * Z = max(mu (c_M + max(delta,0)), mu (a_H + sqrt(a_H^2 + c_M^2))), a_H = (g+1)/8 |delta|; heat split 1/2 : 1/2.
 * Integrated as split modes (two-owner exact flow). Requires GFS_HG_SPLIT (HG itself optional via GFS_HG_OFF). */
static int gfs_bprime(void){ static int on = -1; if(on < 0){ const char *s_ = getenv("GFS_BPRIME"); on = (s_ && s_[0] && atoi(s_) > 0) ? atoi(s_) : 0; } return on; }
/* GFS_BPRIME=2 (Astra design 7): C_f = A mu (g+1)/8 |delta| (no acoustic term); exact flow delta1 = delta0/(1+beta) */
/* GFS_HG_QUAD=1 (Astra audit 21 (a)): project HG modes against [1, y, y^2] (6 columns) instead of [1, y],
 * so that a quadratic (curved-shear) velocity field leaves no HG residual. GFS_HG_QUAD_TEST=1 logs max |r| for a quadratic field. */
static int gfs_hg_quad(void){ static int on = -1; if(on < 0){ const char *s_ = getenv("GFS_HG_QUAD"); on = (s_ && s_[0] && atoi(s_) > 0) ? 1 : 0; } return on; }
static long hg_nqdef = 0; static double hg_qtest = 0;
static int gfs_hg_off(void){ static int on = -1; if(on < 0){ const char *s_ = getenv("GFS_HG_OFF"); on = (s_ && s_[0] && atoi(s_) > 0) ? 1 : 0; } return on; }
#define BP_MAXF 32
static long *bp_sl_start = NULL; static int *bp_sl_len = NULL, *bp_nf = NULL; static size_t bp_cap = 0;
static long *bp_sid = NULL; static double *bp_s = NULL; static long bp_pn = 0, bp_pcap = 0;   /* s: yx, yy, bx, by */
static long *bp_fj = NULL; static double *bp_f = NULL;   /* per face: ex, ey, d, A, rf, kap, ok */
static long bp_nmode = 0, bp_nskip = 0; static double bp_heat = 0, bp_test_err = 0;
static void bp_reset(void){
	if(hg_cap > bp_cap){
		bp_sl_start = (long*)realloc(bp_sl_start, hg_cap*sizeof(long)); bp_sl_len = (int*)realloc(bp_sl_len, hg_cap*sizeof(int));
		bp_nf = (int*)realloc(bp_nf, hg_cap*sizeof(int));
		bp_fj = (long*)realloc(bp_fj, hg_cap*BP_MAXF*sizeof(long)); bp_f = (double*)realloc(bp_f, hg_cap*BP_MAXF*7*sizeof(double));
		bp_cap = hg_cap;
	}
	if(bp_pcap < (long)hg_cap*64){ bp_pcap = (long)hg_cap*64; bp_sid = (long*)realloc(bp_sid, bp_pcap*sizeof(long)); bp_s = (double*)realloc(bp_s, 4*bp_pcap*sizeof(double)); }
	size_t k; for(k = 0; k < bp_cap; k++){ bp_sl_len[k] = -1; bp_nf[k] = 0; }
	bp_pn = 0;
}   /* GFS_BULK_TEST: max |L v - tr A| for random affine v = A y (pre-merge stencil) */   /* 1: hg_cell stores modes instead of forces */
typedef struct { long owner, key, start, len; double A0, b, kap; int onesided; long owner2; double mu, kap2; int kind; } hg_mode;   /* kind 0 HG, 1 bulk, 2 b' */
static hg_mode *hgm = NULL; static long hgm_n = 0, hgm_cap = 0;
static long *hgm_id = NULL; static double *hgm_l = NULL; static long hgm_on = 0, hgm_ocap = 0;
static long hgm_nfail = 0;
static void hgm_reset(long np){
	long nf = np*16 + 64, no = nf*48;
	if(nf > hgm_cap){ hgm = (hg_mode*)realloc(hgm, nf*sizeof(hg_mode)); hgm_cap = nf; }
	if(no > hgm_ocap){ hgm_id = (long*)realloc(hgm_id, no*sizeof(long)); hgm_l = (double*)realloc(hgm_l, 2*no*sizeof(double)); hgm_ocap = no; }
	hgm_n = 0; hgm_on = 0; hgm_nfail = 0;
}
static int hgm_cmp(const void *a, const void *b){ long ka = ((const hg_mode*)a)->key, kb = ((const hg_mode*)b)->key; return (ka > kb) - (ka < kb); }
#define HG_MAXS 512
static void hg_cell(const sub_poly *S, treevorork4particletype *ibp, double rho, double P, double gam){
	int n = S->n, k, q;
	double ix = (double)ibp->x, iy = (double)ibp->y;
	/* cell area, centroid and its corner Jacobians */
	int has_arc = 0; for(k = 0; k < n; k++) if(!S->up[k]) has_arc = 1;
	if(has_arc && (!gfs_hg_arc() || !hgm_building || !(hg_arc_R > 0))) return;
	double Rarc = hg_arc_R;
	double V = 0, Jx = 0, Jy = 0;
	for(k = 0; k < n; k++){
		int k1 = (k+1) % n;
		if(S->up[k]){ double cr = S->ax[k]*S->ay[k1] - S->ay[k]*S->ax[k1]; V += 0.5*cr; Jx += cr*(S->ax[k] + S->ax[k1])/6; Jy += cr*(S->ay[k] + S->ay[k1])/6; }
		else {
			double ta = atan2(S->ay[k], S->ax[k]), tb = atan2(S->ay[k1], S->ax[k1]);
			double dth = atan2(S->ax[k]*S->ay[k1] - S->ay[k]*S->ax[k1], S->ax[k]*S->ax[k1] + S->ay[k]*S->ay[k1]); if(dth <= 0) dth += 2*M_PI;
			V += 0.5*Rarc*Rarc*dth; Jx += Rarc*Rarc*Rarc/3*(sin(tb) - sin(ta)); Jy += Rarc*Rarc*Rarc/3*(cos(ta) - cos(tb));
		}
	}
	if(!(V > 0)) return;
	double cx = Jx/V, cy = Jy/V;
	/* support: index 0 = i, 1..n = neighbour of edge k (k+1) */
	double dc[SUB_MAXC + 1][4];
	for(q = 0; q <= n; q++){ dc[q][0] = dc[q][1] = dc[q][2] = dc[q][3] = 0; }
	static __thread int ac_q[SUB_MAXC][3], ac_n[SUB_MAXC]; static __thread double ac_M[SUB_MAXC][3][4];
	int arc_check = (!has_arc && hgm_building && getenv("GFS_HG_ARC_CHECK") != NULL);
	static __thread double dcg[SUB_MAXC + 1][4];
	if(has_arc || arc_check){
		double (*dcT)[4] = has_arc ? dc : dcg;
		if(arc_check) for(q = 0; q <= n; q++){ dcg[q][0] = dcg[q][1] = dcg[q][2] = dcg[q][3] = 0; }
		/* corner Jacobians da_k/dx_q (relative corners): Voronoi vertex (hg_vjac, J_i - I), bisector-circle point
		 * da/dx_j = Jn (Delta - R n)^T / (Jn.Delta), da/dx_i = -da/dx_j; arc-arc corners are rigid (0) */
		for(k = 0; k < n; k++){
			int kp = (k+n-1) % n; treevorork4particletype *Lo = S->up[kp], *Up = S->up[k];
			ac_n[k] = 0;
			if(Lo && Up){
				double p1x = (double)Lo->x - ix, p1y = (double)Lo->y - iy, p2x = (double)Up->x - ix, p2y = (double)Up->y - iy, Ji[4], J1[4], J2[4];
				if(!hg_vjac(p1x, p1y, p2x, p2y, S->ax[k], S->ay[k], Ji, J1, J2)){
#ifdef _OPENMP
#pragma omp atomic
#endif
					hg_arc_fail++;
					return;
				}
				ac_q[k][0] = 0; ac_M[k][0][0] = Ji[0] - 1; ac_M[k][0][1] = Ji[1]; ac_M[k][0][2] = Ji[2]; ac_M[k][0][3] = Ji[3] - 1;
				ac_q[k][1] = kp + 1; memcpy(ac_M[k][1], J1, sizeof(J1));
				ac_q[k][2] = k + 1; memcpy(ac_M[k][2], J2, sizeof(J2));
				ac_n[k] = 3;
			} else if(Lo || Up){
				treevorork4particletype *jb = Lo ? Lo : Up; int qj = Lo ? kp + 1 : k + 1;
				double Dx = (double)jb->x - ix, Dy = (double)jb->y - iy, ra = hypot(S->ax[k], S->ay[k]);
				double nx_ = S->ax[k]/ra, ny_ = S->ay[k]/ra, jnx = -ny_, jny = nx_, den = jnx*Dx + jny*Dy;
				if(!(fabs(den) > 0)){
#ifdef _OPENMP
#pragma omp atomic
#endif
					hg_arc_fail++;
					return;
				}
				double wx = (Dx - Rarc*nx_)/den, wy = (Dy - Rarc*ny_)/den;
				ac_q[k][0] = qj; ac_M[k][0][0] = jnx*wx; ac_M[k][0][1] = jnx*wy; ac_M[k][0][2] = jny*wx; ac_M[k][0][3] = jny*wy;
				ac_q[k][1] = 0; ac_M[k][1][0] = -jnx*wx; ac_M[k][1][1] = -jnx*wy; ac_M[k][1][2] = -jny*wx; ac_M[k][1][3] = -jny*wy;
				ac_n[k] = 2;
			}
		}
		/* dV/da_k, dM/da_k from the boundary pieces */
		static __thread double gV[SUB_MAXC][2], gM[SUB_MAXC][4];
		for(k = 0; k < n; k++){ gV[k][0] = gV[k][1] = 0; gM[k][0] = gM[k][1] = gM[k][2] = gM[k][3] = 0; }
		for(k = 0; k < n; k++){
			int k1 = (k+1) % n; double axk = S->ax[k], ayk = S->ay[k], bx_ = S->ax[k1], by_ = S->ay[k1];
			if(S->up[k]){
				double cr = axk*by_ - ayk*bx_, sx = axk + bx_, sy_ = ayk + by_;
				gV[k][0] += 0.5*by_; gV[k][1] += -0.5*bx_; gV[k1][0] += -0.5*ayk; gV[k1][1] += 0.5*axk;
				gM[k][0] += (cr + sx*by_)/6; gM[k][1] += (sx*(-bx_))/6; gM[k][2] += (sy_*by_)/6; gM[k][3] += (cr + sy_*(-bx_))/6;
				gM[k1][0] += (cr + sx*(-ayk))/6; gM[k1][1] += (sx*axk)/6; gM[k1][2] += (sy_*(-ayk))/6; gM[k1][3] += (cr + sy_*axk)/6;
			} else {
				double ra = hypot(axk, ayk), rb = hypot(bx_, by_), nax = axk/ra, nay = ayk/ra, nbx = bx_/rb, nby = by_/rb;
				double jax = -nay, jay = nax, jbx = -nby, jby = nbx, c3 = Rarc*Rarc/3;
				gV[k][0] += -0.5*Rarc*jax; gV[k][1] += -0.5*Rarc*jay; gV[k1][0] += 0.5*Rarc*jbx; gV[k1][1] += 0.5*Rarc*jby;
				gM[k][0] += -c3*nax*jax; gM[k][1] += -c3*nax*jay; gM[k][2] += -c3*nay*jax; gM[k][3] += -c3*nay*jay;
				gM[k1][0] += c3*nbx*jbx; gM[k1][1] += c3*nbx*jby; gM[k1][2] += c3*nby*jbx; gM[k1][3] += c3*nby*jby;
			}
		}
		for(k = 0; k < n; k++){
			double C00 = (gM[k][0] - cx*gV[k][0])/V, C01 = (gM[k][1] - cx*gV[k][1])/V, C10 = (gM[k][2] - cy*gV[k][0])/V, C11 = (gM[k][3] - cy*gV[k][1])/V;
			int t;
			for(t = 0; t < ac_n[k]; t++){
				double *M = ac_M[k][t], *D = dcT[ac_q[k][t]];
				D[0] += C00*M[0] + C01*M[2]; D[1] += C00*M[1] + C01*M[3];
				D[2] += C10*M[0] + C11*M[2]; D[3] += C10*M[1] + C11*M[3];
			}
		}
		/* face-mode convention: dc[0] holds d(absolute centroid)/dx_i = dc_rel/dx_i + I */
		dcT[0][0] += 1; dcT[0][3] += 1;
	}
	if(!has_arc)
	for(k = 0; k < n; k++){
		int kn = (k+1) % n, kp = (k+n-1) % n;
		double ax = S->ax[k], ay = S->ay[k], nx_ = S->ax[kn], ny_ = S->ay[kn], px = S->ax[kp], py = S->ay[kp];
		double dVx = 0.5*(ny_ - py), dVy = 0.5*(-(nx_ - px));
		double ck = ax*ny_ - ay*nx_, cp = px*ay - py*ax;
		double C00 = ((ax + nx_)*ny_ + ck + (px + ax)*(-py) + cp)/6.0, C01 = ((ax + nx_)*(-nx_) + (px + ax)*px)/6.0;
		double C10 = ((ay + ny_)*ny_ + (py + ay)*(-py))/6.0,             C11 = ((ay + ny_)*(-nx_) + ck + (py + ay)*px + cp)/6.0;
		C00 = (C00 - cx*dVx)/V; C01 = (C01 - cx*dVy)/V; C10 = (C10 - cy*dVx)/V; C11 = (C11 - cy*dVy)/V;
		treevorork4particletype *lo = S->up[kp], *up = S->up[k];
		double p1x = (double)lo->x - ix, p1y = (double)lo->y - iy, p2x = (double)up->x - ix, p2y = (double)up->y - iy;
		double Ji[4], J1[4], J2[4];
		if(!hg_vjac(p1x, p1y, p2x, p2y, ax, ay, Ji, J1, J2)) continue;
		/* dc/dx_b += Jc_k J_{k,b}; lo is the neighbour of edge kp (index kp+1), up of edge k (index k+1) */
		double *Jb[3] = {Ji, J1, J2}; int ib[3] = {0, kp + 1, k + 1}, t;
		for(t = 0; t < 3; t++){
			double *M = Jb[t], *D = dc[ib[t]];
			D[0] += C00*M[0] + C01*M[2]; D[1] += C00*M[1] + C01*M[3];
			D[2] += C10*M[0] + C11*M[2]; D[3] += C10*M[1] + C11*M[3];
		}
	}
	if(arc_check){
		double mx_ = 0, nr_ = 0;
		for(q = 0; q <= n; q++){ int t_; for(t_ = 0; t_ < 4; t_++){ double e_ = fabs(dcg[q][t_] - dc[q][t_]); if(e_ > mx_) mx_ = e_; if(fabs(dc[q][t_]) > nr_) nr_ = fabs(dc[q][t_]); } }
		if(nr_ > 0 && mx_/nr_ > hg_arc_chk){
#ifdef _OPENMP
#pragma omp critical (hg_arc_chk)
#endif
			{ if(mx_/nr_ > hg_arc_chk) hg_arc_chk = mx_/nr_; }
		}
	}
	/* 2-hop stencil: self, edge neighbours (w = A/d^2), their recorded neighbours (h_i w_ij w_jk) */
	int ns = 0; double sw[HG_MAXS], sy[HG_MAXS][2]; treevorork4particletype *sp[HG_MAXS]; long sid[HG_MAXS];
	double wsum = 0, hi = sqrt(V);
	#define HG_ADD(ID_, PT_, YX_, YY_, W_) do { int z_, f_ = -1; \
		for(z_ = 0; z_ < ns; z_++) if(sid[z_] == (ID_) && fabs(sy[z_][0] - (YX_)) + fabs(sy[z_][1] - (YY_)) < 1e-9*hi){ f_ = z_; break; } \
		if(f_ >= 0) sw[f_] += (W_); else if(ns < HG_MAXS){ sid[ns] = (ID_); sp[ns] = (PT_); sy[ns][0] = (YX_); sy[ns][1] = (YY_); sw[ns] = (W_); ns++; } } while(0)
	int nbase[SUB_MAXC];
	for(k = 0; k < n; k++){
		int k1 = (k+1) % n; treevorork4particletype *j = S->up[k];
		nbase[k] = 0;
		if(!j) continue;   /* arc edge (GFS_HG_ARC): no fluid neighbour */
		double yx = (double)j->x - ix, yy = (double)j->y - iy, A = hypot(S->ax[k1] - S->ax[k], S->ay[k1] - S->ay[k]);
		double w = A/(yx*yx + yy*yy); wsum += w;
		HG_ADD((long)PINDX(j), j, yx, yy, w);
		nbase[k] = 0;
		size_t jid = (size_t)PINDX(j);
		if(jid < shf_cap){
			int m2;
			for(m2 = 0; m2 < shf_n[jid]; m2++){
				long qid = shf_j[jid*SHF_K + m2]; double *g = shf_g + (jid*SHF_K + m2)*6;
				double y2x = yx + g[0], y2y = yy + g[1];
				if(fabs(y2x) + fabs(y2y) < 1e-12) continue;
				double Aq = hypot(g[4] - g[2], g[5] - g[3]), wq = Aq/(g[0]*g[0] + g[1]*g[1]);
				treevorork4particletype *qp = ((size_t)qid < hg_cap) ? hg_ptr[qid] : NULL;
				if(qp) HG_ADD(qid, qp, y2x, y2y, hi*w*wq);
			}
		}
	}
	HG_ADD((long)PINDX(ibp), ibp, 0.0, 0.0, wsum);
	#undef HG_ADD
	if(hgm_building && gfs_bprime() && !has_arc){
		size_t iid = (size_t)PINDX(ibp);
		if(iid < bp_cap){
			int nf = 0;
			for(k = 0; k < n && nf < BP_MAXF; k++){
				int k1 = (k+1) % n; treevorork4particletype *j = S->up[k];
				double Dx = (double)j->x - ix, Dy = (double)j->y - iy, d2 = Dx*Dx + Dy*Dy, d = sqrt(d2), u = 2*(cx*Dx + cy*Dy)/d2;
				double *F_ = bp_f + (iid*BP_MAXF + nf)*7;
				bp_fj[iid*BP_MAXF + nf] = (long)PINDX(j);
				F_[0] = Dx/d; F_[1] = Dy/d; F_[2] = d; F_[3] = hypot(S->ax[k1] - S->ax[k], S->ay[k1] - S->ay[k]);
				if(u < 1 && ibp->ie > 0){ double rf = rho*(1 - u); F_[4] = rf; F_[5] = gam*P*pow(1 - u, gam)/rf/(double)ibp->ie; F_[6] = 1; }
				else { F_[4] = F_[5] = 0; F_[6] = 0; }
				nf++;
			}
			bp_nf[iid] = nf;
			if(n > BP_MAXF){ fprintf(stderr, "[BPRIME] cell with %d > %d faces\n", n, BP_MAXF); exit(3); }
		}
	}
	/* weighted QR of A = W^1/2 [1, y_x, y_y] (modified Gram-Schmidt, twice); projection
	 * P = W^-1/2 Q Q^T W^1/2, l = g (I - P) (Astra audit 2) */
	static __thread double Qm[HG_MAXS][6];
	int ncolp = gfs_hg_quad() ? 6 : 3;
	double nrm0[6] = {0, 0, 0, 0, 0, 0}, Rd[6] = {0, 0, 0, 0, 0, 0}, Rf[6][6]; int rank = 0, col, rank3 = 0;
	memset(Rf, 0, sizeof(Rf));
	for(col = 0; col < ncolp; col++){
		double v_[HG_MAXS]; int z;
		for(z = 0; z < ns; z++){
			double sq = sqrt(sw[z]), Yx = sy[z][0]/hi, Yy = sy[z][1]/hi;
			v_[z] = (col == 0) ? sq : (col < 3) ? sq*sy[z][col-1] : (col == 3) ? sq*Yx*Yx : (col == 4) ? sq*Yx*Yy : sq*Yy*Yy;
			nrm0[col] += v_[z]*v_[z];
		}
		nrm0[col] = sqrt(nrm0[col]);
		int pass, kk;
		for(pass = 0; pass < 2; pass++) for(kk = 0; kk < rank; kk++){
			double dot = 0; for(z = 0; z < ns; z++) dot += Qm[z][kk]*v_[z];
			for(z = 0; z < ns; z++) v_[z] -= dot*Qm[z][kk];
			if(rank < ncolp) Rf[kk][rank] += dot;   /* column 'col' becomes Q column 'rank' if accepted */
		}
		double nn = 0; for(z = 0; z < ns; z++) nn += v_[z]*v_[z]; nn = sqrt(nn);
		if(nn > 64*2.2e-16*sqrt((double)ns)*nrm0[col]){ for(z = 0; z < ns; z++) Qm[z][rank] = v_[z]/nn; Rd[rank] = nn; Rf[rank][rank] = nn; rank++; if(col < 3) rank3++; }
		else { int kk3; for(kk3 = 0; kk3 < ncolp; kk3++) Rf[kk3][rank < ncolp ? rank : ncolp - 1] = 0; }
	}
	if(rank < ncolp && rank3 == 3){
#ifdef _OPENMP
#pragma omp atomic
#endif
		hg_nqdef++;   /* quadratic columns rank-deficient (diagnosed, projection uses the accepted columns) */
	}
	if(rank3 < 3){
#ifdef _OPENMP
#pragma omp atomic
#endif
		hg_npinv++;
	}
	double condR = (rank > 0) ? Rd[0]/Rd[rank-1] : 0;
	if(hgm_building && (gfs_bulk_affine() || gfs_bprime()) && !has_arc){
		if(rank3 == 3 && !gfs_hg_quad()){
			/* P = R^-1 (upper triangular); slope rows: beta_x,z = sqrt(w_z) sum_k P[1][k] Q[z][k], beta_y,z with P[2][k] */
			double Pm[3][3] = {{0,0,0},{0,0,0},{0,0,0}}; int a_, b_;
			for(a_ = 2; a_ >= 0; a_--){
				Pm[a_][a_] = 1/Rf[a_][a_];
				for(b_ = a_ + 1; b_ < 3; b_++){ double sum = 0; int c_; for(c_ = a_ + 1; c_ <= b_; c_++) sum += Rf[a_][c_]*Pm[c_][b_]; Pm[a_][b_] = -sum/Rf[a_][a_]; }
			}
			long ids_[HG_MAXS]; double mlx[HG_MAXS], mly[HG_MAXS]; int nm = 0, z1, z2;
			for(z1 = 0; z1 < ns; z1++){
				double sq = sqrt(sw[z1]), bx = 0, by = 0; int kk4;
				for(kk4 = 0; kk4 < 3; kk4++){ bx += Pm[1][kk4]*Qm[z1][kk4]; by += Pm[2][kk4]*Qm[z1][kk4]; }
				bx *= sq; by *= sq;
				int f2 = -1; for(z2 = 0; z2 < nm; z2++) if(ids_[z2] == sid[z1]){ f2 = z2; break; }
				if(f2 < 0){ ids_[nm] = sid[z1]; mlx[nm] = bx; mly[nm] = by; nm++; } else { mlx[f2] += bx; mly[f2] += by; }
			}
			if(gfs_bprime()){
				size_t iid = (size_t)PINDX(ibp); long ps0;
#ifdef _OPENMP
#pragma omp critical (bp_alloc)
#endif
				{ ps0 = bp_pn; if(bp_pn + ns <= bp_pcap) bp_pn += ns; else ps0 = -1; }
				if(ps0 < 0){ fprintf(stderr, "[BPRIME] slope pool overflow\n"); exit(3); }
				for(z1 = 0; z1 < ns; z1++){
					double sq = sqrt(sw[z1]), bx = 0, by = 0; int kk4;
					for(kk4 = 0; kk4 < 3; kk4++){ bx += Pm[1][kk4]*Qm[z1][kk4]; by += Pm[2][kk4]*Qm[z1][kk4]; }
					bp_sid[ps0 + z1] = sid[z1]; bp_s[4*(ps0 + z1)] = sy[z1][0]; bp_s[4*(ps0 + z1) + 1] = sy[z1][1];
					bp_s[4*(ps0 + z1) + 2] = sq*bx; bp_s[4*(ps0 + z1) + 3] = sq*by;
				}
				if(iid < bp_cap){ bp_sl_start[iid] = ps0; bp_sl_len[iid] = ns; }
			}
			if(gfs_bulk_affine()){
			if(getenv("GFS_BULK_TEST")){
				const double A_[4] = {0.37, -1.21, 0.83, 0.52};   /* v = b0 + A y, tr A = 0.89 */
				double dd = 0;
				for(z1 = 0; z1 < ns; z1++){
					double sq = sqrt(sw[z1]), bx = 0, by = 0; int kk4;
					for(kk4 = 0; kk4 < 3; kk4++){ bx += Pm[1][kk4]*Qm[z1][kk4]; by += Pm[2][kk4]*Qm[z1][kk4]; }
					dd += sq*bx*(0.3 + A_[0]*sy[z1][0] + A_[1]*sy[z1][1]) + sq*by*(-0.7 + A_[2]*sy[z1][0] + A_[3]*sy[z1][1]);
				}
				double er = fabs(dd - 0.89);
#ifdef _OPENMP
#pragma omp critical (bulk_test)
#endif
				{ if(er > bulk_test_err) bulk_test_err = er; }
			}
			long ms, os;
#ifdef _OPENMP
#pragma omp critical (hgm_alloc)
#endif
			{ ms = hgm_n; os = hgm_on; if(hgm_n < hgm_cap && hgm_on + nm <= hgm_ocap){ hgm_n++; hgm_on += nm; } else ms = -1; }
			if(ms < 0){ fprintf(stderr, "[BULK] mode store overflow\n"); exit(3); }
			double hcell = sqrt(V), mi = (double)ibp->mass;
			hgm[ms].owner = (long)PINDX(ibp); hgm[ms].key = (long)PINDX(ibp)*256 + 255; hgm[ms].start = os; hgm[ms].len = nm;
			hgm[ms].A0 = mi*hcell; hgm[ms].b = (gam + 1)/4*hcell; hgm[ms].kap = gam*(gam - 1)/mi; hgm[ms].onesided = 1; hgm[ms].kind = 1; hgm[ms].owner2 = -1;
			for(z1 = 0; z1 < nm; z1++){ hgm_id[os + z1] = ids_[z1]; hgm_l[2*(os + z1)] = mlx[z1]; hgm_l[2*(os + z1) + 1] = mly[z1]; }
#ifdef _OPENMP
#pragma omp atomic
#endif
			bulk_nmode++;
			}
		} else {
#ifdef _OPENMP
#pragma omp atomic
#endif
			bulk_nskip++;
		}
	}
	/* stencil index of each support member (i and edge neighbours) */
	int sidx[SUB_MAXC + 1];
	for(q = 0; q <= n; q++){
		treevorork4particletype *pt = (q == 0) ? ibp : S->up[q-1];
		if(!pt){ sidx[q] = -2; continue; }   /* arc edge: no support point, g = 0 there */
		double yx = (q == 0) ? 0 : (double)pt->x - ix, yy = (q == 0) ? 0 : (double)pt->y - iy; int z; sidx[q] = -1;
		for(z = 0; z < ns; z++) if(sid[z] == (long)PINDX(pt) && fabs(sy[z][0] - yx) + fabs(sy[z][1] - yy) < 1e-9*hi){ sidx[q] = z; break; }
	}
	double Hi = 0, lx[HG_MAXS], ly[HG_MAXS];
	for(k = 0; k < n; k++){
		treevorork4particletype *j = S->up[k]; int k1 = (k+1) % n;
		if(!j) continue;
		double Dx = (double)j->x - ix, Dy = (double)j->y - iy, d2 = Dx*Dx + Dy*Dy, d = sqrt(d2);
		double cD = cx*Dx + cy*Dy, u = 2*cD/d2;
		if(!(u < 1)) continue;
		double sc = -1.0/(1 - u);
		/* g_b over the support */
		double gx[SUB_MAXC + 1], gy[SUB_MAXC + 1];
		for(q = 0; q <= n; q++){
			double *D = dc[q];
			gx[q] = D[0]*2*Dx/d2 + D[2]*2*Dy/d2; gy[q] = D[1]*2*Dx/d2 + D[3]*2*Dy/d2;   /* (dc/dx_b)^T (2 Delta/d^2) */
		}
		gx[0] += -2*Dx/d2 - 2*cx/d2 + 4*cD*Dx/(d2*d2); gy[0] += -2*Dy/d2 - 2*cy/d2 + 4*cD*Dy/(d2*d2);
		gx[k+1] += 2*cx/d2 - 4*cD*Dx/(d2*d2);          gy[k+1] += 2*cy/d2 - 4*cD*Dy/(d2*d2);
		for(q = 0; q <= n; q++){ gx[q] *= sc; gy[q] *= sc; }
		/* l_b = g_b - sqrt(w_b) sum_k Q_bk t_k,  t_k = sum_a Q_ak g_a / sqrt(w_a) over the support */
		double tkx[6] = {0, 0, 0, 0, 0, 0}, tky[6] = {0, 0, 0, 0, 0, 0}; int kk2;
		for(q = 0; q <= n; q++){
			int z = sidx[q]; if(z < 0) continue;
			double iw = 1.0/sqrt(sw[z]);
			for(kk2 = 0; kk2 < rank; kk2++){ tkx[kk2] += Qm[z][kk2]*gx[q]*iw; tky[kk2] += Qm[z][kk2]*gy[q]*iw; }
		}
		double r = 0;
		for(q = 0; q < ns; q++){
			double sq = sqrt(sw[q]), ax_ = 0, ay_ = 0;
			for(kk2 = 0; kk2 < rank; kk2++){ ax_ += Qm[q][kk2]*tkx[kk2]; ay_ += Qm[q][kk2]*tky[kk2]; }
			lx[q] = -sq*ax_; ly[q] = -sq*ay_;
		}
		int ok = 1;
		for(q = 0; q <= n; q++){ if(sidx[q] == -2) continue; if(sidx[q] < 0){ ok = 0; break; } lx[sidx[q]] += gx[q]; ly[sidx[q]] += gy[q]; }
		if(!ok){
			if(hgm_building){
#ifdef _OPENMP
#pragma omp atomic
#endif
				hgm_nfail++;
			}
			continue;
		}
		for(q = 0; q < ns; q++) r += lx[q]*(double)sp[q]->vx + ly[q]*(double)sp[q]->vy;
		if(hgm_building && gfs_hg_quad() && getenv("GFS_HG_QUAD_TEST")){
			double rq = 0, ra = 0;
			for(q = 0; q < ns; q++){
				double Yx = sy[q][0]/hi, Yy = sy[q][1]/hi;
				double ux = 0.3 + 0.7*Yx - 1.1*Yy + 0.9*Yx*Yx - 0.4*Yx*Yy + 1.3*Yy*Yy, uy = -0.2 + 0.5*Yx + 0.6*Yy - 0.8*Yx*Yx + 1.2*Yx*Yy + 0.3*Yy*Yy;
				rq += lx[q]*ux + ly[q]*uy; ra += fabs(lx[q]*ux) + fabs(ly[q]*uy);
			}
			if(ra > 0 && rank == 6 && fabs(rq)/ra > hg_qtest){
#ifdef _OPENMP
#pragma omp critical (hg_qtest)
#endif
				{ if(fabs(rq)/ra > hg_qtest) hg_qtest = fabs(rq)/ra; }
			}
		}
		double A = hypot(S->ax[k1] - S->ax[k], S->ay[k1] - S->ay[k]);
		double hh = 0.5*d, Wf = 0.25*A*d, rf = rho*(1 - u), cf = sqrt(gam*P*pow(1 - u, gam)/rf), bb = (gam + 1)/4*hh;
		if(hgm_building && gfs_hg_off()) continue;
		if(hgm_building){
			/* merge periodic images by global id, store the mode */
			long ids_[HG_MAXS]; double mlx[HG_MAXS], mly[HG_MAXS]; int nm = 0, z1, z2;
			for(z1 = 0; z1 < ns; z1++){
				int f2 = -1; for(z2 = 0; z2 < nm; z2++) if(ids_[z2] == sid[z1]){ f2 = z2; break; }
				if(f2 < 0){ ids_[nm] = sid[z1]; mlx[nm] = lx[z1]; mly[nm] = ly[z1]; nm++; } else { mlx[f2] += lx[z1]; mly[f2] += ly[z1]; }
			}
			long ms, os;
#ifdef _OPENMP
#pragma omp critical (hgm_alloc)
#endif
			{ ms = hgm_n; os = hgm_on; if(hgm_n < hgm_cap && hgm_on + nm <= hgm_ocap){ hgm_n++; hgm_on += nm; } else ms = -1; }
			if(ms < 0){ fprintf(stderr, "[HGSPLIT] mode store overflow\n"); exit(3); }
			hgm[ms].owner = (long)PINDX(ibp); hgm[ms].key = (long)PINDX(ibp)*256 + k; hgm[ms].start = os; hgm[ms].len = nm;
			hgm[ms].A0 = rf*hh*Wf; hgm[ms].b = bb; hgm[ms].kap = gam*(gam - 1)*pow(1 - u, gam - 1)/(double)ibp->mass; hgm[ms].onesided = 0; hgm[ms].kind = 0; hgm[ms].owner2 = -1;
			for(z1 = 0; z1 < nm; z1++){ hgm_id[os + z1] = ids_[z1]; hgm_l[2*(os + z1)] = mlx[z1]; hgm_l[2*(os + z1) + 1] = mly[z1]; }
			continue;
		}
		double sq_ = sqrt(bb*bb*r*r + cf*cf);
		double C = rf*hh*Wf*(bb*fabs(r) + sq_);
		double Cp = C + fabs(r)*rf*hh*Wf*(bb + bb*bb*fabs(r)/sq_);   /* C + q C' (linearized damping) */
		Hi += C*r*r;
		{
			double Ssum = 0;
			double lmx = 0, zl2 = 0, zg2 = 0;
			for(q = 0; q < ns; q++){ double mq = (double)sp[q]->mass, l2 = lx[q]*lx[q] + ly[q]*ly[q]; Ssum += sqrt(l2/mq); zl2 += l2/mq; if(sqrt(l2) > lmx) lmx = sqrt(l2); }
			for(q = 0; q <= n; q++){ treevorork4particletype *pt = (q == 0) ? ibp : S->up[q-1]; zg2 += (gx[q]*gx[q] + gy[q]*gy[q])/(double)pt->mass; }
			double lamf = Cp*zl2;   /* exact nonzero eigenvalue of this face alone */
			if(lamf > hg_worst[0]){
#ifdef _OPENMP
#pragma omp critical (hg_worst)
#endif
				{ if(lamf > hg_worst[0]){ hg_worst[0] = lamf; hg_worst[1] = d/sqrt(V); hg_worst[2] = 1 - u; hg_worst[3] = sqrt(zg2*(double)ibp->mass);
				    hg_worst[4] = sqrt(zl2*(double)ibp->mass); hg_worst[5] = condR; hg_worst[6] = fabs(r)*hh/cf; hg_worst[7] = ns; hg_worst[8] = (double)PINDX(ibp); hg_worst[9] = C/(rf*cf*hh*Wf); } }
			}
			if(gfs_hg_spec()){
				long fs, os;
#ifdef _OPENMP
#pragma omp critical (hgo_alloc)
#endif
				{ fs = hgf_n; os = hgo_n; if(hgf_n < hgf_cap && hgo_n + ns <= hgo_cap){ hgf_n++; hgo_n += ns; } else fs = -1; }
				if(fs < 0){ fprintf(stderr, "[HG] face-operator store overflow\n"); exit(3); }
				hgf_start[fs] = os; hgf_k[fs] = Cp;
				for(q = 0; q < ns; q++){ hgo_id[os + q] = sid[q]; hgo_l[2*(os + q)] = lx[q]; hgo_l[2*(os + q) + 1] = ly[q]; }
				hgf_len[fs] = ns;
			}
			for(q = 0; q < ns; q++){
				size_t id = (size_t)sid[q];
				if(id >= hg_cap) continue;
				double add = Cp*hypot(lx[q], ly[q])/sqrt((double)sp[q]->mass)*Ssum;
#ifdef _OPENMP
#pragma omp atomic
#endif
				hg_lam[id] += add;
			}
		}
		for(q = 0; q < ns; q++){
			size_t id = (size_t)sid[q];
			if(id >= hg_cap) continue;
			double fx_ = -C*r*lx[q], fy_ = -C*r*ly[q];
#ifdef _OPENMP
#pragma omp atomic
#endif
			hg_fx[id] += fx_;
#ifdef _OPENMP
#pragma omp atomic
#endif
			hg_fy[id] += fy_;
		}
	}
	if(has_arc && hgm_building && !gfs_hg_off()){
		/* arc-component modes: xi = log(B(c)/W), W = R^2 dth/2, N = R (sin te - sin ts, cos ts - cos te),
		 * B = W - c.N/2 (relative coordinates); g_q from the corner Jacobians of the two end corners and dc */
		int k0 = -1;
		for(k = 0; k < n; k++) if(S->up[k] && !S->up[(k+1) % n]){ k0 = (k+1) % n; break; }   /* first arc edge after a fluid edge */
		if(k0 >= 0){
			int comp = 0, kk = k0, steps = 0;
			while(steps < n){
				if(S->up[kk]){ kk = (kk+1) % n; steps++; continue; }
				int ks = kk; double dth = 0;
				while(!S->up[kk] && steps < n){
					int k1 = (kk+1) % n;
					double d_ = atan2(S->ax[kk]*S->ay[k1] - S->ay[kk]*S->ax[k1], S->ax[kk]*S->ax[k1] + S->ay[kk]*S->ay[k1]); if(d_ <= 0) d_ += 2*M_PI;
					dth += d_; kk = k1; steps++;
				}
				int ce = kk, cs = ks;   /* start corner cs, end corner ce */
				double ts = atan2(S->ay[cs], S->ax[cs]), te = atan2(S->ay[ce], S->ax[ce]);
				double W = 0.5*Rarc*Rarc*dth, Nx = Rarc*(sin(te) - sin(ts)), Ny = Rarc*(cos(ts) - cos(te));
				double B = W - 0.5*(cx*Nx + cy*Ny), eta = B/W;
				if(!(B > 0) || !(W > 0)){
#ifdef _OPENMP
#pragma omp atomic
#endif
					hg_arc_fail++;
					comp++; continue;
				}
				double gx[SUB_MAXC + 1], gy[SUB_MAXC + 1];
				for(q = 0; q <= n; q++){
					/* dth_e/dx_q, dth_s/dx_q (2-vectors) */
					double dte[2] = {0, 0}, dts[2] = {0, 0}; int t;
					double jex = -sin(te), jey = cos(te), jsx = -sin(ts), jsy = cos(ts);
					for(t = 0; t < ac_n[ce]; t++) if(ac_q[ce][t] == q){ double *M = ac_M[ce][t]; dte[0] += (jex*M[0] + jey*M[2])/Rarc; dte[1] += (jex*M[1] + jey*M[3])/Rarc; }
					for(t = 0; t < ac_n[cs]; t++) if(ac_q[cs][t] == q){ double *M = ac_M[cs][t]; dts[0] += (jsx*M[0] + jsy*M[2])/Rarc; dts[1] += (jsx*M[1] + jsy*M[3])/Rarc; }
					double dW[2] = {0.5*Rarc*Rarc*(dte[0] - dts[0]), 0.5*Rarc*Rarc*(dte[1] - dts[1])};
					/* dN (row: N comp, col: x comp) */
					double dN00 = Rarc*(cos(te)*dte[0] - cos(ts)*dts[0]), dN01 = Rarc*(cos(te)*dte[1] - cos(ts)*dts[1]);
					double dN10 = Rarc*(sin(te)*dte[0] - sin(ts)*dts[0]), dN11 = Rarc*(sin(te)*dte[1] - sin(ts)*dts[1]);
					/* relative-centroid Jacobian */
					double D0 = dc[q][0] - (q == 0 ? 1 : 0), D1 = dc[q][1], D2 = dc[q][2], D3 = dc[q][3] - (q == 0 ? 1 : 0);
					double dBx = dW[0] - 0.5*(Nx*D0 + Ny*D2 + cx*dN00 + cy*dN10);
					double dBy = dW[1] - 0.5*(Nx*D1 + Ny*D3 + cx*dN01 + cy*dN11);
					gx[q] = dBx/B - dW[0]/W; gy[q] = dBy/B - dW[1]/W;
				}
				{
					const char *df_ = getenv("GFS_HG_ARC_DUMP");
					static long ndump_ = 0;
					if(df_ && ndump_ < 400){
#ifdef _OPENMP
#pragma omp critical (hg_arc_dump)
#endif
						{
							FILE *fa_ = fopen(df_, "a");
							if(fa_){
								fprintf(fa_, "CELL %ld R %.17g n %d xi %.17g\n", (long)PINDX(ibp), Rarc, n, log(eta));
								for(q = 0; q < n; q++){ treevorork4particletype *u_ = S->up[q]; fprintf(fa_, "E %d %d %.17g %.17g\n", q, u_ ? 1 : 0, u_ ? (double)u_->x - ix : 0.0, u_ ? (double)u_->y - iy : 0.0); }
								for(q = 0; q <= n; q++) fprintf(fa_, "G %d %.17g %.17g\n", q, gx[q], gy[q]);
								fclose(fa_); ndump_++;
							}
						}
					}
				}
				double tkx[6] = {0, 0, 0, 0, 0, 0}, tky[6] = {0, 0, 0, 0, 0, 0}; int kk2;
				for(q = 0; q <= n; q++){
					int z = sidx[q]; if(z < 0) continue;
					double iw = 1.0/sqrt(sw[z]);
					for(kk2 = 0; kk2 < rank; kk2++){ tkx[kk2] += Qm[z][kk2]*gx[q]*iw; tky[kk2] += Qm[z][kk2]*gy[q]*iw; }
				}
				for(q = 0; q < ns; q++){
					double sq = sqrt(sw[q]), ax_ = 0, ay_ = 0;
					for(kk2 = 0; kk2 < rank; kk2++){ ax_ += Qm[q][kk2]*tkx[kk2]; ay_ += Qm[q][kk2]*tky[kk2]; }
					lx[q] = -sq*ax_; ly[q] = -sq*ay_;
				}
				int ok = 1;
				for(q = 0; q <= n; q++){ if(sidx[q] == -2) continue; if(sidx[q] < 0){ ok = 0; break; } lx[sidx[q]] += gx[q]; ly[sidx[q]] += gy[q]; }
				if(!ok){
#ifdef _OPENMP
#pragma omp atomic
#endif
					hgm_nfail++;
					comp++; continue;
				}
				long ids_[HG_MAXS]; double mlx[HG_MAXS], mly[HG_MAXS]; int nm = 0, z1, z2;
				for(z1 = 0; z1 < ns; z1++){
					int f2 = -1; for(z2 = 0; z2 < nm; z2++) if(ids_[z2] == sid[z1]){ f2 = z2; break; }
					if(f2 < 0){ ids_[nm] = sid[z1]; mlx[nm] = lx[z1]; mly[nm] = ly[z1]; nm++; } else { mlx[f2] += lx[z1]; mly[f2] += ly[z1]; }
				}
				long ms, os;
#ifdef _OPENMP
#pragma omp critical (hgm_alloc)
#endif
				{ ms = hgm_n; os = hgm_on; if(hgm_n < hgm_cap && hgm_on + nm <= hgm_ocap){ hgm_n++; hgm_on += nm; } else ms = -1; }
				if(ms < 0){ fprintf(stderr, "[HGARC] mode store overflow\n"); exit(3); }
				double rf = rho*eta;
				hgm[ms].owner = (long)PINDX(ibp); hgm[ms].key = (long)PINDX(ibp)*256 + 200 + (comp < 50 ? comp : 50); hgm[ms].start = os; hgm[ms].len = nm;
				hgm[ms].A0 = rf*Rarc*W; hgm[ms].b = (gam + 1)/4*Rarc; hgm[ms].kap = gam*(gam - 1)*pow(eta, gam - 1)/(double)ibp->mass;
				hgm[ms].onesided = 0; hgm[ms].kind = 0; hgm[ms].owner2 = -1;
				for(z1 = 0; z1 < nm; z1++){ hgm_id[os + z1] = ids_[z1]; hgm_l[2*(os + z1)] = mlx[z1]; hgm_l[2*(os + z1) + 1] = mly[z1]; }
#ifdef _OPENMP
#pragma omp atomic
#endif
				hg_arc_nmode++;
				comp++;
			}
		}
#ifdef _OPENMP
#pragma omp atomic
#endif
		hg_arc_ncell++;
	}
	{ size_t idi = (size_t)PINDX(ibp); if(idi < hg_cap) hg_H[idi] = Hi; }
#ifdef _OPENMP
#pragma omp atomic
#endif
	hg_ncell++;
}
static long tcf_faces = 0, tcf_zneg = 0, tcf_mm0 = 0, tcf_cellneg = 0, tcf_small = 0;
static double tcf_heat = 0, tcf_zmin = 0;
static double hq_dpos = 0, hq_dneg = 0, hq_d1p = 0, hq_d1n = 0, hq_d2p = 0, hq_d2n = 0;   /* sum of D_f > 0 and of -D_f for D_f < 0 */
static int gfs_schemeH_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_SCHEME_H"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static double *dyn_wd = NULL, *dyn_k[4] = {NULL,NULL,NULL,NULL};
static size_t dyn_cap = 0;
static double dyn_fchk_s2 = 0; static long dyn_fchk_n = 0;
static int gfs_dynw2_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_DYN_W2"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static double gfs_env_d(const char *k, double def){
	const char *s = getenv(k); return (s && s[0]) ? atof(s) : def;
}
static void dyn_ensure(size_t max_id){
	size_t k;
	if(max_id < dyn_cap) return;
	size_t nc = max_id + 1 + max_id/4;
	dyn_wd = (double*)realloc(dyn_wd, nc*sizeof(double));
	for(k=dyn_cap;k<nc;k++) dyn_wd[k] = 0;
	int s;
	for(s=0;s<4;s++){
		dyn_k[s] = (double*)realloc(dyn_k[s], nc*sizeof(double));
		for(k=dyn_cap;k<nc;k++) dyn_k[s][k] = 0;
	}
	dyn_cap = nc;
}
static inline double dyn_wdot(treevorork4particletype *p){
	size_t id = (size_t)PINDX(p);
	return (PINDX(p) != MAX_INDEX && id < dyn_cap) ? dyn_wd[id] : 0;
}
/* GFS_VOL_AUDIT=1 (default 0, diagnostic, 1 rank): per step, compare the
 * geometric volume change of each cell (tessellation at the start of step
 * n+1 minus that at the start of step n, so end-of-step weight jumps are
 * included) with the volume change booked by the face velocities used in
 * the energy equation (RK4-weighted sum of face u·dS). Prints [VAUD] with
 * sum|res|/sumV and the unbooked work sum P|res| per step and cumulative. */
static double gfs_gcl_diag_de = 0;
static int gfs_gcl_diag_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_GCL_DIAG"); on = (s && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static int gfs_vaud_on(void){
	static int on = -1;
	if(on < 0){
		const char *s = getenv("GFS_VOL_AUDIT");
		on = (s && s[0] && atoi(s) > 0) ? 1 : 0;
	}
	return on;
}
static size_t vaud_cap = 0;
static double *vaud_rate = NULL, *vaud_book = NULL, *vaud_vprev = NULL;
static char *vaud_have = NULL;
/* heating decomposition for cells without vacuum arcs (vaud_arc = 0):
 * vaud_diss = RK sum of (dI/dt + P dV_booked/dt) (dissipation), vaud_adb =
 * RK sum of -P dV_booked/dt, vaud_ie0 / vaud_p0 = ie, P at the last close */
static double *vaud_diss = NULL, *vaud_adb = NULL, *vaud_ie0 = NULL, *vaud_p0 = NULL;
static char *vaud_arc = NULL, *vaud_arcp = NULL;
static void vaud_ensure(size_t id){
	if(id < vaud_cap) return;
	size_t nc = id + id/4 + 64, k;
	vaud_rate  = (double*)realloc(vaud_rate,  nc*sizeof(double));
	vaud_book  = (double*)realloc(vaud_book,  nc*sizeof(double));
	vaud_vprev = (double*)realloc(vaud_vprev, nc*sizeof(double));
	vaud_have  = (char*)realloc(vaud_have, nc);
	vaud_diss = (double*)realloc(vaud_diss, nc*sizeof(double));
	vaud_adb  = (double*)realloc(vaud_adb,  nc*sizeof(double));
	vaud_ie0  = (double*)realloc(vaud_ie0,  nc*sizeof(double));
	vaud_p0   = (double*)realloc(vaud_p0,   nc*sizeof(double));
	vaud_arc  = (char*)realloc(vaud_arc, nc);
	vaud_arcp = (char*)realloc(vaud_arcp, nc);
	for(k=vaud_cap;k<nc;k++){ vaud_rate[k] = vaud_book[k] = vaud_vprev[k] = 0; vaud_have[k] = 0;
		vaud_diss[k] = vaud_adb[k] = vaud_ie0[k] = vaud_p0[k] = 0; vaud_arc[k] = vaud_arcp[k] = 0; }
	vaud_cap = nc;
}

#ifdef USE_CUDA
#include "exam_gpu.h"
/* Forward declaration — full signature matches exam_gpu_extract.c */
double getAccVoro2DBlend_GPU(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    postype OrderOfAccuracy, postype Courant, postype Gamma,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
    treevorork4particletype *(*find2DCellBP)(SimParameters *, int, int, int *),
    void (*mkLinkedList2D)(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)));
double getAccVoro2DBlend_GPU_validate(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    postype OrderOfAccuracy, postype Courant, postype Gamma,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
    treevorork4particletype *(*find2DCellBP)(SimParameters *, int, int, int *),
    void (*mkLinkedList2D)(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)));
/* GPU nearest-neighbor (replaces det2d_dpqRK4 for Voronoi av_modes) */
void det2d_dpqRK4_GPU(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    void mkLinkedList2D(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)));
/* LagMFM GPU wrappers (exam_gpu_extract.c) */
void updateDenW2Pressure2D_LagMFM_GPU(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    postype Gamma,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    void mkLinkedList2D(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)),
    postype Dtime);
double getAccVoro2D_LagMFM_GPU(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    postype OrderOfAccuracy, postype Courant, postype Gamma,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    void mkLinkedList2D(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)));
#endif


inline postype getw2forHydroParticle(SimParameters *simpar, treevorork4particletype *bp,
		postype Dtime){
	/* Kappa = 0 is pure Voronoi. A NaN from pow() must not leak into w2,
	 * or the neighbor clip turns it into d_nn^2 and the face leaves the midpoint. */
	if(GAS_Kappa(simpar) == 0) return (postype)0;
//	postype w1old = sqrt(bp->w2old);
	postype gamma = GAS_GAMMA(simpar);
	int ndim = NDIM(simpar);
	if(1){
		postype base = GAS_Kappa(simpar)*GAS_Kappa(simpar)*
			GAS_dMean(simpar)*GAS_dMean(simpar);
		postype q;
		int mode = GAS_W2MODE(simpar);
		if(mode == 1){
			/* Entropy-A power-law: w2 ∝ A^w2Power, A = P/rho^gamma. */
			postype A = bp->pressure / pow((postype)bp->den, gamma);
			q = pow(A * GAS_invw2Scale(simpar), GAS_w2Power(simpar));
		} else if(mode == 2){
			/* Specific-entropy additive: w2 = base * (1 + beta*(s - s_ref)/s_scale).
			 * s = ln(P/rho^gamma). Bounded perturbation, no pow. */
			postype s_local = log(bp->pressure / pow((postype)bp->den, gamma));
			postype s_scale = GAS_W2SSCALE(simpar);
			if(s_scale <= 0) s_scale = 1.0;  /* safety */
			q = 1.0 + GAS_W2SBETA(simpar) * (s_local - GAS_W2SREF(simpar)) / s_scale;
			if(q < 0.1) q = 0.1;  /* clamp to avoid pathological w2 near zero */
		} else if(mode == 3){
			/* Pressure + volume-equalization hybrid (Arepo-Laguerre).
			 * q = q_P · q_V where q_V drives V_i toward Vmean = dMean^ndim.
			 * γ = 0 recovers mode 0 exactly. Multiplicative form keeps the
			 * pressure-tracking baseline; q_V is a bounded mesh-quality knob. */
			postype q_P = pow(bp->pressure * GAS_invw2Scale(simpar), GAS_w2Power(simpar));
			postype Vmean = pow((postype)GAS_dMean(simpar), (postype)ndim);
			postype Vi = (postype)bp->volume;
			postype q_V = 1.0;
			if(Vmean > 0 && Vi > 0){
				q_V = 1.0 + GAS_W2VOLGAMMA(simpar) * (Vmean - Vi) / Vmean;
				if(q_V < 0.1) q_V = 0.1;   /* lower clamp */
				if(q_V > 10.0) q_V = 10.0; /* upper clamp */
			}
			q = q_P * q_V;
		} else if(mode == 4 || mode == 5 || mode == 6){
			/* cs-based adaptive Laguerre weight (per-particle).
			 * Asymmetric response: cells with high cs (active/post-shock) grow
			 * their Laguerre weight, cells with low cs (ambient) stay near-Voronoi.
			 * Cascading propagation: shock arrival raises cs, raises w2,
			 * shifts Laguerre boundary toward downstream neighbor.
			 * Mode 4 (A): linear         length = α·cs·dMean
			 * Mode 5 (B): Padé saturate   length = α·cs·dMean / (1 + α·cs/β)  (→ β·dMean)
			 * Mode 6 (C): power-law       length = α·cs^p·dMean
			 * Note: 'base' (Kappa²·dMean²) is replaced — these formulas
			 * compute length² directly. We bypass Kappa here. */
			postype alpha = (postype)GAS_W2ALPHA(simpar);
			postype cs    = (postype)bp->csound;
			postype dM    = (postype)GAS_dMean(simpar);
			postype length;
			if(mode == 4){
				length = alpha * cs * dM;
			} else if(mode == 5){
				postype beta = (postype)GAS_W2BETA(simpar);
				postype f    = alpha * cs;
				if(beta <= 0) beta = 0.5;  /* safety */
				length = (f * dM) / (1.0 + f / beta);
			} else { /* mode == 6 */
				postype p = (postype)GAS_W2CSEXP(simpar);
				if(p <= 0) p = 1.0;
				postype cs_pow = (cs > 0) ? pow(cs, p) : 0.0;
				length = alpha * cs_pow * dM;
			}
			return length * length;
		} else if(mode == 7){
			/* dK-based source-driven Laguerre weight (mode D).
			 * Heating-expansion equilibrium: cells with positive entropy source
			 * (NS-stress + conduction inflow > 0) expand their Laguerre weight,
			 * pushing boundary outward → V↑ → ρ↓ → P=K·ρ^γ ↓ → relief.
			 *   length² = (α·dM)² · (1 + β·max(0,η))
			 *   η = max(0, dK/K) · (dM/cs)        [relative heating per sound crossing]
			 * Asymmetric: cooling (dK<0) gives baseline (α·dM)²; only heating expands.
			 * Requires av_mode>=1 entropy_mode=1 (stress fields). Otherwise baseline only.
			 */
			postype alpha = (postype)GAS_W2ALPHA(simpar);
			postype dM    = (postype)GAS_dMean(simpar);
			postype length0 = alpha * dM;
			postype eta = 0.0;
			postype beta = 0.0;
			if(GAS_AVMODE(simpar) >= 1 && GAS_ENTROPY_MODE(simpar) == 1){
				treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)bp;
				beta            = (postype)GAS_W2DKBETA(simpar);
				postype etamx   = (postype)GAS_W2DKETAMAX(simpar);
				postype dKv     = (postype)sbp->stress.dK;
				postype K       = (postype)sbp->stress.K;
				postype K_floor = (postype)GAS_K_FLOOR(simpar);
				postype cs      = (postype)bp->csound;
				if(dKv > 0 && cs > 0){
					postype Ksafe = (K > K_floor) ? K : K_floor;
					eta = (dKv / Ksafe) * (dM / cs);
					if(etamx > 0 && eta > etamx) eta = etamx;
				}
			}
			return length0 * length0 * (1.0 + beta * eta);
		} else if(mode == 8){
			/* Local pressure contrast, bounded (default off; w2_mode 8).
			 *   x = P_i / Pbar_i, Pbar_i = face-length-weighted mean of the
			 *   neighbour pressures (avgNeighboringPressure, last tessellation),
			 *   q = x^a / (1 + x^a), a = w2Power, so 0 < q < 1.
			 * Uniform pressure gives q = 1/2 on every cell: W is uniform and
			 * the mesh is exactly Voronoi, with no dependence on a global
			 * pressure scale. Across a jump the high-pressure side tends to
			 * q = 1 and the low side to 0, so |W_p - W_q| < base = kappa^2 dMean^2
			 * and the face moves toward the low-pressure (upstream) side, as
			 * the HLL shock speed does. */
			postype Pn = bp->avgNeighboringPressure;
			postype x = (Pn > 0 && bp->pressure > 0) ? bp->pressure / Pn : (postype)1;
			postype xa = pow(x, GAS_w2Power(simpar));
			q = xa / (1 + xa);
			if(!(q >= 0)) q = 0.5;   /* NaN guard */
		} else {
			/* Mode 0 (default): pressure-based power-law. */
			q = pow(bp->pressure * GAS_invw2Scale(simpar), GAS_w2Power(simpar));
		}
		return base * q;
	}
	else {
		postype w1 = sqrt(bp->w2);
		postype exponent = 1./(ndim*gamma);
		postype dw = GAS_w2Power(simpar)*(pow(bp->pressure/bp->avgNeighboringPressure, exponent)-1)*w1;

		postype shift = bp->csound*Dtime;
		if(PINDX(bp)==16384) 
			DEBUGPRINT("P%d has w1= %g, exponent= %g, dw= %g, shift= %g, Cs= %g Dtime= %g\n", 
					MYID(simpar), w1, exponent, dw, shift, bp->csound, Dtime); 
		if(dw > 0) dw = MIN(dw, shift);
		else dw = MAX(dw, -shift);
		w1 = w1+dw;
		postype w2 = w1*w1;
		return w2;
	}
}

/* Apply w2 control mechanisms: floor, relaxation, rate limiter, ceiling.
 * In 2D, w2ceil is a hard geometric constraint (cell degenerates if exceeded),
 * so ceiling MUST be applied last — rate limiter must not push w2 above w2ceil.
 * (This differs from 1D where ceiling is soft and rate limiter comes after.)
 * Call AFTER getw2forHydroParticle sets bp->w2hydro. */
static inline void applyW2Controls(SimParameters *simpar, treevorork4particletype *bp, postype Dtime){
	postype w2 = bp->w2hydro;

	/* 1. Floor: w >= w2_floor_frac * dMean */
	postype w2ff = GAS_W2FLOORFRAC(simpar);
	if(w2ff > 0){
		postype wf = w2ff * GAS_dMean(simpar);
		postype w2floor = wf * wf;
		if(w2 < w2floor) w2 = w2floor;
	}

	/* 2. Relaxation: exponential smoothing toward target */
	postype w2rt = GAS_W2RELAXTAU(simpar);
	if(w2rt > 0 && Dtime > 0){
		postype ci = fmax(bp->csound, 1e-10);
		postype tau = w2rt * GAS_dMean(simpar) / ci;
		postype alpha_r = fmin(1.0, Dtime / tau);
		postype w2_old = (bp->w2old > 0) ? bp->w2old : w2;
		w2 = w2_old + alpha_r * (w2 - w2_old);
	}

	/* 3. Rate limiter: |Δw2/w2_old| ≤ w2_rate_max */
	postype w2rm = GAS_W2RATEMAX(simpar);
	if(w2rm > 0 && bp->w2old > 0){
		postype w2_hi = bp->w2old * (1.0 + w2rm);
		postype w2_lo = bp->w2old * (1.0 - w2rm);
		if(w2_lo < 0) w2_lo = 0;
		if(w2 > w2_hi) w2 = w2_hi;
		if(w2 < w2_lo) w2 = w2_lo;
	}

	bp->w2hydro = w2;

	/* 4. Ceiling: w2 ≤ w2ceil (MUST be last — hard geometric constraint in 2D) */
	bp->w2 = MIN(w2, bp->w2ceil);
}

// *_w2Masure2D are obsolete and instead use getw2forHydroParticle.
double kh_w2Measure2D(SimParameters *simpar, double dmean, double pressure,double rho){
#define p0_w2 2.5
#define rho0_w2 1.
    double res = GAS_Kappa(simpar)*dmean*pow(pressure/p0_w2*rho0_w2/rho,0.5);
    double res2 = res*res;
    return res2;
#undef p0_w2
#undef rho0_w2
}

double kp_w2Measure2D(SimParameters *simpar, double dmean, double pressure,double rho){
#define p0_w2 2.5
#define rho0_w2 1.
    double res = GAS_Kappa(simpar)*dmean*pow(pressure/p0_w2*rho0_w2/rho,0.5);
    double res2 = res*res;
    return res2;
#undef p0_w2
#undef rho0_w2
}

double rt_w2Measure2D(SimParameters *simpar, double dmean, double pressure,double rho){
#define p0_w2 2.5
#define rho0_w2 1.
    double res = GAS_Kappa(simpar)*dmean*pow(pressure/p0_w2*rho0_w2/rho,0.5);
    double res2 = res*res;
    return res2;
#undef p0_w2
#undef rho0_w2
}

double cyl_w2Measure2D(SimParameters *simpar, double dmean, double pressure,double rho){
#define p0_w2 71.42857
#define rho0_w2 1.
    double res = GAS_Kappa(simpar)*dmean*pow(pressure/p0_w2*rho0_w2/rho,0.5);
    double res2 = res*res;
    return res2;
#undef p0_w2
#undef rho0_w2
}

double gl2d_w2Measure(SimParameters *simpar, double dmean, double pressure,double rho){
#define p0_w2 1.
#define rho0_w2 1.
    double res = GAS_Kappa(simpar)*dmean*pow(pressure/p0_w2*rho0_w2/rho,0.5);
    double res2 = res*res;
    return res2;
#undef p0_w2
#undef rho0_w2
}




void ex2d_findCentroid(TStruct *thisNode){
    TPtlStruct *nodeparticles = (TPtlStruct*)thisNode->daughter;
    thisNode->monox = thisNode->monoy = 0;
    thisNode->mass = 0;
    TPtlStruct *p2ptl;
    for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling){
		postype mass = p2ptl->mass;
        thisNode->mass += mass;
        thisNode->monox += mass* p2ptl->x;
        thisNode->monoy += mass* p2ptl->y;
    }
    thisNode->monox = thisNode->monox/thisNode->mass;
    thisNode->monoy = thisNode->monoy/thisNode->mass;
}

void ex3d_findCentroid(TStruct *thisNode){
    TPtlStruct *nodeparticles = (TPtlStruct*)thisNode->daughter;
    thisNode->monox = thisNode->monoy = thisNode->monoz = 0;
    thisNode->mass = 0;
    TPtlStruct *p2ptl;
    for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling){
		postype mass = p2ptl->mass;
        thisNode->mass += mass;
        thisNode->monox += mass * p2ptl->x;
        thisNode->monoy += mass * p2ptl->y;
        thisNode->monoz += mass * p2ptl->z;
    }
    thisNode->monox = thisNode->monox/thisNode->mass;
    thisNode->monoy = thisNode->monoy/thisNode->mass;
    thisNode->monoz = thisNode->monoz/thisNode->mass;
}


void ex2d_findCellSize(TStruct *thisNode){
    TPtlStruct *nodeparticles = (TPtlStruct *) thisNode->daughter;
    PosType distmax = -1.e20;
    TPtlStruct *p2ptl;
    for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling){
        PosType tmpx = p2ptl->x - thisNode->monox;
        PosType tmpy = p2ptl->y - thisNode->monoy;
        PosType tmpdist2 = tmpx*tmpx + tmpy*tmpy;
        distmax = MAX(distmax, tmpdist2);
    }
    thisNode->nodesize = sqrt(distmax);
}

void ex3d_findCellSize(TStruct *thisNode){
    TPtlStruct *nodeparticles = (TPtlStruct *) thisNode->daughter;
    PosType distmax = -1.e20;
    TPtlStruct *p2ptl;
    for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling){
        PosType tmpx = p2ptl->x - thisNode->monox;
        PosType tmpy = p2ptl->y - thisNode->monoy;
        PosType tmpz = p2ptl->z - thisNode->monoz;
        PosType tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
        distmax = MAX(distmax, tmpdist2);
    }
    thisNode->nodesize = sqrt(distmax);
}



void ex2d_idivision(TStruct *thisNode, TStruct *tmpNode){
    TPtlStruct *nodeparticles = thisNode->daughter;
    TPtlStruct *p2ptl;
    for(p2ptl=nodeparticles;p2ptl;) {
        int mx = (p2ptl->x>thisNode->monox ? 1:0);
        int my = (p2ptl->y>thisNode->monoy ? 1:0);
        int mnode = mx + 2*my;
        tmpNode[mnode].Nparticle ++;
        TPtlStruct *tmpptr = (TPtlStruct *)tmpNode[mnode].daughter;
        TPtlStruct *tmpptr2 = p2ptl->sibling;
        tmpNode[mnode].daughter = (void *)p2ptl;
        p2ptl->sibling = tmpptr;
        p2ptl = tmpptr2;
    }
}

void ex3d_idivision(TStruct *thisNode, TStruct *tmpNode){
    TPtlStruct *nodeparticles = thisNode->daughter;
    TPtlStruct *p2ptl;
    for(p2ptl=nodeparticles;p2ptl;) {
        int mx = (p2ptl->x>thisNode->monox ? 1:0);
        int my = (p2ptl->y>thisNode->monoy ? 1:0);
        int mz = (p2ptl->z>thisNode->monoz ? 1:0);
        int mnode = mx + 2*(my + 2*mz);
        tmpNode[mnode].Nparticle ++;
        TPtlStruct *tmpptr = (TPtlStruct *)tmpNode[mnode].daughter;
        TPtlStruct *tmpptr2 = p2ptl->sibling;
        tmpNode[mnode].daughter = (void *)p2ptl;
        p2ptl->sibling = tmpptr;
        p2ptl = tmpptr2;
    }
}

int nearest2dOpen(void *p, TStruct *tree, postype mindist){
	TPtlStruct *point = (TPtlStruct*)p;
	postype tmpx = point->x - tree->monox;
	postype tmpy = point->y - tree->monoy;
	postype dist2 = tmpx*tmpx + tmpy * tmpy;
	postype dist = sqrt(dist2);
	if( (dist- tree->nodesize) > mindist) return NO;
	else return YES;
}

/* 5-argument nearOpen callback for find_GNear k-NN search */
int near2dOpen_kNN(void *p, TStruct *tree, int npneigh, PosType maxdist, int K){
	if(npneigh >= K){
		TPtlStruct *pt = (TPtlStruct*)p;
		PosType dx = pt->x - tree->monox, dy = pt->y - tree->monoy;
		if(sqrt(dx*dx + dy*dy) - tree->nodesize > maxdist) return NO;
	}
	return YES;
}

int nearest3dOpen(void *p, TStruct *tree, postype mindist){
    TPtlStruct *point = (TPtlStruct*)p;
    postype tmpx = point->x - tree->monox;
    postype tmpy = point->y - tree->monoy;
    postype tmpz = point->z - tree->monoz;
    postype dist2 = tmpx*tmpx + tmpy * tmpy + tmpz *tmpz;
    postype dist = sqrt(dist2);
    if( dist- tree->nodesize > mindist) return NO;
    else return YES;
}


postype ex2d_dist(void *p, void *pp){
	TPtlStruct *point = (TPtlStruct*)p;
	TPtlStruct *neigh = (TPtlStruct *)pp;
	postype tmpx = point->x - neigh->x;
	postype tmpy = point->y - neigh->y;
	postype d2 = tmpx*tmpx + tmpy*tmpy;
	return d2;
}

postype ex3d_dist(void *p, void *pp){
    TPtlStruct *point = (TPtlStruct*)p;
    TPtlStruct *neigh = (TPtlStruct *)pp;
    postype tmpx = point->x - neigh->x;
    postype tmpy = point->y - neigh->y;
    postype tmpz = point->z - neigh->z;
    postype d2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
    return d2;
}



void det2d_dpq(
		SimParameters *simpar, 
		void (*paddingAllTreeParticles)(SimParameters *, postype) 
		){
	treevoroparticletype *bp = VORO_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;

	{
		postype cellsize = HydroGridSize(simpar);
    	paddingAllTreeParticles(simpar, cellsize);
	}
	int mp = VORO_NPAD(simpar);

	TPtlStruct *ptl = (TPtlStruct*)my_malloc(sizeof(TPtlStruct)*(np+mp));
	int nnode = np+mp;
	TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode*2);

	postype xmin,ymin,xmax,ymax;
	xmin = ymin = 1.e20;
	xmax = ymax = -1.e20;
	// this is for the local-domain particles.
	for(i=0;i<np;i++){
		ptl[i].x = bp[i].x;
		ptl[i].y = bp[i].y;
		ptl[i].mass = bp[i].mass; /* in this case mass */
		ptl[i].sibling = ptl + i+1;
		ptl[i].type = TYPE_PTL;
		ptl[i].indx = i;
		xmin = MIN(xmin, bp[i].x);
		ymin = MIN(ymin, bp[i].y);
		xmax = MAX(xmax, bp[i].x);
		ymax = MAX(ymax, bp[i].y);
	}

	// this is for the boundary ghost particles.
	bp = VORO_TBPP(simpar);
	for(i=0;i<mp;i++){
		ptl[np+i].x = bp[i].x;
		ptl[np+i].y = bp[i].y;
		ptl[np+i].mass = bp[i].mass; /* in this case mass */
		ptl[np+i].sibling = ptl + np+i+1;
		ptl[np+i].type = TYPE_PTL;
		ptl[np+i].indx = np+i;
	}
	// terminate the particle sibling pointer
	ptl[np+mp-1].sibling = NULL;
	// tree node has a pointer (daughter) linked to the first particle.
	tree->daughter = ptl;
	// total number of particles linked to the tree node is (np+mp).
	tree->Nparticle = np+mp;

	Make_GNN_Tree(tree, nnode, ex2d_idivision, ex2d_findCentroid,
			ex2d_findCellSize, RECURSIVE);

	// back to the local-domain particles
	bp = VORO_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<np;i++){ 
		bp[i].w2ceil = find_GNearest(ptl+i,  tree, nearest2dOpen, ex2d_dist);
		bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
	}

	my_free(ptl);
	my_free(tree);
	my_free(VORO_TBPP(simpar));
	if(0){
		postype mindpq,maxdpq;
		mindpq = 1.e20;
		maxdpq = -1.e20;
		for(i=0;i<np;i++){
			mindpq = MIN(mindpq, bp[i].w2ceil);
			maxdpq = MAX(maxdpq, bp[i].w2ceil);
		}
		DEBUGPRINT("P%d passed the finding of dpq: min/max= %g %g\n",MYID(simpar), mindpq,maxdpq);
	}
}

void det3d_dpq(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype)
		){
    treevoroparticletype *bp = VORO_TBP(simpar);
    int np = VORO_NP(simpar);

	{
		postype cellsize = HydroGridSize(simpar);
    	paddingAllTreeParticles(simpar, cellsize);

	}
	int mp = VORO_NPAD(simpar);

    TPtlStruct *ptl = (TPtlStruct*)my_malloc(sizeof(TPtlStruct)*(np+mp));
    int nnode = np+mp;
    TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode*2);
    int i;

    for(i=0;i<np;i++){
        ptl[i].x = bp[i].x;
        ptl[i].y = bp[i].y;
        ptl[i].z = bp[i].z;
		ptl[i].mass = bp[i].mass;
        ptl[i].sibling = ptl + i+1;
        ptl[i].type = TYPE_PTL;
        ptl[i].indx = i;
    }

	bp = VORO_TBPP(simpar);
	for(i=0;i<mp;i++){
		ptl[np+i].x = bp[i].x;
		ptl[np+i].y = bp[i].y;
		ptl[np+i].z = bp[i].z;
		ptl[np+i].mass = bp[i].mass; /* in this case mass */
		ptl[np+i].sibling = ptl + np+i+1;
		ptl[np+i].type = TYPE_PTL;
		ptl[np+i].indx = np+i;
	}
	ptl[np+mp-1].sibling = NULL;
    tree->daughter = ptl;
    tree->Nparticle = np+mp;

    Make_GNN_Tree(tree, nnode, ex3d_idivision, ex3d_findCentroid,
            ex3d_findCellSize, RECURSIVE);

	// back to the local-domain particles
	bp = VORO_TBP(simpar);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(i=0;i<np;i++){
        bp[i].w2ceil = find_GNearest(ptl+i,  tree, nearest3dOpen, ex3d_dist);
    }
    my_free(ptl);
    my_free(tree);
	my_free(VORO_TBPP(simpar));
}
void det2d_dpqRK4(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype) ){
	/* Stride-aware: p_size may differ from sizeof(treevorork4particletype)
	   when stress particles are active (av_mode >= 1). */
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;

	int mp;
	{
		postype cellsize = HydroGridSize(simpar);
    	paddingAllTreeParticles(simpar, cellsize);
		mp = VORO_NPAD(simpar);
	}

	TPtlStruct *ptl = (TPtlStruct*)my_malloc(sizeof(TPtlStruct)*(np+mp));
	int nnode = np+mp;
	TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode*2);

	/* Refresh after padding (may realloc) */
	bp_raw = (char*)VORORK4_TBP(simpar);

	postype xmin,ymin,xmax,ymax;
	xmin = ymin = 1.e20;
	xmax = ymax = -1.e20;
	// this is for the local-domain particles.
	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		ptl[i].x = bpi->x;
		ptl[i].y = bpi->y;
		ptl[i].mass = bpi->mass; /* in this case mass */
		ptl[i].sibling = ptl + i+1;
		ptl[i].type = TYPE_PTL;
		ptl[i].indx = i;
		xmin = MIN(xmin, bpi->x);
		ymin = MIN(ymin, bpi->y);
		xmax = MAX(xmax, bpi->x);
		ymax = MAX(ymax, bpi->y);
	}

	// this is for the boundary ghost particles.
	char *pp_raw = (char*)VORORK4_TBPP(simpar);
	for(i=0;i<mp;i++){
		treevorork4particletype *ppi = (treevorork4particletype*)(pp_raw + i*p_size);
		ptl[np+i].x = ppi->x;
		ptl[np+i].y = ppi->y;
		ptl[np+i].mass = ppi->mass; /* in this case mass */
		ptl[np+i].sibling = ptl + np+i+1;
		ptl[np+i].type = TYPE_PTL;
		ptl[np+i].indx = np+i;
	}
	// terminate the particle sibling pointer
	ptl[np+mp-1].sibling = NULL;
	// tree node has a pointer (daughter) linked to the first particle.
	tree->daughter = ptl;
	// total number of particles linked to the tree node is (np+mp).
	tree->Nparticle = np+mp;

	Make_GNN_Tree(tree, nnode, ex2d_idivision, ex2d_findCentroid,
			ex2d_findCellSize, RECURSIVE);

	// back to the local-domain particles
	bp_raw = (char*)VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		bpi->w2ceil = find_GNearest(ptl+i,  tree, nearest2dOpen, ex2d_dist);
		if(gfs_dynw2_active){
			/* GFS_DYN_W2: W is a state variable, no clamp */
		} else if(GAS_Kappa(simpar) == 0)
			bpi->w2 = 0;
		else if(GAS_Kappa(simpar) > 0)
			bpi->w2 = MIN(bpi->w2, bpi->w2ceil);
	}

	my_free(ptl);
	my_free(tree);
	my_free(VORORK4_TBPP(simpar));
	if(0){
		postype mindpq,maxdpq;
		mindpq = 1.e20;
		maxdpq = -1.e20;
		for(i=0;i<np;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			mindpq = MIN(mindpq, bpi->w2ceil);
			maxdpq = MAX(maxdpq, bpi->w2ceil);
		}
		DEBUGPRINT("P%d passed the finding of dpq: min/max= %g %g\n",MYID(simpar), mindpq,maxdpq);
	}
}

/* Build OST tree for GPU NN search.  Same as buildTree2D except it does
   NOT run find_GNearest (the GPU will compute w2ceil instead).  Returns
   tree, ptl, total node count, and the allocated tree-slot capacity.
   Caller frees tree_out, ptl_out, and VORORK4_TBPP. */
void buildTree2D_for_GPU(SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		TStruct **tree_out, TPtlStruct **ptl_out,
		int *nptl_out, int *max_tree_slots_out){
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;

	int mp;
	{
		postype cellsize = HydroGridSize(simpar);
		paddingAllTreeParticles(simpar, cellsize);
		mp = VORO_NPAD(simpar);
	}

	bp_raw = (char*)VORORK4_TBP(simpar);

	TPtlStruct *ptl = (TPtlStruct*)my_malloc(sizeof(TPtlStruct)*(np+mp));
	int nnode = np+mp;
	int max_tree_slots = nnode*2;
	TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*max_tree_slots);

	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		ptl[i].x = bpi->x;
		ptl[i].y = bpi->y;
		ptl[i].mass = bpi->mass;
		ptl[i].sibling = ptl + i+1;
		ptl[i].type = TYPE_PTL;
		ptl[i].indx = i;
		ptl[i].bp = (linkedlisttype*)bpi;
	}

	char *pp_raw = (char*)VORORK4_TBPP(simpar);
	for(i=0;i<mp;i++){
		treevorork4particletype *ppi = (treevorork4particletype*)(pp_raw + i*p_size);
		ptl[np+i].x = ppi->x;
		ptl[np+i].y = ppi->y;
		ptl[np+i].mass = ppi->mass;
		ptl[np+i].sibling = ptl + np+i+1;
		ptl[np+i].type = TYPE_PTL;
		ptl[np+i].indx = np+i;
		ptl[np+i].bp = (linkedlisttype*)ppi;
	}
	ptl[np+mp-1].sibling = NULL;
	tree->daughter = ptl;
	tree->Nparticle = np+mp;

	Make_GNN_Tree(tree, nnode, ex2d_idivision, ex2d_findCentroid,
			ex2d_findCellSize, RECURSIVE);

	*tree_out = tree;
	*ptl_out = ptl;
	*nptl_out = np + mp;
	*max_tree_slots_out = max_tree_slots;
}

/* Build k-d tree for k-NN search. Tree, ptl, and padding are kept alive
   for the caller (updateDen_kNN + getAcc_kNN) to reuse.
   Caller is responsible for freeing tree_out, ptl_out, and VORORK4_TBPP. */
void buildTree2D(SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		int (*skipBP)(treevorork4particletype *),
		TStruct **tree_out, TPtlStruct **ptl_out, int *nptl_out){
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;

	int mp;
	{
		postype cellsize = HydroGridSize(simpar);
		paddingAllTreeParticles(simpar, cellsize);
		mp = VORO_NPAD(simpar);
	}

	DEBUGPRINT("P%d buildTree2D_kNN: after padding np=%d mp=%d\n",
		MYID(simpar), np, mp);

	/* Refresh after padding (may realloc) */
	bp_raw = (char*)VORORK4_TBP(simpar);

	TPtlStruct *ptl = (TPtlStruct*)my_malloc(sizeof(TPtlStruct)*(np+mp));
	int nnode = np+mp;
	TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode*2);

	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		ptl[i].x = bpi->x;
		ptl[i].y = bpi->y;
		ptl[i].mass = bpi->mass;
		ptl[i].sibling = ptl + i+1;
		ptl[i].type = TYPE_PTL;
		ptl[i].indx = i;
		ptl[i].bp = (linkedlisttype*)bpi;
	}

	/* Boundary ghost particles */
	char *pp_raw = (char*)VORORK4_TBPP(simpar);
	for(i=0;i<mp;i++){
		treevorork4particletype *ppi = (treevorork4particletype*)(pp_raw + i*p_size);
		ptl[np+i].x = ppi->x;
		ptl[np+i].y = ppi->y;
		ptl[np+i].mass = ppi->mass;
		ptl[np+i].sibling = ptl + np+i+1;
		ptl[np+i].type = TYPE_PTL;
		ptl[np+i].indx = np+i;
		ptl[np+i].bp = (linkedlisttype*)ppi; /* back-reference to padding particle */
	}
	ptl[np+mp-1].sibling = NULL;
	tree->daughter = ptl;
	tree->Nparticle = np+mp;

	Make_GNN_Tree(tree, nnode, ex2d_idivision, ex2d_findCentroid,
			ex2d_findCellSize, RECURSIVE);

	/* Compute w2ceil for local particles */
	bp_raw = (char*)VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		if(skipBP && skipBP(bpi)) continue;
		bpi->w2ceil = find_GNearest(ptl+i, tree, nearest2dOpen, ex2d_dist);
		if(GAS_Kappa(simpar) >= 0)
			bpi->w2 = MIN(bpi->w2, bpi->w2ceil);
	}

	*tree_out = tree;
	*ptl_out = ptl;
	*nptl_out = np + mp;
}

/* k-NN search around a single query point.
   Returns Voro2D_point neighbor array (malloc'd) and count.
   maxr is set to the distance of the K-th neighbor. */
int searchKNN2D(TPtlStruct *query, int K, TStruct *tree,
		Voro2D_point *neighbors, PosType *maxr){
	void *bpneighbor[MAX_NUM_NEAR];
	int nfound;

	nfound = find_GNear(query, K, tree, maxr, bpneighbor,
			near2dOpen_kNN, ex2d_dist, insGnear);

	/* Debug: if maxr=0, do a manual tree traversal with counters */
	{
		static int dbg_maxr0 = 0;
		if(*maxr < 1e-10 && nfound >= K && dbg_maxr0 < 1){
			dbg_maxr0++;
			fprintf(stderr, "searchKNN2D maxr=0: query=(%g,%g) K=%d nfound=%d maxr=%g\n",
				query->x, query->y, K, nfound, *maxr);
			/* Manual traversal mimicking find_GNear */
			int n_tree=0, n_ptl=0, n_open=0, n_skip=0, n_insert=0;
			int n_dist0=0, n_same_ptr=0;
			void *dbg_ptr = (void*)tree;
			PosType dbg_maxdist2 = 1.E23, dbg_maxdist = 1.E23;
			int dbg_npneigh = 0;
			void *first_ptl_ptr = NULL;
			while(dbg_ptr != NULL){
				if(((TYPE*)dbg_ptr)->type == TYPE_TREE){
					n_tree++;
					/* Always open (same as near2dOpen_kNN with nodesize=3e10) */
					int open_res = near2dOpen_kNN(query, (TStruct*)dbg_ptr,
						dbg_npneigh, dbg_maxdist, K);
					if(open_res == YES){
						n_open++;
						dbg_ptr = (void*)(((TStruct*)dbg_ptr)->daughter);
					} else {
						n_skip++;
						dbg_ptr = (void*)(((TStruct*)dbg_ptr)->sibling);
					}
				} else {
					n_ptl++;
					PosType d2 = ex2d_dist(query, dbg_ptr);
					if(d2 < 1e-20) n_dist0++;
					if(dbg_ptr == (void*)query) n_same_ptr++;
					if(first_ptl_ptr == NULL) first_ptl_ptr = dbg_ptr;
					/* Simulate insertion logic */
					if(dbg_npneigh < K || d2 < dbg_maxdist2){
						n_insert++;
						dbg_npneigh++;
						if(dbg_npneigh > K) dbg_npneigh = K;
						/* Track worst case maxdist2 */
						if(d2 > 0 && d2 < dbg_maxdist2) dbg_maxdist2 = d2;
						if(d2 == 0 && dbg_npneigh == 1) dbg_maxdist2 = 0;
					}
					dbg_maxdist = sqrt(dbg_maxdist2);
					dbg_ptr = (void*)(((TPtlStruct*)dbg_ptr)->sibling);
				}
				if(n_ptl + n_tree > 50000) {
					fprintf(stderr, "  ABORT: too many visits\n");
					break;
				}
			}
			fprintf(stderr, "  Manual traversal: tree_nodes=%d ptl_nodes=%d opens=%d skips=%d\n",
				n_tree, n_ptl, n_open, n_skip);
			fprintf(stderr, "  inserts=%d dist0=%d same_ptr_as_query=%d first_ptl=%p query=%p\n",
				n_insert, n_dist0, n_same_ptr, first_ptl_ptr, (void*)query);
			/* Print first few bpneighbor from actual find_GNear */
			int di;
			for(di=0;di<3 && di<nfound;di++){
				TPtlStruct *nb2 = (TPtlStruct*)bpneighbor[di];
				fprintf(stderr, "  bpneigh[%d]: pos=(%g,%g) bp=%p\n",
					di, nb2->x, nb2->y, (void*)nb2);
			}
		}
	}

	int i;
	for(i=0;i<nfound;i++){
		TPtlStruct *nb = (TPtlStruct*)bpneighbor[i];
		treevorork4particletype *tt = (treevorork4particletype*)(nb->bp);
		neighbors[i].x = tt->x;
		neighbors[i].y = tt->y;
		neighbors[i].indx = PINDX(tt);
		neighbors[i].csound = tt->csound;
		neighbors[i].pressure = tt->pressure;
		neighbors[i].w2 = tt->w2;
		neighbors[i].bp = (void*)tt;
		neighbors[i].u4if.indx = MAX_INDEX;
	}
	return nfound;
}


void det3d_dpqRK4(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype)
		){
    treevorork4particletype *bp = VORORK4_TBP(simpar);
    int np = VORO_NP(simpar);

	{
		postype cellsize = HydroGridSize(simpar);
    	paddingAllTreeParticles(simpar, cellsize);

	}
	int mp = VORO_NPAD(simpar);

    TPtlStruct *ptl = (TPtlStruct*)my_malloc(sizeof(TPtlStruct)*(np+mp));
    int nnode = np+mp;
    TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode*2);
    int i;

    for(i=0;i<np;i++){
        ptl[i].x = bp[i].x;
        ptl[i].y = bp[i].y;
        ptl[i].z = bp[i].z;
		ptl[i].mass = bp[i].mass;
        ptl[i].sibling = ptl + i+1;
        ptl[i].type = TYPE_PTL;
        ptl[i].indx = i;
    }

	bp = VORORK4_TBPP(simpar);
	for(i=0;i<mp;i++){
		ptl[np+i].x = bp[i].x;
		ptl[np+i].y = bp[i].y;
		ptl[np+i].z = bp[i].z;
		ptl[np+i].mass = bp[i].mass; /* in this case mass */
		ptl[np+i].sibling = ptl + np+i+1;
		ptl[np+i].type = TYPE_PTL;
		ptl[np+i].indx = np+i;
	}
	ptl[np+mp-1].sibling = NULL;
    tree->daughter = ptl;
    tree->Nparticle = np+mp;

    Make_GNN_Tree(tree, nnode, ex3d_idivision, ex3d_findCentroid,
            ex3d_findCellSize, RECURSIVE);

	// back to the local-domain particles
	bp = VORORK4_TBP(simpar);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(i=0;i<np;i++){
        bp[i].w2ceil = find_GNearest(ptl+i,  tree, nearest3dOpen, ex3d_dist);
    }
    my_free(ptl);
    my_free(tree);
	my_free(VORORK4_TBPP(simpar));
}

Voro2D_point *searchCellRk4Neighbors2D(
		SimParameters *simpar,
        int ix, 
		int iy, 
		int *nneigh
		){
    int i,j,k;
    int np,mp;
    Voro2D_point *neigh;
    int mx,my;
    mx = BASICCELL_MX(simpar);
    my = BASICCELL_MY(simpar);
    CellType *cells = VORO_BASICCELL(simpar);
    float Lx = KH_SIMBOX(simpar).x.max - KH_SIMBOX(simpar).x.min;
    float Ly = KH_SIMBOX(simpar).y.max - KH_SIMBOX(simpar).y.min;

    np = 0;
    for(j=iy-1;j<=iy+1;j++){
		if(j<0 || j >= my) continue;
        for(i=ix-1;i<=ix+1;i++){
			if(i<0 || i >= mx) continue;
            np += cells[i+mx*j].nmem;
        }
    }
	// twice the size for the boundary conditions
    neigh = (Voro2D_point*)malloc(sizeof(Voro2D_point)*np);
    np = 0;
    for(j=iy-1;j<=iy+1;j++){
		if( j<0 || j >= my) continue;
        for(i=ix-1;i<=ix+1;i++){
			if(i<0 || i >= mx) continue;
            struct linkedlisttype *tmp = cells[i+mx*j].link;
            while(tmp){
                treevorork4particletype *tt = (treevorork4particletype*)tmp;
                neigh[np].x = tt->x;
                neigh[np].y = tt->y;
                neigh[np].indx = PINDX(tt);
                neigh[np].csound = tt->csound;
                neigh[np].pressure = tt->pressure;
                neigh[np].w2 = tt->w2;
                neigh[np].bp = (void*)tt;
                neigh[np].u4if = tt->u4if;
                np++;
                tmp = tmp->next;
            }

        }
    }
    *nneigh = np;
    neigh = realloc(neigh,sizeof(treevorork4particletype)*np);
    return neigh;
}

/* ================================================================
 *  Adaptive neighbor search helpers (Fix B, May 2026)
 *
 *  Root-cause fix for SIGSEGV in M(n,m) face-stencil macros
 *  (VoroStressRK4_Pressure2D etc.).  The macros dereference
 *      neighwork[tmp->upperlink->upperrelated].bp
 *      neighwork[tmp->lowerrelated].bp
 *  When Voro2D_FindVC leaves a corner with negative upperrelated
 *  (= sentinel from Voro2D_InitializeCorner that was never replaced
 *  by a real candidate during cell construction), the dereference
 *  underflows the array and crashes.
 *
 *  The sentinel survives when no candidate particle exists in that
 *  direction within the search radius.  In Sedov-like radial outflow,
 *  particles ahead of the shock front move beyond the fixed 3x3
 *  basic-cell stencil (search reach ~1.5 * basic_cellsize), leaving
 *  the angular sector empty.
 *
 *  Fix: detect surviving negative-relation corners post-FindVC and
 *  retry with progressively larger candidate stencils (R=2 -> 5x5,
 *  R=3 -> 7x7, R=4 -> 9x9) until cell closes or R_MAX is reached.
 * ================================================================ */

/* Radius-aware variant of searchCellRk4Neighbors2D.  R=1 reproduces
 * the original 3x3 stencil exactly. */
Voro2D_point *searchCellRk4Neighbors2D_R(
        SimParameters *simpar, int ix, int iy, int R, int *nneigh){
    int i,j,np;
    Voro2D_point *neigh;
    int mx = BASICCELL_MX(simpar);
    int my = BASICCELL_MY(simpar);
    CellType *cells = VORO_BASICCELL(simpar);

    np = 0;
    for(j=iy-R;j<=iy+R;j++){
        if(j<0 || j>=my) continue;
        for(i=ix-R;i<=ix+R;i++){
            if(i<0 || i>=mx) continue;
            np += cells[i+mx*j].nmem;
        }
    }
    if(np == 0){ *nneigh = 0; return NULL; }
    neigh = (Voro2D_point*)malloc(sizeof(Voro2D_point)*np);
    np = 0;
    for(j=iy-R;j<=iy+R;j++){
        if(j<0 || j>=my) continue;
        for(i=ix-R;i<=ix+R;i++){
            if(i<0 || i>=mx) continue;
            struct linkedlisttype *tmp = cells[i+mx*j].link;
            while(tmp){
                treevorork4particletype *tt = (treevorork4particletype*)tmp;
                neigh[np].x = tt->x;
                neigh[np].y = tt->y;
                neigh[np].indx = PINDX(tt);
                neigh[np].csound = tt->csound;
                neigh[np].pressure = tt->pressure;
                neigh[np].w2 = tt->w2;
                neigh[np].bp = (void*)tt;
                neigh[np].u4if = tt->u4if;
                np++;
                tmp = tmp->next;
            }
        }
    }
    *nneigh = np;
    return neigh;
}

/* Returns 1 if any active corner has a negative upperrelated or
 * lowerrelated index — i.e., FindVC left a sentinel uncut in some
 * angular direction.  Returns 0 if all relations are valid. */
int Voro2D_HasUnassignedCorner(Voro2D_Corner *vc, int ip){
    int i;
    for(i=0;i<ip;i++){
        if(vc[i].status == Active){
            if(vc[i].upperrelated < 0 || vc[i].lowerrelated < 0) return 1;
        }
    }
    return 0;
}

treevorork4particletype *findCellRk4BP2D(SimParameters *simpar, int ix, int iy, int *mp){
    int i,j,k;
    int np;
    treevorork4particletype *res;
    int mx,my;
    mx = BASICCELL_MX(simpar);
    my = BASICCELL_MY(simpar);
    CellType *cells = VORO_BASICCELL(simpar);

    np = cells[ix+mx*iy].nmem;
	// twice the size for the boundary conditions
    res = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*np);
    np = 0;
    int ipixel = ix+mx*iy;
    struct linkedlisttype *tmp = cells[ipixel].link;
    while(tmp){
        if(!IS_FLAG_ON(tmp,BoundaryGhostflag)){
            treevorork4particletype *tt = (treevorork4particletype*)tmp;
            {
                res[np] = *tt;
                res[np].bp = tt;
                np++;
            }
        }
        tmp = tmp->next;
    }
    *mp = np;
    res = realloc(res,sizeof(treevorork4particletype)*np);
    return res;
}

/*
void new_mkLinkedList2D_rt(
		SimParameters *simpar, 
		postype cellsize, 
		postype xmin, 
		postype ymin,
		postype xmax, 
		postype ymax, 
		void (*paddingTreeAllParticle)(SimParameters *, postype)
		){ 
	treevorork4particletype *bp = VORORK4_TBP(simpar); 
	int np = VORO_NP(simpar);
    int i,j,k;
    int mx,my;
	postype Gamma = GAS_GAMMA(simpar);
	postype Ly = SIMBOX(simpar).y.max;

    mx = BASICCELL_MX(simpar);
    my = BASICCELL_MY(simpar);

	CellType *cells = VORO_BASICCELL(simpar); 

	for(i=0;i<mx*my;i++) {
        cells[i].link = NULL;
        cells[i].nmem = 0;
    }
	for(i=0;i<VORO_NP(simpar);i++){
        int ix,iy;
        ix = ((bp[i].x-xmin)/cellsize);
        iy = ((bp[i].y-ymin)/cellsize);
        int index = ix+mx*iy;
        if(index >=mx*my || index <0){
            DEBUGPRINT("error detected: %d : %d %d : %g %g %g %g\n",
                    index, mx,my, bp[i].x, bp[i].y, xmin,ymin);
            continue;
        }
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)(bp+i);
        cells[index].nmem ++;
        bp[i].next = tmp;
    }
	// now for the particle padding from other ranks
    paddingTreeAllParticle(simpar, cellsize);

	{
	    treevorork4particletype *bpp = VORORK4_TBPP(simpar); 
		for(i=0;i<VORO_NPAD(simpar);i++){ 
			int ix,iy; 
			ix = ((bpp[i].x-xmin)/cellsize); 
			iy = ((bpp[i].y-ymin)/cellsize); 
			int index = ix+mx*iy; 
			struct linkedlisttype *tmp = cells[index].link; 
			cells[index].link = (struct linkedlisttype*)(bpp+i); 
			cells[index].nmem ++; 
			bpp[i].next = tmp; 
		}
	}



	{
		postype rho1 = RT_DEN1(simpar);
		postype rho2 = RT_DEN2(simpar);
		postype deltay = RT_Deltay(simpar);
		postype ycen = 0.5*Ly;
		postype getDen(postype , postype , postype , postype ,postype );
		int ny = NY(simpar);
		treevorork4particletype *bpp = VORORK4_TBPP(simpar);
		postype dy = Ly/(postype)ny;
//		DEBUGPRINT("P%d has %g %g deltay= %g pcenter = %g\n", MYID(simpar), rho1, rho2, deltay, RT_Phalf(simpar));
		// top pressure
		postype pup = RT_Phalf(simpar);
		for(j=ny/2+1;j<ny+1;j++){
			postype y = (postype)(j+0.5)*Ly/(postype)ny; 
			postype rho = getDen(rho1, rho2, deltay, ycen, y); 
			postype dr = dy;
			pup += GAS_ACCY(simpar) * rho * dr;
		}
		for(i=0;i<VORO_NPAD(simpar);i++){
			postype getDen(postype , postype , postype , postype ,postype );
			if(bpp[i].y >= Ly) {
				bpp[i].den = getDen(rho1, rho2, deltay, ycen, bpp[i].y);
				postype rho = getDen(rho1, rho2, deltay, ycen, 0.5*(bpp[i].y+Ly));
				postype dr = (bpp[i].y- Ly);
				bpp[i].pressure = pup + GAS_ACCY(simpar)* rho * dr;
				bpp[i].ie = bpp[i].pressure*bpp[i].volume/(Gamma-1);
        		bpp[i].csound = sqrt(Gamma*bpp[i].pressure/bpp[i].den);
		        if(GAS_Kappa(simpar) <0) { 
					bpp[i].w2 = -GAS_Kappa(simpar);
		        }
	       		else if (GAS_Kappa(simpar) >0){ 
					bpp[i].w2 = getw2forHydroParticle(simpar,(bpp+i));
        		}
				bpp[i].w2old = bpp[i].w2; // static boundary conditions.
				bpp[i].vx = bpp[i].vy = 0;
			}
		}
		// bottom pressure
		postype pbottom = RT_Phalf(simpar);
		for(j=ny/2-1;j>=0;j--){
			postype y = (postype)(j+0.5)*Ly/(postype)ny;
			postype rho = getDen(rho1,rho2, deltay,ycen,y);
			postype dr = -dy;
			pbottom += GAS_ACCY(simpar)*rho*dr;
		}
		for(i=0;i<VORO_NPAD(simpar);i++){
			postype getDen(postype , postype , postype , postype ,postype );
			if(bpp[i].y < 0) {
				bpp[i].den = getDen(rho1, rho2, deltay, ycen, bpp[i].y);
				postype rho = getDen(rho1, rho2, deltay, ycen, bpp[i].y*0.5);
				postype dr = bpp[i].y;
				bpp[i].pressure =  pbottom +  GAS_ACCY(simpar)* rho * dr;
				bpp[i].ie = bpp[i].pressure*bpp[i].volume/(Gamma-1);
        		bpp[i].csound = sqrt(Gamma*bpp[i].pressure/bpp[i].den);
		        if(GAS_Kappa(simpar) <0) { 
					bpp[i].w2 = -GAS_Kappa(simpar);
		        }
	       		else if (GAS_Kappa(simpar) >0){ 
					bpp[i].w2 = getw2forHydroParticle(simpar,(bpp+i)); 
        		}
				bpp[i].w2old = bpp[i].w2;// static boundary conditions
				bpp[i].vx = bpp[i].vy = 0;
			}
		}
	}
}
*/

void mkLinkedList2D_rt(
		SimParameters *simpar,
		postype cellsize,
		postype xmin,
		postype ymin,
		postype xmax,
		postype ymax,
		void (*paddingTreeAllParticle)(SimParameters *, postype)
		){
	/* p_size: actual byte size per particle (handles treevorostressrk4particletype
	   when av_mode >= 1, where n_size is overridden in RunRT) */
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
    int i,j,k;
    int mx,my;
	postype Gamma = GAS_GAMMA(simpar);
	postype Ly = SIMBOX(simpar).y.max;
	postype rho1 = RT_DEN1(simpar), rho2 = RT_DEN2(simpar);
	postype deltay = RT_Deltay(simpar);
	postype ycen = 0.5*Ly;

    mx = BASICCELL_MX(simpar);
    my = BASICCELL_MY(simpar);

	CellType *cells = VORO_BASICCELL(simpar);

	for(i=0;i<mx*my;i++) {
        cells[i].link = NULL;
        cells[i].nmem = 0;
    }
	/* Build linked list for local particles (stride-aware) */
	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
        int ix,iy;
        ix = ((bpi->x-xmin)/cellsize);
        iy = ((bpi->y-ymin)/cellsize);
        int index = ix+mx*iy;
        if(index >=mx*my || index <0){
            DEBUGPRINT("error detected: %d : %d %d : %g %g %g %g\n",
                    index, mx,my, bpi->x, bpi->y, xmin,ymin);
            continue;
        }
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)bpi;
        cells[index].nmem ++;
        bpi->next = tmp;
    }
	// now for the particle padding from other ranks
    paddingTreeAllParticle(simpar, cellsize);

	postype getRT_Den(postype , postype , postype , postype ,postype );
	int ny = NY(simpar);

	/* Top wall mirror (y = Ly): compact padding, create mirrors
	   Stride-aware version of the original algorithm. */
	if(Ymax_HydroExam(simpar) >= (Ly - EPS)){
		char *bpp_raw = (char*)VORORK4_TBPP(simpar);
		int mpad = 0;
		for(i=0;i<VORO_NPAD(simpar);i++){
			treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
			if(bppi->y < Ly) {
				if(mpad != i) memcpy(bpp_raw + mpad*p_size, bpp_raw + i*p_size, p_size);
				mpad ++;
			}
		}
		int nmirror = 0;
		for(i=0;i<np;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			if((Ly-bpi->y) < cellsize) nmirror++;
		}
		for(i=0;i<mpad;i++){
			treevorork4particletype *bppi = (treevorork4particletype*)((char*)VORORK4_TBPP(simpar) + i*p_size);
			if((Ly-bppi->y) < cellsize) nmirror++;
		}
		if(VORO_NPAD(simpar) != (mpad+nmirror)) {
			VORORK4_TBPP(simpar) = (treevorork4particletype *)
					realloc(VORORK4_TBPP(simpar), p_size*(mpad+nmirror));
		}
		bpp_raw = (char*)VORORK4_TBPP(simpar);
		char *cursor = bpp_raw + mpad*p_size;
		for(i=0;i<np;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			if((Ly-bpi->y) < cellsize) {
				memcpy(cursor, bpi, p_size);
				treevorork4particletype *gp = (treevorork4particletype*)cursor;
				postype dy = 2*(Ly-bpi->y);
				gp->u4if.indx = MAX_INDEX;
				gp->y += dy;
				gp->vy = -bpi->vy;
				gp->den = getRT_Den(rho1, rho2, deltay, ycen, gp->y);
				gp->pressure += GAS_ACCY(simpar)* getRT_Den(rho1,rho2,deltay,ycen,Ly) * dy;
				gp->ie = gp->pressure*gp->volume/(Gamma-1);
				CLEAR_FLAG(gp);
				SET_FLAG(gp,BoundaryGhostflag);
				SET_FLAG(gp,VOROflag);
				cursor += p_size;
			}
		}
		for(i=0;i<mpad;i++){
			treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
			if((Ly-bppi->y) < cellsize) {
				memcpy(cursor, bppi, p_size);
				treevorork4particletype *gp = (treevorork4particletype*)cursor;
				postype dy = 2*(Ly-bppi->y);
				gp->u4if.indx = MAX_INDEX;
				gp->y += dy;
				gp->vy = -bppi->vy;
				gp->den = getRT_Den(rho1, rho2, deltay, ycen, gp->y);
				gp->pressure += GAS_ACCY(simpar)* getRT_Den(rho1,rho2,deltay,ycen,Ly) * dy;
				gp->ie = gp->pressure*gp->volume/(Gamma-1);
				CLEAR_FLAG(gp);
				SET_FLAG(gp,BoundaryGhostflag);
				SET_FLAG(gp,VOROflag);
				cursor += p_size;
			}
		}
		VORO_NPAD(simpar) = mpad + nmirror;
	}
	/* Bottom wall mirror (y = 0) */
	if(Ymin_HydroExam(simpar)< EPS)
	{
		char *bpp_raw = (char*)VORORK4_TBPP(simpar);
		int mpad = 0;
		for(i=0;i<VORO_NPAD(simpar);i++){
			treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
			if(bppi->y >= 0.) {
				if(mpad != i) memcpy(bpp_raw + mpad*p_size, bpp_raw + i*p_size, p_size);
				mpad ++;
			}
		}
		int nmirror = 0;
		for(i=0;i<np;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			if(bpi->y <= cellsize) nmirror++;
		}
		for(i=0;i<mpad;i++){
			treevorork4particletype *bppi = (treevorork4particletype*)((char*)VORORK4_TBPP(simpar) + i*p_size);
			if(bppi->y <= cellsize) nmirror++;
		}
		if(VORO_NPAD(simpar) != (mpad+nmirror)) {
			VORORK4_TBPP(simpar) = (treevorork4particletype *) realloc(VORORK4_TBPP(simpar),
						p_size*(mpad+nmirror));
		}
		bpp_raw = (char*)VORORK4_TBPP(simpar);
		char *cursor = bpp_raw + mpad*p_size;
		for(i=0;i<np;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			if(bpi->y <= cellsize) {
				memcpy(cursor, bpi, p_size);
				treevorork4particletype *gp = (treevorork4particletype*)cursor;
				gp->u4if.indx = MAX_INDEX;
				postype dy = -bpi->y;
				gp->y = dy;
				gp->vy = -bpi->vy;
				gp->den = getRT_Den(rho1, rho2, deltay, ycen, gp->y);
				gp->pressure -= GAS_ACCY(simpar)* getRT_Den(rho1,rho2,deltay,ycen,0.) * 2*(-dy);
				gp->ie = gp->pressure*gp->volume/(Gamma-1);
				CLEAR_FLAG(gp);
				SET_FLAG(gp,BoundaryGhostflag);
				SET_FLAG(gp,VOROflag);
				cursor += p_size;
			}
		}
		for(i=0;i<mpad;i++){
			treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
			if(bppi->y <= cellsize) {
				memcpy(cursor, bppi, p_size);
				treevorork4particletype *gp = (treevorork4particletype*)cursor;
				gp->u4if.indx = MAX_INDEX;
				postype dy = -bppi->y;
				gp->y = dy;
				gp->vy = -bppi->vy;
				gp->den = getRT_Den(rho1, rho2, deltay, ycen, gp->y);
				gp->pressure -= GAS_ACCY(simpar)* getRT_Den(rho1,rho2,deltay,ycen,0.) * 2*(-dy);
				gp->ie = gp->pressure*gp->volume/(Gamma-1);
				CLEAR_FLAG(gp);
				SET_FLAG(gp,BoundaryGhostflag);
				SET_FLAG(gp,VOROflag);
				cursor += p_size;
			}
		}
		VORO_NPAD(simpar) = mpad + nmirror;
	}
	/* Insert padding + mirror particles into linked list (stride-aware) */
	{
		char *bpp_raw = (char*)VORORK4_TBPP(simpar);
		for(i=0;i<VORO_NPAD(simpar);i++){
			treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
			int ix,iy;
			ix = ((bppi->x-xmin)/cellsize);
			iy = ((bppi->y-ymin)/cellsize);
			int index = ix+mx*iy;
			struct linkedlisttype *tmp = cells[index].link;
			cells[index].link = (struct linkedlisttype*)bppi;
			cells[index].nmem ++;
			bppi->next = tmp;
		}
	}
}

/* ================================================================
 * mkLinkedList2D_sod:
 *   Reflecting walls at x=0 and x=Lx, periodic in y.  Used by the
 *   IC_SOD2D test so that periodic x-wraparound does not pollute
 *   the analysis region.  After calling paddingTreeAllParticle
 *   (fully periodic), x-wraparound padding particles (those with
 *   x outside the physical box) are compacted out and replaced by
 *   wall-mirror particles (vx reversed, flagged MAX_INDEX).
 * ================================================================ */
void mkLinkedList2D_sod(
		SimParameters *simpar,
		postype cellsize,
		postype xmin,
		postype ymin,
		postype xmax,
		postype ymax,
		void (*paddingTreeAllParticle)(SimParameters *, postype)
		){
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;
	int mx = BASICCELL_MX(simpar);
	int my = BASICCELL_MY(simpar);
	CellType *cells = VORO_BASICCELL(simpar);
	postype Lx_lo = SIMBOX(simpar).x.min;
	postype Lx_hi = SIMBOX(simpar).x.max;

	for(i=0;i<mx*my;i++){ cells[i].link = NULL; cells[i].nmem = 0; }

	/* Bin local particles */
	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		int ix = ((bpi->x-xmin)/cellsize);
		int iy = ((bpi->y-ymin)/cellsize);
		int index = ix+mx*iy;
		if(index >= mx*my || index < 0) continue;
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)bpi;
		cells[index].nmem++;
		bpi->next = tmp;
	}

	/* Periodic padding (X and Y) */
	paddingTreeAllParticle(simpar, cellsize);

	/* Compact out X-wraparound padding (x < Lx_lo or x >= Lx_hi) */
	char *bpp_raw = (char*)VORORK4_TBPP(simpar);
	int mpad = 0;
	for(i=0;i<VORO_NPAD(simpar);i++){
		treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
		if(bppi->x >= Lx_lo && bppi->x < Lx_hi){
			if(mpad != i) memcpy(bpp_raw + mpad*p_size, bpp_raw + i*p_size, p_size);
			mpad++;
		}
	}

	/* Count wall mirrors for both walls (local + kept-padding) */
	int nmirror = 0;
	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		if(bpi->x - Lx_lo < cellsize) nmirror++;
		if(Lx_hi - bpi->x < cellsize) nmirror++;
	}
	for(i=0;i<mpad;i++){
		treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
		if(bppi->x - Lx_lo < cellsize) nmirror++;
		if(Lx_hi - bppi->x < cellsize) nmirror++;
	}

	if(VORO_NPAD(simpar) != (mpad + nmirror)){
		VORORK4_TBPP(simpar) = (treevorork4particletype*)
			realloc(VORORK4_TBPP(simpar), p_size * (mpad + nmirror));
	}
	bpp_raw = (char*)VORORK4_TBPP(simpar);
	char *cursor = bpp_raw + mpad*p_size;

	/* Left-wall mirrors from local */
	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		if(bpi->x - Lx_lo < cellsize){
			memcpy(cursor, bpi, p_size);
			treevorork4particletype *gp = (treevorork4particletype*)cursor;
			gp->u4if.indx = MAX_INDEX;
			gp->x = 2.0*Lx_lo - bpi->x;
			gp->vx = -bpi->vx;
			CLEAR_FLAG(gp);
			SET_FLAG(gp, BoundaryGhostflag);
			SET_FLAG(gp, VOROflag);
			cursor += p_size;
		}
	}
	/* Right-wall mirrors from local */
	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		if(Lx_hi - bpi->x < cellsize){
			memcpy(cursor, bpi, p_size);
			treevorork4particletype *gp = (treevorork4particletype*)cursor;
			gp->u4if.indx = MAX_INDEX;
			gp->x = 2.0*Lx_hi - bpi->x;
			gp->vx = -bpi->vx;
			CLEAR_FLAG(gp);
			SET_FLAG(gp, BoundaryGhostflag);
			SET_FLAG(gp, VOROflag);
			cursor += p_size;
		}
	}
	/* Left-wall mirrors from kept padding (y-periodic images near wall) */
	for(i=0;i<mpad;i++){
		treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
		if(bppi->x - Lx_lo < cellsize){
			memcpy(cursor, bppi, p_size);
			treevorork4particletype *gp = (treevorork4particletype*)cursor;
			gp->u4if.indx = MAX_INDEX;
			gp->x = 2.0*Lx_lo - bppi->x;
			gp->vx = -bppi->vx;
			CLEAR_FLAG(gp);
			SET_FLAG(gp, BoundaryGhostflag);
			SET_FLAG(gp, VOROflag);
			cursor += p_size;
		}
	}
	/* Right-wall mirrors from kept padding */
	for(i=0;i<mpad;i++){
		treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
		if(Lx_hi - bppi->x < cellsize){
			memcpy(cursor, bppi, p_size);
			treevorork4particletype *gp = (treevorork4particletype*)cursor;
			gp->u4if.indx = MAX_INDEX;
			gp->x = 2.0*Lx_hi - bppi->x;
			gp->vx = -bppi->vx;
			CLEAR_FLAG(gp);
			SET_FLAG(gp, BoundaryGhostflag);
			SET_FLAG(gp, VOROflag);
			cursor += p_size;
		}
	}
	VORO_NPAD(simpar) = mpad + nmirror;

	/* Bin all padding (y-periodic + wall mirrors) into cells */
	bpp_raw = (char*)VORORK4_TBPP(simpar);
	for(i=0;i<VORO_NPAD(simpar);i++){
		treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
		int ix = ((bppi->x-xmin)/cellsize);
		int iy = ((bppi->y-ymin)/cellsize);
		int index = ix+mx*iy;
		if(index >= mx*my || index < 0) continue;
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)bppi;
		cells[index].nmem++;
		bppi->next = tmp;
	}
}

/* ================================================================
 * mkLinkedList2D_sedov2d:
 *   Reflecting walls on all four sides (x=Lx_lo, Lx_hi, y=Ly_lo, Ly_hi).
 *   Used for the 2D Sedov blast wave. Periodic padding from
 *   paddingTreeAllParticle() is compacted to drop both x- and y-
 *   wraparound images, then 4-wall mirror particles are appended
 *   (vx,vy reversed at the corresponding wall, indx=MAX_INDEX).
 * ================================================================ */
void mkLinkedList2D_sedov2d(
		SimParameters *simpar,
		postype cellsize,
		postype xmin,
		postype ymin,
		postype xmax,
		postype ymax,
		void (*paddingTreeAllParticle)(SimParameters *, postype)
		){
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;
	int mx = BASICCELL_MX(simpar);
	int my = BASICCELL_MY(simpar);
	CellType *cells = VORO_BASICCELL(simpar);
	postype Lx_lo = SIMBOX(simpar).x.min;
	postype Lx_hi = SIMBOX(simpar).x.max;
	postype Ly_lo = SIMBOX(simpar).y.min;
	postype Ly_hi = SIMBOX(simpar).y.max;

	for(i=0;i<mx*my;i++){ cells[i].link = NULL; cells[i].nmem = 0; }

	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		int ix = ((bpi->x-xmin)/cellsize);
		int iy = ((bpi->y-ymin)/cellsize);
		int index = ix+mx*iy;
		if(index >= mx*my || index < 0) continue;
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)bpi;
		cells[index].nmem++;
		bpi->next = tmp;
	}

	paddingTreeAllParticle(simpar, cellsize);

	/* Drop wraparound padding in BOTH x and y: keep only padding inside box */
	char *bpp_raw = (char*)VORORK4_TBPP(simpar);
	int mpad = 0;
	for(i=0;i<VORO_NPAD(simpar);i++){
		treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
		if(bppi->x >= Lx_lo && bppi->x < Lx_hi &&
		   bppi->y >= Ly_lo && bppi->y < Ly_hi){
			if(mpad != i) memcpy(bpp_raw + mpad*p_size, bpp_raw + i*p_size, p_size);
			mpad++;
		}
	}

	/* Count wall mirrors over local + kept padding for all 4 walls */
	int nmirror = 0;
	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		if(bpi->x - Lx_lo < cellsize) nmirror++;
		if(Lx_hi - bpi->x < cellsize) nmirror++;
		if(bpi->y - Ly_lo < cellsize) nmirror++;
		if(Ly_hi - bpi->y < cellsize) nmirror++;
	}
	for(i=0;i<mpad;i++){
		treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
		if(bppi->x - Lx_lo < cellsize) nmirror++;
		if(Lx_hi - bppi->x < cellsize) nmirror++;
		if(bppi->y - Ly_lo < cellsize) nmirror++;
		if(Ly_hi - bppi->y < cellsize) nmirror++;
	}

	if(VORO_NPAD(simpar) != (mpad + nmirror)){
		VORORK4_TBPP(simpar) = (treevorork4particletype*)
			realloc(VORORK4_TBPP(simpar), p_size * (mpad + nmirror));
	}
	bpp_raw = (char*)VORORK4_TBPP(simpar);
	char *cursor = bpp_raw + mpad*p_size;

	#define _SEDOV_MIRROR_X(SRC,XW,SIGN) do { \
		treevorork4particletype *gp = (treevorork4particletype*)cursor; \
		memcpy(cursor, (SRC), p_size); \
		gp->u4if.indx = MAX_INDEX; \
		gp->x = 2.0*(XW) - (SRC)->x; \
		gp->vx = -(SRC)->vx; \
		CLEAR_FLAG(gp); SET_FLAG(gp, BoundaryGhostflag); SET_FLAG(gp, VOROflag); \
		cursor += p_size; \
	} while(0)
	#define _SEDOV_MIRROR_Y(SRC,YW,SIGN) do { \
		treevorork4particletype *gp = (treevorork4particletype*)cursor; \
		memcpy(cursor, (SRC), p_size); \
		gp->u4if.indx = MAX_INDEX; \
		gp->y = 2.0*(YW) - (SRC)->y; \
		gp->vy = -(SRC)->vy; \
		CLEAR_FLAG(gp); SET_FLAG(gp, BoundaryGhostflag); SET_FLAG(gp, VOROflag); \
		cursor += p_size; \
	} while(0)

	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		if(bpi->x - Lx_lo < cellsize) _SEDOV_MIRROR_X(bpi, Lx_lo, -1);
		if(Lx_hi - bpi->x < cellsize) _SEDOV_MIRROR_X(bpi, Lx_hi, +1);
		if(bpi->y - Ly_lo < cellsize) _SEDOV_MIRROR_Y(bpi, Ly_lo, -1);
		if(Ly_hi - bpi->y < cellsize) _SEDOV_MIRROR_Y(bpi, Ly_hi, +1);
	}
	for(i=0;i<mpad;i++){
		treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
		if(bppi->x - Lx_lo < cellsize) _SEDOV_MIRROR_X(bppi, Lx_lo, -1);
		if(Lx_hi - bppi->x < cellsize) _SEDOV_MIRROR_X(bppi, Lx_hi, +1);
		if(bppi->y - Ly_lo < cellsize) _SEDOV_MIRROR_Y(bppi, Ly_lo, -1);
		if(Ly_hi - bppi->y < cellsize) _SEDOV_MIRROR_Y(bppi, Ly_hi, +1);
	}
	#undef _SEDOV_MIRROR_X
	#undef _SEDOV_MIRROR_Y
	VORO_NPAD(simpar) = mpad + nmirror;

	bpp_raw = (char*)VORORK4_TBPP(simpar);
	for(i=0;i<VORO_NPAD(simpar);i++){
		treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
		int ix = ((bppi->x-xmin)/cellsize);
		int iy = ((bppi->y-ymin)/cellsize);
		int index = ix+mx*iy;
		if(index >= mx*my || index < 0) continue;
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)bppi;
		cells[index].nmem++;
		bppi->next = tmp;
	}
}

void mkLinkedList2D_oExam(
		SimParameters *simpar,
		postype cellsize,
		postype xmin,
		postype ymin,
		postype xmax,
		postype ymax,
		void (*paddingTreeAllParticle)(SimParameters *, postype)
		){
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	int np = VORO_NP(simpar);
    int i,j,k;
    int mx,my;

    mx = BASICCELL_MX(simpar);
    my = BASICCELL_MY(simpar);

	CellType *cells = VORO_BASICCELL(simpar);

	for(i=0;i<mx*my;i++) {
        cells[i].link = NULL;
        cells[i].nmem = 0;
    }
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
        treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
        int ix,iy;
        ix = ((bpi->x-xmin)/cellsize);
        iy = ((bpi->y-ymin)/cellsize);
        int index = ix+mx*iy;
        if(index >=mx*my || index < 0){
            DEBUGPRINT("error detected: %d : %d %d : %g %g %g %g\n",
                    index, mx,my, bpi->x, bpi->y, xmin,ymin);
            continue;
        }
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)bpi;
        cells[index].nmem ++;
        bpi->next = tmp;
    }
	// now for the particle padding from other ranks
    paddingTreeAllParticle(simpar, cellsize);
    char *bpp_raw = (char*)VORORK4_TBPP(simpar);
    for(i=0;i<VORO_NPAD(simpar);i++){
        treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
        int ix,iy;
        ix = ((bppi->x-xmin)/cellsize);
        iy = ((bppi->y-ymin)/cellsize);
        int index = ix+mx*iy;
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)bppi;
        cells[index].nmem ++;
        bppi->next = tmp;
    }
}



void updateDenW2Pressure2D(
		SimParameters *simpar, 
		postype xmin, postype ymin, postype xmax, postype ymax,
		 postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		postype Dtime
		){
	postype boxsize = VORO_CELLBOX(simpar);
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int isave = -1;


	det2d_dpqRK4(simpar, paddingAllTreeParticles);
	int i;
	/*
	for(i=0;i<nbp;i++){
		bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
	}
	*/
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)bp;
	for(i=0;i<nbp;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		if(GAS_Kappa(simpar) <0) { 
			bpi->w2 = -GAS_Kappa(simpar); 
		} 
		else if (GAS_Kappa(simpar) >0){ 
			bpi->w2hydro = getw2forHydroParticle(simpar, bpi, Dtime);
			applyW2Controls(simpar, bpi, Dtime);
		}
	}
	/*
	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	*/
	postype cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar);
    int mx, my;
    BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
    BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
//	DEBUGPRINT("P%d has mxy = %d %d cellsize= %g\n", MYID(simpar), mx,my,cellsize);


	CellType *cells = (VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*mx*my));
	mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

    int iy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(iy=0;iy<my;iy++){
    	int mp=1000;
    	Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
    	postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
        for(ix=0;ix<mx;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i); 
                center.csound = p[i].csound;
                center.w2 = p[i].w2;


                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
				treevorork4particletype *ibp = p[i].bp;
//                ibp->volume = Area2DPolygon(vorocorner, mp);
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);
                ibp->den = ibp->mass/ibp->volume;
				/*
                ibp->csound = sqrt(Gamma*ibp->pressure/ibp->den);
                ibp->pressure = ibp->ie/ibp->volume*(Gamma-1);
				*/
            }
            if(np>0) free(p);
			if(nneigh>0) {
				free(neighbors); free(neighwork);
			}

        }
        free(vorocorner);
    }
	for(i=0;i<nbp;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		bpi->pressure = bpi->ie/bpi->volume*(Gamma-1);
		bpi->csound = sqrt(Gamma*bpi->pressure/bpi->den);
	}
	// finalize mkLinkedList2D by my_free memory spaces (cell & boundary ghost particles)
	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
}

/* ================================================================
 *  updateDenW2Pressure2D_kNN:
 *  k-NN based density/pressure computation. No cell grid needed.
 *  Uses pre-built tree from buildTree2D.
 * ================================================================ */
void updateDenW2Pressure2D_kNN(SimParameters *simpar,
		postype xmin, postype ymin, postype xmax, postype ymax, postype Gamma,
		TStruct *tree, TPtlStruct *ptl, int nptl, postype Dtime,
		int (*skipBP)(treevorork4particletype *)){
	postype boxsize = VORO_CELLBOX(simpar);
	int nbp = VORO_NP(simpar);
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int i;

	/* Debug: check ie at function entry (disabled for performance) */

	/* Update w2 from kappa setting */
	for(i=0;i<nbp;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		if(GAS_Kappa(simpar) <0){
			bpi->w2 = -GAS_Kappa(simpar);
		}
		else if(GAS_Kappa(simpar) >0){
			bpi->w2hydro = getw2forHydroParticle(simpar, bpi, Dtime);
			applyW2Controls(simpar, bpi, Dtime);
		}
	}

	/* Debug: ie after w2 loop (disabled for performance) */

	int K_START = 32;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<nbp;i++){
		int skip_this = 0;
		if(skipBP){
			treevorork4particletype *bpi_chk = (treevorork4particletype*)(bp_raw + i*p_size);
			if(skipBP(bpi_chk)) skip_this = 1;
		}
		if(!skip_this){

		int mp = 1000;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		Voro2D_point neighbors[MAX_NUM_NEAR];
		Voro2D_point neighwork[MAX_NUM_NEAR];
		memset(neighwork, 0, sizeof(neighwork));
		PosType maxr;
		int K = K_START;
		int nfound, ip;

		/* Adaptive kNN with retry (no goto/continue — Intel OpenMP safe) */
		int knn_ok = 0;
		while(!knn_ok){
			nfound = searchKNN2D(ptl+i, K, tree, neighbors, &maxr);
			if(nfound < 3 || maxr < 1e-20) break;

			Voro2D_point center;
			treevorork4particletype *ibp_tmp = (treevorork4particletype*)(ptl[i].bp);
			center.x = ibp_tmp->x;
			center.y = ibp_tmp->y;
			center.indx = PINDX(ibp_tmp);
			center.csound = ibp_tmp->csound;
			center.w2 = ibp_tmp->w2;

			ip = Voro2D_FindVC(&center, neighbors, neighwork, nfound, vorocorner, mp, boxsize);

			/* Completeness check */
			Voro2D_Corner *vtmp = vorocorner;
			postype maxvert2 = 0;
			do {
				postype dx = vtmp->x - center.x;
				postype dy = vtmp->y - center.y;
				postype d2 = dx*dx + dy*dy;
				if(d2 > maxvert2) maxvert2 = d2;
				vtmp = vtmp->upperlink;
			} while(vtmp != vorocorner);
			if(sqrt(maxvert2) > 0.95*maxr && K < MAX_NUM_NEAR){
				K = MIN(K*2, MAX_NUM_NEAR);
			} else {
				knn_ok = 1;
			}
		}

		if(knn_ok){
		treevorork4particletype *ibp = (treevorork4particletype*)(ptl[i].bp);

		get2dAreaAvgNeighorPressure(ibp, vorocorner, neighwork, bp);
		ibp->den = ibp->mass / ibp->volume;

		/* Debug: check particles with large volume */
		if(ibp->volume > 0.01){
			/* Count corners and big-box faces */
			int ncorners = 0, nbigbox = 0;
			Voro2D_Corner *dtmp = vorocorner;
			do {
				ncorners++;
				if(dtmp->upperrelated < 0) nbigbox++;
				dtmp = dtmp->upperlink;
			} while(dtmp != vorocorner);
/*			DEBUGPRINT("P%d kNN big cell: PINDX=%ld pos=(%g,%g) vol=%g nfound=%d K=%d maxr=%g ncorners=%d nbigbox=%d ip=%d\n",
				MYID(simpar), PINDX(ibp), ibp->x, ibp->y, ibp->volume,
				nfound, K, maxr, ncorners, nbigbox, ip); */
		}
		} /* if(knn_ok) */

		free(vorocorner);
		} /* if(!skip_this) */
	}

	/* Debug: check ie after kNN Voronoi loop */
	{
		int bad = 0;
		for(i=0;i<nbp;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			if(bpi->ie < 0){
				if(bad < 3)
					DEBUGPRINT("P%d kNN POST-VORO bad ie: i=%d PINDX=%ld ie=%g vol=%g pos=(%g,%g)\n",
						MYID(simpar), i, PINDX(bpi), bpi->ie, bpi->volume, bpi->x, bpi->y);
				bad++;
			}
		}
		if(bad > 0) DEBUGPRINT("P%d kNN POST-VORO: %d particles with negative ie\n", MYID(simpar), bad);
	}

	/* Global pressure/csound update (stride-aware for av_mode >= 1) */
	{
		int neg_ie = 0, neg_vol = 0, has_wall = 0;
		for(i=0;i<nbp;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			if(bpi->ie < 0 && neg_ie < 3){
				DEBUGPRINT("P%d kNN pre-update: i=%d PINDX=%ld ie=%g vol=%g mass=%g pos=(%g,%g)\n",
						MYID(simpar), i, PINDX(bpi), bpi->ie, bpi->volume, bpi->mass,
						bpi->x, bpi->y);
				neg_ie++;
			}
			if(bpi->volume <= 0 && neg_vol < 3){
				DEBUGPRINT("P%d kNN bad vol: i=%d PINDX=%ld vol=%g\n",
						MYID(simpar), i, PINDX(bpi), bpi->volume);
				neg_vol++;
			}
		}
	}
	for(i=0;i<nbp;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		bpi->pressure = bpi->ie/bpi->volume*(Gamma-1);
		if(bpi->pressure <= 0){
			static int pfloor_warn = 0;
			if(pfloor_warn < 10)
				DEBUGPRINT("P%d P-floor: i=%d PINDX=%ld ie=%g vol=%g P=%g pos=(%g,%g)\n",
					MYID(simpar), i, PINDX(bpi), bpi->ie, bpi->volume, bpi->pressure,
					bpi->x, bpi->y);
			pfloor_warn++;
			postype p_old = bpi->csound * bpi->csound * bpi->den / Gamma;
			postype p_floor = (p_old > 0) ? 0.01 * p_old : 1e-6;
			if(p_floor < 1e-6) p_floor = 1e-6;
			bpi->pressure = p_floor;
			bpi->ie = p_floor * bpi->volume / (Gamma-1);
		}
		bpi->csound = sqrt(Gamma*bpi->pressure/bpi->den);
		if(isnan(bpi->csound)){
			DEBUGPRINT("P%d kNN_den NaN at i=%d PINDX=%ld: ie=%g vol=%g den=%g P=%g mass=%g pos=(%g,%g)\n",
					MYID(simpar), i, PINDX(bpi), bpi->ie, bpi->volume, bpi->den, bpi->pressure,
					bpi->mass, bpi->x, bpi->y);
		}
	}
}

/* ================================================================
 *  Voro2D_SutherlandHodgman_CPU:
 *  Sutherland-Hodgman polygon clipping tessellation (array-based).
 *  Same algorithm as the GPU kernel, ported to plain C for benchmarking.
 *  Input: center (x,y,w2), neighbors (relative coords, sorted by drad2).
 *  Output: volume.  out_nfaces = number of polygon faces.
 * ================================================================ */
#ifdef TESS_BENCHMARK
#define SH_MAX_CORNERS 64
static double Voro2D_SutherlandHodgman_CPU(
	double w2c,
	Voro2D_point *neighbors, int nneigh,
	double boxsize, int *out_nfaces)
{
	double halfbox = 0.5 * boxsize;
	double poly_x[SH_MAX_CORNERS], poly_y[SH_MAX_CORNERS];
	double tmp_x[SH_MAX_CORNERS], tmp_y[SH_MAX_CORNERS];
	int nc = 4;
	poly_x[0] = -halfbox; poly_y[0] = -halfbox;
	poly_x[1] =  halfbox; poly_y[1] = -halfbox;
	poly_x[2] =  halfbox; poly_y[2] =  halfbox;
	poly_x[3] = -halfbox; poly_y[3] =  halfbox;

	double maxdist2 = halfbox * sqrt(3.0);

	for (int ni = 0; ni < nneigh; ni++) {
		if (neighbors[ni].drad2 > maxdist2) break;

		double nx_n = neighbors[ni].x, ny_n = neighbors[ni].y;
		double dist2 = nx_n*nx_n + ny_n*ny_n;
		double wfrac = 0.5 + 0.5*(w2c - neighbors[ni].w2)/dist2;
		if(wfrac < 0.05) wfrac = 0.05; else if(wfrac > 0.95) wfrac = 0.95;
		double threshold = wfrac * dist2;

		int out_count = 0;
		for (int e = 0; e < nc && out_count < SH_MAX_CORNERS-2; e++) {
			int e_next = (e+1) % nc;
			double ax = poly_x[e],    ay = poly_y[e];
			double bx = poly_x[e_next], by = poly_y[e_next];
			double da = ax*nx_n + ay*ny_n;
			double db = bx*nx_n + by*ny_n;
			int a_in = (da <= threshold);
			int b_in = (db <= threshold);

			if (a_in && b_in) {
				tmp_x[out_count] = ax;
				tmp_y[out_count] = ay;
				out_count++;
			} else if (a_in && !b_in) {
				tmp_x[out_count] = ax;
				tmp_y[out_count] = ay;
				out_count++;
				double denom = db - da;
				double t = (fabs(denom) > 1e-30) ? (threshold - da)/denom : wfrac;
				if (t < 1e-9) t = 1e-9;
				if (t > 1.0 - 1e-9) t = 1.0 - 1e-9;
				tmp_x[out_count] = ax + t*(bx-ax);
				tmp_y[out_count] = ay + t*(by-ay);
				out_count++;
			} else if (!a_in && b_in) {
				double denom = db - da;
				double t = (fabs(denom) > 1e-30) ? (threshold - da)/denom : wfrac;
				if (t < 1e-9) t = 1e-9;
				if (t > 1.0 - 1e-9) t = 1.0 - 1e-9;
				tmp_x[out_count] = ax + t*(bx-ax);
				tmp_y[out_count] = ay + t*(by-ay);
				out_count++;
			}
		}
		nc = out_count;
		for (int e = 0; e < nc; e++) {
			poly_x[e] = tmp_x[e]; poly_y[e] = tmp_y[e];
		}
		if (nc < 3) break;

		double new_max = 0;
		for (int e = 0; e < nc; e++) {
			double cd = poly_x[e]*poly_x[e] + poly_y[e]*poly_y[e];
			if (cd > new_max) new_max = cd;
		}
		maxdist2 = new_max;
	}

	/* Vertex merging (same as GPU kernel) */
	{
		double max_edge2 = 0;
		for (int e = 0; e < nc; e++) {
			int en = (e+1) % nc;
			double dx = poly_x[en]-poly_x[e], dy = poly_y[en]-poly_y[e];
			double e2 = dx*dx + dy*dy;
			if (e2 > max_edge2) max_edge2 = e2;
		}
		double merge_eps2 = max_edge2 * 1e-6;
		if (merge_eps2 < 1e-30) merge_eps2 = 1e-30;
		int new_nc = 0;
		for (int e = 0; e < nc; e++) {
			int en = (e+1) % nc;
			double dx = poly_x[en]-poly_x[e], dy = poly_y[en]-poly_y[e];
			if (dx*dx + dy*dy > merge_eps2) {
				tmp_x[new_nc] = poly_x[e];
				tmp_y[new_nc] = poly_y[e];
				new_nc++;
			}
		}
		if (new_nc >= 3 && new_nc < nc) {
			nc = new_nc;
			for (int e = 0; e < nc; e++) {
				poly_x[e] = tmp_x[e]; poly_y[e] = tmp_y[e];
			}
		}
	}

	/* Shoelace volume */
	double area = 0;
	for (int e = 0; e < nc; e++) {
		int en = (e+1) % nc;
		area += poly_x[e]*poly_y[en] - poly_y[e]*poly_x[en];
	}
	*out_nfaces = nc;
	return 0.5 * fabs(area);
}
#undef SH_MAX_CORNERS
#endif /* TESS_BENCHMARK */

/* ================================================================
 *  updateDenW2Pressure2DBlend:
 *  Same as updateDenW2Pressure2D but for treevorostressrk4particletype.
 *  Also computes velocity gradient (∇⊗v) via Gauss divergence theorem
 *  and NS stress tensor τ per cell when av_mode >= 1.
 * ================================================================ */
/* Volume-constrained Laguerre weights, no Riemann solve.
 * Target areas follow dV/dt = -V div v from the previous geometric
 * velocity gradient. Weights are a damped diagonal Newton step so the
 * Laguerre cells realize those areas. A constant added to every w2 is
 * removed; it does not move the radical planes. */
static int sedov_lagvol_on(void){
	static int cached = -1;
	if(cached < 0){
		const char *e = getenv("SEDOV_LAGVOL");
		cached = (e && e[0]=='1') ? 1 : 0;
	}
	return cached;
}

static double *lag_vt, *lag_sw;
static size_t lag_vt_cap;
static int lag_vt_ready, lag_sw_ready;

static double *lag_vt_at(size_t indx){
	if(indx >= lag_vt_cap){
		size_t ncap = lag_vt_cap ? lag_vt_cap : 1024;
		while(ncap <= indx) ncap *= 2;
		double *grown = (double*)realloc(lag_vt, ncap*sizeof(double));
		double *sgrown = (double*)realloc(lag_sw, ncap*sizeof(double));
		if(!grown || !sgrown) return NULL;
		memset(grown + lag_vt_cap, 0, (ncap - lag_vt_cap)*sizeof(double));
		memset(sgrown + lag_vt_cap, 0, (ncap - lag_vt_cap)*sizeof(double));
		lag_vt = grown;
		lag_sw = sgrown;
		lag_vt_cap = ncap;
	}
	return lag_vt + indx;
}

static void laguerre_volume_newton(
		SimParameters *simpar,
		postype xmin, postype ymin, postype xmax, postype ymax,
		postype boxsize, postype cellsize, int mx, int my,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int, int, int *),
		void mkLinkedList2D(SimParameters *, postype, postype, postype, postype, postype,
			void (*)(SimParameters *, postype)),
		postype Dtime, int advance)
{
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype dx = SIMBOX(simpar).x.max / NX(simpar);
	(void)dx;
	/* One corrector: ω_pred = ω^n, then a single Laplacian step. */
	int iter, niter = 1;
	size_t max_indx = 0;
	for(int i=0;i<nbp;i++){
		size_t id = (size_t)PINDX(bp+i);
		if(id > max_indx) max_indx = id;
	}
	if(!lag_vt_at(max_indx)) return;

	if(lag_vt_ready && advance && lag_sw_ready){
		double before = 0, after = 0;
		for(int i=0;i<nbp;i++){
			size_t id = (size_t)PINDX(bp+i);
			double *vt = lag_vt_at(id);
			if(!vt || *vt <= 0) continue;
			before += *vt;
			double dV = lag_sw[id] * (double)Dtime;
			if(!isfinite(dV)) dV = 0;
			*vt += dV;
			if(*vt < 1e-8) *vt = 1e-8;
			after += *vt;
		}
		if(after > 0 && before > 0){
			double scale = before / after;
			for(int i=0;i<nbp;i++){
				double *vt = lag_vt_at((size_t)PINDX(bp+i));
				if(vt && *vt > 0) *vt *= scale;
			}
		}
		lag_sw_ready = 0;
	}

	int *loc_of = (int*)malloc(lag_vt_cap * sizeof(int));
	for(size_t k=0;k<lag_vt_cap;k++) loc_of[k] = -1;
	for(int i=0;i<nbp;i++) loc_of[(size_t)PINDX(bp+i)] = i;
	typedef struct { int i, j; double h; } LagEdge;
	LagEdge *edges = (LagEdge*)malloc((size_t)nbp * 24 * sizeof(LagEdge));
	double *rhs = (double*)calloc((size_t)nbp, sizeof(double));
	double *dw = (double*)calloc((size_t)nbp, sizeof(double));
	double *cg_r = (double*)calloc((size_t)nbp, sizeof(double));
	double *cg_p = (double*)calloc((size_t)nbp, sizeof(double));
	double *cg_ap = (double*)calloc((size_t)nbp, sizeof(double));
	double *dmin2 = (double*)malloc((size_t)nbp * sizeof(double));

	for(iter=0; iter<niter; iter++){
		CellType *cells = (VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*mx*my));
		mkLinkedList2D(simpar, cellsize, xmin, ymin, xmax, ymax, paddingAllTreeParticles);
		double err = 0, sum_v = 0, sum_t = 0;
		int nedge = 0, nbad = 0;
		double max_ratio = 0;
		memset(rhs, 0, (size_t)nbp*sizeof(double));
		memset(dmin2, 0, (size_t)nbp*sizeof(double));
		for(int iy=0; iy<my; iy++){
			int mp = 1000;
			Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
			for(int ix=0; ix<mx; ix++){
				int np;
				treevorork4particletype *p = find2DCellBP(simpar, ix, iy, &np);
				int nneigh;
				Voro2D_point *neighbors = find2DNeighborBP(simpar, ix, iy, &nneigh);
				Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
				for(int i=0;i<np;i++){
					Voro2D_point center;
					center.x = p[i].x; center.y = p[i].y;
					center.indx = PINDX(p+i);
					center.csound = p[i].csound;
					center.w2 = p[i].w2;
					Voro2D_FindVC(&center, neighbors, neighwork, nneigh, vorocorner, mp, boxsize);
					treevorostressrk4particletype *ibp =
						(treevorostressrk4particletype*)p[i].bp;
					postype area = Area2DPolygon(vorocorner, mp);
					if(!(area > 0) || !isfinite(area)){
						nbad++;
						ibp->w2 = 0;
						ibp->rk4.w2backup = 0;
						continue;
					}
					int li = loc_of[(size_t)PINDX(ibp)];
					postype dmin = 1e30;
					Voro2D_Corner *tmp = vorocorner;
					do {
						if(tmp->upperrelated >= 0 && li >= 0){
							postype dxe = tmp->x - tmp->upperlink->x;
							postype dye = tmp->y - tmp->upperlink->y;
							postype ds = sqrt(dxe*dxe + dye*dye);
							treevorork4particletype *jbp = neighwork[tmp->upperrelated].bp;
							size_t jid = (size_t)PINDX(jbp);
							postype ddx = jbp->x - ibp->x;
							postype ddy = jbp->y - ibp->y;
							postype dist = sqrt(ddx*ddx + ddy*ddy);
							if(dist > 1e-12){
								if(dist < dmin) dmin = dist;
								double h = (double)ds / (2.0*(double)dist);
								int lj = (jid < lag_vt_cap) ? loc_of[jid] : -1;
								if(lj < 0){
									if(nedge < nbp*24)
										edges[nedge++] = (LagEdge){li, -1, h};
								} else if(lj > li){
									if(nedge < nbp*24)
										edges[nedge++] = (LagEdge){li, lj, h};
								}
							}
						}
						tmp = tmp->upperlink;
					} while(tmp != vorocorner);
					if(li >= 0) dmin2[li] = (dmin < 1e29) ? (double)(dmin*dmin) : 0;
					double *vt = lag_vt + (size_t)PINDX(ibp);
					if(!lag_vt_ready){
						if(vt) *vt = (double)area;
					} else if(li >= 0 && vt){
						rhs[li] = *vt - (double)area;
						err += fabs(rhs[li]);
						sum_t += *vt;
					}
					ibp->volume = area;
					sum_v += (double)area;
				}
				free(neighwork); free(neighbors); free(p);
			}
			free(vorocorner);
		}
		my_free(VORO_BASICCELL(simpar));
		VORO_BASICCELL(simpar) = NULL;
		my_free(VORORK4_TBPP(simpar));
		VORORK4_TBPP(simpar) = NULL;
		VORO_NPAD(simpar) = 0;
		if(nbad > 0){
			for(int i=0;i<nbp;i++){
				bp[i].w2 = 0;
				bp[i].rk4.w2backup = 0;
			}
			if(MYID(simpar)==0){
				fprintf(stderr, "[LAGVOL] reject nbad=%d, back to Voronoi\n", nbad);
				fflush(stderr);
			}
			break;
		}
		if(!lag_vt_ready){
			lag_vt_ready = 1;
			if(MYID(simpar)==0){
				fprintf(stderr, "[LAGVOL] init sumV=%g\n", sum_v);
				fflush(stderr);
			}
			break;
		}
		/* H δw = V_target - V. H is the face Laplacian, kernel = constants. */
		double mean = 0;
		for(int i=0;i<nbp;i++) mean += rhs[i];
		mean /= (nbp > 0 ? nbp : 1);
		for(int i=0;i<nbp;i++){ rhs[i] -= mean; dw[i] = 0; cg_r[i] = rhs[i]; cg_p[i] = rhs[i]; }
		double rr = 0;
		for(int i=0;i<nbp;i++) rr += cg_r[i]*cg_r[i];
		double rr0 = rr;
		for(int k=0; k<40 && rr > 1e-18 && rr > 1e-10*rr0; k++){
			for(int i=0;i<nbp;i++) cg_ap[i] = 0;
			for(int e=0;e<nedge;e++){
				int a = edges[e].i, b = edges[e].j;
				double h = edges[e].h;
				cg_ap[a] += h * cg_p[a];
				if(b >= 0){
					cg_ap[b] += h * cg_p[b];
					cg_ap[a] -= h * cg_p[b];
					cg_ap[b] -= h * cg_p[a];
				}
			}
			double denom = 0;
			for(int i=0;i<nbp;i++) denom += cg_p[i]*cg_ap[i];
			if(!(fabs(denom) > 1e-30)) break;
			double alpha = rr / denom;
			for(int i=0;i<nbp;i++){
				dw[i] += alpha * cg_p[i];
				cg_r[i] -= alpha * cg_ap[i];
			}
			double rr2 = 0;
			for(int i=0;i<nbp;i++) rr2 += cg_r[i]*cg_r[i];
			double beta = rr2 / rr;
			for(int i=0;i<nbp;i++) cg_p[i] = cg_r[i] + beta * cg_p[i];
			rr = rr2;
		}
		mean = 0;
		for(int i=0;i<nbp;i++) mean += dw[i];
		mean /= (nbp > 0 ? nbp : 1);
		int ncap = 0;
		for(int i=0;i<nbp;i++){
			double wnew = (double)bp[i].w2 + (dw[i] - mean);
			double cap = dmin2[i] > 0 ? dmin2[i] : 0;
			if(wnew > cap){ wnew = cap; ncap++; }
			if(wnew < 0) wnew = 0;
			if(!isfinite(wnew)) wnew = 0;
			bp[i].w2 = (postype)wnew;
			bp[i].rk4.w2backup = bp[i].w2;
			if(cap > 0){
				double ratio = wnew / cap;
				if(ratio > max_ratio) max_ratio = ratio;
			}
		}
		if(MYID(simpar)==0){
			fprintf(stderr, "[LAGVOL] iter=%d L1=%g sumV=%g sumT=%g ncap=%d max|w2|/dmin2=%g\n",
				iter, err, sum_v, sum_t, ncap, max_ratio);
			fflush(stderr);
		}
	}
	free(loc_of); free(edges); free(rhs); free(dw); free(cg_r); free(cg_p); free(cg_ap); free(dmin2);
}

static int sedov_wsmooth_on(double *kappa)
{
	const char *e = getenv("SEDOV_WSMOOTH");
	if(!e || !e[0]) return 0;
	*kappa = atof(e);
	if(*kappa <= 0) *kappa = 1.0;
	return 1;
}

/* Voronoi volumes are already on the particles. Shift a face only when the
 * two cells are opposite volume extrema (checkerboard). A one-sided jump,
 * such as a shock, is left at the midpoint. Weights are the least-squares
 * Laguerre diagram of those face shifts. */
static void sedov_wsmooth_apply(
		SimParameters *simpar,
		postype xmin, postype ymin, postype xmax, postype ymax,
		postype boxsize, postype cellsize, int mx, int my,
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		postype Gamma)
{
	double kappa;
	if(!sedov_wsmooth_on(&kappa)) return;
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	const int M = 24;
	int *nn = (int*)calloc((size_t)nbp, sizeof(int));
	int *nei = (int*)malloc((size_t)nbp * M * sizeof(int));
	double *flen = (double*)malloc((size_t)nbp * M * sizeof(double));
	double *dd = (double*)malloc((size_t)nbp * M * sizeof(double));
	double *d2 = (double*)malloc((size_t)nbp * M * sizeof(double));
	double *bb = (double*)calloc((size_t)nbp * M, sizeof(double));
	signed char *flag = (signed char*)calloc((size_t)nbp, 1);
	double *w = (double*)calloc((size_t)nbp, sizeof(double));
	double *V = (double*)malloc((size_t)nbp * sizeof(double));
	if(!nn || !nei || !flen || !dd || !d2 || !bb || !flag || !w || !V){
		free(nn); free(nei); free(flen); free(dd); free(d2); free(bb); free(flag); free(w); free(V);
		return;
	}
	for(int i=0;i<nbp;i++) V[i] = bp[i].volume;

	for(int iy=0; iy<my; iy++){
		int mp = 4000;
		Voro2D_Corner *vc = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		for(int ix=0; ix<mx; ix++){
			int np, nneigh;
			treevorork4particletype *p = find2DCellBP(simpar, ix, iy, &np);
			Voro2D_point *neighbors = find2DNeighboringBP(simpar, ix, iy, &nneigh);
			Voro2D_point *nw = (nneigh>0) ? (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh) : NULL;
			for(int i=0;i<np;i++){
				int ii = (int)((treevorostressrk4particletype*)p[i].bp - bp);
				if(ii < 0 || ii >= nbp) continue;
				Voro2D_point center;
				center.x = p[i].x; center.y = p[i].y;
				center.indx = PINDX(p+i);
				center.csound = p[i].csound;
				center.w2 = 0;
				int ip = Voro2D_FindVC(&center, neighbors, nw, nneigh, vc, mp, boxsize);
				if(ip <= 0) continue;
				Voro2D_Corner *c = vc;
				int guard = 0;
				do {
					int ur = c->upperrelated;
					if(ur >= 0 && ur < nneigh && nn[ii] < M){
						treevorostressrk4particletype *jp = (treevorostressrk4particletype*)nw[ur].bp;
						int jj = (int)(jp - bp);
						if(jj >= 0 && jj < nbp && jj != ii){
							int seen = 0;
							for(int t=0;t<nn[ii];t++) if(nei[ii*M+t]==jj) seen = 1;
							if(seen) { c = c->upperlink; if(++guard > mp) break; continue; }
							Voro2D_Corner *c2 = c->upperlink;
							double L = hypot(c2->x - c->x, c2->y - c->y);
							double dx = jp->x - bp[ii].x, dy = jp->y - bp[ii].y;
							double dist = hypot(dx, dy);
							int slot = ii*M + nn[ii];
							nei[slot] = jj;
							flen[slot] = L;
							dd[slot] = dist;
							d2[slot] = dist*dist;
							nn[ii]++;
						}
					}
					c = c->upperlink;
					if(++guard > mp) break;
				} while(c && c != vc);
			}
			if(np>0) free(p);
			if(nneigh>0){ free(neighbors); if(nw) free(nw); }
		}
		free(vc);
	}

	for(int i=0;i<nbp;i++){
		if(nn[i] <= 0) continue;
		int is_max = 1, is_min = 1;
		for(int k=0;k<nn[i];k++){
			double Vn = V[nei[i*M+k]];
			if(!(V[i] > Vn)) is_max = 0;
			if(!(V[i] < Vn)) is_min = 0;
		}
		if(is_max) flag[i] = 1;
		else if(is_min) flag[i] = -1;
	}

	char *tagged = (char*)calloc((size_t)nbp * M, 1);
	int *tdeg = (int*)calloc((size_t)nbp, sizeof(int));
	int ncheck = 0;
	double rsum = 0;
	if(!tagged || !tdeg){
		free(tagged); free(tdeg);
		free(nn); free(nei); free(flen); free(dd); free(d2); free(bb);
		free(flag); free(w); free(V);
		return;
	}
	for(int i=0;i<nbp;i++){
		for(int k=0;k<nn[i];k++){
			int j = nei[i*M+k];
			if(j <= i) continue;
			double dist = dd[i*M+k];
			double L = flen[i*M+k];
			if(L < 1e-15 || dist < 1e-15) continue;
			double nx = (bp[j].x - bp[i].x) / dist;
			double ny = (bp[j].y - bp[i].y) / dist;
			double wm = 0, numm = 0, wp = 0, nump = 0;
			for(int k2=0;k2<nn[i];k2++){
				int kk = nei[i*M+k2];
				if(kk == j) continue;
				double dk = dd[i*M+k2];
				if(dk < 1e-15) continue;
				double xh = (bp[kk].x - bp[i].x) / dk;
				double yh = (bp[kk].y - bp[i].y) / dk;
				double wt = -(nx*xh + ny*yh);
				if(wt < 0) wt = 0;
				wm += wt;
				numm += wt * (V[i] - V[kk]);
			}
			for(int k2=0;k2<nn[j];k2++){
				int mm = nei[j*M+k2];
				if(mm == i) continue;
				double dk = dd[j*M+k2];
				if(dk < 1e-15) continue;
				double xh = (bp[mm].x - bp[j].x) / dk;
				double yh = (bp[mm].y - bp[j].y) / dk;
				double wt = nx*xh + ny*yh;
				if(wt < 0) wt = 0;
				wp += wt;
				nump += wt * (V[mm] - V[j]);
			}
			if(wm < 0.2 || wp < 0.2) continue;
			double gm = numm / wm, gpv = nump / wp;
			double dV = V[j] - V[i];
			if(!(dV*gm < 0 && dV*gpv < 0)) continue;
			/* A shock has a large pressure ratio. A contact is nearly isobaric.
			 * Only the isobaric oscillation may move the face. */
			{
				double pi = bp[i].pressure, pj = bp[j].pressure;
				double pmin = pi < pj ? pi : pj;
				double pmax = pi > pj ? pi : pj;
				if(!(pmin > 0) || pmax > 2.0 * pmin) continue;
			}
			double mag = fabs(dV);
			if(mag > fabs(gm)) mag = fabs(gm);
			if(mag > fabs(gpv)) mag = fabs(gpv);
			double vsmall = V[i] < V[j] ? V[i] : V[j];
			if(mag > 0.25 * vsmall) mag = 0.25 * vsmall;
			double shift = 0.5 * kappa * (dV > 0 ? mag : -mag) / L;
			double smax = 0.45 * 0.5 * dist;
			if(shift > smax) shift = smax;
			if(shift < -smax) shift = -smax;
			double bij = 2.0 * dist * shift;
			bb[i*M+k] = bij;
			tagged[i*M+k] = 1;
			tdeg[i]++; tdeg[j]++;
			for(int k2=0;k2<nn[j];k2++){
				if(nei[j*M+k2] == i){
					bb[j*M+k2] = -bij;
					tagged[j*M+k2] = 1;
					break;
				}
			}
			ncheck++;
			rsum += hypot(0.5*(bp[i].x+bp[j].x) - 0.5, 0.5*(bp[i].y+bp[j].y) - 0.5);
		}
	}

	for(int it=0; it<40; it++){
		for(int i=0;i<nbp;i++){
			if(tdeg[i] <= 0) continue;
			double sum = 0;
			for(int k=0;k<nn[i];k++){
				if(!tagged[i*M+k]) continue;
				sum += w[nei[i*M+k]] + bb[i*M+k];
			}
			w[i] = sum / tdeg[i];
		}
	}
	for(int sweep=0; sweep<8; sweep++){
		for(int i=0;i<nbp;i++){
			for(int k=0;k<nn[i];k++){
				if(!tagged[i*M+k]) continue;
				int j = nei[i*M+k];
				if(j <= i) continue;
				double lim = 0.95 * d2[i*M+k];
				double diff = w[i] - w[j];
				if(diff > lim || diff < -lim){
					double mid = 0.5*(w[i] + w[j]);
					double h = 0.5 * lim;
					if(diff < 0) h = -h;
					w[i] = mid + h;
					w[j] = mid - h;
				}
			}
		}
	}
	double wmin = 0;
	for(int i=0;i<nbp;i++) if(w[i] < wmin) wmin = w[i];
	double wmax = 0;
	for(int i=0;i<nbp;i++){
		w[i] -= wmin;
		if(w[i] < 0) w[i] = 0;
		bp[i].w2 = w[i];
		if(w[i] > wmax) wmax = w[i];
	}
	if(MYID(simpar)==0){
		static int nlog = 0;
		if(nlog < 8 || (nlog % 50) == 0)
			fprintf(stderr, "[WSMOOTH] kappa=%g ncheck=%d maxw2=%g rmean=%g\n",
				kappa, ncheck, wmax, ncheck ? rsum/ncheck : 0);
		nlog++;
	}

	for(int iy=0; iy<my; iy++){
		int mp = 4000;
		Voro2D_Corner *vc = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		for(int ix=0; ix<mx; ix++){
			int np, nneigh;
			treevorork4particletype *p = find2DCellBP(simpar, ix, iy, &np);
			Voro2D_point *neighbors = find2DNeighboringBP(simpar, ix, iy, &nneigh);
			Voro2D_point *nw = (nneigh>0) ? (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh) : NULL;
			for(int k=0;k<nneigh;k++){
				treevorostressrk4particletype *jp = (treevorostressrk4particletype*)neighbors[k].bp;
				int jj = (int)(jp - bp);
				if(jj >= 0 && jj < nbp) neighbors[k].w2 = bp[jj].w2;
			}
			for(int i=0;i<np;i++){
				treevorork4particletype *ibp0 = p[i].bp;
				int ii = (int)((treevorostressrk4particletype*)ibp0 - bp);
				Voro2D_point center;
				center.x = p[i].x; center.y = p[i].y;
				center.indx = PINDX(p+i);
				center.csound = p[i].csound;
				center.w2 = (ii >= 0 && ii < nbp) ? bp[ii].w2 : p[i].w2;
				Voro2D_FindVC(&center, neighbors, nw, nneigh, vc, mp, boxsize);
				treevorork4particletype *ibp = p[i].bp;
				get2dAreaAvgNeighorPressure(ibp, vc, nw, (treevorork4particletype*)bp);
				if(ibp->volume > 1e-15){
					ibp->den = ibp->mass / ibp->volume;
					ibp->pressure = ibp->ie / ibp->volume * (Gamma - 1);
					if(ibp->pressure < 1e-6) ibp->pressure = 1e-6;
					ibp->csound = sqrt(Gamma * ibp->pressure / ibp->den);
				}
			}
			if(np>0) free(p);
			if(nneigh>0){ free(neighbors); if(nw) free(nw); }
		}
		free(vc);
	}
	free(nn); free(nei); free(flen); free(dd); free(d2);
	free(flag); free(w); free(V); free(bb); free(tagged); free(tdeg);
}

/* CPU stage-floor ledger (see the pressure loop of
 * updateDenW2Pressure2DBlend). Reset each RK4 step, printed in [RK4E]. */
static double gfs_stage_floor_l = 0;
static int gfs_stage_floor_n_l = 0;

/* Widest neighbour search (in linked-list cells) used to close a Voronoi/
 * Laguerre cell. A cell lives inside its cap square of side VORO_CELLBOX,
 * so only generators within twice the cap half-diagonal can own a face of
 * it; searching ceil(2 h_c / cellsize) cells fixes the cell, and one more
 * cell covers the weight shift of a Laguerre face. Directions with no
 * generator inside that radius (e.g. a void) keep cap edges with a
 * negative neighbour id, which every pass treats as vacuum (skipped).
 * EUNHA_VORO_RMAX overrides. Every pass (volume, gradient/limiter, force,
 * scheme-H solves) uses this same bound so they see one tessellation. */
/* EUNHA_VORO_RHO = rho/dx > 0 (default 0 = off): the fluid domain is the
 * union of discs B_i of radius R = rho (W = 0), and each cell is C_i ∩ B_i
 * (exact partition of the union for a common rho; Docs/consult_vacuum_boundary*.md,
 * option (B), stage 1). Faces keep only their part inside the disc; the
 * exposed arcs become single arc edges (VORO_ARC_EDGE, a negative neighbour
 * id, i.e. vacuum); Area2DPolygon adds the exact circular segments, and the
 * normal integral of an arc equals that of its chord. The ring is rebuilt
 * in place, so every CPU pass (volume, gradient, limiter, force) sees the
 * same restricted cell. Stage 1: W = 0 and CPU hydro only. */
double eunha_cellbox_fac(void){
	static double f = -1;
	if(f < 0){ const char *e = getenv("EUNHA_CELLBOX_FAC"); f = (e && e[0]) ? atof(e) : 5.0; if(!(f > 0)) f = 5.0; }
	return f;
}
double voro_disc_radius(SimParameters *simpar){
	static double f = -1;
	if(f < 0){ const char *e = getenv("EUNHA_VORO_RHO"); f = (e && e[0]) ? atof(e) : 0; if(!(f > 0)) f = 0; }
	if(f == 0) return 0;
	if(GAS_GPU_ENABLED(simpar)){
		fprintf(stderr, "EUNHA_VORO_RHO: the GPU path computes volumes without the disc; use gpu_enabled = 0\n");
		exit(1);
	}
	if(getenv("GFS_DYN_W2") && atoi(getenv("GFS_DYN_W2")) > 0){
		fprintf(stderr, "EUNHA_VORO_RHO: stage 1 supports W = 0 only (GFS_DYN_W2 is set)\n");
		exit(1);
	}
	double dx = (SIMMODEL(simpar) == KH) ? (double)(SIMBOX(simpar).x.max - SIMBOX(simpar).x.min)/NX(simpar)
		: (double)BOXSIZE(simpar)/NX(simpar);
	return f*dx;
}
void voro2d_restrict_disc(Voro2D_Corner *vc, int mp, double R){
	/* collect the ring */
	int n = 0, i, k;
	Voro2D_Corner *t = vc;
	double ax[512], ay[512]; int id[512];
	do { if(n < 512){ ax[n] = t->x; ay[n] = t->y; id[n] = t->upperrelated; n++; } t = t->upperlink; } while(t != vc);
	double area2 = 0;
	for(i=0;i<n;i++){ int j = (i+1)%n; area2 += ax[i]*ay[j] - ax[j]*ay[i]; }
	double orient = (area2 >= 0) ? 1.0 : -1.0;
	double R2 = R*R;
	/* kept pieces of each edge */
	double sx0[512], sy0[512], sx1[512], sy1[512]; int kid[512], nk = 0;
	for(i=0;i<n;i++){
		int j = (i+1)%n;
		double ex = ax[j]-ax[i], ey = ay[j]-ay[i];
		double qa = ex*ex + ey*ey, qb = 2*(ax[i]*ex + ay[i]*ey), qc = ax[i]*ax[i] + ay[i]*ay[i] - R2;
		double s0 = 0, s1 = 1;
		if(qa <= 0) continue;
		double disc = qb*qb - 4*qa*qc;
		if(disc <= 0) continue;
		double sq = sqrt(disc), r0 = (-qb - sq)/(2*qa), r1 = (-qb + sq)/(2*qa);
		if(r0 > s0) s0 = r0;
		if(r1 < s1) s1 = r1;
		if(s1 - s0 <= 1e-12) continue;
		sx0[nk] = ax[i] + s0*ex; sy0[nk] = ay[i] + s0*ey;
		sx1[nk] = ax[i] + s1*ex; sy1[nk] = ay[i] + s1*ey;
		kid[nk] = id[i]; nk++;
	}
	/* rebuild: each kept segment, then (if the next one does not start where
	 * it ends) one arc edge along the circle to the next segment's start */
	double ox[1100], oy[1100]; int oid[1100], m = 0;
	if(nk == 0){
		/* the disc lies inside the polygon: four arc edges */
		for(k=0;k<4;k++){ double a = orient*k*0.5*M_PI; ox[m] = R*cos(a); oy[m] = R*sin(a); oid[m] = VORO_ARC_EDGE; m++; }
	} else {
		for(k=0;k<nk && m < 1098;k++){
			int kn = (k+1)%nk;
			ox[m] = sx0[k]; oy[m] = sy0[k]; oid[m] = kid[k]; m++;
			double gx = sx0[kn] - sx1[k], gy = sy0[kn] - sy1[k];
			if(gx*gx + gy*gy < 1e-24*R2) continue;
			/* arc to the next start, in pieces of at most pi/2 so that the
			 * chord polygon always carries the ring orientation */
			double a0 = atan2(sy1[k], sx1[k]), a1 = atan2(sy0[kn], sx0[kn]);
			double sweep = orient*(a1 - a0);
			while(sweep < 0) sweep += 2*M_PI;
			while(sweep >= 2*M_PI) sweep -= 2*M_PI;
			int np_ = (int)ceil(sweep/(0.5*M_PI)); if(np_ < 1) np_ = 1;
			for(i=0;i<np_ && m < 1098;i++){
				double a = a0 + orient*sweep*(double)i/np_;
				ox[m] = (i == 0) ? sx1[k] : R*cos(a); oy[m] = (i == 0) ? sy1[k] : R*sin(a); oid[m] = VORO_ARC_EDGE; m++;
			}
		}
	}
	if(m > mp) m = mp;
	for(k=0;k<m;k++){
		vc[k].x = (postype)ox[k]; vc[k].y = (postype)oy[k];
		vc[k].status = Active; vc[k].indx = k;
		vc[k].upperrelated = oid[k];
		vc[k].lowerrelated = oid[(k-1+m)%m];
		vc[k].upperlink = vc + (k+1)%m; vc[k].lowerlink = vc + (k-1+m)%m;
		vc[k].neighlowerlink = NULL; vc[k].neighupperlink = NULL;
	}
}
int voro_expand_rmax(SimParameters *simpar){
	const char *e = getenv("EUNHA_VORO_RMAX");
	if(e && e[0] && atoi(e) >= 1) return atoi(e);
	double hc = 0.5*sqrt(2.0)*(double)VORO_CELLBOX(simpar);
	double cw = (double)BASICCELL_CELLWIDTH(simpar);
	int r = (cw > 0) ? (int)ceil(2*hc/cw) + 1 : 4;
	return r < 2 ? 2 : r;
}

/* Security-radius completeness (geometry gate, Astra design 7 item 3): a generator farther than
 * 2 R_v from x_i (R_v = max corner radius, W = 0) cannot cut cell i. The searched block of
 * linked cells [ix-R, ix+R] (clipped to the padded grid) is complete only if 2 R_v <= reach,
 * the distance from x_i to the block boundary. Active with GFS_PERIODIC_COMPLETE. */
static long voro_sec_fail = 0;
static int voro_sec_short(SimParameters *simpar, Voro2D_Corner *vc, int ip, double cx, double cy, int ix, int iy, int R){
	if(!gfs_periodic_complete()) return 0;
	double Rv = 0; int k;
	for(k = 0; k < ip; k++) if(vc[k].status == Active){ double r_ = hypot((double)vc[k].x, (double)vc[k].y); if(r_ > Rv) Rv = r_; }
	double cs = (double)BASICCELL_CELLWIDTH(simpar), x0 = (double)Xmin_HydroExam(simpar) - cs, y0 = (double)Ymin_HydroExam(simpar) - cs;
	int mx = BASICCELL_MX(simpar), my = BASICCELL_MY(simpar);
	double lox = x0 + (ix - R > 0 ? ix - R : 0)*cs, hix = x0 + (ix + R + 1 < mx ? ix + R + 1 : mx)*cs;
	double loy = y0 + (iy - R > 0 ? iy - R : 0)*cs, hiy = y0 + (iy + R + 1 < my ? iy + R + 1 : my)*cs;
	double reach = fmin(fmin(cx - lox, hix - cx), fmin(cy - loy, hiy - cy));
	return 2*Rv > reach;
}
void updateDenW2Pressure2DBlend(
		SimParameters *simpar,
		postype xmin, postype ymin, postype xmax, postype ymax,
		postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		postype Dtime
		){
	if(gfs_subzone_on()){
		if(NID(simpar) > 1 || voro_disc_radius(simpar) > 0 || GAS_Kappa(simpar) != 0 || gfs_dynw2_on()){
			fprintf(stderr, "[SUBZ] GFS_SUBZONE: one MPI rank, no EUNHA_VORO_RHO, W=0 only\n"); exit(1);
		}
		sub_ensure(simpar);
	}
	if(gfs_shared_geom_on()){
		if(NID(simpar) > 1){ fprintf(stderr, "[SHGEOM] GFS_SHARED_GEOM: one MPI rank only\n"); exit(1); }
		shf_ensure(simpar);
	}
	if(gfs_tcf_face_z() || gfs_periodic_complete()){ cen_ensure(simpar); pc_ncap = 0; pc_vsum = 0; voro_sec_fail = 0; }
	if(gfs_tcf_qbp()){
		if(gfs_lsf_order() != 2){ fprintf(stderr, "[QBP] GFS_TCF_QBP requires EUNHA_LSF_ORDER=2\n"); exit(1); }
		qbp_ensure(simpar);
	}
	postype boxsize = VORO_CELLBOX(simpar);
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int av_mode = GAS_AVMODE(simpar);
	postype nu_phys_val = GAS_VISCOSITY(simpar);
	int use_stress = (av_mode >= 1 || nu_phys_val > 0);

	postype cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar);
#ifdef USE_CUDA
	/* The GPU kernel applies the kappa rule to w2 (0 for kappa = 0); under
	 * GFS_DYN_W2 W is a state variable, which only det2d_dpqRK4 honours. */
	if(gfs_dynw2_active) det2d_dpqRK4(simpar, paddingAllTreeParticles);
	else det2d_dpqRK4_GPU(simpar, xmin, ymin, xmax, ymax,
		paddingAllTreeParticles, mkLinkedList2D);
#else
	det2d_dpqRK4(simpar, paddingAllTreeParticles);
#endif
	int i;
	if(sedov_lagvol_on()){
		static int lag_calls = 0;
		lag_calls++;
		int mx0 = ceil((xmax-xmin)/cellsize);
		int my0 = ceil((ymax-ymin)/cellsize);
		laguerre_volume_newton(simpar, xmin, ymin, xmax, ymax, boxsize, cellsize,
			mx0, my0, paddingAllTreeParticles, find2DNeighboringBP, find2DCellBP,
			mkLinkedList2D, Dtime, lag_vt_ready && (lag_calls % 2 == 0));
	} else for(i=0;i<nbp;i++){
		if(GAS_Kappa(simpar) <0) {
			bp[i].w2 = -GAS_Kappa(simpar);
		}
		else if (GAS_Kappa(simpar) >0 && !gfs_w2lin_active){
			bp[i].w2hydro = getw2forHydroParticle(simpar,(treevorork4particletype*)(bp+i),Dtime);
			applyW2Controls(simpar, (treevorork4particletype*)(bp+i), Dtime);
		}
	}
	int mx, my;
	BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
	BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	CellType *cells = (VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*mx*my));
	mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax, paddingAllTreeParticles);

#ifdef TESS_BENCHMARK
	/* ---- Tess Benchmark: Incremental vs Sutherland-Hodgman on CPU ---- */
	{
		static int _bench_done = 0;
		if (!_bench_done && MYID(simpar) == 0) {
			_bench_done = 1;
			postype boxsize_b = VORO_CELLBOX(simpar);
			double vol_inc_sum = 0, vol_sh_sum = 0;
			int nf_inc_sum = 0, nf_sh_sum = 0;
			int nbp_b = nbp;
			int my_b = my, mx_b = mx;

			/* --- Pass 1: Incremental (Voro2D_FindVC) --- */
			double t0_inc = MPI_Wtime();
			for (int iy_b = 0; iy_b < my_b; iy_b++) {
				int mp_b = 1000;
				Voro2D_Corner *vc = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp_b);
				for (int ix_b = 0; ix_b < mx_b; ix_b++) {
					int np_b;
					treevorork4particletype *p_b = find2DCellBP(simpar, ix_b, iy_b, &np_b);
					int nn_b;
					Voro2D_point *neigh_b = find2DNeighboringBP(simpar, ix_b, iy_b, &nn_b);
					Voro2D_point *nw_b = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nn_b);
					for (int ii = 0; ii < np_b; ii++) {
						Voro2D_point ctr;
						ctr.x = p_b[ii].x; ctr.y = p_b[ii].y;
						ctr.indx = PINDX(p_b+ii);
						ctr.csound = p_b[ii].csound;
						ctr.w2 = p_b[ii].w2;
						int ip_b = Voro2D_FindVC(&ctr, neigh_b, nw_b, nn_b, vc, mp_b, boxsize_b);
						/* Count faces and accumulate volume */
						Voro2D_Corner *et = vc;
						int fc = 0;
						double vol_b = 0;
						do {
							if (et->status == Active) {
								Voro2D_Corner *et2 = et->upperlink;
								vol_b += et->x * et2->y - et->y * et2->x;
								fc++;
							}
							et = et->upperlink;
						} while (et != vc);
						vol_inc_sum += 0.5 * fabs(vol_b);
						nf_inc_sum += fc;
					}
					if (np_b > 0) free(p_b);
					if (nn_b > 0) { free(neigh_b); free(nw_b); }
				}
				free(vc);
			}
			double t1_inc = MPI_Wtime();

			/* --- Pass 2: Sutherland-Hodgman --- */
			double t0_sh = MPI_Wtime();
			for (int iy_b = 0; iy_b < my_b; iy_b++) {
				for (int ix_b = 0; ix_b < mx_b; ix_b++) {
					int np_b;
					treevorork4particletype *p_b = find2DCellBP(simpar, ix_b, iy_b, &np_b);
					int nn_b;
					Voro2D_point *neigh_b = find2DNeighboringBP(simpar, ix_b, iy_b, &nn_b);
					Voro2D_point *nw_b = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nn_b);
					for (int ii = 0; ii < np_b; ii++) {
						/* Copy & subtract center (relative coords, same as Voro2D_FindVC) */
						postype w2c = p_b[ii].w2;
						int j_b = 0;
						for (int k = 0; k < nn_b; k++) {
							if (neigh_b[k].indx != PINDX(p_b+ii)) {
								nw_b[j_b].x = neigh_b[k].x - p_b[ii].x;
								nw_b[j_b].y = neigh_b[k].y - p_b[ii].y;
								nw_b[j_b].dist2 = nw_b[j_b].x*nw_b[j_b].x + nw_b[j_b].y*nw_b[j_b].y;
								nw_b[j_b].w2 = neigh_b[k].w2;
								double dd = nw_b[j_b].dist2;
								double ww = (w2c - nw_b[j_b].w2)/dd;
								nw_b[j_b].drad2 = 0.25*dd*(1+ww)*(1+ww);
								j_b++;
							}
						}
						qsort(nw_b, j_b, sizeof(Voro2D_point), sort2Ddrad2);
						int nf_b;
						double vol_b = Voro2D_SutherlandHodgman_CPU(w2c, nw_b, j_b, boxsize_b, &nf_b);
						vol_sh_sum += vol_b;
						nf_sh_sum += nf_b;
					}
					if (np_b > 0) free(p_b);
					if (nn_b > 0) { free(neigh_b); free(nw_b); }
				}
			}
			double t1_sh = MPI_Wtime();

			double ms_inc = (t1_inc - t0_inc)*1e3;
			double ms_sh  = (t1_sh  - t0_sh )*1e3;
			fprintf(stderr, "[TESS BENCH] nbp=%d  Incremental=%.1f ms  S-H=%.1f ms  ratio=%.2fx\n",
				nbp_b, ms_inc, ms_sh, ms_inc/ms_sh);
			fprintf(stderr, "[TESS BENCH] faces: inc=%d sh=%d  vol: inc=%.6e sh=%.6e  vol_diff=%.2e\n",
				nf_inc_sum, nf_sh_sum, vol_inc_sum, vol_sh_sum, fabs(vol_inc_sum-vol_sh_sum));
		}
	}
#endif /* TESS_BENCHMARK */

#ifdef USE_CUDA
	size_t _fused_p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *_fused_real_base = (char*)VORORK4_TBP(simpar);
	char *_fused_pad_base  = (char*)VORORK4_TBPP(simpar);
	int _fused_npad = VORO_NPAD(simpar);
	int _fused_nthreads = 1;
#ifdef _OPENMP
	_fused_nthreads = omp_get_max_threads();
#endif

  if (GAS_GPU_ENABLED(simpar)) {
	/* === GPU Tessellation Path === */
	{
		GPUContext *gctx = gpu_get_context();
		int n_total = nbp + _fused_npad;

		/* Lazy init */
		if (!gctx->initialized) {
			int max_p = (int)(n_total * 1.2) + 1024;
			int max_f = max_p * 20;
			gpu_init(gctx, max_p, max_f, MYID(simpar));
		}

		/* Realloc if needed */
		if (n_total > gctx->max_particles) {
			gpu_free(gctx);
			int max_p = (int)(n_total * 1.2) + 1024;
			int max_f = max_p * 20;
			gpu_init(gctx, max_p, max_f, MYID(simpar));
		}

		double _tess_t0 = MPI_Wtime();

		/* Build CellCSR from linked list */
		gpu_build_cell_csr(gctx, (void*)simpar, mx, my, _fused_p_size);

		/* Fill SoA (particles + padding) */
		{
			size_t p_size = _fused_p_size;
			char *real_base = _fused_real_base;
			char *pad_base  = _fused_pad_base;
			ParticleSoA *h = &gctx->h_parts;
			h->n_total = n_total;
			int ii;
			for (ii = 0; ii < n_total; ii++) {
				char *raw;
				if (ii < nbp)
					raw = real_base + ii * p_size;
				else
					raw = pad_base + (ii - nbp) * p_size;
				treevorork4particletype *bpi = (treevorork4particletype *)raw;
				h->x[ii]    = bpi->x;    h->y[ii]    = bpi->y;
				h->vx[ii]   = bpi->vx;   h->vy[ii]   = bpi->vy;
				h->mass[ii] = bpi->mass;
				h->indx[ii] = (long long)PINDX(bpi);
				h->den[ii]      = bpi->den;
				h->pressure[ii] = bpi->pressure;
				h->csound[ii]   = bpi->csound;
				h->volume[ii]   = bpi->volume;
				h->ie[ii]       = bpi->ie;
				h->w2[ii]       = bpi->w2;
				h->w2old[ii]    = bpi->w2old;
				if (use_stress) {
					treevorostressrk4particletype *sbpi =
						(treevorostressrk4particletype *)raw;
					h->alpha_cd[ii] = sbpi->stress.alpha_cd;
					h->K[ii]        = sbpi->stress.K;
				}
			}
		}

		/* Upload particles + cells to GPU */
		gctx->d_parts.n_total = n_total;
		gpu_upload_particles(gctx, n_total, use_stress);
		gpu_upload_cells(gctx, mx * my, gctx->h_cells.n_entries);

		postype OoA = VoroAccuracyOrder(simpar);

		/* Run GPU tessellation: tessellate + gradients + P/cs/stress + face CSR */
		gpu_tessellate_and_build_faces(gctx, nbp, mx, my,
			cellsize, xmin, ymin, boxsize, OoA, av_mode, Gamma,
			nu_phys_val, GAS_PRANDTL(simpar), GAS_CDAMAX(simpar),
			GAS_GRADIENT_METHOD(simpar),
			GAS_ENTROPY_MODE(simpar), GAS_K_FLOOR(simpar));

		/* Download results back to host */
		gpu_download_tess_results(gctx, nbp, use_stress);

		/* Write back to AoS particle structs */
		gpu_writeback_tess_results((void*)simpar, gctx, nbp, use_stress);

		double _tess_t1 = MPI_Wtime();
		if (MYID(simpar) == 0) {
			fprintf(stderr, "[GPU tess] nbp=%d npad=%d faces=%d time=%.3f ms\n",
				nbp, _fused_npad, gctx->d_faces.n_faces_total,
				(_tess_t1 - _tess_t0) * 1000.0);
		}
	}

#ifdef GPU_TESS_COMPARE
	/* === GPU FACE COMPARE: download GPU FaceCSR, let CPU run, compare === */
	{
		static int _face_cmp_call = 0;
		_face_cmp_call++;
		GPUContext *_fcctx = gpu_get_context();
		int _fc_nfaces = _fcctx->d_faces.n_faces_total;

		/* Allocate host buffers for GPU face CSR */
		int *_gf_offset = NULL;
		double *_gf_c1x=NULL, *_gf_c1y=NULL, *_gf_c2x=NULL, *_gf_c2y=NULL;
		int *_gf_neigh=NULL, *_gf_kp=NULL, *_gf_km=NULL, *_gf_ghost=NULL;

		if (_face_cmp_call <= 2 && MYID(simpar) == 0) {
			_gf_offset = (int*)malloc((nbp+1)*sizeof(int));
			_gf_c1x = (double*)malloc(_fc_nfaces*sizeof(double));
			_gf_c1y = (double*)malloc(_fc_nfaces*sizeof(double));
			_gf_c2x = (double*)malloc(_fc_nfaces*sizeof(double));
			_gf_c2y = (double*)malloc(_fc_nfaces*sizeof(double));
			_gf_neigh = (int*)malloc(_fc_nfaces*sizeof(int));
			_gf_kp = (int*)malloc(_fc_nfaces*sizeof(int));
			_gf_km = (int*)malloc(_fc_nfaces*sizeof(int));
			_gf_ghost = (int*)malloc(_fc_nfaces*sizeof(int));
			gpu_download_face_csr(_fcctx, nbp, _fc_nfaces,
				_gf_offset, _gf_c1x, _gf_c1y, _gf_c2x, _gf_c2y,
				_gf_neigh, _gf_kp, _gf_km, _gf_ghost);
			fprintf(stderr, "[FACE_CMP] call=%d gpu_nfaces=%d\n", _face_cmp_call, _fc_nfaces);
		}
		/* Reset tess_faces_on_device so CPU fused path and force use CPU faces */
		_fcctx->tess_faces_on_device = 0;
#else
	/* GPU did tessellation + P/cs/stress + face CSR. Skip CPU path. */
	goto _gpu_tess_skip_cpu;
#endif

	/* Fused face extraction init (for CPU fallback path) */
	gpu_fused_begin(nbp, _fused_npad, _fused_p_size, _fused_nthreads,
	                _fused_real_base, _fused_pad_base);
  } /* end if (GAS_GPU_ENABLED) */
#endif

	int iy;
	/* Two passes (EUNHA_TWOPASS_GRAD, default 1, CPU path): pass 0 builds the
	 * cells and sets volume and density; then pressure and sound speed are
	 * updated and the linked list (neighbour copies) is rebuilt; pass 1
	 * computes gradients and the limiter from the updated states. With one
	 * pass the gradients and limiter bounds used the previous stage's
	 * pressures while the face loop reconstructs around the updated ones,
	 * which produced negative face pressures (1-100% of the cell value). */
	int two_pass = gfs_twopass_on() && !GAS_GPU_ENABLED(simpar);
	int gpass;
	for(gpass = two_pass ? 0 : 1; gpass < 2; gpass++){
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(iy=0;iy<my;iy++){
#ifdef USE_CUDA
		int _fused_tid = 0;
#ifdef _OPENMP
		_fused_tid = omp_get_thread_num();
#endif
#endif
		int mp=4000;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		int ix;
		for(ix=0;ix<mx;ix++){
			int np;
			treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
			int nneigh;
			Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
			Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			int i;
			for(i=0;i<np;i++){
				Voro2D_point center;
				center.x = p[i].x;
				center.y = p[i].y;
				center.indx = PINDX(p+i);
				center.csound = p[i].csound;
				center.w2 = p[i].w2;

				int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
				/* Same cell as the face loop: when the 3x3 stencil leaves a
				 * corner unassigned, widen the search (R = 2..4) exactly as
				 * getAccVoro2DBlend_impl does. Volume, gradients and the
				 * limiter then see the faces the force and fluxes use. */
				Voro2D_point *neighbors_e = NULL, *neighwork_e = NULL;
				int nneigh_e = 0, R_searched = 1;
				Voro2D_point *neighbors_save = neighbors;
				Voro2D_point *neighwork_save = neighwork;
				int nneigh_save = nneigh;
				if(Voro2D_HasUnassignedCorner(vorocorner, ip) || voro_sec_short(simpar, vorocorner, ip, (double)center.x, (double)center.y, ix, iy, 1)){
					int R;
					for(R=2; R<=voro_expand_rmax(simpar); R++){
						if(neighbors_e){ free(neighbors_e); free(neighwork_e); neighbors_e=NULL; neighwork_e=NULL; }
						neighbors_e = searchCellRk4Neighbors2D_R(simpar, ix, iy, R, &nneigh_e);
						neighwork_e = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh_e);
						ip = Voro2D_FindVC(&center, neighbors_e, neighwork_e, nneigh_e, vorocorner, mp, boxsize);
						R_searched = R;
						if(!Voro2D_HasUnassignedCorner(vorocorner, ip) && !voro_sec_short(simpar, vorocorner, ip, (double)center.x, (double)center.y, ix, iy, R)) break;
					}
					if(neighbors_e){
						neighbors = neighbors_e;
						neighwork = neighwork_e;
						nneigh = nneigh_e;
					}
				}
				if(voro_sec_short(simpar, vorocorner, ip, (double)center.x, (double)center.y, ix, iy, R_searched)){
#ifdef _OPENMP
#pragma omp atomic
#endif
					voro_sec_fail++;
				}
				{ double vrR = voro_disc_radius(simpar); if(vrR > 0) voro2d_restrict_disc(vorocorner, mp, vrR); }
				treevorork4particletype *ibp_rk4 = p[i].bp;
				treevorostressrk4particletype *ibp = (treevorostressrk4particletype*)ibp_rk4;
				if(gpass == 0 || !two_pass){
					get2dAreaAvgNeighorPressure(ibp_rk4, vorocorner, neighwork, (treevorork4particletype*)bp);
					ibp_rk4->den = ibp_rk4->mass/ibp_rk4->volume;
					if(gfs_shared_geom_on()) shf_record(vorocorner, neighwork, ibp_rk4);
					if(gfs_tcf_face_z() || gfs_periodic_complete()){
						double cdx_, cdy_; int nc_ = cen_poly(vorocorner, &cdx_, &cdy_);
						size_t id_ = (size_t)PINDX(ibp_rk4);
						if(id_ < cen_cap){ cen_dx[id_] = cdx_; cen_dy[id_] = cdy_; }
						if(nc_ > 0){
#ifdef _OPENMP
#pragma omp atomic
#endif
							pc_ncap++;
						}
#ifdef _OPENMP
#pragma omp atomic
#endif
						pc_vsum += (double)ibp_rk4->volume;
					}
					if(gfs_subzone_on()){
						size_t id_ = (size_t)PINDX(ibp_rk4); sub_poly S_; double g_ = 1;
						if(!(sub_build(vorocorner, neighwork, ibp_rk4, &S_) && sub_eval(&S_, 0, 0, (double)Gamma, &g_, NULL, NULL, NULL, NULL))) g_ = 1;
						if(id_ < sub_cap) sub_g[id_] = g_;
					}
				}
				else {
					/* pass 1: keep the pass-0 volume and density (other threads read
					 * them as neighbour values); refresh only the face-length-weighted
					 * neighbour pressure (kappa > 0 weights) from the updated copies. */
					postype aP = 0, L = 0;
					Voro2D_Corner *ta = vorocorner;
					do {
						if(ta->upperrelated >= 0){
							postype dxa = ta->x - ta->upperlink->x, dya = ta->y - ta->upperlink->y, ds = sqrt(dxa*dxa + dya*dya);
							L += ds; aP += ds*((treevorork4particletype*)neighwork[ta->upperrelated].bp)->pressure;
						}
						ta = ta->upperlink;
					} while(ta != vorocorner);
					ibp_rk4->avgNeighboringPressure = (L > 0) ? aP/L : ibp_rk4->pressure;
				}
				if(sedov_lagvol_on() && lag_vt_ready){
					size_t id = (size_t)PINDX(ibp_rk4);
					if(id < lag_vt_cap && lag_vt[id] > 0){
						ibp_rk4->volume = (postype)lag_vt[id];
						ibp_rk4->den = ibp_rk4->mass / ibp_rk4->volume;
					}
				}
				{
					/* A generator farther than 2 max_v |v - x_i| cannot cut the cell; the
					 * stencil searched so far covers at least R_used * cellwidth. */
					double vmax2 = 0; Voro2D_Corner *tv = vorocorner;
					do { double q = (double)tv->x*(double)tv->x + (double)tv->y*(double)tv->y; if(q > vmax2) vmax2 = q; tv = tv->upperlink; } while(tv != vorocorner);
					if(2*sqrt(vmax2) > R_searched*(double)cellsize){
#ifdef _OPENMP
#pragma omp atomic
#endif
						gfs_cellcheck_fail++;
					}
#ifdef _OPENMP
#pragma omp atomic
#endif
					gfs_cellcheck_n++;
				}
				if(gpass == 0 && gfs_topo_on() && gfs_rk_stage >= 1){
					long ids[16]; int nid = 0, ok16 = 1;
					Voro2D_Corner *tv = vorocorner;
					do {
						if(tv->upperrelated >= 0){
							if(nid < 16) ids[nid++] = (long)PINDX((treevorork4particletype*)(neighwork[tv->upperrelated].bp));
							else ok16 = 0;
						}
						tv = tv->upperlink;
					} while(tv != vorocorner);
					long me = (long)PINDX((treevorork4particletype*)ibp);
					if(ok16 && me >= 0){
						qsort(ids, nid, sizeof(long), gfs_cmp_long);
#ifdef _OPENMP
#pragma omp critical (topo_audit)
#endif
						{
							if(me >= gfs_topo_cap){
								long nc = (me + 1)*2;
								gfs_topo_ids = (long*)realloc(gfs_topo_ids, sizeof(long)*16*nc);
								gfs_topo_cnt = (unsigned char*)realloc(gfs_topo_cnt, nc);
								memset(gfs_topo_cnt + gfs_topo_cap, 255, nc - gfs_topo_cap);
								gfs_topo_cap = nc;
							}
							if(gfs_rk_stage == 1){
								memcpy(gfs_topo_ids + 16*me, ids, sizeof(long)*nid);
								gfs_topo_cnt[me] = (unsigned char)nid;
							} else if(gfs_topo_cnt[me] != 255){
								int diff = (gfs_topo_cnt[me] != nid) || memcmp(gfs_topo_ids + 16*me, ids, sizeof(long)*nid);
								gfs_topo_chg += diff; gfs_topo_n++;
							}
						}
					}
				}
				if(gpass == 0){
					if(neighbors_e){
						free(neighbors_e); free(neighwork_e);
						neighbors = neighbors_save;
						neighwork = neighwork_save;
						nneigh = nneigh_save;
					}
					continue;
				}

#ifdef USE_CUDA
				/* Fallback face extraction when use_stress==false (no gradient traversal).
				   When use_stress==true, extraction is piggybacked on the gradient loop below. */
				if(!use_stress){
					int _fi = gpu_ptr_to_soa_index(ibp_rk4, _fused_real_base,
					            _fused_pad_base, nbp, _fused_npad, _fused_p_size);
					if(_fi >= 0 && _fi < nbp){
						Voro2D_Corner *_et = vorocorner;
						do {
							if(_et->upperrelated >= 0){
								Voro2D_Corner *_et2 = _et->upperlink;
								void *_jptr = neighwork[_et->upperrelated].bp;
								int _fj = gpu_ptr_to_soa_index(_jptr, _fused_real_base,
								            _fused_pad_base, nbp, _fused_npad, _fused_p_size);
								int _ghost = (PINDX((treevorork4particletype*)_jptr) == MAX_INDEX);
								int _fkp = -1, _fkm = -1;
								if(_et2->upperrelated >= 0)
									_fkp = gpu_ptr_to_soa_index(neighwork[_et2->upperrelated].bp,
									         _fused_real_base, _fused_pad_base, nbp, _fused_npad, _fused_p_size);
								if(_et->lowerrelated >= 0)
									_fkm = gpu_ptr_to_soa_index(neighwork[_et->lowerrelated].bp,
									         _fused_real_base, _fused_pad_base, nbp, _fused_npad, _fused_p_size);
								gpu_fused_add_face(_fused_tid, _fi,
									_et->x, _et->y, _et2->x, _et2->y,
									_fj, _fkp, _fkm, _ghost);
							}
							_et = _et->upperlink;
						} while(_et != vorocorner);
					}
				}
#endif

				/* Compute cell-averaged gradients via Green-Gauss:
				 *   ⟨∇Q⟩_i = (1/V_i) Σ_faces Q_face·dS
				 * Face-averaged values use M(n,m) stencil.
				 * Computes: ∇⊗v (velocity gradient), ∇P (MUSCL), ∇ρ (Arepo MUSCL).
				 *
				 * For av_mode==5, a second pass computes per-cell Arepo-style
				 * (Springel 2010, Barth-Jespersen) TVD slope limiter alpha_i
				 * for each of {rho, vx, vy, P} over ALL faces of cell i, and
				 * rescales the stored gradients in-place so that the MUSCL
				 * reconstruction at any face cannot create a new extremum
				 * relative to the cell-averaged values of i's neighbors.
				 * Other av_modes use unlimited gradients (for NS stress). */
				if(use_stress){
					int use_lsf = (GAS_GRADIENT_METHOD(simpar) == 1);
					double q2H[4][3], q2g[4][2]; int q2smooth[4] = {0, 0, 0, 0}, q2ok = 0;   /* EUNHA_LSF_ORDER=2: vx, vy, P, rho */
					/* Green-Gauss accumulators */
					postype sum_gUxx=0, sum_gUxy=0, sum_gUyx=0, sum_gUyy=0;
					postype sum_dPdx=0, sum_dPdy=0;
					postype sum_dRhodx=0, sum_dRhody=0;
					/* LSF (Pakmor 2016) accumulators: 2x2 normal-equation matrix
					 * M = Σ w (Δr Δr^T), and per-quantity RHS b_Q = Σ w Δr ΔQ.
					 * Weight w = A_face / |Δr|^2 (face-area-weighted). */
					postype Mxx=0, Mxy=0, Myy=0;
					postype bxVx=0, byVx=0, bxVy=0, byVy=0;
					postype bxP=0,  byP=0,  bxR=0,  byR=0;
					/* Hyperviscosity Pass 1: lap_v_i = (1/V_i) Σ A_ij (v_j-v_i)/d_ij
					 * Real-neighbor faces only; ghost mirrors excluded. */
					postype sum_lap_vx=0, sum_lap_vy=0;
					postype OoA = VoroAccuracyOrder(simpar);
					Voro2D_Corner *tmp = vorocorner;
					do {
						if(tmp->upperrelated >= 0){
							treevorostressrk4particletype *jbp = (treevorostressrk4particletype*)(neighwork[tmp->upperrelated].bp);
							int jbp_grad_ghost = (PINDX((treevorork4particletype*)jbp) == MAX_INDEX);
							Voro2D_Corner *tmp2 = tmp->upperlink;
							/* Face area vector dS (outward normal * area) */
							postype dSx = tmp2->y - tmp->y;
							postype dSy = -(tmp2->x - tmp->x);

							/* Laguerre-aware face interpolation weight */
							postype dx_ij = jbp->x - ibp->x;
							postype dy_ij = jbp->y - ibp->y;
							postype d2_ij = dx_ij*dx_ij + dy_ij*dy_ij;

							/* Hyperviscosity Pass 1: accumulate (v_j-v_i) weighted
							 * by A_face/d_ij; ghost faces excluded. */
							if(!jbp_grad_ghost && d2_ij > 0){
								postype d_ij = sqrt(d2_ij);
								postype A_face_lap = sqrt(dSx*dSx + dSy*dSy);
								postype invd_lap = A_face_lap / d_ij;
								sum_lap_vx += invd_lap * (jbp->vx - ibp->vx);
								sum_lap_vy += invd_lap * (jbp->vy - ibp->vy);
							}

							if(use_lsf){
								/* Pakmor 2016 LSF: weight = A_face / |Δr|^2 */
								postype A_face = sqrt(dSx*dSx + dSy*dSy);
								postype w = (d2_ij > 0) ? A_face / d2_ij : 0.0;
								postype dvx = jbp->vx - ibp->vx;
								postype dvy = jbp->vy - ibp->vy;
								postype dP  = jbp->pressure - ibp->pressure;
								postype dR  = jbp->den - ibp->den;
								Mxx += w*dx_ij*dx_ij; Mxy += w*dx_ij*dy_ij; Myy += w*dy_ij*dy_ij;
								bxVx += w*dx_ij*dvx; byVx += w*dy_ij*dvx;
								bxVy += w*dx_ij*dvy; byVy += w*dy_ij*dvy;
								bxP  += w*dx_ij*dP;  byP  += w*dy_ij*dP;
								bxR  += w*dx_ij*dR;  byR  += w*dy_ij*dR;
							} else {
								postype wfrac = (d2_ij > 0) ?
									0.5 + 0.5*(ibp_rk4->w2 - ((treevorork4particletype*)jbp)->w2)/d2_ij : 0.5;

								postype vx_face, vy_face, P_face, Rho_face;
								if(jbp_grad_ghost){
									/* Wall ghost: simple weighted average */
									vx_face = ibp->vx + wfrac*(jbp->vx - ibp->vx);
									vy_face = ibp->vy + wfrac*(jbp->vy - ibp->vy);
									P_face  = ibp->pressure + wfrac*(jbp->pressure - ibp->pressure);
									Rho_face= ibp->den + wfrac*(jbp->den - ibp->den);
								} else if(OoA > 0 && (tmp->upperlink)->upperrelated >= 0 && tmp->lowerrelated >= 0){
									/* M(n,m) 4-point stencil: face-averaged value */
									treevorostressrk4particletype *kp = (treevorostressrk4particletype*)(neighwork[(tmp->upperlink)->upperrelated].bp);
									treevorostressrk4particletype *km = (treevorostressrk4particletype*)(neighwork[tmp->lowerrelated].bp);
									postype OoA3 = OoA/3.;
									vx_face = ibp->vx + wfrac*(jbp->vx - ibp->vx)
									        + OoA3*(kp->vx + km->vx - ibp->vx - jbp->vx);
									vy_face = ibp->vy + wfrac*(jbp->vy - ibp->vy)
									        + OoA3*(kp->vy + km->vy - ibp->vy - jbp->vy);
									P_face  = ibp->pressure + wfrac*(jbp->pressure - ibp->pressure)
									        + OoA3*(kp->pressure + km->pressure - ibp->pressure - jbp->pressure);
									Rho_face= ibp->den + wfrac*(jbp->den - ibp->den)
									        + OoA3*(kp->den + km->den - ibp->den - jbp->den);
								} else {
									vx_face = ibp->vx + wfrac*(jbp->vx - ibp->vx);
									vy_face = ibp->vy + wfrac*(jbp->vy - ibp->vy);
									P_face  = ibp->pressure + wfrac*(jbp->pressure - ibp->pressure);
									Rho_face= ibp->den + wfrac*(jbp->den - ibp->den);
								}

								/* Accumulate Green-Gauss: ∮ Q_face · dS */
								sum_gUxx += vx_face * dSx;
								sum_gUxy += vx_face * dSy;
								sum_gUyx += vy_face * dSx;
								sum_gUyy += vy_face * dSy;
								sum_dPdx += P_face * dSx;
								sum_dPdy += P_face * dSy;
								sum_dRhodx += Rho_face * dSx;
								sum_dRhody += Rho_face * dSy;
							}

#ifdef USE_CUDA
							/* Fused CSR face extraction (piggyback on gradient traversal).
							 * Guard: only when GPU is enabled at runtime — the fused buffers
							 * are only initialized via gpu_fused_begin inside the
							 * GAS_GPU_ENABLED branch above. */
							if(GAS_GPU_ENABLED(simpar)){
								int _fi = gpu_ptr_to_soa_index(ibp_rk4, _fused_real_base,
								            _fused_pad_base, nbp, _fused_npad, _fused_p_size);
								if(_fi >= 0 && _fi < nbp){
									int _fj = gpu_ptr_to_soa_index((void*)jbp, _fused_real_base,
									            _fused_pad_base, nbp, _fused_npad, _fused_p_size);
									int _fkp = -1, _fkm = -1;
									if(tmp2->upperrelated >= 0)
										_fkp = gpu_ptr_to_soa_index(neighwork[tmp2->upperrelated].bp,
										         _fused_real_base, _fused_pad_base, nbp, _fused_npad, _fused_p_size);
									if(tmp->lowerrelated >= 0)
										_fkm = gpu_ptr_to_soa_index(neighwork[tmp->lowerrelated].bp,
										         _fused_real_base, _fused_pad_base, nbp, _fused_npad, _fused_p_size);
									gpu_fused_add_face(_fused_tid, _fi,
										tmp->x, tmp->y, tmp2->x, tmp2->y,
										_fj, _fkp, _fkm, jbp_grad_ghost);
								}
							}
#endif
						}
						tmp = tmp->upperlink;
					} while(tmp != vorocorner);

					/* vacuum arcs: boundary data for rho and P only (separate normal matrix) */
					double Mbxx = 0, Mbxy = 0, Mbyy = 0, cbx = 0, cby = 0;
					if(use_lsf && gfs_vacbc_on()){
						double pa = 0;
						Voro2D_Corner *ta = vorocorner;
						do { pa += (double)ta->x*(double)ta->upperlink->y - (double)ta->upperlink->x*(double)ta->y; ta = ta->upperlink; } while(ta != vorocorner);
						double orient = (pa >= 0) ? 1.0 : -1.0;
						ta = vorocorner;
						do {
							if(ta->upperrelated == VORO_ARC_EDGE){
								Voro2D_Corner *tb = ta->upperlink;
								double R = sqrt((double)ta->x*(double)ta->x + (double)ta->y*(double)ta->y);
								double th0 = atan2((double)ta->y, (double)ta->x);
								double sw = orient*(atan2((double)tb->y, (double)tb->x) - th0);
								while(sw < 0) sw += 2*M_PI;
								while(sw >= 2*M_PI) sw -= 2*M_PI;
								double a0 = (orient > 0) ? th0 : th0 - sw, a1 = a0 + sw;   /* increasing angle */
								/* w = 1/R^2, ds = R dtheta */
								Mbxx += R*(0.5*sw + 0.25*(sin(2*a1) - sin(2*a0)));
								Mbyy += R*(0.5*sw - 0.25*(sin(2*a1) - sin(2*a0)));
								Mbxy += R*0.5*(sin(a1)*sin(a1) - sin(a0)*sin(a0));
								cbx += sin(a1) - sin(a0);          /* int cos dtheta */
								cby += -(cos(a1) - cos(a0));       /* int sin dtheta */
							}
							ta = ta->upperlink;
						} while(ta != vorocorner);
					}
					if(gfs_tcf_qbp() && gfs_qbp_edge_inner()){
						Voro2D_Corner *ta = vorocorner; int has_ = 0;
						do { if(ta->upperrelated == VORO_ARC_EDGE){ has_ = 1; break; } ta = ta->upperlink; } while(ta != vorocorner);
						size_t id_ = (size_t)PINDX(ibp_rk4);
						if(id_ < qbp_cap) qbp_arc[id_] = (char)has_;
					}
					if(use_lsf){
						postype det = Mxx*Myy - Mxy*Mxy;
						if(fabs(det) > 1e-30){
							postype invDet = 1.0/det;
							ibp->stress.gUxx   = ( Myy*bxVx - Mxy*byVx) * invDet;
							ibp->stress.gUxy   = (-Mxy*bxVx + Mxx*byVx) * invDet;
							ibp->stress.gUyx   = ( Myy*bxVy - Mxy*byVy) * invDet;
							ibp->stress.gUyy   = (-Mxy*bxVy + Mxx*byVy) * invDet;
							ibp->stress.dPdx   = ( Myy*bxP  - Mxy*byP ) * invDet;
							ibp->stress.dPdy   = (-Mxy*bxP  + Mxx*byP ) * invDet;
							ibp->stress.dRhodx = ( Myy*bxR  - Mxy*byR ) * invDet;
							ibp->stress.dRhody = (-Mxy*bxR  + Mxx*byR ) * invDet;
							if(gfs_tcf_qbp() && (double)det > 64*2.2e-16*(double)Mxx*(double)Myy){
								size_t id_ = (size_t)PINDX(ibp_rk4);
								if(id_ < qbp_cap){ qbp_g[4*id_] = (double)ibp->stress.gUxx; qbp_g[4*id_+1] = (double)ibp->stress.gUxy; qbp_g[4*id_+2] = (double)ibp->stress.gUyx; qbp_g[4*id_+3] = (double)ibp->stress.gUyy; qbp_ok[id_] = 1; }
							}
						}
						if(Mbxx + Mbyy > 0){
							/* rho, P with the vacuum arc data: b += int r (0 - Q_i)/R^2 ds = -Q_i (cbx, cby) */
							double Axx = Mxx + Mbxx, Axy = Mxy + Mbxy, Ayy = Myy + Mbyy, dt2 = Axx*Ayy - Axy*Axy;
							if(fabs(dt2) > 1e-30){
								double bPx = bxP - (double)ibp->pressure*cbx, bPy = byP - (double)ibp->pressure*cby;
								double bRx = bxR - (double)ibp->den*cbx,      bRy = byR - (double)ibp->den*cby;
								ibp->stress.dPdx   = (postype)(( Ayy*bPx - Axy*bPy)/dt2);
								ibp->stress.dPdy   = (postype)((-Axy*bPx + Axx*bPy)/dt2);
								ibp->stress.dRhodx = (postype)(( Ayy*bRx - Axy*bRy)/dt2);
								ibp->stress.dRhody = (postype)((-Axy*bRx + Axx*bRy)/dt2);
							}
						}
						if(!(fabs(det) > 1e-30)){
							ibp->stress.gUxx = ibp->stress.gUxy = 0;
							ibp->stress.gUyx = ibp->stress.gUyy = 0;
							if(!(Mbxx + Mbyy > 0)){
								ibp->stress.dPdx = ibp->stress.dPdy = 0;
								ibp->stress.dRhodx = ibp->stress.dRhody = 0;
							}
						}
						ibp->stress.divv = ibp->stress.gUxx + ibp->stress.gUyy;
						if(gfs_lsf_order() == 2 && gfs_q2_rank()){
							/* GFS_Q2_RANK=1 (Astra audit 13 (a)): weighted LSQ on the design matrix with the
							 * rotation-orthonormal basis (X, Y, X^2/2, XY/sqrt2, Y^2/2), X = dx/hs, w = hs^2/d^2;
							 * modified Gram-Schmidt QR (twice); rank 5 decided at machine precision
							 * R_kk > 64 eps sqrt(n) |col_k|. Start from the points within 2.5 hs and add
							 * whole distance shells of the searched candidates until rank 5; else no fit. */
							double hs = sqrt((double)ibp_rk4->volume);
							int nwork = 0, kk;
							for(kk = 0; kk < nneigh; kk++) if(neighbors[kk].indx != center.indx) nwork++;
							int nc = 0;
							for(kk = 0; kk < nwork && nc < 4096; kk++){
								if(!neighwork[kk].bp || neighwork[kk].bp == (void*)ibp) continue;
								double ex = (double)neighwork[kk].x, ey = (double)neighwork[kk].y, d2 = ex*ex + ey*ey;
								if(!(d2 > 0)) continue;
								treevorostressrk4particletype *jq = (treevorostressrk4particletype*)(neighwork[kk].bp);
								double X = ex/hs, Y = ey/hs, sw = hs/sqrt(d2);
								q2r_cd2[nc] = d2;
								q2r_cph[nc][0] = sw*X; q2r_cph[nc][1] = sw*Y; q2r_cph[nc][2] = sw*0.5*X*X; q2r_cph[nc][3] = sw*X*Y/sqrt(2.0); q2r_cph[nc][4] = sw*0.5*Y*Y;
								q2r_cph[nc][5] = sw*X*X*X/6; q2r_cph[nc][6] = sw*X*X*Y/(2*sqrt(3.0)); q2r_cph[nc][7] = sw*X*Y*Y/(2*sqrt(3.0)); q2r_cph[nc][8] = sw*Y*Y*Y/6;
								q2r_cdq[nc][0] = sw*(double)(jq->vx - ibp->vx); q2r_cdq[nc][1] = sw*(double)(jq->vy - ibp->vy);
								q2r_cdq[nc][2] = sw*(double)(jq->pressure - ibp->pressure); q2r_cdq[nc][3] = sw*(double)(jq->den - ibp->den);
								q2r_cord[nc] = nc; nc++;
							}
							{ int a_, b_; for(a_ = 1; a_ < nc; a_++){ int t_ = q2r_cord[a_]; double v_ = q2r_cd2[t_]; b_ = a_ - 1; while(b_ >= 0 && q2r_cd2[q2r_cord[b_]] > v_){ q2r_cord[b_+1] = q2r_cord[b_]; b_--; } q2r_cord[b_+1] = t_; } }
							int nuse = 0; while(nuse < nc && q2r_cd2[q2r_cord[nuse]] <= 6.25*hs*hs) nuse++;
							int ok = 0, rank = 0, ncol = gfs_q3() ? 9 : 5, nuse0 = nuse;
							q3_retry:;
							double Rm[9][9], xs[4][9];
							for(;;){
								if(nuse >= ncol){
									int col, z, kk2, pass; rank = 0;
									for(col = 0; col < ncol; col++){
										double nrm0 = 0;
										for(z = 0; z < nuse; z++){ q2r_Qc[z][col] = q2r_cph[q2r_cord[z]][col]; nrm0 += q2r_Qc[z][col]*q2r_Qc[z][col]; }
										nrm0 = sqrt(nrm0);
										for(kk2 = 0; kk2 < ncol; kk2++) Rm[kk2][col] = 0;
										for(pass = 0; pass < 2; pass++) for(kk2 = 0; kk2 < col; kk2++){
											double dot = 0; for(z = 0; z < nuse; z++) dot += q2r_Qc[z][kk2]*q2r_Qc[z][col];
											for(z = 0; z < nuse; z++) q2r_Qc[z][col] -= dot*q2r_Qc[z][kk2];
											Rm[kk2][col] += dot;
										}
										double nn = 0; for(z = 0; z < nuse; z++) nn += q2r_Qc[z][col]*q2r_Qc[z][col]; nn = sqrt(nn);
										if(!(nn > 64*2.2e-16*sqrt((double)nuse)*nrm0)) break;
										for(z = 0; z < nuse; z++) q2r_Qc[z][col] /= nn;
										Rm[col][col] = nn; rank++;
									}
									if(rank == ncol){ ok = 1; break; }
								}
								if(nuse >= nc) break;
								double dnext = q2r_cd2[q2r_cord[nuse]];
								while(nuse < nc && q2r_cd2[q2r_cord[nuse]] <= dnext*(1 + 1e-12)) nuse++;
							}
							if(!ok && ncol == 9){
#ifdef _OPENMP
#pragma omp atomic
#endif
								q3_fb++;
								ncol = 5; nuse = nuse0; goto q3_retry;
							}
							if(ok && ncol == 9){
#ifdef _OPENMP
#pragma omp atomic
#endif
								q3_n++;
							}
#ifdef _OPENMP
#pragma omp atomic
#endif
							gfs_q2_n++;
							if(!ok){
#ifdef _OPENMP
#pragma omp atomic
#endif
								gfs_q2_fail++;
								if(gfs_tcf_qbp()){
#ifdef _OPENMP
#pragma omp atomic
#endif
									qbp_q2fail_piv++;
								}
							} else {
								int q_, a_, c_, z;
								if(nuse > (int)(0) && q2r_cd2[q2r_cord[nuse-1]] > 6.25*hs*hs && gfs_tcf_qbp()){
#ifdef _OPENMP
#pragma omp atomic
#endif
									qbp_q2fail_npt++;   /* here: fits that needed a wider shell */
								}
								for(q_ = 0; q_ < 4; q_++){
									double qb[9], S2_ = 0;
									for(a_ = 0; a_ < ncol; a_++){ double v_ = 0; for(z = 0; z < nuse; z++) v_ += q2r_Qc[z][a_]*q2r_cdq[q2r_cord[z]][q_]; qb[a_] = v_; }
									for(z = 0; z < nuse; z++) S2_ += q2r_cdq[q2r_cord[z]][q_]*q2r_cdq[q2r_cord[z]][q_];
									for(a_ = ncol - 1; a_ >= 0; a_--){ double v_ = qb[a_]; for(c_ = a_+1; c_ < ncol; c_++) v_ -= Rm[a_][c_]*xs[q_][c_]; xs[q_][a_] = v_/Rm[a_][a_]; }
									double qq = 0; for(a_ = 0; a_ < ncol; a_++) qq += qb[a_]*qb[a_];
									double res2 = S2_ - qq; if(res2 < 0) res2 = 0;
									q2smooth[q_] = (res2 <= 0.01*S2_);
									q2H[q_][0] = xs[q_][2]/(hs*hs); q2H[q_][1] = xs[q_][3]/(sqrt(2.0)*hs*hs); q2H[q_][2] = xs[q_][4]/(hs*hs);
									q2g[q_][0] = xs[q_][0]/hs; q2g[q_][1] = xs[q_][1]/hs;
								}
								ibp->stress.gUxx = (postype)q2g[0][0]; ibp->stress.gUxy = (postype)q2g[0][1];
								ibp->stress.gUyx = (postype)q2g[1][0]; ibp->stress.gUyy = (postype)q2g[1][1];
								if(gfs_tcf_qbp()){ size_t id_ = (size_t)PINDX(ibp_rk4); if(id_ < qbp_cap){ qbp_g[4*id_] = q2g[0][0]; qbp_g[4*id_+1] = q2g[0][1]; qbp_g[4*id_+2] = q2g[1][0]; qbp_g[4*id_+3] = q2g[1][1]; qbp_ok[id_] = 2; int h_; for(h_ = 0; h_ < 3; h_++){ qbp_h[6*id_+h_] = q2H[0][h_]; qbp_h[6*id_+3+h_] = q2H[1][h_]; }
									if(ncol == 9){ double h3 = hs*hs*hs; for(h_ = 0; h_ < 2; h_++){ qbp_t[8*id_+4*h_] = xs[h_][5]/h3; qbp_t[8*id_+4*h_+1] = xs[h_][6]/(sqrt(3.0)*h3); qbp_t[8*id_+4*h_+2] = xs[h_][7]/(sqrt(3.0)*h3); qbp_t[8*id_+4*h_+3] = xs[h_][8]/h3; } } } }
								ibp->stress.dPdx = (postype)q2g[2][0]; ibp->stress.dPdy = (postype)q2g[2][1];
								ibp->stress.dRhodx = (postype)q2g[3][0]; ibp->stress.dRhody = (postype)q2g[3][1];
								ibp->stress.divv = ibp->stress.gUxx + ibp->stress.gUyy;
								q2ok = 1;
							}
						} else
						if(gfs_lsf_order() == 2){
							/* columns scaled by hs: a = (gx hs, gy hs, hxx hs^2, hxy hs^2, hyy hs^2),
							 * dq = a0 X + a1 Y + a2 X^2/2 + a3 X Y + a4 Y^2/2, X = dx/hs */
							double hs = sqrt((double)ibp_rk4->volume);
							double rs = 2.5*hs, rg = (double)R_searched*(double)cellsize;
							if(rs > rg) rs = rg;
							double N[5][5] = {{0}}, Bq[4][5] = {{0}}, S2[4] = {0, 0, 0, 0}, W = 0;
							/* Voro2D_FindVC fills neighwork with the candidates other than the centre */
							int npt = 0, kk, nwork = 0;
							for(kk = 0; kk < nneigh; kk++) if(neighbors[kk].indx != center.indx) nwork++;
							for(kk = 0; kk < nwork; kk++){
								if(!neighwork[kk].bp || neighwork[kk].bp == (void*)ibp) continue;
								double ex = (double)neighwork[kk].x, ey = (double)neighwork[kk].y, d2 = ex*ex + ey*ey;
								if(d2 <= 0 || d2 > rs*rs) continue;
								treevorostressrk4particletype *jq = (treevorostressrk4particletype*)(neighwork[kk].bp);
								double X = ex/hs, Y = ey/hs, ph[5] = {X, Y, 0.5*X*X, X*Y, 0.5*Y*Y}, w = hs*hs/d2;
								double dq[4] = {(double)(jq->vx - ibp->vx), (double)(jq->vy - ibp->vy), (double)(jq->pressure - ibp->pressure), (double)(jq->den - ibp->den)};
								int a_, b_, q_;
								for(a_ = 0; a_ < 5; a_++){
									for(b_ = a_; b_ < 5; b_++) N[a_][b_] += w*ph[a_]*ph[b_];
									for(q_ = 0; q_ < 4; q_++) Bq[q_][a_] += w*ph[a_]*dq[q_];
								}
								for(q_ = 0; q_ < 4; q_++) S2[q_] += w*dq[q_]*dq[q_];
								W += w; npt++;
							}
							int ok = (npt >= 8);
							if(!ok && gfs_tcf_qbp()){
#ifdef _OPENMP
#pragma omp atomic
#endif
								qbp_q2fail_npt++;
							}
							int ok_npt = ok;
							double L[5][5] = {{0}};
							if(ok){
								int a_, b_, c_;
								double tr = 0;
								for(a_ = 0; a_ < 5; a_++) tr += N[a_][a_];
								for(a_ = 0; a_ < 5 && ok; a_++){
									for(b_ = 0; b_ <= a_; b_++){
										double v = N[b_][a_];
										for(c_ = 0; c_ < b_; c_++) v -= L[a_][c_]*L[b_][c_];
										if(a_ == b_){ if(v <= 1e-8*tr){ ok = 0; break; } L[a_][a_] = sqrt(v); }
										else L[a_][b_] = v/L[b_][b_];
									}
								}
							}
#ifdef _OPENMP
#pragma omp atomic
#endif
							gfs_q2_n++;
							if(!ok){
#ifdef _OPENMP
#pragma omp atomic
#endif
								gfs_q2_fail++;
								if(ok_npt && gfs_tcf_qbp()){
#ifdef _OPENMP
#pragma omp atomic
#endif
									qbp_q2fail_piv++;
								}
							} else {
								int q_, a_, c_;
								for(q_ = 0; q_ < 4; q_++){
									double z[5], x_[5];
									for(a_ = 0; a_ < 5; a_++){ double v = Bq[q_][a_]; for(c_ = 0; c_ < a_; c_++) v -= L[a_][c_]*z[c_]; z[a_] = v/L[a_][a_]; }
									for(a_ = 4; a_ >= 0; a_--){ double v = z[a_]; for(c_ = a_+1; c_ < 5; c_++) v -= L[c_][a_]*x_[c_]; x_[a_] = v/L[a_][a_]; }
									/* weighted residual: S2 - 2 x.B + x.N.x = S2 - x.B at the optimum */
									double xb = 0; for(a_ = 0; a_ < 5; a_++) xb += x_[a_]*Bq[q_][a_];
									double res2 = S2[q_] - xb; if(res2 < 0) res2 = 0;
									q2smooth[q_] = (res2 <= 0.01*S2[q_]);
									q2H[q_][0] = x_[2]/(hs*hs); q2H[q_][1] = x_[3]/(hs*hs); q2H[q_][2] = x_[4]/(hs*hs);
									q2g[q_][0] = x_[0]/hs; q2g[q_][1] = x_[1]/hs;
									if(q2smooth[q_]){
#ifdef _OPENMP
#pragma omp atomic
#endif
										gfs_q2_smooth[q_]++;
									}
								}
								ibp->stress.gUxx = (postype)q2g[0][0]; ibp->stress.gUxy = (postype)q2g[0][1];
								ibp->stress.gUyx = (postype)q2g[1][0]; ibp->stress.gUyy = (postype)q2g[1][1];
								if(gfs_tcf_qbp()){ size_t id_ = (size_t)PINDX(ibp_rk4); if(id_ < qbp_cap){ qbp_g[4*id_] = q2g[0][0]; qbp_g[4*id_+1] = q2g[0][1]; qbp_g[4*id_+2] = q2g[1][0]; qbp_g[4*id_+3] = q2g[1][1]; qbp_ok[id_] = 2; int h_; for(h_ = 0; h_ < 3; h_++){ qbp_h[6*id_+h_] = q2H[0][h_]; qbp_h[6*id_+3+h_] = q2H[1][h_]; } } }
								ibp->stress.dPdx = (postype)q2g[2][0]; ibp->stress.dPdy = (postype)q2g[2][1];
								ibp->stress.dRhodx = (postype)q2g[3][0]; ibp->stress.dRhody = (postype)q2g[3][1];
								ibp->stress.divv = ibp->stress.gUxx + ibp->stress.gUyy;
								q2ok = 1;
							}
						}
					} else {
						postype invVol = 1.0/ibp_rk4->volume;
						ibp->stress.gUxx = sum_gUxx * invVol;
						ibp->stress.gUxy = sum_gUxy * invVol;
						ibp->stress.gUyx = sum_gUyx * invVol;
						ibp->stress.gUyy = sum_gUyy * invVol;
						ibp->stress.divv = ibp->stress.gUxx + ibp->stress.gUyy;
						ibp->stress.dPdx = sum_dPdx * invVol;
						ibp->stress.dPdy = sum_dPdy * invVol;
						ibp->stress.dRhodx = sum_dRhodx * invVol;
						ibp->stress.dRhody = sum_dRhody * invVol;
					}

					/* Hyperviscosity Pass 1 finalize: divide by V_i */
					{
						postype invVol_lap = 1.0/ibp_rk4->volume;
						ibp->stress.lap_vx = sum_lap_vx * invVol_lap;
						ibp->stress.lap_vy = sum_lap_vy * invVol_lap;
					}

					/* ====== Springel 2010 TVD slope limiter (av_mode==5) ======
					 * Pass A: find min/max of {rho, vx, vy, P} over all real
					 *         neighbors of cell i (wall ghosts excluded).
					 * Pass B: for each face centroid x_f, compute trial delta =
					 *         grad Q . (x_f - x_i) and tighten alpha_Q so the
					 *         reconstructed value Q_i + alpha_Q*delta stays in
					 *         [Q_min, Q_max].  Finally rescale the stored
					 *         gradients in-place. */
					int oracle_cell = 0;
					if(gfs_oracle_on()){
						Voro2D_Corner *tv = vorocorner;
						do { if(tv->upperrelated == VORO_ARC_EDGE) oracle_cell = 1; tv = tv->upperlink; } while(tv != vorocorner);
						if(oracle_cell && (gfs_oracle_on() == 1 || gfs_oracle_on() == 2)){
							double cx0 = gfs_env_d("EUNHA_KEPLER_CX", 0), cy0 = gfs_env_d("EUNHA_KEPLER_CY", 0), eps = gfs_env_d("EUNHA_KEPLER_EPS", 0.01);
							double x = (double)ibp->x - cx0, y = (double)ibp->y - cy0, r = sqrt(x*x + y*y);
							double Om = pow(r*r + eps*eps, -0.75), dOm = -1.5*r*pow(r*r + eps*eps, -1.75);
							ibp->stress.gUxx = (postype)(-y*dOm*x/r);
							ibp->stress.gUxy = (postype)(-Om - y*dOm*y/r);
							ibp->stress.gUyx = (postype)(Om + x*dOm*x/r);
							ibp->stress.gUyy = (postype)(x*dOm*y/r);
							ibp->stress.divv = ibp->stress.gUxx + ibp->stress.gUyy;
						}
					}
					if(av_mode == 5){
						postype rho_i = ibp->den, vx_i = ibp->vx, vy_i = ibp->vy, P_i = ibp->pressure;
						postype rho_max=rho_i, rho_min=rho_i;
						postype vx_max=vx_i, vx_min=vx_i;
						postype vy_max=vy_i, vy_min=vy_i;
						postype P_max=P_i, P_min=P_i;
						Voro2D_Corner *tmpA = vorocorner;
						do {
							if(tmpA->upperrelated >= 0){
								treevorostressrk4particletype *jL = (treevorostressrk4particletype*)(neighwork[tmpA->upperrelated].bp);
								int ghostL = (PINDX((treevorork4particletype*)jL) == MAX_INDEX);
								if(!ghostL){
									if(jL->den      > rho_max) rho_max = jL->den;
									if(jL->den      < rho_min) rho_min = jL->den;
									if(jL->vx       > vx_max)  vx_max  = jL->vx;
									if(jL->vx       < vx_min)  vx_min  = jL->vx;
									if(jL->vy       > vy_max)  vy_max  = jL->vy;
									if(jL->vy       < vy_min)  vy_min  = jL->vy;
									if(jL->pressure > P_max)   P_max   = jL->pressure;
									if(jL->pressure < P_min)   P_min   = jL->pressure;
								}
							}
							tmpA = tmpA->upperlink;
						} while(tmpA != vorocorner);

						if(gfs_vacbc_on()){
							/* a cell touching the vacuum arc: the boundary value 0 is in range */
							Voro2D_Corner *tv = vorocorner; int arc = 0;
							do { if(tv->upperrelated == VORO_ARC_EDGE) arc = 1; tv = tv->upperlink; } while(tv != vorocorner);
							if(arc){ rho_min = 0; P_min = 0; }
						}
						if(q2ok && gfs_lim_smooth_on()){
							/* smooth cells: allow the curvature over the stencil, max_j |e_j^T H e_j|/2 */
							double dlt[4] = {0, 0, 0, 0};
							Voro2D_Corner *tc = vorocorner;
							do {
								if(tc->upperrelated >= 0){
									double ex = (double)neighwork[tc->upperrelated].x, ey = (double)neighwork[tc->upperrelated].y;
									int q_;
									for(q_ = 0; q_ < 4; q_++){
										double c2 = 0.5*fabs(q2H[q_][0]*ex*ex + 2*q2H[q_][1]*ex*ey + q2H[q_][2]*ey*ey);
										if(c2 > dlt[q_]) dlt[q_] = c2;
									}
								}
								tc = tc->upperlink;
							} while(tc != vorocorner);
							if(q2smooth[0]){ vx_max += (postype)dlt[0]; vx_min -= (postype)dlt[0]; }
							if(q2smooth[1]){ vy_max += (postype)dlt[1]; vy_min -= (postype)dlt[1]; }
							if(q2smooth[2]){ P_max += (postype)dlt[2]; P_min -= (postype)dlt[2]; }
							if(q2smooth[3]){ rho_max += (postype)dlt[3]; rho_min -= (postype)dlt[3]; }
						}
						postype alpha_rho=1.0, alpha_vx=1.0, alpha_vy=1.0, alpha_P=1.0;
						/* GFS_LIM_PROJ=1 (Astra audit 11): rotation-invariant velocity limiter. One scalar
						 * alpha_v for the whole gradient tensor from the face-normal projections:
						 * delta = e.G r_f, range of v_k.e over S_f = {i, j, k-, k+} (real generators only). */
						double alpha_vp = 1.0; int lim_proj = gfs_lim_proj();
						Voro2D_Corner *tmpB = vorocorner;
						do {
							if(tmpB->upperrelated >= 0){
								Voro2D_Corner *tmpB2 = tmpB->upperlink;
								/* tmpB->x,tmpB->y are already RELATIVE to ibp
								   (Voro2D_FindVC uses Vec2DSub), so the face
								   midpoint IS the displacement from cell center. */
								postype dxf = 0.5*(tmpB->x + tmpB2->x);
								postype dyf = 0.5*(tmpB->y + tmpB2->y);
								if(gfs_lim_anchor_on()){
									/* positivity of rho and P where the flux is evaluated (the
									 * face centroid), separate from the monotonicity check below */
									double dRc = (double)ibp->stress.dRhodx*dxf + (double)ibp->stress.dRhody*dyf;
									double dPc = (double)ibp->stress.dPdx*dxf + (double)ibp->stress.dPdy*dyf;
									if(dRc < -1e-300){ double r_ = -(double)rho_i/dRc; if(r_ < 0) r_ = 0; if(r_ < alpha_rho) alpha_rho = (postype)r_; }
									if(dPc < -1e-300){ double r_ = -(double)P_i/dPc; if(r_ < 0) r_ = 0; if(r_ < alpha_P) alpha_P = (postype)r_; }
									/* the face point on the generator line: phi d e */
									double ex_ = (double)neighwork[tmpB->upperrelated].x, ey_ = (double)neighwork[tmpB->upperrelated].y;
									double d2_ = ex_*ex_ + ey_*ey_;
									if(d2_ > 0){
										double phi_ = 0.5 + ((double)ibp_rk4->w2 - (double)neighwork[tmpB->upperrelated].w2)/(2*d2_);
										dxf = (postype)(phi_*ex_); dyf = (postype)(phi_*ey_);
									}
								}

								postype drhoF = ibp->stress.dRhodx*dxf + ibp->stress.dRhody*dyf;
								postype dvxF  = ibp->stress.gUxx*dxf + ibp->stress.gUxy*dyf;
								postype dvyF  = ibp->stress.gUyx*dxf + ibp->stress.gUyy*dyf;
								postype dPF   = ibp->stress.dPdx*dxf + ibp->stress.dPdy*dyf;

								#define TVD_LIMIT(alpha, delta, Vi, Vmax, Vmin) do { \
									if(delta > 1e-30){ \
										postype r = (Vmax - Vi)/delta; \
										if(r < 0) r = 0; \
										if(r < alpha) alpha = r; \
									} else if(delta < -1e-30){ \
										postype r = (Vmin - Vi)/delta; \
										if(r < 0) r = 0; \
										if(r < alpha) alpha = r; \
									} \
								} while(0)

								if(!gfs_lim_pair_on()){
								TVD_LIMIT(alpha_rho, drhoF, rho_i, rho_max, rho_min);
								if(!lim_proj){
								TVD_LIMIT(alpha_vx,  dvxF,  vx_i,  vx_max,  vx_min);
								TVD_LIMIT(alpha_vy,  dvyF,  vy_i,  vy_max,  vy_min);
								}
								TVD_LIMIT(alpha_P,   dPF,   P_i,   P_max,   P_min);
								}
								if(lim_proj){
									double ex_ = (double)neighwork[tmpB->upperrelated].x, ey_ = (double)neighwork[tmpB->upperrelated].y, dn_ = sqrt(ex_*ex_ + ey_*ey_);
									if(dn_ > 0){
										ex_ /= dn_; ey_ /= dn_;
										double ui_ = (double)vx_i*ex_ + (double)vy_i*ey_, del_ = (double)dvxF*ex_ + (double)dvyF*ey_;
										double Um = ui_, Up = ui_; int kk_[3] = {tmpB->upperrelated, tmpB->lowerrelated, tmpB2->upperrelated}, q_;
										for(q_ = 0; q_ < 3; q_++){
											if(kk_[q_] < 0) continue;
											treevorork4particletype *nb_ = (treevorork4particletype*)neighwork[kk_[q_]].bp;
											if(!nb_) continue;
											double uk_ = (double)nb_->vx*ex_ + (double)nb_->vy*ey_;
											if(uk_ < Um) Um = uk_; if(uk_ > Up) Up = uk_;
										}
										double a_ = 1.0;
										if(del_ > 0) a_ = (Up - ui_)/del_; else if(del_ < 0) a_ = (Um - ui_)/del_;
										if(a_ < 0) a_ = 0; if(a_ < alpha_vp) alpha_vp = a_;
									}
								}
								#undef TVD_LIMIT
							}
							tmpB = tmpB->upperlink;
						} while(tmpB != vorocorner);

						if(lim_proj){ alpha_vx = (postype)alpha_vp; alpha_vy = (postype)alpha_vp; }
						/* Apply limiter in-place (safe: av_mode==5 does not use
						   stress.gUxx/gUxy/... for NS stress computation). */
						if(gfs_tq_on()){
							double cx, cy; tq_center(&cx, &cy);
							double rx = (double)ibp->x - cx, ry = (double)ibp->y - cy;
							int kb = tq_bin(sqrt(rx*rx + ry*ry));
							double am = fmin((double)alpha_vx, (double)alpha_vy);
#ifdef _OPENMP
#pragma omp atomic
#endif
							lq_n[kb] += 1;
#ifdef _OPENMP
#pragma omp atomic
#endif
							lq_asum[kb] += am;
							if(am < 0.5){
#ifdef _OPENMP
#pragma omp atomic
#endif
								lq_low[kb] += 1;
							}
						}
						if(gfs_lim_dump_on()){
							/* GFS_LIM_DUMP=1 (diagnostic): the code's per-component
							 * alphas next to a frame-invariant alternative alpha_nt
							 * that limits, at every face, the reconstructed normal and
							 * tangential velocity by the neighbours' values projected
							 * on that face's normal and tangent. */
							double W_ = gfs_env_d("GFS_EDGE_DUMP_W", 0.05), T1 = gfs_env_d("GFS_EDGE_DUMP_T1", 1.0), T2 = gfs_env_d("GFS_EDGE_DUMP_T2", 5.0);
							double cx = gfs_env_d("EUNHA_KEPLER_CX", 0), cy = gfs_env_d("EUNHA_KEPLER_CY", 0);
							double rx = (double)ibp->x - cx, ry = (double)ibp->y - cy, ri = sqrt(rx*rx + ry*ry);
							if(((gfs_tnow >= T1 && gfs_tnow < T1 + W_) || (gfs_tnow >= T2 && gfs_tnow < T2 + W_))
									&& ri < gfs_env_d("GFS_EDGE_DUMP_R", 0.56)){
								double a_nt = 1.0; int nlim_n = 0, nlim_t = 0, nf = 0;
								Voro2D_Corner *tf = vorocorner;
								do {
									if(tf->upperrelated >= 0){
										Voro2D_Corner *tf2 = tf->upperlink;
										double lx = (double)(tf2->x - tf->x), ly = (double)(tf2->y - tf->y), L = sqrt(lx*lx + ly*ly);
										if(L > 0){
											double tx = lx/L, ty = ly/L, nx = ty, ny = -tx;
											double xf = 0.5*((double)tf->x + (double)tf2->x), yf = 0.5*((double)tf->y + (double)tf2->y);
											double dvx = (double)ibp->stress.gUxx*xf + (double)ibp->stress.gUxy*yf;
											double dvy = (double)ibp->stress.gUyx*xf + (double)ibp->stress.gUyy*yf;
											double dn = dvx*nx + dvy*ny, dt_ = dvx*tx + dvy*ty;
											double vin = (double)vx_i*nx + (double)vy_i*ny, vit = (double)vx_i*tx + (double)vy_i*ty;
											double nmin = vin, nmax = vin, tmin = vit, tmax = vit;
											Voro2D_Corner *tn = vorocorner;
											do {
												if(tn->upperrelated >= 0){
													treevorostressrk4particletype *jn = (treevorostressrk4particletype*)(neighwork[tn->upperrelated].bp);
													if(PINDX((treevorork4particletype*)jn) != MAX_INDEX){
														double a = (double)jn->vx*nx + (double)jn->vy*ny, b = (double)jn->vx*tx + (double)jn->vy*ty;
														if(a < nmin) nmin = a;
														if(a > nmax) nmax = a;
														if(b < tmin) tmin = b;
														if(b > tmax) tmax = b;
													}
												}
												tn = tn->upperlink;
											} while(tn != vorocorner);
											double r_;
											if(dn > 1e-30){ r_ = (nmax - vin)/dn; if(r_ < a_nt){ a_nt = r_ < 0 ? 0 : r_; nlim_n++; } }
											else if(dn < -1e-30){ r_ = (nmin - vin)/dn; if(r_ < a_nt){ a_nt = r_ < 0 ? 0 : r_; nlim_n++; } }
											if(dt_ > 1e-30){ r_ = (tmax - vit)/dt_; if(r_ < a_nt){ a_nt = r_ < 0 ? 0 : r_; nlim_t++; } }
											else if(dt_ < -1e-30){ r_ = (tmin - vit)/dt_; if(r_ < a_nt){ a_nt = r_ < 0 ? 0 : r_; nlim_t++; } }
											nf++;
										}
									}
									tf = tf->upperlink;
								} while(tf != vorocorner);
#ifdef _OPENMP
#pragma omp critical (lim_dump)
#endif
								{
									static FILE *fl = NULL;
									if(!fl){ fl = fopen("limiter_edge.txt", "a"); if(fl) fprintf(fl, "# t id r alpha_rho alpha_vx alpha_vy alpha_P alpha_nt nlim_n nlim_t nfaces gUxx gUxy gUyx gUyy vx vy vx_min vx_max vy_min vy_max\n"); }
									if(fl) fprintf(fl, "%.5f %ld %.5f %.4f %.4f %.4f %.4f %.4f %d %d %d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
										gfs_tnow, (long)PINDX((treevorork4particletype*)ibp), ri, (double)alpha_rho, (double)alpha_vx, (double)alpha_vy, (double)alpha_P, a_nt,
										nlim_n, nlim_t, nf, (double)ibp->stress.gUxx, (double)ibp->stress.gUxy, (double)ibp->stress.gUyx, (double)ibp->stress.gUyy,
										(double)vx_i, (double)vy_i, (double)vx_min, (double)vx_max, (double)vy_min, (double)vy_max);
								}
							}
						}
						if(oracle_cell && (gfs_oracle_on() == 1 || gfs_oracle_on() == 3)){ alpha_vx = 1; alpha_vy = 1; }
						ibp->stress.dRhodx *= alpha_rho;
						ibp->stress.dRhody *= alpha_rho;
						ibp->stress.gUxx   *= alpha_vx;
						ibp->stress.gUxy   *= alpha_vx;
						ibp->stress.gUyx   *= alpha_vy;
						ibp->stress.gUyy   *= alpha_vy;
						ibp->stress.dPdx   *= alpha_P;
						ibp->stress.dPdy   *= alpha_P;
					}
				}
				if(neighbors_e){
					free(neighbors_e); free(neighwork_e);
					neighbors = neighbors_save;
					neighwork = neighwork_save;
					nneigh = nneigh_save;
				}
			}
			if(np>0) free(p);
			if(nneigh>0) {
				free(neighbors); free(neighwork);
			}
		}
		free(vorocorner);
	}
	if(gpass == 0 && gfs_periodic_complete()){
		double Lx_ = (double)(SIMBOX(simpar).x.max - SIMBOX(simpar).x.min), Ly_ = (double)(SIMBOX(simpar).y.max - SIMBOX(simpar).y.min);
		double ev = pc_vsum/(Lx_*Ly_) - 1;
		if(voro_sec_fail > 0){ fprintf(stderr, "[PCOMP] %ld cells exceed the searchable security radius (2 R_v > reach at R_max)\n", voro_sec_fail); }
		if(pc_ncap > 0 || fabs(ev) > 1e-12 || voro_sec_fail > 0){
			fprintf(stderr, "[PCOMP] geometry pass: cap cells %ld, volume sum/box - 1 = %.3e: gate failure\n", pc_ncap, ev);
			{	/* geometry-failure snapshot for offline reconstruction (Astra design 7, item 3) */
				FILE *fp_ = fopen("pcomp_fail.dat", "w");
				if(fp_){
					size_t ps_ = TVORORK4_DDINFO(simpar)[0].n_size; char *b_ = (char*)VORORK4_TBP(simpar); int i_;
					fprintf(fp_, "# box %.17g %.17g %.17g %.17g np %d vsum %.17g\n", (double)SIMBOX(simpar).x.min, (double)SIMBOX(simpar).x.max,
						(double)SIMBOX(simpar).y.min, (double)SIMBOX(simpar).y.max, (int)VORO_NP(simpar), pc_vsum);
					for(i_ = 0; i_ < VORO_NP(simpar); i_++){
						treevorork4particletype *q_ = (treevorork4particletype*)(b_ + i_*ps_);
						fprintf(fp_, "%ld %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", (long)PINDX(q_), (double)q_->x, (double)q_->y,
							(double)q_->volume, (double)q_->mass, (double)q_->vx, (double)q_->vy, (double)q_->ie);
					}
					fclose(fp_);
				}
			}
			exit(2);
		}
	}
	if(gpass == 0){
		/* pressure and sound speed from the new volumes (same rules as below) */
		int em0 = GAS_ENTROPY_MODE(simpar);
		postype Kf0 = GAS_K_FLOOR(simpar);
		double de0 = gfs_dual_energy_eta();
		for(i=0;i<nbp;i++){
			if(em0 == 1){
				if(!(bp[i].stress.K > 0) && bp[i].ie > 0 && bp[i].den > 0 && getenv("GFS_K_INIT"))
					bp[i].stress.K = (Gamma-1)*bp[i].ie/bp[i].volume/pow((double)bp[i].den, (double)Gamma);
				if(bp[i].stress.K < Kf0) bp[i].stress.K = Kf0;
				bp[i].pressure = bp[i].stress.K * pow((double)bp[i].den, (double)Gamma);
				bp[i].ie = bp[i].pressure * bp[i].volume / (Gamma-1);
			} else {
				bp[i].pressure = bp[i].ie/(bp[i].volume*sub_gfac((treevorork4particletype*)(bp+i)))*(Gamma-1);
				if(bp[i].pressure <= 0){
					double ie_before = (double)bp[i].ie;
					if(de0 > 0 && bp[i].stress.K > 0 && bp[i].den > 0)
						bp[i].pressure = bp[i].stress.K * pow((double)bp[i].den, (double)Gamma);
					else
						bp[i].pressure = 1e-6;
					bp[i].ie = bp[i].pressure * bp[i].volume*sub_gfac((treevorork4particletype*)(bp+i)) / (Gamma-1);
					if(i < VORO_NP(simpar)){
						gfs_stage_floor_l += (double)bp[i].ie - ie_before;
						gfs_stage_floor_n_l++; gfs_attempt_failed = 1; gfs_attempt_failed = 1;
						if(gfs_tq_on()){ double cx_, cy_; tq_center(&cx_, &cy_); sfl_bin[tq_bin(sqrt(pow((double)bp[i].x - cx_, 2) + pow((double)bp[i].y - cy_, 2)))] += 1; }
						gfs_sfl_log((treevorork4particletype*)(bp+i), ie_before);
					}
				}
			}
			bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
		}
		/* neighbour copies with the updated states */
		my_free(VORO_BASICCELL(simpar));
		my_free(VORORK4_TBPP(simpar));
		cells = (VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*mx*my));
		mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax, paddingAllTreeParticles);
	}
	}   /* gpass */

	/* Update pressure, csound, and NS stress tensor */
	int entropy_mode = GAS_ENTROPY_MODE(simpar);
	postype K_floor = GAS_K_FLOOR(simpar);
	double de_eta_stage = gfs_dual_energy_eta();
	for(i=0;i<nbp;i++){
		if(entropy_mode == 1){
			/* P recovered from entropy variable K = P/rho^gamma.
			 * Bounded by mass conservation -> structurally avoids ie heating runaway. */
			if(!(bp[i].stress.K > 0) && bp[i].ie > 0 && bp[i].den > 0 && getenv("GFS_K_INIT"))
				bp[i].stress.K = (Gamma-1)*bp[i].ie/bp[i].volume/pow((double)bp[i].den, (double)Gamma);
			if(bp[i].stress.K < K_floor) bp[i].stress.K = K_floor;
			bp[i].pressure = bp[i].stress.K * pow((double)bp[i].den, (double)Gamma);
			bp[i].ie = bp[i].pressure * bp[i].volume / (Gamma-1);
		} else {
			bp[i].pressure = bp[i].ie/(bp[i].volume*sub_gfac((treevorork4particletype*)(bp+i)))*(Gamma-1);
			if(bp[i].pressure <= 0){
				double ie_before = (double)bp[i].ie;
				/* GFS_DUAL_ENERGY: the adiabatic value, not 1e-6. */
				if(de_eta_stage > 0 && bp[i].stress.K > 0 && bp[i].den > 0)
					bp[i].pressure = bp[i].stress.K * pow((double)bp[i].den, (double)Gamma);
				else
					bp[i].pressure = 1e-6;
				bp[i].ie = bp[i].pressure * bp[i].volume*sub_gfac((treevorork4particletype*)(bp+i)) / (Gamma-1);
				/* Stage floor. This reset stays in the final ie (the RK4
				 * combination undoes the k's, not the reset), so it is an
				 * energy source. Booked in [RK4E] as sfl_cum. Local cells
				 * only; the ghosts are copies. */
				if(i < VORO_NP(simpar)){
					gfs_stage_floor_l += (double)bp[i].ie - ie_before;
					gfs_stage_floor_n_l++; gfs_attempt_failed = 1;
						if(gfs_tq_on()){ double cx_, cy_; tq_center(&cx_, &cy_); sfl_bin[tq_bin(sqrt(pow((double)bp[i].x - cx_, 2) + pow((double)bp[i].y - cy_, 2)))] += 1; }
						gfs_sfl_log((treevorork4particletype*)(bp+i), ie_before);
				}
			}
		}
		bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);

		/* Deviatoric NS stress convention: τ_code = -ν ρ S_dev
		   with S_dev = ∇v + ∇v^T - (2/d)·(∇·v) I.
		   We use d = 3 (3D slab / quasi-3D) — the (2/3)·divv trace
		   subtraction treats the 2D plane as a slab of a 3D fluid with
		   ∂_z = 0.  For truly 2D Navier-Stokes (d = 2) the coefficient
		   would be (2/2)·divv = divv.  Impact on nearly-incompressible
		   KH (divv ≈ 0) is negligible; for compressible RT/Sod-2D it
		   matters at the ~ν·divv level.  Do NOT change without
		   updating the documented convention in Docs/cfd_comp.tex.

		   Sign: τ_code is defined with the opposite sign of the
		   physical τ_NS so that the force loop accumulates +τ_NS·dS
		   (viscosity damps shear) via tau_dot_dS = -(τ_code·dS). */
		if(av_mode == 0 && nu_phys_val > 0){
			/* av_mode=0 + nu_phys: NS stress with fixed nu = nu_phys (no CD10) */
			postype divv = bp[i].stress.divv;
			bp[i].stress.tauxx = -nu_phys_val * bp[i].den * (2.0*bp[i].stress.gUxx - (2.0/3.0)*divv);
			bp[i].stress.tauxy = -nu_phys_val * bp[i].den * (bp[i].stress.gUxy + bp[i].stress.gUyx);
			bp[i].stress.tauyy = -nu_phys_val * bp[i].den * (2.0*bp[i].stress.gUyy - (2.0/3.0)*divv);
		} else if(av_mode == 1){
			/* NS-stress uses ONLY physical viscosity. CD10 enters as viscous *pressure*
			 * Π_CD10 = ½ α vsig ρ̄ (-Δv_n) in the face loop, not as NS-stress traction.
			 * Old form ν = max(nu_phys, α_CD·h·c_s) injected 1/Nx via h=sqrt(V), which
			 * combined with 1/m ~ Nx² in dK to drive K-runaway at Nx≥256. Matches
			 * av_mode>=2 design intent (CD10 as viscous pressure). */
			postype nu_phys = GAS_VISCOSITY(simpar);
			postype nu = (nu_phys > 0) ? nu_phys : 0;
			postype divv = bp[i].stress.divv;
			bp[i].stress.tauxx = -nu * bp[i].den * (2.0*bp[i].stress.gUxx - (2.0/3.0)*divv);
			bp[i].stress.tauxy = -nu * bp[i].den * (bp[i].stress.gUxy + bp[i].stress.gUyx);
			bp[i].stress.tauyy = -nu * bp[i].den * (2.0*bp[i].stress.gUyy - (2.0/3.0)*divv);
		} else if(av_mode >= 2){
			/* Design intent — NOT a bug.  av_mode>=2 uses CD10 (α_CD)
			   as the dissipation mechanism in three places:
			     (a) blend factor f_pq = α_max/cd_amax (getAccVoro2DBlend)
			     (b) HLLC viscous pressure Π_CD10 = ½α vsig ρ̄ (-Δv_n)
			     (c) weight (1-f_pq) for M(n,m) path
			   NS stress would double-count α_CD dissipation, so τ is
			   intentionally zeroed here as defensive cleanup.  The
			   force loop gate at getAccVoro2DBlend already excludes
			   av_mode>=2 from applying tau_dot_dS — this zeroing is
			   redundant with that gate but guards against stale reads
			   from diagnostics. */
			bp[i].stress.tauxx = 0;
			bp[i].stress.tauxy = 0;
			bp[i].stress.tauyy = 0;
		}
	}

	sedov_wsmooth_apply(simpar, xmin, ymin, xmax, ymax, boxsize, cellsize, mx, my,
		find2DNeighboringBP, find2DCellBP, Gamma);

#ifdef GPU_TESS_COMPARE
	/* === Compare GPU face CSR vs CPU fused face CSR === */
	if (_face_cmp_call <= 2 && MYID(simpar) == 0 && _gf_offset) {
		/* CPU fused faces are now in gpu_get_context()->h_faces */
		FaceCSR *cf = &_fcctx->h_faces;
		int cpu_nfaces = cf->n_faces_total;
		fprintf(stderr, "[FACE_CMP] cpu_nfaces=%d gpu_nfaces=%d\n", cpu_nfaces, _fc_nfaces);

		/* Compare per-particle face counts */
		int count_mismatch = 0, first_mm = -1;
		for(int ii=0; ii<nbp; ii++){
			int gpu_fc = _gf_offset[ii+1] - _gf_offset[ii];
			int cpu_fc = cf->face_offset[ii+1] - cf->face_offset[ii];
			if(gpu_fc != cpu_fc){
				if(first_mm < 0) first_mm = ii;
				count_mismatch++;
			}
		}
		fprintf(stderr, "[FACE_CMP] face_count_mismatch=%d (of %d particles) first=%d\n",
			count_mismatch, nbp, first_mm);

		/* For each particle, compare faces by matching neighbor_idx */
		int n_neigh_err=0, n_kp_err=0, n_km_err=0, n_ghost_err=0;
		int n_c1_err=0, n_c2_err=0;
		int first_neigh_err=-1;
		for(int ii=0; ii<nbp && ii<nbp; ii++){
			int gbase = _gf_offset[ii], gend = _gf_offset[ii+1];
			int cbase = cf->face_offset[ii], cend = cf->face_offset[ii+1];
			if(gend - gbase != cend - cbase) continue; /* skip count mismatches */
			int nf = gend - gbase;
			/* Match GPU faces to CPU faces by neighbor_idx */
			for(int gf=0; gf<nf; gf++){
				int gn = _gf_neigh[gbase+gf];
				/* Find matching CPU face */
				int found = 0;
				for(int cf2=0; cf2<nf; cf2++){
					if(cf->neighbor_idx[cbase+cf2] == gn){
						found = 1;
						int gi = gbase+gf, ci = cbase+cf2;
						/* Compare kp/km */
						if(_gf_kp[gi] != cf->kp_idx[ci]) n_kp_err++;
						if(_gf_km[gi] != cf->km_idx[ci]) n_km_err++;
						if(_gf_ghost[gi] != cf->is_ghost[ci]) n_ghost_err++;
						/* Compare corners (relative, should be close) */
						double dc1 = fabs(_gf_c1x[gi]-cf->c1x[ci])+fabs(_gf_c1y[gi]-cf->c1y[ci]);
						double dc2 = fabs(_gf_c2x[gi]-cf->c2x[ci])+fabs(_gf_c2y[gi]-cf->c2y[ci]);
						if(dc1 > 1e-10) n_c1_err++;
						if(dc2 > 1e-10) n_c2_err++;
						/* Print first few errors */
						if(n_kp_err+n_km_err<=5 && (_gf_kp[gi]!=cf->kp_idx[ci] || _gf_km[gi]!=cf->km_idx[ci])){
							fprintf(stderr, "[FACE_CMP] i=%d neigh=%d: gpu_kp=%d cpu_kp=%d gpu_km=%d cpu_km=%d\n",
								ii, gn, _gf_kp[gi], cf->kp_idx[ci], _gf_km[gi], cf->km_idx[ci]);
						}
						if(n_c1_err+n_c2_err<=3 && (dc1>1e-10 || dc2>1e-10)){
							fprintf(stderr, "[FACE_CMP] i=%d neigh=%d: corner err dc1=%g dc2=%g\n"
								"  gpu c1=(%g,%g) c2=(%g,%g)\n  cpu c1=(%g,%g) c2=(%g,%g)\n",
								ii, gn, dc1, dc2,
								_gf_c1x[gi],_gf_c1y[gi],_gf_c2x[gi],_gf_c2y[gi],
								cf->c1x[ci],cf->c1y[ci],cf->c2x[ci],cf->c2y[ci]);
						}
						break;
					}
				}
				if(!found){
					n_neigh_err++;
					if(first_neigh_err<0) first_neigh_err = ii;
				}
			}
		}
		fprintf(stderr, "[FACE_CMP] neigh_unmatched=%d kp_err=%d km_err=%d ghost_err=%d c1_err=%d c2_err=%d\n",
			n_neigh_err, n_kp_err, n_km_err, n_ghost_err, n_c1_err, n_c2_err);
		if(first_neigh_err>=0)
			fprintf(stderr, "[FACE_CMP] first_neigh_unmatched at particle %d\n", first_neigh_err);

		/* Print first particle's faces for eyeball check */
		if(nbp > 0){
			int gbase = _gf_offset[0], gend = _gf_offset[1];
			int cbase = cf->face_offset[0], cend = cf->face_offset[1];
			fprintf(stderr, "[FACE_CMP] particle 0: gpu_nf=%d cpu_nf=%d\n", gend-gbase, cend-cbase);
			for(int f=gbase;f<gend && f<gbase+4;f++)
				fprintf(stderr, "[FACE_CMP]   gpu f=%d: neigh=%d kp=%d km=%d ghost=%d c1=(%g,%g) c2=(%g,%g)\n",
					f-gbase, _gf_neigh[f], _gf_kp[f], _gf_km[f], _gf_ghost[f],
					_gf_c1x[f], _gf_c1y[f], _gf_c2x[f], _gf_c2y[f]);
			for(int f=cbase;f<cend && f<cbase+4;f++)
				fprintf(stderr, "[FACE_CMP]   cpu f=%d: neigh=%d kp=%d km=%d ghost=%d c1=(%g,%g) c2=(%g,%g)\n",
					f-cbase, cf->neighbor_idx[f], cf->kp_idx[f], cf->km_idx[f], cf->is_ghost[f],
					cf->c1x[f], cf->c1y[f], cf->c2x[f], cf->c2y[f]);
		}

		free(_gf_offset); free(_gf_c1x); free(_gf_c1y); free(_gf_c2x); free(_gf_c2y);
		free(_gf_neigh); free(_gf_kp); free(_gf_km); free(_gf_ghost);
	}
	} /* close GPU_TESS_COMPARE block */
#endif

#ifdef USE_CUDA
_gpu_tess_skip_cpu:
#endif
	my_free(VORO_BASICCELL(simpar));
#ifdef USE_CUDA
	if (GAS_GPU_ENABLED(simpar)) {
		if (!gpu_get_context()->tess_faces_on_device) {
			gpu_fused_end();
		}
		/* Keep VORORK4_TBPP alive — GPU wrapper will free it after SoA fill */
	} else {
		my_free(VORORK4_TBPP(simpar));
	}
#else
	my_free(VORORK4_TBPP(simpar));
#endif
}

/* ================================================================
 *  LagMFM helpers (av_mode=4): Wendland C2 kernel (2D) + 2x2 inverse
 *  Used by meshless path: kernel density, gradients, effective faces.
 *  Kernel support: r < 2h (q = r/h in [0,2])
 * ================================================================ */
static inline postype mfm_W_wendland2d(postype r, postype h)
{
	/* 2D Wendland C2 kernel, compact support 2h, normalization 7/(4 pi h^2). */
	if(h <= 0) return 0;
	postype q = r / h;
	if(q >= 2.0) return 0;
	postype one_minus_half_q = 1.0 - 0.5*q;
	postype t2 = one_minus_half_q*one_minus_half_q;
	postype t4 = t2*t2;
	return (7.0/(4.0*M_PI*h*h)) * t4 * (1.0 + 2.0*q);
}

static inline postype mfm_dWdr_wendland2d(postype r, postype h)
{
	if(h <= 0) return 0;
	postype q = r / h;
	if(q >= 2.0) return 0;
	postype u = 0.5 * q;
	postype v = 1.0 - u;
	postype v3 = v*v*v;
	postype C = 7.0/(4.0*M_PI*h*h);
	return -10.0 * C * u * v3 / h;
}

/* 2D cubic spline (M4), Hopkins (2015) fiducial kernel.
 * W(q) = σ · f(q), σ = 10/(7π h²), support q ∈ [0,2].
 *   f(q) = 1 - 1.5 q² + 0.75 q³     for 0 ≤ q ≤ 1
 *        = 0.25 (2 - q)³            for 1 ≤ q ≤ 2
 * With h = η·sqrt(V), η ≈ 1.128 gives N_nb ≈ 4π η² ≈ 16.
 */
static inline postype mfm_W_cspline2d(postype r, postype h)
{
	if(h <= 0) return 0;
	postype q = r / h;
	if(q >= 2.0) return 0;
	postype sigma = 10.0/(7.0*M_PI*h*h);
	if(q < 1.0) return sigma*(1.0 - 1.5*q*q + 0.75*q*q*q);
	postype two_q = 2.0 - q;
	return sigma * 0.25 * two_q*two_q*two_q;
}

static inline postype mfm_dWdr_cspline2d(postype r, postype h)
{
	/* dW/dr = (σ/h) · df/dq
	 *   df/dq = -3q + 2.25 q²       (q ≤ 1)
	 *         = -0.75 (2-q)²        (1 ≤ q ≤ 2)
	 */
	if(h <= 0) return 0;
	postype q = r / h;
	if(q >= 2.0) return 0;
	postype sigma = 10.0/(7.0*M_PI*h*h);
	if(q < 1.0) return (sigma/h)*(-3.0*q + 2.25*q*q);
	postype two_q = 2.0 - q;
	return -(sigma/h) * 0.75 * two_q*two_q;
}

static inline int mfm_invert2x2(postype Exx, postype Exy, postype Eyx, postype Eyy,
		postype *Ixx, postype *Ixy, postype *Iyx, postype *Iyy)
{
	/* D-B5 fix: Invert a 2x2 matrix [[Exx,Exy],[Eyx,Eyy]] -> store into Ixx..Iyy.
	 * Returns 1 on success, 0 if near-singular.
	 * On failure, writes a ZERO matrix (not identity).  Callers then get
	 *   ∇Q = (E_inv) · S = 0,
	 * which degrades the particle to first-order (central state) — a safe
	 * and self-degrading fallback.  Identity would have given badly scaled
	 * gradients (~ O(1) instead of O(1/h²)) and corrupted the force loop.
	 * We also use a relative tolerance ‖E‖² · eps² to avoid triggering on
	 * legitimately small but non-singular E matrices in low-density regions.
	 */
	postype det = Exx*Eyy - Exy*Eyx;
	postype adet = det >= 0 ? det : -det;
	postype nrm2 = Exx*Exx + Exy*Exy + Eyx*Eyx + Eyy*Eyy;
	postype tol  = 1.0e-14 * nrm2 + 1.0e-30;
	if(adet < tol){
		*Ixx = 0.0; *Ixy = 0.0;
		*Iyx = 0.0; *Iyy = 0.0;
		return 0;
	}
	postype invdet = 1.0/det;
	*Ixx =  Eyy*invdet;
	*Ixy = -Exy*invdet;
	*Iyx = -Eyx*invdet;
	*Iyy =  Exx*invdet;
	return 1;
}

/* LagMFM tuning knobs (aligned with Hopkins 2015, MNRAS 450, 53).
 *
 *  LAGMFM_ETA         : h_i = eta * sqrt(V_i).  Cubic spline M4 kernel on
 *                       glass IC (Hopkins 2015 fiducial).  eta = 1.128 gives
 *                       N_nb ≈ 4π·eta² ≈ 16.
 *  LAGMFM_H_ITER      : max Newton iterations to converge h_i with
 *                       V_i = 1/n_i(h_i) self-consistent.  50 follows
 *                       Hopkins; bisection fallback triggers on oscillation.
 *  LAGMFM_H_TOL       : relative convergence tolerance on h_i.
 *  LAGMFM_PURE_GIZMO  : 1 = disable non-GIZMO extensions (Monaghan AV,
 *                       NS viscosity, CD10 shock switch) in av_mode=4;
 *                       0 = enable hybrid add-ons (LagMFM+NS).
 *  LAGMFM_WALL_NOSLIP : 0 = free-slip wall, 1 = no-slip wall.
 */
#define LAGMFM_ETA           ((postype)1.4)
#define LAGMFM_H_ITER        50
#define LAGMFM_H_TOL         ((postype)0.001)
#define LAGMFM_PURE_GIZMO    1
#define LAGMFM_WALL_NOSLIP   0

/* ================================================================
 *  updateDenW2Pressure2D_LagMFM (av_mode=4):
 *    Kernel-based density and matrix-weighted gradients (Hopkins 2015).
 *    Replaces Voronoi tessellation with a direct neighbor walk on the
 *    existing cell linked list.  Produces:
 *      - ibp->den, ibp->volume, ibp->csound, ibp->pressure
 *      - stress.gUxx..gUyy, divv, dPdx, dPdy
 *      - stress.E_inv_xx..E_inv_yy, stress.h_mfm
 *    NS stress tau is filled the same way as av_mode=1.
 * ================================================================ */
#ifdef USE_CUDA
/* w2 refresh helper — called from integrator AFTER GPU density kernel,
   which already computed w2ceil via nearest-neighbor in the neighbor loop.
   This handles w2 = f(P) from density results + applyW2Controls. */
static void lagmfm_w2_post_density(SimParameters *simpar, postype Dtime)
{
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int i;
	for(i=0;i<nbp;i++){
		if(GAS_Kappa(simpar) < 0){
			bp[i].w2 = -GAS_Kappa(simpar);
		} else if(GAS_Kappa(simpar) > 0){
			bp[i].w2hydro = getw2forHydroParticle(simpar,(treevorork4particletype*)(bp+i),Dtime);
			applyW2Controls(simpar, (treevorork4particletype*)(bp+i), Dtime);
		} else {
			/* Kappa=0: w2ceil from GPU density kernel, clip w2 */
			bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
		}
	}
}
#endif
void updateDenW2Pressure2D_LagMFM(
		SimParameters *simpar,
		postype xmin, postype ymin, postype xmax, postype ymax,
		postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		postype Dtime
		){
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype nu_phys_val = GAS_VISCOSITY(simpar);
#if LAGMFM_PURE_GIZMO
	(void)nu_phys_val;   /* inviscid Euler → unused */
#endif

	/* Refresh w2 for Laguerre weight (reused for dpq + w2ceil tracking even
	 * though LagMFM uses h_mfm, not w2, as kernel bandwidth). */
	det2d_dpqRK4(simpar, paddingAllTreeParticles);
	int i;
	for(i=0;i<nbp;i++){
		if(GAS_Kappa(simpar) < 0){
			bp[i].w2 = -GAS_Kappa(simpar);
		} else if(GAS_Kappa(simpar) > 0){
			bp[i].w2hydro = getw2forHydroParticle(simpar,(treevorork4particletype*)(bp+i),Dtime);
			applyW2Controls(simpar, (treevorork4particletype*)(bp+i), Dtime);
		}
	}

	/* Build cell linked list */
	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	int mx, my;
	BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
	BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	CellType *cells = (VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*mx*my));
	mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax, paddingAllTreeParticles);

	postype invcs = 1.0/cellsize;
	postype eta   = LAGMFM_ETA;

	/* -----------------------------------------------------------------
	 * Pass 1: kernel density + matrix-weighted gradients
	 * -----------------------------------------------------------------
	 *   n_i        = Σ_j  W(r_ij, h_i)
	 *   V_i        = 1 / n_i
	 *   ρ_i        = m_i · n_i
	 *   E_i^{αβ}   = Σ_j W_ij (x_j-x_i)^α (x_j-x_i)^β
	 *   ψ_j(x_i)   = W_ij · Σ_β E_inv^{αβ} (x_j-x_i)^β   (used in force)
	 *   ∇Q|_i^α    = Σ_β E_inv^{αβ} · Σ_j W_ij (Q_j-Q_i) (x_j-x_i)^β
	 *
	 * h_i is set from the previous volume:   h_i = eta * sqrt(V_i_prev).
	 * On first call V_i comes from kh_findVol (Voronoi) so the
	 * bootstrap is well behaved.
	 * ----------------------------------------------------------------- */
	int iy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(iy=0; iy<my; iy++){
		int ix;
		for(ix=0; ix<mx; ix++){
			struct linkedlisttype *li = cells[ix + mx*iy].link;
			while(li){
				treevorostressrk4particletype *ibp =
						(treevorostressrk4particletype*)li;
				int iwall = (PINDX((treevorork4particletype*)ibp) == MAX_INDEX);
				if(iwall || IS_FLAG_ON(li,BoundaryGhostflag)){
					li = li->next;
					continue;
				}

				/* D-B1 fix: adaptive h_i via Newton iteration.
				 * Start with h_i = eta * sqrt(V_prev), then iterate:
				 *   1) compute n_i = Σ_j W(r_ij, h_i) and moments E, S
				 *   2) V_new = 1/n_i, h_new = eta * sqrt(V_new)
				 *   3) stop if |h_new - h_i| < tol * h_i (converged)
				 * Only the LAST iteration's E, S values are used downstream,
				 * so on non-convergence we accept the final iteration's state.
				 * nbr (cell walk radius) is recomputed each iteration since
				 * a grown h_i needs a wider stencil. */
				postype Vprev = ibp->volume > 0 ? ibp->volume : cellsize*cellsize;
				postype h_i = eta * sqrt(Vprev);

				postype n_i = 0;
				postype Exx=0, Exy=0, Eyy=0;
				postype Sxx_vx=0, Sxy_vx=0;
				postype Sxx_vy=0, Sxy_vy=0;
				postype Sxx_P =0, Sxy_P =0;
				postype Sxx_Rho=0, Sxy_Rho=0;
				/* Per-particle neighbor min/max for Hopkins (2015) B31 limiter. */
				postype P_min = ibp->pressure, P_max = ibp->pressure;
				postype vx_min = ibp->vx, vx_max = ibp->vx;
				postype vy_min = ibp->vy, vy_max = ibp->vy;
				postype rho_min = ibp->den, rho_max = ibp->den;

				/* Newton iteration with bisection fallback (D-B1 / Hopkins).
				 * Newton: h_new = eta*sqrt(1/n(h)).  If Newton oscillates,
				 * fall back to bisection on the residual R(h) = h - eta/√n(h).
				 * h_lo, h_hi bracket is maintained from the residual sign. */
				postype h_lo = 0, h_hi = 0;
				postype h_prev = h_i, dh_prev = 0;
				int h_iter_done = 0;
				int use_bisect = 0;
				int h_it;
				postype rho_ker = 0;
				for(h_it = 0; h_it < LAGMFM_H_ITER; h_it++){
					/* Reset accumulators each iteration. */
					n_i = 0;
					rho_ker = 0;
					Exx = Exy = Eyy = 0;
					Sxx_vx = Sxy_vx = 0;
					Sxx_vy = Sxy_vy = 0;
					Sxx_P  = Sxy_P  = 0;
					Sxx_Rho = Sxy_Rho = 0;
					P_min = P_max = ibp->pressure;
					vx_min = vx_max = ibp->vx;
					vy_min = vy_max = ibp->vy;
					rho_min = rho_max = ibp->den;

					postype two_h = 2.0 * h_i;
					postype two_h2 = two_h * two_h;
					int nbr = (int)ceil(two_h * invcs);
					if(nbr < 1) nbr = 1;
					if(nbr > 10) nbr = 10;

					int jy;
					for(jy = iy-nbr; jy <= iy+nbr; jy++){
						if(jy < 0 || jy >= my) continue;
						int jx;
						for(jx = ix-nbr; jx <= ix+nbr; jx++){
							if(jx < 0 || jx >= mx) continue;
							struct linkedlisttype *lj = cells[jx + mx*jy].link;
							while(lj){
								treevorostressrk4particletype *jbp =
										(treevorostressrk4particletype*)lj;
								postype dx = jbp->x - ibp->x;
								postype dy = jbp->y - ibp->y;
								postype r2 = dx*dx + dy*dy;
								if(r2 < two_h2){
									postype r = sqrt(r2);
									postype Wij = mfm_W_wendland2d(r, h_i);
									n_i += Wij;
									rho_ker += Wij * jbp->mass;
									/* E matrix is symmetric in 2D */
									Exx += Wij * dx * dx;
									Exy += Wij * dx * dy;
									Eyy += Wij * dy * dy;
									/* First moments of (Q_j - Q_i) */
									postype dvx = jbp->vx - ibp->vx;
									postype dvy = jbp->vy - ibp->vy;
									postype dP   = jbp->pressure - ibp->pressure;
									postype dRho = jbp->den      - ibp->den;
									Sxx_vx += Wij * dvx * dx;
									Sxy_vx += Wij * dvx * dy;
									Sxx_vy += Wij * dvy * dx;
									Sxy_vy += Wij * dvy * dy;
									Sxx_P  += Wij * dP  * dx;
									Sxy_P  += Wij * dP  * dy;
									Sxx_Rho+= Wij * dRho* dx;
									Sxy_Rho+= Wij * dRho* dy;
									/* Hopkins B31 neighbor extrema */
									if(jbp->pressure < P_min) P_min = jbp->pressure;
									if(jbp->pressure > P_max) P_max = jbp->pressure;
									if(jbp->vx < vx_min) vx_min = jbp->vx;
									if(jbp->vx > vx_max) vx_max = jbp->vx;
									if(jbp->vy < vy_min) vy_min = jbp->vy;
									if(jbp->vy > vy_max) vy_max = jbp->vy;
									if(jbp->den < rho_min) rho_min = jbp->den;
									if(jbp->den > rho_max) rho_max = jbp->den;
								}
								lj = lj->next;
							}
						}
					}

					if(n_i <= 0) break;  /* isolated */

					postype V_new = 1.0 / n_i;
					postype h_target = eta * sqrt(V_new);
					postype R = h_i - h_target;   /* residual */
					/* Maintain bisection bracket from residual sign. */
					if(R > 0){ h_hi = h_i; }
					else     { h_lo = h_i; }

					postype h_new;
					if(use_bisect && h_lo > 0 && h_hi > 0){
						h_new = 0.5*(h_lo + h_hi);
					} else {
						h_new = h_target;
						/* Detect Newton oscillation: step reversal of similar
						 * magnitude → switch to bisection for the rest. */
						postype dh = h_new - h_i;
						if(h_it > 0 && dh*dh_prev < 0 &&
						   fabs(dh) > 0.5*fabs(dh_prev) &&
						   h_lo > 0 && h_hi > 0){
							use_bisect = 1;
							h_new = 0.5*(h_lo + h_hi);
						}
						dh_prev = dh;
					}
					postype adh = fabs(h_new - h_i);
					h_prev = h_i;
					h_i = h_new;
					if(adh < LAGMFM_H_TOL * h_i){
						h_iter_done = 1;
						break;
					}
				}
				(void)h_prev;
				(void)h_iter_done;  /* reserved for optional diagnostics */

				/* Finalize density / volume */
				if(n_i <= 0){
					/* Isolated particle — retain previous state */
					li = li->next;
					continue;
				}
				postype Vi = 1.0 / n_i;
				ibp->volume = Vi;
				/* Mass-weighted density (correct for unequal-mass IC).
				 * ρ_i = Σ_j m_j · W_ij; reduces to m·n_i when all m_j equal. */
				ibp->den    = rho_ker;

				/* Invert E to get the matrix-weighted gradient operator */
				postype Einv_xx, Einv_xy, Einv_yx, Einv_yy;
				mfm_invert2x2(Exx, Exy, Exy, Eyy,
						&Einv_xx, &Einv_xy, &Einv_yx, &Einv_yy);

				ibp->stress.E_inv_xx = Einv_xx;
				ibp->stress.E_inv_xy = Einv_xy;
				ibp->stress.E_inv_yx = Einv_yx;
				ibp->stress.E_inv_yy = Einv_yy;
				ibp->stress.h_mfm    = h_i;

				/* Gradients:  ∇Q|_i^α = Σ_β (E_inv)^{αβ} Σ_j W_ij dQ_ij dx_ij^β
				 * (no 1/n_i factor — it cancels between E and S because we
				 *  dropped the common 1/n_i prefactor from both.) */
				postype gUxx = Einv_xx*Sxx_vx + Einv_xy*Sxy_vx;
				postype gUxy = Einv_yx*Sxx_vx + Einv_yy*Sxy_vx;
				postype gUyx = Einv_xx*Sxx_vy + Einv_xy*Sxy_vy;
				postype gUyy = Einv_yx*Sxx_vy + Einv_yy*Sxy_vy;
				postype dPdx = Einv_xx*Sxx_P  + Einv_xy*Sxy_P;
				postype dPdy = Einv_yx*Sxx_P  + Einv_yy*Sxy_P;
				postype dRhodx = Einv_xx*Sxx_Rho + Einv_xy*Sxy_Rho;
				postype dRhody = Einv_yx*Sxx_Rho + Einv_yy*Sxy_Rho;

				/* Hopkins (2015) B31 gradient-scaling limiter with ADAPTIVE β_i.
				 *
				 *   β_i = max(β_min, β_max · min(1, N_crit/N_cond))
				 *   α_Q = min(1, β_i·(Q_max-Q_i)/(|∇Q|·h_i), β_i·(Q_i-Q_min)/(|∇Q|·h_i))
				 *
				 * N_cond = λ_max/λ_min of E (2×2 symmetric PD): large when
				 * neighborhood is degenerate (e.g. rolled-up vortex edges).
				 * In that regime β drops toward β_min=1 → tight limiter
				 * → suppresses grid-noise/sub-mode amplification.
				 * Matches Hopkins Appendix B (β_min=1, β_max=2, N_crit=10).
				 */
				postype Ncond_i;
				{
					postype tr  = Exx + Eyy;
					postype det = Exx*Eyy - Exy*Exy;
					postype disc = (Exx-Eyy)*(Exx-Eyy) + 4.0*Exy*Exy;
					postype sqrtdisc = sqrt(disc > 0 ? disc : 0);
					postype lam_max = 0.5*(tr + sqrtdisc);
					postype lam_min = 0.5*(tr - sqrtdisc);
					if(det <= 0 || lam_min <= 0)  Ncond_i = 1.0e30;
					else                          Ncond_i = lam_max / lam_min;
				}
				{
					postype recon_scale = h_i;
					postype N_crit = (postype)10.0;
					postype beta_lim = fmax((postype)1.0,
					                        (postype)2.0 * fmin((postype)1.0, N_crit / Ncond_i));
					/* Pressure */
					postype gmag_P = sqrt(dPdx*dPdx + dPdy*dPdy) * recon_scale;
					if(gmag_P > 0){
						postype a_up = beta_lim * (P_max - ibp->pressure) / gmag_P;
						postype a_dn = beta_lim * (ibp->pressure - P_min) / gmag_P;
						postype alpha = fmin(1.0, fmin(a_up, a_dn));
						if(alpha < 0) alpha = 0;
						dPdx *= alpha; dPdy *= alpha;
					}
					/* vx */
					postype gmag_vx = sqrt(gUxx*gUxx + gUxy*gUxy) * recon_scale;
					if(gmag_vx > 0){
						postype a_up = beta_lim * (vx_max - ibp->vx) / gmag_vx;
						postype a_dn = beta_lim * (ibp->vx - vx_min) / gmag_vx;
						postype alpha = fmin(1.0, fmin(a_up, a_dn));
						if(alpha < 0) alpha = 0;
						gUxx *= alpha; gUxy *= alpha;
					}
					/* vy */
					postype gmag_vy = sqrt(gUyx*gUyx + gUyy*gUyy) * recon_scale;
					if(gmag_vy > 0){
						postype a_up = beta_lim * (vy_max - ibp->vy) / gmag_vy;
						postype a_dn = beta_lim * (ibp->vy - vy_min) / gmag_vy;
						postype alpha = fmin(1.0, fmin(a_up, a_dn));
						if(alpha < 0) alpha = 0;
						gUyx *= alpha; gUyy *= alpha;
					}
					/* ρ */
					postype gmag_rho = sqrt(dRhodx*dRhodx + dRhody*dRhody) * recon_scale;
					if(gmag_rho > 0){
						postype a_up = beta_lim * (rho_max - ibp->den) / gmag_rho;
						postype a_dn = beta_lim * (ibp->den - rho_min) / gmag_rho;
						postype alpha = fmin(1.0, fmin(a_up, a_dn));
						if(alpha < 0) alpha = 0;
						dRhodx *= alpha; dRhody *= alpha;
					}
				}

				ibp->stress.gUxx = gUxx;
				ibp->stress.gUxy = gUxy;
				ibp->stress.gUyx = gUyx;
				ibp->stress.gUyy = gUyy;
				ibp->stress.divv = gUxx + gUyy;
				ibp->stress.dPdx = dPdx;
				ibp->stress.dPdy = dPdy;
				ibp->stress.dRhodx = dRhodx;
				ibp->stress.dRhody = dRhody;

				li = li->next;
			}
		}
	}

	/* -----------------------------------------------------------------
	 * Pass 2: pressure, sound speed, NS stress tensor (τ = -ν ρ S_dev)
	 *
	 * A-M1 analog: deviatoric NS stress with (2/d)·divv trace subtraction.
	 *   We use d = 3 (3D slab / quasi-3D convention).  For truly 2D NS
	 *   the coefficient would be (2/2)·divv.  Impact on KH (divv ≈ 0) is
	 *   negligible, but matters at the ~ν·divv level for compressible
	 *   RT/Sod-2D.  Convention matches av_mode=1 (see updateDenW2Pressure2DBlend).
	 * Sign: τ_code = -ν ρ S_dev (negative of physical τ_NS).  The force
	 *   loop computes τ_NS·A via tau_A = -(τ_code·A).
	 *
	 * A-M2 analog: stride-aware access via TVORORK4_DDINFO.  bp is typed
	 *   treevorostressrk4particletype (av_mode >= 1 allocates this), so a
	 *   plain bp[i] works, but we honor n_size for future-proofing.
	 * ----------------------------------------------------------------- */
	{
		size_t p_size2 = TVORORK4_DDINFO(simpar)[0].n_size;
		char *bp_raw2  = (char*)VORORK4_TBP(simpar);
		for(i=0;i<nbp;i++){
			treevorostressrk4particletype *bpi =
				(treevorostressrk4particletype*)(bp_raw2 + i*p_size2);
			bpi->pressure = bpi->ie/bpi->volume*(Gamma-1);
			if(bpi->pressure <= 0){
				bpi->pressure = 1e-6;
				bpi->ie = bpi->pressure * bpi->volume / (Gamma-1);
			}
			bpi->csound = sqrt(Gamma*bpi->pressure/bpi->den);

#if LAGMFM_PURE_GIZMO
			/* Pure Hopkins (2015): inviscid Euler — τ = 0. */
			bpi->stress.tauxx = 0;
			bpi->stress.tauxy = 0;
			bpi->stress.tauyy = 0;
#else
			/* NS-stress: physical viscosity only. CD10 NS-stress path removed
			 * (caused Nx² K-runaway via h·c_s); CD10 acts as viscous pressure. */
			postype nu = (nu_phys_val > 0) ? nu_phys_val : 0;
			postype divv = bpi->stress.divv;
			bpi->stress.tauxx = -nu * bpi->den * (2.0*bpi->stress.gUxx - (2.0/3.0)*divv);
			bpi->stress.tauxy = -nu * bpi->den * (bpi->stress.gUxy + bpi->stress.gUyx);
			bpi->stress.tauyy = -nu * bpi->den * (2.0*bpi->stress.gUyy - (2.0/3.0)*divv);
#endif
		}
	}

	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
}

/* ================================================================
 *  hllc_face_2d: TRUE HLLC approximate Riemann solver
 *    Ports GIZMO reimann.h:527-593 (get_wavespeeds_and_pressure_star).
 *    Returns star-region pressure P_M and normal velocity S_M (contact).
 *
 *    Three-stage fallback chain:
 *      1) Primary: Gaburov/Davis wave speeds + HLLC star state
 *      2) Fallback: Roe-averaged wave speeds (better accuracy near shocks)
 *      3) Last resort: Primitive-variable Rusanov-like bounds
 *
 *    Previous impl was TSRS (Toro 1999, Eq.10.72) — less dissipative at
 *    contacts and produces faster-than-Hopkins KH growth.  Proper HLLC
 *    uses S_L/S_R wave-speed upwinding which adds correct contact
 *    dissipation matching GIZMO behavior.
 * ================================================================ */
/* Phase-1 Sedov probe (env SEDOV_PHASE1=1). Fixed particle mass.
 * Extreme faces, including an ambient-pressure floor, take P* and v*
 * from the HLL average state. Dissipative dK heating is off, and the
 * kick stores that face work in the total energy, then sets
 * ie = E − kinetic. Gated so other tests keep the existing path. */
int cyl_wall_query(postype x, postype y, postype *nx, postype *ny);

static int sedov_phase1_on(void){
	static int cached = -1;
	if(cached < 0){
		const char *e = getenv("SEDOV_PHASE1");
		cached = (e && e[0] == '1') ? 1 : 0;
	}
	return cached;
}

static inline void hll_star_state(
		postype rhoL, postype pL, postype vnL, postype cL,
		postype rhoR, postype pR, postype vnR, postype cR,
		postype Gamma,
		postype *pstar, postype *vnstar)
{
	postype SL = fmin(vnL - cL, vnR - cR);
	postype SR = fmax(vnL + cL, vnR + cR);
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
	postype gm1 = Gamma - (postype)1;
	if(gm1 < (postype)1e-8) gm1 = (postype)1e-8;
	postype rhoLs = fmax(rhoL, (postype)1e-30);
	postype rhoRs = fmax(rhoR, (postype)1e-30);
	/* EL, ER are specific. The HLL average is on the conserved densities. */
	postype EL = pL/(gm1*rhoLs) + (postype)0.5*vnL*vnL;
	postype ER = pR/(gm1*rhoRs) + (postype)0.5*vnR*vnR;
	postype EdenL = rhoL * EL, EdenR = rhoR * ER;
	postype denom = SR - SL;
	if(fabs(denom) < (postype)1e-30){
		*pstar = (postype)0.5*(pL+pR);
		*vnstar = (postype)0.5*(vnL+vnR);
		return;
	}
	postype FrL = rhoL*vnL, FrR = rhoR*vnR;
	postype FmL = rhoL*vnL*vnL + pL, FmR = rhoR*vnR*vnR + pR;
	postype FeL = vnL*(EdenL + pL), FeR = vnR*(EdenR + pR);
	postype rhoS = (SR*rhoR - SL*rhoL - (FrR - FrL)) / denom;
	postype momS = (SR*(rhoR*vnR) - SL*(rhoL*vnL) - (FmR - FmL)) / denom;
	postype ES   = (SR*EdenR - SL*EdenL - (FeR - FeL)) / denom;
	if(!(rhoS > (postype)1e-30) || isnan(rhoS) || isinf(rhoS)
			|| isnan(momS) || isinf(momS) || isnan(ES) || isinf(ES)){
		*pstar = (postype)0.5*(pL+pR);
		*vnstar = (postype)0.5*(vnL+vnR);
		return;
	}
	postype vnS = momS / rhoS;
	postype eint = ES/rhoS - (postype)0.5*vnS*vnS;
	postype PS = gm1 * rhoS * eint;
	if(!(PS > 0) || isnan(PS) || isinf(PS) || isnan(vnS) || isinf(vnS)){
		*pstar = (postype)0.5*(pL+pR);
		*vnstar = (postype)0.5*(vnL+vnR);
		return;
	}
	*pstar = PS;
	*vnstar = vnS;
}

static int phase1_ie_clips = 0;
/* KDK clip bookkeeping, summed over the kicks of one step on this rank.
 * clip_dke is the kinetic energy removed by the velocity rescale.
 * clip_inj is the energy added when E < 0 is raised to 0.
 * clip_dL is the change of angular momentum about the Kepler centre. */
static double phase1_clip_dke = 0, phase1_clip_inj = 0, phase1_clip_dL = 0;
static long phase1_clip_events = 0;
/* Central mass GM=1. Active only for EUNHA_IC=kepler.
 * Default is the GIZMO point-mass acceleration -r_hat/R^2.
 * EUNHA_KEPLER_EPS>0 switches on the softened force, which must
 * match the circular speed set in the Kepler initial condition. */
/* EUNHA_KEPLER_EVOLVE_R = R_e (Kepler IC): particles at R >= R_e from the
 * Kepler centre are not evolved (all RK4 increments zero); they stay in
 * the tessellation as neighbours. As the original LagEunha Kepler code
 * ("evolve only for R < Lx/2", R_e = 2.5). */
/* EUNHA_KEPLER_SOFT=spline (with EUNHA_KEPLER_RC, default 0.25): the point mass is spread with the cubic-spline
 * kernel of radius r_c (GADGET/AREPO softening): g = M(<r)/r^2, exactly 1/r^2 for r >= r_c, C^2-continuous,
 * g(0) = 0. Enclosed mass fraction of the 3D cubic spline (u = r/r_c). */
static double kepler_spline_mass(double u){
	if(u >= 1) return 1.0;
	if(u <= 0.5) return u*u*u*(32.0/3 - 192.0/5*u*u + 32*u*u*u);
	return -1.0/15 + u*u*u*(64.0/3 - 48*u + 192.0/5*u*u - 32.0/3*u*u*u);
}
static int gfs_kepler_soft_spline(void){ static int on = -1; if(on < 0){ const char *s = getenv("EUNHA_KEPLER_SOFT"); on = (s && strcmp(s, "spline") == 0) ? 1 : 0; } return on; }
static double gfs_kepler_rc(void){ static double rc = -1; if(rc < 0){ const char *s = getenv("EUNHA_KEPLER_RC"); rc = (s && s[0]) ? atof(s) : 0.25; } return rc; }
/* EUNHA_KEPLER_FREEZE_RIN = R_f: particles at R < R_f are frozen too (static inner background, Pakmor et al. 2016 App. B). */
static int gfs_kepler_frozen(postype x, postype y){
	static double re = -1, cx = 0, cy = 0, rf = 0;
	if(re < 0){ const char *s = getenv("EUNHA_KEPLER_EVOLVE_R"), *a = getenv("EUNHA_KEPLER_CX"), *b = getenv("EUNHA_KEPLER_CY"), *f = getenv("EUNHA_KEPLER_FREEZE_RIN");
		re = (s && s[0]) ? atof(s) : 0; cx = (a && a[0]) ? atof(a) : 2; cy = (b && b[0]) ? atof(b) : 2; rf = (f && f[0]) ? atof(f) : 0; }
	if(re <= 0 && rf <= 0) return 0;
	double dx = (double)x - cx, dy = (double)y - cy, r2 = dx*dx + dy*dy;
	if(re > 0 && r2 >= re*re) return 1;
	if(rf > 0 && r2 < rf*rf) return 1;
	return 0;
}
/* GFS_HYDRO_OFF=1 (test): RK4 stages use gravity only; hydro accelerations
 * and internal-energy rates are multiplied by 0 (tessellation still runs). */
static double gfs_hydro_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_HYDRO_OFF"); on = (s && atoi(s) > 0) ? 0 : 1; }
	return (double)on;
}
/* EUNHA_KEPLER_FORCE=gizmo: the GIZMO Keplerian test-problem force */
static int gfs_kepler_force_gizmo(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("EUNHA_KEPLER_FORCE"); on = (s && strcmp(s, "gizmo") == 0) ? 1 : 0; }
	return on;
}
/* EUNHA_HILL_OMEGA = Omega0 > 0: local shearing-sheet (Hill) frame about
 * x_c = EUNHA_HILL_XC (default 3). Tidal acceleration 3 Omega0^2 (x - x_c)
 * here (potential -3/2 Omega0^2 (x - x_c)^2 in kepler_phi) and the
 * Coriolis acceleration -2 Omega0 z x v in hill_coriolis(); equilibrium
 * v_y = -3/2 Omega0 (x - x_c). Test problem only (cold sheared slab). */
static postype hill_omega(postype *xc){
	static int cached = 0; static postype om = 0, xc0 = 3;
	if(!cached){ const char *s = getenv("EUNHA_HILL_OMEGA"), *c = getenv("EUNHA_HILL_XC");
		om = (s && s[0]) ? atof(s) : 0; xc0 = (c && c[0]) ? atof(c) : 3; cached = 1; }
	*xc = xc0; return om;
}
static void hill_coriolis(postype vx, postype vy, postype *ax, postype *ay){
	postype xc, om = hill_omega(&xc);
	if(om <= 0) return;
	*ax += 2*om*vy;
	*ay -= 2*om*vx;
}
static void kepler_accel(postype x, postype y, postype *ax, postype *ay){
	{ postype xc, om = hill_omega(&xc); if(om > 0) *ax += 3*om*om*(x - xc); }
	const char *e = getenv("EUNHA_IC");
	if(!e || strcmp(e, "kepler") != 0) return;
	const char *cxe = getenv("EUNHA_KEPLER_CX");
	const char *cye = getenv("EUNHA_KEPLER_CY");
	const char *ee = getenv("EUNHA_KEPLER_EPS");
	postype cx = cxe && cxe[0] ? atof(cxe) : (postype)2.0;
	postype cy = cye && cye[0] ? atof(cye) : (postype)2.0;
	postype eps = ee && ee[0] ? atof(ee) : (postype)0.0;
	postype dx = x - cx, dy = y - cy;
	postype r2 = dx*dx + dy*dy;
	postype f;
	if(gfs_kepler_force_gizmo()){
		/* GIZMO GravAccel_KeplerianTestProblem (gravity/analytic_gravity.h):
		 * a = -rhat/r^2 for 0.35 < r < 2.1; times (r/0.35)^2 - (0.35-r)/0.35
		 * for r <= 0.35 (repulsive for r < 0.216); times 1 + (r-2.1)/0.1 for
		 * r >= 2.1. No softening. */
		double r = sqrt((double)r2), g;
		if(!(r > 1e-12)) return;
		if(r <= 0.35) g = pow(r/0.35, 2) - (0.35 - r)/0.35;
		else if(r < 2.1) g = 1;
		else g = 1 + (r - 2.1)/0.1;
		f = (postype)(-g/(r*r*r));
		*ax += dx*f; *ay += dy*f;
		return;
	}
	if(gfs_kepler_soft_spline()){
		double r = sqrt((double)r2), rc = gfs_kepler_rc();
		if(!(r > 1e-12)) return;
		f = (postype)(-kepler_spline_mass(r/rc)/(r*r*r));
		*ax += dx*f; *ay += dy*f;
		return;
	}
	if(eps > 0){
		postype R2 = r2 + eps*eps;
		f = -1.0 / (R2 * sqrt(R2));
	}else if(r2 > (postype)1e-16){
		postype invr = 1.0 / sqrt(r2);
		f = -invr * invr * invr;
	}else{
		return;
	}
	*ax += dx * f;
	*ay += dy * f;
}


/* GFS_E_AUDIT=1: split the RK4 energy change of one step into
 *   hyd  = sum_k w_k sum_i die_i dt   (face fluxes; 0 if every face pair
 *          uses one pressure and one face velocity on both ends)
 *   wg   = sum_k w_k sum_i m v_i.g_i dt (gravity work seen by RK4)
 * and compare with the actual change of E_hyd and E_pot at the end of the
 * step: rk_hyd = dE_hyd - hyd - wg - floors, rk_pot = dE_pot + wg. These two
 * are the RK4 time-discretisation residuals (KE is quadratic in v; the
 * stage rates cancel pairwise but the final kinetic energy is not their
 * weighted sum). Diagnostic only. */
static int gfs_eaud_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_E_AUDIT"); on = (s && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static double gfs_eaud_hyd = 0, gfs_eaud_wg = 0, gfs_eaud_eh0 = 0, gfs_eaud_ep0 = 0;
static double gfs_eaud_hyd_abs = 0;
static postype kepler_phi(postype x, postype y);
static void gfs_eaud_stage(treevorostressrk4particletype *sbp, int np,
		int (*targetBP)(treevorork4particletype*, postype, postype),
		postype Lx, postype Ly, postype dt, double w);
static void gfs_eaud_start(treevorostressrk4particletype *sbp, int np,
		int (*targetBP)(treevorork4particletype*, postype, postype),
		postype Lx, postype Ly);

/* GFS mode (av_mode 5) and Laguerre weights from GAS kappa do not mix.
 * kappa > 0: hydro-derived weights move the faces and the HLL face flux is
 * applied on top, so the same physics is counted twice. kappa < 0: a
 * uniform nonzero w2 silently turns off the face-rotation term. Force
 * kappa = 0 (w2 = 0 exactly) unless GFS_ALLOW_KAPPA=1 (default 0, for
 * future experiments). Called after the parameters are read (eunha2.c) and
 * at the first RK4 blend step. Returns 1 when kappa was reset. */
int gfs_kappa_guard(SimParameters *simpar, const char *where){
	const char *e = getenv("GFS_ALLOW_KAPPA");
	int allow = (e && atoi(e) > 0);
	if(GAS_AVMODE(simpar) != 5 || GAS_Kappa(simpar) == 0) return 0;
	if(allow){
		if(MYID(simpar) == 0)
			fprintf(stderr, "[KAPPA_GUARD] %s: av_mode=5 with GAS kappa=%g kept (GFS_ALLOW_KAPPA=1). Experimental.\n",
				where, (double)GAS_Kappa(simpar));
		return 0;
	}
	if(MYID(simpar) == 0)
		fprintf(stderr, "[KAPPA_GUARD] WARNING %s: av_mode=5 (GFS) with GAS kappa=%g. %s Forcing kappa=0 (w2=0). Set GFS_ALLOW_KAPPA=1 to keep it.\n",
			where, (double)GAS_Kappa(simpar),
			GAS_Kappa(simpar) > 0 ? "Hydro Laguerre weights plus the HLL face flux count the same physics twice."
			                      : "A uniform w2 disables the face-rotation term.");
	GAS_Kappa(simpar) = 0;
	return 1;
}

static void kepler_center(postype *cx, postype *cy){
	const char *cxe = getenv("EUNHA_KEPLER_CX");
	const char *cye = getenv("EUNHA_KEPLER_CY");
	*cx = cxe && cxe[0] ? atof(cxe) : (postype)2.0;
	*cy = cye && cye[0] ? atof(cye) : (postype)2.0;
}

/* Potential of the same central mass, for the energy diagnostics only.
 * Plummer when EUNHA_KEPLER_EPS>0, which is the potential of the
 * softened force above. Returns 0 when EUNHA_IC is not kepler. */
static postype kepler_phi(postype x, postype y){
	{ postype xc, om = hill_omega(&xc); if(om > 0) return -(postype)1.5*om*om*(x - xc)*(x - xc); }
	const char *e = getenv("EUNHA_IC");
	if(!e || strcmp(e, "kepler") != 0) return 0;
	const char *cxe = getenv("EUNHA_KEPLER_CX");
	const char *cye = getenv("EUNHA_KEPLER_CY");
	const char *ee = getenv("EUNHA_KEPLER_EPS");
	postype cx = cxe && cxe[0] ? atof(cxe) : (postype)2.0;
	postype cy = cye && cye[0] ? atof(cye) : (postype)2.0;
	postype eps = ee && ee[0] ? atof(ee) : (postype)0.0;
	postype dx = x - cx, dy = y - cy;
	if(gfs_kepler_force_gizmo()){
		/* potential of the GIZMO test-problem force, continuous at 0.35, 2.1:
		 * -1/r (middle); r/R0^2 + 1/r + ln(r)/R0 + C_in (r <= R0 = 0.35);
		 * 10 ln r + 20/r + C_out (r >= 2.1) */
		double r = sqrt((double)(dx*dx + dy*dy)), R0 = 0.35, R1 = 2.1;
		if(!(r > 1e-12)) r = 1e-12;
		if(r <= R0) return (postype)(r/(R0*R0) + 1/r + log(r)/R0 - 3/R0 - log(R0)/R0);
		if(r < R1) return (postype)(-1/r);
		return (postype)(10*log(r) + 20/r - 21/R1 - 10*log(R1));
	}
	if(gfs_kepler_soft_spline()){
		/* potential of the cubic-spline softened point mass (GADGET-2 W2 kernel), continuous with -1/r at r_c */
		double r = sqrt((double)(dx*dx + dy*dy)), rc = gfs_kepler_rc(), u = r/rc;
		if(u >= 1) return (postype)(-1.0/r);
		if(u < 0.5) return (postype)((16.0/3*u*u - 48.0/5*u*u*u*u + 32.0/5*u*u*u*u*u - 14.0/5)/rc);
		return (postype)((1.0/(15*u) + 32.0/3*u*u - 16*u*u*u + 48.0/5*u*u*u*u - 32.0/15*u*u*u*u*u - 16.0/5)/rc);
	}
	postype r2 = dx*dx + dy*dy + eps*eps;
	if(!(r2 > (postype)1e-32)) return 0;
	return -(postype)1.0 / sqrt(r2);
}

static void gfs_eaud_stage(treevorostressrk4particletype *sbp, int np,
		int (*targetBP)(treevorork4particletype*, postype, postype),
		postype Lx, postype Ly, postype dt, double w){
	int i;
	double h = 0, ha = 0, g = 0;
	for(i=0;i<np;i++){
		if(!targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)) continue;
		postype gx = 0, gy = 0;
		kepler_accel(sbp[i].x, sbp[i].y, &gx, &gy);
		h  += (double)sbp[i].die * (double)dt;
		ha += fabs((double)sbp[i].die * (double)dt);
		g  += (double)sbp[i].mass * ((double)sbp[i].vx*gx + (double)sbp[i].vy*gy) * (double)dt;
	}
	gfs_eaud_hyd += w*h; gfs_eaud_hyd_abs += w*ha; gfs_eaud_wg += w*g;
}
static void gfs_eaud_start(treevorostressrk4particletype *sbp, int np,
		int (*targetBP)(treevorork4particletype*, postype, postype),
		postype Lx, postype Ly){
	int i;
	double eh = 0, ep = 0;
	for(i=0;i<np;i++){
		if(!targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)) continue;
		double m = (double)sbp[i].mass;
		eh += (double)sbp[i].ie + 0.5*m*((double)sbp[i].vx*sbp[i].vx + (double)sbp[i].vy*sbp[i].vy);
		ep += m*(double)kepler_phi(sbp[i].x, sbp[i].y);
	}
	gfs_eaud_eh0 = eh; gfs_eaud_ep0 = ep;
	gfs_eaud_hyd = gfs_eaud_wg = gfs_eaud_hyd_abs = 0;
}

/* Orbital time step for the central mass. EUNHA_KEPLER_DTETA = eta > 0
 * limits dt to eta * (R^2+eps^2)^(3/4), i.e. eta/(2 pi) of the local
 * orbital period. Default 0 leaves dt alone. Returns a global minimum. */
static postype kepler_dt_limit(SimParameters *simpar, void *base, size_t stride, int np){
	static int cached = 0;
	static postype eta = 0;
	if(!cached){
		const char *s = getenv("EUNHA_KEPLER_DTETA");
		eta = (s && s[0]) ? atof(s) : (postype)0;
		cached = 1;
	}
	const char *e = getenv("EUNHA_IC");
	postype dtl = (postype)1e30, dtg;
	if(eta > 0 && e && strcmp(e, "kepler") == 0){
		const char *cxe = getenv("EUNHA_KEPLER_CX");
		const char *cye = getenv("EUNHA_KEPLER_CY");
		const char *ee = getenv("EUNHA_KEPLER_EPS");
		postype cx = cxe && cxe[0] ? atof(cxe) : (postype)2.0;
		postype cy = cye && cye[0] ? atof(cye) : (postype)2.0;
		postype eps = ee && ee[0] ? atof(ee) : (postype)0.0;
		int i;
		for(i=0;i<np;i++){
			treevorork4particletype *b = (treevorork4particletype*)((char*)base + (size_t)i*stride);
			postype dx = b->x - cx, dy = b->y - cy;
			postype r2 = dx*dx + dy*dy + eps*eps;
			postype d = eta * pow((double)r2, 0.75);
			if(d < dtl) dtl = d;
		}
	}
	MPI_Allreduce(&dtl, &dtg, 1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
	return dtg;
}

/* Voronoi bisector normal speed at the face centroid.
 * w = 0: u_n = n·ū − (u_q−u_p)·(c_f−m)/d. The second term is
 * zero when the centroid sits on the midpoint.
 * (cx, cy) is the face centroid relative to particle i (corner coordinates
 * are local to the generator, see Voro2D_FindVC), so m = (x_j - x_i)/2.
 * Laguerre faces (either w2 nonzero): only with laguerre = 1
 * (GFS_LAGUERRE_ROTATION). The reference point is then the anchor
 * a - x_i = fact1 (x_j - x_i), fact1 = (1 + (w2i - w2j)/d^2)/2, the point
 * whose velocity get2dUpqradRk4 returns. The w = 0 branch is unchanged. */
static void voronoi_face_rotation(Voro2D_point *urad,
		postype cx, postype cy,
		postype xi, postype yi, postype xj, postype yj,
		postype vxi, postype vyi, postype vxj, postype vyj,
		postype erx, postype ery, postype dist,
		postype w2i, postype w2j, int laguerre)
{
	postype dx, dy, delta;
	if(!(dist > (postype)0)) return;
	if(w2i == (postype)0 && w2j == (postype)0){
		dx = cx - (postype)0.5 * (xj - xi);
		dy = cy - (postype)0.5 * (yj - yi);
	} else {
		postype rx = xj - xi, ry = yj - yi;
		postype fact1;
		if(!laguerre) return;
		/* Same d^2 and fact1 as get2dUpqradRk4. */
		fact1 = (postype)0.5 * ((postype)1 + (w2i - w2j) / (rx*rx + ry*ry));
		dx = cx - fact1 * rx;
		dy = cy - fact1 * ry;
	}
	delta = -((vxj - vxi) * dx + (vyj - vyi) * dy) / dist;
	urad->x += delta * erx;
	urad->y += delta * ery;
}

/* get2dUpqradRk4 plus, under GFS_W2_LINEAR, the exact weight term
 * (Wdot_p - Wdot_q)/(2d) along e (Wdot in w2hydro). */
static inline Voro2D_point gfs_upqrad_rk4(treevorork4particletype *p,
		treevorork4particletype *q, postype dtold){
	Voro2D_point u = get2dUpqradRk4(p, q, dtold);
	if(gfs_dynw2_active){
		/* Wall mirror (MAX_INDEX): use the cell's own rate. */
		double wp = dyn_wdot(p);
		double wq = (PINDX(q) == MAX_INDEX) ? wp : dyn_wdot(q);
		postype ex = q->x - p->x, ey = q->y - p->y;
		postype d2 = ex*ex + ey*ey;
		if(d2 > 0){
			postype f = (postype)((wp - wq) / (2*d2));
			u.x += f*ex; u.y += f*ey;
		}
	} else if(gfs_w2lin_active){
		postype ex = q->x - p->x, ey = q->y - p->y;
		postype d2 = ex*ex + ey*ey;
		if(d2 > 0){
			postype f = (p->w2hydro - q->w2hydro) / (2*d2);
			u.x += f*ex; u.y += f*ey;
		}
	}
	return u;
}

/* SEDOV_PHASE1 stores dE/dt in die. RK4 must not add that to ie. */
static postype phase1_ie_stage(treevorostressrk4particletype *p, postype dt)
{
	postype dke;
	if(!sedov_phase1_on()) return p->die * dt;
	dke = p->mass * (p->vx * p->ax + p->vy * p->ay);
	return (p->die - dke) * dt;
}

static inline void phase1_half_kick(treevorostressrk4particletype *p,
		postype half_dt, postype ax_ext, postype ay_ext){
	postype m = p->mass;
	/* Hydro acceleration only. External gravity is applied after the
	 * energy fix so a potential, which is not in E, cannot zero ie. */
	postype ax = p->ax, ay = p->ay;
	postype vx0 = p->vx, vy0 = p->vy;
	postype ke0 = (postype)0.5 * m * (vx0*vx0 + vy0*vy0);
	postype E = p->ie + ke0 + p->die * half_dt;
	p->vx = vx0 + ax * half_dt;
	p->vy = vy0 + ay * half_dt;
	postype ke1 = (postype)0.5 * m * (p->vx*p->vx + p->vy*p->vy);
	postype ie1 = E - ke1;
	if(ie1 > (postype)0){
		p->ie = ie1;
	} else {
		/* (1/2) m |a Δt|^2 is not carried by the face flux. Subtracting
		 * the whole new kinetic energy therefore wipes the PdV heat and
		 * the cold shock never builds pressure. Keep the first-order
		 * internal energy and shed only the surplus kinetic energy. */
		postype work = m * (vx0*ax + vy0*ay) * half_dt;
		postype ie_pdv = p->ie + p->die * half_dt - work;
		if(ie_pdv < (postype)0) ie_pdv = (postype)0;
		postype ie_keep = (E > ie_pdv) ? ie_pdv : ((E > (postype)0) ? E : (postype)0);
		postype ke_allow = E - ie_keep;
		if(ke_allow < (postype)0) ke_allow = (postype)0;
		if(ie_keep <= (postype)0) phase1_ie_clips++;
		phase1_clip_events++;
		{
			postype vxb = p->vx, vyb = p->vy;
			if(ke1 > (postype)0 && ke_allow < ke1){
				postype sc = sqrt(ke_allow / ke1);
				phase1_clip_dke += (double)(ke1 - ke_allow);
				postype kcx, kcy;
				kepler_center(&kcx, &kcy);
				phase1_clip_dL += (double)(m * (sc - (postype)1)
						* ((p->x - kcx)*vyb - (p->y - kcy)*vxb));
			}
			if(E < (postype)0) phase1_clip_inj += (double)(-E);
		}
		if(ke1 > (postype)0 && ke_allow < ke1){
			postype s = sqrt(ke_allow / ke1);
			p->vx *= s;
			p->vy *= s;
		} else if(!(ke_allow > (postype)0)){
			p->vx = 0;
			p->vy = 0;
		}
		p->ie = ie_keep;
	}
	p->vx += ax_ext * half_dt;
	p->vy += ay_ext * half_dt;
}

/* Fix C v2 (2026-05-12): face-level near-vacuum guard with floor discriminator.
 *   v1 (engage on any ratio>100) also caught shock-front faces where pmin is at
 *   the ambient pressure floor (P_amb~1e-5). At a real shock front, pmin sits
 *   near the floor while pmax is post-shock — ratio is millions but the face
 *   is physical, and HLLC must drive shock propagation. Engaging Fix C there
 *   collapsed pstar to floor and killed outward force → shock front instability.
 *
 *   v2 adds pmin > HLLC_VACUUM_PMIN_GUARD: only engage when BOTH sides are
 *   well above floor (i.e., interior post-shock cells, one rarefaction-cooled).
 *   The Sedov Mode 2 pathological case had pmin=0.74 — well above 0.01. */
#ifndef HLLC_VACUUM_PRATIO
#define HLLC_VACUUM_PRATIO 100.0
#endif
#ifndef HLLC_VACUUM_PMIN_GUARD
#define HLLC_VACUUM_PMIN_GUARD 0.01
#endif

/* Acceleration CFL (Mode 2 mitigation):
 * Constrain dt < ACC_CFL_FRAC · sqrt(h_i / |a_i|), so cell displacement
 * under current acceleration stays a small fraction of cell size.  Catches
 * Sedov shock-front overshoot where ∇P/ρ jumps 4+ orders in one step
 * (acceleration explosion that K-based CFL alone cannot anticipate). */
#ifndef ACC_CFL_FRAC
#define ACC_CFL_FRAC ((postype)0.25)
#endif

/* dK-aware CFL (Mode 2 mitigation):
 * Shock-front cells can develop |dK| ~ 1e+9 spikes in Sedov 2D, which
 * overshoot the cell's K in a single substep. Rather than clipping dK
 * (breaks energy/entropy conservation), constrain dt so that the relative
 * change |dK|·dt / K stays below DKCFL_FRAC.  DKCFL_KFLOOR prevents
 * degeneracy for cold cells (K≈0) where dt would otherwise → 0. */
#ifndef DKCFL_FRAC
#define DKCFL_FRAC ((postype)0.1)
#endif
#ifndef DKCFL_KFLOOR
#define DKCFL_KFLOOR ((postype)1.0e-10)
#endif

static inline void hllc_face_2d(
		postype rhoL, postype pL, postype vnL, postype cL,
		postype rhoR, postype pR, postype vnR, postype cR,
		postype Gamma,
		postype *pstar, postype *vnstar)
{
	const postype tiny = (postype)1.0e-30;
	postype S_L, S_R, S_M, P_M;
	postype cmax = cL > cR ? cL : cR;

	/* Phase 1: HLL average only when the two pressures differ by the extreme
	 * ratio. A uniformly cold face, as in the Noh inflow, stays on HLLC. */
	if(sedov_phase1_on()){
		postype pmin = pL < pR ? pL : pR;
		postype pmax = pL > pR ? pL : pR;
		if(pmax > (postype)100.0 * fmax(pmin, (postype)1.0e-30)){
			hll_star_state(rhoL, pL, vnL, cL, rhoR, pR, vnR, cR, Gamma, pstar, vnstar);
			return;
		}
	}

	/* Fix C v2: extreme P-ratio guard (HLLE fallback), only for interior cells */
	{
		postype pmin = pL < pR ? pL : pR;
		postype pmax = pL > pR ? pL : pR;
		if(pmin > HLLC_VACUUM_PMIN_GUARD && pmax > HLLC_VACUUM_PRATIO * pmin){
			*pstar  = pmin;
			*vnstar = 0.5*(vnL + vnR);
			return;
		}
	}

	/* Vacuum check: large rarefaction → no valid Riemann solution */
	if((vnR - vnL) > cmax){
		*pstar = tiny; *vnstar = 0.5*(vnL + vnR);
		return;
	}

	/* --- Primary: Gaburov 'simplest' HLLC (Davis wave speeds) --- */
	{
		postype vmin = vnL < vnR ? vnL : vnR;
		postype vmax = vnL > vnR ? vnL : vnR;
		S_L = vmin - cmax;
		S_R = vmax + cmax;
		postype rho_wt_L = rhoL * (S_L - vnL);   /* ≤ 0 */
		postype rho_wt_R = rhoR * (S_R - vnR);   /* ≥ 0 */
		postype denom = rho_wt_L - rho_wt_R;     /* ≤ 0 */
		if(fabs(denom) > tiny){
			S_M = ((pR - pL) + rho_wt_L*vnL - rho_wt_R*vnR) / denom;
			P_M = (pL*rho_wt_R - pR*rho_wt_L + rho_wt_L*rho_wt_R*(vnR - vnL))
			      / (rho_wt_R - rho_wt_L);
			if(P_M > 0 && !isnan(P_M) && !isinf(P_M)){
				*pstar = P_M; *vnstar = S_M; return;
			}
		}
	}

	/* --- Fallback 1: Roe-averaged wave speeds + HLLC --- */
	{
		postype sqL = sqrt(rhoL), sqR = sqrt(rhoR);
		postype sq_inv = 1.0/(sqL + sqR + tiny);
		postype vn_roe = (sqL*vnL + sqR*vnR) * sq_inv;
		/* Enthalpy per unit mass (normal-direction only — 1D Roe slice):
		 *   h = γp/((γ-1)ρ) + 0.5·vₙ² */
		postype h_L = Gamma*pL/((Gamma - 1.0)*rhoL) + 0.5*vnL*vnL;
		postype h_R = Gamma*pR/((Gamma - 1.0)*rhoR) + 0.5*vnR*vnR;
		postype h_roe = (sqL*h_L + sqR*h_R) * sq_inv;
		postype c2_roe = (Gamma - 1.0) * (h_roe - 0.5*vn_roe*vn_roe);
		if(c2_roe < tiny) c2_roe = tiny;
		postype c_roe = sqrt(c2_roe);

		postype srA = vnR + cR, srB = vn_roe + c_roe;
		postype slA = vnL - cL, slB = vn_roe - c_roe;
		S_R = srA > srB ? srA : srB;
		S_L = slA < slB ? slA : slB;

		postype rho_wt_R =  rhoR * (S_R - vnR);   /* ≥ 0 */
		postype rho_wt_L = -rhoL * (S_L - vnL);   /* ≥ 0 (sign flip) */
		postype denom = rho_wt_R + rho_wt_L;
		if(fabs(denom) > tiny){
			S_M = (rho_wt_R*vnR + rho_wt_L*vnL + (pL - pR)) / denom;
			P_M = rhoL*(vnL - S_L)*(vnL - S_M) + pL;
			if(P_M > 0 && !isnan(P_M) && !isinf(P_M)){
				*pstar = P_M; *vnstar = S_M; return;
			}
		}
	}

	/* --- Last resort: Primitive-variable + Rusanov-like bounds --- */
	{
		P_M = 0.5*((pL + pR) + (vnL - vnR)*0.25*(rhoL + rhoR)*(cL + cR));
		S_M = 0.5*(vnR + vnL) + 2.0*(pL - pR)/((rhoL + rhoR)*(cL + cR) + tiny);
		postype a1 = fabs(vnL - cL), a2 = fabs(vnR - cR);
		postype a3 = fabs(vnL + cL), a4 = fabs(vnR + cR);
		postype S_plus = a1;
		if(a2 > S_plus) S_plus = a2;
		if(a3 > S_plus) S_plus = a3;
		if(a4 > S_plus) S_plus = a4;
		if(S_M < -S_plus) S_M = -S_plus;
		if(S_M >  S_plus) S_M =  S_plus;
		if(P_M < tiny)    P_M = tiny;
		*pstar = P_M; *vnstar = S_M;
	}
}

/* ================================================================
 *  hllc_face_2d_rest_frame: HLLC in the face rest frame (Arepo-style)
 *
 *  Boost the normal velocities into the face frame w_n (face moves
 *  with velocity w in the lab frame), solve HLLC there, then unboost.
 *  Pressure is Galilean-invariant, so p_star is returned unchanged.
 *  The normal star velocity is returned in the lab frame.
 *
 *  Key property (Springel 2010 §32): for Lagrangian particles where
 *  w_ij = (v_i + v_j)/2 is close to both v_i and v_j, the face-frame
 *  relative velocities are small, so HLLC dissipation is small and
 *  contact discontinuities do not suffer from the "fake compression"
 *  artifact seen with lab-frame Godunov on moving meshes.
 * ================================================================ */
static inline void hllc_face_2d_rest_frame(
		postype rhoL, postype pL, postype vnL_lab, postype cL,
		postype rhoR, postype pR, postype vnR_lab, postype cR,
		postype wn,
		postype Gamma,
		postype *pstar, postype *vnstar_lab)
{
	postype vnL = vnL_lab - wn;
	postype vnR = vnR_lab - wn;
	postype pst, vnst;
	hllc_face_2d(rhoL, pL, vnL, cL,
	             rhoR, pR, vnR, cR,
	             Gamma, &pst, &vnst);
	*pstar       = pst;        /* pressure is frame-invariant */
	*vnstar_lab  = vnst + wn;  /* boost back to lab frame */
}

/* TCF face pressure P*(rho, p, u on both sides): Davis-speed HLLC with the
 * reconstructed sound speeds, the exact vacuum condition and no pressure-ratio
 * or fallback branches, so that P* is one function of du = uL - uR. Negative
 * HLLC values (strong rarefaction) are set to 0, which keeps P* monotone. */
static double tcf_pstar(double rhoL, double pL, double uL,
		double rhoR, double pR, double uR, double Gamma)
{
	double cL = sqrt(Gamma*pL/rhoL), cR = sqrt(Gamma*pR/rhoR);
	if(uR - uL >= 2*(cL + cR)/(Gamma - 1)) return 0;
	double cmax = cL > cR ? cL : cR;
	double SL = (uL < uR ? uL : uR) - cmax, SR = (uL > uR ? uL : uR) + cmax;
	double aL = rhoL*(SL - uL), aR = rhoR*(SR - uR);
	double P = (pL*aR - pR*aL + aL*aR*(uR - uL))/(aR - aL);
	return P > 0 ? P : 0;
}

/* ================================================================
 *  update_alpha_cd_2d: Cullen-Dehnen (2010) viscosity switch
 *  Called once per full RK4 timestep after final state combination.
 * ================================================================ */
static void update_alpha_cd_2d(SimParameters *simpar, postype dt){
	/* Stride-aware access: honor TVORORK4_DDINFO n_size in case future
	   extensions append fields beyond Stress. */
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype cd_amax = GAS_CDAMAX(simpar);
	postype cd_ell  = GAS_CDELL(simpar);
	postype cd_amin = GAS_CDAMIN(simpar);

	/* Revision 4 defensive guard: skip if CD10 is disabled or dt is
	 * invalid (happens on first integrator call before dt estimation
	 * or if parameters are corrupt). Prevents NaN from divisions. */
	if(cd_amax <= 0 || cd_ell <= 0 || dt <= 0) return;

	int i;
	postype d_amin = 1e30, d_amax = -1e30, d_asum = 0.0;
	postype d_Amax = 0.0, d_vsmax = 0.0;
	int d_vs_idx = -1, d_A_idx = -1;
	long d_nfloor = 0, d_ncap = 0, d_nact = 0;
	postype d_Vmin = 1e30, d_rhomax = 0.0, d_w2min = 1e30, d_csmax = 0.0;
	int d_V_idx = -1, d_rho_idx = -1, d_w2_idx = -1, d_cs_idx = -1;
	postype d_vmax = 0.0, d_accmax = 0.0;  /* SHOCK_DIAG_VMAX / SHOCK_DIAG_AMAX tracking */
	int d_v_idx = -1, d_acc_idx = -1;
	for(i=0;i<nbp;i++){
		treevorostressrk4particletype *bpi =
			(treevorostressrk4particletype*)(bp_raw + i*p_size);
		postype h = sqrt(bpi->volume);
		postype A = fmax(0.0, -(bpi->stress.divv - bpi->stress.divv_old) / dt);
		postype vs = bpi->stress.vsig_max;
		postype alpha_loc = cd_amax * h*h * A / (vs*vs + h*h * A + 1e-30);
		postype tau = h / (2.0 * cd_ell * vs + 1e-30);
		if(alpha_loc > bpi->stress.alpha_cd)
			bpi->stress.alpha_cd = alpha_loc;
		else
			bpi->stress.alpha_cd = alpha_loc + (bpi->stress.alpha_cd - alpha_loc) * exp(-dt/tau);
		if(bpi->stress.alpha_cd < cd_amin) bpi->stress.alpha_cd = cd_amin;
		if(bpi->stress.alpha_cd > cd_amax) bpi->stress.alpha_cd = cd_amax;
		postype a = bpi->stress.alpha_cd;
		if(a < d_amin) d_amin = a;
		if(a > d_amax) d_amax = a;
		d_asum += a;
		if(a <= cd_amin + 1e-9*(cd_amax+1e-30)) d_nfloor++;
		else if(a >= cd_amax - 1e-9*(cd_amax+1e-30)) d_ncap++;
		else d_nact++;
		if(A > d_Amax) { d_Amax = A; d_A_idx = i; }
		if(vs > d_vsmax) { d_vsmax = vs; d_vs_idx = i; }
		if(bpi->volume < d_Vmin) { d_Vmin = bpi->volume; d_V_idx = i; }
		if(bpi->den    > d_rhomax){ d_rhomax = bpi->den; d_rho_idx = i; }
		if(bpi->w2     < d_w2min){ d_w2min = bpi->w2; d_w2_idx = i; }
		if(bpi->csound > d_csmax){ d_csmax = bpi->csound; d_cs_idx = i; }
		/* Velocity / acceleration runaway diagnostics */
		{
			postype vsq = (postype)bpi->vx*(postype)bpi->vx + (postype)bpi->vy*(postype)bpi->vy;
			postype accsq = (postype)bpi->ax*(postype)bpi->ax + (postype)bpi->ay*(postype)bpi->ay;
			if(vsq > d_vmax){ d_vmax = vsq; d_v_idx = i; }
			if(accsq > d_accmax){ d_accmax = accsq; d_acc_idx = i; }
		}
		bpi->stress.divv_old = bpi->stress.divv;
		bpi->stress.vsig_max = 0;
	}
	d_vmax = sqrt(d_vmax);
	d_accmax = sqrt(d_accmax);
	{
		postype g_amin, g_amax, g_asum, g_Amax, g_vsmax;
		long g_nfloor, g_ncap, g_nact, g_nbp = nbp;
		MPI_Allreduce(&d_amin, &g_amin, 1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
		MPI_Allreduce(&d_amax, &g_amax, 1, MPI_POSTYPE, MPI_MAX, MPI_COMM(simpar));
		MPI_Allreduce(&d_asum, &g_asum, 1, MPI_POSTYPE, MPI_SUM, MPI_COMM(simpar));
		MPI_Allreduce(&d_Amax, &g_Amax, 1, MPI_POSTYPE, MPI_MAX, MPI_COMM(simpar));
		MPI_Allreduce(&d_vsmax, &g_vsmax, 1, MPI_POSTYPE, MPI_MAX, MPI_COMM(simpar));
		MPI_Allreduce(&d_nfloor, &g_nfloor, 1, MPI_LONG, MPI_SUM, MPI_COMM(simpar));
		MPI_Allreduce(&d_ncap, &g_ncap, 1, MPI_LONG, MPI_SUM, MPI_COMM(simpar));
		MPI_Allreduce(&d_nact, &g_nact, 1, MPI_LONG, MPI_SUM, MPI_COMM(simpar));
		MPI_Allreduce(MPI_IN_PLACE, &g_nbp, 1, MPI_LONG, MPI_SUM, MPI_COMM(simpar));
		if(MYID(simpar)==0)
			fprintf(stderr,"[ALPHA_CD] dt=%.4e mean=%.4e min=%.4e max=%.4e Amax=%.4e vsmax=%.4e floor=%ld cap=%ld act=%ld N=%ld\n",
				(double)dt,(double)(g_asum/(g_nbp+1e-30)),(double)g_amin,(double)g_amax,
				(double)g_Amax,(double)g_vsmax,g_nfloor,g_ncap,g_nact,g_nbp);
	}
	/* Identify particles achieving global max vsig and max A.
	 * Print their full state to localize the dt-collapse pathology. */
	{
		struct { double val; int rnk; } loc, glo;
		loc.val = (double)d_vsmax; loc.rnk = MYID(simpar);
		MPI_Allreduce(&loc, &glo, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM(simpar));
		if(MYID(simpar) == glo.rnk && d_vs_idx >= 0) {
			treevorostressrk4particletype *bpw =
				(treevorostressrk4particletype*)(bp_raw + d_vs_idx*p_size);
			fprintf(stderr,"[ALPHA_CD_TOPV] rnk=%d idx=%d x=%.6e y=%.6e "
				"den=%.4e P=%.4e V=%.4e cs=%.4e m=%.4e divv=%.4e divvo=%.4e a=%.4e\n",
				MYID(simpar), d_vs_idx, (double)bpw->x, (double)bpw->y,
				(double)bpw->den, (double)bpw->pressure, (double)bpw->volume,
				(double)bpw->csound, (double)bpw->mass,
				(double)bpw->stress.divv, (double)bpw->stress.divv_old,
				(double)bpw->stress.alpha_cd);
		}
		loc.val = (double)d_Amax; loc.rnk = MYID(simpar);
		MPI_Allreduce(&loc, &glo, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM(simpar));
		if(MYID(simpar) == glo.rnk && d_A_idx >= 0) {
			treevorostressrk4particletype *bpw =
				(treevorostressrk4particletype*)(bp_raw + d_A_idx*p_size);
			fprintf(stderr,"[ALPHA_CD_TOPA] rnk=%d idx=%d x=%.6e y=%.6e "
				"den=%.4e P=%.4e V=%.4e cs=%.4e m=%.4e divv=%.4e divvo=%.4e a=%.4e\n",
				MYID(simpar), d_A_idx, (double)bpw->x, (double)bpw->y,
				(double)bpw->den, (double)bpw->pressure, (double)bpw->volume,
				(double)bpw->csound, (double)bpw->mass,
				(double)bpw->stress.divv, (double)bpw->stress.divv_old,
				(double)bpw->stress.alpha_cd);
		}
	}
	/* SHOCK_DIAG: track runaway candidates (V_min, rho_max, w2_min, cs_max).
	 * Used to diagnose dt-collapse mechanism in strong-shock tests. */
	{
		postype g_Vmin, g_rhomax, g_w2min, g_csmax;
		MPI_Allreduce(&d_Vmin,   &g_Vmin,   1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
		MPI_Allreduce(&d_rhomax, &g_rhomax, 1, MPI_POSTYPE, MPI_MAX, MPI_COMM(simpar));
		MPI_Allreduce(&d_w2min,  &g_w2min,  1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
		MPI_Allreduce(&d_csmax,  &g_csmax,  1, MPI_POSTYPE, MPI_MAX, MPI_COMM(simpar));
		if(MYID(simpar)==0)
			fprintf(stderr,"[SHOCK_DIAG] dt=%.4e Vmin=%.4e rho_max=%.4e w2_min=%.4e cs_max=%.4e\n",
				(double)dt,(double)g_Vmin,(double)g_rhomax,(double)g_w2min,(double)g_csmax);

		/* Print state of V_min particle (smallest cell — runaway compression) */
		struct { double val; int rnk; } loc, glo;
		loc.val = -((double)d_Vmin); loc.rnk = MYID(simpar);
		MPI_Allreduce(&loc, &glo, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM(simpar));
		if(MYID(simpar) == glo.rnk && d_V_idx >= 0) {
			treevorostressrk4particletype *bpw =
				(treevorostressrk4particletype*)(bp_raw + d_V_idx*p_size);
			fprintf(stderr,"[SHOCK_DIAG_VMIN] rnk=%d idx=%d x=%.6e y=%.6e "
				"V=%.4e den=%.4e P=%.4e cs=%.4e w2=%.4e w2ceil=%.4e divv=%.4e a_cd=%.4e\n",
				MYID(simpar), d_V_idx, (double)bpw->x, (double)bpw->y,
				(double)bpw->volume, (double)bpw->den, (double)bpw->pressure,
				(double)bpw->csound, (double)bpw->w2, (double)bpw->w2ceil,
				(double)bpw->stress.divv, (double)bpw->stress.alpha_cd);
		}
		/* Print state of rho_max particle */
		loc.val = (double)d_rhomax; loc.rnk = MYID(simpar);
		MPI_Allreduce(&loc, &glo, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM(simpar));
		if(MYID(simpar) == glo.rnk && d_rho_idx >= 0) {
			treevorostressrk4particletype *bpw =
				(treevorostressrk4particletype*)(bp_raw + d_rho_idx*p_size);
			fprintf(stderr,"[SHOCK_DIAG_RHOMAX] rnk=%d idx=%d x=%.6e y=%.6e "
				"V=%.4e den=%.4e P=%.4e cs=%.4e w2=%.4e w2ceil=%.4e divv=%.4e a_cd=%.4e\n",
				MYID(simpar), d_rho_idx, (double)bpw->x, (double)bpw->y,
				(double)bpw->volume, (double)bpw->den, (double)bpw->pressure,
				(double)bpw->csound, (double)bpw->w2, (double)bpw->w2ceil,
				(double)bpw->stress.divv, (double)bpw->stress.alpha_cd);
		}
		/* Velocity / acceleration global maxima — root-cause for vsig CFL collapse */
		postype g_vmax, g_accmax;
		MPI_Allreduce(&d_vmax,   &g_vmax,   1, MPI_POSTYPE, MPI_MAX, MPI_COMM(simpar));
		MPI_Allreduce(&d_accmax, &g_accmax, 1, MPI_POSTYPE, MPI_MAX, MPI_COMM(simpar));
		if(MYID(simpar)==0)
			fprintf(stderr,"[SHOCK_DIAG_KIN] dt=%.4e vmax=%.4e accmax=%.4e\n",
				(double)dt, (double)g_vmax, (double)g_accmax);

		/* Print state of v_max particle (highest |v| globally) */
		loc.val = (double)d_vmax; loc.rnk = MYID(simpar);
		MPI_Allreduce(&loc, &glo, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM(simpar));
		if(MYID(simpar) == glo.rnk && d_v_idx >= 0) {
			treevorostressrk4particletype *bpw =
				(treevorostressrk4particletype*)(bp_raw + d_v_idx*p_size);
			fprintf(stderr,"[SHOCK_DIAG_VMAX] rnk=%d idx=%d x=%.6e y=%.6e "
				"vx=%.4e vy=%.4e ax=%.4e ay=%.4e V=%.4e den=%.4e P=%.4e cs=%.4e m=%.4e\n",
				MYID(simpar), d_v_idx, (double)bpw->x, (double)bpw->y,
				(double)bpw->vx, (double)bpw->vy,
				(double)bpw->ax, (double)bpw->ay,
				(double)bpw->volume, (double)bpw->den, (double)bpw->pressure,
				(double)bpw->csound, (double)bpw->mass);
		}
		/* Print state of acc_max particle (highest |a|) — predicts next-step v jump */
		loc.val = (double)d_accmax; loc.rnk = MYID(simpar);
		MPI_Allreduce(&loc, &glo, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM(simpar));
		if(MYID(simpar) == glo.rnk && d_acc_idx >= 0) {
			treevorostressrk4particletype *bpw =
				(treevorostressrk4particletype*)(bp_raw + d_acc_idx*p_size);
			fprintf(stderr,"[SHOCK_DIAG_AMAX] rnk=%d idx=%d x=%.6e y=%.6e "
				"vx=%.4e vy=%.4e ax=%.4e ay=%.4e V=%.4e den=%.4e P=%.4e cs=%.4e m=%.4e\n",
				MYID(simpar), d_acc_idx, (double)bpw->x, (double)bpw->y,
				(double)bpw->vx, (double)bpw->vy,
				(double)bpw->ax, (double)bpw->ay,
				(double)bpw->volume, (double)bpw->den, (double)bpw->pressure,
				(double)bpw->csound, (double)bpw->mass);
		}
	}
}

/* ================================================================
 *  getAccVoro2DBlend_impl: Inner implementation of blended force.
 * ================================================================ */
/* Face-ends with d < d_c on the latest force evaluation. Each interior
 * face is seen by both owners, so a shared pair contributes two. */
long long gfs_pair_face_ends = 0;

/* GFS_FACE_AUDIT=1 (diagnostic, 1 rank CPU blend): record, for every real
 * face end, the energy flux -p_f u_face.dS that enters dte, and after the
 * force loop pair the two ends of each face (keyed by particle index).
 * Energy is conserved only if p_f and u_face are single-valued, i.e. the two
 * ends sum to zero. Prints [FAUD] per force call: faces paired, ends with
 * no partner, sum and sum|.| of the mismatch, split by face type. */
#define FAUD_MAX 48
static int gfs_faud_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_FACE_AUDIT"); on = (s && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static size_t faud_cap = 0;
static int *faud_n = NULL;
static long *faud_j = NULL;
static double *faud_f = NULL, *faud_p = NULL, *faud_u = NULL;
static int *faud_fl = NULL;
static void faud_prepare(treevorork4particletype *bp0, size_t stride, int nbp){
	size_t mx = 0; int i;
	for(i=0;i<nbp;i++){
		treevorork4particletype *b = (treevorork4particletype*)((char*)bp0 + (size_t)i*stride);
		size_t k = (size_t)PINDX(b);
		if(k + 1 > mx) mx = k + 1;
	}
	if(mx > faud_cap){
		faud_cap = mx + mx/4 + 16;
		faud_n  = (int*)realloc(faud_n, faud_cap*sizeof(int));
		faud_j  = (long*)realloc(faud_j, faud_cap*FAUD_MAX*sizeof(long));
		faud_f  = (double*)realloc(faud_f, faud_cap*FAUD_MAX*sizeof(double));
		faud_p  = (double*)realloc(faud_p, faud_cap*FAUD_MAX*sizeof(double));
		faud_u  = (double*)realloc(faud_u, faud_cap*FAUD_MAX*sizeof(double));
		faud_fl = (int*)realloc(faud_fl, faud_cap*FAUD_MAX*sizeof(int));
	}
	memset(faud_n, 0, faud_cap*sizeof(int));
}
/* flags: 1 riemann_vstar (extreme face), 2 phase1 cell-centred HLL,
 * 4 pair pressure on, 8 pair charge cap bound, 16 face charge cap bound */
static void faud_match(int call){
	size_t s; int k, m;
	long npair = 0, nmiss = 0, nover = 0;
	double sum = 0, sabs = 0, msum = 0, dmax = 0, fabs_tot = 0;
	double sabs_fl[5] = {0,0,0,0,0}, sabs_plain = 0;
	double pmx_i = 0, pmx_j = 0, umx_i = 0, umx_j = 0; int flmx = 0; long imx = -1, jmx = -1;
	for(s=0;s<faud_cap;s++){
		if(faud_n[s] > FAUD_MAX){ nover++; }
		int ns = faud_n[s] < FAUD_MAX ? faud_n[s] : FAUD_MAX;
		for(k=0;k<ns;k++){
			size_t a = s*FAUD_MAX + k;
			long j = faud_j[a];
			fabs_tot += fabs(faud_f[a]);
			if(j < 0 || (size_t)j >= faud_cap){ nmiss++; msum += faud_f[a]; continue; }
			int nj = faud_n[j] < FAUD_MAX ? faud_n[j] : FAUD_MAX, found = -1;
			for(m=0;m<nj;m++) if(faud_j[(size_t)j*FAUD_MAX + m] == (long)s){ found = m; break; }
			if(found < 0){ nmiss++; msum += faud_f[a]; continue; }
			if((long)s > j) continue;
			size_t b = (size_t)j*FAUD_MAX + found;
			double d = faud_f[a] + faud_f[b];
			int fl = faud_fl[a] | faud_fl[b];
			npair++; sum += d; sabs += fabs(d);
			int q, any = 0;
			for(q=0;q<5;q++) if(fl & (1<<q)){ sabs_fl[q] += fabs(d); any = 1; }
			if(!any) sabs_plain += fabs(d);
			if(fabs(d) > dmax){ dmax = fabs(d); pmx_i = faud_p[a]; pmx_j = faud_p[b];
				umx_i = faud_u[a]; umx_j = faud_u[b]; flmx = fl; imx = (long)s; jmx = j; }
		}
	}
	fprintf(stderr, "[FAUD] call=%d faces=%ld sum|f|=%.3e mismatch_sum=%.3e mismatch_abs=%.3e unpaired=%ld unpaired_sum=%.3e over=%ld abs[vstar]=%.2e abs[p1]=%.2e abs[pair]=%.2e abs[paircap]=%.2e abs[facecap]=%.2e abs[plain]=%.2e max=%.3e i=%ld j=%ld fl=%d p=%.6e/%.6e uA=%.6e/%.6e\n",
		call, npair, fabs_tot, sum, sabs, nmiss, msum, nover,
		sabs_fl[0], sabs_fl[1], sabs_fl[2], sabs_fl[3], sabs_fl[4], sabs_plain,
		dmax, imx, jmx, flmx, pmx_i, pmx_j, umx_i, umx_j);
}

static double getAccVoro2DBlend_impl(SimParameters *simpar, postype xmin, postype ymin,
		postype xmax, postype ymax,
		postype OrderOfAccuracy, postype Courant, postype Gamma,
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *)){
	static int s_force_call = 0;
	s_force_call++;
	tcf_dump_switch(s_force_call);   /* serial entry, outside the threaded cell loop */
	gfs_pair_face_ends = 0;
	/* Diagnostic call thresholds.  Defaults catch NP=2 OMP=8 (call=752);
	 * env vars EUNHA_DBG_LO/HI let solo (call=464) reuse same code without
	 * spamming the larger-NP logs. */
	static int s_dbg_lo = -1, s_dbg_hi = -1;
	if(s_dbg_lo < 0){
		const char *lo = getenv("EUNHA_DBG_LO");
		const char *hi = getenv("EUNHA_DBG_HI");
		s_dbg_lo = lo ? atoi(lo) : 750;
		s_dbg_hi = hi ? atoi(hi) : 999999;
	}
	int do_iy_dbg = (s_force_call >= s_dbg_lo && s_force_call <= s_dbg_hi);
	int do_i_dbg  = (s_force_call >= s_dbg_lo+2 && s_force_call <= s_dbg_hi);
	if(MYID(simpar)<=1){ fprintf(stderr,"[FORCE_DBG] r%d gamma impl entry call=%d\n", MYID(simpar), s_force_call); fflush(stderr); }
	postype boxsize = VORO_CELLBOX(simpar);
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype Dtime = 1.e10;
	int nx = NX(simpar);
	int av_mode = GAS_AVMODE(simpar);
	int use_muscl = GAS_USEMUSCL(simpar);
	postype cd_amax = GAS_CDAMAX(simpar);
	postype blend_theta = GAS_BLENDTHETA(simpar);
	postype nu_phys = GAS_VISCOSITY(simpar);
	postype prandtl = GAS_PRANDTL(simpar);

	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	int mx = BASICCELL_MX(simpar);
	int my = BASICCELL_MY(simpar);

	postype dtold = GAS_dtold(simpar);
	/* Read before the OpenMP region. */
	int lag_rot = gfs_laguerre_rotation_on() || gfs_dynw2_active;
	int vaud_on_l = gfs_vaud_on();
	if(vaud_on_l){
		size_t max_id = 0;
		for(int ii=0; ii<nbp; ii++){
			size_t id = (size_t)PINDX(bp+ii);
			if(id > max_id) max_id = id;
		}
		vaud_ensure(max_id);
		for(int ii=0; ii<nbp; ii++) vaud_rate[(size_t)PINDX(bp+ii)] = 0;
	}
	if(sedov_lagvol_on() && sedov_phase1_on() && use_muscl){
		size_t max_id = 0;
		for(int ii=0; ii<nbp; ii++){
			size_t id = (size_t)PINDX(bp+ii);
			if(id > max_id) max_id = id;
		}
		if(lag_vt_at(max_id)){
			for(int ii=0; ii<nbp; ii++)
				lag_sw[(size_t)PINDX(bp+ii)] = 0;
			lag_sw_ready = 1;
		}
	}

	/* Monaghan AV parameters (used when av_mode==0) */
	float alphavis = GAS_AlphaVis(simpar);
	float betavis = GAS_BetaVis(simpar);
	float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
	float epsvis = GAS_EPSVIS(simpar);

	int iy;
	if(gfs_faud_on()) faud_prepare(VORORK4_TBP(simpar), TVORORK4_DDINFO(simpar)[0].n_size, VORO_NP(simpar));
	if(MYID(simpar)<=1){ fprintf(stderr,"[FORCE_DBG] r%d gamma1 before OMP loop mx=%d my=%d call=%d\n",MYID(simpar),mx,my,s_force_call); fflush(stderr); }
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
	for(iy=1;iy<my-1;iy++){
		if(do_iy_dbg){
#pragma omp critical (iy_dbg_in)
			{ fprintf(stderr,"[FORCE_DBG] r%d iy=%d in call=%d\n", MYID(simpar), iy, s_force_call); fflush(stderr); }
		}
		int mp=4000;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		int ix;
		for(ix=1;ix<mx-1;ix++){
			if(do_iy_dbg){
#pragma omp critical (ix_dbg)
				{ fprintf(stderr,"[FORCE_DBG] r%d ix=%d iy=%d call=%d\n", MYID(simpar), ix, iy, s_force_call); fflush(stderr); }
			}
			int np;
			treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
			int nneigh;
			Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
			Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			int i;
			for(i=0;i<np;i++){
				if(do_i_dbg){
#pragma omp critical (i_dbg_in)
					{ fprintf(stderr,"[FORCE_DBG] r%d ix=%d iy=%d i=%d/%d call=%d enter\n", MYID(simpar), ix, iy, i, np, s_force_call); fflush(stderr); }
				}
				treevorork4particletype *ibp_rk4 = p[i].bp;

				Voro2D_point center;
				center.x = p[i].x;
				center.y = p[i].y;
				center.indx = PINDX(p+i);
				center.csound = p[i].csound;
				center.w2 = p[i].w2;

				treevorostressrk4particletype *ibp = (treevorostressrk4particletype*)ibp_rk4;
				ibp_rk4->dt = 1.e10;

				/* Universal pre-FindVC sentinel dump (rate-limited).
				 * Captures exact center+pressure/den state for every entered
				 * cell at call=752, so we know what was about to be tessellated
				 * if FindVC SIGSEGVs. */
				if(do_i_dbg){
#pragma omp critical (i_dbg_dump)
					{
						static int s_dump_n = 0;
						s_dump_n++;
						if(s_dump_n <= 5000){
							fprintf(stderr,"[DUMP] n=%d r%d ix=%d iy=%d i=%d call=%d c=(%g,%g) indx=%lld w2=%g cs=%g P=%g den=%g K=%g nn=%d\n",
								s_dump_n, MYID(simpar), ix, iy, i, s_force_call,
								(double)center.x, (double)center.y, (long long)center.indx,
								(double)center.w2, (double)center.csound,
								(double)ibp->pressure, (double)ibp->den,
								(double)ibp->stress.K, nneigh);
							fflush(stderr);
						}
					}
				}

				int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
				/* === Adaptive expansion (Fix B for negative-upperrelated SIGSEGV) ===
				 * If FindVC left any active corner with negative upperrelated
				 * or lowerrelated (= sentinel from Voro2D_InitializeCorner that
				 * was never replaced by a real candidate during cell construction
				 * because no particle exists in that angular sector within the
				 * 3x3 stencil), retry with progressively larger radii.  Save
				 * cell-level pointers so the next i-iteration restores them. */
				Voro2D_point *neighbors_e = NULL, *neighwork_e = NULL;
				int nneigh_e = 0;
				Voro2D_point *neighbors_save = neighbors;
				Voro2D_point *neighwork_save = neighwork;
				int nneigh_save = nneigh;
				int adapt_R = 1;
				if(Voro2D_HasUnassignedCorner(vorocorner, ip) || voro_sec_short(simpar, vorocorner, ip, (double)center.x, (double)center.y, ix, iy, 1)){
					int R;
					for(R=2; R<=voro_expand_rmax(simpar); R++){
						if(neighbors_e){ free(neighbors_e); free(neighwork_e); neighbors_e=NULL; neighwork_e=NULL; }
						neighbors_e = searchCellRk4Neighbors2D_R(simpar, ix, iy, R, &nneigh_e);
						neighwork_e = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh_e);
						ip = Voro2D_FindVC(&center, neighbors_e, neighwork_e, nneigh_e, vorocorner, mp, boxsize);
						if(!Voro2D_HasUnassignedCorner(vorocorner, ip) && !voro_sec_short(simpar, vorocorner, ip, (double)center.x, (double)center.y, ix, iy, R)){
							neighbors = neighbors_e;
							neighwork = neighwork_e;
							nneigh = nneigh_e;
							adapt_R = R;
							break;
						}
					}
					if(R > voro_expand_rmax(simpar)){
						/* No generator within the bound in some direction:
						 * the cell keeps cap edges (negative neighbour id,
						 * vacuum), as in every other pass. */
						if(neighbors_e){
							neighbors = neighbors_e;
							neighwork = neighwork_e;
							nneigh = nneigh_e;
							adapt_R = R;
						}
					}
				}
				{ double vrR = voro_disc_radius(simpar); if(vrR > 0) voro2d_restrict_disc(vorocorner, mp, vrR); }
				int cell_has_arc = 0;
				{ Voro2D_Corner *tv = vorocorner; do { if(tv->upperrelated == VORO_ARC_EDGE) cell_has_arc = 1; tv = tv->upperlink; } while(tv != vorocorner); }
				if(do_i_dbg){
#pragma omp critical (i_dbg_vc)
					{ fprintf(stderr,"[FORCE_DBG] r%d ix=%d iy=%d i=%d call=%d post-FindVC ip=%d nneigh=%d R=%d\n",
						MYID(simpar), ix, iy, i, s_force_call, ip, nneigh, adapt_R); fflush(stderr); }
				}
				postype ibp_vx = ibp->vx;
				postype ibp_vy = ibp->vy;
				postype ibp_csound = ibp->csound;
				postype ibp_pressure = ibp->pressure;
				postype ibp_den = ibp->den;

				Voro2D_Corner *tmp,*tmp2;
				tmp = vorocorner;
				double die, dte, dke, fx, fy, dvb = 0;
				die = dte = dke = fx = fy = 0;
				/* Dissipative-only energy rate for entropy variable K=P/rho^gamma:
				 * accumulates ONLY viscous heating and heat conduction, NOT
				 * adiabatic pressure work (which preserves K).  See Springel &
				 * Hernquist 2002 MNRAS 333, 649. */
				double dK_diss = 0;
				double tcf_tds = 0, tcf_dV = 0;
				double tg_nx = 0, tg_ny = 0, tg_ex = 0, tg_ey = 0, tg_bis = 0, tg_en = 0, tg_amin = 1e300; int tg_nf = 0;   /* GFS_TCF_DUMP geometry */
				postype dK_pdv_factor = (postype)GAS_DKPDVFACTOR(simpar);
				/* XSPH face-area-weighted velocity offset accumulator
				 * (real-neighbor faces only; ghost mirrors excluded). */
				double xsph_vx_acc = 0, xsph_vy_acc = 0, xsph_w_acc = 0;
				/* Hyperviscosity Pass 2: biharm_v_i = (1/V_i) Σ A_ij (lap_v_j-lap_v_i)/d_ij */
				double biharm_vx_acc = 0, biharm_vy_acc = 0;
				int corner_n = 0;
				int do_face_dbg = (do_i_dbg && ix==31 && iy==32 && i==4 && s_force_call==464);
				do {
					corner_n++;
					if(do_face_dbg){
#pragma omp critical (face_dbg)
						{ fprintf(stderr,"[FACE] r%d call=%d corner=%d tmp=%p upperrel=%d upperlink=%p nneigh=%d\n",
							MYID(simpar), s_force_call, corner_n, (void*)tmp, tmp->upperrelated, (void*)tmp->upperlink, nneigh); fflush(stderr); }
					}
					if(vaud_on_l && tmp->upperrelated == VORO_ARC_EDGE && (size_t)PINDX(ibp_rk4) < vaud_cap)
						vaud_arc[(size_t)PINDX(ibp_rk4)] = 1;
					if(tmp->upperrelated >= 0){
						if(tmp->upperrelated >= nneigh){
							fprintf(stderr,"[GUARD] r%d ix=%d iy=%d i=%d call=%d upperrelated=%d nneigh=%d OOB!\n",
								MYID(simpar), ix, iy, i, s_force_call, tmp->upperrelated, nneigh);
							fflush(stderr);
							tmp = tmp->upperlink;
							continue;
						}
						if(neighwork[tmp->upperrelated].bp == NULL){
							fprintf(stderr,"[GUARD] r%d ix=%d iy=%d i=%d call=%d upperrelated=%d bp=NULL!\n",
								MYID(simpar), ix, iy, i, s_force_call, tmp->upperrelated);
							fflush(stderr);
							tmp = tmp->upperlink;
							continue;
						}
						treevorostressrk4particletype *jbp = (treevorostressrk4particletype*)(neighwork[tmp->upperrelated].bp);
						if(do_face_dbg){
#pragma omp critical (face_dbg)
							{ fprintf(stderr,"[FACE] r%d call=%d corner=%d jbp=%p pre-deref upperrel=%d\n",
								MYID(simpar), s_force_call, corner_n, (void*)jbp, tmp->upperrelated); fflush(stderr); }
						}
						int jbp_is_ghost = (PINDX((treevorork4particletype*)jbp) == MAX_INDEX);
						postype jbp_vx = jbp->vx;
						postype jbp_vy = jbp->vy;
						postype jbp_csound = jbp->csound;
						postype jbp_pressure = jbp->pressure;
						postype jbp_den = jbp->den;
						if(do_face_dbg){
#pragma omp critical (face_dbg)
							{ fprintf(stderr,"[FACE] r%d call=%d corner=%d jbp loaded ghost=%d vx=%g vy=%g cs=%g P=%g den=%g K=%g\n",
								MYID(simpar), s_force_call, corner_n, jbp_is_ghost,
								(double)jbp_vx, (double)jbp_vy, (double)jbp_csound,
								(double)jbp_pressure, (double)jbp_den, (double)jbp->stress.K); fflush(stderr); }
						}

						/* Face geometry */
						tmp2 = tmp->upperlink;
						Voro2D_point line;
						line.x = tmp2->x - tmp->x;
						line.y = tmp2->y - tmp->y;
						Voro2D_point dS = voro2D_norm(&line);
						postype facearea = Vec2DLength(tmp,tmp2);
						double sh_cx = 0.5*((double)tmp->x + (double)tmp2->x), sh_cy = 0.5*((double)tmp->y + (double)tmp2->y);
						if(gfs_tcf_on() && gfs_shared_geom_on() && PINDX((treevorork4particletype*)jbp) != MAX_INDEX){
							facearea = (postype)shf_common((long)PINDX(ibp_rk4), (long)PINDX((treevorork4particletype*)jbp),
								(double)jbp->x - (double)ibp_rk4->x, (double)jbp->y - (double)ibp_rk4->y,
								(double)tmp->x, (double)tmp->y, (double)tmp2->x, (double)tmp2->y, &sh_cx, &sh_cy);
						}
						dS.x = facearea*dS.x;
						dS.y = facearea*dS.y;

						/* XSPH accumulation: face-area weighted velocity
						 * difference. Skip ghost faces — wall mirrors hold
						 * reflected v that would pull particles off-axis. */
						if(!jbp_is_ghost){
							xsph_vx_acc += facearea * (jbp_vx - ibp_vx);
							xsph_vy_acc += facearea * (jbp_vy - ibp_vy);
							xsph_w_acc  += facearea;
						}

						/* Inter-particle direction */
						Voro2D_point dr = EunhaVec2DSub((treevorork4particletype*)jbp, ibp_rk4);
						postype dramp = sqrt(Vec2DDotP(&dr, &dr));
						Voro2D_point er;
						er.x = dr.x/dramp;
						er.y = dr.y/dramp;

						/* Hyperviscosity Pass 2: accumulate A_face/d_ij × Δlap_v
						 * Real-neighbor faces only; ghost mirrors excluded. */
						if(!jbp_is_ghost && dramp > 0){
							postype invd_bh = facearea / dramp;
							biharm_vx_acc += invd_bh * (jbp->stress.lap_vx - ibp->stress.lap_vx);
							biharm_vy_acc += invd_bh * (jbp->stress.lap_vy - ibp->stress.lap_vy);
						}

						postype pi_total;
						postype tau_dot_dS_x = 0, tau_dot_dS_y = 0;
						int phase1_extreme = 0;
						double hq_pL = 0, hq_pR = 0, hq_uL = 0, hq_uR = 0, hq_us = 0, hq_rhoL = 0, hq_rhoR = 0; int hq_ok = 0; double gq_cov = 0;
						double tcf_durec = 0, tcf_Z = 0; int tcf_ok = 0;   /* EUNHA_FACE_GAUSS2: (1/A) int (P* - Pbar)(w_n - wbar) ds */
						postype phase1_vn = 0, phase1_nx = 0, phase1_ny = 0;
						int riemann_vstar = 0;
						int faud_flags = 0;
						const int geom_fv = gfs_geom_fv_on() || gfs_dynw2_active;
						postype riemann_vn = 0, riemann_nx = 0, riemann_ny = 0;
						int hard_wall = 0;
						postype hard_nx = 0, hard_ny = 0;
						/* Fix A+B: dissipative AV pressure extracted from pi_total.
						 * Sum of (Monaghan AV + CD10 + HLLC-jump-excess) ≥ 0 when active.
						 * Used downstream to compute pair-symmetric ε_AV ≥ 0 heating. */
						postype p_av_for_heat = 0;

						if(av_mode == 0 && nu_phys <= 0){
							/* Original Monaghan AV path (no NS stress) */
							postype pi = VoroRK4_Pressure2D(
									ibp_pressure, jbp_pressure, tmp, neighwork, simpar, RT);

							Voro2D_point uij;
							uij.x = jbp_vx - ibp_vx;
							uij.y = jbp_vy - ibp_vy;
							postype rvel = Vec2DDotP(&er, &uij);
							if(rvel < 0){
								postype wcomp = sqrt(ibp_rk4->w2)+sqrt(((treevorork4particletype*)jbp)->w2);
								postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
								postype drampScale = dramp/scaleFactor;
								postype mu = rvel /* Voronoi-AV-reform: no length regularization */;
								postype meanden = 0.5*(ibp_den + jbp_den);
								postype meanCsound = 0.5*(ibp_csound + jbp_csound);
								pi = pi + (-alphavis*meanCsound*mu + betavis*mu*mu)*meanden;
							}
							pi_total = pi;

						} else if(jbp_is_ghost) {
							/* Cylinder mirror: reflective Riemann. The ghost is the
							 * mirror state, so the contact sits on the wall and p*
							 * is the wall pressure. Flat-wall ghosts stay on M(n,m);
							 * HLLC there sees a fake collision when the face normal
							 * is not the wall normal. */
							postype gnx, gny;
							if(cyl_wall_query(jbp->x, jbp->y, &gnx, &gny)){
								postype pst, vnst;
								postype vnL = ibp_vx*er.x + ibp_vy*er.y;
								postype vnR = jbp_vx*er.x + jbp_vy*er.y;
								hllc_face_2d_rest_frame(
										ibp_den, ibp_pressure, vnL, ibp_csound,
										jbp_den, jbp_pressure, vnR, jbp_csound,
										(postype)0, Gamma, &pst, &vnst);
								if(!(pst > 0) || isnan(pst))
									pst = 0.5*(ibp_pressure + jbp_pressure);
								pi_total = pst;
								hard_wall = 1;
								hard_nx = -gnx;
								hard_ny = -gny;
							} else {
							/* Ghost face in blend path: use M(n,m) pressure
							   + Monaghan AV for wall damping (same as av_mode=0).
							   Skip NS stress and HLLC — ghost stress fields are
							   not properly mirrored and HLLC sees artificial collision.
							   Guard M(n,m) stencil indices like Tier 1 path
							   (neighwork[-N] would SIGSEGV at corner faces). */
							postype pi;
							if(OrderOfAccuracy > 0 && (tmp->upperlink)->upperrelated >= 0 && tmp->lowerrelated >= 0){
								pi = VoroStressRK4_Pressure2D(
										ibp_pressure, jbp_pressure, tmp, neighwork, simpar, RT);
							} else {
								pi = 0.5*(ibp_pressure + jbp_pressure);
							}
							Voro2D_point uij;
							uij.x = jbp_vx - ibp_vx;
							uij.y = jbp_vy - ibp_vy;
							postype rvel = Vec2DDotP(&er, &uij);
							if(rvel < 0){
								postype wcomp = sqrt(ibp_rk4->w2)+sqrt(((treevorork4particletype*)jbp)->w2);
								postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
								postype drampScale = dramp/scaleFactor;
								postype mu = rvel /* Voronoi-AV-reform: no length regularization */;
								postype meanden = 0.5*(ibp_den + jbp_den);
								postype meanCsound = 0.5*(ibp_csound + jbp_csound);
								pi = pi + (-alphavis*meanCsound*mu + betavis*mu*mu)*meanden;
							}
							pi_total = pi;
							}

						} else if(av_mode == 3) {
							/* ========== DEPRECATED: Pure HLLC (lab frame) ==========
							 * WARNING (Revision 4, expert review):
							 *   This mode is mathematically guaranteed to fail on
							 *   Lagrangian Voronoi/Laguerre meshes in flows with
							 *   contact discontinuities (e.g. KH, RT). The contact
							 *   wave is a linearly degenerate characteristic, so a
							 *   Riemann solver alone carries zero entropy flux. In
							 *   an Eulerian code the upwind advection term provides
							 *   implicit numerical dissipation, but Lagrangian cells
							 *   move with the fluid and have no advective term.
							 *   Consequence: grid-scale noise at contacts grows
							 *   unbounded → dt collapse (confirmed by KH Step 4/5).
							 *   Proof: Toro 2009 §10.2 (HLLC contact wave), LeVeque
							 *   2002 §9.5 (Godunov on linearly degenerate fields),
							 *   Springel 2010 §32 (Arepo face rest-frame motivation).
							 * USE av_mode=5 INSTEAD (face-rest-frame HLLC + mesh
							 * regularization, which is the proven Arepo recipe).
							 * No M(n,m), no NS stress, no CD10 blending.
							 * HLLC Riemann solver is the sole dissipation mechanism. */
							postype ds_mag_inv = 1.0/(facearea + 1e-30);
							postype nx_hat = dS.x * ds_mag_inv;
							postype ny_hat = dS.y * ds_mag_inv;
							postype vnL = ibp_vx*nx_hat + ibp_vy*ny_hat;
							postype vnR = jbp_vx*nx_hat + jbp_vy*ny_hat;
							postype pL = ibp_pressure;
							postype pR = jbp_pressure;

							/* MUSCL reconstruction at face midpoint.
							   tmp->x,tmp->y are stored RELATIVE to generator ibp
							   (Voro2D_FindVC uses Vec2DSub). */
							postype xf_rel_i = 0.5*(tmp->x + tmp2->x);
							postype yf_rel_i = 0.5*(tmp->y + tmp2->y);
							postype dx_ij_rel = jbp->x - ibp->x;
							postype dy_ij_rel = jbp->y - ibp->y;
							postype dx_iF = xf_rel_i;
							postype dy_iF = yf_rel_i;
							postype dx_jF = xf_rel_i - dx_ij_rel;
							postype dy_jF = yf_rel_i - dy_ij_rel;

							/* Pressure reconstruction */
							postype dp_i = ibp->stress.dPdx*dx_iF + ibp->stress.dPdy*dy_iF;
							postype dp_j = jbp->stress.dPdx*dx_jF + jbp->stress.dPdy*dy_jF;
							pL += dp_i;
							pR += dp_j;

							/* Velocity reconstruction */
							postype dvnL = (ibp->stress.gUxx*nx_hat + ibp->stress.gUxy*ny_hat)*dx_iF
							             + (ibp->stress.gUyx*nx_hat + ibp->stress.gUyy*ny_hat)*dy_iF;
							postype dvnR = (jbp->stress.gUxx*nx_hat + jbp->stress.gUxy*ny_hat)*dx_jF
							             + (jbp->stress.gUyx*nx_hat + jbp->stress.gUyy*ny_hat)*dy_jF;
							vnL += dvnL;
							vnR += dvnR;

							/* Slope limiter: clamp to prevent new extrema */
							postype pmin = fmin(ibp_pressure, jbp_pressure);
							postype pmax = fmax(ibp_pressure, jbp_pressure);
							if(pL < pmin) pL = pmin; if(pL > pmax) pL = pmax;
							if(pR < pmin) pR = pmin; if(pR > pmax) pR = pmax;
							postype vnmin = fmin(vnL - dvnL, vnR - dvnR);
							postype vnmax = fmax(vnL - dvnL, vnR - dvnR);
							if(vnL < vnmin) vnL = vnmin; if(vnL > vnmax) vnL = vnmax;
							if(vnR < vnmin) vnR = vnmin; if(vnR > vnmax) vnR = vnmax;

							/* Positivity: floor on reconstructed pressure */
							if(pL < 1e-10) pL = 1e-10;
							if(pR < 1e-10) pR = 1e-10;

							postype pst, vnst;
							hllc_face_2d(ibp_den, pL, vnL, ibp_csound,
							             jbp_den, pR, vnR, jbp_csound,
							             Gamma, &pst, &vnst);
							pi_total = pst;

							/* CD10 viscous pressure (compression only) */
							if(cd_amax > 0){
								postype dvn_hllc = vnR - vnL;
								if(dvn_hllc < 0){
									postype alpha_face = 0.5*(ibp->stress.alpha_cd + jbp->stress.alpha_cd);
									postype rho_mean = 0.5*(ibp_den + jbp_den);
									postype vsig_cd = ibp_csound + jbp_csound - fmin(0.0, dvn_hllc);
									pi_total += 0.5 * alpha_face * vsig_cd * rho_mean * (-dvn_hllc);
								}
							}

							/* Optional Monaghan AV for contact noise control */
							if(alphavis > 0){
								Voro2D_point uij;
								uij.x = jbp_vx - ibp_vx;
								uij.y = jbp_vy - ibp_vy;
								postype rvel = Vec2DDotP(&er, &uij);
								if(rvel < 0){
									postype wcomp = sqrt(ibp_rk4->w2)+sqrt(((treevorork4particletype*)jbp)->w2);
									postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
									postype drampScale = dramp/scaleFactor;
									postype mu = rvel /* Voronoi-AV-reform: no length regularization */;
									postype meanden = 0.5*(ibp_den + jbp_den);
									postype meanCsound = 0.5*(ibp_csound + jbp_csound);
									pi_total += (-alphavis*meanCsound*mu + betavis*mu*mu)*meanden;
								}
							}
							/* tau_dot_dS_x/y remain 0 — no NS stress */

						} else if(av_mode == 5) {
							/* ========== AREPO-LAGUERRE: face rest frame HLLC + MUSCL ========== */
							/* Pure HLLC on every (real-real) face, solved in the
							   face rest frame.  MUSCL reconstruction of P and
							   v_n using cell-averaged gradients from
							   updateDenW2Pressure2DBlend, with a pair-wise
							   clamp limiter (prevents new extrema between i and
							   j).  Pressure floor for positivity.  No M(n,m)
							   blending, no NS stress, no explicit AV — the face
							   rest frame boost is what makes contact
							   discontinuities survive without the "fake
							   compression" artefact of lab-frame Godunov.

							   E-C1 fix: face velocity is taken from
							   get2dUpqradRk4, which is the same Laguerre-aware
							   (wfrac + AA-corrected) value used by the energy
							   path below.  This keeps momentum and internal
							   energy updates consistent with each other on
							   Laguerre meshes (w2_i ≠ w2_j). */
							postype ds_mag_inv = 1.0/(facearea + 1e-30);
							postype nx_hat = dS.x * ds_mag_inv;
							postype ny_hat = dS.y * ds_mag_inv;

							/* Face velocity in lab frame: v_face = v_i + uradix_ui
							   with uradix_ui = (v_face - v_i) from get2dUpqradRk4. */
							/* EUNHA_FACE_GAUSS2=1: two-point Gauss quadrature along the face
							 * (Zier & Springel 2022 sec. 4.2). Each point gets its own face
							 * velocity, reconstruction and HLLC state; the force uses the mean
							 * P* and the work adds the covariance of P* and the face normal
							 * velocity along the face (gq_cov). */
							int gq_n = gfs_gauss2_on() ? 2 : 1, gq_k;
							double gq_P = 0, gq_w[2] = {0, 0}, gq_p[2] = {0, 0};
							postype pst, vnst_lab;
							for(gq_k = 0; gq_k < gq_n; gq_k++){
							postype gq_px = (postype)sh_cx, gq_py = (postype)sh_cy;   /* face centroid (common segment under GFS_SHARED_GEOM) */
							if(gq_n == 2){
								postype sg = (gq_k == 0 ? -1 : 1)*(postype)0.28867513459481287;   /* 1/(2 sqrt 3) */
								gq_px += sg*(tmp2->x - tmp->x); gq_py += sg*(tmp2->y - tmp->y);
							}
							Voro2D_point uradix_ui_boost =
								gfs_upqrad_rk4(ibp_rk4,
									(treevorork4particletype*)jbp, dtold);
							voronoi_face_rotation(&uradix_ui_boost,
									gq_px,
									gq_py,
									ibp->x, ibp->y, jbp->x, jbp->y,
									ibp_vx, ibp_vy, jbp_vx, jbp_vy,
									er.x, er.y, dramp,
									ibp_rk4->w2,
									((treevorork4particletype*)jbp)->w2, lag_rot);
							postype wx = ibp_vx + uradix_ui_boost.x;
							postype wy = ibp_vy + uradix_ui_boost.y;
							postype wn = wx*nx_hat + wy*ny_hat;

							postype vnL_lab = ibp_vx*nx_hat + ibp_vy*ny_hat;
							postype vnR_lab = jbp_vx*nx_hat + jbp_vy*ny_hat;

							postype pL = ibp_pressure;
							postype pR = jbp_pressure;
							postype rhoL = ibp_den;
							postype rhoR = jbp_den;

							if(use_muscl){
								/* Springel-Arepo MUSCL: reconstruct {rho, vn, P}
								   at the face centroid using Green-Gauss
								   gradients that were already scaled by the
								   per-cell Barth-Jespersen/Springel TVD limiter
								   in updateDenW2Pressure2DBlend.  NOTE: Voronoi
								   corners tmp->x,tmp->y are stored in LOCAL
								   coordinates relative to the generator ibp
								   (see Voro2D_FindVC which calls Vec2DSub).
								   So xf is already (face - ibp), and the
								   face-relative-to-j is xf - (jbp-ibp). */
								postype xf_rel_i = gq_px;
								postype yf_rel_i = gq_py;
								postype dx_ij_rel = jbp->x - ibp->x;
								postype dy_ij_rel = jbp->y - ibp->y;
								postype dx_iF = xf_rel_i;
								postype dy_iF = yf_rel_i;
								postype dx_jF = xf_rel_i - dx_ij_rel;
								postype dy_jF = yf_rel_i - dy_ij_rel;

								postype drho_i = ibp->stress.dRhodx*dx_iF + ibp->stress.dRhody*dy_iF;
								postype drho_j = jbp->stress.dRhodx*dx_jF + jbp->stress.dRhody*dy_jF;
								postype dp_i = ibp->stress.dPdx*dx_iF + ibp->stress.dPdy*dy_iF;
								postype dp_j = jbp->stress.dPdx*dx_jF + jbp->stress.dPdy*dy_jF;
								postype dvxL = ibp->stress.gUxx*dx_iF + ibp->stress.gUxy*dy_iF;
								postype dvyL = ibp->stress.gUyx*dx_iF + ibp->stress.gUyy*dy_iF;
								postype dvnL0 = nx_hat*dvxL + ny_hat*dvyL;
								postype dvxR = jbp->stress.gUxx*dx_jF + jbp->stress.gUxy*dy_jF;
								postype dvyR = jbp->stress.gUyx*dx_jF + jbp->stress.gUyy*dy_jF;
								postype dvnR0 = nx_hat*dvxR + ny_hat*dvyR;
								if(gfs_rec_tcommon_on()){
									/* split the extrapolation at the anchor a = phi e: each side
									 * goes to a with its own gradient, a -> c with the pair mean, so
									 * the face jump equals the (limited) anchor jump and
									 * (g_j - g_i).t adds nothing; exact for linear fields */
									double ex_ = (double)dx_ij_rel, ey_ = (double)dy_ij_rel, d2_ = ex_*ex_ + ey_*ey_;
									double phi_ = 0.5 + ((double)ibp_rk4->w2 - (double)((treevorork4particletype*)jbp)->w2)/(2*d2_);
									double ax = phi_*ex_, ay = phi_*ey_;               /* anchor rel. to i */
									double tx = (double)dx_iF - ax, ty = (double)dy_iF - ay;   /* anchor -> centroid */
									double bx = ax - ex_, by = ay - ey_;               /* anchor rel. to j */
									/* GFS_TC(...): increments of side i, j at the centroid; with
									 * EUNHA_LIM_PAIR the anchor parts are clamped to the pair interval
									 * and the common tangential part is returned in tau */
									int pair_ = gfs_lim_pair_on();
#define GFS_TC(qi_, qj_, gix, giy, gjx, gjy, outI, outJ, tau) do { \
										double mx_ = 0.5*((double)(gix) + (double)(gjx)), my_ = 0.5*((double)(giy) + (double)(gjy)); \
										double aI = (double)(gix)*ax + (double)(giy)*ay, aJ = (double)(gjx)*bx + (double)(gjy)*by; \
										if(pair_){ double qa = (double)(qi_), qb = (double)(qj_), lo = fmin(qa, qb), hi = fmax(qa, qb); \
											if(qa + aI < lo) aI = lo - qa; if(qa + aI > hi) aI = hi - qa; \
											if(qb + aJ < lo) aJ = lo - qb; if(qb + aJ > hi) aJ = hi - qb; } \
										tau = mx_*tx + my_*ty; \
										outI = (postype)(aI + tau); outJ = (postype)(aJ + tau); } while(0)
									double tR_, tP_, tVx_, tVy_;
									GFS_TC(ibp_den, jbp_den, ibp->stress.dRhodx, ibp->stress.dRhody, jbp->stress.dRhodx, jbp->stress.dRhody, drho_i, drho_j, tR_);
									GFS_TC(ibp_pressure, jbp_pressure, ibp->stress.dPdx, ibp->stress.dPdy, jbp->stress.dPdx, jbp->stress.dPdy, dp_i, dp_j, tP_);
									GFS_TC(ibp_vx, jbp_vx, ibp->stress.gUxx, ibp->stress.gUxy, jbp->stress.gUxx, jbp->stress.gUxy, dvxL, dvxR, tVx_);
									GFS_TC(ibp_vy, jbp_vy, ibp->stress.gUyx, ibp->stress.gUyy, jbp->stress.gUyx, jbp->stress.gUyy, dvyL, dvyR, tVy_);
#undef GFS_TC
									if(pair_){
										/* bound both sides at the evaluation point by the values of the
										 * generators of the two Delaunay triangles sharing the face
										 * (i, j, k+, k-): linear interpolation on the face stays in this
										 * range, rho and P stay positive, and an ill-conditioned gradient
										 * cannot carry both sides to a common unphysical state */
										double lo4[4], hi4[4], qi4[4] = {(double)ibp_den, (double)ibp_pressure, (double)ibp_vx, (double)ibp_vy},
										       qj4[4] = {(double)jbp_den, (double)jbp_pressure, (double)jbp_vx, (double)jbp_vy};
										int q4;
										for(q4 = 0; q4 < 4; q4++){ lo4[q4] = fmin(qi4[q4], qj4[q4]); hi4[q4] = fmax(qi4[q4], qj4[q4]); }
										int kk4[2] = {tmp2->upperrelated, tmp->lowerrelated}, k4;
										for(k4 = 0; k4 < 2; k4++){
											if(kk4[k4] < 0) continue;
											treevorostressrk4particletype *kb = (treevorostressrk4particletype*)(neighwork[kk4[k4]].bp);
											double qk[4] = {(double)kb->den, (double)kb->pressure, (double)kb->vx, (double)kb->vy};
											for(q4 = 0; q4 < 4; q4++){ if(qk[q4] < lo4[q4]) lo4[q4] = qk[q4]; if(qk[q4] > hi4[q4]) hi4[q4] = qk[q4]; }
										}
#define GFS_B4(qv, d, q4_) do { double v_ = (double)(qv) + (double)(d); if(v_ < lo4[q4_]) d = (postype)(lo4[q4_] - (double)(qv)); else if(v_ > hi4[q4_]) d = (postype)(hi4[q4_] - (double)(qv)); } while(0)
										GFS_B4(ibp_den, drho_i, 0); GFS_B4(jbp_den, drho_j, 0);
										GFS_B4(ibp_pressure, dp_i, 1); GFS_B4(jbp_pressure, dp_j, 1);
										GFS_B4(ibp_vx, dvxL, 2); GFS_B4(jbp_vx, dvxR, 2);
										GFS_B4(ibp_vy, dvyL, 3); GFS_B4(jbp_vy, dvyR, 3);
#undef GFS_B4
										(void)tR_; (void)tP_; (void)tVx_; (void)tVy_;
									}
									dvnL0 = nx_hat*dvxL + ny_hat*dvyL;
									dvnR0 = nx_hat*dvxR + ny_hat*dvyR;
								}

								rhoL += drho_i;
								rhoR += drho_j;
								pL += dp_i;
								pR += dp_j;
								vnL_lab += dvnL0;
								vnR_lab += dvnR0;
							}
							if(gfs_tq_on() && !jbp_is_ghost){
								double cx, cy; tq_center(&cx, &cy);
								double rx = (double)ibp_rk4->x - cx, ry = (double)ibp_rk4->y - cy;
								int kb = tq_bin(sqrt(rx*rx + ry*ry));
								double jc = fabs((double)((jbp_vx - ibp_vx)*nx_hat + (jbp_vy - ibp_vy)*ny_hat));
								double jr = fabs((double)(vnR_lab - vnL_lab));
#ifdef _OPENMP
#pragma omp atomic
#endif
								jq_n[kb] += 1;
#ifdef _OPENMP
#pragma omp atomic
#endif
								jq_cell[kb] += jc;
#ifdef _OPENMP
#pragma omp atomic
#endif
								jq_rec[kb] += jr;
								if(jc > 1e-12 && jr > 0.5*jc){
#ifdef _OPENMP
#pragma omp atomic
#endif
									jq_keep[kb] += 1;
								}
							}

							/* Positivity floor */
							int rf_floor = (pL < 1e-10 || pR < 1e-10 || rhoL < 1e-10 || rhoR < 1e-10);
							{
								double qs[4] = {(double)pL, (double)pR, (double)rhoL, (double)rhoR};
								double qc[4] = {(double)ibp_pressure, (double)jbp_pressure, (double)ibp_den, (double)jbp_den};
								int q;
								for(q=0;q<4;q++) if(qs[q] < 0){
									double rel = -qs[q]/(fabs(qc[q]) + 1e-300);
									int b = rel < 1e-12 ? 0 : rel < 1e-6 ? 1 : rel < 1e-2 ? 2 : rel < 1 ? 3 : 4;
#ifdef _OPENMP
#pragma omp atomic
#endif
									gfs_rf_mag[q/2][b]++;
								}
							}
							if(rf_floor){
								int nP = (pL < 0) + (pR < 0), tP = (pL >= 0 && pL < 1e-10) + (pR >= 0 && pR < 1e-10);
								int nR = (rhoL < 0) + (rhoR < 0), tR = (rhoL >= 0 && rhoL < 1e-10) + (rhoR >= 0 && rhoR < 1e-10);
#ifdef _OPENMP
#pragma omp critical (rf_count)
#endif
								{
									gfs_recfloor_n++;
									gfs_rf_negP += nP; gfs_rf_tinyP += tP; gfs_rf_negR += nR; gfs_rf_tinyR += tR;
									gfs_rf_arc += cell_has_arc;
								}
							}
							if(pL < 1e-10) pL = 1e-10;
							if(pR < 1e-10) pR = 1e-10;
							if(rhoL < 1e-10) rhoL = 1e-10;
							if(rhoR < 1e-10) rhoR = 1e-10;

							/* Face rest frame HLLC */
							hllc_face_2d_rest_frame(rhoL, pL, vnL_lab, ibp_csound,
							                        rhoR, pR, vnR_lab, jbp_csound,
							                        wn, Gamma, &pst, &vnst_lab);
							gq_p[gq_k] = (double)pst;
							gq_w[gq_k] = (double)(uradix_ui_boost.x*nx_hat + uradix_ui_boost.y*ny_hat);
							gq_P += (double)pst/(double)gq_n;
							{
								double aps = (double)facearea*(double)pst;
#ifdef _OPENMP
#pragma omp atomic
#endif
								gfs_rf_faces++;
#ifdef _OPENMP
#pragma omp atomic
#endif
								gfs_rf_APs_all += aps;
								if(rf_floor){
#ifdef _OPENMP
#pragma omp atomic
#endif
									gfs_rf_APs_floor += aps;
								}
							}
							hq_pL = pL; hq_pR = pR; hq_uL = vnL_lab; hq_uR = vnR_lab; hq_us = vnst_lab; hq_ok = 1; hq_rhoL = rhoL; hq_rhoR = rhoR;
							if(gfs_tcf_on()){
								if(gq_n != 1 || !sedov_phase1_on() || gfs_schemeH_on() || gfs_dynw2_active){
									fprintf(stderr, "[TCF] requires one face point, SEDOV_PHASE1 total energy, no GFS_SCHEME_H, no GFS_DYN_W2\n");
									exit(1);
								}
								double d2w = (double)dramp*(double)dramp;
								if(fabs((double)ibp_rk4->w2 - (double)((treevorork4particletype*)jbp)->w2) > 1e-12*d2w){
									fprintf(stderr, "[TCF] W_i != W_j on a face: stage 2 is W=0 only\n");
									exit(1);
								}
								/* Z_eff = (P*(du) - P*(0))/du; for |du| below 1e-8 c the
								 * difference quotient at du = 1e-6 c (its limit) */
								double g_ = (double)Gamma, uL_ = vnL_lab, uR_ = vnR_lab, um_ = 0.5*(uL_ + uR_);
								double cs_ = sqrt(g_*pL/rhoL) + sqrt(g_*pR/rhoR);
								double P0_ = tcf_pstar(rhoL, pL, um_, rhoR, pR, um_, g_);
								tcf_durec = uL_ - uR_;
								if(fabs(tcf_durec) > 1e-8*cs_){
									tcf_Z = (tcf_pstar(rhoL, pL, uL_, rhoR, pR, uR_, g_) - P0_)/tcf_durec;
								} else {
									double dl_ = 1e-6*cs_;
									tcf_Z = (tcf_pstar(rhoL, pL, um_ + 0.5*dl_, rhoR, pR, um_ - 0.5*dl_, g_) - P0_)/dl_;
#ifdef _OPENMP
#pragma omp atomic
#endif
									tcf_small++;
								}
								tcf_ok = 1;
							}
							}   /* gq_k */
							pst = (postype)gq_P;
							pi_total = pst;
							if(gq_n == 2) gq_cov = 0.25*(gq_p[1] - gq_p[0])*(gq_w[1] - gq_w[0]);
							/* Extreme faces use one HLL state built from the
							 * cell-centered pair along i→j. MUSCL values differ
							 * on the two ends of the same face, so one end was
							 * taking the shock speed and the other the
							 * geometric speed. The pair then does not cancel
							 * and the blast creates energy. */
							if(sedov_phase1_on() && dramp > 0){
								postype pmin_c = fmin(ibp_pressure, jbp_pressure);
								postype pmax_c = fmax(ibp_pressure, jbp_pressure);
								int face_extreme = (pmax_c > (postype)100.0 * fmax(pmin_c, (postype)1.0e-30));
								if(face_extreme){
									postype ps, vns;
									hll_star_state(ibp_den, ibp_pressure,
										ibp_vx*er.x + ibp_vy*er.y, ibp_csound,
										jbp_den, jbp_pressure,
										jbp_vx*er.x + jbp_vy*er.y, jbp_csound,
										Gamma, &ps, &vns);
									pi_total = ps;
									riemann_vstar = 1;
									riemann_vn = vns;
									riemann_nx = er.x;
									riemann_ny = er.y;
									hq_pL = ibp_pressure; hq_pR = jbp_pressure;
									hq_uL = ibp_vx*er.x + ibp_vy*er.y; hq_uR = jbp_vx*er.x + jbp_vy*er.y;
									hq_us = vns; hq_ok = 1;
								}
							}

							/* Optional Monaghan AV for grid-scale noise control
							   (uses GAS_AlphaVis; disabled when alphavis==0).
							   Applied only on compression to damp cell jitter. */
							if(alphavis > 0){
								Voro2D_point uij;
								uij.x = jbp_vx - ibp_vx;
								uij.y = jbp_vy - ibp_vy;
								postype rvel = Vec2DDotP(&er, &uij);
								if(rvel < 0){
									postype wcomp = sqrt(ibp_rk4->w2)+sqrt(((treevorork4particletype*)jbp)->w2);
									postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
									postype drampScale = dramp/scaleFactor;
									postype mu = rvel /* Voronoi-AV-reform: no length regularization */;
									postype meanden = 0.5*(ibp_den + jbp_den);
									postype meanCsound = 0.5*(ibp_csound + jbp_csound);
									pi_total += (-alphavis*meanCsound*mu + betavis*mu*mu)*meanden;
								}
							}
							/* tau_dot_dS_x/y remain 0 — no NS stress */

						} else {
							/* Two-tier blended path (real-real faces only) */
							postype OoA = OrderOfAccuracy;

							if(do_face_dbg){
#pragma omp critical (face_dbg)
								{
									int ulr = tmp->upperlink ? tmp->upperlink->upperrelated : -999;
									int llr = tmp->lowerrelated;
									void *bp_u = (ulr >= 0 && ulr < nneigh) ? neighwork[ulr].bp : (void*)0xbaadf00d;
									void *bp_l = (llr >= 0 && llr < nneigh) ? neighwork[llr].bp : (void*)0xbaadf00d;
									fprintf(stderr,"[FACE] r%d call=%d corner=%d pre-Mnm OoA=%g ulink=%p ulr=%d llr=%d nneigh=%d bp_u=%p bp_l=%p\n",
										MYID(simpar), s_force_call, corner_n, (double)OoA,
										(void*)tmp->upperlink, ulr, llr, nneigh, bp_u, bp_l);
									fflush(stderr);
								}
							}

							/* === Tier 1: M(n,m) pressure + NS stress ===
							 * BUG FIX (May 2026, Sedov2D step08_K diag):
							 * VoroStressRK4_Pressure2D macro dereferences
							 *   neighwork[tmp->upperlink->upperrelated].bp
							 *   neighwork[tmp->lowerrelated].bp
							 * unconditionally.  When either index is < 0
							 * (no stencil neighbor at corner faces),
							 * neighwork[-4] underflows and SIGSEGVs.
							 * Guard like the tau stencil at the next block:
							 * fall back to simple average. */
							postype p_mnm;
							if(OoA > 0 && (tmp->upperlink)->upperrelated >= 0 && tmp->lowerrelated >= 0){
								p_mnm = VoroStressRK4_Pressure2D(
										ibp_pressure, jbp_pressure, tmp, neighwork, simpar, RT);
							} else {
								p_mnm = 0.5*(ibp_pressure + jbp_pressure);
							}

							if(do_face_dbg){
#pragma omp critical (face_dbg)
								{ fprintf(stderr,"[FACE] r%d call=%d corner=%d post-Mnm p_mnm=%g\n",
									MYID(simpar), s_force_call, corner_n, (double)p_mnm); fflush(stderr); }
							}

							/* Viscous traction vector: τ·dS (full tensor, not just normal projection) */
							tau_dot_dS_x = 0; tau_dot_dS_y = 0;
							if(av_mode == 1 || (av_mode == 0 && nu_phys > 0)){
								/* M(n,m) face-average of τ components */
								postype txx_face, txy_face, tyy_face;
								if(OoA > 0 && (tmp->upperlink)->upperrelated >= 0 && tmp->lowerrelated >= 0){
									txx_face = VoroStressRK4_Scalar2D(
										ibp->stress.tauxx, jbp->stress.tauxx, tmp, neighwork, simpar, RT, stress.tauxx);
									txy_face = VoroStressRK4_Scalar2D(
										ibp->stress.tauxy, jbp->stress.tauxy, tmp, neighwork, simpar, RT, stress.tauxy);
									tyy_face = VoroStressRK4_Scalar2D(
										ibp->stress.tauyy, jbp->stress.tauyy, tmp, neighwork, simpar, RT, stress.tauyy);
								} else {
									txx_face = 0.5*(ibp->stress.tauxx + jbp->stress.tauxx);
									txy_face = 0.5*(ibp->stress.tauxy + jbp->stress.tauxy);
									tyy_face = 0.5*(ibp->stress.tauyy + jbp->stress.tauyy);
								}
								/* Viscous traction: τ_code = -τ_NS, so
								   τ_NS·dS = -(τ_code·dS).  Negate here so that
								   subsequent "+tau_dot_dS" gives correct
								   +τ_NS·dS direction (viscosity damps shear). */
								tau_dot_dS_x = -(txx_face * dS.x + txy_face * dS.y);
								tau_dot_dS_y = -(txy_face * dS.x + tyy_face * dS.y);
							}

							if(av_mode == 0 || av_mode == 1){
								/* M(n,m) pressure + NS stress + optional Monaghan AV
								 * Balsara (1995) switch suppresses Monaghan AV in shear-dominated
								 * regions: f = |divv| / (|divv| + |curlv| + xi). Without it,
								 * Monaghan damps rotation (KH rolls) as well as compressive noise. */
								pi_total = p_mnm;
								if(alphavis > 0){
									Voro2D_point uij;
									uij.x = jbp_vx - ibp_vx;
									uij.y = jbp_vy - ibp_vy;
									postype rvel = Vec2DDotP(&er, &uij);
									if(rvel < 0){
										postype wcomp = sqrt(ibp_rk4->w2)+sqrt(((treevorork4particletype*)jbp)->w2);
										postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
										postype drampScale = dramp/scaleFactor;
										postype mu = rvel /* Voronoi-AV-reform: no length regularization */;
										postype meanden = 0.5*(ibp_den + jbp_den);
										postype meanCsound = 0.5*(ibp_csound + jbp_csound);
										postype divv_i  = ibp->stress.divv;
										postype curlv_i = ibp->stress.gUyx - ibp->stress.gUxy;
										postype divv_j  = jbp->stress.divv;
										postype curlv_j = jbp->stress.gUyx - jbp->stress.gUxy;
										postype xi_b    = 1.0e-4 * meanCsound / dramp;
										postype f_i = fabs(divv_i) / (fabs(divv_i) + fabs(curlv_i) + xi_b);
										postype f_j = fabs(divv_j) / (fabs(divv_j) + fabs(curlv_j) + xi_b);
										postype f_balsara = 0.5*(f_i + f_j);
										postype mon_av_pi = f_balsara * (-alphavis*meanCsound*mu + betavis*mu*mu)*meanden;
										pi_total += mon_av_pi;
										if(mon_av_pi > 0) p_av_for_heat += mon_av_pi;
									}
								}
								/* CD10 viscous pressure (Cullen-Dehnen 2010) — mirror
								 * of av_mode=3 path. Strong-shock detector that
								 * survives Voronoi cell collapse: alpha_cd→1 forces
								 * extra dissipation regardless of w2 size. */
								if(av_mode == 1 && cd_amax > 0){
									Voro2D_point uij_cd;
									uij_cd.x = jbp_vx - ibp_vx;
									uij_cd.y = jbp_vy - ibp_vy;
									postype rvel_cd = Vec2DDotP(&er, &uij_cd);
									if(rvel_cd < 0){
										postype alpha_face = 0.5*(ibp->stress.alpha_cd + jbp->stress.alpha_cd);
										postype rho_mean   = 0.5*(ibp_den + jbp_den);
										postype vsig_cd    = ibp_csound + jbp_csound - rvel_cd;
										postype cd10_pi = 0.5 * alpha_face * vsig_cd * rho_mean * (-rvel_cd);
										pi_total += cd10_pi;
										p_av_for_heat += cd10_pi;
									}
								}
								/* tau_dot_dS already computed above */
							} else {
								/* === av_mode == 2: Two-tier blend ===
								 * Revision 4 (expert review):
								 *   FIX1: use hllc_face_2d_rest_frame (Springel 2010 §32).
								 *         Lab-frame HLLC at Lagrangian contacts gives ~zero
								 *         dissipation and produces surface-tension artifacts.
								 *   FIX2: gate CD10 viscous pressure when f_pq > 0.5
								 *         (HLLC alone provides dissipation; avoid double
                                 *         dissipation positive feedback).
								 *   FIX3: contact-aware MUSCL — skip velocity reconstruction
								 *         when density jump > 10% but pressure continuous
								 *         (< 2%). Shear gradients at smooth contacts otherwise
								 *         create spurious relative v_n for HLLC to dissipate.
								 */

								/* Blending factor: symmetric f_pq = f_qp */
								postype alpha_max_pq = fmax(ibp->stress.alpha_cd, jbp->stress.alpha_cd);
								postype theta_pq = 0;
								if(ibp->stress.divv < -1e-10 || jbp->stress.divv < -1e-10)
									theta_pq = blend_theta;
								postype f_pq = fmin(1.0, alpha_max_pq/(cd_amax + 1e-30) + theta_pq);

								postype p_hllc = p_mnm;  /* fallback */
								postype Pi_cd10 = 0;

								if(f_pq > 1e-6){
									/* HLLC path */
									postype vnL, vnR, pL, pR, rhoL, rhoR, cL, cR;

									/* Project velocities onto face normal */
									postype ds_mag_inv = 1.0/(facearea + 1e-30);
									postype nx_hat = dS.x * ds_mag_inv;
									postype ny_hat = dS.y * ds_mag_inv;
									vnL = ibp_vx*nx_hat + ibp_vy*ny_hat;
									vnR = jbp_vx*nx_hat + jbp_vy*ny_hat;

									/* FIX3: detect pure contact (ρ jump, P continuous).
									 *       At such faces, velocity MUSCL must be
									 *       suppressed to avoid spurious dvn from
									 *       transverse shear gradient. */
									postype rho_mean_det = 0.5*(ibp_den + jbp_den);
									postype p_mean_det   = 0.5*(ibp_pressure + jbp_pressure);
									postype drho_rel = fabs(ibp_den - jbp_den) /
									                   (rho_mean_det + 1e-30);
									postype dp_rel   = fabs(ibp_pressure - jbp_pressure) /
									                   (p_mean_det + 1e-30);
									#define CONTACT_DRHO_THRESH 0.1
									#define CONTACT_DP_THRESH   0.02
									int is_pure_contact = (drho_rel > CONTACT_DRHO_THRESH) && (dp_rel < CONTACT_DP_THRESH);

									if(use_muscl){
										/* MUSCL reconstruction using cell-averaged gradients.
										 * ∇P, ∇v computed via Green-Gauss in updateDenW2Pressure2DBlend.
										 * Reconstruct at face midpoint: Q_L = Q̄_i + ⟨∇Q⟩_i · (x_f - x_i).
										 * tmp->x,tmp->y are stored RELATIVE to generator ibp
										 * (Voro2D_FindVC uses Vec2DSub). */
										postype xf_rel_i = 0.5*(tmp->x + tmp2->x);
										postype yf_rel_i = 0.5*(tmp->y + tmp2->y);
										postype dx_ij_rel = jbp->x - ibp->x;
										postype dy_ij_rel = jbp->y - ibp->y;

										/* Displacement: cell center → face midpoint */
										postype dx_iF = xf_rel_i;
										postype dy_iF = yf_rel_i;
										postype dx_jF = xf_rel_i - dx_ij_rel;
										postype dy_jF = yf_rel_i - dy_ij_rel;

										/* Pressure reconstruction */
										postype dp_i = ibp->stress.dPdx*dx_iF + ibp->stress.dPdy*dy_iF;
										postype dp_j = jbp->stress.dPdx*dx_jF + jbp->stress.dPdy*dy_jF;
										pL = ibp_pressure + dp_i;
										pR = jbp_pressure + dp_j;

										/* Velocity reconstruction (project gradient onto face normal).
										 * FIX3: skip at pure contacts — shear gradient
										 * ∇v is nonzero even with v·n̂ = 0, which fabricates
										 * relative v_n that HLLC would then dissipate. */
										postype dvnL = 0, dvnR = 0;
										if(!is_pure_contact){
											dvnL = (ibp->stress.gUxx*nx_hat + ibp->stress.gUxy*ny_hat)*dx_iF
											     + (ibp->stress.gUyx*nx_hat + ibp->stress.gUyy*ny_hat)*dy_iF;
											dvnR = (jbp->stress.gUxx*nx_hat + jbp->stress.gUxy*ny_hat)*dx_jF
											     + (jbp->stress.gUyx*nx_hat + jbp->stress.gUyy*ny_hat)*dy_jF;
											vnL += dvnL;
											vnR += dvnR;
										}

										/* Simple limiter: clamp to prevent new extrema */
										postype pmin = fmin(ibp_pressure, jbp_pressure);
										postype pmax = fmax(ibp_pressure, jbp_pressure);
										if(pL < pmin) pL = pmin; if(pL > pmax) pL = pmax;
										if(pR < pmin) pR = pmin; if(pR > pmax) pR = pmax;
										postype vnmin = fmin(vnL - dvnL, vnR - dvnR);
										postype vnmax = fmax(vnL - dvnL, vnR - dvnR);
										if(vnL < vnmin) vnL = vnmin; if(vnL > vnmax) vnL = vnmax;
										if(vnR < vnmin) vnR = vnmin; if(vnR > vnmax) vnR = vnmax;
									}

									pL = use_muscl ? pL : ibp_pressure;
									pR = use_muscl ? pR : jbp_pressure;
									rhoL = ibp_den;
									rhoR = jbp_den;
									cL = ibp_csound;
									cR = jbp_csound;

									/* FIX1: Face rest frame HLLC (Arepo-style).
									 *       w_n = ((v_i + v_j)/2) · n̂ boosts the Riemann
									 *       problem into the face's rest frame, where
									 *       contact discontinuities show proper relative
									 *       velocity. Unboost pstar (invariant) and
									 *       vnstar (lab = rest + w_n). */
									postype wx = 0.5*(ibp_vx + jbp_vx);
									postype wy = 0.5*(ibp_vy + jbp_vy);
									postype wn = wx*nx_hat + wy*ny_hat;

									postype pst, vnst_lab;
									hllc_face_2d_rest_frame(rhoL, pL, vnL, cL,
									                        rhoR, pR, vnR, cR,
									                        wn, Gamma,
									                        &pst, &vnst_lab);
									p_hllc = pst;

									/* FIX2: CD10 viscous pressure is active ONLY when
									 * HLLC weight is low (f_pq <= 0.5). When HLLC
									 * dominates the face it provides its own
									 * dissipation; stacking CD10 on top causes
									 * double heating → dt collapse positive feedback. */
									if(f_pq <= 0.5){
										postype dvn = vnR - vnL;
										if(dvn < 0){ /* approaching (in face rest frame same sign) */
											postype alpha_face = 0.5*(ibp->stress.alpha_cd + jbp->stress.alpha_cd);
											postype rho_mean = 0.5*(rhoL + rhoR);
											postype vsig = cL + cR - fmin(0.0, dvn);
											Pi_cd10 = 0.5 * alpha_face * vsig * rho_mean * (-dvn);
										}
									}
								}

								/* Blend: pressure is blended, NS traction scaled by (1-f) */
								pi_total = (1.0-f_pq)*p_mnm + f_pq*(p_hllc + Pi_cd10);
								tau_dot_dS_x *= (1.0-f_pq);
								tau_dot_dS_y *= (1.0-f_pq);
								/* Fix B: extract dissipative AV part — HLLC star excess + CD10.
								 *   p_hllc - p_mnm  is the Riemann jump dissipation
								 *   (Toro 2009 §10.5).  Clamped at 0 so rarefaction
								 *   stars (p_hllc < p_mnm) don't subtract — entropy
								 *   non-decrease is enforced at the source. */
								postype Pi_hllc_extra = (p_hllc > p_mnm ? (p_hllc - p_mnm) : 0);
								p_av_for_heat = f_pq * (Pi_hllc_extra + Pi_cd10);
							}
						}

						/* === Accumulate forces and energy rates === */
						/* Internal energy: reversible PdV + AV heating + NS viscous + heat conduction.
						 *
						 * Fix A+B redesign (entropy non-decrease + energy conservation joint):
						 *   pi_total = p_rev + p_av_for_heat   (decomposition)
						 *     p_rev          := p_mnm    (face-symmetric, reversible)
						 *     p_av_for_heat  := Monaghan AV + CD10 + HLLC-jump-excess   (dissipative)
						 *
						 *   FACE GATE: if u_n_compress > 0 (approaching cells):
						 *     - F_AV acts (pi_total includes p_av_for_heat) → KE drains
						 *     - ε_AV_face = p_av · u_n · |dS|  (≥0 by gate)
						 *     - die_AV_i += ½·ε_AV   (pair-symmetric: j gets the other ½ in its loop)
						 *     - dK_diss += ½·ε_AV    (AV → K via dissipative channel)
						 *   else (expansion / p_av≤0):
						 *     - pi_total reverted to p_rev only (F_AV gated off)
						 *     - ε_AV = 0
						 *   This jointly satisfies pair-energy conservation AND per-face entropy non-decrease.
						 *
						 *   Reversible PdV (-p_rev · u_rad · dS):
						 *     - enters die (em=0 ie integration)
						 *     - does NOT enter dK_diss (K invariant under adiabatic PdV;
						 *       Springel & Hernquist 2002 MNRAS 333, 649)
						 *
						 *   NS viscous stress (face-pair exact, τ·dS·u_rad):
						 *     - enters both die and dK_diss
						 *     - per-cell sum ≥ 0 in well-resolved limit (τ:∇v ≥ 0)
						 *     - per-face can be either sign; cell-level clamp applied after loop.
						 */
						Voro2D_point uradix_ui = gfs_upqrad_rk4(ibp_rk4, (treevorork4particletype*)jbp, dtold);
						voronoi_face_rotation(&uradix_ui,
								(postype)0.5*(tmp->x + tmp2->x),
								(postype)0.5*(tmp->y + tmp2->y),
								ibp->x, ibp->y, jbp->x, jbp->y,
								ibp_vx, ibp_vy, jbp_vx, jbp_vy,
								er.x, er.y, dramp,
								ibp_rk4->w2,
								((treevorork4particletype*)jbp)->w2, lag_rot);
						/* Stationary hard wall: the face velocity has no normal part. */
						if(hard_wall){
							postype vn_i = ibp_vx*hard_nx + ibp_vy*hard_ny;
							postype urn = uradix_ui.x*hard_nx + uradix_ui.y*hard_ny;
							postype corr = urn + vn_i;
							uradix_ui.x -= corr*hard_nx;
							uradix_ui.y -= corr*hard_ny;
						}

						/* Cell-centered HLL. With MUSCL on, av_mode 5 already built
						 * P* from limited face states; do not overwrite it.
						 * Ghost faces stay on the wall path. */
						if(sedov_phase1_on() && !use_muscl && !jbp_is_ghost && facearea > 0){
							postype ps, vns;
							postype invA = (postype)1.0 / facearea;
							postype nx = dS.x * invA;
							postype ny = dS.y * invA;
							p_av_for_heat = 0;
							tau_dot_dS_x = 0;
							tau_dot_dS_y = 0;
							hll_star_state(ibp_den, ibp_pressure, ibp_vx*nx + ibp_vy*ny, ibp_csound,
									jbp_den, jbp_pressure, jbp_vx*nx + jbp_vy*ny, jbp_csound,
									Gamma, &ps, &vns);
							pi_total = ps;
							phase1_extreme = 1;
							phase1_vn = vns;
							phase1_nx = nx;
							phase1_ny = ny;
						}

						/* Compression gate on (v_i - v_j)·n̂.  er already points i→j. */
						/* Restoring pressure along this face. Zero at and above the
						 * local cutoff, and still present when Δv = 0. It is
						 * part of pi_total, so the pair work stays in dE. */
						postype tq_ppair = 0;
						int tq_ext = (riemann_vstar || phase1_extreme);
						if(!jbp_is_ghost){
							postype vol_j = ((treevorork4particletype*)jbp)->volume;
							double vclose = -((jbp_vx - ibp_vx)*er.x + (jbp_vy - ibp_vy)*er.y);
							if(vclose < 0) vclose = 0;
							double qpair = (jbp_vx - ibp_vx)*er.x
									+ (jbp_vy - ibp_vy)*er.y;
							postype p_pair = (postype)gfs_pair_pressure(
									(double)dramp,
									(double)ibp_rk4->volume, (double)vol_j,
									(double)ibp_den, (double)jbp_den,
									(double)ibp_csound, (double)jbp_csound,
									vclose);
							p_pair = (postype)gfs_pair_work_limit(
									(double)p_pair, (double)facearea,
									(double)GAS_dtold(simpar),
									(double)ibp->mass, (double)jbp->mass,
									qpair,
									(double)ibp->ie, (double)jbp->ie, 6.0);
							/* Charge cap: s_i from the face velocity this end
							 * will use for the work (ua below), s_j = qA - s_i. */
							if(p_pair > 0){
								double uax, uay;
								if(phase1_extreme && !geom_fv){
									uax = phase1_vn*phase1_nx; uay = phase1_vn*phase1_ny;
								} else if(riemann_vstar && !geom_fv){
									uax = riemann_vn*riemann_nx; uay = riemann_vn*riemann_ny;
								} else {
									uax = ibp_vx + uradix_ui.x; uay = ibp_vy + uradix_ui.y;
								}
								double s_i = (uax - ibp_vx)*dS.x + (uay - ibp_vy)*dS.y;
								double qA  = (jbp_vx - ibp_vx)*dS.x + (jbp_vy - ibp_vy)*dS.y;
								{
								postype pp0 = p_pair;
								p_pair = (postype)gfs_pair_charge_limit(
										(double)p_pair, s_i, qA - s_i,
										(double)GAS_dtold(simpar),
										(double)ibp->ie, (double)jbp->ie, 6.0);
								if(p_pair < pp0) faud_flags |= 8;
								}
							}
							if(p_pair > 0){
								pi_total += p_pair;
								tq_ppair = p_pair;
								faud_flags |= 4;
#ifdef _OPENMP
#pragma omp atomic
#endif
								gfs_pair_face_ends++;
							}
						}

						/* GFS_FACE_CHARGE_LIMIT (default off). */
						if(gfs_face_charge_on() && !jbp_is_ghost && pi_total > 0){
							double uax, uay;
							if(phase1_extreme && !geom_fv){
								uax = phase1_vn*phase1_nx; uay = phase1_vn*phase1_ny;
							} else if(riemann_vstar && !geom_fv){
								uax = riemann_vn*riemann_nx; uay = riemann_vn*riemann_ny;
							} else {
								uax = ibp_vx + uradix_ui.x; uay = ibp_vy + uradix_ui.y;
							}
							double s_i = (uax - ibp_vx)*dS.x + (uay - ibp_vy)*dS.y;
							double qA  = (jbp_vx - ibp_vx)*dS.x + (jbp_vy - ibp_vy)*dS.y;
							double pc = gfs_pair_charge_limit((double)pi_total, s_i, qA - s_i,
									(double)GAS_dtold(simpar),
									(double)ibp->ie, (double)jbp->ie, 6.0);
							if(pc < (double)pi_total){
								double f = pc / (double)pi_total;
								faud_flags |= 16;
								p_av_for_heat *= f;
								pi_total = (postype)pc;
							}
						}

						postype u_n_compress;
						{
							Voro2D_point dv_av;
							dv_av.x = jbp_vx - ibp_vx;
							dv_av.y = jbp_vy - ibp_vy;
							u_n_compress = -Vec2DDotP(&dv_av, &er);
						}

						/* Apply AV gate: when not compressing, strip p_av_for_heat from pi_total.
						 * This is the key change that makes BOTH energy and entropy hold. */
						postype eps_av_face = 0;
						if(p_av_for_heat > 0 && u_n_compress > 0){
							eps_av_face = p_av_for_heat * u_n_compress * facearea;
						} else if(p_av_for_heat != 0){
							pi_total -= p_av_for_heat;   /* gate dissipative force off */
							p_av_for_heat = 0;
						}

						/* Reversible PdV from face-symmetric reference pressure.
						 *   p_rev_face = pi_total - p_av_for_heat
						 * In av_mode==0/1/2:  p_rev_face = p_mnm (M(n,m) symmetric centered avg)
						 * In av_mode==3/5/ghost:  p_av_for_heat=0 → p_rev_face = pi_total
						 *   (preserves OLD behavior for those paths)
						 * This term does NOT touch dK_diss (Springel & Hernquist 2002). */
						{
							postype p_rev_face = pi_total - p_av_for_heat;
							postype die_rev_face = -p_rev_face * Vec2DDotP(&uradix_ui, &dS) - (postype)gq_cov*facearea;
							die += die_rev_face;
							(void)dK_pdv_factor;  /* deprecated knob — superseded by ε_AV routing */
						}
						/* Pair-symmetric AV heating (½ here, ½ in j's loop) — ≥0 by construction */
						if(eps_av_face > 0){
							postype half_eps = 0.5 * eps_av_face;
							die     += half_eps;
							dK_diss += half_eps;
						}
						/* NS viscous heating: face-pair exact */
						{
							postype die_visc_face = tau_dot_dS_x * uradix_ui.x + tau_dot_dS_y * uradix_ui.y;
							die     += die_visc_face;
							dK_diss += die_visc_face;
						}

						/* Heat conduction: Q = χ ρ_face (Tj-Ti)/d_ij * facearea — dissipative.
						 * dramp is floored at heff (= cell radius) — same floor used by CFL —
						 * to avoid divide-by-near-zero when Voronoi tessellation produces
						 * degenerate face geometry (closely-spaced generators / S-H clipping
						 * residual / ring-IC ring-boundary near-duplicate). Without this,
						 * a single bad face can inflate dK_diss by ~10²⁰ and blow up K. */
						if(nu_phys > 0 && prandtl > 0 && !jbp_is_ghost && !sedov_phase1_on()){
							postype chi = nu_phys / prandtl;
							postype Ti = ibp_pressure / ibp_den;
							postype Tj = jbp_pressure / jbp_den;
							postype rho_face = 0.5*(ibp_den + jbp_den);
							postype dramp_cond = fmax(dramp, 0.25*sqrt(ibp_rk4->volume));
							postype die_cond_face = chi * rho_face * (Tj - Ti) / dramp_cond * facearea;
							die += die_cond_face;
							dK_diss += die_cond_face;
						}

						/* Total energy: ua = v_face = v_i + (v_face - v_i)
						   ensures dte = die + dke consistency */
						Voro2D_point ua;
						if(phase1_extreme && !geom_fv){
							ua.x = phase1_vn * phase1_nx;
							ua.y = phase1_vn * phase1_ny;
						} else if(riemann_vstar && !geom_fv){
							ua.x = riemann_vn * riemann_nx;
							ua.y = riemann_vn * riemann_ny;
						} else {
							ua.x = ibp_vx + uradix_ui.x;
							ua.y = ibp_vy + uradix_ui.y;
						}
						double hq_ftx = 0, hq_fty = 0;
						if(tcf_ok){
							/* TCF (W=0): F_i = -A[(Pbar + Q) e + (P_j - P_i) t/d], t = c - m;
							 * dE_i = -A[(P_j u_i + P_i u_j)/2 + (P_j v_i - P_i v_j).t/d + Q (u_i + u_j)/2].
							 * Arc and cap edges add nothing: their force and work are absorbed here. */
							double d = (double)dramp, ex = (double)er.x, ey = (double)er.y, A = (double)facearea;
							double tx = sh_cx - 0.5*(double)dr.x;
							double ty = sh_cy - 0.5*(double)dr.y;
							double te = tx*ex + ty*ey; tx -= te*ex; ty -= te*ey;
							double Pi_ = (double)ibp_pressure, Pj_ = (double)jbp_pressure;
							double vix = (double)ibp_vx, viy = (double)ibp_vy, vjx = (double)jbp_vx, vjy = (double)jbp_vy;
							double ui = vix*ex + viy*ey, uj = vjx*ex + vjy*ey, ducell = ui - uj;
							double qin = tcf_durec, qin_del = 0, qin_chi = 0;   /* qin_del, qin_chi: face dump only */
							if(gfs_tcf_qbp()){
								size_t ii_ = (size_t)PINDX(ibp_rk4), jj_ = (size_t)PINDX((treevorork4particletype*)jbp);
								if(ii_ < qbp_cap && jj_ < qbp_cap && qbp_ok[ii_] && qbp_ok[jj_]){
									const double *Gi = qbp_g + 4*ii_, *Gj = qbp_g + 4*jj_;
									double Gx_[4];
									if(gfs_qbp_edge_inner() && qbp_arc[ii_] != qbp_arc[jj_]){
										/* x_j - x_i = d e: the arc cell sits at r = -d e (arc i) or +d e (arc j) from the arc-free cell */
										size_t in_ = qbp_arc[ii_] ? jj_ : ii_;
										const double *Gn = qbp_g + 4*in_, *Hn = qbp_h + 6*in_;
										double rx = (qbp_arc[ii_] ? -d : d)*ex, ry = (qbp_arc[ii_] ? -d : d)*ey;
										if(gfs_qbp_edge_inner() >= 2){
											Gx_[0] = Gn[0] + Hn[0]*rx + Hn[1]*ry; Gx_[1] = Gn[1] + Hn[1]*rx + Hn[2]*ry;
											Gx_[2] = Gn[2] + Hn[3]*rx + Hn[4]*ry; Gx_[3] = Gn[3] + Hn[4]*rx + Hn[5]*ry;
											if(gfs_qbp_edge_inner() == 3){
												/* + (1/2) T r r: one-sided quadratic extrapolation leaves an O(d^3) error of fixed sign in curved shear */
												const double *Tn = qbp_t + 8*in_; int a_;
												for(a_ = 0; a_ < 2; a_++){
													const double *T_ = Tn + 4*a_;
													Gx_[2*a_]   += 0.5*(T_[0]*rx*rx + 2*T_[1]*rx*ry + T_[2]*ry*ry);
													Gx_[2*a_+1] += 0.5*(T_[1]*rx*rx + 2*T_[2]*rx*ry + T_[3]*ry*ry);
												}
											}
										} else { Gx_[0] = Gn[0]; Gx_[1] = Gn[1]; Gx_[2] = Gn[2]; Gx_[3] = Gn[3]; }
										if(qbp_arc[ii_]){ Gi = Gx_; Gj = Gn; } else { Gi = Gn; Gj = Gx_; }
#ifdef _OPENMP
#pragma omp atomic
#endif
										qbp_nmixed++;
									}
									double eGe = ex*(Gi[0] + Gj[0])*ex + ex*(Gi[1] + Gj[1])*ey + ey*(Gi[2] + Gj[2])*ex + ey*(Gi[3] + Gj[3])*ey;
									qin = ducell + 0.5*d*eGe; qin_del = qin;
									if(gfs_qbp_chi()){
										/* Astra audit 15: rotation-invariant compression share chi = min(1, [-div v]_+ / |S|_F)
										 * of each cell's unlimited gradient; q_in = (1-chi) delta_QBP + chi du_cell */
										double chi_f = 0; int c_;
										for(c_ = 0; c_ < 2; c_++){
											const double *G_ = c_ ? Gj : Gi;
											double th_ = G_[0] + G_[3], s01 = 0.5*(G_[1] + G_[2]);
											double sF = sqrt(G_[0]*G_[0] + G_[3]*G_[3] + 2*s01*s01), ch_ = 0;
											if(sF > 0){ ch_ = (-th_ > 0 ? -th_ : 0)/sF; if(ch_ > 1) ch_ = 1; }
											if(ch_ > chi_f) chi_f = ch_;
										}
										qin = (1 - chi_f)*qin + chi_f*ducell; qin_chi = chi_f;
									}
								} else {
									qin = 0;   /* invalid gradient: gate failure (counted), no Q on this face */
#ifdef _OPENMP
#pragma omp atomic
#endif
									qbp_nfail++;
								}
#ifdef _OPENMP
#pragma omp atomic
#endif
								qbp_nface++;
							}
							double mm = (qin*ducell <= 0) ? 0 : (fabs(qin) < fabs(ducell) ? qin : ducell);
							double Zq = tcf_Z, fz_fi = 1, fz_fj = 1;
							if(gfs_tcf_face_z()){
								size_t ii_ = (size_t)PINDX(ibp_rk4), jj_ = (size_t)PINDX((treevorork4particletype*)jbp);
								double dli = (ii_ < cen_cap) ? (cen_dx[ii_]*ex + cen_dy[ii_]*ey) : 0;
								double dlj = (jj_ < cen_cap) ? (cen_dx[jj_]*ex + cen_dy[jj_]*ey) : 0;
								double fi = 1 - 2*dli/d, fj = 1 + 2*dlj/d;
								if(fi < 1e-12 || fj < 1e-12){
#ifdef _OPENMP
#pragma omp atomic
#endif
									fz_nclamp++;
								}
								if(fi < 1e-12) fi = 1e-12;
								if(fj < 1e-12) fj = 1e-12;
								double g_ = (double)Gamma;
								double rLf = (double)ibp_den*fi, rRf = (double)jbp_den*fj;
								double pLf = Pi_*pow(fi, g_), pRf = Pj_*pow(fj, g_);
								double cL = sqrt(g_*pLf/rLf), cR = sqrt(g_*pRf/rRf), cM = cL > cR ? cL : cR;
								double mu = rLf*rRf/(rLf + rRf);
								double ZD = mu*(cM + (mm > 0 ? mm : 0));
								double aH = (g_ + 1)*fabs(mm)/8, ZH = mu*(aH + sqrt(aH*aH + cM*cM));
								Zq = ZD > ZH ? ZD : ZH; fz_fi = fi; fz_fj = fj;
							} else
							if(gfs_tcf_zh()){
								double ri = (double)ibp_den, rj = (double)jbp_den, rh = 2*ri*rj/(ri + rj + 1e-300);
								double cmx = (double)(ibp_csound > jbp_csound ? ibp_csound : jbp_csound);
								double aH = ((double)Gamma + 1)*fabs(mm)/8;
								double ZH = 0.5*rh*(aH + sqrt(aH*aH + cmx*cmx));
								if(ZH > Zq) Zq = ZH;
							}
							double Q = gfs_tcf_q0() ? 0
								: (gfs_tcf_qlaw() == 2 ? (ducell != 0 ? Zq*mm*mm/ducell : 0) : Zq*mm);
							double Pb = 0.5*(Pi_ + Pj_);
							double ftx_ = -A*(Pj_ - Pi_)*tx/d, fty_ = -A*(Pj_ - Pi_)*ty/d;
							fx += -A*(Pb + Q)*ex;
							fy += -A*(Pb + Q)*ey;
							if(gfs_tcf_notan()){
#ifdef _OPENMP
#pragma omp atomic
#endif
								tcf_tanwork += vix*ftx_ + viy*fty_;
							} else { fx += ftx_; fy += fty_; }
							dte += -A*(0.5*(Pj_*ui + Pi_*uj) + (Pj_*(vix*tx + viy*ty) - Pi_*(vjx*tx + vjy*ty))/d + 0.5*Q*(ui + uj));
							tcf_tds += 0.5*A*Q*ducell;
							tcf_dV += A*(0.5*(ui + uj) + ((vix - vjx)*tx + (viy - vjy)*ty)/d);
							if(tcf_dump_fp()){
								/* closure and bisector residuals of this face */
								double nx_ = (double)dS.x/A, ny_ = (double)dS.y/A;
								tg_nx += (double)dS.x; tg_ny += (double)dS.y; tg_ex += A*ex; tg_ey += A*ey; tg_nf++;
								double en_ = fabs(ex*ny_ - ey*nx_); if(en_ > tg_en) tg_en = en_;
								double b1 = fabs((double)tmp->x*ex + (double)tmp->y*ey - 0.5*d), b2 = fabs((double)tmp2->x*ex + (double)tmp2->y*ey - 0.5*d);
								if(b1/d > tg_bis) tg_bis = b1/d; if(b2/d > tg_bis) tg_bis = b2/d;
								if(A/d < tg_amin) tg_amin = A/d;
								if(tcf_fdump){
#pragma omp critical (tcf_fdump)
									fprintf(tcf_fdump, "%ld %ld %.10g %.10g %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
										(long)PINDX(ibp_rk4), (long)PINDX((treevorork4particletype*)jbp), (double)ibp_rk4->x, (double)ibp_rk4->y,
										A, d, ex, ey, tcf_durec, ducell, tcf_Z, (double)ibp_csound, hq_uL, hq_uR, Q, Zq, fz_fi, fz_fj, qin, qin_del, qin_chi);   /* cols 16-19: applied Zq, face-state factors, Q input; 20-21: QBP delta before chi, chi */
								}
							}
#ifdef _OPENMP
#pragma omp critical (tcf_count)
#endif
							{
								tcf_faces++;
								tcf_heat += 0.5*A*Q*ducell;
								if(tcf_Z < 0){ tcf_zneg++; if(tcf_Z < tcf_zmin) tcf_zmin = tcf_Z; }
								if(mm == 0 && tcf_durec != 0) tcf_mm0++;
							}
						} else
						if(gfs_schemeH_on() && gfs_dynw2_active){
							/* (H): anchor a = x_i + phi d e, t = c_f - a (both relative to x_i). */
							double d = (double)dramp, ex = (double)er.x, ey = (double)er.y;
							double Wi = (double)ibp_rk4->w2, Wj = (double)((treevorork4particletype*)jbp)->w2;
							double phi = 0.5 + (Wi - Wj)/(2*d*d);
							double tx = 0.5*((double)tmp->x + (double)tmp2->x) - phi*d*ex;
							double ty = 0.5*((double)tmp->y + (double)tmp2->y) - phi*d*ey;
							double Pi_ = (double)ibp_pressure, Pj_ = (double)jbp_pressure;
							double A = (double)facearea;
							hq_ftx = -A*(Pj_ - Pi_)*tx/d;
							hq_fty = -A*(Pj_ - Pi_)*ty/d;
							double fxf = -(double)pi_total*(double)dS.x + hq_ftx;
							double fyf = -(double)pi_total*(double)dS.y + hq_fty;
							double wi = dyn_wdot((treevorork4particletype*)ibp_rk4);
							double wj = (PINDX((treevorork4particletype*)jbp) == MAX_INDEX) ? wi
								: dyn_wdot((treevorork4particletype*)jbp);
							double delta = (wi - wj)/(2*d);
							double ui = (double)ibp_vx*ex + (double)ibp_vy*ey;
							double uj = (double)jbp_vx*ex + (double)jbp_vy*ey;
							double Ptil = phi*Pj_ + (1 - phi)*Pi_;
							double Df = A*((Pi_ - Pj_)*delta + ((double)pi_total - Ptil)*(ui - uj));
							double th = 0.5;
							if(hq_ok){
								double gi = -((double)pi_total - hq_pL)*(hq_us - hq_uL); if(gi < 0) gi = 0;
								double gj =  ((double)pi_total - hq_pR)*(hq_us - hq_uR); if(gj < 0) gj = 0;
								if(gi + gj > 0) th = gi/(gi + gj);
							}
							{
								double d1 = A*(Pi_ - Pj_)*delta, d2 = A*((double)pi_total - Ptil)*(ui - uj);
#ifdef _OPENMP
#pragma omp atomic
#endif
								hq_dpos += (Df > 0 ? Df : 0);
#ifdef _OPENMP
#pragma omp atomic
#endif
								hq_dneg += (Df < 0 ? -Df : 0);
#ifdef _OPENMP
#pragma omp atomic
#endif
								hq_d1p += (d1 > 0 ? d1 : 0);
#ifdef _OPENMP
#pragma omp atomic
#endif
								hq_d1n += (d1 < 0 ? -d1 : 0);
#ifdef _OPENMP
#pragma omp atomic
#endif
								hq_d2p += (d2 > 0 ? d2 : 0);
#ifdef _OPENMP
#pragma omp atomic
#endif
								hq_d2n += (d2 < 0 ? -d2 : 0);
							}
							dte += (double)ibp_vx*fxf + (double)ibp_vy*fyf
							     - Pi_*(double)Vec2DDotP(&ua, &dS) + th*Df
							     + tau_dot_dS_x * ua.x + tau_dot_dS_y * ua.y;
						} else
						dte += -pi_total * Vec2DDotP(&ua, &dS)
						     + tau_dot_dS_x * ua.x + tau_dot_dS_y * ua.y;
						dvb += Vec2DDotP(&ua, &dS);
						if(gfs_dynw2_active && riemann_vstar && !jbp_is_ghost){
							/* realised face speed minus HLL v* (independent check) */
							double dv = (double)(ua.x*riemann_nx + ua.y*riemann_ny) - (double)riemann_vn;
#ifdef _OPENMP
#pragma omp atomic
#endif
							dyn_fchk_s2 += dv*dv;
#ifdef _OPENMP
#pragma omp atomic
#endif
							dyn_fchk_n += 1;
						}
						if(gfs_faud_on() && !jbp_is_ghost){
							size_t fs = (size_t)PINDX(ibp_rk4);
							if(fs < faud_cap){
								int fk = faud_n[fs]++;
								if(fk < FAUD_MAX){
									size_t fa = fs*FAUD_MAX + fk;
									faud_j[fa] = (long)PINDX((treevorork4particletype*)jbp);
									faud_f[fa] = -(double)pi_total * (double)Vec2DDotP(&ua, &dS);
									faud_p[fa] = (double)pi_total;
									faud_u[fa] = (double)Vec2DDotP(&ua, &dS);
									faud_fl[fa] = faud_flags | (riemann_vstar ? 1 : 0) | (phase1_extreme ? 2 : 0);
								}
							}
						}
						/* Volume swept by the Riemann contact. Interior faces cancel. */
						if(riemann_vstar && sedov_lagvol_on() && lag_sw &&
								(size_t)PINDX(ibp_rk4) < lag_vt_cap){
							lag_sw[(size_t)PINDX(ibp_rk4)] += Vec2DDotP(&ua, &dS);
						}

						/* Kinetic energy */
						Voro2D_point ub;
						ub.x = ibp_vx;
						ub.y = ibp_vy;
						dke += -pi_total * Vec2DDotP(&ub, &dS)
						     + tau_dot_dS_x * ub.x + tau_dot_dS_y * ub.y;

						/* Force: -p·dS + τ·dS */
						if(!tcf_ok){
						fx += -pi_total * dS.x + tau_dot_dS_x;
						fy += -pi_total * dS.y + tau_dot_dS_y;
						}
						fx += hq_ftx; fy += hq_fty;   /* (H) VPH tangential term; 0 otherwise */
						if(gfs_edge_on() && !jbp_is_ghost){
							double W_ = gfs_env_d("GFS_EDGE_DUMP_W", 0.05), T1 = gfs_env_d("GFS_EDGE_DUMP_T1", 0.1), T2 = gfs_env_d("GFS_EDGE_DUMP_T2", 5.0);
							if((gfs_tnow >= T1 && gfs_tnow < T1 + W_) || (gfs_tnow >= T2 && gfs_tnow < T2 + W_)){
								double cx = gfs_env_d("EUNHA_KEPLER_CX", 0), cy = gfs_env_d("EUNHA_KEPLER_CY", 0);
								double rx = (double)ibp_rk4->x - cx, ry = (double)ibp_rk4->y - cy, ri = sqrt(rx*rx + ry*ry);
								int slab = (int)gfs_env_d("GFS_EDGE_DUMP_SLAB", 0);   /* select |x - CX| > R instead of r < R */
								if(slab) ri = fabs(rx);
								if(slab ? (ri > gfs_env_d("GFS_EDGE_DUMP_R", 0.56)) : (ri < gfs_env_d("GFS_EDGE_DUMP_R", 0.56))){
									double jx = (double)((treevorork4particletype*)jbp)->x - cx, jy = (double)((treevorork4particletype*)jbp)->y - cy;
									double A = (double)facearea, nx = (double)dS.x/(A + 1e-300), ny = (double)dS.y/(A + 1e-300);
									double cxf = 0.5*((double)tmp->x + (double)tmp2->x), cyf = 0.5*((double)tmp->y + (double)tmp2->y);
									double toff = fabs(-cxf*ny + cyf*nx);   /* tangential offset of the face centre from the pair axis */
									double du_cell = ((double)jbp_vx - (double)ibp_vx)*nx + ((double)jbp_vy - (double)ibp_vy)*ny;
									double du_rec = hq_ok ? (hq_uR - hq_uL) : du_cell;
									double tqn = -(double)pi_total*(rx*(double)dS.y - ry*(double)dS.x);
									double tqt = rx*hq_fty - ry*hq_ftx;
#ifdef _OPENMP
#pragma omp critical (edge_dump)
#endif
									{
										static FILE *fe = NULL;
										if(!fe){ fe = fopen("edge_faces.txt", "a"); if(fe) fprintf(fe, "# t call id_i r_i id_j r_j A d toff du_n_cell du_n_rec cs_i P_i Pstar extreme tq_normal tq_tangential rhoL rhoR pL pR uL uR ua_n flux_E force_x vx_i\n"); }
										if(fe) fprintf(fe, "%.5f %ld %ld %.5f %ld %.5f %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %d %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",
											gfs_tnow, gfs_edge_call, (long)PINDX(ibp_rk4), ri, (long)PINDX((treevorork4particletype*)jbp), sqrt(jx*jx + jy*jy),
											A, (double)dramp, toff, du_cell, du_rec, (double)ibp_csound, (double)ibp_pressure, (double)pi_total,
											(int)(phase1_extreme != 0 || riemann_vstar != 0), tqn, tqt,
											hq_rhoL, hq_rhoR, hq_pL, hq_pR, hq_uL, hq_uR,
											((double)ua.x*nx + (double)ua.y*ny), -(double)pi_total*(double)Vec2DDotP(&ua, &dS), -(double)pi_total*(double)dS.x, (double)ibp_vx);
									}
								}
							}
						}
						if(gfs_tq_on()){
							double cx, cy; tq_center(&cx, &cy);
							double rx = (double)ibp_rk4->x - cx, ry = (double)ibp_rk4->y - cy;
							int kb = tq_bin(sqrt(rx*rx + ry*ry));
							double arm = rx*(double)dS.y - ry*(double)dS.x;   /* (r x dS)_z */
							double tq_pair = -(double)tq_ppair * arm;
							double tq_riem = -((double)pi_total - (double)tq_ppair) * arm;
#ifdef _OPENMP
#pragma omp atomic
#endif
							tq_acc[kb][tq_ext ? 1 : 0] += tq_riem;
#ifdef _OPENMP
#pragma omp atomic
#endif
							tq_acc[kb][2] += tq_pair;
							{
								double dun = hq_ok ? (hq_uR - hq_uL)
									: ((double)(jbp_vx - ibp_vx)*(double)dS.x + (double)(jbp_vy - ibp_vy)*(double)dS.y)/((double)facearea + 1e-300);
								int k4 = 2*(dun < 0) + (arm > 0);
								double cs_ = (double)ibp_csound > 0 ? (double)ibp_csound : 1e-30;
#ifdef _OPENMP
#pragma omp atomic
#endif
								tqs_acc[kb][k4] += tq_riem;
#ifdef _OPENMP
#pragma omp atomic
#endif
								tqs_n[kb][k4] += 1;
#ifdef _OPENMP
#pragma omp atomic
#endif
								tqs_pp[kb][k4] += (double)pi_total/((double)ibp_pressure + 1e-300);
#ifdef _OPENMP
#pragma omp atomic
#endif
								tqs_du[kb][k4] += fabs(dun)/cs_;
							}
						}

						/* Signal velocity for CFL and CD10 */
						Voro2D_point dv;
						dv.x = jbp_vx - ibp_vx;
						dv.y = jbp_vy - ibp_vy;
						postype VdotR = Vec2DDotP(&dv, &er);
						postype vsig = jbp_csound + ibp_csound - MIN(0, VdotR);
						/* A flattened real cell must use its true thickness.
						 * Flooring dramp up to sqrt(volume) let a Sedov shell
						 * particle step across its neighbor. Ghost mirrors
						 * keep the floor so a particle on the wall does not
						 * drive dt to zero. */
						postype heff = 0.25*sqrt(ibp_rk4->volume);
						postype dramp_cfl = (jbp_is_ghost) ? fmax(dramp, heff) : dramp;
						postype dt = 2*Courant*dramp_cfl/vsig;
						/* faces between two frozen Kepler-background particles never move: exclude them from the dt limits */
						if(gfs_kepler_frozen(ibp_rk4->x, ibp_rk4->y) && gfs_kepler_frozen(jbp->x, jbp->y)) dt = (postype)1e300;
						if((double)dt < gfs_dtb[0]){
#ifdef _OPENMP
#pragma omp critical (gfs_dtb)
#endif
							{ if((double)dt < gfs_dtb[0]) gfs_dtb[0] = (double)dt; }
						}
						if(gfs_tcf_on() && (gfs_hg_affine() || gfs_tcf_face_z()) && !jbp_is_ghost){
							/* design 10: dt_geom = 0.25 h_f/(c_f + |dv.e|), h_f = d/2, c_f the face subzone sound speed */
							double cfi = (double)ibp_csound, cfj = (double)jbp_csound;
							if(gfs_tcf_face_z()){
								size_t ii_ = (size_t)PINDX(ibp_rk4), jj_ = (size_t)PINDX((treevorork4particletype*)jbp);
								double dli = (ii_ < cen_cap) ? (cen_dx[ii_]*er.x + cen_dy[ii_]*er.y) : 0;
								double dlj = (jj_ < cen_cap) ? (cen_dx[jj_]*er.x + cen_dy[jj_]*er.y) : 0;
								double fi = 1 - 2*dli/(double)dramp, fj = 1 + 2*dlj/(double)dramp;
								if(fi > 0) cfi *= pow(fi, 0.5*((double)Gamma - 1));
								if(fj > 0) cfj *= pow(fj, 0.5*((double)Gamma - 1));
							}
							double dtg = 0.25*0.5*(double)dramp/((cfi > cfj ? cfi : cfj) + fabs((double)VdotR));
							if(dt >= (postype)1e300) dtg = 1e300;   /* frozen-frozen face */
							if(dtg < gfs_dtb[1]){
#ifdef _OPENMP
#pragma omp critical (gfs_dtb)
#endif
								{ if(dtg < gfs_dtb[1]) gfs_dtb[1] = dtg; }
							}
							if(dtg < (double)dt) dt = (postype)dtg;
						}

						/* dt3: 0.1 * d / |Δv|. The old 0.1*d/|Δv/2|^2 has units time^2/length. */
						postype uij_x_half = 0.5*(jbp_vx - ibp_vx);
						postype uij_y_half = 0.5*(jbp_vy - ibp_vy);
						postype uij2 = uij_x_half*uij_x_half + uij_y_half*uij_y_half;
						if(uij2 > 0 && dt < (postype)1e300){   /* skip frozen-frozen faces */
							/* GFS_DT3_COEF (default 0.1): for time-convergence tests */
							static double dt3c = -1;
							if(dt3c < 0){ const char *e_ = getenv("GFS_DT3_COEF"); dt3c = (e_ && e_[0]) ? atof(e_) : 0.1; }
							postype dt3 = (postype)dt3c*dramp_cfl / (2.0*sqrt(uij2));
							if((double)dt3 < gfs_dtb[2]){
#ifdef _OPENMP
#pragma omp critical (gfs_dtb)
#endif
								{ if((double)dt3 < gfs_dtb[2]) gfs_dtb[2] = (double)dt3; }
							}
							if(dt3 < dt) dt = dt3;
						}

						if(isnan(dt)){
							DEBUGPRINT("P%d blend: nan dt %d p%ld xy= %g %g : j=%ld jxy= %g %g dramp= %g vsig= %g cs_i= %g cs_j= %g pi_tot= %g\n",
									MYID(simpar), i, (long)PINDX(p+i), p[i].x, p[i].y,
									(long)PINDX((treevorork4particletype*)jbp), jbp->x, jbp->y,
									dramp, vsig, ibp_csound, jbp_csound, pi_total);
							dt = 1e-10;
						}
						if(dt < 1e-6){
							fprintf(stderr,"[HYDRO_CFL] P%d i=%d dt=%g dramp=%g dramp_cfl=%g vsig=%g csound_i=%g csound_j=%g vol=%g x=%g y=%g uij2=%g\n",
								MYID(simpar), i, dt, dramp, dramp_cfl, vsig, ibp_csound, jbp_csound, ibp_rk4->volume, ibp_rk4->x, ibp_rk4->y, uij2);
						}
						ibp_rk4->dt = MIN(ibp_rk4->dt, dt);
						if(dt < Dtime) Dtime = dt;

						/* Track vsig_max for CD10 (home cell only).
						   The face i-j also appears in j's face list where j
						   is home; by face-symmetry of vsig, j's vsig_max
						   will be set there. Writing only ibp->stress.vsig_max
						   removes the OMP critical-section hot spot. */
						if(av_mode >= 1){
							/* For ghost faces, use sound-speed-only signal velocity
							   (reversed velocity inflates the full vsig artificially) */
							postype vsig_cd = jbp_is_ghost ?
								(ibp_csound + jbp_csound) : vsig;
							if(vsig_cd > ibp->stress.vsig_max) ibp->stress.vsig_max = vsig_cd;
						}
					}
					if(tmp->upperrelated < 0 && gfs_tcf_on() && tcf_dump_fp()){
						fprintf(stderr, "[TCFGEO] cap edge in cell %ld\n", (long)PINDX(ibp_rk4));
					}
					if(tmp->upperrelated < 0 && gfs_tcf_on()){
						/* geometric volume rate of a cap/arc edge moving with the generator (diagnostic) */
						Voro2D_point cl; cl.x = tmp->upperlink->x - tmp->x; cl.y = tmp->upperlink->y - tmp->y;
						Voro2D_point cS = voro2D_norm(&cl);
						tcf_dV += (double)Vec2DLength(tmp, tmp->upperlink)*((double)ibp_vx*(double)cS.x + (double)ibp_vy*(double)cS.y);
					}
					if(tmp->upperrelated < 0 && gfs_schemeH_on() && gfs_dynw2_active){
						/* (H): a cap edge moves rigidly with the generator and
						 * borders vacuum (P* = 0): no force, work -P_i A (v_i.n). */
						Voro2D_point cl; cl.x = tmp->upperlink->x - tmp->x; cl.y = tmp->upperlink->y - tmp->y;
						Voro2D_point cS = voro2D_norm(&cl);   /* unit normal */
						postype cA = Vec2DLength(tmp, tmp->upperlink);
						dte += -ibp_pressure*cA*(ibp_vx*cS.x + ibp_vy*cS.y);
					}
					if(tmp->upperrelated < 0 && gfs_tq_on()){
						double cx, cy; tq_center(&cx, &cy);
						double rx = (double)ibp_rk4->x - cx, ry = (double)ibp_rk4->y - cy;
						int kb = tq_bin(sqrt(rx*rx + ry*ry));
						double len = (double)Vec2DLength(tmp, tmp->upperlink);
#ifdef _OPENMP
#pragma omp atomic
#endif
						tq_cap_n[kb] += 1;
#ifdef _OPENMP
#pragma omp atomic
#endif
						tq_cap_len[kb] += len;
					}
					tmp = tmp->upperlink;
				} while(tmp != vorocorner);
				if(gfs_hg_affine() && (!gfs_hg_split() || hgm_building)){
					sub_poly S_;
					int capfree = 1; { Voro2D_Corner *t_ = vorocorner; do { if(t_->upperrelated < 0 && !(gfs_hg_arc() && t_->upperrelated == VORO_ARC_EDGE)) capfree = 0; t_ = t_->upperlink; } while(t_ != vorocorner); }
					if(gfs_hg_arc()) hg_arc_R = voro_disc_radius(simpar);
					if(capfree && sub_build(vorocorner, neighwork, ibp_rk4, &S_))
						hg_cell(&S_, ibp_rk4, (double)ibp_den, (double)ibp_pressure, (double)Gamma);
					else {
#ifdef _OPENMP
#pragma omp atomic
#endif
						hg_nskip++;
					}
				}
				if(gfs_subzone_on()){
					/* Phi gradients of this cell, scattered to the generator and to the two
					 * other generators of every corner (circumcentre Jacobians) */
					sub_poly S_; double g_, phi_, gP_[2], gvx_[SUB_MAXC], gvy_[SUB_MAXC];
					size_t idi = (size_t)PINDX(ibp_rk4);
					if(sub_build(vorocorner, neighwork, ibp_rk4, &S_)
							&& sub_eval(&S_, (double)ibp_rk4->mass, (double)ibp_pressure, (double)Gamma, &g_, &phi_, gP_, gvx_, gvy_)){
						int n_ = S_.n, k_;
						double dphi = gP_[0]*(double)ibp_vx + gP_[1]*(double)ibp_vy;
						double fix = -gP_[0], fiy = -gP_[1];
						for(k_ = 0; k_ < n_; k_++){
							int kp_ = (k_+n_-1) % n_;
							treevorork4particletype *lo = S_.up[kp_], *up = S_.up[k_];
							double vx_ = S_.ax[k_], vy_ = S_.ay[k_], gx = gvx_[k_], gy = gvy_[k_];
							double Ji00, Ji01, Ji10, Ji11, Jl00 = 0, Jl01 = 0, Jl10 = 0, Jl11 = 0, Ju00 = 0, Ju01 = 0, Ju10 = 0, Ju11 = 0;
							if(lo && up){
								/* corner = circumcentre of (i, lo, up): rows 2(x_lo - x_i), 2(x_up - x_i) */
								double p1x = (double)lo->x - (double)ibp_rk4->x, p1y = (double)lo->y - (double)ibp_rk4->y;
								double p2x = (double)up->x - (double)ibp_rk4->x, p2y = (double)up->y - (double)ibp_rk4->y;
								double a00 = 2*p1x, a01 = 2*p1y, a10 = 2*p2x, a11 = 2*p2y, det = a00*a11 - a01*a10;
								if(!(fabs(det) > 0)) continue;
								double i00 = a11/det, i01 = -a01/det, i10 = -a10/det, i11 = a00/det;
								double ri0x = 2*vx_, ri0y = 2*vy_;
								Ji00 = (i00 + i01)*ri0x; Ji01 = (i00 + i01)*ri0y; Ji10 = (i10 + i11)*ri0x; Ji11 = (i10 + i11)*ri0y;
								double ql_x = 2*(p1x - vx_), ql_y = 2*(p1y - vy_), qu_x = 2*(p2x - vx_), qu_y = 2*(p2y - vy_);
								Jl00 = i00*ql_x; Jl01 = i00*ql_y; Jl10 = i10*ql_x; Jl11 = i10*ql_y;
								Ju00 = i01*qu_x; Ju01 = i01*qu_y; Ju10 = i11*qu_x; Ju11 = i11*qu_y;
							} else if(lo || up){
								/* corner = bisector(i, j) x box side (normal nc, moves with x_i):
								 * rows d^T dv = (v - x_i).dx_i + (x_j - v).dx_j, nc^T dv = nc.dx_i */
								treevorork4particletype *jb = lo ? lo : up;
								int kc = lo ? k_ : kp_;   /* the cap edge at this corner */
								double dx_ = (double)jb->x - (double)ibp_rk4->x, dy_ = (double)jb->y - (double)ibp_rk4->y;
								double ncx = S_.nx[kc], ncy = S_.ny[kc];
								double det = dx_*ncy - dy_*ncx;
								if(!(fabs(det) > 0)) continue;
								double i00 = ncy/det, i01 = -dy_/det, i10 = -ncx/det, i11 = dx_/det;
								/* J_i = Ai [ (v)^T ; nc^T ],  J_j = Ai [ (x_j - v)^T ; 0 ] */
								Ji00 = i00*vx_ + i01*ncx; Ji01 = i00*vy_ + i01*ncy; Ji10 = i10*vx_ + i11*ncx; Ji11 = i10*vy_ + i11*ncy;
								double qx = dx_ - vx_, qy = dy_ - vy_;
								if(lo){ Jl00 = i00*qx; Jl01 = i00*qy; Jl10 = i10*qx; Jl11 = i10*qy; }
								else  { Ju00 = i00*qx; Ju01 = i00*qy; Ju10 = i10*qx; Ju11 = i10*qy; }
							} else {
								/* box corner: moves rigidly with x_i */
								Ji00 = 1; Ji01 = 0; Ji10 = 0; Ji11 = 1;
							}
							fix -= Ji00*gx + Ji10*gy; fiy -= Ji01*gx + Ji11*gy;
							double dvx = Ji00*ibp_vx + Ji01*ibp_vy, dvy = Ji10*ibp_vx + Ji11*ibp_vy;
							if(lo){
								double flx = -(Jl00*gx + Jl10*gy), fly = -(Jl01*gx + Jl11*gy);
								dvx += Jl00*(double)lo->vx + Jl01*(double)lo->vy; dvy += Jl10*(double)lo->vx + Jl11*(double)lo->vy;
								size_t idl = (size_t)PINDX(lo);
								if(idl < sub_cap){
#ifdef _OPENMP
#pragma omp atomic
#endif
									sub_fx[idl] += flx;
#ifdef _OPENMP
#pragma omp atomic
#endif
									sub_fy[idl] += fly;
								}
							}
							if(up){
								double fux = -(Ju00*gx + Ju10*gy), fuy = -(Ju01*gx + Ju11*gy);
								dvx += Ju00*(double)up->vx + Ju01*(double)up->vy; dvy += Ju10*(double)up->vx + Ju11*(double)up->vy;
								size_t idu = (size_t)PINDX(up);
								if(idu < sub_cap){
#ifdef _OPENMP
#pragma omp atomic
#endif
									sub_fx[idu] += fux;
#ifdef _OPENMP
#pragma omp atomic
#endif
									sub_fy[idu] += fuy;
								}
							}
							dphi += gx*dvx + gy*dvy;
						}
						if(idi < sub_cap){
#ifdef _OPENMP
#pragma omp atomic
#endif
							sub_fx[idi] += fix;
#ifdef _OPENMP
#pragma omp atomic
#endif
							sub_fy[idi] += fiy;
							sub_dphi[idi] = dphi; sub_phi[idi] = phi_;
						}
#ifdef _OPENMP
#pragma omp atomic
#endif
						sub_ncell++;
					} else {
						if(idi < sub_cap){ sub_dphi[idi] = 0; sub_phi[idi] = 0; }
#ifdef _OPENMP
#pragma omp atomic
#endif
						sub_ncap++;
					}
				}
				if(do_i_dbg){
#pragma omp critical (i_dbg_pf)
					{ fprintf(stderr,"[FORCE_DBG] r%d ix=%d iy=%d i=%d call=%d post-face fx=%g fy=%g die=%g dKd=%g vol=%g\n", MYID(simpar), ix, iy, i, s_force_call, (double)fx, (double)fy, (double)die, (double)dK_diss, (double)ibp_rk4->volume); fflush(stderr); }
				}

				/* Viscous CFL: dt_visc = 0.5 h^2 / max(nu_phys, nu_cd, chi)
				   Evaluated once per particle (depends only on i's state). */
				if(av_mode >= 1 || nu_phys > 0){
					postype h_i = sqrt(ibp_rk4->volume);
					postype nu_cd_i = ibp->stress.alpha_cd * h_i * ibp_csound;
					postype chi = (nu_phys > 0 && prandtl > 0) ? nu_phys / prandtl : 0;
					postype nu_eff = fmax(fmax(nu_phys, nu_cd_i), chi);
					if(nu_eff > 0){
						postype dt_visc = 0.5 * h_i * h_i / nu_eff;
						if(dt_visc < 1e-6){
							fprintf(stderr,"[VISC_CFL] P%d i=%d dt_visc=%g h_i=%g nu_eff=%g alpha_cd=%g csound=%g den=%g P=%g ie=%g vol=%g x=%g y=%g\n",
								MYID(simpar), i, dt_visc, h_i, nu_eff, ibp->stress.alpha_cd, ibp_csound, ibp_den, ibp_pressure, ibp_rk4->ie, ibp_rk4->volume, ibp_rk4->x, ibp_rk4->y);
						}
						ibp_rk4->dt = MIN(ibp_rk4->dt, dt_visc);
						if(dt_visc < Dtime) Dtime = dt_visc;
					}
				}

				if(gfs_tcf_on() && tcf_sacc_file()){
					size_t id_ = (size_t)PINDX(ibp_rk4);
					if(id_ < tcf_cap) tcf_rate[id_] = tcf_tds/(double)ibp_rk4->ie;
				}
				if(gfs_tcf_on() && tcf_tds < 0){
#ifdef _OPENMP
#pragma omp atomic
#endif
					tcf_cellneg++;
				}
				if(vaud_on_l && (size_t)PINDX(ibp_rk4) < vaud_cap)
					vaud_rate[(size_t)PINDX(ibp_rk4)] = dvb;
				if(sedov_phase1_on()){
					/* die holds dE/dt. The kick rebuilds ie from E − kinetic. */
					ibp_rk4->die = (postype)dte;
					ibp->stress.dK = 0;
				} else {
				ibp_rk4->die = die;
				/* Fix A: per-cell entropy gate — dK_diss must be ≥ 0.
				 * The face-pair NS viscous term Σ_faces τ·dS·(v_f-v_i) is positive
				 * in the well-resolved limit (= q_NS·V ≥ 0), but individual face
				 * terms can be negative, and discretization noise/boundary effects
				 * can drive the cell sum slightly negative.  Clamping enforces the
				 * 2nd law per cell, costing only that small noise (well below the
				 * AV heating signal). */
				if(dK_diss < 0) dK_diss = 0;
				/* Convert dissipative die to entropy rate dK/dt = (γ-1) * q_diss / ρ^γ
				 * where q_diss = dK_diss/V_cell (volumetric rate), so per cell:
				 * dK/dt = (γ-1) * dK_diss / (V_cell * ρ^γ) = (γ-1) * dK_diss / (m * ρ^(γ-1)) */
				if(!isfinite(dK_diss)){
					fprintf(stderr,"[GUARD] r%d ix=%d iy=%d i=%d call=%d dK_diss=%g (NaN/Inf) den=%g P=%g cs=%g K=%g\n",
						MYID(simpar), ix, iy, i, s_force_call, (double)dK_diss,
						(double)ibp_den, (double)ibp_pressure, (double)ibp_csound, (double)ibp->stress.K);
					fflush(stderr);
					dK_diss = 0;
				}
				if(ibp_rk4->mass > 0 && ibp_den > 0){
					postype rho_pow = pow((double)ibp_den, (double)(Gamma-1));
					ibp->stress.dK = (Gamma-1) * dK_diss / (ibp_rk4->mass * rho_pow);
					if(!isfinite(ibp->stress.dK)){
						fprintf(stderr,"[GUARD] r%d ix=%d iy=%d i=%d call=%d dK=%g (NaN/Inf) rho_pow=%g mass=%g dK_diss=%g\n",
							MYID(simpar), ix, iy, i, s_force_call, (double)ibp->stress.dK,
							(double)rho_pow, (double)ibp_rk4->mass, (double)dK_diss);
						fflush(stderr);
						ibp->stress.dK = 0;
					}
				} else {
					ibp->stress.dK = 0;
				}
				/* dK-aware CFL: dt ≤ DKCFL_FRAC · K / |dK|.
				 * Mode 2 (Sedov shock-front, P-ratio 5700:1 face) produced
				 * |dK|≈1e+9 vs K≈1450 → ΔK/K~40% in one substep. Shrinking dt
				 * preserves conservation; no direct dK clipping.
				 * K derived locally from P/ρ^γ so the constraint works in
				 * both entropy_mode=0 (ie-mode) and entropy_mode=1. */
				{
					postype dK_abs = fabs(ibp->stress.dK);
					if(dK_abs > 0 && ibp_den > 0){
						postype rho_pow_K = pow((double)fmax(ibp_den,(postype)1e-30), (double)Gamma);
						postype K_eff = ibp_pressure / rho_pow_K;
						if(K_eff < DKCFL_KFLOOR) K_eff = DKCFL_KFLOOR;
						postype dt_dK = DKCFL_FRAC * K_eff / dK_abs;
						if(dt_dK < 1e-9){
							fprintf(stderr,"[DKCFL] P%d i=%d dt_dK=%g K_eff=%g dK=%g den=%g P=%g cs=%g x=%g y=%g\n",
								MYID(simpar), i, (double)dt_dK, (double)K_eff, (double)ibp->stress.dK,
								(double)ibp_den, (double)ibp_pressure, (double)ibp_csound,
								(double)ibp_rk4->x, (double)ibp_rk4->y);
						}
						ibp_rk4->dt = MIN(ibp_rk4->dt, dt_dK);
						if(dt_dK < Dtime) Dtime = dt_dK;
					}
				}
				} /* end !sedov_phase1 dK path */
				if(do_i_dbg){
#pragma omp critical (i_dbg_pdk)
					{ fprintf(stderr,"[FORCE_DBG] r%d ix=%d iy=%d i=%d call=%d post-dK dK=%g mass=%g vol=%g\n", MYID(simpar), ix, iy, i, s_force_call, (double)ibp->stress.dK, (double)ibp_rk4->mass, (double)ibp_rk4->volume); fflush(stderr); }
				}
				ibp_rk4->ax = fx/ibp_rk4->mass;
				ibp_rk4->ay = fy/ibp_rk4->mass;
				if(tcf_dump_fp()){
#pragma omp critical (tcf_dump)
					fprintf(tcf_dump_fp(), "%ld %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.6e %.6e %.6e %.6e %.6e %.6e %.6e %d %.6e %.6e %.6e %.6e\n",
						(long)PINDX(ibp_rk4), (double)ibp_rk4->x, (double)ibp_rk4->y, (double)ibp_vx, (double)ibp_vy,
						(double)ibp_pressure, (double)ibp_rk4->volume, (double)ibp_rk4->mass, fx, fy, dte, tcf_tds, tcf_dV,
						tg_nx, tg_ny, tg_ex, tg_ey, tg_bis, tg_en, tg_amin, tg_nf,
						(double)ibp->stress.gUxx, (double)ibp->stress.gUxy, (double)ibp->stress.gUyx, (double)ibp->stress.gUyy);
				}

				/* Hyperviscosity force: dv/dt -= ν₄ × biharm_v_i
				 * ν₄ = α₄ × cs × Δx³, with Δx = sqrt(V_i) (cell-local scale).
				 * Selectively damps high-k modes (rate ∝ ν₄ k⁴).
				 * Optional cap: |a_HV| ≤ cap·cs²/Δx (CFL acceleration) — prevents
				 * lap_v runaway from feeding back into a destabilizing force. */
				{
					postype hyperv_alpha = GAS_HYPERVALPHA(simpar);
					if(hyperv_alpha > 0 && ibp_rk4->volume > 0){
						postype invVol_hv = 1.0/ibp_rk4->volume;
						postype dx_eff = sqrt(ibp_rk4->volume);
						postype nu4 = hyperv_alpha * ibp_csound * dx_eff*dx_eff*dx_eff;
						postype a_hv_x = -nu4 * biharm_vx_acc * invVol_hv;
						postype a_hv_y = -nu4 * biharm_vy_acc * invVol_hv;
						postype hv_cap = GAS_HYPERVFORCECAP(simpar);
						if(hv_cap > 0){
							postype a_hv_mag = sqrt(a_hv_x*a_hv_x + a_hv_y*a_hv_y);
							postype a_hv_max = hv_cap * ibp_csound * ibp_csound / dx_eff;
							if(a_hv_mag > a_hv_max && a_hv_mag > 0){
								postype s = a_hv_max / a_hv_mag;
								a_hv_x *= s; a_hv_y *= s;
							}
						}
						ibp_rk4->ax += a_hv_x;
						ibp_rk4->ay += a_hv_y;
					}
				}

				/* Acceleration CFL: dt ≤ ACC_CFL_FRAC · sqrt(h_i / |a_i|).
				 * Mode 2 (Sedov shock-front): ∇P/ρ jumps ~4 orders across an
				 * extreme-P-ratio face → |a|~1e+10 → free-fall time
				 * τ_a = sqrt(h/|a|) collapses to ns. This CFL keeps the cell
				 * from traversing more than ACC_CFL_FRAC² · h_i per substep,
				 * preventing position+velocity overshoot that the K-based CFL
				 * cannot anticipate (since K-CFL acts after-the-fact on the
				 * resulting dK). Preserves conservation. */
				{
					postype ax_i = ibp_rk4->ax;
					postype ay_i = ibp_rk4->ay;
					postype amag2 = ax_i*ax_i + ay_i*ay_i;
					if(amag2 > 0 && ibp_rk4->volume > 0){
						postype amag = sqrt(amag2);
						postype h_i = sqrt(ibp_rk4->volume);
						postype dt_acc = ACC_CFL_FRAC * sqrt(h_i / amag);
						if(dt_acc < 1e-9){
							fprintf(stderr,"[ACC_CFL] P%d i=%d dt_acc=%g |a|=%g h_i=%g x=%g y=%g den=%g P=%g\n",
								MYID(simpar), i, (double)dt_acc, (double)amag, (double)h_i,
								(double)ibp_rk4->x, (double)ibp_rk4->y,
								(double)ibp_den, (double)ibp_pressure);
						}
						ibp_rk4->dt = MIN(ibp_rk4->dt, dt_acc);
						if(dt_acc < Dtime) Dtime = dt_acc;
					}
				}

				/* Store face-averaged velocity offset for XSPH position drift.
				 * Used by integrator only when GAS_XSPHEPS > 0. */
				if(xsph_w_acc > 0){
					ibp->stress.v_smooth_x = xsph_vx_acc / xsph_w_acc;
					ibp->stress.v_smooth_y = xsph_vy_acc / xsph_w_acc;
				} else {
					ibp->stress.v_smooth_x = 0;
					ibp->stress.v_smooth_y = 0;
				}

				if(isnan(fx)){
					DEBUGPRINT("P%d blend: nan fx %d : %g %g %ld p= %g\n",
							MYID(simpar), i, ibp->x, ibp->y, PINDX(ibp), ibp_pressure);
					if(gfs_dynw2_active){ fprintf(stderr,"[DYNW2] nan fx: den=%g vol=%g w2=%g w2old=%g cs=%g dtold=%g fy=%g die=%g\n",
						(double)ibp_den,(double)ibp_rk4->volume,(double)ibp_rk4->w2,(double)ibp_rk4->w2old,(double)ibp_csound,(double)dtold,(double)fy,(double)die); fflush(stderr); }
					exit(0);
				}
				if(do_i_dbg){
#pragma omp critical (i_dbg_out)
					{ fprintf(stderr,"[FORCE_DBG] r%d ix=%d iy=%d i=%d call=%d exit fx=%g fy=%g die=%g dK_diss=%g\n",
						MYID(simpar), ix, iy, i, s_force_call, (double)fx, (double)fy, (double)die, (double)dK_diss); fflush(stderr); }
				}
				/* Adaptive expansion cleanup: free per-particle expansion buffers (if used)
				 * and restore cell-level pointers for the next i-iteration. */
				if(neighbors_e){
					free(neighbors_e); free(neighwork_e);
					neighbors = neighbors_save;
					neighwork = neighwork_save;
					nneigh = nneigh_save;
				}
			}
			free(p); free(neighbors); free(neighwork);
		}
		free(vorocorner);
		if(do_iy_dbg){
#pragma omp critical (iy_dbg_out)
			{ fprintf(stderr,"[FORCE_DBG] r%d iy=%d out call=%d\n", MYID(simpar), iy, s_force_call); fflush(stderr); }
		}
	}
	if(gfs_faud_on()) faud_match(s_force_call);
	if(MYID(simpar)<=1){ fprintf(stderr,"[FORCE_DBG] r%d delta impl exit Dtime=%g call=%d\n", MYID(simpar), (double)Dtime, s_force_call); fflush(stderr); }
	return Dtime;
}

/* ================================================================
 *  getAccVoro2DBlend: Original entry point (calls _impl).
 * ================================================================ */
/* GFS_DYN_W2 (see gfs_dynw2_on): face collection, PCG solve, diagnostics. */
typedef struct { int i, j, ext; double a, b, d, dd, dW, vstar_rel; } dyn_face_t;
static dyn_face_t *dyn_f = NULL; static size_t dyn_fn = 0, dyn_fcap = 0;
/* step accumulators (sums over the 4 stages) */
static int dyn_stage_cur = 0; static double dyn_dtlim_last = 1e30;
static double dyn_r2e, dyn_b2e, dyn_r2o, dyn_b2o, dyn_dv2e, dyn_vs2e, dyn_dv2o;
static long dyn_ne, dyn_no, dyn_nbad, dyn_nbar; static int dyn_itmax; static double dyn_wratio_max;
static void dyn_step_reset(void){
	dyn_r2e = dyn_b2e = dyn_r2o = dyn_b2o = dyn_dv2e = dyn_vs2e = dyn_dv2o = 0;
	dyn_ne = dyn_no = dyn_nbad = dyn_nbar = 0; dyn_itmax = 0; dyn_wratio_max = 0;
	dyn_fchk_s2 = 0; dyn_fchk_n = 0;
}
static void dyn_face_push(dyn_face_t *buf, size_t n){
	if(dyn_fn + n > dyn_fcap){
		dyn_fcap = (dyn_fn + n)*2 + 1024;
		dyn_f = (dyn_face_t*)realloc(dyn_f, dyn_fcap*sizeof(dyn_face_t));
	}
	memcpy(dyn_f + dyn_fn, buf, n*sizeof(dyn_face_t));
	dyn_fn += n;
}
/* Returns the Dt limit from cell validity (1e30 if none). */
static double gfs_dynw2_solve(SimParameters *simpar, int stage, postype Gamma,
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *)){
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype boxsize = VORO_CELLBOX(simpar);
	int mx = BASICCELL_MX(simpar), my = BASICCELL_MY(simpar);
	double ratio = gfs_env_d("GFS_DYN_W2_RATIO", 100.0);
	double px = gfs_env_d("GFS_DYN_W2_PX", 100.0);
	int mid = (int)gfs_env_d("GFS_DYN_W2_MID", 0);
	double delta = gfs_env_d("GFS_DYN_W2_DELTA", 0.2);
	double tol = gfs_env_d("GFS_DYN_W2_TOL", 1e-10);
	int maxit = (int)gfs_env_d("GFS_DYN_W2_MAXIT", 2000);
	double vbw = gfs_env_d("GFS_DYN_W2_BARRIER", 1000.0);
	double tau = (double)GAS_dtold(simpar);
	int i, iy;
	size_t max_id = 0;
	for(i=0;i<nbp;i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(bp+i)); if(id > max_id) max_id = id; }
	dyn_ensure(max_id);
	int *imap = (int*)malloc((max_id+1)*sizeof(int));
	for(size_t k=0;k<=max_id;k++) imap[k] = -1;
	for(i=0;i<nbp;i++) imap[(size_t)PINDX((treevorork4particletype*)(bp+i))] = i;
	dyn_fn = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if((int)gfs_env_d("GFS_DYN_W2_OMP", 0))
#endif
	for(iy=1;iy<my-1;iy++){
		int mp = 4000, ix;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		size_t nb = 0, cb = 256;
		dyn_face_t *buf = (dyn_face_t*)malloc(cb*sizeof(dyn_face_t));
		for(ix=1;ix<mx-1;ix++){
			int np, nneigh, k;
			treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
			Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
			Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			for(k=0;k<np;k++){
				treevorork4particletype *ib = p[k].bp;
				treevorostressrk4particletype *ibs = (treevorostressrk4particletype*)ib;
				Voro2D_point center;
				center.x = p[k].x; center.y = p[k].y;
				center.indx = PINDX(p+k); center.csound = p[k].csound; center.w2 = p[k].w2;
				Voro2D_point *nb_use = neighbors, *nw_use = neighwork, *ne = NULL, *nwe = NULL;
				int nn_use = nneigh;
				int ip = Voro2D_FindVC(&center, nb_use, nw_use, nn_use, vorocorner, mp, boxsize);
				if(Voro2D_HasUnassignedCorner(vorocorner, ip) || voro_sec_short(simpar, vorocorner, ip, (double)center.x, (double)center.y, ix, iy, 1)){
					int R;
					for(R=2;R<=voro_expand_rmax(simpar);R++){
						int nne;
						if(ne){ free(ne); free(nwe); }
						ne = searchCellRk4Neighbors2D_R(simpar, ix, iy, R, &nne);
						nwe = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nne);
						ip = Voro2D_FindVC(&center, ne, nwe, nne, vorocorner, mp, boxsize);
						nb_use = ne; nw_use = nwe; nn_use = nne;
						if(!Voro2D_HasUnassignedCorner(vorocorner, ip) && !voro_sec_short(simpar, vorocorner, ip, (double)center.x, (double)center.y, ix, iy, R)) break;
					}
				}
				size_t idi = (size_t)PINDX(ib);
				Voro2D_Corner *tmp = vorocorner;
				do {
					if(tmp->upperrelated >= 0 && tmp->upperrelated < nn_use && nw_use[tmp->upperrelated].bp){
						treevorork4particletype *jb = (treevorork4particletype*)nw_use[tmp->upperrelated].bp;
						treevorostressrk4particletype *jbs = (treevorostressrk4particletype*)jb;
						if(PINDX(jb) != MAX_INDEX && (size_t)PINDX(jb) > idi && (size_t)PINDX(jb) <= max_id
								&& imap[(size_t)PINDX(jb)] >= 0){
							Voro2D_Corner *tmp2 = tmp->upperlink;
							double area = Vec2DLength(tmp,tmp2);
							double rx = jb->x - ib->x, ry = jb->y - ib->y;
							double d = sqrt(rx*rx + ry*ry);
							if(d > 0 && area > 0){
								double ex = rx/d, ey = ry/d;
								double vni = ib->vx*ex + ib->vy*ey, vnj = jb->vx*ex + jb->vy*ey;
								double dd = vnj - vni, nub = 0.5*(vni + vnj);
								double dW = (double)ib->w2 - (double)jb->w2;
								Voro2D_point ur; ur.x = 0; ur.y = 0;
								voronoi_face_rotation(&ur, (postype)0.5*(tmp->x + tmp2->x), (postype)0.5*(tmp->y + tmp2->y),
									ib->x, ib->y, jb->x, jb->y, ib->vx, ib->vy, jb->vx, jb->vy,
									(postype)ex, (postype)ey, (postype)d, ib->w2, jb->w2, 1);
								double rot = ur.x*ex + ur.y*ey;
								postype ps, vns;
								hll_star_state(ibs->den, ibs->pressure, (postype)vni, ibs->csound,
									jbs->den, jbs->pressure, (postype)vnj, jbs->csound, Gamma, &ps, &vns);
								double pmn = fmin(ibs->pressure, jbs->pressure), pmxv = fmax(ibs->pressure, jbs->pressure);
								int ext = (pmxv > px*fmax(pmn, 1e-30));
								double vt = (mid && !ext) ? (nub + rot) : (double)vns;
								dyn_face_t f;
								f.i = imap[idi]; f.j = imap[(size_t)PINDX(jb)]; f.ext = ext;
								f.a = area*(ext ? ratio : 1.0);
								f.b = 2*d*(vt - nub - rot) + dW*dd/d;
								/* Validity barrier: past (1-2 delta) d^2 the target keeps
								 * dW/d^2 fixed (dW' = 2 dW ddot/d) and relaxes |dW| back over
								 * tau = dt_old, with weight a_f*vbw. Counted in [DYNW2] nbar. */
								if(vbw > 0 && fabs(dW) > (1-2*delta)*d*d && tau > 0){
									double ex_ = fabs(dW) - (1-2*delta)*d*d;
									f.b = 2*dW*dd/d - (dW > 0 ? 1 : -1)*ex_/tau;
									f.a = area*vbw*(ext ? ratio : 1.0);
									f.ext = ext | 2;
								}
								f.d = d; f.dd = dd; f.dW = dW;
								f.vstar_rel = (double)vns - nub - rot;
								if(nb == cb){ cb *= 2; buf = (dyn_face_t*)realloc(buf, cb*sizeof(dyn_face_t)); }
								buf[nb++] = f;
							}
						}
					}
					tmp = tmp->upperlink;
				} while(tmp != vorocorner);
				if(ne){ free(ne); free(nwe); }
			}
			free(p); free(neighbors); free(neighwork);
		}
#ifdef _OPENMP
#pragma omp critical (dynw2_push)
#endif
		dyn_face_push(buf, nb);
		free(buf); free(vorocorner);
	}
	/* PCG on L_a x = B^T A b */
	size_t f; int n = nbp;
	double *x = (double*)calloc(n, sizeof(double)), *g = (double*)calloc(n, sizeof(double));
	double *r = (double*)malloc(n*sizeof(double)), *z = (double*)malloc(n*sizeof(double));
	double *pp = (double*)malloc(n*sizeof(double)), *q = (double*)malloc(n*sizeof(double));
	double *dg = (double*)calloc(n, sizeof(double));
	for(f=0;f<dyn_fn;f++){
		dyn_face_t *F = dyn_f+f;
		g[F->i] += F->a*F->b; g[F->j] -= F->a*F->b;
		dg[F->i] += F->a; dg[F->j] += F->a;
	}
	for(i=0;i<n;i++){ x[i] = dyn_wd[(size_t)PINDX((treevorork4particletype*)(bp+i))]; if(!(dg[i] > 0)) dg[i] = 1; }
	{ double m = 0; for(i=0;i<n;i++) m += x[i]; m /= n; for(i=0;i<n;i++) x[i] -= m; }
#define DYN_APPLY(v, out) do{ for(i=0;i<n;i++) out[i] = 0; \
	for(f=0;f<dyn_fn;f++){ dyn_face_t *F = dyn_f+f; double t = F->a*(v[F->i]-v[F->j]); out[F->i] += t; out[F->j] -= t; } }while(0)
	DYN_APPLY(x, q);
	double gn = 0, rz = 0, rn = 0;
	for(i=0;i<n;i++){ r[i] = g[i] - q[i]; gn += g[i]*g[i]; }
	gn = sqrt(gn);
	for(i=0;i<n;i++){ z[i] = r[i]/dg[i]; pp[i] = z[i]; rz += r[i]*z[i]; rn += r[i]*r[i]; }
	int it = 0;
	while(gn > 0 && sqrt(rn) > tol*gn && it < maxit){
		DYN_APPLY(pp, q);
		double pq = 0; for(i=0;i<n;i++) pq += pp[i]*q[i];
		if(!(pq > 0)) break;
		double al = rz/pq, rz2 = 0; rn = 0;
		for(i=0;i<n;i++){ x[i] += al*pp[i]; r[i] -= al*q[i]; z[i] = r[i]/dg[i]; rz2 += r[i]*z[i]; rn += r[i]*r[i]; }
		double be = rz2/rz; rz = rz2;
		for(i=0;i<n;i++) pp[i] = z[i] + be*pp[i];
		it++;
	}
#undef DYN_APPLY
	{ double m = 0; for(i=0;i<n;i++) m += x[i]; m /= n; for(i=0;i<n;i++) x[i] -= m; }
	if(it > dyn_itmax) dyn_itmax = it;
	for(i=0;i<n;i++){
		size_t id = (size_t)PINDX((treevorork4particletype*)(bp+i));
		if(!isfinite(x[i])) x[i] = 0;
		dyn_wd[id] = x[i]; dyn_k[stage][id] = x[i];
	}
	/* residuals, validity */
	double dtlim = 1e30;
	for(f=0;f<dyn_fn;f++){
		dyn_face_t *F = dyn_f+f;
		double dx = x[F->i] - x[F->j], rr = dx - F->b;
		if(F->ext & 2) dyn_nbar++;
		if(F->ext & 1){
			dyn_r2e += rr*rr; dyn_b2e += F->b*F->b; dyn_ne++;
			dyn_dv2e += rr*rr/(4*F->d*F->d);
			dyn_vs2e += F->vstar_rel*F->vstar_rel;
		} else { dyn_r2o += rr*rr; dyn_b2o += F->b*F->b; dyn_no++; dyn_dv2o += rr*rr/(4*F->d*F->d); }
		double d2 = F->d*F->d, lim = (1-delta)*d2;
		double wr = fabs(F->dW)/d2; if(wr > dyn_wratio_max) dyn_wratio_max = wr;
		double m = lim - fabs(F->dW);
		/* d|dW|/dt - (1-delta) d(d^2)/dt */
		double grow = (F->dW > 0 ? dx : (F->dW < 0 ? -dx : fabs(dx))) - (1-delta)*2*F->d*F->dd;
		if(m <= 0) dyn_nbad++;
		else if(grow > 1e-30 && m/grow < dtlim) dtlim = m/grow;
	}
	free(x); free(g); free(r); free(z); free(pp); free(q); free(dg); free(imap);
	return dtlim;
}
/* GFS_SCHEME_H: per-particle target volume V^T (an RK4 state like W). */
static double *hq_vt0 = NULL, *hq_vts = NULL, *hq_kvt[4] = {NULL,NULL,NULL,NULL};
static char *hq_init = NULL; static size_t hq_cap = 0;
static double hq_rmax, hq_vrel2, hq_vtrel2, hq_proj; static long hq_nbad, hq_ncell, hq_ncomp;
static void hq_ensure(size_t max_id){
	size_t k, nc; int s_;
	if(max_id < hq_cap) return;
	nc = max_id + 1 + max_id/4;
	hq_vt0 = (double*)realloc(hq_vt0, nc*sizeof(double));
	hq_vts = (double*)realloc(hq_vts, nc*sizeof(double));
	hq_init = (char*)realloc(hq_init, nc);
	for(s_=0;s_<4;s_++) hq_kvt[s_] = (double*)realloc(hq_kvt[s_], nc*sizeof(double));
	for(k=hq_cap;k<nc;k++){ hq_vt0[k] = hq_vts[k] = 0; hq_init[k] = 0; for(s_=0;s_<4;s_++) hq_kvt[s_][k] = 0; }
	hq_cap = nc;
}
static int hq_find(int *par, int a){ while(par[a] != a){ par[a] = par[par[a]]; a = par[a]; } return a; }
static double hq_polyarea(Voro2D_Corner *vc){
	return (double)Area2DPolygon(vc, 1);   /* includes arc segments */
}
/* GFS_H_OMEGA=1 (with GFS_SCHEME_H): W_i = h_i^2 omega_i with |omega_i| <=
 * omega_max < 1/2 and the smooth lower bound (Docs/consult_scheme_H_round4_fable.md 5a)
 *   h_i = (R_1^-p + sum_{d_ij < R_c} phi(d_ij/R_c) d_ij^-p)^(-1/p) <= min_j d_ij,
 * so |W_i - W_j| <= 2 omega_max d_ij^2 < d_ij^2 for every pair: every
 * generator stays inside its own Laguerre cell. omega is the RK4 state; its
 * rate solves the W-law least squares in the tangent cone of the box. */
static double *hq_om0 = NULL, *hq_oms = NULL, *hq_kom[4] = {NULL,NULL,NULL,NULL};
static double *hq_h = NULL, *hq_hd = NULL; static size_t hq_ocap = 0;
static double hq_Lx = 0, hq_Ly = 0, hq_R1 = 0;
static long hq_nsat; static double hq_ommax_seen;
static int gfs_h_omega_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_H_OMEGA"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static void hq_oensure(size_t max_id){
	size_t k, nc; int s_;
	if(max_id < hq_ocap) return;
	nc = max_id + 1 + max_id/4;
	hq_om0 = (double*)realloc(hq_om0, nc*sizeof(double));
	hq_oms = (double*)realloc(hq_oms, nc*sizeof(double));
	hq_h = (double*)realloc(hq_h, nc*sizeof(double));
	hq_hd = (double*)realloc(hq_hd, nc*sizeof(double));
	for(s_=0;s_<4;s_++) hq_kom[s_] = (double*)realloc(hq_kom[s_], nc*sizeof(double));
	for(k=hq_ocap;k<nc;k++){ hq_om0[k] = hq_oms[k] = 0; hq_h[k] = hq_hd[k] = 0; for(s_=0;s_<4;s_++) hq_kom[s_][k] = 0; }
	hq_ocap = nc;
}
/* h_i and dh_i/dt for all local particles (periodic box hq_Lx x hq_Ly,
 * 1 rank), cubic-Hermite cutoff between R_1 and R_c, power p. */
static void hq_compute_h(treevorork4particletype *bp0, size_t psz, int np){
	double p = gfs_env_d("GFS_H_P", 8.0), R1 = hq_R1, Rc = gfs_env_d("GFS_H_RC", 2.5)*hq_R1;
	double s1 = R1/Rc;
	int nx = (int)floor(hq_Lx/Rc), ny = (int)floor(hq_Ly/Rc), i;
	if(nx < 1) nx = 1;
	if(ny < 1) ny = 1;
	int *head = (int*)malloc((size_t)nx*ny*sizeof(int)), *next = (int*)malloc(np*sizeof(int));
	for(i=0;i<nx*ny;i++) head[i] = -1;
#define HQ_P(k) ((treevorork4particletype*)((char*)bp0 + (size_t)(k)*psz))
	for(i=0;i<np;i++){
		double x = fmod((double)HQ_P(i)->x, hq_Lx), y = fmod((double)HQ_P(i)->y, hq_Ly);
		if(x < 0) x += hq_Lx;
		if(y < 0) y += hq_Ly;
		int cx = (int)(x/hq_Lx*nx), cy = (int)(y/hq_Ly*ny);
		if(cx >= nx) cx = nx-1;
		if(cy >= ny) cy = ny-1;
		next[i] = head[cx + nx*cy]; head[cx + nx*cy] = i;
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<np;i++){
		treevorork4particletype *a = HQ_P(i);
		double x = fmod((double)a->x, hq_Lx), y = fmod((double)a->y, hq_Ly);
		if(x < 0) x += hq_Lx;
		if(y < 0) y += hq_Ly;
		int cx = (int)(x/hq_Lx*nx), cy = (int)(y/hq_Ly*ny), ox, oy;
		if(cx >= nx) cx = nx-1;
		if(cy >= ny) cy = ny-1;
		double S = pow(R1, -p), Sd = 0;
		for(oy=-1;oy<=1;oy++) for(ox=-1;ox<=1;ox++){
			int qx = (cx+ox+nx)%nx, qy = (cy+oy+ny)%ny, k;
			if((nx < 3 && ox != 0) || (ny < 3 && oy != 0)) continue;
			for(k=head[qx+nx*qy];k>=0;k=next[k]){
				if(k == i) continue;
				treevorork4particletype *b = HQ_P(k);
				double rx = (double)b->x - (double)a->x, ry = (double)b->y - (double)a->y;
				rx -= hq_Lx*floor(rx/hq_Lx + 0.5); ry -= hq_Ly*floor(ry/hq_Ly + 0.5);
				double d = sqrt(rx*rx + ry*ry);
				if(!(d > 0) || d >= Rc) continue;
				double sd = d/Rc, ph, php;
				if(sd <= s1){ ph = 1; php = 0; }
				else { double u = (sd - s1)/(1 - s1); ph = 1 - 3*u*u + 2*u*u*u; php = (-6*u + 6*u*u)/(1 - s1); }
				double dd = (rx*((double)b->vx - (double)a->vx) + ry*((double)b->vy - (double)a->vy))/d;
				S += ph*pow(d, -p);
				Sd += (ph*pow(d, -p-1) - php*pow(d, -p)/(p*Rc))*dd;
			}
		}
		double h = pow(S, -1.0/p);
		size_t id = (size_t)PINDX(a);
		if(id < hq_ocap){ hq_h[id] = h; hq_hd[id] = pow(h, p+1)*Sd; }
	}
#undef HQ_P
	free(head); free(next);
}
/* GFS_H_DAE=1 (with GFS_H_OMEGA): W is algebraic. At the start of every RK
 * stage, before densities are computed, W = h^2 omega is solved from
 * V(x, W) = V^T inside the omega box (box-constrained least squares,
 * projected Newton with backtracking; Jacobian dV/dW = 1/2 L_{A/d}).
 * The thermodynamic volume stays the geometric V. A cell whose target is
 * out of reach (omega on its bound, residual pointing outward) gets
 * V^T := V at the start of a step. No tau_C feedback, no omega step limit.
 * Docs/consult_schemeH_omega_stiff_reply_astra.md (a). */
static int gfs_h_dae_on(void){
	static int on = -1;
	if(on < 0){ const char *s = getenv("GFS_H_DAE"); on = (s && s[0] && atoi(s) > 0) ? 1 : 0; }
	return on;
}
static int hq_dae_it, hq_dae_nsat, hq_dae_nclip; static double hq_dae_res0, hq_dae_res, hq_dae_ommax_pred;
/* Volumes and A/(2d) face graph of the current (x, W) on a freshly built
 * linked list; faces go to dyn_f (i < j by particle index map). */
static void hq_geom(SimParameters *simpar, postype xmin, postype ymin, postype xmax, postype ymax,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		int *imap, size_t max_id, double *vgeo){
	int nbp = VORO_NP(simpar), iy;
	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	int mx, my;
	BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
	BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*mx*my);
	mkLinkedList2D(simpar, cellsize, xmin, ymin, xmax, ymax, paddingAllTreeParticles);
	postype boxsize = VORO_CELLBOX(simpar);
	int i;
	for(i=0;i<nbp;i++) vgeo[i] = 0;
	dyn_fn = 0;
	for(iy=1;iy<my-1;iy++){
		int mp = 4000, ix;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		size_t nb = 0, cb = 256;
		dyn_face_t *buf = (dyn_face_t*)malloc(cb*sizeof(dyn_face_t));
		for(ix=1;ix<mx-1;ix++){
			int np, nneigh, k;
			treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
			Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
			Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			for(k=0;k<np;k++){
				treevorork4particletype *ib = p[k].bp;
				Voro2D_point center;
				center.x = p[k].x; center.y = p[k].y;
				center.indx = PINDX(p+k); center.csound = p[k].csound; center.w2 = p[k].w2;
				Voro2D_point *nw_use = neighwork, *ne = NULL, *nwe = NULL;
				int nn_use = nneigh;
				int ip = Voro2D_FindVC(&center, neighbors, nw_use, nn_use, vorocorner, mp, boxsize);
				if(Voro2D_HasUnassignedCorner(vorocorner, ip) || voro_sec_short(simpar, vorocorner, ip, (double)center.x, (double)center.y, ix, iy, 1)){
					int R;
					for(R=2;R<=voro_expand_rmax(simpar);R++){
						int nne;
						if(ne){ free(ne); free(nwe); }
						ne = searchCellRk4Neighbors2D_R(simpar, ix, iy, R, &nne);
						nwe = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nne);
						ip = Voro2D_FindVC(&center, ne, nwe, nne, vorocorner, mp, boxsize);
						nw_use = nwe; nn_use = nne;
						if(!Voro2D_HasUnassignedCorner(vorocorner, ip) && !voro_sec_short(simpar, vorocorner, ip, (double)center.x, (double)center.y, ix, iy, R)) break;
					}
				}
				size_t idi = (size_t)PINDX(ib);
				int ii = (idi <= max_id) ? imap[idi] : -1;
				if(ii >= 0){
					vgeo[ii] = hq_polyarea(vorocorner);
					Voro2D_Corner *tmp = vorocorner;
					do {
						if(tmp->upperrelated >= 0 && tmp->upperrelated < nn_use && nw_use[tmp->upperrelated].bp){
							treevorork4particletype *jb = (treevorork4particletype*)nw_use[tmp->upperrelated].bp;
							if(PINDX(jb) != MAX_INDEX && (size_t)PINDX(jb) > idi && (size_t)PINDX(jb) <= max_id
									&& imap[(size_t)PINDX(jb)] >= 0){
								double area = Vec2DLength(tmp, tmp->upperlink);
								double rx = jb->x - ib->x, ry = jb->y - ib->y, d = sqrt(rx*rx + ry*ry);
								if(d > 0 && area > 0){
									dyn_face_t f;
									memset(&f, 0, sizeof(f));
									f.i = ii; f.j = imap[(size_t)PINDX(jb)]; f.a = area/(2*d); f.d = d;
									if(nb == cb){ cb *= 2; buf = (dyn_face_t*)realloc(buf, cb*sizeof(dyn_face_t)); }
									buf[nb++] = f;
								}
							}
						}
						tmp = tmp->upperlink;
					} while(tmp != vorocorner);
				}
				if(ne){ free(ne); free(nwe); }
			}
			free(p); free(neighbors); free(neighwork);
		}
		dyn_face_push(buf, nb);
		free(buf); free(vorocorner);
	}
	my_free(VORO_BASICCELL(simpar));
	VORO_BASICCELL(simpar) = NULL;
	my_free(VORORK4_TBPP(simpar));
	VORORK4_TBPP(simpar) = NULL;
	VORO_NPAD(simpar) = 0;
}
static void hq_setW(treevorostressrk4particletype *bp, int n, const double *om){
	int i; double wmin = 1e300;
	for(i=0;i<n;i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(bp+i)); double W = hq_h[id]*hq_h[id]*om[i]; if(W < wmin) wmin = W; }
	for(i=0;i<n;i++){
		size_t id = (size_t)PINDX((treevorork4particletype*)(bp+i));
		double W = hq_h[id]*hq_h[id]*om[i] - wmin;
		if(!(W > 0)) W = 0;
		bp[i].w2 = (postype)W; bp[i].w2old = bp[i].w2;
		hq_oms[id] = om[i];
	}
}
static void hq_project(SimParameters *simpar, postype xmin, postype ymin, postype xmax, postype ymax,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		int step_start){
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int n = VORO_NP(simpar), i, it;
	size_t f, max_id = 0;
	double om_max = gfs_env_d("GFS_H_OMAX", 0.4), tol = gfs_env_d("GFS_H_DAE_TOL", 1e-9);
	int maxit = (int)gfs_env_d("GFS_H_DAE_MAXIT", 8);
	for(i=0;i<n;i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(bp+i)); if(id > max_id) max_id = id; }
	dyn_ensure(max_id); hq_ensure(max_id); hq_oensure(max_id);
	hq_Lx = (double)(SIMBOX(simpar).x.max - SIMBOX(simpar).x.min);
	hq_Ly = (double)(SIMBOX(simpar).y.max - SIMBOX(simpar).y.min);
	hq_R1 = gfs_env_d("GFS_H_R1", 1.0)*hq_Lx/NX(simpar);
	hq_compute_h((treevorork4particletype*)bp, sizeof(treevorostressrk4particletype), n);
	int *imap = (int*)malloc((max_id+1)*sizeof(int));
	for(size_t k=0;k<=max_id;k++) imap[k] = -1;
	for(i=0;i<n;i++) imap[(size_t)PINDX((treevorork4particletype*)(bp+i))] = i;
	double *V = (double*)malloc(n*sizeof(double)), *r = (double*)malloc(n*sizeof(double));
	double *om = (double*)malloc(n*sizeof(double)), *om_old = (double*)malloc(n*sizeof(double));
	double *h2 = (double*)malloc(n*sizeof(double)), *y = (double*)malloc(n*sizeof(double));
	double *rr = (double*)malloc(n*sizeof(double)), *z = (double*)malloc(n*sizeof(double));
	double *pp = (double*)malloc(n*sizeof(double)), *q = (double*)malloc(n*sizeof(double));
	double *tq = (double*)malloc(n*sizeof(double)), *dh = (double*)malloc(n*sizeof(double));
	int *par = (int*)malloc(n*sizeof(int)); char *fix = (char*)malloc(n);
	/* The unknown lives in the box: the RK predictor is only a first guess. */
	hq_dae_nclip = 0; hq_dae_ommax_pred = 0;
	for(i=0;i<n;i++){
		size_t id = (size_t)PINDX((treevorork4particletype*)(bp+i));
		om[i] = hq_oms[id]; h2[i] = hq_h[id]*hq_h[id];
		if(fabs(om[i]) > hq_dae_ommax_pred) hq_dae_ommax_pred = fabs(om[i]);
		if(om[i] > om_max){ om[i] = om_max; hq_dae_nclip++; }
		if(om[i] < -om_max){ om[i] = -om_max; hq_dae_nclip++; }
	}
	hq_setW(bp, n, om);
	double res = 0, vs = 0, alpha = 1;
	hq_dae_it = 0;
#define HQ_RES() do{ \
		hq_geom(simpar, xmin, ymin, xmax, ymax, paddingAllTreeParticles, find2DNeighboringBP, find2DCellBP, mkLinkedList2D, imap, max_id, V); \
		for(i=0;i<n;i++) par[i] = i; \
		for(f=0;f<dyn_fn;f++){ int a_ = hq_find(par, dyn_f[f].i), b_ = hq_find(par, dyn_f[f].j); if(a_ != b_) par[a_] = b_; } \
		res = 0; vs = 0; \
		for(i=0;i<n;i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(bp+i)); \
			if(!hq_init[id]){ hq_vt0[id] = hq_vts[id] = V[i]; hq_init[id] = 1; } \
			r[i] = hq_vts[id] - V[i]; vs += V[i]; } \
		{ double *cs = (double*)calloc(n, sizeof(double)); int *cn = (int*)calloc(n, sizeof(int)); \
		  for(i=0;i<n;i++){ int c = hq_find(par, i); cs[c] += r[i]; cn[c]++; } \
		  for(i=0;i<n;i++){ int c = hq_find(par, i); r[i] -= cs[c]/cn[c]; } free(cs); free(cn); } \
		for(i=0;i<n;i++) res += r[i]*r[i]; res = sqrt(res); }while(0)
#define HQ_APPLY(v, out) do{ for(i=0;i<n;i++) out[i] = 0; \
	for(f=0;f<dyn_fn;f++){ dyn_face_t *F = dyn_f+f; double t = F->a*(v[F->i]-v[F->j]); out[F->i] += t; out[F->j] -= t; } }while(0)
	HQ_RES();
	hq_dae_res0 = res/(vs > 0 ? vs : 1);
	for(it=0; it<maxit && res > tol*vs; it++){
		/* Newton direction in omega: (D K D) y = D r, D = diag h^2, active set on the box */
		for(i=0;i<n;i++){ y[i] = 0; fix[i] = (om[i] >= om_max && r[i] > 0) || (om[i] <= -om_max && r[i] < 0); }
		for(i=0;i<n;i++) dh[i] = 0;
		for(f=0;f<dyn_fn;f++){ dyn_face_t *F = dyn_f+f; dh[F->i] += F->a*h2[F->i]*h2[F->i]; dh[F->j] += F->a*h2[F->j]*h2[F->j]; }
		for(i=0;i<n;i++) if(!(dh[i] > 0)) dh[i] = 1;
		double bn = 0, rz = 0, rn = 0;
		for(i=0;i<n;i++){ rr[i] = fix[i] ? 0 : h2[i]*r[i]; bn += rr[i]*rr[i]; }
		bn = sqrt(bn);
		for(i=0;i<n;i++){ z[i] = rr[i]/dh[i]; pp[i] = z[i]; rz += rr[i]*z[i]; rn += rr[i]*rr[i]; }
		int it2 = 0;
		while(bn > 0 && sqrt(rn) > 1e-10*bn && it2 < 2000){
			for(i=0;i<n;i++) tq[i] = h2[i]*pp[i];
			HQ_APPLY(tq, q);
			double pq = 0; for(i=0;i<n;i++){ q[i] = fix[i] ? 0 : h2[i]*q[i]; pq += pp[i]*q[i]; }
			if(!(pq > 0)) break;
			double al = rz/pq, rz2 = 0; rn = 0;
			for(i=0;i<n;i++){ y[i] += al*pp[i]; rr[i] -= al*q[i]; z[i] = rr[i]/dh[i]; rz2 += rr[i]*z[i]; rn += rr[i]*rr[i]; }
			double be = rz2/rz; rz = rz2;
			for(i=0;i<n;i++) pp[i] = z[i] + be*pp[i];
			it2++;
		}
		/* projected step with backtracking on |r| */
		double res_old = res; int ls;
		for(i=0;i<n;i++) om_old[i] = om[i];
		alpha = 1;
		for(ls=0; ls<6; ls++){
			for(i=0;i<n;i++){
				double v = om_old[i] + alpha*y[i];
				if(v > om_max) v = om_max;
				if(v < -om_max) v = -om_max;
				om[i] = v;
			}
			hq_setW(bp, n, om);
			HQ_RES();
			if(res < res_old) break;
			alpha *= 0.5;
		}
		hq_dae_it++;
		if(!(res < res_old)){
			for(i=0;i<n;i++) om[i] = om_old[i];
			hq_setW(bp, n, om);
			HQ_RES();
			break;
		}
	}
#undef HQ_RES
#undef HQ_APPLY
	hq_dae_res = res/(vs > 0 ? vs : 1);
	hq_dae_nsat = 0;
	for(i=0;i<n;i++){
		size_t id = (size_t)PINDX((treevorork4particletype*)(bp+i));
		int out = (om[i] >= om_max && r[i] > 0) || (om[i] <= -om_max && r[i] < 0);
		if(out){
			hq_dae_nsat++;
			if(step_start){ hq_vts[id] = hq_vt0[id] = V[i]; }   /* geometry decides */
		}
		if(step_start) hq_om0[id] = om[i];
	}
	free(V); free(r); free(om); free(om_old); free(h2); free(y); free(rr); free(z); free(pp); free(q); free(tq); free(dh); free(par); free(fix); free(imap);
}
/* Scheme (H) W law (Docs/gfs_scheme3_formulation.md §19), one stage:
 *   dV^T/dt  = sum_f A u*_f + (V_vor - V^T)/tau_V
 *   1/2 L_{A/d} Wdot = sum_f A(u*_f - u~_f - rho_f) + (V_vor - V^T)/tau_V + (V^T - V)/tau_C
 * i.e. the geometric volume follows V^T. u* is the cell-state HLL contact
 * speed; V is the Laguerre area, V_vor the W = 0 area of the same
 * generators (second tessellation). The right side is projected onto the
 * range of L (zero mean per connected face graph). Returns no dt limit:
 * validity |dW| < d^2 is only counted ([SCHH]). */
static double gfs_schemeH_solve(SimParameters *simpar, int stage, postype Gamma,
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *)){
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype boxsize = VORO_CELLBOX(simpar);
	int mx = BASICCELL_MX(simpar), my = BASICCELL_MY(simpar);
	double tauV = gfs_env_d("GFS_H_TAUV", 1.0), tauC = gfs_env_d("GFS_H_TAUC", 0.02);
	int use_muscl_h = GAS_USEMUSCL(simpar) && (int)gfs_env_d("GFS_H_USTAR_MUSCL", 1);
	double tol = gfs_env_d("GFS_DYN_W2_TOL", 1e-10);
	int maxit = (int)gfs_env_d("GFS_DYN_W2_MAXIT", 2000);
	int i, iy;
	size_t max_id = 0;
	for(i=0;i<nbp;i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(bp+i)); if(id > max_id) max_id = id; }
	dyn_ensure(max_id); hq_ensure(max_id);
	int omega = gfs_h_omega_on();
	double om_max = gfs_env_d("GFS_H_OMAX", 0.4);
	if(omega){
		hq_oensure(max_id);
		hq_Lx = (double)(SIMBOX(simpar).x.max - SIMBOX(simpar).x.min);
		hq_Ly = (double)(SIMBOX(simpar).y.max - SIMBOX(simpar).y.min);
		hq_R1 = gfs_env_d("GFS_H_R1", 1.0)*hq_Lx/NX(simpar);
		hq_compute_h((treevorork4particletype*)bp, sizeof(treevorostressrk4particletype), nbp);
	}
	int *imap = (int*)malloc((max_id+1)*sizeof(int));
	for(size_t k=0;k<=max_id;k++) imap[k] = -1;
	for(i=0;i<nbp;i++) imap[(size_t)PINDX((treevorork4particletype*)(bp+i))] = i;
	double *vgeo = (double*)calloc(nbp, sizeof(double)), *vvor = (double*)calloc(nbp, sizeof(double));
	char *seen = (char*)calloc(nbp, 1), *seen2 = (char*)calloc(nbp, 1);
	double *capr = (double*)calloc(nbp, sizeof(double));   /* v_i . sum_cap A n */
	dyn_fn = 0;
	for(iy=1;iy<my-1;iy++){
		int mp = 4000, ix;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		Voro2D_Corner *vc0 = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		size_t nb = 0, cb = 256;
		dyn_face_t *buf = (dyn_face_t*)malloc(cb*sizeof(dyn_face_t));
		for(ix=1;ix<mx-1;ix++){
			int np, nneigh, k;
			treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
			Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
			Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			for(k=0;k<np;k++){
				treevorork4particletype *ib = p[k].bp;
				treevorostressrk4particletype *ibs = (treevorostressrk4particletype*)ib;
				Voro2D_point center;
				center.x = p[k].x; center.y = p[k].y;
				center.indx = PINDX(p+k); center.csound = p[k].csound; center.w2 = p[k].w2;
				Voro2D_point *nb_use = neighbors, *nw_use = neighwork, *ne = NULL, *nwe = NULL;
				int nn_use = nneigh;
				int ip = Voro2D_FindVC(&center, nb_use, nw_use, nn_use, vorocorner, mp, boxsize);
				if(Voro2D_HasUnassignedCorner(vorocorner, ip) || voro_sec_short(simpar, vorocorner, ip, (double)center.x, (double)center.y, ix, iy, 1)){
					int R;
					for(R=2;R<=voro_expand_rmax(simpar);R++){
						int nne;
						if(ne){ free(ne); free(nwe); }
						ne = searchCellRk4Neighbors2D_R(simpar, ix, iy, R, &nne);
						nwe = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nne);
						ip = Voro2D_FindVC(&center, ne, nwe, nne, vorocorner, mp, boxsize);
						nb_use = ne; nw_use = nwe; nn_use = nne;
						if(!Voro2D_HasUnassignedCorner(vorocorner, ip) && !voro_sec_short(simpar, vorocorner, ip, (double)center.x, (double)center.y, ix, iy, R)) break;
					}
				}
				size_t idi = (size_t)PINDX(ib);
				int ii = (idi <= max_id) ? imap[idi] : -1;
				if(ii >= 0 && !seen[ii]){
					/* Voronoi (W = 0) cell of the same generators. */
					Voro2D_point *z = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nn_use);
					Voro2D_point *zw = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nn_use);
					Voro2D_point c0 = center; int q;
					for(q=0;q<nn_use;q++){ z[q] = nb_use[q]; z[q].w2 = 0; }
					c0.w2 = 0;
					Voro2D_FindVC(&c0, z, zw, nn_use, vc0, mp, boxsize);
					vvor[ii] = hq_polyarea(vc0);
					vgeo[ii] = hq_polyarea(vorocorner);
					seen[ii] = 1;
					free(z); free(zw);
				}
				Voro2D_Corner *tmp = vorocorner;
				do {
					if(tmp->upperrelated < 0 && ii >= 0 && !seen2[ii]){
						Voro2D_point cl; cl.x = tmp->upperlink->x - tmp->x; cl.y = tmp->upperlink->y - tmp->y;
						Voro2D_point cS = voro2D_norm(&cl);   /* unit normal */
						capr[ii] += (double)Vec2DLength(tmp, tmp->upperlink)*((double)ib->vx*cS.x + (double)ib->vy*cS.y);
					}
					if(tmp->upperrelated >= 0 && tmp->upperrelated < nn_use && nw_use[tmp->upperrelated].bp){
						treevorork4particletype *jb = (treevorork4particletype*)nw_use[tmp->upperrelated].bp;
						treevorostressrk4particletype *jbs = (treevorostressrk4particletype*)jb;
						if(PINDX(jb) != MAX_INDEX && (size_t)PINDX(jb) > idi && (size_t)PINDX(jb) <= max_id
								&& imap[(size_t)PINDX(jb)] >= 0 && ii >= 0){
							Voro2D_Corner *tmp2 = tmp->upperlink;
							double area = Vec2DLength(tmp,tmp2);
							double rx = jb->x - ib->x, ry = jb->y - ib->y;
							double d = sqrt(rx*rx + ry*ry);
							if(d > 0 && area > 0){
								double ex = rx/d, ey = ry/d;
								double vni = ib->vx*ex + ib->vy*ey, vnj = jb->vx*ex + jb->vy*ey;
								double dd = vnj - vni, nub = 0.5*(vni + vnj);
								double dW = (double)ib->w2 - (double)jb->w2;
								Voro2D_point ur; ur.x = 0; ur.y = 0;
								voronoi_face_rotation(&ur, (postype)0.5*(tmp->x + tmp2->x), (postype)0.5*(tmp->y + tmp2->y),
									ib->x, ib->y, jb->x, jb->y, ib->vx, ib->vy, jb->vx, jb->vy,
									(postype)ex, (postype)ey, (postype)d, ib->w2, jb->w2, 1);
								double rot = ur.x*ex + ur.y*ey;
								/* u*: the Riemann contact speed of the face states the
								 * force uses (face-loop MUSCL at the face centroid +
								 * rest-frame HLLC; cell-state HLL on extreme faces). */
								postype ps, vns;
								double pmn_ = fmin((double)ibs->pressure, (double)jbs->pressure);
								double pmx_ = fmax((double)ibs->pressure, (double)jbs->pressure);
								if(pmx_ > 100.0*fmax(pmn_, 1e-30)){
									hll_star_state(ibs->den, ibs->pressure, (postype)vni, ibs->csound,
										jbs->den, jbs->pressure, (postype)vnj, jbs->csound, Gamma, &ps, &vns);
								} else {
									Voro2D_point fl; fl.x = tmp2->x - tmp->x; fl.y = tmp2->y - tmp->y;
									Voro2D_point nh = voro2D_norm(&fl);
									postype xf = 0.5*(tmp->x + tmp2->x), yf = 0.5*(tmp->y + tmp2->y);
									postype xjf = xf - (jb->x - ib->x), yjf = yf - (jb->y - ib->y);
									postype rL = ibs->den, rR = jbs->den, pL_ = ibs->pressure, pR_ = jbs->pressure;
									postype uL = ib->vx*nh.x + ib->vy*nh.y, uR = jb->vx*nh.x + jb->vy*nh.y;
									if(use_muscl_h){
										rL += ibs->stress.dRhodx*xf + ibs->stress.dRhody*yf;
										rR += jbs->stress.dRhodx*xjf + jbs->stress.dRhody*yjf;
										pL_ += ibs->stress.dPdx*xf + ibs->stress.dPdy*yf;
										pR_ += jbs->stress.dPdx*xjf + jbs->stress.dPdy*yjf;
										uL += nh.x*(ibs->stress.gUxx*xf + ibs->stress.gUxy*yf) + nh.y*(ibs->stress.gUyx*xf + ibs->stress.gUyy*yf);
										uR += nh.x*(jbs->stress.gUxx*xjf + jbs->stress.gUxy*yjf) + nh.y*(jbs->stress.gUyx*xjf + jbs->stress.gUyy*yjf);
									}
									double rL_raw = rL, rR_raw = rR;
									if(pL_ < 1e-10) pL_ = 1e-10;
									if(pR_ < 1e-10) pR_ = 1e-10;
									if(rL < 1e-10) rL = 1e-10;
									if(rR < 1e-10) rR = 1e-10;
									hllc_face_2d_rest_frame(rL, pL_, uL, ibs->csound, rR, pR_, uR, jbs->csound,
										(postype)0.5*(uL + uR), Gamma, &ps, &vns);
									if((int)gfs_env_d("GFS_H_DEBUG", 0) == 2){
										double sc = fabs((double)uL - (double)uR) + (double)ibs->csound + (double)jbs->csound + 1e-30;
										if(fabs((double)vns - 0.5*((double)uL + (double)uR)) > 10*sc){
											double fdist = sqrt((double)(xf*xf + yf*yf)), fdj = sqrt((double)(xjf*xjf + yjf*yjf));
											fprintf(stderr, "[SCHH_US] stage=%d u*=%.3e uL=%.3e uR=%.3e rho_raw L=%.3e R=%.3e (cell %.3e %.3e) gradrho L=(%.3e,%.3e) R=(%.3e,%.3e) |xf-xi|=%.3e |xf-xj|=%.3e pL=%.3e pR=%.3e (cell %.3e %.3e) d=%.3e expanded=%d nn=%d\n",
												stage, (double)vns, (double)uL, (double)uR, rL_raw, rR_raw, (double)ibs->den, (double)jbs->den,
												(double)ibs->stress.dRhodx, (double)ibs->stress.dRhody, (double)jbs->stress.dRhodx, (double)jbs->stress.dRhody,
												fdist, fdj, (double)pL_, (double)pR_, (double)ibs->pressure, (double)jbs->pressure, d, (int)(ne != NULL), nn_use);
										}
									}
									/* along e (nh = e for a Laguerre face up to rounding) */
									vns = vns*(postype)(nh.x*ex + nh.y*ey);
								}
								dyn_face_t f;
								f.i = ii; f.j = imap[(size_t)PINDX(jb)]; f.ext = 0;
								f.a = area/(2*d);
								/* 2d(u* - u~ - rho); u~ = phi u_i + (1-phi) u_j */
								f.b = 2*d*((double)vns - nub - rot) + dW*dd/d;
								f.d = d; f.dd = dd; f.dW = dW;
								f.vstar_rel = area*(double)vns;   /* A u*, outward for i */
								if(nb == cb){ cb *= 2; buf = (dyn_face_t*)realloc(buf, cb*sizeof(dyn_face_t)); }
								buf[nb++] = f;
							}
						}
					}
					tmp = tmp->upperlink;
				} while(tmp != vorocorner);
				if(ii >= 0) seen2[ii] = 1;
				if(ne){ free(ne); free(nwe); }
			}
			free(p); free(neighbors); free(neighwork);
		}
		dyn_face_push(buf, nb);
		free(buf); free(vorocorner); free(vc0);
	}
	size_t f; int n = nbp;
	double *x = (double*)calloc(n, sizeof(double)), *g = (double*)calloc(n, sizeof(double));
	double *r = (double*)malloc(n*sizeof(double)), *z = (double*)malloc(n*sizeof(double));
	double *pp = (double*)malloc(n*sizeof(double)), *q = (double*)malloc(n*sizeof(double));
	double *dg = (double*)calloc(n, sizeof(double)), *au = (double*)calloc(n, sizeof(double));
	int *par = (int*)malloc(n*sizeof(int));
	for(i=0;i<n;i++) par[i] = i;
	for(f=0;f<dyn_fn;f++){
		dyn_face_t *F = dyn_f+f;
		g[F->i] += F->a*F->b; g[F->j] -= F->a*F->b;
		dg[F->i] += F->a; dg[F->j] += F->a;
		au[F->i] += F->vstar_rel; au[F->j] -= F->vstar_rel;
		int a_ = hq_find(par, F->i), b_ = hq_find(par, F->j); if(a_ != b_) par[a_] = b_;
	}
	hq_vrel2 = hq_vtrel2 = 0; hq_ncell = 0;
	for(i=0;i<n;i++){
		size_t id = (size_t)PINDX((treevorork4particletype*)(bp+i));
		if(!hq_init[id]){ hq_vt0[id] = hq_vts[id] = vgeo[i]; hq_init[id] = 1; }
		double vt = hq_vts[id];
		double relax = (vvor[i] - vt)/tauV;
		hq_kvt[stage][id] = au[i] + capr[i] + relax;
		g[i] += relax + (gfs_h_dae_on() ? 0.0 : (vt - vgeo[i])/tauC);
		if(vgeo[i] > 0){ hq_vrel2 += pow((vt - vgeo[i])/vgeo[i], 2); hq_vtrel2 += pow((vvor[i] - vt)/vgeo[i], 2); hq_ncell++; }
	}
	if((int)gfs_env_d("GFS_H_DEBUG", 0)){
		/* top |g| cells: face term, relaxation, feedback, cap rate */
		static int dbg_calls = 0;
		if(dbg_calls++ < 8){
			double cx0 = gfs_env_d("EUNHA_KEPLER_CX", 0), cy0 = gfs_env_d("EUNHA_KEPLER_CY", 0);
			int top[5] = {-1,-1,-1,-1,-1}, a_, b_;
			for(i=0;i<n;i++){
				for(a_=0;a_<5;a_++) if(top[a_] < 0 || fabs(g[i]) > fabs(g[top[a_]])) break;
				if(a_ < 5){ for(b_=4;b_>a_;b_--) top[b_] = top[b_-1]; top[a_] = i; }
			}
			for(a_=0;a_<5 && top[a_]>=0;a_++){
				int k_ = top[a_]; size_t id = (size_t)PINDX((treevorork4particletype*)(bp+k_));
				double vt = hq_vts[id];
				double rr = sqrt(pow((double)bp[k_].x - cx0, 2) + pow((double)bp[k_].y - cy0, 2));
				fprintf(stderr, "[SCHH_DBG] call=%d stage=%d r=%.4f g=%.3e face=%.3e relax=%.3e fb=%.3e capr=%.3e au=%.3e V=%.3e Vvor=%.3e VT=%.3e diag=%.3e W=%.3e P=%.3e\n",
					dbg_calls, stage, rr, g[k_], g[k_] - (vvor[k_] - vt)/tauV - (vt - vgeo[k_])/tauC,
					(vvor[k_] - vt)/tauV, (vt - vgeo[k_])/tauC, capr[k_], au[k_], vgeo[k_], vvor[k_], vt, dg[k_],
					(double)bp[k_].w2, (double)bp[k_].pressure);
			}
		}
	}
	/* project g onto range(L): zero mean on each connected face graph */
	{
		double *cs = (double*)calloc(n, sizeof(double)); int *cn = (int*)calloc(n, sizeof(int));
		hq_ncomp = 0;
		for(i=0;i<n;i++){ int c = hq_find(par, i); cs[c] += g[i]; cn[c]++; }
		for(i=0;i<n;i++) if(hq_find(par, i) == i) hq_ncomp++;
		double gp2 = 0, gg2 = 0;
		for(i=0;i<n;i++){ int c = hq_find(par, i); gg2 += g[i]*g[i]; gp2 += pow(cs[c]/cn[c], 2); g[i] -= cs[c]/cn[c]; }
		hq_proj = gg2 > 0 ? sqrt(gp2/gg2) : 0;
		free(cs); free(cn);
	}
	for(i=0;i<n;i++){ x[i] = dyn_wd[(size_t)PINDX((treevorork4particletype*)(bp+i))]; if(!(dg[i] > 0)) dg[i] = 1; }
#define DYN_APPLY(v, out) do{ for(i=0;i<n;i++) out[i] = 0; \
	for(f=0;f<dyn_fn;f++){ dyn_face_t *F = dyn_f+f; double t = F->a*(v[F->i]-v[F->j]); out[F->i] += t; out[F->j] -= t; } }while(0)
	DYN_APPLY(x, q);
	double gn = 0, rz = 0, rn = 0;
	for(i=0;i<n;i++){ r[i] = g[i] - q[i]; gn += g[i]*g[i]; }
	gn = sqrt(gn);
	for(i=0;i<n;i++){ z[i] = r[i]/dg[i]; pp[i] = z[i]; rz += r[i]*z[i]; rn += r[i]*r[i]; }
	int it = 0;
	while(gn > 0 && sqrt(rn) > tol*gn && it < maxit){
		DYN_APPLY(pp, q);
		double pq = 0; for(i=0;i<n;i++) pq += pp[i]*q[i];
		if(!(pq > 0)) break;
		double al = rz/pq, rz2 = 0; rn = 0;
		for(i=0;i<n;i++){ x[i] += al*pp[i]; r[i] -= al*q[i]; z[i] = r[i]/dg[i]; rz2 += r[i]*z[i]; rn += r[i]*r[i]; }
		double be = rz2/rz; rz = rz2;
		for(i=0;i<n;i++) pp[i] = z[i] + be*pp[i];
		it++;
	}
	if(it > dyn_itmax) dyn_itmax = it;
	dyn_r2o += rn; dyn_b2o += gn*gn;
	double dtlim = 1e30;
	if(omega){
		/* min_y 1/2 Wd^T K Wd - g^T Wd, Wd = h^2 y + c, c = 2 h hdot omega,
		 * y_i <= 0 where omega_i >= om_max, y_i >= 0 where omega_i <= -om_max.
		 * Active set: fixed y = 0 on blocked components; PCG on the rest of
		 * H y = b, H = D K D (D = diag h^2), b = D (g - K c). */
		double *h2 = (double*)malloc(n*sizeof(double)), *cc = (double*)malloc(n*sizeof(double));
		double *om = (double*)malloc(n*sizeof(double)), *b = (double*)malloc(n*sizeof(double));
		double *y = (double*)calloc(n, sizeof(double)), *tq = (double*)malloc(n*sizeof(double));
		double *dh = (double*)calloc(n, sizeof(double));
		char *fix = (char*)calloc(n, 1);
		for(i=0;i<n;i++){
			size_t id = (size_t)PINDX((treevorork4particletype*)(bp+i));
			h2[i] = hq_h[id]*hq_h[id];
			om[i] = hq_oms[id];
			cc[i] = 2*hq_h[id]*hq_hd[id]*om[i];
		}
		DYN_APPLY(cc, tq);
		for(i=0;i<n;i++) b[i] = h2[i]*(g[i] - tq[i]);
		for(f=0;f<dyn_fn;f++){ dyn_face_t *F = dyn_f+f; dh[F->i] += F->a*h2[F->i]*h2[F->i]; dh[F->j] += F->a*h2[F->j]*h2[F->j]; }
		for(i=0;i<n;i++) if(!(dh[i] > 0)) dh[i] = 1;
		int pass, npass = 0;
		for(i=0;i<n;i++) fix[i] = (om[i] >= om_max || om[i] <= -om_max);
		for(pass=0;pass<20;pass++){
			npass++;
			/* PCG on free components, y_fixed = 0 */
			double bn = 0, rz_ = 0, rn_ = 0;
			for(i=0;i<n;i++){ if(fix[i]) y[i] = 0; }
			for(i=0;i<n;i++) z[i] = h2[i]*y[i];
			DYN_APPLY(z, q);
			for(i=0;i<n;i++){ r[i] = fix[i] ? 0 : b[i] - h2[i]*q[i]; bn += fix[i] ? 0 : b[i]*b[i]; }
			bn = sqrt(bn);
			for(i=0;i<n;i++){ z[i] = r[i]/dh[i]; pp[i] = z[i]; rz_ += r[i]*z[i]; rn_ += r[i]*r[i]; }
			int it2 = 0;
			while(bn > 0 && sqrt(rn_) > tol*bn && it2 < maxit){
				for(i=0;i<n;i++) tq[i] = h2[i]*pp[i];
				DYN_APPLY(tq, q);
				double pq = 0; for(i=0;i<n;i++){ q[i] = fix[i] ? 0 : h2[i]*q[i]; pq += pp[i]*q[i]; }
				if(!(pq > 0)) break;
				double al = rz_/pq, rz2 = 0; rn_ = 0;
				for(i=0;i<n;i++){ y[i] += al*pp[i]; r[i] -= al*q[i]; z[i] = r[i]/dh[i]; rz2 += r[i]*z[i]; rn_ += r[i]*r[i]; }
				double be = rz2/rz_; rz_ = rz2;
				for(i=0;i<n;i++) pp[i] = z[i] + be*pp[i];
				it2++;
			}
			if(it2 > dyn_itmax) dyn_itmax = it2;
			/* KKT: gradient G = H y - b on fixed components; release a blocked
			 * component whose gradient points into the box; block a free
			 * saturated component moving outward. */
			for(i=0;i<n;i++) tq[i] = h2[i]*y[i];
			DYN_APPLY(tq, q);
			int changed = 0;
			for(i=0;i<n;i++){
				double G = h2[i]*q[i] - b[i];
				if(fix[i]){
					if((om[i] >= om_max && G > 0) || (om[i] <= -om_max && G < 0)){ fix[i] = 0; changed = 1; }
				} else if((om[i] >= om_max && y[i] > 0) || (om[i] <= -om_max && y[i] < 0)){ fix[i] = 1; changed = 1; }
			}
			if(!changed) break;
		}
		if((int)gfs_env_d("GFS_H_DEBUG", 0)){
			static int odbg = 0;
			if(odbg++ < 12){
				int im = -1, nviol = 0; double ym = 0;
				for(i=0;i<n;i++){
					if(fabs(y[i]) > ym){ ym = fabs(y[i]); im = i; }
					if((om[i] >= om_max && y[i] > 0) || (om[i] <= -om_max && y[i] < 0)) nviol++;
				}
				if(im >= 0){
					size_t idm = (size_t)PINDX((treevorork4particletype*)(bp+im));
					fprintf(stderr, "[SCHH_OMDBG] stage=%d passes=%d nviol=%d max|y|=%.3e at omega=%.4f h=%.4e hdot=%.3e fix=%d g=%.3e b=%.3e c=%.3e\n",
						stage, npass, nviol, ym, om[im], hq_h[idm], hq_hd[idm], (int)fix[im], g[im], b[im], cc[im]);
				}
			}
		}
		if((int)gfs_env_d("GFS_H_DEBUG", 0) == 2){
			/* largest actual Wdot = h^2 y + c this stage, with its neighbourhood */
			int im = -1; double wm = 0, cx0 = gfs_env_d("EUNHA_KEPLER_CX", 0), cy0 = gfs_env_d("EUNHA_KEPLER_CY", 0);
			for(i=0;i<n;i++){ double wd = fabs(h2[i]*y[i] + cc[i]); if(wd > wm){ wm = wd; im = i; } }
			if(im >= 0){
				size_t idm = (size_t)PINDX((treevorork4particletype*)(bp+im));
				double dmin = 1e30, amax = 0; int nf = 0;
				for(f=0;f<dyn_fn;f++){ dyn_face_t *F = dyn_f+f; if(F->i == im || F->j == im){ nf++; if(F->d < dmin) dmin = F->d; if(F->a > amax) amax = F->a; } }
				fprintf(stderr, "[SCHH_WD] stage=%d max|Wdot|=%.3e r=%.4f h=%.3e omega=%.4f y=%.3e c=%.3e fix=%d g=%.3e V=%.3e VT=%.3e nfaces=%d dmin=%.3e amax=%.3e\n",
					stage, wm, sqrt(pow((double)bp[im].x - cx0, 2) + pow((double)bp[im].y - cy0, 2)), hq_h[idm], om[im], y[im], cc[im], (int)fix[im],
					g[im], vgeo[im], hq_vts[idm], nf, dmin, amax);
			}
		}
		hq_nsat = 0; hq_ommax_seen = 0;
		for(i=0;i<n;i++){
			size_t id = (size_t)PINDX((treevorork4particletype*)(bp+i));
			x[i] = h2[i]*y[i] + cc[i];
			hq_kom[stage][id] = y[i];
			if(fix[i]) hq_nsat++;
			if(fabs(om[i]) > hq_ommax_seen) hq_ommax_seen = fabs(om[i]);
			/* stay inside |omega| < om_ok over this step */
		}
		free(h2); free(cc); free(om); free(b); free(y); free(tq); free(dh); free(fix);
	}
#undef DYN_APPLY
	for(i=0;i<n;i++){
		size_t id = (size_t)PINDX((treevorork4particletype*)(bp+i));
		dyn_wd[id] = x[i]; dyn_k[stage][id] = x[i];
	}
	hq_rmax = 0; hq_nbad = 0;
	for(f=0;f<dyn_fn;f++){
		dyn_face_t *F = dyn_f+f;
		double wr = fabs(F->dW)/(F->d*F->d);
		if(wr > hq_rmax) hq_rmax = wr;
		if(wr >= 1) hq_nbad++;
	}
	if(hq_rmax > dyn_wratio_max) dyn_wratio_max = hq_rmax;
	dyn_nbad += hq_nbad; dyn_no += (long)dyn_fn;
	free(x); free(g); free(r); free(z); free(pp); free(q); free(dg); free(au); free(par);
	free(imap); free(vgeo); free(vvor); free(seen); free(seen2); free(capr);
	return dtlim;
}
/* Stage W = W_n + cdt*k[kst] (kst < 0: final RK4 combination with Dt = cdt),
 * shifted so min W = 0; w2old = w2. */
static void gfs_dynw2_set(treevorostressrk4particletype *sbp, int np, postype cdt, int kst){
	int i; double wmin = 1e300;
	if(gfs_schemeH_on() && gfs_h_omega_on() && hq_ocap > 0){
		/* omega = omega_n + cdt k; W = h(x)^2 omega at the stage positions */
		hq_compute_h((treevorork4particletype*)sbp, sizeof(treevorostressrk4particletype), np);
		for(i=0;i<np;i++){
			size_t id = (size_t)PINDX((treevorork4particletype*)(sbp+i));
			if(id >= hq_ocap) continue;
			double k = (kst >= 0) ? hq_kom[kst][id]
				: (hq_kom[0][id] + 2*hq_kom[1][id] + 2*hq_kom[2][id] + hq_kom[3][id])/6.0;
			hq_oms[id] = hq_om0[id] + (double)cdt*k;
			if(kst < 0) hq_om0[id] = hq_oms[id];
			double W = hq_h[id]*hq_h[id]*hq_oms[id];
			if(W < wmin) wmin = W;
		}
		for(i=0;i<np;i++){
			size_t id = (size_t)PINDX((treevorork4particletype*)(sbp+i));
			double W = (id < hq_ocap) ? hq_h[id]*hq_h[id]*hq_oms[id] - wmin : 0;
			if(!(W > 0)) W = 0;
			sbp[i].w2 = (postype)W; sbp[i].w2old = sbp[i].w2;
		}
		if(gfs_schemeH_on()){
			for(i=0;i<np;i++){
				size_t id = (size_t)PINDX((treevorork4particletype*)(sbp+i));
				if(id >= hq_cap || !hq_init[id]) continue;
				double k = (kst >= 0) ? hq_kvt[kst][id]
					: (hq_kvt[0][id] + 2*hq_kvt[1][id] + 2*hq_kvt[2][id] + hq_kvt[3][id])/6.0;
				hq_vts[id] = hq_vt0[id] + (double)cdt*k;
				if(kst < 0) hq_vt0[id] = hq_vts[id];
			}
		}
		return;
	}
	for(i=0;i<np;i++){
		size_t id = (size_t)PINDX((treevorork4particletype*)(sbp+i));
		double k = (kst >= 0) ? dyn_k[kst][id]
			: (dyn_k[0][id] + 2*dyn_k[1][id] + 2*dyn_k[2][id] + dyn_k[3][id])/6.0;
		double W = (double)sbp[i].rk4.w2backup + (double)cdt*k;
		sbp[i].w2 = (postype)W;
		if(W < wmin) wmin = W;
	}
	for(i=0;i<np;i++){
		size_t id = (size_t)PINDX((treevorork4particletype*)(sbp+i));
		double k = (kst >= 0) ? dyn_k[kst][id]
			: (dyn_k[0][id] + 2*dyn_k[1][id] + 2*dyn_k[2][id] + dyn_k[3][id])/6.0;
		double W = (double)sbp[i].rk4.w2backup + (double)cdt*k - wmin;
		if(!(W > 0)) W = 0;   /* rounding only; sqrt(w2) is taken elsewhere */
		sbp[i].w2 = (postype)W; sbp[i].w2old = sbp[i].w2;
	}	if(gfs_schemeH_on()){
		for(i=0;i<np;i++){
			size_t id = (size_t)PINDX((treevorork4particletype*)(sbp+i));
			if(id >= hq_cap || !hq_init[id]) continue;
			double k = (kst >= 0) ? hq_kvt[kst][id]
				: (hq_kvt[0][id] + 2*hq_kvt[1][id] + 2*hq_kvt[2][id] + hq_kvt[3][id])/6.0;
			hq_vts[id] = hq_vt0[id] + (double)cdt*k;
			if(kst < 0) hq_vt0[id] = hq_vts[id];
		}
	}
}
static void gfs_dynw2_report(SimParameters *simpar, postype Dt, double dtlim, int bound){
	static int st = 0; static long nbind = 0;
	st++; if(bound) nbind++;
	if(MYID(simpar) != 0) return;
	fprintf(stderr, "[DYNW2] step=%d Dt=%.4e faces/stage=%ld ext=%ld relres_ext=%.3e relres_ord=%.3e relres_all=%.3e rms_dv_ext=%.3e rms_vstar_rel_ext=%.3e rms_dv_ext_force=%.3e rms_dv_ord=%.3e cg_it_max=%d maxdW/d2=%.3f invalid=%ld nbar=%ld rms_b_ord=%.3e dtlim=%.3e bound=%d nbind=%ld\n",
		st, (double)Dt, (dyn_ne + dyn_no)/4, dyn_ne/4,
		dyn_b2e > 0 ? sqrt(dyn_r2e/dyn_b2e) : 0.0,
		dyn_b2o > 0 ? sqrt(dyn_r2o/dyn_b2o) : 0.0,
		(dyn_b2e + dyn_b2o) > 0 ? sqrt((dyn_r2e + dyn_r2o)/(dyn_b2e + dyn_b2o)) : 0.0,
		dyn_ne > 0 ? sqrt(dyn_dv2e/dyn_ne) : 0.0,
		dyn_ne > 0 ? sqrt(dyn_vs2e/dyn_ne) : 0.0,
		dyn_fchk_n > 0 ? sqrt(dyn_fchk_s2/dyn_fchk_n) : 0.0,
		dyn_no > 0 ? sqrt(dyn_dv2o/dyn_no) : 0.0,
		dyn_itmax, dyn_wratio_max, dyn_nbad/4, dyn_nbar/4, dyn_no > 0 ? sqrt(dyn_b2o/dyn_no) : 0.0, dtlim, bound, nbind);	if(gfs_schemeH_on())
		fprintf(stderr, "[SCHH] step=%d maxdW/d2(last stage)=%.4f nfaces>=1=%ld proj=%.2e rms(VT-V)/V=%.3e rms(Vvor-VT)/V=%.3e comps=%ld cells=%ld\n",
			st, hq_rmax, hq_nbad, hq_proj, hq_ncell > 0 ? sqrt(hq_vrel2/hq_ncell) : 0.0,
			hq_ncell > 0 ? sqrt(hq_vtrel2/hq_ncell) : 0.0, hq_ncomp, hq_ncell);
	if(gfs_schemeH_on() && gfs_h_omega_on())
		fprintf(stderr, "[SCHH_OM] step=%d max|omega|=%.4f saturated=%ld dtlim=%.3e\n", st, hq_ommax_seen, hq_nsat, dtlim);
	if(gfs_schemeH_on() && gfs_h_dae_on())
		fprintf(stderr, "[SCHH_DAE] step=%d newton_it=%d |V^T-V|/sumV before=%.3e after=%.3e out_of_reach=%d predictor_outside_box=%d max|omega_pred|=%.3e\n",
			st, hq_dae_it, hq_dae_res0, hq_dae_res, hq_dae_nsat, hq_dae_nclip, hq_dae_ommax_pred);
	if(gfs_schemeH_on()){
		fprintf(stderr, "[SCHH_D] step=%d sumD+=%.4e sumD-=%.4e  dPdelta+=%.4e -=%.4e  dPstar+=%.4e -=%.4e (both sides, 4 stages)\n",
			st, hq_dpos, hq_dneg, hq_d1p, hq_d1n, hq_d2p, hq_d2n);
		hq_dpos = hq_dneg = hq_d1p = hq_d1n = hq_d2p = hq_d2n = 0;
	}
}
double getAccVoro2DBlend(SimParameters *simpar, postype xmin, postype ymin,
		postype xmax, postype ymax,
		postype OrderOfAccuracy, postype Courant, postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype))
		){
	if(MYID(simpar)<=1){ fprintf(stderr,"[FORCE_DBG] r%d alpha wrapper entry\n", MYID(simpar)); fflush(stderr); }
	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	int mx, my;
	BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
	BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*mx*my);
	mkLinkedList2D(simpar, cellsize, xmin, ymin, xmax, ymax, paddingAllTreeParticles);
	if(MYID(simpar)<=1){ fprintf(stderr,"[FORCE_DBG] r%d beta after mkLinkedList2D\n", MYID(simpar)); fflush(stderr); }
	/* GFS_DYN_W2: solve Wdot for this stage on the linked list just built. */
	if(gfs_dynw2_active){
		dyn_dtlim_last = gfs_schemeH_on()
			? gfs_schemeH_solve(simpar, dyn_stage_cur, Gamma, find2DNeighboringBP, find2DCellBP)
			: gfs_dynw2_solve(simpar, dyn_stage_cur, Gamma, find2DNeighboringBP, find2DCellBP);
	}

	if(gfs_subzone_on()){
		size_t k_; sub_ensure(simpar);
		for(k_ = 0; k_ < sub_cap; k_++){ sub_fx[k_] = sub_fy[k_] = 0; }
	}
	if((gfs_hg_affine() || gfs_subzone_on()) && !sedov_phase1_on()){
		fprintf(stderr, "[HG] GFS_HG_AFFINE / GFS_SUBZONE require SEDOV_PHASE1=1 (die = dE/dt path)\n"); exit(1);
	}
	if(gfs_hg_affine() && gfs_subzone_on()){ fprintf(stderr, "[HG] GFS_HG_AFFINE and GFS_SUBZONE are exclusive (design: Phi = 0)\n"); exit(1); }
	if((gfs_bprime() || gfs_hg_off()) && !(gfs_hg_affine() && gfs_hg_split())){ fprintf(stderr, "[BPRIME] GFS_BPRIME/GFS_HG_OFF require GFS_HG_AFFINE=1 GFS_HG_SPLIT=1\n"); exit(1); }
	if(gfs_bulk_affine() && !(gfs_hg_affine() && gfs_hg_split())){ fprintf(stderr, "[BULK] GFS_BULK_AFFINE requires GFS_HG_AFFINE=1 GFS_HG_SPLIT=1\n"); exit(1); }
	if(gfs_hg_affine()){
		if(NID(simpar) > 1 || !gfs_shared_geom_on()){ fprintf(stderr, "[HG] GFS_HG_AFFINE needs GFS_SHARED_GEOM and one MPI rank\n"); exit(1); }
		hg_ensure(simpar);
		hgo_reset((long)VORO_NP(simpar));
	}
	postype Dtime = getAccVoro2DBlend_impl(simpar, xmin, ymin, xmax, ymax,
		OrderOfAccuracy, Courant, Gamma,
		find2DNeighboringBP, find2DCellBP);
	if(gfs_subzone_on()){
		/* gathered Phi forces; total-energy rate of cell i gets v_i.F_i + dPhi_i/dt
		 * (sum over cells of both is zero identically) */
		size_t ps = TVORORK4_DDINFO(simpar)[0].n_size; char *b = (char*)VORORK4_TBP(simpar); int i_;
		for(i_ = 0; i_ < VORO_NP(simpar); i_++){
			treevorork4particletype *q = (treevorork4particletype*)(b + i_*ps);
			size_t id = (size_t)PINDX(q);
			if(id >= sub_cap) continue;
			q->ax += sub_fx[id]/q->mass; q->ay += sub_fy[id]/q->mass;
			q->die += (postype)((double)q->vx*sub_fx[id] + (double)q->vy*sub_fy[id] + sub_dphi[id]);
		}
		sub_vf_last = 0;
		(void)0;
		for(i_ = 0; i_ < VORO_NP(simpar); i_++){
			treevorork4particletype *q = (treevorork4particletype*)(b + i_*ps); size_t id = (size_t)PINDX(q);
			if(id < sub_cap) sub_vf_last += (double)q->vx*sub_fx[id] + (double)q->vy*sub_fy[id];
		}
	}
	if(gfs_hg_affine() && !gfs_hg_split()){
		/* gathered hourglass forces; total-energy rate of cell i gets v_i.F_i + H_i */
		size_t ps = TVORORK4_DDINFO(simpar)[0].n_size; char *b = (char*)VORORK4_TBP(simpar); int i_;
		hg_vf_last = 0; hg_H_last = 0;
		double lmax = 0;
		for(i_ = 0; i_ < VORO_NP(simpar); i_++){ size_t id = (size_t)PINDX((treevorork4particletype*)(b + i_*ps)); if(id < hg_cap && hg_lam[id] > lmax) lmax = hg_lam[id]; }
		hg_dtlim = (lmax > 0) ? 1.4/lmax : 1e300;   /* half the RK4 real-axis limit 2.785 */
		if(gfs_hg_spec() && hgf_n > 0){
			/* power iteration on B = sum_f k_f z_f z_f^T, z = l/sqrt(m) (diagnostic lower estimate) */
			size_t nc = hg_cap, ii; int it;
			double *xv = (double*)calloc(2*nc, sizeof(double)), *yv = (double*)calloc(2*nc, sizeof(double)), lam = 0, res = 0;
			unsigned long long sd = 12345;
			for(ii = 0; ii < 2*nc; ii++){ sd = sd*6364136223846793005ULL + 1442695040888963407ULL; xv[ii] = (hg_ptr[ii/2] ? ((double)(sd >> 11)/9007199254740992.0 - 0.5) : 0); }
			for(it = 0; it < 60; it++){
				double nx = 0; for(ii = 0; ii < 2*nc; ii++) nx += xv[ii]*xv[ii]; nx = sqrt(nx); if(!(nx > 0)) break;
				for(ii = 0; ii < 2*nc; ii++){ xv[ii] /= nx; yv[ii] = 0; }
				long f_;
				for(f_ = 0; f_ < hgf_n; f_++){
					long o = hgf_start[f_], L_ = hgf_len[f_], z; double dot = 0;
					for(z = 0; z < L_; z++){ long id = hgo_id[o + z]; if(id < 0 || (size_t)id >= nc || !hg_ptr[id]) continue; double im = 1/sqrt((double)hg_ptr[id]->mass);
						dot += im*(hgo_l[2*(o+z)]*xv[2*id] + hgo_l[2*(o+z)+1]*xv[2*id+1]); }
					dot *= hgf_k[f_];
					for(z = 0; z < L_; z++){ long id = hgo_id[o + z]; if(id < 0 || (size_t)id >= nc || !hg_ptr[id]) continue; double im = 1/sqrt((double)hg_ptr[id]->mass);
						yv[2*id] += dot*im*hgo_l[2*(o+z)]; yv[2*id+1] += dot*im*hgo_l[2*(o+z)+1]; }
				}
				double xy = 0, yy = 0; for(ii = 0; ii < 2*nc; ii++){ xy += xv[ii]*yv[ii]; yy += yv[ii]*yv[ii]; }
				lam = xy; res = sqrt(fabs(yy - xy*xy));
				double *t_ = xv; xv = yv; yv = t_;
			}
			free(xv); free(yv);
			hg_spec_lam = lam; hg_spec_res = res; hg_spec_bound = lmax;
		}
		if((double)Dtime > hg_dtlim) Dtime = (postype)hg_dtlim;
		for(i_ = 0; i_ < VORO_NP(simpar); i_++){
			treevorork4particletype *q = (treevorork4particletype*)(b + i_*ps);
			size_t id = (size_t)PINDX(q);
			if(id >= hg_cap) continue;
			double vf = (double)q->vx*hg_fx[id] + (double)q->vy*hg_fy[id];
			q->ax += hg_fx[id]/q->mass; q->ay += hg_fy[id]/q->mass;
			q->die += (postype)(vf + hg_H[id]);
			if(tcf_sacc_file() && id < tcf_cap && q->ie > 0) tcf_rate[id] += hg_H[id]/(double)q->ie;   /* entropy integral includes HG heat */
			hg_vf_last += vf; hg_H_last += hg_H[id];
		}
	}
	if((double)Dtime < gfs_step_dtmin) gfs_step_dtmin = (double)Dtime;
	if(gfs_tq_on()) tq_ncalls++;
	if(gfs_edge_on()) gfs_edge_call++;
	if(MYID(simpar)<=1){ fprintf(stderr,"[FORCE_DBG] r%d eta wrapper after impl\n", MYID(simpar)); fflush(stderr); }

	if(MYID(simpar)<=1){ fprintf(stderr,"[FORCE_DBG] r%d theta-A before BASICCELL free\n", MYID(simpar)); fflush(stderr); }
	my_free(VORO_BASICCELL(simpar));
	if(MYID(simpar)<=1){ fprintf(stderr,"[FORCE_DBG] r%d theta-B after BASICCELL free\n", MYID(simpar)); fflush(stderr); }
	my_free(VORORK4_TBPP(simpar));
	if(MYID(simpar)<=1){ fprintf(stderr,"[FORCE_DBG] r%d theta wrapper after free\n", MYID(simpar)); fflush(stderr); }
	{
		postype TDtime;
		MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
		Dtime = TDtime;
	}
	if(MYID(simpar)==0){ fprintf(stderr,"[FORCE_DBG] iota wrapper after allreduce\n"); fflush(stderr); }
	return Dtime;
}


/* ================================================================
 *  getAccVoro2D_LagMFM (av_mode=4):
 *    Meshless force computation using MFM effective faces
 *    (Hopkins 2015 GIZMO).  Lagrangian limit: mass per particle is
 *    fixed, face flux carries pressure work only.
 *
 *    Effective face:
 *      A_ij^α = V_i² W_ij · (Ẽ_inv_i^{αβ} dx^β)
 *             + V_j² W_ji · (Ẽ_inv_j^{αβ} dx^β)
 *    (where Ẽ_inv is the inverse of the unnormalized E = Σ W dx⊗dx
 *     already stored in stress.E_inv_xx..yy by updateDenW2Pressure2D_LagMFM.)
 *
 *    Riemann:   pstar, vnstar from hllc_face_2d (first-order).
 *    Force on i: fx += -pstar * A_ij_x + τ·A_ij_x (viscous)
 *    Energy on i: die += -pstar * (vnstar - vn_i) * |A_ij|
 * ================================================================ */
double getAccVoro2D_LagMFM(SimParameters *simpar, postype xmin, postype ymin,
		postype xmax, postype ymax,
		postype OrderOfAccuracy, postype Courant, postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype))
		){
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype Dtime = 1.e10;
	postype nu_phys = GAS_VISCOSITY(simpar);
	float alphavis = GAS_AlphaVis(simpar);
	float betavis  = GAS_BetaVis(simpar);
	float etavis   = GAS_ETAVIS(simpar) * SIMBOX(simpar).x.max / NX(simpar);
	float epsvis   = GAS_EPSVIS(simpar);
#if LAGMFM_PURE_GIZMO
	/* Pure Hopkins(2015): no Monaghan AV / NS viscosity. */
	(void)nu_phys; (void)alphavis; (void)betavis; (void)etavis; (void)epsvis;
#endif

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	int mx, my;
	BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
	BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	CellType *cells = (VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*mx*my));
	mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax, paddingAllTreeParticles);

	postype invcs = 1.0/cellsize;

#if defined(LAGMFV)
	/* Reset MFV mass-flux accumulator at start of every force eval. */
	{ int k; for(k=0; k<nbp; k++) bp[k].stress.mass_flux_accum = 0.0; }
#endif

	int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
	for(iy=0; iy<my; iy++){
		int ix;
		for(ix=0; ix<mx; ix++){
			struct linkedlisttype *li = cells[ix + mx*iy].link;
			while(li){
				treevorostressrk4particletype *ibp = (treevorostressrk4particletype*)li;
				int iwall = (PINDX((treevorork4particletype*)ibp) == MAX_INDEX);
				if(iwall || IS_FLAG_ON(li,BoundaryGhostflag)){
					li = li->next;
					continue;
				}

				postype ibp_vx  = ibp->vx;
				postype ibp_vy  = ibp->vy;
				postype ibp_cs  = ibp->csound;
				postype ibp_P   = ibp->pressure;
				postype ibp_rho = ibp->den;
				postype Vi      = ibp->volume;
				postype h_i     = ibp->stress.h_mfm;
				if(h_i <= 0) h_i = LAGMFM_ETA * sqrt(Vi > 0 ? Vi : cellsize*cellsize);
				postype Ei_xx = ibp->stress.E_inv_xx;
				postype Ei_xy = ibp->stress.E_inv_xy;
				postype Ei_yx = ibp->stress.E_inv_yx;
				postype Ei_yy = ibp->stress.E_inv_yy;

				postype fx=0, fy=0, die=0;
				postype dt_cell = 1.e10;

				/* Neighbor search: up to 2 max(h_i, h_neighbor_max).
				 * We don't know h_j a priori, so use 2*h_i and validate via
				 * per-pair 2*max(h_i,h_j) cutoff inside the inner loop. */
				postype two_h_i = 2.0 * h_i;
				int nbr = (int)ceil(two_h_i * invcs) + 1;
				if(nbr < 1) nbr = 1;
				if(nbr > 10) nbr = 10;

				int jy;
				for(jy = iy-nbr; jy <= iy+nbr; jy++){
					if(jy < 0 || jy >= my) continue;
					int jx;
					for(jx = ix-nbr; jx <= ix+nbr; jx++){
						if(jx < 0 || jx >= mx) continue;
						struct linkedlisttype *lj = cells[jx + mx*jy].link;
						while(lj){
							treevorostressrk4particletype *jbp =
									(treevorostressrk4particletype*)lj;
							if((treevorork4particletype*)jbp == (treevorork4particletype*)ibp){
								lj = lj->next;
								continue;
							}
							int jwall = (PINDX((treevorork4particletype*)jbp) == MAX_INDEX);

							postype dx = jbp->x - ibp->x;
							postype dy = jbp->y - ibp->y;
							postype r2 = dx*dx + dy*dy;

							postype h_j = jbp->stress.h_mfm;
							if(h_j <= 0) h_j = h_i;
							postype h_max = h_i > h_j ? h_i : h_j;
							postype two_h_max = 2.0*h_max;
							if(r2 > two_h_max*two_h_max){
								lj = lj->next;
								continue;
							}
							postype r = sqrt(r2);

							/* Effective face area vector (antisymmetric in i↔j).
							 * Wall-ghost: cannot use j's E_inv (unmirrored), so
							 * use only the i-side contribution.
							 *
							 * GIZMO compute_finitevol_faces.h: when |V_i-V_j|/min(V)/NUMDIMS > 1.25
							 * (2D threshold = 2.5), switch from (V_i, V_j) weights to a
							 * symmetric volume-centered weight to keep the face well-behaved
							 * at sharp density jumps.  At KH ρ=1↔2 contact, V ratio ≈ 2 —
							 * on the edge, triggers intermittently. */
							postype W_ij = mfm_W_wendland2d(r, h_i);
							postype W_ji = mfm_W_wendland2d(r, h_j);
							postype Vj   = jbp->volume;
							postype wt_i, wt_j;
							{
								postype Vmin = (Vi < Vj) ? Vi : Vj;
								postype Vabs = Vi - Vj; if(Vabs < 0) Vabs = -Vabs;
								/* 2D: NUMDIMS=2, threshold 2*1.25 = 2.5 */
								if(Vmin > 0 && (Vabs/Vmin) > (postype)2.5){
									postype denom = Vi*W_ij + Vj*W_ji;
									postype wcomm = (denom > 0) ? (Vi*Vj*(W_ij+W_ji)/denom) : Vi;
									wt_i = wt_j = wcomm;
								} else {
									wt_i = Vi;
									wt_j = Vj;
								}
							}
							postype Vi2 = wt_i;
							postype Vj2 = wt_j;
							postype Ej_xx = jbp->stress.E_inv_xx;
							postype Ej_xy = jbp->stress.E_inv_xy;
							postype Ej_yx = jbp->stress.E_inv_yx;
							postype Ej_yy = jbp->stress.E_inv_yy;

							postype Ax_i = Vi2 * W_ij * (Ei_xx*dx + Ei_xy*dy);
							postype Ay_i = Vi2 * W_ij * (Ei_yx*dx + Ei_yy*dy);
							postype Ax_j, Ay_j;
							if(jwall){
								Ax_j = 0;
								Ay_j = 0;
							} else {
								Ax_j = Vj2 * W_ji * (Ej_xx*dx + Ej_xy*dy);
								Ay_j = Vj2 * W_ji * (Ej_yx*dx + Ej_yy*dy);
							}
							postype Ax = Ax_i + Ax_j;
							postype Ay = Ay_i + Ay_j;
							postype Amag2 = Ax*Ax + Ay*Ay;
							if(Amag2 < 1e-30){
								lj = lj->next;
								continue;
							}
							postype Amag = sqrt(Amag2);
							postype nx_hat = Ax / Amag;
							postype ny_hat = Ay / Amag;

							postype jbp_vx = jbp->vx;
							postype jbp_vy = jbp->vy;
							postype jbp_cs = jbp->csound;
							postype jbp_P  = jbp->pressure;
							postype jbp_rho = jbp->den;

							/* GIZMO default (hydro_core_meshless.h): equal-weight midpoint,
							 * s_star_ij = 0 → s_i = s_j = 0.5·r, face at arithmetic mid.
							 * (Hopkins 2015 Eq. 20 h-weighted form is only used in special
							 * branches; main MFM loop uses 0.5.) */
							postype hdx = 0.5*dx, hdy = 0.5*dy;

							postype pL_rec  = ibp_P   + ibp->stress.dPdx*hdx   + ibp->stress.dPdy*hdy;
							postype pR_rec  = jbp_P   - jbp->stress.dPdx*hdx   - jbp->stress.dPdy*hdy;
							postype vxL_rec = ibp_vx  + ibp->stress.gUxx*hdx   + ibp->stress.gUxy*hdy;
							postype vxR_rec = jbp_vx  - jbp->stress.gUxx*hdx   - jbp->stress.gUxy*hdy;
							postype vyL_rec = ibp_vy  + ibp->stress.gUyx*hdx   + ibp->stress.gUyy*hdy;
							postype vyR_rec = jbp_vy  - jbp->stress.gUyx*hdx   - jbp->stress.gUyy*hdy;
							postype rhoL_rec= ibp_rho + ibp->stress.dRhodx*hdx + ibp->stress.dRhody*hdy;
							postype rhoR_rec= jbp_rho - jbp->stress.dRhodx*hdx - jbp->stress.dRhody*hdy;

							/* GIZMO reconstruct_face_states (reimann.h:188-237): asymmetric
							 * pair-wise limiter with fac_minmax=0.5, fac_meddev=0.375.
							 * Our convention: L=from i, R=from j.  GIZMO convention opposite
							 * (Q_R from i, Q_L from j) — role-swap the clamp logic below.
							 *
							 *   Qmed = 0.5*(Q_i + Q_j)
							 *   Qmax_eff = Qmax + 0.5·(Qmax-Qmin)  Qmin_eff = Qmin - 0.5·(Qmax-Qmin)
							 *   Qmed_max = min(Qmed + 0.375·(Qmax-Qmin), Qmax_eff)
							 *   Qmed_min = max(Qmed - 0.375·(Qmax-Qmin), Qmin_eff)
							 *
							 * For Q_i<Q_j:
							 *   L (from i, low side): clamp to [Qmin_eff, Qmed_max]
							 *   R (from j, high side): clamp to [Qmed_min, Qmax_eff]
							 * Q_i>Q_j: swap roles. */
#define GIZMO_PAIRLIM(QL, QR, Qi, Qj) do {                                   \
    if((Qi) != (Qj)) {                                                        \
        postype _Qmed = (postype)0.5 * ((Qi) + (Qj));                         \
        postype _Qmin, _Qmax;                                                 \
        if((Qi) < (Qj)) { _Qmin = (Qi); _Qmax = (Qj); }                       \
        else            { _Qmin = (Qj); _Qmax = (Qi); }                       \
        postype _dQ  = _Qmax - _Qmin;                                         \
        postype _f1  = (postype)0.5   * _dQ;                                  \
        postype _f2  = (postype)0.375 * _dQ;                                  \
        postype _Qmax_eff = _Qmax + _f1;                                      \
        postype _Qmin_eff = _Qmin - _f1;                                      \
        postype _Qmed_max = _Qmed + _f2;                                      \
        postype _Qmed_min = _Qmed - _f2;                                      \
        if(_Qmed_max > _Qmax_eff) _Qmed_max = _Qmax_eff;                      \
        if(_Qmed_min < _Qmin_eff) _Qmed_min = _Qmin_eff;                      \
        if((Qi) < (Qj)) {                                                      \
            if((QL) < _Qmin_eff) (QL) = _Qmin_eff;                            \
            if((QL) > _Qmed_max) (QL) = _Qmed_max;                            \
            if((QR) > _Qmax_eff) (QR) = _Qmax_eff;                            \
            if((QR) < _Qmed_min) (QR) = _Qmed_min;                            \
        } else {                                                               \
            if((QL) > _Qmax_eff) (QL) = _Qmax_eff;                            \
            if((QL) < _Qmed_min) (QL) = _Qmed_min;                            \
            if((QR) < _Qmin_eff) (QR) = _Qmin_eff;                            \
            if((QR) > _Qmed_max) (QR) = _Qmed_max;                            \
        }                                                                      \
    }                                                                          \
} while(0)
							if(!jwall){
								GIZMO_PAIRLIM(pL_rec,   pR_rec,   ibp_P,   jbp_P);
								GIZMO_PAIRLIM(vxL_rec,  vxR_rec,  ibp_vx,  jbp_vx);
								GIZMO_PAIRLIM(vyL_rec,  vyR_rec,  ibp_vy,  jbp_vy);
								GIZMO_PAIRLIM(rhoL_rec, rhoR_rec, ibp_rho, jbp_rho);
							}
#undef GIZMO_PAIRLIM

							/* Wall ghost: j-side gradients unmirrored → use
							 * central state on both sides to avoid artifacts. */
							if(jwall){
								pL_rec   = ibp_P;
								vxL_rec  = ibp_vx;
								vyL_rec  = ibp_vy;
								rhoL_rec = ibp_rho;
								pR_rec   = jbp_P;
								vxR_rec  = jbp_vx;
								vyR_rec  = jbp_vy;
								rhoR_rec = jbp_rho;
							}
							/* Safety floor for reconstructed ρ, P (HLLC requires >0). */
							if(rhoL_rec <= 0) rhoL_rec = ibp_rho;
							if(rhoR_rec <= 0) rhoR_rec = jbp_rho;
							if(pL_rec   <= 0) pL_rec   = ibp_P;
							if(pR_rec   <= 0) pR_rec   = jbp_P;

							/* Project reconstructed velocities to the normal */
							postype vnL_lab = vxL_rec*nx_hat + vyL_rec*ny_hat;
							postype vnR_lab = vxR_rec*nx_hat + vyR_rec*ny_hat;

							postype pL = pL_rec;
							postype pR = pR_rec;
							postype rhoL = rhoL_rec;
							postype rhoR = rhoR_rec;
							/* Sound speed uses reconstructed ρ, P */
							postype cL = sqrt(Gamma * pL / rhoL);
							postype cR = sqrt(Gamma * pR / rhoR);

							/* GIZMO default (hydro_core_meshless.h:86-102):
							 * equal-weight face velocity w = 0.5·(v_i + v_j). */
							postype wx = 0.5*(ibp_vx + jbp_vx);
							postype wy = 0.5*(ibp_vy + jbp_vy);
							postype wn = wx*nx_hat + wy*ny_hat;

							postype pstar, vnstar_lab;
							if(jwall){
								/* Wall face: bypass Riemann, use central
								 * pressure and v_n = 0 (wall enforces no
								 * normal flux). */
								pstar       = 0.5*(pL + pR);
								vnstar_lab  = 0;
							} else {
								hllc_face_2d_rest_frame(rhoL, pL, vnL_lab, cL,
								                        rhoR, pR, vnR_lab, cR,
								                        wn, Gamma,
								                        &pstar, &vnstar_lab);
							}

#if !LAGMFM_PURE_GIZMO
							/* NON-GIZMO: Monaghan AV (SPH-style) add-on.
							 * Disabled in pure Hopkins(2015) mode. */
							if(alphavis > 0 && !jwall){
								postype uij_x = jbp_vx - ibp_vx;
								postype uij_y = jbp_vy - ibp_vy;
								postype er_x  = dx/r;
								postype er_y  = dy/r;
								postype rvel  = er_x*uij_x + er_y*uij_y;
								if(rvel < 0){
									postype wcomp = sqrt(ibp->w2) +
									                sqrt(((treevorork4particletype*)jbp)->w2);
									postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
									postype drampScale  = r/scaleFactor;
									postype mu = rvel /* Voronoi-AV-reform: no length regularization */;
									postype meanden = 0.5*(ibp_rho + jbp_rho);
									postype meanCs  = 0.5*(ibp_cs  + jbp_cs);
									pstar += (-alphavis*meanCs*mu + betavis*mu*mu)*meanden;
								}
							}
#endif
							/* NS viscous traction (GIZMO-absent; gated).
							 * τ in stress.tauxx..tauyy is -νρS. */
							postype tau_A_x = 0, tau_A_y = 0;
#if !LAGMFM_PURE_GIZMO
							if(!jwall){
								postype txx_f = 0.5*(ibp->stress.tauxx + jbp->stress.tauxx);
								postype txy_f = 0.5*(ibp->stress.tauxy + jbp->stress.tauxy);
								postype tyy_f = 0.5*(ibp->stress.tauyy + jbp->stress.tauyy);
								/* τ_NS = -τ_code */
								tau_A_x = -(txx_f*Ax + txy_f*Ay);
								tau_A_y = -(txy_f*Ax + tyy_f*Ay);
							}
#endif

							/* Force: -pstar * A  (pure MFM); viscous adds zero when pure. */
							fx += -pstar * Ax + tau_A_x;
							fy += -pstar * Ay + tau_A_y;

#if defined(LAGMFV)
							/* MFV (Hopkins 2015 Sec 2.5, GIZMO hydro_core_meshless.h):
							 *   m_flux = rho_up * (v_n_star - w_n) * Adotn    (+ out of i)
							 *   rho_up, vt_up selected upwind by contact direction.
							 *   Force: -m_flux * v_up   (momentum advection)
							 *   Energy: -m_flux * (e_tot + P_up/rho_up)  (enthalpy flux)
							 * Mass change accumulated in stress.mass_flux_accum;
							 * integrator updates particle mass each half-step.   */
							if(!jwall){
								postype vn_rel = vnstar_lab - wn;
								int upL_mfv = (vn_rel >= 0);
								postype rho_up = upL_mfv ? rhoL_rec : rhoR_rec;
								postype p_up   = upL_mfv ? pL : pR;
								postype vx_up_mfv = upL_mfv ? vxL_rec : vxR_rec;
								postype vy_up_mfv = upL_mfv ? vyL_rec : vyR_rec;

								postype Adotn = Ax*nx_hat + Ay*ny_hat;
								postype m_flux = rho_up * vn_rel * Adotn;

								fx += -m_flux * vx_up_mfv;
								fy += -m_flux * vy_up_mfv;

								postype e_tot = 0.5*(vx_up_mfv*vx_up_mfv + vy_up_mfv*vy_up_mfv)
								               + p_up/((Gamma - 1.0)*rho_up);
								die += -m_flux * (e_tot + p_up/rho_up);

								ibp->stress.mass_flux_accum += -m_flux;
							}
#endif

							/* ============================================
							 * D-B4 Full face velocity + viscous heating.
							 *
							 * Newton 3rd law & energy conservation check:
							 *
							 *   Per face,  A_ij is antisymmetric (A_ji = -A_ij)
							 *   and τ_face, pstar are symmetric in (i,j).
							 *   Hence  F_ij + F_ji = 0  (checked).
							 *
							 *   Total viscous heating per face:
							 *     die_i + die_j = (τ_NS·A_ij)·(v_face - v_i)
							 *                   + (τ_NS·A_ji)·(v_face - v_j)
							 *                   = (τ_NS·A_ij)·(v_j - v_i).
							 *   Total kinetic-energy change per face:
							 *     dKE_i + dKE_j = F_ij·v_i + F_ji·v_j
							 *                   = (τ_NS·A_ij)·(v_i - v_j).
							 *   Sum: 0  ⇒ energy conserved.
							 *   Since τ_NS is symmetric negative semidefinite
							 *   on the shear it dissipates v_j-v_i in the
							 *   correct direction — viscous heating ≥ 0.
							 *
							 * Full face velocity (normal + tangential):
							 *   v_face = vnstar · n̂ + v_t · t̂
							 * Hopkins (2015) §2.5: tangential velocity is
							 * **passively advected** across the contact wave.
							 * Upwind rule based on the sign of (vnstar - wn):
							 *   (vnstar - wn) > 0  ⇒ contact moves toward R,
							 *                       L-state passes through → v_t_L
							 *   (vnstar - wn) < 0  ⇒ contact moves toward L,
							 *                       R-state passes through → v_t_R
							 * ============================================ */
							postype tx_hat = -ny_hat;
							postype ty_hat =  nx_hat;
							postype vtL = vxL_rec*tx_hat + vyL_rec*ty_hat;
							postype vtR = vxR_rec*tx_hat + vyR_rec*ty_hat;
							/* wn already defined above; upwind by contact direction
							 * in the face rest frame. */
							postype vt_up = (vnstar_lab - wn) >= 0 ? vtL : vtR;
							postype vfx = vnstar_lab * nx_hat + vt_up * tx_hat;
							postype vfy = vnstar_lab * ny_hat + vt_up * ty_hat;

							/* D-B6 Wall slip option:
							 *   Free-slip (default):  v_face keeps the
							 *     particle's tangential velocity
							 *     (approximated by the avg above, which for
							 *      a wall pair is 0.5*v_i_tan in the limit
							 *      v_j = mirrored v_i_tan).
							 *   No-slip (LAGMFM_WALL_NOSLIP=1): zero the
							 *     face velocity entirely at walls. */
							if(jwall){
#if LAGMFM_WALL_NOSLIP
								vfx = 0;
								vfy = 0;
#else
								/* Free-slip: tangential = v_i, normal = 0 */
								postype vti = ibp_vx*tx_hat + ibp_vy*ty_hat;
								vfx = vti*tx_hat;
								vfy = vti*ty_hat;
#endif
							}

							die += -pstar*((vfx - ibp_vx)*Ax + (vfy - ibp_vy)*Ay);
							/* Viscous heating: (τ_NS·A) · (v_face - v_i) */
							die += tau_A_x*(vfx - ibp_vx) + tau_A_y*(vfy - ibp_vy);

							/* Signal velocity for CFL (Hopkins 2015, eq. 9).
							 *   dt_i = C · h_i / vsig,
							 *   vsig = c_L + c_R - min(0, (v_j-v_i)·r̂)
							 * Uses kernel bandwidth h_i as length scale. */
							postype dvx = jbp_vx - ibp_vx;
							postype dvy = jbp_vy - ibp_vy;
							postype er_x = dx/r;
							postype er_y = dy/r;
							postype VdotR = dvx*er_x + dvy*er_y;
							postype vsig = cL + cR - MIN(0, VdotR);
							postype dt_face = Courant * h_i / vsig;
							if(dt_face < dt_cell) dt_cell = dt_face;

#if !LAGMFM_PURE_GIZMO
							/* CD10 shock switch tracking (non-GIZMO). */
							if(!jwall){
								if(vsig > ibp->stress.vsig_max)
									ibp->stress.vsig_max = vsig;
							}
#endif

							lj = lj->next;
						}
					}
				}

#if !LAGMFM_PURE_GIZMO
				/* Viscous CFL (NS + CD10; GIZMO-absent, gated). */
				{
					postype h_eff = sqrt(Vi);
					postype nu_cd = ibp->stress.alpha_cd * h_eff * ibp_cs;
					postype nu_eff = fmax(nu_phys, nu_cd);
					if(nu_eff > 0){
						postype dt_visc = 0.5*h_eff*h_eff/nu_eff;
						if(dt_visc < dt_cell) dt_cell = dt_visc;
					}
				}
#endif

				ibp->die = die / ibp->mass;  /* store as rate per unit mass? */
				/* Actually die accumulated above is m_i * de/dt.  For consistency
				 * with the blend code convention (ibp_rk4->die = dE/dt), store
				 * die as the total rate, not per-mass.  We must match the
				 * integrator which does: ie += die*Dtime. */
				ibp->die = die;
				ibp->ax  = fx / ibp->mass;
				ibp->ay  = fy / ibp->mass;
				ibp->dt  = dt_cell;
				if(dt_cell < Dtime) Dtime = dt_cell;

				li = li->next;
			}
		}
	}

	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
	{
		postype TDtime;
		MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
		Dtime = TDtime;
	}
	return Dtime;
}

double getAccVoro2D_rt(SimParameters *simpar, postype xmin, postype ymin,
		postype xmax, postype ymax,
		postype OrderOfAccuracy, postype Courant, postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype))
		){
	postype boxsize = VORO_CELLBOX(simpar);
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int isave = -1;
    postype Dtime = 1.e10;
	int nx=NX(simpar);


	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	postype cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar);
    int mx, my;
    BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
    BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
	mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);


	postype dtold = GAS_dtold(simpar);



	float alphavis = GAS_AlphaVis(simpar);
	float betavis = GAS_BetaVis(simpar);
	float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
	float epsvis = GAS_EPSVIS(simpar);


//	DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);

    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=0;iy<my;iy++){
    	int mp=1000;
    	Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
    	postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
        for(ix=0;ix<mx;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			int i;
            for(i=0;i<np;i++){

                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i); 
                center.csound = p[i].csound;
                center.w2 = p[i].w2;

				treevorork4particletype *ibp = p[i].bp;
				ibp->dt = 1.e10; // initialization of dt

                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
//                ibp->volume = Area2DPolygon(vorocorner, mp);
				// to use registers
				postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
				ibp_vx = ibp->vx;
				ibp_vy = ibp->vy;
				ibp_csound = ibp->csound;
				ibp_pressure = ibp->pressure;
				ibp_den = ibp->den;


				Voro2D_Corner *tmp,*tmp2;
				tmp = vorocorner;
				int j;
				double die, dte,dke,fx,fy;
				die = dte = dke = fx = fy = 0;
				do {
					/*
					if(tmp->upperrelated <0) {
						tmp = vorocorner;
						do { 
							treevorork4particletype *jbp = (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
							DEBUGPRINT("P%d has p%d :neighbor %d : %g %g : %g %g\n", 
									MYID(simpar), center.indx, tmp->upperrelated, ibp->x,ibp->y,jbp->x, jbp->y);
						} while(tmp);
						exit(99);
					}
					*/
					if(tmp->upperrelated >=0)
					{
						// to use registers
						treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
						postype jbp_vx = jbp->vx;
						postype jbp_vy = jbp->vy;
						postype jbp_csound = jbp->csound;
                        postype jbp_pressure = jbp->pressure;
//						postype jbp_pressure = getNeighborPressure(jbp,ibp); // to consider the wallflag
						postype jbp_den = jbp->den;

						Voro2D_point line; 
	                    tmp2 = tmp->upperlink; 
						line.x = tmp2->x - tmp->x; 
						line.y = tmp2->y - tmp->y;
	                    Voro2D_point dS = voro2D_norm(&line); 
						postype facearea = Vec2DLength(tmp,tmp2); 
						dS.x = facearea*dS.x;
						dS.y = facearea*dS.y;

						/*
						postype pi = getPressure2D(
								ibp_pressure,jbp_pressure, tmp, neighwork, OrderOfAccuracy);
								*/
						postype pi = VoroRK4_Pressure2D(
								ibp_pressure,jbp_pressure, tmp, neighwork, simpar, LAG);

						Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

						postype dramp = sqrt(Vec2DDotP(&dr, &dr));
						Voro2D_point er;
						er.x = dr.x/dramp;
						er.y = dr.y/dramp;

						Voro2D_point uij,ua,ub;


						uij.x = (jbp_vx - ibp_vx); 
						uij.y = (jbp_vy - ibp_vy);

						// artificial bulk viscosity
						{
							postype rvel = Vec2DDotP(&er, &uij);
							if(rvel<0){
								postype wcomp = 0.5*(sqrt(ibp->w2hydro)+sqrt(jbp->w2hydro));
								postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
								postype drampScale = dramp/wcomp;
								postype mu = rvel /* Voronoi-AV-reform: no length regularization */;
								postype meanden = 0.5*(ibp_den + jbp_den);
								postype meanCsound = 0.5*(ibp_csound + jbp_csound);

								/*
								if(center.indx == 58443){
									DEBUGPRINT("p(%ld,%ld): has viscos: wcomp= %g drampScale= %g mu= %g meanden= %g meanCsound= %g :: pi= %g  dpi= %g\n", 
											PINDX(ibp)%256, PINDX(ibp)/256, wcomp, drampScale, mu, meanden, meanCsound, 
											pi, (-alphavis * meanCsound*mu + betavis*mu*mu)*meanden);
								}
								*/

								pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
							}
						}

						// for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);

						die += -pi * Vec2DDotP(&uradix_ui,&dS);

						// for the total energy
						ua.x = ibp_vx + uradix_ui.x;
						ua.y = ibp_vy + uradix_ui.y;
						dte += -pi*Vec2DDotP(&ua, &dS);

						// for the kinetic energy 
						ub.x = ibp_vx;
						ub.y = ibp_vy;
						dke += -pi*Vec2DDotP(&ub, &dS);

						// for the force 
						fx += -pi * dS.x;
						fy += -pi * dS.y;

						/*
						if(center.indx == 58443)
						{
							DEBUGPRINT("p(%ld,%ld): xy= %g %g has ax/y= %g %g pij= %g dS= %g %g pi= %g : p%lu rj= %g %g pj= %g :: jindx= (%ld,%ld)\n", 
									PINDX(ibp)%256, PINDX(ibp)/256,
									ibp->x,ibp->y,
									fx/ibp->mass,fy/ibp->mass,
									pi,dS.x,dS.y, ibp_pressure, PINDX(jbp),jbp->x, jbp->y, jbp_pressure, PINDX(jbp)%256, PINDX(jbp)/256);
						}
						*/

						Voro2D_point dv;
						dv.x = (jbp_vx - ibp_vx); 
						dv.y = (jbp_vy - ibp_vy); 

						postype VdotR = Vec2DDotP(&dv,&er);
						postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
						postype dramp_cfl = fmax(dramp, 0.25*sqrt(ibp->volume));
						postype dt = 2*Courant*dramp_cfl/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d p(%ld,%ld) x/y= %g %g : %g %g : %g %g : %g : dv.xy= %g %g for p(%ld,%ld)\n",
                                    MYID(simpar), i, PINDX(p+i)%nx, PINDX(p+i)/nx,
									p[i].x, p[i].y, dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y,
									PINDX(jbp)%nx, PINDX(jbp)/nx);
//							DEBUGGING(simpar);
							dt = 1e-10;
                        }
                        ibp->dt = MIN(ibp->dt,dt);
                        if(dt < Dtime){
                            Dtime = dt;
                            // isave = id;
                        }
                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
               	ibp->die = die;

                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;

				/*
				if(center.indx == 58443) 
				{
					DEBUGPRINT("+p%ld has ax/y= %g %g xy= %g %g : vxy= %g %g ::: rho= %g\n", PINDX(ibp),ibp->ax,ibp->ay, ibp->x, ibp->y,
							ibp->vx, ibp->vy, ibp->den);
					exit(0);
				}
				*/

                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g with p= %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin,
							ibp->pressure);
					exit(0);
                }

            }
            free(p);free(neighbors); free(neighwork);

        }
        free(vorocorner);
    }
	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
	{
		postype TDtime;
        MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
        Dtime = TDtime;
	}
	return Dtime;
}
double getAccVoro2D(SimParameters *simpar, postype xmin, postype ymin, 
		postype xmax, postype ymax,
		postype OrderOfAccuracy, postype Courant, postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	postype boxsize = VORO_CELLBOX(simpar);
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int isave = -1;
    postype Dtime = 1.e10;
	int nx = NX(simpar);


	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	postype cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar);
    int mx, my;
    BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
    BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
	mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);


	postype dtold = GAS_dtold(simpar);



	float alphavis = GAS_AlphaVis(simpar);
	float betavis = GAS_BetaVis(simpar);
	float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
	float epsvis = GAS_EPSVIS(simpar);


//	DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);

    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=1;iy<my-1;iy++){
    	int mp=1000;
    	Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
    	postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			int i;
            for(i=0;i<np;i++){

                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i); 
                center.csound = p[i].csound;
                center.w2 = p[i].w2;

				treevorork4particletype *ibp = p[i].bp;
				ibp->dt = 1.e10; // initialization of dt

                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
//                ibp->volume = Area2DPolygon(vorocorner, mp);
				// to use registers
				postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
				ibp_vx = ibp->vx;
				ibp_vy = ibp->vy;
				ibp_csound = ibp->csound;
				ibp_pressure = ibp->pressure;
				ibp_den = ibp->den;


				Voro2D_Corner *tmp,*tmp2;
				tmp = vorocorner;
				int j;
				double die, dte,dke,fx,fy;
				die = dte = dke = fx = fy = 0;
				do {
					if(tmp->upperrelated >=0)
					{
						// to use registers
						treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
						postype jbp_vx = jbp->vx;
						postype jbp_vy = jbp->vy;
						postype jbp_csound = jbp->csound;
						postype jbp_pressure = jbp->pressure;
						postype jbp_den = jbp->den;

						Voro2D_point line; 
	                    tmp2 = tmp->upperlink; 
						line.x = tmp2->x - tmp->x; 
						line.y = tmp2->y - tmp->y;
	                    Voro2D_point dS = voro2D_norm(&line); 
						postype facearea = Vec2DLength(tmp,tmp2); 
						dS.x = facearea*dS.x;
						dS.y = facearea*dS.y;

						/*
						postype pi = getPressure2D(
								ibp_pressure,jbp_pressure, tmp, neighwork, OrderOfAccuracy);
								*/
						postype pi = VoroRK4_Pressure2D(
								ibp_pressure,jbp_pressure, tmp, neighwork, simpar, RT);

						Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

						postype dramp = sqrt(Vec2DDotP(&dr, &dr));
						Voro2D_point er;
						er.x = dr.x/dramp;
						er.y = dr.y/dramp;

						Voro2D_point uij,ua,ub;
						uij.x = (jbp_vx - ibp_vx); 
						uij.y = (jbp_vy - ibp_vy);

						{
							postype rvel = Vec2DDotP(&er, &uij);
							if(rvel<0){
								postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
								postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
								postype drampScale = dramp/scaleFactor;
								postype mu = rvel /* Voronoi-AV-reform: no length regularization */;
								postype meanden = 0.5*(ibp_den + jbp_den);
								postype meanCsound = 0.5*(ibp_csound + jbp_csound);
								pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
							}
						}

						// for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

						// for the total energy
						ua.x = ibp_vx + uradix_ui.x;
						ua.y = ibp_vy + uradix_ui.y;
						dte += -pi*Vec2DDotP(&ua, &dS);

						// for the kinetic energy
						ub.x = ibp_vx;
						ub.y = ibp_vy;
						dke += -pi*Vec2DDotP(&ub, &dS);

						// for the force
						fx += -pi * dS.x;
						fy += -pi * dS.y;

						/*
						if(center.indx == 130816)
						{
							DEBUGPRINT("p(%ld,%ld): xy= %g %g has ax/y= %g %g pij= %g dS= %g %g pi= %g : p%lu rj= %g %g pj= %g :: jindx= (%ld,%ld)\n", 
									PINDX(ibp)/256, PINDX(ibp)%256,
									ibp->x,ibp->y,
									fx/ibp->mass,fy/ibp->mass,
									pi,dS.x,dS.y, ibp_pressure, PINDX(jbp),jbp->x, jbp->y, jbp_pressure, PINDX(jbp)/256, PINDX(jbp)%256);
						}
						*/

						Voro2D_point dv;
						dv.x = (jbp_vx - ibp_vx); 
						dv.y = (jbp_vy - ibp_vy); 

						postype VdotR = Vec2DDotP(&dv,&er);
						postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
						postype dramp_cfl = fmax(dramp, 0.25*sqrt(ibp->volume));
						postype dt = 2*Courant*dramp_cfl/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d %ld x/y= %g %g : %g %g : %g %g : %g : dv.xy= %g %g\n",
                                    MYID(simpar), i, PINDX(p+i), p[i].x, p[i].y, dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y);
                            dt = 1e-10;
                        }
                        ibp->dt = MIN(ibp->dt,dt);
                        if(dt < Dtime){
                            Dtime = dt;
                            // isave = id;
                        }
                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
               	ibp->die = die;
                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;

				/*
				if(center.indx == 130816) 
				{
					DEBUGPRINT("+p%ld has ax/y= %g %g xy= %g %g : vxy= %g %g ::: rho= %g\n", PINDX(ibp),ibp->ax,ibp->ay, ibp->x, ibp->y,
							ibp->vx, ibp->vy, ibp->den);
							
				}
				*/

                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g with p= %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin,
							ibp->pressure);
					exit(0);
                }

            }
            free(p);free(neighbors); free(neighwork);

        }
        free(vorocorner);
    }
	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
	{
		postype TDtime;
        MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN,  MPI_COMM(simpar));
        Dtime = TDtime;
	}
	return Dtime;
}

/* ================================================================
 *  getAccVoro2D_kNN:
 *  k-NN based force computation. Reuses pre-built tree.
 *  No cell grid needed.
 * ================================================================ */
double getAccVoro2D_kNN(SimParameters *simpar,
		postype xmin, postype ymin, postype xmax, postype ymax,
		postype OrderOfAccuracy, postype Courant, postype Gamma,
		TStruct *tree, TPtlStruct *ptl, int nptl,
		int (*skipBP)(treevorork4particletype *)){
	postype boxsize = VORO_CELLBOX(simpar);
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype Dtime = 1.e10;

	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	postype dtold = GAS_dtold(simpar);

	float alphavis = GAS_AlphaVis(simpar);
	float betavis = GAS_BetaVis(simpar);
	float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
	float epsvis = GAS_EPSVIS(simpar);

	int K_START = 32;
	size_t p_size_a = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw_a = (char*)VORORK4_TBP(simpar);

#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
	for(int ii=0;ii<nbp;ii++){
		int skip_this = 0;
		if(skipBP){
			treevorork4particletype *bpi_chk = (treevorork4particletype*)(bp_raw_a + ii*p_size_a);
			if(skipBP(bpi_chk)) skip_this = 1;
		}
		if(!skip_this){

		int mp = 1000;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		Voro2D_point neighbors[MAX_NUM_NEAR];
		Voro2D_point neighwork[MAX_NUM_NEAR];
		memset(neighwork, 0, sizeof(neighwork));
		PosType maxr;
		int K = K_START;
		int nfound, ip;

		/* Adaptive kNN with retry (no goto/continue — Intel OpenMP safe) */
		int knn_ok = 0;
		while(!knn_ok){
			nfound = searchKNN2D(ptl+ii, K, tree, neighbors, &maxr);
			if(nfound < 3 || maxr < 1e-20) break;

			treevorork4particletype *ibp_tmp = (treevorork4particletype*)(ptl[ii].bp);
			Voro2D_point center;
			center.x = ibp_tmp->x;
			center.y = ibp_tmp->y;
			center.indx = PINDX(ibp_tmp);
			center.csound = ibp_tmp->csound;
			center.w2 = ibp_tmp->w2;

			ibp_tmp->dt = 1.e10;
			ip = Voro2D_FindVC(&center, neighbors, neighwork, nfound, vorocorner, mp, boxsize);

			/* Completeness check */
			Voro2D_Corner *vtmp = vorocorner;
			postype maxvert2 = 0;
			do {
				postype dx = vtmp->x - center.x;
				postype dy = vtmp->y - center.y;
				postype d2 = dx*dx + dy*dy;
				if(d2 > maxvert2) maxvert2 = d2;
				vtmp = vtmp->upperlink;
			} while(vtmp != vorocorner);
			if(sqrt(maxvert2) > 0.95*maxr && K < MAX_NUM_NEAR){
				K = MIN(K*2, MAX_NUM_NEAR);
			} else {
				knn_ok = 1;
			}
		}

		if(knn_ok){
		treevorork4particletype *ibp = (treevorork4particletype*)(ptl[ii].bp);

		postype ibp_vx = ibp->vx;
		postype ibp_vy = ibp->vy;
		postype ibp_csound = ibp->csound;
		postype ibp_pressure = ibp->pressure;
		postype ibp_den = ibp->den;

		Voro2D_Corner *tmp, *tmp2;
		tmp = vorocorner;
		double die, dte, dke, fx, fy;
		die = dte = dke = fx = fy = 0;
		do {
			if(tmp->upperrelated >= 0){
				treevorork4particletype *jbp = (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
				postype jbp_vx = jbp->vx;
				postype jbp_vy = jbp->vy;
				postype jbp_csound = jbp->csound;
				postype jbp_pressure = jbp->pressure;
				postype jbp_den = jbp->den;

				Voro2D_point line;
				tmp2 = tmp->upperlink;
				line.x = tmp2->x - tmp->x;
				line.y = tmp2->y - tmp->y;
				Voro2D_point dS = voro2D_norm(&line);
				postype facearea = Vec2DLength(tmp, tmp2);
				dS.x = facearea*dS.x;
				dS.y = facearea*dS.y;

				postype pi = VoroRK4_Pressure2D(
						ibp_pressure, jbp_pressure, tmp, neighwork, simpar, RT);

				Voro2D_point dr = EunhaVec2DSub(jbp, ibp);
				postype dramp = sqrt(Vec2DDotP(&dr, &dr));
				/* Skip degenerate face: particles at same position */
				if(dramp < 1e-20){
					tmp = tmp->upperlink;
					continue;
				}
				Voro2D_point er;
				er.x = dr.x/dramp;
				er.y = dr.y/dramp;

				Voro2D_point uij;
				uij.x = (jbp_vx - ibp_vx);
				uij.y = (jbp_vy - ibp_vy);

				/* Skip AV for boundary ghost neighbors — the velocity
				   difference from wall reflection is not a physical shock */
				if(PINDX(jbp) != MAX_INDEX){
					postype rvel = Vec2DDotP(&er, &uij);
					if(rvel < 0){
						postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
						postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
						postype drampScale = dramp/scaleFactor;
						postype mu = rvel /* Voronoi-AV-reform: no length regularization */;
						postype meanden = 0.5*(ibp_den + jbp_den);
						postype meanCsound = 0.5*(ibp_csound + jbp_csound);
						pi = pi + (-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
					}
				}

				/* internal energy */
				Voro2D_point uradix_ui = get2dUpqradRk4(ibp, jbp, dtold);
				die += -pi * Vec2DDotP(&uradix_ui, &dS);

				/* total energy */
				Voro2D_point ua;
				ua.x = ibp_vx + uradix_ui.x;
				ua.y = ibp_vy + uradix_ui.y;
				dte += -pi*Vec2DDotP(&ua, &dS);

				/* kinetic energy */
				Voro2D_point ub;
				ub.x = ibp_vx;
				ub.y = ibp_vy;
				dke += -pi*Vec2DDotP(&ub, &dS);

				/* force */
				fx += -pi * dS.x;
				fy += -pi * dS.y;

				/* CFL timestep */
				Voro2D_point dv;
				dv.x = (jbp_vx - ibp_vx);
				dv.y = (jbp_vy - ibp_vy);
				postype VdotR = Vec2DDotP(&dv, &er);
				postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
				/* CFL floor: prevent dt→0 from particle clustering
				   or geometric proximity of ghost mirrors */
				postype dramp_cfl = fmax(dramp, 0.25*sqrt(ibp->volume));
				postype dt = 2*Courant*dramp_cfl/vsig;
				if(isnan(dt)){
					DEBUGPRINT("P%d kNN has error dt %d %ld : %g %g : %g %g : %g\n",
							MYID(simpar), ii, PINDX(ibp), dramp, vsig, jbp_csound, ibp_csound, VdotR);
					dt = 1e-10;  /* Fallback: dont crash, use tiny dt */
				}
				ibp->dt = MIN(ibp->dt, dt);
				if(dt < Dtime) Dtime = dt;
			}
			tmp = tmp->upperlink;
		} while(tmp != vorocorner);
		ibp->die = die;
		ibp->ax = fx/ibp->mass;
		ibp->ay = fy/ibp->mass;

		if(isnan(fx)){
			static int nan_warn = 0;
			if(nan_warn < 10)
				DEBUGPRINT("P%d kNN has nan %d : %g %g %ld\n",
						MYID(simpar), ii, ibp->x, ibp->y, PINDX(ibp));
			nan_warn++;
			ibp->ax = 0; ibp->ay = 0; ibp->die = 0;
		}
		} /* if(knn_ok) */
		free(vorocorner);
		} /* if(!skip_this) */
	}
	{
		postype TDtime;
		MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
		Dtime = TDtime;
	}
	return Dtime;
}

void exam2d_centroidShift(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	/* Stride-aware: p_size may differ from sizeof(treevorork4particletype)
	   when stress particles are active (av_mode >= 1). */
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	int nbp = VORO_NP(simpar);
	postype boxsize = VORO_CELLBOX(simpar);

	postype xmin,ymin,zmin,xmax,ymax,zmax, cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
	// xmin,ymin,xmax,ymax are the boundaries of the local domain.
	xmin = Xmin_HydroExam(simpar)-cellsize;
    ymin = Ymin_HydroExam(simpar)-cellsize;
    xmax = Xmax_HydroExam(simpar)+cellsize;
    ymax = Ymax_HydroExam(simpar)+cellsize;

	int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);
#ifdef USE_CUDA
	if(GAS_Kappa(simpar)>0) det2d_dpqRK4_GPU(simpar, xmin, ymin, xmax, ymax,
		paddingAllTreeParticles, mkLinkedList2D);
#else
	if(GAS_Kappa(simpar)>0 || gfs_dynw2_active) det2d_dpqRK4(simpar,paddingAllTreeParticles);
#endif
	// prepare the linked-list cells for mkLinkedList2D()  and tree findings below
    CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
    // building the linked list with "paddingAllTreeParticles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
    mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);
	// repositioning the memory space
	char *bp_raw = (char*)VORORK4_TBP(simpar);
    int i;
    int iy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(iy=1;iy<my-1;iy++){
		int mp=2000;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
		for(ix=1;ix<mx-1;ix++){
			int np;
			treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
			int nneigh;
			Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
			Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			int i;
			for(i=0;i<np;i++){
				Voro2D_point center;
				center.x = p[i].x;
				center.y = p[i].y;
				center.indx = PINDX(p+i);
				center.w2 = p[i].w2;
				treevorork4particletype *ibp = p[i].bp;
				int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
                // in this case ax, ay are used for the centroid vector
                Voro2D_Corner centroid = Centroid2DPolygon(vorocorner);
				ibp->ax = centroid.x;
				ibp->ay = centroid.y;
			}
			free(p);free(neighbors); free(neighwork);
		}
		free(vorocorner);
	}
	my_free(VORORK4_TBPP(simpar));
	my_free(VORO_BASICCELL(simpar));
    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;
	float fshift = GAS_FCENTROID(simpar);
	bp_raw = (char*)VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		bpi->x += fshift*bpi->ax; // move it to the centroid  with a factor of shift
		bpi->y += fshift*bpi->ay;
		bpi->x = fmod(bpi->x+Lx,Lx);
		bpi->y = fmod(bpi->y+Ly,Ly);
	}
    // migrate particles between mpi ranks
	migrateTreeVorork4Particles(simpar);
}
void exam2dUpdateVol(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype boxsize = VORO_CELLBOX(simpar);
	fprintf(stderr,"[UPDV] P%d enter np=%d boxsize=%g\n", MYID(simpar), nbp, boxsize); fflush(stderr);

	postype xmin,ymin,zmin,xmax,ymax,zmax, cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
    ymin = Ymin_HydroExam(simpar)-cellsize;
    xmax = Xmax_HydroExam(simpar)+cellsize;
    ymax = Ymax_HydroExam(simpar)+cellsize;
	fprintf(stderr,"[UPDV] P%d cellsize=%g xrng=(%g,%g) yrng=(%g,%g) Kappa=%g\n",
		MYID(simpar), cellsize, xmin, xmax, ymin, ymax, GAS_Kappa(simpar)); fflush(stderr);
//  determine the minimum dpq for each particle by updating the w2ceil & w2
#ifdef USE_CUDA
	if(GAS_Kappa(simpar)>0){
		fprintf(stderr,"[UPDV] P%d before det2d_dpqRK4_GPU\n", MYID(simpar)); fflush(stderr);
		det2d_dpqRK4_GPU(simpar, xmin, ymin, xmax, ymax,
			paddingAllTreeParticles, mkLinkedList2D);
		fprintf(stderr,"[UPDV] P%d after det2d_dpqRK4_GPU\n", MYID(simpar)); fflush(stderr);
	}
#else
	if(GAS_Kappa(simpar)>0 || gfs_dynw2_active) det2d_dpqRK4(simpar,paddingAllTreeParticles);
#endif

	int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);
	fprintf(stderr,"[UPDV] P%d mx=%d my=%d\n", MYID(simpar), mx, my); fflush(stderr);
	// prepare the linked-list cells for mkLinkedList2D()  and tree findings below
	CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
	fprintf(stderr,"[UPDV] P%d basiccell allocated, calling mkLinkedList2D\n", MYID(simpar)); fflush(stderr);
	// building the linked list with "paddingAllTreeParticles()" which pads the domain
	// with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
	mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);
	fprintf(stderr,"[UPDV] P%d after mkLinkedList2D npad=%d\n", MYID(simpar), VORO_NPAD(simpar)); fflush(stderr);

	int i;
	int iy;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(iy=0;iy<my;iy++){
        int mp=1000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=0;ix<mx;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.csound = p[i].csound;
                center.w2 = p[i].w2;
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
				/* same cell as the stage passes: widened search, then the disc */
				Voro2D_point *ne_u = NULL, *nwe_u = NULL, *nw_u = neighwork;
				if(Voro2D_HasUnassignedCorner(vorocorner, ip) || voro_sec_short(simpar, vorocorner, ip, (double)center.x, (double)center.y, ix, iy, 1)){
					int R, nne;
					for(R=2; R<=voro_expand_rmax(simpar); R++){
						if(ne_u){ free(ne_u); free(nwe_u); }
						ne_u = searchCellRk4Neighbors2D_R(simpar, ix, iy, R, &nne);
						nwe_u = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nne);
						ip = Voro2D_FindVC(&center, ne_u, nwe_u, nne, vorocorner, mp, boxsize);
						nw_u = nwe_u;
						if(!Voro2D_HasUnassignedCorner(vorocorner, ip) && !voro_sec_short(simpar, vorocorner, ip, (double)center.x, (double)center.y, ix, iy, R)) break;
					}
				}
				{ double vrR = voro_disc_radius(simpar); if(vrR > 0) voro2d_restrict_disc(vorocorner, mp, vrR); }
				treevorork4particletype *ibp = p[i].bp;
//				ibp->volume = Area2DPolygon(vorocorner,mp);
                get2dAreaAvgNeighorPressure(ibp,vorocorner, nw_u,bp);
				ibp->den = ibp->mass / ibp->volume;
				if(ne_u){ free(ne_u); free(nwe_u); }
			}
			free(neighwork);
			free(neighbors); 
			free(p);
		}
		free(vorocorner);
	
	}
	my_free(VORO_BASICCELL(simpar));

	VORO_NPAD(simpar) = 0;
	my_free(VORORK4_TBPP(simpar));
}
double exam2d_vph_rk4(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	postype xmin,ymin,xmax,ymax,OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);
//	DEBUGPRINT("P%d has x/y minmax= %g %g : %g %g ::: %g\n", MYID(simpar), xmin,ymin,xmax,ymax, cellsize);


    int i;
    postype Dtime,dt;
    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) bp[i].rk4.w2backup = bp[i].w2;
	// Runge-Kutta 4-th order time evolution of r and vr
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, 1);
    Dtime = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){ 
		bp[i].rk4.k1x = bp[i].vx*Dtime;
        bp[i].rk4.k1y = bp[i].vy*Dtime;
        bp[i].rk4.k1vx = bp[i].ax*Dtime;
        bp[i].rk4.k1vy = bp[i].ay*Dtime;
        bp[i].rk4.k1ie = bp[i].die*Dtime;
    }
    // RK4 second
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += bp[i].rk4.k1x*0.5;
        bp[i].y += bp[i].rk4.k1y*0.5;
        bp[i].vx += bp[i].rk4.k1vx*0.5;
        bp[i].vy += bp[i].rk4.k1vy*0.5;
        bp[i].ie += bp[i].rk4.k1ie*0.5;
        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
	migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
    dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
        paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k2x = bp[i].vx*Dtime;
        bp[i].rk4.k2y = bp[i].vy*Dtime;
        bp[i].rk4.k2vx = bp[i].ax*Dtime;
        bp[i].rk4.k2vy = bp[i].ay*Dtime;
        bp[i].rk4.k2ie = bp[i].die*Dtime;
    }
    // RK4 third
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += (bp[i].rk4.k2x-bp[i].rk4.k1x)*0.5;
        bp[i].y += (bp[i].rk4.k2y-bp[i].rk4.k1y)*0.5;
        bp[i].vx += (bp[i].rk4.k2vx-bp[i].rk4.k1vx)*0.5;
        bp[i].vy += (bp[i].rk4.k2vy-bp[i].rk4.k1vy)*0.5;
        bp[i].ie += (bp[i].rk4.k2ie-bp[i].rk4.k1ie)*0.5;

        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
	dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k3x = bp[i].vx*Dtime;
        bp[i].rk4.k3y = bp[i].vy*Dtime;
        bp[i].rk4.k3vx = bp[i].ax*Dtime;
        bp[i].rk4.k3vy = bp[i].ay*Dtime;
        bp[i].rk4.k3ie = bp[i].die*Dtime;
    }
    // RK4 fourth
    for(i=0;i<VORO_NP(simpar);i++){

        bp[i].x += bp[i].rk4.k3x - bp[i].rk4.k2x*0.5;
        bp[i].y += bp[i].rk4.k3y - bp[i].rk4.k2y*0.5;
        bp[i].vx += bp[i].rk4.k3vx - bp[i].rk4.k2vx*0.5;
        bp[i].vy += bp[i].rk4.k3vy - bp[i].rk4.k2vy*0.5;
        bp[i].ie += bp[i].rk4.k3ie - bp[i].rk4.k2ie*0.5;

        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
    dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k4x = bp[i].vx*Dtime;
        bp[i].rk4.k4y = bp[i].vy*Dtime;
        bp[i].rk4.k4vx = bp[i].ax*Dtime;
        bp[i].rk4.k4vy = bp[i].ay*Dtime;
        bp[i].rk4.k4ie = bp[i].die*Dtime;
    }
    // migrate particles between mpi ranks
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x -= bp[i].rk4.k3x;
        bp[i].y -= bp[i].rk4.k3y;
        bp[i].vx -= bp[i].rk4.k3vx;
        bp[i].vy -= bp[i].rk4.k3vy;
        bp[i].ie -= bp[i].rk4.k3ie;
        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
	// migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

    // to finalize ke & ie updates 
	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
    // Update position and velocity using RK4
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].x  += (bp[i].rk4.k1x +2*bp[i].rk4.k2x +2*bp[i].rk4.k3x +bp[i].rk4.k4x )/6.;
        bp[i].y  += (bp[i].rk4.k1y +2*bp[i].rk4.k2y +2*bp[i].rk4.k3y +bp[i].rk4.k4y )/6.;
        bp[i].vx += (bp[i].rk4.k1vx+2*bp[i].rk4.k2vx+2*bp[i].rk4.k3vx+bp[i].rk4.k4vx)/6.;
        bp[i].vy += (bp[i].rk4.k1vy+2*bp[i].rk4.k2vy+2*bp[i].rk4.k3vy+bp[i].rk4.k4vy)/6.;
        bp[i].ie += (bp[i].rk4.k1ie+2*bp[i].rk4.k2ie+2*bp[i].rk4.k3ie+bp[i].rk4.k4ie)/6.;

        bp[i].x = fmod(bp[i].x+Lx,Lx);
        bp[i].y = fmod(bp[i].y+Ly,Ly);
    }

    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);

	//centroid shift
	if(GAS_FCENTROID(simpar)>0) 
		exam2d_centroidShift(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	// update volume & density with one-step evolved positions 
    exam2dUpdateVol(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	bp = VORORK4_TBP(simpar);
	postype dmean = SIMBOX(simpar).x.max/NX(simpar);

	// split getw2forHydroParticle and pressure update because the former uses the latter.
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(GAS_Kappa(simpar) <0) { 
			bp[i].w2 = -GAS_Kappa(simpar); 
		} 
		else if (GAS_Kappa(simpar) >0){ 
			bp[i].w2hydro = getw2forHydroParticle(simpar,(bp+i),Dtime);
			applyW2Controls(simpar, bp+i, Dtime);
		}
    }
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].den = bp[i].mass/bp[i].volume;
        bp[i].pressure = bp[i].ie/bp[i].volume * (Gamma-1);
        bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
    }

    return Dtime;
}

void dumpRk4Data2D(SimParameters *simpar, int nstep, postype t, postype dt){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int myid, nid;
	myid = MYID(simpar);
	nid = NID(simpar);
	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);


	char outfile[190];
	sprintf(outfile,"glout.%.6d.dat",nstep);
	int i;
	int snp=0;
	FILE *wp;
	int tnp = VORO_TNP(simpar);
	for(i=0;i<nid;i++){
		if(myid ==i){
			if(myid==0){
				wp = fopen(outfile,"w"); 
				fwrite(&tnp, sizeof(int),1,wp);
			}
			else wp = fopen(outfile,"a");
			fwrite(bp, sizeof(treevorork4particletype), np, wp);
			fclose(wp);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if(myid==0){
		wp = fopen(outfile,"a");
		fwrite(&t, sizeof(postype), 1, wp);
		fwrite(&nx, sizeof(int), 1, wp);
		fwrite(&ny, sizeof(int), 1, wp);
		fwrite(&(SIMBOX(simpar).x.max), sizeof(postype), 1, wp);
		fwrite(&(SIMBOX(simpar).y.max), sizeof(postype), 1, wp);
		fwrite(&GAS_AlphaVis(simpar), sizeof(float), 1, wp);
		fwrite(&GAS_BetaVis(simpar), sizeof(float), 1, wp);
		fwrite(&GAS_ETAVIS(simpar), sizeof(float), 1, wp);
		fwrite(&GAS_EPSVIS(simpar), sizeof(float), 1, wp);
		fclose(wp);
	}
	/*
	MPI_Barrier(MPI_COMM_WORLD);
	DEBUGPRINT("P%d is now exiting the dump\n", myid);
	*/
}


void readRk4Data2D(SimParameters *simpar, postype *t, postype *dt, int nstep){
	float dtold;
	char infile[190];
	sprintf(infile,"glout.%.6d.dat", nstep);
	int myid = MYID(simpar);
	if(myid ==0){
		printf("P0: is now reading starting file: %s\n", infile);
		FILE *fp = fopen(infile,"r");
		int np;
		fread(&np, sizeof(int), 1, fp);
		VORO_TNP(simpar) = VORO_NP(simpar) = np;
		VORORK4_TBP(simpar) = (treevorork4particletype*)my_malloc(sizeof(treevorork4particletype)*np);
		fread(VORORK4_TBP(simpar),sizeof(treevorork4particletype), np,fp);
		fread(t,sizeof(postype), 1,fp);
		fread(dt,sizeof(postype), 1,fp);
		fread(&GAS_dtold(simpar),sizeof(float), 1,fp);
		fclose(fp);
	}
	else {
		VORO_NP(simpar) = 0;
		VORORK4_TBP(simpar) = (treevorork4particletype*)my_malloc(sizeof(treevorork4particletype)*100);;

	}
	DEBUGPRINT("P%d: after reading iput data\n", myid);
	migrateTreeVorork4Particles(simpar);
	MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(dt, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(&GAS_dtold(simpar), 1, MPI_FLOAT,0, MPI_COMM(simpar));
	MPI_Barrier(MPI_COMM(simpar));
	BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
}

double exam2d_vph(SimParameters *simpar, 
		double (*measureW2Exam2D)(SimParameters *, double, double, double),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		)
{
    postype boxsize = VORO_CELLBOX(simpar);
    treevorork4particletype *bp = VORORK4_TBP(simpar);
    int nbp = VORO_NP(simpar);
    int isave = -1;
    postype Dtime = 1.e10;
    postype Courant = GAS_COURANT(simpar);
    postype cellsize;


    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;
    postype xmin,ymin,zmin,xmax,ymax,zmax, pwidth;
    pwidth = BASICCELL_CELLWIDTH(simpar) = cellsize = HydroGridSize(simpar);
    xmin = Xmin_HydroExam(simpar)-pwidth;
    ymin = Ymin_HydroExam(simpar)-pwidth;
    xmax = Xmax_HydroExam(simpar)+pwidth;
    ymax = Ymax_HydroExam(simpar)+pwidth;

    postype dtold = GAS_dtold(simpar);

//  DEBUGPRINT("P%d has dtold = %g\n", MYID(simpar), dtold);

    int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);

    // prepare the linked-list cells for mkLinkedList2D()  and tree findings below
    CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
    // building the linked list with "paddingTreeVorork4Particles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingTreeVorork4Particles);

    postype Gamma = GAS_GAMMA(simpar);

    float alphavis = GAS_AlphaVis(simpar);
    float betavis = GAS_BetaVis(simpar);
    float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
    float epsvis = GAS_EPSVIS(simpar);


//  DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);


    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = findCellRk4BP2D(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = searchCellRk4Neighbors2D(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.csound = p[i].csound;
                center.w2 = p[i].w2;

                treevorork4particletype *ibp = p[i].bp;
                ibp->dt = 1.e10; // initialization of dt
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
//                ibp->volume = Area2DPolygon(vorocorner, mp);
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);
                /* These are to use register */
                postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
                ibp_vx = ibp->vx;
                ibp_vy = ibp->vy;
                ibp_csound = ibp->csound;
                ibp_pressure = ibp->pressure;
                ibp_den = ibp->den;


                Voro2D_Corner *tmp,*tmp2;
                tmp = vorocorner;
                int j;
                double die, dte,dke,fx,fy;
                die = dte = dke = fx = fy = 0;
                do {
                    if(tmp->upperrelated >=0){
                        /* These are to use register */
                        treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
                        postype jbp_vx = jbp->vx;
                        postype jbp_vy = jbp->vy;
                        postype jbp_csound = jbp->csound;
//                        postype jbp_pressure = jbp->pressure;
						postype jbp_pressure = getNeighborPressure(jbp,ibp);

                        postype jbp_den = jbp->den;

                        Voro2D_point line;
                        tmp2 = tmp->upperlink;
                        line.x = tmp2->x - tmp->x;
                        line.y = tmp2->y - tmp->y;
                        Voro2D_point dS = voro2D_norm(&line);
                        postype facearea = Vec2DLength(tmp,tmp2);
                        dS.x = facearea*dS.x;
                        dS.y = facearea*dS.y;

                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,RT);

                        Voro2D_point uij,ua,ub;
                        uij.x = (jbp_vx - ibp_vx); 
                        uij.y = (jbp_vy - ibp_vy);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &uij);
                        postype mu, meanden, meanCsound;
                        if(rvel<0){
							postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
							postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
							postype drampScale = dramp/scaleFactor;
                            mu = rvel /* Voronoi-AV-reform: no length regularization */;
                            meanden = 0.5*(ibp_den + jbp_den);
                            meanCsound = 0.5*(ibp_csound + jbp_csound);
                            pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        /* for the internal energy */
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

                        /* for the total energy */
                        ua.x = ibp_vx + uradix_ui.x;
                        ua.y = ibp_vy + uradix_ui.y;
                        dte += -pi*Vec2DDotP(&ua, &dS);

                        /* for the kinetic energy */
                        ub.x = ibp_vx;
                        ub.y = ibp_vy;
                        dke += -pi*Vec2DDotP(&ub, &dS);

                        /* for the force */
                        fx += -pi * dS.x;
                        fy += -pi * dS.y;

                        Voro2D_point dv;
                        dv.x = (jbp_vx - ibp_vx);
                        dv.y = (jbp_vy - ibp_vy);

                        postype VdotR = Vec2DDotP(&dv,&er);
                        postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
                        postype dramp_cfl = fmax(dramp, 0.25*sqrt(ibp->volume));
                        postype dt = 2*Courant*dramp_cfl/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d %ld : %g %g : %g %g : %g : dv.xy= %g %g\n",
                                    MYID(simpar), i, PINDX(p+i), dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y);
                            dt = 1e-10;
                        }
                        ibp->dt = MIN(ibp->dt,dt);
                        if(dt < Dtime){
                            Dtime = dt;
                            // isave = id;
                        }
                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
				/*
                ibp->dte = dte;
                ibp->dke = dke;
				*/
                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;
                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin);
                }

            }
            free(p);free(neighbors); free(neighwork);

        }
        free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));

//  DEBUGPRINT("P%d is exiting here\n", MYID(simpar)); exit(9);

    postype TDtime;
    MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
    Dtime = TDtime;

    int i;
    postype dmax=0;
//  postype eps2 = KH_EPS(simpar)*KH_EPS(simpar);
#ifdef _OPENMP
#pragma omp parallel for reduction(max: dmax)
#endif
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].vx += bp[i].ax *Dtime;
        bp[i].vy += bp[i].ay *Dtime;
        bp[i].x += bp[i].vx *Dtime;
        bp[i].y += bp[i].vy *Dtime;
        dmax = MAX(dmax, fabs(bp[i].ay));
        bp[i].x = fmod(bp[i].x+Lx,Lx);
        bp[i].y = fmod(bp[i].y+Ly,Ly);
        bp[i].w2old = bp[i].w2;
    }
    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);
    if(0){
        postype xmin,ymin,xmax,ymax;
        xmin = ymin = 1.e20;
        xmax = ymax =-1.e20;
        for(i=0;i<VORO_NP(simpar);i++){
            xmin = MIN(xmin, bp[i].x);
            ymin = MIN(ymin, bp[i].y);
            xmax = MAX(xmax, bp[i].x);
            ymax = MAX(ymax, bp[i].y);
        }
    }

// This is for centroid shift.
	if(GAS_FCENTROID(simpar) > 0 ){
	    if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingTreeVorork4Particles);
    // Update the linked list 
		mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingTreeVorork4Particles);

    // repositioning the memory space 
	    bp = VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	    for(iy=1;iy<my-1;iy++){ 
			int mp=2000; 
			Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp); 
			postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a; 
			int ix; 
			for(ix=1;ix<mx-1;ix++){ 
				int np; 
				treevorork4particletype *p = findCellRk4BP2D(simpar,ix,iy,&np); 
				int nneigh; 
				Voro2D_point *neighbors = searchCellRk4Neighbors2D(simpar,ix,iy,&nneigh); 
				Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh); 
				int i; 
				for(i=0;i<np;i++){ 
					Voro2D_point center; 
					center.x = p[i].x; 
					center.y = p[i].y; 
					center.indx = PINDX(p+i); 
					center.w2 = p[i].w2; 
					
					treevorork4particletype *ibp = p[i].bp; 
					int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
					// in this case ax, ay are used for the centroid vector
					Voro2D_Corner centroid = Centroid2DPolygon(vorocorner);
					ibp->ax = centroid.x;
					ibp->ay = centroid.y; 
				} 
				free(p);free(neighbors); free(neighwork); 
			} 
			free(vorocorner); 
		} 
		my_free(VORORK4_TBPP(simpar));
		float fshift = GAS_FCENTROID(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	    for(i=0;i<VORO_NP(simpar);i++){ 
			bp[i].x += fshift*bp[i].ax; // move it to the centroid  with a factor of shift
			bp[i].y += fshift*bp[i].ay; 
        	bp[i].x = fmod(bp[i].x+Lx,Lx);
	        bp[i].y = fmod(bp[i].y+Ly,Ly);
		} 
		// migrate particles between mpi ranks 
		migrateTreeVorork4Particles(simpar);
	}
//
//
	if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingTreeVorork4Particles);
    // Update the linked list
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingTreeVorork4Particles);


    // repositioning the memory space
    bp = VORORK4_TBP(simpar);



#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = findCellRk4BP2D(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = searchCellRk4Neighbors2D(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.w2 = p[i].w2;

                treevorork4particletype *ibp = p[i].bp;
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
//              ibp->w2 = ibp->w2ceil = center.w2ceil;
                // These are to use register 
                postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
                ibp_vx = ibp->vx;
                ibp_vy = ibp->vy;
                ibp_csound = ibp->csound;
                ibp_pressure = ibp->pressure;
                ibp_den = ibp->den;

//                ibp->volume = Area2DPolygon(vorocorner, mp);
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);




                Voro2D_Corner *tmp,*tmp2;
                tmp = vorocorner;
                int j;
//                double die, dte,dke,fx,fy;
//               die = dte = dke = fx = fy = 0;
				double die, cx,cy;
				die = cx = cy = 0;
                do {
                    if(tmp->upperrelated >=0){
                        // These are to use register 
                        treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
                        postype jbp_vx = jbp->vx;
                        postype jbp_vy = jbp->vy;
                        postype jbp_csound = jbp->csound;
//                        postype jbp_pressure = jbp->pressure;
						postype jbp_pressure = getNeighborPressure(jbp,ibp);
                        postype jbp_den = jbp->den;

                        Voro2D_point line;
                        tmp2 = tmp->upperlink;
                        line.x = tmp2->x - tmp->x;
                        line.y = tmp2->y - tmp->y;
                        Voro2D_point dS = voro2D_norm(&line);
                        postype facearea = Vec2DLength(tmp,tmp2);
                        dS.x = facearea*dS.x;
                        dS.y = facearea*dS.y;
//                      postype pi = Half*(jbp_pressure + ibp_pressure);
                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,RT);

                        Voro2D_point uij,ua,ub;
                        uij.x = (jbp_vx - ibp_vx);
                        uij.y = (jbp_vy - ibp_vy);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &uij);
                        if(rvel<0){
							 postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
							 postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
							 postype drampScale = dramp/scaleFactor;
                             postype mu = rvel /* Voronoi-AV-reform: no length regularization */;
                             postype meanden = 0.5*(ibp_den + jbp_den);
                             postype meanCsound = 0.5*(ibp_csound + jbp_csound);
                             pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
//                ibp->ax = fx/ibp->mass;
//                ibp->ay = fy/ibp->mass;
            }
            free(p);free(neighbors); free(neighwork);
        }
        free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));
    my_free(VORO_BASICCELL(simpar));

    // update the w2ceil
    postype dmean = SIMBOX(simpar).x.max/NX(simpar);

//  DEBUGPRINT("P%d is now updating the internal energy\n", MYID(simpar));
	// split getw2forHydroParticle and pressure update because the former uses the latter.
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
        if(GAS_Kappa(simpar) <0) {
            bp[i].w2 = -GAS_Kappa(simpar);
        }
        else if (GAS_Kappa(simpar) >0){
          	bp[i].w2hydro = getw2forHydroParticle(simpar,(bp+i),Dtime);
            applyW2Controls(simpar, bp+i, Dtime);
        }
    }
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].ie += bp[i].die * Dtime;
        bp[i].den  = bp[i].mass/bp[i].volume;
        bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1);
        bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return Dtime;
}
//
//
//
double exam2d_vph_int_rt(SimParameters *simpar, 
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		) {
    postype boxsize = VORO_CELLBOX(simpar);
    treevorork4particletype *bp = VORORK4_TBP(simpar);
    int nbp = VORO_NP(simpar);
    int isave = -1;
    postype Dtime = 1.e10;
    postype Courant = GAS_COURANT(simpar);
    postype cellsize;


    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;
    postype xmin,ymin,zmin,xmax,ymax,zmax, pwidth;
    pwidth = BASICCELL_CELLWIDTH(simpar) = cellsize = HydroGridSize(simpar);
    xmin = Xmin_HydroExam(simpar)-pwidth;
    ymin = Ymin_HydroExam(simpar)-pwidth;
    xmax = Xmax_HydroExam(simpar)+pwidth;
    ymax = Ymax_HydroExam(simpar)+pwidth;

    postype dtold = GAS_dtold(simpar);

//  DEBUGPRINT("P%d has dtold = %g\n", MYID(simpar), dtold);

    int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);

    // prepare the linked-list cells for mkLinkedList2D()  and tree findings below
    CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
    // building the linked list with "paddingAllTreeParticles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

    postype Gamma = GAS_GAMMA(simpar);

    float alphavis = GAS_AlphaVis(simpar);
    float betavis = GAS_BetaVis(simpar);
    float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
    float epsvis = GAS_EPSVIS(simpar);


//  DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);


    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.csound = p[i].csound;
                center.w2 = p[i].w2;

                treevorork4particletype *ibp = p[i].bp;
                ibp->dt = 1.e10; // initialization of dt
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
//                ibp->volume = Area2DPolygon(vorocorner, mp);
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);
                // These are to use register 
                postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
                ibp_vx = ibp->vx;
                ibp_vy = ibp->vy;
                ibp_csound = ibp->csound;
                ibp_pressure = ibp->pressure;
                ibp_den = ibp->den;


                Voro2D_Corner *tmp,*tmp2;
                tmp = vorocorner;
                int j;
                double die, dte,dke,fx,fy;
                die = dte = dke = fx = fy = 0;
                do {
                    if(tmp->upperrelated >=0){
                        // These are to use register 
                        treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
                        postype jbp_vx = jbp->vx;
                        postype jbp_vy = jbp->vy;
                        postype jbp_csound = jbp->csound;
                        postype jbp_pressure = jbp->pressure; // for RT, the wall particle doesnot mirror the pressure
//						postype jbp_pressure = getNeighborPressure(jbp,ibp);
                        postype jbp_den = jbp->den;

                        Voro2D_point line;
                        tmp2 = tmp->upperlink;
                        line.x = tmp2->x - tmp->x;
                        line.y = tmp2->y - tmp->y;
                        Voro2D_point dS = voro2D_norm(&line);
                        postype facearea = Vec2DLength(tmp,tmp2);
                        dS.x = facearea*dS.x;
                        dS.y = facearea*dS.y;

                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,RT);

                        Voro2D_point uij,ua,ub;
                        uij.x = (jbp_vx - ibp_vx); 
                        uij.y = (jbp_vy - ibp_vy);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &uij);
                        postype mu, meanden, meanCsound;
                        if(rvel<0){ 
							postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
							postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
							postype drampScale = dramp/scaleFactor;
                            mu = rvel /* Voronoi-AV-reform: no length regularization */;
                            meanden = 0.5*(ibp_den + jbp_den);
                            meanCsound = 0.5*(ibp_csound + jbp_csound);
                            pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

                        // for the total energy
                        ua.x = ibp_vx + uradix_ui.x;
                        ua.y = ibp_vy + uradix_ui.y;
                        dte += -pi*Vec2DDotP(&ua, &dS);

                        // for the kinetic energy 
                        ub.x = ibp_vx;
                        ub.y = ibp_vy;
                        dke += -pi*Vec2DDotP(&ub, &dS);

                        // for the force 
                        fx += -pi * dS.x;
                        fy += -pi * dS.y;

                        Voro2D_point dv;
                        dv.x = (jbp_vx - ibp_vx);
                        dv.y = (jbp_vy - ibp_vy);

                        postype VdotR = Vec2DDotP(&dv,&er);
                        postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
                        postype dramp_cfl = fmax(dramp, 0.25*sqrt(ibp->volume));
                        postype dt = 2*Courant*dramp_cfl/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d %ld : %g %g : %g %g : %g : dv.xy= %g %g\n",
                                    MYID(simpar), i, PINDX(p+i), dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y);
                            dt = 1e-10;
                        }
                        ibp->dt = MIN(ibp->dt,dt);
                        if(dt < Dtime){
                            Dtime = dt;
                            // isave = id;
                        }
                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
				/*
                ibp->dte = dte;
                ibp->dke = dke;
				*/
                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;
                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin);
                }

            }
            free(p);free(neighbors); free(neighwork);

        }
        free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));

//  DEBUGPRINT("P%d is exiting here\n", MYID(simpar)); exit(9);

    postype TDtime;
    MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN,  MPI_COMM(simpar));
    Dtime = TDtime;

    int i;
    postype dmax=0;
//  postype eps2 = KH_EPS(simpar)*KH_EPS(simpar);
#ifdef _OPENMP
#pragma omp parallel for reduction(max: dmax)
#endif
    for(i=0;i<VORO_NP(simpar);i++){
//		if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
		if(targetBP(bp+i, Lx, Ly))
		{
			bp[i].ax += GAS_ACCX(simpar);
			bp[i].ay += GAS_ACCY(simpar);
	        bp[i].vx += bp[i].ax *Dtime; 
			bp[i].vy += bp[i].ay *Dtime; 
			bp[i].x += bp[i].vx *Dtime; 
			bp[i].y += bp[i].vy *Dtime; 
			dmax = MAX(dmax, fabs(bp[i].ay)); 
			bp[i].x = fmod(bp[i].x+Lx,Lx); 
			bp[i].y = fmod(bp[i].y+Ly,Ly); 
			bp[i].w2old = bp[i].w2;
		}
		else {
			bp[i].vx = bp[i].vy = bp[i].ax = bp[i].ay = 0;
		}
    }
    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);

// This is for centroid shift.
	if(GAS_FCENTROID(simpar) > 0 ){
	    if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingAllTreeParticles);
    // Update the linked list 
		mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

    // repositioning the memory space 
	    bp = VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	    for(iy=1;iy<my-1;iy++){ 
			int mp=2000; 
			Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp); 
			postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a; 
			int ix; 
			for(ix=1;ix<mx-1;ix++){ 
				int np; 
				treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np); 
				int nneigh; 
				Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh); 
				Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh); 
				int i; 
				for(i=0;i<np;i++){ 
					Voro2D_point center; 
					center.x = p[i].x; 
					center.y = p[i].y; 
					center.indx = PINDX(p+i); 
					center.w2 = p[i].w2; 
					
					treevorork4particletype *ibp = p[i].bp; 
					int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
					// in this case ax, ay are used for the centroid vector
					Voro2D_Corner centroid = Centroid2DPolygon(vorocorner);
					ibp->ax = centroid.x;
					ibp->ay = centroid.y; 
				} 
				free(p);free(neighbors); free(neighwork); 
			} 
			free(vorocorner); 
		} 
		my_free(VORORK4_TBPP(simpar));
		float fshift = GAS_FCENTROID(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	    for(i=0;i<VORO_NP(simpar);i++){ 
//			if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
			if(targetBP(bp+i, Lx, Ly))
			{
				bp[i].x += fshift*bp[i].ax; // move it to the centroid  with a factor of shift
				bp[i].y += fshift*bp[i].ay; 
   		     	bp[i].x = fmod(bp[i].x+Lx,Lx);
		        bp[i].y = fmod(bp[i].y+Ly,Ly);
			}
		} 
		// migrate particles between mpi ranks 
		migrateTreeVorork4Particles(simpar);
	}
//
//
	if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingAllTreeParticles);
    // Update the linked list
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);


    // repositioning the memory space
    bp = VORORK4_TBP(simpar);



#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.w2 = p[i].w2;

                treevorork4particletype *ibp = p[i].bp;
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
//              ibp->w2 = ibp->w2ceil = center.w2ceil;
                // These are to use register 
                postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
                ibp_vx = ibp->vx;
                ibp_vy = ibp->vy;
                ibp_csound = ibp->csound;
                ibp_pressure = ibp->pressure;
                ibp_den = ibp->den;

//                ibp->volume = Area2DPolygon(vorocorner, mp);
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);




                Voro2D_Corner *tmp,*tmp2;
                tmp = vorocorner;
                int j;
//                double die, dte,dke,fx,fy;
//               die = dte = dke = fx = fy = 0;
				double die, cx,cy;
				die = cx = cy = 0;
                do {
                    if(tmp->upperrelated >=0){
                        // These are to use register 
                        treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
                        postype jbp_vx = jbp->vx;
                        postype jbp_vy = jbp->vy;
                        postype jbp_csound = jbp->csound;
                        postype jbp_pressure = jbp->pressure; // for RT, the wall particle doesnot mirror the pressure
//						postype jbp_pressure = getNeighborPressure(jbp,ibp);
                        postype jbp_den = jbp->den;

                        Voro2D_point line;
                        tmp2 = tmp->upperlink;
                        line.x = tmp2->x - tmp->x;
                        line.y = tmp2->y - tmp->y;
                        Voro2D_point dS = voro2D_norm(&line);
                        postype facearea = Vec2DLength(tmp,tmp2);
                        dS.x = facearea*dS.x;
                        dS.y = facearea*dS.y;
//                      postype pi = Half*(jbp_pressure + ibp_pressure);
                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,RT);

                        Voro2D_point uij,ua,ub;
                        uij.x = (jbp_vx - ibp_vx);
                        uij.y = (jbp_vy - ibp_vy);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &uij);
                        if(rvel<0){
							 postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
							 postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
							 postype drampScale = dramp/scaleFactor;
                             postype mu = rvel /* Voronoi-AV-reform: no length regularization */;
                             postype meanden = 0.5*(ibp_den + jbp_den);
                             postype meanCsound = 0.5*(ibp_csound + jbp_csound);
                             pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
//                ibp->ax = fx/ibp->mass;
//                ibp->ay = fy/ibp->mass;
            }
            free(p);free(neighbors); free(neighwork);
        }
        free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));
    my_free(VORO_BASICCELL(simpar));

    // update the w2ceil
    postype dmean = SIMBOX(simpar).x.max/NX(simpar);
//  DEBUGPRINT("P%d is now updating the internal energy\n", MYID(simpar));
	// split getw2forHydroParticle and pressure update because the former uses the latter.
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
			if(GAS_Kappa(simpar) <0) { 
				bp[i].w2 = -GAS_Kappa(simpar); 
			} 
			else if (GAS_Kappa(simpar) >0){ 
          		bp[i].w2hydro = getw2forHydroParticle(simpar,(bp+i),Dtime);
				applyW2Controls(simpar, bp+i, Dtime);
			} 
		}
	}
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
	        bp[i].ie += bp[i].die * Dtime; 
			bp[i].den  = bp[i].mass/bp[i].volume; 
			bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1); 
			bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den); 
			if(bp[i].volume ==0) 
				DEBUGPRINT("P%d has wrong volume at i= %d with den= %g at xy= %g %g\n", MYID(simpar), 
						i,bp[i].den, bp[i].x, bp[i].y); 
			if(isnan(bp[i].csound)){ 
				DEBUGPRINT("P%d has nan at the sound speed = %d with" 
						"id= %ld: pressure= %g den= %g Delta IE=%g dIE= %g Dt= %g\n", 
						MYID(simpar), i,PINDX(bp+i), bp[i].pressure, bp[i].den, bp[i].ie,bp[i].die,Dtime); 
				exit(9);
			}
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return Dtime;
}
double exam2d_vph_int(SimParameters *simpar, 
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		) {
    postype boxsize = VORO_CELLBOX(simpar);
    treevorork4particletype *bp = VORORK4_TBP(simpar);
    int nbp = VORO_NP(simpar);
    int isave = -1;
    postype Dtime = 1.e10;
    postype Courant = GAS_COURANT(simpar);
    postype cellsize;


    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;
    postype xmin,ymin,zmin,xmax,ymax,zmax, pwidth;
    pwidth = BASICCELL_CELLWIDTH(simpar) = cellsize = HydroGridSize(simpar);
    xmin = Xmin_HydroExam(simpar)-pwidth;
    ymin = Ymin_HydroExam(simpar)-pwidth;
    xmax = Xmax_HydroExam(simpar)+pwidth;
    ymax = Ymax_HydroExam(simpar)+pwidth;

    postype dtold = GAS_dtold(simpar);

//  DEBUGPRINT("P%d has dtold = %g\n", MYID(simpar), dtold);

    int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);

    // prepare the linked-list cells for mkLinkedList2D()  and tree findings below
    CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
    // building the linked list with "paddingAllTreeParticles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

    postype Gamma = GAS_GAMMA(simpar);

    float alphavis = GAS_AlphaVis(simpar);
    float betavis = GAS_BetaVis(simpar);
    float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
    float epsvis = GAS_EPSVIS(simpar);


//  DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);


    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.csound = p[i].csound;
                center.w2 = p[i].w2;

                treevorork4particletype *ibp = p[i].bp;
                ibp->dt = 1.e10; // initialization of dt
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
//                ibp->volume = Area2DPolygon(vorocorner, mp);
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);
                // These are to use register 
                postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
                ibp_vx = ibp->vx;
                ibp_vy = ibp->vy;
                ibp_csound = ibp->csound;
                ibp_pressure = ibp->pressure;
                ibp_den = ibp->den;


                Voro2D_Corner *tmp,*tmp2;
                tmp = vorocorner;
                int j;
                double die, dte,dke,fx,fy;
                die = dte = dke = fx = fy = 0;
                do {
                    if(tmp->upperrelated >=0){
                        // These are to use register 
                        treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
                        postype jbp_vx = jbp->vx;
                        postype jbp_vy = jbp->vy;
                        postype jbp_csound = jbp->csound;
//                        postype jbp_pressure = jbp->pressure;
						postype jbp_pressure = getNeighborPressure(jbp,ibp);
                        postype jbp_den = jbp->den;

                        Voro2D_point line;
                        tmp2 = tmp->upperlink;
                        line.x = tmp2->x - tmp->x;
                        line.y = tmp2->y - tmp->y;
                        Voro2D_point dS = voro2D_norm(&line);
                        postype facearea = Vec2DLength(tmp,tmp2);
                        dS.x = facearea*dS.x;
                        dS.y = facearea*dS.y;

                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,RT);

                        Voro2D_point uij,ua,ub;
                        //
                        uij.x = (jbp_vx - ibp_vx); 
                        uij.y = (jbp_vy - ibp_vy);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &uij);
                        postype mu, meanden, meanCsound;
                        if(rvel<0){ 
							postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
							postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
							postype drampScale = dramp/scaleFactor;
                            mu = rvel /* Voronoi-AV-reform: no length regularization */;
                            meanden = 0.5*(ibp_den + jbp_den);
                            meanCsound = 0.5*(ibp_csound + jbp_csound);
                            pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

                        // for the total energy
                        ua.x = ibp_vx + uradix_ui.x;
                        ua.y = ibp_vy + uradix_ui.y;
                        dte += -pi*Vec2DDotP(&ua, &dS);

                        // for the kinetic energy 
                        ub.x = ibp_vx;
                        ub.y = ibp_vy;
                        dke += -pi*Vec2DDotP(&ub, &dS);

                        // for the force 
                        fx += -pi * dS.x;
                        fy += -pi * dS.y;

                        Voro2D_point dv;
                        dv.x = (jbp_vx - ibp_vx);
                        dv.y = (jbp_vy - ibp_vy);

                        postype VdotR = Vec2DDotP(&dv,&er);
                        postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
                        postype dramp_cfl = fmax(dramp, 0.25*sqrt(ibp->volume));
                        postype dt = 2*Courant*dramp_cfl/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d %ld : %g %g : %g %g : %g : dv.xy= %g %g\n",
                                    MYID(simpar), i, PINDX(p+i), dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y);
                            dt = 1e-10;
                        }
                        ibp->dt = MIN(ibp->dt,dt);
                        if(dt < Dtime){
                            Dtime = dt;
                            // isave = id;
                        }
                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
				/*
                ibp->dte = dte;
                ibp->dke = dke;
				*/
                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;
                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin);
                }

            }
            free(p);free(neighbors); free(neighwork);

        }
        free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));

//  DEBUGPRINT("P%d is exiting here\n", MYID(simpar)); exit(9);

    postype TDtime;
    MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN,  MPI_COMM(simpar));
    Dtime = TDtime;

    int i;
    postype dmax=0;
//  postype eps2 = KH_EPS(simpar)*KH_EPS(simpar);
#ifdef _OPENMP
#pragma omp parallel for reduction(max: dmax)
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
			bp[i].ax += GAS_ACCX(simpar);
			bp[i].ay += GAS_ACCY(simpar);
	        bp[i].vx += bp[i].ax *Dtime; 
			bp[i].vy += bp[i].ay *Dtime; 
			bp[i].x += bp[i].vx *Dtime; 
			bp[i].y += bp[i].vy *Dtime; 
			dmax = MAX(dmax, fabs(bp[i].ay)); 
			bp[i].x = fmod(bp[i].x+Lx,Lx); 
			bp[i].y = fmod(bp[i].y+Ly,Ly); 
			bp[i].w2old = bp[i].w2;
		}
		else {
			bp[i].vx = bp[i].vy = bp[i].ax = bp[i].ay = 0;
		}
    }
    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);

// This is for centroid shift.
	if(GAS_FCENTROID(simpar) > 0 ){
	    if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingAllTreeParticles);
    // Update the linked list 
		mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

    // repositioning the memory space 
	    bp = VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	    for(iy=1;iy<my-1;iy++){ 
			int mp=2000; 
			Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp); 
			postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a; 
			int ix; 
			for(ix=1;ix<mx-1;ix++){ 
				int np; 
				treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np); 
				int nneigh; 
				Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh); 
				Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh); 
				int i; 
				for(i=0;i<np;i++){ 
					Voro2D_point center; 
					center.x = p[i].x; 
					center.y = p[i].y; 
					center.indx = PINDX(p+i); 
					center.w2 = p[i].w2; 
					
					treevorork4particletype *ibp = p[i].bp; 
					int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
					// in this case ax, ay are used for the centroid vector
					Voro2D_Corner centroid = Centroid2DPolygon(vorocorner);
					ibp->ax = centroid.x;
					ibp->ay = centroid.y; 
				} 
				free(p);free(neighbors); free(neighwork); 
			} 
			free(vorocorner); 
		} 
		my_free(VORORK4_TBPP(simpar));
		float fshift = GAS_FCENTROID(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	    for(i=0;i<VORO_NP(simpar);i++){ 
//			if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
			if(targetBP(bp+i, Lx, Ly))
			{
				bp[i].x += fshift*bp[i].ax; // move it to the centroid  with a factor of shift
				bp[i].y += fshift*bp[i].ay; 
   		     	bp[i].x = fmod(bp[i].x+Lx,Lx);
		        bp[i].y = fmod(bp[i].y+Ly,Ly);
			}
		} 
		// migrate particles between mpi ranks 
		migrateTreeVorork4Particles(simpar);
	}
//
//
	if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingAllTreeParticles);
    // Update the linked list
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);


    // repositioning the memory space
    bp = VORORK4_TBP(simpar);



#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.w2 = p[i].w2;

                treevorork4particletype *ibp = p[i].bp;
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
//              ibp->w2 = ibp->w2ceil = center.w2ceil;
                // These are to use register 
                postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
                ibp_vx = ibp->vx;
                ibp_vy = ibp->vy;
                ibp_csound = ibp->csound;
                ibp_pressure = ibp->pressure;
                ibp_den = ibp->den;

//                ibp->volume = Area2DPolygon(vorocorner, mp);
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);




                Voro2D_Corner *tmp,*tmp2;
                tmp = vorocorner;
                int j;
//                double die, dte,dke,fx,fy;
//               die = dte = dke = fx = fy = 0;
				double die, cx,cy;
				die = cx = cy = 0;
                do {
                    if(tmp->upperrelated >=0){
                        // These are to use register 
                        treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
                        postype jbp_vx = jbp->vx;
                        postype jbp_vy = jbp->vy;
                        postype jbp_csound = jbp->csound;
//                        postype jbp_pressure = jbp->pressure;
						postype jbp_pressure = getNeighborPressure(jbp,ibp);
                        postype jbp_den = jbp->den;

                        Voro2D_point line;
                        tmp2 = tmp->upperlink;
                        line.x = tmp2->x - tmp->x;
                        line.y = tmp2->y - tmp->y;
                        Voro2D_point dS = voro2D_norm(&line);
                        postype facearea = Vec2DLength(tmp,tmp2);
                        dS.x = facearea*dS.x;
                        dS.y = facearea*dS.y;
//                      postype pi = Half*(jbp_pressure + ibp_pressure);
                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,RT);

                        Voro2D_point uij,ua,ub;
                        uij.x = (jbp_vx - ibp_vx);
                        uij.y = (jbp_vy - ibp_vy);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &uij);
                        if(rvel<0){ 
							postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
							postype scaleFactor = (wcomp > etavis ? wcomp : etavis);
							postype drampScale = dramp/scaleFactor; 
							postype mu = rvel /* Voronoi-AV-reform: no length regularization */; 
							postype meanden = 0.5*(ibp_den + jbp_den); 
							postype meanCsound = 0.5*(ibp_csound + jbp_csound); 
							pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
//                ibp->ax = fx/ibp->mass;
//                ibp->ay = fy/ibp->mass;
            }
            free(p);free(neighbors); free(neighwork);
        }
        free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));
    my_free(VORO_BASICCELL(simpar));

    // update the w2ceil
    postype dmean = SIMBOX(simpar).x.max/NX(simpar);
//  DEBUGPRINT("P%d is now updating the internal energy\n", MYID(simpar));
	// split getw2forHydroParticle and pressure update because the former uses the latter.
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
//		if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
		if(targetBP(bp+i, Lx, Ly))
		{
			if(GAS_Kappa(simpar) <0) { 
				bp[i].w2 = -GAS_Kappa(simpar); 
			} 
			else if (GAS_Kappa(simpar) >0){ 
          		bp[i].w2hydro = getw2forHydroParticle(simpar,(bp+i),Dtime);
				applyW2Controls(simpar, bp+i, Dtime);
			} 
		}
	}
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
//		if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
		if(targetBP(bp+i, Lx, Ly))
		{
	        bp[i].ie += bp[i].die * Dtime; 
			bp[i].den  = bp[i].mass/bp[i].volume; 
			bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1); 
			bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den); 
			if(bp[i].volume ==0) 
				DEBUGPRINT("P%d has wrong volume at i= %d with den= %g at xy= %g %g\n", MYID(simpar), 
						i,bp[i].den, bp[i].x, bp[i].y); 
			if(isnan(bp[i].csound)){ 
				DEBUGPRINT("P%d has nan at the sound speed = %d with" 
						"id= %ld: pressure= %g den= %g Delta IE=%g dIE= %g Dt= %g\n", 
						MYID(simpar), i,PINDX(bp+i), bp[i].pressure, bp[i].den, bp[i].ie,bp[i].die,Dtime); 
				exit(9);
			}
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return Dtime;
}

double exam2d_vph_rk4_int(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	postype xmin,ymin,xmax,ymax,OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);
//	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has x/y minmax= %g %g : %g %g ::: %g\n", MYID(simpar), xmin,ymin,xmax,ymax, cellsize);

	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);


    int i;
    postype Dtime,dt;
    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) {
		bp[i].w2old = bp[i].w2;
		bp[i].rk4.w2backup = bp[i].w2;
	}
	// Runge-Kutta 4-th order time evolution of r and vr
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, 1);
//	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed updateDenW2P...()\n", MYID(simpar));
    Dtime = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){ 
		bp[i].rk4.k1x = bp[i].vx*Dtime;
        bp[i].rk4.k1y = bp[i].vy*Dtime;
        bp[i].rk4.k1vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k1vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k1ie = bp[i].die*Dtime;
    }
    // RK4 second
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += bp[i].rk4.k1x*0.5;
        bp[i].y += bp[i].rk4.k1y*0.5;
        bp[i].vx += bp[i].rk4.k1vx*0.5;
        bp[i].vy += bp[i].rk4.k1vy*0.5;
        bp[i].ie += bp[i].rk4.k1ie*0.5;
        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
	migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed updateDenW2P...()\n", MYID(simpar));
    dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
        paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k2x = bp[i].vx*Dtime;
        bp[i].rk4.k2y = bp[i].vy*Dtime;
        bp[i].rk4.k2vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k2vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k2ie = bp[i].die*Dtime;
    }
    // RK4 third
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += (bp[i].rk4.k2x-bp[i].rk4.k1x)*0.5;
        bp[i].y += (bp[i].rk4.k2y-bp[i].rk4.k1y)*0.5;
        bp[i].vx += (bp[i].rk4.k2vx-bp[i].rk4.k1vx)*0.5;
        bp[i].vy += (bp[i].rk4.k2vy-bp[i].rk4.k1vy)*0.5;
        bp[i].ie += (bp[i].rk4.k2ie-bp[i].rk4.k1ie)*0.5;

        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed updateDenW2P...()\n", MYID(simpar));
	dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k3x = bp[i].vx*Dtime;
        bp[i].rk4.k3y = bp[i].vy*Dtime;
        bp[i].rk4.k3vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k3vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k3ie = bp[i].die*Dtime;
    }
    // RK4 fourth
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += bp[i].rk4.k3x - bp[i].rk4.k2x*0.5; 
		bp[i].y += bp[i].rk4.k3y - bp[i].rk4.k2y*0.5; 
		bp[i].vx += bp[i].rk4.k3vx - bp[i].rk4.k2vx*0.5; 
		bp[i].vy += bp[i].rk4.k3vy - bp[i].rk4.k2vy*0.5; 
		bp[i].ie += bp[i].rk4.k3ie - bp[i].rk4.k2ie*0.5;

        bp[i].x = fmod(bp[i].x+Lx, Lx); 
		bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed updateDenW2P...()\n", MYID(simpar));
    dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k4x = bp[i].vx*Dtime;
        bp[i].rk4.k4y = bp[i].vy*Dtime;
        bp[i].rk4.k4vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k4vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k4ie = bp[i].die*Dtime;
    }
    // migrate particles between mpi ranks
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x -= bp[i].rk4.k3x;
        bp[i].y -= bp[i].rk4.k3y;
        bp[i].vx -= bp[i].rk4.k3vx;
        bp[i].vy -= bp[i].rk4.k3vy;
        bp[i].ie -= bp[i].rk4.k3ie;
        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
	// migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

    // to finalize ke & ie updates 
	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
    // Update position and velocity using RK4
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
//		if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
		if(targetBP(bp+i, Lx, Ly))
		{
			bp[i].x  += (bp[i].rk4.k1x +2*bp[i].rk4.k2x +2*bp[i].rk4.k3x +bp[i].rk4.k4x )/6.; 
			bp[i].y  += (bp[i].rk4.k1y +2*bp[i].rk4.k2y +2*bp[i].rk4.k3y +bp[i].rk4.k4y )/6.; 
			bp[i].vx += (bp[i].rk4.k1vx+2*bp[i].rk4.k2vx+2*bp[i].rk4.k3vx+bp[i].rk4.k4vx)/6.; 
			bp[i].vy += (bp[i].rk4.k1vy+2*bp[i].rk4.k2vy+2*bp[i].rk4.k3vy+bp[i].rk4.k4vy)/6.; 
			bp[i].ie += (bp[i].rk4.k1ie+2*bp[i].rk4.k2ie+2*bp[i].rk4.k3ie+bp[i].rk4.k4ie)/6.;

	        bp[i].x = fmod(bp[i].x+Lx,Lx); 
			bp[i].y = fmod(bp[i].y+Ly,Ly);
		}
    }

    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);

	//centroid shift
	if(GAS_FCENTROID(simpar)>0) 
		exam2d_centroidShift(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	// update volume & density with one-step evolved positions 
    exam2dUpdateVol(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);
	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed exam2dUpdateVol...()\n", MYID(simpar));

	bp = VORORK4_TBP(simpar);
	postype dmean = SIMBOX(simpar).x.max/NX(simpar);
	// split getw2forHydroParticle and pressure update because the former uses the latter.
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
			if(GAS_Kappa(simpar) <0) { 
				bp[i].w2 = -GAS_Kappa(simpar); 
			} 
			else if (GAS_Kappa(simpar) >0){ 
          		bp[i].w2hydro = getw2forHydroParticle(simpar,(bp+i),Dtime);
				applyW2Controls(simpar, bp+i, Dtime);
			} 
		}
	}
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
			bp[i].den = bp[i].mass/bp[i].volume; 
			bp[i].pressure = bp[i].ie/bp[i].volume * (Gamma-1); 
			bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den); 
			if(bp[i].volume ==0) 
				DEBUGPRINT("P%d has wrong volume at i= %d with den= %g at xy= %g %g\n", MYID(simpar), 
						i,bp[i].den, bp[i].x, bp[i].y); 
			if(isnan(bp[i].csound)){ 
				DEBUGPRINT("P%d has nan at the sound speed = %d with" 
						"id= %zu: pressure= %g den= %g Delta IE=%g dIE= %g Dt= %g\n", 
						MYID(simpar), i,PINDX(bp+i), bp[i].pressure, bp[i].den, bp[i].ie,bp[i].die,Dtime); 
				exit(9);
			}
		}
		else {
			bp[i].vx = bp[i].vy = 0;
			bp[i].w2= bp[i].rk4.w2backup;
		}
    }
    return Dtime;
}

/* ================================================================
 *  periodic_postStage_kNN / periodic_postStage_blend:
 *  Default postStage callbacks for periodic boundary models (KH etc.).
 *  Wrap x,y and migrate particles between MPI ranks.
 * ================================================================ */
void periodic_postStage_kNN(SimParameters *simpar){
	postype Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	postype Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *raw = (char*)VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;
	for(i=0;i<np;i++){
		treevorork4particletype *bp = (treevorork4particletype*)(raw + i*p_size);
		bp->x = fmod(bp->x + Lx, Lx);
		bp->y = fmod(bp->y + Ly, Ly);
	}
	migrateTreeVorork4Particles(simpar);
}

void periodic_postStage_blend(SimParameters *simpar){
	postype Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	postype Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;
	for(i=0;i<np;i++){
		sbp[i].x = fmod(sbp[i].x + Lx, Lx);
		sbp[i].y = fmod(sbp[i].y + Ly, Ly);
	}
	migrateTreeVorork4Particles(simpar);
}

/* walls_xy_postStage_blend: reflecting walls on all four sides.
 * Used by the 2D Sedov blast wave. Mirror logic identical to wallx but
 * applied to both x and y; both axes also clamp near-wall particles
 * to half-IC-spacing inside the box to prevent density pile-up.
 *
 * Sedov2D-specific: a particle with indx == nx*ny is the central anchor
 * (sentinel value, never produced by regular IC layouts). It is pinned
 * to the box center with zero velocity to keep the blast geometric center
 * fixed during the rebound. */
void walls_xy_postStage_blend(SimParameters *simpar){
	postype xlo = SIMBOX(simpar).x.min, xhi = SIMBOX(simpar).x.max;
	postype ylo = SIMBOX(simpar).y.min, yhi = SIMBOX(simpar).y.max;
	postype dx_min = 0.5 * (xhi - xlo) / (postype)NX(simpar);
	postype dy_min = 0.5 * (yhi - ylo) / (postype)NY(simpar);
	postype cx = 0.5 * (xlo + xhi);
	postype cy = 0.5 * (ylo + yhi);
	size_t anchor_idx = (size_t)NX(simpar) * (size_t)NY(simpar);
	treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;
	for(i=0;i<np;i++){
		if(PINDX(sbp+i) == anchor_idx){
			/* Phase 1 rebuilds ie from E−kinetic. Zeroing the pin velocity
			 * would drop that kinetic energy out of the budget. */
			if(sedov_phase1_on()){
				postype m = sbp[i].mass;
				sbp[i].ie += (postype)0.5 * m * (sbp[i].vx*sbp[i].vx + sbp[i].vy*sbp[i].vy);
			}
			sbp[i].x = cx; sbp[i].y = cy;
			sbp[i].vx = 0; sbp[i].vy = 0;
			sbp[i].ax = 0; sbp[i].ay = 0;
			continue;
		}
		if(sbp[i].x < xlo){ sbp[i].x = 2.0*xlo - sbp[i].x; sbp[i].vx = -sbp[i].vx; }
		if(sbp[i].x > xhi){ sbp[i].x = 2.0*xhi - sbp[i].x; sbp[i].vx = -sbp[i].vx; }
		if(sbp[i].y < ylo){ sbp[i].y = 2.0*ylo - sbp[i].y; sbp[i].vy = -sbp[i].vy; }
		if(sbp[i].y > yhi){ sbp[i].y = 2.0*yhi - sbp[i].y; sbp[i].vy = -sbp[i].vy; }
		if(sbp[i].x < xlo + dx_min){ sbp[i].x = xlo + dx_min; if(sbp[i].vx < 0) sbp[i].vx = 0; }
		if(sbp[i].x > xhi - dx_min){ sbp[i].x = xhi - dx_min; if(sbp[i].vx > 0) sbp[i].vx = 0; }
		if(sbp[i].y < ylo + dy_min){ sbp[i].y = ylo + dy_min; if(sbp[i].vy < 0) sbp[i].vy = 0; }
		if(sbp[i].y > yhi - dy_min){ sbp[i].y = yhi - dy_min; if(sbp[i].vy > 0) sbp[i].vy = 0; }
	}
	migrateTreeVorork4Particles(simpar);
}

/* wallx_postStage_blend: reflecting walls in x, periodic in y.
 * Reflect any overshoot, then clamp near-wall particles so they keep at least
 * dx_min away from each wall — otherwise the wall mirror sits ~2·x away and
 * halves the Voronoi cell width, producing a density pile-up at the wall.
 * dx_min is taken to be 0.5 * (Lx / NX), i.e. half the IC grid spacing. */
void wallx_postStage_blend(SimParameters *simpar){
	postype xlo = SIMBOX(simpar).x.min, xhi = SIMBOX(simpar).x.max;
	postype Ly  = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	postype dx_ic = (xhi - xlo) / (postype)NX(simpar);
	postype dx_min = 0.5 * dx_ic;
	treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;
	for(i=0;i<np;i++){
		if(sbp[i].x < xlo){ sbp[i].x = 2.0*xlo - sbp[i].x; sbp[i].vx = -sbp[i].vx; }
		if(sbp[i].x > xhi){ sbp[i].x = 2.0*xhi - sbp[i].x; sbp[i].vx = -sbp[i].vx; }
		if(sbp[i].x < xlo + dx_min){ sbp[i].x = xlo + dx_min; if(sbp[i].vx < 0) sbp[i].vx = 0; }
		if(sbp[i].x > xhi - dx_min){ sbp[i].x = xhi - dx_min; if(sbp[i].vx > 0) sbp[i].vx = 0; }
		sbp[i].y = fmod(sbp[i].y + Ly, Ly);
	}
	migrateTreeVorork4Particles(simpar);
}

/* dK rate-limiter (entropy_mode=1 only): clamp |dK|*dt_ref/K <= rate_max.
 * Uses first-step Dtime as dt_ref so clamp does NOT loosen when dt collapses
 * during a runaway. Bounds K growth absolutely (per fixed reference time).
 * Off when GAS_DKRATEMAX==0. */
static void clamp_dK_rate_limiter(SimParameters *simpar, postype Dtime)
{
	static postype dt_ref = 0;
	static int entry_diag_done = 0;
	postype rate_max = (postype)GAS_DKRATEMAX(simpar);
	if(!entry_diag_done){
		if(MYID(simpar) == 0){
			fprintf(stderr, "[DK_LIMIT_ENTRY] rate_max=%.6f em=%d Dtime=%.3e K_floor=%.3e\n",
				(double)rate_max, GAS_ENTROPY_MODE(simpar),
				(double)Dtime, (double)GAS_K_FLOOR(simpar));
			fflush(stderr);
		}
		entry_diag_done = 1;
	}
	if(rate_max <= 0 || GAS_ENTROPY_MODE(simpar) != 1) return;
	if(dt_ref <= 0){
		if(Dtime > 0){
			dt_ref = Dtime;
			if(MYID(simpar) == 0){
				fprintf(stderr, "[DK_LIMIT] dt_ref initialized = %.3e (rate_max=%.3f)\n",
					(double)dt_ref, (double)rate_max);
				fflush(stderr);
			}
		} else {
			return;
		}
	}
	postype K_floor_rl = (postype)GAS_K_FLOOR(simpar);
	treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int n = VORO_NP(simpar);
	int i;
	int nclamp = 0;
	postype absdK_pre = 0, absdK_post = 0;
	for(i=0;i<n;i++){
		postype Ki = (postype)sbp[i].stress.K;
		postype Ksafe = (Ki > K_floor_rl) ? Ki : K_floor_rl;
		postype dKmax = rate_max * Ksafe / dt_ref;
		postype dKi = (postype)sbp[i].stress.dK;
		postype adKi = (dKi >= 0) ? dKi : -dKi;
		if(adKi > absdK_pre) absdK_pre = adKi;
		if(dKi >  dKmax){ sbp[i].stress.dK = (float)( dKmax); nclamp++; }
		else if(dKi < -dKmax){ sbp[i].stress.dK = (float)(-dKmax); nclamp++; }
		postype dKf = (postype)sbp[i].stress.dK;
		postype adKf = (dKf >= 0) ? dKf : -dKf;
		if(adKf > absdK_post) absdK_post = adKf;
	}
	{
		static int rate_diag_count = 0;
		if(rate_diag_count < 20 && MYID(simpar) == 0){
			fprintf(stderr, "[DK_LIMIT_TICK] nclamp=%d/%d absdK pre=%.3e post=%.3e dt_ref=%.3e\n",
				nclamp, n, (double)absdK_pre, (double)absdK_post, (double)dt_ref);
			fflush(stderr);
			rate_diag_count++;
		}
	}
	(void)nclamp;
}

/* ================================================================
 *  exam2d_vph_rk4_int_blend:
 *  RK4 integrator with NS stress + HLLC + CD10 two-tier blending.
 *  Uses treevorostressrk4particletype for stress fields.
 *  Same RK4 structure as exam2d_vph_rk4_int, but calls
 *  updateDenW2Pressure2DBlend and getAccVoro2DBlend.
 * ================================================================ */
/* GFS_W2_LINEAR helpers (see gfs_w2lin_on). */
static long gfs_w2lin_nclip = 0;
static void gfs_w2lin_start(SimParameters *simpar, treevorostressrk4particletype *sbp, int np,
		int (*targetBP)(treevorork4particletype*, postype, postype), postype Lx, postype Ly){
	postype dto = GAS_dtold(simpar);
	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<np;i++){
		treevorork4particletype *b = (treevorork4particletype*)(sbp+i);
		postype Wn = b->w2;
		postype rate = 0;
		if(targetBP(b, Lx, Ly) && dto > 0){
			b->w2hydro = getw2forHydroParticle(simpar, b, dto);
			applyW2Controls(simpar, b, dto);   /* sets b->w2 = target (<= w2ceil) */
			postype Wt = b->w2;
			/* Margin: target <= w2ceil/2, so the ceiling clip in det2d_dpqRK4
			 * (a weight change that P dV does not see) needs the nearest
			 * distance to shrink by ~30% within one step. */
			if(Wt > 0.5*b->w2ceil) Wt = 0.5*b->w2ceil;
			rate = (Wt - Wn)/dto;
			postype wm = sqrt((double)(Wn > Wt ? Wn : Wt));
			postype rmax = 2*b->csound*wm;
			if(rate >  rmax) rate =  rmax;
			if(rate < -rmax) rate = -rmax;
		}
		b->w2 = Wn;
		b->w2old = Wn;
		b->w2hydro = rate;
	}
}
/* W = W_n + cdt*Wdot (W_n in rk4.w2backup), w2old = w2. */
static void gfs_w2lin_set(treevorostressrk4particletype *sbp, int np, postype cdt){
	int i; long nc = 0;
	for(i=0;i<np;i++){
		postype W = sbp[i].rk4.w2backup + cdt*sbp[i].w2hydro;
		if(W < 0){ W = 0; nc++; }
		if(W > sbp[i].w2ceil) nc++;   /* det2d_dpqRK4 will clip it */
		sbp[i].w2 = W;
		sbp[i].w2old = W;
	}
	gfs_w2lin_nclip += nc;
}
/* GFS_VOL_AUDIT: add w*dt*(booked dV/dt) of this stage. */
static void vaud_stage(treevorostressrk4particletype *sbp, int np, postype wdt){
	int i;
	for(i=0;i<np;i++){
		size_t id = (size_t)PINDX((treevorork4particletype*)(sbp+i));
		if(id < vaud_cap){
			vaud_book[id] += (double)wdt*vaud_rate[id];
			double ier = (double)phase1_ie_stage(sbp+i, (postype)1.0), P = (double)sbp[i].pressure;
			vaud_diss[id] += (double)wdt*(ier + P*vaud_rate[id]);
			vaud_adb[id]  += (double)wdt*(-P*vaud_rate[id]);
		}
	}
}
/* GFS_VOL_AUDIT: close the previous step against the current tessellation. */
static void vaud_close(SimParameters *simpar, treevorostressrk4particletype *sbp, int np,
		int (*targetBP)(treevorork4particletype*, postype, postype), postype Lx, postype Ly){
	static int st = 0; static double pw_cum = 0, pws_cum = 0;
	static double c_dI = 0, c_diss = 0, c_adb = 0, c_ag = 0, c_absgcl = 0;
	long vdec_n = 0; double vdec_dI = 0, vdec_diss = 0, vdec_adb = 0, vdec_ag = 0, vdec_absgcl = 0;
	double ra = 0, vs = 0, pw = 0, pws = 0, wb = 0; int i, nres = 0;
	size_t max_id = 0;
	for(i=0;i<np;i++){ size_t id = (size_t)PINDX((treevorork4particletype*)(sbp+i)); if(id > max_id) max_id = id; }
	vaud_ensure(max_id);
	for(i=0;i<np;i++){
		size_t id = (size_t)PINDX((treevorork4particletype*)(sbp+i));
		double V = (double)sbp[i].volume;
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly) && vaud_have[id]){
			double r = V - vaud_vprev[id] - vaud_book[id];
			ra += fabs(r); vs += V; nres++;
			pw += (double)sbp[i].pressure*fabs(r);
			pws += (double)sbp[i].pressure*r;
			wb += fabs(vaud_book[id]);
		}
		if(gfs_gcl_diag_on() && targetBP((treevorork4particletype*)(sbp+i), Lx, Ly) && vaud_have[id] && V > 0){
			/* GFS_GCL_DIAG (diagnostic, not conservative): make the adiabatic
			 * work of the last step match the geometric volume change,
			 * dI = -P (dV_geom - dV_booked), then refresh P and c. */
			double r = V - vaud_vprev[id] - vaud_book[id];
			double ie_new = (double)sbp[i].ie - (double)sbp[i].pressure*r;
			if(ie_new > 0){
				gfs_gcl_diag_de += ie_new - (double)sbp[i].ie;
				sbp[i].ie = (postype)ie_new;
				sbp[i].pressure = sbp[i].ie/sbp[i].volume*(GAS_GAMMA(simpar)-1);
				sbp[i].csound = sqrt(GAS_GAMMA(simpar)*sbp[i].pressure/sbp[i].den);
			}
		}
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly) && vaud_have[id] && !vaud_arcp[id] && !vaud_arc[id]){
			double dI = (double)sbp[i].ie - vaud_ie0[id], Pm = 0.5*((double)sbp[i].pressure + vaud_p0[id]);
			double ag = -Pm*(V - vaud_vprev[id]);
			vdec_n++; vdec_dI += dI; vdec_diss += vaud_diss[id]; vdec_adb += vaud_adb[id]; vdec_ag += ag;
			vdec_absgcl += fabs(vaud_adb[id] - ag);
		}
		vaud_ie0[id] = (double)sbp[i].ie; vaud_p0[id] = (double)sbp[i].pressure;
		vaud_diss[id] = vaud_adb[id] = 0; vaud_arcp[id] = vaud_arc[id]; vaud_arc[id] = 0;
		vaud_vprev[id] = V; vaud_book[id] = 0; vaud_have[id] = 1;
	}
	if(nres > 0 && MYID(simpar) == 0){
		st++; pw_cum += pw; pws_cum += pws;
		fprintf(stderr, "[VAUD] step=%d n=%d res_abs/V=%.3e booked_abs/V=%.3e unbooked_PdV_abs=%.3e signed=%.3e cum_abs=%.3e cum_signed=%.3e w2clip=%ld\n",
			st, nres, ra/vs, wb/vs, pw, pws, pw_cum, pws_cum, gfs_w2lin_nclip);
		if(gfs_gcl_diag_on()) fprintf(stderr, "[GCLDIAG] step=%d cum_dI=%.3e\n", st, gfs_gcl_diag_de);
		c_dI += vdec_dI; c_diss += vdec_diss; c_adb += vdec_adb; c_ag += vdec_ag; c_absgcl += vdec_absgcl;
		fprintf(stderr, "[VDEC] step=%d n_int=%ld cum: dI=%.4e diss=%.4e adb_booked=%.4e adb_geom=%.4e gcl(booked-geom)=%.4e |gcl|=%.4e resid(dI-adb_b-diss)=%.4e\n",
			st, vdec_n, c_dI, c_diss, c_adb, c_ag, c_adb - c_ag, c_absgcl, c_dI - c_adb - c_diss);
	}
}

/* HG split step of length tau: all modes forward with tau/2, then backward with tau/2.
 * One mode: r0 = sum l.v, alpha = sum |l|^2/m; with s = ln(q0/q), ds/dt = alpha C(q, U),
 * C = A0 (b q + sqrt(b^2 q^2 + kap U)), U = U0 + q0^2 (1 - e^{-2s})/(2 alpha); then
 * v_a -= (r0 - r1) l_a/(alpha m_a), U += (r0^2 - r1^2)/(2 alpha) (exact energy transfer). */
static double hgs_heat = 0; static long hgs_nmode = 0, hgs_bad = 0;
static void bp_apply_mode(hg_mode *M, double tau);
static void hg_apply_mode(long mi, double tau){
	hg_mode *M = hgm + mi; long z; double r0 = 0, al = 0;
	if(M->kind == 2){ bp_apply_mode(M, tau); return; }
	treevorork4particletype *ow = ((size_t)M->owner < hg_cap) ? hg_ptr[M->owner] : NULL;
	if(!ow){ hgs_bad++; return; }
	for(z = 0; z < M->len; z++){
		long id = hgm_id[M->start + z]; treevorork4particletype *q = ((size_t)id < hg_cap) ? hg_ptr[id] : NULL;
		if(!q){ hgs_bad++; return; }
		double lx_ = hgm_l[2*(M->start + z)], ly_ = hgm_l[2*(M->start + z) + 1];
		r0 += lx_*(double)q->vx + ly_*(double)q->vy; al += (lx_*lx_ + ly_*ly_)/(double)q->mass;
	}
	if(!(al > 0) || r0 == 0) return;
	if(M->onesided && r0 > 0) return;   /* bulk acts on compression only */
	double U0 = (double)ow->ie, q0 = fabs(r0);
	if(!(U0 > 0)){ hgs_bad++; return; }
	#define HGS_F(s_) ({ double qq_ = q0*exp(-(s_)), UU_ = U0 + q0*q0*(-expm1(-2*(s_)))/(2*al); \
		al*M->A0*(M->b*qq_ + (M->onesided ? sqrt(M->kap*UU_) : sqrt(M->b*M->b*qq_*qq_ + M->kap*UU_))); })
	double f0 = HGS_F(0.0);
	int N = (int)ceil(f0*tau/0.25); if(N < 1) N = 1; if(N > 4096) N = 4096;
	double h = tau/N, sv = 0; int it;
	for(it = 0; it < N; it++){
		double k1 = HGS_F(sv), k2 = HGS_F(sv + 0.5*h*k1), k3 = HGS_F(sv + 0.5*h*k2), k4 = HGS_F(sv + h*k3);
		sv += h*(k1 + 2*k2 + 2*k3 + k4)/6;
	}
	#undef HGS_F
	double q1 = q0*exp(-sv), r1 = (r0 > 0) ? q1 : -q1;
	double dU = q0*q0*(-expm1(-2*sv))/(2*al), D = r0 - r1;
	for(z = 0; z < M->len; z++){
		long id = hgm_id[M->start + z]; treevorork4particletype *q = hg_ptr[id];
		double lx_ = hgm_l[2*(M->start + z)], ly_ = hgm_l[2*(M->start + z) + 1], c_ = D/(al*(double)q->mass);
		q->vx -= (postype)(c_*lx_); q->vy -= (postype)(c_*ly_);
	}
	ow->ie += (postype)dU;
	if(tcf_sacc_file() && (size_t)M->owner < tcf_cap) tcf_sacc[M->owner] += log1p(dU/U0);
	hgs_heat += dU; hgs_nmode++; if(M->onesided) bulk_heat += dU;
}
/* b' mode: two owners, exact scalar flow in s = ln(q0/q) */
static void bp_apply_mode(hg_mode *M, double tau){
	long z; double r0 = 0, al = 0;
	treevorork4particletype *o1 = ((size_t)M->owner < hg_cap) ? hg_ptr[M->owner] : NULL, *o2 = ((size_t)M->owner2 < hg_cap) ? hg_ptr[M->owner2] : NULL;
	if(!o1 || !o2){ hgs_bad++; return; }
	for(z = 0; z < M->len; z++){
		long id = hgm_id[M->start + z]; treevorork4particletype *q = ((size_t)id < hg_cap) ? hg_ptr[id] : NULL;
		if(!q){ hgs_bad++; return; }
		double lx_ = hgm_l[2*(M->start + z)], ly_ = hgm_l[2*(M->start + z) + 1];
		r0 += lx_*(double)q->vx + ly_*(double)q->vy; al += (lx_*lx_ + ly_*ly_)/(double)q->mass;
	}
	if(!(al > 0) || r0 == 0) return;
	double U1 = (double)o1->ie, U2 = (double)o2->ie, q0 = fabs(r0), sg = (r0 > 0) ? 1 : -1;
	if(!(U1 > 0) || !(U2 > 0)){ hgs_bad++; return; }
	if(gfs_bprime() == 2){
		double kf = M->A0*M->mu*M->b, beta = al*kf*q0*tau, sig = beta/(1 + beta);
		double Dd = r0*sig, dUq = q0*q0/(4*al)*sig*(2 - sig);
		for(z = 0; z < M->len; z++){
			long id = hgm_id[M->start + z]; treevorork4particletype *q = hg_ptr[id];
			double lx_ = hgm_l[2*(M->start + z)], ly_ = hgm_l[2*(M->start + z) + 1], c_ = Dd/(al*(double)q->mass);
			q->vx -= (postype)(c_*lx_); q->vy -= (postype)(c_*ly_);
		}
		o1->ie += (postype)dUq; o2->ie += (postype)dUq;
		if(tcf_sacc_file()){ if((size_t)M->owner < tcf_cap) tcf_sacc[M->owner] += log1p(dUq/U1); if((size_t)M->owner2 < tcf_cap) tcf_sacc[M->owner2] += log1p(dUq/U2); }
		hgs_heat += 2*dUq; hgs_nmode++; bp_heat += 2*dUq;
		return;
	}
	#define BPS_F(s_) ({ double qq_ = q0*exp(-(s_)), dU_ = q0*q0*(-expm1(-2*(s_)))/(4*al); \
		double cM_ = fmax(sqrt(M->kap*(U1 + dU_)), sqrt(M->kap2*(U2 + dU_))), aH_ = M->b*qq_; \
		al*M->A0*fmax(M->mu*(cM_ + fmax(sg*qq_, 0.0)), M->mu*(aH_ + sqrt(aH_*aH_ + cM_*cM_))); })
	double f0 = BPS_F(0.0);
	int N = (int)ceil(f0*tau/0.25); if(N < 1) N = 1; if(N > 4096) N = 4096;
	double h = tau/N, sv = 0; int it;
	for(it = 0; it < N; it++){
		double k1 = BPS_F(sv), k2 = BPS_F(sv + 0.5*h*k1), k3 = BPS_F(sv + 0.5*h*k2), k4 = BPS_F(sv + h*k3);
		sv += h*(k1 + 2*k2 + 2*k3 + k4)/6;
	}
	#undef BPS_F
	double q1 = q0*exp(-sv), r1 = sg*q1, dU = q0*q0*(-expm1(-2*sv))/(4*al), D = r0 - r1;
	for(z = 0; z < M->len; z++){
		long id = hgm_id[M->start + z]; treevorork4particletype *q = hg_ptr[id];
		double lx_ = hgm_l[2*(M->start + z)], ly_ = hgm_l[2*(M->start + z) + 1], c_ = D/(al*(double)q->mass);
		q->vx -= (postype)(c_*lx_); q->vy -= (postype)(c_*ly_);
	}
	o1->ie += (postype)dU; o2->ie += (postype)dU;
	if(tcf_sacc_file()){ if((size_t)M->owner < tcf_cap) tcf_sacc[M->owner] += log1p(dU/U1); if((size_t)M->owner2 < tcf_cap) tcf_sacc[M->owner2] += log1p(dU/U2); }
	hgs_heat += 2*dU; hgs_nmode++; bp_heat += 2*dU;
}
/* assemble b' face modes from the per-cell slope operators and face records (sequential, after the build pass) */
static void bp_assemble(double gam){
	static long *mid = NULL; static double *mlx = NULL, *mly = NULL; static long mcap = 0;
	int test = getenv("GFS_BP_TEST") != NULL; size_t i;
	const double A_[4] = {0.37, -1.21, 0.83, 0.52}, b0x = 0.3, b0y = -0.7;
	for(i = 0; i < bp_cap; i++){
		int k;
		for(k = 0; k < bp_nf[i]; k++){
			long j = bp_fj[i*BP_MAXF + k]; double *Fi = bp_f + (i*BP_MAXF + k)*7;
			if(j <= (long)i || (size_t)j >= bp_cap) continue;
			int kj, mk = -1;
			for(kj = 0; kj < bp_nf[j]; kj++){
				double *Fj = bp_f + (j*BP_MAXF + kj)*7;
				if(bp_fj[j*BP_MAXF + kj] == (long)i && Fi[0]*Fj[0] + Fi[1]*Fj[1] < -0.999 && fabs(Fi[2] - Fj[2]) < 1e-6*Fi[2]){ mk = kj; break; }
			}
			double *Fj = (mk >= 0) ? bp_f + (j*BP_MAXF + mk)*7 : NULL;
			if(!Fj || !(Fi[6] > 0) || !(Fj[6] > 0) || bp_sl_len[i] < 0 || bp_sl_len[j] < 0){ bp_nskip++; continue; }
			long need = 2 + bp_sl_len[i] + bp_sl_len[j];
			if(need > mcap){ mcap = need; mid = (long*)realloc(mid, mcap*sizeof(long)); mlx = (double*)realloc(mlx, mcap*sizeof(double)); mly = (double*)realloc(mly, mcap*sizeof(double)); }
			double ex = Fi[0], ey = Fi[1], d = Fi[2]; int nm = 0, side;
			double tdel = 0;
			#define BP_PUT(ID_, CX_, CY_) do { int z_, f_ = -1; for(z_ = 0; z_ < nm; z_++) if(mid[z_] == (ID_)){ f_ = z_; break; } \
				if(f_ < 0){ mid[nm] = (ID_); mlx[nm] = (CX_); mly[nm] = (CY_); nm++; } else { mlx[f_] += (CX_); mly[f_] += (CY_); } } while(0)
			BP_PUT((long)i, ex, ey); BP_PUT(j, -ex, -ey);
			if(test){ double vjx = b0x + A_[0]*d*ex + A_[1]*d*ey, vjy = b0y + A_[2]*d*ex + A_[3]*d*ey; tdel += (b0x - vjx)*ex + (b0y - vjy)*ey; }
			for(side = 0; side < 2; side++){
				size_t c = side ? (size_t)j : i; long st = bp_sl_start[c]; int ln = bp_sl_len[c], z;
				double ox = side ? d*ex : 0, oy = side ? d*ey : 0;
				for(z = 0; z < ln; z++){
					double *sz = bp_s + 4*(st + z), w = sz[2]*ex + sz[3]*ey, cc = 0.5*d*w;
					BP_PUT(bp_sid[st + z], cc*ex, cc*ey);
					if(test){ double yx = ox + sz[0], yy = oy + sz[1]; tdel += cc*(ex*(b0x + A_[0]*yx + A_[1]*yy) + ey*(b0y + A_[2]*yx + A_[3]*yy)); }
				}
			}
			#undef BP_PUT
			if(test && fabs(tdel) > bp_test_err) bp_test_err = fabs(tdel);
			if(hgm_n >= hgm_cap || hgm_on + nm > hgm_ocap){ fprintf(stderr, "[BPRIME] mode store overflow\n"); exit(3); }
			long ms = hgm_n++, os = hgm_on; hgm_on += nm;
			hgm[ms].owner = (long)i; hgm[ms].owner2 = j; hgm[ms].key = (long)i*256 + 128 + k; hgm[ms].start = os; hgm[ms].len = nm;
			hgm[ms].A0 = 0.5*(Fi[3] + Fj[3]); hgm[ms].b = (gam + 1)/8; hgm[ms].mu = Fi[4]*Fj[4]/(Fi[4] + Fj[4]);
			hgm[ms].kap = Fi[5]; hgm[ms].kap2 = Fj[5]; hgm[ms].kind = 2; hgm[ms].onesided = 0;
			int z1; for(z1 = 0; z1 < nm; z1++){ hgm_id[os + z1] = mid[z1]; hgm_l[2*(os + z1)] = mlx[z1]; hgm_l[2*(os + z1) + 1] = mly[z1]; }
			bp_nmode++;
		}
	}
}
static void hg_apply(SimParameters *simpar, double tau){
	long m;
	for(m = 0; m < hgm_n; m++) hg_apply_mode(m, 0.5*tau);
	for(m = hgm_n - 1; m >= 0; m--) hg_apply_mode(m, 0.5*tau);
	/* pressure and sound speed from the new internal energies */
	size_t ps = TVORORK4_DDINFO(simpar)[0].n_size; char *b = (char*)VORORK4_TBP(simpar); int i; double gam = (double)GAS_GAMMA(simpar);
	for(i = 0; i < VORO_NP(simpar); i++){
		treevorork4particletype *q = (treevorork4particletype*)(b + i*ps);
		if(q->volume > 0 && q->ie > 0){ q->pressure = (postype)((gam - 1)*(double)q->ie/(double)q->volume); q->csound = (postype)sqrt(gam*(double)q->pressure/(double)q->den); }
	}
}
/* Step rejection and halving (Astra audit 2): the step is repeated from its starting state
 * with half the time step when any force evaluation of the step returned a smaller limit
 * than the step used, or a stage floor / non-positive state occurred. Active with
 * GFS_HG_AFFINE or GFS_STEP_RETRY=1. */
static int gfs_step_retry_on(void){ static int on = -1; if(on < 0){ const char *s_ = getenv("GFS_STEP_RETRY"); on = (s_ && s_[0] && atoi(s_) > 0) ? 1 : 0; } return on; }
static double exam2d_vph_rk4_int_blend_inner(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		void (*postStage)(SimParameters *)
		);
double exam2d_vph_rk4_int_blend(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		void (*postStage)(SimParameters *)
		){
	if(!(gfs_hg_affine() || gfs_step_retry_on()))
		return exam2d_vph_rk4_int_blend_inner(simpar, paddingAllTreeParticles, measureW2, find2DNeighborBP, find2DCellBP, targetBP, mkLinkedList2D, postStage);
	if(NID(simpar) > 1){ fprintf(stderr, "[RETRY] step rejection supports one MPI rank\n"); exit(1); }
	int split = gfs_hg_affine() && gfs_hg_split();
	postype cellsize = BASICCELL_CELLWIDTH(simpar), Gamma = GAS_GAMMA(simpar), Courant = GAS_COURANT(simpar);
	postype xmin = Xmin_HydroExam(simpar)-cellsize, ymin = Ymin_HydroExam(simpar)-cellsize;
	postype xmax = Xmax_HydroExam(simpar)+cellsize, ymax = Ymax_HydroExam(simpar)+cellsize;
	postype OrderOfAccuracy = VoroAccuracyOrder(simpar);
	size_t ps = TVORORK4_DDINFO(simpar)[0].n_size; int np0 = VORO_NP(simpar);
	if(gfs_tcf_on() && tcf_sacc_file()) tcf_sacc_ensure((treevorostressrk4particletype*)VORORK4_TBP(simpar), np0);
	char *save = (char*)malloc(ps*(size_t)np0);
	memcpy(save, (char*)VORORK4_TBP(simpar), ps*(size_t)np0);
	double *ssave = NULL; size_t scap = tcf_cap;
	if(tcf_sacc && scap){ ssave = (double*)malloc(scap*sizeof(double)); memcpy(ssave, tcf_sacc, scap*sizeof(double)); }
	static double fdt = -1;
	if(fdt < 0){ const char *e_ = getenv("GFS_FIXED_DT"); fdt = (e_ && e_[0]) ? atof(e_) : 0; }
	int tries;
	double cap = 0, dt = 0;
	gfs_in_retry = 1;
	for(tries = 0; tries < 40; tries++){
		gfs_attempt_failed = 0;
		double dtH = 0;
		hgs_heat = 0; hgs_nmode = 0; hgs_bad = 0;
		if(split){
			/* modes and step at x_n, then HG(dt/2) */
			hgm_reset((long)np0); hgm_building = 1; gfs_step_dtmin = 1e300;
			updateDenW2Pressure2DBlend(simpar, xmin, ymin, xmax, ymax, Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP, mkLinkedList2D, 1);
			double dpr = getAccVoro2DBlend(simpar, xmin, ymin, xmax, ymax, OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP, mkLinkedList2D);
			hgm_building = 0;
			if(gfs_bprime()) bp_assemble((double)GAS_GAMMA(simpar));
			qsort(hgm, hgm_n, sizeof(hg_mode), hgm_cmp);
			dtH = dpr*(1 - 1e-3);   /* margin for the change of the RK4 limit across HG(dt/2) (audit item) */
			if(fdt > 0){ if(fdt > dpr){ fprintf(stderr, "[HGSPLIT] GFS_FIXED_DT=%g exceeds the limit %g: test invalid\n", fdt, dpr); exit(5); } dtH = fdt; }
			if(cap > 0 && dtH > cap) dtH = cap;
			hg_apply(simpar, 0.5*dtH);
			gfs_dt_cap = dtH;
		} else gfs_dt_cap = cap;
		gfs_step_dtmin = 1e300;
		dt = exam2d_vph_rk4_int_blend_inner(simpar, paddingAllTreeParticles, measureW2, find2DNeighborBP, find2DCellBP, targetBP, mkLinkedList2D, postStage);
		int bad = gfs_attempt_failed ? 1 : 0, i, over = 0;
		if(split && dt < dtH*(1 - 1e-9)) over = 1;   /* the RK4 part could not use the HG step */
		if(!bad && !over && split){
			/* modes at x_{n+1}, then HG(dt/2) */
			hgm_reset((long)np0); hgm_building = 1;
			updateDenW2Pressure2DBlend(simpar, xmin, ymin, xmax, ymax, Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP, mkLinkedList2D, 1);
			(void)getAccVoro2DBlend(simpar, xmin, ymin, xmax, ymax, OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP, mkLinkedList2D);
			hgm_building = 0;
			if(gfs_bprime()) bp_assemble((double)GAS_GAMMA(simpar));
			qsort(hgm, hgm_n, sizeof(hg_mode), hgm_cmp);
			hg_apply(simpar, 0.5*dtH);
			if(hgs_bad > 0) bad = 4;
		}
		char *b = (char*)VORORK4_TBP(simpar);
		if(VORO_NP(simpar) != np0) bad = 2;
		else for(i = 0; i < np0 && !bad; i++){
			treevorork4particletype *q = (treevorork4particletype*)(b + i*ps);
			if(!(q->ie > 0) || !(q->volume > 0) || !isfinite((double)q->vx) || !isfinite((double)q->vy) || !isfinite((double)q->x) || !isfinite((double)q->y)) bad = 3;
		}
		if(!split && dt > gfs_step_dtmin*(1 + 1e-9)) over = 1;
		if(!bad && !over) break;
		if(fdt > 0){ fprintf(stderr, "[RETRY] fixed dt %g rejected (bad=%d over=%d): test invalid\n", fdt, bad, over); exit(5); }
		double ncap = bad ? 0.5*(split ? dtH : dt) : (split ? dt : gfs_step_dtmin)*(1 - 1e-9);
		fprintf(stderr, "[RETRY] step rejected: dt=%.4e limit=%.4e bad=%d over=%d -> retry with dt <= %.4e\n", split ? dtH : dt, split ? dt : gfs_step_dtmin, bad, over, ncap);
		memcpy((char*)VORORK4_TBP(simpar), save, ps*(size_t)np0); VORO_NP(simpar) = np0;
		if(ssave) memcpy(tcf_sacc, ssave, scap*sizeof(double));
		cap = ncap;
	}
	gfs_in_retry = 0;
	if(tries >= 40){ fprintf(stderr, "[RETRY] step could not be accepted after 40 retries\n"); exit(4); }
	gfs_dt_cap = 0;
	if(split && gfs_hg_arc()){ fprintf(stderr, "[HGARC] cells=%ld modes=%ld fail=%ld dc_check=%.2e\n", hg_arc_ncell, hg_arc_nmode, hg_arc_fail, hg_arc_chk); hg_arc_ncell = hg_arc_nmode = hg_arc_fail = 0; hg_arc_chk = 0; }
	if(split) fprintf(stderr, "[HGSPLIT] dt=%.4e modes=%ld heat=%.4e skipped=%ld bulk_modes=%ld bulk_skip=%ld bulk_heat=%.4e bulk_test=%.2e bp_modes=%ld bp_skip=%ld bp_heat=%.4e bp_test=%.2e\n", dt, hgs_nmode, hgs_heat, hgm_nfail, bulk_nmode, bulk_nskip, bulk_heat, bulk_test_err, bp_nmode, bp_nskip, bp_heat, bp_test_err);
	bulk_nmode = bulk_nskip = 0; bulk_heat = 0; bp_nmode = bp_nskip = 0; bp_heat = 0;
	if(gfs_tcf_on() && tcf_sacc_file()) tcf_sacc_write((treevorostressrk4particletype*)VORORK4_TBP(simpar), VORO_NP(simpar));
	free(save); if(ssave) free(ssave);
	return dt;
}
static double exam2d_vph_rk4_int_blend_inner(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		void (*postStage)(SimParameters *)
		){
	/* All array indexing uses treevorostressrk4particletype* (sbp) for correct
	   stride when particles are stress-sized.  Common fields (x, y, vx, vy,
	   rk4, VoroQ, etc.) are at the same offsets in both struct types. */
	treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	postype xmin,ymin,xmax,ymax,OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);
	int av_mode = GAS_AVMODE(simpar);
	postype xsph_eps = GAS_XSPHEPS(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);

	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);

	int i;
	postype Dtime,dt;
	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) {
		sbp[i].w2old = sbp[i].w2;
		sbp[i].rk4.w2backup = sbp[i].w2;
	}
	/* before the stage-1 force call, which already stores the rates */
	if(gfs_tcf_on() && tcf_sacc_file()) tcf_sacc_ensure(sbp, VORO_NP(simpar));

	/* GFS_DUAL_ENERGY: give every particle an entropy on the first call.
	 * The KH/Kepler IC leaves stress.K = 0. */
	{
		static int de_k_init = 0;
		if(!de_k_init && (gfs_dual_energy_eta() > 0 || (GAS_ENTROPY_MODE(simpar) == 1 && getenv("GFS_K_INIT")))){   /* GFS_K_INIT: entropy_mode 1 from an ie-mode state */
			int nk = 0;
			for(i=0;i<VORO_NP(simpar);i++){
				if(sbp[i].stress.K > 0) continue;
				if(sbp[i].den > 0 && sbp[i].pressure > 0){
					sbp[i].stress.K = sbp[i].pressure / pow((double)sbp[i].den, (double)Gamma);
					nk++;
				}
			}
			if(MYID(simpar) == 0)
				fprintf(stderr, "[DUALE] eta=%g K set from P/rho^gamma on %d particles of rank 0\n",
						gfs_dual_energy_eta(), nk);
			de_k_init = 1;
		}
	}

	/* Reset vsig_max at start */
	if(av_mode >= 1){
		for(i=0;i<VORO_NP(simpar);i++)
			sbp[i].stress.vsig_max = 0;
	}

	/* Revision 4 defensive init: first-call guarantee that CD10 state
	 * fields (divv_old, alpha_cd) are not garbage, regardless of
	 * allocator behaviour.  my_malloc may or may not zero-init; this
	 * removes the dependency. */
	{
		static int first_call = 1;
		if(first_call && av_mode >= 1){
			postype cd_amin_init = GAS_CDAMIN(simpar);
			for(i=0;i<VORO_NP(simpar);i++){
				sbp[i].stress.divv_old = 0;
				sbp[i].stress.alpha_cd = cd_amin_init;
			}
			first_call = 0;
		}
	}

	double _t0, _t_update=0, _t_force=0, _t_post=0, _t_fin=0;

	{
		/* av_mode 5 / kappa guard, and clear any w2 the IC or a restart
		 * already set from the old kappa (CPU and GPU both read w2 here). */
		static int kg_done = 0;
		if(!kg_done){
			kg_done = 1;
			gfs_kappa_guard(simpar, "RK4 start");
			if(GAS_AVMODE(simpar) == 5 && GAS_Kappa(simpar) == 0){
				int nz = 0;
				for(i=0;i<VORO_NP(simpar);i++){
					if(sbp[i].w2 != 0 || sbp[i].w2old != 0 || sbp[i].rk4.w2backup != 0) nz++;
					sbp[i].w2 = 0; sbp[i].w2old = 0; sbp[i].rk4.w2backup = 0;
				}
				if(nz > 0)
					fprintf(stderr, "[KAPPA_GUARD] r%d: cleared nonzero w2 on %d particles (av_mode=5, kappa=0)\n", MYID(simpar), nz);
			}
		}
	}
	/* GFS_W2_LINEAR (default off): one Wdot per particle for this step. */
	int w2lin = gfs_w2lin_on() && GAS_Kappa(simpar) > 0 && !sedov_lagvol_on();
	{
		static int wl_banner = 0;
		if(!wl_banner && gfs_w2lin_on() && MYID(simpar) == 0){
			fprintf(stderr, "[W2LIN] GFS_W2_LINEAR=1: active=%d (kappa>0, no lagvol). W = W_n + c Dt Wdot in stages, exact weight term in the face speed\n", w2lin);
			wl_banner = 1;
		}
	}
	gfs_w2lin_active = w2lin;
	if(w2lin) gfs_w2lin_start(simpar, sbp, VORO_NP(simpar), targetBP, Lx, Ly);
	/* GFS_DYN_W2 (default off): scheme 3, see gfs_dynw2_on. */
	int dynw2 = gfs_dynw2_on() && av_mode == 5 && GAS_Kappa(simpar) == 0
		&& NID(simpar) == 1 && !w2lin && !sedov_lagvol_on();
#ifdef USE_CUDA
	if(GAS_GPU_ENABLED(simpar)) dynw2 = 0;   /* CPU only */
#endif
	{
		static int dw_banner = 0;
		if(!dw_banner && gfs_dynw2_on() && MYID(simpar) == 0){
			fprintf(stderr, "[DYNW2] GFS_DYN_W2=1: active=%d (needs av_mode 5, kappa 0, 1 rank, CPU, no lagvol/W2_LINEAR) ratio=%g px=%g mid=%d delta=%g barrier=%g\n",
				dynw2, gfs_env_d("GFS_DYN_W2_RATIO", 100.0), gfs_env_d("GFS_DYN_W2_PX", 100.0),
				(int)gfs_env_d("GFS_DYN_W2_MID", 0), gfs_env_d("GFS_DYN_W2_DELTA", 0.2),
				gfs_env_d("GFS_DYN_W2_BARRIER", 1000.0));
			dw_banner = 1;
		}
	}
	gfs_dynw2_active = dynw2;
	double dyn_dtlim = 1e30; int dyn_bound = 0;
	if(dynw2) dyn_step_reset();
	if(dynw2 && gfs_h_dae_on()) hq_project(simpar, xmin, ymin, xmax, ymax,
			paddingAllTreeParticles, find2DNeighborBP, find2DCellBP, mkLinkedList2D, 1);
	if(gfs_eaud_on()) gfs_eaud_start(sbp, VORO_NP(simpar), targetBP, Lx, Ly);
	/* === K1 evaluation === */
	gfs_rk_stage = 1;
	_t0 = MPI_Wtime();
	updateDenW2Pressure2DBlend(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, 1);
	_t_update += MPI_Wtime() - _t0;
	dyn_stage_cur = 0;
	if(gfs_vaud_on()) vaud_close(simpar, sbp, VORO_NP(simpar), targetBP, Lx, Ly);
	/* GFS_ENTROPY_SWITCH / GFS_HALF_LIMIT: remember ie_n and V_n at the
	 * start of the step. Stage arithmetic never touches these fields. */
	gfs_es_config();
	int gfs_es_step = (gfs_es_on || gfs_half_on) && av_mode != 4
		&& GAS_ENTROPY_MODE(simpar) != 1;
	{
		static int es_banner = 0;
		if(!es_banner && (gfs_es_on || gfs_half_on) && MYID(simpar) == 0){
			fprintf(stderr, "[ESW] entropy_switch=%d coef=%g half_limit=%d active=%d (off for av_mode 4 / entropy_mode 1)\n",
					gfs_es_on, gfs_es_coef, gfs_half_on, gfs_es_step);
			es_banner = 1;
		}
		static int lr_banner = 0;
		if(!lr_banner && gfs_laguerre_rotation_on() && MYID(simpar) == 0){
			fprintf(stderr, "[LAGROT] GFS_LAGUERRE_ROTATION=1: face-rotation correction on Laguerre faces (anchor fact1 = (1+(w2i-w2j)/d^2)/2)\n");
			lr_banner = 1;
		}
	}
	if(gfs_es_step){
		sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
		for(i=0;i<VORO_NP(simpar);i++){
			GFS_ES_IE0(sbp[i]) = sbp[i].ie;
			GFS_ES_V0(sbp[i])  = sbp[i].volume;
		}
	}
	_t0 = MPI_Wtime();
#if defined(USE_CUDA) && defined(GPU_VALIDATE)
	if (GAS_GPU_ENABLED(simpar))
		Dtime = getAccVoro2DBlend_GPU_validate(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	else
#elif defined(USE_CUDA)
	if (GAS_GPU_ENABLED(simpar))
		Dtime = getAccVoro2DBlend_GPU(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	else
#endif
	Dtime = getAccVoro2DBlend(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	_t_force += MPI_Wtime() - _t0;
	{
		postype dtk = kepler_dt_limit(simpar, VORORK4_TBP(simpar),
				TVORORK4_DDINFO(simpar)[0].n_size, VORO_NP(simpar));
		if(dtk < Dtime) Dtime = dtk;
	}
	if(dynw2) dyn_dtlim = dyn_dtlim_last;
	if(dynw2 && dyn_dtlim < Dtime){ Dtime = dyn_dtlim; dyn_bound = 1; }
	{
		/* GFS_FIXED_DT (time-convergence tests): constant step, must not exceed the CFL step */
		static double fdt = -1;
		if(fdt < 0){ const char *e_ = getenv("GFS_FIXED_DT"); fdt = (e_ && e_[0]) ? atof(e_) : 0; }
		if(fdt > 0){
			if(fdt > (double)Dtime){ fprintf(stderr, "[TCF] GFS_FIXED_DT=%g exceeds the stability limit %g: test invalid\n", fdt, (double)Dtime); exit(5); }
			Dtime = (postype)fdt;
		}
		if(gfs_dt_cap > 0 && (double)Dtime > gfs_dt_cap) Dtime = (postype)gfs_dt_cap;
	}

	clamp_dK_rate_limiter(simpar, Dtime);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].rk4.k1x  = (sbp[i].vx + xsph_eps*sbp[i].stress.v_smooth_x)*Dtime;
		sbp[i].rk4.k1y  = (sbp[i].vy + xsph_eps*sbp[i].stress.v_smooth_y)*Dtime;
		{
			/* Central mass at the stage position (EUNHA_IC=kepler only).
			 * It enters v only. phase1_ie_stage uses the hydro ax, ay. */
			postype gx = 0, gy = 0;
			kepler_accel(sbp[i].x, sbp[i].y, &gx, &gy);
			hill_coriolis(sbp[i].vx, sbp[i].vy, &gx, &gy);
			sbp[i].rk4.k1vx = (gfs_hydro_on()*sbp[i].ax + GAS_ACCX(simpar) + gx)*Dtime;
			sbp[i].rk4.k1vy = (gfs_hydro_on()*sbp[i].ay + GAS_ACCY(simpar) + gy)*Dtime;
		}
		sbp[i].rk4.k1ie = gfs_hydro_on()*phase1_ie_stage(sbp+i, Dtime);
		sbp[i].rk4.k1K  = sbp[i].stress.dK*Dtime;
		if(gfs_kepler_frozen(sbp[i].x, sbp[i].y)){
			sbp[i].rk4.k1x = sbp[i].rk4.k1y = 0; sbp[i].rk4.k1vx = sbp[i].rk4.k1vy = 0;
			sbp[i].rk4.k1ie = 0; sbp[i].rk4.k1K = 0;
		}
	}
	if(gfs_eaud_on()) gfs_eaud_stage(sbp, VORO_NP(simpar), targetBP, Lx, Ly, Dtime, 1.0/6.0);
	if(gfs_vaud_on()) vaud_stage(sbp, VORO_NP(simpar), Dtime*(1.0/6.0));
	if(gfs_tcf_on() && tcf_sacc_file()) tcf_sacc_stage(sbp, VORO_NP(simpar), (double)Dtime/6.0);

	/* === K2 evaluation === */
	gfs_rk_stage = 2;
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].x  += sbp[i].rk4.k1x*0.5;
		sbp[i].y  += sbp[i].rk4.k1y*0.5;
		sbp[i].vx += sbp[i].rk4.k1vx*0.5;
		sbp[i].vy += sbp[i].rk4.k1vy*0.5;
		sbp[i].ie += sbp[i].rk4.k1ie*0.5;
		sbp[i].stress.K += sbp[i].rk4.k1K*0.5;
	}
	_t0 = MPI_Wtime();
	postStage(simpar);
	_t_post += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;
	if(w2lin) gfs_w2lin_set(sbp, VORO_NP(simpar), (postype)0.5*Dtime);
	if(dynw2) gfs_dynw2_set(sbp, VORO_NP(simpar), (postype)0.5*Dtime, 0);
	if(dynw2 && gfs_h_dae_on()) hq_project(simpar, xmin, ymin, xmax, ymax,
			paddingAllTreeParticles, find2DNeighborBP, find2DCellBP, mkLinkedList2D, 0);
	if(av_mode >= 1)
		for(i=0;i<VORO_NP(simpar);i++) sbp[i].stress.vsig_max = 0;
	_t0 = MPI_Wtime();
	updateDenW2Pressure2DBlend(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
	_t_update += MPI_Wtime() - _t0;
	dyn_stage_cur = 1;
	_t0 = MPI_Wtime();
#if defined(USE_CUDA) && defined(GPU_VALIDATE)
	if (GAS_GPU_ENABLED(simpar))
		dt = getAccVoro2DBlend_GPU_validate(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	else
#elif defined(USE_CUDA)
	if (GAS_GPU_ENABLED(simpar))
		dt = getAccVoro2DBlend_GPU(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	else
#endif
	dt = getAccVoro2DBlend(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	_t_force += MPI_Wtime() - _t0;
	clamp_dK_rate_limiter(simpar, Dtime);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].rk4.k2x  = (sbp[i].vx + xsph_eps*sbp[i].stress.v_smooth_x)*Dtime;
		sbp[i].rk4.k2y  = (sbp[i].vy + xsph_eps*sbp[i].stress.v_smooth_y)*Dtime;
		{
			/* Central mass at the stage position (EUNHA_IC=kepler only).
			 * It enters v only. phase1_ie_stage uses the hydro ax, ay. */
			postype gx = 0, gy = 0;
			kepler_accel(sbp[i].x, sbp[i].y, &gx, &gy);
			hill_coriolis(sbp[i].vx, sbp[i].vy, &gx, &gy);
			sbp[i].rk4.k2vx = (gfs_hydro_on()*sbp[i].ax + GAS_ACCX(simpar) + gx)*Dtime;
			sbp[i].rk4.k2vy = (gfs_hydro_on()*sbp[i].ay + GAS_ACCY(simpar) + gy)*Dtime;
		}
		sbp[i].rk4.k2ie = gfs_hydro_on()*phase1_ie_stage(sbp+i, Dtime);
		sbp[i].rk4.k2K  = sbp[i].stress.dK*Dtime;
		if(gfs_kepler_frozen(sbp[i].x, sbp[i].y)){
			sbp[i].rk4.k2x = sbp[i].rk4.k2y = 0; sbp[i].rk4.k2vx = sbp[i].rk4.k2vy = 0;
			sbp[i].rk4.k2ie = 0; sbp[i].rk4.k2K = 0;
		}
	}
	if(gfs_eaud_on()) gfs_eaud_stage(sbp, VORO_NP(simpar), targetBP, Lx, Ly, Dtime, 2.0/6.0);
	if(gfs_vaud_on()) vaud_stage(sbp, VORO_NP(simpar), Dtime*(2.0/6.0));
	if(gfs_tcf_on() && tcf_sacc_file()) tcf_sacc_stage(sbp, VORO_NP(simpar), (double)Dtime/3.0);

	/* === K3 evaluation === */
	gfs_rk_stage = 3;
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].x  += (sbp[i].rk4.k2x -sbp[i].rk4.k1x )*0.5;
		sbp[i].y  += (sbp[i].rk4.k2y -sbp[i].rk4.k1y )*0.5;
		sbp[i].vx += (sbp[i].rk4.k2vx-sbp[i].rk4.k1vx)*0.5;
		sbp[i].vy += (sbp[i].rk4.k2vy-sbp[i].rk4.k1vy)*0.5;
		sbp[i].ie += (sbp[i].rk4.k2ie-sbp[i].rk4.k1ie)*0.5;
		sbp[i].stress.K += (sbp[i].rk4.k2K-sbp[i].rk4.k1K)*0.5;
	}
	_t0 = MPI_Wtime();
	postStage(simpar);
	_t_post += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;
	if(w2lin) gfs_w2lin_set(sbp, VORO_NP(simpar), (postype)0.5*Dtime);
	if(dynw2) gfs_dynw2_set(sbp, VORO_NP(simpar), (postype)0.5*Dtime, 1);
	if(dynw2 && gfs_h_dae_on()) hq_project(simpar, xmin, ymin, xmax, ymax,
			paddingAllTreeParticles, find2DNeighborBP, find2DCellBP, mkLinkedList2D, 0);
	if(av_mode >= 1)
		for(i=0;i<VORO_NP(simpar);i++) sbp[i].stress.vsig_max = 0;
	_t0 = MPI_Wtime();
	updateDenW2Pressure2DBlend(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
	_t_update += MPI_Wtime() - _t0;
	dyn_stage_cur = 2;
	_t0 = MPI_Wtime();
#if defined(USE_CUDA) && defined(GPU_VALIDATE)
	if (GAS_GPU_ENABLED(simpar))
		dt = getAccVoro2DBlend_GPU_validate(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	else
#elif defined(USE_CUDA)
	if (GAS_GPU_ENABLED(simpar))
		dt = getAccVoro2DBlend_GPU(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	else
#endif
	dt = getAccVoro2DBlend(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	_t_force += MPI_Wtime() - _t0;
	clamp_dK_rate_limiter(simpar, Dtime);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].rk4.k3x  = (sbp[i].vx + xsph_eps*sbp[i].stress.v_smooth_x)*Dtime;
		sbp[i].rk4.k3y  = (sbp[i].vy + xsph_eps*sbp[i].stress.v_smooth_y)*Dtime;
		{
			/* Central mass at the stage position (EUNHA_IC=kepler only).
			 * It enters v only. phase1_ie_stage uses the hydro ax, ay. */
			postype gx = 0, gy = 0;
			kepler_accel(sbp[i].x, sbp[i].y, &gx, &gy);
			hill_coriolis(sbp[i].vx, sbp[i].vy, &gx, &gy);
			sbp[i].rk4.k3vx = (gfs_hydro_on()*sbp[i].ax + GAS_ACCX(simpar) + gx)*Dtime;
			sbp[i].rk4.k3vy = (gfs_hydro_on()*sbp[i].ay + GAS_ACCY(simpar) + gy)*Dtime;
		}
		sbp[i].rk4.k3ie = gfs_hydro_on()*phase1_ie_stage(sbp+i, Dtime);
		sbp[i].rk4.k3K  = sbp[i].stress.dK*Dtime;
		if(gfs_kepler_frozen(sbp[i].x, sbp[i].y)){
			sbp[i].rk4.k3x = sbp[i].rk4.k3y = 0; sbp[i].rk4.k3vx = sbp[i].rk4.k3vy = 0;
			sbp[i].rk4.k3ie = 0; sbp[i].rk4.k3K = 0;
		}
	}
	if(gfs_eaud_on()) gfs_eaud_stage(sbp, VORO_NP(simpar), targetBP, Lx, Ly, Dtime, 2.0/6.0);
	if(gfs_vaud_on()) vaud_stage(sbp, VORO_NP(simpar), Dtime*(2.0/6.0));
	if(gfs_tcf_on() && tcf_sacc_file()) tcf_sacc_stage(sbp, VORO_NP(simpar), (double)Dtime/3.0);

	/* === K4 evaluation === */
	gfs_rk_stage = 4;
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].x  += sbp[i].rk4.k3x  - sbp[i].rk4.k2x*0.5;
		sbp[i].y  += sbp[i].rk4.k3y  - sbp[i].rk4.k2y*0.5;
		sbp[i].vx += sbp[i].rk4.k3vx - sbp[i].rk4.k2vx*0.5;
		sbp[i].vy += sbp[i].rk4.k3vy - sbp[i].rk4.k2vy*0.5;
		sbp[i].ie += sbp[i].rk4.k3ie - sbp[i].rk4.k2ie*0.5;
		sbp[i].stress.K += sbp[i].rk4.k3K - sbp[i].rk4.k2K*0.5;
	}
	_t0 = MPI_Wtime();
	postStage(simpar);
	_t_post += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;
	if(w2lin) gfs_w2lin_set(sbp, VORO_NP(simpar), (postype)1.0*Dtime);
	if(dynw2) gfs_dynw2_set(sbp, VORO_NP(simpar), (postype)1.0*Dtime, 2);
	if(dynw2 && gfs_h_dae_on()) hq_project(simpar, xmin, ymin, xmax, ymax,
			paddingAllTreeParticles, find2DNeighborBP, find2DCellBP, mkLinkedList2D, 0);
	if(av_mode >= 1)
		for(i=0;i<VORO_NP(simpar);i++) sbp[i].stress.vsig_max = 0;
	_t0 = MPI_Wtime();
	updateDenW2Pressure2DBlend(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
	_t_update += MPI_Wtime() - _t0;
	dyn_stage_cur = 3;
	_t0 = MPI_Wtime();
#if defined(USE_CUDA) && defined(GPU_VALIDATE)
	if (GAS_GPU_ENABLED(simpar))
		dt = getAccVoro2DBlend_GPU_validate(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	else
#elif defined(USE_CUDA)
	if (GAS_GPU_ENABLED(simpar))
		dt = getAccVoro2DBlend_GPU(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	else
#endif
	dt = getAccVoro2DBlend(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	_t_force += MPI_Wtime() - _t0;
	clamp_dK_rate_limiter(simpar, Dtime);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].rk4.k4x  = (sbp[i].vx + xsph_eps*sbp[i].stress.v_smooth_x)*Dtime;
		sbp[i].rk4.k4y  = (sbp[i].vy + xsph_eps*sbp[i].stress.v_smooth_y)*Dtime;
		{
			/* Central mass at the stage position (EUNHA_IC=kepler only).
			 * It enters v only. phase1_ie_stage uses the hydro ax, ay. */
			postype gx = 0, gy = 0;
			kepler_accel(sbp[i].x, sbp[i].y, &gx, &gy);
			hill_coriolis(sbp[i].vx, sbp[i].vy, &gx, &gy);
			sbp[i].rk4.k4vx = (gfs_hydro_on()*sbp[i].ax + GAS_ACCX(simpar) + gx)*Dtime;
			sbp[i].rk4.k4vy = (gfs_hydro_on()*sbp[i].ay + GAS_ACCY(simpar) + gy)*Dtime;
		}
		sbp[i].rk4.k4ie = gfs_hydro_on()*phase1_ie_stage(sbp+i, Dtime);
		sbp[i].rk4.k4K  = sbp[i].stress.dK*Dtime;
		if(gfs_kepler_frozen(sbp[i].x, sbp[i].y)){
			sbp[i].rk4.k4x = sbp[i].rk4.k4y = 0; sbp[i].rk4.k4vx = sbp[i].rk4.k4vy = 0;
			sbp[i].rk4.k4ie = 0; sbp[i].rk4.k4K = 0;
		}
	}
	if(gfs_eaud_on()) gfs_eaud_stage(sbp, VORO_NP(simpar), targetBP, Lx, Ly, Dtime, 1.0/6.0);
	if(gfs_vaud_on()) vaud_stage(sbp, VORO_NP(simpar), Dtime*(1.0/6.0));
	if(gfs_tcf_on() && tcf_sacc_file()) tcf_sacc_stage(sbp, VORO_NP(simpar), (double)Dtime/6.0);
	_t0 = MPI_Wtime();

	/* === Undo K4 shift and prepare for final combination === */
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].x  -= sbp[i].rk4.k3x;
		sbp[i].y  -= sbp[i].rk4.k3y;
		sbp[i].vx -= sbp[i].rk4.k3vx;
		sbp[i].vy -= sbp[i].rk4.k3vy;
		sbp[i].ie -= sbp[i].rk4.k3ie;
		sbp[i].stress.K -= sbp[i].rk4.k3K;
	}
	_t0 = MPI_Wtime();
	postStage(simpar);
	_t_post += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;

	/* === Final RK4 combination === */
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	{
		int entropy_mode = GAS_ENTROPY_MODE(simpar);
		postype K_floor = GAS_K_FLOOR(simpar);
		for(i=0;i<VORO_NP(simpar);i++){
			if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
				sbp[i].x  += (sbp[i].rk4.k1x +2*sbp[i].rk4.k2x +2*sbp[i].rk4.k3x +sbp[i].rk4.k4x )/6.;
				sbp[i].y  += (sbp[i].rk4.k1y +2*sbp[i].rk4.k2y +2*sbp[i].rk4.k3y +sbp[i].rk4.k4y )/6.;
				sbp[i].vx += (sbp[i].rk4.k1vx+2*sbp[i].rk4.k2vx+2*sbp[i].rk4.k3vx+sbp[i].rk4.k4vx)/6.;
				sbp[i].vy += (sbp[i].rk4.k1vy+2*sbp[i].rk4.k2vy+2*sbp[i].rk4.k3vy+sbp[i].rk4.k4vy)/6.;
				sbp[i].ie += (sbp[i].rk4.k1ie+2*sbp[i].rk4.k2ie+2*sbp[i].rk4.k3ie+sbp[i].rk4.k4ie)/6.;
				sbp[i].stress.K += (sbp[i].rk4.k1K+2*sbp[i].rk4.k2K+2*sbp[i].rk4.k3K+sbp[i].rk4.k4K)/6.;
			}
		}
	}
	if(gfs_tcf_on() && tcf_sacc_file() && !gfs_in_retry) tcf_sacc_write(sbp, VORO_NP(simpar));
	postStage(simpar);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	migrateTreeVorork4Particles(simpar);

	/* Arepo-Laguerre (av_mode=5): optional centroidal steering for
	   mesh regularization.  Nudges each generator toward its Voronoi
	   cell centroid by fraction GAS_FCENTROID, preventing cell
	   degeneracy in strong shear flows (KH).  Safe no-op when
	   GAS_FCENTROID == 0. */
	if(av_mode == 5 && GAS_FCENTROID(simpar) > 0){
		exam2d_centroidShift(simpar, paddingAllTreeParticles,
				find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	}

	/* GFS_W2_LINEAR: end-of-step weights before the volume update. */
	if(w2lin){
		sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
		gfs_w2lin_set(sbp, VORO_NP(simpar), Dtime);
	}
	/* GFS_DYN_W2: W_{n+1} = W_n + Dt (k1 + 2k2 + 2k3 + k4)/6. */
	if(dynw2){
		sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
		gfs_dynw2_set(sbp, VORO_NP(simpar), Dtime, -1);
		gfs_dynw2_report(simpar, Dtime, dyn_dtlim, dyn_bound);
	}
	/* Update volume & density */
	exam2dUpdateVol(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	/* Update w2 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
			if(GAS_Kappa(simpar) <0){
				sbp[i].w2 = -GAS_Kappa(simpar);
			}
			else if (GAS_Kappa(simpar) >0 && !w2lin){
				sbp[i].w2hydro = getw2forHydroParticle(simpar,(treevorork4particletype*)(sbp+i),Dtime);
				applyW2Controls(simpar, (treevorork4particletype*)(sbp+i), Dtime);
			}
		}
	}

	/* Energy bookkeeping for the disk (SEDOV_PHASE1=1 or EUNHA_IC=kepler).
	 * ie <= 0 here is what the refresh below floors to P = 1e-6. */
	double rk4e_floor_l = 0; int rk4e_neg_l = 0;
	int rk4e_on = sedov_phase1_on()
		|| (getenv("EUNHA_IC") && strcmp(getenv("EUNHA_IC"), "kepler") == 0);

	/* GFS_ENTROPY_SWITCH / GFS_HALF_LIMIT, once per step on the final
	 * state (see gfs_es_config). ie_E is the energy-equation value,
	 * including any CPU stage floor. The ledger books ie - ie_E:
	 * es_de from the switch and hl_de from the half-limit. */
	double es_de_l = 0, hl_de_l = 0; int es_n_l = 0, hl_n_l = 0;
	if(gfs_es_step){
		postype gax = GAS_ACCX(simpar), gay = GAS_ACCY(simpar);
		for(i=0;i<VORO_NP(simpar);i++){
			if(!targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)) continue;
			double V = (double)sbp[i].volume, V0 = (double)GFS_ES_V0(sbp[i]);
			double ie0 = (double)GFS_ES_IE0(sbp[i]);
			double ieE = (double)sbp[i].ie, ie1 = ieE;
			if(!(V > 0) || !(V0 > 0) || !(ie0 > 0)) continue;
			if(gfs_es_on){
				postype gx = gax, gy = gay;
				kepler_accel(sbp[i].x, sbp[i].y, &gx, &gy);
				double eg = (double)sbp[i].mass * sqrt((double)(gx*gx + gy*gy)) * sqrt(V);
				double eth = (0.5*ie0 > ieE) ? 0.5*ie0 : ieE;
				if(gfs_es_coef*eg > eth){
					ie1 = ie0 * pow(V0/V, (double)(Gamma-1));
					es_de_l += ie1 - ieE;
					es_n_l++;
				}
			}
			if(gfs_half_on && ie1 < 0.5*ie0){
				hl_de_l += 0.5*ie0 - ie1;
				hl_n_l++;
				ie1 = 0.5*ie0;
			}
			sbp[i].ie = (postype)ie1;
		}
	}

	/* GFS_DUAL_ENERGY: ie < eta*ie_K takes ie_K before the floor below.
	 * GFS_FLOOR_LOG=1 prints up to 8 such cells per rank per step (before
	 * the change), with R from the Kepler centre. Both off by default. */
	double de_inj_l = 0; int de_n_l = 0;
	{
		double de_eta = gfs_dual_energy_eta();
		static int flog = -1;
		if(flog < 0){
			const char *fl = getenv("GFS_FLOOR_LOG");
			flog = (fl && fl[0] == '1') ? 1 : 0;
		}
		if((de_eta > 0 || flog) && GAS_ENTROPY_MODE(simpar) != 1){
			postype kcx, kcy;
			int nlog = 0;
			kepler_center(&kcx, &kcy);
			for(i=0;i<VORO_NP(simpar);i++){
				if(!targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)) continue;
				if(!(sbp[i].volume > 0)) continue;
				double den_i = (double)sbp[i].mass/(double)sbp[i].volume;
				double ieK = (sbp[i].stress.K > 0)
					? (double)sbp[i].stress.K * pow(den_i, (double)Gamma)
					  * (double)sbp[i].volume / (double)(Gamma-1)
					: 0.0;
				int reset = (de_eta > 0 && ieK > 0 && (double)sbp[i].ie < de_eta*ieK);
				if(flog && (reset || sbp[i].ie <= 0) && nlog < 8){
					double dx = (double)(sbp[i].x - kcx), dy = (double)(sbp[i].y - kcy);
					fprintf(stderr, "[RK4F] r%d x=%.5f y=%.5f R=%.4f den=%.4e ie=%.4e ie_K=%.4e ie/ie_K=%.3e vol=%.4e reset=%d\n",
							MYID(simpar), (double)sbp[i].x, (double)sbp[i].y, sqrt(dx*dx+dy*dy),
							den_i, (double)sbp[i].ie, ieK,
							(ieK > 0 ? (double)sbp[i].ie/ieK : 0.0), (double)sbp[i].volume, reset);
					nlog++;
				}
				if(reset){
					de_inj_l += ieK - (double)sbp[i].ie;
					de_n_l++;
					sbp[i].ie = (postype)ieK;
				}
			}
		}
	}

	if(rk4e_on && GAS_ENTROPY_MODE(simpar) != 1){
		for(i=0;i<VORO_NP(simpar);i++){
			if(!targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)) continue;
			if(sbp[i].ie <= 0){
				rk4e_neg_l++;
				rk4e_floor_l += 1e-6*(double)sbp[i].volume/(double)(Gamma-1) - (double)sbp[i].ie;
			}
		}
	}

	/* Update pressure, csound, NaN check */
	{
		int entropy_mode = GAS_ENTROPY_MODE(simpar);
		postype K_floor = GAS_K_FLOOR(simpar);
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(i=0;i<VORO_NP(simpar);i++){
			if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
				sbp[i].den = sbp[i].mass/sbp[i].volume;
				if(entropy_mode == 1){
					if(sbp[i].stress.K < K_floor) sbp[i].stress.K = K_floor;
					sbp[i].pressure = sbp[i].stress.K * pow((double)sbp[i].den, (double)Gamma);
					sbp[i].ie = sbp[i].pressure * sbp[i].volume / (Gamma-1);
				} else {
					sbp[i].pressure = sbp[i].ie/(sbp[i].volume*sub_gfac((treevorork4particletype*)(sbp+i))) * (Gamma-1);   /* g of the last stage geometry (diagnostic) */
					if(sbp[i].pressure <= 0){
						gfs_attempt_failed = 1;
						sbp[i].pressure = 1e-6;
						sbp[i].ie = sbp[i].pressure * sbp[i].volume*sub_gfac((treevorork4particletype*)(sbp+i)) / (Gamma-1);
					}
				}
				sbp[i].csound = sqrt(Gamma*sbp[i].pressure/sbp[i].den);
			}
			else {
				sbp[i].vx = sbp[i].vy = 0;
				sbp[i].w2 = sbp[i].rk4.w2backup;
			}
		}
	}

	if(rk4e_on){
		double eh_l = 0, ep_l = 0, eh_g, ep_g, fl_g;
		double ei_l = 0, ei_g = 0, dei_g = 0;
		int ng_g, dn_g = 0;
		for(i=0;i<VORO_NP(simpar);i++){
			if(!targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)) continue;
			double m = (double)sbp[i].mass;
			eh_l += (double)sbp[i].ie + 0.5*m*((double)sbp[i].vx*(double)sbp[i].vx
					+ (double)sbp[i].vy*(double)sbp[i].vy);
			ei_l += (double)sbp[i].ie;
			ep_l += m*(double)kepler_phi(sbp[i].x, sbp[i].y);
		}
		MPI_Allreduce(&eh_l, &eh_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
		MPI_Allreduce(&ep_l, &ep_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
		MPI_Allreduce(&ei_l, &ei_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
		MPI_Allreduce(&de_inj_l, &dei_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
		MPI_Allreduce(&de_n_l, &dn_g, 1, MPI_INT, MPI_SUM, MPI_COMM(simpar));
		MPI_Allreduce(&rk4e_floor_l, &fl_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
		MPI_Allreduce(&rk4e_neg_l, &ng_g, 1, MPI_INT, MPI_SUM, MPI_COMM(simpar));
		double sfl_g = 0; int sfl_n_g = 0;
		MPI_Allreduce(&gfs_stage_floor_l, &sfl_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
		MPI_Allreduce(&gfs_stage_floor_n_l, &sfl_n_g, 1, MPI_INT, MPI_SUM, MPI_COMM(simpar));
		gfs_stage_floor_l = 0; gfs_stage_floor_n_l = 0;
		double es_de_g = 0, hl_de_g = 0; int es_n_g = 0, hl_n_g = 0;
		if(gfs_es_step){
			MPI_Allreduce(&es_de_l, &es_de_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
			MPI_Allreduce(&hl_de_l, &hl_de_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
			MPI_Allreduce(&es_n_l, &es_n_g, 1, MPI_INT, MPI_SUM, MPI_COMM(simpar));
			MPI_Allreduce(&hl_n_l, &hl_n_g, 1, MPI_INT, MPI_SUM, MPI_COMM(simpar));
		}
		long long npair_l = gfs_pair_face_ends, npair_g = 0;
		MPI_Allreduce(&npair_l, &npair_g, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM(simpar));
		/* [RK4E] energies above are global; audit stage sums must have the
		 * same scope before they are compared. Keep the reduction diagnostic
		 * only so ordinary runs incur no additional MPI collective. */
		double eaud_g[5] = {0, 0, 0, 0, 0};
		if(gfs_eaud_on()){
			double eaud_l[5] = {gfs_eaud_hyd, gfs_eaud_hyd_abs,
				gfs_eaud_wg, gfs_eaud_eh0, gfs_eaud_ep0};
			MPI_Allreduce(eaud_l, eaud_g, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
		}
		if(gfs_tq_on()){
			double acc_g[TQ_NB*3], cap_g[2*TQ_NB], loc[2*TQ_NB]; long nc_g = 0;
			int kk, jj;
			MPI_Allreduce(&tq_acc[0][0], acc_g, TQ_NB*3, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
			for(kk=0;kk<TQ_NB;kk++){ loc[kk] = tq_cap_n[kk]; loc[TQ_NB+kk] = tq_cap_len[kk]; }
			MPI_Allreduce(loc, cap_g, 2*TQ_NB, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
			MPI_Allreduce(&tq_ncalls, &nc_g, 1, MPI_LONG, MPI_SUM, MPI_COMM(simpar));
			long ncr = nc_g / NID(simpar); if(ncr < 1) ncr = 1;
			for(kk=0;kk<TQ_NB;kk++) for(jj=0;jj<3;jj++) tq_cum[kk][jj] += acc_g[kk*3+jj]/(double)ncr*(double)Dtime;
			if(MYID(simpar)==0){
				fprintf(stderr, "[TORQ] dt=%.3e calls=%ld cumLz(bin: HLLC,extHLL,pair)", (double)Dtime, ncr);
				for(kk=0;kk<TQ_NB;kk++) fprintf(stderr, " b%d:%.3e,%.3e,%.3e", kk, tq_cum[kk][0], tq_cum[kk][1], tq_cum[kk][2]);
				fprintf(stderr, " capfaces(n,len per call)");
				for(kk=0;kk<TQ_NB;kk++) fprintf(stderr, " b%d:%.1f,%.3f", kk, cap_g[kk]/(double)ncr, cap_g[TQ_NB+kk]/(double)ncr);
				fprintf(stderr, "\n");
			}
			{
				double lj[7*TQ_NB], lj_g[7*TQ_NB];
				for(kk=0;kk<TQ_NB;kk++){ lj[kk]=lq_n[kk]; lj[TQ_NB+kk]=lq_asum[kk]; lj[2*TQ_NB+kk]=lq_low[kk];
					lj[3*TQ_NB+kk]=jq_n[kk]; lj[4*TQ_NB+kk]=jq_cell[kk]; lj[5*TQ_NB+kk]=jq_rec[kk]; lj[6*TQ_NB+kk]=jq_keep[kk]; }
				MPI_Allreduce(lj, lj_g, 7*TQ_NB, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
				if(MYID(simpar)==0){
					fprintf(stderr, "[LIMQ] bin:<alpha_v>,frac(alpha<0.5),<|du_n|cell>,<|du_n|rec>,frac(rec>0.5cell)");
					for(kk=0;kk<TQ_NB;kk++){
						double nl = lj_g[kk] > 0 ? lj_g[kk] : 1, nf = lj_g[3*TQ_NB+kk] > 0 ? lj_g[3*TQ_NB+kk] : 1;
						fprintf(stderr, " b%d:%.3f,%.3f,%.3e,%.3e,%.3f", kk, lj_g[TQ_NB+kk]/nl, lj_g[2*TQ_NB+kk]/nl,
							lj_g[4*TQ_NB+kk]/nf, lj_g[5*TQ_NB+kk]/nf, lj_g[6*TQ_NB+kk]/nf);
					}
					fprintf(stderr, "\n");
				}
				for(kk=0;kk<TQ_NB;kk++){ lq_n[kk]=lq_asum[kk]=lq_low[kk]=0; jq_n[kk]=jq_cell[kk]=jq_rec[kk]=jq_keep[kk]=0; }
			}
			{
				double a4[TQ_NB*4], n4[TQ_NB*4], p4[TQ_NB*4], u4[TQ_NB*4];
				MPI_Allreduce(&tqs_acc[0][0], a4, TQ_NB*4, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
				MPI_Allreduce(&tqs_n[0][0], n4, TQ_NB*4, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
				MPI_Allreduce(&tqs_pp[0][0], p4, TQ_NB*4, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
				MPI_Allreduce(&tqs_du[0][0], u4, TQ_NB*4, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
				for(kk=0;kk<TQ_NB;kk++) for(jj=0;jj<4;jj++) tqs_cum[kk][jj] += a4[kk*4+jj]/(double)ncr*(double)Dtime;
				if(MYID(simpar)==0){
					fprintf(stderr, "[TQSIGN] cumLz(exp-trail,exp-lead,comp-trail,comp-lead) | mean P*/P | mean |du_n|/c");
					for(kk=0;kk<5;kk++){
						fprintf(stderr, " b%d:", kk);
						for(jj=0;jj<4;jj++) fprintf(stderr, "%s%.3e", jj ? "," : "", tqs_cum[kk][jj]);
						fprintf(stderr, "|");
						for(jj=0;jj<4;jj++) fprintf(stderr, "%s%.2e", jj ? "," : "", n4[kk*4+jj] > 0 ? p4[kk*4+jj]/n4[kk*4+jj] : 0.0);
						fprintf(stderr, "|");
						for(jj=0;jj<4;jj++) fprintf(stderr, "%s%.2e", jj ? "," : "", n4[kk*4+jj] > 0 ? u4[kk*4+jj]/n4[kk*4+jj] : 0.0);
					}
					fprintf(stderr, "\n");
				}
				memset(tqs_acc, 0, sizeof(tqs_acc)); memset(tqs_n, 0, sizeof(tqs_n));
				memset(tqs_pp, 0, sizeof(tqs_pp)); memset(tqs_du, 0, sizeof(tqs_du));
			}
			{
				double sg[TQ_NB];
				MPI_Allreduce(sfl_bin, sg, TQ_NB, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
				if(MYID(simpar)==0){
					fprintf(stderr, "[SFLB] stage floor resets per bin (edges 0.5,0.55,0.6,0.75,1,1.5,1.9,2.1)");
					for(kk=0;kk<TQ_NB;kk++) fprintf(stderr, " b%d:%.0f", kk, sg[kk]);
					fprintf(stderr, "\n");
				}
				for(kk=0;kk<TQ_NB;kk++) sfl_bin[kk] = 0;
			}
			for(kk=0;kk<TQ_NB;kk++){ for(jj=0;jj<3;jj++) tq_acc[kk][jj] = 0; tq_cap_n[kk] = 0; tq_cap_len[kk] = 0; }
			tq_ncalls = 0;
		}
		if(MYID(simpar)==0){
			static double Et0 = 0, Eh0 = 0, fl_cum = 0, de_cum = 0;
			static double es_cum = 0, hl_cum = 0, sfl_cum = 0;
			static int st = 0, first = 1;
			st++;
			fl_cum += fl_g;
			de_cum += dei_g;
			double et = eh_g + ep_g;
			if(first){ Et0 = et; Eh0 = eh_g; first = 0; }
			fprintf(stderr,
				"[RK4E] step=%d dt=%.3e E_hyd=%.6e E_pot=%.6e E_tot=%.6e dEtot/|Etot0|=%.3e dEhyd/Ehyd0=%.3e n_ie_le0=%d floor_inj=%.3e floor_cum=%.3e n_pair=%lld E_int=%.6e n_de=%d de_inj=%.3e de_cum=%.3e",
				st, (double)Dtime, eh_g, ep_g, et,
				(Et0 != 0 ? (et-Et0)/fabs(Et0) : 0), (Eh0 != 0 ? (eh_g-Eh0)/Eh0 : 0),
				ng_g, fl_g, fl_cum, npair_g, ei_g, dn_g, dei_g, de_cum);
			fprintf(stderr, " cellcheck=%ld/%ld", gfs_cellcheck_fail, gfs_cellcheck_n);
			if(gfs_topo_on()){ fprintf(stderr, " topo_chg=%ld/%ld", gfs_topo_chg, gfs_topo_n); gfs_topo_chg = gfs_topo_n = 0; }
			if(gfs_lsf_order() == 2){
				fprintf(stderr, " q2fail=%ld/%ld q2smooth=%ld,%ld,%ld,%ld", gfs_q2_fail, gfs_q2_n, gfs_q2_smooth[0], gfs_q2_smooth[1], gfs_q2_smooth[2], gfs_q2_smooth[3]);
				gfs_q2_fail = gfs_q2_n = 0; gfs_q2_smooth[0] = gfs_q2_smooth[1] = gfs_q2_smooth[2] = gfs_q2_smooth[3] = 0;
			}
			gfs_cellcheck_fail = gfs_cellcheck_n = 0;
			fprintf(stderr, " recfloor=%ld rf_faces=%ld negP=%ld tinyP=%ld negR=%ld tinyR=%ld rf_arccell=%ld APs_floor/all=%.3e",
				gfs_recfloor_n, gfs_rf_faces, gfs_rf_negP, gfs_rf_tinyP, gfs_rf_negR, gfs_rf_tinyR, gfs_rf_arc,
				gfs_rf_APs_all > 0 ? gfs_rf_APs_floor/gfs_rf_APs_all : 0.0);
			fprintf(stderr, " negmagP=%ld,%ld,%ld,%ld,%ld negmagR=%ld,%ld,%ld,%ld,%ld",
				gfs_rf_mag[0][0], gfs_rf_mag[0][1], gfs_rf_mag[0][2], gfs_rf_mag[0][3], gfs_rf_mag[0][4],
				gfs_rf_mag[1][0], gfs_rf_mag[1][1], gfs_rf_mag[1][2], gfs_rf_mag[1][3], gfs_rf_mag[1][4]);
			memset(gfs_rf_mag, 0, sizeof(gfs_rf_mag));
			if(gfs_hg_affine()){
				fprintf(stderr, " [HG] H=%.4e vF=%.4e sum=%.2e cells=%ld skip=%ld rankdef=%ld dtlim=%.3e worst: lamf=%.3e d/sqrtV=%.3f 1-u=%.3f |gM|=%.3e |lM|=%.3e cond=%.2e |r|h/c=%.2e ns=%.0f id=%.0f C/acou=%.3f",
					hg_H_last, hg_vf_last, hg_H_last + hg_vf_last, hg_ncell, hg_nskip, hg_npinv, hg_dtlim,
					hg_worst[0], hg_worst[1], hg_worst[2], hg_worst[3], hg_worst[4], hg_worst[5], hg_worst[6], hg_worst[7], hg_worst[8], hg_worst[9]);
				{ double m_ = gfs_dtb[0]; int c_ = 0, z_; double hgl = hg_dtlim; for(z_ = 1; z_ < 3; z_++) if(gfs_dtb[z_] < m_){ m_ = gfs_dtb[z_]; c_ = z_; } if(hgl < m_){ m_ = hgl; c_ = 3; }
				  gfs_dtb_cnt[c_]++;
				  fprintf(stderr, " dtlims: cfl=%.3e geom=%.3e dt3=%.3e hg=%.3e bind(cfl,geom,dt3,hg)=%ld,%ld,%ld,%ld", gfs_dtb[0], gfs_dtb[1], gfs_dtb[2], hgl, gfs_dtb_cnt[0], gfs_dtb_cnt[1], gfs_dtb_cnt[2], gfs_dtb_cnt[3]);
				  for(z_ = 0; z_ < 3; z_++) gfs_dtb[z_] = 1e300; }
				if(gfs_hg_spec()) fprintf(stderr, " spec: lam_power=%.4e res=%.2e bound=%.4e bound/lam=%.2f", hg_spec_lam, hg_spec_res, hg_spec_bound, hg_spec_lam > 0 ? hg_spec_bound/hg_spec_lam : 0);
				{ int z_; for(z_ = 0; z_ < 12; z_++) hg_worst[z_] = 0; }
				hg_ncell = hg_nskip = hg_npinv = 0;
				if(gfs_hg_quad()){ fprintf(stderr, " [HGQUAD] qrankdef=%ld qtest=%.2e", hg_nqdef, hg_qtest); hg_nqdef = 0; hg_qtest = 0; }
			}
			if(gfs_periodic_complete()){
				double Lx_ = (double)(SIMBOX(simpar).x.max - SIMBOX(simpar).x.min), Ly_ = (double)(SIMBOX(simpar).y.max - SIMBOX(simpar).y.min);
				fprintf(stderr, " [PCOMP] cap_cells=%ld vsum/box-1=%.3e", pc_ncap, pc_vsum/(Lx_*Ly_) - 1);
				if(pc_ncap > 0){ fprintf(stderr, "\n[PCOMP] active cap edges in a periodic run: gate failure\n"); exit(2); }
			}
			if(gfs_shared_geom_on()){
				fprintf(stderr, " [SHGEOM] look=%ld cut=%ld miss=%ld maxcut(A>1e-6d)=%.2e fzclamp=%ld", shf_nlook, shf_ncut, shf_nmiss, shf_maxcut, fz_nclamp);
				shf_nlook = shf_ncut = shf_nmiss = 0; shf_maxcut = 0; fz_nclamp = 0;
			}
			if(gfs_subzone_on()){
				double sphi = 0, sdp = 0, svf = sub_vf_last; size_t k_;
				for(k_ = 0; k_ < sub_cap; k_++){ sphi += sub_phi[k_]; sdp += sub_dphi[k_]; }
				fprintf(stderr, " [SUBZ] Phi=%.6e dPhi=%.4e vF=%.4e sum=%.2e cells=%ld cap=%ld", sphi, sdp, svf, sdp + svf, sub_ncell, sub_ncap);
				sub_ncell = sub_ncap = 0;
			}
			if(gfs_tcf_on()){
				fprintf(stderr, " [TCF] faces=%ld heat=%.4e zneg=%ld zmin=%.3e mm0=%ld small=%ld cellneg=%ld tanwork=%.4e",
					tcf_faces, tcf_heat, tcf_zneg, tcf_zmin, tcf_mm0, tcf_small, tcf_cellneg, tcf_tanwork);
				tcf_tanwork = 0;
				tcf_faces = tcf_zneg = tcf_mm0 = tcf_small = tcf_cellneg = 0; tcf_heat = tcf_zmin = 0;
				if(gfs_tcf_qbp()){
					long ns_[3] = {0, 0, 0}; size_t k_; for(k_ = 0; k_ < qbp_cap; k_++) if(qbp_ok[k_] >= 0 && qbp_ok[k_] <= 2) ns_[(int)qbp_ok[k_]]++;
					if(gfs_q3()){ fprintf(stderr, " [Q3] cubic=%ld fallback=%ld", q3_n, q3_fb); q3_n = q3_fb = 0; }
					fprintf(stderr, " [QBP] faces=%ld invalidface=%ld cells(q2,lin,none)=%ld,%ld,%ld q2fail(npt<8|widened,pivot|norank)=%ld,%ld", qbp_nface, qbp_nfail, ns_[2], ns_[1], ns_[0], qbp_q2fail_npt, qbp_q2fail_piv);
					qbp_nface = qbp_nfail = 0; qbp_q2fail_npt = qbp_q2fail_piv = 0;
					if(gfs_qbp_edge_inner()){ long na_ = 0; for(k_ = 0; k_ < qbp_cap; k_++) na_ += qbp_arc[k_]; fprintf(stderr, " [QBPEDGE] arccells=%ld mixedfaces=%ld", na_, qbp_nmixed); qbp_nmixed = 0; }
				}
			}
			gfs_recfloor_n = 0; gfs_rf_negP = gfs_rf_tinyP = gfs_rf_negR = gfs_rf_tinyR = gfs_rf_arc = gfs_rf_faces = 0;
			gfs_rf_APs_floor = gfs_rf_APs_all = 0;
			/* Extra columns at the end of [RK4E], only when the switch or
			 * the half-limit is on. The ledger becomes
			 * dEtot*|Etot0| ~ floor_cum + de_cum + es_cum + hl_cum. */
			if(gfs_es_step){
				es_cum += es_de_g;
				hl_cum += hl_de_g;
				fprintf(stderr,
					" n_es=%d es_de=%.3e es_cum=%.3e n_half=%d hl_de=%.3e hl_cum=%.3e",
					es_n_g, es_de_g, es_cum, hl_n_g, hl_de_g, hl_cum);
			}
			/* CPU stage floor. The GPU path floors only P in the stages, so
			 * a negative ie reaches the end-of-step floor and is counted in
			 * floor_cum there. Ledger:
			 * dEtot*|Etot0| ~ floor_cum + sfl_cum + de_cum (+ es_cum, hl_cum). */
			sfl_cum += sfl_g;
			fprintf(stderr, " n_sfl=%d sfl_inj=%.3e sfl_cum=%.3e", sfl_n_g, sfl_g, sfl_cum);
			if(gfs_eaud_on()){
				static double rkh_cum = 0, rkp_cum = 0, hyd_cum = 0;
				double dEh = eh_g - eaud_g[3], dEp = ep_g - eaud_g[4];
				double rkh = dEh - eaud_g[0] - eaud_g[2] - fl_g - sfl_g;
				double rkp = dEp + eaud_g[2];
				rkh_cum += rkh; rkp_cum += rkp; hyd_cum += eaud_g[0];
				fprintf(stderr, " aud_hyd=%.3e aud_habs=%.3e aud_rkh=%.3e aud_rkp=%.3e hyd_cum=%.3e rkh_cum=%.3e rkp_cum=%.3e",
					eaud_g[0], eaud_g[1], rkh, rkp, hyd_cum, rkh_cum, rkp_cum);
			}
			fprintf(stderr, "\n");
			fflush(stderr);
		}
#ifdef USE_CUDA
		/* Trigger all four force stages of the next step after three
		 * consecutive large imbalances. This rejects isolated spikes.
		 * All MPI ranks see the same globally reduced eaud_g. */
		if(GAS_GPU_ENABLED(simpar) && gfs_eaud_on()){
			static double dump_hyd_threshold = -1;
			static int dump_requested = 0;
			static int consecutive = 0;
			if(dump_hyd_threshold < 0){
				const char *s = getenv("GFS_GPU_FACE_DUMP_HYD_THRESHOLD");
				dump_hyd_threshold = (s && s[0]) ? atof(s) : 0;
				if(!isfinite(dump_hyd_threshold) || dump_hyd_threshold <= 0)
					dump_hyd_threshold = 0;
			}
			if(dump_hyd_threshold > 0 && fabs(eaud_g[0]) >= dump_hyd_threshold)
				consecutive++;
			else
				consecutive = 0;
			if(!dump_requested && consecutive >= 3){
				gpu_request_face_dumps(4);
				dump_requested = 1;
				if(MYID(simpar)==0){
					fprintf(stderr,
						"[FACE_DUMP_TRIGGER] aud_hyd=%.6e threshold=%.6e next_force_calls=4\n",
						eaud_g[0], dump_hyd_threshold);
					fflush(stderr);
				}
			}
		}
#endif
	}

	/* CD10 alpha update after full RK4 step */
	if(av_mode >= 1){
		update_alpha_cd_2d(simpar, Dtime);
	}

	gfs_tnow += (double)Dtime;
	gfs_w2lin_active = 0;
	gfs_dynw2_active = 0;
	_t_fin = MPI_Wtime() - _t0;
	{
		static int _step_timer = 0;
		_step_timer++;
		double _t_total = _t_update + _t_force + _t_post + _t_fin;
		if(MYID(simpar) == 0)
			printf("[CPU step %d] update=%.1fs force=%.1fs post=%.1fs final=%.1fs | total=%.1fs\n",
				_step_timer, _t_update, _t_force, _t_post, _t_fin, _t_total);
		fflush(stdout);
	}

	return Dtime;
}

/* ================================================================
 *  exam2d_vph_kdk_int_blend (EVOLMETHOD=3):
 *    KDK leapfrog integrator for the *blend* (Voronoi) hydro path.
 *    2 force evals per step vs RK4's 4 → ~2x faster, 2nd-order time
 *    accuracy (matches 2nd-order spatial reconstruction).
 *
 *    CFL: Athena++-style FIXED Δx_uniform = Lx/Nx.
 *      Dtime = Courant * Δx_uniform / max_global(stress.vsig_max)
 *
 *    Flow (mirrors exam2d_vph_rk4_int_lagmfm but calls blend routines):
 *      1) update density/pressure + force eval at (x_n, v_n)
 *         → (a_n, du_n) ; compute Dtime from fixed-Δx CFL
 *      2) half-kick: v += 0.5·Dtime·a_n;  ie += 0.5·Dtime·du_n
 *      3) drift:     x += Dtime·v
 *      4) postStage + migrate handled below
 *      5) update density/pressure + force at (x_{n+1}, v_{n+1/2})
 *         → (a_{n+1}, du_{n+1})
 *      6) half-kick: v += 0.5·Dtime·a_{n+1}; ie += 0.5·Dtime·du_{n+1}
 *
 *    Same Dtime is used throughout the step (consistent leapfrog).
 * ================================================================ */
double exam2d_vph_kdk_int_blend(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		void (*postStage)(SimParameters *)
		){
	treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	postype xmin,ymin,xmax,ymax,OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);
	int av_mode = GAS_AVMODE(simpar);
	postype xsph_eps = GAS_XSPHEPS(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);

	int i;
	postype Dtime, dt_dummy;
	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	/* Athena++-style fixed Δx_uniform = Lx/Nx (initial grid spacing) */
	postype dx_uniform = Lx / NX(simpar);
	postype accx_ext = GAS_ACCX(simpar), accy_ext = GAS_ACCY(simpar);

	if(MYID(simpar)==0){ fprintf(stderr,"[KDK_DBG] H entry np=%td sbp=%p\n", VORO_NP(simpar), (void*)sbp); fflush(stderr); }
	for(i=0;i<VORO_NP(simpar);i++) {
		sbp[i].w2old = sbp[i].w2;
		sbp[i].rk4.w2backup = sbp[i].w2;
	}
	if(MYID(simpar)==0){ fprintf(stderr,"[KDK_DBG] I after w2backup loop\n"); fflush(stderr); }

	if(av_mode >= 1)
		for(i=0;i<VORO_NP(simpar);i++)
			sbp[i].stress.vsig_max = 0;
	if(MYID(simpar)==0){ fprintf(stderr,"[KDK_DBG] J after vsig_max reset\n"); fflush(stderr); }

	/* CD10 first-call defensive init (mirror of RK4 path) */
	{
		static int first_call = 1;
		if(first_call && av_mode >= 1){
			postype cd_amin_init = GAS_CDAMIN(simpar);
			for(i=0;i<VORO_NP(simpar);i++){
				sbp[i].stress.divv_old = 0;
				sbp[i].stress.alpha_cd = cd_amin_init;
			}
			first_call = 0;
		}
	}
	if(MYID(simpar)==0){ fprintf(stderr,"[KDK_DBG] K before updateDen\n"); fflush(stderr); }

	double _t0, _t_update=0, _t_force=0, _t_post=0, _t_fin=0;

	/* === Stage 1: force eval at (x_n, v_n) === */
	_t0 = MPI_Wtime();
	updateDenW2Pressure2DBlend(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, 1);
	if(MYID(simpar)==0){ fprintf(stderr,"[KDK_DBG] L after stage1 updateDen\n"); fflush(stderr); }
	_t_update += MPI_Wtime() - _t0;
	_t0 = MPI_Wtime();
#ifdef USE_CUDA
	if (GAS_GPU_ENABLED(simpar))
		dt_dummy = getAccVoro2DBlend_GPU(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	else
#endif
	dt_dummy = getAccVoro2DBlend(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	(void)dt_dummy;
	if(MYID(simpar)==0){ fprintf(stderr,"[KDK_DBG] M after stage1 force\n"); fflush(stderr); }
	_t_force += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	/* Athena++-style CFL: Dtime = Courant * Δx_uniform / max_global(vsig_max) */
	{
		postype vsig_max_local = 1e-30;
		for(i=0;i<VORO_NP(simpar);i++){
			postype v = (av_mode >= 1) ? sbp[i].stress.vsig_max : 0;
			/* Fallback for av_mode==0: vsig ~ 2*csound */
			if(v <= 0) v = 2.0 * sbp[i].csound;
			if(v > vsig_max_local) vsig_max_local = v;
		}
		postype vsig_max_global;
		MPI_Allreduce(&vsig_max_local, &vsig_max_global, 1, MPI_POSTYPE, MPI_MAX, MPI_COMM(simpar));
		Dtime = Courant * dx_uniform / vsig_max_global;
		/* Face dt uses the true generator separation. The uniform-dx
		 * step above stays larger than a gap that has already narrowed,
		 * so a converging pair crosses in one drift. */
		if(dt_dummy > (postype)0 && dt_dummy < Dtime) Dtime = dt_dummy;
		{
			postype dtk = kepler_dt_limit(simpar, VORORK4_TBP(simpar),
					TVORORK4_DDINFO(simpar)[0].n_size, VORO_NP(simpar));
			if(dtk < Dtime) Dtime = dtk;
		}
	}

	/* dK rate-limiter: clamp |dK|*dt_ref/K <= rate_max before stage-1 half-kick */
	clamp_dK_rate_limiter(simpar, Dtime);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	/* Half-kick with a_n, du_n */
	postype half_dt = 0.5 * Dtime;
	for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
			postype axk = accx_ext, ayk = accy_ext;
			kepler_accel(sbp[i].x, sbp[i].y, &axk, &ayk);
			if(sedov_phase1_on()){
				phase1_half_kick(sbp+i, half_dt, axk, ayk);
			} else {
				sbp[i].vx += (sbp[i].ax + axk) * half_dt;
				sbp[i].vy += (sbp[i].ay + ayk) * half_dt;
				sbp[i].ie += sbp[i].die             * half_dt;
				sbp[i].stress.K += sbp[i].stress.dK * half_dt;
			}
		}
	}

	/* Full drift with half-kicked v (XSPH offset added when xsph_eps>0) */
	for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
			sbp[i].x += (sbp[i].vx + xsph_eps*sbp[i].stress.v_smooth_x) * Dtime;
			sbp[i].y += (sbp[i].vy + xsph_eps*sbp[i].stress.v_smooth_y) * Dtime;
		}
	}

	_t0 = MPI_Wtime();
	postStage(simpar);
	if(MYID(simpar)==0){ fprintf(stderr,"[KDK_DBG] N after stage1 drift+postStage\n"); fflush(stderr); }
	_t_post += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;
	if(av_mode >= 1)
		for(i=0;i<VORO_NP(simpar);i++) sbp[i].stress.vsig_max = 0;

	/* === Stage 2: force eval at drifted state, then second half-kick === */
	_t0 = MPI_Wtime();
	updateDenW2Pressure2DBlend(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
	if(MYID(simpar)==0){ fprintf(stderr,"[KDK_DBG] O after stage2 updateDen\n"); fflush(stderr); }
	_t_update += MPI_Wtime() - _t0;
	_t0 = MPI_Wtime();
#ifdef USE_CUDA
	if (GAS_GPU_ENABLED(simpar))
		dt_dummy = getAccVoro2DBlend_GPU(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	else
#endif
	dt_dummy = getAccVoro2DBlend(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	(void)dt_dummy;
	if(MYID(simpar)==0){ fprintf(stderr,"[KDK_DBG] P after stage2 force\n"); fflush(stderr); }
	_t_force += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	/* dK rate-limiter: clamp |dK|*dt_ref/K <= rate_max before stage-2 half-kick */
	clamp_dK_rate_limiter(simpar, Dtime);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	/* Second half-kick with a_{n+1}, du_{n+1} */
	{
		for(i=0;i<VORO_NP(simpar);i++){
			if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
				postype axk = accx_ext, ayk = accy_ext;
				kepler_accel(sbp[i].x, sbp[i].y, &axk, &ayk);
				if(sedov_phase1_on()){
					phase1_half_kick(sbp+i, half_dt, axk, ayk);
				} else {
					sbp[i].vx += (sbp[i].ax + axk) * half_dt;
					sbp[i].vy += (sbp[i].ay + ayk) * half_dt;
					sbp[i].ie += sbp[i].die             * half_dt;
					sbp[i].stress.K += sbp[i].stress.dK * half_dt;
				}
			} else {
				sbp[i].vx = sbp[i].vy = 0;
			}
		}
	}

	_t0 = MPI_Wtime();
	postStage(simpar);
	if(MYID(simpar)==0){ fprintf(stderr,"[KDK_DBG] Q after 2nd postStage\n"); fflush(stderr); }
	_t_post += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	migrateTreeVorork4Particles(simpar);
	if(MYID(simpar)==0){ fprintf(stderr,"[KDK_DBG] R after migrate\n"); fflush(stderr); }

	/* Centroidal steering (av_mode=5 mesh regularization) */
	if(av_mode == 5 && GAS_FCENTROID(simpar) > 0){
		exam2d_centroidShift(simpar, paddingAllTreeParticles,
				find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	}

	/* Refresh volume & density after migration/centroidShift */
	exam2dUpdateVol(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);
	if(MYID(simpar)==0){ fprintf(stderr,"[KDK_DBG] S after end-of-step exam2dUpdateVol\n"); fflush(stderr); }

	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	/* Update w2 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
			if(GAS_Kappa(simpar) < 0){
				sbp[i].w2 = -GAS_Kappa(simpar);
			}
			else if(GAS_Kappa(simpar) > 0){
				sbp[i].w2hydro = getw2forHydroParticle(simpar,(treevorork4particletype*)(sbp+i),Dtime);
				applyW2Controls(simpar, (treevorork4particletype*)(sbp+i), Dtime);
			}
		}
	}

	if(sedov_phase1_on()){
		double vmin_l = 1e300, rhomax_l = 0, etot_l = 0;
		int nneg_l = 0, ncell_l = 0;
		for(i=0;i<VORO_NP(simpar);i++){
			if(!targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)) continue;
			double V = (double)sbp[i].volume;
			double m = (double)sbp[i].mass;
			double rho = (V > 0) ? m/V : 0;
			double ke = 0.5*m*((double)sbp[i].vx*(double)sbp[i].vx
					+ (double)sbp[i].vy*(double)sbp[i].vy);
			if(V < vmin_l) vmin_l = V;
			if(rho > rhomax_l) rhomax_l = rho;
			etot_l += (double)sbp[i].ie + ke;
			if(sbp[i].ie <= 0) nneg_l++;
			ncell_l++;
		}
		double vmin_g, rhomax_g, etot_g;
		int nneg_g, ncell_g;
		MPI_Allreduce(&vmin_l, &vmin_g, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM(simpar));
		MPI_Allreduce(&rhomax_l, &rhomax_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM(simpar));
		MPI_Allreduce(&etot_l, &etot_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
		MPI_Allreduce(&nneg_l, &nneg_g, 1, MPI_INT, MPI_SUM, MPI_COMM(simpar));
		MPI_Allreduce(&ncell_l, &ncell_g, 1, MPI_INT, MPI_SUM, MPI_COMM(simpar));
		long long npair_l = gfs_pair_face_ends, npair_g = 0;
		MPI_Allreduce(&npair_l, &npair_g, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM(simpar));
		if(MYID(simpar)==0){
			static double E0 = -1;
			static int p1step = 0;
			p1step++;
			if(E0 < 0) E0 = etot_g;
			double nx = (double)NX(simpar);
			double vnorm = vmin_g * nx * nx;
			fprintf(stderr,
				"[PHASE1] step=%d dt=%.3e Vmin=%.3e Vmin*Nx2=%.3e rho_max=%.3e E=%.6e dE/E0=%.3e n_neg=%d/%d n_clip=%d n_pair=%lld\n",
				p1step, (double)Dtime, vmin_g, vnorm, rhomax_g, etot_g,
				(E0 != 0 ? (etot_g-E0)/E0 : 0), nneg_g, ncell_g, phase1_ie_clips,
				(long long)npair_g);
			fflush(stderr);
		}
	}

	/* Floor energy of the refresh below, and the clip ledger of this step. */
	double kdke_floor_l = 0;
	if(sedov_phase1_on() && GAS_ENTROPY_MODE(simpar) != 1){
		for(i=0;i<VORO_NP(simpar);i++){
			if(!targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)) continue;
			if(sbp[i].ie <= 0)
				kdke_floor_l += 1e-6*(double)sbp[i].volume/(double)(Gamma-1) - (double)sbp[i].ie;
		}
	}
	if(sedov_phase1_on()){
		double loc[5], glb[5];
		double ep_l = 0;
		for(i=0;i<VORO_NP(simpar);i++){
			if(!targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)) continue;
			ep_l += (double)sbp[i].mass*(double)kepler_phi(sbp[i].x, sbp[i].y);
		}
		loc[0] = phase1_clip_dke; loc[1] = phase1_clip_inj; loc[2] = phase1_clip_dL;
		loc[3] = (double)phase1_clip_events; loc[4] = kdke_floor_l;
		MPI_Allreduce(loc, glb, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
		double ep_g;
		MPI_Allreduce(&ep_l, &ep_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM(simpar));
		int nclip_g = 0;
		MPI_Allreduce(&phase1_ie_clips, &nclip_g, 1, MPI_INT, MPI_SUM, MPI_COMM(simpar));
		phase1_clip_dke = phase1_clip_inj = phase1_clip_dL = 0;
		phase1_clip_events = 0;
		if(MYID(simpar)==0){
			static double inj_cum = 0, dke_cum = 0, fl_cum = 0;
			static int st = 0;
			st++;
			inj_cum += glb[1]; dke_cum += glb[0]; fl_cum += glb[4];
			/* E_pot is the central-mass potential energy. E in [PHASE1]
			 * excludes it, so E + E_pot is the conserved total. */
			fprintf(stderr,
				"[PHASE1E] step=%d E_pot=%.6e clip_ev=%.0f n_clip_all=%d clip_dKE=%.3e clip_inj=%.3e clip_dL=%.3e floor_inj=%.3e cum_dKE=%.3e cum_inj=%.3e cum_floor=%.3e\n",
				st, ep_g, glb[3], nclip_g, glb[0], glb[1], glb[2], glb[4],
				dke_cum, inj_cum, fl_cum);
			fflush(stderr);
		}
	}

	/* Pressure/csound refresh */
	{
		int entropy_mode = GAS_ENTROPY_MODE(simpar);
		postype K_floor = GAS_K_FLOOR(simpar);
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(i=0;i<VORO_NP(simpar);i++){
			if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
				sbp[i].den = sbp[i].mass/sbp[i].volume;
				if(entropy_mode == 1){
					if(sbp[i].stress.K < K_floor) sbp[i].stress.K = K_floor;
					sbp[i].pressure = sbp[i].stress.K * pow((double)sbp[i].den, (double)Gamma);
					sbp[i].ie = sbp[i].pressure * sbp[i].volume / (Gamma-1);
				} else {
					sbp[i].pressure = sbp[i].ie/(sbp[i].volume*sub_gfac((treevorork4particletype*)(sbp+i))) * (Gamma-1);   /* g of the last stage geometry (diagnostic) */
					if(sbp[i].pressure <= 0){
						gfs_attempt_failed = 1;
						sbp[i].pressure = 1e-6;
						sbp[i].ie = sbp[i].pressure * sbp[i].volume*sub_gfac((treevorork4particletype*)(sbp+i)) / (Gamma-1);
					}
				}
				sbp[i].csound = sqrt(Gamma*sbp[i].pressure/sbp[i].den);
			}
			else {
				sbp[i].vx = sbp[i].vy = 0;
				sbp[i].w2 = sbp[i].rk4.w2backup;
			}
		}
	}

	/* CD10 alpha update */
	if(av_mode >= 1){
		update_alpha_cd_2d(simpar, Dtime);
	}

	_t_fin = MPI_Wtime() - _t0;
	{
		static int _step_timer = 0;
		_step_timer++;

		/* [K_DIAG] entropy_mode per-step diagnostic. Tracks K floor
		 * saturation, |dK| spike, vsig spike — to localize where the
		 * GPU+entropy_mode dt-collapse originates (Sedov2D Nx=512). */
		{
			int em_diag = GAS_ENTROPY_MODE(simpar);
			postype Kf_diag = GAS_K_FLOOR(simpar);
			double Kmin_loc = 1e300, Kmax_loc = 0;
			double absdKmax_loc = 0, vsigmax_loc = 0;
			double Pmin_loc = 1e300, Pmax_loc = 0;
			double csmax_loc = 0;
			int nfloor_loc = 0, ntot_loc = 0;
			ptrdiff_t ii;
			for(ii = 0; ii < VORO_NP(simpar); ii++){
				if(!targetBP((treevorork4particletype*)(sbp+ii), Lx, Ly)) continue;
				double Kv = (double)sbp[ii].stress.K;
				double dKv = (double)sbp[ii].stress.dK;
				double vs = (double)sbp[ii].stress.vsig_max;
				double Pv = (double)sbp[ii].pressure;
				double cs = (double)sbp[ii].csound;
				if(Kv < Kmin_loc) Kmin_loc = Kv;
				if(Kv > Kmax_loc) Kmax_loc = Kv;
				if(Pv < Pmin_loc) Pmin_loc = Pv;
				if(Pv > Pmax_loc) Pmax_loc = Pv;
				if(cs > csmax_loc) csmax_loc = cs;
				double adk = fabs(dKv);
				if(adk > absdKmax_loc) absdKmax_loc = adk;
				if(vs > vsigmax_loc) vsigmax_loc = vs;
				if(em_diag == 1 && Kv <= Kf_diag*1.0001) nfloor_loc++;
				ntot_loc++;
			}
			double Kmin_g, Kmax_g, absdKmax_g, vsigmax_g;
			double Pmin_g, Pmax_g, csmax_g;
			int nfloor_g, ntot_g;
			MPI_Allreduce(&Kmin_loc,    &Kmin_g,    1, MPI_DOUBLE, MPI_MIN, MPI_COMM(simpar));
			MPI_Allreduce(&Kmax_loc,    &Kmax_g,    1, MPI_DOUBLE, MPI_MAX, MPI_COMM(simpar));
			MPI_Allreduce(&Pmin_loc,    &Pmin_g,    1, MPI_DOUBLE, MPI_MIN, MPI_COMM(simpar));
			MPI_Allreduce(&Pmax_loc,    &Pmax_g,    1, MPI_DOUBLE, MPI_MAX, MPI_COMM(simpar));
			MPI_Allreduce(&csmax_loc,   &csmax_g,   1, MPI_DOUBLE, MPI_MAX, MPI_COMM(simpar));
			MPI_Allreduce(&absdKmax_loc,&absdKmax_g,1, MPI_DOUBLE, MPI_MAX, MPI_COMM(simpar));
			MPI_Allreduce(&vsigmax_loc, &vsigmax_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM(simpar));
			MPI_Allreduce(&nfloor_loc,  &nfloor_g,  1, MPI_INT,    MPI_SUM, MPI_COMM(simpar));
			MPI_Allreduce(&ntot_loc,    &ntot_g,    1, MPI_INT,    MPI_SUM, MPI_COMM(simpar));
			if(MYID(simpar) == 0){
				fprintf(stderr,
					"[K_DIAG] step=%d em=%d K=[%.3e..%.3e] P=[%.3e..%.3e] csmax=%.3e |dK|max=%.3e vsig_max=%.3e floor=%d/%d Dtime=%.3e\n",
					_step_timer, em_diag, Kmin_g, Kmax_g, Pmin_g, Pmax_g, csmax_g,
					absdKmax_g, vsigmax_g, nfloor_g, ntot_g, Dtime);
				fflush(stderr);
			}
		}

		double _t_total = _t_update + _t_force + _t_post + _t_fin;
		if(MYID(simpar) == 0)
			printf("[CPU step %d KDK] update=%.1fs force=%.1fs post=%.1fs final=%.1fs | total=%.1fs Dtime=%.3e\n",
				_step_timer, _t_update, _t_force, _t_post, _t_fin, _t_total, Dtime);
		fflush(stdout);
	}

	return Dtime;
}

/* ================================================================
 *  exam2d_vph_rk4_int_lagmfm (av_mode=4):
 *    **KDK leapfrog** integrator (replaces RK4 — name kept for compat).
 *    Matches GIZMO hydro_toplevel.c structure: 2 density+force evals per
 *    step vs RK4's 4 → half the work, Hopkins-reference scheme.
 *
 *    Flow:
 *      1) force eval → (a_n, du_n, Dtime_CFL)
 *      2) half-kick:  v += 0.5·Dtime·a_n;  ie += 0.5·Dtime·du_n
 *      3) drift:      x += Dtime·v  (using half-kicked v)
 *      4) postStage + migrate + re-pad
 *      5) force eval at new state → (a_{n+1}, du_{n+1})
 *      6) half-kick:  v += 0.5·Dtime·a_{n+1};  ie += 0.5·Dtime·du_{n+1}
 *
 *    Uses the SAME Dtime from step 1 throughout (second eval's CFL suggestion
 *    would be used by next call — consistent with single-step leapfrog).
 * ================================================================ */
double exam2d_vph_rk4_int_lagmfm(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		void (*postStage)(SimParameters *)
		){
	treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	postype xmin,ymin,xmax,ymax,OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);
	int av_mode = GAS_AVMODE(simpar);
	(void)av_mode;

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);

	int i;
	postype Dtime, dt;
	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) {
		sbp[i].w2old = sbp[i].w2;
		sbp[i].rk4.w2backup = sbp[i].w2;
	}
	for(i=0;i<VORO_NP(simpar);i++)
		sbp[i].stress.vsig_max = 0;

	/* === Stage 1: force eval at current state === */
#ifdef USE_CUDA
	if (GAS_GPU_ENABLED(simpar)) {
		updateDenW2Pressure2D_LagMFM_GPU(simpar, xmin,ymin,xmax,ymax,
				Gamma, paddingAllTreeParticles, mkLinkedList2D, 1);
		lagmfm_w2_post_density(simpar, 1);
		Dtime = getAccVoro2D_LagMFM_GPU(simpar, xmin,ymin,xmax,ymax,
				OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, mkLinkedList2D);
	} else
#endif
	{
		updateDenW2Pressure2D_LagMFM(simpar, xmin,ymin,xmax,ymax,
				Gamma, paddingAllTreeParticles, mkLinkedList2D, 1);
		Dtime = getAccVoro2D_LagMFM(simpar, xmin,ymin,xmax,ymax,
				OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, mkLinkedList2D);
	}
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	/* Half-kick with a_n, du_n */
	postype half_dt = 0.5 * Dtime;
	postype accx_ext = GAS_ACCX(simpar), accy_ext = GAS_ACCY(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].vx += (sbp[i].ax + accx_ext) * half_dt;
		sbp[i].vy += (sbp[i].ay + accy_ext) * half_dt;
		sbp[i].ie += sbp[i].die             * half_dt;
#if defined(LAGMFV)
		/* MFV: update particle mass from accumulated flux (half-step). */
		sbp[i].mass += sbp[i].stress.mass_flux_accum * half_dt;
		if(sbp[i].mass < 1e-30) sbp[i].mass = 1e-30;
#endif
	}

	/* Full drift with half-kicked v */
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].x += sbp[i].vx * Dtime;
		sbp[i].y += sbp[i].vy * Dtime;
	}

	postStage(simpar);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;
	for(i=0;i<VORO_NP(simpar);i++) sbp[i].stress.vsig_max = 0;

	/* === Stage 2: force eval at drifted state, then second half-kick === */
#ifdef USE_CUDA
	if (GAS_GPU_ENABLED(simpar)) {
		updateDenW2Pressure2D_LagMFM_GPU(simpar, xmin,ymin,xmax,ymax,
				Gamma, paddingAllTreeParticles, mkLinkedList2D, Dtime);
		lagmfm_w2_post_density(simpar, Dtime);
		dt = getAccVoro2D_LagMFM_GPU(simpar, xmin,ymin,xmax,ymax,
				OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, mkLinkedList2D);
	} else
#endif
	{
		updateDenW2Pressure2D_LagMFM(simpar, xmin,ymin,xmax,ymax,
				Gamma, paddingAllTreeParticles, mkLinkedList2D, Dtime);
		dt = getAccVoro2D_LagMFM(simpar, xmin,ymin,xmax,ymax,
				OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, mkLinkedList2D);
	}
	(void)dt;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	/* Second half-kick with a_{n+1}, du_{n+1} */
	for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
			sbp[i].vx += (sbp[i].ax + accx_ext) * half_dt;
			sbp[i].vy += (sbp[i].ay + accy_ext) * half_dt;
			sbp[i].ie += sbp[i].die             * half_dt;
#if defined(LAGMFV)
			sbp[i].mass += sbp[i].stress.mass_flux_accum * half_dt;
			if(sbp[i].mass < 1e-30) sbp[i].mass = 1e-30;
#endif
		} else {
			sbp[i].vx = 0;
			sbp[i].vy = 0;
		}
	}

	postStage(simpar);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	migrateTreeVorork4Particles(simpar);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	/* Refresh w2 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
			if(GAS_Kappa(simpar) < 0){
				sbp[i].w2 = -GAS_Kappa(simpar);
			}
			else if(GAS_Kappa(simpar) > 0){
				sbp[i].w2hydro = getw2forHydroParticle(simpar,(treevorork4particletype*)(sbp+i),Dtime);
				applyW2Controls(simpar, (treevorork4particletype*)(sbp+i), Dtime);
			}
		}
	}

	/* Pressure/csound refresh (volume already updated by last density pass) */
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
			if(sbp[i].volume <= 0) sbp[i].volume = cellsize*cellsize;
			sbp[i].den = sbp[i].mass/sbp[i].volume;
			sbp[i].pressure = sbp[i].ie/sbp[i].volume * (Gamma-1);
			if(sbp[i].pressure <= 0){
				sbp[i].pressure = 1e-6;
				sbp[i].ie = sbp[i].pressure * sbp[i].volume / (Gamma-1);
			}
			sbp[i].csound = sqrt(Gamma*sbp[i].pressure/sbp[i].den);
		}
		else {
			sbp[i].vx = sbp[i].vy = 0;
			sbp[i].w2 = sbp[i].rk4.w2backup;
		}
	}

	/* CD10 alpha update */
	update_alpha_cd_2d(simpar, Dtime);

	return Dtime;
}

/* ================================================================
 *  exam2d_vph_rk4_int_kNN:
 *  RK4 integrator using k-NN tree search instead of cell-based search.
 *  No cell callbacks needed — uses buildTree2D + kNN functions.
 * ================================================================ */
double exam2d_vph_rk4_int_kNN(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		int (*skipBP)(treevorork4particletype *),
		void (*postStage)(SimParameters *)
		){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	size_t p_sz = TVORORK4_DDINFO(simpar)[0].n_size;
	postype xmin, ymin, xmax, ymax, OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);

	int i;
	postype Dtime, dt;
	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;

	/* Debug: check ie at very start of RK4 */
	DEBUGPRINT("P%d RK4 start: np=%d bp=%p bp[0] ie=%g vol=%g pos=(%g,%g) PINDX=%ld\n",
		MYID(simpar), VORO_NP(simpar), (void*)bp, bp->ie, bp->volume, bp->x, bp->y, PINDX(bp));

	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->w2old = bpi->w2;
		bpi->rk4.w2backup = bpi->w2;
	}

	/* ---- RK4 Stage 1 ---- */
	{
		TStruct *tree; TPtlStruct *ptl; int nptl;
		buildTree2D(simpar, paddingAllTreeParticles, skipBP, &tree, &ptl, &nptl);
		bp = VORORK4_TBP(simpar);
		updateDenW2Pressure2D_kNN(simpar, xmin, ymin, xmax, ymax, Gamma,
				tree, ptl, nptl, 1, skipBP);
		Dtime = getAccVoro2D_kNN(simpar, xmin, ymin, xmax, ymax,
				OrderOfAccuracy, Courant, Gamma, tree, ptl, nptl, skipBP);
		my_free(ptl); my_free(tree); my_free(VORORK4_TBPP(simpar));
	}
	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->rk4.k1x = bpi->vx*Dtime;
		bpi->rk4.k1y = bpi->vy*Dtime;
		bpi->rk4.k1vx = (bpi->ax + GAS_ACCX(simpar))*Dtime;
		bpi->rk4.k1vy = (bpi->ay + GAS_ACCY(simpar))*Dtime;
		bpi->rk4.k1ie = bpi->die*Dtime;
	}
	/* RK4 second half-step */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->x += bpi->rk4.k1x*0.5;
		bpi->y += bpi->rk4.k1y*0.5;
		bpi->vx += bpi->rk4.k1vx*0.5;
		bpi->vy += bpi->rk4.k1vy*0.5;
		bpi->ie += bpi->rk4.k1ie*0.5;
	}
	postStage(simpar);
	bp = VORORK4_TBP(simpar);
	/* ---- RK4 Stage 2 ---- */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->w2 = bpi->rk4.w2backup;
	}
	{
		TStruct *tree; TPtlStruct *ptl; int nptl;
		buildTree2D(simpar, paddingAllTreeParticles, skipBP, &tree, &ptl, &nptl);
		updateDenW2Pressure2D_kNN(simpar, xmin, ymin, xmax, ymax, Gamma,
				tree, ptl, nptl, Dtime, skipBP);
		dt = getAccVoro2D_kNN(simpar, xmin, ymin, xmax, ymax,
				OrderOfAccuracy, Courant, Gamma, tree, ptl, nptl, skipBP);
		my_free(ptl); my_free(tree); my_free(VORORK4_TBPP(simpar));
	}

	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->rk4.k2x = bpi->vx*Dtime;
		bpi->rk4.k2y = bpi->vy*Dtime;
		bpi->rk4.k2vx = (bpi->ax + GAS_ACCX(simpar))*Dtime;
		bpi->rk4.k2vy = (bpi->ay + GAS_ACCY(simpar))*Dtime;
		bpi->rk4.k2ie = bpi->die*Dtime;
	}
	/* RK4 third half-step */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->x += (bpi->rk4.k2x-bpi->rk4.k1x)*0.5;
		bpi->y += (bpi->rk4.k2y-bpi->rk4.k1y)*0.5;
		bpi->vx += (bpi->rk4.k2vx-bpi->rk4.k1vx)*0.5;
		bpi->vy += (bpi->rk4.k2vy-bpi->rk4.k1vy)*0.5;
		bpi->ie += (bpi->rk4.k2ie-bpi->rk4.k1ie)*0.5;
	}
	postStage(simpar);
	bp = VORORK4_TBP(simpar);

	/* ---- RK4 Stage 3 ---- */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->w2 = bpi->rk4.w2backup;
	}
	{
		TStruct *tree; TPtlStruct *ptl; int nptl;
		buildTree2D(simpar, paddingAllTreeParticles, skipBP, &tree, &ptl, &nptl);
		updateDenW2Pressure2D_kNN(simpar, xmin, ymin, xmax, ymax, Gamma,
				tree, ptl, nptl, Dtime, skipBP);
		dt = getAccVoro2D_kNN(simpar, xmin, ymin, xmax, ymax,
				OrderOfAccuracy, Courant, Gamma, tree, ptl, nptl, skipBP);
		my_free(ptl); my_free(tree); my_free(VORORK4_TBPP(simpar));
	}

	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->rk4.k3x = bpi->vx*Dtime;
		bpi->rk4.k3y = bpi->vy*Dtime;
		bpi->rk4.k3vx = (bpi->ax + GAS_ACCX(simpar))*Dtime;
		bpi->rk4.k3vy = (bpi->ay + GAS_ACCY(simpar))*Dtime;
		bpi->rk4.k3ie = bpi->die*Dtime;
	}
	/* RK4 fourth full-step */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->x += bpi->rk4.k3x - bpi->rk4.k2x*0.5;
		bpi->y += bpi->rk4.k3y - bpi->rk4.k2y*0.5;
		bpi->vx += bpi->rk4.k3vx - bpi->rk4.k2vx*0.5;
		bpi->vy += bpi->rk4.k3vy - bpi->rk4.k2vy*0.5;
		bpi->ie += bpi->rk4.k3ie - bpi->rk4.k2ie*0.5;
	}
	postStage(simpar);
	bp = VORORK4_TBP(simpar);

	/* ---- RK4 Stage 4 ---- */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->w2 = bpi->rk4.w2backup;
	}
	{
		TStruct *tree; TPtlStruct *ptl; int nptl;
		buildTree2D(simpar, paddingAllTreeParticles, skipBP, &tree, &ptl, &nptl);
		updateDenW2Pressure2D_kNN(simpar, xmin, ymin, xmax, ymax, Gamma,
				tree, ptl, nptl, Dtime, skipBP);
		dt = getAccVoro2D_kNN(simpar, xmin, ymin, xmax, ymax,
				OrderOfAccuracy, Courant, Gamma, tree, ptl, nptl, skipBP);
		my_free(ptl); my_free(tree); my_free(VORORK4_TBPP(simpar));
	}

	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->rk4.k4x = bpi->vx*Dtime;
		bpi->rk4.k4y = bpi->vy*Dtime;
		bpi->rk4.k4vx = (bpi->ax + GAS_ACCX(simpar))*Dtime;
		bpi->rk4.k4vy = (bpi->ay + GAS_ACCY(simpar))*Dtime;
		bpi->rk4.k4ie = bpi->die*Dtime;
	}
	/* Rollback k3 step */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->x -= bpi->rk4.k3x;
		bpi->y -= bpi->rk4.k3y;
		bpi->vx -= bpi->rk4.k3vx;
		bpi->vy -= bpi->rk4.k3vy;
		bpi->ie -= bpi->rk4.k3ie;
	}
	postStage(simpar);
	bp = VORORK4_TBP(simpar);

	/* Final RK4 combination */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->w2 = bpi->rk4.w2backup;
	}
	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		if(targetBP(bpi, Lx, Ly)){
			bpi->x  += (bpi->rk4.k1x +2*bpi->rk4.k2x +2*bpi->rk4.k3x +bpi->rk4.k4x )/6.;
			bpi->y  += (bpi->rk4.k1y +2*bpi->rk4.k2y +2*bpi->rk4.k3y +bpi->rk4.k4y )/6.;
			bpi->vx += (bpi->rk4.k1vx+2*bpi->rk4.k2vx+2*bpi->rk4.k3vx+bpi->rk4.k4vx)/6.;
			bpi->vy += (bpi->rk4.k1vy+2*bpi->rk4.k2vy+2*bpi->rk4.k3vy+bpi->rk4.k4vy)/6.;
			bpi->ie += (bpi->rk4.k1ie+2*bpi->rk4.k2ie+2*bpi->rk4.k3ie+bpi->rk4.k4ie)/6.;
		}
	}
	postStage(simpar);
	bp = VORORK4_TBP(simpar);
	migrateTreeVorork4Particles(simpar);

	/* Update volume & density with final positions — uses cell-based method */
	exam2dUpdateVol(simpar, paddingAllTreeParticles,
			searchCellRk4Neighbors2D, findCellRk4BP2D, mkLinkedList2D_oExam);

	bp = VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		if(targetBP(bpi, Lx, Ly)){
			if(GAS_Kappa(simpar) < 0){
				bpi->w2 = -GAS_Kappa(simpar);
			}
			else if(GAS_Kappa(simpar) > 0){
				bpi->w2hydro = getw2forHydroParticle(simpar, bpi, Dtime);
				applyW2Controls(simpar, bpi, Dtime);
			}
		}
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		if(targetBP(bpi, Lx, Ly)){
			bpi->den = bpi->mass/bpi->volume;
			bpi->pressure = bpi->ie/bpi->volume * (Gamma-1);
			if(bpi->pressure <= 0){
				bpi->pressure = 1e-6;
				bpi->ie = bpi->pressure * bpi->volume / (Gamma-1);
			}
			bpi->csound = sqrt(Gamma*bpi->pressure/bpi->den);
		}
	}

	return Dtime;
}

double exam2d_vph_rk4_int_rt(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	postype xmin,ymin,xmax,ymax,OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);
//	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has x/y minmax= %g %g : %g %g ::: %g\n", MYID(simpar), xmin,ymin,xmax,ymax, cellsize);
//
	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);


    int i;
    postype Dtime,dt;
    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) {
		bp[i].w2old = bp[i].w2;
		bp[i].rk4.w2backup = bp[i].w2;
	}
	// Runge-Kutta 4-th order time evolution of r and vr
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D,1);
//	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed updateDenW2P...()\n", MYID(simpar));
    Dtime = getAccVoro2D_rt(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);


	bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){ 
		bp[i].rk4.k1x = bp[i].vx*Dtime;
        bp[i].rk4.k1y = bp[i].vy*Dtime;
        bp[i].rk4.k1vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k1vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k1ie = bp[i].die*Dtime;
    }
    // RK4 second
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += bp[i].rk4.k1x*0.5;
        bp[i].y += bp[i].rk4.k1y*0.5;
        bp[i].vx += bp[i].rk4.k1vx*0.5;
        bp[i].vy += bp[i].rk4.k1vy*0.5;
        bp[i].ie += bp[i].rk4.k1ie*0.5;
		if(bp[i].ie <0) DEBUGPRINT("P%d has p(%ld,%ld) : ie= %g, die= %g :: x,y= %g %g vx,y = %g %g\n", MYID(simpar), PINDX(bp+i)%nx, PINDX(bp+i)/nx, bp[i].ie, bp[i].die, bp[i].x, bp[i].y, bp[i].vx, bp[i].vy);
        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }

	migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
    dt = getAccVoro2D_rt(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
        paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k2x = bp[i].vx*Dtime;
        bp[i].rk4.k2y = bp[i].vy*Dtime;
        bp[i].rk4.k2vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k2vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k2ie = bp[i].die*Dtime;
    }
    // RK4 third
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += (bp[i].rk4.k2x-bp[i].rk4.k1x)*0.5;
        bp[i].y += (bp[i].rk4.k2y-bp[i].rk4.k1y)*0.5;
        bp[i].vx += (bp[i].rk4.k2vx-bp[i].rk4.k1vx)*0.5;
        bp[i].vy += (bp[i].rk4.k2vy-bp[i].rk4.k1vy)*0.5;
        bp[i].ie += (bp[i].rk4.k2ie-bp[i].rk4.k1ie)*0.5;
		if(bp[i].ie <0) DEBUGPRINT("P%d has p(%ld,%ld) : ie= %g, die= %g %g :: x,y= %g %g vx,y = %g %g\n", MYID(simpar), PINDX(bp+i)%nx, PINDX(bp+i)/nx, bp[i].ie, bp[i].rk4.k2ie, bp[i].rk4.k1ie, bp[i].x, bp[i].y, bp[i].vx, bp[i].vy);

        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
	dt = getAccVoro2D_rt(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k3x = bp[i].vx*Dtime;
        bp[i].rk4.k3y = bp[i].vy*Dtime;
        bp[i].rk4.k3vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k3vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k3ie = bp[i].die*Dtime;
    }
    // RK4 fourth
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += bp[i].rk4.k3x - bp[i].rk4.k2x*0.5; 
		bp[i].y += bp[i].rk4.k3y - bp[i].rk4.k2y*0.5; 
		bp[i].vx += bp[i].rk4.k3vx - bp[i].rk4.k2vx*0.5; 
		bp[i].vy += bp[i].rk4.k3vy - bp[i].rk4.k2vy*0.5; 
		bp[i].ie += bp[i].rk4.k3ie - bp[i].rk4.k2ie*0.5;
		if(bp[i].ie <0) DEBUGPRINT("P%d has p(%ld,%ld) : ie= %g, die= %g %g :: x,y= %g %g vx,y = %g %g\n", MYID(simpar), PINDX(bp+i)%nx, PINDX(bp+i)/nx, bp[i].ie, bp[i].rk4.k3ie, bp[i].rk4.k2ie, bp[i].x, bp[i].y, bp[i].vx, bp[i].vy);

        bp[i].x = fmod(bp[i].x+Lx, Lx);
		bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);



    dt = getAccVoro2D_rt(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k4x = bp[i].vx*Dtime;
        bp[i].rk4.k4y = bp[i].vy*Dtime;
        bp[i].rk4.k4vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k4vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k4ie = bp[i].die*Dtime;
    }
    // migrate particles between mpi ranks
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x -= bp[i].rk4.k3x;
        bp[i].y -= bp[i].rk4.k3y;
        bp[i].vx -= bp[i].rk4.k3vx;
        bp[i].vy -= bp[i].rk4.k3vy;
        bp[i].ie -= bp[i].rk4.k3ie;
        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
	// migrate particles between mpi ranks
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

    // to finalize ke & ie updates 
	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
    // Update position and velocity using RK4
	//
	//

    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly)) {
			bp[i].x  += (bp[i].rk4.k1x +2*bp[i].rk4.k2x +2*bp[i].rk4.k3x +bp[i].rk4.k4x )/6.; 
			bp[i].y  += (bp[i].rk4.k1y +2*bp[i].rk4.k2y +2*bp[i].rk4.k3y +bp[i].rk4.k4y )/6.; 
			bp[i].vx += (bp[i].rk4.k1vx+2*bp[i].rk4.k2vx+2*bp[i].rk4.k3vx+bp[i].rk4.k4vx)/6.; 
			bp[i].vy += (bp[i].rk4.k1vy+2*bp[i].rk4.k2vy+2*bp[i].rk4.k3vy+bp[i].rk4.k4vy)/6.; 
			bp[i].ie += (bp[i].rk4.k1ie+2*bp[i].rk4.k2ie+2*bp[i].rk4.k3ie+bp[i].rk4.k4ie)/6.;
		   if(bp[i].ie <0) DEBUGPRINT("P%d has p(%ld,%ld) : ie= %g, die= %g %g %g %g :: x,y= %g %g vx,y = %g %g\n", MYID(simpar), PINDX(bp+i)%nx, PINDX(bp+i)/nx, bp[i].ie, 
				   bp[i].rk4.k1ie, bp[i].rk4.k2ie, 
				   bp[i].rk4.k3ie, bp[i].rk4.k4ie, 
				   bp[i].x, bp[i].y, bp[i].vx, bp[i].vy);
	        bp[i].x = fmod(bp[i].x+Lx,Lx);
			bp[i].y = fmod(bp[i].y+Ly,Ly);
		}
    }

    // migrate particles between mpi ranks
    migrateTreeVorork4Particles(simpar);

	//centroid shift
	if(GAS_FCENTROID(simpar)>0) 
		exam2d_centroidShift(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	// update volume & density with one-step evolved positions 
    exam2dUpdateVol(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	bp = VORORK4_TBP(simpar);
	postype dmean = SIMBOX(simpar).x.max/NX(simpar);
	// split getw2forHydroParticle and pressure update because the former uses the latter.
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
			if(GAS_Kappa(simpar) <0) { 
				bp[i].w2 = -GAS_Kappa(simpar); 
			} 
			else if (GAS_Kappa(simpar) >0){ 
          		bp[i].w2hydro = getw2forHydroParticle(simpar,(bp+i),Dtime);
				applyW2Controls(simpar, bp+i, Dtime);
			} 
		}
	}
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
			bp[i].den = bp[i].mass/bp[i].volume; 
			bp[i].pressure = bp[i].ie/bp[i].volume * (Gamma-1); 
			bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den); 
			if(bp[i].volume ==0) 
				DEBUGPRINT("P%d has wrong volume at i= %d with den= %g at xy= %g %g\n", MYID(simpar), 
						i,bp[i].den, bp[i].x, bp[i].y); 
			if(isnan(bp[i].csound)){ 
				DEBUGPRINT("P%d has nan at the sound speed = %d with" 
						"id= %ld: pressure= %g den= %g Delta IE=%g dIE= %g Dt= %g\n", 
						MYID(simpar), i,PINDX(bp+i), bp[i].pressure, bp[i].den, bp[i].ie,bp[i].die,Dtime); 
				exit(9);
			}
		}
		else {
			bp[i].vx = bp[i].vy = 0;
			bp[i].w2= bp[i].rk4.w2backup;
		}
    }
    return Dtime;
}
