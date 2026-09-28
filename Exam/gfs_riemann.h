#ifndef GFS_RIEMANN_H
#define GFS_RIEMANN_H
/*
 * Face Riemann solvers for the 3D GFS path (Exam/exam3d_gfs.c).
 *
 * These are copies of the 2D production solvers in Exam/exam.c:
 *   gfs_hll_star_state        = hll_star_state          (exam.c, "static inline void hll_star_state")
 *   gfs_hllc_face             = hllc_face_2d            (exam.c)
 *   gfs_hllc_face_rest_frame  = hllc_face_2d_rest_frame (exam.c)
 * The only change is that the SEDOV_PHASE1 switch is a parameter instead of
 * the cached getenv inside exam.c. The solvers are 1D along the face normal,
 * so nothing in them depends on the dimension. exam.c is deliberately not
 * changed to include this header (the 2D objects stay bit-identical); if one
 * of the two copies is edited, edit the other one too.
 */
#include <math.h>

#ifndef GFS_HLLC_VACUUM_PRATIO
#define GFS_HLLC_VACUUM_PRATIO 100.0      /* = HLLC_VACUUM_PRATIO in exam.c */
#endif
#ifndef GFS_HLLC_VACUUM_PMIN_GUARD
#define GFS_HLLC_VACUUM_PMIN_GUARD 0.01   /* = HLLC_VACUUM_PMIN_GUARD in exam.c */
#endif

static inline void gfs_hll_star_state(
		double rhoL, double pL, double vnL, double cL,
		double rhoR, double pR, double vnR, double cR,
		double Gamma, double *pstar, double *vnstar)
{
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
	if(gm1 < 1e-8) gm1 = 1e-8;
	double rhoLs = fmax(rhoL, 1e-30);
	double rhoRs = fmax(rhoR, 1e-30);
	double EL = pL/(gm1*rhoLs) + 0.5*vnL*vnL;
	double ER = pR/(gm1*rhoRs) + 0.5*vnR*vnR;
	double EdenL = rhoL * EL, EdenR = rhoR * ER;
	double denom = SR - SL;
	if(fabs(denom) < 1e-30){
		*pstar = 0.5*(pL+pR); *vnstar = 0.5*(vnL+vnR); return;
	}
	double FrL = rhoL*vnL, FrR = rhoR*vnR;
	double FmL = rhoL*vnL*vnL + pL, FmR = rhoR*vnR*vnR + pR;
	double FeL = vnL*(EdenL + pL), FeR = vnR*(EdenR + pR);
	double rhoS = (SR*rhoR - SL*rhoL - (FrR - FrL)) / denom;
	double momS = (SR*(rhoR*vnR) - SL*(rhoL*vnL) - (FmR - FmL)) / denom;
	double ES   = (SR*EdenR - SL*EdenL - (FeR - FeL)) / denom;
	if(!(rhoS > 1e-30) || isnan(rhoS) || isinf(rhoS)
			|| isnan(momS) || isinf(momS) || isnan(ES) || isinf(ES)){
		*pstar = 0.5*(pL+pR); *vnstar = 0.5*(vnL+vnR); return;
	}
	double vnS = momS / rhoS;
	double eint = ES/rhoS - 0.5*vnS*vnS;
	double PS = gm1 * rhoS * eint;
	if(!(PS > 0) || isnan(PS) || isinf(PS) || isnan(vnS) || isinf(vnS)){
		*pstar = 0.5*(pL+pR); *vnstar = 0.5*(vnL+vnR); return;
	}
	*pstar = PS; *vnstar = vnS;
}

static inline void gfs_hllc_face(
		double rhoL, double pL, double vnL, double cL,
		double rhoR, double pR, double vnR, double cR,
		double Gamma, int phase1, double *pstar, double *vnstar)
{
	const double tiny = 1.0e-30;
	double S_L, S_R, S_M, P_M;
	double cmax = cL > cR ? cL : cR;

	if(phase1){
		double pmin = pL < pR ? pL : pR;
		double pmax = pL > pR ? pL : pR;
		if(pmax > 100.0 * fmax(pmin, 1.0e-30)){
			gfs_hll_star_state(rhoL, pL, vnL, cL, rhoR, pR, vnR, cR, Gamma, pstar, vnstar);
			return;
		}
	}
	{	/* Fix C v2 */
		double pmin = pL < pR ? pL : pR;
		double pmax = pL > pR ? pL : pR;
		if(pmin > GFS_HLLC_VACUUM_PMIN_GUARD && pmax > GFS_HLLC_VACUUM_PRATIO * pmin){
			*pstar = pmin; *vnstar = 0.5*(vnL + vnR); return;
		}
	}
	if((vnR - vnL) > cmax){ *pstar = tiny; *vnstar = 0.5*(vnL + vnR); return; }
	{	/* Gaburov simplest HLLC */
		double vmin = vnL < vnR ? vnL : vnR;
		double vmax = vnL > vnR ? vnL : vnR;
		S_L = vmin - cmax;
		S_R = vmax + cmax;
		double rho_wt_L = rhoL * (S_L - vnL);
		double rho_wt_R = rhoR * (S_R - vnR);
		double denom = rho_wt_L - rho_wt_R;
		if(fabs(denom) > tiny){
			S_M = ((pR - pL) + rho_wt_L*vnL - rho_wt_R*vnR) / denom;
			P_M = (pL*rho_wt_R - pR*rho_wt_L + rho_wt_L*rho_wt_R*(vnR - vnL))
			      / (rho_wt_R - rho_wt_L);
			if(P_M > 0 && !isnan(P_M) && !isinf(P_M)){ *pstar = P_M; *vnstar = S_M; return; }
		}
	}
	{	/* Roe fallback */
		double sqL = sqrt(rhoL), sqR = sqrt(rhoR);
		double sq_inv = 1.0/(sqL + sqR + tiny);
		double vn_roe = (sqL*vnL + sqR*vnR) * sq_inv;
		double h_L = Gamma*pL/((Gamma - 1.0)*rhoL) + 0.5*vnL*vnL;
		double h_R = Gamma*pR/((Gamma - 1.0)*rhoR) + 0.5*vnR*vnR;
		double h_roe = (sqL*h_L + sqR*h_R) * sq_inv;
		double c2_roe = (Gamma - 1.0) * (h_roe - 0.5*vn_roe*vn_roe);
		if(c2_roe < tiny) c2_roe = tiny;
		double c_roe = sqrt(c2_roe);
		double srA = vnR + cR, srB = vn_roe + c_roe;
		double slA = vnL - cL, slB = vn_roe - c_roe;
		S_R = srA > srB ? srA : srB;
		S_L = slA < slB ? slA : slB;
		double rho_wt_R =  rhoR * (S_R - vnR);
		double rho_wt_L = -rhoL * (S_L - vnL);
		double denom = rho_wt_R + rho_wt_L;
		if(fabs(denom) > tiny){
			S_M = (rho_wt_R*vnR + rho_wt_L*vnL + (pL - pR)) / denom;
			P_M = rhoL*(vnL - S_L)*(vnL - S_M) + pL;
			if(P_M > 0 && !isnan(P_M) && !isinf(P_M)){ *pstar = P_M; *vnstar = S_M; return; }
		}
	}
	{	/* Rusanov-like last resort */
		P_M = 0.5*((pL + pR) + (vnL - vnR)*0.25*(rhoL + rhoR)*(cL + cR));
		S_M = 0.5*(vnR + vnL) + 2.0*(pL - pR)/((rhoL + rhoR)*(cL + cR) + tiny);
		double a1 = fabs(vnL - cL), a2 = fabs(vnR - cR);
		double a3 = fabs(vnL + cL), a4 = fabs(vnR + cR);
		double S_plus = a1;
		if(a2 > S_plus) S_plus = a2;
		if(a3 > S_plus) S_plus = a3;
		if(a4 > S_plus) S_plus = a4;
		if(S_M < -S_plus) S_M = -S_plus;
		if(S_M >  S_plus) S_M =  S_plus;
		if(P_M < tiny)    P_M = tiny;
		*pstar = P_M; *vnstar = S_M;
	}
}

static inline void gfs_hllc_face_rest_frame(
		double rhoL, double pL, double vnL_lab, double cL,
		double rhoR, double pR, double vnR_lab, double cR,
		double wn, double Gamma, int phase1,
		double *pstar, double *vnstar_lab)
{
	double pst, vnst;
	gfs_hllc_face(rhoL, pL, vnL_lab - wn, cL, rhoR, pR, vnR_lab - wn, cR,
			Gamma, phase1, &pst, &vnst);
	*pstar = pst;
	*vnstar_lab = vnst + wn;
}
#endif
