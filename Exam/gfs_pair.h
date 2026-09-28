#ifndef GFS_PAIR_H
#define GFS_PAIR_H
/*
 * Short-range pair restoring pressure for the geometric-face scheme.
 * Identically zero for d >= d_c. The length entering d_c is sqrt(area)
 * in 2D and the cell length in 1D:
 *   d_c = 0.2 * min(length_i, length_j)
 *   P   = 16 * rho_bar * (c_bar^2 + min(v_close^2, c_bar^2)) * (1 - d/d_c)^2
 * v_close is the approaching speed along the pair (zero when separating
 * or already at rest). The value and the first derivative vanish at d_c.
 * As d -> 0 the pressure stays finite. A stationary pair is still pushed
 * by the thermal piece 16 rho c^2 (1-d/d_c)^2. The ram piece is capped at
 * the thermal piece, so P never exceeds 32 rho_bar c_bar^2.
 * Why the cap. The uncapped 16 rho v_close^2 is large on ordinary shear
 * pairs in a cold disk (Kepler job 406510 gained energy through it).
 * Thermal only (cap 16) does not stop the Mach-2 approach in the 1D
 * Woodward-Colella blast, and neither does 20. 24 and above finish.
 * 32 is the smallest round value with margin, and it gives the same blast
 * L1 and energy error as the uncapped law. Callers add P to the face
 * pressure.
 */
#include <math.h>

#ifdef __CUDACC__
__host__ __device__
#endif
static inline double gfs_pair_pressure_len(
		double dist, double len_i, double len_j,
		double rho_i, double rho_j, double cs_i, double cs_j,
		double vclose)
{
	double dc, q, rho, c2, s2, p;
	if(!(dist > 0.0) || !(len_i > 0.0) || !(len_j > 0.0)) return 0.0;
	dc = 0.2 * (len_i < len_j ? len_i : len_j);
	if(!(dist < dc)) return 0.0;
	q = 1.0 - dist / dc;
	rho = 0.5 * (rho_i + rho_j);
	c2 = 0.5 * (cs_i * cs_i + cs_j * cs_j);
	if(vclose < 0.0) vclose = 0.0;
	if(vclose > 0.0 && vclose * vclose > c2)
		s2 = 2.0 * c2;          /* ram piece capped at the thermal piece */
	else
		s2 = c2 + vclose * vclose;
	if(!(rho > 0.0) || !(s2 > 0.0) || !(q > 0.0)) return 0.0;
	p = 16.0 * rho * s2 * q * q;
	if(!isfinite(p) || !(p > 0.0)) return 0.0;
	return p;
}

/* Cap the pair impulse so its discrete work cannot exceed the spendable
 * internal energy B. q = (v_j - v_i)·n with q < 0 approaching.
 * nshare is the number of faces that share one cell's budget (6 in 2D, 2 in 1D).
 * W(J) = q J + J^2/(2 mu). The returned pressure is the one actually added. */
#ifdef __CUDACC__
__host__ __device__
#endif
static inline double gfs_pair_work_limit(
		double p0, double area, double dt,
		double m_i, double m_j, double q,
		double ie_i, double ie_j, double nshare)
{
	double mu, B, disc, jmax, pmax;
	if(!(p0 > 0.0)) return 0.0;
	if(!(area > 0.0) || !(dt > 0.0) || !(m_i > 0.0) || !(m_j > 0.0))
		return p0;
	if(!(nshare > 0.0)) nshare = 1.0;
	mu = m_i * m_j / (m_i + m_j);
	/* Budget from the poorer cell: B = 2 min(ie_i, ie_j) / nshare.
	 * The old budget (ie_i + ie_j) / nshare pooled both cells, but the
	 * work is not charged to the pair as a whole. Each end pays P dV of
	 * its own volume change through this face, and on a separating pair
	 * next to a Laguerre face almost all of it lands on one cell. With a
	 * hot j and a cold i the pooled budget let j's energy pay for a pair
	 * pressure that i alone had to do work against: in the 64^2 Kepler
	 * run cells at R ~ 0.02 lost 500 times their ie in one stage through
	 * the pair term, went to ie < 0 and were refilled by the stage floor
	 * (the whole energy error, see code_review_grokbot.md). For equal ie
	 * the budget is the same as before; it is symmetric in i and j, so
	 * both ends of the face still get the same pressure. */
	B = 0.0;
	if(ie_i > 0.0 && ie_j > 0.0)
		B = 2.0 * (ie_i < ie_j ? ie_i : ie_j) / nshare;
	disc = q * q + 2.0 * B / mu;
	if(!(disc > 0.0)) return 0.0;
	jmax = mu * (sqrt(disc) - q);
	if(!(jmax > 0.0)) return 0.0;
	pmax = jmax / (area * dt);
	return (p0 < pmax) ? p0 : pmax;
}

/* Second cap, on what each end actually pays. The ie of cell i changes by
 * -P s_i dt on this face, s_i = (u_face - v_i)·dS the rate at which the
 * face sweeps i's volume. s_i + s_j = (v_j - v_i)·dS, but the split is set
 * by the face velocity (Laguerre weights, face rotation), not by q: near
 * the Kepler centre a cold cell next to a hot one had s_i several times
 * q A, so the pair pressure allowed by gfs_pair_work_limit still took
 * about twice its ie per stage (stage floor, energy error). Only an end
 * that expands (s > 0) pays; cap P so that P s dt <= ie / nshare there.
 * s_j is computed from s_i on the i end and vice versa, so both ends get
 * the same pressure as long as they see the same face velocity. */
#ifdef __CUDACC__
__host__ __device__
#endif
static inline double gfs_pair_charge_limit(
		double p0, double s_i, double s_j, double dt,
		double ie_i, double ie_j, double nshare)
{
	double p = p0, pm;
	if(!(p0 > 0.0)) return 0.0;
	if(!(dt > 0.0)) return p0;
	if(!(nshare > 0.0)) nshare = 1.0;
	if(s_i > 0.0){
		pm = (ie_i > 0.0 ? ie_i : 0.0) / (nshare * s_i * dt);
		if(pm < p) p = pm;
	}
	if(s_j > 0.0){
		pm = (ie_j > 0.0 ? ie_j : 0.0) / (nshare * s_j * dt);
		if(pm < p) p = pm;
	}
	if(!isfinite(p) || !(p > 0.0)) return 0.0;
	return p;
}

#ifdef __CUDACC__
__host__ __device__
#endif
static inline double gfs_pair_pressure(
		double dramp, double vol_i, double vol_j,
		double rho_i, double rho_j, double cs_i, double cs_j,
		double vclose)
{
	if(!(vol_i > 0.0) || !(vol_j > 0.0)) return 0.0;
	return gfs_pair_pressure_len(dramp, sqrt(vol_i), sqrt(vol_j),
			rho_i, rho_j, cs_i, cs_j, vclose);
}
#endif
