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
