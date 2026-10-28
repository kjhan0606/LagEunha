#ifndef LAGFORCE_WEIGHT_REMAP_POINT_ENERGY_H
#define LAGFORCE_WEIGHT_REMAP_POINT_ENERGY_H
#include "weight_remap_flux.h"
#ifdef __CUDACC__
#define WP_INLINE static __host__ __device__ __forceinline__
#else
#define WP_INLINE static inline
#endif
/* Diagnostic ONLY: frozen-generator weight remap with point potential m*phi.
 * Input: dm(j->i), vi[2], vj[2], ei, ej, phi_i, phi_j.
 * Output: shared m,px,py,hydro-energy flux plus local energy sources si,sj.
 * E_i'=flux.E+si, E_j'=-flux.E+sj; phi is static during this remap.
 * Centered face potential gives si=sj=-dm*(phi_i-phi_j)/2.
 * Thus sum(E_hydro+m*phi)'=0. With the kinetic-compatible flux,
 * sum(K)'=0, so sum(U)'=-sum(phi*m'). This is NOT heat-free relabeling
 * and does NOT repair point angular momentum. Do not enable in production.
 */
WP_INLINE void wp_pair(const double *a, double *b) {
    WrFlux f=wr_pair_flux(WR_KINETIC_THERMAL_DONOR,a[0],a[1],a[2],a[3],a[4],a[5],a[6]);
    double correction=-0.5*a[0]*(a[7]-a[8]);
    b[0]=f.mass; b[1]=f.px; b[2]=f.py; b[3]=f.energy;
    b[4]=correction; b[5]=correction;
}
#undef WP_INLINE
#endif
