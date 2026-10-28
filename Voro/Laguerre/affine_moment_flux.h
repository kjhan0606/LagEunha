#ifndef LAGFORCE_AFFINE_MOMENT_FLUX_H
#define LAGFORCE_AFFINE_MOMENT_FLUX_H
#include "weight_remap_moments.h"
#ifdef __CUDACC__
#define AMF_INLINE static __host__ __device__ __forceinline__
#else
#define AMF_INLINE static inline
#endif
/* Experimental transport-only rate in the existing 17-moment layout.
 * Input first 11 fields match affine_face_transport.h; fields 11,12 are
 * specific thermal energy at midpoint and its tangential derivative.
 * Signed quadrature mass is a FLUX, not a positive material atom.
 * Two Gauss nodes exactly integrate this affine trace model (degree <=3).
 * Pressure forces/work, volume sources and gravity are NOT included.
 * Potential-energy flux is explicitly zero; not valid as a gravity closure.
 */
AMF_INLINE void affine_moment_flux(const double *a,double *out) {
    for(int k=0;k<17;k++)out[k]=0.0;
    for(int sign=-1;sign<=1;sign+=2){
        const double z=sign*a[0]*0.28867513459481288225;
        const double atom[7]={0.5*a[0]*(a[1]+a[2]*z),
            a[3]+a[5]*z,a[4]+a[6]*z,a[7]+a[9]*z,a[8]+a[10]*z,
            a[11]+a[12]*z,0.0};
        double contribution[17];wm_atom(atom,contribution);
        for(int k=0;k<17;k++)out[k]+=contribution[k];
    }
}
#undef AMF_INLINE
#endif
