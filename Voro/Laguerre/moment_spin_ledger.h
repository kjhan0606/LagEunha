#ifndef LAGFORCE_MOMENT_SPIN_LEDGER_H
#define LAGFORCE_MOMENT_SPIN_LEDGER_H
#include <math.h>
#ifdef __CUDACC__
#define MSL_INLINE static __host__ __device__ __forceinline__
#else
#define MSL_INLINE static inline
#endif
/* Diagnostic chain rule for the existing 17 raw-moment layout.
 * q: state; d: COMPLETE moment rate, including ALE fluxes/volume sources.
 * out: I,L_internal,omega,Krot,Idot,Ldot,Krot_dot,thermal_dot.
 * No forces, heat source, relaxation law or realizability repair are added.
 * Krot is a projection of existing kinetic energy, NOT extra total energy.
 * Return 1 for invalid/nonfinite input, 2 for nonpositive central inertia.
 * Positive inertia alone does not establish full moment realizability.
 */
MSL_INLINE int msl_rate(const double *q,const double *d,double *out) {
    for(int k=0;k<17;k++)if(!isfinite(q[k]) || !isfinite(d[k]))return 1;
    const double m=q[0];if(!(m>0.0))return 1;
    const double cx=q[1]/m,cy=q[2]/m;
    const double I=q[3]+q[5]-cx*q[1]-cy*q[2];
    if(!(I>0.0) || !isfinite(I))return 2;
    const double orbital=cx*q[7]-cy*q[6];
    const double L=q[9]-q[10]-orbital,omega=L/I;
    const double Id=d[3]+d[5]-2.0*(cx*d[1]+cy*d[2])+(cx*cx+cy*cy)*d[0];
    const double Ld=d[9]-d[10]-(d[1]*q[7]-d[2]*q[6])/m
                    -(cx*d[7]-cy*d[6])+orbital*d[0]/m;
    const double v[8]={I,L,omega,0.5*L*omega,Id,Ld,
                       omega*Ld-0.5*omega*omega*Id,
                       d[15]-0.5*(d[12]+d[14])};
    for(int k=0;k<8;k++)if(!isfinite(v[k]))return 1;
    for(int k=0;k<8;k++)out[k]=v[k];
    return 0;
}
#undef MSL_INLINE
#endif
