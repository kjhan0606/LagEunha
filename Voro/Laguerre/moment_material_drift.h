#ifndef LAGFORCE_MOMENT_MATERIAL_DRIFT_H
#define LAGFORCE_MOMENT_MATERIAL_DRIFT_H
#ifdef __CUDACC__
#define MD_INLINE static __host__ __device__ __forceinline__
#else
#define MD_INLINE static inline
#endif
/* Instantaneous force-free material source, 17-moment layout.
 * This is the dt=0 derivative of md_drift, with zero potential source.
 * Requires moving-face fluxes and force/thermal sources for a complete RHS.
 */
MD_INLINE void md_material_rate(const double *a,double *b) {
    for(int k=0;k<17;k++)b[k]=0.0;
    b[1]=a[6];b[2]=a[7];
    b[3]=2*a[8];b[4]=a[9]+a[10];b[5]=2*a[11];
    b[8]=a[12];b[9]=a[13];b[10]=a[13];b[11]=a[14];
}
/* Exact force-free MATERIAL moment drift x'=x+dt*v, v'=v.
 * Input first16 hydro moments (weight_remap_moments.h layout), then dt.
 * Output16 hydro moments. Potential energy is deliberately NOT represented:
 * external gravity/potential work must be closed by the caller separately.
 * Includes full Pi, hence unresolved velocity covariance in Q's dt^2 term.
 * This follows the old material, NOT new power cells; moving-face moment
 * flux/remap is still required. Never use this as a complete ALE update.
 */
MD_INLINE int md_drift(const double *a,double *b) {
    if(!(a[0]>0.0))return 1;
    double dt=a[16],dt2=dt*dt;
    for(int k=0;k<16;k++)b[k]=a[k];
    b[1]+=dt*a[6];b[2]+=dt*a[7];
    b[3]+=2*dt*a[8]+dt2*a[12];
    b[4]+=dt*(a[9]+a[10])+dt2*a[13];
    b[5]+=2*dt*a[11]+dt2*a[14];
    b[8]+=dt*a[12];b[9]+=dt*a[13];
    b[10]+=dt*a[13];b[11]+=dt*a[14];
    return 0;
}
#undef MD_INLINE
#endif
