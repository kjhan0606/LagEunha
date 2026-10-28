#ifndef LAGFORCE_MOMENT_PRESSURE_RATE_H
#define LAGFORCE_MOMENT_PRESSURE_RATE_H
#ifdef __CUDACC__
#define MPR_INLINE static __host__ __device__ __forceinline__
#else
#define MPR_INLINE static inline
#endif
/* Experimental smooth affine pressure/velocity moments, not a Riemann solver.
 * Face input: A,P0,Pt,cx,cy,tx,ty,nx,ny,vx,vy,gtx,gty.
 * P=P0+Pt*z, v=v0+gt*z, position=c+t*z. Outward unit normal n.
 * Output uses the existing17moment layout; no advective flux or gravity.
 */
MPR_INLINE void mpr_face(const double *a,double *o) {
    for(int k=0;k<17;k++)o[k]=0.0;
    const double p=a[0]*a[1],s=a[0]*a[0]*a[0]*a[2]/12.0;
    const double xp[2]={p*a[3]+s*a[5],p*a[4]+s*a[6]};
    const double vp[2]={p*a[9]+s*a[11],p*a[10]+s*a[12]};
    o[6]=-p*a[7];o[7]=-p*a[8];
    for(int i=0;i<2;i++)for(int j=0;j<2;j++)o[8+2*i+j]=-xp[i]*a[7+j];
    o[12]=-2*vp[0]*a[7];o[13]=-vp[0]*a[8]-vp[1]*a[7];o[14]=-2*vp[1]*a[8];
    o[15]=-vp[0]*a[7]-vp[1]*a[8];
}
/* Input: integral_cell P dV, Gxx,Gxy,Gyx,Gyy, with G_ab=partial_b v_a.
 * Valid for spatially constant G. Necessary volume terms from integration by
 * parts; adding only face pressure terms produces an incorrect moment RHS.
 */
MPR_INLINE void mpr_volume(const double *a,double *o) {
    for(int k=0;k<17;k++)o[k]=0.0;
    o[8]=o[11]=a[0];
    o[12]=2*a[0]*a[1];o[13]=a[0]*(a[2]+a[3]);o[14]=2*a[0]*a[4];
}
#undef MPR_INLINE
#endif
