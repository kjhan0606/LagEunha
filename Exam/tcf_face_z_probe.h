#ifndef TCF_FACE_Z_PROBE_H
#define TCF_FACE_Z_PROBE_H
#include <math.h>
#ifdef __CUDACC__
#define TF_INLINE static __host__ __device__ __forceinline__
#else
#define TF_INLINE static inline
#endif
/* Experimental isolated active FACE_Z / QLAW=1 algebra, no production callers.
 * Inputs: A,ex,ey,Pi,Pj,rhoi,rhoj,fi,fj,qin,du_cell,mean_un,gamma.
 * fi/fj and qin are upstream recorded data, not reconstructed here.
 * Outputs: Zq,Q,FQx,FQy,EQ,heat_i. No reversible pressure/HG/volume term.
 */
TF_INLINE void tcf_face_z_probe(const double *a,double *o){
 double mm=(a[9]*a[10]<=0)?0:(fabs(a[9])<fabs(a[10])?a[9]:a[10]);
 double fi=a[7]<1e-12?1e-12:a[7],fj=a[8]<1e-12?1e-12:a[8];
 double rl=a[5]*fi,rr=a[6]*fj,pl=a[3]*pow(fi,a[12]),pr=a[4]*pow(fj,a[12]);
 double cl=sqrt(a[12]*pl/rl),cr=sqrt(a[12]*pr/rr),cm=cl>cr?cl:cr;
 double mu=rl*rr/(rl+rr),zd=mu*(cm+(mm>0?mm:0));
 double ah=(a[12]+1)*fabs(mm)/8,zh=mu*(ah+sqrt(ah*ah+cm*cm));
 double z=zd>zh?zd:zh,Q=z*mm;
 o[0]=z;o[1]=Q;o[2]=-a[0]*Q*a[1];o[3]=-a[0]*Q*a[2];
 o[4]=-a[0]*Q*a[11];o[5]=.5*a[0]*Q*a[10];
}
#undef TF_INLINE
#endif
