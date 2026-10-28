#ifndef TCF_FULL_FACE_PROBE_H
#define TCF_FULL_FACE_PROBE_H
/* Diagnostic wrapper of the shared default TCF assembly.
 * Inputs A,d,ex,ey,tx,ty,Pi,Pj,vix,viy,vjx,vjy,Q.
 * Outputs Fx,Fy,Edot,Vdot,heat,0. Production uses tcf_face_terms directly. */
#include "tcf_face_terms.h"
#ifdef __CUDACC__
#define TFF_INLINE static __host__ __device__ __forceinline__
#else
#define TFF_INLINE static inline
#endif
TFF_INLINE void tcf_full_face_probe(const double*a,double*o){
 double A=a[0],d=a[1],ex=a[2],ey=a[3],tx=a[4],ty=a[5];
 double Pi_=a[6],Pj_=a[7],vix=a[8],viy=a[9],vjx=a[10],vjy=a[11],Q=a[12];
 double ui=vix*ex+viy*ey,uj=vjx*ex+vjy*ey;
 TcfFaceTerms q=tcf_face_terms(A,d,ex,ey,tx,ty,Pi_,Pj_,vix,viy,vjx,vjy,ui,uj,Q);
 o[0]=q.normal_x+q.tangent_x;o[1]=q.normal_y+q.tangent_y;
 o[2]=q.energy;o[3]=q.volume;o[4]=q.heat;o[5]=0;
}
#undef TFF_INLINE
#endif
