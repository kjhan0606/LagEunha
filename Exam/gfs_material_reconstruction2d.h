#ifndef GFS_MATERIAL_RECONSTRUCTION2D_H
#define GFS_MATERIAL_RECONSTRUCTION2D_H
#include "gfs_velocity_limiter.h"
#ifdef __CUDACC__
#define GFS_REC_HD static __host__ __device__ inline
#else
#define GFS_REC_HD static inline
#endif
/* Unweighted angular least squares, using actual cell-centroid offsets.
 * Rank failure is explicit; it is not replaced by a zero gradient. */
GFS_REC_HD int gfs_material_fit2d(int n,const double *dx,const double *dy,
 const double *dq,double gradient[2],double *condition){
 double xx=0,xy=0,yy=0,bx=0,by=0;
 for(int i=0;i<n;i++){
  double r2=dx[i]*dx[i]+dy[i]*dy[i];if(!(r2>0)||!isfinite(dq[i]))return 0;
  xx+=dx[i]*dx[i]/r2;xy+=dx[i]*dy[i]/r2;yy+=dy[i]*dy[i]/r2;
  bx+=dx[i]*dq[i]/r2;by+=dy[i]*dq[i]/r2;
 }
 *condition=gfs_velocity_matrix_condition(xx,xy,yy);
 if(!isfinite(*condition)||*condition*DBL_EPSILON>1e-6)return 0;
 double det=xx*yy-xy*xy;
 gradient[0]=(yy*bx-xy*by)/det;gradient[1]=(xx*by-xy*bx)/det;
 return isfinite(gradient[0])&&isfinite(gradient[1]);
}
/* A single scalar limiter bounds every reconstructed vertex value by the
 * declared neighboring cell extrema. No coefficient or post-update repair. */
GFS_REC_HD double gfs_material_scalar_alpha2d(int n,const double *dx,const double *dy,
 double gx,double gy,double q,double low,double high){
 double alpha=1;
 for(int i=0;i<n;i++){
  double d=gx*dx[i]+gy*dy[i];
  if(d>0)alpha=fmin(alpha,(high-q)/d);
  else if(d<0)alpha=fmin(alpha,(low-q)/d);
 }
 return fmax(0,fmin(1,alpha));
}
#undef GFS_REC_HD
#endif
