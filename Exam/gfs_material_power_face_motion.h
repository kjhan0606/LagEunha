#ifndef GFS_MATERIAL_POWER_FACE_MOTION_H
#define GFS_MATERIAL_POWER_FACE_MOTION_H
#include <math.h>
#include "gfs_material_swept_measure.h"
#ifdef __CUDACC__
#define GFS_POWER_MOTION_HD static __host__ __device__ inline
#else
#define GFS_POWER_MOTION_HD static inline
#endif
/* One directed power-cell face, all dimensions. separation = xj-xi;
 * centroid is the face centroid relative to xi; measure is length/area in
 * 2D/3D (one for a 1D endpoint). The exact affine plane velocity gives
 * Vdot_face = left.vi + right.vj + weight*(Wi_dot-Wj_dot).
 * For an image xj=g(x_owner), use J*v_owner here and J^T*right for
 * the transpose force. A fluid reflection is not generally that J.
 * This is geometry only, not a pressure law or a mass flux. */
typedef struct {
    double left[3], right[3], weight;
} GfsMaterialPowerFaceMotion;
GFS_POWER_MOTION_HD int gfs_material_power_face_motion(int dim,
 const double *separation,const double *centroid,double measure,
 GfsMaterialPowerFaceMotion *out){
 if(dim<1||dim>3||!separation||!centroid||!out||!(measure>0)||!isfinite(measure))return 0;
 for(int a=0;a<dim;a++)
  if(!isfinite(separation[a])||!isfinite(centroid[a]))return 0;
 /* hypot(0,x) is exactly |x|; avoid an extra library norm per face.
  * Subsequent hypot calls retain overflow/underflow-safe 2D/3D norms. */
 double d=fabs(separation[0]);
 for(int a=1;a<dim;a++)d=hypot(d,separation[a]);
 if(!(d>0)||!isfinite(d))return 0;
 out->weight=measure/(2*d);
 if(!(out->weight>0)||!isfinite(out->weight))return 0;
 for(int a=0;a<3;a++){
  out->left[a]=a<dim?measure*(centroid[a]/d):0;
  out->right[a]=a<dim?measure*((separation[a]-centroid[a])/d):0;
  if(!isfinite(out->left[a])||!isfinite(out->right[a]))return 0;
 }
 return 1;
}
GFS_POWER_MOTION_HD int gfs_material_power_face_motion_rate(int dim,
 const GfsMaterialPowerFaceMotion *f,const double *vi,const double *vj,
 double wi_dot,double wj_dot,double *rate){
 if(dim<1||dim>3||!f||!vi||!vj||!rate||!isfinite(wi_dot)||!isfinite(wj_dot))return 0;
 double v=f->weight*(wi_dot-wj_dot);
 for(int a=0;a<dim;a++)v+=f->left[a]*vi[a]+f->right[a]*vj[a];
 if(!isfinite(v))return 0;
 *rate=v;return 1;
}
/* Integral of |(w-vi).n| over a face simplex: one point, a segment, or a
 * triangle in 1D/2D/3D. Vertices are owner-relative, row-major dim*dim;
 * measure is supplied by the actual geometry. A polygonal 3D face sums
 * its triangle contributions. vi/vj are GENERATOR velocities (including
 * the coordinate Jacobian for an image), not reflected fluid traces.
 * Use velocity differences before multiplication: a uniform boost must
 * not introduce sweep by subtracting large laboratory-frame fluxes.
 * This is a time-resolution measure; it neither transports mass nor adds
 * dissipation. The caller supplies genuine face vertices and measure. */
/* Internal shared evaluation with a distance computed from this separation.
 * The segment wrapper already needs that distance to construct endpoints;
 * passing it here avoids repeating robust hypot calls on every face. */
GFS_POWER_MOTION_HD int gfs_material_power_face_sweep_distance(int dim,
 const double*separation,const double*vertices,double measure,double d,
 const double*vi,const double*vj,double wi_dot,double wj_dot,double*out){
 if(dim<1||dim>3||!separation||!vertices||!vi||!vj||!out||
    !(measure>0)||!isfinite(measure)||!(d>0)||!isfinite(d)||
    !isfinite(wi_dot)||!isfinite(wj_dot))return 0;
 double dv[3]={0},value[3]={0};
 for(int a=0;a<dim;a++){
  if(!isfinite(separation[a])||!isfinite(vi[a])||!isfinite(vj[a]))return 0;
  dv[a]=vj[a]-vi[a];
  if(!isfinite(dv[a]))return 0;
 }
 for(int k=0;k<dim;k++){
  double v=.5*(wi_dot-wj_dot);
  for(int a=0;a<dim;a++){
   double x=vertices[k*dim+a];if(!isfinite(x))return 0;
   v+=(separation[a]-x)*dv[a];
  }
  value[k]=v/d;if(!isfinite(value[k]))return 0;
 }
 double mean;if(!gfs_material_affine_abs_mean(dim,value,&mean))return 0;
 *out=measure*mean;return *out>=0&&isfinite(*out);
}
GFS_POWER_MOTION_HD int gfs_material_power_face_sweep(int dim,
 const double*separation,const double*vertices,double measure,
 const double*vi,const double*vj,double wi_dot,double wj_dot,double*out){
 if(dim<1||dim>3||!separation)return 0;
 double d=fabs(separation[0]);
 for(int a=1;a<dim;a++)d=hypot(d,separation[a]);
 return gfs_material_power_face_sweep_distance(dim,separation,vertices,measure,d,
    vi,vj,wi_dot,wj_dot,out);
}
/* A straight power edge is perpendicular to separation, so centroid and
 * length recover its endpoints without extending the existing 2D ABI. */
GFS_POWER_MOTION_HD int gfs_material_power_segment_sweep(
 const double*separation,const double*centroid,double length,
 const double*vi,const double*vj,double wi_dot,double wj_dot,double*out){
 if(!separation||!centroid||!(length>0)||!isfinite(length))return 0;
 double d=hypot(separation[0],separation[1]);if(!(d>0)||!isfinite(d))return 0;
 double tx=-.5*length*(separation[1]/d),ty=.5*length*(separation[0]/d);
 double v[4]={centroid[0]-tx,centroid[1]-ty,centroid[0]+tx,centroid[1]+ty};
 return gfs_material_power_face_sweep_distance(2,separation,v,length,d,
    vi,vj,wi_dot,wj_dot,out);
}
#undef GFS_POWER_MOTION_HD
#endif
