#ifndef GFS_MATERIAL_CYLINDER_SUPPORT_H
#define GFS_MATERIAL_CYLINDER_SUPPORT_H
#include <math.h>
#include <float.h>
#ifdef __CUDACC__
#define GFS_CYL_SUPPORT_HD static __host__ __device__ inline
#else
#define GFS_CYL_SUPPORT_HD static inline
#endif
/* Equal-weight self reflection imposes n.(X-center)>=R. Certifying this
 * linear inequality at every vertex certifies the entire convex cell,
 * including edges which could cross the disk with both endpoints outside.
 * dx,dy locate the generator relative to the cylinder; vx,vy are the
 * local vertex coordinates. This is a support certificate, not a force,
 * clipping operation, physical length tolerance or radius-only check. */
GFS_CYL_SUPPORT_HD int gfs_material_cylinder_vertex_support(
 double dx,double dy,double R,double vx,double vy){
 double r=hypot(dx,dy);
 if(!(R>0)||!(r>0)||!isfinite(R)||!isfinite(r)||
    !isfinite(vx)||!isfinite(vy))return -1;
 double nx=dx/r,ny=dy/r;
 double residual=(r-R)+nx*vx+ny*vy;
 double scale=r+R+fabs(nx*vx)+fabs(ny*vy);
 return residual>=-64*DBL_EPSILON*scale?1:0;
}
/* At r=2R the near-side radial image reaches the obstacle centre. Its
 * removal is continuous only if that limiting site is inactive. The power
 * gap to owner i at X=x_i+vertex is |offset|^2+2 offset.vertex+Wi-Wmax.
 * Wmax bounds every reflected owner's weight. Testing all convex-cell
 * vertices proves inactivity throughout the fluid partition. Strict
 * roundoff separation is required; uncertain/active cuts fail visibly.
 * The identity is dimension-independent (1--3); it adds no force or scale. */
GFS_CYL_SUPPORT_HD int gfs_material_center_image_inactive(int dim,
 const double*offset,const double*vertex,double Wi,double Wmax){
 if(dim<1||dim>3||!offset||!vertex||!isfinite(Wi)||!isfinite(Wmax)||Wi>Wmax)return -1;
 double gap=Wi-Wmax,scale=fabs(Wi)+fabs(Wmax);
 for(int k=0;k<dim;k++){
  if(!isfinite(offset[k])||!isfinite(vertex[k]))return -1;
  double square=offset[k]*offset[k],cross=2*offset[k]*vertex[k];
  gap+=square+cross;scale+=square+fabs(cross);
 }
 if(!isfinite(gap)||!isfinite(scale))return -1;
 return gap>64*DBL_EPSILON*scale?1:0;
}
#undef GFS_CYL_SUPPORT_HD
#endif
