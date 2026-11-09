#ifndef GFS_MATERIAL_PLANE_H
#define GFS_MATERIAL_PLANE_H
#include <math.h>
#include <float.h>
#ifdef __CUDACC__
#define GFS_PLANE_HD static __host__ __device__ inline
#else
#define GFS_PLANE_HD static inline
#endif
/* Geometry of a generator and its equal-weight reflection in a fixed plane.
 * d is image minus generator; c is any point on that plane relative to the
 * generator. The same formulas apply in 1, 2 and 3 dimensions. */
GFS_PLANE_HD int gfs_material_plane_bisector(int dim,const double*d,const double*c){
 if(dim<1||dim>3||!d||!c)return 0;
 double length=0;int axis=-1,nonzero=0;
 for(int a=0;a<dim;a++){
  if(!isfinite(d[a])||!isfinite(c[a]))return 0;
  length=hypot(length,d[a]);if(d[a]!=0){axis=a;nonzero++;}
 }
 if(!(length>0)||!isfinite(length))return 0;
 /* Preserve the arithmetic of the existing axis-aligned wall. */
 if(nonzero==1)return fabs(c[axis]-.5*d[axis])<=256*DBL_EPSILON*length;
 double normal_offset=0;
 for(int a=0;a<dim;a++)normal_offset+=c[a]*(d[a]/length);
 return isfinite(normal_offset)&&fabs(normal_offset-.5*length)<=256*DBL_EPSILON*length;
}
GFS_PLANE_HD int gfs_material_plane_reflect(int dim,const double*d,double*v){
 if(dim<1||dim>3||!d||!v)return 0;
 double length=0;int axis=-1,nonzero=0;
 for(int a=0;a<dim;a++){
  if(!isfinite(d[a])||!isfinite(v[a]))return 0;
  length=hypot(length,d[a]);if(d[a]!=0){axis=a;nonzero++;}
 }
 if(!(length>0)||!isfinite(length))return 0;
 if(nonzero==1){v[axis]=-v[axis];return 1;}
 double vn=0;for(int a=0;a<dim;a++)vn+=v[a]*(d[a]/length);
 for(int a=0;a<dim;a++)v[a]-=2*vn*(d[a]/length);
 for(int a=0;a<dim;a++)if(!isfinite(v[a]))return 0;
 return 1;
}
#undef GFS_PLANE_HD
#endif
