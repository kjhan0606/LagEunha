#ifndef GFS_MATERIAL_ORDERED_JUMP_H
#define GFS_MATERIAL_ORDERED_JUMP_H
#include <math.h>
#ifdef __CUDACC__
#define GFS_JUMP_HD static __host__ __device__ inline
#else
#define GFS_JUMP_HD static inline
#endif
/* Scalar face-normal increments. Geometry projects velocity before this call;
 * the same ordering rule therefore applies in every spatial dimension.
 * No physical parameter, velocity scale, or dimensional tolerance is added. */
GFS_JUMP_HD int gfs_material_ordered_jump_alpha(double mean_jump,double change,double *alpha){
 if(!alpha||!isfinite(mean_jump)||!isfinite(change))return 0;
 double a=1;
 if(change!=0){if(mean_jump==0||((mean_jump>0)==(change>0)))a=0;else a=fmin(1.0,-mean_jump/change);}
 *alpha=a;return 1;
}
#undef GFS_JUMP_HD
#endif
