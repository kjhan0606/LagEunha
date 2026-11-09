#ifndef GFS_MATERIAL_PRESSURE_OFFSET_H
#define GFS_MATERIAL_PRESSURE_OFFSET_H
#include <math.h>
#ifdef __CUDACC__
#define GFS_PRESSURE_HD static __host__ __device__ inline
#else
#define GFS_PRESSURE_HD static inline
#endif
/* Pressure continuity at a shared face, assuming locally common grad(P)/rho:
 * Pj-Pi = grad(P)i dot [ri - (rhoj/rhoi) rj], rj=ri-dij.
 * No gravity value, timestep, sensor, or fitted physical coefficient enters.
 * For equal density this is the original centroid separation, exactly.
 * This is a reconstruction hypothesis, not a new force or constitutive law.
 */
GFS_PRESSURE_HD int gfs_material_pressure_offset(int dim,double rhoi,double rhoj,
 const double*ri,const double*dij,double*out){
 if(dim<1||dim>3||!(rhoi>0)||!(rhoj>0)||!isfinite(rhoi)||!isfinite(rhoj)||!ri||!dij||!out)return 0;
 double ratio=rhoj/rhoi;if(!isfinite(ratio))return 0;
 for(int a=0;a<dim;a++){
  if(!isfinite(ri[a])||!isfinite(dij[a]))return 0;
  out[a]=rhoi==rhoj?dij[a]:ri[a]-ratio*(ri[a]-dij[a]);
  if(!isfinite(out[a]))return 0;
 }
 return 1;
}
#undef GFS_PRESSURE_HD
#endif
