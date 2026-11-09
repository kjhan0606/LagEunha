#ifndef GFS_MATERIAL_SWEPT_MEASURE_H
#define GFS_MATERIAL_SWEPT_MEASURE_H
#include <math.h>
#ifdef __CUDACC__
#define GFS_SWEPT_HD static __host__ __device__ inline
#else
#define GFS_SWEPT_HD static inline
#endif
/* Exact mean |affine f| on a face simplex: a point in 1D, segment in 2D,
 * triangle in 3D. Multiplying by its measure gives the unsigned sweep.
 * This is geometric time resolution, not viscosity or a material flux.
 * General planar 3D faces can sum the measured triangle contributions. */
GFS_SWEPT_HD int gfs_material_affine_abs_mean(int dim,const double*f,double*out){
 if(dim<1||dim>3||!f||!out)return 0;
 double scale=0;int positive=0,negative=0;
 for(int i=0;i<dim;i++){
  if(!isfinite(f[i]))return 0;
  scale=fmax(scale,fabs(f[i]));positive+=f[i]>0;negative+=f[i]<0;
 }
 if(scale==0){*out=0;return 1;}
 double v[3]={0,0,0};for(int i=0;i<dim;i++)v[i]=f[i]/scale;
 double mean=0;
 if(!positive||!negative){for(int i=0;i<dim;i++)mean+=fabs(v[i])/dim;}
 else if(dim==2){
  double a=fabs(v[0]),b=fabs(v[1]);mean=(a*a+b*b)/(2*(a+b));
 }else{
  /* Flip the sign if necessary so only the isolated corner is negative.
   * Its negative subtriangle has relative area a^2/((a+b)(a+c)). */
  int isolated=-1;
  for(int i=0;i<3;i++)if((positive==1&&v[i]>0)||(positive!=1&&v[i]<0))isolated=i;
  if(isolated<0)return 0;
  double a=fabs(v[isolated]),b=fabs(v[(isolated+1)%3]),c=fabs(v[(isolated+2)%3]);
  mean=(b+c-a)/3+(2*a/3)*(a/(a+b))*(a/(a+c));
 }
 *out=scale*mean;return *out>=0&&isfinite(*out);
}
#undef GFS_SWEPT_HD
#endif
