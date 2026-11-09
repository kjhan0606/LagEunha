#ifndef GFS_POWER_SEARCH_EXTENT_H
#define GFS_POWER_SEARCH_EXTENT_H
#include <math.h>
#include <limits.h>
/* Computational bins only. Deficits are lower[dim],upper[dim].
 * Preserve bin width and extend by whole bins; never alter physical sites.
 * Return -1 invalid/capacity overflow, 0 unchanged, 1 extended. */
static inline int gfs_power_expand_search_box(int dim,double cell,
 const double*deficit,double*lower,double*upper){
 if(dim<1||dim>3||!(cell>0)||!isfinite(cell)||!deficit||!lower||!upper)return -1;
 double lo[3],hi[3],count=1.;int changed=0;
 for(int d=0;d<dim;d++){
  if(!isfinite(lower[d])||!isfinite(upper[d])||!(upper[d]>lower[d])||
     !isfinite(deficit[d])||deficit[d]<0||
     !isfinite(deficit[dim+d])||deficit[dim+d]<0)return -1;
  lo[d]=lower[d];hi[d]=upper[d];
  if(deficit[d]>0){
   double bins=ceil(nextafter(deficit[d]/cell,INFINITY));
   lo[d]-=bins*cell;changed=1;
   /* A sub-ULP move at a large coordinate must not trigger an endless retry. */
   if(!(lo[d]<lower[d])||lower[d]-lo[d]<deficit[d])return -1;
  }
  if(deficit[dim+d]>0){
   double bins=ceil(nextafter(deficit[dim+d]/cell,INFINITY));
   hi[d]+=bins*cell;changed=1;
   if(!(hi[d]>upper[d])||hi[d]-upper[d]<deficit[dim+d])return -1;
  }
  double bins=ceil((hi[d]-lo[d])/cell);
  if(!isfinite(lo[d])||!isfinite(hi[d])||!(hi[d]>lo[d])||
     !(bins>=1)||!isfinite(bins)||bins>INT_MAX/count)return -1;
  count*=bins;
 }
 for(int d=0;d<dim;d++){lower[d]=lo[d];upper[d]=hi[d];}
 return changed;
}
#endif
