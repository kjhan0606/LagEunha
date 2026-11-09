#ifndef GFS_POWER_CHART2D_H
#define GFS_POWER_CHART2D_H
#include "gfs_power_vertex2d.h"
#include "gfs_power_predicate2d.h"
#ifdef __CUDACC__
#define GFS_CHART_HD static __host__ __device__ inline
#else
#define GFS_CHART_HD static inline
#endif
/* The first three sites define a vertex. A fourth is an optional query.
 * Select the chart from the defining sites only, so a side test uses exactly
 * the chart used by the final vertex. Restore original owner coordinates
 * before translating them: translated image slots may already have rounded.
 * A zero period is a physical, nonperiodic coordinate. */
GFS_CHART_HD int gfs_power_chart_sites2d(const double base[][3],
 const double image[][3],const double period[2],int count,double site[][3],double low[][2]){
 if(count<3||count>4)return 0;
 for(int i=0;i<count;i++)for(int q=0;q<3;q++){
  if(!isfinite(base[i][q])||!isfinite(image[i][q]))return 0;
  site[i][q]=base[i][q];
  if(q<2)low[i][q]=0;
 }
 for(int q=0;q<2;q++){
  double L=period[q];if(!isfinite(L)||L<0)return 0;
  if(L==0){for(int i=0;i<count;i++)site[i][q]=image[i][q];continue;}
  double lo=image[0][q],hi=lo;
  for(int i=1;i<3;i++){lo=fmin(lo,image[i][q]);hi=fmax(hi,image[i][q]);}
  double origin=floor((lo/L+hi/L)*.5+.5);
  for(int i=0;i<count;i++){
   double shift=nearbyint((image[i][q]-base[i][q])/L)-origin;
   site[i][q]=fma(shift,L,base[i][q]);
   if(!isfinite(site[i][q]))return 0;
   /* A chart translation must not change the physical power planes. Near
    * a binade boundary, even subtracting one exact period can lose bits of
    * the original coordinate. Retain that low part; most shifts are exact
    * and keep the old double path bit for bit. */
   if(shift!=0){
    GfsPowerDD translated=gfs_pdd_add(gfs_pdd_mul(gfs_pdd_value(shift),
      gfs_pdd_value(L)),gfs_pdd_value(base[i][q]));
    GfsPowerDD error=gfs_pdd_sub(translated,gfs_pdd_value(site[i][q]));
    low[i][q]=error.hi+error.lo;
    if(!isfinite(low[i][q]))return 0;
   }
  }
 }
 return 1;
}
GFS_CHART_HD int gfs_power_chart_vertex2d(const double base[3][3],
 const double image[3][3],const double period[2],double out[2]){
 double site[3][3],low[3][2];
 if(!gfs_power_chart_sites2d(base,image,period,3,site,low)||
    !gfs_power_vertex2d(site,out))return 0; /* Preserve the conditioning gate. */
 int rounded=0;for(int i=0;i<3;i++)for(int q=0;q<2;q++)rounded|=low[i][q]!=0;
 if(!rounded)return 1;
 GfsPowerDD exact[4][3];
 for(int i=0;i<4;i++)for(int q=0;q<3;q++){
  int j=i==3?1:i;exact[i][q].hi=site[j][q];exact[i][q].lo=q<2?low[j][q]:0;
 }
 double residual;int sign;
 return gfs_power_vertex_side_dd2d(exact,out,&residual,&sign);
}
#undef GFS_CHART_HD
#endif
