#ifndef GFS_POWER_COORDINATES2D_H
#define GFS_POWER_COORDINATES2D_H
#include "gfs_power_predicate2d.h"
#ifdef __CUDACC__
#define GFS_COORD_HD static __host__ __device__ inline
#else
#define GFS_COORD_HD static inline
#endif
/* Optional low coordinate parts retain local image displacements when an
 * absolute coordinate cannot represent them accurately. They change neither
 * the intended site nor its weight. Ordinary sites use a zero low part. */
GFS_COORD_HD GfsPowerDD gfs_power_coordinate(double high,const double*low,int i){
 GfsPowerDD a={high,low?low[i]:0};return a;
}
GFS_COORD_HD double gfs_power_coordinate_delta(const double*high,
 const double*low,int i,int j){
 if(!low||(low[i]==0&&low[j]==0))return high[j]-high[i];
 GfsPowerDD d=gfs_pdd_sub(gfs_power_coordinate(high[j],low,j),gfs_power_coordinate(high[i],low,i));
 return d.hi+d.lo;
}
GFS_COORD_HD int gfs_power_coordinates_active(int count,const int*id,
 const double*lowx,const double*lowy){
 if(!lowx||!lowy)return 0;
 for(int k=0;k<count;k++)if(lowx[id[k]]!=0||lowy[id[k]]!=0)return 1;
 return 0;
}
GFS_COORD_HD int gfs_power_coordinates_side2d(const int id[4],const double*x,
 const double*y,const double*w,const double*lowx,const double*lowy,
 double vertex[2],double*residual,int*sign){
 GfsPowerDD s[4][3];
 for(int k=0;k<4;k++){
  int j=id[k];s[k][0]=gfs_power_coordinate(x[j],lowx,j);
  s[k][1]=gfs_power_coordinate(y[j],lowy,j);s[k][2]=gfs_pdd_value(w[j]);
 }
 return gfs_power_vertex_side_dd2d(s,vertex,residual,sign);
}
GFS_COORD_HD int gfs_power_coordinates_vertex2d(const int id[3],const double*x,
 const double*y,const double*w,const double*lowx,const double*lowy,double vertex[2]){
 int q[4]={id[0],id[1],id[2],id[1]},sign;double residual;
 return gfs_power_coordinates_side2d(q,x,y,w,lowx,lowy,vertex,&residual,&sign);
}
#undef GFS_COORD_HD
#endif
