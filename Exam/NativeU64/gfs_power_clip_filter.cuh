#ifndef GFS_POWER_CLIP_FILTER_CUH
#define GFS_POWER_CLIP_FILTER_CUH
#include "../gfs_power_predicate2d.h"
#include "../gfs_power_chart2d.h"
#include "../gfs_power_coordinates2d.h"
/* Included after weight_geometry_value/fraction. Only opt-in raw, unbounded
 * power planes use this numerical predicate; legacy clamped planes do not. */
__device__ __forceinline__
int gpu_power_clip_side(int owner,int left,int right,int query,
 const double*px,const double*py,const double*W,
 double nx,double ny,double threshold,double*x,double*y,double*residual,
 const int*alias=NULL,double lx=0,double ly=0,const double*lowx=NULL,const double*lowy=NULL){
 if(!px||!py||!W||owner<0||query<0)return 0;
 /* A defining plane contains its symbolic vertex, including a vertex
  * shared with an artificial seed-box edge. Establish this identity before
  * excluding negative supports from the two-real-plane intersection solve.
  * Otherwise a repeated cut can discard its own seed intersection. */
 if(query==left||query==right){*residual=0;return 1;}
 if(left<0||right<0||left==right)return 0;
 double scale=fabs(nx*(*x))+fabs(ny*(*y))+fabs(threshold);
 if(fabs(*residual)>1024*DBL_EPSILON*scale)return 0;
 int id[4]={owner,left,right,query};double s[4][3];
 if(gfs_power_coordinates_active(4,id,lowx,lowy)){
  double v[2],r;int sign;
  if(!gfs_power_coordinates_side2d(id,px,py,W,lowx,lowy,v,&r,&sign))return -1;
  *x=v[0];*y=v[1];*residual=r;return 1;
 }
 for(int k=0;k<4;k++){
  s[k][0]=px[id[k]];s[k][1]=py[id[k]];s[k][2]=weight_geometry_value(W,id[k]);
  if(k){
   double dx=s[k][0]-s[0][0],dy=s[k][1]-s[0][1],d2=dx*dx+dy*dy;
   if(!(d2>0))return 0;
   double wf=.5+.5*(s[0][2]-s[k][2])/d2;
   if(!isfinite(wf)||weight_geometry_fraction(wf)!=wf)return 0;
  }
 }
 /* Preserve the reflected-support convention. Otherwise the side test and
  * final vertex share the original-owner chart, including the query plane. */
 if(alias&&alias[owner]>=0&&alias[left]>=0&&alias[right]>=0){
  double base[4][3],image[4][3],period[2]={lx,ly},low[4][2];
  for(int k=0;k<4;k++){
   for(int q=0;q<3;q++)base[k][q]=image[k][q]=s[k][q];
   int a=alias[id[k]];
   if(a>=0){base[k][0]=px[a];base[k][1]=py[a];}
  }
  if(!gfs_power_chart_sites2d(base,image,period,4,s,low))return -1;
  int rounded=0;for(int k=0;k<4;k++)for(int q=0;q<2;q++)rounded|=low[k][q]!=0;
  if(rounded){
   GfsPowerDD exact[4][3];
   for(int k=0;k<4;k++)for(int q=0;q<3;q++){
    exact[k][q].hi=s[k][q];exact[k][q].lo=q<2?low[k][q]:0;
   }
   double v[2],r;int sign;
   if(!gfs_power_vertex_side_dd2d(exact,v,&r,&sign))return -1;
   *x=v[0];*y=v[1];*residual=r;return 1;
  }
 }
 double v[2],r;int sign;
 if(!gfs_power_vertex_side2d(s,v,&r,&sign))return -1;
 *x=v[0];*y=v[1];
 *residual=r;
 return 1;
}
#endif
