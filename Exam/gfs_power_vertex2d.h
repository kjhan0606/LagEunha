#ifndef GFS_POWER_VERTEX2D_H
#define GFS_POWER_VERTEX2D_H
#include <math.h>
#include <float.h>
#ifdef __CUDACC__
#define GFS_VERTEX_HD static __host__ __device__ inline
#else
#define GFS_VERTEX_HD static inline
#endif
/* Exact power-plane intersection, evaluated from the best-conditioned of
 * the three equivalent anchors. Input: x,y,W for three distinct sites.
 * Output relative to input site0. This does not define or change topology,
 * apply weight clipping, or repair degenerate cells. Caller certifies planes.
 * Lexicographic order removes caller-order dependence in anchor tie breaks. */
GFS_VERTEX_HD int gfs_power_vertex2d(const double site[3][3],double out[2]){
 int order[3]={0,1,2};
 for(int i=0;i<3;i++)for(int a=0;a<3;a++)if(!isfinite(site[i][a]))return 0;
 for(int i=0;i<3;i++)for(int j=i+1;j<3;j++){
  int a=order[i],b=order[j];
  if(site[b][0]<site[a][0]||(site[b][0]==site[a][0]&&site[b][1]<site[a][1])){order[i]=b;order[j]=a;}
 }
 double best=-1,dx[2]={0},dy[2]={0},det=0;int anchor=-1,other[2]={0};
 for(int q=0;q<3;q++){
  int a=order[q],b=order[(q+1)%3],c=order[(q+2)%3];
  double x=site[b][0]-site[a][0],y=site[b][1]-site[a][1];
  double u=site[c][0]-site[a][0],v=site[c][1]-site[a][1];
  double z=fma(x,v,-y*u),den=hypot(x,y)*hypot(u,v);
  if(!(den>0)||!isfinite(den)||!isfinite(z))continue;
  double score=fabs(z)/den;
  if(score>best){best=score;anchor=a;other[0]=b;other[1]=c;dx[0]=x;dy[0]=y;dx[1]=u;dy[1]=v;det=z;}
 }
 if(anchor<0||!(best>64*DBL_EPSILON))return 0;
 double h[2];
 for(int k=0;k<2;k++)h[k]=.5*(fma(dx[k],dx[k],dy[k]*dy[k])+(site[anchor][2]-site[other[k]][2]));
 double x=fma(h[0],dy[1],-dy[0]*h[1])/det;
 double y=fma(dx[0],h[1],-h[0]*dx[1])/det;
 out[0]=x+(site[anchor][0]-site[0][0]);out[1]=y+(site[anchor][1]-site[0][1]);
 return isfinite(out[0])&&isfinite(out[1]);
}
#undef GFS_VERTEX_HD
#endif
