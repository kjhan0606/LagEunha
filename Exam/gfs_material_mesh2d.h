#ifndef GFS_MATERIAL_MESH2D_H
#define GFS_MATERIAL_MESH2D_H
#include <math.h>
#include "gfs_face_work.h"
#ifdef __CUDACC__
#define GFS_MM_HD __host__ __device__
#else
#define GFS_MM_HD
#endif
/* Cell connectivity indexes shared GLOBAL vertices, never owner-local copies.
 * The caller supplies one velocity per vertex and immutable cell masses.
 * These kinematics do not select a hydrodynamic vertex velocity or perform a
 * remap. The caller must reject a failed trial before committing any vertex.
 */
static GFS_MM_HD inline double gfs_mm_cross(double ax,double ay,double bx,double by){return ax*by-ay*bx;}
static GFS_MM_HD inline int gfs_mm_convex_step(int n,const int *ids,
 const double (*x)[2],const double (*v)[2],double dt){
 if(n<3 || !(dt>=0) || !isfinite(dt))return 0;
 for(int k=0;k<n;k++){
  int i=ids[k],j=ids[(k+1)%n];
  double ax=x[j][0]-x[i][0],ay=x[j][1]-x[i][1];
  double ux=v[j][0]-v[i][0],uy=v[j][1]-v[i][1];
  for(int l=0;l<n;l++){
   if(l==k || l==(k+1)%n)continue;
   int q=ids[l];
   double bx=x[q][0]-x[i][0],by=x[q][1]-x[i][1];
   double wx=v[q][0]-v[i][0],wy=v[q][1]-v[i][1];
   double a=gfs_mm_cross(ux,uy,wx,wy);
   double b=gfs_mm_cross(ux,uy,bx,by)+gfs_mm_cross(ax,ay,wx,wy);
   double c=gfs_mm_cross(ax,ay,bx,by),end=(a*dt+b)*dt+c;
   if(!isfinite(c)||!isfinite(end)||c<=0||end<=0)return 0;
   if(a>0){double t=-b/(2*a);if(t>0 && t<dt && (a*t+b)*t+c<=0)return 0;}
  }
 }
 return 1;
}
static GFS_MM_HD inline double gfs_mm_area(int n,const int *ids,const double (*x)[2]){
 double a=0;int o=ids[0];
 for(int k=1;k<n-1;k++){
  int i=ids[k],j=ids[k+1];
  a+=gfs_mm_cross(x[i][0]-x[o][0],x[i][1]-x[o][1],x[j][0]-x[o][0],x[j][1]-x[o][1]);
 }
 return .5*a;
}
/* Exact area sweep for straight-line shared-vertex motion over one step.
 * Use midpoint edge geometry for both volume and traction work. Pressure
 * temporal accuracy is the responsibility of the hydrodynamic integrator.
 */
static GFS_MM_HD inline double gfs_mm_edge_sweep(int i,int j,
 const double (*x)[2],const double (*v)[2],double dt){
 double sx=x[j][1]-x[i][1]+.5*dt*(v[j][1]-v[i][1]);
 double sy=x[i][0]-x[j][0]+.5*dt*(v[i][0]-v[j][0]);
 return dt*gfs_face_sweep2d(sx,sy,.5*(v[i][0]+v[j][0]),.5*(v[i][1]+v[j][1]));
}
#undef GFS_MM_HD
#endif
