#ifndef GFS_MATERIAL_CENTROID2D_H
#define GFS_MATERIAL_CENTROID2D_H
#include "gfs_material_plane.h"
/* Included after GfsTraceDirected2D is declared. */
/* Exact polygon area/centroid from owner-relative straight-face midpoints.
 * Directed edges must cover each closed cell boundary exactly once. */
static inline int gfs_material_volume_centroid2d(int n,int nf,
 const GfsTraceDirected2D*edge,const double*expected_volume,double*center){
 if(n<1||nf<1||!edge||!expected_volume||!center)return 0;
 double *work=(double*)calloc((size_t)n*4,sizeof(double));if(!work)return 0;
 memset(center,0,(size_t)n*2*sizeof(double));int ok=0;
 for(int k=0;k<nf;k++){
  const GfsTraceDirected2D*e=edge+k;int i=e->i;double d=hypot(e->dx,e->dy);
  if(i<0||i>=n||!(d>0)||!(e->area>0)||!isfinite(d)||!isfinite(e->area)||!isfinite(e->cx)||!isfinite(e->cy))goto done;
  double nx=e->dx/d,ny=e->dy/d,cross=e->area*(e->cx*nx+e->cy*ny);
  work[4*i]+=cross;work[4*i+1]+=e->area*nx;work[4*i+2]+=e->area*ny;work[4*i+3]+=e->area;
  center[2*i]+=e->cx*cross;center[2*i+1]+=e->cy*cross;
 }
 for(int i=0;i<n;i++){
  double V=.5*work[4*i];
  if(!(V>0)||!isfinite(V)||!(expected_volume[i]>0)||!isfinite(expected_volume[i])||
   fabs(V-expected_volume[i])>1e-10*fmax(V,expected_volume[i])||
   hypot(work[4*i+1],work[4*i+2])>1e-10*work[4*i+3])goto done;
  center[2*i]/=3*V;center[2*i+1]/=3*V;
  if(!isfinite(center[2*i])||!isfinite(center[2*i+1]))goto done;
 }
 ok=1;
 done:free(work);return ok;
}
/* Fitting coordinates only: do not use these edges for physical normals,
 * volume derivatives, gravity, or weight constraints. */
static inline GfsTraceDirected2D *gfs_trace_centroid_coordinates2d(int n,int nf,
 const GfsTraceDirected2D *edge,const unsigned char *reflected){
 if(n<1||nf<1||!edge)return NULL;
 double *V=(double*)calloc(n,sizeof(double)),*c=(double*)calloc((size_t)n*2,sizeof(double));
 GfsTraceDirected2D *fit=(GfsTraceDirected2D*)malloc((size_t)nf*sizeof(*fit));
 if(!V||!c||!fit)goto fail;
 for(int k=0;k<nf;k++){const GfsTraceDirected2D *e=edge+k;
  if(e->i<0||e->i>=n||e->j<0||e->j>=n)goto fail;
  double d=hypot(e->dx,e->dy);if(!(d>0))goto fail;
  V[e->i]+=.5*e->area*(e->cx*e->dx+e->cy*e->dy)/d;
 }
 if(!gfs_material_volume_centroid2d(n,nf,edge,V,c))goto fail;
 for(int k=0;k<nf;k++){const GfsTraceDirected2D *e=edge+k;fit[k]=*e;
  double cx=c[2*e->j],cy=c[2*e->j+1];if(reflected&&reflected[k]){
#ifdef GFS_MATERIAL_DIRECT_WALLS
   double d[2]={e->dx,e->dy},v[2]={cx,cy};
   if(!gfs_material_plane_reflect(2,d,v))goto fail;
   cx=v[0];cy=v[1];
#else
   cy=-cy;
#endif
  }
  fit[k].dx+=cx-c[2*e->i];fit[k].dy+=cy-c[2*e->i+1];
  fit[k].cx-=c[2*e->i];fit[k].cy-=c[2*e->i+1];
 }free(V);free(c);return fit;
 fail:free(V);free(c);free(fit);return NULL;
}
#endif
