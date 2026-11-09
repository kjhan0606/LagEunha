#ifndef GFS_MATERIAL_FIXED_SOLID2D_H
#define GFS_MATERIAL_FIXED_SOLID2D_H
#include <math.h>
#include <float.h>
#include "gfs_power_predicate2d.h"
#ifdef __CUDACC__
#define GFS_SOLID_HD __host__ __device__
#else
#define GFS_SOLID_HD
#endif
/* P and S are convex CCW polygons. The material cell is P minus the fixed
 * solid S: one state, possibly several convex integration pieces. Only the
 * exposed boundary is emitted. Artificial piece interfaces never enter the
 * Riemann solver or the weight Jacobian. Coordinates have the same origin.
 * Negative tags encode stationary solid facets, never power neighbours. */
#define GFS_SOLID_TAG_BASE (-2000)
typedef struct { double ax,ay,bx,by; int neighbor; } GfsSolidSegment2D;
typedef struct { int edge,part,solid; } GfsSolidIterator2D;
typedef struct { int count,walls; double volume,cx,cy,nx,ny; } GfsSolidMeasure2D;
static GFS_SOLID_HD inline double gfs_solid_dd_value(GfsPowerDD a){return a.hi+a.lo;}
static GFS_SOLID_HD inline double gfs_solid_cross(double ax,double ay,double bx,double by){
 return gfs_solid_dd_value(gfs_pdd_sub(gfs_pdd_mul((GfsPowerDD){ax,0},(GfsPowerDD){by,0}),
                                gfs_pdd_mul((GfsPowerDD){ay,0},(GfsPowerDD){bx,0})));
}
/* Convex clipping of a line segment. Empty and point intersections have no
 * boundary measure. There is no geometric floor or adjustable tolerance. */
static GFS_SOLID_HD inline int gfs_solid_interval(int n,const double*x,const double*y,
 double ax,double ay,double bx,double by,double*lo,double*hi){
 double l=0,h=1;
 for(int k=0;k<n;k++){
  int j=(k+1)%n;double ex=x[j]-x[k],ey=y[j]-y[k];
  double a=gfs_solid_cross(ex,ey,ax-x[k],ay-y[k]);
  double b=gfs_solid_cross(ex,ey,bx-x[k],by-y[k]);
  if(a<0&&b<0)return 0;
  if(a<0||b<0){double t=a/(a-b);if(a<0)l=fmax(l,t);else h=fmin(h,t);}
  if(!(h>l))return 0;
 }
 *lo=l;*hi=h;return h>l;
}
/* Power-cell support normals are authoritative. Inferring a halfplane
 * from a rounded, almost zero-length edge can reverse its orientation. */
typedef struct {const double *nx,*ny,*h;} GfsSolidSupports2D;
static GFS_SOLID_HD inline double gfs_solid_plane_value(double nx,double ny,double h,double x,double y){
 GfsPowerDD z=gfs_pdd_sub((GfsPowerDD){h,0},gfs_pdd_add(gfs_pdd_mul((GfsPowerDD){nx,0},(GfsPowerDD){x,0}),gfs_pdd_mul((GfsPowerDD){ny,0},(GfsPowerDD){y,0})));
 return z.hi+z.lo;
}
static GFS_SOLID_HD inline int gfs_solid_support_interval(int n,const GfsSolidSupports2D*p,
 double ax,double ay,double bx,double by,double*lo,double*hi){
 double l=0,h=1;
 for(int k=0;k<n;k++){
  double a=gfs_solid_plane_value(p->nx[k],p->ny[k],p->h[k],ax,ay);
  double b=gfs_solid_plane_value(p->nx[k],p->ny[k],p->h[k],bx,by);
  if(a<0&&b<0)return 0;
  if(a<0||b<0){double t=a/(a-b);if(a<0)l=fmax(l,t);else h=fmin(h,t);}
  if(!(h>l))return 0;
 }
 *lo=l;*hi=h;return h>l;
}
static GFS_SOLID_HD inline double gfs_solid_lerp(double a,double b,double t){
 if(t==0)return a;if(t==1)return b;return fma(t,b-a,a);
}
static GFS_SOLID_HD inline int gfs_solid_next_supports(int np,const double*px,const double*py,const int*pn,
 int ns,const double*sx,const double*sy,GfsSolidIterator2D*it,GfsSolidSegment2D*out,const GfsSolidSupports2D*supports){
 while(!it->solid){
  if(it->edge==np){it->solid=1;it->edge=0;it->part=0;break;}
  int k=it->edge,j=(k+1)%np;double l,h,ta=0,tb=1;
  int cut=gfs_solid_interval(ns,sx,sy,px[k],py[k],px[j],py[j],&l,&h);
  if(!cut){it->edge++;it->part=0;}
  else if(it->part==0){ta=0;tb=l;it->part=1;}
  else {ta=h;tb=1;it->edge++;it->part=0;}
  if(!(tb>ta))continue;
  out->ax=gfs_solid_lerp(px[k],px[j],ta);out->ay=gfs_solid_lerp(py[k],py[j],ta);
  out->bx=gfs_solid_lerp(px[k],px[j],tb);out->by=gfs_solid_lerp(py[k],py[j],tb);
  if(out->ax==out->bx&&out->ay==out->by)continue;
  out->neighbor=pn?pn[k]:k;return 1;
 }
 while(it->edge<ns){
  int k=it->edge++,j=(k+1)%ns;double l,h;
  int clipped=supports?gfs_solid_support_interval(np,supports,sx[k],sy[k],sx[j],sy[j],&l,&h):gfs_solid_interval(np,px,py,sx[k],sy[k],sx[j],sy[j],&l,&h);
  if(!clipped)continue;
  double ax=gfs_solid_lerp(sx[k],sx[j],l),ay=gfs_solid_lerp(sy[k],sy[j],l);
  double bx=gfs_solid_lerp(sx[k],sx[j],h),by=gfs_solid_lerp(sy[k],sy[j],h);
  /* A coincident facet contributes only when the fluid side is in P.
   * Two polygons touching with disjoint interiors still have a wall; a
   * coincident facet bounding only solid has no fluid face. */
  int fluid=1;
  for(int e=0;e<np;e++){
   int f=(e+1)%np;double ex=px[f]-px[e],ey=py[f]-py[e];
   double a=gfs_solid_cross(ex,ey,ax-px[e],ay-py[e]);
   double b=gfs_solid_cross(ex,ey,bx-px[e],by-py[e]);
   double outward=-gfs_solid_cross(ex,ey,sy[j]-sy[k],sx[k]-sx[j]);
   if(supports){a=gfs_solid_plane_value(supports->nx[e],supports->ny[e],supports->h[e],ax,ay);b=gfs_solid_plane_value(supports->nx[e],supports->ny[e],supports->h[e],bx,by);outward=supports->nx[e]*(sy[j]-sy[k])+supports->ny[e]*(sx[k]-sx[j]);}
   if(a==0&&b==0&&outward>0){fluid=0;break;}
  }
  if(!fluid)continue;
  /* Reverse S's orientation: the fluid remains on the left. */
  out->ax=bx;out->ay=by;out->bx=ax;out->by=ay;
  if(out->ax==out->bx&&out->ay==out->by)continue;
  out->neighbor=GFS_SOLID_TAG_BASE-k;return 1;
 }
 return 0;
}
static GFS_SOLID_HD inline int gfs_solid_measure_supports(int np,const double*px,const double*py,const int*pn,
 int ns,const double*sx,const double*sy,GfsSolidMeasure2D*out,const GfsSolidSupports2D*supports){
 GfsSolidIterator2D it={0,0,0};GfsSolidSegment2D e;
 GfsPowerDD a={0,0},mx={0,0},my={0,0},nx={0,0},ny={0,0};int count=0,walls=0;
 if(np<3||ns<3)return 0;
 while(gfs_solid_next_supports(np,px,py,pn,ns,sx,sy,&it,&e,supports)){
  double z=gfs_solid_cross(e.ax,e.ay,e.bx,e.by);
  if(!isfinite(z))return 0;
  a=gfs_pdd_add(a,(GfsPowerDD){z,0});
  mx=gfs_pdd_add(mx,gfs_pdd_mul((GfsPowerDD){e.ax+e.bx,0},(GfsPowerDD){z,0}));
  my=gfs_pdd_add(my,gfs_pdd_mul((GfsPowerDD){e.ay+e.by,0},(GfsPowerDD){z,0}));
  nx=gfs_pdd_add(nx,(GfsPowerDD){e.by-e.ay,0});ny=gfs_pdd_add(ny,(GfsPowerDD){e.ax-e.bx,0});
  count++;walls+=e.neighbor<=GFS_SOLID_TAG_BASE;
 }
 double area=gfs_solid_dd_value(a);if(!(area>0)||!isfinite(area))return 0;
 out->count=count;out->walls=walls;out->volume=.5*area;
 out->cx=gfs_solid_dd_value(mx)/(3*area);out->cy=gfs_solid_dd_value(my)/(3*area);
 out->nx=gfs_solid_dd_value(nx);out->ny=gfs_solid_dd_value(ny);
 return isfinite(out->cx)&&isfinite(out->cy);
}
static GFS_SOLID_HD inline int gfs_solid_next(int np,const double*px,const double*py,const int*pn,
 int ns,const double*sx,const double*sy,GfsSolidIterator2D*it,GfsSolidSegment2D*out){
 return gfs_solid_next_supports(np,px,py,pn,ns,sx,sy,it,out,0);
}
static GFS_SOLID_HD inline int gfs_solid_measure(int np,const double*px,const double*py,const int*pn,
 int ns,const double*sx,const double*sy,GfsSolidMeasure2D*out){
 return gfs_solid_measure_supports(np,px,py,pn,ns,sx,sy,out,0);
}
#undef GFS_SOLID_HD
#endif
