#ifndef GFS_MATERIAL_LIMIT_MAPS2D_H
#define GFS_MATERIAL_LIMIT_MAPS2D_H
#include "gfs_material_rhs2d.h"
#include "gfs_material_reconstruction2d.h"
/* Fit with exactly the metric used to construct the CSR face traces. */
static inline int gfs_material_trace_fit2d(int n,const GfsTraceDirected2D*edge,
 const double*dq,double g[2],double*condition){
 double xx=0,xy=0,yy=0,bx=0,by=0;
 for(int k=0;k<n;k++){
  const GfsTraceDirected2D*e=edge+k;double d2=e->dx*e->dx+e->dy*e->dy;
  if(!(d2>0)||!isfinite(dq[k]))return 0;
  double w=gfs_trace_metric2d(e,d2);if(!(w>0)||!isfinite(w))return 0;
  xx+=w*e->dx*e->dx/d2;xy+=w*e->dx*e->dy/d2;yy+=w*e->dy*e->dy/d2;
  bx+=w*e->dx*dq[k]/d2;by+=w*e->dy*dq[k]/d2;
 }
 *condition=gfs_velocity_matrix_condition(xx,xy,yy);
 if(!isfinite(*condition)||*condition*DBL_EPSILON>1e-6)return 0;
 double det=xx*yy-xy*xy;g[0]=(yy*bx-xy*by)/det;g[1]=(xx*by-xy*bx)/det;
 return isfinite(g[0])&&isfinite(g[1]);
}
/* Candidate coefficient selection, no state clipping or physical parameter.
 * Uses the validated owner-sorted geometry that built t. Scalar extrema are
 * checked at actual face traces; velocity uses the experimental affine-preserving
 * geometric-direction midpoint limiter. Return 0 invalidates the entire alpha output. */
static inline int gfs_material_limit_maps_coordinates2d(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const GfsFvPrimitive*c,double*alpha3,const unsigned char*reflected_y){
 if(!t||!edge||!c||!alpha3||nd<1||t->cells<1)return 0;
 int n=t->cells,ok=0,maxdegree=0;int*off=calloc(n+1,sizeof(int));
 double*bounds=malloc((size_t)n*4*sizeof(double)),*work=NULL;
 if(!off||!bounds)goto done;
 for(int i=0;i<n;i++){
  if(!(c[i].rho>0)||!(c[i].p>0)||!isfinite(c[i].rho)||!isfinite(c[i].p)||!isfinite(c[i].u[0])||!isfinite(c[i].u[1]))goto done;
  bounds[4*i]=bounds[4*i+1]=c[i].rho;bounds[4*i+2]=bounds[4*i+3]=c[i].p;
  alpha3[3*i]=alpha3[3*i+1]=alpha3[3*i+2]=1;
 }
 for(int k=0;k<nd;k++){
  int i=edge[k].i,j=edge[k].j;
  int wall=reflected_y?reflected_y[k]:0;
  if(wall>1||i<0||i>=n||j<0||j>=n||(!wall&&i==j)||(k&&edge[k-1].i>i))goto done;
  if(wall&&!gfs_material_stationary_wall2d(edge+k))goto done;
  off[i+1]++;bounds[4*i]=fmin(bounds[4*i],c[j].rho);bounds[4*i+1]=fmax(bounds[4*i+1],c[j].rho);
  bounds[4*i+2]=fmin(bounds[4*i+2],c[j].p);bounds[4*i+3]=fmax(bounds[4*i+3],c[j].p);
 }
 for(int i=0;i<n;i++){if(off[i+1]>maxdegree)maxdegree=off[i+1];off[i+1]+=off[i];}
#ifndef GFS_MATERIAL_THERMO_LIMIT_ONLY
 work=malloc((size_t)maxdegree*4*sizeof(double));if(!work)goto done;
 for(int i=0;i<n;i++){
  int m=off[i+1]-off[i];double*dx=work,*dy=dx+maxdegree,*dvx=dy+maxdegree,*dvy=dvx+maxdegree;
  double scale=hypot(c[i].u[0],c[i].u[1]);
  for(int j=0;j<m;j++){const GfsTraceDirected2D*e=edge+off[i]+j;
   double v[2]={c[e->j].u[0],c[e->j].u[1]};
   if(reflected_y&&reflected_y[off[i]+j]){
#ifdef GFS_MATERIAL_DIRECT_WALLS
    double d[2]={e->dx,e->dy};if(!gfs_material_plane_reflect(2,d,v))goto done;
#else
    v[1]=-v[1];
#endif
   }
   dx[j]=e->dx;dy[j]=e->dy;dvx[j]=v[0]-c[i].u[0];dvy[j]=v[1]-c[i].u[1];scale=fmax(scale,hypot(c[e->j].u[0],c[e->j].u[1]));}
  double gx[2],gy[2],condition,c2;
  if(!gfs_material_trace_fit2d(m,edge+off[i],dvx,gx,&condition)||!gfs_material_trace_fit2d(m,edge+off[i],dvy,gy,&c2))goto done;
  double a=gfs_velocity_geometric_alpha(m,dx,dy,dvx,dvy,gx[0],gx[1],gy[0],gy[1],condition,scale);
  if(!(a>=0)||!isfinite(a))goto done;alpha3[3*i+2]=a;
 }
#endif
 for(int k=0;k<t->faces;k++){
  const GfsTraceFace2D*f=t->face+k;int owner[2]={f->i,f->j},start[2]={f->sl,f->sr},len[2]={f->nl,f->nr};
  for(int s=0;s<2;s++){
   int i=owner[s];if(i<0||i>=n||start[s]<0||len[s]<1||start[s]+len[s]>t->entries)goto done;
   for(int a=0;a<2;a++){
    double center=a?c[i].p:c[i].rho,delta=0;
    for(int j=start[s];j<start[s]+len[s];j++){
     const GfsTraceEntry*e=t->entry+j;if(e->cell<0||e->cell>=n||!isfinite(e->weight))goto done;
     delta+=e->weight*((a?c[e->cell].p:c[e->cell].rho)-center);
    }
    if(!isfinite(delta))goto done;double zero=0;
    double lower=bounds[4*i+2*a],upper=bounds[4*i+2*a+1];
#ifdef GFS_MATERIAL_DIRECT_WALLS
    /* No exterior scalar datum at a slip wall. Preserve positivity while
     * allowing one-sided affine extrapolation beyond interior extrema. */
    if(f->i==f->j){lower=0;upper=DBL_MAX;}
#endif
    double bound=gfs_material_scalar_alpha2d(1,&delta,&zero,1,0,center,lower,upper);
    alpha3[3*i+a]=fmin(alpha3[3*i+a],bound);
   }
  }
 }
 ok=1;
 done:free(off);free(bounds);free(work);return ok;
}
static inline int gfs_material_limit_maps_boundaries2d(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const GfsFvPrimitive*c,double*alpha3,const unsigned char*reflected){
#ifdef GFS_MATERIAL_CENTROID_TRACE
 if(!t)return 0;GfsTraceDirected2D *fit=gfs_trace_centroid_coordinates2d(t->cells,nd,edge,reflected);
 if(!fit)return 0;int ok=gfs_material_limit_maps_coordinates2d(t,nd,fit,c,alpha3,reflected);free(fit);return ok;
#else
 return gfs_material_limit_maps_coordinates2d(t,nd,edge,c,alpha3,reflected);
#endif
}
/* Preserve the periodic entrypoint and its arithmetic. Wall scalar ghosts
 * are even; the normal velocity is odd. No hydrostatic pressure offset. */
static inline int gfs_material_limit_maps2d(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const GfsFvPrimitive*c,double*alpha3){
 return gfs_material_limit_maps_boundaries2d(t,nd,edge,c,alpha3,NULL);
}
#endif
