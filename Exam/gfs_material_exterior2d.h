#ifndef GFS_MATERIAL_EXTERIOR2D_H
#define GFS_MATERIAL_EXTERIOR2D_H
#include "gfs_material_trace_csr2d.h"
#include "gfs_material_exterior.h"
#include "gfs_material_face_response.h"
#include "gfs_power_coordinates2d.h"
#include <stdio.h>
#ifdef GFS_MATERIAL_FIXED_SOLID
#include "gfs_material_fixed_solid2d.h"
#include "gfs_material_fixed_cylinder_data.h"
#endif
/* Adapter for native FaceCSR arrays with explicit nonperiodic exterior sites.
 * Slots [0,n) are owned; [n,nt) are prescribed sites, never periodic aliases.
 * All coordinates below are physical coordinates; native corners are relative
 * to their owner. This is not the periodic/reflected-wall adapter.
 * Interior faces remain reciprocal. Exterior faces have one owned side.
 */
typedef struct {
 int n,nt,directed_count,pair_count;
 GfsTraceDirected2D *directed;
 GfsMaterialExteriorEdge *pair;
 int *pair_directed;
 double *center;
} GfsMaterialExterior2D;
static inline void gfs_material_exterior2d_free(GfsMaterialExterior2D*g){
 if(!g)return;
 free(g->directed);free(g->pair);free(g->pair_directed);free(g->center);
 memset(g,0,sizeof(*g));
}
static inline int gfs_material_exterior2d_build_coordinates(int n,int nt,int nf,
 const int*offset,const int*neighbor,const double*x,const double*y,const double*w,
 const double*c1x,const double*c1y,const double*c2x,const double*c2y,
 const double*volume,const double*lowx,const double*lowy,GfsMaterialExterior2D*out){
 if(!out)return 0;memset(out,0,sizeof(*out));
 if((lowx==NULL)!=(lowy==NULL))return 0;
 if(n<1||nt<=n||nf<1||!offset||!neighbor||!x||!y||!w||!volume||
    !c1x||!c1y||!c2x||!c2y||offset[0]!=0||offset[n]!=nf)return 0;
 GfsMaterialExterior2D g={0};int ok=0;
 const int diagnostics=getenv("GFS_MATERIAL_EXTERIOR_DIAGNOSTICS")!=NULL;
 int*off=(int*)calloc((size_t)n+1,sizeof(int));
 g.n=n;g.nt=nt;g.directed=(GfsTraceDirected2D*)malloc((size_t)nf*sizeof(*g.directed));
 g.pair=(GfsMaterialExteriorEdge*)malloc((size_t)nf*sizeof(*g.pair));
 g.pair_directed=(int*)malloc((size_t)nf*sizeof(int));g.center=(double*)calloc((size_t)2*n,sizeof(double));
 if(!off||!g.directed||!g.pair||!g.pair_directed||!g.center)goto done;
 for(int i=0;i<nt;i++)if(!isfinite(x[i])||!isfinite(y[i])||!isfinite(w[i])||
  (lowx&&(!isfinite(lowx[i])||!isfinite(lowy[i]))))goto done;
 for(int i=0;i<n;i++)if(offset[i]<0||offset[i+1]<offset[i]||offset[i+1]>nf)goto done;
 for(int i=0;i<n;i++){
  off[i]=g.directed_count;
  for(int k=offset[i];k<offset[i+1];k++){
   int j=neighbor[k];
#ifdef GFS_MATERIAL_FIXED_SOLID
   if(j<=GFS_SOLID_TAG_BASE){
    int side=GFS_SOLID_TAG_BASE-j;if(side>=GFS_FIXED_CYLINDER_N)goto done;
    int next=(side+1)%GFS_FIXED_CYLINDER_N;
    double dx=c2x[k]-c1x[k],dy=c2y[k]-c1y[k],A=hypot(dx,dy);
    if(!(A>0)||!isfinite(A))goto done;
    double sx=gfs_fixed_cylinder_x[side]-x[i],sy=gfs_fixed_cylinder_y[side]-y[i];
    double ex=gfs_fixed_cylinder_x[next]-gfs_fixed_cylinder_x[side],ey=gfs_fixed_cylinder_y[next]-gfs_fixed_cylinder_y[side];
    double tol=512*DBL_EPSILON*hypot(ex,ey)*(hypot(sx,sy)+hypot(c1x[k],c1y[k])+hypot(c2x[k],c2y[k]));
    if(fabs(gfs_solid_cross(ex,ey,c1x[k]-sx,c1y[k]-sy))>tol||fabs(gfs_solid_cross(ex,ey,c2x[k]-sx,c2y[k]-sy))>tol||dx*ex+dy*ey>0)goto done;
    g.directed[g.directed_count++]=(GfsTraceDirected2D){i,i,-ey,ex,.5*(c1x[k]+c2x[k]),.5*(c1y[k]+c2y[k]),A};
    continue;
   }
#endif
   if(j<0||j>=nt||j==i)goto done;
   double dx=gfs_power_coordinate_delta(x,lowx,i,j),dy=gfs_power_coordinate_delta(y,lowy,i,j),d=hypot(dx,dy);
   double A=hypot(c2x[k]-c1x[k],c2y[k]-c1y[k]);
   if(!(d>0)||!isfinite(d)||!isfinite(A)||!isfinite(c1x[k])||!isfinite(c1y[k])||!isfinite(c2x[k])||!isfinite(c2y[k]))goto done;
   /* Same resolved-face threshold as the native direct-FV path. No absolute
    * length floor and no new adjustable geometry tolerance. */
   if(A<256*DBL_EPSILON*d)continue;
   double cx=.5*(c1x[k]+c2x[k]),cy=.5*(c1y[k]+c2y[k]);
   /* Same paired-cutoff policy as gfs_material_trace_probe.inc. An
    * unresolved interior face must not survive on only one side because
    * its two independently rounded lengths straddle the identical cutoff.
    * Require the raw reverse face and its matching centroid; missing or
    * geometrically inconsistent faces are still rejected below. */
   double cut=256*DBL_EPSILON*d;
   if(j<n&&A<=2*cut){
    int omit=0;
    for(int h=offset[j];h<offset[j+1];h++){
     if(neighbor[h]!=i)continue;
     double other=hypot(c2x[h]-c1x[h],c2y[h]-c1y[h]);
     double bx=.5*(c1x[h]+c2x[h]),by=.5*(c1y[h]+c2y[h]);
     if(other<cut&&fabs(other-A)<=cut&&fabs(bx-(cx-dx))<=cut&&fabs(by-(cy-dy))<=cut)omit=1;
    }
    if(omit)continue;
   }
   /* Subtract the gauge-dependent weights first. Adding their common offset
    * to d^2 before cancellation destroys the small-separation plane value. */
   double plane=.5*(dx*dx+dy*dy+(w[i]-w[j]));
   /* Validate the computed endpoints, whose rounding scale does not shrink
    * when their midpoint cancels. Testing only that midpoint both understates
    * endpoint roundoff and can hide opposite endpoint errors. The epsilon
    * multiplier is unchanged; no absolute-coordinate or length floor. */
   double ex[2]={c1x[k],c2x[k]},ey[2]={c1y[k],c2y[k]};
   for(int endpoint=0;endpoint<2;endpoint++){
    double scale=fabs(dx*ex[endpoint])+fabs(dy*ey[endpoint])+.5*(dx*dx+dy*dy+fabs(w[i]-w[j]));
    double residual=fabs(dx*ex[endpoint]+dy*ey[endpoint]-plane);
    if(!isfinite(scale)||!isfinite(residual)||residual>256*DBL_EPSILON*scale){
     if(getenv("GFS_MATERIAL_EXTERIOR_DIAGNOSTICS")){
      fprintf(stderr,"[EXTERIOR_GEOMETRY] plane i=%d j=%d endpoint=%d error=%.17g tolerance=%.17g\n",i,j,endpoint,residual,256*DBL_EPSILON*scale);
      fprintf(stderr,"[EXTERIOR_PLANE_INPUT] x_i=%.17g y_i=%.17g x_j=%.17g y_j=%.17g w_i=%.17g w_j=%.17g c1x=%.17g c1y=%.17g c2x=%.17g c2y=%.17g d=%.17g area=%.17g\n",x[i],y[i],x[j],y[j],w[i],w[j],c1x[k],c1y[k],c2x[k],c2y[k],d,A);
     }
     goto done;
    }
   }
   g.directed[g.directed_count++]=(GfsTraceDirected2D){i,j,dx,dy,cx,cy,A};
  }
 }off[n]=g.directed_count;
 for(int k=0;k<g.directed_count;k++){
  const GfsTraceDirected2D*a=g.directed+k;double d=hypot(a->dx,a->dy);
#ifdef GFS_MATERIAL_FIXED_SOLID
  if(a->i==a->j){int z=g.pair_count++;g.pair[z]=(GfsMaterialExteriorEdge){a->i,a->j,0};g.pair_directed[z]=k;continue;}
#endif
  if(a->j<n){
   int reverse=0;
   for(int h=off[a->j];h<off[a->j+1];h++){
    const GfsTraceDirected2D*b=g.directed+h;if(b->j!=a->i)continue;
    double tol=256*DBL_EPSILON*d;
    double at=256*DBL_EPSILON*fmax(d,fmax(a->area,b->area));
    double sa[2],sb[2];
    if(!gfs_trace_endpoint_scales2d(a,d,sa)||!gfs_trace_endpoint_scales2d(b,d,sb))goto done;
    double xt=256*DBL_EPSILON*fmax(d,fmax(sa[0],sb[0]));
    double yt=256*DBL_EPSILON*fmax(d,fmax(sa[1],sb[1]));
    if(diagnostics&&
       (fabs(b->cx-(a->cx-a->dx))>256*DBL_EPSILON*fmax(d,fmax(fabs(a->cx),fabs(b->cx)))||
        fabs(b->cy-(a->cy-a->dy))>256*DBL_EPSILON*fmax(d,fmax(fabs(a->cy),fabs(b->cy)))))
     fprintf(stderr,"[EXTERIOR_CENTROID_SCALE] i=%d j=%d d=%.17g area_a=%.17g area_b=%.17g dx=%.17g dy=%.17g cx_a=%.17g cy_a=%.17g cx_b=%.17g cy_b=%.17g error_x=%.17g error_y=%.17g tolerance_x=%.17g tolerance_y=%.17g\n",
      a->i,a->j,d,a->area,b->area,a->dx,a->dy,a->cx,a->cy,b->cx,b->cy,
      fabs(b->cx-(a->cx-a->dx)),fabs(b->cy-(a->cy-a->dy)),xt,yt);
#ifdef GFS_MATERIAL_FIXED_SOLID
    if(fabs(b->dx+a->dx)>tol||fabs(b->dy+a->dy)>tol||fabs(b->cx-(a->cx-a->dx))>xt||fabs(b->cy-(a->cy-a->dy))>yt||fabs(b->area-a->area)>at)continue;
#endif
    if(fabs(b->dx+a->dx)>tol||fabs(b->dy+a->dy)>tol||fabs(b->cx-(a->cx-a->dx))>xt||fabs(b->cy-(a->cy-a->dy))>yt||fabs(b->area-a->area)>at){if(getenv("GFS_MATERIAL_EXTERIOR_DIAGNOSTICS"))fprintf(stderr,"[EXTERIOR_GEOMETRY] reverse i=%d j=%d dx=%.17g dy=%.17g cx=%.17g cy=%.17g area=%.17g tolerances=%.17g,%.17g,%.17g,%.17g\n",a->i,a->j,b->dx+a->dx,b->dy+a->dy,b->cx-(a->cx-a->dx),b->cy-(a->cy-a->dy),b->area-a->area,tol,xt,yt,at);goto done;}
    reverse++;
   }
   if(reverse!=1){if(getenv("GFS_MATERIAL_EXTERIOR_DIAGNOSTICS"))fprintf(stderr,"[EXTERIOR_GEOMETRY] missing_reverse i=%d j=%d count=%d area=%.17g\n",a->i,a->j,reverse,a->area);goto done;}
  }
  if(a->i<a->j){
   double delta[3]={a->dx,a->dy,0},centroid[3]={a->cx,a->cy,0},zero[3]={0},base,coefficient;
   if(!gfs_material_power_face_rate(2,delta,centroid,zero,zero,a->area,&base,&coefficient))goto done;
   int z=g.pair_count++;g.pair[z]=(GfsMaterialExteriorEdge){a->i,a->j,coefficient};g.pair_directed[z]=k;
  }
 }
 if(!gfs_material_volume_centroid2d(n,g.directed_count,g.directed,volume,g.center))goto done;
 *out=g;memset(&g,0,sizeof(g));ok=1;
 done:free(off);gfs_material_exterior2d_free(&g);return ok;
}
static inline int gfs_material_exterior2d_build(int n,int nt,int nf,
 const int*offset,const int*neighbor,const double*x,const double*y,const double*w,
 const double*c1x,const double*c1y,const double*c2x,const double*c2y,
 const double*volume,GfsMaterialExterior2D*out){
 return gfs_material_exterior2d_build_coordinates(n,nt,nf,offset,neighbor,x,y,w,
  c1x,c1y,c2x,c2y,volume,NULL,NULL,out);
}
/* Existing RHS layout Fx,Fy,Vdot,H. Only Vdot enters this geometric solve.
 * Exterior Wdot is prescribed; its contribution stays on the RHS. */
static inline int gfs_material_exterior2d_weight_rate(const GfsMaterialExterior2D*g,
 const GfsFvPrimitive*cell,const double*rhs4,const double*known,
 double*q,double*realized,int*iterations,double*relative){
 if(!g||g->n<1||!cell||!rhs4||!known||!q||!realized||!iterations||!relative)return 0;
 int n=g->n,ok=0;double*buf=(double*)calloc((size_t)3*n,sizeof(double));if(!buf)return 0;
 double*rhs=buf,*base=rhs+n,*external=base+n;
 for(int k=0;k<g->pair_count;k++){
  const GfsTraceDirected2D*e=g->directed+g->pair_directed[k];
  double d[3]={e->dx,e->dy,0},c[3]={e->cx,e->cy,0},b,coefficient;
  if(!gfs_material_power_face_rate(2,d,c,cell[e->i].u,cell[e->j].u,e->area,&b,&coefficient))goto done;
  base[e->i]+=b;
  if(e->j<n)base[e->j]-=b;
  else {if(!isfinite(known[e->j]))goto done;external[e->i]+=coefficient*known[e->j];}
 }
 for(int i=0;i<n;i++)rhs[i]=rhs4[4*i+2]-base[i]+external[i];
 if(!gfs_material_exterior_solve_selected(n,g->nt,g->pair_count,g->pair,rhs,q,iterations,relative))goto done;
 gfs_material_exterior_matvec(n,g->pair_count,g->pair,q,realized);
 for(int i=0;i<n;i++)realized[i]+=base[i]-external[i];
 ok=1;
 done:free(buf);return ok;
}
#endif
