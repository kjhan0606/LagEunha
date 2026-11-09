#ifndef GFS_MATERIAL_IMAGE_TRACE2D_H
#define GFS_MATERIAL_IMAGE_TRACE2D_H
#include "gfs_material_exterior_trace2d.h"
#include "gfs_material_image_weight2d.h"
#if !defined(GFS_MATERIAL_DIRECT_WALLS) || !defined(GFS_MATERIAL_GRADIENT_PROJECT)
#error Image trace requires the existing direct-wall and scalar-reconstruction method
#endif
/* Mixture of interior faces, physical image walls, and prescribed exterior
 * samples. Image geometry and mesh velocity remain in g, unchanged. A wall
 * has only an interior fluid trace; its ghost is not a physical state used
 * to estimate gradients. Physical sample/face slots use the existing exterior
 * trace representation, so both existing C and CUDA RHS consume this object.
 * Periodic aliases must have been resolved before calling this adapter.
 */
typedef struct {
 GfsMaterialExteriorTrace2D fluid;
 GfsTraceDirected2D *physical;
 unsigned char *wall;
} GfsMaterialImageTrace2D;
static inline void gfs_material_image_trace2d_free(GfsMaterialImageTrace2D*t){
 if(!t)return;gfs_material_exterior_trace2d_free(&t->fluid);
 free(t->physical);free(t->wall);memset(t,0,sizeof(*t));
}
static inline int gfs_material_image_trace2d_build(const GfsMaterialExterior2D*g,
 const int*image_owner,const GfsFvPrimitive*sample,const double*sample_offset,
 const GfsFvPrimitive*exterior_face,GfsMaterialImageTrace2D*out){
 if(!out)return 0;memset(out,0,sizeof(*out));
 if(!g||g->n<1||g->nt<=g->n||!image_owner||!sample||!sample_offset||!exterior_face)return 0;
 GfsMaterialImageTrace2D a={0};GfsMaterialExteriorTrace2D*t=&a.fluid;int ok=0,nb=0;
 for(int j=g->n;j<g->nt;j++)if(image_owner[j]<-1||image_owner[j]>=g->n)return 0;
 for(int k=0;k<g->pair_count;k++){int j=g->pair[k].j;if(j>=g->n&&image_owner[j]<0)nb++;}
 if(g->nt>INT_MAX-nb)return 0;
 a.physical=malloc((size_t)g->directed_count*sizeof(*a.physical));a.wall=calloc(g->directed_count,1);
 t->fit=malloc((size_t)g->directed_count*sizeof(*t->fit));t->boundary_slot=malloc((size_t)g->directed_count*sizeof(int));
 t->state=calloc((size_t)g->nt+nb,sizeof(*t->state));
 if(!a.physical||!a.wall||!t->fit||!t->boundary_slot||!t->state)goto done;
 memcpy(t->state,sample,(size_t)g->nt*sizeof(*sample));
 for(int k=0;k<g->directed_count;k++){
  const GfsTraceDirected2D*e=g->directed+k;int i=e->i,j=e->j;
  double cx=j<g->n?g->center[2*j]:sample_offset[2*j],cy=j<g->n?g->center[2*j+1]:sample_offset[2*j+1];
  if(!isfinite(cx)||!isfinite(cy))goto done;
  a.physical[k]=*e;t->fit[k]=*e;t->boundary_slot[k]=-1;
  t->fit[k].dx+=cx-g->center[2*i];t->fit[k].dy+=cy-g->center[2*i+1];
  t->fit[k].cx-=g->center[2*i];t->fit[k].cy-=g->center[2*i+1];
  if((j>=g->n&&image_owner[j]>=0)
#ifdef GFS_MATERIAL_FIXED_SOLID
      || i==j
#endif
  ){
   a.wall[k]=1;a.physical[k].j=i;t->fit[k].j=i;
  }
 }
 int slot=g->nt;
 for(int k=0;k<g->pair_count;k++){
  int z=g->pair_directed[k],j=g->pair[k].j;
  if(j<g->n||image_owner[j]>=0)continue;
  GfsFvPrimitive s=exterior_face[k];
  if(!(s.rho>0)||!(s.p>0)||!isfinite(s.rho)||!isfinite(s.p)||!isfinite(s.u[0])||!isfinite(s.u[1]))goto done;
  t->boundary_slot[z]=slot;t->state[slot++]=s;
 }
 if(slot!=g->nt+nb||!gfs_trace_csr2d_build_coordinates_domain_planes(g->n,g->nt,slot,g->directed_count,a.physical,&t->trace,a.wall,t->fit,t->boundary_slot,1))goto done;
 if(t->trace.faces!=g->pair_count)goto done;
 for(int k=0;k<g->pair_count;k++){
  int z=g->pair_directed[k],j=g->pair[k].j;GfsTraceFace2D*f=t->trace.face+k;
  int expected=a.wall[z]?g->pair[k].i:j<g->n?j:t->boundary_slot[z];
  if(f->i!=g->pair[k].i||f->j!=expected)goto done;
 }
 *out=a;memset(&a,0,sizeof(a));ok=1;
 done:gfs_material_image_trace2d_free(&a);return ok;
}
/* Scalar projection is independent of the wall-normal velocity. The image
 * RHS prepares velocity traces before each Riemann evaluation, so callers
 * handing this object directly to that RHS need only this scalar pass. */
static inline int gfs_material_image_trace2d_prepare_scalar(const GfsMaterialExterior2D*g,
 GfsMaterialImageTrace2D*a,GfsMaterialProjectScalar2D project){
 if(!g||!a||!project)return 0;
 GfsMaterialExteriorTrace2D*t=&a->fluid;int n=g->n,ok=0,maxdegree=0;
 int*off=calloc((size_t)n+1,sizeof(int));double*grad=calloc((size_t)4*n,sizeof(double));
 if(!off||!grad)goto done;
 for(int k=0;k<g->directed_count;k++)off[g->directed[k].i+1]++;
 for(int i=0;i<n;i++){if(off[i+1]>maxdegree)maxdegree=off[i+1];off[i+1]+=off[i];}
 if(!gfs_material_project_scalar_batch2d(n,g->nt,maxdegree,off,
    a->physical,t->fit,t->state,a->wall,grad,project))goto done;
 for(int k=0;k<g->pair_count;k++){
  int z=g->pair_directed[k];const GfsTraceDirected2D*e=g->directed+z;GfsTraceFace2D*f=t->trace.face+k;
  for(int side=0;side<2;side++){
   int i=side&&!a->wall[z]?e->j:e->i;
   if(i>=n){f->thermo[2]=t->state[f->j].rho;f->thermo[3]=t->state[f->j].p;continue;}
   int other=side&&!a->wall[z];double x=t->fit[z].cx-(other?t->fit[z].dx:0),y=t->fit[z].cy-(other?t->fit[z].dy:0);
   for(int s=0;s<2;s++){
    double value=(s?t->state[i].p:t->state[i].rho)+grad[4*i+2*s]*x+grad[4*i+2*s+1]*y;
    if(!(value>0)||!isfinite(value))goto done;f->thermo[2*side+s]=value;
   }
  }
 }
 ok=1;
 done:free(off);free(grad);return ok;
}
/* Preserve the full preparation contract for callers that consume traces
 * directly. Scalar-only preparation leaves velocity_alpha unprepared. */
static inline int gfs_material_image_trace2d_prepare(const GfsMaterialExterior2D*g,
 GfsMaterialImageTrace2D*a,GfsMaterialProjectScalar2D project){
 if(!gfs_material_image_trace2d_prepare_scalar(g,a,project))return 0;
 return gfs_material_face_normal_prepare2d(&a->fluid.trace,a->fluid.state);
}
#endif
