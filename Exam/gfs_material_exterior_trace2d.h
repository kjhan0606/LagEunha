#ifndef GFS_MATERIAL_EXTERIOR_TRACE2D_H
#define GFS_MATERIAL_EXTERIOR_TRACE2D_H
#include "gfs_material_exterior2d.h"
#include "gfs_material_face_normal2d.h"
/* Prescribed exterior data live at two explicitly distinguished locations:
 * physical samples used by the owned-cell gradient fit, and boundary face
 * states used by the Riemann problem. A boundary datum is not a new fluid
 * degree of freedom. The same existing C/CUDA face kernel consumes this CSR.
 */
typedef struct {
 GfsTraceCsr2D trace;
 GfsTraceDirected2D *fit;
 GfsFvPrimitive *state;
 int *boundary_slot;
} GfsMaterialExteriorTrace2D;
static inline void gfs_material_exterior_trace2d_free(GfsMaterialExteriorTrace2D*t){
 if(!t)return;gfs_trace_csr2d_free(&t->trace);free(t->fit);free(t->state);free(t->boundary_slot);memset(t,0,sizeof(*t));
}
/* sample_offset[2*j:2*j+2] is the prescribed physical sample location
 * relative to generator j for j>=n. Owned sample locations use the verified
 * volume centroid. exterior_face has pair_count entries; only external ones
 * are read. State and geometry must belong to this same immutable RHS. */
static inline int gfs_material_exterior_trace2d_build(const GfsMaterialExterior2D*g,
 const GfsFvPrimitive*sample,const double*sample_offset,
 const GfsFvPrimitive*exterior_face,GfsMaterialExteriorTrace2D*out){
 if(!out)return 0;memset(out,0,sizeof(*out));
 if(!g||g->n<1||!sample||!sample_offset||!exterior_face)return 0;
 GfsMaterialExteriorTrace2D t={0};int ok=0,nb=0;
 for(int k=0;k<g->pair_count;k++)if(g->pair[k].j>=g->n)nb++;
 if(nb<1||g->nt>INT_MAX-nb)return 0;
 t.fit=(GfsTraceDirected2D*)malloc((size_t)g->directed_count*sizeof(*t.fit));
 t.boundary_slot=(int*)malloc((size_t)g->directed_count*sizeof(int));
 t.state=(GfsFvPrimitive*)calloc((size_t)g->nt+nb,sizeof(*t.state));
 if(!t.fit||!t.boundary_slot||!t.state)goto done;
 memcpy(t.state,sample,(size_t)g->nt*sizeof(*sample));
 for(int k=0;k<g->directed_count;k++){
  const GfsTraceDirected2D*e=g->directed+k;int i=e->i,j=e->j;
  double cx=j<g->n?g->center[2*j]:sample_offset[2*j],cy=j<g->n?g->center[2*j+1]:sample_offset[2*j+1];
  if(!isfinite(cx)||!isfinite(cy))goto done;
  t.fit[k]=*e;t.fit[k].dx+=cx-g->center[2*i];t.fit[k].dy+=cy-g->center[2*i+1];
  t.fit[k].cx-=g->center[2*i];t.fit[k].cy-=g->center[2*i+1];t.boundary_slot[k]=-1;
 }
 int slot=g->nt;
 for(int k=0;k<g->pair_count;k++)if(g->pair[k].j>=g->n){
  GfsFvPrimitive s=exterior_face[k];
  if(!(s.rho>0)||!(s.p>0)||!isfinite(s.rho)||!isfinite(s.p)||!isfinite(s.u[0])||!isfinite(s.u[1]))goto done;
  t.boundary_slot[g->pair_directed[k]]=slot;t.state[slot++]=s;
 }
 if(slot!=g->nt+nb||!gfs_trace_csr2d_build_coordinates_domain(g->n,g->nt,slot,g->directed_count,g->directed,&t.trace,NULL,t.fit,t.boundary_slot))goto done;
 *out=t;memset(&t,0,sizeof(t));ok=1;
 done:gfs_material_exterior_trace2d_free(&t);return ok;
}
#ifdef GFS_MATERIAL_GRADIENT_PROJECT
#include "gfs_material_scalar_batch2d.h"
/* Reuse the selected native one-cell gradient implementation by passing its
 * function. This adapter introduces no alternative slope limiter or physics.
 */
static inline int gfs_material_exterior_trace2d_prepare(const GfsMaterialExterior2D*g,
 GfsMaterialExteriorTrace2D*t,GfsMaterialProjectScalar2D project){
 if(!g||!t||!project)return 0;
 int n=g->n,ok=0,maxdegree=0,*off=(int*)calloc((size_t)n+1,sizeof(int));
 double*grad=(double*)calloc((size_t)4*n,sizeof(double));
 if(!off||!grad)goto done;
 for(int k=0;k<g->directed_count;k++)off[g->directed[k].i+1]++;
 for(int i=0;i<n;i++){if(off[i+1]>maxdegree)maxdegree=off[i+1];off[i+1]+=off[i];}
 if(!gfs_material_project_scalar_batch2d(n,g->nt,maxdegree,off,
    g->directed,t->fit,t->state,NULL,grad,project))goto done;
 if(t->trace.faces!=g->pair_count)goto done;
 for(int k=0;k<g->pair_count;k++){
  int z=g->pair_directed[k];const GfsTraceDirected2D*e=g->directed+z;GfsTraceFace2D*f=t->trace.face+k;
  if(f->i!=e->i||f->j!=(e->j<n?e->j:t->boundary_slot[z]))goto done;
  for(int side=0;side<2;side++){
   int i=side?e->j:e->i;
   if(i>=n){f->thermo[2]=t->state[f->j].rho;f->thermo[3]=t->state[f->j].p;continue;}
   double x=t->fit[z].cx-(side?t->fit[z].dx:0),y=t->fit[z].cy-(side?t->fit[z].dy:0);
   for(int scalar=0;scalar<2;scalar++){
    double v=(scalar?t->state[i].p:t->state[i].rho)+grad[4*i+2*scalar]*x+grad[4*i+2*scalar+1]*y;
    if(!(v>0)||!isfinite(v))goto done;f->thermo[2*side+scalar]=v;
   }
  }
 }
 if(!gfs_material_face_normal_prepare2d(&t->trace,t->state))goto done;
 ok=1;
 done:free(off);free(grad);return ok;
}
#endif
/* RHS convention H=Edot-v.F+P.Vdot. Return the exterior-to-owned ledger;
 * callers evolve only [0,n), never the prescribed-data slots. */
static inline int gfs_material_exterior_trace2d_ledger(int n,const GfsMaterialExteriorTrace2D*t,
 const double*rhs,double ledger[4]){
 if(!t||!rhs||!ledger||n<1||n>=t->trace.cells)return 0;
 for(int a=0;a<4;a++)ledger[a]=0;
 for(int i=n;i<t->trace.cells;i++){
  double energy=rhs[4*i+3]-t->state[i].p*rhs[4*i+2]+t->state[i].u[0]*rhs[4*i]+t->state[i].u[1]*rhs[4*i+1];
  for(int a=0;a<3;a++){if(!isfinite(rhs[4*i+a]))return 0;ledger[a]-=rhs[4*i+a];}
  if(!isfinite(energy))return 0;ledger[3]-=energy;
 }return 1;
}
#endif
