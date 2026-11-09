#ifndef GFS_MATERIAL_EXTERIOR_NATIVE2D_H
#define GFS_MATERIAL_EXTERIOR_NATIVE2D_H
#include "exam_gpu.h"
#include "../gfs_material_rhs2d.h"
#include "../gfs_material_cuda_rhs.h"
#include "../gfs_material_exterior_trace2d.h"
#if !defined(GFS_MATERIAL_DIRECT_FV) || !defined(GFS_MATERIAL_GRADIENT_PROJECT)
#error Exterior native bridge requires the direct pressure-work and projected-gradient method
#endif
/* The driver supplies exterior fields at the current RK-stage time and
 * positions. This bridge never selects a boundary physics model, moves a site,
 * changes mass, or advances prescribed samples as fluid unknowns. */
typedef int (*GfsMaterialExteriorState2D)(double,double,void*,GfsFvPrimitive*);
typedef struct {
 int faces,exterior_faces,weight_iterations;
 double weight_residual,volume_rate_error,boundary_ledger[4],energy_rate_defect,face_crossing_time;
 double thermal_ledger[3]; /* heat, reversible pressure work, kinetic work */
 double geometry_timestep,geometry_sweep_coefficient;
 int geometry_limiting_cell;
} GfsMaterialExteriorNativeStats2D;
static inline int gfs_material_exterior_native_rhs2d(GPUContext*ctx,int n,int nt,
 double gamma,GfsMaterialExteriorState2D prescribe,void*data,
 const double*sample_offset,const double*known_weight_rate,
 GfsMaterialProjectScalar2D project,int use_cuda,double*realized_volume_rate,GfsMaterialExteriorNativeStats2D*stats){
 if(!ctx||n<1||nt<=n||!(gamma>1)||!isfinite(gamma)||!prescribe||!sample_offset||!known_weight_rate||!project||!stats)return 0;
 GfsMaterialExterior2D g={0};GfsMaterialExteriorTrace2D t={0};GfsMaterialExteriorNativeStats2D report={0};
 ParticleSoA*p=&ctx->h_parts;FaceCSR*f=&ctx->h_faces;
 GfsFvPrimitive*sample=NULL,*face=NULL;double*rhs=NULL,*rates=NULL,*sweep=NULL;int ok=0;const char*phase="geometry";
 if(!gfs_material_exterior2d_build(n,nt,f->face_offset[n],f->face_offset,f->neighbor_idx,p->x,p->y,p->weight_W,f->c1x,f->c1y,f->c2x,f->c2y,p->weight_volume64,&g))goto done;
 phase="samples";
 sample=(GfsFvPrimitive*)calloc(nt,sizeof(*sample));face=(GfsFvPrimitive*)calloc(g.pair_count,sizeof(*face));
 rates=(double*)calloc((size_t)3*n,sizeof(double));if(!sample||!face||!rates)goto done;
 for(int i=0;i<n;i++){
  if(!(p->mass[i]>0)||!isfinite(p->mass[i])||!(p->weight_U64[i]>0)||!isfinite(p->weight_U64[i]))goto done;
  sample[i]=(GfsFvPrimitive){p->mass[i]/p->weight_volume64[i],(gamma-1)*p->weight_U64[i]/p->weight_volume64[i],{p->vx[i],p->vy[i],0}};
 }
 for(int i=n;i<nt;i++)if(!prescribe(p->x[i]+sample_offset[2*i],p->y[i]+sample_offset[2*i+1],data,sample+i))goto done;
 for(int k=0;k<g.pair_count;k++)if(g.pair[k].j>=n){
  const GfsTraceDirected2D*e=g.directed+g.pair_directed[k];
  if(!prescribe(p->x[e->i]+e->cx,p->y[e->i]+e->cy,data,face+k))goto done;report.exterior_faces++;
 }
 phase="trace_build";if(!gfs_material_exterior_trace2d_build(&g,sample,sample_offset,face,&t))goto done;
 phase="trace_prepare";if(!gfs_material_exterior_trace2d_prepare(&g,&t,project))goto done;
 rhs=(double*)calloc((size_t)4*t.trace.cells,sizeof(double));if(!rhs)goto done;
 phase="rhs";if((use_cuda?gfs_material_cuda_rhs2d(&t.trace,t.state,gamma,NULL,rhs):gfs_material_rhs_maps2d(&t.trace,t.state,gamma,NULL,rhs))!=1)goto done;
 phase="weight_solve";
 if(!gfs_material_exterior2d_weight_rate(&g,sample,rhs,known_weight_rate,rates,rates+n,&report.weight_iterations,&report.weight_residual))goto done;
 phase="ledger";if(!gfs_material_exterior_trace2d_ledger(n,&t,rhs,report.boundary_ledger))goto done;
 report.faces=g.pair_count;report.energy_rate_defect=-report.boundary_ledger[3];
 for(int i=0;i<n;i++){
  double udot=-sample[i].p*rates[n+i]+rhs[4*i+3];rates[2*n+i]=udot;
  if(!isfinite(udot)||!isfinite(rates[i])||!isfinite(rhs[4*i]/p->mass[i])||!isfinite(rhs[4*i+1]/p->mass[i]))goto done;
  report.energy_rate_defect+=udot+sample[i].u[0]*rhs[4*i]+sample[i].u[1]*rhs[4*i+1];
  report.thermal_ledger[0]+=rhs[4*i+3];
  report.thermal_ledger[1]+=-sample[i].p*rates[n+i];
  report.thermal_ledger[2]+=sample[i].u[0]*rhs[4*i]+sample[i].u[1]*rhs[4*i+1];
  report.volume_rate_error=fmax(report.volume_rate_error,fabs(rates[n+i]-rhs[4*i+2]));
 }
 if(!isfinite(report.energy_rate_defect))goto done;
 /* Same material face-crossing CFL as the periodic native adapter, with
  * known exterior Wdot substituted for the non-owned endpoint. */
 phase="crossing";memset(rates+2*n,0,(size_t)n*sizeof(double));
#ifdef GFS_MATERIAL_SWEPT_CFL
 const char*ge=getenv("GFS_VARPRESS_GEOM_CFL");char*end=NULL;
 double cg=ge?strtod(ge,&end):0;
 if(!isfinite(cg)||cg<0||cg>1||(ge&&(!*ge||*end)))goto done;
 report.geometry_limiting_cell=-1;
 if(cg>0){sweep=calloc((size_t)n,sizeof(double));if(!sweep)goto done;}
#endif
 for(int k=0;k<g.directed_count;k++){
  const GfsTraceDirected2D*e=g.directed+k;int i=e->i,j=e->j;
  double d[3]={e->dx,e->dy,0},c[3]={e->cx,e->cy,0},base,coef,dist=hypot(e->dx,e->dy);
  if(!gfs_material_power_face_rate(2,d,c,sample[i].u,sample[j].u,e->area,&base,&coef))goto done;
  double wn=(base+coef*(rates[i]-(j<n?rates[j]:known_weight_rate[j])))/e->area;
  double vi=(d[0]*sample[i].u[0]+d[1]*sample[i].u[1])/dist,vj=(d[0]*sample[j].u[0]+d[1]*sample[j].u[1])/dist;
  double ci=sqrt(gamma*sample[i].p/sample[i].rho),cj=sqrt(gamma*sample[j].p/sample[j].rho);
  double signal=fmax(ci+fabs(vi-wn),cj+fabs(vj-wn));if(!(signal>0)||!isfinite(signal))goto done;
  rates[2*n+i]+=e->area*signal;
#ifdef GFS_MATERIAL_SWEPT_CFL
  if(sweep){
   double absolute;
   if(!gfs_material_power_segment_sweep(d,c,e->area,sample[i].u,sample[j].u,
      rates[i],j<n?rates[j]:known_weight_rate[j],&absolute))goto done;
   sweep[i]+=absolute;if(!isfinite(sweep[i]))goto done;
  }
#endif
 }
 report.face_crossing_time=HUGE_VAL;
 for(int i=0;i<n;i++){
  if(!(rates[2*n+i]>0)||!isfinite(rates[2*n+i]))goto done;
  report.face_crossing_time=fmin(report.face_crossing_time,p->weight_volume64[i]/rates[2*n+i]);
  rates[2*n+i]=-sample[i].p*rates[n+i]+rhs[4*i+3];
 }
 if(!(report.face_crossing_time>0)||!isfinite(report.face_crossing_time))goto done;
#ifdef GFS_MATERIAL_SWEPT_CFL
 if(sweep){
  report.geometry_timestep=HUGE_VAL;report.geometry_sweep_coefficient=cg;
  for(int i=0;i<n;i++)if(sweep[i]>0){
   double h=cg*p->weight_volume64[i]/sweep[i];if(!(h>0)||!isfinite(h))goto done;
   if(h<report.geometry_timestep){report.geometry_timestep=h;report.geometry_limiting_cell=i;}
  }
 }
#endif
 /* Commit native rates only after every cell/face/solve has succeeded. Same
  * field conventions as the existing periodic material adapter. */
 for(int i=0;i<n;i++){
  p->ax_out[i]=rhs[4*i]/p->mass[i];p->ay_out[i]=rhs[4*i+1]/p->mass[i];
  p->weight_Udot64[i]=rates[2*n+i];p->die_out[i]=rates[2*n+i];
  p->weight_rate[i]=rates[i];p->weight_force[i]=0;
 }
 /* Optional observation of the SAME realized rate used in pressure work.
  * The caller owns this output; exterior samples are not material unknowns. */
 if(realized_volume_rate)memcpy(realized_volume_rate,rates+n,(size_t)n*sizeof(double));
 *stats=report;ok=1;
 done:
 if(!ok&&getenv("GFS_MATERIAL_EXTERIOR_DIAGNOSTICS")){
  fprintf(stderr,"[EXTERIOR_RHS_FAIL] phase=%s\n",phase);
  if(!strcmp(phase,"geometry")){
   FILE*dump=fopen("exterior_failed_geometry.txt","wx");
   if(dump){fprintf(dump,"%d %d %d\n",n,nt,f->face_offset[n]);
    for(int i=0;i<nt;i++)fprintf(dump,"%d %.17g %.17g %.17g %.17g\n",i,p->x[i],p->y[i],p->weight_W[i],p->weight_volume64[i]);
    for(int i=0;i<n;i++)for(int k=f->face_offset[i];k<f->face_offset[i+1];k++)fprintf(dump,"%d %d %.17g %.17g %.17g %.17g\n",i,f->neighbor_idx[k],f->c1x[k],f->c1y[k],f->c2x[k],f->c2y[k]);
    fclose(dump);
   }
  }
 }
 free(sweep);free(sample);free(face);free(rhs);free(rates);gfs_material_exterior_trace2d_free(&t);gfs_material_exterior2d_free(&g);return ok;
}
#endif
