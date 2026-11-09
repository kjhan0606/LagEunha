#ifndef GFS_MATERIAL_CYLINDER_RHS2D_H
#define GFS_MATERIAL_CYLINDER_RHS2D_H
#include <time.h>
#include "exam_gpu.h"
#include "../gfs_material_radial_image.h"
#include "../gfs_material_exterior_alias2d.h"
#include "../gfs_material_image_rhs2d.h"
#ifdef GFS_MATERIAL_SWEPT_CFL
#include "../gfs_material_power_face_motion.h"
#include "../gfs_material_swept_measure.h"
#endif

/* Native SoA adapter for the existing mixed periodic/image/open operator.
 * The caller owns padding, geometric precision sidecars and stage time.
 * alias identifies ONLY certified periodic copies; image identifies ONLY
 * radial wall copies. Other exterior slots carry prescribed m,U,V,v,W.
 * No particle motion, ownership changes, remapping or wall impulse here.
 * Native RunCylinder must connect these separately before enabling this RHS.
 */
typedef struct {
 GfsMaterialImageRhsStats2D exchange;
 int faces;
 double thermal_ledger[3];
 double laboratory_crossing_time;
 double geometry_timestep,geometry_sweep_coefficient;
 int geometry_limiting_cell;
} GfsMaterialCylinderNativeStats2D;

static inline double gfs_material_cylinder_profile_seconds(void){
 struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);
 return (double)t.tv_sec+1e-9*t.tv_nsec;
}

static inline int gfs_material_cylinder_native_rhs2d(GPUContext*ctx,int n,int nt,
 double gamma,const double*center,double radius,const int*alias,const int*image,
 const double*lowx,const double*lowy,GfsMaterialProjectScalar2D project,
 GfsMaterialImageRhsEval2D evaluate,double*volume_rate,
 GfsMaterialCylinderNativeStats2D*stats){
 if(!ctx||n<1||nt<=n||!(gamma>1)||!isfinite(gamma)||!center||!(radius>0)||
    !isfinite(radius)||!alias||!image||!lowx||!lowy||!project||!evaluate||
    !volume_rate||!stats)return 0;
 ParticleSoA*p=&ctx->h_parts;FaceCSR*f=&ctx->h_faces;
 GfsMaterialExterior2D g={0};GfsMaterialImageTrace2D trace={0};
 GfsMaterialCylinderNativeStats2D report={0};int ok=0;
 int profile=getenv("GFS_MATERIAL_CYLINDER_PROFILE")!=NULL;
 double times[7]={0};
 if(profile)times[0]=gfs_material_cylinder_profile_seconds();
 GfsFvPrimitive *state=NULL,*mesh=NULL,*face=NULL;
 double *offset=NULL,*known=NULL,*q=NULL,*realized=NULL,*speed=NULL;
 double *self_speed=NULL,*rhs=NULL,*area=NULL,*sweep=NULL;
 if(!gfs_material_exterior2d_build_coordinates(n,nt,f->face_offset[n],
    f->face_offset,f->neighbor_idx,p->x,p->y,p->weight_W,
    f->c1x,f->c1y,f->c2x,f->c2y,p->weight_volume64,lowx,lowy,&g)||
    !gfs_material_exterior_alias2d(&g,alias))goto done;
 if(profile)times[1]=gfs_material_cylinder_profile_seconds();
 state=calloc(nt,sizeof(*state));mesh=calloc(nt,sizeof(*mesh));
 face=calloc(g.pair_count,sizeof(*face));offset=calloc((size_t)2*nt,sizeof(double));
 known=calloc(nt,sizeof(double));q=calloc(n,sizeof(double));
 realized=calloc(n,sizeof(double));speed=calloc(g.pair_count,sizeof(double));
 self_speed=calloc(g.pair_count,sizeof(double));area=calloc(n,sizeof(double));
 if(!state||!mesh||!face||!offset||!known||!q||!realized||!speed||!self_speed||!area)goto done;
 for(int i=0;i<n;i++){
  double V=p->weight_volume64[i];
  if(!(V>0)||!isfinite(V)||!(p->mass[i]>0)||!isfinite(p->mass[i])||
     !(p->weight_U64[i]>0)||!isfinite(p->weight_U64[i]))goto done;
  state[i]=(GfsFvPrimitive){p->mass[i]/V,(gamma-1)*p->weight_U64[i]/V,
                          {p->vx[i],p->vy[i],0}};
  mesh[i]=state[i];
 }
 for(int i=n;i<nt;i++){
  if(image[i]<-1||image[i]>=n||(image[i]>=0&&alias[i]>=0))goto done;
  int j=image[i]>=0?image[i]:alias[i];
  if(j>=0){
   state[i]=state[j];mesh[i]=mesh[j];
   double off[2]={g.center[2*j],g.center[2*j+1]};
   if(image[i]>=0){
    GfsMaterialRadialImage a;double position[2]={p->x[j],p->y[j]};
    if(!gfs_material_radial_image(2,center,radius,position,&a)||
       !gfs_material_radial_vector(2,&a,state[j].u,0,state[i].u)||
       !gfs_material_radial_vector(2,&a,mesh[j].u,1,mesh[i].u)||
       !gfs_material_radial_vector(2,&a,off,0,offset+2*i))goto done;
   }else{offset[2*i]=off[0];offset[2*i+1]=off[1];}
  }else{
   double V=p->weight_volume64[i];
   if(!(V>0)||!isfinite(V)||!(p->mass[i]>0)||!isfinite(p->mass[i])||
      !(p->weight_U64[i]>0)||!isfinite(p->weight_U64[i]))goto done;
   state[i]=(GfsFvPrimitive){p->mass[i]/V,(gamma-1)*p->weight_U64[i]/V,
                           {p->vx[i],p->vy[i],0}};
   mesh[i]=state[i]; /* ballistic/prescribed samples: known Wdot = 0 */
  }
 }
 for(int k=0;k<g.pair_count;k++){
  const GfsTraceDirected2D*e=g.directed+g.pair_directed[k];
  if(e->j>=n&&image[e->j]<0)face[k]=state[e->j];
  if(e->j>=n&&image[e->j]==e->i){
   double position[2]={p->x[e->i],p->y[e->i]},c[2]={e->cx,e->cy};
   if(!gfs_material_radial_self_face_speed(2,center,radius,position,c,
                                         state[e->i].u,self_speed+k))goto done;
  }
 }
 if(profile)times[2]=gfs_material_cylinder_profile_seconds();
 if(!gfs_material_image_trace2d_build(&g,image,state,offset,face,&trace))goto done;
 if(profile)times[3]=gfs_material_cylinder_profile_seconds();
 /* image_rhs prepares normal velocities at zero and then recovered wall
  * speed. No consumer needs a third, immediately overwritten normal pass. */
 if(!gfs_material_image_trace2d_prepare_scalar(&g,&trace,project))goto done;
 if(profile)times[4]=gfs_material_cylinder_profile_seconds();
 rhs=calloc((size_t)4*trace.fluid.trace.cells,sizeof(double));if(!rhs)goto done;
 if(!gfs_material_image_rhs2d_with_self_speeds(&g,&trace,mesh,image,known,gamma,
    evaluate,q,realized,speed,rhs,&report.exchange,self_speed))goto done;
 if(profile)times[5]=gfs_material_cylinder_profile_seconds();
#ifdef GFS_MATERIAL_SWEPT_CFL
 /* Isolated candidate: reuse the existing Kepler coefficient, but integrate
  * |relative normal speed| over the whole straight face. A rotating face
  * can have zero midpoint speed and still sweep a finite area. */
 const char*sweep_env=getenv("GFS_VARPRESS_GEOM_CFL");char*parsed=NULL;
 double sweep_coefficient=sweep_env?strtod(sweep_env,&parsed):0;
 if(!isfinite(sweep_coefficient)||sweep_coefficient<0||sweep_coefficient>1||
    (sweep_env&&(!*sweep_env||*parsed)))goto done;
 report.geometry_limiting_cell=-1;
 if(sweep_coefficient>0){
  sweep=calloc(n,sizeof(double));if(!sweep)goto done;
  for(int k=0;k<g.directed_count;k++){
   const GfsTraceDirected2D*e=g.directed+k;int i=e->i,j=e->j;
   double A=e->area;
   double absolute;
#ifdef GFS_MATERIAL_FIXED_SOLID
   if(i==j){double d=hypot(e->dx,e->dy);absolute=A*fabs((e->dx*mesh[i].u[0]+e->dy*mesh[i].u[1])/d);}
   else
#endif
   if(j>=n&&image[j]==i){
    double d=hypot(e->dx,e->dy);
    double own=A*(e->dx*mesh[i].u[0]+e->dy*mesh[i].u[1])/d;
    double position[2]={p->x[i],p->y[i]},c[2]={e->cx,e->cy},velocity;
    if(!gfs_material_radial_self_face_speed(2,center,radius,position,c,mesh[i].u,&velocity))goto done;
    double mean=A*velocity-own;
    double x=position[0]-center[0],y=position[1]-center[1],r2=x*x+y*y;
    double span=.5*A*A*(x*mesh[i].u[1]-y*mesh[i].u[0])/r2;
    double endpoints[2]={mean-span,mean+span};
    if(!gfs_material_affine_abs_mean(2,endpoints,&absolute))goto done;
   }else{
    double delta[2]={e->dx,e->dy},c[2]={e->cx,e->cy};
    double qj=j<n?q[j]:(image[j]>=0?q[image[j]]:known[j]);
    if(!gfs_material_power_segment_sweep(delta,c,A,mesh[i].u,mesh[j].u,q[i],qj,&absolute))goto done;
   }
   sweep[i]+=absolute;if(!isfinite(sweep[i]))goto done;
  }
  double bound=HUGE_VAL;
  for(int i=0;i<n;i++)if(sweep[i]>0){
   double dt=sweep_coefficient*p->weight_volume64[i]/sweep[i];
   if(!(dt>0)||!isfinite(dt))goto done;
   if(dt<bound){bound=dt;report.geometry_limiting_cell=i;}
  }
  report.geometry_timestep=bound;report.geometry_sweep_coefficient=sweep_coefficient;
 }
#endif
 for(int k=0;k<g.directed_count;k++)area[g.directed[k].i]+=g.directed[k].area;
 report.laboratory_crossing_time=HUGE_VAL;report.faces=g.pair_count;
 for(int i=0;i<n;i++){
  double udot=rhs[4*i+3]-state[i].p*rhs[4*i+2];
  double bound=p->weight_volume64[i]/(area[i]*(sqrt(gamma*state[i].p/state[i].rho)+
                                          hypot(state[i].u[0],state[i].u[1])));
  if(!isfinite(udot)||!isfinite(q[i])||!(bound>0)||!isfinite(bound)||
     !isfinite(rhs[4*i]/p->mass[i])||!isfinite(rhs[4*i+1]/p->mass[i]))goto done;
  report.laboratory_crossing_time=fmin(report.laboratory_crossing_time,bound);
  report.thermal_ledger[0]+=rhs[4*i+3];
  report.thermal_ledger[1]-=state[i].p*rhs[4*i+2];
  report.thermal_ledger[2]+=state[i].u[0]*rhs[4*i]+state[i].u[1]*rhs[4*i+1];
 }
 /* Publish rates only after all cells succeed. The caller applies its same
  * volume/energy tolerances and records wall power separately from physical
  * stationary-wall conservation. This bound is the existing comparison's
  * laboratory-frame bound; it is not a new native face-CFL qualification. */
 for(int i=0;i<n;i++){
  p->ax_out[i]=rhs[4*i]/p->mass[i];p->ay_out[i]=rhs[4*i+1]/p->mass[i];
  p->weight_Udot64[i]=rhs[4*i+3]-state[i].p*rhs[4*i+2];
  p->die_out[i]=p->weight_Udot64[i];p->weight_rate[i]=q[i];p->weight_force[i]=0;
  volume_rate[i]=rhs[4*i+2];
 }
 *stats=report;ok=1;
 if(profile)times[6]=gfs_material_cylinder_profile_seconds();
 done:
 free(state);free(mesh);free(face);free(offset);free(known);free(q);free(realized);
 free(speed);free(self_speed);free(rhs);free(area);free(sweep);
 gfs_material_image_trace2d_free(&trace);gfs_material_exterior2d_free(&g);
 if(ok&&profile)fprintf(stderr,"[CYLINDER_RHS_PROFILE] n=%d geometry_alias_ms=%.9g state_ms=%.9g trace_build_ms=%.9g scalar_prepare_ms=%.9g exchange_ms=%.9g ledger_publish_ms=%.9g cleanup_ms=%.9g total_ms=%.9g\n",n,
  1e3*(times[1]-times[0]),1e3*(times[2]-times[1]),1e3*(times[3]-times[2]),
  1e3*(times[4]-times[3]),1e3*(times[5]-times[4]),1e3*(times[6]-times[5]),
  1e3*(gfs_material_cylinder_profile_seconds()-times[6]),1e3*(gfs_material_cylinder_profile_seconds()-times[0]));
 return ok;
}
#endif
