#ifndef GFS_MATERIAL_SSPRK2_H
#define GFS_MATERIAL_SSPRK2_H
#include <math.h>
#ifdef __CUDACC__
#define GFS_SSP_HD static __host__ __device__ inline
#else
#define GFS_SSP_HD static inline
#endif
/* Fixed-mass SSPRK2 state. W is constrained geometry, not a canonical
 * coordinate with a pi/inertia kinetic energy. No clipping or repair.
 * Periodic length 0 means a nonperiodic coordinate; otherwise displacements
 * must stay below half a period so the material trajectory is unambiguous.
 */
typedef struct {double mass,x[3],v[3],U,W;} GfsMaterialStepState;
typedef struct {double force[3],Udot,Wdot;} GfsMaterialStepRate;
/* Non-autonomous boundary/source data must follow the SSPRK2 abscissae:
 * RHS 1 at t, RHS 2 at t+h. The final geometry is also evaluated at t+h.
 * h is the current trial duration, not a preceding accepted/rejected step.
 * Reject an unrepresentable endpoint instead of advancing state at fixed t.
 */
GFS_SSP_HD int gfs_material_ssprk2_rhs_time(double t,double h,int stage,double *out){
 if(!out||!isfinite(t)||!(h>0)||!isfinite(h)||(stage!=1&&stage!=2))return 0;
 double end=t+h;
 if(!isfinite(end)||!(end>t))return 0;
 *out=stage==1?t:end;return 1;
}
/* External-force substep at fixed position: update momentum and its exact
 * kinetic work without subtracting bulk kinetic energy from internal energy.
 * For a second-order split use half kicks around the hydro step, evaluating
 * acceleration at each endpoint. This helper does not supply wall treatment,
 * hydrostatic balance, a timestep criterion, or native-driver integration.
 */
GFS_SSP_HD int gfs_material_external_kick(int dim,double dt,
 const double acceleration[3],const GfsMaterialStepState *in,
 GfsMaterialStepState *out,double *work){
 if(!in||!out||!acceleration||!work||dim<1||dim>3||
    !isfinite(dt)||dt<0||!(in->mass>0)||!isfinite(in->mass)||
    !(in->U>0)||!isfinite(in->U))return 0;
 GfsMaterialStepState q=*in;double w=0;
 for(int a=0;a<dim;a++){
  if(!isfinite(in->v[a])||!isfinite(acceleration[a]))return 0;
  double dv=dt*acceleration[a];q.v[a]=in->v[a]+dv;
  w+=in->mass*dv*(in->v[a]+.5*dv);
  if(!isfinite(q.v[a])||!isfinite(w))return 0;
 }
 *out=q;*work=w;return 1;
}
GFS_SSP_HD int gfs_material_ssprk2_stage(int dim,int stage,double dt,
 const double period[3],const GfsMaterialStepState *base,
 const GfsMaterialStepState *current,const GfsMaterialStepRate *rate,
 GfsMaterialStepState *out){
 if(!base||!current||!rate||!out||!period||dim<1||dim>3||
    (stage!=1&&stage!=2)||!(dt>0)||!isfinite(dt)||!(current->mass>0)||
    !isfinite(current->mass)||base->mass!=current->mass||
    !(base->U>0)||!(current->U>0)||!isfinite(rate->Udot)||!isfinite(rate->Wdot))return 0;
 GfsMaterialStepState q=*current;double E=current->U,Ebase=base->U,Edot=rate->Udot;
 for(int a=0;a<dim;a++){
  if(!isfinite(current->x[a])||!isfinite(current->v[a])||!isfinite(base->x[a])||
     !isfinite(base->v[a])||!isfinite(rate->force[a])||period[a]<0||!isfinite(period[a]))return 0;
  E+=.5*current->mass*current->v[a]*current->v[a];
  Ebase+=.5*base->mass*base->v[a]*base->v[a];Edot+=current->v[a]*rate->force[a];
  double displacement=dt*current->v[a];
  if(period[a]>0&&fabs(displacement)>=.5*period[a])return 0;
  q.x[a]=current->x[a]+displacement;q.v[a]=current->v[a]+dt*rate->force[a]/current->mass;
  if(stage==2){
   double delta=current->x[a]-base->x[a];
   if(period[a]>0)delta-=period[a]*nearbyint(delta/period[a]);
   q.x[a]=base->x[a]+.5*(delta+displacement);q.v[a]=.5*(base->v[a]+q.v[a]);
  }
  if(period[a]>0)q.x[a]-=period[a]*floor(q.x[a]/period[a]);
  if(!isfinite(q.x[a])||!isfinite(q.v[a]))return 0;
 }
 E+=dt*Edot;q.W=current->W+dt*rate->Wdot;
 if(stage==2){E=.5*(E+Ebase);q.W=.5*(base->W+q.W);}
 q.U=E;for(int a=0;a<dim;a++)q.U-=.5*q.mass*q.v[a]*q.v[a];
 if(!(q.U>0)||!isfinite(q.U)||!isfinite(q.W))return 0;
 *out=q;return 1;
}
#undef GFS_SSP_HD
#endif
