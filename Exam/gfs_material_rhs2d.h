#include <stdio.h>
#ifndef GFS_MATERIAL_RHS2D_H
#define GFS_MATERIAL_RHS2D_H
#include "gfs_material_trace_csr2d.h"
#include "gfs_material_face_response.h"
#include "gfs_material_face_normal2d.h"
#include "gfs_material_power_face_motion.h"
#ifdef GFS_MATERIAL_DIRECT_FV
#include "gfs_material_direct_flux.h"
#endif
#ifdef GFS_MATERIAL_STAR_ENTROPY
#ifndef GFS_MATERIAL_DIRECT_FV
#error Material star entropy candidate requires direct force/volume/work
#endif
#include "gfs_material_star_entropy.h"
#endif
#ifdef GFS_MATERIAL_CUDA_WEIGHT
#include "gfs_material_cuda_weight.h"
#endif
/* In-memory CPU reference for the same transpose assembly used on CUDA.
 * Output layout per cell: Fx,Fy,dV/dt,H. Caller owns zeroing/output transaction.
 * Failure invalidates the entire output, never a partial force update.
 */
/* Optional per-cell alpha[rho,P,velocity] scales the reconstructed increment.
 * Velocity uses one scalar for both components. Each field's matching adjoint
 * must use its own coefficients; density limiting must not limit velocity.
 * NULL recovers the existing unrestricted operator without extra arithmetic.
 * This kernel alone does not select/qualify a limiter or constrain weights. */
static inline int gfs_material_rhs_shifted_owner_boundary2d(const GfsTraceCsr2D *t,
 const GfsFvPrimitive *cell,double gamma,const double *alpha3,double *out,
 const double *pressure_shift,double *source_power,const double *owner_shift,
 const double *wall_normal,double *boundary_power){
 if(!t||!cell||!out||(pressure_shift&&!source_power))return 0;
 if(wall_normal&&!boundary_power)return 0;
#if !defined(GFS_MATERIAL_DIRECT_FV) || !defined(GFS_MATERIAL_DIRECT_WALLS)
 if(wall_normal)return 0;
#endif
 if(boundary_power)*boundary_power=0;
 #ifdef GFS_MATERIAL_DIRECT_FV
 if(pressure_shift)return 0; /* Affine wall/gravity work not qualified here. */
#ifndef GFS_MATERIAL_DIRECT_WALLS
 for(int j=0;j<t->entries;j++)if(t->entry[j].pad)return 0;
#endif
#else
 if(owner_shift)return 0;
#endif
 if(source_power)*source_power=0;
 if(alpha3)for(int i=0;i<3*t->cells;i++)if(!isfinite(alpha3[i])||alpha3[i]<0||alpha3[i]>1)return 0;
 memset(out,0,(size_t)t->cells*4*sizeof(double));
 for(int k=0;k<t->faces;k++){
  const GfsTraceFace2D *f=t->face+k;GfsFvPrimitive state[2]={{0}};
  int start[2]={f->sl,f->sr},len[2]={f->nl,f->nr};double reversible_speed[2]={0},offset[2]={0};
  for(int side=0;side<2;side++)for(int j=start[side];j<start[side]+len[side];j++){
   const GfsTraceEntry *e=t->entry+j;GfsFvPrimitive c=cell[e->cell];
   int owner=side?f->j:f->i;
   GfsMaterialTraceWeights w=gfs_material_entry_weights(e,owner,alpha3?alpha3+3*owner:NULL);
   state[side].rho+=w.rho*c.rho;state[side].p+=w.pressure*c.p;
   if(pressure_shift){
    if(!isfinite(pressure_shift[j]))return 0;
    double shift=w.pressure*pressure_shift[j];state[side].p+=shift;offset[side]+=shift;
   }
   for(int a=0;a<2;a++)state[side].u[a]+=w.velocity*gfs_trace_velocity_sign(e,a)*c.u[a];
   reversible_speed[side]+=w.reversible*(f->nx*c.u[0]+f->ny*gfs_trace_velocity_sign(e,1)*c.u[1]);
  }
#ifdef GFS_MATERIAL_GRADIENT_PROJECT
  for(int side=0;side<2;side++){state[side].rho=f->thermo[2*side];state[side].p=f->thermo[2*side+1];}
#endif
  if(owner_shift){if(!isfinite(owner_shift[f->i])||!isfinite(owner_shift[f->j]))return 0;state[0].p-=owner_shift[f->i];state[1].p-=owner_shift[f->j];}
  double normal[3]={f->nx,f->ny,0};GfsMaterialFaceResponse response;GfsMaterialFaceExchange x;
#ifdef GFS_MATERIAL_DIRECT_WALLS
  if(f->i==f->j){
   if(owner_shift)return 0;
   GfsMaterialDirectFlux wall;
   if(!gfs_material_direct_moving_wall_flux(2,state[0],normal,f->area,gamma,
       wall_normal?wall_normal[k]:0,&wall))return 0;
   out[4*f->i]+=wall.force[0];out[4*f->i+1]+=wall.force[1];
   if(wall_normal){
    out[4*f->i+2]+=wall.volume;out[4*f->i+3]+=wall.energy;
    *boundary_power+=wall.energy;if(!isfinite(*boundary_power))return 0;
   }
   continue;
  }
#endif
  int status=gfs_material_face_response(2,state[0],state[1],normal,gamma,&response);
  if(status!=1)return status==2?2:0;
#ifdef GFS_MATERIAL_DIRECT_FV
  #ifdef GFS_MATERIAL_STAR_ENTROPY
  double entropy_theta;
  if(!gfs_material_star_entropy_limit(2,cell[f->i],cell[f->j],normal,gamma,
       &response.pstar,&response.ustar,&entropy_theta))return 0;
  #endif
  GfsMaterialDirectFlux flux;if(!gfs_material_direct_flux(2,normal,f->area,&response,&flux))return 0;
  for(int side=0;side<2;side++){int i=side?f->j:f->i;double sign=side?-1:1;
   out[4*i]+=sign*flux.force[0];out[4*i+1]+=sign*flux.force[1];out[4*i+2]+=sign*flux.volume;out[4*i+3]+=sign*flux.energy;}
  continue;
#endif
  #ifdef GFS_MATERIAL_SPLIT_REVERSIBLE
  response.mean_speed=response.alpha_left*reversible_speed[0]+response.alpha_right*reversible_speed[1];
#endif
  if(!gfs_material_face_exchange(2,normal,f->area,cell[f->i].p,cell[f->j].p,&response,&x))return 0;
  /* Affine pressure traces J P+b: the quadratic heat minus P-weighted
   * volume work leaves b^T A J_full P/Z. Keep this explicit; it is not
   * automatically physical heating or moving-wall work. */
  if(source_power)*source_power+=(offset[0]-offset[1])*f->area*response.pressure_speed;
  for(int side=0;side<2;side++)for(int j=start[side];j<start[side]+len[side];j++){
   const GfsTraceEntry *e=t->entry+j;double *r=out+4*e->cell;
   int owner=side?f->j:f->i;
   GfsMaterialTraceWeights w=gfs_material_entry_weights(e,owner,alpha3?alpha3+3*owner:NULL);
   for(int a=0;a<2;a++)r[a]+=gfs_trace_velocity_sign(e,a)*gfs_material_scatter_force(w,&x,side,a);
   r[2]+=w.pressure*x.side_volume[side];
  }
  out[4*f->i+2]+=x.owner_volume[0];out[4*f->j+2]+=x.owner_volume[1];
  out[4*f->i+3]+=x.owner_heat[0];out[4*f->j+3]+=x.owner_heat[1];
 }
#ifdef GFS_MATERIAL_DIRECT_FV
 /* Adapter convention H=Edot-v.F+P.Vdot. Signed reversible subcell work is
  * included; this is NOT the adjoint model's positive entropy-production H. */
 for(int i=0;i<t->cells;i++)out[4*i+3]+=cell[i].p*out[4*i+2]-cell[i].u[0]*out[4*i]-cell[i].u[1]*out[4*i+1];
#endif
 return 1;
}
/* Same wall speed must already have prepared the trace limiter and must
 * come from the geometric volume derivative. Ledger is total wall power. */
static inline int gfs_material_rhs_moving_walls2d(const GfsTraceCsr2D*t,
 const GfsFvPrimitive*cell,double gamma,const double*alpha3,
 const double*wall_normal,double*out,double*boundary_power){
 if(!wall_normal||!boundary_power)return 0;
 return gfs_material_rhs_shifted_owner_boundary2d(t,cell,gamma,alpha3,out,
   NULL,NULL,NULL,wall_normal,boundary_power);
}
static inline int gfs_material_rhs_shifted_owner2d(const GfsTraceCsr2D*t,
 const GfsFvPrimitive*cell,double gamma,const double*alpha3,double*out,
 const double*pressure_shift,double*source_power,const double*owner_shift){
 return gfs_material_rhs_shifted_owner_boundary2d(t,cell,gamma,alpha3,out,
   pressure_shift,source_power,owner_shift,NULL,NULL);
}
/* Owner shifts act AFTER reconstruction; distinct from entry wall/gravity offsets.
 * Opt-in direct FV only. The caller must solve and certify geometric constraints. */
static inline int gfs_material_direct_owner_rhs2d(const GfsTraceCsr2D *t,
 const GfsFvPrimitive *cell,double gamma,const double *alpha3,double *out,const double *owner_shift){
 return gfs_material_rhs_shifted_owner2d(t,cell,gamma,alpha3,out,NULL,NULL,owner_shift);
}
static inline int gfs_material_rhs_shifted2d(const GfsTraceCsr2D *t,
 const GfsFvPrimitive *cell,double gamma,const double *alpha3,double *out,
 const double *pressure_shift,double *source_power){
 return gfs_material_rhs_shifted_owner2d(t,cell,gamma,alpha3,out,pressure_shift,source_power,NULL);
}
static inline int gfs_material_rhs_maps2d(const GfsTraceCsr2D *t,
 const GfsFvPrimitive *cell,double gamma,const double *alpha3,double *out){
 return gfs_material_rhs_shifted2d(t,cell,gamma,alpha3,out,NULL,NULL);
}
static inline int gfs_material_rhs2d(const GfsTraceCsr2D *t,
 const GfsFvPrimitive *cell,double gamma,double *out){
 return gfs_material_rhs_maps2d(t,cell,gamma,NULL,out);
}
static inline void gfs_material_weight_matvec(int n,int nf,
 const GfsTraceDirected2D *edge,const double *coef,const double *x,double *y){
 memset(y,0,(size_t)n*sizeof(double));
 for(int k=0;k<nf;k++)y[edge[k].i]+=coef[k]*(x[edge[k].i]-x[edge[k].j]);
}
/* Pack the repeatedly applied graph separately from face geometry. Preserve
 * directed edge order and arithmetic, including the reference summation order. */
typedef struct { int i,j; double coefficient; } GfsMaterialWeightEdge;
static inline void gfs_material_weight_packed_matvec(int n,int nf,
 const GfsMaterialWeightEdge *edge,const double *x,double *y){
 memset(y,0,(size_t)n*sizeof(double));
 for(int k=0;k<nf;k++)y[edge[k].i]+=edge[k].coefficient*(x[edge[k].i]-x[edge[k].j]);
}
/* Projected Jacobi-PCG; certify the 1e-10 true residual. The work budget
 * scales with graph size (at least the legacy 256 iterations). No global variables, file I/O, clipping, or physics
 * relaxation. Returns 1 success, 0 failure. Wdot is a kinematic constraint,
 * not the previous canonical pi/inertia evolution. */
static inline int gfs_material_weight_rate_boundaries2d(int n,int nf,
 const GfsTraceDirected2D *edge,const GfsFvPrimitive *cell,const double *rhs4,
 double *q,double *realized,int *iterations,double *relative,
 const unsigned char *reflected_y){
 if(n<1||nf<1||!edge||!cell||!rhs4||!q||!realized||!iterations||!relative)return 0;
 double *buf=(double*)calloc((size_t)n*7,sizeof(double));
 GfsMaterialWeightEdge *packed=malloc((size_t)nf*sizeof(*packed));
 if(!buf||!packed){free(buf);free(packed);return 0;}
 double *rhs=buf,*r=rhs+n,*z=r+n,*p=z+n,*Ap=p+n,*diag=Ap+n,*Bv=diag+n;
 int ok=0;double transport_scale=0;*iterations=0;*relative=1;
 for(int k=0;k<nf;k++){
  const GfsTraceDirected2D *e=edge+k;
  if(e->i<0||e->i>=n||e->j<0||e->j>=n)goto done;
  int wall=reflected_y?reflected_y[k]:0;
  if(wall>1)goto done;
  if(wall){
   /* Equal weights at a real generator and its own reflected image fix
    * their bisector at the physical wall. Reflection reverses normal velocity, so
    * its normal geometric volume flux is identically zero. */
   if(!gfs_material_stationary_wall2d(e))goto done;
   packed[k]=(GfsMaterialWeightEdge){e->i,e->i,0};
   continue;
  }
  double delta[2]={e->dx,e->dy},centroid[2]={e->cx,e->cy},volume_rate;
  GfsMaterialPowerFaceMotion face;
  if(!gfs_material_power_face_motion(2,delta,centroid,e->area,&face)||
     !gfs_material_power_face_motion_rate(2,&face,cell[e->i].u,cell[e->j].u,0,0,&volume_rate))goto done;
  packed[k]=(GfsMaterialWeightEdge){e->i,e->j,face.weight};diag[e->i]+=face.weight;
  Bv[e->i]+=volume_rate;
  transport_scale+=fabs(volume_rate);
 }
 double mean=0,norm=0,sumabs=0;
 for(int i=0;i<n;i++){rhs[i]=rhs4[4*i+2]-Bv[i];mean+=rhs[i];sumabs+=fabs(rhs[i]);q[i]=0;if(!(diag[i]>0))goto done;}
 /* Cancellation must be measured against transported face volume, not only
  * its small divergence. No dimensionful absolute floor is introduced. */
 double compatibility_scale=sumabs+transport_scale;
 if(!isfinite(compatibility_scale)||!isfinite(mean)||fabs(mean)>1024*DBL_EPSILON*compatibility_scale){fprintf(stderr,"[MATERIAL_WEIGHT_COMPAT_FAIL] sum=%.17g sumabs=%.17g tolerance=%.17g\n",mean,sumabs,1024*DBL_EPSILON*compatibility_scale);goto done;}
 mean/=n;
 for(int i=0;i<n;i++){rhs[i]-=mean;r[i]=rhs[i];norm+=r[i]*r[i];z[i]=r[i]/diag[i];}
 if(norm==0){for(int i=0;i<n;i++){if(rhs[i]!=0)goto done;realized[i]=Bv[i];}*relative=0;ok=1;goto done;}
 if(!(norm>0)||!isfinite(norm))goto done;
#ifdef GFS_MATERIAL_CUDA_WEIGHT
 if(!gfs_material_cuda_weight_pcg(n,nf,edge,reflected_y,rhs,diag,q,iterations,relative))goto done;
 ok=1;
#else
 mean=0;for(int i=0;i<n;i++)mean+=z[i]/n;
 double rz=0;for(int i=0;i<n;i++){p[i]=z[i]-mean;rz+=r[i]*p[i];}
 const int max_iterations=n>256?n:256;
 for(int it=0;it<max_iterations;it++){
  gfs_material_weight_packed_matvec(n,nf,packed,p,Ap);double curvature=0;
  for(int i=0;i<n;i++)curvature+=p[i]*Ap[i];if(!(curvature>0)||!isfinite(curvature))goto done;
  double alpha=rz/curvature,res=0;for(int i=0;i<n;i++){q[i]+=alpha*p[i];r[i]-=alpha*Ap[i];res+=r[i]*r[i];}
  *iterations=it+1;
  if(sqrt(res/norm)<=1e-10){
   /* Recursive residual can drift from rhs-Lq. Certify before accepting;
    * otherwise replace it and restart PCG without weakening tolerance. */
   gfs_material_weight_packed_matvec(n,nf,packed,q,Ap);res=0;
   for(int i=0;i<n;i++){r[i]=rhs[i]-Ap[i];res+=r[i]*r[i];}
   if(sqrt(res/norm)<=1e-10){ok=1;break;}
   mean=0;for(int i=0;i<n;i++){z[i]=r[i]/diag[i];mean+=z[i]/n;}
   rz=0;for(int i=0;i<n;i++){p[i]=z[i]-mean;rz+=r[i]*p[i];}
   continue;
  }
  mean=0;for(int i=0;i<n;i++){z[i]=r[i]/diag[i];mean+=z[i]/n;}
  double next=0;for(int i=0;i<n;i++){z[i]-=mean;next+=r[i]*z[i];}
  double beta=next/rz;rz=next;for(int i=0;i<n;i++)p[i]=z[i]+beta*p[i];
 }
#endif
 mean=0;for(int i=0;i<n;i++)mean+=q[i]/n;for(int i=0;i<n;i++)q[i]-=mean;
 gfs_material_weight_packed_matvec(n,nf,packed,q,Ap);double res=0;
 for(int i=0;i<n;i++){double error=Ap[i]-rhs[i];res+=error*error;realized[i]=Bv[i]+Ap[i];}
 *relative=sqrt(res/norm);if(!isfinite(*relative)||*relative>1e-10)ok=0;
 done:free(packed);free(buf);return ok;
}
/* Preserve the periodic entry point and arithmetic for existing callers. */
static inline int gfs_material_weight_rate2d(int n,int nf,
 const GfsTraceDirected2D *edge,const GfsFvPrimitive *cell,const double *rhs4,
 double *q,double *realized,int *iterations,double *relative){
 return gfs_material_weight_rate_boundaries2d(n,nf,edge,cell,rhs4,q,realized,iterations,relative,NULL);
}
#endif
