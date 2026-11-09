#ifndef GFS_MATERIAL_CONSTRAINT_PROJECT2D_H
#define GFS_MATERIAL_CONSTRAINT_PROJECT2D_H
#include <stdio.h>
#include <time.h>
#ifdef GFS_MATERIAL_BLAS_SCHUR
#include <mkl_cblas.h>
#endif
static inline double gfs_active_profile_clock(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static inline void gfs_active_profile_mark(int enabled,double *last,double *bucket){if(enabled){double now=gfs_active_profile_clock();*bucket+=now-*last;*last=now;}}
#ifdef GFS_MATERIAL_CUDA_ACTIONS
#include "gfs_material_cuda_rhs.h"
#ifndef GFS_MATERIAL_DIRECT_BASIS
#error CUDA actions require batched direct basis
#endif
#endif
#include "gfs_material_constraint_ops2d.h"
typedef struct {int i,j;double sign,dx,dy;} GfsMaterialActiveFace2D;
#ifdef GFS_MATERIAL_DIRECT_BASIS
#include "gfs_material_direct_basis2d.h"
#endif
/* Select current neighbor pairs whose explicit predictor would violate the
 * existing geometric bound. No fitted near-boundary threshold. New neighbors
 * formed during the step must still be checked by the native tessellator.
 * Returns selected count, or -1 for invalid input/capacity. */
static inline int gfs_material_predict_constraints_support2d(int n,int nd,
 const GfsTraceDirected2D*edges,const GfsFvPrimitive*cell,const double*W,
 const double*q,double dt,double beta,int capacity,GfsMaterialActiveFace2D*out,const GfsTraceDirected2D*support){
 if(n<1||nd<1||!edges||!cell||!W||!q||!(dt>0)||!isfinite(dt)||
    !(beta>0&&beta<1)||capacity<0||!out)return -1;
 int m=0;
 for(int k=0;k<nd;k++){
  const GfsTraceDirected2D*e=edges+k;int i=e->i,j=e->j;
  if(i<0||i>=n||j<0||j>=n)return -1;if(i>=j)continue;
  double dx=e->dx+dt*(cell[j].u[0]-cell[i].u[0]);
  double dy=e->dy+dt*(cell[j].u[1]-cell[i].u[1]);
  double dw=W[i]-W[j]+dt*(q[i]-q[j]),d2=dx*dx+dy*dy;
  double sx=e->dx,sy=e->dy;
  if(support){
   sx=support[k].dx;sy=support[k].dy;
   /* Tangent lower bound to squared separation: affine in endpoint x.
    * Freeze support over the whole SSP step, including rejected retries. */
   d2=2*(sx*dx+sy*dy)-sx*sx-sy*sy;
  }
  if(!(d2>0)||!isfinite(d2)||!isfinite(dw))return -1;
  if(fabs(dw)>beta*d2){
   if(m==capacity)return -1;
   out[m++]=(GfsMaterialActiveFace2D){i,j,dw>0?1:-1,sx,sy};
  }
 }return m;
}
static inline int gfs_material_predict_constraints2d(int n,int nd,
 const GfsTraceDirected2D*edges,const GfsFvPrimitive*cell,const double*W,
 const double*q,double dt,double beta,int capacity,GfsMaterialActiveFace2D*out){
 return gfs_material_predict_constraints_support2d(n,nd,edges,cell,W,q,dt,beta,capacity,out,NULL);
}
/* Small SPD active system. No regularization of singular constraints. */
static inline int gfs_material_small_spd(int m,const double*S,const double*b,double*x){
 double*L=calloc((size_t)m*m,sizeof(double));if(!L)return 0;int ok=0;
 for(int i=0;i<m;i++)for(int j=0;j<=i;j++){
  double a=S[i*m+j];for(int k=0;k<j;k++)a-=L[i*m+k]*L[j*m+k];
  if(i==j){if(!(a>64*DBL_EPSILON*fabs(S[i*m+i]))||!isfinite(a)){fprintf(stderr,"[MATERIAL_SPD_FAIL] phase=pivot count=%d index=%d pivot=%.17g diagonal=%.17g\n",m,i,a,S[i*m+i]);goto done;}L[i*m+j]=sqrt(a);}
  else L[i*m+j]=a/L[j*m+j];
 }
 for(int i=0;i<m;i++){double a=b[i];for(int k=0;k<i;k++)a-=L[i*m+k]*x[k];x[i]=a/L[i*m+i];}
 for(int i=m-1;i>=0;i--){double a=x[i];for(int k=i+1;k<m;k++)a-=L[k*m+i]*x[k];x[i]=a/L[i*m+i];}
 double err=0,scale=0;for(int i=0;i<m;i++){double a=-b[i];scale+=fabs(b[i]);for(int j=0;j<m;j++){a+=S[i*m+j]*x[j];scale+=fabs(S[i*m+j]*x[j]);}err+=fabs(a);}
 ok=isfinite(err)&&err<=1e-10*scale;
 if(!ok)fprintf(stderr,"[MATERIAL_SPD_FAIL] phase=residual count=%d error=%.17g scale=%.17g\n",m,err,scale);
 done:free(L);return ok;
}
#include "gfs_material_dependent_release.h"
/* Per-RHS workspace ONLY: destroy before any state, geometry, boundary metadata, trace maps, beta or gamma
 * change, including an RK stage or retry. Never persist across time steps.
 * Reuse only a verified unchanged active prefix during union expansion. */
#ifdef GFS_MATERIAL_REUSE_RESPONSE
typedef struct {int n,nf,setup_ready,rhs_ready;GfsMaterialFaceResponse*face;double*DP,*rhs;} GfsMaterialFrozenResponse2D;
static inline void gfs_material_frozen_clear2d(GfsMaterialFrozenResponse2D*f){if(f){free(f->face);free(f->DP);free(f->rhs);free(f);}}
#endif
typedef struct {
 int n,m,iterations;
 GfsMaterialActiveFace2D *active;
 double *G,*DG,*Q,*S;
#ifdef GFS_MATERIAL_REUSE_RESPONSE
 GfsMaterialFrozenResponse2D*frozen;
#endif
#ifdef GFS_MATERIAL_REUSE_FACTOR
 GfsMaterialDirectFactor2D factor;
#endif
} GfsMaterialConstraintReuse2D;
static inline void gfs_material_constraint_reuse_clear2d(GfsMaterialConstraintReuse2D*r){
 if(!r)return;
#ifdef GFS_MATERIAL_REUSE_FACTOR
 gfs_material_direct_factor_clear2d(&r->factor);
#endif
 #ifdef GFS_MATERIAL_REUSE_RESPONSE
 gfs_material_frozen_clear2d(r->frozen);
#endif
 free(r->active);free(r->G);free(r->DG);free(r->Q);free(r->S);memset(r,0,sizeof(*r));
}
/* Explicit active-set candidate. Caller supplies the geometrically selected
 * inequalities; negative multipliers fail rather than silently forcing them.
 * Output layout per cell is Fx,Fy,Vtarget,H,Wdot,Vrealized,Udot.
 * No time integration, state clipping or Python dependency. */
static inline int gfs_material_project_active_reuse_drive_maps2d_status(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const GfsFvPrimitive*cell,double gamma,
 int m,const GfsMaterialActiveFace2D*active,double beta,double*out7,int*max_iterations,GfsMaterialConstraintReuse2D*reuse,const double*alpha3,const unsigned char*reflected_y,const double*drive,const double*extra_force,int *release_index){
 if(release_index)*release_index=-1;
 if(!t||nd<1||!edge||!cell||m<1||!active||!(beta>0&&beta<1)||!out7||!max_iterations)return 0;
 int n=t->cells,ok=0;*max_iterations=0;
 int profile=getenv("GFS_MATERIAL_ACTIVE_DETAIL")!=NULL;double stamps[6]={0},last=profile?gfs_active_profile_clock():0;int profile_phase=0;
 GfsMaterialFaceResponse*cache=malloc((size_t)t->faces*sizeof(*cache));GfsFvPrimitive*zero=calloc(n,sizeof(*zero));
 double*G=calloc((size_t)m*n,sizeof(double)),*DG=calloc((size_t)m*n,sizeof(double)),*Q=calloc((size_t)m*2*n,sizeof(double));
 double*S=calloc((size_t)m*m,sizeof(double)),*b=calloc((size_t)m*4,sizeof(double));
 int*gram_map=NULL;
#ifdef GFS_MATERIAL_REUSE_GRAM
 gram_map=malloc((size_t)m*sizeof(int));
#endif
 double*work=calloc((size_t)n*15,sizeof(double));
 if(!cache||!zero||!G||!DG||!Q||!S||!b||!work)goto done;
 double*rhs=work,*realized=rhs+4*n,*ct=realized+n,*bt=ct+2*n,*P=bt+2*n,*peff=P+n,*heat=peff+n,*q=heat+n,*DP=q+n,*dummy=DP+n;
 double*bp=b+m,*lambda=bp+m,*lp=lambda+m;
 for(int i=0;i<n;i++){P[i]=cell[i].p;peff[i]=drive?drive[i]:P[i];if(!isfinite(peff[i]))goto done;}
#ifdef GFS_MATERIAL_REUSE_RESPONSE
 /* Owned by the caller's immutable RHS workspace; never cross a stage. */
 GfsMaterialFrozenResponse2D*frozen=reuse?reuse->frozen:NULL;
 if(reuse&&!frozen){
  frozen=calloc(1,sizeof(*frozen));
  if(frozen){frozen->n=n;frozen->nf=t->faces;frozen->face=malloc((size_t)t->faces*sizeof(*cache));frozen->DP=malloc((size_t)n*sizeof(double));frozen->rhs=malloc((size_t)4*n*sizeof(double));
   if(!frozen->face||!frozen->DP||!frozen->rhs){gfs_material_frozen_clear2d(frozen);frozen=NULL;}
  }reuse->frozen=frozen;
 }
 if(frozen&&(frozen->n!=n||frozen->nf!=t->faces))goto done;
 if(frozen&&frozen->setup_ready){memcpy(cache,frozen->face,(size_t)t->faces*sizeof(*cache));memcpy(DP,frozen->DP,(size_t)n*sizeof(double));if(profile)fprintf(stderr,"[MATERIAL_RESPONSE_REUSE] setup=1\n");}
 else
#endif
 {
  if(!gfs_material_constraint_cache_maps2d(t,cell,gamma,alpha3,cache))goto done;
  if(!gfs_material_constraint_actions_maps2d(t,cache,alpha3,peff,DP,ct))goto done;
#ifdef GFS_MATERIAL_REUSE_RESPONSE
  if(frozen){memcpy(frozen->face,cache,(size_t)t->faces*sizeof(*cache));memcpy(frozen->DP,DP,(size_t)n*sizeof(double));frozen->setup_ready=1;}
#endif
 }
 int cached=0;
#ifdef GFS_MATERIAL_REUSE_RELEASE
 /* Same RHS only: match retained rows after removal/reordering. All appended
  * rows follow this matched prefix; unmatched suffix is recomputed. */
 if(reuse&&reuse->n==n){
  for(int k=0;k<m;k++){
   int found=-1;
   for(int j=0;j<reuse->m;j++){
    const GfsMaterialActiveFace2D*a=active+k,*b=reuse->active+j;
    if(a->i==b->i&&a->j==b->j&&a->sign==b->sign&&a->dx==b->dx&&a->dy==b->dy){found=j;break;}
   }
   if(found<0)break;
   if(gram_map)gram_map[k]=found;
   memcpy(G+(size_t)k*n,reuse->G+(size_t)found*n,(size_t)n*sizeof(double));
   memcpy(DG+(size_t)k*n,reuse->DG+(size_t)found*n,(size_t)n*sizeof(double));
   memcpy(Q+(size_t)k*2*n,reuse->Q+(size_t)found*2*n,(size_t)2*n*sizeof(double));
   cached++;
  }
  if(cached)*max_iterations=reuse->iterations;
 }
#else
 if(reuse&&reuse->n==n&&reuse->m<=m){
  cached=reuse->m;
  for(int k=0;k<cached;k++){
   const GfsMaterialActiveFace2D*a=active+k,*b=reuse->active+k;
   if(a->i!=b->i||a->j!=b->j||a->sign!=b->sign||a->dx!=b->dx||a->dy!=b->dy){cached=0;break;}
  }
  if(cached){
   if(gram_map)for(int k=0;k<cached;k++)gram_map[k]=k;
   memcpy(G,reuse->G,(size_t)cached*n*sizeof(double));
   memcpy(DG,reuse->DG,(size_t)cached*n*sizeof(double));
   memcpy(Q,reuse->Q,(size_t)cached*2*n*sizeof(double));
   *max_iterations=reuse->iterations;
  }
 }
#endif
gfs_active_profile_mark(profile,&last,&stamps[0]);profile_phase=1;
#ifdef GFS_MATERIAL_DIRECT_BASIS
#ifdef GFS_MATERIAL_REUSE_FACTOR
 if(cached<m&&!gfs_material_direct_basis_cached2d(n,nd,edge,reflected_y,m-cached,active+cached,G+(size_t)cached*n,reuse?&reuse->factor:NULL))goto done;
#else
 if(cached<m&&!gfs_material_direct_basis2d(n,nd,edge,reflected_y,m-cached,active+cached,G+(size_t)cached*n))goto done;
#endif
#endif
gfs_active_profile_mark(profile,&last,&stamps[1]);profile_phase=2;
#ifdef GFS_MATERIAL_CUDA_ACTIONS
 if(cached<m&&!gfs_material_cuda_actions2d(t,cache,alpha3,m-cached,G+(size_t)cached*n,DG+(size_t)cached*n,Q+(size_t)cached*2*n))goto done;
#endif
 for(int k=0;k<m;k++){
  const GfsMaterialActiveFace2D*a=active+k;
  if(a->i<0||a->i>=n||a->j<0||a->j>=n||a->i==a->j||fabs(a->sign)!=1)goto done;
  memset(rhs,0,(size_t)n*4*sizeof(double));rhs[4*a->i+2]=a->sign;rhs[4*a->j+2]=-a->sign;
  int iterations=0;double rel=0;double*g=G+k*n,*dg=DG+k*n,*qq=Q+2*k*n;
  if(k>=cached){
#ifndef GFS_MATERIAL_DIRECT_BASIS
  if(!gfs_material_weight_rate_boundaries2d(n,nd,edge,zero,rhs,g,realized,&iterations,&rel,reflected_y))goto done;
#endif
  if(iterations>*max_iterations)*max_iterations=iterations;
  if(!gfs_material_geometry_transpose_boundaries2d(n,nd,edge,g,bt,reflected_y))goto done;
#ifdef GFS_MATERIAL_CUDA_ACTIONS
  for(int i=0;i<2*n;i++)qq[i]-=bt[i];
#else
  if(!gfs_material_constraint_actions_maps2d(t,cache,alpha3,g,dg,ct))goto done;
  for(int i=0;i<2*n;i++)qq[i]=ct[i]-bt[i];
#endif
  qq[2*a->i]+=2*beta*a->dx;qq[2*a->i+1]+=2*beta*a->dy;
  qq[2*a->j]-=2*beta*a->dx;qq[2*a->j+1]-=2*beta*a->dy;
  }
  for(int i=0;i<n;i++){bp[k]+=g[i]*DP[i];b[k]+=qq[2*i]*cell[i].u[0]+qq[2*i+1]*cell[i].u[1];}
  b[k]+=bp[k];
 }
 gfs_active_profile_mark(profile,&last,&stamps[2]);profile_phase=3;
#ifdef GFS_MATERIAL_REUSE_GRAM
 if(gram_map&&reuse&&reuse->S&&cached==m){
  for(int i=0;i<m;i++)for(int j=0;j<m;j++)S[i*m+j]=reuse->S[gram_map[i]*reuse->m+gram_map[j]];
  if(profile)fprintf(stderr,"[MATERIAL_GRAM_REUSE] rows=%d previous_rows=%d\n",m,reuse->m);
 }else
#endif
 {
#ifdef GFS_MATERIAL_BLAS_SCHUR
 /* Same Gram matrix, batched BLAS; keep the original symmetry/SPD checks. */
 cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,m,m,n,1.0,G,n,DG,n,0.0,S,m);
#else
 for(int i=0;i<m;i++)for(int j=0;j<m;j++)for(int k=0;k<n;k++)S[i*m+j]+=G[i*n+k]*DG[j*n+k];
#endif
 }
 /* Optional first real Gram operands, for backend precision qualification. */
 if(getenv("GFS_MATERIAL_EXPORT_GRAM")){
  static int exported=0;
  if(!exported){const char*path=getenv("GFS_MATERIAL_EXPORT_GRAM");FILE*f=fopen(path,"wx");
   if(f){int shape[2]={m,n};size_t count=0;count+=fwrite(shape,sizeof(int),2,f);count+=fwrite(G,sizeof(double),(size_t)m*n,f);count+=fwrite(DG,sizeof(double),(size_t)m*n,f);count+=fwrite(S,sizeof(double),(size_t)m*m,f);int close_status=fclose(f);
    if(count!=2+(size_t)2*m*n+(size_t)m*m||close_status)fprintf(stderr,"[MATERIAL_GRAM_EXPORT_FAIL] path=%s\n",path);
    else exported=1;
   }
  }
 }
 /* Symmetry is required, not repaired. */
 for(int i=0;i<m;i++)for(int j=0;j<i;j++)if(fabs(S[i*m+j]-S[j*m+i])>1e-10*sqrt(fabs(S[i*m+i]*S[j*m+j])))goto done;
 if(!gfs_material_small_spd(m,S,b,lambda)||!gfs_material_small_spd(m,S,bp,lp)){fprintf(stderr,"[MATERIAL_ACTIVE_FAIL] phase=spd count=%d\n",m);
  if(release_index){
   int dependent=gfs_material_dependent_release(m,S,b,bp,gfs_material_small_spd);
   if(dependent>=0){*release_index=dependent;
    fprintf(stderr,"[MATERIAL_DEPENDENT_RELEASE] count=%d index=%d certified_all_inequalities=1\n",m,dependent);goto done;}
  }
  if(getenv("GFS_MATERIAL_EXPORT_SPD")){
   FILE*f=fopen("failed_active_system.txt","wx");
   if(f){fprintf(f,"%d\n",m);for(int i=0;i<m;i++){
    fprintf(f,"%d %d %.17g %.17g %.17g",active[i].i,active[i].j,active[i].sign,b[i],bp[i]);
    for(int j=0;j<m;j++)fprintf(f," %.17g",S[i*m+j]);fprintf(f,"\n");}fclose(f);}
  }goto done;}
 if(getenv("GFS_MATERIAL_ACTIVE_WORK")){
  for(int k=0;k<m;k++)fprintf(stderr,"[MATERIAL_ACTIVE_WORK] row=%d i=%d j=%d b=%.17g bp=%.17g lambda=%.17g lp=%.17g\n",k,active[k].i,active[k].j,b[k],bp[k],lambda[k],lp[k]);
 }
 int worst=-1;
 for(int k=0;k<m;k++){
  if(!isfinite(lambda[k]))goto done;
  if(lambda[k]<0&&(worst<0||lambda[k]<lambda[worst]))worst=k;
 }
 if(worst>=0){
  if(release_index)*release_index=worst;
  fprintf(stderr,"[MATERIAL_ACTIVE_FAIL] phase=multiplier count=%d index=%d value=%.17g\n",m,worst,lambda[worst]);goto done;
 }
 gfs_active_profile_mark(profile,&last,&stamps[3]);profile_phase=4;
 for(int k=0;k<m;k++)for(int i=0;i<n;i++)peff[i]-=G[k*n+i]*lp[k];
#ifdef GFS_MATERIAL_REUSE_RESPONSE
 if(frozen&&frozen->rhs_ready){memcpy(rhs,frozen->rhs,(size_t)4*n*sizeof(double));if(profile)fprintf(stderr,"[MATERIAL_RESPONSE_REUSE] rhs=1\n");}
 else
#endif
 {
 if(drive){if(!gfs_material_cached_drive_rhs2d(t,cache,cell,alpha3,drive,rhs))goto done;}
 else if(gfs_material_rhs_maps2d(t,cell,gamma,alpha3,rhs)!=1)goto done;
#ifdef GFS_MATERIAL_REUSE_RESPONSE
  if(frozen){memcpy(frozen->rhs,rhs,(size_t)4*n*sizeof(double));frozen->rhs_ready=1;}
#endif
 }
 if(!gfs_material_projected_heat_maps2d(t,cache,cell,alpha3,peff,heat))goto done;
 if(extra_force)for(int i=0;i<n;i++){
  if(!isfinite(extra_force[2*i])||!isfinite(extra_force[2*i+1]))goto done;
  rhs[4*i]+=extra_force[2*i];rhs[4*i+1]+=extra_force[2*i+1];
 }
 for(int i=0;i<n;i++){
  for(int k=0;k<m;k++){rhs[4*i]-=Q[2*k*n+2*i]*lp[k];rhs[4*i+1]-=Q[2*k*n+2*i+1]*lp[k];rhs[4*i+2]-=DG[k*n+i]*lambda[k];}
  rhs[4*i+3]=heat[i];
 }
 gfs_active_profile_mark(profile,&last,&stamps[4]);profile_phase=5;
 int iterations=0;double rel=0;
 if(!gfs_material_weight_rate_boundaries2d(n,nd,edge,cell,rhs,q,realized,&iterations,&rel,reflected_y))goto done;
 if(iterations>*max_iterations)*max_iterations=iterations;
 for(int i=0;i<n;i++){double*r=out7+7*i;for(int j=0;j<4;j++)r[j]=rhs[4*i+j];r[4]=q[i];r[5]=realized[i];r[6]=-P[i]*realized[i]+heat[i];}
 ok=1;
#ifdef GFS_MATERIAL_REUSE_RELEASE
 done:
 if(reuse&&(ok||(release_index&&*release_index>=0))){
#else
 if(reuse){
#endif
  GfsMaterialActiveFace2D*saved=malloc((size_t)m*sizeof(*saved));
  if(saved){
#ifdef GFS_MATERIAL_REUSE_FACTOR
   GfsMaterialDirectFactor2D kept_factor=reuse->factor;memset(&reuse->factor,0,sizeof(reuse->factor));
#endif
#ifdef GFS_MATERIAL_REUSE_RESPONSE
   GfsMaterialFrozenResponse2D*kept_frozen=reuse->frozen;reuse->frozen=NULL;
#endif
   memcpy(saved,active,(size_t)m*sizeof(*saved));gfs_material_constraint_reuse_clear2d(reuse);
#ifdef GFS_MATERIAL_REUSE_RESPONSE
   reuse->frozen=kept_frozen;
#endif
#ifdef GFS_MATERIAL_REUSE_FACTOR
   reuse->factor=kept_factor;
#endif
   reuse->n=n;reuse->m=m;reuse->iterations=*max_iterations;reuse->active=saved;
   reuse->G=G;reuse->DG=DG;reuse->Q=Q;G=DG=Q=NULL;
#ifdef GFS_MATERIAL_REUSE_GRAM
   reuse->S=S;S=NULL;
#endif
  }
 }
#ifndef GFS_MATERIAL_REUSE_RELEASE
 done:
#endif
 if(profile){gfs_active_profile_mark(profile,&last,&stamps[profile_phase]);fprintf(stderr,"[MATERIAL_ACTIVE_DETAIL] rows=%d ok=%d setup_ms=%.9g basis_ms=%.9g actions_ms=%.9g schur_ms=%.9g physical_ms=%.9g final_solve_ms=%.9g\n",m,ok,1000*stamps[0],1000*stamps[1],1000*stamps[2],1000*stamps[3],1000*stamps[4],1000*stamps[5]);}
 free(cache);free(zero);free(G);free(DG);free(Q);free(S);free(b);free(work);free(gram_map);return ok;
}
static inline int gfs_material_project_active_reuse_drive_maps2d(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const GfsFvPrimitive*cell,double gamma,
 int m,const GfsMaterialActiveFace2D*active,double beta,double*out7,int*max_iterations,GfsMaterialConstraintReuse2D*reuse,const double*alpha3,const unsigned char*reflected_y,const double*drive,const double*extra_force){
 return gfs_material_project_active_reuse_drive_maps2d_status(t,nd,edge,cell,gamma,m,active,beta,out7,max_iterations,reuse,alpha3,reflected_y,drive,extra_force,NULL);
}
/* Generalized drive requires matching external-potential work in the caller;
 * it is not an EOS pressure and does not replace physical sound speeds. */
static inline int gfs_material_project_active_reuse_boundaries_maps2d(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const GfsFvPrimitive*cell,double gamma,
 int m,const GfsMaterialActiveFace2D*active,double beta,double*out7,int*max_iterations,GfsMaterialConstraintReuse2D*reuse,const double*alpha3,const unsigned char*reflected_y){
 return gfs_material_project_active_reuse_drive_maps2d(t,nd,edge,cell,gamma,m,active,beta,out7,max_iterations,reuse,alpha3,reflected_y,NULL,NULL);
}
/* Stationary slip walls with even pressure ghosts only. Affine hydrostatic
 * pressure offsets need their own force/work treatment and are not implied. */
static inline int gfs_material_project_active_reuse_maps2d(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const GfsFvPrimitive*cell,double gamma,
 int m,const GfsMaterialActiveFace2D*active,double beta,double*out7,int*max_iterations,GfsMaterialConstraintReuse2D*reuse,const double*alpha3){
 return gfs_material_project_active_reuse_boundaries_maps2d(t,nd,edge,cell,gamma,m,active,beta,out7,max_iterations,reuse,alpha3,NULL);
}
static inline int gfs_material_project_active_reuse2d(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const GfsFvPrimitive*cell,double gamma,
 int m,const GfsMaterialActiveFace2D*active,double beta,double*out7,int*max_iterations,GfsMaterialConstraintReuse2D*reuse){
 return gfs_material_project_active_reuse_maps2d(t,nd,edge,cell,gamma,m,active,beta,out7,max_iterations,reuse,NULL);
}
static inline int gfs_material_project_active2d(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const GfsFvPrimitive*cell,double gamma,
 int m,const GfsMaterialActiveFace2D*active,double beta,double*out7,int*max_iterations){
 return gfs_material_project_active_reuse2d(t,nd,edge,cell,gamma,m,active,beta,out7,max_iterations,NULL);
}
#endif
