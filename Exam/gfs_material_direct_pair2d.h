#ifndef GFS_MATERIAL_DIRECT_PAIR2D_H
#define GFS_MATERIAL_DIRECT_PAIR2D_H
#ifndef GFS_MATERIAL_DIRECT_FV
#error Direct pair pressure correction requires direct FV
#endif
#include "gfs_material_constraint_project2d.h"
#include "gfs_material_pressure_jacobian.h"
#ifdef GFS_MATERIAL_STAR_ENTROPY
#include "gfs_material_star_entropy.h"
#endif
#include "gfs_material_small_general.h"
/* Experimental periodic frozen-active-set solve. Caller rechecks ALL bounds
 * and retries failed steps. No complementarity or entropy proof is asserted.
 * pressure shifts are endpoint incidence times pair amplitudes, not L^-1 A.
 * out7: force x/y, Vdot, signed excess work, Wdot, realized Vdot, reserved.
 */
static inline int gfs_material_direct_pair2d(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const GfsFvPrimitive*cell,double gamma,
 const double*alpha,int m,const GfsMaterialActiveFace2D*active,double beta,
 const double*W,double dt,const GfsTraceDirected2D*support,double*out7,int*iterations){
 if(!t||m<1||!active||!out7||!iterations||!W||!(dt>0))return 0;
 int n=t->cells,nf=t->faces,ok=0;*iterations=0;
 for(int j=0;j<t->entries;j++)if(t->entry[j].pad)return 0;
 double*G=calloc((size_t)m*n,sizeof(double)),*buf=calloc((size_t)n*12,sizeof(double));
 double*S=calloc((size_t)m*m,sizeof(double)),*r=calloc(m,sizeof(double)),*dl=calloc(m,sizeof(double));
 double*lambda=calloc(m,sizeof(double)),*trial=calloc(m,sizeof(double));
 GfsMaterialActiveFace2D*certificate=malloc((size_t)nd*sizeof(*certificate));
 GfsFvPrimitive*still=calloc(n,sizeof(*still));GfsFvPrimitive*face=calloc((size_t)2*nf,sizeof(*face));
 if(!G||!buf||!S||!r||!dl||!lambda||!trial||!still||!face||!certificate)goto done;
 {double*rhs=buf,*q=rhs+4*n,*real=q+n,*chi=real+n,*candidate=chi+n;int its=0;double rel=0;
#ifdef GFS_MATERIAL_DIRECT_PAIR_BATCH_BASIS
#ifndef GFS_MATERIAL_DIRECT_BASIS
#error Pair batched basis requires existing direct graph backend
#endif
 if(!gfs_material_direct_basis2d(n,nd,edge,NULL,m,active,G))goto done;
#else
 for(int a=0;a<m;a++){
  memset(rhs,0,(size_t)4*n*sizeof(double));rhs[4*active[a].i+2]=active[a].sign;rhs[4*active[a].j+2]=-active[a].sign;
  if(!gfs_material_weight_rate_boundaries2d(n,nd,edge,still,rhs,G+(size_t)a*n,real,&its,&rel,NULL))goto done;*iterations+=its;
 }
#endif
 for(int k=0;k<nf;k++){
  const GfsTraceFace2D*f=t->face+k;int start[2]={f->sl,f->sr},len[2]={f->nl,f->nr};
  for(int side=0;side<2;side++)for(int h=start[side];h<start[side]+len[side];h++){
   const GfsTraceEntry*e=t->entry+h;int owner=side?f->j:f->i;GfsMaterialTraceWeights w=gfs_material_entry_weights(e,owner,alpha?alpha+3*owner:NULL);GfsFvPrimitive*z=face+2*k+side;
   z->rho+=w.rho*cell[e->cell].rho;z->p+=w.pressure*cell[e->cell].p;for(int d=0;d<2;d++)z->u[d]+=w.velocity*cell[e->cell].u[d];
  }
 }
 if(gfs_material_direct_owner_rhs2d(t,cell,gamma,alpha,rhs,chi)!=1||!gfs_material_weight_rate_boundaries2d(n,nd,edge,cell,rhs,q,real,&its,&rel,NULL))goto done;
 for(int step=0;step<8;step++){
  double worst=-INFINITY;
  for(int a=0;a<m;a++){int i=active[a].i,j=active[a].j;r[a]=active[a].sign*(q[i]-q[j])-2*beta*(active[a].dx*(cell[j].u[0]-cell[i].u[0])+active[a].dy*(cell[j].u[1]-cell[i].u[1]));double roundoff=128*DBL_EPSILON*(fabs(q[i])+fabs(q[j])+2*beta*(fabs(active[a].dx*(cell[j].u[0]-cell[i].u[0]))+fabs(active[a].dy*(cell[j].u[1]-cell[i].u[1]))));worst=fmax(worst,r[a]-roundoff);}
  /* The physical acceptance criterion is the unchanged whole-step geometric
   * predictor. A tangent equality need not be resolved below roundoff once
   * every actual finite-step inequality is already satisfied. */
  int feasible=gfs_material_predict_constraints_support2d(n,nd,edge,cell,W,q,dt,beta,nd,certificate,support);
  if(feasible<0)goto done;
  if(worst<=0||feasible==0){ok=1;break;}
  memset(S,0,(size_t)m*m*sizeof(double));
  for(int k=0;k<nf;k++){
   const GfsTraceFace2D*f=t->face+k;GfsFvPrimitive l=face[2*k],rr=face[2*k+1];l.p-=chi[f->i];rr.p-=chi[f->j];double normal[3]={f->nx,f->ny,0};GfsMaterialFaceResponse fr;GfsMaterialPressureJacobian jac;
   #ifdef GFS_MATERIAL_STAR_ENTROPY
   double entropy_theta;
   if(!gfs_material_star_entropy_jacobian(2,cell[f->i],cell[f->j],l,rr,normal,gamma,&fr,&jac,&entropy_theta))goto done;
   #else
   if(gfs_material_face_response(2,l,rr,normal,gamma,&fr)!=1||!gfs_material_pressure_jacobian(l.rho,l.p,rr.rho,rr.p,fr.pstar,gamma,&jac))goto done;
   #endif
   for(int b=0;b<m;b++){
    double left=active[b].sign*((f->i==active[b].i)-(f->i==active[b].j));double right=active[b].sign*((f->j==active[b].i)-(f->j==active[b].j));double change=f->area*(jac.du_left*left+jac.du_right*right);
    if(change==0)continue;for(int a=0;a<m;a++)S[a*m+b]+=(G[(size_t)a*n+f->i]-G[(size_t)a*n+f->j])*change;
   }
  }
  if(!gfs_material_small_general(m,S,r,dl)){fprintf(stderr,"[MATERIAL_PAIR_FAIL] phase=small_system active=%d step=%d\n",m,step);goto done;}
  int accepted=0;
  for(int back=0;back<16;back++){
   double factor=ldexp(1.,-back);memset(candidate,0,(size_t)n*sizeof(double));
   for(int a=0;a<m;a++){trial[a]=lambda[a]+factor*dl[a];candidate[active[a].i]+=active[a].sign*trial[a];candidate[active[a].j]-=active[a].sign*trial[a];}
   if(gfs_material_direct_owner_rhs2d(t,cell,gamma,alpha,rhs,candidate)!=1)continue;
   if(!gfs_material_weight_rate_boundaries2d(n,nd,edge,cell,rhs,q,real,&its,&rel,NULL))continue;*iterations+=its;
   double next=-INFINITY;for(int a=0;a<m;a++){int i=active[a].i,j=active[a].j;double residual=active[a].sign*(q[i]-q[j])-2*beta*(active[a].dx*(cell[j].u[0]-cell[i].u[0])+active[a].dy*(cell[j].u[1]-cell[i].u[1]));double roundoff=128*DBL_EPSILON*(fabs(q[i])+fabs(q[j])+2*beta*(fabs(active[a].dx*(cell[j].u[0]-cell[i].u[0]))+fabs(active[a].dy*(cell[j].u[1]-cell[i].u[1]))));next=fmax(next,residual-roundoff);}
   int feasible_trial=gfs_material_predict_constraints_support2d(n,nd,edge,cell,W,q,dt,beta,nd,certificate,support);
   if(feasible_trial<0)continue;
   if(next<worst||feasible_trial==0){memcpy(lambda,trial,m*sizeof(double));memcpy(chi,candidate,n*sizeof(double));accepted=1;if(next<=0||feasible_trial==0)ok=1;break;}
  }
  if(!accepted){fprintf(stderr,"[MATERIAL_PAIR_FAIL] phase=line_search active=%d step=%d worst=%.17g\n",m,step,worst);goto done;}if(ok)break;
 }
 if(!ok){double worst=-INFINITY;for(int a=0;a<m;a++){int i=active[a].i,j=active[a].j;double residual=active[a].sign*(q[i]-q[j])-2*beta*(active[a].dx*(cell[j].u[0]-cell[i].u[0])+active[a].dy*(cell[j].u[1]-cell[i].u[1]));worst=fmax(worst,residual);}fprintf(stderr,"[MATERIAL_PAIR_FAIL] phase=iteration_budget active=%d worst=%.17g\n",m,worst);}
 if(ok)for(int i=0;i<n;i++){for(int d=0;d<4;d++)out7[7*i+d]=rhs[4*i+d];out7[7*i+4]=q[i];out7[7*i+5]=real[i];out7[7*i+6]=0;}
 }
 done:free(certificate);free(G);free(buf);free(S);free(r);free(dl);free(lambda);free(trial);free(still);free(face);return ok;
}
#endif
