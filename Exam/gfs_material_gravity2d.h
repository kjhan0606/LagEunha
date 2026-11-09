#ifndef GFS_MATERIAL_GRAVITY2D_H
#define GFS_MATERIAL_GRAVITY2D_H
#include "gfs_material_potential2d.h"
#include "gfs_material_constraint_project2d.h"
/* Experimental instantaneous linear-potential coupling, NOT a time driver.
 * Uses physical rho/P for impedances, P-eta only for the pressure drive.
 * L eta=dEpot/dW; F=C^T(P-eta)-Kv+B^T eta-dEpot/dx.
 * Requires stationary horizontal reflected walls / consistent periodicity.
 * No affine ghost pressure offsets or extra gravity kick. Active constraints
 * use the same effective drive and geometric potential work.
 * Output layout Fx,Fy,Vtarget,H,Wdot,Vrealized,Udot; potential_power is direct
 * gx.v+gW.Wdot. Failure invalidates outputs. Caller owns geometry validation.
 */
static inline int gfs_material_gravity_active_rhs2d_reuse_status(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const unsigned char*reflected_y,
 const GfsFvPrimitive*cell,const double*centroid_relative,const double grad_phi[2],
 double gamma,const double*alpha3,double*out7,double*potential_power,int*iterations,
 int active_count,const GfsMaterialActiveFace2D*active,double beta,int *release_index,GfsMaterialConstraintReuse2D *reuse){
 if(release_index)*release_index=-1;
 if(active_count<0||(active_count&&(!active||!(beta>0&&beta<1))))return 0;
 if(!t||nd<1||!edge||!cell||!centroid_relative||!grad_phi||!out7||!potential_power||!iterations)return 0;
 const char *phase="allocation";
 int n=t->cells,ok=0;*potential_power=0;*iterations=0;
 double*work=calloc((size_t)n*15,sizeof(double));
 GfsFvPrimitive*zero=calloc(n,sizeof(*zero));
 GfsMaterialFaceResponse*cache=malloc((size_t)t->faces*sizeof(*cache));
 if(!work||!zero||!cache)goto done;
 double*gx=work,*gw=gx+2*n,*eta=gw+n,*bt=eta+n;
 double*pe=bt+2*n,*q=pe+n,*real=q+n,*rhs=real+n,*heat=rhs+4*n,*rho=heat+n;
 for(int i=0;i<n;i++)rho[i]=cell[i].rho;
 phase="potential_gradient";
 if(!gfs_material_potential_gradient2d(n,nd,edge,reflected_y,rho,centroid_relative,grad_phi,gx,gw))goto done;
 for(int i=0;i<n;i++)rhs[4*i+2]=gw[i];
 int its=0;double residual=0;
 phase="potential_solve";
 if(!gfs_material_weight_rate_boundaries2d(n,nd,edge,zero,rhs,eta,real,&its,&residual,reflected_y))goto done;
 *iterations=its;
 phase="geometry_transpose";
 if(!gfs_material_geometry_transpose_boundaries2d(n,nd,edge,eta,bt,reflected_y))goto done;
 for(int i=0;i<n;i++)pe[i]=cell[i].p-eta[i];
 if(active_count){
  for(int i=0;i<2*n;i++)bt[i]-=gx[i];
  int projected_iterations=0;
  phase="active_projection";
  if(!gfs_material_project_active_reuse_drive_maps2d_status(t,nd,edge,cell,gamma,active_count,active,beta,out7,&projected_iterations,reuse,alpha3,reflected_y,pe,bt,release_index))goto done;
  *iterations+=projected_iterations;
  for(int i=0;i<n;i++)*potential_power+=gx[2*i]*cell[i].u[0]+gx[2*i+1]*cell[i].u[1]+gw[i]*out7[7*i+4];
  ok=isfinite(*potential_power);goto done;
 }
 phase="face_cache";
 if(!gfs_material_constraint_cache_maps2d(t,cell,gamma,alpha3,cache))goto done;
 phase="pressure_drive";
 if(!gfs_material_cached_drive_rhs2d(t,cache,cell,alpha3,pe,rhs))goto done;
 for(int i=0;i<n;i++){
  rhs[4*i]+=bt[2*i]-gx[2*i];rhs[4*i+1]+=bt[2*i+1]-gx[2*i+1];
  heat[i]=rhs[4*i+3];
 }
 phase="weight_solve";
 if(!gfs_material_weight_rate_boundaries2d(n,nd,edge,cell,rhs,q,real,&its,&residual,reflected_y))goto done;
 *iterations+=its;
 for(int i=0;i<n;i++){
  double*r=out7+7*i;r[0]=rhs[4*i];r[1]=rhs[4*i+1];r[2]=rhs[4*i+2];
  r[3]=heat[i];r[4]=q[i];r[5]=real[i];r[6]=-cell[i].p*real[i]+heat[i];
  *potential_power+=gx[2*i]*cell[i].u[0]+gx[2*i+1]*cell[i].u[1]+gw[i]*q[i];
 }
 ok=isfinite(*potential_power);
 done:if(!ok)fprintf(stderr,"[MATERIAL_GRAVITY_FAIL] phase=%s prior_iterations=%d\n",phase,*iterations);
 free(work);free(zero);free(cache);return ok;
}
/* Reuse is scoped to one unchanged RHS; legacy callers retain no-cache behavior. */
static inline int gfs_material_gravity_active_rhs2d_status(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const unsigned char*reflected_y,
 const GfsFvPrimitive*cell,const double*centroid_relative,const double grad_phi[2],
 double gamma,const double*alpha3,double*out7,double*potential_power,int*iterations,
 int active_count,const GfsMaterialActiveFace2D*active,double beta,int *release_index){
 return gfs_material_gravity_active_rhs2d_reuse_status(t,nd,edge,reflected_y,cell,centroid_relative,grad_phi,gamma,alpha3,out7,potential_power,iterations,active_count,active,beta,release_index,NULL);
}
static inline int gfs_material_gravity_active_rhs2d(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const unsigned char*reflected_y,
 const GfsFvPrimitive*cell,const double*centroid_relative,const double grad_phi[2],
 double gamma,const double*alpha3,double*out7,double*potential_power,int*iterations,
 int active_count,const GfsMaterialActiveFace2D*active,double beta){
 return gfs_material_gravity_active_rhs2d_status(t,nd,edge,reflected_y,cell,centroid_relative,grad_phi,gamma,alpha3,out7,potential_power,iterations,active_count,active,beta,NULL);
}
static inline int gfs_material_gravity_rhs2d(const GfsTraceCsr2D*t,int nd,
 const GfsTraceDirected2D*edge,const unsigned char*reflected_y,
 const GfsFvPrimitive*cell,const double*centroid_relative,const double grad_phi[2],
 double gamma,const double*alpha3,double*out7,double*potential_power,int*iterations){
 return gfs_material_gravity_active_rhs2d(t,nd,edge,reflected_y,cell,centroid_relative,grad_phi,gamma,alpha3,out7,potential_power,iterations,0,NULL,.9);
}
#endif
