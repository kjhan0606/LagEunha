#ifndef GFS_MATERIAL_CONSTRAINT_OPS2D_H
#define GFS_MATERIAL_CONSTRAINT_OPS2D_H
#include <stdio.h>
#include "gfs_material_rhs2d.h"
/* Frozen-state linear operator actions for the compatible constraint audit.
 * No active-set selection, dense global matrix, clipping or time evolution.
 * Acoustic impedances are evaluated from the ORIGINAL primitive traces. */
static inline int gfs_material_constraint_cache_maps2d(const GfsTraceCsr2D*t,
 const GfsFvPrimitive*c,double gamma,const double*alpha3,GfsMaterialFaceResponse*cache){
 if(!t||!c||!cache)return 0;
 if(alpha3)for(int i=0;i<3*t->cells;i++)if(!isfinite(alpha3[i])||alpha3[i]<0||alpha3[i]>1)return 0;
 for(int k=0;k<t->faces;k++){
  const GfsTraceFace2D*f=t->face+k;GfsFvPrimitive q[2]={{0}};
  int starts[2]={f->sl,f->sr},counts[2]={f->nl,f->nr};
  for(int s=0;s<2;s++)for(int j=starts[s];j<starts[s]+counts[s];j++){
   const GfsTraceEntry*e=t->entry+j;int owner=s?f->j:f->i;
   GfsMaterialTraceWeights w=gfs_material_entry_weights(e,owner,alpha3?alpha3+3*owner:NULL);
   q[s].rho+=w.rho*c[e->cell].rho;q[s].p+=w.pressure*c[e->cell].p;
   for(int a=0;a<2;a++)q[s].u[a]+=w.velocity*gfs_trace_velocity_sign(e,a)*c[e->cell].u[a];
  }
  double normal[3]={f->nx,f->ny,0};
  int response=gfs_material_face_response(2,q[0],q[1],normal,gamma,cache+k);
  if(response!=1){
   fprintf(stderr,"[MATERIAL_FACE_FAIL] face=%d i=%d j=%d status=%d area=%.17g rhoL=%.17g rhoR=%.17g pL=%.17g pR=%.17g vnL=%.17g vnR=%.17g\n",k,f->i,f->j,response,f->area,q[0].rho,q[1].rho,q[0].p,q[1].p,q[0].u[0]*f->nx+q[0].u[1]*f->ny,q[1].u[0]*f->nx+q[1].u[1]*f->ny);return 0;
  }
 }return 1;
}
/* Dg=Jp^T diag(A/(ZL+ZR)) Jp g and C^T g, using one face traversal. */
static inline int gfs_material_constraint_actions_maps2d(const GfsTraceCsr2D*t,
 const GfsMaterialFaceResponse*cache,const double*alpha3,const double*g,double*Dg,double*Ctg){
 if(!t||!cache||!g||!Dg||!Ctg)return 0;
 memset(Dg,0,(size_t)t->cells*sizeof(double));memset(Ctg,0,(size_t)t->cells*2*sizeof(double));
 for(int k=0;k<t->faces;k++){
  const GfsTraceFace2D*f=t->face+k;const GfsMaterialFaceResponse*r=cache+k;
  int starts[2]={f->sl,f->sr},counts[2]={f->nl,f->nr};double jump=0;
  for(int s=0;s<2;s++)for(int j=starts[s];j<starts[s]+counts[s];j++){
   const GfsTraceEntry*e=t->entry+j;int owner=s?f->j:f->i;
   GfsMaterialTraceWeights w=gfs_material_entry_weights(e,owner,alpha3?alpha3+3*owner:NULL);
   jump+=(s?-1:1)*w.pressure*g[e->cell];
  }
  double z=r->zleft+r->zright,dp=g[f->i]-g[f->j];if(!(z>0)||!isfinite(jump)||!isfinite(dp))return 0;
  double rate=f->area*jump/z;
  for(int s=0;s<2;s++)for(int j=starts[s];j<starts[s]+counts[s];j++){
   const GfsTraceEntry*e=t->entry+j;int owner=s?f->j:f->i;
   GfsMaterialTraceWeights w=gfs_material_entry_weights(e,owner,alpha3?alpha3+3*owner:NULL);
   double force=w.reversible*f->area*(s?r->alpha_right:r->alpha_left)*dp;
   Dg[e->cell]+=(s?-1:1)*w.pressure*rate;Ctg[2*e->cell]+=force*f->nx;Ctg[2*e->cell+1]+=force*f->ny*gfs_trace_velocity_sign(e,1);
  }
 }return 1;
}
/* Geometry transpose B^T g. Directed faces must be the same validated
 * owner-relative, periodic-image-preserving geometry as the weight solve. */
static inline int gfs_material_geometry_transpose_boundaries2d(int n,int nf,
 const GfsTraceDirected2D*edges,const double*g,double*out,const unsigned char*reflected_y){
 if(n<1||nf<1||!edges||!g||!out)return 0;memset(out,0,(size_t)n*2*sizeof(double));
 for(int k=0;k<nf;k++){
  const GfsTraceDirected2D*e=edges+k;
  if(e->i<0||e->i>=n||e->j<0||e->j>=n)return 0;
  int wall=reflected_y?reflected_y[k]:0;
  if(wall>1)return 0;
  if(wall){
   /* Own reflected generator follows the real generator; its fixed-wall
    * bisector has identically zero volume derivative, as in the weight solve.
    * Do not pull back this image as an untranslated periodic owner. */
   if(!gfs_material_stationary_wall2d(e))return 0;
   continue;
  }
  double delta[2]={e->dx,e->dy},centroid[2]={e->cx,e->cy};
  GfsMaterialPowerFaceMotion face;
  if(!gfs_material_power_face_motion(2,delta,centroid,e->area,&face))return 0;
  for(int a=0;a<2;a++){
   out[2*e->i+a]+=g[e->i]*face.left[a];
   out[2*e->j+a]+=g[e->i]*face.right[a];
  }
 }return 1;
}
static inline int gfs_material_geometry_transpose2d(int n,int nf,
 const GfsTraceDirected2D*edges,const double*g,double*out){
 return gfs_material_geometry_transpose_boundaries2d(n,nf,edges,g,out,NULL);
}
/* Assemble the effective drive directly. Subtracting two separately formed
 * volume RHS vectors loses compatibility near hydrostatic cancellation. */
static inline int gfs_material_cached_drive_rhs2d(const GfsTraceCsr2D*t,
 const GfsMaterialFaceResponse*cache,const GfsFvPrimitive*cell,
 const double*alpha3,const double*drive,double*out){
 memset(out,0,(size_t)t->cells*4*sizeof(double));
 for(int k=0;k<t->faces;k++){
  const GfsTraceFace2D*f=t->face+k;GfsMaterialFaceResponse r=cache[k];
  int starts[2]={f->sl,f->sr},counts[2]={f->nl,f->nr};double jp=0,jv=0,mean=0;
  for(int s=0;s<2;s++)for(int j=starts[s];j<starts[s]+counts[s];j++){
   const GfsTraceEntry*e=t->entry+j;int owner=s?f->j:f->i;
   GfsMaterialTraceWeights w=gfs_material_entry_weights(e,owner,alpha3?alpha3+3*owner:NULL);
   double vn=f->nx*cell[e->cell].u[0]+f->ny*gfs_trace_velocity_sign(e,1)*cell[e->cell].u[1];
   jp+=(s?-1:1)*w.pressure*drive[e->cell];jv+=(s?1:-1)*w.velocity*vn;
   mean+=(s?r.alpha_right:r.alpha_left)*w.reversible*vn;
  }
  double z=r.zleft+r.zright;if(!(z>0))return 0;
  r.pressure_speed=jp/z;r.mean_speed=mean;
  double dl=r.pressure_speed+r.alpha_right*jv,dr=r.pressure_speed-r.alpha_left*jv;
  r.heat_left=r.zleft*dl*dl;r.heat_right=r.zright*dr*dr;
  double normal[3]={f->nx,f->ny,0};GfsMaterialFaceExchange x;
  if(!gfs_material_face_exchange_drive(2,normal,f->area,drive[f->i]-drive[f->j],&r,&x))return 0;
  for(int s=0;s<2;s++)for(int j=starts[s];j<starts[s]+counts[s];j++){
   const GfsTraceEntry*e=t->entry+j;int owner=s?f->j:f->i;
   GfsMaterialTraceWeights w=gfs_material_entry_weights(e,owner,alpha3?alpha3+3*owner:NULL);
   out[4*e->cell]+=gfs_material_scatter_force(w,&x,s,0);
   out[4*e->cell+1]+=gfs_trace_velocity_sign(e,1)*gfs_material_scatter_force(w,&x,s,1);
   out[4*e->cell+2]+=w.pressure*x.side_volume[s];
  }
  out[4*f->i+2]+=x.owner_volume[0];out[4*f->j+2]+=x.owner_volume[1];
  out[4*f->i+3]+=x.owner_heat[0];out[4*f->j+3]+=x.owner_heat[1];
 }
 return 1;
}
/* Positive local heat after pressure projection: P_eff=P-G lambda_P.
 * Original trace velocities and frozen impedances are retained. */
static inline int gfs_material_projected_heat_maps2d(const GfsTraceCsr2D*t,
 const GfsMaterialFaceResponse*cache,const GfsFvPrimitive*c,
 const double*alpha3,const double*effective_pressure,double*heat){
 if(!t||!cache||!c||!effective_pressure||!heat)return 0;memset(heat,0,(size_t)t->cells*sizeof(double));
 for(int k=0;k<t->faces;k++){
  const GfsTraceFace2D*f=t->face+k;const GfsMaterialFaceResponse*r=cache+k;
  int starts[2]={f->sl,f->sr},counts[2]={f->nl,f->nr};double jp=0,jv=0;
  for(int s=0;s<2;s++)for(int j=starts[s];j<starts[s]+counts[s];j++){
   const GfsTraceEntry*e=t->entry+j;int owner=s?f->j:f->i;
   GfsMaterialTraceWeights w=gfs_material_entry_weights(e,owner,alpha3?alpha3+3*owner:NULL);
   jp+=(s?-1:1)*w.pressure*effective_pressure[e->cell];jv+=(s?1:-1)*w.velocity*(f->nx*c[e->cell].u[0]+f->ny*gfs_trace_velocity_sign(e,1)*c[e->cell].u[1]);
  }
  double wp=jp/(r->zleft+r->zright),dl=wp+r->alpha_right*jv,dr=wp-r->alpha_left*jv;
  double hl=f->area*r->zleft*dl*dl,hr=f->area*r->zright*dr*dr;
  if(!isfinite(hl)||!isfinite(hr)||hl<0||hr<0)return 0;
  heat[f->i]+=hl;heat[f->j]+=hr;
 }return 1;
}
static inline int gfs_material_constraint_cache2d(const GfsTraceCsr2D*t,const GfsFvPrimitive*c,double gamma,GfsMaterialFaceResponse*cache){return gfs_material_constraint_cache_maps2d(t,c,gamma,NULL,cache);}
static inline int gfs_material_constraint_actions2d(const GfsTraceCsr2D*t,const GfsMaterialFaceResponse*cache,const double*g,double*Dg,double*Ctg){return gfs_material_constraint_actions_maps2d(t,cache,NULL,g,Dg,Ctg);}
static inline int gfs_material_projected_heat2d(const GfsTraceCsr2D*t,const GfsMaterialFaceResponse*cache,const GfsFvPrimitive*c,const double*effective_pressure,double*heat){return gfs_material_projected_heat_maps2d(t,cache,c,NULL,effective_pressure,heat);}
#endif
