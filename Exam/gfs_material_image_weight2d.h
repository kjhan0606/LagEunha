#ifndef GFS_MATERIAL_IMAGE_WEIGHT2D_H
#define GFS_MATERIAL_IMAGE_WEIGHT2D_H
#include "gfs_material_exterior2d.h"
/* Curved image boundary: its flux volume equals its ACTUAL geometric rate,
 * including Wdot of the image owner. That identical term cancels from the
 * volume constraint. Solve only interior/prescribed-exterior edges, then
 * recover all image-face speeds. No stationary-wall residual projection.
 *
 * image_owner[j]>=0 means external site j has Wdot_j=Wdot_owner; its supplied
 * mesh velocity must be the derivative of the coordinate image, not fluid
 * reflection. All other external sites have prescribed known Wdot. Periodic
 * aliases must be resolved before this adapter. The existing anchored PCG
 * is used, so every remaining graph component needs a prescribed exterior.
 * target_without_images contains only non-image volume fluxes, per owned cell.
 * Outputs are invalid on failure; caller commits no partial physical state.
 */
static inline int gfs_material_image_weight_rate2d_with_self_speeds(const GfsMaterialExterior2D*g,
 const GfsFvPrimitive*mesh_velocity,const int*image_owner,const double*known,
 const double*target_without_images,double*q,double*realized,
 double*image_normal_speed,int*iterations,double*relative,const double*self_speed){
 if(!g||g->n<1||g->nt<=g->n||!mesh_velocity||!image_owner||!known||
    !target_without_images||!q||!realized||!image_normal_speed||!iterations||!relative)return 0;
 int n=g->n,ok=0,count=0;
 GfsMaterialExteriorEdge*edge=malloc((size_t)g->pair_count*sizeof(*edge));
 double*work=calloc((size_t)2*n,sizeof(double));
 if(!edge||!work)goto done;
 {
 double*rhs=work,*base=work+n;
 for(int i=0;i<n;i++){if(!isfinite(target_without_images[i]))goto done;rhs[i]=target_without_images[i];}
 for(int j=n;j<g->nt;j++)if(image_owner[j]<-1||image_owner[j]>=n)goto done;
 for(int k=0;k<g->pair_count;k++){
  const GfsTraceDirected2D*e=g->directed+g->pair_directed[k];int i=e->i,j=e->j;
  image_normal_speed[k]=0;
  if(j>=n&&image_owner[j]>=0)continue;
#ifdef GFS_MATERIAL_FIXED_SOLID
  if(i==j)continue;
#endif
  double d[2]={e->dx,e->dy},c[2]={e->cx,e->cy},b,coefficient;
  if(!gfs_material_power_face_rate(2,d,c,mesh_velocity[i].u,mesh_velocity[j].u,e->area,&b,&coefficient))goto done;
  edge[count++]=(GfsMaterialExteriorEdge){i,j,coefficient};base[i]+=b;
  if(j<n)base[j]-=b;
  else {if(!isfinite(known[j]))goto done;rhs[i]+=coefficient*known[j];}
 }
 for(int i=0;i<n;i++)rhs[i]-=base[i];
 if(!gfs_material_exterior_solve_selected(n,g->nt,count,edge,rhs,q,iterations,relative))goto done;
 memset(realized,0,(size_t)n*sizeof(double));
 for(int k=0;k<g->pair_count;k++){
  const GfsTraceDirected2D*e=g->directed+g->pair_directed[k];int i=e->i,j=e->j;
#ifdef GFS_MATERIAL_FIXED_SOLID
  if(i==j)continue;
#endif
  double d[2]={e->dx,e->dy},c[2]={e->cx,e->cy},rate,neighbor_q;
  GfsMaterialPowerFaceMotion f;
  if(j<n)neighbor_q=q[j];else if(image_owner[j]>=0)neighbor_q=q[image_owner[j]];else neighbor_q=known[j];
  if(self_speed&&j>=n&&image_owner[j]==i){
   /* The equal-weight self-image has no Wdot difference. Its caller supplies
    * the factored derivative of the same image coordinates used by geometry. */
   rate=e->area*self_speed[k];if(!isfinite(rate))goto done;
  }else if(!gfs_material_power_face_motion(2,d,c,e->area,&f)||
     !gfs_material_power_face_motion_rate(2,&f,mesh_velocity[i].u,mesh_velocity[j].u,q[i],neighbor_q,&rate))goto done;
  realized[i]+=rate;if(j<n)realized[j]-=rate;
  else if(image_owner[j]>=0)image_normal_speed[k]=rate/e->area;
 }
 for(int i=0;i<n;i++)if(!isfinite(realized[i]))goto done;
 ok=1;
 }
 done:free(edge);free(work);return ok;
}
static inline int gfs_material_image_weight_rate2d(const GfsMaterialExterior2D*g,
 const GfsFvPrimitive*mesh_velocity,const int*image_owner,const double*known,
 const double*target_without_images,double*q,double*realized,
 double*image_normal_speed,int*iterations,double*relative){
 return gfs_material_image_weight_rate2d_with_self_speeds(g,mesh_velocity,image_owner,
  known,target_without_images,q,realized,image_normal_speed,iterations,relative,NULL);
}
#endif
