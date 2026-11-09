#ifndef GFS_MATERIAL_CYLINDER_BOUNDARY2D_H
#define GFS_MATERIAL_CYLINDER_BOUNDARY2D_H
#include "Cylinder/cylinder_padding_owner.h"
#include "gfs_material_radial_image.h"

/* Kinematics of the ORIGINAL Cylinder padding, not a new boundary law.
 * The caller supplies the physical ghost velocity made by that padding.
 * In particular, an open ghost copies/prescribes fluid velocity although its
 * reflected position moves with (-vx,vy). Every generated ghost copies W.
 * Only radial images are walls; only pure y translations are owned aliases.
 * Keeping these maps separate is essential: copied W does not imply a wall.
 * This helper alone does not qualify a material RHS at a stationary opening.
 */
typedef struct {
    double position[2], low[2], mesh_velocity[2], fluid_velocity[2];
    double weight, weight_rate;
    int weight_owner, wall_owner, periodic_owner;
} GfsMaterialCylinderImage2D;
#ifdef __CUDACC__
#define GFS_CYL_IMAGE_HD static __host__ __device__ inline
#else
#define GFS_CYL_IMAGE_HD static inline
#endif
GFS_CYL_IMAGE_HD int gfs_material_cylinder_image2d(
 const CylinderPaddingOwner *map,int owned_count,long long source_id,
 const double source_position[2],const double source_velocity[2],
 double source_weight,double source_weight_rate,const double ghost_velocity[2],
 const double center[2],double radius,double lx,double ly,
 GfsMaterialCylinderImage2D *out){
 if(!map||!source_position||!source_velocity||!ghost_velocity||!center||!out||
    owned_count<1||map->source_index<0||map->source_index>=owned_count||
    map->source_id!=source_id||map->kind<CYL_IMAGE_RADIAL||map->kind>CYL_IMAGE_OUTFLOW||
    !(lx>0)||!(ly>0)||!isfinite(lx)||!isfinite(ly)||
    !isfinite(source_weight)||!isfinite(source_weight_rate))return 0;
 GfsMaterialCylinderImage2D a={{0},{0},{0},{0},source_weight,source_weight_rate,
                              map->source_index,-1,-1};
 for(int k=0;k<2;k++){
  if(!isfinite(source_position[k])||!isfinite(source_velocity[k])||!isfinite(ghost_velocity[k]))return 0;
  a.position[k]=source_position[k];a.mesh_velocity[k]=source_velocity[k];
  a.fluid_velocity[k]=ghost_velocity[k];
 }
 if(map->kind==CYL_IMAGE_RADIAL){
  GfsMaterialRadialImage radial;
  if(map->y_shift!=0||!gfs_material_radial_image(2,center,radius,source_position,&radial)||
     !gfs_material_radial_split_point(2,center,radius,source_position,a.position,a.low)||
     !gfs_material_radial_vector(2,&radial,source_velocity,1,a.mesh_velocity))return 0;
  a.wall_owner=map->source_index;
 }else{
  if(map->y_shift < -1||map->y_shift > 1)return 0;
  GfsPowerDD y=gfs_pdd_sum(source_position[1],map->y_shift*ly);
  a.position[1]=y.hi;a.low[1]=y.lo;
  if(map->kind==CYL_IMAGE_PERIODIC){
   if(map->y_shift==0)return 0;
   a.periodic_owner=map->source_index;
  }else{
   GfsPowerDD x=gfs_pdd_sum(map->kind==CYL_IMAGE_OUTFLOW?2*lx:0,-source_position[0]);
   a.position[0]=x.hi;a.low[0]=x.lo;a.mesh_velocity[0]=-source_velocity[0];
  }
 }
 for(int k=0;k<2;k++)if(!isfinite(a.position[k])||!isfinite(a.low[k])||!isfinite(a.mesh_velocity[k]))return 0;
 *out=a;return 1;
}
#undef GFS_CYL_IMAGE_HD
#endif
