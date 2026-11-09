#ifndef GFS_MATERIAL_DIRECT_FLUX_H
#define GFS_MATERIAL_DIRECT_FLUX_H
#include "gfs_material_face_response.h"
/* Candidate shared material-face exchange. Apply to left owner and its
 * negative to right owner. Energy is TOTAL energy, not internal heating.
 * No mass flux. Does not by itself specify weight constraints or a limiter. */
typedef struct { double force[3],volume,energy; } GfsMaterialDirectFlux;
#ifdef __CUDACC__
#define GFS_MATERIAL_DIRECT_HD static __host__ __device__ inline
#else
#define GFS_MATERIAL_DIRECT_HD static inline
#endif
GFS_MATERIAL_DIRECT_HD int gfs_material_direct_flux(int dim,const double normal[3],
 double area,const GfsMaterialFaceResponse *r,GfsMaterialDirectFlux *out){
 if(dim<1||dim>3||!normal||!r||!out||!(area>0)||!isfinite(area)||
    !isfinite(r->pstar)||!isfinite(r->ustar)||r->pstar<0)return 0;
 GfsMaterialDirectFlux q={{0},area*r->ustar,-area*r->pstar*r->ustar};
 double norm=0;
 for(int a=0;a<dim;a++){norm+=normal[a]*normal[a];q.force[a]=-area*r->pstar*normal[a];if(!isfinite(q.force[a]))return 0;}
 if(!isfinite(norm)||fabs(norm-1)>128*DBL_EPSILON||!isfinite(q.volume)||!isfinite(q.energy))return 0;
 *out=q;return 1;
}
/* Prescribed normal wall motion: the SAME wall speed enters reflection,
 * volume change and total-energy work. The caller obtains that speed from
 * its geometry and records the external wall force/work. This does not make
 * an approximate moving boundary equivalent to a stationary physical wall.
 * A symmetric vacuum fan exerts zero pressure while the boundary still moves.
 */
GFS_MATERIAL_DIRECT_HD int gfs_material_direct_moving_wall_flux(int dim,
 GfsFvPrimitive interior,const double normal[3],double area,double gamma,
 double wall_speed,GfsMaterialDirectFlux *out){
 if(dim<1||dim>3||!normal||!isfinite(wall_speed))return 0;
 GfsFvPrimitive ghost=interior;double wall[3]={0};
 for(int a=0;a<dim;a++)wall[a]=wall_speed*normal[a];
 if(!gfs_material_slip_velocity(dim,normal,wall,interior.u,ghost.u))return 0;
 GfsMaterialFaceResponse r;int status=gfs_material_face_response(dim,interior,ghost,normal,gamma,&r);
 if(!status)return 0;
 if(status==2){r=(GfsMaterialFaceResponse){0};}
 else{
  double scale=sqrt(gamma*interior.p/interior.rho);
  for(int a=0;a<dim;a++)scale+=fabs(interior.u[a]);
  scale+=fabs(wall_speed);
  if(!isfinite(r.ustar)||fabs(r.ustar-wall_speed)>256*DBL_EPSILON*scale)return 0;
 }
 r.ustar=wall_speed; /* boundary datum, not a fluid-state repair */
 return gfs_material_direct_flux(dim,normal,area,&r,out);
}
/* Existing stationary-wall entry point, identical physical boundary condition. */
GFS_MATERIAL_DIRECT_HD int gfs_material_direct_wall_flux(int dim,
 GfsFvPrimitive interior,const double normal[3],double area,double gamma,
 GfsMaterialDirectFlux *out){
 return gfs_material_direct_moving_wall_flux(dim,interior,normal,area,gamma,0,out);
}
#undef GFS_MATERIAL_DIRECT_HD
#endif
