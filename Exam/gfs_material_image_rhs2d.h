#ifndef GFS_MATERIAL_IMAGE_RHS2D_H
#define GFS_MATERIAL_IMAGE_RHS2D_H
#include "gfs_material_image_trace2d.h"
/* Existing CPU or CUDA moving-wall RHS signature. The same kernel is used
 * before and after recovering the wall speed; no separate flux model. */
typedef int (*GfsMaterialImageRhsEval2D)(const GfsTraceCsr2D*,const GfsFvPrimitive*,
 double,const double*,const double*,double*,double*);
typedef struct {
 int weight_iterations;
 double weight_residual,volume_error,wall_power,exterior_ledger[4],energy_defect;
} GfsMaterialImageRhsStats2D;
/* This is an immutable-state RHS adapter, not a time integrator. Only cells
 * [0,g->n) may be advanced; remaining slots are prescribed sample/face data.
 * Geometry, sample and mesh velocities must describe this same state.
 * Caller must discard ALL outputs on failure. Does not certify a moving
 * image boundary as a stationary obstacle, nor supply periodic aliasing.
 */
static inline int gfs_material_image_rhs2d_with_self_speeds(const GfsMaterialExterior2D*g,
 GfsMaterialImageTrace2D*t,const GfsFvPrimitive*mesh_velocity,const int*image_owner,
 const double*known_weight_rate,double gamma,GfsMaterialImageRhsEval2D evaluate,
 double*q,double*realized,double*wall_normal,double*rhs,GfsMaterialImageRhsStats2D*stats,const double*self_speed){
 if(!g||!t||!mesh_velocity||!image_owner||!known_weight_rate||!evaluate||!q||!realized||!wall_normal||!rhs||!stats)return 0;
 GfsMaterialImageRhsStats2D s={0};int n=g->n,ok=0;
 if(t->fluid.trace.faces!=g->pair_count||n<1)return 0;
 double*target=calloc(n,sizeof(double));if(!target)return 0;
 memset(wall_normal,0,(size_t)g->pair_count*sizeof(double));
 if(!gfs_material_face_normal_prepare2d_wall_velocity(&t->fluid.trace,t->fluid.state,wall_normal))goto done;
 if(evaluate(&t->fluid.trace,t->fluid.state,gamma,NULL,wall_normal,rhs,&s.wall_power)!=1)goto done;
 for(int i=0;i<n;i++)target[i]=rhs[4*i+2];
 if(!gfs_material_image_weight_rate2d_with_self_speeds(g,mesh_velocity,image_owner,known_weight_rate,target,q,realized,wall_normal,&s.weight_iterations,&s.weight_residual,self_speed))goto done;
 if(!gfs_material_face_normal_refresh_walls2d(&t->fluid.trace,t->fluid.state,wall_normal))goto done;
 if(evaluate(&t->fluid.trace,t->fluid.state,gamma,NULL,wall_normal,rhs,&s.wall_power)!=1)goto done;
 if(!gfs_material_exterior_trace2d_ledger(n,&t->fluid,rhs,s.exterior_ledger))goto done;
 s.energy_defect=-s.wall_power-s.exterior_ledger[3];
 for(int i=0;i<n;i++){
  GfsFvPrimitive c=t->fluid.state[i];
  double e=rhs[4*i+3]-c.p*rhs[4*i+2]+c.u[0]*rhs[4*i]+c.u[1]*rhs[4*i+1];
  if(!isfinite(e)||!isfinite(realized[i])||!isfinite(rhs[4*i+2]))goto done;
  s.energy_defect+=e;s.volume_error=fmax(s.volume_error,fabs(realized[i]-rhs[4*i+2]));
 }
 if(!isfinite(s.energy_defect))goto done;
 *stats=s;ok=1;
 done:free(target);return ok;
}
static inline int gfs_material_image_rhs2d(const GfsMaterialExterior2D*g,
 GfsMaterialImageTrace2D*t,const GfsFvPrimitive*mesh_velocity,const int*image_owner,
 const double*known_weight_rate,double gamma,GfsMaterialImageRhsEval2D evaluate,
 double*q,double*realized,double*wall_normal,double*rhs,GfsMaterialImageRhsStats2D*stats){
 return gfs_material_image_rhs2d_with_self_speeds(g,t,mesh_velocity,image_owner,
  known_weight_rate,gamma,evaluate,q,realized,wall_normal,rhs,stats,NULL);
}
#endif
