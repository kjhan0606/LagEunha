#ifndef GFS_MATERIAL_POINT_GRAVITY_H
#define GFS_MATERIAL_POINT_GRAVITY_H
#include <math.h>
/* Fixed-mass particle quadrature for a prescribed linear potential.
 * xdot=v, Phi=grad_phi.x: F=-m grad_phi, Edot=v.F, Udot=0.
 * This is NOT the exact integral over an off-centre uniform-density cell.
 * Caller must use the same quadrature in its endpoint energy ledger.
 * No mass exchange, weight force, thermal repair or timestep rule. */
typedef struct {double force[3],potential,potential_rate,energy_rate;} GfsMaterialPointGravity;
#ifdef __CUDACC__
#define GFS_POINT_HD static __host__ __device__ inline
#else
#define GFS_POINT_HD static inline
#endif
GFS_POINT_HD int gfs_material_point_gravity(int dim,double mass,
 const double x[3],const double v[3],const double grad_phi[3],GfsMaterialPointGravity*out){
 if(dim<1||dim>3||!(mass>0)||!isfinite(mass)||!x||!v||!grad_phi||!out)return 0;
 GfsMaterialPointGravity s={{0},0,0,0};
 for(int a=0;a<dim;a++){
  if(!isfinite(x[a])||!isfinite(v[a])||!isfinite(grad_phi[a]))return 0;
  s.force[a]=-mass*grad_phi[a];
  s.potential+=mass*grad_phi[a]*x[a];
  s.potential_rate+=mass*grad_phi[a]*v[a];
  if(!isfinite(s.force[a]))return 0;
 }
 s.energy_rate=-s.potential_rate;
 if(!isfinite(s.potential)||!isfinite(s.potential_rate))return 0;
 *out=s;return 1;
}
#undef GFS_POINT_HD
#endif
