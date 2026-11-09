#ifndef GFS_NOH_EXTERIOR_H
#define GFS_NOH_EXTERIOR_H
#include "gfs_material_face_response.h"
/* Prescribed upstream Noh reservoir. The density and radial trajectories are
 * the cold Noh solution in dimension d. p0>0 is the stated finite-pressure
 * boundary approximation, held fixed as in the existing test convention;
 * this is not an exact finite-pressure solution of the exterior Euler flow.
 * Never use this provider inside the expanding cold-solution shock.
 */
typedef struct {int dim;double center[3],rho0,speed,p0,gamma,time;} GfsNohExterior;
static inline int gfs_noh_exterior_state(const GfsNohExterior*b,const double*x,GfsFvPrimitive*out){
 if(!b||!x||!out||b->dim<1||b->dim>3||!(b->rho0>0)||!(b->speed>0)||!(b->p0>0)||!(b->gamma>1)||b->time<0||!isfinite(b->time))return 0;
 double r2=0;for(int a=0;a<b->dim;a++)r2+=(x[a]-b->center[a])*(x[a]-b->center[a]);
 double r=sqrt(r2),shock=.5*(b->gamma-1)*b->speed*b->time;
 if(!(r>shock)||!isfinite(r))return 0;
 GfsFvPrimitive q={b->rho0*pow(1+b->speed*b->time/r,b->dim-1),b->p0,{0}};
 for(int a=0;a<b->dim;a++)q.u[a]=-b->speed*(x[a]-b->center[a])/r;
 if(!isfinite(q.rho)||!isfinite(q.p))return 0;*out=q;return 1;
}
static inline int gfs_noh_exterior_position(const GfsNohExterior*b,const double*x0,double*x){
 if(!b||!x0||!x||b->dim<1||b->dim>3||!(b->speed>0)||b->time<0||!isfinite(b->time))return 0;
 double r2=0;for(int a=0;a<b->dim;a++)r2+=(x0[a]-b->center[a])*(x0[a]-b->center[a]);
 double r=sqrt(r2);if(!(r>b->speed*b->time)||!isfinite(r))return 0;
 for(int a=0;a<3;a++)x[a]=a<b->dim?b->center[a]+(1-b->speed*b->time/r)*(x0[a]-b->center[a]):0;
 return 1;
}
#endif
