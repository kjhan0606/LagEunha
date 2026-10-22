#ifndef MATERIAL_P3_FINITE_BATCH_H
#define MATERIAL_P3_FINITE_BATCH_H
#include <stddef.h>
#include "gfs_material_map12.h"
#ifdef __CUDACC__
#define P3_FINITE_HD __host__ __device__
#else
#define P3_FINITE_HD
#endif
typedef struct {GfsMaterialMap12 map;double origin[2],mass[12][12];} P3FiniteMap;
typedef struct {int triangle;double N[12],Nx[12],Ny[12],rho0,weight;} P3FinitePoint;
typedef struct {double gamma,pressure0,eps;int source_mode;} P3FinitePhysics;
typedef struct {
    GfsMaterialPoint12 point;
    double force[2][12],Udot,gravity_power,gravity_torque,gravity_force[2],thermal,pressure;
    double energy_identity,angular_identity,energy_scale,angular_scale;
    int valid;
} P3FiniteOutput;
typedef struct {GfsMaterialMoments12 moment;int valid;} P3FiniteMoment;
/* Pressure is the ADIABATIC GEOMETRY CONTROL P0*J^-gamma, not a production
 * update or a reset of transported thermal state. No physical viscosity.
 */
P3_FINITE_HD static inline P3FiniteOutput p3_finite_evaluate(
    const P3FiniteMap *m,const P3FinitePoint *q,const P3FinitePhysics *physics){
    P3FiniteOutput out={};if(!gfs_material12_point(&m->map,q->N,q->Nx,q->Ny,&out.point))return out;
    out.pressure=physics->pressure0*pow(out.point.J,-physics->gamma);
    out.thermal=q->weight*physics->pressure0/(physics->gamma-1)*pow(out.point.J,1-physics->gamma);
    if(!(out.pressure>0) || !isfinite(out.pressure) || !(out.thermal>=0))return out;
    double x=out.point.position[0]+m->origin[0],y=out.point.position[1]+m->origin[1],gx=0,gy=0;
    if(physics->source_mode==1){gx=-x;gy=-y;}
    if(physics->source_mode==2){double inv=pow(x*x+y*y+physics->eps*physics->eps,-1.5);gx=-x*inv;gy=-y*inv;}
    double sigma[3]={-out.pressure,0,-out.pressure};
    if(!gfs_material12_stress_at_point(m->map.n,&out.point,q->N,q->Nx,q->Ny,q->rho0,q->weight,sigma,gx,gy,out.force,&out.Udot,&out.gravity_power))return out;
    out.gravity_torque=q->weight*q->rho0*(x*gy-y*gx);out.energy_identity=out.Udot-out.gravity_power;out.angular_identity=-out.gravity_torque;
    out.gravity_force[0]=q->weight*q->rho0*gx;out.gravity_force[1]=q->weight*q->rho0*gy;
    out.energy_scale=fabs(out.Udot)+fabs(out.gravity_power);out.angular_scale=fabs(out.gravity_torque);
    for(int a=0;a<m->map.n;a++){
        double constant=a<10?1:0,X=m->map.x[0][a]+m->origin[0]*constant,Y=m->map.x[1][a]+m->origin[1]*constant;
        double ex=m->map.v[0][a]*out.force[0][a],ey=m->map.v[1][a]*out.force[1][a],ax=X*out.force[1][a],ay=-Y*out.force[0][a];
        out.energy_identity+=ex+ey;out.angular_identity+=ax+ay;out.energy_scale+=fabs(ex)+fabs(ey);out.angular_scale+=fabs(ax)+fabs(ay);
    }
    out.valid=1;return out;
}
#ifdef __cplusplus
extern "C" {
#endif
int p3_finite_cuda_create(const P3FinitePoint *,size_t,size_t,void **,double *);
int p3_finite_cuda_eval(void *,const P3FiniteMap *,const P3FinitePhysics *,P3FiniteOutput *,P3FiniteMoment *,double *,double *);
int p3_finite_cuda_destroy(void *);
const char *p3_finite_cuda_error(int);
#ifdef __cplusplus
}
#endif
#undef P3_FINITE_HD
#endif
