#ifndef MATERIAL_GEOMETRIC_PRESSURE_BATCH_H
#define MATERIAL_GEOMETRIC_PRESSURE_BATCH_H
#include "material_shared_trace_batch.h"
#ifdef __CUDACC__
#define GP_HD __host__ __device__
#else
#define GP_HD
#endif
// U is the physical point energy supplied by the caller for this state.
// Saved CONTROL audits keep OLD U immutable. Smooth evolution supplies
// stage U from immutable generating entropy data; this kernel never resets U.
// deltaJ and deltaJdot are same-stage material-cell volume corrections.
struct GPInput {double U,deltaJ,mean_deltaP,deltaJdot;};
struct GPOutput {STOutput native;double nu,rho,deltaP,nudot;};
GP_HD static inline GPOutput gp_evaluate(const STTriangle*t,const STCoefficient*c,const STPoint&q,const STPhysics&physics,const GPInput&input){
    GPOutput out={};auto&o=out.native;const auto&T=t[q.triangle];STJet N[ST_MAX];
    if(!st_shape(T,q,N)||!st_point(T,c,N,&o.point))return out;
    out.nu=q.weight*(o.point.J+input.deltaJ);if(!(out.nu>0)||!(input.U>0))return out;
    out.rho=q.weight*q.rho/out.nu;
    // Compensated subtraction of the fixed reference pressure. This is
    // arithmetic only: U is not projected, reset, floored, or corrected.
    double target=physics.P0*out.nu;
    out.deltaP=(fma(physics.gamma-1,input.U,-target)+fma(-physics.P0,out.nu,target))/out.nu;
    o.pressure=physics.P0+out.deltaP;if(!(o.pressure>0))return out;
    o.U=input.U;double residual=out.deltaP-input.mean_deltaP;
    double x=o.point.position[0],y=o.point.position[1],gx=0,gy=0;
    if(physics.gravity==1){double z=pow(x*x+y*y+physics.eps*physics.eps,-1.5);gx=-x*z;gy=-y*z;}
    if(physics.gravity==2){gx=-x;gy=-y;}
    double sigma[3]={-residual,0,-residual};
    for(int begin=0;begin<T.n;begin+=12){int n=T.n-begin;if(n>12)n=12;double A[12],Ax[12],Ay[12],force[2][12]={},u=0,g=0;
        for(int j=0;j<n;j++){A[j]=N[begin+j].v;Ax[j]=N[begin+j].x;Ay[j]=N[begin+j].y;}
        if(!gfs_material12_stress_at_point(n,&o.point,A,Ax,Ay,q.rho,q.weight,sigma,gx,gy,force,&u,&g))return out;
        if(begin==0)o.gpower=g;for(int j=0;j<n;j++)for(int d=0;d<2;d++)o.force[d][begin+j]=force[d][j];
    }
    out.nudot=q.weight*(o.point.J*(o.point.grad_v[0]+o.point.grad_v[3])+input.deltaJdot);
    o.Udot=-o.pressure*out.nudot;o.gtorque=q.weight*q.rho*(x*gy-y*gx);o.valid=isfinite(o.Udot);return out;
}
extern "C" int gp_cuda_evaluate(const STTriangle*,size_t,const STCoefficient*,size_t,const STPoint*,size_t,const STPhysics*,const GPInput*,GPOutput*,double*);
#undef GP_HD
#endif
