#ifndef GFS_MOMENT_INITIAL_H
#define GFS_MOMENT_INITIAL_H
#include "gfs_moment_density.h"
typedef struct {
    int kind; /* 1: exact affine patch, 2: existing smooth native Kepler IC */
    double gamma,rho,pressure,cx,cy,eps,rin,rout,vx,vy,B[4];
} GfsMomentIC;
typedef struct {
    GfsMomentPayload payload;
    GfsMomentState coarse;
    double thermal,thermal_coarse,nonaffine_kinetic,min_density;
    double moment_residual,energy_identity_error;
    int valid,iterations,restart;
} GfsMomentInitialAudit;

/* Same continuous field as Exam/KH/util.c's default smooth Kepler branch,
 * with GM=1 and constant P. Hopkins-profile mode is deliberately unsupported
 * here until its separate IC contract is supplied. This is fresh integration,
 * never energy added to an evolved checkpoint. */
GFS_MOM_HD static inline void gfs_moment_ic_value(const GfsMomentIC *ic,
    double x,double y,double *rho,double *vx,double *vy){
    double dx=x-ic->cx,dy=y-ic->cy;
    if(ic->kind==1){
        *rho=ic->rho;*vx=ic->vx+ic->B[0]*dx+ic->B[1]*dy;
        *vy=ic->vy+ic->B[2]*dx+ic->B[3]*dy;
    }else{
        double r=hypot(dx,dy);
        double inner=.5*(1+tanh((r-ic->rin)/.25));
        double outer=.5*(1+tanh((ic->rout-r)/.25));
        *rho=.001+ic->rin/(r+1.e-6)*inner*outer;
        double omega=pow(r*r+ic->eps*ic->eps,-.75);
        *vx=-dy*omega;*vy=dx*omega;
    }
}

GFS_MOM_HD static inline void gfs_moment_integrate_ic(GfsMomentPolygon poly,
    const GfsMomentDensity *geom,const GfsMomentIC *ic,double gx,double gy,
    GfsMomentState *s,double *thermal){
    gfs_moment_zero(s);*thermal=0;
    for(int k=0;k<gfs_moment_nquad(poly);k++){
        double x,y,w,rho,vx,vy;gfs_moment_quad(poly,geom,k,&x,&y,&w);
        if(w==0)continue;
        gfs_moment_ic_value(ic,gx+x,gy+y,&rho,&vx,&vy);
        double u=w*ic->pressure/(ic->gamma-1);
        gfs_moment_accumulate(s,w*rho,x,y,vx,vy,u);*thermal+=u;
    }
}
#endif
