#ifndef GFS_KEP_ES_ALE_H
#define GFS_KEP_ES_ALE_H
/*
 * One complete kinetic-energy-preserving, entropy-stable ALE face flux.
 * The centered part is Chandrashekar's cell-mean KEP/EC flux. A symmetric
 * rank-two acoustic entropy-variable dissipation replaces the experimental
 * compression-stress sensor. Its coefficient is frozen from the two fluid
 * states in their arithmetic-mean normal frame, so the mass flux remains
 * affine in the actual ALE face speed and retains the scalar SPD weight solve.
 * This is an algebraic flux candidate, not a qualified evolution method.
 */
#include "gfs_kep_ale.h"
#ifdef __CUDACC__
#define GFS_KES_HD static __host__ __device__ inline
#else
#define GFS_KES_HD static inline
#endif

GFS_KES_HD void gfs_kep_entropy_variables(int dim,GfsFvPrimitive p,
        double gamma,double v[5]){
    for(int c=0;c<5;c++)v[c]=0;
    double beta=.5*p.rho/p.p,s=log(p.p)-gamma*log(p.rho),u2=0;
    for(int a=0;a<dim;a++){
        u2+=p.u[a]*p.u[a];
        v[a+1]=2*beta*p.u[a];
    }
    v[0]=(gamma-s)/(gamma-1)-beta*u2;
    v[4]=-2*beta;
}

GFS_KES_HD int gfs_kep_es_flux(int dim,GfsFvPrimitive l,GfsFvPrimitive r,
        const double normal[3],double gamma,double wn,double flux[5],
        double *density,double *dissipation_speed){
    if(!gfs_kep_flux(dim,l,r,normal,gamma,wn,0,flux,density))return 0;
    double beta_l=.5*l.rho/l.p,beta_r=.5*r.rho/r.p;
    double beta=gfs_kep_logmean(beta_l,beta_r);
    if(!(beta>0)||!isfinite(beta))return 0;
    double ubar[3]={0},un=0,ul=0,ur=0;
    for(int a=0;a<dim;a++){
        ubar[a]=.5*l.u[a]+.5*r.u[a];
        un+=ubar[a]*normal[a];
        ul+=l.u[a]*normal[a];ur+=r.u[a]*normal[a];
    }
    double cl=sqrt(gamma*l.p/l.rho),cr=sqrt(gamma*r.p/r.rho);
    double lambda=fmax(fabs(ul-un)+cl,fabs(ur-un)+cr);
    double a=sqrt(gamma/(2*beta)),u2=0;
    for(int k=0;k<dim;k++)u2+=ubar[k]*ubar[k];
    double H=a*a/(gamma-1)+.5*u2;
    if(!(lambda>0)||!(a>0)||!isfinite(lambda)||!isfinite(H))return 0;

    double vl[5],vr[5],dv[5],rm[5]={0},rp[5]={0};
    gfs_kep_entropy_variables(dim,l,gamma,vl);
    gfs_kep_entropy_variables(dim,r,gamma,vr);
    for(int c=0;c<5;c++)dv[c]=vr[c]-vl[c];
    double am=dv[0],ap=dv[0];
    for(int k=0;k<dim;k++){
        rm[k+1]=ubar[k]-a*normal[k];
        rp[k+1]=ubar[k]+a*normal[k];
        am+=rm[k+1]*dv[k+1];
        ap+=rp[k+1]*dv[k+1];
    }
    double unbar=un;
    rm[0]=rp[0]=1;
    rm[4]=H-a*unbar;rp[4]=H+a*unbar;
    am+=rm[4]*dv[4];ap+=rp[4]*dv[4];
    double coefficient=gfs_kep_logmean(l.rho,r.rho)*lambda/(2*gamma);
    if(!(coefficient>0)||!isfinite(coefficient))return 0;
    for(int c=0;c<5;c++){
        double qdv=coefficient*(rm[c]*am+rp[c]*ap);
        flux[c]-=.5*qdv;
        if(!isfinite(flux[c]))return 0;
    }
    if(dissipation_speed)*dissipation_speed=lambda;
    return 1;
}
#undef GFS_KES_HD
#endif
