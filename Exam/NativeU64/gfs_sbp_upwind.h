#ifndef GFS_SBP_UPWIND_H
#define GFS_SBP_UPWIND_H
#include "gfs_sbp_flux.h"
#ifdef __CUDACC__
#define GFS_UP_HD __host__ __device__
#else
#define GFS_UP_HD
#endif
/* Candidate selected by ASTRA_VOLUME_UPWIND_PLAN.md. Complete Euler HLL
 * fluctuation evaluated on two midpoint reconstructed primitive traces.
 * Stored nodal conserved states are not replaced by those traces. Wave bounds
 * include both endpoint and Roe acoustic speeds; no general positivity proof
 * is asserted for this high-order corrected operator. Invalid traces or HLL
 * intermediate states fail explicitly. No tunable dissipation coefficient.
 */
GFS_UP_HD static inline int gfs_sbp_hll_flux(GfsSbpPrimitive l,GfsSbpPrimitive r,
    double nx,double ny,double wn,double gamma,double f[4]){
    double check[4];if(!gfs_sbp_volume_flux(l,l,nx,ny,wn,gamma,check) ||
        !gfs_sbp_volume_flux(r,r,nx,ny,wn,gamma,check))return 0;
    double norm=hypot(nx,ny);if(norm==0){for(int c=0;c<4;c++)f[c]=0;return wn==0;}
    nx/=norm;ny/=norm;wn/=norm;
    double ul[4]={l.rho,l.rho*l.vx,l.rho*l.vy,l.p/(gamma-1)+.5*l.rho*(l.vx*l.vx+l.vy*l.vy)};
    double ur[4]={r.rho,r.rho*r.vx,r.rho*r.vy,r.p/(gamma-1)+.5*r.rho*(r.vx*r.vx+r.vy*r.vy)};
    double vl=l.vx*nx+l.vy*ny,vr=r.vx*nx+r.vy*ny,al=vl-wn,ar=vr-wn;
    double fl[4]={l.rho*al,ul[1]*al+l.p*nx,ul[2]*al+l.p*ny,ul[3]*al+l.p*vl};
    double fr[4]={r.rho*ar,ur[1]*ar+r.p*nx,ur[2]*ar+r.p*ny,ur[3]*ar+r.p*vr};
    double cl=sqrt(gamma*l.p/l.rho),cr=sqrt(gamma*r.p/r.rho),rl=sqrt(l.rho),rr=sqrt(r.rho);
    double vx=(rl*l.vx+rr*r.vx)/(rl+rr),vy=(rl*l.vy+rr*r.vy)/(rl+rr);
    double H=(rl*(ul[3]+l.p)/l.rho+rr*(ur[3]+r.p)/r.rho)/(rl+rr);
    double cs2=(gamma-1)*(H-.5*(vx*vx+vy*vy));if(!(cs2>0) || !isfinite(cs2))return 0;
    double roe=vx*nx+vy*ny-wn,cs=sqrt(cs2);
    double sl=fmin(fmin(al-cl,ar-cr),roe-cs),sr=fmax(fmax(al+cl,ar+cr),roe+cs);
    if(sl>=0){for(int c=0;c<4;c++)f[c]=norm*fl[c];}
    else if(sr<=0){for(int c=0;c<4;c++)f[c]=norm*fr[c];}
    else{
        double middle[4];for(int c=0;c<4;c++){
            middle[c]=(sr*ur[c]-sl*ul[c]-(fr[c]-fl[c]))/(sr-sl);
            f[c]=norm*(sr*fl[c]-sl*fr[c]+sl*sr*(ur[c]-ul[c]))/(sr-sl);
        }
        if(!(middle[0]>0) || !(middle[3]-.5*(middle[1]*middle[1]+middle[2]*middle[2])/middle[0]>0))return 0;
    }
    for(int c=0;c<4;c++)if(!isfinite(f[c]))return 0;return 1;
}
GFS_UP_HD static inline int gfs_sbp_upwind_correction(GfsSbpPrimitive l_trace,GfsSbpPrimitive r_trace,
    double nx,double ny,double wn,double gamma,double correction[4]){
    double smooth[4],hll[4];if(!gfs_sbp_volume_flux(l_trace,r_trace,nx,ny,wn,gamma,smooth) ||
        !gfs_sbp_hll_flux(l_trace,r_trace,nx,ny,wn,gamma,hll))return 0;
    for(int c=0;c<4;c++)correction[c]=hll[c]-smooth[c];return 1;
}
#undef GFS_UP_HD
#endif
