#ifndef GFS_SBP_FLUX_H
#define GFS_SBP_FLUX_H
#include <math.h>
#ifdef __CUDACC__
#define GFS_SBP_HD __host__ __device__
#else
#define GFS_SBP_HD
#endif
/* Symmetric kinetic/pressure-compatible smooth VOLUME flux selected in
 * ASTRA_NODAL_SBP_PLAN.md. Actual nodal diagonal mass/SBP/GCL is required.
 * This is not an entropy-stable shock flux or a time integrator. Boundary
 * SAT uses the canonical complete Euler HLLC/HLL numerical flux separately.
 * n may be a unit direction or a signed metric vector; wn has matching units.
 */
typedef struct {double rho,p,vx,vy;} GfsSbpPrimitive;
GFS_SBP_HD static inline int gfs_sbp_volume_flux(GfsSbpPrimitive l,GfsSbpPrimitive r,
    double nx,double ny,double wn,double gamma,double f[4]){
    if(!(l.rho>0) || !(r.rho>0) || !(l.p>0) || !(r.p>0) || !(gamma>1) ||
       !isfinite(l.rho) || !isfinite(r.rho) || !isfinite(l.p) || !isfinite(r.p) ||
       !isfinite(l.vx) || !isfinite(r.vx) || !isfinite(l.vy) || !isfinite(r.vy) ||
       !isfinite(nx) || !isfinite(ny) || !isfinite(wn))return 0;
    double rho=.5*(l.rho+r.rho),p=.5*(l.p+r.p),vx=.5*(l.vx+r.vx),vy=.5*(l.vy+r.vy);
    double relative=vx*nx+vy*ny-wn,mass=rho*relative;
    f[0]=mass;f[1]=mass*vx+p*nx;f[2]=mass*vy+p*ny;
    f[3]=.5*mass*(l.vx*r.vx+l.vy*r.vy)+p/(gamma-1)*relative+
        .5*((l.p*r.vx+r.p*l.vx)*nx+(l.p*r.vy+r.p*l.vy)*ny);
    for(int c=0;c<4;c++)if(!isfinite(f[c]))return 0;return 1;
}
/* Exact directional derivative at fixed spatial metrics/mesh speed. Primitive
 * increments use the same layout but are signed, not admissible states. This
 * diagnoses the actual volume flux's local stability; it is not extra damping.
 */
GFS_SBP_HD static inline int gfs_sbp_volume_jvp(GfsSbpPrimitive l,GfsSbpPrimitive r,
    GfsSbpPrimitive dl,GfsSbpPrimitive dr,double nx,double ny,double wn,double gamma,double df[4]){
    double f[4];if(!gfs_sbp_volume_flux(l,r,nx,ny,wn,gamma,f))return 0;
    double rho=.5*(l.rho+r.rho),p=.5*(l.p+r.p),vx=.5*(l.vx+r.vx),vy=.5*(l.vy+r.vy);
    double d_rho=.5*(dl.rho+dr.rho),dp=.5*(dl.p+dr.p),dvx=.5*(dl.vx+dr.vx),dvy=.5*(dl.vy+dr.vy);
    double a=vx*nx+vy*ny-wn,da=dvx*nx+dvy*ny,dm=d_rho*a+rho*da;
    df[0]=dm;df[1]=dm*vx+f[0]*dvx+dp*nx;df[2]=dm*vy+f[0]*dvy+dp*ny;
    df[3]=.5*dm*(l.vx*r.vx+l.vy*r.vy)+.5*f[0]*(dl.vx*r.vx+l.vx*dr.vx+dl.vy*r.vy+l.vy*dr.vy)+
        dp/(gamma-1)*a+p/(gamma-1)*da+
        .5*((dl.p*r.vx+l.p*dr.vx+dr.p*l.vx+r.p*dl.vx)*nx+
            (dl.p*r.vy+l.p*dr.vy+dr.p*l.vy+r.p*dl.vy)*ny);
    for(int c=0;c<4;c++)if(!isfinite(df[c]))return 0;return 1;
}
#undef GFS_SBP_HD
#endif
