#ifndef GFS_KEP_ALE_H
#define GFS_KEP_ALE_H
/* Cell-pair kinetic-compatible, entropy-conservative ALE Euler core.
 * Chandrashekar's spatial flux plus its temporal entropy state. These
 * are CELL MEANS, not independently reconstructed face states. This
 * primitive does not qualify a timestep or production PDE solver.
 * Optional normal compression stress has nonnegative physical entropy
 * production and changes momentum and energy together; mass is unchanged.
 */
#include "gfs_fv_ale.h"
#ifdef __CUDACC__
#define GFS_KEP_HD static __host__ __device__ inline
#else
#define GFS_KEP_HD static inline
#endif
GFS_KEP_HD double gfs_kep_logmean(double a,double b){
    if(!(a>0)||!(b>0)||!isfinite(a)||!isfinite(b))return -1;
    if(a==b)return a;
    double mean=.5*a+.5*b,f=(.5*b-.5*a)/mean,s=f*f;
    if(s<1e-4){
        double ratio=1+s*(1./3+s*(1./5+s*(1./7+s*(1./9+s/11))));
        return mean/ratio;
    }
    return (b-a)/(log(b)-log(a));
}
GFS_KEP_HD int gfs_kep_flux(int dim,GfsFvPrimitive l,GfsFvPrimitive r,
        const double normal[3],double gamma,double wn,double compression,
        double flux[5],double *density){
    for(int c=0;c<5;c++)flux[c]=0;
    *density=0;
    if(!gfs_fv_input(dim,l)||!gfs_fv_input(dim,r)||!(gamma>1)||!isfinite(gamma)
       ||!isfinite(wn)||!(compression>=0)||!isfinite(compression))return 0;
    double rho=gfs_kep_logmean(l.rho,r.rho),bl=.5*(l.rho/l.p),br=.5*(r.rho/r.p);
    double beta=gfs_kep_logmean(bl,br),mean_beta=.5*bl+.5*br;
    if(!(rho>0)||!(beta>0)||!isfinite(mean_beta))return 0;
    double mean[3]={0},un=0,cross=0,du=0,norm=0;
    for(int a=0;a<dim;a++){
        if(!isfinite(normal[a]))return 0;
        norm+=normal[a]*normal[a];mean[a]=.5*l.u[a]+.5*r.u[a];
        un+=mean[a]*normal[a];cross+=l.u[a]*r.u[a];du+=(r.u[a]-l.u[a])*normal[a];
    }
    if(fabs(norm-1)>128*DBL_EPSILON||(compression>0&&du>=0))return 0;
    double pressure=(.5*l.rho+.5*r.rho)/(2*mean_beta);
    flux[0]=rho*(un-wn);
    for(int a=0;a<dim;a++)flux[a+1]=flux[0]*mean[a]+(pressure+compression)*normal[a];
    flux[4]=flux[0]*(1/(2*(gamma-1)*beta)+.5*cross)+(pressure+compression)*un;
    for(int c=0;c<5;c++)if(!isfinite(flux[c]))return 0;
    *density=rho;return 1;
}
#undef GFS_KEP_HD
#endif
