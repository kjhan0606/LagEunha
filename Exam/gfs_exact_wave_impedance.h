#ifndef GFS_EXACT_WAVE_IMPEDANCE_H
#define GFS_EXACT_WAVE_IMPEDANCE_H
#include <math.h>
#ifdef __CUDACC__
#define GFS_WAVE_HD static __host__ __device__ inline
#else
#define GFS_WAVE_HD static inline
#endif
/* Secant of the exact ideal-gas wave curve at a supplied nonnegative star
 * pressure. Does not solve the Riemann problem or select a vacuum face speed.
 * Callers must reject a vacuum fan before assuming a shared material face.
 * This helper is not enabled in a production force/integrator path. */
GFS_WAVE_HD int gfs_exact_wave_impedance(double rho,double p,double star,
        double gamma,double *out)
{
    if(!out||!(rho>0)||!(p>0)||!(star>=0)||!(gamma>1)||
       !isfinite(rho)||!isfinite(p)||!isfinite(star)||!isfinite(gamma))return 0;
    double c=sqrt(gamma*p/rho),z;
    if(!(c>0)||!isfinite(c))return 0;
    if(star>p){
        double A=2/((gamma+1)*rho),B=(gamma-1)/(gamma+1)*p;
        z=sqrt((star+B)/A);
    }else if(star==p){
        z=rho*c;
    }else{
        double a=(gamma-1)/(2*gamma);
        if(star==0)z=rho*c*a;
        else{
            /* log1p/expm1 avoid cancellation for a weak rarefaction. */
            double r=(star-p)/p;
            double ell=r>-.5?log1p(r):log(star)-log(p);
            z=rho*c*a*expm1(ell)/expm1(a*ell);
        }
    }
    if(!(z>0)||!isfinite(z))return 0;
    *out=z;return 1;
}
#undef GFS_WAVE_HD
#endif
