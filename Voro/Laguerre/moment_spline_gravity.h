#ifndef LAGFORCE_MOMENT_SPLINE_GRAVITY_H
#define LAGFORCE_MOMENT_SPLINE_GRAVITY_H
#include <math.h>
#ifdef __CUDACC__
#define WG_INLINE static __host__ __device__ __forceinline__
#else
#define WG_INLINE static inline
#endif

/* Experimental material-quadrature source, not production cell evolution.
 * Input dm,x,y,vx,vy,rc (coordinates relative to the central mass, GM=1).
 * Output: time derivative of the 17 raw moments in weight_remap_moments.h,
 * followed by gx,gy,phi. Pressure, boundaries and weight remap are absent.
 * Spline formula matches Exam/exam.c kepler_spline_mass/kepler_phi.
 */
WG_INLINE int wg_material_source(const double *a,double *b) {
    double m=a[0],x=a[1],y=a[2],vx=a[3],vy=a[4],rc=a[5];
    if (!(m>0.0 && rc>0.0)) return 1;
    double r=sqrt(x*x+y*y),u=r/rc,enclosed,phi;
    if(u>=1.0) {
        enclosed=1.0; phi=-1.0/r;
    } else if(u<=0.5) {
        enclosed=u*u*u*(32.0/3-192.0/5*u*u+32*u*u*u);
        phi=(16.0/3*u*u-48.0/5*u*u*u*u+32.0/5*u*u*u*u*u-14.0/5)/rc;
    } else {
        enclosed=-1.0/15+u*u*u*(64.0/3-48*u+192.0/5*u*u-32.0/3*u*u*u);
        phi=(1.0/(15*u)+32.0/3*u*u-16*u*u*u+48.0/5*u*u*u*u-32.0/15*u*u*u*u*u-16.0/5)/rc;
    }
    double factor=r>0.0 ? -enclosed/(r*r*r) : 0.0;
    double gx=factor*x,gy=factor*y,work=m*(vx*gx+vy*gy);
    b[0]=0.0; b[1]=m*vx; b[2]=m*vy;
    b[3]=2*m*x*vx; b[4]=m*(vx*y+x*vy); b[5]=2*m*y*vy;
    b[6]=m*gx; b[7]=m*gy;
    b[8]=m*(vx*vx+x*gx); b[9]=m*(vx*vy+x*gy);
    b[10]=m*(vy*vx+y*gx); b[11]=m*(vy*vy+y*gy);
    b[12]=2*m*vx*gx; b[13]=m*(gx*vy+vx*gy); b[14]=2*m*vy*gy;
    b[15]=work; b[16]=-work; b[17]=gx; b[18]=gy; b[19]=phi;
    return 0;
}
#undef WG_INLINE
#endif
