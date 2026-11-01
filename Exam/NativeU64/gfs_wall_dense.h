#ifndef GFS_WALL_DENSE_H
#define GFS_WALL_DENSE_H
#include <math.h>
#include <float.h>

/* Cubic continuous extension of classical RK4. This proposes the first
 * outward event; native RK trials must subsequently refine its bracket.
 * Its dense accuracy is third order, not fourth-order event certification.
 */
typedef struct { double a0,a1,a2,a3; } GfsWallCubic;

static inline GfsWallCubic gfs_wall_cubic(double gap,double k1,double k2,double k3,double k4) {
    GfsWallCubic p={gap,k1,-1.5*k1+k2+k3-.5*k4,
                   (2./3.)*(k1-k2-k3+k4)};
    return p;
}
static inline double gfs_wall_value(GfsWallCubic p,double x) {
    return ((p.a3*x+p.a2)*x+p.a1)*x+p.a0;
}
static inline double gfs_wall_derivative(GfsWallCubic p,double x) {
    return (3.*p.a3*x+2.*p.a2)*x+p.a1;
}
static inline void gfs_wall_add_break(double *t,int *n,double x) {
    if(isfinite(x) && x>0. && x<1.) t[(*n)++]=x;
}
/* Return 1 for a first outward crossing, 0 for no resolved crossing, -1 for
 * invalid coefficients. Monotonic subdivisions find out-and-back events.
 * A small initial positive roundoff offset moving inward is not a new hit.
 */
static inline int gfs_wall_first_outward(GfsWallCubic p,double *root) {
    if(!isfinite(p.a0)||!isfinite(p.a1)||!isfinite(p.a2)||!isfinite(p.a3)||!root) return -1;
    double t[4]={0.,1.,0.,0.};int n=2;
    const double a=3.*p.a3,b=2.*p.a2,c=p.a1;
    if(a==0.) {
        if(b!=0.) gfs_wall_add_break(t,&n,-c/b);
    } else {
        const double d=b*b-4.*a*c;
        if(!isfinite(d)) return -1;
        if(d>=0.) {
            const double q=-.5*(b+copysign(sqrt(d),b));
            if(q==0.) gfs_wall_add_break(t,&n,-b/(2.*a));
            else {gfs_wall_add_break(t,&n,q/a);gfs_wall_add_break(t,&n,c/q);}
        }
    }
    for(int i=1;i<n;i++) for(int j=i;j>0 && t[j]<t[j-1];j--) {
        double q=t[j];t[j]=t[j-1];t[j-1]=q;
    }
    if(p.a0>=0. && p.a1>0.) {*root=0.;return 1;}
    for(int i=0;i+1<n;i++) {
        double lo=t[i],hi=t[i+1];
        const double left=gfs_wall_value(p,lo),right=gfs_wall_value(p,hi);
        if(left<=0. && (right>0. || (right==0. && gfs_wall_derivative(p,hi)>0.))) {
            if(left==0. && gfs_wall_derivative(p,lo)>0.) {*root=lo;return 1;}
            if(right==0.) {*root=hi;return 1;}
            for(int k=0;k<64;k++) {
                const double mid=.5*(lo+hi);
                if(mid==lo || mid==hi) break;
                if(gfs_wall_value(p,mid)>0.) hi=mid;else lo=mid;
            }
            *root=.5*(lo+hi);return 1;
        }
    }
    return 0;
}
#endif
