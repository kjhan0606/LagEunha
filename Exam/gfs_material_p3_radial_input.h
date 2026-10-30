#ifndef GFS_MATERIAL_P3_RADIAL_INPUT_H
#define GFS_MATERIAL_P3_RADIAL_INPUT_H
#include <math.h>
#include <float.h>
#ifdef __CUDACC__
#define GFS_P3_INPUT_HD __host__ __device__
#else
#define GFS_P3_INPUT_HD
#endif
/* Accurate INITIAL v0-T3(v0), v0=J X (r^2+eps^2)^(-3/4).
 * This is the SAME initial velocity interpolation residual as direct
 * subtraction: no gravity basis, force correction, future orbit or reset.
 */
GFS_P3_INPUT_HD static inline int gfs_material_p3_radial_remainder(
    const double anchor[2],double x,double y,double eps,double q[2],double Dq[4]){
    double dx=x-anchor[0],dy=y-anchor[1],z=anchor[0]*anchor[0]+anchor[1]*anchor[1]+eps*eps;
    if(!(z>0))return 0;
    double L=2*(anchor[0]*dx+anchor[1]*dy)/z,Q=(dx*dx+dy*dy)/z,s=L+Q,abs_s=fabs(s);
    if(!(abs_s<.75))return 0;
    const double c2=.65625,c3=-.6015625,c4=.56396484375;
    double tail=0,dtail=0;
    if(s!=0){
        double term=c4*s*s*s*s;tail=term;dtail=4*term/s;int converged=0;
        for(int n=4;n<256;n++){
            double value_bound=fabs(term)*abs_s/(1-abs_s);
            double derivative_bound=fabs(term)*((n+1)/(1-abs_s)+abs_s/((1-abs_s)*(1-abs_s)));
            if(value_bound<=DBL_EPSILON*fabs(tail) && derivative_bound<=DBL_EPSILON*fabs(dtail)){converged=1;break;}
            term*=(-.75-n)/(n+1)*s;tail+=term;dtail+=(n+1)*term/s;
        }
        if(!converged)return 0;
    }
    double A=c2*Q*Q+c3*(3*L*L*Q+3*L*Q*Q+Q*Q*Q)+tail;
    double B=c2*(2*L*Q+Q*Q)+c3*s*s*s+tail;
    double AL=c3*(6*L*Q+3*Q*Q)+dtail,AQ=2*c2*Q+c3*(3*L*L+6*L*Q+3*Q*Q)+dtail;
    double BL=2*c2*Q+3*c3*s*s+dtail,BQ=2*c2*(L+Q)+3*c3*s*s+dtail;
    double omega=pow(z,-.75),ja[2]={-anchor[1],anchor[0]},jd[2]={-dy,dx};
    for(int k=0;k<2;k++){
        q[k]=omega*(ja[k]*A+jd[k]*B);
        for(int d=0;d<2;d++){
            double ld=2*anchor[d]/z,qd=2*(d?dy:dx)/z,dA=AL*ld+AQ*qd,dB=BL*ld+BQ*qd,dJ=(k==0 && d==1)?-1:(k==1 && d==0)?1:0;
            Dq[2*k+d]=omega*(ja[k]*dA+jd[k]*dB+dJ*B);
        }
    }return 1;
}
#undef GFS_P3_INPUT_HD
#endif
