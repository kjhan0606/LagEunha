#ifndef LAGFORCE_MOMENT_KINETIC_FACE_FLUX_H
#define LAGFORCE_MOMENT_KINETIC_FACE_FLUX_H
#include "moment_state_layout.h"

/* One planar face quadrature node, fixed four-sample kinetic closure.
 * Explicit input order (33 doubles): point[2],normal[2],scale,
 * vl[4][2],vr[4][2],wl[4],wr[4],eL,eR,Pstar,ustar_normal.
 * wl>=0,wr<=0 include density/4 and relative-normal split, NOT scale.
 * scale=face length*positive quadrature weight. e is thermal energy/mass.
 * Output is outward raw18 flux, applied once with opposite signs to cells.
 * Pi here contains transport only. Compatible pressure Pi/volume sources,
 * geometry, reconstruction, massroot, gravity and Psi are separate operators.
 * Caller validates unit normal/physical thermodynamics. No extra tolerance or
 * positivity assumption on e/Pstar is introduced in this algebra helper.
 */
enum LagKineticNodeIndex {
    LKN_X=0,LKN_N=2,LKN_SCALE=4,LKN_VL=5,LKN_VR=13,
    LKN_WL=21,LKN_WR=25,LKN_EL=29,LKN_ER=30,LKN_PSTAR=31,
    LKN_USTAR=32,LKN_COUNT=33
};
#ifdef __CUDACC__
#define LKF_INLINE static __host__ __device__ __forceinline__
#else
#define LKF_INLINE static inline
#endif
/* 0 success,1 invalid input,2 overflow. Failures leave output untouched. */
LKF_INLINE int lmk_face_node_flux(const double *a,double *out) {
    if(!a || !out)return 1;
    for(int k=0;k<LKN_COUNT;k++)if(!isfinite(a[k]))return 1;
    if(a[LKN_SCALE]<0.)return 1;
    double ml=0.,mr=0.,pl[2]={0.,0.},pr[2]={0.,0.};
    double pixx=0.,pixy=0.,piyy=0.,el=0.,er=0.;
    for(int k=0;k<4;k++){
        const double wl=a[LKN_WL+k],wr=a[LKN_WR+k];
        if(wl<0. || wr>0.)return 1;
        const double lx=a[LKN_VL+2*k],ly=a[LKN_VL+2*k+1];
        const double rx=a[LKN_VR+2*k],ry=a[LKN_VR+2*k+1];
        ml+=wl;mr+=wr;pl[0]+=wl*lx;pl[1]+=wl*ly;
        pr[0]+=wr*rx;pr[1]+=wr*ry;
        pixx+=wl*lx*lx+wr*rx*rx;
        pixy+=wl*lx*ly+wr*rx*ry;
        piyy+=wl*ly*ly+wr*ry*ry;
        el+=wl*(.5*(lx*lx+ly*ly)+a[LKN_EL]);
        er+=wr*(.5*(rx*rx+ry*ry)+a[LKN_ER]);
    }
    const double fm=ml+mr,x=a[0],y=a[1],scale=a[LKN_SCALE];
    const double px=pl[0]+pr[0]+a[LKN_PSTAR]*a[2];
    const double py=pl[1]+pr[1]+a[LKN_PSTAR]*a[3];
    double tmp[18]={fm,x*fm,y*fm,x*x*fm,x*y*fm,x*y*fm,y*y*fm,
        px,py,x*px,x*py,y*px,y*py,pixx,pixy,pixy,piyy,
        el+er+a[LKN_PSTAR]*a[LKN_USTAR]};
    for(int k=0;k<LMR_COUNT;k++){
        tmp[k]*=scale;if(!isfinite(tmp[k]))return 2;
    }
    for(int k=0;k<LMR_COUNT;k++)out[k]=tmp[k];
    return 0;
}
#undef LKF_INLINE
#endif
