#ifndef GFS_MOMENT_H
#define GFS_MOMENT_H
/* Opt-in affine cell moments. All positions below are measured from a fixed
 * per-label reference origin O, never from the evolving generator/centroid.
 * D_ab = integral rho (x-O)_a v_b. The state therefore mixes linearly in SSP
 * stages and migrates without a moving-origin source. No physical damping.
 */
#include <math.h>
#include <float.h>
#ifdef __CUDACC__
#define GFS_MOM_HD __host__ __device__
#else
#define GFS_MOM_HD
#endif
enum { GFM_M, GFM_QX, GFM_QY, GFM_QXX, GFM_QXY, GFM_QYY,
       GFM_PX, GFM_PY, GFM_DXX, GFM_DXY, GFM_DYX, GFM_DYY, GFM_E,
       GFM_NSTATE };
typedef struct { double a[GFM_NSTATE]; } GfsMomentState;
/* Exactly 24 doubles: the old RK4 scratch prefix, excluding w2backup. This
 * tagged payload is used only by the new mode; legacy RK/SSP cannot touch it.
 * No native particle ABI/stride change: halos, pmigrate and KH full-stride
 * snapshots carry the complete state. Access through memcpy, not type punning.
 */
#define GFS_MOMENT_TAG 2026102001.0
typedef struct {
    GfsMomentState s;
    double m0, mu, w, ox, oy, tag, lambda[5];
} GfsMomentPayload;
typedef struct {
    double c[2], S[3], C[4], B[4], v[2];
    double Kmean, Kaffine, U, subtraction_scale;
} GfsMomentPrimitive;

/* Adjacent exact power supports n.x=b. Recomputing an intersection avoids
 * the one-ulp backwards edge left by repeated segment interpolation at a
 * cocircular vertex. No edge-length threshold or weight clamp is used. */
GFS_MOM_HD static inline int gfs_moment_support_intersection(
    double ax,double ay,double ab,double bx,double by,double bb,double *x,double *y){
    double det=fma(ax,by,-ay*bx);
    if(det==0 || !isfinite(det))return 0;
    *x=fma(ab,by,-ay*bb)/det;*y=fma(ax,bb,-ab*bx)/det;
    return isfinite(*x) && isfinite(*y);
}

GFS_MOM_HD static inline void gfs_moment_zero(GfsMomentState *s){
    for(int k=0;k<GFM_NSTATE;k++) s->a[k]=0;
}

GFS_MOM_HD static inline int gfs_moment_primitive(const GfsMomentState *s,
    GfsMomentPrimitive *r){
    const double *a=s->a; double m=a[GFM_M];
    for(int k=0;k<GFM_NSTATE;k++) if(!isfinite(a[k])) return 0;
    if(!(m>0)) return 0;
    r->c[0]=a[GFM_QX]/m; r->c[1]=a[GFM_QY]/m;
    r->v[0]=a[GFM_PX]/m; r->v[1]=a[GFM_PY]/m;
    r->S[0]=fma(-a[GFM_QX],r->c[0],a[GFM_QXX]);
    r->S[1]=fma(-a[GFM_QX],r->c[1],a[GFM_QXY]);
    r->S[2]=fma(-a[GFM_QY],r->c[1],a[GFM_QYY]);
    double det=fma(r->S[0],r->S[2],-r->S[1]*r->S[1]);
    if(!(r->S[0]>0) || !(r->S[2]>0) || !(det>0) || !isfinite(det)) return 0;
    r->C[0]=fma(-a[GFM_QX],r->v[0],a[GFM_DXX]);
    r->C[1]=fma(-a[GFM_QX],r->v[1],a[GFM_DXY]);
    r->C[2]=fma(-a[GFM_QY],r->v[0],a[GFM_DYX]);
    r->C[3]=fma(-a[GFM_QY],r->v[1],a[GFM_DYY]);
    double ixx=r->S[2]/det,ixy=-r->S[1]/det,iyy=r->S[0]/det;
    r->B[0]=fma(r->C[0],ixx,r->C[2]*ixy);
    r->B[1]=fma(r->C[0],ixy,r->C[2]*iyy);
    r->B[2]=fma(r->C[1],ixx,r->C[3]*ixy);
    r->B[3]=fma(r->C[1],ixy,r->C[3]*iyy);
    r->Kmean=.5*(a[GFM_PX]*r->v[0]+a[GFM_PY]*r->v[1]);
    r->Kaffine=.5*(r->B[0]*r->C[0]+r->B[1]*r->C[2]+
                        r->B[2]*r->C[1]+r->B[3]*r->C[3]);
    r->U=(a[GFM_E]-r->Kmean)-r->Kaffine;
    r->subtraction_scale=DBL_EPSILON*(fabs(a[GFM_E])+r->Kmean+fabs(r->Kaffine));
    return r->Kaffine>=0 && r->U>0 && isfinite(r->U);
}

GFS_MOM_HD static inline void gfs_moment_accumulate(GfsMomentState *s,
    double dm,double x,double y,double vx,double vy,double eint){
    double *a=s->a;
    a[GFM_M]+=dm; a[GFM_QX]+=dm*x; a[GFM_QY]+=dm*y;
    a[GFM_QXX]+=dm*x*x; a[GFM_QXY]+=dm*x*y; a[GFM_QYY]+=dm*y*y;
    a[GFM_PX]+=dm*vx; a[GFM_PY]+=dm*vy;
    a[GFM_DXX]+=dm*x*vx; a[GFM_DXY]+=dm*x*vy;
    a[GFM_DYX]+=dm*y*vx; a[GFM_DYY]+=dm*y*vy;
    a[GFM_E]+=eint+.5*dm*(vx*vx+vy*vy);
}

/* O' = O + delta. Physical angular momentum and energy are unchanged. */
GFS_MOM_HD static inline GfsMomentState gfs_moment_shift(GfsMomentState s,
    double dx,double dy){
    double *a=s.a,m=a[GFM_M],qx=a[GFM_QX],qy=a[GFM_QY];
    a[GFM_QXX]+= -2*dx*qx+m*dx*dx;
    a[GFM_QXY]+= -dx*qy-dy*qx+m*dx*dy;
    a[GFM_QYY]+= -2*dy*qy+m*dy*dy;
    a[GFM_QX]-=m*dx; a[GFM_QY]-=m*dy;
    a[GFM_DXX]-=dx*a[GFM_PX]; a[GFM_DXY]-=dx*a[GFM_PY];
    a[GFM_DYX]-=dy*a[GFM_PX]; a[GFM_DYY]-=dy*a[GFM_PY];
    return s;
}

GFS_MOM_HD static inline GfsMomentState gfs_moment_boost(GfsMomentState s,
    double ux,double uy){
    double *a=s.a;
    a[GFM_E]+=ux*a[GFM_PX]+uy*a[GFM_PY]+.5*a[GFM_M]*(ux*ux+uy*uy);
    a[GFM_PX]+=a[GFM_M]*ux; a[GFM_PY]+=a[GFM_M]*uy;
    a[GFM_DXX]+=a[GFM_QX]*ux; a[GFM_DXY]+=a[GFM_QX]*uy;
    a[GFM_DYX]+=a[GFM_QY]*ux; a[GFM_DYY]+=a[GFM_QY]*uy;
    return s;
}

GFS_MOM_HD static inline double gfs_moment_Lz(const GfsMomentState *s,
    double ox,double oy){
    return s->a[GFM_DXY]-s->a[GFM_DYX]+ox*s->a[GFM_PY]-oy*s->a[GFM_PX];
}

/* One actual quadrature point, outward Euler ALE flux [mass,px,py,E].
 * Neighbours use the SAME canonical flux and physical point, opposite signs,
 * but express its coordinates in their own fixed reference origins.
 */
GFS_MOM_HD static inline void gfs_moment_face_rhs(GfsMomentState *rhs,
    double x,double y,const double flux[4],double signed_area){
    double fm=signed_area*flux[0],fx=signed_area*flux[1];
    double fy=signed_area*flux[2],fe=signed_area*flux[3],*a=rhs->a;
    a[GFM_M]-=fm; a[GFM_QX]-=x*fm; a[GFM_QY]-=y*fm;
    a[GFM_QXX]-=x*x*fm; a[GFM_QXY]-=x*y*fm; a[GFM_QYY]-=y*y*fm;
    a[GFM_PX]-=fx; a[GFM_PY]-=fy;
    a[GFM_DXX]-=x*fx; a[GFM_DXY]-=x*fy;
    a[GFM_DYX]-=y*fx; a[GFM_DYY]-=y*fy; a[GFM_E]-=fe;
}

/* Exact affine Euler volume terms; pressure integral=(gamma-1) U. Gravity
 * is added separately with the same positive density/velocity quadrature.
 */
GFS_MOM_HD static inline void gfs_moment_euler_source(const GfsMomentState *s,
    const GfsMomentPrimitive *r,double gamma,GfsMomentState *rhs){
    double *a=rhs->a,px=s->a[GFM_PX],py=s->a[GFM_PY];
    a[GFM_QX]+=px; a[GFM_QY]+=py;
    a[GFM_QXX]+=2*s->a[GFM_DXX];
    a[GFM_QXY]+=s->a[GFM_DXY]+s->a[GFM_DYX];
    a[GFM_QYY]+=2*s->a[GFM_DYY];
    double pii=(gamma-1)*r->U;
    a[GFM_DXX]+=px*r->v[0]+r->B[0]*r->C[0]+r->B[1]*r->C[2]+pii;
    a[GFM_DXY]+=px*r->v[1]+r->B[0]*r->C[1]+r->B[1]*r->C[3];
    a[GFM_DYX]+=py*r->v[0]+r->B[2]*r->C[0]+r->B[3]*r->C[2];
    a[GFM_DYY]+=py*r->v[1]+r->B[2]*r->C[1]+r->B[3]*r->C[3]+pii;
}

GFS_MOM_HD static inline void gfs_moment_gravity_source(GfsMomentState *rhs,
    double dm,double x,double y,double vx,double vy,double gx,double gy){
    double fx=dm*gx,fy=dm*gy,*a=rhs->a;
    a[GFM_PX]+=fx; a[GFM_PY]+=fy;
    a[GFM_DXX]+=x*fx; a[GFM_DXY]+=x*fy;
    a[GFM_DYX]+=y*fx; a[GFM_DYY]+=y*fy;
    a[GFM_E]+=vx*fx+vy*fy;
}
#endif
