#ifndef GFS_POWER_PREDICATE2D_H
#define GFS_POWER_PREDICATE2D_H
#include <math.h>
/* Intel's default fast FP model can fold the error-free sum's low part to
 * zero. Preserve the specified rounding operations inside these helpers;
 * changing only the surrounding caller's precision is insufficient. */
#if !defined(__CUDACC__) && (defined(__clang__) || defined(__INTEL_COMPILER) || defined(_MSC_VER))
#pragma float_control(precise, on, push)
#define GFS_PRED_PRECISE_SCOPE 1
#endif
#ifdef __CUDACC__
#define GFS_PRED_HD static __host__ __device__ inline
#else
#define GFS_PRED_HD static inline
#endif
/* Compensated arithmetic for a FILTERED geometric predicate. This is not
 * arbitrary-precision arithmetic or an exact-sign certificate for all inputs.
 * The normal clipping path calls it only near its floating-point uncertainty
 * band. It introduces no physical tolerance, displacement or weight clamp. */
typedef struct {double hi,lo;} GfsPowerDD;
GFS_PRED_HD GfsPowerDD gfs_pdd_sum(double a,double b){
 double s=a+b,z=s-a;GfsPowerDD r={s,(a-(s-z))+(b-z)};return r;
}
GFS_PRED_HD GfsPowerDD gfs_pdd_value(double x){GfsPowerDD r={x,0};return r;}
GFS_PRED_HD GfsPowerDD gfs_pdd_neg(GfsPowerDD a){GfsPowerDD r={-a.hi,-a.lo};return r;}
GFS_PRED_HD GfsPowerDD gfs_pdd_add(GfsPowerDD a,GfsPowerDD b){
 GfsPowerDD s=gfs_pdd_sum(a.hi,b.hi),t=gfs_pdd_sum(a.lo,b.lo);
 s=gfs_pdd_sum(s.hi,s.lo+t.hi);return gfs_pdd_sum(s.hi,s.lo+t.lo);
}
GFS_PRED_HD GfsPowerDD gfs_pdd_sub(GfsPowerDD a,GfsPowerDD b){return gfs_pdd_add(a,gfs_pdd_neg(b));}
GFS_PRED_HD GfsPowerDD gfs_pdd_mul(GfsPowerDD a,GfsPowerDD b){
 double p=a.hi*b.hi,e=fma(a.hi,b.hi,-p);
 e+=(a.hi*b.lo+a.lo*b.hi);GfsPowerDD s=gfs_pdd_sum(p,e);
 return gfs_pdd_sum(s.hi,s.lo+a.lo*b.lo);
}
GFS_PRED_HD GfsPowerDD gfs_pdd_div(GfsPowerDD a,GfsPowerDD b){
 double q=a.hi/b.hi;
 GfsPowerDD r=gfs_pdd_sub(a,gfs_pdd_mul(b,gfs_pdd_value(q)));
 return gfs_pdd_sum(q,(r.hi+r.lo)/b.hi);
}
GFS_PRED_HD GfsPowerDD gfs_pdd_cross(GfsPowerDD a,GfsPowerDD b,GfsPowerDD c,GfsPowerDD d){
 return gfs_pdd_sub(gfs_pdd_mul(a,b),gfs_pdd_mul(c,d));
}
/* Three support sites define a vertex. Site0 is the owner, sites1/2 its
 * two incident supports. Site3 supplies the halfplane being tested.
 * Outputs are owner-relative vertex, signed power-plane residual, and sign.
 * Computing the residual from the support equations avoids testing rounded
 * Cartesian vertex coordinates against another rounded plane. */
GFS_PRED_HD int gfs_power_vertex_side_dd2d(const GfsPowerDD s[4][3],double vertex[2],double*residual,int*sign){
 GfsPowerDD x[3],y[3],h[3];
 for(int k=0;k<4;k++)for(int a=0;a<3;a++)if(!isfinite(s[k][a].hi)||!isfinite(s[k][a].lo))return 0;
 for(int k=0;k<3;k++){
  x[k]=gfs_pdd_sub(s[k+1][0],s[0][0]);y[k]=gfs_pdd_sub(s[k+1][1],s[0][1]);
  GfsPowerDD z=gfs_pdd_add(gfs_pdd_mul(x[k],x[k]),gfs_pdd_mul(y[k],y[k]));
  z=gfs_pdd_add(z,gfs_pdd_sub(s[0][2],s[k+1][2]));
  h[k]=gfs_pdd_mul(gfs_pdd_value(.5),z);
 }
 GfsPowerDD det=gfs_pdd_cross(x[0],y[1],y[0],x[1]);
 if(!isfinite(det.hi)||det.hi==0)return 0;
 GfsPowerDD nx=gfs_pdd_cross(h[0],y[1],y[0],h[1]);
 GfsPowerDD ny=gfs_pdd_cross(x[0],h[1],h[0],x[1]);
 GfsPowerDD num=gfs_pdd_sub(gfs_pdd_add(gfs_pdd_mul(x[2],nx),gfs_pdd_mul(y[2],ny)),gfs_pdd_mul(h[2],det));
 GfsPowerDD vx=gfs_pdd_div(nx,det),vy=gfs_pdd_div(ny,det),v=gfs_pdd_div(num,det);
 vertex[0]=vx.hi+vx.lo;vertex[1]=vy.hi+vy.lo;*residual=v.hi+v.lo;
 double q=num.hi!=0?num.hi:num.lo;*sign=(q>0?1:(q<0?-1:0))*(det.hi>0?1:-1);
 return isfinite(vertex[0])&&isfinite(vertex[1])&&isfinite(*residual);
}
GFS_PRED_HD int gfs_power_vertex_side2d(const double s[4][3],double vertex[2],double*residual,int*sign){
 GfsPowerDD d[4][3];
 for(int k=0;k<4;k++)for(int a=0;a<3;a++)d[k][a]=gfs_pdd_value(s[k][a]);
 return gfs_power_vertex_side_dd2d(d,vertex,residual,sign);
}
#undef GFS_PRED_HD
#ifdef GFS_PRED_PRECISE_SCOPE
#pragma float_control(pop)
#undef GFS_PRED_PRECISE_SCOPE
#endif
#endif
