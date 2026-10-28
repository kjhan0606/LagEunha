#ifndef EXAM_MOMENT_GEOMETRY_INTERVAL_H
#define EXAM_MOMENT_GEOMETRY_INTERVAL_H
#include <math.h>
#include <stdint.h>
#ifdef __CUDACC__
#define MG2_HD __host__ __device__
#else
#define MG2_HD
#endif
/* Requires IEEE binary64 and no fast-math. Every elementary operation is
 * rounded separately and expanded outward. NaN/overflow returns whole line.
 * This is a sign filter; a zero-containing result is never assigned a sign. */
struct MG2Interval{double lo,hi;};
MG2_HD inline MG2Interval mg2_whole(){return {-INFINITY,INFINITY};}
MG2_HD inline MG2Interval mg2_point(double x){return {x,x};}
MG2_HD inline MG2Interval mg2_add(MG2Interval a,MG2Interval b){
 volatile double lo=a.lo+b.lo,hi=a.hi+b.hi;
 if(!isfinite(lo)||!isfinite(hi))return mg2_whole();
 return {::nextafter(lo,-HUGE_VAL),::nextafter(hi,HUGE_VAL)};
}
MG2_HD inline MG2Interval mg2_neg(MG2Interval a){return {-a.hi,-a.lo};}
MG2_HD inline MG2Interval mg2_sub(MG2Interval a,MG2Interval b){return mg2_add(a,mg2_neg(b));}
MG2_HD inline MG2Interval mg2_mul(MG2Interval a,MG2Interval b){
 volatile double p0=a.lo*b.lo,p1=a.lo*b.hi,p2=a.hi*b.lo,p3=a.hi*b.hi;
 if(!isfinite(p0)||!isfinite(p1)||!isfinite(p2)||!isfinite(p3))return mg2_whole();
 double lo=fmin(fmin(p0,p1),fmin(p2,p3)),hi=fmax(fmax(p0,p1),fmax(p2,p3));
 return {::nextafter(lo,-HUGE_VAL),::nextafter(hi,HUGE_VAL)};
}
MG2_HD inline int mg2_sign(MG2Interval a){return a.lo>0.?1:(a.hi<0.?-1:0);}
struct MG2PlaneInterval{MG2Interval x,y,b;};
MG2_HD inline MG2PlaneInterval mg2_plane(int64_t i,int64_t j,const double *xy,
 const double *W,const double *bounds){
 MG2Interval zero=mg2_point(0.),one=mg2_point(1.),minus=mg2_point(-1.);
 if(j==-1)return {zero,minus,mg2_sub(mg2_point(xy[2*i+1]),mg2_point(bounds[2]))};
 if(j==-2)return {one,zero,mg2_sub(mg2_point(bounds[1]),mg2_point(xy[2*i]))};
 if(j==-3)return {zero,one,mg2_sub(mg2_point(bounds[3]),mg2_point(xy[2*i+1]))};
 if(j==-4)return {minus,zero,mg2_sub(mg2_point(xy[2*i]),mg2_point(bounds[0]))};
 MG2Interval dx=mg2_sub(mg2_point(xy[2*j]),mg2_point(xy[2*i]));
 MG2Interval dy=mg2_sub(mg2_point(xy[2*j+1]),mg2_point(xy[2*i+1]));
 MG2Interval b=mg2_mul(mg2_point(.5),mg2_add(mg2_add(mg2_mul(dx,dx),mg2_mul(dy,dy)),mg2_sub(mg2_point(W[i]),mg2_point(W[j]))));
 return {dx,dy,b};
}
/* Sign of n3 dot intersection(p1,p2)-b3, determined from support planes.
 * Zero means unresolved, including exact boundary equality and singular D. */
MG2_HD inline int mg2_support_sign(MG2PlaneInterval p,MG2PlaneInterval q,MG2PlaneInterval r){
 MG2Interval D=mg2_sub(mg2_mul(p.x,q.y),mg2_mul(p.y,q.x));
 MG2Interval X=mg2_sub(mg2_mul(p.b,q.y),mg2_mul(p.y,q.b));
 MG2Interval Y=mg2_sub(mg2_mul(p.x,q.b),mg2_mul(p.b,q.x));
 MG2Interval N=mg2_sub(mg2_add(mg2_mul(r.x,X),mg2_mul(r.y,Y)),mg2_mul(r.b,D));
 int ds=mg2_sign(D),ns=mg2_sign(N);return ds&&ns?ds*ns:0;
}
#endif
