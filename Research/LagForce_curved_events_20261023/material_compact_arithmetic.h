#ifndef MATERIAL_COMPACT_ARITHMETIC_H
#define MATERIAL_COMPACT_ARITHMETIC_H
// Same two-double arithmetic as the user's original CUDA mass operator,
// namespaced and available on host/device for tiny geometry differences.
#include <cmath>
#ifdef __CUDACC__
#define COMPACT_ARITH_HD __host__ __device__
#else
#define COMPACT_ARITH_HD
#endif
namespace compact_arithmetic {
struct Pair {double hi,lo;COMPACT_ARITH_HD Pair(double h=0,double l=0):hi(h),lo(l){}};
COMPACT_ARITH_HD static inline Pair operator+(Pair a,Pair b){double s=a.hi+b.hi,z=s-a.hi,e=(a.hi-(s-z))+(b.hi-z);e+=a.lo+b.lo;double h=s+e;return {h,e-(h-s)};}
COMPACT_ARITH_HD static inline Pair operator-(Pair a){return {-a.hi,-a.lo};}
COMPACT_ARITH_HD static inline Pair operator-(Pair a,Pair b){return a+(-b);}
COMPACT_ARITH_HD static inline Pair operator*(Pair a,Pair b){double p=a.hi*b.hi,e=fma(a.hi,b.hi,-p);e+=a.hi*b.lo+a.lo*b.hi;e+=a.lo*b.lo;double h=p+e;return {h,e-(h-p)};}
COMPACT_ARITH_HD static inline Pair operator/(Pair a,Pair b){double q=a.hi/b.hi;Pair r=a-Pair(q)*b;return Pair(q)+Pair((r.hi+r.lo)/b.hi);}
COMPACT_ARITH_HD static inline double rounded(Pair a){return a.hi+a.lo;}
struct Jet {
    Pair v,x,y;
    COMPACT_ARITH_HD Jet(Pair a=Pair(),Pair b=Pair(),Pair c=Pair()):v(a),x(b),y(c){}
};
COMPACT_ARITH_HD static inline Jet operator+(Jet a,Jet b){return {a.v+b.v,a.x+b.x,a.y+b.y};}
COMPACT_ARITH_HD static inline Jet operator*(Jet a,Pair b){return {a.v*b,a.x*b,a.y*b};}
} // namespace compact_arithmetic
#undef COMPACT_ARITH_HD
#endif
