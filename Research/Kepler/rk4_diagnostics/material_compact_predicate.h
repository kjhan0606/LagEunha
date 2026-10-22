#ifndef MATERIAL_COMPACT_PREDICATE_H
#define MATERIAL_COMPACT_PREDICATE_H
// Exact sign of a binary64 orientation, by error-free sums/FMA products.
// Bounded metadata predicate; no tolerance merge or approximate zero test.
#include "material_compact_arithmetic.h"
#include <array>
#include <stdexcept>
namespace compact_predicate {
inline void coordinate(double x){
    if(!std::isfinite(x)||fabs(x)>16||(x!=0&&fabs(x)<ldexp(1.,-128)))throw std::runtime_error("orientation input outside normal-product engineering range");
}
struct Expansion {
    std::array<double,32>value={};size_t size=0;
    void add(double x){
        if(x==0)return;std::array<double,32>out={};size_t n=0;double q=x;
        for(size_t j=0;j<size;j++){
            double a=q,b=value[j],s=a+b,z=s-a,error=(a-(s-z))+(b-z);
            if(error!=0){if(n>=out.size())throw std::runtime_error("bounded exact predicate expansion");out[n++]=error;}q=s;
        }
        if(q!=0){if(n>=out.size())throw std::runtime_error("bounded exact predicate expansion");out[n++]=q;}
        value=out;size=n;
    }
    void product(compact_arithmetic::Pair a,compact_arithmetic::Pair b,double sign){
        for(double x:{a.lo,a.hi})for(double y:{b.lo,b.hi}){
            double p=x*y,e=fma(x,y,-p);add(sign*e);add(sign*p);
        }
    }
    int sign()const{return size?(value[size-1]>0?1:-1):0;}
};
inline int orientation(const double*A,const double*B,const double*P){
    for(int k=0;k<2;k++){coordinate(A[k]);coordinate(B[k]);coordinate(P[k]);}
    using compact_arithmetic::Pair;
    Pair ax=Pair(B[0])-Pair(A[0]),ay=Pair(B[1])-Pair(A[1]);
    Pair bx=Pair(P[0])-Pair(A[0]),by=Pair(P[1])-Pair(A[1]);
    Expansion result;result.product(ax,by,1);result.product(ay,bx,-1);return result.sign();
}
} // namespace compact_predicate
#endif
