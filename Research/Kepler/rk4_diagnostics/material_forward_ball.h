#ifndef MATERIAL_FORWARD_BALL_H
#define MATERIAL_FORWARD_BALL_H
#include <quadmath.h>
/* Absolute forward errors through the ACTUAL scalar operation graph.
 * Inputs are exact unless a radius is explicitly supplied. Each radius
 * operation is rounded outward. No observed-error multiplier is used.
 * Elementary pow/tanh/sqrt are an explicit <=8ulp engineering model,
 * independently checked by the original-basis binary128 oracle below.
 * This is not a formal proof of any CPU/GPU libm implementation. */
namespace material_forward {
template<class T> struct FP;
template<> struct FP<double>{
    static double abs(double x){return fabs(x);}static double up(double x){return nextafter(x,INFINITY);}static double down(double x){return nextafter(x,-INFINITY);}
    static double eps(){return DBL_EPSILON;}static double tiny(){return std::numeric_limits<double>::denorm_min();}
    static double power(double x,double p){return pow(x,p);}static double root(double x){return sqrt(x);}static double tanh_(double x){return tanh(x);}static double fma_(double a,double b,double c){return fma(a,b,c);}
};
template<> struct FP<long double>{
    static long double abs(long double x){return fabsl(x);}static long double up(long double x){return nextafterl(x,INFINITY);}static long double down(long double x){return nextafterl(x,-INFINITY);}
    static long double eps(){return LDBL_EPSILON;}static long double tiny(){return std::numeric_limits<long double>::denorm_min();}
    static long double power(long double x,long double p){return powl(x,p);}static long double root(long double x){return sqrtl(x);}static long double tanh_(long double x){return tanhl(x);}static long double fma_(long double a,long double b,long double c){return fmal(a,b,c);}
};
template<> struct FP<__float128>{
    using T=__float128;static T abs(T x){return fabsq(x);}static T up(T x){return nextafterq(x,__builtin_huge_valq());}static T down(T x){return nextafterq(x,-__builtin_huge_valq());}
    static T eps(){return ldexpq(1,-112);}static T tiny(){return ldexpq(1,-16494);}
    static T power(T x,T p){return powq(x,p);}static T root(T x){return sqrtq(x);}static T tanh_(T x){return tanhq(x);}static T fma_(T a,T b,T c){return fmaq(a,b,c);}
};
template<class T> static T plus_up(T a,T b){if(a==0)return b;if(b==0)return a;return FP<T>::up(a+b);}
template<class T> static T times_up(T a,T b){if(a==0 || b==0)return 0;return FP<T>::up(a*b);}
template<class T> static T round_error(T z){
    T u=FP<T>::eps()/2,d=FP<T>::down(T(1)-u);
    return plus_up(FP<T>::up(times_up(FP<T>::abs(z),u)/d),FP<T>::tiny());
}
template<class T> static T libm_error(T z){T hi=z,lo=z;for(int k=0;k<8;k++){hi=FP<T>::up(hi);lo=FP<T>::down(lo);}return std::max(FP<T>::abs(hi-z),FP<T>::abs(z-lo));}
template<class T> struct B {T n,e;B(T value=0,T error=0):n(value),e(error){need(error>=0,"nonnegative scalar forward radius");}};
template<class T> static B<T> operator+(B<T> a,B<T> b){T z=a.n+b.n;if(a.n==0 && a.e==0)return b;if(b.n==0 && b.e==0)return a;return {z,plus_up(plus_up(a.e,b.e),round_error(z))};}
template<class T> static B<T> operator-(B<T> a){return {-a.n,a.e};}
template<class T> static B<T> operator-(B<T> a,B<T> b){return a+(-b);}
template<class T> static B<T> operator*(B<T> a,B<T> b){if((a.n==0 && a.e==0)||(b.n==0 && b.e==0))return B<T>();T z=a.n*b.n,e=plus_up(times_up(FP<T>::abs(a.n),b.e),times_up(FP<T>::abs(b.n),a.e));e=plus_up(e,times_up(a.e,b.e));return {z,plus_up(e,round_error(z))};}
template<class T> static B<T> operator/(B<T> a,B<T> b){
    T z=a.n/b.n,lower=FP<T>::down(FP<T>::abs(b.n)-b.e);need(lower>0,"forward quotient denominator separated");T local=round_error(z),qmax=plus_up(FP<T>::abs(z),local),e=plus_up(a.e,times_up(qmax,b.e));e=FP<T>::up(e/lower);return {z,plus_up(e,local)};
}
template<class T> static B<T> fused(B<T> a,B<T> b,B<T> c){T z=FP<T>::fma_(a.n,b.n,c.n),e=plus_up(times_up(FP<T>::abs(a.n),b.e),times_up(FP<T>::abs(b.n),a.e));e=plus_up(e,times_up(a.e,b.e));e=plus_up(e,c.e);return {z,plus_up(e,round_error(z))};}
template<class T> static B<T> root(B<T> a){
    T low=FP<T>::down(a.n-a.e),high=FP<T>::up(a.n+a.e);need(low>0,"positive sqrt argument in forward model");T n=FP<T>::root(a.n),l=FP<T>::root(low),h=FP<T>::root(high);l=FP<T>::down(l-libm_error(l));h=FP<T>::up(h+libm_error(h));return {n,FP<T>::up(std::max(FP<T>::abs(n-l),FP<T>::abs(h-n)))};
}
template<class T> static B<T> power(B<T> a,T exponent){
    T low=FP<T>::down(a.n-a.e),high=FP<T>::up(a.n+a.e);need(low>0,"positive pow argument in forward model");T n=FP<T>::power(a.n,exponent),x=FP<T>::power(low,exponent),y=FP<T>::power(high,exponent);T l=std::min(FP<T>::down(x-libm_error(x)),FP<T>::down(y-libm_error(y))),h=std::max(FP<T>::up(x+libm_error(x)),FP<T>::up(y+libm_error(y)));return {n,FP<T>::up(std::max(FP<T>::abs(n-l),FP<T>::abs(h-n)))};
}
template<class T> static B<T> tanh_(B<T> a){T n=FP<T>::tanh_(a.n),l=FP<T>::tanh_(FP<T>::down(a.n-a.e)),h=FP<T>::tanh_(FP<T>::up(a.n+a.e));l=FP<T>::down(l-libm_error(l));h=FP<T>::up(h+libm_error(h));return {n,FP<T>::up(std::max(FP<T>::abs(n-l),FP<T>::abs(h-n)))};}
template<class T> struct D {B<T> v,x,y;D(B<T> value=B<T>(),B<T> dx=B<T>(),B<T> dy=B<T>()):v(value),x(dx),y(dy){}};
template<class T> static D<T> operator+(D<T>a,D<T>b){return {a.v+b.v,a.x+b.x,a.y+b.y};}
template<class T> static D<T> operator-(D<T>a){return {-a.v,-a.x,-a.y};}
template<class T> static D<T> operator-(D<T>a,D<T>b){return a+(-b);}
template<class T> static D<T> operator*(D<T>a,D<T>b){return {a.v*b.v,a.x*b.v+a.v*b.x,a.y*b.v+a.v*b.y};}
template<class T> static D<T> operator/(D<T>a,D<T>b){B<T> z=a.v/b.v;return {z,(a.x-z*b.x)/b.v,(a.y-z*b.y)/b.v};}
template<class T> static D<T> fused(D<T>a,D<T>b,D<T>c){return {fused(a.v,b.v,c.v),a.x*b.v+a.v*b.x+c.x,a.y*b.v+a.v*b.y+c.y};}
template<class T> static D<T> power(D<T>a,T p){B<T> n=power(a.v,p),factor=B<T>(p)*n/a.v;return {n,factor*a.x,factor*a.y};}
template<class T> struct Point {B<T> position[2],velocity[2],F[4],V[4],J,grad_v[4];};
template<class T> static void finish(Point<T> &p){
    p.J=p.F[0]*p.F[3]-p.F[1]*p.F[2];need(p.J.n-p.J.e>0,"forward map Jacobian positive");B<T> inv[4]={p.F[3]/p.J,-p.F[1]/p.J,-p.F[2]/p.J,p.F[0]/p.J};
    for(int k=0;k<2;k++){p.grad_v[2*k]=p.V[2*k]*inv[0]+p.V[2*k+1]*inv[2];p.grad_v[2*k+1]=p.V[2*k]*inv[1]+p.V[2*k+1]*inv[3];}
}
template<class T> static T cast_radius(T x){T a=FP<T>::up(x)-x,b=x-FP<T>::down(x);return FP<T>::up(std::max(a,b)/2);}
}
#endif
