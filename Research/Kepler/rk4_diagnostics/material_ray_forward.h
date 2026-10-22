#ifndef MATERIAL_RAY_FORWARD_H
#define MATERIAL_RAY_FORWARD_H
#include "material_forward_ball.h"
/* Range of the actually rounded scalar and its absolute error relative to
 * the exact midpoint representative. No operation-count multiplier. The
 * nominal range encloses every radius in the unchanged accepted panel.
 * All libm endpoints retain the explicit eight-ulp engineering assumption. */
namespace material_ray_forward {
using R=long double;using I=curved_bounds::I;namespace mf=material_forward;
struct Range {I n;R e;Range(R x=0):n(x),e(0){}Range(I x,R error=0):n(x),e(error){}};
static bool zero(Range a){return curved_bounds::exact_zero(a.n)&&a.e==0;}
static R mag(I a){return curved_bounds::absmax(a);}
static Range operator+(Range a,Range b){if(zero(a))return b;if(zero(b))return a;I n=a.n+b.n;return {n,mf::plus_up(mf::plus_up(a.e,b.e),mf::round_error(mag(n)))};}
static Range operator-(Range a){return {-a.n,a.e};}
static Range operator-(Range a,Range b){return a+(-b);}
static Range operator*(Range a,Range b){if(zero(a)||zero(b))return {};I n=a.n*b.n;R e=mf::plus_up(mf::times_up(mag(a.n),b.e),mf::times_up(mag(b.n),a.e));e=mf::plus_up(e,mf::times_up(a.e,b.e));return {n,mf::plus_up(e,mf::round_error(mag(n)))};}
static Range operator/(Range a,Range b){I n=a.n/b.n;R lower=std::min(fabsl(b.n.lo),fabsl(b.n.hi));lower=curved_bounds::down(lower-b.e);need((b.n.lo>0||b.n.hi<0)&&lower>0,"ray forward denominator separated");R local=mf::round_error(mag(n)),e=mf::plus_up(a.e,mf::times_up(mf::plus_up(mag(n),local),b.e));return {n,mf::plus_up(curved_bounds::up(e/lower),local)};}
static Range power(Range a,R exponent){
    I input=curved_bounds::widen(a.n,a.e);need(input.lo>0,"ray power argument positive");
    R x=powl(a.n.lo,exponent),y=powl(a.n.hi,exponent),l=powl(input.lo,exponent),h=powl(input.hi,exponent);
    I n(std::min(curved_bounds::down(x-mf::libm_error(x)),curved_bounds::down(y-mf::libm_error(y))),std::max(curved_bounds::up(x+mf::libm_error(x)),curved_bounds::up(y+mf::libm_error(y))));
    R derivative=fabsl(exponent)*std::max(powl(input.lo,exponent-1),powl(input.hi,exponent-1));derivative=curved_bounds::up(derivative+mf::libm_error(derivative));
    R round=std::max({mf::libm_error(x),mf::libm_error(y),mf::libm_error(l),mf::libm_error(h)});
    return {n,mf::plus_up(round,mf::times_up(derivative,a.e))};
}
template<class S> struct Poly {int degree;std::array<S,7> c;explicit Poly(int d=0):degree(d),c{} {}};
template<class S> static Poly<S> constant(S a){Poly<S> p;p.c[0]=a;return p;}
template<class S> static Poly<S> add(const Poly<S>&a,const Poly<S>&b){Poly<S> p(std::max(a.degree,b.degree));for(int i=0;i<=p.degree;i++)p.c[i]=a.c[i]+b.c[i];return p;}
template<class S> static Poly<S> scale(const Poly<S>&a,S b){Poly<S> p(a.degree);for(int i=0;i<=a.degree;i++)p.c[i]=a.c[i]*b;return p;}
template<class S> static Poly<S> sub(const Poly<S>&a,const Poly<S>&b){return add(a,scale(b,S(-1)));}
template<class S> static Poly<S> mul(const Poly<S>&a,const Poly<S>&b){need(a.degree+b.degree<=6,"unchanged degree-six ray graph");Poly<S> p(a.degree+b.degree);for(int i=0;i<=a.degree;i++)for(int j=0;j<=b.degree;j++)p.c[i+j]=p.c[i+j]+a.c[i]*b.c[j];return p;}
template<class S> static Poly<S> power_poly(Poly<S>a,int n){auto p=constant(S(1));for(int j=0;j<n;j++)p=mul(p,a);return p;}
template<class S> static S value(const Poly<S>&p,S u){S v=p.c[p.degree];for(int i=p.degree-1;i>=0;i--)v=v*u+p.c[i];return v;}
template<class S> struct Ray {Poly<S> D,position[2],bary[3];};
template<class S> static Ray<S> ray(const curved_branches::FixedMap&m,S rad,int chart,R eps){
    Ray<S> out;Poly<S> D(2);D.c[0]=S(1);D.c[2]=S(1);out.D=D;auto D2=mul(D,D);Poly<S> X(2),Y(2);X.c[0]=rad*S(chart);X.c[2]=rad*S(-chart);Y.c[1]=rad*S(2*chart);
    auto x=scale(sub(X,scale(D,S(m.anchor[0]))),S(1)/S(m.h)),y=scale(sub(Y,scale(D,S(m.anchor[1]))),S(1)/S(m.h));
    S z=rad*rad+S(eps)*S(eps),omega=power(z,R(-.75L));Poly<S> initial[2]={scale(mul(scale(Y,S(-1)),D2),omega),scale(mul(X,D2),omega)};
    for(int k=0;k<2;k++){auto p=constant(S(0));for(int a=0;a<10;a++)if(m.position[k][a]!=0){auto e=curved_compiled::exponents(a);p=add(p,scale(mul(mul(power_poly(x,e[0]),power_poly(y,e[1])),power_poly(D,3-e[0]-e[1])),S(m.position[k][a])));}for(int j=0;j<2;j++)if(m.B[k][j]!=0)p=add(p,scale(initial[j],S(m.B[k][j])));out.position[k]=p;}
    for(int k=0;k<3;k++)out.bary[k]=sub(add(scale(X,S(m.source_n[k][0])),scale(Y,S(m.source_n[k][1]))),scale(D,S(m.source_c[k])));return out;
}
template<class S> static Poly<S> polynomial(const Ray<S>&r,const curved_branches::Curve&c,const std::vector<curved_partition::Site>&site,const std::vector<R>&weight){
    if(c.source>=0)return r.bary[c.source];const auto&a=site[c.i],&b=site[c.j];
    /* These rounded power-plane scalars are EXACT inputs to this graph.
     * Their construction error is already in conversion(). Charging it
     * again here would duplicate the same dependency path. */
    R nx=b.x-a.x,ny=b.y-a.y,constant_=(b.x*b.x+b.y*b.y-a.x*a.x-a.y*a.y+weight[c.i]-weight[c.j])/2;
    return sub(add(scale(r.position[0],S(nx)),scale(r.position[1],S(ny))),scale(power_poly(r.D,3),S(constant_)));
}
}
#endif
