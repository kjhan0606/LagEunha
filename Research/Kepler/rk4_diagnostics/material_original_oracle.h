#ifndef MATERIAL_ORIGINAL_ORACLE_H
#define MATERIAL_ORIGINAL_ORACLE_H
#include "material_forward_ball.h"
namespace material_forward {
using Q=__float128;
/* Independent ORIGINAL basis, not the rounded compiled cubic. Stored
 * doubles are exact input constants. Analytic differentiation uses
 * l0=1-l1-l2, including any stored-gradient compatibility discrepancy. */
static Point<Q> original_at(const curved::Record &r,B<Q> px,B<Q> py,double eps){
    using S=B<Q>;using J=D<Q>;J X(px,S(1),S(0)),Y(py,S(0),S(1));
    S ax{Q(r.triangle.geometry.x[0][0])},ay{Q(r.triangle.geometry.x[0][1])},ep{Q(eps)};J dx=X-J(ax),dy=Y-J(ay),l[3],N[12];
    l[1]=J(S(Q(r.triangle.geometry.grad[1][0])))*dx+J(S(Q(r.triangle.geometry.grad[1][1])))*dy;
    l[2]=J(S(Q(r.triangle.geometry.grad[2][0])))*dx+J(S(Q(r.triangle.geometry.grad[2][1])))*dy;l[0]=J(S(1))-l[1]-l[2];
    for(int a=0;a<10;a++){int id[3];gfs_material_p3_index(a,id);N[a]=J(S(1));for(int k=0;k<3;k++)for(int j=0;j<id[k];j++)N[a]=N[a]*(J(S(3))*l[k]-J(S(j)))/J(S(id[k]-j));}
    S za=ax*ax+ay*ay+ep*ep,oa=power(za,Q(-.75L));J L=(J(S(2)*ax)*dx+J(S(2)*ay)*dy)/J(za),R=(dx*dx+dy*dy)/J(za),L2=L*L;
    J A=J(S(1))-J(S(.75L))*L-J(S(.75L))*R+J(S(.65625L))*L2+J(S(1.3125L))*L*R-J(S(.6015625L))*L2*L;
    J C=J(S(1))-J(S(.75L))*L-J(S(.75L))*R+J(S(.65625L))*L2;
    J omega=power(X*X+Y*Y+J(ep*ep),Q(-.75L));J q[2]={(-Y)*omega-J(oa)*(-J(ay)*A-dy*C),X*omega-J(oa)*(J(ax)*A+dx*C)};
    for(int k=0;k<2;k++){J v=q[k]-J(S(Q(r.triangle.remainder[0][k])));for(int a=1;a<10;a++)v=v-N[a]*J(S(Q(r.triangle.remainder[a][k]))-S(Q(r.triangle.remainder[0][k])));N[10+k]=v/J(S(Q(r.scale[k])));}
    Point<Q> p={};for(int k=0;k<2;k++){J x(S(Q(r.origin[k]))),v;for(int a=0;a<12;a++){x=x+J(S(Q(r.map.x[k][a])))*N[a];v=v+J(S(Q(r.map.v[k][a])))*N[a];}p.position[k]=x.v;p.velocity[k]=v.v;p.F[2*k]=x.x;p.F[2*k+1]=x.y;p.V[2*k]=v.x;p.V[2*k+1]=v.y;}finish(p);return p;
}
static Point<Q> original(const curved::Record&r,double px,double py,double eps){return original_at(r,B<Q>(Q(px)),B<Q>(Q(py)),eps);}
/* Stable local cubic + B v0, evaluated by explicit FMA Horner. The
 * representative uses midpoint coefficients; conversion=true includes
 * every original->compiled coefficient interval. Neither modifies state. */
static Point<long double> stable(const curved_compiled::Map &m,double px,double py,double eps,bool conversion){
    using T=long double;using S=B<T>;using J=D<T>;J X(S(px),S(1),S(0)),Y(S(py),S(0),S(1)),x=(X-J(S(m.anchor[0])))/J(S(m.h)),y=(Y-J(S(m.anchor[1])))/J(S(m.h));J omega=power(X*X+Y*Y+J(S(eps)*S(eps)),T(-.75L)),v0[2]={-Y*omega,X*omega};
    auto coefficient=[&](curved_bounds::I z){T n=z.lo+(z.hi-z.lo)/2,e=conversion?FP<T>::up(std::max(n-z.lo,z.hi-n)):0;return J(S(n,e));};
    auto cubic=[&](const curved_compiled::C &c){
        auto at=[&](int a,int b){return coefficient(c.a[curved_compiled::index(a,b)]);};
        J a0=fused(y,fused(y,fused(y,at(0,3),at(0,2)),at(0,1)),at(0,0));
        J a1=fused(y,fused(y,at(1,2),at(1,1)),at(1,0)),a2=fused(y,at(2,1),at(2,0));
        return fused(x,fused(x,fused(x,at(3,0),a2),a1),a0);
    };
    Point<T> p={};for(int k=0;k<2;k++){J z=cubic(m.position[k]),v=cubic(m.velocity[k]);for(int j=0;j<2;j++){z=z+coefficient(m.position_B[k][j])*v0[j];v=v+coefficient(m.velocity_B[k][j])*v0[j];}p.position[k]=z.v;p.velocity[k]=v.v;p.F[2*k]=z.x;p.F[2*k+1]=z.y;p.V[2*k]=v.x;p.V[2*k+1]=v.y;}finish(p);return p;
}
/* Independent monomial evaluation of the EXACT midpoint representative,
 * used only to check FMA-Horner evaluation error, not as physical reference. */
static Point<Q> representative(const curved_compiled::Map &m,double px,double py,double eps){
    using S=B<Q>;using J=D<Q>;J X(S(Q(px)),S(1),S(0)),Y(S(Q(py)),S(0),S(1)),x=(X-J(S(Q(m.anchor[0]))))/J(S(Q(m.h))),y=(Y-J(S(Q(m.anchor[1]))))/J(S(Q(m.h)));
    auto c=[&](curved_bounds::I a){return J(S(Q(a.lo+(a.hi-a.lo)/2)));};
    auto ipow=[&](J a,int n){J z(S(1));for(int j=0;j<n;j++)z=z*a;return z;};
    J omega=power(X*X+Y*Y+J(S(Q(eps))*S(Q(eps))),Q(-.75L)),v0[2]={-Y*omega,X*omega};Point<Q> p={};
    for(int k=0;k<2;k++){J z,v;for(int a=0;a<10;a++){auto e=curved_compiled::exponents(a);J n=ipow(x,e[0])*ipow(y,e[1]);z=z+c(m.position[k].a[a])*n;v=v+c(m.velocity[k].a[a])*n;}for(int j=0;j<2;j++){z=z+c(m.position_B[k][j])*v0[j];v=v+c(m.velocity_B[k][j])*v0[j];}p.position[k]=z.v;p.velocity[k]=v.v;p.F[2*k]=z.x;p.F[2*k+1]=z.y;p.V[2*k]=v.x;p.V[2*k+1]=v.y;}finish(p);return p;
}
}
#endif
