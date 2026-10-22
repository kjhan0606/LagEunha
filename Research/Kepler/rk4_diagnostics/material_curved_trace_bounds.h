#ifndef MATERIAL_CURVED_TRACE_BOUNDS_H
#define MATERIAL_CURVED_TRACE_BOUNDS_H
#include "material_curved_roots.h"
namespace curved_trace {
using namespace curved_roots;
struct EdgeForm {P polynomial[2],radial[2],X,Y,T3[2];};
static std::array<P,2> taylor(const P &X,const P &Y,I ax,I ay,R eps){
    I z=curved_bounds::square(ax)+curved_bounds::square(ay)+curved_bounds::square(I(eps)),w=curved_bounds::root(curved_bounds::root(z)),omega=I(1)/(w*w*w);P dx=sub(X,constant(ax)),dy=sub(Y,constant(ay)),L=scale(add(scale(dx,I(2)*ax),scale(dy,I(2)*ay)),I(1)/z),Q=scale(add(mul(dx,dx),mul(dy,dy)),I(1)/z),L2=mul(L,L);
    P A=add(add(add(constant(I(1)),scale(L,I(-.75L))),scale(Q,I(-.75L))),add(add(scale(L2,I(.65625L)),scale(mul(L,Q),I(1.3125L))),scale(mul(L2,L),I(-.6015625L))));
    P B=add(add(constant(I(1)),scale(L,I(-.75L))),add(scale(Q,I(-.75L)),scale(L2,I(.65625L))));return {scale(add(scale(A,-ay),mul(scale(dy,I(-1)),B)),omega),scale(add(scale(A,ax),mul(dx,B)),omega)};
}
static EdgeForm edge_form(const curved::Record &r,const std::array<double,2> &a,const std::array<double,2> &b,R eps,bool velocity=false){
    EdgeForm out;out.X=P(1);out.Y=P(1);out.X.c[0]=I(a[0]);out.X.c[1]=I(b[0])-I(a[0]);out.Y.c[0]=I(a[1]);out.Y.c[1]=I(b[1])-I(a[1]);I ax(r.triangle.geometry.x[0][0]),ay(r.triangle.geometry.x[0][1]);P dx=sub(out.X,constant(ax)),dy=sub(out.Y,constant(ay)),bary[3],N[10];
    bary[1]=add(scale(dx,I(r.triangle.geometry.grad[1][0])),scale(dy,I(r.triangle.geometry.grad[1][1])));bary[2]=add(scale(dx,I(r.triangle.geometry.grad[2][0])),scale(dy,I(r.triangle.geometry.grad[2][1])));bary[0]=sub(sub(constant(I(1)),bary[1]),bary[2]);
    for(int j=0;j<10;j++){N[j]=constant(I(1));int index[3];gfs_material_p3_index(j,index);for(int k=0;k<3;k++)for(int l=0;l<index[k];l++)N[j]=mul(N[j],scale(sub(scale(bary[k],I(3)),constant(I(l))),I(1)/I(index[k]-l)));}
    auto T=taylor(out.X,out.Y,ax,ay,eps);P vbase[2]={scale(out.Y,I(-1)),out.X};const auto &coefficient=velocity?r.map.v:r.map.x;
    for(int k=0;k<2;k++){
        out.polynomial[k]=constant(I(velocity?0:r.origin[k]));out.radial[k]=constant(I(0));out.T3[k]=T[k];for(int j=0;j<10;j++)out.polynomial[k]=add(out.polynomial[k],scale(N[j],I(coefficient[k][j])));
        for(int d=0;d<2;d++){I B=I(coefficient[k][10+d])/I(r.scale[d]);P initial_poly=T[d];for(int j=0;j<10;j++)initial_poly=add(initial_poly,scale(N[j],I(r.triangle.remainder[j][d])));out.polynomial[k]=sub(out.polynomial[k],scale(initial_poly,B));out.radial[k]=add(out.radial[k],scale(vbase[d],B));}
    }return out;
}
static I hull(const P &p,R a,R b){auto c=bernstein(p,a,b);I h=c[0];for(int j=1;j<=p.degree;j++){h.lo=std::min(h.lo,c[j].lo);h.hi=std::max(h.hi,c[j].hi);}return h;}
struct Gap {R position=0,velocity=0;};
static R difference_bound(const EdgeForm &a,const EdgeForm &b,R eps){
    R gap[2]={};for(int part=0;part<16;part++){
        R lo=part/16.L,hi=(part+1)/16.L;I X=hull(a.X,lo,hi),Y=hull(a.Y,lo,hi),z=curved_bounds::square(X)+curved_bounds::square(Y)+curved_bounds::square(I(eps)),w=curved_bounds::root(curved_bounds::root(z)),omega=I(1)/(w*w*w);
        for(int k=0;k<2;k++){P dpoly=sub(a.polynomial[k],b.polynomial[k]),dradial=sub(a.radial[k],b.radial[k]);I d=hull(dpoly,lo,hi)+omega*hull(dradial,lo,hi);gap[k]=std::max(gap[k],absmax(d));}
    }return curved_bounds::up(hypotl(gap[0],gap[1]));
}
}
#endif
