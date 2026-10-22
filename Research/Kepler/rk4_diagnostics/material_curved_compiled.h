#ifndef MATERIAL_CURVED_COMPILED_H
#define MATERIAL_CURVED_COMPILED_H
#include "material_curved_bounds.h"
namespace curved_compiled {
using R=long double;using I=curved_bounds::I;
static int index(int x,int y){need(x>=0 && y>=0 && x+y<=3,"complete bivariate cubic index");int k=0;for(int n=0;n<=3;n++)for(int a=0;a<=n;a++,k++)if(a==x && n-a==y)return k;need(false,"cubic monomial index");return -1;}
static std::array<int,2> exponents(int k){int j=0;for(int n=0;n<=3;n++)for(int a=0;a<=n;a++,j++)if(j==k)return {a,n-a};need(false,"cubic monomial exponent");return {};}
static bool zero(I a){return a.lo==0 && a.hi==0;}
struct C {int degree;std::array<I,10> a;explicit C(int n=0):degree(n),a{} {}};
static C constant(I x){C c;c.a[0]=x;return c;}
static C add(const C &a,const C &b){C c(std::max(a.degree,b.degree));for(int i=0;i<10;i++){if(zero(a.a[i]))c.a[i]=b.a[i];else if(zero(b.a[i]))c.a[i]=a.a[i];else c.a[i]=a.a[i]+b.a[i];}return c;}
static C scale(const C &a,I x){if(zero(x))return C();C c(a.degree);for(int i=0;i<10;i++)if(!zero(a.a[i]))c.a[i]=a.a[i]*x;return c;}
static C sub(const C &a,const C &b){return add(a,scale(b,I(-1)));}
static C mul(const C &a,const C &b){need(a.degree+b.degree<=3,"cubic static conversion, no truncation");C c(a.degree+b.degree);for(int i=0;i<10;i++)if(!zero(a.a[i]))for(int j=0;j<10;j++)if(!zero(b.a[j])){auto p=exponents(i),q=exponents(j);int k=index(p[0]+q[0],p[1]+q[1]);I value=a.a[i]*b.a[j];c.a[k]=zero(c.a[k])?value:c.a[k]+value;}return c;}
struct Map {R h,anchor[2];C position[2],velocity[2];I position_B[2][2],velocity_B[2][2];};
static Map compile(const curved::Record &r,R eps){
    Map out={};out.h=0;for(int a=0;a<3;a++)for(int k=0;k<2;k++)out.h=std::max(out.h,fabsl((R)r.triangle.geometry.x[a][k]-r.triangle.geometry.x[0][k]));need(out.h>0,"positive static polynomial scale");for(int k=0;k<2;k++)out.anchor[k]=r.triangle.geometry.x[0][k];C dx(1),dy(1);dx.a[index(1,0)]=I(out.h);dy.a[index(0,1)]=I(out.h);C bary[3],N[10];
    bary[1]=add(scale(dx,I(r.triangle.geometry.grad[1][0])),scale(dy,I(r.triangle.geometry.grad[1][1])));bary[2]=add(scale(dx,I(r.triangle.geometry.grad[2][0])),scale(dy,I(r.triangle.geometry.grad[2][1])));bary[0]=sub(sub(constant(I(1)),bary[1]),bary[2]);
    for(int a=0;a<10;a++){int id[3];gfs_material_p3_index(a,id);N[a]=constant(I(1));for(int k=0;k<3;k++)for(int j=0;j<id[k];j++)N[a]=mul(N[a],scale(sub(scale(bary[k],I(3)),constant(I(j))),I(1)/I(id[k]-j)));}
    auto interpolate=[&](const std::array<I,10> &value){C p=add(constant(value[0]),add(scale(bary[1],value[1]-value[0]),scale(bary[2],value[2]-value[0])));for(int a=3;a<10;a++){int id[3];gfs_material_p3_index(a,id);I linear(0);for(int k=0;k<3;k++)linear=linear+I(id[k])*value[k]/I(3);p=add(p,scale(N[a],value[a]-linear));}return p;};
    I ax(out.anchor[0]),ay(out.anchor[1]),z=curved_bounds::square(ax)+curved_bounds::square(ay)+curved_bounds::square(I(eps)),w=curved_bounds::root(curved_bounds::root(z)),omega=I(1)/(w*w*w);C L=scale(add(scale(dx,I(2)*ax),scale(dy,I(2)*ay)),I(1)/z),Q=scale(add(mul(dx,dx),mul(dy,dy)),I(1)/z),L2=mul(L,L);
    C A=add(add(add(constant(I(1)),scale(L,I(-.75L))),scale(Q,I(-.75L))),add(add(scale(L2,I(.65625L)),scale(mul(L,Q),I(1.3125L))),scale(mul(L2,L),I(-.6015625L))));C B=add(add(constant(I(1)),scale(L,I(-.75L))),add(scale(Q,I(-.75L)),scale(L2,I(.65625L))));C initial_poly[2]={scale(add(scale(A,-ay),mul(scale(dy,I(-1)),B)),omega),scale(add(scale(A,ax),mul(dx,B)),omega)};
    for(int j=0;j<2;j++){std::array<I,10> values;for(int a=0;a<10;a++)values[a]=I(r.triangle.remainder[a][j]);initial_poly[j]=add(initial_poly[j],interpolate(values));}
    for(int k=0;k<2;k++){
        std::array<I,10> x,v;for(int a=0;a<10;a++){x[a]=I(r.map.x[k][a]);v[a]=I(r.map.v[k][a]);}out.position[k]=add(constant(I(r.origin[k])),interpolate(x));out.velocity[k]=interpolate(v);
        for(int j=0;j<2;j++){out.position_B[k][j]=r.map.x[k][10+j]==0?I(0):I(r.map.x[k][10+j])/I(r.scale[j]);out.velocity_B[k][j]=r.map.v[k][10+j]==0?I(0):I(r.map.v[k][10+j])/I(r.scale[j]);out.position[k]=sub(out.position[k],scale(initial_poly[j],out.position_B[k][j]));out.velocity[k]=sub(out.velocity[k],scale(initial_poly[j],out.velocity_B[k][j]));}
    }return out;
}
static I power(I x,int n){I p(1);for(int j=0;j<n;j++)p=p*x;return p;}
static I value(const C &p,I x,I y){I z(0);for(int a=0;a<10;a++)if(!zero(p.a[a])){auto n=exponents(a);z=z+p.a[a]*power(x,n[0])*power(y,n[1]);}return z;}
static I gradient(const C &p,I x,I y,int d){I z(0);for(int a=0;a<10;a++)if(!zero(p.a[a])){auto n=exponents(a);if(n[d]){int c=n[d];n[d]--;z=z+p.a[a]*I(c)*power(x,n[0])*power(y,n[1]);}}return z;}
}
#endif
