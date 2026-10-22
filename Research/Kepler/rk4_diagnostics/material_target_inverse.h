#ifndef MATERIAL_TARGET_INVERSE_H
#define MATERIAL_TARGET_INVERSE_H
#include "material_original_oracle.h"
/* The original basis with zero residual coefficients is a genuine cubic.
 * This helper certifies that exact stored function, in local coordinates.
 * No sample fit, trace average or coefficient replacement is performed. */
namespace material_target_inverse {
namespace mf=material_forward;using Q=mf::Q;using S=mf::B<Q>;using R=long double;
static Q mag(S a){return mf::plus_up(fabsq(a.n),a.e);}
static Q lower(S a){return mf::FP<Q>::down(a.n-a.e);}
static Q upper(S a){return mf::FP<Q>::up(a.n+a.e);}
static Q add(Q a,Q b){return mf::plus_up(a,b);}static Q mul(Q a,Q b){return mf::times_up(a,b);}
struct P {int degree=0;S c[7][7]={};};
static P constant(S x){P p;p.c[0][0]=x;return p;}
static P plus(const P&a,const P&b){P p;p.degree=std::max(a.degree,b.degree);for(int i=0;i<=p.degree;i++)for(int j=0;j+i<=p.degree;j++)p.c[i][j]=a.c[i][j]+b.c[i][j];return p;}
static P scale(const P&a,S b){P p;p.degree=a.degree;for(int i=0;i<=p.degree;i++)for(int j=0;j+i<=p.degree;j++)p.c[i][j]=a.c[i][j]*b;return p;}
static P minus(const P&a,const P&b){return plus(a,scale(b,S(-1)));}
static P product(const P&a,const P&b){P p;p.degree=a.degree+b.degree;need(p.degree<=6,"finite original polynomial degree");for(int i=0;i<=a.degree;i++)for(int j=0;j+i<=a.degree;j++)for(int k=0;k<=b.degree;k++)for(int l=0;l+k<=b.degree;l++)p.c[i+k][j+l]=p.c[i+k][j+l]+a.c[i][j]*b.c[k][l];return p;}
static P derivative(const P&a,int d){P p;p.degree=std::max(0,a.degree-1);for(int i=0;i<=a.degree;i++)for(int j=0;j+i<=a.degree;j++)if(d?j:i)p.c[i-(d==0)][j-(d==1)]=a.c[i][j]*S(d?j:i);return p;}
static Q ipow(Q x,int n){Q z=1;for(int k=0;k<n;k++)z=mul(z,x);return z;}
static Q norm(const P&p,Q radius){Q n=0;for(int i=0;i<=p.degree;i++)for(int j=0;j+i<=p.degree;j++)n=add(n,mul(mag(p.c[i][j]),ipow(radius,i+j)));return n;}
static S value(const P&p,S x,S y){S result;for(int i=p.degree;i>=0;i--){S row;for(int j=p.degree-i;j>=0;j--)row=row*y+p.c[i][j];result=result*x+row;}return result;}
struct PolynomialMap {P x[2],v[2],F[4],V[4],J;Q anchor[2];};
static PolynomialMap compile(const curved::Record&r,const Q*anchor){
    for(int k=0;k<2;k++)for(int j=10;j<12;j++)need(r.map.x[k][j]==0&&r.map.v[k][j]==0,"this first lens gate requires actual zero residual coefficients");
    PolynomialMap out;out.anchor[0]=anchor[0];out.anchor[1]=anchor[1];P X=constant(S(anchor[0])-S(Q(r.triangle.geometry.x[0][0]))),Y=constant(S(anchor[1])-S(Q(r.triangle.geometry.x[0][1])));X.degree=Y.degree=1;X.c[1][0]=S(1);Y.c[0][1]=S(1);P l[3],N[10];l[1]=plus(scale(X,S(Q(r.triangle.geometry.grad[1][0]))),scale(Y,S(Q(r.triangle.geometry.grad[1][1]))));l[2]=plus(scale(X,S(Q(r.triangle.geometry.grad[2][0]))),scale(Y,S(Q(r.triangle.geometry.grad[2][1]))));l[0]=minus(minus(constant(S(1)),l[1]),l[2]);
    for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);N[a]=constant(S(1));for(int k=0;k<3;k++)for(int j=0;j<index[k];j++)N[a]=product(N[a],scale(minus(scale(l[k],S(3)),constant(S(j))),S(1)/S(index[k]-j)));}
    for(int k=0;k<2;k++){out.x[k]=constant(S(Q(r.origin[k])));for(int a=0;a<10;a++){out.x[k]=plus(out.x[k],scale(N[a],S(Q(r.map.x[k][a]))));out.v[k]=plus(out.v[k],scale(N[a],S(Q(r.map.v[k][a]))));}for(int d=0;d<2;d++){out.F[2*k+d]=derivative(out.x[k],d);out.V[2*k+d]=derivative(out.v[k],d);}}
    out.J=minus(product(out.F[0],out.F[3]),product(out.F[1],out.F[2]));return out;
}
struct InverseCertificate {Q A[4]={},radius=0,contraction=0,invariance=0,global_injectivity=0,J_lower=0,J_upper=0,source_margin=0;};
static InverseCertificate certify(const curved::Record&r,const PolynomialMap&m,const Q*physical_center,const Q*complex_radius){
    InverseCertificate c;c.radius=Q(.01L);S F[4];for(int k=0;k<4;k++)F[k]=m.F[k].c[0][0];S det=F[0]*F[3]-F[1]*F[2];S A[4]={F[3]/det,-F[1]/det,-F[2]/det,F[0]/det};for(int k=0;k<4;k++)c.A[k]=A[k].n;
    for(int i=0;i<2;i++){P G=constant(S(0));for(int k=0;k<2;k++)G=plus(G,scale(m.x[k],S(-c.A[2*i+k])));G.c[i==0][i==1]=G.c[i==0][i==1]+S(1);Q q=add(norm(derivative(G,0),c.radius),norm(derivative(G,1),c.radius));c.contraction=std::max(c.contraction,q);S offset;Q physical=0;for(int k=0;k<2;k++){offset=offset+S(c.A[2*i+k])*(m.x[k].c[0][0]-S(physical_center[k]));physical=add(physical,mul(fabsq(c.A[2*i+k]),complex_radius[k]));}c.invariance=std::max(c.invariance,add(mag(offset),add(physical,mul(q,c.radius))));}
    need(c.contraction<1&&c.invariance<c.radius,"one complex contraction box proves inverse for full target and Cauchy discs");
    // chi-I Lipschitz in infinity norm on the whole convex source triangle.
    Q full_radius=0;for(int a=0;a<3;a++)for(int k=0;k<2;k++)full_radius=std::max(full_radius,mag(S(Q(r.triangle.geometry.x[a][k]))-S(m.anchor[k])));full_radius=mf::FP<Q>::up(full_radius);
    for(int i=0;i<2;i++){Q row=0;for(int k=0;k<2;k++){P d=m.F[2*i+k];d.c[0][0]=d.c[0][0]-S(i==k?1:0);row=add(row,norm(d,full_radius));}c.global_injectivity=std::max(c.global_injectivity,row);}need(c.global_injectivity<1,"original map globally injective on the convex source");
    P remainder=m.J;remainder.c[0][0]=S(0);Q variation=norm(remainder,c.radius);c.J_lower=mf::FP<Q>::down(lower(m.J.c[0][0])-variation);c.J_upper=add(mag(m.J.c[0][0]),variation);need(c.J_lower>0,"positive real J and complex zero exclusion");
    // The identical lens uses the exact unit reference triangle.
    need(r.triangle.geometry.x[0][0]==0&&r.triangle.geometry.x[0][1]==0&&r.triangle.geometry.x[1][0]==1&&r.triangle.geometry.x[1][1]==0&&r.triangle.geometry.x[2][0]==0&&r.triangle.geometry.x[2][1]==1,"retained lens source domain");
    c.source_margin=std::min({lower(S(m.anchor[0])-S(c.radius)),lower(S(m.anchor[1])-S(c.radius)),lower(S(1)-S(m.anchor[0])-S(m.anchor[1])-S(2)*S(c.radius))});need(c.source_margin>0,"entire inverse box inside original source; no gluing or source baseline for this target");return c;
}
static std::array<S,2> inverse(const curved::Record&r,const PolynomialMap&m,const InverseCertificate&c,S x,S y,double eps){
    Q X=m.anchor[0],Y=m.anchor[1];for(int it=0;it<12;it++){auto p=mf::original_at(r,S(X),S(Y),eps);Q dx=p.position[0].n-x.n,dy=p.position[1].n-y.n;Q a=c.A[0]*dx+c.A[1]*dy,b=c.A[2]*dx+c.A[3]*dy;X-=a;Y-=b;if(std::max(fabsq(a),fabsq(b))<1e-29Q)break;}
    auto p=mf::original_at(r,S(X),S(Y),eps);S residual[2]={p.position[0]-x,p.position[1]-y};Q e=0;for(int i=0;i<2;i++){S z=S(c.A[2*i])*residual[0]+S(c.A[2*i+1])*residual[1];e=std::max(e,mag(z));}e=mf::FP<Q>::up(e/mf::FP<Q>::down(1-c.contraction));Q distance=std::max(mag(S(X)-S(m.anchor[0])),mag(S(Y)-S(m.anchor[1])));need(add(distance,e)<c.radius,"certified inverse point remains inside proof box");return {S(X,e),S(Y,e)};
}
static mf::B<R> cast(S q){R n=R(q.n);Q e=add(fabsq(Q(n)-q.n),q.e);R er=nextafterl(R(e),INFINITY);return {n,er};}
static mf::Point<R> evaluate(const PolynomialMap&m,S X,S Y){
    using B=mf::B<R>;using D=mf::D<R>;D x(cast(X-S(m.anchor[0])),B(1),B(0)),y(cast(Y-S(m.anchor[1])),B(0),B(1));
    auto poly=[&](const P&p){D z;for(int i=p.degree;i>=0;i--){D row;for(int j=p.degree-i;j>=0;j--)row=mf::fused(row,y,D(cast(p.c[i][j])));z=mf::fused(z,x,row);}return z;};mf::Point<R> out;for(int k=0;k<2;k++){D p=poly(m.x[k]),v=poly(m.v[k]);out.position[k]=p.v;out.velocity[k]=v.v;out.F[2*k]=p.x;out.F[2*k+1]=p.y;out.V[2*k]=v.x;out.V[2*k+1]=v.y;}mf::finish(out);return out;
}
struct Gauss {std::array<S,16> x,w;};
static std::pair<S,S> legendre(S z){S a(1),b=z,da,db(1);for(int k=2;k<=16;k++){S c=(S(2*k-1)*z*b-S(k-1)*a)/S(k),dc=(S(2*k-1)*(b+z*db)-S(k-1)*da)/S(k);a=b;b=c;da=db;db=dc;}return {b,db};}
static Gauss gauss16(){Gauss g;Q pi=acosq(-1);for(int i=0;i<8;i++){Q z=cosq(pi*(Q(i)+Q(.75L))/Q(16.5L));for(int it=0;it<24;it++){auto p=legendre(S(z));Q next=z-p.first.n/p.second.n;if(next==z)break;z=next;}Q e=1e-29Q;auto l=legendre(S(z-e)),r=legendre(S(z+e));need(upper(l.first)<0&&lower(r.first)>0||lower(l.first)>0&&upper(r.first)<0,"Gauss root sign bracket");auto p=legendre(S(z,e));need(lower(p.second)>0||upper(p.second)<0,"Gauss root unique in bracket");S node(z,e),weight=S(1)/((S(1)-node*node)*p.second*p.second);g.x[i]=(S(1)-node)/S(2);g.x[15-i]=(S(1)+node)/S(2);g.w[i]=g.w[15-i]=weight;need(lower(weight)>0,"positive certified Gauss weight");}for(int i=1;i<16;i++)need(upper(g.x[i-1])<lower(g.x[i]),"all16 distinct roots cover degree16");return g;}
}
#endif
