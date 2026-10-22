#ifndef MATERIAL_TARGET_GENERAL_H
#define MATERIAL_TARGET_GENERAL_H
#include "material_target_inverse.h"
#include <type_traits>
/* Exact original P3 + initial-v0 residual map, expanded algebraically at
 * a local reference anchor. Q128 intervals enclose conversion; neither
 * stored coefficients nor source state are changed. */
namespace material_target_general {
namespace ti=material_target_inverse;namespace mf=material_forward;
using Q=ti::Q;using S=ti::S;using P=ti::P;using R=long double;
struct Map {P x[2],v[2];S Bx[2][2],Bv[2][2];Q anchor[2];double eps;};
static Map compile(const curved::Record&r,const Q*anchor,double eps){
    Map out;out.anchor[0]=anchor[0];out.anchor[1]=anchor[1];out.eps=eps;
    P dx=ti::constant(S(anchor[0])-S(Q(r.triangle.geometry.x[0][0]))),dy=ti::constant(S(anchor[1])-S(Q(r.triangle.geometry.x[0][1])));dx.degree=dy.degree=1;dx.c[1][0]=S(1);dy.c[0][1]=S(1);P l[3],N[10];
    l[1]=ti::plus(ti::scale(dx,S(Q(r.triangle.geometry.grad[1][0]))),ti::scale(dy,S(Q(r.triangle.geometry.grad[1][1]))));l[2]=ti::plus(ti::scale(dx,S(Q(r.triangle.geometry.grad[2][0]))),ti::scale(dy,S(Q(r.triangle.geometry.grad[2][1]))));l[0]=ti::minus(ti::minus(ti::constant(S(1)),l[1]),l[2]);
    for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);N[a]=ti::constant(S(1));for(int k=0;k<3;k++)for(int j=0;j<index[k];j++)N[a]=ti::product(N[a],ti::scale(ti::minus(ti::scale(l[k],S(3)),ti::constant(S(j))),S(1)/S(index[k]-j)));}
    S ax(Q(r.triangle.geometry.x[0][0])),ay(Q(r.triangle.geometry.x[0][1])),z=ax*ax+ay*ay+S(Q(eps))*S(Q(eps)),omega=mf::power(z,Q(-.75L));
    P L=ti::scale(ti::plus(ti::scale(dx,S(2)*ax),ti::scale(dy,S(2)*ay)),S(1)/z),W=ti::scale(ti::plus(ti::product(dx,dx),ti::product(dy,dy)),S(1)/z),L2=ti::product(L,L);
    P A=ti::plus(ti::minus(ti::minus(ti::constant(S(1)),ti::scale(L,S(.75L))),ti::scale(W,S(.75L))),ti::minus(ti::plus(ti::scale(L2,S(.65625L)),ti::scale(ti::product(L,W),S(1.3125L))),ti::scale(ti::product(L2,L),S(.6015625L))));
    P C=ti::plus(ti::minus(ti::minus(ti::constant(S(1)),ti::scale(L,S(.75L))),ti::scale(W,S(.75L))),ti::scale(L2,S(.65625L)));
    P initial[2]={ti::scale(ti::minus(ti::scale(A,-ay),ti::product(dy,C)),omega),ti::scale(ti::plus(ti::scale(A,ax),ti::product(dx,C)),omega)};
    for(int j=0;j<2;j++)for(int a=0;a<10;a++)initial[j]=ti::plus(initial[j],ti::scale(N[a],S(Q(r.triangle.remainder[a][j]))));
    for(int k=0;k<2;k++){out.x[k]=ti::constant(S(Q(r.origin[k])));for(int a=0;a<10;a++){out.x[k]=ti::plus(out.x[k],ti::scale(N[a],S(Q(r.map.x[k][a]))));out.v[k]=ti::plus(out.v[k],ti::scale(N[a],S(Q(r.map.v[k][a]))));}for(int j=0;j<2;j++){out.Bx[k][j]=S(Q(r.map.x[k][10+j]))/S(Q(r.scale[j]));out.Bv[k][j]=S(Q(r.map.v[k][10+j]))/S(Q(r.scale[j]));out.x[k]=ti::minus(out.x[k],ti::scale(initial[j],out.Bx[k][j]));out.v[k]=ti::minus(out.v[k],ti::scale(initial[j],out.Bv[k][j]));}}
    return out;
}
template<class T> static mf::B<T> adapt(S q){if constexpr(std::is_same<T,Q>::value)return q;else return ti::cast(q);}
template<class T> static mf::Point<T> evaluate(const Map&m,S X,S Y){
    using B=mf::B<T>;using D=mf::D<T>;D dx(adapt<T>(X-S(m.anchor[0])),B(1),B(0)),dy(adapt<T>(Y-S(m.anchor[1])),B(0),B(1));
    auto poly=[&](const P&p){D z;for(int i=p.degree;i>=0;i--){D row;for(int j=p.degree-i;j>=0;j--)row=mf::fused(row,dy,D(adapt<T>(p.c[i][j])));z=mf::fused(z,dx,row);}return z;};
    D wx(adapt<T>(X),B(1),B(0)),wy(adapt<T>(Y),B(0),B(1)),omega=mf::power(wx*wx+wy*wy+D(B(T(m.eps))*B(T(m.eps))),T(-.75L)),base[2]={-wy*omega,wx*omega};mf::Point<T> out;
    for(int k=0;k<2;k++){D x=poly(m.x[k]),v=poly(m.v[k]);for(int j=0;j<2;j++){x=x+D(adapt<T>(m.Bx[k][j]))*base[j];v=v+D(adapt<T>(m.Bv[k][j]))*base[j];}out.position[k]=x.v;out.velocity[k]=v.v;out.F[2*k]=x.x;out.F[2*k+1]=x.y;out.V[2*k]=v.x;out.V[2*k+1]=v.y;}mf::finish(out);return out;
}
struct Jet {mf::Point<Q> point;S H[2][2][2];};
static Jet jet(const Map&m,S X,S Y){Jet out;out.point=evaluate<Q>(m,X,Y);S c[2]={X,Y},rot[2]={-Y,X},z=X*X+Y*Y+S(Q(m.eps))*S(Q(m.eps)),omega=mf::power(z,Q(-.75L)),od[2],odd[2][2];int E[2][2]={{0,-1},{1,0}};for(int i=0;i<2;i++){od[i]=S(-1.5L)*c[i]*omega/z;for(int j=0;j<2;j++)odd[i][j]=S(-1.5L*(i==j))*omega/z+S(5.25L)*c[i]*c[j]*omega/(z*z);}
    S dx=X-S(m.anchor[0]),dy=Y-S(m.anchor[1]);for(int k=0;k<2;k++)for(int i=0;i<2;i++)for(int j=0;j<2;j++){out.H[k][i][j]=ti::value(ti::derivative(ti::derivative(m.x[k],i),j),dx,dy);for(int a=0;a<2;a++)out.H[k][i][j]=out.H[k][i][j]+m.Bx[k][a]*(rot[a]*odd[i][j]+S(E[a][i])*od[j]+S(E[a][j])*od[i]);}return out;
}
struct Inverse {S X,Y;Q radius=0,contraction=0,global_L=0;bool inside=false,outside=false;};
static Inverse inverse(const Map&m,const curved::Record&r,S x,S y){
    Q X=m.anchor[0],Y=m.anchor[1];for(int it=0;it<16;it++){auto p=evaluate<Q>(m,S(X),S(Y));Q det=p.F[0].n*p.F[3].n-p.F[1].n*p.F[2].n,dx=p.position[0].n-x.n,dy=p.position[1].n-y.n,a=(p.F[3].n*dx-p.F[1].n*dy)/det,b=(-p.F[2].n*dx+p.F[0].n*dy)/det;X-=a;Y-=b;if(std::max(fabsq(a),fabsq(b))<1e-29Q)break;}
    auto center=evaluate<Q>(m,S(X),S(Y));S F[4];for(int k=0;k<4;k++)F[k]=center.F[k];S det=F[0]*F[3]-F[1]*F[2];Q A[4]={(F[3]/det).n,(-F[1]/det).n,(-F[2]/det).n,(F[0]/det).n};Q proof_radius=1e-5Q;auto box=evaluate<Q>(m,S(X,proof_radius),S(Y,proof_radius));Inverse out;Q residual=0;S f[2]={center.position[0]-x,center.position[1]-y};for(int i=0;i<2;i++){Q row=0;for(int j=0;j<2;j++){S z(i==j?1:0);for(int k=0;k<2;k++)z=z-S(A[2*i+k])*box.F[2*k+j];row=ti::add(row,ti::mag(z));}out.contraction=std::max(out.contraction,row);S r0=S(A[2*i])*f[0]+S(A[2*i+1])*f[1];residual=std::max(residual,ti::mag(r0));}need(out.contraction<1&&ti::add(residual,ti::mul(out.contraction,proof_radius))<proof_radius,"uniform original-map inverse contraction");out.radius=mf::FP<Q>::up(residual/mf::FP<Q>::down(1-out.contraction));out.X=S(X,out.radius);out.Y=S(Y,out.radius);
    // Exclude a second source preimage on the convex hull of T and the box.
    S global[2];for(int k=0;k<2;k++){Q lo=(k?Y:X)-proof_radius,hi=(k?Y:X)+proof_radius;for(int a=0;a<3;a++){lo=std::min(lo,Q(r.triangle.geometry.x[a][k]));hi=std::max(hi,Q(r.triangle.geometry.x[a][k]));}S mid=(S(lo)+S(hi))/S(2);global[k]=S(mid.n,ti::add((S(hi)-S(lo)).n/2,mid.e));}auto whole=evaluate<Q>(m,global[0],global[1]);for(int i=0;i<2;i++){Q row=0;for(int j=0;j<2;j++)row=ti::add(row,ti::mag(whole.F[2*i+j]-S(i==j?1:0)));out.global_L=std::max(out.global_L,row);}need(out.global_L<1,"whole original source plus inverse box is injective");
    S ax(Q(r.triangle.geometry.x[0][0])),ay(Q(r.triangle.geometry.x[0][1])),ux=S(Q(r.triangle.geometry.x[1][0]))-ax,uy=S(Q(r.triangle.geometry.x[1][1]))-ay,vx=S(Q(r.triangle.geometry.x[2][0]))-ax,vy=S(Q(r.triangle.geometry.x[2][1]))-ay,delta=ux*vy-uy*vx,dx=out.X-ax,dy=out.Y-ay,l1=(dx*vy-dy*vx)/delta,l2=(ux*dy-uy*dx)/delta,l0=S(1)-l1-l2;out.inside=ti::lower(l0)>0&&ti::lower(l1)>0&&ti::lower(l2)>0;out.outside=ti::upper(l0)<0||ti::upper(l1)<0||ti::upper(l2)<0;return out;
}
struct CurvePoint {S x[2],d[2],dd[2],X,Y;};
static CurvePoint edge(const Map&m,const curved::Record&r,int a,S t){int b=(a+1)%3;S A[2]={S(Q(r.triangle.geometry.x[a][0])),S(Q(r.triangle.geometry.x[a][1]))},D[2]={S(Q(r.triangle.geometry.x[b][0]))-A[0],S(Q(r.triangle.geometry.x[b][1]))-A[1]};CurvePoint out;out.X=A[0]+t*D[0];out.Y=A[1]+t*D[1];auto q=jet(m,out.X,out.Y);for(int k=0;k<2;k++){out.x[k]=q.point.position[k];out.d[k]=q.point.F[2*k]*D[0]+q.point.F[2*k+1]*D[1];for(int i=0;i<2;i++)for(int j=0;j<2;j++)out.dd[k]=out.dd[k]+q.H[k][i][j]*D[i]*D[j];}return out;
}
struct Plane {S n[2],d;int competitor=-1;};
static S exact_difference(Q a,Q b){Q s=a-b,z=s-a,e=(a-(s-z))+(-b-z);if(e==0)return S(s);return S(a)-S(b);}
static Plane plane(const std::vector<curved_partition::Site>&sites,const std::vector<R>&weights,int own,int other,const Q*origin){Plane p;for(int k=0;k<2;k++)p.n[k]=exact_difference(Q(k?sites[other].y:sites[other].x),Q(k?sites[own].y:sites[own].x));S dx=p.n[0],dy=p.n[1];p.d=(dx*dx+dy*dy+S(Q(weights[own]))-S(Q(weights[other])))/S(2);p.d=p.d+p.n[0]*(S(Q(sites[own].x))-S(origin[0]))+p.n[1]*(S(Q(sites[own].y))-S(origin[1]));p.competitor=other;return p;}
static S value(const Plane&p,S x,S y,const Q*origin){return p.n[0]*(x-S(origin[0]))+p.n[1]*(y-S(origin[1]))-p.d;}
struct Root {S t;int plane=-1;};
struct Roots {std::vector<Root> root;uint64_t evaluations=0;bool complete=true;};
static Roots roots(const Map&m,const curved::Record&r,int a,const Plane&p,const Q*origin,int plane_id){
    struct Window {Q lo,hi;int depth;bool exists=false;};std::vector<Window> stack={{0,1,0,false}};Roots out;while(!stack.empty()){auto w=stack.back();stack.pop_back();if(out.evaluations++>=256){out.complete=false;return out;}S mid=(S(w.lo)+S(w.hi))/S(2);Q radius=ti::mag(S(w.hi)-S(w.lo))/2;auto q=edge(m,r,a,S(mid.n,ti::add(radius,mid.e)));S f=value(p,q.x[0],q.x[1],origin),df=p.n[0]*q.d[0]+p.n[1]*q.d[1];if(ti::lower(f)>0||ti::upper(f)<0)continue;
        auto left=edge(m,r,a,S(w.lo)),right=edge(m,r,a,S(w.hi));S fl=value(p,left.x[0],left.x[1],origin),fr=value(p,right.x[0],right.x[1],origin);bool separated=ti::lower(df)>0||ti::upper(df)<0;
        if(separated){bool opposite=ti::upper(fl)<0&&ti::lower(fr)>0||ti::lower(fl)>0&&ti::upper(fr)<0;if(!opposite&&((ti::lower(fl)>0&&ti::lower(fr)>0)||(ti::upper(fl)<0&&ti::upper(fr)<0)))continue;
            auto center=edge(m,r,a,S(mid.n));S fc=value(p,center.x[0],center.x[1],origin),N=S(mid.n)-fc/df;Q lo=std::max(w.lo,ti::lower(N)),hi=std::min(w.hi,ti::upper(N));if(lo>hi)continue;
            if((opposite||w.exists)&&hi-lo<1e-27Q){S t=(S(lo)+S(hi))/S(2);t.e=ti::add(t.e,ti::mag(S(hi)-S(lo))/2);out.root.push_back({t,plane_id});continue;}
            if(hi-lo<Q(.75L)*(w.hi-w.lo)){stack.push_back({lo,hi,w.depth+1,opposite||w.exists});continue;}}
        if(w.depth>=64){out.complete=false;return out;}Q c=mid.n;need(c>w.lo&&c<w.hi,"finite boundary parameter chart");stack.push_back({c,w.hi,w.depth+1});stack.push_back({w.lo,c,w.depth+1});}
    std::sort(out.root.begin(),out.root.end(),[](const Root&a,const Root&b){return a.t.n<b.t.n;});for(size_t k=1;k<out.root.size();k++)need(ti::upper(out.root[k-1].t)<ti::lower(out.root[k].t),"distinct certified source-edge intersections");return out;
}
}
#endif
