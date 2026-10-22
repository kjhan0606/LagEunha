#ifndef MATERIAL_SHARED_FUNCTION_H
#define MATERIAL_SHARED_FUNCTION_H
#include <limits>
#include "material_forward_ball.h"
// Canonical function specification, independent of native qnode/gradient
// caches. Q128 outward balls use the retained <=8ulp libm engineering model.
namespace shared_function {
namespace mf=material_forward;using Q=__float128;using S=mf::B<Q>;using D=mf::D<Q>;
static Q lower(S x){return mf::FP<Q>::down(x.n-x.e);}static Q upper(S x){return mf::FP<Q>::up(x.n+x.e);}static Q mag(S x){return mf::plus_up(fabsq(x.n),x.e);}
static mf::B<double> cache(double native,S original){Q delta=mf::plus_up(fabsq(Q(native)-original.n),original.e);double e=nextafter(double(delta),INFINITY);return {native,e};}
static S exact(double x){return S(Q(x));}
static D ipow(D x,int n){D out(S(1));for(int j=0;j<n;j++)out=out*x;return out;}
static D lagrange(int j,D t){D out(S(1));for(int k=0;k<4;k++)if(j!=k)out=out*(D(S(3))*t-D(S(k)))/D(S(j-k));return out;}
static std::array<D,2> velocity(D X,D Y,double eps){auto w=mf::power(X*X+Y*Y+D(exact(eps)*exact(eps)),Q(-.75L));return {-Y*w,X*w};}
// The Taylor polynomial is used only to certify a rounded native cache.
// Canonical functions below use v0 and exact rational interpolation nodes.
static std::array<S,2> remainder(const double*anchor,S X,S Y,double eps){
    S a=exact(anchor[0]),b=exact(anchor[1]),dx=X-a,dy=Y-b,z=a*a+b*b+exact(eps)*exact(eps),L=S(2)*(a*dx+b*dy)/z,W=(dx*dx+dy*dy)/z;
    S A=S(1)-S(.75L)*L-S(.75L)*W+S(.65625L)*L*L+S(1.3125L)*L*W-S(.6015625L)*L*L*L,C=S(1)-S(.75L)*L-S(.75L)*W+S(.65625L)*L*L,wa=mf::power(z,Q(-.75L)),w=mf::power(X*X+Y*Y+exact(eps)*exact(eps),Q(-.75L));return {-Y*w-wa*(-b*A-dy*C),X*w-wa*(a*A+dx*C)};
}
struct Edge {
    STEdge native;S A[2],Dref[2],node[4][2];mf::B<double> qcache[4][2],Dcache[2];
};
static Edge edge(const STEdge&e,const STTriangle&t,int local){
    Edge out={};out.native=e;int a=t.endpoint[local][0],b=t.endpoint[local][1];for(int k=0;k<2;k++){out.A[k]=exact(t.vertex[a][k]);out.Dref[k]=exact(t.vertex[b][k])-out.A[k];out.Dcache[k]=cache(e.D[k],out.Dref[k]);}
    for(int j=0;j<4;j++){S X=out.A[0]+S(j)*out.Dref[0]/S(3),Y=out.A[1]+S(j)*out.Dref[1]/S(3);auto v=velocity(D(X),D(Y),e.eps);auto q=remainder(e.A,X,Y,e.eps);for(int k=0;k<2;k++){out.node[j][k]=v[k].v;out.qcache[j][k]=cache(e.qnode[j][k],q[k]);}}return out;
}
struct Map {
    const st::Factory*state=nullptr;size_t source=0;S vertex[3][2],grad[3][2],node[10][2],nodal[2][2][10],bulk[2][2][2],extension[2][2][3][2];
    Edge edge[3];mf::B<double> qcache[10][2],gradcache[3][2];double eps=0;Q anchor[2]={};
};
static Map compile(const st::Factory&f,size_t source,const std::vector<STCoefficient>&coefficient){
    need(source<f.triangle.size()&&coefficient.size()==f.coefficient.size(),
        "source map uses the same complete native basis with selected current or generating coefficients");
    Map m={};m.state=&f;m.source=source;const auto&t=f.triangle[source];m.eps=t.eps;
    for(int a=0;a<3;a++)for(int k=0;k<2;k++)m.vertex[a][k]=exact(t.vertex[a][k]);S dx=m.vertex[1][0]-m.vertex[0][0],dy=m.vertex[1][1]-m.vertex[0][1],ex=m.vertex[2][0]-m.vertex[0][0],ey=m.vertex[2][1]-m.vertex[0][1],det=dx*ey-dy*ex;need(lower(det)>0,"canonical exact positive reference geometry");m.grad[1][0]=ey/det;m.grad[1][1]=-ex/det;m.grad[2][0]=-dy/det;m.grad[2][1]=dx/det;for(int k=0;k<2;k++){m.grad[0][k]=-m.grad[1][k]-m.grad[2][k];m.anchor[k]=((m.vertex[0][k]+m.vertex[1][k]+m.vertex[2][k])/S(3)).n;for(int a=0;a<3;a++)m.gradcache[a][k]=cache(t.grad[a][k],m.grad[a][k]);}
    for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);S x,y;for(int j=0;j<3;j++){x=x+S(index[j])*m.vertex[j][0]/S(3);y=y+S(index[j])*m.vertex[j][1]/S(3);}auto v=velocity(D(x),D(y),m.eps);auto q=remainder(t.vertex[0],x,y,m.eps);for(int k=0;k<2;k++){m.node[a][k]=v[k].v;m.qcache[a][k]=cache(t.qnode[a][k],q[k]);}}
    for(int e=0;e<3;e++)m.edge[e]=edge(t.edge[e],t,e);
    // Algebraic condensation of the SAME stored state. No sampled fit,
    // coefficient replacement, component tie, or altered enrichment.
    for(int mode=0;mode<2;mode++)for(int k=0;k<2;k++){
        for(int j=0;j<2;j++)for(int c=0;c<2;c++){const auto&q=coefficient[t.id[t.bubble_offset+j]];m.bulk[mode][k][c]=m.bulk[mode][k][c]+exact(mode?q.v[k]:q.x[k])*exact(t.bubble_encode[j][c]);}
        for(int a=0;a<10;a++){const auto&q=coefficient[t.id[a]];S v=exact(mode?q.v[k]:q.x[k]);for(int c=0;c<2;c++)v=v-m.bulk[mode][k][c]*m.node[a][c];if(a==9)for(int c=0;c<2;c++)v=v-m.bulk[mode][k][c]*exact(t.bubble_shift[c]);m.nodal[mode][k][a]=v;}
        for(int e=0;e<3;e++)for(int c=0;c<2;c++){S value=-m.bulk[mode][k][c];for(int j=0;j<t.edge[e].rank;j++){const auto&q=coefficient[t.id[t.offset[e]+j]];value=value+exact(mode?q.v[k]:q.x[k])*exact(t.edge[e].encode[j][c]);}m.extension[mode][k][e][c]=value;}
    }return m;
}
static Map compile(const st::Factory&f,size_t source){return compile(f,source,f.coefficient);}
static std::array<D,2> edge_raw(const Edge&e,D t){D X=D(e.A[0])+t*D(e.Dref[0]),Y=D(e.A[1])+t*D(e.Dref[1]);auto r=velocity(X,Y,e.native.eps);for(int j=0;j<4;j++){auto L=lagrange(j,t);for(int k=0;k<2;k++)r[k]=r[k]-L*D(e.node[j][k]);}return r;}
static std::array<D,3> barycentric(const Map&m,S X,S Y){D dx(X-m.vertex[0][0],S(1),S(0)),dy(Y-m.vertex[0][1],S(0),S(1));std::array<D,3>l;l[1]=D(m.grad[1][0])*dx+D(m.grad[1][1])*dy;l[2]=D(m.grad[2][0])*dx+D(m.grad[2][1])*dy;l[0]=D(S(1))-l[1]-l[2];return l;}
static mf::Point<Q> evaluate(const Map&m,S X,S Y){
    const auto&t=m.state->triangle[m.source];auto l=barycentric(m,X,Y);D N[10];for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);N[a]=D(S(1));for(int k=0;k<3;k++)for(int j=0;j<index[k];j++)N[a]=N[a]*(D(S(3))*l[k]-D(S(j)))/D(S(index[k]-j));}
    auto v0=velocity(D(X,S(1),S(0)),D(Y,S(0),S(1)),m.eps);D E[3][2];for(int e=0;e<3;e++){auto s=l[t.endpoint[e][0]]+l[t.endpoint[e][1]];need(lower(s.v)>0,"canonical extension chart excludes opposite vertex; use its canonical edge/vertex chart there");auto r=edge_raw(m.edge[e],l[t.endpoint[e][1]]/s);for(int k=0;k<2;k++)E[e][k]=ipow(s,4)*r[k];}
    mf::Point<Q>out;for(int mode=0;mode<2;mode++)for(int k=0;k<2;k++){D value;for(int a=0;a<10;a++)value=value+D(m.nodal[mode][k][a])*N[a];for(int c=0;c<2;c++)value=value+D(m.bulk[mode][k][c])*v0[c];for(int e=0;e<3;e++)for(int c=0;c<2;c++)value=value+D(m.extension[mode][k][e][c])*E[e][c];if(mode){out.velocity[k]=value.v;out.V[2*k]=value.x;out.V[2*k+1]=value.y;}else{out.position[k]=value.v;out.F[2*k]=value.x;out.F[2*k+1]=value.y;}}mf::finish(out);return out;
}
struct RealBounds {Q L=0,displacement[2]={},FminusI[4]={};};
static S range(Q lo,Q hi){S mid=(S(lo)+S(hi))/S(2),half=(S(hi)-S(lo))/S(2);return {mid.n,mf::plus_up(mid.e,mag(half))};}
static Q add(Q a,Q b){return mf::plus_up(a,b);}static Q mul(Q a,Q b){return mf::times_up(a,b);}
static std::array<S,2> reference_node(const Map&m,int a){int index[3];gfs_material_p3_index(a,index);std::array<S,2>out;for(int j=0;j<3;j++)for(int k=0;k<2;k++)out[k]=out[k]+S(index[j])*m.vertex[j][k]/S(3);return out;}
static std::array<Q,2> coordinate_bound(const Map&m){std::array<Q,2>out={};for(int a=0;a<3;a++)for(int k=0;k<2;k++)out[k]=std::max(out[k],mag(m.vertex[a][k]));return out;}
static RealBounds real_bounds(const Map&m){
    // Real triangle bound, including the C3 opposite-vertex extension.
    // Independent s,t in[0,1] bound its factored first derivative WITHOUT
    // division by s. This is not an analytic complex rectangle certificate.
    RealBounds out;const auto&t=m.state->triangle[m.source];D l[3],N[10];for(int a=0;a<3;a++)l[a]=D(S(Q(.5L),Q(.5L)),m.grad[a][0],m.grad[a][1]);for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);N[a]=D(S(1));for(int k=0;k<3;k++)for(int j=0;j<index[k];j++)N[a]=N[a]*(D(S(3))*l[k]-D(S(j)))/D(S(index[k]-j));}
    Q zlo=lower(exact(m.eps)*exact(m.eps));auto coord=coordinate_bound(m);Q omega=mag(mf::power(S(zlo),Q(-.75L))),vbound[2]={mul(coord[1],omega),mul(coord[0],omega)},Dvelocity[2][2];for(int k=0;k<2;k++)for(int j=0;j<2;j++)Dvelocity[k][j]=add(k!=j?omega:Q(0),mul(Q(1.5L),mul(mul(coord[1-k],coord[j]),mf::FP<Q>::up(omega/zlo))));
    // Improve the radial denominator with the exact triangle's axis-aligned
    // distance bound. It is a geometric bound, not an accuracy scan.
    Q radius2lo=0;for(int k=0;k<2;k++){Q lo=lower(m.vertex[0][k]),hi=upper(m.vertex[0][k]);for(int a=1;a<3;a++){lo=std::min(lo,lower(m.vertex[a][k]));hi=std::max(hi,upper(m.vertex[a][k]));}Q d=lo>0?lo:(hi<0?-hi:Q(0));radius2lo=mf::FP<Q>::down(radius2lo+mf::FP<Q>::down(d*d));}zlo=mf::FP<Q>::down(zlo+radius2lo);omega=mag(mf::power(S(zlo),Q(-.75L)));for(int k=0;k<2;k++){vbound[k]=mul(coord[1-k],omega);for(int j=0;j<2;j++)Dvelocity[k][j]=add(k!=j?omega:Q(0),mul(Q(1.5L),mul(mul(coord[1-k],coord[j]),mf::FP<Q>::up(omega/zlo))));}
    for(int k=0;k<2;k++){S reference0=reference_node(m,0)[k],p0=m.nodal[0][k][0]-reference0;S poly[3]={p0,S(),S()};for(int a=1;a<10;a++){S difference=(m.nodal[0][k][a]-reference_node(m,a)[k])-p0;poly[0]=poly[0]+difference*N[a].v;poly[1]=poly[1]+difference*N[a].x;poly[2]=poly[2]+difference*N[a].y;}out.displacement[k]=mag(poly[0]);for(int j=0;j<2;j++)out.FminusI[2*k+j]=mag(poly[1+j]);for(int c=0;c<2;c++){out.displacement[k]=add(out.displacement[k],mul(mag(m.bulk[0][k][c]),vbound[c]));for(int j=0;j<2;j++)out.FminusI[2*k+j]=add(out.FminusI[2*k+j],mul(mag(m.bulk[0][k][c]),Dvelocity[c][j]));}}
    for(int e=0;e<3;e++){const auto&edge=m.edge[e];Q edgecoord[2]={},edgezlo=lower(exact(m.eps)*exact(m.eps));for(int k=0;k<2;k++){S a=edge.A[k],b=a+edge.Dref[k];Q lo=std::min(lower(a),lower(b)),hi=std::max(upper(a),upper(b)),d=lo>0?lo:(hi<0?-hi:Q(0));edgecoord[k]=std::max(mag(a),mag(b));edgezlo=mf::FP<Q>::down(edgezlo+mf::FP<Q>::down(d*d));}Q eo=mag(mf::power(S(edgezlo),Q(-.75L)));D parameter(S(Q(.5L),Q(.5L)),S(1),S(0));Q residual[2],derivative[2];for(int c=0;c<2;c++){residual[c]=add(mul(edgecoord[1-c],eo),mag(edge.node[0][c]));derivative[c]=0;for(int j=0;j<2;j++){Q dg=add(c!=j?eo:Q(0),mul(Q(1.5L),mul(mul(edgecoord[1-c],edgecoord[j]),mf::FP<Q>::up(eo/edgezlo))));derivative[c]=add(derivative[c],mul(dg,mag(edge.Dref[j])));}for(int j=1;j<4;j++){auto L=lagrange(j,parameter);Q delta=mag(edge.node[j][c]-edge.node[0][c]);residual[c]=add(residual[c],mul(mag(L.v),delta));derivative[c]=add(derivative[c],mul(mag(L.x),delta));}}
        for(int j=0;j<2;j++){Q ds=mag(m.grad[t.endpoint[e][0]][j]+m.grad[t.endpoint[e][1]][j]),b=mag(m.grad[t.endpoint[e][1]][j]);for(int c=0;c<2;c++){Q bound=add(mul(mul(Q(4),ds),residual[c]),mul(add(b,ds),derivative[c]));for(int k=0;k<2;k++)out.FminusI[2*k+j]=add(out.FminusI[2*k+j],mul(mag(m.extension[0][k][e][c]),bound));}}
        for(int k=0;k<2;k++)for(int c=0;c<2;c++)out.displacement[k]=add(out.displacement[k],mul(mag(m.extension[0][k][e][c]),residual[c]));
    }for(int k=0;k<2;k++)out.L=std::max(out.L,add(out.FminusI[2*k],out.FminusI[2*k+1]));return out;
}
struct Inverse {S X,Y;Q radius=0,L=0;bool inside=false;};
static Inverse inverse_inside(const Map&m,const RealBounds&bound,S x,S y){
    need(bound.L<1,"canonical real whole-source injectivity/contraction bound");Q X=m.anchor[0],Y=m.anchor[1];for(int it=0;it<16;it++){auto p=evaluate(m,S(X),S(Y));Q det=p.F[0].n*p.F[3].n-p.F[1].n*p.F[2].n,fx=p.position[0].n-x.n,fy=p.position[1].n-y.n,a=(p.F[3].n*fx-p.F[1].n*fy)/det,b=(-p.F[2].n*fx+p.F[0].n*fy)/det;X-=a;Y-=b;if(std::max(fabsq(a),fabsq(b))<1e-29Q)break;}auto p=evaluate(m,S(X),S(Y));Q residual=std::max(mag(p.position[0]-x),mag(p.position[1]-y));Inverse out;out.L=bound.L;out.radius=mf::FP<Q>::up(Q(2)*residual/mf::FP<Q>::down(Q(1)-bound.L));out.X=S(X,out.radius);out.Y=S(Y,out.radius);auto l=barycentric(m,out.X,out.Y);out.inside=lower(l[0].v)>0&&lower(l[1].v)>0&&lower(l[2].v)>0;need(out.inside,"data-derived inverse contraction box remains strictly inside canonical source");return out;
}
// Immutable thermal generating geometry is copied once. The recipe applies
// ONLY to these saved control states, never to later current geometry.
struct FrozenThermal {
    const st::Factory generating;const double P0,gamma;const uint64_t epoch;
    FrozenThermal(const st::Factory&f,double pressure,double g,uint64_t e):generating(f),P0(pressure),gamma(g),epoch(e){}
    S per_reference_volume(size_t source,S X,S Y)const{auto m=compile(generating,source);auto p=evaluate(m,X,Y);return exact(P0)/(exact(gamma)-S(1))*mf::power(p.J,Q(1)-Q(gamma));}
};
}
#endif
