#ifndef MATERIAL_SHARED_OBLIQUE_POINT_H
#define MATERIAL_SHARED_OBLIQUE_POINT_H
#include "material_shared_oblique_slab.h"
namespace oblique_point {
namespace os=oblique_slab;
struct Linear {S s[2],q[2];int pivot=0;};
static Linear linear(const sd::Compiled&m,const Frame29&f){Linear out;int i=(f.original.apex+1)%3,j=(f.original.apex+2)%3;S Bi=f.original.vertex[1][i],Bj=f.original.vertex[1][j],Ci=f.original.vertex[2][i],Cj=f.original.vertex[2][j],det=Bi*Cj-Bj*Ci;need(sf::lower(det)>0,"original linear chart determinant");for(int axis=0;axis<2;axis++){const S*D=axis?f.direction:f.transverse;S li,lj;for(int k=0;k<2;k++){li=li+m.base.grad[i][k]*D[k];lj=lj+m.base.grad[j][k]*D[k];}S a=(li*Cj-lj*Ci)/det,b=(Bi*lj-Bj*li)/det;out.s[axis]=a+b;out.q[axis]=b;}out.pivot=fabsq(out.s[1].n)>fabsq(out.s[0].n)?1:0;need(fabsq(out.s[out.pivot].n)>out.s[out.pivot].e,"fixed oblique algebraic pivot excludes zero");return out;}
struct Coordinate {S s,t;sd::AffineAngular angular;};
static Coordinate coordinate(const Linear&l,S u,S v,const S*collected_s=nullptr){Coordinate out;S input[2]={u,v};out.s=collected_s?*collected_s:l.s[0]*u+l.s[1]*v;int p=l.pivot,other=1-p;out.angular.slope=l.q[p]/l.s[p];out.angular.intercept=(l.q[other]-out.angular.slope*l.s[other])*input[other];if(!(sf::lower(out.s)>0))throw sd::DomainFailure("collected oblique regular-Duffy radial pole excluded",out.s,Q(0));out.t=out.angular.slope+out.angular.intercept/out.s;return out;}
static Coordinate affine_coordinate(const Linear&l,S u,S beta,S gamma){
 Coordinate out;auto&a=out.angular;a.oblique_collected=true;a.common=u;a.beta=beta;a.gamma=gamma;a.su=l.s[0];a.sv=l.s[1];a.qu=l.q[0];a.qv=l.q[1];
 out.s=(a.su+a.sv*beta)*u+a.sv*gamma;
 if(!(sf::lower(out.s)>0))throw sd::DomainFailure("collected oblique regular-Duffy radial pole excluded",out.s,Q(0));
 a.intercept=(a.qu+a.qv*beta)*u+a.qv*gamma;out.t=a.intercept/out.s;return out;
}
static sd::Jet face(const target::Face&p,const Frame29&f,const Linear&l,S u,S v,const Coordinate*collected=nullptr){auto c=collected?*collected:coordinate(l,u,v);auto a=target::face_q(p,f.original,c.s,c.angular);sd::Jet out(a.v);S D[2][2]={{l.s[0],l.q[0]},{l.s[1],l.q[1]}};for(int i=0;i<2;i++)for(int k=0;k<2;k++){out.g[i]=out.g[i]+a.g[k]*D[i][k];for(int j=0;j<2;j++)for(int h=0;h<2;h++)out.H[i][j]=out.H[i][j]+a.H[k][h]*D[i][k]*D[j][h];}return out;}
struct Graph {target::Implicit certificate;S slope,intercept,correction;};
#include "material_shared_oblique_v_domain.h"
static Graph implicit(const sd::Compiled&m,const target::Face&p,const Frame29&f,const Linear&l,const os::Slab&slab,int index,S u){
 S u0(u.n),lo=os::edge_v(f,slab.geometric_lower,u0),hi=os::edge_v(f,slab.geometric_upper,u0);os::Proof real={};real.derivative=f.derivative[index];S root=os::graph(m,p,f,u0,lo,hi,real);auto at=face(p,f,l,u0,S(root.n));Q A=Q(1)/at.g[1].n,slope=(-at.g[0]/at.g[1]).n;S delta=u-u0,predictor=S(root.n)+S(slope)*delta;Q prior=[&]{sd::TracePhase phase(65);return geometric_v_prior(m,f,l,index,u,root,slope);}();auto collected=affine_coordinate(l,u,S(slope),S(root.n)-S(slope)*u0+S(0,prior));auto bound=[&]{sd::TracePhase phase(66);return face(p,f,l,u,S(predictor.n,sf::add(predictor.e,prior)),&collected);}();
 S first=at.g[0]+S(slope)*at.g[1],second=bound.H[0][0]+S(2)*S(slope)*bound.H[0][1]+S(slope)*S(slope)*bound.H[1][1];Q D=sf::mag(delta),aa=fabsq(A),residual=sf::mul(aa,sf::add(sf::mag(at.v),sf::add(sf::mul(sf::mag(first),D),sf::mul(Q(.5L),sf::mul(sf::mag(second),sf::mul(D,D))))));S mixed=bound.H[1][0]+S(slope)*bound.H[1][1];Q g=sf::add(sf::mag(S(1)-S(A)*at.g[1]),sf::mul(aa,sf::mul(sf::mag(mixed),D))),h=sf::mul(aa,sf::mag(bound.H[1][1])),gap=mf::FP<Q>::down(Q(1)-g);
 if(!(gap>0))throw target::GraphFailure("oblique affine derivative cone unproved",g,0,residual,prior,prior,root.e);Q R=std::max(root.e,sf::quotient(sf::mul(Q(2),residual),gap));
 // Same existing generic graph_disk sufficient-prior construction.
 if(R>prior){Q q=sf::add(g,sf::mul(h,prior)),inv=sf::quotient(sf::add(residual,sf::mul(q,prior)),prior);if(q<1&&inv<1)R=prior;}
 Q q=sf::add(g,sf::mul(h,R)),inv=R>0?sf::quotient(sf::add(residual,sf::mul(q,R)),R):__builtin_huge_valq();if(!(R>0&&R<=prior))throw target::GraphFailure("oblique correction exceeds original Hessian neighborhood",q,inv,residual,R,prior,root.e);if(!(q<1&&inv<1))throw target::GraphFailure("oblique affine contraction/invariance unproved",q,inv,residual,R,prior,root.e);Graph out;out.certificate={S(predictor.n,sf::add(predictor.e,R)),q,inv,residual,R,prior,root.e};out.slope=S(slope);out.intercept=S(root.n)-S(slope)*u0;out.correction=S(0,R);return out;
}
struct Point {sd::Point value;S u,v,s,t,width;Q contraction=0,invariance=0,forcing=0,radius=0,prior=0;};
struct Branches {Graph lo,hi;};
static Branches prepare(const sd::Compiled&m,const std::vector<target::Face>&faces,const Frame29&f,const Linear&l,const os::Slab&slab,S u){auto branch=[&](int index,int side){Graph g;if(index<0){int b=(side+1)%4;g.slope=(f.projected[b][1]-f.projected[side][1])/(f.projected[b][0]-f.projected[side][0]);g.intercept=f.projected[side][1]-g.slope*f.projected[side][0];return g;}return implicit(m,faces[index],f,l,slab,index,u);};return {branch(slab.lower,slab.geometric_lower),branch(slab.upper,slab.geometric_upper)};}
static Point point(const sd::Compiled&m,const Frame29&f,const Linear&l,S u,S eta,const Branches&branches){Point out;out.u=u;const auto&lo=branches.lo;const auto&hi=branches.hi;for(const auto*g:{&lo,&hi}){const auto&p=g->certificate;out.contraction=std::max(out.contraction,p.contraction);out.invariance=std::max(out.invariance,p.invariance);out.forcing=std::max(out.forcing,p.forcing);out.radius=std::max(out.radius,p.correction_radius);out.prior=std::max(out.prior,p.prior_radius);}S left=S(1)-S(eta.n),right(eta.n),z(0,eta.e),db=hi.slope-lo.slope,dg=(hi.intercept+hi.correction)-(lo.intercept+lo.correction);out.width=db*u+dg;
 // Same weighted affine identity used in the previous accepted ray proof.
 S beta=left*lo.slope+right*hi.slope+z*db,gamma=left*(lo.intercept+lo.correction)+right*(hi.intercept+hi.correction)+z*dg;out.v=beta*u+gamma;auto c=affine_coordinate(l,u,beta,gamma);out.s=c.s;out.t=c.t;out.value=sd::evaluate(m,f.original,c.s,c.t,&c.angular);return out;}

}
#endif
