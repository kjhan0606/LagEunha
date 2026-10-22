#ifndef MATERIAL_SHARED_RUNTIME_REFERENCE_H
#define MATERIAL_SHARED_RUNTIME_REFERENCE_H
namespace runtime_reference {
struct Point {S X[2],lambda[3],s,t,measure;sd::AffineAngular angular;bool affine=false;};
static Point coordinates(const STTriangle&triangle,const sd::Chart&chart,S s,S t,const sd::AffineAngular*angular=nullptr){Point out;out.s=s;out.t=t;out.affine=bool(angular);if(angular)out.angular=*angular;
 for(int a=0;a<3;a++)out.lambda[a]=angular?sd::affine_coordinate(chart.vertex[0][a],chart.vertex[1][a],chart.vertex[2][a],s,*angular):chart.vertex[0][a]+s*((chart.vertex[1][a]-chart.vertex[0][a])+t*(chart.vertex[2][a]-chart.vertex[1][a]));
 for(int k=0;k<2;k++){S B,C;for(int a=0;a<3;a++){B=B+chart.vertex[1][a]*sf::exact(triangle.vertex[a][k]);C=C+chart.vertex[2][a]*sf::exact(triangle.vertex[a][k]);}S A=sf::exact(triangle.vertex[chart.apex][k]);out.X[k]=angular?sd::affine_coordinate(A,B,C,s,*angular):A+s*((B-A)+t*(C-B));}return out;
}
static Point whole(const STTriangle&t,const TypedPanel&panel,S u,S v){const auto&p=panel.geometry;need(p.whole,"runtime whole reference panel");const auto&g=p.ordinary;S ds=g.shi-g.slo,dt=g.thi-g.tlo,s=g.slo+ds*u,tau=g.tlo+dt*v;auto out=coordinates(t,p.reference,s,tau);out.measure=(s*p.reference.det)*ds*dt;return out;}
template<class Real>static Point regular(const STTriangle&t,const TypedPanel&panel,const Real*real,S free,S eta){const auto&p=panel.geometry;need(!p.whole,"runtime clipped reference panel");const auto&strip=p.strip;S df=strip.hi-strip.lo,de=p.eta_hi-p.eta_lo;Point out;
 if(strip.dependent==1){
  if(strip.lower>=0||strip.upper>=0){need(real&&real->passed,"saved angular actual Gauss branch");const auto&lo=real->angular_lo;const auto&hi=real->angular_hi;auto mixed=lo;mixed.interpolated=true;mixed.eta=eta;mixed.width_slope=hi.slope-lo.slope;mixed.width_intercept=hi.intercept-lo.intercept;mixed.width_correction=hi.correction-lo.correction;S width=free*mixed.width_slope+mixed.width_intercept+mixed.width_correction;S tau=target::weighted_angular(lo,hi,free,eta);out=coordinates(t,p.reference,free,tau,&mixed);out.measure=p.reference.det*width*df*de;}
  else{out=coordinates(t,p.reference,free,eta);out.measure=(free*p.reference.det)*df*de;}
 }else{S lo(0),hi(1);if(strip.lower>=0||strip.upper>=0){need(real&&real->passed,"saved radial actual Gauss branch");lo=real->radial_lo;hi=real->radial_hi;}S s=lo+eta*(hi-lo);out=coordinates(t,p.reference,s,free);out.measure=(s*p.reference.det)*(hi-lo)*df*de;}
 return out;
}
static Point oblique(const STTriangle&t,const Frame29&frame,const op::Linear&linear,const Panel&panel,S u,S eta,const op::Branches&branches){const auto&lo=branches.lo;const auto&hi=branches.hi;S left=S(1)-S(eta.n),right(eta.n),z(0,eta.e),db=hi.slope-lo.slope,dg=(hi.intercept+hi.correction)-(lo.intercept+lo.correction),width=db*u+dg;S beta=left*lo.slope+right*hi.slope+z*db,gamma=left*(lo.intercept+lo.correction)+right*(hi.intercept+hi.correction)+z*dg;auto c=op::affine_coordinate(linear,u,beta,gamma);auto out=coordinates(t,frame.original,c.s,c.t,&c.angular);out.measure=frame.determinant*width*(panel.u1-panel.u0)*(panel.eta1-panel.eta0);return out;}
}
#endif
