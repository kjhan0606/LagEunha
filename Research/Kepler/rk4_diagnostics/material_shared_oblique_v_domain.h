#ifndef MATERIAL_SHARED_OBLIQUE_V_DOMAIN_H
#define MATERIAL_SHARED_OBLIQUE_V_DOMAIN_H
// Included inside oblique_point after its original chart types.
struct VDomainProof {
 uint64_t slot,source,target,branch,panel,axis,stage,complete;
 S u,alpha,radial,radial_derivative,lambda[3],lambda_derivative[3],X[2],denominator[3],parameter[3],parameter_determinant[3];
 Q beta,radius,apex_radius,parameter_radius[3],parameter_lipschitz[3];
};
static thread_local std::vector<VDomainProof>*v_domain_observer=nullptr;
static thread_local uint64_t v_domain_panel=UINT64_MAX,v_domain_axis=UINT64_MAX;
struct VDomainTrace {VDomainProof p={};~VDomainTrace(){if(v_domain_observer)v_domain_observer->push_back(p);}};
static Q geometric_v_prior(const sd::Compiled&m,const Frame29&f,const Linear&l,int index,S u,S root,Q beta){
 VDomainTrace trace;auto&z=trace.p;z.slot=f.slot;z.source=f.source;z.target=f.target;z.branch=index;z.panel=v_domain_panel;z.axis=v_domain_axis;z.u=u;z.beta=beta;z.alpha=S(root.n)-S(beta)*S(u.n);z.stage=1;
 z.radial=(l.s[0]+l.s[1]*S(beta))*u+l.s[1]*z.alpha;z.radial_derivative=l.s[1];Q floor=sf::lower(z.radial);if(!(floor>0))throw sd::DomainFailure("geometric oblique predictor radial domain",z.radial,Q(0));Q radius=Q(1),change=sf::mag(z.radial_derivative);if(change>0)radius=std::min(radius,mf::FP<Q>::down(floor/sf::mul(Q(2),change)));z.radius=radius;
 S gu[3],gv[3];for(int a=0;a<3;a++)for(int k=0;k<2;k++){gu[a]=gu[a]+m.base.grad[a][k]*f.transverse[k];gv[a]=gv[a]+m.base.grad[a][k]*f.direction[k];}
 auto combination=[&](const S*weight){S A,du,dv;for(int a=0;a<3;a++){if(a==f.original.apex)A=weight[a];du=du+weight[a]*gu[a];dv=dv+weight[a]*gv[a];}return A+(du+S(beta)*dv)*u+dv*z.alpha;};
 for(int a=0;a<3;a++){S weight[3];weight[a]=S(1);z.lambda[a]=combination(weight);z.lambda_derivative[a]=gv[a];}
 for(int k=0;k<2;k++)z.X[k]=f.apex[k]+(f.transverse[k]+S(beta)*f.direction[k])*u+f.direction[k]*z.alpha;
 const auto&t=m.base.state->triangle[m.base.source];z.stage=2;z.apex_radius=target::binomial_clearance(t.vertex[f.original.apex],z.X[0],z.X[1],f.direction[0],f.direction[1],m.base.eps);radius=std::min(radius,z.apex_radius);z.radius=radius;
 for(int e=0;e<3;e++){z.stage=3+e;int a=t.endpoint[e][0],b=t.endpoint[e][1];S weight[3];weight[a]=weight[b]=S(1);S D=combination(weight);z.denominator[e]=D;floor=sf::lower(D);if(!(floor>0))throw sd::DomainFailure("geometric oblique predictor edge denominator",D,Q(0));S derivative=gv[a]+gv[b];change=sf::mag(derivative);if(change>0)radius=std::min(radius,mf::FP<Q>::down(floor/sf::mul(Q(2),change)));
  Q center=(z.lambda[b]/D).n;S rw[3];rw[a]=S(-center);rw[b]=S(1)-S(center);S parameter=S(center)+combination(rw)/D;z.parameter[e]=parameter;const auto&edge=m.base.edge[e];S X=edge.A[0]+parameter*edge.Dref[0],Y=edge.A[1]+parameter*edge.Dref[1];const double*anchor=parameter.n>Q(.5L)?t.vertex[t.endpoint[e][1]]:edge.native.A;
  try{z.parameter_radius[e]=target::binomial_clearance(anchor,X,Y,edge.Dref[0],edge.Dref[1],m.base.eps);}catch(sd::DomainFailure&failure){failure.edge=e;failure.edge_parameter=parameter;throw;}
  // Finite Mobius difference: delta-v*(gb*D-N*(ga+gb))/(D*D_new).
  // D_new>=D_floor/2 gives the factor2 and no parameter-radius scan.
  S dw[3];dw[a]=gv[b];dw[b]=-gv[a];z.parameter_determinant[e]=combination(dw);Q squared=mf::FP<Q>::down(floor*floor);z.parameter_lipschitz[e]=sf::quotient(sf::mul(Q(2),sf::mag(z.parameter_determinant[e])),squared);if(z.parameter_lipschitz[e]>0)radius=std::min(radius,mf::FP<Q>::down(z.parameter_radius[e]/z.parameter_lipschitz[e]));z.radius=radius;
 }
 need(radius>0,"positive geometric oblique correction neighborhood");z.complete=1;z.stage=6;return radius;
}
#endif
