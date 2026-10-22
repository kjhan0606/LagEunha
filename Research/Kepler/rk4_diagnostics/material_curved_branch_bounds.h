#ifndef MATERIAL_CURVED_BRANCH_BOUNDS_H
#define MATERIAL_CURVED_BRANCH_BOUNDS_H
#include "material_curved_branches.h"
/* Whole-radius-panel inclusion. Root intervals are not integration cuts.
 * Direct endpoint evaluation avoids the dependency growth of repeatedly
 * subdividing uncertain Bernstein coefficients of two nearby roots. */
namespace curved_branch_bounds {
using R=long double;using I=curved_bounds::I;using P=curved_roots::P;using W=curved_roots::Window;
static curved_roots::Ray representative_ray(const curved::Record &r,const curved_branches::FixedMap &m,I radius,int chart,R eps){
    auto fixed=m.interval;for(int k=0;k<2;k++){for(int a=0;a<10;a++)fixed.position[k].a[a]=I(m.position[k][a]);for(int j=0;j<2;j++)fixed.position_B[k][j]=I(m.B[k][j]);}
    auto ray=curved_roots::ray(r,radius,chart,eps,&fixed);P X(2),Y(2);X.c[0]=radius*I(chart);X.c[2]=-radius*I(chart);Y.c[1]=radius*I(2*chart);
    for(int k=0;k<3;k++)ray.bary[k]=curved_roots::sub(curved_roots::add(curved_roots::scale(X,I(m.source_n[k][0])),curved_roots::scale(Y,I(m.source_n[k][1]))),curved_roots::scale(ray.denominator,I(m.source_c[k])));return ray;
}
static P unite(P a,const P &b){a.degree=std::max(a.degree,b.degree);for(int k=0;k<=a.degree;k++){a.c[k].lo=std::min(a.c[k].lo,b.c[k].lo);a.c[k].hi=std::max(a.c[k].hi,b.c[k].hi);}return a;}
struct Tube {W interval={-1,1};R min_dp=0,left_margin=0,right_margin=0;int queries=0;bool complete=false;};
static Tube axis_tube(const curved_branches::FixedMap &m,const curved_branches::Root &root,R lo,R hi){Tube out;int k=root.curve.source;if(k<0 || m.source_axis[k]<0 || !(lo>fabsl(m.source_axis_coordinate[k])))return out;R d=fabsl(m.source_axis_coordinate[k]);I r(lo,hi),minus=curved_branches::tight_add(r,I(-d)),plus=curved_branches::tight_add(r,I(d)),product=minus*plus;product.lo=std::max(0.L,product.lo);I radial=curved_bounds::root(product),u;
    if(m.source_axis[k]==0){u=curved_bounds::root(minus/plus);if(root.u<0)u=-u;out.min_dp=curved_bounds::down(2*fabsl(m.source_n[k][0])*radial.lo);}else{u=I(root.chart*m.source_axis_coordinate[k])/(r+radial);out.min_dp=curved_bounds::down(2*fabsl(m.source_n[k][1])*radial.lo);}out.interval={std::max(-1.L,u.lo),std::min(1.L,u.hi)};out.complete=out.min_dp>0;out.queries=1;return out;}
static Tube certify(const P &p,const curved_branches::Root &root,R lo,R hi,W domain={-1,1}){
    Tube out;auto dp=curved_branches::derivative(p);R pad=std::max({1e-12L,4*fabsl(root.du_dr)*(hi-lo),4*(root.band_hi-root.band_lo)});
    for(int attempt=0;attempt<12;attempt++){
        W w={std::max(domain.lo,root.u-pad),std::min(domain.hi,root.u+pad)};out.interval=w;out.queries++;I a=curved_roots::value(p,I(w.lo)),b=curved_roots::value(p,I(w.hi)),d=curved_roots::value(dp,I(w.lo,w.hi));auto bern=curved_branches::hull(dp,w);d.lo=std::max(d.lo,bern.lo);d.hi=std::min(d.hi,bern.hi);need(d.lo<=d.hi,"two valid derivative enclosures intersect");bool opposite=(a.hi<0 && b.lo>0)||(a.lo>0 && b.hi<0);bool monotone=d.lo>0 || d.hi<0;
        if(opposite && monotone){out.complete=true;out.min_dp=std::min(fabsl(d.lo),fabsl(d.hi));out.left_margin=std::min(fabsl(a.lo),fabsl(a.hi));out.right_margin=std::min(fabsl(b.lo),fabsl(b.hi));return out;}if(w.lo==domain.lo && w.hi==domain.hi)break;pad*=2;
    }return out;
}
/* Analytic conversion component only: full cubic AND B, exact source
 * inverse, and exact power-plane scalar arithmetic. Numerical ray/root,
 * field evaluation and gluing remain distinct terms in the caller. */
struct Conversion {R analytic_numerator=0,native_numerator=0;};
struct SourceFold {bool complete=false;I critical;R C=0,numerical=0;};
static SourceFold source_fold(const curved::Record &r,const curved_branches::FixedMap &m,int k,R maximum){
    SourceFold out;if(k<0 || m.source_axis[k]>=0)return out;for(int d=0;d<2;d++)if(m.source_exact_n[k][d].lo!=m.source_exact_n[k][d].hi)return out;if(m.source_exact_c[k].lo!=m.source_exact_c[k].hi)return out;
    int a=(k+1)%3,b=(k+2)%3;I ax(r.triangle.geometry.x[a][0]),ay(r.triangle.geometry.x[a][1]),dx=curved_branches::tight_add(I(r.triangle.geometry.x[b][0]),-ax),dy=curved_branches::tight_add(I(r.triangle.geometry.x[b][1]),-ay),length2=curved_branches::tight_add(curved_branches::tight_mul(dx,dx),curved_branches::tight_mul(dy,dy)),dot=curved_branches::tight_add(curved_branches::tight_mul(ax,dx),curved_branches::tight_mul(ay,dy)),t=curved_branches::tight_div(-dot,length2);if(t.lo<0 || t.hi>1)return out;I x=curved_branches::tight_add(ax,curved_branches::tight_mul(t,dx)),y=curved_branches::tight_add(ay,curved_branches::tight_mul(t,dy)),r2=curved_branches::tight_add(curved_branches::tight_mul(x,x),curved_branches::tight_mul(y,y));if(r2.lo<=0)return out;I radius=curved_bounds::root(r2),norm=curved_bounds::root(length2);R nx=fabsl(m.source_n[k][0]),ny=fabsl(m.source_n[k][1]),c=fabsl(m.source_c[k]),scale=2*maximum*(nx+ny)+2*c,nu=256*LDBL_EPSILON/2,dp=2*(nx*maximum+c)+2*ny*maximum;out.numerical=curved_bounds::up(nu/(1-nu)*scale+256*LDBL_EPSILON*dp);R margin=curved_bounds::up(2*out.numerical/norm.lo);out.critical=I(curved_bounds::down(radius.lo-margin),curved_bounds::up(radius.hi+margin));if(out.critical.lo<=0)return out;out.C=(norm*curved_bounds::root(I(2)*I(out.critical.lo))).lo;out.complete=out.C>0;return out;
}
/* Conservative arithmetic model, separate from analytic conversion:
 * <16384 long-double operations in the degree-six ray construction, normal
 * round-to-nearest arithmetic, and native powl within8 ulps. The L1 norm is
 * taken BEFORE cancellation, not from the cancelled final coefficients.
 * This is an engineering envelope, not a formal libm proof. */
static R arithmetic(const curved_branches::FixedMap &m,const curved_branches::Curve &curve,R radius,R eps,const std::vector<curved_partition::Site> &site,const std::vector<R> &weight){
    I r(radius),norm(0);int operations=16384;
    if(curve.source>=0){int k=curve.source;norm=I(2)*r*(I(fabsl(m.source_n[k][0]))+I(fabsl(m.source_n[k][1])))+I(2*fabsl(m.source_c[k]));operations=256;}
    else{I X=(I(2)*r+I(2*fabsl(m.anchor[0])))/I(m.h),Y=(I(2)*r+I(2*fabsl(m.anchor[1])))/I(m.h),root=curved_bounds::root(curved_bounds::root(curved_bounds::square(I(eps)))),omega=I(1)/(root*root*root),position[2]={};for(int k=0;k<2;k++){for(int a=0;a<10;a++){auto e=curved_compiled::exponents(a);position[k]=position[k]+I(fabsl(m.position[k][a]))*curved_roots::power(X,e[0])*curved_roots::power(Y,e[1])*curved_roots::power(I(2),3-e[0]-e[1]);}for(int j=0;j<2;j++)position[k]=position[k]+I(8)*r*omega*I(fabsl(m.B[k][j]));}const auto &a=site[curve.i],&b=site[curve.j];I nx=I(fabsl(a.x))+I(fabsl(b.x)),ny=I(fabsl(a.y))+I(fabsl(b.y)),c=(curved_bounds::square(I(a.x))+curved_bounds::square(I(a.y))+curved_bounds::square(I(b.x))+curved_bounds::square(I(b.y))+I(fabsl(weight[curve.i]))+I(fabsl(weight[curve.j])))/I(2);norm=nx*position[0]+ny*position[1]+I(8)*c;}
    I nu=I(operations)*I(LDBL_EPSILON/2),gamma=nu/(I(1)-nu);return curved_bounds::absmax(gamma*norm);
}
static Conversion conversion(const curved::Record &r,const curved_branches::FixedMap &m,const curved_branches::Curve &curve,I radius,W tube,int chart,R eps,const std::vector<curved_partition::Site> &site,const std::vector<R> &weight,R native_position_error){
    I u(tube.lo,tube.hi),D=I(1)+curved_bounds::square(u),X=I(chart)*radius*(I(1)-curved_bounds::square(u))/D,Y=I(2*chart)*radius*u/D;X.lo=std::max(X.lo,-radius.hi);X.hi=std::min(X.hi,radius.hi);Y.lo=std::max(Y.lo,-radius.hi);Y.hi=std::min(Y.hi,radius.hi);I local[2]={(X-I(m.anchor[0]))/I(m.h),(Y-I(m.anchor[1]))/I(m.h)},z=curved_bounds::square(radius)+curved_bounds::square(I(eps)),q=curved_bounds::root(curved_bounds::root(z)),omega=I(1)/(q*q*q),v0[2]={-Y*omega,X*omega};Conversion out;
    if(curve.source>=0){
        int k=curve.source;I dn0=curved_branches::tight_add(m.source_exact_n[k][0],I(-m.source_n[k][0])),dn1=curved_branches::tight_add(m.source_exact_n[k][1],I(-m.source_n[k][1])),dc=curved_branches::tight_add(m.source_exact_c[k],I(-m.source_c[k]));I delta=dn0*X+dn1*Y-dc;out.analytic_numerator=curved_bounds::absmax(delta*D);return out;
    }
    I position[2]={},delta[2]={};for(int k=0;k<2;k++){for(int a=0;a<10;a++){auto e=curved_compiled::exponents(a);I basis=curved_roots::power(local[0],e[0])*curved_roots::power(local[1],e[1]);position[k]=position[k]+I(m.position[k][a])*basis;delta[k]=delta[k]+(m.interval.position[k].a[a]-I(m.position[k][a]))*basis;}for(int j=0;j<2;j++){position[k]=position[k]+I(m.B[k][j])*v0[j];delta[k]=delta[k]+(m.interval.position_B[k][j]-I(m.B[k][j]))*v0[j];}}
    const auto &a=site[curve.i],&b=site[curve.j];I n[2]={I(b.x)-I(a.x),I(b.y)-I(a.y)},c=(curved_bounds::square(I(b.x))+curved_bounds::square(I(b.y))-curved_bounds::square(I(a.x))-curved_bounds::square(I(a.y))+I(weight[curve.i])-I(weight[curve.j]))/I(2);R scalar_n[2]={b.x-a.x,b.y-a.y},scalar_c=(b.x*b.x+b.y*b.y-a.x*a.x-a.y*a.y+weight[curve.i]-weight[curve.j])/2;I df=-(c-I(scalar_c));
    /* n_exact*delta includes the dn*delta cross term. */
    for(int k=0;k<2;k++)df=df+n[k]*delta[k]+(n[k]-I(scalar_n[k]))*position[k];I D3=D*D*D;out.analytic_numerator=curved_bounds::absmax(df*D3);out.native_numerator=curved_bounds::absmax(D3)*native_position_error*(curved_bounds::absmax(n[0])+curved_bounds::absmax(n[1]));return out;
}
}
#endif
