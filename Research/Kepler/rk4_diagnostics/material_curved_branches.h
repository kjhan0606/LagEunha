#ifndef MATERIAL_CURVED_BRANCHES_H
#define MATERIAL_CURVED_BRANCHES_H
#include "material_curved_events.h"
namespace curved_branches {
using R=long double;using I=curved_bounds::I;using P=curved::Poly;using IP=curved_roots::P;using W=curved_roots::Window;
static R mid(I a){return a.lo+(a.hi-a.lo)/2;}
static int sign(R a){return (a>0)-(a<0);}
static P pow_poly(P p,int n){P out=curved::constant(1);for(int j=0;j<n;j++)out=curved::mul(out,p);return out;}
/* Exactness checks are local to binary source-edge predicates. No global
 * interval rule or physical map coefficient is changed. TwoSum verifies an
 * exact addition; FMA verifies a product when residual underflow is excluded. */
static I tight_add(I a,I b){if(a.lo==a.hi && b.lo==b.hi){R s=a.lo+b.lo;if(std::isfinite(s)){R z=s-a.lo,e=(a.lo-(s-z))+(b.lo-z);if(e==0)return I(s);}}return a+b;}
static I tight_mul(I a,I b){if(a.lo==a.hi && b.lo==b.hi){R p=a.lo*b.lo;R minimum=16*std::numeric_limits<R>::min()/(LDBL_EPSILON*LDBL_EPSILON);if((a.lo==0 || b.lo==0) || (std::isfinite(p) && fabsl(p)>=minimum && fmal(a.lo,b.lo,-p)==0))return I(p);}return a*b;}
static I tight_div(I a,I b){if(a.lo==a.hi && b.lo==b.hi && b.lo!=0){R q=a.lo/b.lo,minimum=16*std::numeric_limits<R>::min()/(LDBL_EPSILON*LDBL_EPSILON);if(a.lo==0 || (std::isfinite(q) && fabsl(a.lo)>=minimum && fmal(q,b.lo,-a.lo)==0))return I(q);}return a/b;}
struct FixedMap {curved_compiled::Map interval;R position[2][10],B[2][2],anchor[2],h,source_n[3][2],source_c[3],source_axis_coordinate[3];I source_exact_n[3][2],source_exact_c[3];int source_axis[3];};
static FixedMap fixed_map(const curved::Record &r,R eps){
    FixedMap m;m.interval=curved_compiled::compile(r,eps);m.h=m.interval.h;for(int k=0;k<2;k++){m.anchor[k]=m.interval.anchor[k];for(int a=0;a<10;a++)m.position[k][a]=mid(m.interval.position[k].a[a]);for(int j=0;j<2;j++)m.B[k][j]=mid(m.interval.position_B[k][j]);}
    for(int k=0;k<3;k++){int a=(k+1)%3,b=(k+2)%3;I ax(r.triangle.geometry.x[a][0]),ay(r.triangle.geometry.x[a][1]),bx(r.triangle.geometry.x[b][0]),by(r.triangle.geometry.x[b][1]);m.source_exact_n[k][0]=tight_add(ay,-by);m.source_exact_n[k][1]=tight_add(bx,-ax);m.source_exact_c[k]=tight_add(tight_mul(m.source_exact_n[k][0],ax),tight_mul(m.source_exact_n[k][1],ay));bool exact=m.source_exact_c[k].lo==m.source_exact_c[k].hi;for(int d=0;d<2;d++){m.source_n[k][d]=mid(m.source_exact_n[k][d]);exact=exact && m.source_exact_n[k][d].lo==m.source_exact_n[k][d].hi;}m.source_c[k]=mid(m.source_exact_c[k]);m.source_axis[k]=-1;if(exact && ax.lo==bx.lo){m.source_axis[k]=0;m.source_axis_coordinate[k]=ax.lo;}else if(exact && ay.lo==by.lo){m.source_axis[k]=1;m.source_axis_coordinate[k]=ay.lo;}}
    return m;
}
static void source_ray(const FixedMap &m,I radius,int chart,curved_roots::Ray &ray){IP X(2),Y(2);X.c[0]=radius*I(chart);X.c[2]=-radius*I(chart);Y.c[1]=radius*I(2*chart);for(int k=0;k<3;k++)ray.bary[k]=curved_roots::sub(curved_roots::add(curved_roots::scale(X,m.source_exact_n[k][0]),curved_roots::scale(Y,m.source_exact_n[k][1])),curved_roots::scale(ray.denominator,m.source_exact_c[k]));}
/* Stable circle/axis-line intersections of the EXACT reference geometry.
 * No inverse barycentric normal and no nearly equal squared radii. */
static std::vector<R> axis_roots(const FixedMap &m,int edge,R radius,int chart){std::vector<R> out;R coord=m.source_axis_coordinate[edge],d=fabsl(coord);if(radius<d || !(radius>0))return out;if(m.source_axis[edge]==0){if(coord==0){out={-1,1};return out;}if(chart*coord<0)return out;R s=sqrtl((radius-d)/(radius+d));out.push_back(-s);if(s!=0)out.push_back(s);}else{R s=sqrtl((radius-d)*(radius+d));out.push_back(chart*coord/(radius+s));}return out;}
struct Ray {P D,position[2],bary[3];};
static Ray ray(const FixedMap &m,R radius,int chart,R eps){
    Ray out;P D(2);D.c[0]=curved::Dual(1);D.c[2]=curved::Dual(1);out.D=D;P D2=curved::mul(D,D),D3=curved::mul(D2,D),X(2),Y(2);curved::Dual rad(radius,1);X.c[0]=rad*curved::Dual(chart);X.c[2]=rad*curved::Dual(-chart);Y.c[1]=rad*curved::Dual(2*chart);P x=curved::scale(curved::sub(X,curved::scale(D,m.anchor[0])),1/m.h),y=curved::scale(curved::sub(Y,curved::scale(D,m.anchor[1])),1/m.h);R z=radius*radius+eps*eps,o=powl(z,-.75L);curved::Dual omega(o,-1.5L*radius*o/z);P initial[2]={curved::scale_dual(curved::mul(curved::scale(Y,-1),D2),omega),curved::scale_dual(curved::mul(X,D2),omega)};
    for(int k=0;k<2;k++){P p=curved::constant(0);for(int a=0;a<10;a++)if(m.position[k][a]!=0){auto e=curved_compiled::exponents(a);p=curved::add(p,curved::scale(curved::mul(curved::mul(pow_poly(x,e[0]),pow_poly(y,e[1])),pow_poly(D,3-e[0]-e[1])),m.position[k][a]));}for(int j=0;j<2;j++)if(m.B[k][j]!=0)p=curved::add(p,curved::scale(initial[j],m.B[k][j]));out.position[k]=p;}
    for(int k=0;k<3;k++)out.bary[k]=curved::sub(curved::add(curved::scale(X,m.source_n[k][0]),curved::scale(Y,m.source_n[k][1])),curved::scale(D,m.source_c[k]));return out;
}
struct Curve {int id,source=-1,i=-1,j=-1;};
static P polynomial(const Ray &r,const Curve &c,const std::vector<curved_partition::Site> &site,const std::vector<R> &weight){if(c.source>=0)return r.bary[c.source];const auto &a=site[c.i],&b=site[c.j];R constant=(b.x*b.x+b.y*b.y-a.x*a.x-a.y*a.y+weight[c.i]-weight[c.j])/2;return curved::sub(curved::add(curved::scale(r.position[0],b.x-a.x),curved::scale(r.position[1],b.y-a.y)),curved::scale(pow_poly(r.D,3),constant));}
static IP polynomial(const curved_roots::Ray &r,const Curve &c,const std::vector<curved_partition::Site> &site,const std::vector<R> &weight,R error=0){if(c.source>=0)return r.bary[c.source];return curved_partition::plane(r,site[c.i],site[c.j],weight[c.i],weight[c.j],error);}
static IP derivative(const IP &p){IP q(std::max(0,p.degree-1));for(int a=1;a<=p.degree;a++)q.c[a-1]=p.c[a]*I(a);return q;}
static I hull(const IP &p,W w){return curved_partition::hull(p,w);}
struct Root {int id,chart;Curve curve;R radius,u,du_dr,f_theta,band_lo,band_hi;bool simple=false,exists=false;};
struct Trace {std::vector<Root> root;curved_roots::Ray interval_ray;Ray representative_ray;uint64_t root_panels=0,failures=0;};
static bool solve(const P &p,R low,R high,R prediction,R &u){
    R a=curved::value(p,low).v,b=curved::value(p,high).v;if(a==0){u=low;return true;}if(b==0){u=high;return true;}if(sign(a)==sign(b))return false;u=std::max(low,std::min(high,prediction));
    for(int iteration=0;iteration<48;iteration++){R f=curved::value(p,u).v,df=curved::derivative(p,u);if(f==0)return true;if(sign(f)==sign(a)){low=u;a=f;}else{high=u;b=f;}if(high-low<=64*LDBL_EPSILON*(1+fabsl(u))){u=low+(high-low)/2;return true;}R next=df!=0?u-f/df:low+(high-low)/2;if(!(next>low && next<high))next=low+(high-low)/2;if(next==u){u=low+(high-low)/2;return true;}u=next;}return false;
}
static Trace trace(const curved::Record &record,const FixedMap &map,R radius,int chart,R eps,const std::vector<curved_partition::Site> &site,const std::vector<R> &weight,const std::vector<Curve> &curves,const std::vector<int> &candidate,const std::map<int,Root> *previous=nullptr){
    Trace out;out.interval_ray=curved_roots::ray(record,I(radius),chart,eps,&map.interval);source_ray(map,I(radius),chart,out.interval_ray);out.representative_ray=ray(map,radius,chart,eps);for(auto curve:curves){auto p=polynomial(out.interval_ray,curve,site,weight),dp=derivative(p);auto scalar=polynomial(out.representative_ray,curve,site,weight);curved_roots::Isolation isolated;bool axis=curve.source>=0 && map.source_axis[curve.source]>=0;std::vector<R> exact_axis;if(axis){exact_axis=axis_roots(map,curve.source,radius,chart);for(R u:exact_axis){R e=128*LDBL_EPSILON*(1+fabsl(u));isolated.window.push_back({std::max(-1.L,u-e),std::min(1.L,u+e)});}}else isolated=curved_roots::isolate(p,-1,1,1e-12L,46,32768);out.root_panels+=isolated.panels;out.failures+=!isolated.budget_complete;int ordinal=0;
        for(auto window:isolated.window){
            if(curve.source>=0){int a=(curve.source+1)%3,b=(curve.source+2)%3;if(record.triangle.geometry.x[a][0]==0 && record.triangle.geometry.x[b][0]==0 && (window.lo==-1 || window.hi==1))continue;}
            bool possible=true;for(int k=0;k<3;k++)if(k!=curve.source && hull(out.interval_ray.bary[k],window).hi<0)possible=false;if(!possible)continue;
            if(curve.source<0)for(int k:candidate)if(k!=curve.i && k!=curve.j){auto a=curved_partition::plane(out.interval_ray,site[curve.i],site[k],weight[curve.i],weight[k],0),b=curved_partition::plane(out.interval_ray,site[curve.j],site[k],weight[curve.j],weight[k],0);if(hull(a,window).lo>0 && hull(b,window).lo>0)possible=false;}if(!possible)continue;
            Root root;root.id=(chart==1?0:100000)+curve.id*8+ordinal++;root.chart=chart;root.curve=curve;root.radius=radius;R pad=std::max((window.hi-window.lo)*.2L,128*LDBL_EPSILON);W bracket={std::max(-1.L,window.lo-pad),std::min(1.L,window.hi+pad)};R prediction=(bracket.lo+bracket.hi)/2;if(previous){auto prior=previous->find(root.id);if(prior!=previous->end())prediction=prior->second.u+prior->second.du_dr*(radius-prior->second.radius);}
            if(axis){root.u=exact_axis.front();for(R u:exact_axis)if(fabsl(u-(window.lo+window.hi)/2)<fabsl(root.u-(window.lo+window.hi)/2))root.u=u;root.exists=true;root.simple=radius>fabsl(map.source_axis_coordinate[curve.source]);}else root.exists=solve(scalar,bracket.lo,bracket.hi,prediction,root.u);I d=hull(dp,bracket),a=curved_roots::value(p,I(bracket.lo)),b=curved_roots::value(p,I(bracket.hi));if(!axis)root.simple=(d.lo>0 || d.hi<0) && ((a.hi<0 && b.lo>0)||(a.lo>0 && b.hi<0));
            if(!root.exists){out.failures++;continue;}R du=curved::derivative(scalar,root.u),D=1+root.u*root.u;root.du_dr=du!=0?-curved::value(scalar,root.u).d/du:0;root.f_theta=du/(2*powl(D,curve.source>=0?0:2));W band=bracket;
            if(root.simple && !axis)for(int iteration=0;iteration<12;iteration++){I deriv=hull(dp,band);if(!(deriv.lo>0 || deriv.hi<0))break;R center=band.lo+(band.hi-band.lo)/2;I N=I(center)-curved_roots::value(p,I(center))/deriv;R lo=std::max(band.lo,N.lo),hi=std::min(band.hi,N.hi);if(lo>hi){out.failures++;break;}if(hi-lo>=.9L*(band.hi-band.lo))break;band={lo,hi};}
            R numerical=128*LDBL_EPSILON*(1+fabsl(root.u));root.band_lo=std::min(band.lo,root.u-numerical);root.band_hi=std::max(band.hi,root.u+numerical);out.root.push_back(root);
        }
    }std::sort(out.root.begin(),out.root.end(),[](const Root &a,const Root &b){return a.u<b.u;});return out;
}
static std::vector<Curve> curves(const std::vector<int> &candidate){std::vector<Curve> out;for(int k=0;k<3;k++)out.push_back({k,k,-1,-1});int id=3;for(size_t a=0;a<candidate.size();a++)for(size_t b=a+1;b<candidate.size();b++)out.push_back({id++,-1,candidate[a],candidate[b]});return out;}
}
#endif
