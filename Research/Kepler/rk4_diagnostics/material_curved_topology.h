#ifndef MATERIAL_CURVED_TOPOLOGY_H
#define MATERIAL_CURVED_TOPOLOGY_H
#include "material_curved_partition.h"
namespace curved_topology {
using R=long double;using I=curved_bounds::I;using P=curved_roots::P;using W=curved_roots::Window;using Pt=std::array<R,2>;
static I hull(const P &p,R a,R b){auto c=curved_roots::bernstein(p,a,b);I h=c[0];for(int j=1;j<=p.degree;j++){h.lo=std::min(h.lo,c[j].lo);h.hi=std::max(h.hi,c[j].hi);}return h;}
static P derivative(const P &p){P d(std::max(0,p.degree-1));for(int i=1;i<=p.degree;i++)d.c[i-1]=p.c[i]*I(i);return d;}
static bool cone(const curved::Record &r,int chart,W &range){
    std::vector<Pt> vertex;for(int j=0;j<3;j++)vertex.push_back({r.triangle.geometry.x[j][0],r.triangle.geometry.x[j][1]});std::vector<Pt> clipped;for(size_t j=0;j<vertex.size();j++){Pt a=vertex[j],b=vertex[(j+1)%vertex.size()];bool ia=chart*a[0]>=0,ib=chart*b[0]>=0;if(ia)clipped.push_back(a);if(ia!=ib){R t=-a[0]/(b[0]-a[0]);clipped.push_back({0,a[1]+t*(b[1]-a[1])});}}
    if(clipped.size()<3)return false;R lo=2,hi=-2;for(auto a:clipped){R radius=hypotl(a[0],a[1]);if(radius==0)continue;R u=chart*a[1]/(radius+chart*a[0]);lo=std::min(lo,u);hi=std::max(hi,u);}if(lo>=hi)return false;range={std::max(-1.L,lo-128*LDBL_EPSILON),std::min(1.L,hi+128*LDBL_EPSILON)};return true;
}
struct Curve {P p;int source=-1,left=-1,right=-1;};
struct Root {W interval;int curve;bool simple;};
struct Event {std::vector<std::pair<int,W>> unresolved;uint64_t polynomials=0,root_panels=0;bool stable=true;};
static Event events(const curved::Record &r,R low,R high,R eps,const std::vector<curved_partition::Site> &site,const std::vector<R> &weights,R position_error,const curved_compiled::Map *compiled=nullptr){
    Event result;int own=-1;for(size_t i=0;i<site.size();i++)if(site[i].id==r.cell)own=i;need(own>=0,"event guard source site");
    for(int chart:{1,-1}){
        W domain;if(!cone(r,chart,domain))continue;auto ray=curved_roots::ray(r,I(low,high),chart,eps,compiled);bool possible=true;for(int k=0;k<3;k++)possible=possible && hull(ray.bary[k],domain.lo,domain.hi).hi>=0;if(!possible)continue;
        std::vector<int> candidate={own};for(int j=0;j<(int)site.size();j++)if(j!=own){auto p=curved_partition::plane(ray,site[own],site[j],weights[own],weights[j],position_error);if(hull(p,domain.lo,domain.hi).hi>=0)candidate.push_back(j);}
        std::vector<Curve> curves;for(int k=0;k<3;k++)curves.push_back({ray.bary[k],k,-1,-1});for(size_t i=0;i<candidate.size();i++)for(size_t j=i+1;j<candidate.size();j++){int a=candidate[i],b=candidate[j];curves.push_back({curved_partition::plane(ray,site[a],site[b],weights[a],weights[b],position_error),-1,a,b});}
        std::vector<Root> roots;std::vector<W> bad;R tolerance=std::min((domain.hi-domain.lo)*.02L,std::max(1e-11L,.05L*(high-low)/std::max(low,eps)));
        for(size_t ip=0;ip<curves.size();ip++){
            const auto &curve=curves[ip];auto root=curved_roots::isolate(curve.p,domain.lo,domain.hi,tolerance,36,2048);result.polynomials++;result.root_panels+=root.panels;
            for(auto w:root.window){
                bool relevant=true;for(int k=0;k<3;k++)if(k!=curve.source && hull(ray.bary[k],w.lo,w.hi).hi<0)relevant=false;if(!relevant)continue;
                if(curve.left>=0)for(int k:candidate)if(k!=curve.left && k!=curve.right){auto f=curved_partition::plane(ray,site[curve.left],site[k],weights[curve.left],weights[k],position_error);if(hull(f,w.lo,w.hi).lo>0)relevant=false;}if(!relevant)continue;
                /* A reference source edge on x=0 is exactly the common chart
                 * boundary. Its persistent endpoint root has no interior
                 * topology event. This uses original edge coordinates. */
                if(curve.source>=0){int a=(curve.source+1)%3,b=(curve.source+2)%3;if(r.triangle.geometry.x[a][0]==0 && r.triangle.geometry.x[b][0]==0 && (w.lo==domain.lo || w.hi==domain.hi) && (domain.lo==-1 || domain.hi==1))continue;}
                R padding=std::max(1e-13L,(w.hi-w.lo)*.1L);w.lo=std::max(domain.lo,w.lo-padding);w.hi=std::min(domain.hi,w.hi+padding);I d=hull(derivative(curve.p),w.lo,w.hi),a=curved_roots::value(curve.p,I(w.lo)),b=curved_roots::value(curve.p,I(w.hi));bool monotone=d.lo>0 || d.hi<0,opposite=(a.hi<0 && b.lo>0)||(a.lo>0 && b.hi<0);bool simple=root.budget_complete && monotone && opposite;roots.push_back({w,(int)ip,simple});if(!simple)bad.push_back(w);
            }
        }
        for(size_t i=0;i<roots.size();i++)for(size_t j=i+1;j<roots.size();j++)if(roots[i].curve!=roots[j].curve && std::max(roots[i].interval.lo,roots[j].interval.lo)<=std::min(roots[i].interval.hi,roots[j].interval.hi)){bad.push_back(roots[i].interval);bad.push_back(roots[j].interval);}
        auto merged=curved_partition::merge(bad);for(auto w:merged)result.unresolved.push_back({chart,w});if(!merged.empty())result.stable=false;
    }return result;
}
static R density_upper(const curved::Packet &p,R lo,R hi){R a=.5L*(1+tanhl((hi-p.parameters[4])/p.parameters[8])),b=.5L*(1+tanhl((p.parameters[5]-lo)/p.parameters[8]));return curved_bounds::up(p.parameters[6]+p.parameters[9]/(lo+p.parameters[7])*a*b);}
struct Panel {R low,high,K_bound=0,mass_bound=0,angle_bound=0;int depth=0;bool stable=false;uint64_t polynomials=0,root_panels=0;};
static Panel panel(const curved::Record &r,const curved::Packet &p,R low,R high,int depth,const std::vector<curved_partition::Site> &site,const std::vector<R> &weight,const curved_bounds::Enclosure &bounds,const curved_compiled::Map *compiled=nullptr){
    auto e=events(r,low,high,p.parameters[3],site,weight,bounds.native_position_envelope,compiled);Panel out={};out.low=low;out.high=high;out.depth=depth;out.stable=e.stable;out.polynomials=e.polynomials;out.root_panels=e.root_panels;for(auto window:e.unresolved)out.angle_bound+=2*(atanl(window.second.hi)-atanl(window.second.lo));R rho=density_upper(p,low,high),measure=curved_bounds::up((high-low)*(high+low)/2*out.angle_bound),v2=0;for(int k=0;k<2;k++)v2+=curved_bounds::absmax(bounds.velocity[k])*curved_bounds::absmax(bounds.velocity[k]);out.mass_bound=curved_bounds::up(3*measure*rho);out.K_bound=curved_bounds::up(out.mass_bound*v2/2);return out;
}
}
#endif
