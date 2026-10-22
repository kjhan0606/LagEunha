#ifndef MATERIAL_CURVED_PARTITION_H
#define MATERIAL_CURVED_PARTITION_H
#define MATERIAL_CURVED_PARTITION_RAY_CACHE
#include "material_curved_roots.h"
namespace curved_partition {
using R=long double;using I=curved_bounds::I;using W=curved_roots::Window;
struct Site {uint64_t id;R x,y;};
struct Fragment {int target,chart;R left,right;};
struct Partition {std::vector<Fragment> fragment;std::vector<std::pair<int,W>> source_unknown,target_unknown;curved_roots::Ray chart_ray[2];uint64_t root_panels=0,pair_polynomials=0;bool complete=true;};
static std::vector<Site> read_sites(const char *path){FILE *f=fopen(path,"r");need(f,"advected native sites");char header[512];need(fgets(header,sizeof(header),f),"site CSV header");std::vector<Site> v;unsigned long long id;double X,Y,x,y,vx,vy;while(fscanf(f,"%llu,%lf,%lf,%lf,%lf,%lf,%lf\n",&id,&X,&Y,&x,&y,&vx,&vy)==7)v.push_back({id,x,y});need(feof(f) && v.size()==24,"all actual sites");fclose(f);return v;}
static curved_roots::P plane(const curved_roots::Ray &ray,const Site &a,const Site &b,R wa,R wb,R uncertainty){I ax=I(b.x)-I(a.x),ay=I(b.y)-I(a.y),c=(I(b.x)*I(b.x)+I(b.y)*I(b.y)-I(a.x)*I(a.x)-I(a.y)*I(a.y)+I(wa)-I(wb))/I(2);return curved_roots::plane(ray,ax,ay,c,uncertainty);}
static std::vector<W> merge(std::vector<W> input){std::sort(input.begin(),input.end(),[](W a,W b){return a.lo<b.lo;});std::vector<W> v;for(auto w:input){if(!v.empty() && w.lo<=v.back().hi)v.back().hi=std::max(v.back().hi,w.hi);else v.push_back(w);}return v;}
static I hull(const curved_roots::P &p,W w){auto b=curved_roots::bernstein(p,w.lo,w.hi);I out=b[0];for(int i=1;i<=p.degree;i++){out.lo=std::min(out.lo,b[i].lo);out.hi=std::max(out.hi,b[i].hi);}return out;}
static Partition partition(const curved::Record &r,R radius,R eps,const std::vector<Site> &sites,const std::vector<R> &weight,R uncertainty=0,R tolerance=1e-12L,const curved_compiled::Map *compiled=nullptr,const std::vector<int> *whole_source_candidates=nullptr){
    need(radius>0 && sites.size()==weight.size(),"finite radial partition input");int own=-1;for(size_t j=0;j<sites.size();j++)if(sites[j].id==r.cell)own=j;need(own>=0,"source site in target partition");Partition out;
    for(int chart:{1,-1}){
        auto ray=curved_roots::ray(r,I(radius),chart,eps,compiled);out.chart_ray[chart==1?0:1]=ray;std::vector<W> unknown;
        for(int k=0;k<3;k++){auto roots=curved_roots::isolate(ray.bary[k],-1,1,tolerance,46,32768);out.root_panels+=roots.panels;out.complete=out.complete && roots.budget_complete;unknown.insert(unknown.end(),roots.window.begin(),roots.window.end());}auto source_bands=merge(unknown);
        std::vector<W> source;R left=-1;auto inside=[&](R a,R b){if(a>=b)return;bool yes=true;for(int k=0;k<3;k++){I v=curved_roots::value(ray.bary[k],I((a+b)/2));need(v.lo>0 || v.hi<0,"source sign separated from retained bands");yes=yes && v.lo>0;}if(yes)source.push_back({a,b});};for(auto w:source_bands){inside(left,w.lo);bool possible=true;for(int k=0;k<3;k++)possible=possible && hull(ray.bary[k],w).hi>=0;if(possible)out.source_unknown.push_back({chart,w});left=w.hi;}inside(left,1);
        for(auto segment:source){
            std::vector<int> candidate={own};for(int j=0;j<(int)sites.size();j++)if(j!=own && (!whole_source_candidates || std::find(whole_source_candidates->begin(),whole_source_candidates->end(),j)!=whole_source_candidates->end())){auto f=plane(ray,sites[own],sites[j],weight[own],weight[j],uncertainty);if(curved_roots::sign(curved_roots::bernstein(f,segment.lo,segment.hi),f.degree)!=-1)candidate.push_back(j);}
            std::vector<W> target_bands;for(size_t i=0;i<candidate.size();i++)for(size_t j=i+1;j<candidate.size();j++){int ci=candidate[i],cj=candidate[j];auto f=plane(ray,sites[ci],sites[cj],weight[ci],weight[cj],uncertainty);auto roots=curved_roots::isolate(f,segment.lo,segment.hi,tolerance,46,32768);out.pair_polynomials++;out.root_panels+=roots.panels;out.complete=out.complete && roots.budget_complete;for(auto w:roots.window){bool possible=true;for(int k:candidate)if(k!=ci && k!=cj){auto a=plane(ray,sites[ci],sites[k],weight[ci],weight[k],uncertainty),b=plane(ray,sites[cj],sites[k],weight[cj],weight[k],uncertainty);if(hull(a,w).lo>0 && hull(b,w).lo>0){possible=false;break;}}if(possible)target_bands.push_back(w);}}
            auto bands=merge(target_bands);R left=segment.lo;auto label=[&](R a,R b){if(a>=b)return;R mid=(a+b)/2,D=1+mid*mid;I x=curved_roots::value(ray.position[0],I(mid))/I(D*D*D),y=curved_roots::value(ray.position[1],I(mid))/I(D*D*D);int winner=0;R best=1e99L;for(int j:candidate){R score=sites[j].x*sites[j].x+sites[j].y*sites[j].y-sites[j].x*(x.lo+x.hi)-sites[j].y*(y.lo+y.hi)-weight[j];if(score<best){best=score;winner=j;}}
                for(int j:candidate)if(j!=winner){auto f=plane(ray,sites[winner],sites[j],weight[winner],weight[j],uncertainty);need(curved_roots::value(f,I(mid)).hi<0,"target label certified outside bands");}out.fragment.push_back({winner,chart,a,b});};
            for(auto w:bands){label(left,w.lo);out.target_unknown.push_back({chart,w});left=w.hi;}label(left,segment.hi);
        }
    }return out;
}
}
#endif
