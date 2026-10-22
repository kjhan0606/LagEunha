#include <algorithm>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <chrono>
#include <string>
static void need(bool ok,const char *message){if(!ok){fprintf(stderr,"FAIL: %s\n",message);exit(2);}}
#include "material_curved_packet.h"
#include "material_curved_bounds.h"
using R=long double;using Pt=std::array<double,2>;using Key=std::pair<Pt,Pt>;
struct Edge {size_t triangle;int vertex,component;Pt a,b;};
static R dot(Pt a,Pt b){return (R)a[0]*b[0]+(R)a[1]*b[1];}
static Pt difference(Pt a,Pt b){return {a[0]-b[0],a[1]-b[1]};}
static R point_segment(Pt p,Pt a,Pt b){auto d=difference(b,a);R t=std::max(0.L,std::min(1.L,dot(difference(p,a),d)/dot(d,d)));return hypotl((R)p[0]-a[0]-t*d[0],(R)p[1]-a[1]-t*d[1]);}
static R segment_distance(const Edge &a,const Edge &b){return std::min({point_segment(a.a,b.a,b.b),point_segment(a.b,b.a,b.b),point_segment(b.a,a.a,a.b),point_segment(b.b,a.a,a.b)});}
static curved_bounds::Tri triangle(const curved::Record &r){curved_bounds::Tri t;for(int a=0;a<3;a++)for(int k=0;k<2;k++)t[a][k]=r.triangle.geometry.x[a][k];return t;}
static std::map<Key,std::vector<Edge>> edges(const curved::Packet &p){
    std::map<Key,std::vector<Edge>> result;for(size_t it=0;it<p.record.size();it++){const auto &r=p.record[it];for(int e=0;e<3;e++){Pt a={r.triangle.geometry.x[e][0],r.triangle.geometry.x[e][1]},b={r.triangle.geometry.x[(e+1)%3][0],r.triangle.geometry.x[(e+1)%3][1]};Key key=a<b?Key(a,b):Key(b,a);result[key].push_back({it,e,r.component,a,b});}}return result;
}
struct Boundary {bool complete=true;uint64_t exterior_edges=0,interior_edges=0,nonadjacent_pairs=0,adjacent_pairs=0,loops=0;R min_nonadjacent_margin=1e99L,max_adjacent_factor=0,max_sampled_trace_gap=0;};
static Boundary boundary(const curved::Packet &p,const std::map<int,curved_bounds::Certificate> &cert){
    Boundary out;auto all=edges(p);std::vector<Edge> exterior;std::map<std::pair<int,Pt>,std::vector<int>> incoming,outgoing;
    for(const auto &entry:all){const auto &list=entry.second;need(list.size()==1 || list.size()==2,"manifold reference triangle edge");if(list.size()==2){
        out.interior_edges++;need(list[0].a==list[1].b && list[0].b==list[1].a && list[0].component==list[1].component,"oriented interior reference edge");
        for(int q=0;q<=8;q++){double t=q/8.,X=(1-t)*list[0].a[0]+t*list[0].b[0],Y=(1-t)*list[0].a[1]+t*list[0].b[1];auto a=curved::native_point(p.record[list[0].triangle],X,Y,p.parameters[3]),b=curved::native_point(p.record[list[1].triangle],X,Y,p.parameters[3]);out.max_sampled_trace_gap=std::max(out.max_sampled_trace_gap,hypotl((R)a.position[0]-b.position[0],(R)a.position[1]-b.position[1]));}
    }else{int index=exterior.size();exterior.push_back(list[0]);outgoing[{list[0].component,list[0].a}].push_back(index);incoming[{list[0].component,list[0].b}].push_back(index);}}
    out.exterior_edges=exterior.size();for(const auto &entry:outgoing)need(entry.second.size()==1 && incoming[entry.first].size()==1,"degree two simple reference boundary graph");need(incoming.size()==outgoing.size(),"closed reference boundary graph");
    std::set<int> visited;for(int start=0;start<(int)exterior.size();start++)if(!visited.count(start)){out.loops++;int i=start;do{need(!visited.count(i),"reference boundary has no early cycle");visited.insert(i);i=outgoing.at({exterior[i].component,exterior[i].b})[0];}while(i!=start);}
    need(out.loops==cert.size(),"one reference boundary and no holes per selected component");
    for(size_t i=0;i<exterior.size();i++)for(size_t j=i+1;j<exterior.size();j++){
        const auto &a=exterior[i],&b=exterior[j];bool adjacent=false;Pt shared={},pa={},pb={};for(Pt x:{a.a,a.b})for(Pt y:{b.a,b.b})if(x==y){need(a.component==b.component,"components have disjoint reference boundaries");adjacent=true;shared=x;pa=x==a.a?a.b:a.a;pb=y==b.a?b.b:b.a;}
        if(adjacent){out.adjacent_pairs++;auto da=difference(pa,shared),db=difference(pb,shared);R cosine=dot(da,db)/sqrtl(dot(da,da)*dot(db,db));need(cosine<1,"nondegenerate boundary corner");R ratio=sqrtl(2/(1-cosine)),factor=ratio*cert.at(a.component).max_lipschitz;out.max_adjacent_factor=std::max(out.max_adjacent_factor,factor);out.complete=out.complete && factor<1;}
        else{out.nonadjacent_pairs++;R distance=segment_distance(a,b),margin=distance-cert.at(a.component).max_displacement-cert.at(b.component).max_displacement;out.min_nonadjacent_margin=std::min(out.min_nonadjacent_margin,margin);out.complete=out.complete && margin>0;}
    }
    /* This test deliberately reports native continuity as sampled evidence.
     * Global injectivity is conditional on exact shared-trace gluing. It is
     * NOT promoted to a complete injectivity proof by this routine. */
    return out;
}
int main(int argc,char **argv){
    need(argc==3 || argc==4,"actual curved packet, new output directory, optional preserved identical-map interval ledger");auto started=std::chrono::steady_clock::now();auto p=curved::read(argv[1]);std::map<int,curved_bounds::Certificate> cert;
    FILE *resume=argc==4?fopen(argv[3],"r"):nullptr;if(argc==4){need(resume,"preserved terminal interval ledger");char header[1024];need(fgets(header,sizeof(header),resume),"preserved ledger header");}
    std::string path=std::string(argv[2])+"/interval_bounds.csv";FILE *f=fopen(path.c_str(),"w");need(f,"new interval-bound ledger");fprintf(f,"triangle,cell,component,complete,min_J,max_J,max_distortion,max_displacement,panels,leaves,max_depth,native_position_envelope,native_gradient_envelope,max_abs_s\n");
    curved_bounds::Certificate total;R sample_excess=0;uint64_t samples=0;
    for(size_t it=0;it<p.record.size();it++){
        const auto &r=p.record[it];curved_bounds::Certificate c;
        if(resume){size_t index;unsigned long long cell,panels,leaves;int component,complete;need(fscanf(resume,"%zu,%llu,%d,%d,%Lg,%Lg,%Lg,%Lg,%llu,%llu,%d,%Lg,%Lg,%Lg\n",&index,&cell,&component,&complete,&c.min_J,&c.max_J,&c.max_lipschitz,&c.max_displacement,&panels,&leaves,&c.max_depth,&c.native_position_envelope,&c.native_gradient_envelope,&c.max_abs_s)==14 && index==it && cell==r.cell && component==r.component,"preserved interval ledger identity");c.complete=complete;c.panels=panels;c.leaves=leaves;}
        else curved_bounds::certify(r,triangle(r),p.parameters[3],0,5,.24L,c);
        auto &z=cert[r.component];z.complete=z.complete && c.complete;z.min_J=std::min(z.min_J,c.min_J);z.max_J=std::max(z.max_J,c.max_J);z.max_lipschitz=std::max(z.max_lipschitz,c.max_lipschitz);z.max_displacement=std::max(z.max_displacement,c.max_displacement);z.panels+=c.panels;z.leaves+=c.leaves;z.max_depth=std::max(z.max_depth,c.max_depth);
        total.complete=total.complete && c.complete;total.min_J=std::min(total.min_J,c.min_J);total.max_J=std::max(total.max_J,c.max_J);total.max_lipschitz=std::max(total.max_lipschitz,c.max_lipschitz);total.max_displacement=std::max(total.max_displacement,c.max_displacement);total.panels+=c.panels;total.leaves+=c.leaves;total.max_depth=std::max(total.max_depth,c.max_depth);total.native_position_envelope=std::max(total.native_position_envelope,c.native_position_envelope);total.native_gradient_envelope=std::max(total.native_gradient_envelope,c.native_gradient_envelope);total.max_abs_s=std::max(total.max_abs_s,c.max_abs_s);
        fprintf(f,"%zu,%llu,%d,%d,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%llu,%llu,%d,%.21Lg,%.21Lg,%.21Lg\n",it,(unsigned long long)r.cell,r.component,c.complete,c.min_J,c.max_J,c.max_lipschitz,c.max_displacement,(unsigned long long)c.panels,(unsigned long long)c.leaves,c.max_depth,c.native_position_envelope,c.native_gradient_envelope,c.max_abs_s);
        for(int i=0;i<5;i++)for(int j=0;j<5;j++){
            R a=(i+.5L)/5,b=(1-a)*(j+.5L)/5,c0=1-a-b;double X=c0*r.triangle.geometry.x[0][0]+a*r.triangle.geometry.x[1][0]+b*r.triangle.geometry.x[2][0],Y=c0*r.triangle.geometry.x[0][1]+a*r.triangle.geometry.x[1][1]+b*r.triangle.geometry.x[2][1];curved_bounds::Tri point={std::array<R,2>{X,Y},std::array<R,2>{X,Y},std::array<R,2>{X,Y}};auto e=curved_bounds::enclose(r,point,p.parameters[3]);auto q=curved::native_point(r,X,Y,p.parameters[3]);
            for(int k=0;k<2;k++){R d=(R)q.position[k]-(k?Y:X);sample_excess=std::max({sample_excess,e.displacement[k].lo-d,d-e.displacement[k].hi});}
            for(int k=0;k<4;k++){R d=(R)q.F[k]-(k==0 || k==3?1:0);sample_excess=std::max({sample_excess,e.gradient[k].lo-d,d-e.gradient[k].hi});}samples++;
        }
    }fclose(f);if(resume){need(fgetc(resume)==EOF,"preserved interval ledger exact end");fclose(resume);}auto b=boundary(p,cert);
    auto reflected=p.record[0];for(int a=0;a<12;a++)reflected.map.x[0][a]=-reflected.map.x[0][a];reflected.origin[0]=-reflected.origin[0];curved_bounds::Certificate negative;curved_bounds::certify(reflected,triangle(reflected),p.parameters[3],0,1,.24L,negative);need(!negative.complete && negative.max_J<0,"inverted map rejected");
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    printf("{\"scope\":\"outward interval local analytic map bounds with explicit native evaluation envelope; boundary separation conditional on exact shared-trace gluing; sampled trace continuity is not a full injectivity proof; no remap/Newton promotion\",\"triangles\":%zu,\"components\":%zu,\"local_interval_complete\":%s,\"sampled_enclosure_queries\":%llu,\"sample_enclosure_excess\":%.17Lg,\"interval_min_J\":%.17Lg,\"interval_max_J\":%.17Lg,\"max_distortion_norm_bound\":%.17Lg,\"max_displacement_bound\":%.17Lg,\"native_position_envelope\":%.17Lg,\"native_gradient_envelope\":%.17Lg,\"max_abs_basis_s\":%.17Lg,\"panels\":%llu,\"leaves\":%llu,\"max_depth\":%d,\"boundary_geometric_separation_complete\":%s,\"boundary_loops\":%llu,\"exterior_edges\":%llu,\"interior_edges\":%llu,\"nonadjacent_boundary_pairs\":%llu,\"adjacent_boundary_pairs\":%llu,\"min_nonadjacent_separation_margin\":%.17Lg,\"max_adjacent_corner_factor\":%.17Lg,\"sampled_shared_trace_gap_max\":%.17Lg,\"inverted_map_rejected\":true,\"global_injectivity_proved\":false,\"seconds\":%.17g}\n",p.record.size(),cert.size(),total.complete?"true":"false",(unsigned long long)samples,sample_excess,total.min_J,total.max_J,total.max_lipschitz,total.max_displacement,total.native_position_envelope,total.native_gradient_envelope,total.max_abs_s,(unsigned long long)total.panels,(unsigned long long)total.leaves,total.max_depth,b.complete?"true":"false",(unsigned long long)b.loops,(unsigned long long)b.exterior_edges,(unsigned long long)b.interior_edges,(unsigned long long)b.nonadjacent_pairs,(unsigned long long)b.adjacent_pairs,b.min_nonadjacent_margin,b.max_adjacent_factor,b.max_sampled_trace_gap,seconds);
    need(total.complete && sample_excess==0,"local interval enclosures complete at fixed cap");need(b.complete,"boundary separation at actual finite control");return 0;
}
