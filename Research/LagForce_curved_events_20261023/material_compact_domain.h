#ifndef MATERIAL_COMPACT_DOMAIN_H
#define MATERIAL_COMPACT_DOMAIN_H
// Boundary graph of Omega* intersect a bounded plane power cell. All vertices
// keep canonical implicit definitions/enclosures. No rounded-root chart,
// conservative packet, positive projection or compact commit is made here.
#include "material_compact_power_cell.h"
namespace compact_domain {
namespace cp=compact_power;namespace sf=shared_function;
using Q=sf::Q;using S=sf::S;using D=sf::D;
inline void require(bool b,const char*s){compact_p3::require(b,s);}
struct Domain {
    const compact_target::Registry&source;std::vector<cp::Curve>curve;
    bool accepted_generator_anchors_bound=false;
    explicit Domain(const compact_target::Registry&r):source(r){for(const auto&e:r.exterior)curve.push_back(cp::curve(e));}
    Domain(const Domain&)=delete;Domain&operator=(const Domain&)=delete;
    void bind_actual_source(const compact_source::GeometryBound&g){
        require(g.canonical_global_injective&&g.exact_shared_trace&&g.source_steps==source.source_steps&&g.source_epoch==source.source_epoch&&g.source_time==source.source_time&&g.external_edges==source.exterior.size(),"same proven accepted source domain before interior-anchor shortcut");
        require(std::find(source.source_shard_hash.begin(),source.source_shard_hash.end(),g.source_coefficient_hash)!=source.source_shard_hash.end(),"same actual source coefficient shard binding");accepted_generator_anchors_bound=true;
    }
};
struct Key {
    uint64_t kind=0;std::array<uint64_t,5>id={}; // 0 powercorner,1 sourcevertex,2 line/curve root
    bool operator<(const Key&x)const{return std::tie(kind,id)<std::tie(x.kind,x.id);}
    bool operator==(const Key&x)const{return kind==x.kind&&id==x.id;}
};
struct Vertex {
    Key key;std::array<S,2>x;
    cp::LineKey line;size_t source_curve=SIZE_MAX;
    Q parameter_lower=0,parameter_upper=0;
    // For kind2: SAME source curve root of SAME line; ordinal lives in key.
    // For kind0: canonical power-corner function; kind1: exact shared DOF.
};
struct Arc {
    size_t from=0,to=0;bool curved=false;
    cp::LineKey line;size_t source_curve=SIZE_MAX;
    Q parameter_lower[2]={},parameter_upper[2]={};
};
struct Graph {
    uint64_t cell=0,source_registry_hash=0,weight_hash=0,source_epoch=0,source_steps=0;double source_time=0;
    std::vector<Vertex>vertex;std::vector<Arc>arc;
    std::vector<std::vector<size_t>>component;
    size_t curve_candidates=0,root_calls=0,root_boxes=0,winding_calls=0;
    bool generator_anchor_used=false;
    // Topological boundary graph only. Component chart/rank/positive area,
    // source-overlap integration and physical fields still require guards.
};
struct Budget {
    cp::RootBudget root;
    size_t maximum_root_calls=0,maximum_winding_calls=0,maximum_arcs=0;
};
inline Key power_key(const cp::Corner&c){return {0,{c.key.kind,c.key.id[0],c.key.id[1],c.key.id[2],0}};}
inline Key endpoint_key(const compact_source::Exterior&e,int end){return {1,{e.vertex[end],0,0,0,0}};}
inline std::array<D,2> as_point(const std::array<S,2>&x){return {D(x[0]),D(x[1])};}
inline std::array<D,2> curve_box(const compact_source::Exterior&e){return {D(sf::range(Q(e.lower[0]),Q(e.upper[0]))),D(sf::range(Q(e.lower[1]),Q(e.upper[1])))};}
inline cp::Roots roots(const cp::Plane&p,const cp::Curve&curve,const Budget&budget,Graph&out){
    require(out.root_calls++<budget.maximum_root_calls,"finite curved-domain root-call budget");auto r=cp::intersections(p,curve,budget.root);out.root_boxes+=r.boxes;
    require(r.complete,"unresolved canonical domain root: endpoint/tangency/uncertainty/budget requires rejection, not a tie");return r;
}
inline bool inside(const Domain&domain,const std::array<S,2>&point,const Budget&budget,Graph&out){
    require(out.winding_calls++<budget.maximum_winding_calls,"finite curved-domain point-classification budget");
    // Fixed non-axis ray; no retries/tolerance shift to escape a degeneracy.
    const S slope(Q(137)/256);cp::Plane ray;ray.cell[0]=UINT64_MAX-1;ray.cell[1]=UINT64_MAX;ray.normal[0]=-slope;ray.normal[1]=S(1);ray.offset=point[1]-slope*point[0];
    int winding=0;
    for(size_t j=0;j<domain.curve.size();j++){
        const auto&e=domain.source.exterior[j];auto box=curve_box(e);
        if(cp::strict_sign(cp::value(ray,box).v))continue;
        S forward=(box[0].v-point[0])+slope*(box[1].v-point[1]);if(sf::upper(forward)<0)continue;
        auto r=roots(ray,domain.curve[j],budget,out);
        for(const auto&event:r.event){
            auto p=cp::position(domain.curve[j],sf::range(event.lower,event.upper));S dot=(p[0].v-point[0])+slope*(p[1].v-point[1]);int sign=cp::strict_sign(dot);
            require(sign!=0,"point classification unresolved at actual exterior, no boundary snap");if(sign>0)winding+=event.derivative_sign*e.orientation;
        }
    }
    require(winding==0||winding==1,"oriented one-component source winding must be zero or one");return winding==1;
}
inline void close(Graph&out){
    out.component.clear();
    require(out.vertex.size()<=out.arc.size(),"closed curved-cell boundary vertex/arc count");
    std::vector<size_t>incoming(out.vertex.size(),SIZE_MAX),outgoing(out.vertex.size(),SIZE_MAX);
    for(size_t j=0;j<out.arc.size();j++){
        const auto&a=out.arc[j];require(a.from<out.vertex.size()&&a.to<out.vertex.size()&&a.from!=a.to,"nonzero canonical boundary arc");
        require(outgoing[a.from]==SIZE_MAX&&incoming[a.to]==SIZE_MAX,"unique oriented incoming/outgoing domain boundary arcs");outgoing[a.from]=j;incoming[a.to]=j;
    }
    for(size_t j=0;j<incoming.size();j++)require(incoming[j]!=SIZE_MAX&&outgoing[j]!=SIZE_MAX,"no open or unreferenced actual-domain boundary vertex");
    std::vector<bool>seen(out.arc.size());for(size_t start=0;start<out.arc.size();start++)if(!seen[start]){
        out.component.emplace_back();size_t current=start;
        do{require(!seen[current],"domain boundary cycle closes at its own first arc");seen[current]=true;out.component.back().push_back(current);current=outgoing[out.arc[current].to];}while(current!=start);
        require(out.component.back().size()>=2,"at least two distinct arcs per selected domain component");
        if(out.component.back().size()==2)require(out.arc[out.component.back()[0]].curved||out.arc[out.component.back()[1]].curved,"two-arc component needs a true curved boundary, not a degenerate straight polygon");
    }
}
inline Graph intersect(const Domain&domain,const cp::Spatial&tree,size_t index,const cp::BoundedCell&cell,const Budget&budget){
    require(&domain.source==&tree.registry&&index<tree.registry.site.size()&&cell.cell==tree.registry.site[index].cell,"same immutable source registry/weight construction for curved clipping");
    require(cell.source_registry_hash==domain.source.wire_hash&&cell.weight_hash==tree.weight_hash,"same sealed weight trial and source registry for bounded cell");
    require(budget.maximum_root_calls>0&&budget.maximum_root_calls<=1048576&&budget.maximum_winding_calls>0&&budget.maximum_winding_calls<=65536&&budget.maximum_arcs>=4&&budget.maximum_arcs<=65536,"explicit finite whole-cell curved topology budget");
    Graph out;out.cell=cell.cell;out.source_registry_hash=domain.source.wire_hash;out.weight_hash=tree.weight_hash;out.source_epoch=domain.source.source_epoch;out.source_steps=domain.source.source_steps;out.source_time=domain.source.source_time;if(cell.corner.empty())return out;
    std::map<cp::LineKey,cp::Constraint>line;for(const auto&c:cell.corner)for(const auto&l:c.line)line.emplace(l.key,l);
    Q lower[2]={HUGE_VALQ,HUGE_VALQ},upper[2]={-HUGE_VALQ,-HUGE_VALQ};
    for(const auto&v:cell.corner)for(int k=0;k<2;k++){lower[k]=std::min(lower[k],sf::lower(v.x[k]));upper[k]=std::max(upper[k],sf::upper(v.x[k]));}
    std::vector<size_t>candidate;for(size_t j=0;j<domain.curve.size();j++){
        const auto&e=domain.source.exterior[j];bool overlap=true;for(int k=0;k<2;k++)overlap&=Q(e.upper[k])>=lower[k]&&Q(e.lower[k])<=upper[k];if(overlap)candidate.push_back(j);
    }out.curve_candidates=candidate.size();
    std::map<Key,size_t>vertex;
    auto add_vertex=[&](const Vertex&v){
        auto p=vertex.emplace(v.key,out.vertex.size());if(p.second)out.vertex.push_back(v);
        else for(int k=0;k<2;k++)require(out.vertex[p.first->second].x[k].n==v.x[k].n&&out.vertex[p.first->second].x[k].e==v.x[k].e,"identical canonical domain vertex expressions within a cell");return p.first->second;
    };
    auto add_arc=[&](const Vertex&a,const Vertex&b,bool curved,const cp::LineKey&key,size_t source){
        require(out.arc.size()<budget.maximum_arcs,"finite actual-domain boundary arc budget");Arc arc;arc.from=add_vertex(a);arc.to=add_vertex(b);arc.curved=curved;arc.line=key;arc.source_curve=source;
        if(curved){arc.parameter_lower[0]=a.parameter_lower;arc.parameter_upper[0]=a.parameter_upper;arc.parameter_lower[1]=b.parameter_lower;arc.parameter_upper[1]=b.parameter_upper;}out.arc.push_back(arc);
    };
    // Safe interior shortcut ONLY for a generator transported by the same
    // actually certified square-source map, strictly inside this power cell.
    if(candidate.empty()&&domain.accepted_generator_anchors_bound){
        const auto&q=tree.registry.site[index];bool anchor=true;std::array<S,2>x;
        for(int k=0;k<2;k++){anchor&=q.reference[k]>-12&&q.reference[k]<12;x[k]=compact_source::canonical_position(q,k);}
        for(const auto&entry:line)anchor&=cp::strict_sign(cp::value(entry.second.plane,as_point(x)).v)*entry.second.side<0;
        if(anchor){
            out.generator_anchor_used=true;for(size_t j=0;j<cell.corner.size();j++){
                const auto&a=cell.corner[j];const auto&b=cell.corner[(j+1)%cell.corner.size()];Vertex A;A.key=power_key(a);A.x=a.x;Vertex B;B.key=power_key(b);B.x=b.x;
                add_arc(A,B,false,cp::common_line(a,b).key,SIZE_MAX);
            }close(out);return out;
        }
    }
    std::map<cp::LineKey,std::vector<Vertex>>edge_event;std::map<size_t,std::vector<Vertex>>curve_event;
    for(size_t j:candidate){
        const auto&e=domain.source.exterior[j];auto box=curve_box(e);bool whole_outside=false;
        for(const auto&entry:line)whole_outside|=cp::strict_sign(cp::value(entry.second.plane,box).v)*entry.second.side>0;
        if(whole_outside)continue;
        for(const auto&entry:line){
            const auto&constraint=entry.second;if(cp::strict_sign(cp::value(constraint.plane,box).v))continue;
            auto r=roots(constraint.plane,domain.curve[j],budget,out);
            for(const auto&event:r.event){
                auto position=cp::position(domain.curve[j],sf::range(event.lower,event.upper));Vertex v;v.x={position[0].v,position[1].v};bool keep=true;
                for(const auto&other:line)if(!(other.first==entry.first)){
                    int sign=cp::strict_sign(cp::value(other.second.plane,as_point(v.x)).v);
                    require(sign!=0,"unresolved curve intersection at a power corner; canonical multiple-event adapter required");if(sign*other.second.side>0){keep=false;break;}
                }
                if(!keep)continue;v.key={2,{entry.first.kind,entry.first.a,entry.first.b,e.edge,event.ordinal}};v.line=entry.first;v.source_curve=j;v.parameter_lower=event.lower;v.parameter_upper=event.upper;
                edge_event[entry.first].push_back(v);curve_event[j].push_back(v);
            }
        }
    }
    for(size_t j=0;j<cell.corner.size();j++){
        const auto&a=cell.corner[j];const auto&b=cell.corner[(j+1)%cell.corner.size()];auto support=cp::common_line(a,b);auto events=edge_event[support.key];
        int axis=sf::mag(b.x[0]-a.x[0])>=sf::mag(b.x[1]-a.x[1])?0:1;S delta=b.x[axis]-a.x[axis];require(cp::strict_sign(delta)!=0,"nonzero canonical straight boundary edge");
        struct Along {Q lo,hi;Vertex v;};std::vector<Along>point;
        Vertex A;A.key=power_key(a);A.x=a.x;Vertex B;B.key=power_key(b);B.x=b.x;point.push_back({0,0,A});
        for(const auto&v:events){S t=(v.x[axis]-a.x[axis])/delta;require(sf::lower(t)>0&&sf::upper(t)<1,"root strictly inside canonical straight edge, not snapped to corner");point.push_back({sf::lower(t),sf::upper(t),v});}point.push_back({1,1,B});
        std::sort(point.begin(),point.end(),[](const Along&a,const Along&b){return a.lo<b.lo;});
        for(size_t k=0;k+1<point.size();k++){
            require(point[k].hi<point[k+1].lo,"ordered disjoint straight-edge event enclosures");Q t=(point[k].hi+point[k+1].lo)/2;std::array<S,2>sample;
            for(int d=0;d<2;d++)sample[d]=a.x[d]+S(t)*(b.x[d]-a.x[d]);
            if(inside(domain,sample,budget,out))add_arc(point[k].v,point[k+1].v,false,support.key,SIZE_MAX);
        }
    }
    for(size_t j:candidate){
        const auto&e=domain.source.exterior[j];auto events=curve_event[j];
        auto endpoint=[&](int end){Vertex v;v.key=endpoint_key(e,end);v.source_curve=j;v.parameter_lower=v.parameter_upper=end;for(int k=0;k<2;k++)v.x[k]=sf::exact(e.frozen.coefficient[end?3:0].x[k]);return v;};
        events.push_back(endpoint(0));events.push_back(endpoint(1));std::sort(events.begin(),events.end(),[](const Vertex&a,const Vertex&b){return a.parameter_lower<b.parameter_lower;});
        for(size_t k=0;k+1<events.size();k++){
            require(events[k].parameter_upper<events[k+1].parameter_lower,"ordered disjoint canonical curve events");Q t=(events[k].parameter_upper+events[k+1].parameter_lower)/2;auto sample=cp::position(domain.curve[j],S(t));bool keep=true;
            for(const auto&entry:line){int sign=cp::strict_sign(cp::value(entry.second.plane,sample).v);require(sign!=0,"curve arc membership unresolved, no tolerance tie");if(sign*entry.second.side>0){keep=false;break;}}
            if(keep){if(e.orientation==1)add_arc(events[k],events[k+1],true,{},j);else{require(e.orientation==-1,"oriented source exterior");add_arc(events[k+1],events[k],true,{},j);}}
        }
    }
    close(out);return out;
}
} // namespace compact_domain
#endif
