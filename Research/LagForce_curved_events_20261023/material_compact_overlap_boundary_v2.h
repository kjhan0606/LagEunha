#ifndef MATERIAL_COMPACT_OVERLAP_BOUNDARY_H
#define MATERIAL_COMPACT_OVERLAP_BOUNDARY_H
// Complete broad phase for all straight patch edges in one certified global
// target chart, including internal fan/ear edges and physical power faces.
// A source edge/target patch-edge crossing must pass the
// two conservative AABBs, then a certified curve/plane root and a strict
// finite-segment membership test. Exact accepted exterior-trace nodes or
// shared original source vertex IDs may resolve certified vertex ties.
#include "material_compact_chart_degree.h"
#include "material_compact_source_spatial.h"
#include <map>
#include <numeric>
#include <sstream>
#include <iomanip>

namespace compact_overlap_boundary {
namespace cp=compact_power;
namespace cc=compact_chart;
namespace sf=shared_function;
using Q=sf::Q;using S=sf::S;using R=long double;
using Box=compact_source_spatial::Box;
inline void require(bool ok,const char*why){if(!ok)throw std::runtime_error(why);}

struct Budget {
    size_t maximum_segments=0,maximum_nodes=0,maximum_candidates_per_edge=0;
    size_t maximum_root_calls_per_edge=0,maximum_events_per_edge=0;
    size_t leaf_size=8,maximum_depth=64;
};
struct ProductionBudget {
    size_t maximum_global_edges=0,maximum_local_incidences=0,maximum_local_unique_edges=0;
    size_t maximum_total_candidate_segments=0,maximum_total_root_calls=0;
    size_t maximum_total_root_boxes=0,maximum_local_wires=0;
};
struct Segment {
    bool physical_boundary=false;
    cp::LineKey line;
    std::array<S,2> endpoint[2];
    size_t patch[2]={SIZE_MAX,SIZE_MAX},local_edge[2]={SIZE_MAX,SIZE_MAX};
    uint64_t cell[2]={UINT64_MAX,UINT64_MAX};
    size_t uses=0;
    cp::Plane plane;
    Box physical;
};
struct Crossing {
    cp::Event event;
    size_t boundary_segment=SIZE_MAX;
    int target_endpoint=-1; // 0/1 only for a proved canonical target chart node
    uint64_t source_vertex=UINT64_MAX; // exact shared original source vertex, if proved
};
struct Result {
    std::vector<Crossing>crossing;
    size_t candidate_segments=0,root_calls=0,root_boxes=0;
    bool complete=false;
};
inline Box segment_box(const std::array<S,2>&a,const std::array<S,2>&b){
    Box box;for(int k=0;k<2;k++){
        box.lo[k]=nextafterl(R(std::min(sf::lower(a[k]),sf::lower(b[k]))),-INFINITY);
        box.hi[k]=nextafterl(R(std::max(sf::upper(a[k]),sf::upper(b[k]))),INFINITY);
    }
    require(compact_source_spatial::valid(box),"finite enclosing target straight segment box");return box;
}
inline Box source_edge_box(const cp::Curve&curve){
    auto x=cp::position(curve,sf::range(Q(0),Q(1)));
    Box box;for(int k=0;k<2;k++){
        box.lo[k]=nextafterl(R(sf::lower(x[k].v)),-INFINITY);
        box.hi[k]=nextafterl(R(sf::upper(x[k].v)),INFINITY);
    }
    require(compact_source_spatial::valid(box),"finite enclosing accepted current source edge box");return box;
}
inline Box source_edge_box(const cp::Curve&curve,Q first,Q last){
    require(first>=0&&first<last&&last<=1,
        "finite dyadic source-edge broad-phase parameter interval");
    auto x=cp::position(curve,sf::range(first,last));
    Box box;for(int k=0;k<2;k++){
        box.lo[k]=nextafterl(R(sf::lower(x[k].v)),-INFINITY);
        box.hi[k]=nextafterl(R(sf::upper(x[k].v)),INFINITY);
    }
    require(compact_source_spatial::valid(box),
        "finite enclosing accepted source subedge box");return box;
}
enum class SegmentLocation {outside,interior,endpoint0,endpoint1,unresolved};
inline SegmentLocation locate_on_segment(const Segment&segment,const cp::Curve&curve,const cp::Event&event){
    auto x=cp::position(curve,sf::range(event.lower,event.upper));
    int axis=cp::strict_sign(segment.endpoint[1][0]-segment.endpoint[0][0])?0:1;
    S delta=segment.endpoint[1][axis]-segment.endpoint[0][axis];
    require(cp::strict_sign(delta)!=0,"target straight boundary has a certified nonzero projection");
    S along=(x[axis].v-segment.endpoint[0][axis])/delta;
    if(sf::upper(along)<0||sf::lower(along)>1)return SegmentLocation::outside;
    if(sf::lower(along)>0&&sf::upper(along)<1)return SegmentLocation::interior;
    bool at0=sf::lower(along)<=0&&sf::upper(along)>=0;
    bool at1=sf::lower(along)<=1&&sf::upper(along)>=1;
    if(at0&&!at1)return SegmentLocation::endpoint0;
    if(at1&&!at0)return SegmentLocation::endpoint1;
    return SegmentLocation::unresolved;
}
inline std::pair<SegmentLocation,cp::Event> refine_ambiguous_segment_location(
        const Segment&segment,const cp::Curve&curve,cp::Event root){
    auto location=locate_on_segment(segment,curve,root);
    for(unsigned iteration=0;iteration<64&&
            location!=SegmentLocation::interior&&
            location!=SegmentLocation::outside;iteration++){
        Q a=root.lower,b=root.upper,mid=(a+b)/2;
        if(!(a<mid&&mid<b))break;
        auto whole=cp::value_centered(segment.plane,curve,sf::range(a,b));
        if(cp::strict_sign(whole.x)!=root.derivative_sign)break;
        auto center=cp::value_centered(segment.plane,curve,S(mid));
        S newton=S(mid)-center.v/whole.x;
        Q na=std::max(a,sf::lower(newton));
        Q nb=std::min(b,sf::upper(newton));
        require(na<=nb,"nonempty certified local source-edge root refinement");
        if(nb-na<(b-a)*Q(.75L)){root.lower=na;root.upper=nb;}
        else{
            int sign=cp::strict_sign(center.v);
            if(!sign)break;
            if(sign==-root.derivative_sign)root.lower=mid;
            else root.upper=mid;
        }
        location=locate_on_segment(segment,curve,root);
    }
    return {location,root};
}

class Index {
    struct TraceParameter {uint64_t source_hash=0,source_epoch=0;S parameter;int source_endpoint=-1;};
    struct Node {Box box;size_t begin=0,end=0,left=SIZE_MAX,right=SIZE_MAX;};
    std::vector<Segment>segment_;
    std::vector<size_t>order_;
    std::vector<Node>node_;
    std::map<std::tuple<size_t,int,uint64_t>,TraceParameter>segment_trace_;
    std::map<std::pair<size_t,int>,cc::Key>trace_key_;
    std::map<std::pair<uint64_t,cc::Key>,std::vector<std::pair<size_t,int>>>trace_star_;
    std::map<cc::Key,std::vector<std::pair<size_t,int>>>node_star_;
    Budget budget_;
    size_t physical_segments_=0,internal_segments_=0;
    uint64_t source_registry_hash_=0,source_epoch_=0,weight_hash_=0,chart_hash_=0;
    std::array<uint64_t,4>source_shard_hash_={};
    static bool source_vertex_key(const cc::Key&key,uint64_t vertex){
        return key.size()==7&&key[0]==0&&key[1]==1&&key[2]==vertex&&
            key[3]==0&&key[4]==0&&key[5]==0&&key[6]==0;
    }
    bool accepted_shard(uint64_t hash)const{
        return hash&&std::find(source_shard_hash_.begin(),source_shard_hash_.end(),hash)!=
            source_shard_hash_.end();
    }
    int certified_source_vertex_endpoint(size_t id,int endpoint,const cp::Curve&curve)const{
        if(segment_[id].physical_boundary||curve.source_epoch!=source_epoch_||
            !accepted_shard(curve.source_hash))return -1;
        auto found=trace_key_.find({id,endpoint});
        if(found==trace_key_.end())return -1;
        for(int side=0;side<2;side++)if(curve.vertex[side]!=UINT64_MAX&&
            source_vertex_key(found->second,curve.vertex[side])){
            for(int d=0;d<2;d++)require(cc::same_ball(
                segment_[id].endpoint[endpoint][d],
                sf::exact(curve.coefficient[side?3:0].x[d])),
                "shared original source vertex has the identical accepted coefficient in both curves");
            return side;
        }
        return -1;
    }
    size_t build(size_t begin,size_t end,size_t depth){
        require(begin<end&&depth<=budget_.maximum_depth&&node_.size()<budget_.maximum_nodes,
            "finite target straight-boundary BVH node/depth budget");
        size_t at=node_.size();node_.emplace_back();
        Box box=segment_[order_[begin]].physical;
        for(size_t p=begin+1;p<end;p++)box=compact_source_spatial::hull(box,segment_[order_[p]].physical);
        node_[at].box=box;node_[at].begin=begin;node_[at].end=end;
        if(end-begin<=budget_.leaf_size)return at;
        int axis=box.hi[1]-box.lo[1]>box.hi[0]-box.lo[0]?1:0;
        size_t middle=begin+(end-begin)/2;
        std::nth_element(order_.begin()+begin,order_.begin()+middle,order_.begin()+end,
            [&](size_t a,size_t b){
                const auto&A=segment_[a].physical,&B=segment_[b].physical;
                R ca=A.lo[axis]/2+A.hi[axis]/2,cb=B.lo[axis]/2+B.hi[axis]/2;
                return ca<cb||(ca==cb&&a<b);
            });
        size_t left=build(begin,middle,depth+1),right=build(middle,end,depth+1);
        node_[at].left=left;node_[at].right=right;return at;
    }
    bool query(const Box&box,std::vector<size_t>&result)const{
        require(compact_source_spatial::valid(box),"finite current source edge box query");
        result.clear();std::vector<size_t>stack={0};
        while(!stack.empty()){
            size_t at=stack.back();stack.pop_back();const auto&node=node_[at];
            if(!compact_source_spatial::intersects(node.box,box))continue;
            if(node.left==SIZE_MAX){
                for(size_t p=node.begin;p<node.end;p++){
                    size_t id=order_[p];
                    if(compact_source_spatial::intersects(segment_[id].physical,box)){
                        if(result.size()==budget_.maximum_candidates_per_edge)
                            return false;
                        result.push_back(id);
                    }
                }
            }else{stack.push_back(node.right);stack.push_back(node.left);}
        }
        std::sort(result.begin(),result.end());return true;
    }
    bool certified_trace_endpoint(size_t id,int endpoint,const cp::Curve&curve,
            const cp::Event&event)const{
        auto found=segment_trace_.find(std::make_tuple(id,endpoint,curve.edge));
        if(found==segment_trace_.end())return false;
        const auto&trace=found->second;
        require(trace.source_hash==curve.source_hash&&trace.source_epoch==curve.source_epoch,
            "target chart exterior node and current source edge have the same accepted trace seal");
        Q lo=sf::lower(trace.parameter),hi=sf::upper(trace.parameter);
        require(std::max(lo,event.lower)<=std::min(hi,event.upper)&&
            (trace.source_endpoint<0||
                (event.lower==Q(trace.source_endpoint)&&event.upper==event.lower))&&
            (trace.source_endpoint>=0||(lo>=0&&hi<=1)),
            "same exact exterior-trace parameter lies in the isolated source-edge root");
        const auto&node=segment_[id].endpoint[endpoint];
        auto point=cp::position(curve,trace.parameter);
        for(int k=0;k<2;k++)require(sf::lower(point[k].v)<=sf::upper(node[k])&&
            sf::lower(node[k])<=sf::upper(point[k].v),
            "target exterior node and source edge enclose the same physical point");
        S union_parameter=trace.source_endpoint>=0?S(Q(trace.source_endpoint)):
            sf::range(std::min(lo,event.lower),std::max(hi,event.upper));
        auto value=cp::value(segment_[id].plane,cp::position(curve,union_parameter));
        require(cp::strict_sign(value.x)==event.derivative_sign,
            "strictly monotone curve/target line proves the shared endpoint is the unique isolated root");
        return true;
    }
public:
    Index(const cp::Spatial&sites,const std::vector<cc::Candidate>&candidate,
            const cc::Certificate&topology,const Budget&budget):budget_(budget),
            source_registry_hash_(topology.source_registry_hash),source_epoch_(topology.source_epoch),
            weight_hash_(topology.weight_hash),chart_hash_(topology.chart_hash),
            source_shard_hash_(sites.registry.source_shard_hash){
        require(budget_.maximum_segments>0&&budget_.maximum_nodes>0&&
            budget_.maximum_candidates_per_edge>0&&budget_.maximum_root_calls_per_edge>0&&
            budget_.maximum_events_per_edge>0&&budget_.leaf_size>0&&
            budget_.maximum_depth>0&&budget_.maximum_depth<=64,
            "explicit finite target boundary search and source-edge event budgets");
        cc::require_same_chart_complex(sites.registry,candidate,topology);
        require(topology.disk_topology&&topology.exact_C0_connectivity&&
            topology.oriented_source_boundary_complete&&
            topology.source_registry_hash==sites.registry.wire_hash&&
            topology.weight_hash==sites.weight_hash,
            "complete certified global target chart and one sealed trial weight vector");
        std::map<cc::Key,std::map<uint64_t,TraceParameter>>trace_node;
        for(const auto&ref:topology.patch){
            const auto&patch=candidate.at(ref.candidate).patch.at(ref.patch);
            for(int k=0;k<3;k++)if(patch.edge[k].curved){
                const auto&e=patch.edge[k];
                require(e.source<sites.registry.exterior.size(),
                    "curved chart endpoint has one accepted source exterior trace");
                const auto&exterior=sites.registry.exterior[e.source];
                for(int end=0;end<2;end++){
                    const auto&key=patch.vertex[(k+end)%3].key;
                    int source_endpoint=-1;
                    if(key.size()==7&&key[0]==0&&key[1]==1&&
                        key[3]==0&&key[4]==0&&key[5]==0&&key[6]==0){
                        if(key[2]==exterior.vertex[0])source_endpoint=0;
                        else if(key[2]==exterior.vertex[1])source_endpoint=1;
                        require(source_endpoint>=0&&
                            sf::lower(e.parameter[end])<=Q(source_endpoint)&&
                            sf::upper(e.parameter[end])>=Q(source_endpoint),
                            "source exterior vertex key and original curve endpoint parameter agree");
                        for(int d=0;d<2;d++)require(cc::same_ball(
                            patch.vertex[(k+end)%3].x[d],
                            sf::exact(exterior.frozen.coefficient[source_endpoint?3:0].x[d])),
                            "canonical chart source vertex is the same accepted exterior coefficient");
                    }
                    TraceParameter trace{exterior.frozen.source_hash,
                        exterior.frozen.source_epoch,e.parameter[end],source_endpoint};
                    auto inserted=trace_node[key].emplace(exterior.edge,trace);
                    if(!inserted.second)require(
                        inserted.first->second.source_hash==trace.source_hash&&
                        inserted.first->second.source_epoch==trace.source_epoch&&
                        cc::same_ball(inserted.first->second.parameter,trace.parameter)&&
                        inserted.first->second.source_endpoint==trace.source_endpoint,
                        "one canonical exterior parameter at a shared target chart node");
                }
            }
        }
        std::map<std::pair<cc::Key,cc::Key>,size_t>unique;
        for(size_t flat=0;flat<topology.patch.size();flat++){
            const auto&ref=topology.patch[flat];
            const auto&patch=candidate.at(ref.candidate).patch.at(ref.patch);
            for(int k=0;k<3;k++){
                const auto&e=patch.edge[k];if(e.curved)continue;
                if(e.boundary)require(e.line.kind==1&&e.line.a<e.line.b&&
                    (ref.cell==e.line.a||ref.cell==e.line.b),
                    "actual target straight power boundary has one sorted pair of cells");
                else require(e.source==SIZE_MAX,
                    "internal straight patch edge has no exterior curve source");
                cc::Key a=patch.vertex[k].key,b=patch.vertex[(k+1)%3].key;
                std::array<S,2>xa=patch.vertex[k].x,xb=patch.vertex[(k+1)%3].x;
                if(b<a){std::swap(a,b);std::swap(xa,xb);}
                auto inserted=unique.emplace(std::make_pair(a,b),segment_.size());
                if(!inserted.second){
                    auto&prior=segment_[inserted.first->second];
                    require(prior.uses==1&&prior.physical_boundary==e.boundary&&
                        (!e.boundary||prior.line==e.line),
                        "paired target patch edge has the same canonical physical role and line");
                    for(int d=0;d<2;d++)require(cc::same_ball(prior.endpoint[0][d],xa[d])&&
                        cc::same_ball(prior.endpoint[1][d],xb[d]),
                        "paired target patch edge has identical canonical endpoint enclosures");
                    prior.patch[1]=flat;prior.local_edge[1]=k;prior.cell[1]=ref.cell;prior.uses=2;
                    continue;
                }
                require(segment_.size()<budget_.maximum_segments,
                    "finite complete straight target patch-edge registry");
                Segment out;out.physical_boundary=e.boundary;out.line=e.line;
                out.endpoint[0]=xa;out.endpoint[1]=xb;
                out.patch[0]=flat;out.local_edge[0]=k;out.cell[0]=ref.cell;out.uses=1;
                if(e.boundary){
                    auto lo=std::lower_bound(sites.registry.site.begin(),sites.registry.site.end(),e.line.a,
                        [](const compact_source::Site&s,uint64_t id){return s.cell<id;});
                    auto hi=std::lower_bound(sites.registry.site.begin(),sites.registry.site.end(),e.line.b,
                        [](const compact_source::Site&s,uint64_t id){return s.cell<id;});
                    require(lo!=sites.registry.site.end()&&lo->cell==e.line.a&&
                        hi!=sites.registry.site.end()&&hi->cell==e.line.b,
                        "both canonical power plane generators exist in the sealed registry");
                    size_t i=lo-sites.registry.site.begin(),j=hi-sites.registry.site.begin();
                    out.plane=cp::plane(*lo,sites.weight[i],*hi,sites.weight[j]);
                    physical_segments_++;
                }else{
                    S dx=xb[0]-xa[0],dy=xb[1]-xa[1];
                    require(cp::strict_sign(dx*dx+dy*dy)==1,
                        "positive canonical internal target patch edge length");
                    out.plane.cell[0]=ref.cell;out.plane.cell[1]=UINT64_MAX;
                    out.plane.origin[0]=xa[0];out.plane.origin[1]=xa[1];
                    out.plane.normal[0]=-dy;out.plane.normal[1]=dx;
                    out.plane.offset=S(0);internal_segments_++;
                }
                out.physical=segment_box(xa,xb);segment_.push_back(std::move(out));
                for(int end=0;end<2;end++){
                    const auto&key=end?b:a;auto trace=trace_node.find(key);
                    if(trace!=trace_node.end())for(const auto&entry:trace->second){
                        require(segment_trace_.emplace(
                            std::make_tuple(segment_.size()-1,end,entry.first),entry.second).second,
                            "unique original exterior trace parameter at target segment endpoint");
                        auto inserted_key=trace_key_.emplace(
                            std::make_pair(segment_.size()-1,end),key);
                        require(inserted_key.second||inserted_key.first->second==key,
                            "one canonical target chart key per straight segment endpoint");
                        if(inserted_key.second)node_star_[key].push_back({segment_.size()-1,end});
                        trace_star_[{entry.first,key}].push_back({segment_.size()-1,end});
                    }
                }
            }
        }
        require(!segment_.empty()&&physical_segments_==topology.straight_shared_edges&&
            internal_segments_==topology.internal_edges,
            "all certified physical and internal straight target patch edges indexed once");
        for(const auto&s:segment_)require(s.uses==2&&
            (s.physical_boundary?s.cell[0]!=s.cell[1]:s.cell[0]==s.cell[1]),
            "each straight target patch edge has exactly its certified two patch incidences");
        order_.resize(segment_.size());std::iota(order_.begin(),order_.end(),size_t(0));
        node_.reserve(std::min(budget_.maximum_nodes,2*segment_.size()));
        (void)build(0,order_.size(),0);
    }
    size_t segments()const{return segment_.size();}
    size_t physical_segments()const{return physical_segments_;}
    size_t internal_segments()const{return internal_segments_;}
    size_t trace_ties()const{return segment_trace_.size();}
    const cc::Key&trace_key(size_t segment,int endpoint,uint64_t source_edge,
            uint64_t source_vertex=UINT64_MAX)const{
        if(source_vertex==UINT64_MAX)require(
            segment_trace_.count(std::make_tuple(segment,endpoint,source_edge))==1,
            "source exterior trace exists at this canonical target segment endpoint");
        auto found=trace_key_.find({segment,endpoint});
        require(found!=trace_key_.end(),"sparse canonical target chart node identity");
        if(source_vertex!=UINT64_MAX)require(
            source_vertex_key(found->second,source_vertex),
            "target chart endpoint is the same original source vertex");
        return found->second;
    }
    const std::vector<std::pair<size_t,int>>&trace_star(uint64_t source_edge,
            const cc::Key&node,uint64_t source_vertex=UINT64_MAX)const{
        if(source_vertex!=UINT64_MAX){
            require(source_vertex_key(node,source_vertex),
                "same original source vertex defines the target chart node star");
            auto found=node_star_.find(node);
            require(found!=node_star_.end(),"complete target straight-segment star at source vertex");
            return found->second;
        }
        auto found=trace_star_.find({source_edge,node});
        require(found!=trace_star_.end(),"known target chart node on the original exterior source edge");
        return found->second;
    }
    const Segment&segment(size_t id)const{require(id<segment_.size(),"known target patch edge identity");return segment_[id];}
    size_t nodes()const{return node_.size();}
    uint64_t source_registry_hash()const{return source_registry_hash_;}
    uint64_t source_epoch()const{return source_epoch_;}
    uint64_t weight_hash()const{return weight_hash_;}
    uint64_t chart_hash()const{return chart_hash_;}
    Result intersections(const cp::Curve&curve,const cp::RootBudget&root_budget)const{
        require(curve.source_epoch==source_epoch_&&curve.source_hash,
            "source edge and complete target boundary use the same accepted epoch");
        Result out;std::vector<size_t>possible;
        bool complete_candidates=query(source_edge_box(curve),possible);
        if(!complete_candidates){
            for(int parts=2;parts<=128&&!complete_candidates;parts*=2){
                std::vector<size_t>all,piece;bool every_subedge_complete=true;
                for(int j=0;j<parts;j++){
                    Q first=Q(j)/parts,last=Q(j+1)/parts;
                    if(!query(source_edge_box(curve,first,last),piece)){
                        every_subedge_complete=false;break;
                    }
                    all.insert(all.end(),piece.begin(),piece.end());
                }
                if(!every_subedge_complete)continue;
                std::sort(all.begin(),all.end());
                all.erase(std::unique(all.begin(),all.end()),all.end());
                if(all.size()>budget_.maximum_candidates_per_edge)continue;
                possible=std::move(all);complete_candidates=true;
                fprintf(stderr,
                    "CERTIFIED_SUBEDGE_BROADPHASE edge=%llu parts=%d candidates=%zu\n",
                    (unsigned long long)curve.edge,parts,possible.size());
            }
        }
        require(complete_candidates,
            "finite certified source subedge target candidates after 128 dyadic boxes");
        out.candidate_segments=possible.size();
        std::map<std::tuple<uint64_t,uint64_t,uint64_t>,std::vector<size_t>>by_line;
        for(size_t id:possible){
            const auto&s=segment_[id];
            auto key=s.physical_boundary?
                std::make_tuple(uint64_t(0),s.line.a,s.line.b):
                std::make_tuple(uint64_t(1),uint64_t(id),uint64_t(0));
            by_line[key].push_back(id);
        }
        for(const auto&item:by_line){
            require(out.root_calls<budget_.maximum_root_calls_per_edge,
                "finite source-edge/target-plane root-call budget");
            const auto&plane=segment_[item.second.front()].plane;
            bool exact_endpoint[2]={false,false};
            for(size_t id:item.second)for(int end=0;end<2;end++){
                // An internal line is constructed through its chart node.
                // A physical power plane needs a separate algebraic endpoint
                // witness; sharing the vertex key alone does not prove it.
                if(segment_[id].physical_boundary)continue;
                int shared=certified_source_vertex_endpoint(id,end,curve);
                if(shared>=0)exact_endpoint[shared]=true;
                auto found=segment_trace_.find(std::make_tuple(id,end,curve.edge));
                if(found==segment_trace_.end())continue;
                const auto&trace=found->second;
                if(trace.source_endpoint<0||trace.source_hash!=curve.source_hash||
                    trace.source_epoch!=curve.source_epoch)continue;
                exact_endpoint[trace.source_endpoint]=true;
            }
            auto roots=cp::intersections(plane,curve,root_budget,
                exact_endpoint[0],exact_endpoint[1]);out.root_calls++;out.root_boxes+=roots.boxes;
            if(!roots.complete){
                std::string trace="";
                for(size_t id:item.second)for(int end=0;end<2;end++){
                    auto found=segment_trace_.find(std::make_tuple(id,end,curve.edge));
                    if(found!=segment_trace_.end())trace+=" seg"+std::to_string(id)+
                        ":end"+std::to_string(end)+
                        ":u"+std::to_string(double(found->second.parameter.n))+
                        ":source-end"+std::to_string(found->second.source_endpoint);
                }
                throw std::runtime_error(
                "unresolved source-edge/target-boundary root edge="+std::to_string(curve.edge)+
                " segment="+std::to_string(item.second.front())+
                " physical="+std::to_string(int(segment_[item.second.front()].physical_boundary))+
                " reason="+std::to_string(int(roots.unresolved))+
                " witness=["+std::to_string(double(roots.witness_lower))+","+
                    std::to_string(double(roots.witness_upper))+"]"+trace);
            }
            for(const auto&event:roots.event){
                size_t matched=SIZE_MAX;int tie=-1;uint64_t source_vertex=UINT64_MAX;
                cp::Event matched_event=event;
                for(size_t id:item.second){
                    cp::Event local_event=event;
                    auto where=locate_on_segment(segment_[id],curve,local_event);
                    if(where==SegmentLocation::outside)continue;
                    int endpoint=-1;bool endpoint_proved=false;
                    if(where==SegmentLocation::endpoint0||where==SegmentLocation::endpoint1){
                        endpoint=where==SegmentLocation::endpoint0?0:1;
                        int shared=certified_source_vertex_endpoint(id,endpoint,curve);
                        if(shared>=0&&local_event.lower==Q(shared)&&
                                local_event.upper==local_event.lower){
                            auto derivative=cp::value(segment_[id].plane,
                                cp::position(curve,S(Q(shared))));
                            require(cp::strict_sign(derivative.x)==local_event.derivative_sign,
                                "same original source vertex is a unique transverse target-edge root");
                            source_vertex=curve.vertex[shared];endpoint_proved=true;
                        }else if(certified_trace_endpoint(id,endpoint,curve,local_event))
                            endpoint_proved=true;
                        else{
                            auto refined=refine_ambiguous_segment_location(
                                segment_[id],curve,local_event);
                            where=refined.first;local_event=refined.second;
                            if(where==SegmentLocation::outside)continue;
                            if(where==SegmentLocation::interior)endpoint=-1;
                        }
                    }else if(where==SegmentLocation::unresolved){
                        auto refined=refine_ambiguous_segment_location(
                            segment_[id],curve,local_event);
                        where=refined.first;local_event=refined.second;
                        if(where==SegmentLocation::outside)continue;
                    }
                    if(where==SegmentLocation::endpoint0||
                            where==SegmentLocation::endpoint1||
                            where==SegmentLocation::unresolved){
                        if(!endpoint_proved){
                            auto x=cp::position(curve,sf::range(local_event.lower,local_event.upper));
                            int axis=cp::strict_sign(segment_[id].endpoint[1][0]-
                                segment_[id].endpoint[0][0])?0:1;
                            S along=(x[axis].v-segment_[id].endpoint[0][axis])/
                                (segment_[id].endpoint[1][axis]-segment_[id].endpoint[0][axis]);
                            std::ostringstream message;message<<std::setprecision(21)
                                <<"source edge root near an unproved target patch vertex"
                                <<" edge="<<curve.edge<<" segment="<<id
                                <<" endpoint="<<endpoint
                                <<" root=["<<R(local_event.lower)<<","<<R(local_event.upper)<<"]"
                                <<" along=["<<R(sf::lower(along))<<","<<R(sf::upper(along))<<"]";
                            throw std::runtime_error(message.str());
                        }
                    }else require(where==SegmentLocation::interior,
                        "source edge root has strict interior target segment membership");
                    require(matched==SIZE_MAX,"one source-edge root lies on two target boundary segments");
                    matched=id;tie=endpoint;matched_event=local_event;
                }
                if(matched!=SIZE_MAX){
                    require(out.crossing.size()<budget_.maximum_events_per_edge,
                        "finite complete source-edge/target-boundary event budget");
                    out.crossing.push_back({matched_event,matched,tie,source_vertex});
                }
            }
        }
        std::sort(out.crossing.begin(),out.crossing.end(),[](const Crossing&a,const Crossing&b){
            return std::tie(a.event.edge,a.boundary_segment,a.event.ordinal)<
                std::tie(b.event.edge,b.boundary_segment,b.event.ordinal);
        });
        out.complete=true;return out;
    }
};
}
#endif
