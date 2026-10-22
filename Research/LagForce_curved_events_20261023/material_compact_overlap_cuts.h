#ifndef MATERIAL_COMPACT_OVERLAP_CUTS_H
#define MATERIAL_COMPACT_OVERLAP_CUTS_H
// Ordered source-edge cuts from already reconciled target patch-edge events.
// A shared target chart node contributes one cut with its complete straight
// edge star. This is a boundary arrangement input, not an overlap region.
#include "material_compact_overlap_boundary.h"
#include "material_compact_overlap_events.h"
#include <map>
#include <set>

namespace compact_overlap_cuts {
namespace cb=compact_overlap_boundary;
namespace co=compact_overlap_events;
namespace cc=compact_chart;
using Q=compact_power::Q;
inline void require(bool ok,const char*why){if(!ok)throw std::runtime_error(why);}

struct Incidence {
    uint64_t triangle=UINT64_MAX,source_hash=0;
    int rank=-1,local_edge=-1,orientation=0;
    bool operator==(const Incidence&x)const{
        return triangle==x.triangle&&source_hash==x.source_hash&&rank==x.rank&&
            local_edge==x.local_edge&&orientation==x.orientation;
    }
};
struct Cut {
    uint64_t source_edge=UINT64_MAX;
    uint64_t source_vertex=UINT64_MAX;
    int edge_owner=-1;
    Q lower=0,upper=0;
    cc::Key target_node;
    std::vector<std::pair<size_t,int>>target_segments;
    bool at_target_node=false;
    Incidence incidence[2];int incidences=0;
};
struct Budget {
    size_t maximum_events=0,maximum_cuts=0,maximum_segments_per_cut=0;
};
struct Result {
    std::vector<Cut>cut;
    size_t grouped_endpoint_events=0;
    bool complete=false;
};
inline Result order(const cb::Index&boundary,const std::vector<co::Canonical>&canonical,
        const Budget&budget){
    require(budget.maximum_events>0&&budget.maximum_cuts>0&&
        budget.maximum_segments_per_cut>0&&canonical.size()<=budget.maximum_events,
        "finite canonical source-edge cut arrangement budgets");
    struct Item {size_t source=0;cc::Key node;bool at_node=false;};
    std::vector<Item>item;item.reserve(canonical.size());
    for(size_t j=0;j<canonical.size();j++){
        const auto&e=canonical[j];
        require(e.source_registry_hash==boundary.source_registry_hash()&&
            e.source_epoch==boundary.source_epoch()&&
            e.target_chart_hash==boundary.chart_hash()&&
            e.weight_hash==boundary.weight_hash()&&
            e.target_segment<boundary.segments()&&
            e.target_endpoint>=-1&&e.target_endpoint<=1&&
            (e.source_vertex==UINT64_MAX||
                (e.target_endpoint>=0&&e.event.lower==e.event.upper&&
                 (e.event.lower==0||e.event.lower==1)))&&
            e.edge_owner>=0&&
            e.event.lower>=0&&e.event.lower<=e.event.upper&&e.event.upper<=1&&
            (e.incidences==1||e.incidences==2),
            "reconciled source-edge event belongs to one sealed target chart and accepted source");
        for(int k=0;k<e.incidences;k++)require(e.incident_triangle[k]!=UINT64_MAX&&
            e.source_hash[k]&&e.source_rank[k]>=0&&
            e.local_edge[k]>=0&&e.local_edge[k]<3&&
            (e.orientation[k]==1||e.orientation[k]==-1),
            "each source-edge cut retains its original oriented triangle incidence");
        Item v;v.source=j;
        if(e.target_endpoint>=0){
            v.at_node=true;
            v.node=boundary.trace_key(e.target_segment,e.target_endpoint,e.event.edge,
                e.source_vertex);
            require(!v.node.empty(),"canonical target node key at exact source-edge crossing");
        }
        item.push_back(std::move(v));
    }
    std::sort(item.begin(),item.end(),[&](const Item&a,const Item&b){
        const auto&A=canonical[a.source],&B=canonical[b.source];
        return std::tie(A.event.edge,A.event.lower,A.event.upper,A.target_segment)<
            std::tie(B.event.edge,B.event.lower,B.event.upper,B.target_segment);
    });
    Result out;
    for(const auto&i:item){
        const auto&e=canonical[i.source];
        bool overlapping=!out.cut.empty()&&out.cut.back().source_edge==e.event.edge&&
            e.event.lower<=out.cut.back().upper;
        if(overlapping){
            auto&cut=out.cut.back();
            require(cut.at_target_node&&i.at_node&&cut.target_node==i.node,
                "overlapping source-edge root intervals require the same canonical target chart node");
            require(cut.incidences==e.incidences,
                "one target node cut retains the same complete original source edge incidence count");
            require(cut.edge_owner==e.edge_owner,
                "one target node cut has one canonical original source edge owner");
            require(cut.source_vertex==e.source_vertex,
                "one target node cut has one exact original source vertex identity");
            for(int k=0;k<e.incidences;k++)require(cut.incidence[k]==Incidence{
                    e.incident_triangle[k],e.source_hash[k],e.source_rank[k],
                    e.local_edge[k],e.orientation[k]},
                "every target segment at a source edge cut retains identical oriented source triangle incidences");
            cut.lower=std::max(cut.lower,e.event.lower);
            cut.upper=std::min(cut.upper,e.event.upper);
            require(cut.lower<=cut.upper,"one isolated source-edge root shared by target node star");
            require(cut.target_segments.size()<budget.maximum_segments_per_cut,
                "finite target segment star at one source-edge cut");
            cut.target_segments.push_back({size_t(e.target_segment),e.target_endpoint});
            out.grouped_endpoint_events++;
        }else{
            require(out.cut.size()<budget.maximum_cuts,"finite ordered source-edge cut budget");
            Cut cut;cut.source_edge=e.event.edge;cut.source_vertex=e.source_vertex;
            cut.edge_owner=e.edge_owner;
            cut.lower=e.event.lower;cut.upper=e.event.upper;
            cut.at_target_node=i.at_node;cut.target_node=i.node;
            cut.target_segments.push_back({size_t(e.target_segment),e.target_endpoint});
            cut.incidences=e.incidences;
            for(int k=0;k<e.incidences;k++)cut.incidence[k]={
                e.incident_triangle[k],e.source_hash[k],e.source_rank[k],
                e.local_edge[k],e.orientation[k]};
            out.cut.push_back(std::move(cut));
        }
    }
    std::set<std::pair<uint64_t,cc::Key>>seen_node;
    struct Basis {Incidence incidence[2];int count=0,owner=-1;};
    std::map<uint64_t,Basis>edge_basis;
    for(auto&cut:out.cut){
        std::sort(cut.target_segments.begin(),cut.target_segments.end());
        require(std::adjacent_find(cut.target_segments.begin(),cut.target_segments.end())==
            cut.target_segments.end(),"one target segment endpoint appears once per source-edge cut");
        auto basis=edge_basis.emplace(cut.source_edge,Basis{});
        if(basis.second){
            basis.first->second.count=cut.incidences;
            basis.first->second.owner=cut.edge_owner;
            for(int k=0;k<cut.incidences;k++)basis.first->second.incidence[k]=cut.incidence[k];
        }else{
            require(basis.first->second.count==cut.incidences,
                "all cuts of one original source edge retain its complete incidence count");
            require(basis.first->second.owner==cut.edge_owner,
                "all cuts of one original source edge retain one canonical edge owner");
            for(int k=0;k<cut.incidences;k++)require(
                basis.first->second.incidence[k]==cut.incidence[k],
                "all cuts of one original source edge retain the same oriented source triangles");
        }
        if(cut.at_target_node){
            require(seen_node.insert({cut.source_edge,cut.target_node}).second,
                "one canonical target node has one ordered parameter cut on an injective source edge");
            auto star=boundary.trace_star(cut.source_edge,cut.target_node,
                cut.source_vertex);
            std::sort(star.begin(),star.end());
            require(cut.target_segments==star,
                "every registered target straight segment at a source exterior node contributes its event");
        }
    }
    out.complete=true;return out;
}
}
#endif
