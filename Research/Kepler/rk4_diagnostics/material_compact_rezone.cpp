// Components of the one full accepted-state compact-rezoning transaction.
// No standalone partial-run main: the complete transfer/EOS/step controller
// must be linked before any physical transaction is launched.
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <climits>
#include <type_traits>
#include "material_compact_checkpoint.h"
#include "material_compact_geometry.h"
#include "material_compact_volume.h"
#include "material_compact_thermo.h"
using R=long double;
static void need(bool b,const char*s){compact_p3::require(b,s);}
#include "material_full_owned.h"
#include "material_distributed_mass.h"
#include "material_compact_mass.h"
#include "material_compact_operator.h"
namespace st=shared_trace;
#include "material_compact_source.h"
#include "material_compact_source_geometry.h"
#include "material_compact_source_context.h"
#include "material_compact_source_target.h"
#include "material_compact_target_registry.h"
#include "material_compact_mass_law.h"
#include "material_compact_power_curve.h"
#include "material_compact_overlap_events.h"
#include "material_compact_power_cell.h"
#include "material_compact_domain.h"
#include "material_compact_chart_topology.h"
#include "material_compact_chart_certificate.h"
#include "material_compact_chart_jacobian.h"
#include "material_compact_chart_area.h"
#include "material_compact_chart_degree.h"
#include "material_compact_moments.h"
#include "material_compact_source_fields.h"
#include "material_compact_source_spatial.h"
#include "material_compact_overlap_boundary.h"
#include "material_compact_overlap_cuts.h"
#include "material_compact_overlap_triangle_route.h"
namespace compact_source {
Model import_accepted_state(const Inputs&input,MPI_Comm comm){return import_source(input,comm);}
GeometryBound certify_accepted_geometry(const Model&source){return bound_source_geometry(source);}
CurrentContext derive_accepted_context(Model&source,const Inputs&input,const GeometryBound&geometry){return reconstruct_endpoint_context(source,input,geometry);}
TargetSource extract_target_inputs(const Model&source,const Inputs&input,const GeometryBound&geometry,const CurrentContext&context){return target_source_inputs(source,input,geometry,context);}
}
namespace compact_source_fields {
struct PreparedActualTriangle {
    sf::Map current,generating;
    sf::RealBounds geometry_bound;
    Law law;
    uint64_t original_triangle=UINT64_MAX;
    BoundedField at(S X,S Y)const{return evaluate(current,generating,X,Y,law);}
    BoundedField at_physical(S x,S y)const{
        return evaluate_at_physical(current,generating,geometry_bound,x,y,law);
    }
    Containment contain_physical_box(S x,S y)const{
        return enclose_physical_box(current,geometry_bound,x,y);
    }
};
PreparedActualTriangle prepare_actual_source_triangle(const compact_source::Model&source,
        const compact_source::GeometryBound&geometry,size_t local_triangle){
    const auto&o=source.owned;
    require(geometry.canonical_global_injective&&geometry.exact_shared_trace&&
        local_triangle<o.f.triangle.size()&&local_triangle<geometry.complete.size()&&
        geometry.complete[local_triangle]&&source.current_context_bound&&
        source.generating.size()==o.f.coefficient.size()&&
        geometry.source_coefficient_hash==source.accepted.coefficient_hash&&
        geometry.source_epoch==source.accepted.physical.epoch&&
        geometry.source_steps==source.accepted.physical.steps&&
        geometry.source_time==source.accepted.physical.time,
        "known source triangle uses the complete accepted map, generating map and certified source geometry");
    const auto&t=o.f.triangle[local_triangle];
    require(t.cell>=0&&size_t(t.cell)<source.current_context.size()&&
        local_triangle<o.source_triangle.size(),"known original triangle and accepted cell context");
    const auto&cell=source.current_context[t.cell];
    PreparedActualTriangle out;
    out.current=sf::compile(o.f,local_triangle,o.f.coefficient);
    out.generating=sf::compile(o.f,local_triangle,source.generating);
    out.geometry_bound=geometry.triangle[local_triangle];
    out.law.gamma=source.gamma;out.law.P0=source.P0;out.law.eps=source.eps;
    out.law.current_delta=cell.delta;out.law.generating_delta=cell.generating_delta;
    out.original_triangle=o.source_triangle[local_triangle];
    return out;
}
PreparedActualTriangle prepare_owned_spatial_candidate(const compact_source::Model&source,
        const compact_source::GeometryBound&geometry,
        const compact_source_spatial::Entry&candidate){
    require(source.layout&&candidate.owner==source.layout->rank&&
        candidate.local_triangle<source.owned.source_triangle.size()&&
        source.owned.source_triangle[candidate.local_triangle]==candidate.original_triangle,
        "broad-phase candidate belongs to this actual source MPI owner and original triangle");
    return prepare_actual_source_triangle(source,geometry,candidate.local_triangle);
}
}
namespace compact_source_spatial {
namespace sf=shared_function;
static_assert(std::is_trivially_copyable<Entry>::value,"source AABB entries have a fixed MPI byte layout");
Index gather_actual_source_index(const compact_source::Model&source,
        const compact_source::GeometryBound&geometry,const Budget&budget){
    require(source.layout&&source.layout->ranks==budget.ranks&&
        geometry.canonical_global_injective&&geometry.exact_shared_trace&&
        geometry.source_coefficient_hash==source.accepted.coefficient_hash&&
        geometry.source_epoch==source.accepted.physical.epoch&&
        geometry.source_steps==source.accepted.physical.steps&&
        geometry.source_time==source.accepted.physical.time,
        "actual source broad phase uses one complete accepted geometry and MPI ownership");
    const auto&o=source.owned;require(o.f.triangle.size()==geometry.triangle.size()&&
        o.f.triangle.size()==geometry.complete.size()&&o.f.triangle.size()==o.source_triangle.size(),
        "one physical enclosure per accepted source triangle");
    std::vector<Entry>local;local.reserve(o.f.triangle.size());
    for(size_t s=0;s<o.f.triangle.size();s++){
        require(geometry.complete[s],"all source triangle whole-map enclosures complete before broad phase");
        Entry e{};e.original_triangle=o.source_triangle[s];e.local_triangle=s;e.owner=source.layout->rank;
        for(int d=0;d<2;d++){
            sf::Q q=geometry.triangle[s].displacement[d];
            require(q>=0&&::finiteq(q),"finite nonnegative whole-triangle physical displacement bound");
            R displacement=nextafterl(R(q),INFINITY),lo=o.f.triangle[s].vertex[0][d],hi=lo;
            for(int a=1;a<3;a++){
                lo=std::min(lo,R(o.f.triangle[s].vertex[a][d]));
                hi=std::max(hi,R(o.f.triangle[s].vertex[a][d]));
            }
            e.physical.lo[d]=nextafterl(lo-displacement,-INFINITY);
            e.physical.hi[d]=nextafterl(hi+displacement,INFINITY);
        }
        require(valid(e.physical),"finite conservative deformed-source physical AABB");
        local.push_back(e);
    }
    uint64_t local_count=local.size();std::vector<uint64_t>counts(source.layout->ranks);
    MPI_Allgather(&local_count,1,MPI_UINT64_T,counts.data(),1,MPI_UINT64_T,source.layout->comm);
    uint64_t total=0,bytes=0;std::vector<int>byte_counts(counts.size()),offsets(counts.size());
    for(size_t rank=0;rank<counts.size();rank++){
        require(counts[rank]<=budget.maximum_entries&&counts[rank]<=uint64_t(INT_MAX/sizeof(Entry)),
            "finite per-rank source AABB count and MPI byte count");
        uint64_t count_bytes=counts[rank]*sizeof(Entry);
        require(bytes+count_bytes<=INT_MAX,"finite global source AABB MPI byte displacement");
        offsets[rank]=int(bytes);byte_counts[rank]=int(count_bytes);bytes+=count_bytes;total+=counts[rank];
    }
    require(total==o.global_triangles&&total<=budget.maximum_entries,
        "complete bounded original source triangle registry across four owners");
    std::vector<Entry>global(total);
    MPI_Allgatherv(local.data(),int(local_count*sizeof(Entry)),MPI_BYTE,
        global.data(),byte_counts.data(),offsets.data(),MPI_BYTE,source.layout->comm);
    return Index(std::move(global),budget);
}
class TargetChartBoxes {
    const compact_target::Registry&registry_;
    const std::vector<compact_chart::Candidate>&candidate_;
    const compact_chart::Certificate&topology_;
public:
    TargetChartBoxes(const compact_target::Registry&registry,
            const std::vector<compact_chart::Candidate>&candidate,
            const compact_chart::Certificate&topology):registry_(registry),candidate_(candidate),topology_(topology){
        compact_chart::require_same_chart_complex(registry_,candidate_,topology_);
    }
    std::array<sf::S,2> interval(size_t flat_patch,sf::Q s0,sf::Q s1,sf::Q t0,sf::Q t1)const{
        require(flat_patch<topology_.patch.size()&&s0>=0&&s0<s1&&s1<=1&&
            t0>=0&&t0<t1&&t1<=1,"bounded target chart reference box");
        const auto&ref=topology_.patch[flat_patch];
        const auto&patch=candidate_[ref.candidate].patch[ref.patch];
        compact_power::Curve curve;const compact_power::Curve*ptr=nullptr;
        if(patch.edge[1].curved){
            require(patch.edge[1].source<registry_.exterior.size(),"sealed exterior curve in target chart box");
            curve=compact_power::curve(registry_.exterior[patch.edge[1].source]);ptr=&curve;
        }
        return compact_chart::chart_position_box(patch,ptr,sf::range(s0,s1),sf::range(t0,t1));
    }
    static Box enclosing_box(const std::array<sf::S,2>&position){
        Box out;for(int d=0;d<2;d++){
            out.lo[d]=nextafterl(R(sf::lower(position[d])),-INFINITY);
            out.hi[d]=nextafterl(R(sf::upper(position[d])),INFINITY);
        }
        require(valid(out),"finite conservative target physical chart box");return out;
    }
    Box box(size_t flat_patch,sf::Q s0,sf::Q s1,sf::Q t0,sf::Q t1)const{
        return enclosing_box(interval(flat_patch,s0,s1,t0,t1));
    }
};
struct WholePatchSource {
    uint64_t triangle=UINT64_MAX;
    int owner=-1;
    bool certified=false;
};
struct WholePatchResult {
    uint64_t source_registry_hash=0,source_epoch=0,weight_hash=0,chart_hash=0;
    std::vector<WholePatchSource>patch;
    uint64_t certified=0,unresolved=0,candidates=0;
    // A sufficient whole-patch containment result, not an overlap rule or
    // a proof that unresolved patches are outside any source triangle.
};
struct WholePatchBudget {
    uint64_t maximum_patches=0,maximum_total_candidates_per_rank=0;
};
WholePatchResult certify_whole_patch_source_containment(
        const compact_source::Model&source,const compact_source::GeometryBound&geometry,
        const compact_target::Registry&registry,const Budget&source_index_budget,
        const std::vector<compact_chart::Candidate>&candidate,
        const compact_chart::Certificate&topology,
        const compact_chart::DegreeOneCertificate&degree,
        const WholePatchBudget&budget){
    require(source.layout&&source.layout->ranks==4&&
        geometry.canonical_global_injective&&geometry.exact_shared_trace&&
        source.current_context_bound&&source.original_reference_receipt_bound&&
        geometry.source_coefficient_hash==source.accepted.coefficient_hash&&
        geometry.source_epoch==source.accepted.physical.epoch&&
        geometry.source_steps==source.accepted.physical.steps&&
        geometry.source_time==source.accepted.physical.time&&
        std::find(registry.source_shard_hash.begin(),registry.source_shard_hash.end(),
            geometry.source_coefficient_hash)!=registry.source_shard_hash.end()&&
        geometry.source_epoch==registry.source_epoch&&
        geometry.source_steps==registry.source_steps&&
        geometry.source_time==registry.source_time&&
        degree.degree_one_physical_coverage&&
        degree.source_registry_hash==registry.wire_hash&&
        degree.source_epoch==registry.source_epoch&&
        degree.weight_hash==topology.weight_hash&&
        degree.chart_hash==topology.chart_hash&&
        budget.maximum_patches>0&&topology.patch.size()<=budget.maximum_patches&&
        budget.maximum_total_candidates_per_rank>0,
        "same certified actual source and complete target chart for whole-patch containment");
    uint64_t seal[5]={registry.wire_hash,registry.source_epoch,
        topology.weight_hash,topology.chart_hash,uint64_t(topology.patch.size())};
    uint64_t low[5],high[5];
    require(MPI_Allreduce(seal,low,5,MPI_UINT64_T,MPI_MIN,source.layout->comm)==MPI_SUCCESS&&
        MPI_Allreduce(seal,high,5,MPI_UINT64_T,MPI_MAX,source.layout->comm)==MPI_SUCCESS,
        "collective actual whole-patch source/target seal preflight");
    for(int k=0;k<5;k++)require(low[k]==high[k],
        "all source owners classify the same sealed target chart");
    // Construct the BVH from every accepted source owner here. A caller-owned
    // index with the right entry count is not a completeness certificate.
    Index source_index=gather_actual_source_index(source,geometry,source_index_budget);
    require(source_index.size()==source.owned.global_triangles,
        "complete original source triangle BVH for whole-patch containment");
    TargetChartBoxes boxes(registry,candidate,topology);
    const size_t n=topology.patch.size();
    require(n>0&&n<=size_t(INT_MAX),"finite collective whole-patch classification size");
    std::vector<uint64_t>local_id(n,UINT64_MAX),local_owner(n,UINT64_MAX),local_count(n,0);
    WholePatchResult out;out.source_registry_hash=registry.wire_hash;
    out.source_epoch=registry.source_epoch;out.weight_hash=topology.weight_hash;
    out.chart_hash=topology.chart_hash;
    bool valid=true;
    try{
        for(size_t flat=0;flat<n;flat++){
            Box physical=boxes.box(flat,0,1,0,1);
            auto possible=source_index.query(physical);
            require(possible.size()<=budget.maximum_total_candidates_per_rank-out.candidates,
                "finite whole-patch source candidate budget without truncation");
            out.candidates+=possible.size();
            for(const auto&entry:possible){
                if(entry.owner!=source.layout->rank)continue;
                auto triangle=compact_source_fields::prepare_owned_spatial_candidate(
                    source,geometry,entry);
                auto contained=triangle.contain_physical_box(
                    sf::range(sf::Q(physical.lo[0]),sf::Q(physical.hi[0])),
                    sf::range(sf::Q(physical.lo[1]),sf::Q(physical.hi[1])));
                if(!contained.entire_physical_box_inside)continue;
                require(local_count[flat]==0,
                    "two distinct original source triangles claim the same whole target patch");
                local_id[flat]=entry.original_triangle;
                local_owner[flat]=uint64_t(entry.owner);local_count[flat]=1;
            }
        }
    }catch(const std::exception&e){
        fprintf(stderr,"rank%d WHOLE_PATCH_SOURCE_CONTAINMENT_FAIL %s\n",
            source.layout->rank,e.what());valid=false;
    }
    compact_mass::collective_require(valid,source.layout->comm,
        "bounded whole-patch source containment evaluation on every owner");
    std::vector<uint64_t>global_id(n),global_owner(n),global_count(n);
    require(MPI_Allreduce(local_id.data(),global_id.data(),int(n),MPI_UINT64_T,
            MPI_MIN,source.layout->comm)==MPI_SUCCESS&&
        MPI_Allreduce(local_owner.data(),global_owner.data(),int(n),MPI_UINT64_T,
            MPI_MIN,source.layout->comm)==MPI_SUCCESS&&
        MPI_Allreduce(local_count.data(),global_count.data(),int(n),MPI_UINT64_T,
            MPI_SUM,source.layout->comm)==MPI_SUCCESS,
        "complete MPI whole-patch source containment claims");
    out.patch.resize(n);
    for(size_t flat=0;flat<n;flat++){
        require(global_count[flat]<=1,
            "globally injective accepted source gives at most one whole-patch owner");
        if(global_count[flat]){
            require(global_id[flat]!=UINT64_MAX&&global_owner[flat]<4,
                "complete original triangle identity for a certified whole patch");
            out.patch[flat]={global_id[flat],int(global_owner[flat]),true};out.certified++;
        }else out.unresolved++;
    }
    return out;
}
}
namespace compact_domain {
Graph intersect_actual_domain(const Domain&domain,const compact_power::Spatial&sites,size_t index,const compact_power::BoundedCell&cell,const Budget&budget){require(domain.accepted_generator_anchors_bound,"actual source geometry binding required even when no anchor shortcut is used");return intersect(domain,sites,index,cell,budget);}
}
namespace compact_chart {
Candidate connect_actual_domain(const compact_domain::Domain&domain,const compact_power::Spatial&sites,const compact_domain::Graph&graph,const Budget&budget){require(domain.accepted_generator_anchors_bound&& &sites.registry== &domain.source&&graph.weight_hash==sites.weight_hash,"actual source and sealed weight binding before chart connectivity");return connect(domain,graph,budget);}
Certificate certify_actual_chart_connectivity(const compact_target::Registry&registry,const std::vector<Candidate>&candidate,const compact_mass_law::BoundLaw&law,size_t maximum_patches){
    require(law.source_registry_hash()==registry.wire_hash&&law.source_epoch()==registry.source_epoch&&law.source_steps()==registry.source_steps&&law.source_time()==registry.source_time,"same immutable fixed-M law and chart source epoch");
    return certify_global_connectivity(registry,candidate,law.law_hash(),maximum_patches);
}
JacobianCertificate certify_actual_chart_jacobians(const compact_target::Registry&registry,const std::vector<Candidate>&candidate,const Certificate&connectivity,const compact_mass_law::BoundLaw&law,const JacobianBudget&budget){
    require_same_chart_complex(registry,candidate,connectivity);
    require(connectivity.law_hash==law.law_hash()&&connectivity.exact_C0_connectivity&&connectivity.oriented_source_boundary_complete&&
        connectivity.disk_topology&&connectivity.global_components==1,"positive chart-J proof follows the same complete connected disk/source/W/mass-law chart complex");
    return certify_all_chart_jacobians(registry,candidate,budget);
}
DegreeOneCertificate certify_actual_chart_degree_one(const compact_target::Registry&registry,const std::vector<Candidate>&candidate,
        const Certificate&connectivity,const JacobianCertificate&jacobian,const compact_mass_law::BoundLaw&law,
        const compact_source::GeometryBound&source_geometry){
    require(connectivity.law_hash==law.law_hash()&&law.source_registry_hash()==registry.wire_hash&&
        law.source_epoch()==registry.source_epoch&&law.source_steps()==registry.source_steps&&law.source_time()==registry.source_time,
        "degree-one coverage uses the same sealed fixed-M law and accepted source epoch");
    return certify_degree_one(registry,candidate,connectivity,jacobian,source_geometry);
}
AreaBuilder prepare_actual_chart_area(const compact_target::Registry&registry,const std::vector<Candidate>&candidate,
        const Certificate&connectivity,const JacobianCertificate&jacobian,const DegreeOneCertificate&degree,
        const compact_mass_law::BoundLaw&law,const AreaBudget&budget){
    require(degree.degree_one_physical_coverage&&degree.source_boundary_simple&&degree.chart_is_disk&&
        degree.chart_J_positive&&degree.continuous_shared_edges&&degree.boundary_degree==1&&
        degree.chart_hash==connectivity.chart_hash&&degree.chart_hash==jacobian.chart_hash&&
        degree.law_hash==law.law_hash()&&degree.source_registry_hash==registry.wire_hash&&
        degree.weight_hash==connectivity.weight_hash&&degree.source_epoch==registry.source_epoch&&
        degree.source_steps==registry.source_steps&&degree.source_time==registry.source_time,
        "physical-area quadrature requires actual degree-one coverage for the same accepted source, chart and immutable mass law");
    return AreaBuilder(registry,candidate,connectivity,jacobian,budget);
}
}
namespace compact_overlap {
StreamAccumulator prepare_actual_overlap_stream(const compact_target::Registry&registry,
        const compact_chart::Certificate&topology,const compact_chart::DegreeOneCertificate&degree,
        const compact_mass_law::BoundLaw&law,const Seal&expected_seal,const Seal&packet_seal,
        const Seal&field_seal,
        const std::vector<uint64_t>&source_triangle_ids,
        const std::vector<R>&source_area,const std::vector<R>&source_area_error,
        const std::vector<R>&target_area,const std::vector<R>&target_area_error,
        const Budget&budget){
    require(degree.degree_one_physical_coverage&&degree.chart_hash==topology.chart_hash&&
        degree.source_registry_hash==registry.wire_hash&&degree.source_epoch==registry.source_epoch&&
        degree.source_steps==registry.source_steps&&degree.source_time==registry.source_time&&
        degree.weight_hash==topology.weight_hash&&degree.law_hash==law.law_hash()&&
        packet_seal.source_registry_hash==registry.wire_hash&&packet_seal.target_chart_hash==topology.chart_hash&&
        packet_seal.weight_hash==topology.weight_hash&&packet_seal.mass_law_hash==law.law_hash()&&
        packet_seal.source_epoch==registry.source_epoch&&packet_seal.source_steps==registry.source_steps&&
        packet_seal.source_time==registry.source_time&&packet_seal.target_cells==registry.site.size()&&
        packet_seal.target_patches==topology.patch.size(),
        "common physical moment packets bind to the same actual degree-one chart, source, weight and fixed-M law");
    std::vector<uint64_t>target_cell_ids,target_patch_cell_ids;
    target_cell_ids.reserve(registry.site.size());target_patch_cell_ids.reserve(topology.patch.size());
    for(const auto&site:registry.site)target_cell_ids.push_back(site.cell);
    for(const auto&patch:topology.patch)target_patch_cell_ids.push_back(patch.cell);
    return StreamAccumulator(expected_seal,packet_seal,field_seal,source_triangle_ids,
        target_cell_ids,target_patch_cell_ids,source_area,source_area_error,
        target_area,target_area_error,budget);
}
Ledger account_actual_overlap_packets(const compact_target::Registry&registry,
        const compact_chart::Certificate&topology,const compact_chart::DegreeOneCertificate&degree,
        const compact_mass_law::BoundLaw&law,const Seal&expected_seal,const Seal&packet_seal,
        const std::vector<Point>&points,const FieldStream&fields,
        const std::vector<uint64_t>&source_triangle_ids,
        const std::vector<R>&source_area,const std::vector<R>&source_area_error,
        const std::vector<R>&target_area,const std::vector<R>&target_area_error,
        const Budget&budget){
    require(points.size()==fields.sample.size(),"one physical field sample per actual overlap packet");
    auto stream=prepare_actual_overlap_stream(registry,topology,degree,law,
        expected_seal,packet_seal,fields.seal,source_triangle_ids,
        source_area,source_area_error,target_area,target_area_error,budget);
    for(size_t i=0;i<points.size();i++)stream.append(points[i],fields.sample[i]);
    return stream.finish();
}
}
namespace compact_target {
Registry register_target_inputs(const compact_source::TargetSource&source,MPI_Comm comm){return gather(source,comm);}
}
namespace compact_mass_law {
BoundLaw bind_accepted_reference_law(const compact_source::Inputs&input,const compact_source::Model&source,
    const compact_source::GeometryBound&geometry,const compact_source::CurrentContext&context,
    const compact_target::Registry&registry,MPI_Comm comm){
    return bind_actual(input,source,geometry,context,registry,comm);
}
}
namespace compact_power {
Roots canonical_exterior_events(const Plane&plane,const compact_source::Exterior&exterior,const RootBudget&budget){return intersections(plane,curve(exterior),budget);}
Curve canonical_actual_source_edge(const compact_source::Model&source,
        const compact_source::GeometryBound&geometry,size_t local_triangle,int local_edge){
    const auto&o=source.owned;
    require(geometry.canonical_global_injective&&geometry.exact_shared_trace&&
        geometry.source_coefficient_hash==source.accepted.coefficient_hash&&
        geometry.source_epoch==source.accepted.physical.epoch&&
        local_triangle<o.f.triangle.size()&&local_edge>=0&&local_edge<3&&
        o.triangle_edge.size()==3*o.f.triangle.size(),
        "actual current source edge from complete accepted shared-trace geometry");
    const auto&e=o.f.triangle[local_triangle].edge[local_edge];
    require(e.poly[1]>=0&&size_t(e.poly[1])<o.global.size(),"original shared source edge polynomial index");
    uint64_t first=o.global[e.poly[1]],edge=o.triangle_edge[3*local_triangle+local_edge];
    require(first>=o.global_vertices&&(first-o.global_vertices)%2==0&&
        edge==(first-o.global_vertices)/2&&edge<o.global_edges,
        "same canonical original edge identity on both incident triangles");
    require(e.poly[0]>=0&&e.poly[3]>=0&&
        size_t(e.poly[0])<o.global.size()&&size_t(e.poly[3])<o.global.size(),
        "actual original source edge global endpoint identities");
    auto curve=source_edge_curve(o.f,local_triangle,local_edge,edge,
        geometry.source_coefficient_hash,geometry.source_epoch);
    curve.vertex[0]=o.global[e.poly[0]];curve.vertex[1]=o.global[e.poly[3]];
    require(curve.vertex[0]<o.global_vertices&&curve.vertex[1]<o.global_vertices&&
        curve.vertex[0]!=curve.vertex[1],
        "distinct canonical original source edge global vertex IDs");
    return curve;
}
Roots canonical_actual_source_edge_events(const Plane&plane,const Curve&edge,
        const compact_source::GeometryBound&geometry,const RootBudget&budget){
    require(edge.source_hash==geometry.source_coefficient_hash&&
        edge.source_epoch==geometry.source_epoch&&geometry.canonical_global_injective,
        "power-plane/source-edge events use the same accepted coefficient epoch");
    return intersections(plane,edge,budget);
}
std::vector<compact_overlap_events::Wire> actual_complete_source_edge_boundary_wires(
        const compact_source::Model&source,const compact_source::GeometryBound&geometry,
        const compact_target::Registry&registry,const compact_overlap_boundary::Index&boundary,
        size_t local_triangle,int local_edge,const RootBudget&root_budget){
    require(source.layout&&source.layout->rank>=0&&source.layout->rank<4&&
        local_triangle<source.owned.f.triangle.size()&&local_edge>=0&&local_edge<3&&
        registry.wire_hash&&boundary.source_registry_hash()==registry.wire_hash&&
        boundary.source_epoch()==registry.source_epoch&&
        registry.source_shard_hash[source.layout->rank]==geometry.source_coefficient_hash&&
        registry.source_epoch==geometry.source_epoch&&registry.source_steps==geometry.source_steps&&
        registry.source_time==geometry.source_time,
        "complete target boundary and local actual-source shard share one accepted epoch");
    auto curve=canonical_actual_source_edge(source,geometry,local_triangle,local_edge);
    auto crossings=boundary.intersections(curve,root_budget);
    require(crossings.complete,"all target-boundary intersections of this source edge resolved");
    const auto&o=source.owned;const auto&t=o.f.triangle[local_triangle];
    const auto&e=t.edge[local_edge];
    require(e.poly[1]>=0&&size_t(e.poly[1])<o.owner.size()&&
        local_triangle<o.source_triangle.size()&&
        ((t.endpoint[local_edge][0]==local_edge&&
          t.endpoint[local_edge][1]==(local_edge+1)%3)||
         (t.endpoint[local_edge][1]==local_edge&&
          t.endpoint[local_edge][0]==(local_edge+1)%3)),
        "original source edge owner, local orientation and triangle identity");
    std::vector<compact_overlap_events::Wire>out;out.reserve(crossings.crossing.size());
    for(const auto&crossing:crossings.crossing){
        compact_overlap_events::Wire wire;wire.event=crossing.event;
        wire.source_registry_hash=registry.wire_hash;
        wire.source_hash=geometry.source_coefficient_hash;wire.source_epoch=geometry.source_epoch;
        wire.target_chart_hash=boundary.chart_hash();wire.weight_hash=boundary.weight_hash();
        wire.target_segment=crossing.boundary_segment;
        wire.target_endpoint=crossing.target_endpoint;
        wire.source_vertex=crossing.source_vertex;
        wire.source_rank=source.layout->rank;wire.incident_triangle=o.source_triangle[local_triangle];
        wire.local_edge=local_edge;
        wire.orientation=t.endpoint[local_edge][0]==local_edge?1:-1;
        wire.edge_owner=o.owner[e.poly[1]];out.push_back(wire);
    }
    return out;
}
struct ActualSourceEdgeBatch {
    std::vector<compact_overlap_events::Wire>wire;
    size_t incidences=0,unique_edges=0,candidate_segments=0,root_calls=0,root_boxes=0;
    bool complete=false;
};
ActualSourceEdgeBatch produce_all_actual_source_edge_boundary_wires(
        const compact_source::Model&source,const compact_source::GeometryBound&geometry,
        const compact_target::Registry&registry,const compact_overlap_boundary::Index&boundary,
        const compact_chart::DegreeOneCertificate&degree,
        const RootBudget&root_budget,const compact_overlap_boundary::ProductionBudget&budget){
    require(source.layout&&source.layout->rank>=0&&source.layout->rank<4&&
        geometry.canonical_global_injective&&geometry.exact_shared_trace&&
        geometry.source_coefficient_hash==source.accepted.coefficient_hash&&
        geometry.source_epoch==source.accepted.physical.epoch&&
        registry.wire_hash&&boundary.source_registry_hash()==registry.wire_hash&&
        boundary.source_epoch()==registry.source_epoch&&
        registry.source_shard_hash[source.layout->rank]==geometry.source_coefficient_hash&&
        registry.source_epoch==geometry.source_epoch&&registry.source_steps==geometry.source_steps&&
        registry.source_time==geometry.source_time&&degree.degree_one_physical_coverage&&
        degree.source_registry_hash==registry.wire_hash&&
        degree.source_epoch==registry.source_epoch&&degree.source_steps==registry.source_steps&&
        degree.source_time==registry.source_time&&degree.weight_hash==boundary.weight_hash()&&
        degree.chart_hash==boundary.chart_hash()&&
        source.owned.external_edge.size()==source.owned.global_edges,
        "complete accepted source and sealed whole target boundary before all-edge event production");
    const auto&o=source.owned;
    require(budget.maximum_global_edges>0&&o.global_edges<=budget.maximum_global_edges&&
        budget.maximum_local_incidences>0&&o.f.triangle.size()<=budget.maximum_local_incidences/3&&
        budget.maximum_local_unique_edges>0&&budget.maximum_total_candidate_segments>0&&
        budget.maximum_total_root_calls>0&&budget.maximum_total_root_boxes>0&&
        budget.maximum_local_wires>0&&
        o.triangle_edge.size()==3*o.f.triangle.size()&&
        o.source_triangle.size()==o.f.triangle.size(),
        "finite complete original source edge incidence and event production budgets");
    struct Slot {uint32_t triangle[2]={UINT32_MAX,UINT32_MAX};uint8_t local_edge[2]={};uint8_t count=0;};
    require(o.f.triangle.size()<=UINT32_MAX,"bounded local original source triangle index");
    std::vector<Slot>slot(o.global_edges);ActualSourceEdgeBatch out;
    for(size_t s=0;s<o.f.triangle.size();s++)for(int k=0;k<3;k++){
        uint64_t edge=o.triangle_edge[3*s+k];
        require(edge<slot.size()&&slot[edge].count<2,
            "every local original source edge has at most two triangle incidences");
        auto&v=slot[edge];v.triangle[v.count]=uint32_t(s);v.local_edge[v.count]=uint8_t(k);v.count++;
        out.incidences++;
    }
    auto charge=[](size_t&total,size_t add,size_t maximum,const char*why){
        require(total<=maximum&&add<=maximum-total,why);total+=add;
    };
    for(size_t edge=0;edge<slot.size();edge++){
        const auto&v=slot[edge];if(!v.count)continue;
        charge(out.unique_edges,1,budget.maximum_local_unique_edges,
            "finite all-source unique edge event budget");
        auto curve=canonical_actual_source_edge(source,geometry,v.triangle[0],v.local_edge[0]);
        require(curve.edge==edge,"complete original edge index and current source curve agree");
        auto roots=boundary.intersections(curve,root_budget);
        require(roots.complete,"all selected target boundary intersections of each source edge resolved");
        charge(out.candidate_segments,roots.candidate_segments,budget.maximum_total_candidate_segments,
            "finite full-source target boundary candidate budget");
        charge(out.root_calls,roots.root_calls,budget.maximum_total_root_calls,
            "finite full-source target plane root-call budget");
        charge(out.root_boxes,roots.root_boxes,budget.maximum_total_root_boxes,
            "finite full-source target plane root-box budget");
        require(out.wire.size()<=budget.maximum_local_wires,
            "finite full-source edge event wire budget");
        require(roots.crossing.size()<=
            (budget.maximum_local_wires-out.wire.size())/v.count,
            "finite complete original triangle incidence wires, no truncation");
        int canonical_owner=-1;
        for(int i=0;i<v.count;i++){
            size_t s=v.triangle[i];int k=v.local_edge[i];
            require(s<o.f.triangle.size()&&o.triangle_edge[3*s+k]==edge&&
                (i==0||o.source_triangle[s]!=o.source_triangle[v.triangle[0]]),
                "distinct complete local original source edge incidences");
            const auto&t=o.f.triangle[s];const auto&e=t.edge[k];
            require(e.poly[1]>=0&&size_t(e.poly[1])<o.owner.size(),
                "canonical original source edge owner for all incidences");
            require((t.endpoint[k][0]==k&&t.endpoint[k][1]==(k+1)%3)||
                    (t.endpoint[k][1]==k&&t.endpoint[k][0]==(k+1)%3),
                "original triangle boundary edge has canonical local orientation");
            if(i==0)canonical_owner=o.owner[e.poly[1]];
            else require(o.owner[e.poly[1]]==canonical_owner,
                "both local incidences route to the same original edge owner");
            for(const auto&crossing:roots.crossing){
                compact_overlap_events::Wire w;w.event=crossing.event;
                w.source_registry_hash=registry.wire_hash;
                w.source_hash=geometry.source_coefficient_hash;w.source_epoch=geometry.source_epoch;
                w.target_chart_hash=boundary.chart_hash();w.weight_hash=boundary.weight_hash();
                w.target_segment=crossing.boundary_segment;
                w.target_endpoint=crossing.target_endpoint;
                w.source_vertex=crossing.source_vertex;
                w.source_rank=source.layout->rank;w.incident_triangle=o.source_triangle[s];
                w.local_edge=k;w.orientation=t.endpoint[k][0]==k?1:-1;
                w.edge_owner=canonical_owner;out.wire.push_back(w);
            }
        }
    }
    require(out.incidences==3*o.f.triangle.size(),
        "every local source triangle contributes all three original edge incidences");
    out.complete=true;return out;
}
std::vector<compact_overlap_events::Canonical> reconcile_actual_source_edge_events(
        const compact_source::Model&source,const compact_source::GeometryBound&geometry,
        const compact_target::Registry&registry,const compact_overlap_boundary::Index&boundary,
        const std::vector<compact_overlap_events::Wire>&local,
        const compact_overlap_events::Budget&budget){
    require(source.layout&&geometry.canonical_global_injective&&geometry.exact_shared_trace&&
        geometry.source_coefficient_hash==source.accepted.coefficient_hash&&
        geometry.source_epoch==source.accepted.physical.epoch&&
        source.owned.external_edge.size()==source.owned.global_edges&&
        source.layout->rank>=0&&source.layout->rank<4&&registry.wire_hash&&
        registry.source_shard_hash[source.layout->rank]==geometry.source_coefficient_hash&&
        registry.source_epoch==geometry.source_epoch&&registry.source_steps==geometry.source_steps&&
        registry.source_time==geometry.source_time&&
        boundary.source_registry_hash()==registry.wire_hash&&
        boundary.source_epoch()==geometry.source_epoch&&boundary.chart_hash()&&
        boundary.weight_hash()&&boundary.segments()>0,
        "complete accepted source edge manifold before global event reconciliation");
    std::vector<compact_overlap_events::Canonical>canonical;bool valid=true;
    try{
        canonical=compact_overlap_events::route_and_reconcile(local,source.owned.external_edge,
            registry.source_shard_hash,registry.wire_hash,geometry.source_epoch,
            boundary.chart_hash(),boundary.weight_hash(),boundary.segments(),budget,source.layout->comm);
    }catch(const std::exception&e){
        fprintf(stderr,"rank%d SOURCE_EDGE_EVENT_RECONCILIATION_FAIL %s\n",source.layout->rank,e.what());
        valid=false;
    }
    compact_mass::collective_require(valid,source.layout->comm,
        "complete original source edge incidences and canonical events on every MPI owner");
    return canonical;
}
struct AllActualSourceEdgeEvents {
    std::vector<compact_overlap_events::Canonical>canonical;
    std::vector<compact_overlap_cuts::Cut>ordered_cut;
    std::vector<compact_overlap_triangle_route::Wire>local_triangle_cut;
};
AllActualSourceEdgeEvents produce_and_reconcile_all_actual_source_edge_events(
        const compact_source::Model&source,const compact_source::GeometryBound&geometry,
        const compact_target::Registry&registry,const compact_overlap_boundary::Index&boundary,
        const compact_chart::DegreeOneCertificate&degree,const RootBudget&root_budget,
        const compact_overlap_boundary::ProductionBudget&production_budget,
        const compact_overlap_events::Budget&merge_budget,
        const compact_overlap_triangle_route::Budget&route_budget){
    require(source.layout!=nullptr,"actual source MPI ownership exists before complete edge-event exchange");
    ActualSourceEdgeBatch local;bool valid=true;
    try{
        local=produce_all_actual_source_edge_boundary_wires(source,geometry,registry,boundary,
            degree,root_budget,production_budget);
        valid=local.complete&&local.wire.size()<=merge_budget.maximum_local_wires;
    }catch(const std::exception&e){
        fprintf(stderr,"rank%d SOURCE_EDGE_EVENT_PRODUCTION_FAIL %s\n",source.layout->rank,e.what());
        valid=false;
    }
    compact_mass::collective_require(valid,source.layout->comm,
        "every accepted source triangle edge has complete bounded target-boundary events on every rank");
    AllActualSourceEdgeEvents out;
    out.canonical=reconcile_actual_source_edge_events(
        source,geometry,registry,boundary,local.wire,merge_budget);
    valid=true;
    std::vector<compact_overlap_triangle_route::Wire>send_triangle_cut;
    std::vector<compact_overlap_triangle_route::OwnedTriangle>owned_triangle;
    try{
        compact_overlap_cuts::Budget cut_budget{merge_budget.maximum_canonical_events,
            merge_budget.maximum_canonical_events,boundary.segments()};
        auto ordered=compact_overlap_cuts::order(boundary,out.canonical,cut_budget);
        valid=ordered.complete;
        out.ordered_cut=std::move(ordered.cut);
        send_triangle_cut=compact_overlap_triangle_route::pack(boundary,out.ordered_cut,
            source.layout->rank,source.layout->ranks,route_budget.maximum_local_wires);
        const auto&o=source.owned;
        require(o.source_triangle.size()==o.f.triangle.size()&&
            o.triangle_edge.size()==3*o.f.triangle.size(),
            "complete original triangle topology for distributed source-edge cuts");
        owned_triangle.reserve(o.f.triangle.size());
        for(size_t s=0;s<o.f.triangle.size();s++){
            compact_overlap_triangle_route::OwnedTriangle t;t.triangle=o.source_triangle[s];
            for(int k=0;k<3;k++){
                const auto&edge=o.f.triangle[s];
                require((edge.endpoint[k][0]==k&&edge.endpoint[k][1]==(k+1)%3)||
                    (edge.endpoint[k][1]==k&&edge.endpoint[k][0]==(k+1)%3),
                    "same actual triangle local boundary orientation after source-edge event merge");
                t.edge[k]=o.triangle_edge[3*s+k];
                t.orientation[k]=edge.endpoint[k][0]==k?1:-1;
                const auto&trace=edge.edge[k];
                require(trace.poly[0]>=0&&trace.poly[3]>=0&&
                    size_t(trace.poly[0])<o.global.size()&&
                    size_t(trace.poly[3])<o.global.size(),
                    "original triangle source edge endpoint global vertices");
                int a=edge.endpoint[k][0],b=edge.endpoint[k][1];
                uint64_t first=o.global[trace.poly[0]],last=o.global[trace.poly[3]];
                require(first<o.global_vertices&&last<o.global_vertices&&first!=last&&
                    (t.vertex[a]==UINT64_MAX||t.vertex[a]==first)&&
                    (t.vertex[b]==UINT64_MAX||t.vertex[b]==last),
                    "same canonical original vertex ID on both incident triangle edges");
                t.vertex[a]=first;t.vertex[b]=last;
            }
            owned_triangle.push_back(t);
        }
    }catch(const std::exception&e){
        fprintf(stderr,"rank%d SOURCE_EDGE_CUT_ORDER_FAIL %s\n",source.layout->rank,e.what());
        valid=false;
    }
    compact_mass::collective_require(valid,source.layout->comm,
        "ordered source edge cuts and complete exterior target-node stars on every MPI owner");
    auto routed=compact_overlap_triangle_route::route(boundary,send_triangle_cut,owned_triangle,
        registry.source_shard_hash,route_budget,source.layout->comm);
    require(routed.complete,"complete ordered original source triangle boundary cuts");
    out.local_triangle_cut=std::move(routed.triangle_cut);
    return out;
}
BoundedCell enclose_power_cell(const Spatial&sites,size_t index,const CellBudget&budget){return bounded_cell(sites,index,budget);}
}
