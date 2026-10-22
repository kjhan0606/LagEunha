#ifndef MATERIAL_COMPACT_CHART_DEGREE_H
#define MATERIAL_COMPACT_CHART_DEGREE_H
// Degree-one coverage from a certified disk chart, positive local Jacobians,
// and the globally injective, positively oriented accepted-source boundary.
#include "material_compact_chart_jacobian.h"

namespace compact_chart {
struct DegreeOneCertificate {
    uint64_t source_registry_hash=0,weight_hash=0,source_epoch=0,source_steps=0,law_hash=0,chart_hash=0;
    double source_time=0;
    size_t chart_patches=0,vertex_stars=0,boundary_edges=0;
    int64_t euler_characteristic=0;
    unsigned boundary_degree=0;
    bool source_boundary_simple=false,chart_is_disk=false,chart_J_positive=false;
    bool continuous_shared_edges=false,degree_one_physical_coverage=false;
};

inline DegreeOneCertificate certify_degree_one(const compact_target::Registry&registry,
        const std::vector<Candidate>&candidate,const Certificate&topology,
        const JacobianCertificate&jacobian,const compact_source::GeometryBound&source_geometry){
    require_same_chart_complex(registry,candidate,topology);
    require(topology.exact_C0_connectivity&&topology.oriented_source_boundary_complete&&topology.disk_topology&&
        topology.manifold_vertex_stars&&topology.global_components==1&&topology.euler_characteristic==1,
        "global chart is one connected oriented topological disk with a complete original boundary cycle");
    require(source_geometry.canonical_global_injective&&source_geometry.exact_shared_trace&&
        source_geometry.source_epoch==registry.source_epoch&&source_geometry.source_steps==registry.source_steps&&
        source_geometry.source_time==registry.source_time&&source_geometry.external_edges==registry.exterior.size()&&
        source_geometry.maximum_Dchi_minus_I<1&&source_geometry.minimum_canonical_J>0&&
        std::find(registry.source_shard_hash.begin(),registry.source_shard_hash.end(),source_geometry.source_coefficient_hash)!=registry.source_shard_hash.end(),
        "same accepted source geometry has a globally injective, positively oriented, shared-trace boundary map");
    require(jacobian.all_reference_J_positive&&jacobian.source_registry_hash==registry.wire_hash&&
        jacobian.weight_hash==topology.weight_hash&&jacobian.source_epoch==registry.source_epoch&&
        jacobian.source_steps==registry.source_steps&&jacobian.source_time==registry.source_time&&
        jacobian.chart_hash==topology.chart_hash&&jacobian.patch_lower.size()==topology.patch.size()&&jacobian.minimum_lower>0,
        "positive whole-reference Jacobians certify this exact source/W/chart set");

    // The map is continuous across canonical shared straight edges and follows
    // the same accepted source curve around the exterior exactly once. Its
    // domain is a disk, and every patch preserves orientation. The boundary
    // is the simple positive image of the original square. The planar degree
    // is therefore one inside that boundary and zero outside; positive local
    // degrees make the multiplicity exactly one for every regular interior
    // value and exclude overlap or missing open regions.
    DegreeOneCertificate out;out.source_registry_hash=registry.wire_hash;out.weight_hash=topology.weight_hash;
    out.source_epoch=registry.source_epoch;out.source_steps=registry.source_steps;out.source_time=registry.source_time;
    out.law_hash=topology.law_hash;out.chart_hash=topology.chart_hash;out.chart_patches=topology.patch.size();
    out.vertex_stars=topology.vertex_stars;out.boundary_edges=topology.boundary_edges;
    out.euler_characteristic=topology.euler_characteristic;out.boundary_degree=1;
    out.source_boundary_simple=true;out.chart_is_disk=true;out.chart_J_positive=true;
    out.continuous_shared_edges=topology.exact_C0_connectivity;out.degree_one_physical_coverage=true;
    return out;
}
}
#endif
