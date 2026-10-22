#ifndef MATERIAL_COMPACT_SOURCE_GEOMETRY_H
#define MATERIAL_COMPACT_SOURCE_GEOMETRY_H
// Native TU provides original st alias/need, then includes compact_source.h.
#include "material_shared_real_bound.h"
#include <atomic>

namespace compact_source {
namespace sf=shared_function;
struct EdgeWire {
    uint64_t edge=0,owner=0,rank=0,incidence=0,external=0,dof[6]={};
    int64_t orientation=0;
    double parameter[21]={};
};
inline bool same_edge(const EdgeWire&a,const EdgeWire&b){
    if(a.edge!=b.edge||a.owner!=b.owner||a.rank!=b.rank||a.external!=b.external)return false;
    for(int j=0;j<6;j++)if(a.dof[j]!=b.dof[j])return false;
    for(int j=0;j<21;j++)if(a.parameter[j]!=b.parameter[j])return false;
    return true;
}
inline uint64_t check_shared_edges(const Model&m){
    const auto&o=m.owned;const auto&layout=*m.layout;std::map<uint64_t,EdgeWire>local;
    for(const auto&t:o.f.triangle)for(int k=0;k<3;k++){
        const auto&e=t.edge[k];require(e.rank==1||e.rank==2,"original canonical edge rank, no truncated direction");
        for(int j=0;j<4;j++)require(e.poly[j]>=0&&size_t(e.poly[j])<o.global.size(),"original shared polynomial trace indices");
        uint64_t first=o.global[e.poly[1]];
        require(first>=o.global_vertices&&(first-o.global_vertices)%2==0&&o.global[e.poly[2]]==first+1,"canonical global edge polynomial identity");
        EdgeWire w;w.edge=(first-o.global_vertices)/2;require(w.edge<o.global_edges,"canonical source edge range");
        w.owner=o.owner[e.poly[1]];w.rank=e.rank;w.external=o.external_edge.at(w.edge);
        for(int j=0;j<4;j++)w.dof[j]=o.global[e.poly[j]];
        for(int j=0;j<2;j++){
            w.dof[j+4]=UINT64_MAX;
            if(j<e.rank){require(e.id[j]>=0&&size_t(e.id[j])<o.global.size(),"canonical enrichment trace coefficient");w.dof[j+4]=o.global[e.id[j]];}
        }
        require(w.dof[0]<w.dof[3]&&w.dof[3]<o.global_vertices,"canonical oriented shared reference endpoints");
        size_t at=0;for(double z:e.A)w.parameter[at++]=z;for(double z:e.D)w.parameter[at++]=z;
        for(const auto&row:e.qnode)for(double z:row)w.parameter[at++]=z;w.parameter[at++]=e.eps;
        for(const auto&row:e.encode)for(double z:row)w.parameter[at++]=z;
        for(const auto&row:e.decode)for(double z:row)w.parameter[at++]=z;
        require(at==21,"complete canonical trace function descriptor");for(double z:w.parameter)require(std::isfinite(z),"finite actual canonical edge function");
        int a=t.endpoint[k][0],b=t.endpoint[k][1];require((a==k&&b==(k+1)%3)||(b==k&&a==(k+1)%3),"original trace endpoint chart");
        for(int d=0;d<2;d++)require(e.A[d]==t.vertex[a][d]&&e.D[d]==t.vertex[b][d]-t.vertex[a][d],"same reference edge function on each incident triangle");
        w.incidence=1;w.orientation=a==k?1:-1;
        auto entry=local.emplace(w.edge,w);
        if(!entry.second){require(same_edge(entry.first->second,w),"exact local shared trace function, no approximate gluing");entry.first->second.incidence++;entry.first->second.orientation+=w.orientation;}
    }
    std::vector<std::vector<EdgeWire>>send(layout.ranks),receive(layout.ranks);
    for(const auto&entry:local)send.at(entry.second.owner).push_back(entry.second);
    std::vector<uint64_t>counts(layout.ranks),incoming(layout.ranks);
    for(int p=0;p<layout.ranks;p++)counts[p]=send[p].size();
    MPI_Alltoall(counts.data(),1,MPI_UINT64_T,incoming.data(),1,MPI_UINT64_T,layout.comm);
    for(int p=0;p<layout.ranks;p++){require(incoming[p]<=o.global_edges,"bounded actual edge incidence messages");receive[p].resize(incoming[p]);}
    receive[layout.rank]=send[layout.rank];
    MPI_Datatype wire;require(MPI_Type_contiguous(sizeof(EdgeWire),MPI_BYTE,&wire)==MPI_SUCCESS&&MPI_Type_commit(&wire)==MPI_SUCCESS,"canonical edge descriptor MPI type");
    md::exchange(send,receive,wire,layout.comm,layout.rank,340);
    require(MPI_Type_free(&wire)==MPI_SUCCESS,"canonical edge MPI type released");
    std::map<uint64_t,EdgeWire>owned;
    for(const auto&peer:receive)for(const auto&w:peer){
        require(w.owner==uint64_t(layout.rank),"same unique canonical edge owner");auto entry=owned.emplace(w.edge,w);
        if(!entry.second){require(same_edge(entry.first->second,w),"exact cross-rank shared trace function and global DOFs");entry.first->second.incidence+=w.incidence;entry.first->second.orientation+=w.orientation;}
    }
    uint64_t summary[3]={owned.size(),0,0};
    for(const auto&entry:owned){const auto&e=entry.second;require(e.incidence==(e.external?1:2)&&e.orientation==(e.external?(e.orientation>0?1:-1):0),"same oriented reference manifold and exterior registry");summary[1]+=e.incidence;summary[2]+=e.external;}
    MPI_Allreduce(MPI_IN_PLACE,summary,3,MPI_UINT64_T,MPI_SUM,layout.comm);
    require(summary[0]==o.global_edges&&summary[1]==3*o.global_triangles&&summary[2]==2*o.global_edges-3*o.global_triangles,"complete global canonical trace incidence, not sampled C0");
    return summary[2];
}
struct GeometryBound {
    std::vector<sf::RealBounds>triangle;
    std::vector<uint8_t>complete;
    uint64_t unresolved=0,first_unresolved=UINT64_MAX,external_edges=0;
    uint64_t source_coefficient_hash=0,source_epoch=0,source_steps=0;
    double source_time=0;
    R maximum_Dchi_minus_I=0,minimum_canonical_J=0;
    bool exact_shared_trace=false,canonical_global_injective=false;
    // Original convex reference construction + existing <=8ulp elementary
    // function engineering model are explicit premises of this certificate.
    // It does not certify arbitrary native libm or a target domain/remap.
};
inline GeometryBound bound_source_geometry(const Model&m){
    require(m.layout&&m.original_reference_receipt_bound,"actual full source and immutable original square-reference construction required");
    GeometryBound result;result.external_edges=check_shared_edges(m);result.exact_shared_trace=true;
    result.source_coefficient_hash=m.accepted.coefficient_hash;result.source_epoch=m.accepted.physical.epoch;
    result.source_steps=m.accepted.physical.steps;result.source_time=m.accepted.physical.time;
    require(compact_p3::bytes_hash(m.owned.f.coefficient.data(),m.owned.f.coefficient.size()*sizeof(STCoefficient))==result.source_coefficient_hash,"geometry enclosure is bound to the actual accepted coefficient bits");
    const auto&o=m.owned;result.triangle.resize(o.f.triangle.size());result.complete.resize(o.f.triangle.size());
    std::vector<sf::Q>lower_J(o.f.triangle.size());
    #pragma omp parallel for schedule(dynamic,32)
    for(size_t s=0;s<o.f.triangle.size();s++)try{
        const auto&t=o.f.triangle[s];for(const auto&v:t.vertex)for(double x:v)require(std::isfinite(x)&&fabs(x)<=12,"same original convex reference square");
        auto map=sf::compile(o.f,s);sf::RealBounds bound;
        try{bound=sf::tight_real_bounds(map);}catch(const std::exception&){bound=sf::real_bounds(map);}
        require(finiteq(bound.L)&&bound.L>=0,"finite whole-enriched-map derivative enclosure");result.triangle[s]=bound;
        using S=sf::S;
        auto det=(S(1)-S(bound.FminusI[0]))*(S(1)-S(bound.FminusI[3]))-S(bound.FminusI[1])*S(bound.FminusI[2]);
        lower_J[s]=sf::lower(det);
        result.complete[s]=bound.L<1&&lower_J[s]>0;
    }catch(const std::exception&){result.complete[s]=0;}
    R maximum=0,minimum=std::numeric_limits<R>::infinity();
    for(size_t s=0;s<o.f.triangle.size();s++){
        if(!result.complete[s]){result.unresolved++;result.first_unresolved=std::min(result.first_unresolved,o.source_triangle[s]);continue;}
        maximum=std::max(maximum,nextafterl(R(result.triangle[s].L),INFINITY));
        minimum=std::min(minimum,nextafterl(R(lower_J[s]),-INFINITY));
    }
    MPI_Allreduce(MPI_IN_PLACE,&result.unresolved,1,MPI_UINT64_T,MPI_SUM,m.layout->comm);
    MPI_Allreduce(MPI_IN_PLACE,&result.first_unresolved,1,MPI_UINT64_T,MPI_MIN,m.layout->comm);
    MPI_Allreduce(MPI_IN_PLACE,&maximum,1,MPI_LONG_DOUBLE,MPI_MAX,m.layout->comm);
    MPI_Allreduce(MPI_IN_PLACE,&minimum,1,MPI_LONG_DOUBLE,MPI_MIN,m.layout->comm);
    if(result.unresolved){maximum=std::numeric_limits<R>::infinity();minimum=-std::numeric_limits<R>::infinity();}
    result.maximum_Dchi_minus_I=maximum;result.minimum_canonical_J=minimum;
    // On this already verified convex reference disk, exact shared traces and
    // global sup||Dchi-I||_infinity<1 give ||chi(X)-chi(Y)||_infinity >=
    // (1-L)||X-Y||_infinity. This is a sufficient CURRENT-state criterion,
    // not a restriction on permanent rigid rotation or long mesh lifetime.
    result.canonical_global_injective=!result.unresolved&&maximum<1&&minimum>0;
    return result;
}
GeometryBound certify_accepted_geometry(const Model&);
} // namespace compact_source
#endif
