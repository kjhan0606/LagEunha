#ifndef MATERIAL_COMPACT_MASS_H
#define MATERIAL_COMPACT_MASS_H
// Adapter to the ORIGINAL MPI owned/ghost, binary80 outer equation and DD
// GPU/resident-PCG APIs. Local geometry-mode blocks may be semidefinite even
// when the assembled operator is positive. Never regularize the physical M.
#include "material_distributed_mass.h"
#include "material_compact_geometry.h"
#include <memory>
#include <climits>
#include <cstring>
namespace compact_mass {
namespace md=material_distributed;
using R=long double;
struct Component {uint64_t id=0,global_mode[2]={};};
struct LocalElement {
    uint64_t component=UINT64_MAX;
    int n=10,id[12]={};
    std::array<R,144>mass={};
};
struct Prepared {
    std::vector<STTriangle>triangle; // IDs/n only; NEVER an old physics basis.
    std::unique_ptr<md::Mass>exact;
    std::vector<std::array<R,3>>component_mode_metric;
    std::vector<uint64_t>component_support;
    // exact owns references to triangle and external Ownership. Keep both alive
    // and destroy exact before MPI_Finalize. Moving this object would invalidate
    // the triangle reference; use only a heap-owned Prepared, never by value.
    Prepared()=default;Prepared(const Prepared&)=delete;Prepared&operator=(const Prepared&)=delete;
    Prepared(Prepared&&)=delete;Prepared&operator=(Prepared&&)=delete;
};
inline void collective_require(bool valid,MPI_Comm comm,const char*message){
    int bad=valid?0:1,any=0;MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,comm);md::require(!any,message);
}
inline void bind_components(const std::vector<Component>&components,const md::Ownership&layout){
    uint64_t count=components.size(),root_count=count;MPI_Bcast(&root_count,1,MPI_UINT64_T,0,layout.comm);
    collective_require(count==root_count&&count<=size_t(INT_MAX/3),layout.comm,"same bounded actual geometry-component registry");
    std::vector<uint64_t>wire(3*count);if(layout.rank==0)for(size_t j=0;j<count;j++){
        wire[3*j]=components[j].id;wire[3*j+1]=components[j].global_mode[0];wire[3*j+2]=components[j].global_mode[1];
    }
    MPI_Bcast(wire.data(),int(wire.size()),MPI_UINT64_T,0,layout.comm);bool same=true;std::set<uint64_t>ids;
    for(size_t j=0;j<count;j++){
        const auto&c=components[j];same&=c.id==wire[3*j]&&c.global_mode[0]==wire[3*j+1]&&c.global_mode[1]==wire[3*j+2];
        same&=c.global_mode[0]!=c.global_mode[1]&&(j==0||components[j-1].id<c.id);
        same&=ids.insert(c.global_mode[0]).second&&ids.insert(c.global_mode[1]).second;
    }
    collective_require(same,layout.comm,"exact shared two-mode component IDs, not a hash-only registry");
}
inline std::unique_ptr<Prepared> prepare(const md::Ownership&layout,const std::vector<LocalElement>&element,
    const std::vector<Component>&components){
    bind_components(components,layout);std::map<uint64_t,size_t>component;std::set<uint64_t>mode_ids;
    for(size_t j=0;j<components.size();j++){
        component.emplace(components[j].id,j);for(auto id:components[j].global_mode)mode_ids.insert(id);
    }
    auto out=std::make_unique<Prepared>();out->component_mode_metric.resize(components.size());out->component_support.resize(components.size());
    bool valid=element.size()<=UINT32_MAX&&layout.global.size()<=INT_MAX;
    for(const auto&e:element){
        bool local=e.n==10||e.n==12;std::set<int>ids;
        for(int a=0;a<e.n&&a<12;a++)local&=e.id[a]>=0&&size_t(e.id[a])<layout.global.size()&&ids.insert(e.id[a]).second;
        if(!local){valid=false;continue;}
        for(int a=0;a<10;a++)local&=!mode_ids.count(layout.global[e.id[a]]);
        for(int a=0;a<e.n;a++)for(int b=0;b<e.n;b++)local&=std::isfinite(e.mass[a*12+b])&&e.mass[a*12+b]==e.mass[b*12+a];
        for(int a=0;a<e.n;a++)local&=e.mass[a*12+a]>=0;
        if(e.n==12){
            auto c=component.find(e.component);local&=c!=component.end();
            if(c!=component.end()){
                size_t j=c->second;
                for(int k=0;k<2;k++)local&=layout.global[e.id[10+k]]==components[j].global_mode[k];
                out->component_mode_metric[j][0]+=e.mass[10*12+10];
                out->component_mode_metric[j][1]+=e.mass[10*12+11];
                out->component_mode_metric[j][2]+=e.mass[11*12+11];out->component_support[j]++;
            }
        }
        valid&=local;
    }
    collective_require(valid,layout.comm,"finite symmetric physical compact blocks and actual local mode ownership");
    static_assert(sizeof(std::array<R,3>)==3*sizeof(R),"contiguous three-entry component metric");
    MPI_Allreduce(MPI_IN_PLACE,out->component_mode_metric.data(),int(3*components.size()),MPI_LONG_DOUBLE,MPI_SUM,layout.comm);
    MPI_Allreduce(MPI_IN_PLACE,out->component_support.data(),int(components.size()),MPI_UINT64_T,MPI_SUM,layout.comm);
    valid=true;for(size_t j=0;j<components.size();j++){
        const auto&C=out->component_mode_metric[j];
        // Necessary finite-functional rank check only. Actual geometry-function
        // rank and continuous quadrature positivity still need their own proof.
        valid&=out->component_support[j]>0&&C[0]>0&&C[2]>0&&std::isfinite(C[0])&&std::isfinite(C[1])&&std::isfinite(C[2]);
        valid&=C[0]*C[2]-C[1]*C[1]>0;
    }
    collective_require(valid,layout.comm,"positive assembled two-mode metric, no tolerance-based mode removal");
    out->triangle.resize(element.size());md::Vec physical(element.size()*size_t(ST_MAX)*ST_MAX),surrogate(physical.size());
    for(size_t s=0;s<element.size();s++){
        const auto&e=element[s];auto&t=out->triangle[s];t={};t.n=e.n;for(int a=0;a<e.n;a++)t.id[a]=e.id[a];
        size_t at=s*ST_MAX*ST_MAX;
        for(int a=0;a<e.n;a++)for(int b=0;b<e.n;b++)physical[at+a*ST_MAX+b]=e.mass[a*12+b];
        // Only the PRECONDITIONER uses this block diagonal surrogate. Every
        // bulk--mode and mode--mode physical entry remains in physical above.
        for(int a=0;a<10;a++)for(int b=0;b<10;b++)surrogate[at+a*ST_MAX+b]=e.mass[a*12+b];
        if(e.n==12){
            size_t j=component.at(e.component);const auto&C=out->component_mode_metric[j];R n=out->component_support[j];
            surrogate[at+10*ST_MAX+10]=C[0]/n;surrogate[at+10*ST_MAX+11]=surrogate[at+11*ST_MAX+10]=C[1]/n;
            surrogate[at+11*ST_MAX+11]=C[2]/n;
        }
    }
    // Original constructor builds the support registry/partition and SPD local
    // factors from surrogate. Before ANY apply/GPU upload, replace ONLY block
    // with the authoritative physical operator. Factors remain preconditioners.
    out->exact=std::make_unique<md::Mass>(layout,out->triangle,std::move(surrogate));
    out->exact->block=std::move(physical);
    return out;
}
} // namespace compact_mass
#endif
