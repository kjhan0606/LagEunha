#ifndef MATERIAL_COMPACT_OVERLAP_EVENTS_H
#define MATERIAL_COMPACT_OVERLAP_EVENTS_H
// Reconcile duplicate current-source edge / Laguerre-plane roots from the
// one or two original source-triangle incidences. Geometry is never inferred
// from a single incidence or an approximately equal root coordinate.
#include "material_compact_power_curve.h"
#include <mpi.h>
#include <climits>
#include <cstring>
#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

namespace compact_overlap_events {
namespace cp=compact_power;
using Q=cp::Q;
inline void require(bool ok,const char*why){if(!ok)throw std::runtime_error(why);}
struct Wire {
    cp::Event event;
    uint64_t source_registry_hash=0,source_hash=0,source_epoch=0,incident_triangle=UINT64_MAX;
    uint64_t target_chart_hash=0,weight_hash=0,target_segment=UINT64_MAX;
    uint64_t source_vertex=UINT64_MAX;
    int source_rank=-1,edge_owner=-1,local_edge=-1,orientation=0;
    int target_endpoint=-1;
};
struct Canonical {
    cp::Event event;
    uint64_t source_registry_hash=0,source_epoch=0,target_chart_hash=0,weight_hash=0;
    uint64_t target_segment=UINT64_MAX;
    uint64_t source_vertex=UINT64_MAX;
    uint64_t source_hash[2]={},incident_triangle[2]={UINT64_MAX,UINT64_MAX};
    int source_rank[2]={-1,-1};
    int local_edge[2]={-1,-1},orientation[2]={};
    int edge_owner=-1,incidences=0;
    int target_endpoint=-1;
};
struct Budget {
    uint64_t maximum_local_wires=0,maximum_received_wires=0,maximum_canonical_events=0;
};
static_assert(std::is_trivially_copyable<Wire>::value,"fixed bytewise source-edge MPI event wire");
inline auto key(const Wire&w){return std::make_tuple(w.event.edge,w.event.cell[0],w.event.cell[1],w.event.ordinal,w.target_segment);}
inline std::vector<Canonical> reconcile(std::vector<Wire>wires,
        const std::vector<uint8_t>&external_edge,int rank,int ranks,
        const std::array<uint64_t,4>&expected_shard_hash,uint64_t expected_registry_hash,
        uint64_t expected_epoch,uint64_t expected_chart_hash,uint64_t expected_weight_hash,
        uint64_t expected_segments,const Budget&budget){
    bool shards_valid=ranks>0&&ranks<=4;
    for(int p=0;p<ranks&&p<4;p++)shards_valid&=expected_shard_hash[p]!=0;
    require(shards_valid&&rank>=0&&rank<ranks&&expected_registry_hash&&expected_epoch&&
        expected_chart_hash&&expected_weight_hash&&expected_segments>0&&
        budget.maximum_received_wires>0&&budget.maximum_canonical_events>0&&
        wires.size()<=budget.maximum_received_wires&& !external_edge.empty(),
        "finite complete source-edge event owner registry and budget");
    std::sort(wires.begin(),wires.end(),[](const Wire&a,const Wire&b){
        return key(a)<key(b)||(key(a)==key(b)&&a.incident_triangle<b.incident_triangle);
    });
    std::vector<Canonical>out;
    for(size_t i=0;i<wires.size();){
        size_t j=i+1;while(j<wires.size()&&key(wires[j])==key(wires[i]))j++;
        const Wire&first=wires[i];const auto&e=first.event;
        require(e.edge<external_edge.size()&&external_edge[e.edge]<=1&&e.cell[0]<e.cell[1]&&
            first.target_segment<expected_segments&&
            first.target_endpoint>=-1&&first.target_endpoint<=1&&
            (first.source_vertex==UINT64_MAX||
                (first.target_endpoint>=0&&e.lower==e.upper&&
                 (e.lower==0||e.lower==1)))&&
            e.lower>=0&&e.lower<=e.upper&&e.upper<=1&&
            (e.derivative_sign==1||e.derivative_sign==-1),
            "valid oriented original-edge/power-plane root interval");
        size_t need_count=external_edge[e.edge]?1:2;
        require(j-i==need_count,"every accepted source edge event has its complete original triangle incidence");
        require(out.size()<budget.maximum_canonical_events,"finite canonical source-edge event budget");
        Canonical canonical;canonical.event=e;canonical.source_registry_hash=expected_registry_hash;
        canonical.source_epoch=expected_epoch;canonical.target_chart_hash=expected_chart_hash;
        canonical.weight_hash=expected_weight_hash;canonical.target_segment=first.target_segment;
        canonical.target_endpoint=first.target_endpoint;
        canonical.source_vertex=first.source_vertex;
        canonical.edge_owner=rank;
        canonical.incidences=int(need_count);
        for(size_t k=i;k<j;k++){
            const auto&w=wires[k];const auto&a=w.event;
            require(w.source_rank>=0&&w.source_rank<ranks&&
                w.local_edge>=0&&w.local_edge<3&&
                (w.orientation==1||w.orientation==-1)&&
                w.source_hash==expected_shard_hash[w.source_rank]&&
                w.source_registry_hash==expected_registry_hash&&w.source_epoch==expected_epoch&&
                w.target_chart_hash==expected_chart_hash&&w.weight_hash==expected_weight_hash&&
                w.target_segment==first.target_segment&&
                w.target_endpoint==first.target_endpoint&&
                w.source_vertex==first.source_vertex&&
                w.edge_owner==rank&&w.incident_triangle!=UINT64_MAX&&
                a.edge==e.edge&&a.cell[0]==e.cell[0]&&a.cell[1]==e.cell[1]&&
                a.ordinal==e.ordinal&&a.derivative_sign==e.derivative_sign&&
                a.lower<=a.upper&&a.lower>=0&&a.upper<=1,
                "both edge incidences use identical source epoch, owner, plane, ordinal and derivative orientation");
            if(k>i)require(wires[k-1].incident_triangle<w.incident_triangle,
                "distinct ordered original source triangle incidences");
            canonical.incident_triangle[k-i]=w.incident_triangle;
            canonical.source_rank[k-i]=w.source_rank;
            canonical.source_hash[k-i]=w.source_hash;
            canonical.local_edge[k-i]=w.local_edge;
            canonical.orientation[k-i]=w.orientation;
            canonical.event.lower=std::max(canonical.event.lower,a.lower);
            canonical.event.upper=std::min(canonical.event.upper,a.upper);
        }
        require(canonical.event.lower<=canonical.event.upper,
            "incident root enclosures intersect for the same canonical event");
        require(need_count==1||canonical.orientation[0]+canonical.orientation[1]==0,
            "two original source triangles traverse their shared oriented edge oppositely");
        out.push_back(canonical);i=j;
    }
    return out;
}
inline std::vector<Canonical> route_and_reconcile(const std::vector<Wire>&local,
        const std::vector<uint8_t>&external_edge,
        const std::array<uint64_t,4>&expected_shard_hash,uint64_t expected_registry_hash,
        uint64_t expected_epoch,uint64_t expected_chart_hash,uint64_t expected_weight_hash,
        uint64_t expected_segments,const Budget&budget,MPI_Comm comm){
    int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
    int local_valid=ranks>0&&ranks<=4&&budget.maximum_local_wires>0&&
        local.size()<=budget.maximum_local_wires&&expected_registry_hash&&expected_epoch&&
        expected_chart_hash&&expected_weight_hash&&expected_segments>0;
    for(int p=0;p<ranks&&p<4;p++)local_valid&=expected_shard_hash[p]!=0;
    for(const auto&w:local)local_valid&=rank>=0&&rank<4&&w.source_rank==rank&&
        w.source_hash==expected_shard_hash[rank]&&
        w.source_registry_hash==expected_registry_hash&&w.source_epoch==expected_epoch&&
        w.target_chart_hash==expected_chart_hash&&w.weight_hash==expected_weight_hash&&
        w.target_segment<expected_segments&&w.target_endpoint>=-1&&w.target_endpoint<=1&&
        (w.source_vertex==UINT64_MAX||
            (w.target_endpoint>=0&&w.event.lower==w.event.upper&&
             (w.event.lower==0||w.event.lower==1)))&&
        w.local_edge>=0&&w.local_edge<3&&
        (w.orientation==1||w.orientation==-1)&&
        w.edge_owner>=0&&w.edge_owner<ranks;
    int all_valid=0;
    require(MPI_Allreduce(&local_valid,&all_valid,1,MPI_INT,MPI_MIN,comm)==MPI_SUCCESS,
        "collective source-edge event preflight");
    require(all_valid,"finite locally generated source-edge events and valid edge owners on every rank");
    std::array<uint64_t,10> seal={expected_registry_hash,expected_epoch,
        expected_chart_hash,expected_weight_hash,expected_segments,
        expected_shard_hash[0],expected_shard_hash[1],expected_shard_hash[2],
        expected_shard_hash[3],uint64_t(external_edge.size())};
    std::vector<uint64_t>all_seals(seal.size()*ranks);
    require(MPI_Allgather(seal.data(),int(seal.size()),MPI_UINT64_T,
        all_seals.data(),int(seal.size()),MPI_UINT64_T,comm)==MPI_SUCCESS,
        "complete replicated source/target event manifest exchange");
    for(int p=0;p<ranks;p++)for(size_t j=0;j<seal.size();j++)
        require(all_seals[p*seal.size()+j]==seal[j],
            "all ranks use the same source, target chart, trial weights and segment registry");
    std::vector<std::vector<Wire>>send(ranks);
    for(const auto&w:local)send[w.edge_owner].push_back(w);
    std::vector<uint64_t>send_count(ranks),recv_count(ranks);
    for(int p=0;p<ranks;p++)send_count[p]=send[p].size();
    require(MPI_Alltoall(send_count.data(),1,MPI_UINT64_T,recv_count.data(),1,MPI_UINT64_T,comm)==MPI_SUCCESS,
        "complete source-edge event peer counts");
    uint64_t send_total=0,recv_total=0;
    std::vector<int>send_bytes(ranks),recv_bytes(ranks),send_offset(ranks),recv_offset(ranks);
    local_valid=1;
    for(int p=0;p<ranks;p++){
        if(send_count[p]>uint64_t(INT_MAX/sizeof(Wire))||
            recv_count[p]>uint64_t(INT_MAX/sizeof(Wire))){local_valid=0;continue;}
        uint64_t sb=send_count[p]*sizeof(Wire),rb=recv_count[p]*sizeof(Wire);
        if(send_total+sb>INT_MAX||recv_total+rb>INT_MAX){local_valid=0;continue;}
        send_offset[p]=int(send_total);recv_offset[p]=int(recv_total);
        send_bytes[p]=int(sb);recv_bytes[p]=int(rb);
        send_total+=sb;recv_total+=rb;
    }
    local_valid&=recv_total/sizeof(Wire)<=budget.maximum_received_wires;
    require(MPI_Allreduce(&local_valid,&all_valid,1,MPI_INT,MPI_MIN,comm)==MPI_SUCCESS,
        "collective source-edge MPI byte budget preflight");
    require(all_valid,"finite source-edge MPI bytes and received event budget on every rank");
    std::vector<Wire>packed(send_total/sizeof(Wire)),received(recv_total/sizeof(Wire));
    for(int p=0;p<ranks;p++)if(!send[p].empty())
        std::memcpy(reinterpret_cast<char*>(packed.data())+send_offset[p],send[p].data(),send_bytes[p]);
    require(MPI_Alltoallv(packed.data(),send_bytes.data(),send_offset.data(),MPI_BYTE,
        received.data(),recv_bytes.data(),recv_offset.data(),MPI_BYTE,comm)==MPI_SUCCESS,
        "complete source-edge event owner routing");
    return reconcile(std::move(received),external_edge,rank,ranks,
        expected_shard_hash,expected_registry_hash,expected_epoch,
        expected_chart_hash,expected_weight_hash,expected_segments,budget);
}
}
#endif
