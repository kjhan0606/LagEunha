#ifndef MATERIAL_COMPACT_OVERLAP_TRIANGLE_ROUTE_H
#define MATERIAL_COMPACT_OVERLAP_TRIANGLE_ROUTE_H
// Route complete source-edge cuts from canonical edge owners back to the
// original triangle owners. The receiver reconstructs target-node stars from
// its sealed chart and orders cuts along each local triangle boundary edge.
// This supplies boundary incidences; it does not construct overlap regions.
#include "material_compact_overlap_cuts.h"
#include <climits>
#include <cstdio>
#include <cstring>
#include <type_traits>

namespace compact_overlap_triangle_route {
namespace oc=compact_overlap_cuts;
namespace cc=compact_chart;
using Q=compact_power::Q;
inline void require(bool ok,const char*why){if(!ok)throw std::runtime_error(why);}

struct OwnedTriangle {
    uint64_t triangle=UINT64_MAX,edge[3]={UINT64_MAX,UINT64_MAX,UINT64_MAX};
    uint64_t vertex[3]={UINT64_MAX,UINT64_MAX,UINT64_MAX};
    int orientation[3]={};
};
struct Wire {
    uint64_t source_registry_hash=0,source_epoch=0,target_chart_hash=0,weight_hash=0;
    uint64_t source_edge=UINT64_MAX,source_triangle=UINT64_MAX,source_hash=0;
    uint64_t source_vertex=UINT64_MAX;
    uint64_t target_segment=UINT64_MAX,star_count=0,star_hash=0;
    Q lower=0,upper=0;
    int source_rank=-1,edge_owner=-1,local_edge=-1,orientation=0,target_endpoint=-1;
};
static_assert(std::is_trivially_copyable<Wire>::value,"fixed bytewise triangle cut routing wire");
struct Budget {
    uint64_t maximum_local_wires=0,maximum_received_wires=0;
};
struct Result {
    struct EdgeRange {
        uint64_t triangle=UINT64_MAX,source_edge=UINT64_MAX;
        uint64_t start_vertex=UINT64_MAX,end_vertex=UINT64_MAX;
        int local_edge=-1,orientation=0;
        size_t begin=0,end=0;
    };
    std::vector<Wire>triangle_cut;
    std::vector<EdgeRange>edge_range;
    bool complete=false;
};
inline uint64_t star_hash(const std::vector<std::pair<size_t,int>>&star){
    uint64_t h=14695981039346656037ULL;
    auto word=[&](uint64_t v){for(int k=0;k<8;k++){h^=(v>>(8*k))&255U;h*=1099511628211ULL;}};
    word(star.size());for(const auto&v:star){word(v.first);word(uint64_t(int64_t(v.second)));}
    return h;
}
template<class Boundary>inline std::vector<Wire>pack(const Boundary&boundary,
        const std::vector<oc::Cut>&cuts,int rank,int ranks,uint64_t maximum_wires){
    require(rank>=0&&rank<ranks&&ranks<=4&&maximum_wires>0,
        "bounded original source triangle cut fanout rank and budget");
    std::vector<Wire>out;
    for(const auto&cut:cuts){
        require(cut.edge_owner==rank&&cut.source_edge!=UINT64_MAX&&
            cut.lower>=0&&cut.lower<=cut.upper&&cut.upper<=1&&
            !cut.target_segments.empty()&&cut.incidences>=1&&cut.incidences<=2&&
            cut.target_segments.size()<=boundary.segments(),
            "canonical source edge cut on its original edge owner");
        auto member=cut.target_segments.front();
        require(member.first<boundary.segments()&&
            ((cut.at_target_node&&member.second>=0)||
             (!cut.at_target_node&&member.second==-1))&&
            member.second>=-1&&member.second<=1,
            "canonical representative target patch edge of source cut");
        if(cut.at_target_node){
            require(member.second>=0&&
                boundary.trace_key(member.first,member.second,cut.source_edge,
                    cut.source_vertex)==cut.target_node,
                "same canonical target chart node for routed source edge cut");
            auto expected=boundary.trace_star(cut.source_edge,cut.target_node,
                cut.source_vertex);
            std::sort(expected.begin(),expected.end());
            require(expected==cut.target_segments,
                "complete target chart node star before triangle cut routing");
        }else require(cut.target_segments.size()==1&&member.second==-1&&
            cut.source_vertex==UINT64_MAX,
            "one finite target patch segment for a nonvertex source edge cut");
        require(uint64_t(cut.incidences)<=maximum_wires-out.size(),
            "finite local source triangle cut fanout budget");
        for(int k=0;k<cut.incidences;k++){
            const auto&i=cut.incidence[k];
            require(i.triangle!=UINT64_MAX&&i.source_hash&&
                i.rank>=0&&i.rank<ranks&&i.local_edge>=0&&i.local_edge<3&&
                (i.orientation==1||i.orientation==-1),
                "original oriented source triangle incidence for routed cut");
            Wire w;w.source_registry_hash=boundary.source_registry_hash();
            w.source_epoch=boundary.source_epoch();w.target_chart_hash=boundary.chart_hash();
            w.weight_hash=boundary.weight_hash();w.source_edge=cut.source_edge;
            w.source_triangle=i.triangle;w.source_hash=i.source_hash;
            w.source_vertex=cut.source_vertex;
            w.target_segment=member.first;w.target_endpoint=member.second;
            w.star_count=cut.target_segments.size();w.star_hash=star_hash(cut.target_segments);
            w.lower=cut.lower;w.upper=cut.upper;w.source_rank=i.rank;
            w.edge_owner=rank;w.local_edge=i.local_edge;w.orientation=i.orientation;
            out.push_back(w);
        }
    }
    return out;
}
template<class Boundary>inline Result route(const Boundary&boundary,const std::vector<Wire>&local,
        const std::vector<OwnedTriangle>&owned,
        const std::array<uint64_t,4>&expected_shard_hash,const Budget&budget,MPI_Comm comm){
    int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
    int valid=ranks>0&&ranks<=4&&rank>=0&&rank<ranks&&
        budget.maximum_local_wires>0&&budget.maximum_received_wires>0&&
        local.size()<=budget.maximum_local_wires&&boundary.source_registry_hash()&&
        boundary.source_epoch()&&boundary.chart_hash()&&boundary.weight_hash()&&
        boundary.segments()>0;
    for(int p=0;p<ranks&&p<4;p++)valid&=expected_shard_hash[p]!=0;
    for(const auto&w:local)valid&=w.edge_owner==rank&&w.source_rank>=0&&
        w.source_rank<ranks&&w.source_triangle!=UINT64_MAX&&w.source_hash&&
        w.source_registry_hash==boundary.source_registry_hash()&&
        w.source_epoch==boundary.source_epoch()&&
        w.target_chart_hash==boundary.chart_hash()&&w.weight_hash==boundary.weight_hash()&&
        w.target_segment<boundary.segments()&&w.local_edge>=0&&w.local_edge<3&&
        (w.orientation==1||w.orientation==-1)&&
        w.target_endpoint>=-1&&w.target_endpoint<=1&&
        (w.source_vertex==UINT64_MAX||
            (w.target_endpoint>=0&&w.lower==w.upper&&
             (w.lower==0||w.lower==1)))&&
        w.lower>=0&&w.lower<=w.upper&&w.upper<=1&&w.star_count>0;
    int all_valid=0;
    require(MPI_Allreduce(&valid,&all_valid,1,MPI_INT,MPI_MIN,comm)==MPI_SUCCESS,
        "collective original triangle cut route preflight");
    require(all_valid,"all edge owners provide bounded sealed triangle cut wires");
    std::array<uint64_t,9>seal={boundary.source_registry_hash(),boundary.source_epoch(),
        boundary.chart_hash(),boundary.weight_hash(),uint64_t(boundary.segments()),
        expected_shard_hash[0],expected_shard_hash[1],expected_shard_hash[2],
        expected_shard_hash[3]};
    std::vector<uint64_t>all_seals(seal.size()*ranks);
    require(MPI_Allgather(seal.data(),int(seal.size()),MPI_UINT64_T,
        all_seals.data(),int(seal.size()),MPI_UINT64_T,comm)==MPI_SUCCESS,
        "complete source/target triangle cut routing seal exchange");
    for(int p=0;p<ranks;p++)for(size_t j=0;j<seal.size();j++)
        require(all_seals[p*seal.size()+j]==seal[j],
            "all ranks route source edge cuts against the same accepted source and target chart");
    std::vector<std::vector<Wire>>send(ranks);
    for(const auto&w:local)send[w.source_rank].push_back(w);
    std::vector<uint64_t>send_count(ranks),recv_count(ranks);
    for(int p=0;p<ranks;p++)send_count[p]=send[p].size();
    require(MPI_Alltoall(send_count.data(),1,MPI_UINT64_T,
        recv_count.data(),1,MPI_UINT64_T,comm)==MPI_SUCCESS,
        "complete original triangle cut peer counts");
    uint64_t send_total=0,recv_total=0;
    std::vector<int>send_bytes(ranks),recv_bytes(ranks),send_offset(ranks),recv_offset(ranks);
    valid=1;
    for(int p=0;p<ranks;p++){
        if(send_count[p]>uint64_t(INT_MAX/sizeof(Wire))||
            recv_count[p]>uint64_t(INT_MAX/sizeof(Wire))){valid=0;continue;}
        uint64_t sb=send_count[p]*sizeof(Wire),rb=recv_count[p]*sizeof(Wire);
        if(send_total+sb>INT_MAX||recv_total+rb>INT_MAX){valid=0;continue;}
        send_offset[p]=int(send_total);recv_offset[p]=int(recv_total);
        send_bytes[p]=int(sb);recv_bytes[p]=int(rb);
        send_total+=sb;recv_total+=rb;
    }
    valid&=recv_total/sizeof(Wire)<=budget.maximum_received_wires;
    require(MPI_Allreduce(&valid,&all_valid,1,MPI_INT,MPI_MIN,comm)==MPI_SUCCESS,
        "collective original triangle cut routing byte budget");
    require(all_valid,"finite complete original triangle cut MPI bytes and receiver budget");
    std::vector<Wire>packed(send_total/sizeof(Wire)),received(recv_total/sizeof(Wire));
    for(int p=0;p<ranks;p++)if(!send[p].empty())
        std::memcpy(reinterpret_cast<char*>(packed.data())+send_offset[p],
            send[p].data(),send_bytes[p]);
    require(MPI_Alltoallv(packed.data(),send_bytes.data(),send_offset.data(),MPI_BYTE,
        received.data(),recv_bytes.data(),recv_offset.data(),MPI_BYTE,comm)==MPI_SUCCESS,
        "complete original source triangle cut delivery");
    Result out;valid=1;
    try{
        for(int p=0;p<ranks;p++)for(uint64_t j=0;j<recv_count[p];j++){
            const auto&w=received[size_t(recv_offset[p]/sizeof(Wire))+j];
            require(w.edge_owner==p,
                "routed triangle cut names its actual sending original edge owner");
        }
        auto topology=owned;
        std::sort(topology.begin(),topology.end(),[](const auto&a,const auto&b){return a.triangle<b.triangle;});
        for(size_t j=0;j<topology.size();j++)require(topology[j].triangle!=UINT64_MAX&&
            (j==0||topology[j-1].triangle<topology[j].triangle),
            "unique original source triangles on this receiving rank");
        for(const auto&t:topology)require(t.vertex[0]!=UINT64_MAX&&
            t.vertex[1]!=UINT64_MAX&&t.vertex[2]!=UINT64_MAX&&
            t.vertex[0]!=t.vertex[1]&&t.vertex[1]!=t.vertex[2]&&
            t.vertex[2]!=t.vertex[0],
            "three distinct canonical original vertices per owned source triangle");
        for(const auto&w:received){
            auto found=std::lower_bound(topology.begin(),topology.end(),w.source_triangle,
                [](const OwnedTriangle&t,uint64_t id){return t.triangle<id;});
            require(w.source_rank==rank&&w.source_hash==expected_shard_hash[rank]&&
                w.edge_owner>=0&&w.edge_owner<ranks&&
                w.source_registry_hash==boundary.source_registry_hash()&&
                w.source_epoch==boundary.source_epoch()&&
                w.target_chart_hash==boundary.chart_hash()&&w.weight_hash==boundary.weight_hash()&&
                found!=topology.end()&&found->triangle==w.source_triangle&&
                w.local_edge>=0&&w.local_edge<3&&
                found->edge[w.local_edge]==w.source_edge&&
                found->orientation[w.local_edge]==w.orientation&&
                w.target_segment<boundary.segments()&&
                w.target_endpoint>=-1&&w.target_endpoint<=1&&
                w.lower>=0&&w.lower<=w.upper&&w.upper<=1,
                "routed cut belongs to an actually owned original source triangle boundary edge");
            if(w.target_endpoint>=0){
                auto node=boundary.trace_key(w.target_segment,w.target_endpoint,
                    w.source_edge,w.source_vertex);
                auto star=boundary.trace_star(w.source_edge,node,w.source_vertex);
                std::sort(star.begin(),star.end());
                require(w.star_count==star.size()&&w.star_hash==star_hash(star),
                    "same complete target chart node star reconstructed on triangle owner");
            }else require(w.source_vertex==UINT64_MAX&&w.star_count==1&&
                w.star_hash==star_hash({{size_t(w.target_segment),-1}}),
                "one finite target patch segment at nonvertex triangle cut");
        }
        std::sort(received.begin(),received.end(),[](const Wire&a,const Wire&b){
            if(a.source_triangle!=b.source_triangle)return a.source_triangle<b.source_triangle;
            if(a.local_edge!=b.local_edge)return a.local_edge<b.local_edge;
            if(a.orientation!=b.orientation)return a.orientation<b.orientation;
            if(a.orientation==1)return std::tie(a.lower,a.upper,a.target_segment)<
                std::tie(b.lower,b.upper,b.target_segment);
            return std::tie(a.upper,a.lower,a.target_segment)>
                std::tie(b.upper,b.lower,b.target_segment);
        });
        for(size_t j=1;j<received.size();j++){
            const auto&a=received[j-1],&b=received[j];
            if(a.source_triangle!=b.source_triangle||a.local_edge!=b.local_edge)continue;
            require(a.source_edge==b.source_edge&&a.orientation==b.orientation&&
                (a.orientation==1?a.upper<b.lower:b.upper<a.lower),
                "strictly ordered disjoint cuts along each original source triangle boundary edge");
        }
        out.triangle_cut=std::move(received);out.complete=true;
        size_t cursor=0;
        for(const auto&t:topology)for(int k=0;k<3;k++){
            Result::EdgeRange range;range.triangle=t.triangle;
            range.source_edge=t.edge[k];range.local_edge=k;
            range.orientation=t.orientation[k];
            range.start_vertex=t.vertex[k];range.end_vertex=t.vertex[(k+1)%3];
            range.begin=cursor;
            while(cursor<out.triangle_cut.size()&&
                out.triangle_cut[cursor].source_triangle==t.triangle&&
                out.triangle_cut[cursor].local_edge==k){
                const auto&w=out.triangle_cut[cursor];
                if(w.lower==0||w.upper==1){
                    require(w.lower==w.upper&&w.target_endpoint>=0,
                        "original source edge endpoint cut is exact and names a target chart node");
                    auto node=boundary.trace_key(w.target_segment,w.target_endpoint,
                        w.source_edge,w.source_vertex);
                    uint64_t vertex=w.lower==0?
                        (w.orientation==1?range.start_vertex:range.end_vertex):
                        (w.orientation==1?range.end_vertex:range.start_vertex);
                    require(node.size()==7&&node[0]==0&&node[1]==1&&
                        node[2]==vertex&&node[3]==0&&node[4]==0&&
                        node[5]==0&&node[6]==0&&w.source_vertex==vertex,
                        "target chart endpoint cut is the same exact original source triangle vertex");
                }else require(w.lower>0&&w.upper<1,
                    "nonendpoint cut lies strictly inside its original source triangle edge");
                cursor++;
            }
            range.end=cursor;out.edge_range.push_back(range);
        }
        require(cursor==out.triangle_cut.size()&&
            out.edge_range.size()==3*topology.size(),
            "all routed cuts belong to exactly one of three ordered original triangle boundary edges");
    }catch(const std::exception&e){
        fprintf(stderr,"rank%d TRIANGLE_CUT_RECEIVE_FAIL %s\n",rank,e.what());valid=0;
    }
    require(MPI_Allreduce(&valid,&all_valid,1,MPI_INT,MPI_MIN,comm)==MPI_SUCCESS,
        "collective original triangle cut receiving audit");
    require(all_valid,"every routed cut reconstructs the same chart node and original triangle edge");
    return out;
}
}
#endif
