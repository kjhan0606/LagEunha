#include "../material_compact_rezone_parallel.cpp"
#include <quadmath.h>
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace ce=compact_overlap_events;

static std::vector<std::string> fields(const std::string&line){
    std::vector<std::string>out;
    size_t first=0;
    for(size_t next;(next=line.find(',',first))!=std::string::npos;first=next+1)
        out.push_back(line.substr(first,next-first));
    out.push_back(line.substr(first));return out;
}
static __float128 parse_q(const std::string&word){
    char*end=nullptr;__float128 value=strtoflt128(word.c_str(),&end);
    ce::require(end==word.c_str()+word.size()&&finiteq(value),
        "complete finite binary128 event endpoint");
    char roundtrip[128];
    int n=quadmath_snprintf(roundtrip,sizeof(roundtrip),"%.36Qg",value);
    ce::require(n>0&&size_t(n)<sizeof(roundtrip)&&word==roundtrip,
        "lossless binary128 source event decimal");
    return value;
}
static std::string q_text(__float128 value){
    char text[128];int n=quadmath_snprintf(text,sizeof(text),"%.36Qg",value);
    ce::require(n>0&&size_t(n)<sizeof(text),"write binary128 canonical endpoint");
    return text;
}
int main(int argc,char**argv){
    try{
        ce::require(argc==5,"usage: reconcile CHUNK EXTERNAL_FLAGS OUTPUT EXPECTED_SEGMENTS");
        const std::filesystem::path chunk=argv[1],flags=argv[2],output=argv[3];
        uint64_t expected_segments=std::stoull(argv[4]);
        ce::require(expected_segments>0,"finite complete target segment registry");
        std::ifstream exterior_file(flags,std::ios::binary);
        ce::require(bool(exterior_file),"read sealed accepted source exterior flags");
        std::vector<uint8_t>external_edge((std::istreambuf_iterator<char>(exterior_file)),{});
        ce::require(!external_edge.empty()&&external_edge.size()<=4000000,
            "finite complete accepted global source edge registry");
        std::array<uint64_t,4>shard{};
        uint64_t registry=0,epoch=0,chart=0,weight=0,law=0;
        size_t first=SIZE_MAX,last=0,total_rows=0;
        std::array<std::vector<ce::Wire>,4>owned;
        for(int rank=0;rank<4;rank++){
            std::ifstream manifest(chunk/("source_edge_manifest_rank"+
                std::to_string(rank)+".csv"));
            ce::require(bool(manifest),"read sealed source-edge rank manifest");
            std::string header,line;
            ce::require(bool(std::getline(manifest,header))&&
                bool(std::getline(manifest,line))&&
                !std::getline(manifest,header),
                "one complete sealed rank manifest row");
            auto m=fields(line);ce::require(m.size()==18,
                "complete source-edge rank manifest schema");
            ce::require(std::stoi(m[0])==rank&&std::stoull(m[3])==external_edge.size(),
                "same rank and complete original edge registry");
            if(rank==0){first=std::stoull(m[1]);last=std::stoull(m[2]);
                registry=std::stoull(m[4]);epoch=std::stoull(m[6]);
                law=std::stoull(m[9]);chart=std::stoull(m[10]);weight=std::stoull(m[11]);}
            ce::require(first==std::stoull(m[1])&&last==std::stoull(m[2])&&
                registry==std::stoull(m[4])&&epoch==std::stoull(m[6])&&
                law==std::stoull(m[9])&&chart==std::stoull(m[10])&&
                weight==std::stoull(m[11])&&first<last&&last<=external_edge.size(),
                "all ranks have one source, target, mass law and edge range");
            shard[rank]=std::stoull(m[5]);
            std::ifstream csv(chunk/("source_edge_wires_rank"+
                std::to_string(rank)+".csv"));
            ce::require(bool(csv)&&bool(std::getline(csv,header))&&
                header.find("event_lower,event_upper")!=std::string::npos,
                "complete source-edge wire CSV header");
            size_t rows=0;
            while(std::getline(csv,line)){
                auto f=fields(line);ce::require(f.size()==22,
                    "complete source-edge incidence CSV row");
                ce::Wire w;
                w.source_registry_hash=std::stoull(f[0]);
                w.source_hash=std::stoull(f[1]);w.source_epoch=std::stoull(f[2]);
                w.target_chart_hash=std::stoull(f[3]);w.weight_hash=std::stoull(f[4]);
                ce::require(std::stoull(f[5])==first&&std::stoull(f[6])==last,
                    "one explicit checkpointed original edge range");
                w.event.cell[0]=std::stoull(f[7]);w.event.cell[1]=std::stoull(f[8]);
                w.event.edge=std::stoull(f[9]);w.event.ordinal=std::stoull(f[10]);
                w.event.lower=parse_q(f[11]);w.event.upper=parse_q(f[12]);
                w.event.derivative_sign=std::stoi(f[13]);
                w.target_segment=std::stoull(f[14]);w.target_endpoint=std::stoi(f[15]);
                w.source_vertex=std::stoull(f[16]);w.source_rank=std::stoi(f[17]);
                w.incident_triangle=std::stoull(f[18]);w.local_edge=std::stoi(f[19]);
                w.orientation=std::stoi(f[20]);w.edge_owner=std::stoi(f[21]);
                ce::require(w.source_rank==rank&&w.source_hash==shard[rank]&&
                    w.event.edge>=first&&w.event.edge<last&&
                    w.edge_owner>=0&&w.edge_owner<4,
                    "source rank and original edge identity before owner routing");
                owned[w.edge_owner].push_back(w);rows++;
            }
            ce::require(csv.eof()&&rows==std::stoull(m[16]),
                "all sealed source rank event rows read once");
            total_rows+=rows;
        }
        size_t total_events=0;
        for(int owner=0;owner<4;owner++){
            ce::Budget budget;budget.maximum_received_wires=1000000;
            budget.maximum_canonical_events=1000000;
            auto canonical=ce::reconcile(std::move(owned[owner]),external_edge,
                owner,4,shard,registry,epoch,chart,weight,expected_segments,budget);
            std::ofstream csv(output/("canonical_owner"+std::to_string(owner)+".csv"));
            ce::require(bool(csv),"create certified canonical event owner ledger");
            csv<<"event_edge,event_cell_low,event_cell_high,event_ordinal,event_lower,event_upper,event_derivative_sign,target_segment,target_endpoint,source_vertex,edge_owner,incidences,source_rank0,source_rank1,incident_triangle0,incident_triangle1,local_edge0,local_edge1,orientation0,orientation1\n";
            for(const auto&c:canonical)
                csv<<c.event.edge<<","<<c.event.cell[0]<<","<<c.event.cell[1]<<","
                    <<c.event.ordinal<<","<<q_text(c.event.lower)<<","
                    <<q_text(c.event.upper)<<","<<c.event.derivative_sign<<","
                    <<c.target_segment<<","<<c.target_endpoint<<","<<c.source_vertex<<","
                    <<c.edge_owner<<","<<c.incidences<<","<<c.source_rank[0]<<","
                    <<c.source_rank[1]<<","<<c.incident_triangle[0]<<","
                    <<c.incident_triangle[1]<<","<<c.local_edge[0]<<","
                    <<c.local_edge[1]<<","<<c.orientation[0]<<","
                    <<c.orientation[1]<<"\n";
            csv.close();ce::require(bool(csv),"flush canonical event owner ledger");
            total_events+=canonical.size();
            std::cout<<"NATIVE_CANONICAL_RECONCILE_PASS owner="<<owner
                <<" events="<<canonical.size()<<std::endl;
        }
        std::cout<<"NATIVE_SOURCE_EDGE_CHUNK_PASS first="<<first<<" last="<<last
            <<" incidence_wires="<<total_rows<<" canonical_events="<<total_events
            <<" mass_law_solved=0"<<std::endl;
        return 0;
    }catch(const std::exception&e){
        std::cerr<<"NATIVE_SOURCE_EDGE_CHUNK_FAIL "<<e.what()<<std::endl;return 2;
    }
}
