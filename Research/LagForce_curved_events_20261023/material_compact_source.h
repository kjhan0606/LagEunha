#ifndef MATERIAL_COMPACT_SOURCE_H
#define MATERIAL_COMPACT_SOURCE_H
// Included after the original native Owned and distributed Ownership headers.
// Imports a real accepted full field. Does not create ICs, replay an RHS,
// reconstruct thermal state from mean U, or perform a state transfer.
#include "material_compact_checkpoint.h"
#include "material_full_initial_batch.h"
#include <fstream>
#include <charconv>
#include <cctype>
#include <memory>
#include <cstdlib>

namespace compact_source {
using R=compact_p3::R;
namespace md=material_distributed;
inline void require(bool b,const char*s){compact_p3::require(b,s);}
inline std::string small_text(const std::string&path){
    std::ifstream in(path,std::ios::binary|std::ios::ate);require(bool(in),"verified source metadata readable");
    auto size=in.tellg();require(size>0&&size<1048576,"bounded nonempty immutable source metadata");
    std::string out(size_t(size),'\0');in.seekg(0);in.read(out.data(),out.size());require(bool(in),"complete immutable source metadata");return out;
}
inline size_t json_field(const std::string&s,const char*name){
    std::string key="\""+std::string(name)+"\"";size_t at=s.find(key);
    require(at!=std::string::npos&&s.find(key,at+key.size())==std::string::npos,"one explicit saved metadata key");
    at+=key.size();while(at<s.size()&&std::isspace(static_cast<unsigned char>(s[at])))at++;
    require(at<s.size()&&s[at++]==':',"saved metadata field separator");
    while(at<s.size()&&std::isspace(static_cast<unsigned char>(s[at])))at++;return at;
}
inline uint64_t json_u64(const std::string&s,const char*key){
    size_t at=json_field(s,key);uint64_t value=0;auto result=std::from_chars(s.data()+at,s.data()+s.size(),value);
    require(result.ec==std::errc()&&result.ptr>s.data()+at,"exact uint64 saved metadata, not rounded double");
    require(result.ptr==s.data()+s.size()||*result.ptr==','||*result.ptr=='}'||std::isspace(static_cast<unsigned char>(*result.ptr)),"complete uint64 metadata token");return value;
}
inline bool json_bool(const std::string&s,const char*key){
    size_t at=json_field(s,key);bool yes=s.compare(at,4,"true")==0,no=s.compare(at,5,"false")==0;
    require(yes||no,"explicit saved boolean");size_t end=at+(yes?4:5);
    require(end==s.size()||s[end]==','||s[end]=='}'||std::isspace(static_cast<unsigned char>(s[end])),"complete boolean metadata token");return yes;
}
inline double json_double(const std::string&s,const char*key){
    size_t at=json_field(s,key);char*end=nullptr;double value=strtod(s.c_str()+at,&end);
    require(end>s.c_str()+at&&std::isfinite(value),"finite saved physical parameter");
    require(*end==0||*end==','||*end=='}'||std::isspace(static_cast<unsigned char>(*end)),"complete physical parameter token");return value;
}
struct Inputs {
    // Every payload below must be bound by the transaction's prelaunch SHA256
    // manifest. No former step4 k1 or finite-CONTROL packet is an input.
    std::string initial,accepted_prefix,immutable_manifest,positive_directory,whole_audit,reference_receipt;
    uint64_t accepted_steps=0;double accepted_time=0;
};
struct Model {
    material_full::Owned owned;
    std::vector<STCoefficient>generating;
    std::vector<std::array<R,2>>generating_delta;
    compact_p3::ImportedCheckpoint accepted;
    std::unique_ptr<md::Ownership>layout;
    std::string fixed_mass_path,positive_binding;
    double gamma=0,P0=0,eps=0;
    std::vector<MFEvolvedCell>current_context;
    // Missing ENDPOINT context is derived from actual current geometry/rule,
    // never copied from stage4 or fitted to accepted mean internal energy.
    bool current_context_bound=false;
    bool original_reference_receipt_bound=false;
};
inline uint64_t audited_rule_generation(const Inputs&input,uint64_t triangles){
    auto contract=small_text(input.whole_audit+"/contract.json");
    require(json_bool(contract,"whole_rule_checked")&&!json_bool(contract,"prefix_only")&&json_bool(contract,"audit_pass")&&json_u64(contract,"triangles_checked")==triangles,"independent complete stored-rule audit, not prefix or proposal");
    require(json_u64(contract,"original_positive_points_accounted")==2040256512ULL,"same original full reference rule");
    uint64_t streams[4]={},count=0;
    for(uint64_t rank=0;rank<4;rank++){
        auto record=small_text(input.whole_audit+"/run/records_r"+std::to_string(rank)+".json");
        require(json_u64(record,"rank")==rank&&!json_bool(record,"prefix_only")&&json_bool(record,"record_and_contraction_checks_pass"),"actual whole-rule shard identity");
        count+=json_u64(record,"records_checked");streams[rank]=json_u64(record,"stream_FNV");
    }
    require(count==triangles,"all independent rule streams accounted for");
    return compact_p3::bytes_hash(streams,sizeof(streams));
}
inline Model import_source(const Inputs&input,MPI_Comm comm){
    int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
    require(ranks==4&&input.accepted_steps>=5&&input.accepted_time>0,"same actual four-owner accepted source epoch");
    Model out;std::string suffix=std::to_string(rank);
    out.owned=material_full::load_basis(input.initial+"/full_conditioned_basis_r"+suffix+".bin");
    auto&o=out.owned;
    require(o.global_triangles==1138432&&o.global_dofs==10816441&&o.f.triangle.size()==o.source_triangle.size()&&o.f.coefficient.size()==o.global.size(),"full original source, not a small control state");
    require(std::is_sorted(o.f.cell.begin(),o.f.cell.end())&&std::adjacent_find(o.f.cell.begin(),o.f.cell.end())==o.f.cell.end(),"unique original material source cells");
    uint64_t counts[2]={o.f.cell.size(),o.f.triangle.size()};
    MPI_Allreduce(MPI_IN_PLACE,counts,2,MPI_UINT64_T,MPI_SUM,comm);
    require(counts[0]==71080&&counts[1]==o.global_triangles,"all original physical cells and source triangles imported");
    uint64_t generation=audited_rule_generation(input,o.global_triangles);
    out.positive_binding="POSITIVE_ORIGINAL_POINT_SUBSET_V1\n"+input.positive_directory+"\n"+std::to_string(generation);
    compact_p3::AnchorExpectation expected;
    expected.rank=rank;expected.ranks=ranks;expected.local_coefficients=o.global.size();
    expected.steps=input.accepted_steps;expected.time=input.accepted_time;expected.generation=generation;
    expected.exact_binding="FRESH_CONSTANT_PRESSURE_GENERATING_CONTEXT\n"+input.initial+"\n"+small_text(input.immutable_manifest)+"\n"+out.positive_binding;
    out.accepted=compact_p3::read_checkpoint(input.accepted_prefix+suffix+".bin",expected);
    out.generating=std::move(o.f.coefficient);o.f.coefficient=std::move(out.accepted.coefficient);
    out.layout=std::make_unique<md::Ownership>(o.global,o.owner,comm);
    require(out.layout->global_owned==o.global_dofs,"same complete native material DOF ownership");
    // Verify ghost bits rather than silently replacing a mismatched checkpoint.
    md::Vec values(o.global.size());uint64_t mismatch=0;
    for(int d=0;d<4;d++){
        for(size_t j=0;j<values.size();j++)values[j]=d<2?o.f.coefficient[j].x[d]:o.f.coefficient[j].v[d-2];
        out.layout->pull(values);
        for(size_t j=0;j<values.size();j++){
            double actual=double(values[j]),stored=d<2?o.f.coefficient[j].x[d]:o.f.coefficient[j].v[d-2];
            mismatch+=memcmp(&actual,&stored,sizeof(double))!=0;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE,&mismatch,1,MPI_UINT64_T,MPI_SUM,comm);
    require(!mismatch,"all actual accepted owned/ghost position and velocity bits agree");
    FILE*f=fopen((input.initial+"/full_generating_thermal_r"+suffix+".bin").c_str(),"rb");require(f,"immutable original generating thermal state");
    uint64_t header[6];material_full::read(f,header,6);
    require(header[0]==0x4d46544845524d31ULL&&header[1]==1&&header[2]==o.f.cell.size()&&header[3]==sizeof(R)&&header[4]==uint64_t(rank)&&header[5]==1,"same original generating thermal source ABI");
    double law[4];material_full::read(f,law,4);require(std::isfinite(law[0])&&law[0]>1&&std::isfinite(law[1])&&law[1]>0&&std::isfinite(law[2])&&law[2]>0&&law[3]==0,"actual immutable constant-pressure generating law");
    out.gamma=law[0];out.P0=law[1];out.eps=law[2];out.generating_delta.resize(header[2]);
    material_full::read(f,out.generating_delta.data(),header[2]);require(fgetc(f)==EOF&&fclose(f)==0,"complete immutable thermal payload");
    for(const auto&delta:out.generating_delta)for(R value:delta)require(std::isfinite(value),"finite immutable source thermal correction");
    double low[3]={law[0],law[1],law[2]},high[3]={law[0],law[1],law[2]};
    MPI_Allreduce(MPI_IN_PLACE,low,3,MPI_DOUBLE,MPI_MIN,comm);MPI_Allreduce(MPI_IN_PLACE,high,3,MPI_DOUBLE,MPI_MAX,comm);
    for(int d=0;d<3;d++)require(low[d]==high[d],"all owners share the original thermal/gravity law");
    for(size_t s=0;s<o.f.triangle.size();s++){
        const auto&t=o.f.triangle[s];require(t.cell>=0&&size_t(t.cell)<o.f.cell.size()&&o.source_cell[s]==o.f.cell[t.cell]&&t.eps==out.eps,"same source-local topology and softened gravity");
        for(int a=0;a<t.n;a++)require(t.id[a]>=0&&size_t(t.id[a])<o.global.size(),"full accepted basis coefficient map");
    }
    auto receipt=small_text(input.reference_receipt);
    require(json_bool(receipt,"reference_preparation_complete")&&json_bool(receipt,"exact_support_topology_no_proximity_merge")&&json_u64(receipt,"cells")==counts[0]&&json_u64(receipt,"reference_vertices")==o.global_vertices&&json_u64(receipt,"reference_triangles")==o.global_triangles&&json_u64(receipt,"reference_edges")==o.global_edges,"same previously constructed full native reference, bound by input SHA");
    require(o.global_vertices+o.global_triangles==o.global_edges+1&&json_double(receipt,"total_area")==576&&json_double(receipt,"gamma")==out.gamma&&json_double(receipt,"P0")==out.P0&&json_double(receipt,"eps")==out.eps,"original square-reference topology and actual native physical parameters");
    out.original_reference_receipt_bound=true;
    out.fixed_mass_path=input.initial+"/full_fixed_mass_r"+suffix+".bin";
    // Fixed M remains its original immutable payload and is loaded only when
    // actual overlap/projection needs it, after geometry feasibility checks.
    return out;
}
inline void require_context(const Model&m){
    require(m.current_context_bound&&m.current_context.size()==m.owned.f.cell.size(),"actual accepted endpoint thermal context required before physical transfer");
}
Model import_accepted_state(const Inputs&,MPI_Comm);
} // namespace compact_source
#endif
