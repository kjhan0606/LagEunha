#ifndef MATERIAL_COMPACT_SOURCE_CONTEXT_H
#define MATERIAL_COMPACT_SOURCE_CONTEXT_H
// Missing accepted-endpoint context, evaluated within the ONE full transaction.
// Only original CUDA passes0/1 are used: no force pass, mass solve or evolution.
#include <cuda_runtime.h>
#include <filesystem>
#include "material_full_reduced_rule.h"

namespace compact_source {
struct SourceRule {
    std::vector<uint64_t>offset;
    std::vector<uint32_t>index;
    std::vector<double>weight;
    uint64_t global_points=0;
};
inline SourceRule read_context_rule(const Inputs&input,const Model&m){
    const auto&o=m.owned;const auto&layout=*m.layout;
    auto contract=small_text(input.positive_directory+"/contract.json");
    require(json_bool(contract,"all_original_triangles_processed")&&json_bool(contract,"candidate_ready_for_actual_RK_oracle"),"complete immutable positive functional");
    uint64_t expected=json_u64(contract,"retained_points_with_original_rule_fallbacks");
    require(expected>0&&expected*4<=2040256512ULL,"same complete compressed functional, including fallbacks");
    auto progress=small_text(input.positive_directory+"/run/progress_r"+std::to_string(layout.rank)+".json");
    require(json_u64(progress,"completed_triangles")==o.f.triangle.size(),"completed rule shard, no live prefix");
    uint64_t local_limit=json_u64(progress,"retained_points");require(local_limit<=expected,"bounded same-rank rule storage");
    std::unique_ptr<FILE,int(*)(FILE*)>f(fopen((input.positive_directory+"/run/positive_rule_r"+std::to_string(layout.rank)+".bin").c_str(),"rb"),fclose);
    require(bool(f),"verified complete positive-rule input");uint64_t h[12];material_full::read(f.get(),h,12);
    auto immutable=small_text(input.immutable_manifest);
    require(h[0]==0x4d46504f53525531ULL&&h[1]==1&&h[2]==uint64_t(layout.rank)&&h[3]==uint64_t(layout.ranks)&&h[4]==o.f.triangle.size()&&h[5]==o.global_triangles&&h[6]==sizeof(R)&&h[8]==compact_p3::bytes_hash(immutable.data(),immutable.size())&&h[9]==MF_EVOL_METRICS&&h[10]==16&&h[11]==0,"original source/rule/owner/immutable ABI");
    double anchor;material_full::read(f.get(),&anchor);require(anchor>0&&anchor<m.accepted.physical.time,"same earlier fixed-rule anchor, not a new initial condition");
    SourceRule out;out.offset.resize(o.f.triangle.size()+1);out.index.reserve(local_limit);out.weight.reserve(local_limit);
    uint64_t stream=compact_p3::bytes_hash(h,sizeof(h));stream=compact_p3::bytes_hash(&anchor,sizeof(anchor),stream);
    uint64_t original=0,selected=0;
    for(size_t s=0;s<o.f.triangle.size();s++){
        uint64_t record[6];material_full::read(f.get(),record,6);uint64_t n=o.f.triangle[s].n,rows=(n+2)*(n+2),nodes=record[2];
        require(record[0]==o.source_triangle[s]&&record[1]>0&&record[1]<=1048576&&record[4]==rows&&record[5]==0&&record[3]<=1,"same ordered full-source rule record");
        require((record[3]&&nodes>0&&4*nodes<=record[1])||(!record[3]&&nodes==0),"positive subset or explicit unmodified original rule");
        // Every diagnostic/residual/roundoff byte was checked by the completed
        // WHOLE audit. The transaction prelaunch SHA256 binds those same files.
        // Skip those existing records; do not repeat an audit/NNLS proposal.
        require(fseeko(f.get(),off_t((16+2*rows)*sizeof(R)),SEEK_CUR)==0,"skip already audited immutable all-row diagnostics");
        size_t at=out.index.size();require(nodes<=local_limit-at,"bounded complete selected-point prefix");
        out.index.resize(at+nodes);out.weight.resize(at+nodes);
        material_full::read(f.get(),out.index.data()+at,nodes);material_full::read(f.get(),out.weight.data()+at,nodes);
        uint64_t checksum;material_full::read(f.get(),&checksum);stream=compact_p3::bytes_hash(&checksum,sizeof(checksum),stream);
        for(size_t j=0;j<nodes;j++)require(out.index[at+j]<record[1]&&(!j||out.index[at+j-1]<out.index[at+j])&&out.weight[at+j]>0&&std::isfinite(out.weight[at+j]),"same positive unique original-point subset");
        out.offset[s+1]=out.index.size();original+=record[1];selected+=record[3]?nodes:record[1];
    }
    require(fgetc(f.get())==EOF&&fclose(f.release())==0&&selected==local_limit,"complete immutable positive-rule shard");
    uint64_t streams[4],counts[2]={original,selected};
    MPI_Allgather(&stream,1,MPI_UINT64_T,streams,1,MPI_UINT64_T,layout.comm);
    MPI_Allreduce(MPI_IN_PLACE,counts,2,MPI_UINT64_T,MPI_SUM,layout.comm);
    require(compact_p3::bytes_hash(streams,sizeof(streams))==m.accepted.numerical.generation&&counts[0]==2040256512ULL&&counts[1]==expected,"same independently audited accepted numerical functional");
    out.global_points=expected;return out;
}
inline std::vector<std::string> csv_columns(const std::string&s){
    std::vector<std::string>out;size_t begin=0;for(size_t j=0;j<=s.size();j++)if(j==s.size()||s[j]==','){out.push_back(s.substr(begin,j-begin));begin=j+1;}return out;
}
inline R csv_real(const std::string&s){char*end=nullptr;R value=strtold(s.c_str(),&end);require(end==s.c_str()+s.size()&&end!=s.c_str()&&std::isfinite(value),"finite complete saved endpoint scalar");return value;}
struct EndpointBaseline {uint64_t cell=0;R volume=0,vdot=0,mass=0,U=0,Phi=0,K=0,L=0;double mean_deltaP=0;};
inline std::vector<EndpointBaseline> saved_endpoint(const Inputs&input,const Model&m){
    auto directory=std::filesystem::path(input.accepted_prefix).parent_path()/"endpoint";
    std::string suffix=std::to_string(m.layout->rank),line;std::vector<EndpointBaseline>out(m.owned.f.cell.size());
    auto load=[&](const std::string&file,auto consume){
        std::ifstream in(directory/file);require(bool(in),"actual accepted endpoint baseline, not an RK stage");
        require(bool(std::getline(in,line)),"saved endpoint CSV schema");auto header=csv_columns(line);std::map<std::string,size_t>column;
        for(size_t j=0;j<header.size();j++)require(column.emplace(header[j],j).second,"unique saved endpoint column");
        size_t row=0;while(std::getline(in,line)){
            auto fields=csv_columns(line);require(row<out.size()&&fields.size()==header.size(),"complete ordered accepted endpoint cells");
            auto get=[&](const char*key){auto at=column.find(key);require(at!=column.end(),"required actual endpoint field");return csv_real(fields[at->second]);};
            R id=get("cell");require(id==R(m.owned.f.cell[row]),"same accepted endpoint physical cell identity");out[row].cell=m.owned.f.cell[row];consume(out[row],get);row++;
        }require(in.eof()&&row==out.size(),"all accepted endpoint cell rows present");
    };
    load("quadrature_cells_r"+suffix+".csv",[&](EndpointBaseline&b,auto get){
        require(get("acceleration_compared")==0,"endpoint baseline, no stage acceleration response substituted");
        b.volume=get("volume");b.vdot=get("volume_dot");b.mass=get("mean_rho_low")*b.volume;
        b.U=get("U_low");b.Phi=get("Phi_low");b.mean_deltaP=double(get("mean_deltaP_low"));
        require(b.volume>0&&b.mass>0&&b.U>0&&m.P0+b.mean_deltaP>0,"positive actual accepted source baseline");
    });
    load("physical_cells_r"+suffix+".csv",[&](EndpointBaseline&b,auto get){
        require(get("volume")==b.volume&&get("U")==b.U,"same actual physical and quadrature endpoint state");
        R mass=get("mass");require(fabsl(mass-b.mass)<128*LDBL_EPSILON*mass,"same saved endpoint mass accounting");
        b.mass=mass;b.K=get("K_point");b.L=get("L_point");require(b.K>0,"same positive accepted cell kinetic scale");
    });return out;
}
struct SourceDevice {
    void*ptr=nullptr;
    SourceDevice()=default;SourceDevice(const SourceDevice&)=delete;SourceDevice&operator=(const SourceDevice&)=delete;
    ~SourceDevice(){mf_initial_destroy(ptr);}
};
struct CurrentContext {
    std::unique_ptr<SourceDevice>device;
    std::vector<EndpointBaseline>baseline;
    std::vector<std::array<R,MF_EVOL_METRICS>>moment;
    R mass_difference=0,U_difference=0,Phi_difference=0,K_difference=0,L_difference=0,pressure_difference=0;
    R conditional_corrected_J_lower=0;
    double kernel_seconds[2]={},wall_seconds=0;
    uint64_t coefficient_hash=0;
};
inline CurrentContext reconstruct_endpoint_context(Model&m,const Inputs&input,const GeometryBound&geometry){
    require(m.layout&&geometry.canonical_global_injective&&geometry.exact_shared_trace&&!m.current_context_bound,"one missing-context derivation after the actual source geometry prerequisite");
    const auto&o=m.owned;const auto&layout=*m.layout;double begin=MPI_Wtime();CurrentContext out;
    out.coefficient_hash=compact_p3::bytes_hash(o.f.coefficient.data(),o.f.coefficient.size()*sizeof(STCoefficient));
    require(out.coefficient_hash==geometry.source_coefficient_hash&&geometry.source_time==m.accepted.physical.time&&geometry.source_steps==m.accepted.physical.steps,"same accepted source geometry and physical epoch");
    out.baseline=saved_endpoint(input,m);auto rule=read_context_rule(input,m);uint64_t expected_points=rule.global_points;
    MPI_Comm local;MPI_Comm_split_type(layout.comm,MPI_COMM_TYPE_SHARED,layout.rank,MPI_INFO_NULL,&local);int local_rank=0;MPI_Comm_rank(local,&local_rank);MPI_Comm_free(&local);
    int devices=0;require(cudaGetDeviceCount(&devices)==cudaSuccess&&devices>0&&cudaSetDevice(local_rank%devices)==cudaSuccess,"same native MPI local GPU mapping");
    GfsEnrichConfig cfg={};cfg.na=64;cfg.nr=28;cfg.field.ic.kind=2;cfg.field.ic.gamma=m.gamma;cfg.field.ic.pressure=m.P0;
    cfg.field.ic.cx=cfg.field.ic.cy=12;cfg.field.ic.rin=2;cfg.field.ic.rout=8;cfg.field.ic.eps=m.eps;
    std::vector<double>x,w;shared_trace::gauss(64,x,w);for(int j=0;j<64;j++){cfg.an[j]=x[j];cfg.aw[j]=w[j];}
    shared_trace::gauss(28,x,w);for(int j=0;j<28;j++){cfg.rn[j]=x[j];cfg.rw[j]=w[j];}
    out.device=std::make_unique<SourceDevice>();
    require(!mf_initial_create(o.f.triangle.data(),o.f.triangle.size(),m.generating.data(),m.generating.size(),o.rule.data(),&cfg,&out.device->ptr),"original native source device allocation; no initial-force kernel");
    require(!mf_evolved_bind(out.device->ptr,m.generating.data(),m.generating.size(),o.f.cell.size()),"immutable source generating context bound");
    require(!mf_evolved_rule_upload(out.device->ptr,rule.offset.data(),o.f.triangle.size(),rule.index.data(),rule.weight.data(),rule.index.size())&&!mf_evolved_rule_use(out.device->ptr,1),"same accepted positive functional uploaded");
    rule=SourceRule{}; // Device owns the copy; no all-point host atlas retained.
    std::vector<MFEvolvedCell>candidate(o.f.cell.size());for(size_t c=0;c<candidate.size();c++)candidate[c].generating_delta=double(m.generating_delta[c][0]);
    out.moment.resize(o.f.cell.size());std::vector<MFEvolvedBlock>batch(MF_BATCH_TRIANGLES);
    for(int pass=0;pass<2;pass++){
        require(!mf_evolved_stage(out.device->ptr,o.f.coefficient.data(),o.f.coefficient.size(),candidate.data(),candidate.size()),"actual accepted source coefficients/context uploaded");
        for(auto&cell:out.moment)cell={};uint64_t points=0;
        for(size_t at=0;at<o.f.triangle.size();at+=batch.size()){
            size_t count=std::min(batch.size(),o.f.triangle.size()-at);double seconds=0;
            require(!mf_evolved_batch(out.device->ptr,pass,at,count,batch.data(),&seconds),"missing-context original CUDA geometry/thermal pass");out.kernel_seconds[pass]+=seconds;
            for(size_t j=0;j<count;j++){
                const auto&b=batch[j];require(b.bad_point<0&&b.minJ>0,"actual source point validity in context derivation");points+=b.points;
                auto&cell=out.moment.at(o.f.triangle[at+j].cell);
                for(int k=0;k<MF_EVOL_METRICS;k++){R value=R(b.metric[k].hi)+b.metric[k].lo;require(std::isfinite(value),"finite actual source functional summary");cell[k]+=value;}
            }
        }
        MPI_Allreduce(MPI_IN_PLACE,&points,1,MPI_UINT64_T,MPI_SUM,layout.comm);require(points==expected_points,"all accepted source points accounted for");
        for(size_t c=0;c<candidate.size();c++){
            const auto&p=out.moment[c];require(p[0]>0&&p[3]>0,"positive source measure and mass");
            if(pass==0){candidate[c].delta=double((out.baseline[c].volume-p[1])/p[0]);candidate[c].delta_dot=double((out.baseline[c].vdot-p[2])/p[0]);}
            else candidate[c].mean_deltaP=double(p[4]/p[0]);
            require(std::isfinite(candidate[c].delta)&&std::isfinite(candidate[c].delta_dot)&&std::isfinite(candidate[c].mean_deltaP)&&std::isfinite(candidate[c].generating_delta),"finite actual source context, no clamp or heat fit");
        }
    }
    for(size_t c=0;c<out.moment.size();c++){
        const auto&p=out.moment[c];const auto&b=out.baseline[c];require(p[5]>0&&p[8]>=0,"positive reconstructed source heat and kinetic energy");
        out.mass_difference=std::max(out.mass_difference,fabsl(p[3]-b.mass)/b.mass);
        out.U_difference=std::max(out.U_difference,fabsl(p[5]-b.U)/b.U);
        out.Phi_difference=std::max(out.Phi_difference,fabsl(p[10]-b.Phi)/b.U);
        out.K_difference=std::max(out.K_difference,fabsl(p[8]-b.K)/b.U);
        out.L_difference=std::max(out.L_difference,fabsl(p[9]-b.L)/(2*b.K));
        out.pressure_difference=std::max(out.pressure_difference,fabsl(R(candidate[c].mean_deltaP)-b.mean_deltaP)/m.P0);
        require(fabsl(p[7]/b.volume-1)<1e-10L,"same canonical accepted corrected volume");
    }
    R checks[6]={out.mass_difference,out.U_difference,out.Phi_difference,out.K_difference,out.L_difference,out.pressure_difference};
    MPI_Allreduce(MPI_IN_PLACE,checks,6,MPI_LONG_DOUBLE,MPI_MAX,layout.comm);
    require(checks[0]<1e-12L&&checks[1]<1e-8L&&checks[2]<1e-8L&&checks[3]<1e-8L&&checks[4]<1e-11L&&checks[5]<1e-8L,"reconstructed current source matches saved accepted physical functional");
    out.mass_difference=checks[0];out.U_difference=checks[1];out.Phi_difference=checks[2];out.K_difference=checks[3];out.L_difference=checks[4];out.pressure_difference=checks[5];
    R min_delta=std::numeric_limits<R>::infinity();for(const auto&c:candidate)min_delta=std::min(min_delta,R(c.delta));
    MPI_Allreduce(MPI_IN_PLACE,&min_delta,1,MPI_LONG_DOUBLE,MPI_MIN,layout.comm);
    out.conditional_corrected_J_lower=nextafterl(geometry.minimum_canonical_J+min_delta,-INFINITY);
    require(out.conditional_corrected_J_lower>0,"corrected source measure remains positive under the canonical geometry certificate premises");
    require(!mf_evolved_stage(out.device->ptr,o.f.coefficient.data(),o.f.coefficient.size(),candidate.data(),candidate.size()),"verified actual source context kept on device for transfer");
    m.current_context=std::move(candidate);m.current_context_bound=true;out.wall_seconds=MPI_Wtime()-begin;return out;
}
CurrentContext derive_accepted_context(Model&,const Inputs&,const GeometryBound&);
} // namespace compact_source
#endif
