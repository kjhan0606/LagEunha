#define GFS_CURVED_PROOF_EXCEPTIONS
#include <climits>
#define GFS_CURVED_TARGET_INTEGRAL_LIBRARY
#include "test_curved_target_integral_v1.cpp"
#include "material_target_analytic.h"
namespace ta=material_target_analytic;
#include "material_whole_target.h"
namespace wt=whole_target;namespace wg=whole_geometry;

static std::vector<wt::Contribution> read_completed_source_proofs(const char*path,const wt::Input&in,int target,const std::vector<wg::Source>&geometry) {
    std::ifstream f(path);need(bool(f),"retained immutable source proof ledger");std::string line;std::getline(f,line);
    std::vector<wt::Contribution> out(in.packet.record.size());std::vector<std::array<bool,12>> seen(out.size());
    while(std::getline(f,line)){
        std::vector<std::string> c;size_t first=0;for(size_t i=0;i<=line.size();i++)if(i==line.size()||line[i]==','){c.push_back(line.substr(first,i-first));first=i+1;}
        need(c.size()==11,"retained source proof schema");size_t s=std::stoull(c[0]);need(s<out.size(),"retained original source id");int k=-1;for(int j=0;j<12;j++)if(c[7]==moment_names[j])k=j;need(k>=0&&!seen[s][k],"one retained proof per source and moment");
        auto&r=out[s];uint64_t cell=std::stoull(c[1]),points=std::stoull(c[5]),outgoing=std::stoull(c[6]);int classification=std::stoi(c[2]),complete=std::stoi(c[3]);
        bool firstrow=true;for(bool flag:seen[s])firstrow=firstrow&&!flag;
        if(firstrow){r.source=s;r.cell=cell;r.classification=classification;r.complete=complete;r.points=points;r.outgoing=outgoing;std::snprintf(r.reason,sizeof(r.reason),"%s",c[4].c_str());}
        need(r.cell==cell&&r.classification==classification&&r.complete==complete&&r.points==points&&r.outgoing==outgoing&&c[4]==r.reason,"consistent retained source metadata");
        need(cell==in.packet.record[s].cell&&classification==wg::classify(geometry[s],in.sites,in.weights,target),"fresh all384 classification agrees with retained immutable sources");
        if(complete){r.moment[k]=std::stold(c[8]);r.native_error[k]=std::stold(c[9]);char*end=nullptr;r.integral_error[k]=strtoflt128(c[10].c_str(),&end);need(end&&*end==0&&r.native_error[k]>=0&&r.integral_error[k]>=0,"roundtrip exact retained binary80/Q128 values");}
        else need(c[8]=="unproved"&&c[9]=="unproved"&&c[10]=="unproved","missing source proof stays UNPROVED");seen[s][k]=true;
    }
    for(const auto&row:seen)for(bool flag:row)need(flag,"all384 by12 retained proof records");return out;
}

static void write_packet(const std::string&path,const wt::TargetPacket&p) {
    FILE*f=fopen(path.c_str(),"wb");need(f,"canonical target packet output");
    uint64_t h[6]={0x4746575441523031ULL,1,p.epoch,p.target,12,10};need(fwrite(h,sizeof(h),1,f)==1,"target packet header");
    // x86 binary80 VALUE bytes only; six ABI padding bytes are not data.
    need(LDBL_MANT_DIG==64,"frozen binary80 value serialization");
    for(int k=0;k<12;k++)need(fwrite(&p.moment[k],10,1,f)==1&&fwrite(&p.absolute_error[k],10,1,f)==1,"target packet moments and errors");fclose(f);
}
int main(int argc,char**argv) {
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);active_mpi=true;
    int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    int result=2;
    try {
        need(argc==6||argc==7,"packet sites guard target-lower NEW-output [unchanged native source-proof ledger]");
        auto started=std::chrono::steady_clock::now();auto p=curved::read(argv[1]);auto sites=curved_partition::read_sites(argv[2]);auto weights=read_weights(argv[3],sites);std::string dest=argv[5];
        need(p.record.size()==384&&p.cells==24&&sites.size()==24,"same immutable original finite 24-cell packet");int target=-1;for(int j=0;j<int(sites.size());j++)if(sites[j].id==68526)target=j;need(target>=0,"fixed target68526");
        wt::Input input={p,sites,weights,1};R normalizers[3]={};
        {std::ifstream f(argv[4]);std::string line;std::getline(f,line);bool found=false;while(std::getline(f,line)){unsigned long long id;R U,M,A;if(sscanf(line.c_str(),"%llu,%Lg,%Lg,%Lg",&id,&U,&M,&A)==4&&id==68526){normalizers[0]=U;normalizers[1]=M;normalizers[2]=A;found=true;}}need(found&&normalizers[0]>0&&normalizers[1]>0&&normalizers[2]>0,"unchanged positive target lower bounds");}
        std::vector<wg::Source> geometry;geometry.reserve(p.record.size());for(const auto&r:p.record)geometry.push_back(wg::source(r,p.parameters[3]));
        std::string suffix=std::to_string(rank);std::ofstream sectors(dest+"/sectors_rank"+suffix+".csv");sectors<<"source,target,loop,boundary,sector,pass,reason,contraction,invariance\n";
        FILE*points=fopen((dest+"/positive_rank"+suffix+".bin").c_str(),"wb");need(points,"rank positive packet");uint64_t h[6]={0x474657504f533031ULL,1,uint64_t(rank),uint64_t(ranks),sizeof(Positive),uint64_t(p.record.size())};need(fwrite(h,sizeof(h),1,points)==1,"rank positive header");
        wt::Logs logs={sectors,points};std::vector<wt::Contribution> local(p.record.size());
        if(argc==7){auto retained=read_completed_source_proofs(argv[6],input,target,geometry);for(size_t s=rank;s<p.record.size();s+=ranks)local[s]=retained[s];}
        else for(size_t s=rank;s<p.record.size();s+=ranks)local[s]=wt::integrate(input,s,target,geometry[s],logs);
        sectors.close();fclose(points);
        std::vector<wt::Contribution> gathered(rank==0?p.record.size()*ranks:0);
        need(local.size()*sizeof(wt::Contribution)<INT_MAX,"bounded native MPI message");
        MPI_Gather(local.data(),local.size()*sizeof(wt::Contribution),MPI_BYTE,gathered.data(),local.size()*sizeof(wt::Contribution),MPI_BYTE,0,MPI_COMM_WORLD);
        if(rank==0) {
            std::vector<wt::Contribution> rows(p.record.size());for(size_t s=0;s<rows.size();s++)rows[s]=gathered[(s%ranks)*rows.size()+s];
            std::ofstream source(dest+"/source_contributions.csv");source<<std::setprecision(21)<<"source,cell,classification,complete,reason,points,outgoing,moment,native,native_geometry_nodes_sum_error,analytic_remainder\n";
            for(const auto&r:rows)for(int k=0;k<12;k++)source<<r.source<<','<<r.cell<<','<<r.classification<<','<<r.complete<<','<<r.reason<<','<<r.points<<','<<r.outgoing<<','<<moment_names[k]<<','<<(r.complete?qs(Q(r.moment[k])):"unproved")<<','<<(r.complete?qs(Q(r.native_error[k])):"unproved")<<','<<(r.complete?qs(r.integral_error[k]):"unproved")<<'\n';source.close();
            wt::TargetPacket caller={};caller.epoch=777;caller.target=999;for(int k=0;k<12;k++){caller.moment[k]=k+1;caller.absolute_error[k]=0;}wt::TargetPacket before;std::memcpy(&before,&caller,sizeof(caller));
            auto epoch=wt::assemble(input,2,target,rows,geometry,normalizers,dest,caller);
            bool epoch_rollback=!epoch.published&&!std::memcmp(&before,&caller,sizeof(caller));
            auto missing=rows;missing.pop_back();auto absent=wt::assemble(input,1,target,missing,geometry,normalizers,dest,caller);
            bool missing_rollback=!absent.published&&!std::memcmp(&before,&caller,sizeof(caller));
            auto report=wt::assemble(input,1,target,rows,geometry,normalizers,dest,caller);
            bool rollback=report.published||!std::memcmp(&before,&caller,sizeof(caller));
            std::ofstream ledger(dest+"/whole_target_ledger.csv");ledger<<std::setprecision(21)<<"target,moment,candidate,native_geometry_node_source_sum_error,analytic_remainder,canonical_seam_bound,vertex_cap_bound,final_sum_rounding,total_absolute_error,unchanged_physical_tolerance,error_over_tolerance,all_sources_complete,geometry_complete,pass\n";
            for(int k=0;k<12;k++){bool complete=report.sources_complete&&report.geometry_complete;ledger<<68526<<','<<moment_names[k]<<','<<qs(Q(report.candidate.moment[k]))<<','<<qs(report.native[k])<<','<<(report.sources_complete?qs(report.integral[k]):"unproved")<<','<<(report.geometry_complete?qs(report.geometry.seam[k]):"unproved")<<','<<(report.geometry_complete?qs(report.geometry.cap[k]):"unproved")<<','<<qs(report.sum[k])<<','<<(complete?qs(report.error[k]):"unproved")<<','<<qs(report.tolerance[k])<<','<<(complete?qs(report.error[k]/report.tolerance[k]):"unproved")<<','<<report.sources_complete<<','<<report.geometry_complete<<','<<(complete&&report.error[k]<=report.tolerance[k])<<'\n';}ledger.close();
            std::ofstream duplicated(dest+"/canonical_vs_duplicated_seams.csv");duplicated<<"moment,canonical_once,source_local_duplicate_diagnostic\n";if(report.geometry_complete)for(int k=0;k<12;k++)duplicated<<moment_names[k]<<','<<qs(report.geometry.seam[k])<<','<<qs(report.geometry.duplicated_seam[k])<<'\n';duplicated.close();
            write_packet(dest+"/caller_after.bin",caller);if(report.published)write_packet(dest+"/validated_target.bin",caller);
            std::ofstream meta(dest+"/native_contract.json");meta<<"{\n  \"target\":68526,\n  \"epoch\":1,\n  \"MPI_ranks\":"<<ranks<<",\n  \"retained_source_certificates\":"<<(argc==7?"true":"false")<<",\n  \"source_owner\":\"source_id % MPI_ranks\",\n  \"seam_owner\":\"rank0 canonical sorted reference edge then fixed part\",\n  \"final_sum_order\":\"ascending original source id, actual binary80 scalar addition\",\n  \"all384_source_rows\":"<<rows.size()<<",\n  \"sources_complete\":"<<(report.sources_complete?"true":"false")<<",\n  \"failed_sources\":"<<report.failed_sources<<",\n  \"whole_image_exclusions\":"<<report.excluded<<",\n  \"whole_inclusions\":"<<report.inside<<",\n  \"intersections_or_certified_empty\":"<<report.intersected<<",\n  \"positive_points\":"<<report.points<<",\n  \"geometry_complete\":"<<(report.geometry_complete?"true":"false")<<",\n  \"source_pairs_examined\":"<<report.geometry.pairs<<",\n  \"canonical_target_seam_parts\":"<<report.geometry.seam_parts<<",\n  \"external_reference_edges\":"<<report.geometry.external<<",\n  \"relevant_external_edges\":"<<report.geometry.relevant_external<<",\n  \"epoch_reject_preserves_output\":"<<(epoch_rollback?"true":"false")<<",\n  \"missing_row_reject_preserves_output\":"<<(missing_rollback?"true":"false")<<",\n  \"actual_failure_preserves_output\":"<<(rollback?"true":"false")<<",\n  \"physical_budgets_pass\":"<<(report.budgets_pass?"true":"false")<<",\n  \"positive_mass_and_U\":"<<(report.positive?"true":"false")<<",\n  \"published_one_target\":"<<(report.published?"true":"false")<<",\n  \"witness\":\""<<report.witness<<"\",\n  \"backend\":\"CPU binary80 integration, Q128 geometry; IEEE nearest / <=8ulp engineering libm model\",\n  \"CUDA_execution\":false,\n  \"exact_global_C0_proved\":false,\n  \"physical_time_evolution\":false,\n  \"full_remap_proved\":false,\n  \"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"\n}\n";meta.close();
            result=report.published&&epoch_rollback&&missing_rollback&&rollback?0:2;
            printf("{\"target\":68526,\"sources\":%zu,\"failed_sources\":%llu,\"geometry_complete\":%s,\"published\":%s,\"witness\":\"%s\",\"seconds\":%.17g}\n",rows.size(),(unsigned long long)report.failed_sources,report.geometry_complete?"true":"false",report.published?"true":"false",report.witness.c_str(),std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count());
        }
        MPI_Bcast(&result,1,MPI_INT,0,MPI_COMM_WORLD);
    }catch(const std::exception&e){fprintf(stderr,"rank%d FAIL immutable metadata/preflight: %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);}
    active_mpi=false;MPI_Finalize();return result;
}
