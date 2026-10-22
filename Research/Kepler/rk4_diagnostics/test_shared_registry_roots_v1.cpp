#define GFS_SHARED_FUNCTION_LIBRARY
#include "test_shared_function_v1.cpp"
#include "material_shared_retained.h"
#include <tuple>

namespace target=shared_target;namespace sd=shared_duffy;
namespace {
using Key=std::tuple<int,uint64_t,uint64_t>;
struct Job {uint64_t id;int edge;uint64_t first,second;bool retained;};
static std::string quoted(const char*text){std::string s="\"";for(const char*p=text;*p;p++){if(*p=='\"')s+='\"';s+=*p;}return s+'\"';}
}
int main(int argc,char**argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);int result=2;
    try{
        need(argc==8,"original packet saved sites weights audited_registry retained_roots NEWoutput");
        auto begin=std::chrono::steady_clock::now();auto p=curved::read(argv[1]);const auto f=load_saved(std::string(argv[2])+"/shared_state.bin",p);auto coefficients=f.coefficient;
        auto site=sites(argv[3],argv[4]);std::string registry=argv[5],retained=argv[6],dest=argv[7];
        std::ifstream terminal(registry+"/native_exit_code.txt");int exit=-1;need(bool(terminal>>exit)&&exit==0,"registry preflight terminal success");
        std::vector<Job>jobs;std::string line;std::ifstream plan(registry+"/canonical_root_work_registry.csv");need(bool(plan),"immutable deduplicated root work registry");std::getline(plan,line);
        std::set<Key>unique;
        while(std::getline(plan,line)){auto c=target::columns(line);Job job={std::stoull(c[0]),std::stoi(c[1]),std::stoull(c[2]),std::stoull(c[3]),c[4]=="1"};need(job.id==jobs.size()&&job.first<job.second&&unique.emplace(job.edge,job.first,job.second).second&&std::stoi(c[5])==(job.retained?0:256),"fixed original canonical work identity and cap");jobs.push_back(job);}
        need(jobs.size()==2452,"audited 2289 new plus163 reused jobs");
        std::map<Key,target::Roots>old;std::map<Key,uint64_t>old_count;std::ifstream input(retained+"/canonical_event_root_packets.csv");need(bool(input),"same-state original retained root packets");std::getline(input,line);
        while(std::getline(input,line)){
            auto c=target::columns(line);uint64_t other=std::stoull(c[1]);Key key(std::stoi(c[0]),std::min(uint64_t(68526),other),std::max(uint64_t(68526),other));auto&r=old[key];
            need(c[2]=="1","retained exact event complete");r.complete=true;r.evaluations=std::stoull(c[3]);uint64_t count=std::stoull(c[5]);
            if(old_count.count(key))need(old_count[key]==count,"consistent retained count");else old_count[key]=count;
            if(c[4]!="-1"){need(std::stoull(c[4])==r.root.size(),"retained original root ordering");r.root.push_back({S(target::number(c[6]),target::number(c[7])),std::stoi(c[8])});}
            snprintf(r.reason,sizeof(r.reason),"retained_same_function_unordered_plane_sign_invariant");
        }
        for(const auto&entry:old)need(entry.second.root.size()==old_count[entry.first],"all retained roots recovered without evaluation");
        std::vector<sd::Compiled>maps(f.triangle.size());int bad=0;
        #pragma omp parallel for schedule(static) reduction(+:bad)
        for(size_t s=0;s<maps.size();s++)try{maps[s]=sd::compile(f,s);}catch(const std::exception&){bad++;}
        need(bad==0,"same canonical stable Duffy compilation for edge functions");auto edges=target::topology(maps);
        std::map<std::pair<uint64_t,uint64_t>,target::Plane>plane;
        for(const auto&s:site)for(const auto&entry:target::planes(site,s.id))if(s.id<entry.competitor)plane.emplace(std::make_pair(s.id,entry.competitor),entry);
        std::vector<target::Roots>outcome(jobs.size());
        for(const auto&job:jobs)need(edges.count(job.edge)&&plane.count({job.first,job.second})&&(!job.retained||old.count({job.edge,job.first,job.second})),"every job binds one same canonical edge and exact pair-plane");
        #pragma omp parallel for schedule(dynamic,1)
        for(size_t i=0;i<jobs.size();i++)if(i%ranks==size_t(rank)){
            const auto&job=jobs[i];if(job.retained)outcome[i]=old.at({job.edge,job.first,job.second});
            else outcome[i]=target::isolate(edges.at(job.edge),plane.at({job.first,job.second}));
        }
        FILE*summary=fopen((dest+"/canonical_root_results_rank"+std::to_string(rank)+".csv").c_str(),"w");
        FILE*packets=fopen((dest+"/canonical_event_root_packets_rank"+std::to_string(rank)+".csv").c_str(),"w");
        FILE*witness=fopen((dest+"/failed_root_windows_rank"+std::to_string(rank)+".csv").c_str(),"w");need(summary&&packets&&witness,"actual roots and fixed-cap failure witness outputs");
        fprintf(summary,"work_id,canonical_edge,first_target,second_target,retained,complete,root_windows,roots,endpoint_roots,owner_rank,reason\n");
        fprintf(packets,"work_id,canonical_edge,first_target,second_target,retained,complete,root_index,root_count,t,t_radius,endpoint\n");
        fprintf(witness,"work_id,canonical_edge,first_target,second_target,last_lo,last_hi,exists,f,f_radius,df,df_radius,left_f,left_radius,right_f,right_radius,reason\n");
        uint64_t counts[6]={};
        for(size_t i=0;i<jobs.size();i++)if(i%ranks==size_t(rank)){
            const auto&job=jobs[i];const auto&r=outcome[i];uint64_t endpoint=0;for(const auto&root:r.root)endpoint+=root.endpoint>=0;
            counts[0]+=!job.retained;counts[1]+=job.retained;counts[2]+=!r.complete;counts[3]+=job.retained?0:r.evaluations;counts[4]+=r.root.size();counts[5]+=endpoint;
            fprintf(summary,"%llu,%d,%llu,%llu,%d,%d,%llu,%zu,%llu,%d,%s\n",(unsigned long long)job.id,job.edge,(unsigned long long)job.first,(unsigned long long)job.second,job.retained,r.complete,(unsigned long long)r.evaluations,r.root.size(),(unsigned long long)endpoint,rank,quoted(r.reason).c_str());
            if(r.root.empty())fprintf(packets,"%llu,%d,%llu,%llu,%d,%d,-1,0,0,0,-1\n",(unsigned long long)job.id,job.edge,(unsigned long long)job.first,(unsigned long long)job.second,job.retained,r.complete);
            else for(size_t k=0;k<r.root.size();k++)fprintf(packets,"%llu,%d,%llu,%llu,%d,%d,%zu,%zu,%s,%s,%d\n",(unsigned long long)job.id,job.edge,(unsigned long long)job.first,(unsigned long long)job.second,job.retained,r.complete,k,r.root.size(),qs(r.root[k].t.n).c_str(),qs(r.root[k].t.e).c_str(),r.root[k].endpoint);
            if(!r.complete)fprintf(witness,"%llu,%d,%llu,%llu,%s,%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",(unsigned long long)job.id,job.edge,(unsigned long long)job.first,(unsigned long long)job.second,qs(r.last_lo).c_str(),qs(r.last_hi).c_str(),r.last_exists,qs(r.last_f.n).c_str(),qs(r.last_f.e).c_str(),qs(r.last_df.n).c_str(),qs(r.last_df.e).c_str(),qs(r.last_fl.n).c_str(),qs(r.last_fl.e).c_str(),qs(r.last_fr.n).c_str(),qs(r.last_fr.e).c_str(),quoted(r.reason).c_str());
        }
        fclose(summary);fclose(packets);fclose(witness);MPI_Allreduce(MPI_IN_PLACE,counts,6,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
        need(!memcmp(coefficients.data(),f.coefficient.data(),coefficients.size()*sizeof(STCoefficient)),"root geometry does not change state");
        if(rank==0){checkpoint(dest+"/caller_unchanged_state.bin",f,1);printf("{\"scope\":\"ONLY audited deduplicated canonical edge/power-plane roots on same original state, fixed256 budget; no quadrature/evolution\",\"MPI_ranks\":%d,\"OMP_threads\":%d,\"planned_jobs\":%zu,\"new_jobs\":%llu,\"retained_jobs\":%llu,\"incomplete_jobs\":%llu,\"new_root_windows\":%llu,\"retained_and_new_roots\":%llu,\"endpoint_roots\":%llu,\"all_registered_roots_complete\":%s,\"actual_all_target_domains_complete\":false,\"source_partition_complete\":false,\"new_physical_quadrature_queries\":0,\"full_quadrature_allowed\":false,\"accepted_RK4_evolution\":false,\"long_Kepler_validated\":false,\"seconds\":%.17g}\n",ranks,omp_get_max_threads(),jobs.size(),(unsigned long long)counts[0],(unsigned long long)counts[1],(unsigned long long)counts[2],(unsigned long long)counts[3],(unsigned long long)counts[4],(unsigned long long)counts[5],counts[2]?"false":"true",std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count());}result=counts[2]?2:0;
    }catch(const std::exception&e){fprintf(stderr,"rank%d FAIL %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);}MPI_Finalize();return result;
}
