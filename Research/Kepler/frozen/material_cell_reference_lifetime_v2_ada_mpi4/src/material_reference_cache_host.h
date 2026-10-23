#ifndef MATERIAL_REFERENCE_CACHE_HOST_H
#define MATERIAL_REFERENCE_CACHE_HOST_H
// Included after the original full_evolution definitions. The immutable cache
// belongs to one live source/basis/rule/generating context, never to x/v at a stage.
#include <iomanip>
namespace reference_cache {
namespace fe=full_evolution;
using U=uint64_t;
struct Cache {
    fe::Material&material;const md::Ownership&layout;void*physical=nullptr,*ptr=nullptr;
    double*host=nullptr;bool resident=false;U generation=0,points=0,total_doubles=0,device_bytes=0,key_hash=0;
    size_t max_points=0,max_doubles=0;std::vector<fe::PointTile>tiles;std::vector<U>offset;
    std::vector<std::vector<U>>prefix;std::string binding;
    double generation_seconds=0,wall_seconds=0;bool inputs_released=false;
    Cache(fe::Material&m,const md::Ownership&l,Device&gpu,const GfsEnrichConfig&cfg,const std::string&key,const std::string&directory,const cell_native::DevicePlan&solver_plan):material(m),layout(l),physical(gpu.ptr),generation(m.positive.generation){
        double begin=MPI_Wtime();need(!active.self&&m.positive.loaded&&m.positive.active,"one immutable selected reference cache in this native context");
        need(sizeof(GfsMaterialPoint12)==104&&sizeof(STJet)==24&&sizeof(double)==8,"packed reference/current ABI");
        tiles=fe::point_tiles(m,cfg);offset.resize(m.owned.f.triangle.size()+1);prefix.resize(tiles.size());
        for(size_t s=0;s<m.owned.f.triangle.size();s++){
            const auto&t=m.owned.f.triangle[s];U n=m.positive.count[s],stride=6+3*(t.n-10);need(n>0&&stride<=30,"original selected point layout");offset[s+1]=offset[s]+n*stride;points+=n;
        }
        total_doubles=offset.back();for(size_t i=0;i<tiles.size();i++){
            const auto&t=tiles[i];prefix[i].resize(t.count+1);
            for(size_t j=0;j<=t.count;j++)prefix[i][j]=offset[t.first+j]-offset[t.first];
            max_points=std::max<size_t>(max_points,t.prefix.back());max_doubles=std::max<size_t>(max_doubles,prefix[i].back());
        }
        U counts[2]={points,total_doubles};MPI_Allreduce(MPI_IN_PLACE,counts,2,MPI_UINT64_T,MPI_SUM,l.comm);
        need(counts[0]==249691129&&counts[1]*8==59891910888ULL,"same whole original selected point and exact enriched cache bytes");
        std::ostringstream b;b<<"IMMUTABLE_REFERENCE_CACHE_LIFETIME_V2\n"<<key<<"\nABI double8 STJet24 current112 P3-order-original\n";
        b<<"rule_generation "<<generation<<"\nlocal_triangles "<<m.owned.f.triangle.size()<<"\nrank "<<l.rank<<"\ngenerating_coeff_hash "<<fe::hash(m.generating.data(),m.generating.size()*sizeof(STCoefficient))<<"\ngenerating_delta_hash "<<fe::hash(m.generating_delta.data(),m.generating_delta.size()*sizeof(m.generating_delta[0]))<<"\n";
        cudaDeviceProp prop;cell_native::gpu(cudaGetDeviceProperties(&prop,0),"reference cache actual device ABI");b<<"compute_capability "<<prop.major<<" "<<prop.minor<<"\n";binding=b.str();key_hash=fe::hash(binding.data(),binding.size());
        // All four ranks must be charged to the same actual device UUID.
        // No assumption that one MPI rank means one distinct GPU is allowed.
        std::vector<char>uuids(16*l.ranks);
        MPI_Allgather(prop.uuid.bytes,16,MPI_CHAR,uuids.data(),16,MPI_CHAR,l.comm);
        for(int rank=0;rank<l.ranks;rank++)need(memcmp(uuids.data(),uuids.data()+16*rank,16)==0,"one actual GPU UUID for this original MPI4 candidate");
        std::ostringstream uuid;uuid<<std::hex<<std::setfill('0');for(int k=0;k<16;k++)uuid<<std::setw(2)<<unsigned(static_cast<unsigned char>(prop.uuid.bytes[k]));
        U shared=max_points*112+(MF_BATCH_TRIANGLES+1)*12+MF_BATCH_TRIANGLES*4;
        U resident_local=total_doubles*8+shared,tiled_local=max_doubles*8+shared,resident_all=0,tiled_all=0;
        MPI_Allreduce(&resident_local,&resident_all,1,MPI_UINT64_T,MPI_SUM,l.comm);MPI_Allreduce(&tiled_local,&tiled_all,1,MPI_UINT64_T,MPI_SUM,l.comm);
        need(m.positive.index.size()==points&&m.positive.weight.size()==points,"every selected reference generation input accounted for");
        U release_local=points*12,release_all=0;MPI_Allreduce(&release_local,&release_all,1,MPI_UINT64_T,MPI_SUM,l.comm);
        MPI_Barrier(l.comm);size_t free_before,total;cell_native::gpu(cudaMemGetInfo(&free_before,&total),"generation phase actual memory before cache or solver");
        U free_min=free_before;MPI_Allreduce(MPI_IN_PLACE,&free_min,1,MPI_UINT64_T,MPI_MIN,l.comm);free_before=free_min;
        int can_reside=free_before>=total/5&&resident_all<=free_before-total/5&&free_before+release_all>=resident_all+solver_plan.aggregate+total/5,all;
        MPI_Allreduce(&can_reside,&all,1,MPI_INT,MPI_MIN,l.comm);resident=all;
        U actual_all=resident?resident_all:tiled_all;
        U phase1_free=free_before>=actual_all?free_before-actual_all:0;
        U phase2_free=free_before+release_all>=actual_all+solver_plan.aggregate?free_before+release_all-actual_all-solver_plan.aggregate:0;
        int fits=phase1_free>=total/5&&phase2_free>=total/5;
        MPI_Allreduce(&fits,&all,1,MPI_INT,MPI_MIN,l.comm);
        std::ifstream mem("/proc/meminfo");std::string label,rest;U available=0,n;
        while(mem>>label>>n){std::getline(mem,rest);if(label=="MemAvailable:"){available=n*1024;break;}}
        MPI_Allreduce(MPI_IN_PLACE,&available,1,MPI_UINT64_T,MPI_MIN,l.comm);
        U pinned_all=resident?0:counts[1]*8;int host_fits=available>=pinned_all+32ULL*1024*1024*1024,host_all;
        MPI_Allreduce(&host_fits,&host_all,1,MPI_INT,MPI_MIN,l.comm);
        if(!l.rank){FILE*f=fopen((directory+"/two_phase_memory_plan.json").c_str(),"wx");need(f,"new preallocation two-phase device and host accounting");
            fprintf(f,"{\"gpu_uuid\":\"%s\",\"all_MPI_ranks_on_this_UUID\":%d,\"compiled_tile_triangles\":%d,\"point_capacity\":%zu,\"total_device_bytes\":%llu,\"free_after_actual_physical_context_bytes\":%llu,\"planned_cache_bytes_all_ranks\":%llu,\"planned_solver_bytes_all_ranks\":%llu,\"generation_input_bytes_to_release_all_ranks\":%llu,\"predicted_generation_phase_free_bytes\":%llu,\"predicted_normal_stage_phase_free_bytes\":%llu,\"required_free_device_bytes\":%llu,\"host_available_bytes\":%llu,\"host_pinned_reference_bytes\":%llu,\"host_workspace_reserve_bytes\":34359738368,\"both_device_phases_pass\":%s,\"host_capacity_pass\":%s,\"cache_generation_started\":false,\"solver_device_allocation_started\":false,\"physical_trajectory_advanced\":false}\n",uuid.str().c_str(),l.ranks,MF_BATCH_TRIANGLES,m.cache_capacity,(unsigned long long)total,(unsigned long long)free_before,(unsigned long long)actual_all,(unsigned long long)solver_plan.aggregate,(unsigned long long)release_all,(unsigned long long)phase1_free,(unsigned long long)phase2_free,(unsigned long long)(total/5),(unsigned long long)available,(unsigned long long)pinned_all,all?"true":"false",host_all?"true":"false");material_full::durable(f);}
        need(all&&host_all,"both actual UUID device phases keep original twenty percent headroom and host cache/workspace fit");
        if(!resident)cell_native::gpu(cudaMallocHost(reinterpret_cast<void**>(&host),total_doubles*8),"rank-owned immutable pinned host reference cache");
        cell_native::gpu(mfr_create(physical,max_points,max_doubles,total_doubles,resident?1:0,host,&ptr),"original immutable/dynamic split cache allocation");device_bytes=mfr_device_bytes(ptr);need(device_bytes==(resident?resident_local:tiled_local),"exact reference cache CUDA allocation accounting");
        MPI_Barrier(l.comm);size_t free_after;cell_native::gpu(cudaMemGetInfo(&free_after,&total),"cache actual allocated memory");fits=free_after>=total/5;MPI_Allreduce(&fits,&all,1,MPI_INT,MPI_MIN,l.comm);need(all,"actual generation phase before solver keeps twenty percent device memory free");
        std::vector<MFEvolvedCell>context(m.owned.f.cell.size());for(size_t c=0;c<context.size();c++)context[c].generating_delta=double(m.generating_delta[c][0]);
        cell_native::gpu(mf_evolved_stage(physical,m.owned.f.coefficient.data(),m.owned.f.coefficient.size(),context.data(),context.size()),"same immutable generating delta before one reference pass");
        FILE*progress=fopen((directory+"/reference_progress_r"+std::to_string(l.rank)+".tsv").c_str(),"wx");need(progress,"new immutable reference generation progress");fprintf(progress,"tiles_done\ttiles_total\tpoints_done\tCUDA_seconds\n");
        U done=0;for(size_t i=0;i<tiles.size();i++){
            const auto&t=tiles[i];double seconds;cell_native::gpu(mfr_generate(physical,ptr,t.first,t.count,offset[t.first],t.prefix.data(),prefix[i].data(),&seconds),"every immutable original enriched point and generating j0");
            done+=t.prefix.back();generation_seconds+=seconds;fprintf(progress,"%zu\t%zu\t%llu\t%.9g\n",i+1,tiles.size(),(unsigned long long)done,generation_seconds);fflush(progress);
        }
        material_full::durable(progress);need(done==points,"each original selected point cached once");MPI_Barrier(l.comm);wall_seconds=MPI_Wtime()-begin;
        FILE*keyfile=fopen((directory+"/reference_binding_r"+std::to_string(l.rank)+".txt").c_str(),"wx");need(keyfile,"new immutable source/basis/rule/thermal/ABI binding");need(fwrite(binding.data(),1,binding.size(),keyfile)==binding.size(),"reference binding write");material_full::durable(keyfile);
        if(!l.rank){FILE*f=fopen((directory+"/reference_cache.json").c_str(),"wx");need(f,"new actual reference cache report");fprintf(f,"{\"all_original_selected_points_built\":true,\"selected_points\":%llu,\"immutable_double_bytes_all_ranks\":%llu,\"mode\":\"%s\",\"new_GPU_bytes_all_ranks\":%llu,\"device_free_after_allocation_bytes\":%llu,\"device_total_bytes\":%llu,\"twenty_percent_headroom_pass\":true,\"same_original_complete_source_tiles\":true,\"current_state_recomputed_each_stage\":true,\"source_basis_rule_thermal_and_ABI_binding\":true,\"rank0_key_hash\":%llu,\"generation_GPU_seconds_rank0\":%.9g,\"setup_wall_seconds_rank0\":%.9g,\"solver_allocation_delayed\":true,\"generation_input_release_pending\":true,\"cache_direct_comparison_pending\":true,\"physical_trajectory_advanced\":false}\n",(unsigned long long)counts[0],(unsigned long long)(counts[1]*8),resident?"device_resident":"fixed_source_tiles_from_rank_owned_pinned_host",(unsigned long long)actual_all,(unsigned long long)free_after,(unsigned long long)total,(unsigned long long)key_hash,generation_seconds,wall_seconds);material_full::durable(f);}
        active={this,build_bridge,batch_bridge};
    }
    void release_generation_inputs(const std::string&directory){
        valid(physical);need(!inputs_released&&material.cache_checked&&active.self==this,"complete original cache and same-stage direct witness before any input release");
        MPI_Barrier(layout.comm);size_t before,total;cell_native::gpu(cudaMemGetInfo(&before,&total),"memory before generation inputs release");
        U freed=0;cell_native::gpu(mfr_release_generation_inputs(physical,ptr,points,&freed),"enter reference-only mode and release only original selected GPU IDs/weights");
        need(freed==points*12,"exact released generation-only payload bytes");U all_freed=0;MPI_Allreduce(&freed,&all_freed,1,MPI_UINT64_T,MPI_SUM,layout.comm);
        MPI_Barrier(layout.comm);size_t after;cell_native::gpu(cudaMemGetInfo(&after,&total),"memory after every rank generation input release");
        int fits=after>=total/5,all;MPI_Allreduce(&fits,&all,1,MPI_INT,MPI_MIN,layout.comm);need(all,"actual post-release phase retains original headroom");inputs_released=true;
        if(!layout.rank){FILE*f=fopen((directory+"/generation_input_release.json").c_str(),"wx");need(f,"new completed generation/witness/release ownership evidence");
            fprintf(f,"{\"all_original_reference_points_generated\":true,\"same_actual_stage_direct_witness_pass\":true,\"freed_GPU_ID_weight_bytes_all_ranks\":%llu,\"actual_device_free_before_release_bytes\":%llu,\"actual_device_free_after_all_rank_release_bytes\":%llu,\"host_rule_and_GPU_counts_offsets_retained\":true,\"direct_selected_point_entry_disabled\":true,\"current_stage_assembled_F_not_recomputed\":true,\"saved_k1_not_recomputed\":true,\"solver_allocation_next\":true}\n",(unsigned long long)all_freed,(unsigned long long)before,(unsigned long long)after);material_full::durable(f);}
    }
    ~Cache(){if(active.self==this)active={};mfr_destroy(ptr);if(host)cudaFreeHost(host);}
    Cache(const Cache&)=delete;Cache&operator=(const Cache&)=delete;
    void valid(void*ctx)const{need(ctx==physical&&material.positive.loaded&&material.positive.active&&material.positive.generation==generation&&material.owned.f.triangle.size()+1==offset.size(),"immutable reference binding at actual RK stage");}
    static int build_bridge(void*self,void*ctx,size_t first,size_t count,const uint32_t*pp,MFEvolvedBlock*out,double*seconds){
        auto&s=*static_cast<Cache*>(self);s.valid(ctx);auto it=std::lower_bound(s.tiles.begin(),s.tiles.end(),first,[](const fe::PointTile&t,size_t n){return t.first<n;});
        need(it!=s.tiles.end()&&it->first==first&&it->count==count&&memcmp(pp,it->prefix.data(),(count+1)*sizeof(uint32_t))==0,"exact original fixed source tile and positive point order");size_t i=it-s.tiles.begin();
        return mfr_build(ctx,s.ptr,first,count,s.offset[first],pp,s.prefix[i].data(),out,seconds);
    }
    static int batch_bridge(void*self,void*ctx,int pass,size_t first,size_t count,MFEvolvedBlock*out,double*seconds){auto&s=*static_cast<Cache*>(self);s.valid(ctx);return mfr_batch(ctx,s.ptr,pass,first,count,out,seconds);}
};
}
#endif
