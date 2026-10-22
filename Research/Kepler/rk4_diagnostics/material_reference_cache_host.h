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
    double generation_seconds=0,wall_seconds=0;
    Cache(fe::Material&m,const md::Ownership&l,Device&gpu,const GfsEnrichConfig&cfg,const std::string&key,const std::string&directory):material(m),layout(l),physical(gpu.ptr),generation(m.positive.generation){
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
        std::ostringstream b;b<<"IMMUTABLE_REFERENCE_CACHE_V1\n"<<key<<"\nABI double8 STJet24 current112 P3-order-original\n";
        b<<"rule_generation "<<generation<<"\nlocal_triangles "<<m.owned.f.triangle.size()<<"\nrank "<<l.rank<<"\ngenerating_coeff_hash "<<fe::hash(m.generating.data(),m.generating.size()*sizeof(STCoefficient))<<"\ngenerating_delta_hash "<<fe::hash(m.generating_delta.data(),m.generating_delta.size()*sizeof(m.generating_delta[0]))<<"\n";
        cudaDeviceProp prop;cell_native::gpu(cudaGetDeviceProperties(&prop,0),"reference cache actual device ABI");b<<"compute_capability "<<prop.major<<" "<<prop.minor<<"\n";binding=b.str();key_hash=fe::hash(binding.data(),binding.size());
        U shared=max_points*112+(MF_BATCH_TRIANGLES+1)*12+MF_BATCH_TRIANGLES*4;
        U resident_local=total_doubles*8+shared,tiled_local=max_doubles*8+shared,resident_all=0,tiled_all=0;
        MPI_Allreduce(&resident_local,&resident_all,1,MPI_UINT64_T,MPI_SUM,l.comm);MPI_Allreduce(&tiled_local,&tiled_all,1,MPI_UINT64_T,MPI_SUM,l.comm);
        MPI_Barrier(l.comm);size_t free_before,total;cell_native::gpu(cudaMemGetInfo(&free_before,&total),"cache fresh memory capacity");
        int can_reside=free_before>=total/5&&resident_all<=free_before-total/5,all;MPI_Allreduce(&can_reside,&all,1,MPI_INT,MPI_MIN,l.comm);resident=all;
        U actual_all=resident?resident_all:tiled_all;int fits=free_before>=total/5&&actual_all<=free_before-total/5;MPI_Allreduce(&fits,&all,1,MPI_INT,MPI_MIN,l.comm);need(all,"fixed same-source cache fits actual device with twenty percent headroom");
        if(!resident)cell_native::gpu(cudaMallocHost(reinterpret_cast<void**>(&host),total_doubles*8),"rank-owned immutable pinned host reference cache");
        cell_native::gpu(mfr_create(physical,max_points,max_doubles,total_doubles,resident?1:0,host,&ptr),"original immutable/dynamic split cache allocation");device_bytes=mfr_device_bytes(ptr);need(device_bytes==(resident?resident_local:tiled_local),"exact reference cache CUDA allocation accounting");
        MPI_Barrier(l.comm);size_t free_after;cell_native::gpu(cudaMemGetInfo(&free_after,&total),"cache actual allocated memory");fits=free_after>=total/5;MPI_Allreduce(&fits,&all,1,MPI_INT,MPI_MIN,l.comm);need(all,"actual cell and reference cache keeps twenty percent device memory free");
        std::vector<MFEvolvedCell>context(m.owned.f.cell.size());for(size_t c=0;c<context.size();c++)context[c].generating_delta=double(m.generating_delta[c][0]);
        cell_native::gpu(mf_evolved_stage(physical,m.owned.f.coefficient.data(),m.owned.f.coefficient.size(),context.data(),context.size()),"same immutable generating delta before one reference pass");
        FILE*progress=fopen((directory+"/reference_progress_r"+std::to_string(l.rank)+".tsv").c_str(),"wx");need(progress,"new immutable reference generation progress");fprintf(progress,"tiles_done\ttiles_total\tpoints_done\tCUDA_seconds\n");
        U done=0;for(size_t i=0;i<tiles.size();i++){
            const auto&t=tiles[i];double seconds;cell_native::gpu(mfr_generate(physical,ptr,t.first,t.count,offset[t.first],t.prefix.data(),prefix[i].data(),&seconds),"every immutable original enriched point and generating j0");
            done+=t.prefix.back();generation_seconds+=seconds;fprintf(progress,"%zu\t%zu\t%llu\t%.9g\n",i+1,tiles.size(),(unsigned long long)done,generation_seconds);fflush(progress);
        }
        material_full::durable(progress);need(done==points,"each original selected point cached once");MPI_Barrier(l.comm);wall_seconds=MPI_Wtime()-begin;
        FILE*keyfile=fopen((directory+"/reference_binding_r"+std::to_string(l.rank)+".txt").c_str(),"wx");need(keyfile,"new immutable source/basis/rule/thermal/ABI binding");need(fwrite(binding.data(),1,binding.size(),keyfile)==binding.size(),"reference binding write");material_full::durable(keyfile);
        if(!l.rank){FILE*f=fopen((directory+"/reference_cache.json").c_str(),"wx");need(f,"new actual reference cache report");fprintf(f,"{\"all_original_selected_points_built\":true,\"selected_points\":%llu,\"immutable_double_bytes_all_ranks\":%llu,\"mode\":\"%s\",\"new_GPU_bytes_all_ranks\":%llu,\"device_free_after_allocation_bytes\":%llu,\"device_total_bytes\":%llu,\"twenty_percent_headroom_pass\":true,\"same_original_complete_source_tiles\":true,\"current_state_recomputed_each_stage\":true,\"source_basis_rule_thermal_and_ABI_binding\":true,\"rank0_key_hash\":%llu,\"generation_GPU_seconds_rank0\":%.9g,\"setup_wall_seconds_rank0\":%.9g,\"cache_direct_comparison_pending\":true,\"physical_trajectory_advanced\":false}\n",(unsigned long long)counts[0],(unsigned long long)(counts[1]*8),resident?"device_resident":"fixed_source_tiles_from_rank_owned_pinned_host",(unsigned long long)actual_all,(unsigned long long)free_after,(unsigned long long)total,(unsigned long long)key_hash,generation_seconds,wall_seconds);material_full::durable(f);}
        active={this,build_bridge,batch_bridge};
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
