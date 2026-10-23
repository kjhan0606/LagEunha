#!/usr/bin/env python3
"""Derive the one approved lifetime variant from the frozen reference candidate.

This script changes ownership guards and allocation order. It does not change
the numerical rule, tile size, force/EOS formula, mass, tolerances or RK4.
"""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
SRC = HERE / 'src'
changes = []


def edit(name, old, new):
    path = SRC / name
    text = path.read_text()
    assert text.count(old) == 1, (name, old)
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    path.write_text(text.replace(old, new, 1))
    changes.append({'file': name, 'before': before,
                    'after': hashlib.sha256(path.read_bytes()).hexdigest()})


def cuda_ownership():
    edit('material_full_initial_cuda.cu',
         'bool reduced_active=false;};',
         'bool reduced_active=false,reduced_reference_only=false;};')
    edit('material_full_initial_cuda.cu',
         'if(n&&uint64_t(slot)>=n)return 0;',
         'if(n&&(!v.index||!v.weight||uint64_t(slot)>=n))return 0;')
    edit('material_full_evolved_cuda.inc',
         'if(pass<0||pass>2||!count||',
         'if(!s||s->reduced_reference_only||pass<0||pass>2||!count||')
    for old, new in [
        ('if(!s||!s->generating||s->cache||!capacity)',
         'if(!s||s->reduced_reference_only||!s->generating||s->cache||!capacity)'),
        ('if(!s||!s->cache||!count||',
         'if(!s||s->reduced_reference_only||!s->cache||!count||'),
        ('if(!s||pass<1||pass>2||!s->cache||',
         'if(!s||s->reduced_reference_only||pass<1||pass>2||!s->cache||'),
    ]:
        edit('material_full_cached_cuda.inc', old, new)
    edit('material_full_reduced_rule.inc',
         'if(!s||s->reduced_offset||nt!=s->nt||',
         'if(!s||s->reduced_reference_only||s->reduced_offset||nt!=s->nt||')
    edit('material_full_reduced_rule.inc',
         'if(!s||use<0||use>1||(use&&!s->reduced_offset))',
         'if(!s||use<0||use>1||(s->reduced_reference_only&&!use)||(use&&!s->reduced_offset))')
    edit('material_reference_cache.h',
         'extern "C" void mfr_destroy(void*);',
         'extern "C" void mfr_destroy(void*);\nextern "C" int mfr_release_generation_inputs(void*,void*,uint64_t,uint64_t*);')
    edit('material_reference_cache.cu',
         'size_t start=0,count=0;uint64_t offset=0,bytes=0;',
         'size_t start=0,count=0;uint64_t offset=0,bytes=0,generated_doubles=0;')
    edit('material_reference_cache.cu',
         'int error=mfr_prefixes(s,p,start,count,offset,prefix,doubles);if(error)return error;double*payload=s->payload+(s->resident?offset:0);\n    MFR_TRY(cudaEventRecord(p->begin));mfr_generate_kernel',
         'if(!p||!s||p->reduced_reference_only||offset!=s->generated_doubles)return cudaErrorInvalidValue;\n    int error=mfr_prefixes(s,p,start,count,offset,prefix,doubles);if(error)return error;double*payload=s->payload+(s->resident?offset:0);\n    MFR_TRY(cudaEventRecord(p->begin));mfr_generate_kernel')
    edit('material_reference_cache.cu',
         '*seconds=ms/1000.;return 0;\n    #undef MFR_TRY\n}\nextern "C" int mfr_build',
         '*seconds=ms/1000.;s->generated_doubles+=doubles[count];return 0;\n    #undef MFR_TRY\n}\nextern "C" int mfr_build')
    extra = r'''
// Only the completed cache's owning context can enter reference-only mode.
// The caller additionally proves the same-stage direct witness on the host.
extern "C" int mfr_release_generation_inputs(void*ctx,void*ptr,uint64_t expected_nodes,uint64_t*freed){
    auto*p=static_cast<Context*>(ctx);auto*s=static_cast<MFRContext*>(ptr);
    if(!p||!s||s->physical!=p||!freed||!expected_nodes||p->reduced_reference_only||
       !p->reduced_active||!p->reduced_offset||!p->reduced_index||!p->reduced_weight||
       s->generated_doubles!=s->total_doubles)return cudaErrorInvalidValue;
    uint64_t nodes=0;auto e=cudaMemcpy(&nodes,p->reduced_offset+p->nt,sizeof(nodes),cudaMemcpyDeviceToHost);
    if(e!=cudaSuccess)return e;if(nodes!=expected_nodes)return cudaErrorInvalidValue;
    e=cudaDeviceSynchronize();if(e!=cudaSuccess)return e;
    // Block every original direct selected-point entry before releasing inputs.
    p->reduced_reference_only=true;
    e=cudaFree(p->reduced_index);if(e!=cudaSuccess)return e;p->reduced_index=nullptr;
    e=cudaFree(p->reduced_weight);if(e!=cudaSuccess)return e;p->reduced_weight=nullptr;
    *freed=nodes*(sizeof(uint32_t)+sizeof(double));return 0;
}
'''
    path = SRC / 'material_reference_cache.cu'
    old = path.read_text()
    edit(path.name, old, old + extra)


def host_lifetime():
    edit('material_cell_continuation_common.h',
         '#include "material_cell_condensation_solver.h"',
         '#include "material_cell_condensation_solver.h"\n#include "material_cell_device_plan.h"')
    edit('material_reference_cache_host.h',
         'double generation_seconds=0,wall_seconds=0;',
         'double generation_seconds=0,wall_seconds=0;bool inputs_released=false;')
    edit('material_reference_cache_host.h',
         'const std::string&key,const std::string&directory):material(m)',
         'const std::string&key,const std::string&directory,const cell_native::DevicePlan&solver_plan):material(m)')
    edit('material_reference_cache_host.h',
         'IMMUTABLE_REFERENCE_CACHE_V1', 'IMMUTABLE_REFERENCE_CACHE_LIFETIME_V2')
    path = SRC / 'material_reference_cache_host.h'
    text = path.read_text()
    first = text.index('        U shared=max_points*112+')
    last = text.index('        if(!resident)cell_native::gpu(cudaMallocHost', first)
    old = text[first:last]
    new = r'''        // All four ranks must be charged to the same actual device UUID.
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
'''
    edit('material_reference_cache_host.h', old, new)
    edit('material_reference_cache_host.h',
         'actual cell and reference cache keeps twenty percent device memory free',
         'actual generation phase before solver keeps twenty percent device memory free')
    # Preserve an explicit distinction between phase1 allocation and the later
    # phase2 runtime proof in the pre-existing cache summary.
    edit('material_reference_cache_host.h',
         '\\"cache_direct_comparison_pending\\":true,',
         '\\"solver_allocation_delayed\\":true,\\"generation_input_release_pending\\":true,\\"cache_direct_comparison_pending\\":true,')
    extra = r'''    void release_generation_inputs(const std::string&directory){
        valid(physical);need(!inputs_released&&material.cache_checked&&active.self==this,"complete original cache and same-stage direct witness before any input release");
        MPI_Barrier(layout.comm);size_t before,total;cell_native::gpu(cudaMemGetInfo(&before,&total),"memory before generation inputs release");
        U freed=0;cell_native::gpu(mfr_release_generation_inputs(physical,ptr,points,&freed),"enter reference-only mode and release only original selected GPU IDs/weights");
        need(freed==points*12,"exact released generation-only payload bytes");U all_freed=0;MPI_Allreduce(&freed,&all_freed,1,MPI_UINT64_T,MPI_SUM,layout.comm);
        MPI_Barrier(layout.comm);size_t after;cell_native::gpu(cudaMemGetInfo(&after,&total),"memory after every rank generation input release");
        int fits=after>=total/5,all;MPI_Allreduce(&fits,&all,1,MPI_INT,MPI_MIN,layout.comm);need(all,"actual post-release phase retains original headroom");inputs_released=true;
        if(!layout.rank){FILE*f=fopen((directory+"/generation_input_release.json").c_str(),"wx");need(f,"new completed generation/witness/release ownership evidence");
            fprintf(f,"{\"all_original_reference_points_generated\":true,\"same_actual_stage_direct_witness_pass\":true,\"freed_GPU_ID_weight_bytes_all_ranks\":%llu,\"actual_device_free_before_release_bytes\":%llu,\"actual_device_free_after_all_rank_release_bytes\":%llu,\"host_rule_and_GPU_counts_offsets_retained\":true,\"direct_selected_point_entry_disabled\":true,\"current_stage_assembled_F_not_recomputed\":true,\"saved_k1_not_recomputed\":true,\"solver_allocation_next\":true}\n",(unsigned long long)all_freed,(unsigned long long)before,(unsigned long long)after);material_full::durable(f);}
    }
'''
    edit('material_reference_cache_host.h', '    ~Cache(){', extra + '    ~Cache(){')


def driver_lifetime():
    name = 'run_cell_reference_continuation_20261022.cpp'
    edit(name, 'static std::unique_ptr<cell_native::Solver>solver;', r'''static std::unique_ptr<cell_native::Solver>solver;
static std::unique_ptr<cell_native::DevicePlan>device_plan;
static reference_cache::Cache*active_cache=nullptr;
static std::string partition_directory,factor_directory,report_directory;''')
    edit(name,
         '    need(solver&&&solver->original==&M,',
         r'''    if(!solver){
        need(active_cache&&device_plan,"one completed reference cache and exact device plan before lazy solver allocation");
        active_cache->release_generation_inputs(report_directory);
        solver=std::make_unique<cell_native::Solver>(M,partition_directory,factor_directory,report_directory);
        need(solver->planned_device_bytes==device_plan->aggregate&&solver->actual_local_device_bytes==device_plan->local,"actual solver allocation equals both-phase preflight");
    }
    need(solver&&&solver->original==&M,''')
    edit(name,
         '        solver=std::make_unique<cell_native::Solver>(*M,argv[5],argv[6],output);',
         '        partition_directory=argv[5];factor_directory=argv[6];report_directory=output;device_plan=std::make_unique<cell_native::DevicePlan>(*M,argv[5]);')
    edit(name,
         'auto*cache=new reference_cache::Cache(*m,*layout,*gpu,cfg,fe::textfile(argv[8]),output);',
         'auto*cache=new reference_cache::Cache(*m,*layout,*gpu,cfg,fe::textfile(argv[8]),output,*device_plan);active_cache=cache;')
    edit(name, 'delete cache;solver.reset();', 'delete cache;active_cache=nullptr;solver.reset();device_plan.reset();')


def main():
    assert not (HERE / 'lifetime_patch_provenance.json').exists()
    cuda_ownership()
    host_lifetime()
    driver_lifetime()
    (HERE / 'lifetime_patch_provenance.json').write_text(json.dumps({
        'scope': 'allocation order, exact byte planning and explicit reference-only ownership guards',
        'original_fixed_tile_constant': 4096,
        'original_point_capacity': 1048576,
        'physical_expressions_rule_mass_and_RK4_unchanged': True,
        'changes': changes,
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
