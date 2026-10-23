// First actual selected-rule RHS at retained step5; never a fresh IC.
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>
#include <cuda_runtime.h>
static void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
using R=long double;using U=uint64_t;
#include "material_full_owned.h"
#include "material_full_initial_batch.h"
#include "material_distributed_mass.h"
namespace st=shared_trace;namespace md=material_distributed;
namespace shared_trace_rhs {using Vec=std::vector<R>;}
#include "material_geometric_pressure_geometry.h"
namespace gp=geometric_pressure;
static R value(MFDoubleDouble a){return R(a.hi)+a.lo;}
struct Device {void*ptr=nullptr;~Device(){mf_initial_destroy(ptr);}};
static void config_rule(GfsEnrichConfig&cfg,double gamma,double pressure,double eps){
    cfg={};cfg.na=64;cfg.nr=28;cfg.field.ic.kind=2;cfg.field.ic.gamma=gamma;cfg.field.ic.pressure=pressure;cfg.field.ic.cx=cfg.field.ic.cy=12;cfg.field.ic.rin=2;cfg.field.ic.rout=8;cfg.field.ic.eps=eps;
    std::vector<double>x,w;st::gauss(64,x,w);for(int i=0;i<64;i++){cfg.an[i]=x[i];cfg.aw[i]=w[i];}st::gauss(28,x,w);for(int i=0;i<28;i++){cfg.rn[i]=x[i];cfg.rw[i]=w[i];}
}
#include "material_cell_condensation_solver.h"
// This isolated frozen header has one explicit solve hook; all original
// physical force, energy/torque, per-source spin and reaction gates are retained.
#include "material_full_evolution_host.h"
namespace fe=full_evolution;
struct Writer {
    FILE*f;U checksum=14695981039346656037ULL;
    explicit Writer(const std::string&path):f(fopen(path.c_str(),"wbx")){need(f,"new actual RHS artifact");}
    template<class T>void put(const T*p,size_t n){need(fwrite(p,sizeof(T),n,f)==n,"actual RHS artifact write");checksum=fe::hash(p,n*sizeof(T),checksum);}
    void finish(){need(fwrite(&checksum,sizeof(U),1,f)==1,"actual RHS artifact checksum");material_full::durable(f);f=nullptr;}
    ~Writer(){if(f)fclose(f);}
};
static std::string partition_dir,factor_dir;
static std::unique_ptr<cell_native::Solver>selected;
static U state_hash=0,immutable_hash=0,rule_generation=0;static double accepted_time=0,next_h=0;
static std::array<md::Solution,2>actual_solve(const md::Mass&M,const std::array<md::Vec,2>&F,const std::string&path,R Uheat,R kinetic,const std::vector<STCoefficient>&coefficient){
    const auto&layout=M.layout;std::string suffix=std::to_string(layout.rank);double begin=MPI_Wtime();
    need(!selected,"this bounded executable consumes exactly one genuinely new RHS");
    Writer saved(path+"/actual_force_r"+suffix+".bin");U h[10]={0x4d43414354463031ULL,1,U(layout.rank),U(layout.ranks),U(F[0].size()),sizeof(R),state_hash,immutable_hash,rule_generation,5};saved.put(h,10);saved.put(&accepted_time,1);saved.put(&next_h,1);for(auto&f:F)saved.put(f.data(),f.size());saved.finish();
    std::array<full_gpu_mass::Solution,2>reference;double shadow_start=MPI_Wtime(),old_upload=0;
    {
        full_gpu_mass::Operator old(M);old_upload=old.upload_seconds;
        for(int d=0;d<2;d++)reference[d]=full_gpu_mass::solve(old,F[d]);
        need(old.resident,"old full native solve actually uses resident CUDA vectors and DD MPI pairs");
    } // release the old full mass GPU buffers before allocating the cell solver.
    double shadow_seconds=MPI_Wtime()-shadow_start;
    Writer old(path+"/same_force_old_acceleration_r"+suffix+".bin");old.put(h,10);for(auto&a:reference)old.put(a.x.data(),a.x.size());old.finish();
    if(!layout.rank){FILE*f=fopen((path+"/old_full_shadow.json").c_str(),"wx");need(f,"actual original full shadow evidence");fprintf(f,"{\"same_actual_force\":true,\"resident_CUDA_vectors_and_DD_MPI\":true,\"old_full_true_residual\":[%.21Lg,%.21Lg],\"upload_seconds_rank0\":%.9g,\"total_shadow_seconds_rank0\":%.9g,\"old_GPU_buffers_released_before_new_solver\":true}\n",reference[0].true_relative_residual,reference[1].true_relative_residual,old_upload,shadow_seconds);material_full::durable(f);}
    selected=std::make_unique<cell_native::Solver>(M,partition_dir,factor_dir,path);
    auto result=selected->solve(F);R work=0,torque=0,mass_error=0,maximum=0,scale=0,internal=0,load=0;
    for(int d=0;d<2;d++){
        md::Vec diff(F[d].size()),v(diff.size()),x(diff.size());
        for(size_t i=0;i<diff.size();i++){diff[i]=result.component[d].x[i]-reference[d].x[i];v[i]=coefficient[i].v[d];x[i]=d==1?coefficient[i].x[0]:-coefficient[i].x[1];}
        auto mdiff=M.apply(diff),ma=M.apply(result.component[d].x),mref=M.apply(reference[d].x);
        work+=layout.dot(v,mdiff);torque+=layout.dot(x,mdiff);
        R numerator=layout.dot(diff,mdiff),denominator=layout.dot(reference[d].x,mref);need(numerator>=0&&denominator>0,"positive original mass norm of actual acceleration difference");mass_error=std::max(mass_error,sqrtl(numerator/denominator));
        for(auto i:layout.owned){maximum=std::max(maximum,fabsl(diff[i]));scale=std::max(scale,fabsl(reference[d].x[i]));load=std::max(load,fabsl(F[d][i]));}
        for(const auto&cell:selected->partition.cell)for(auto i:cell.I)internal=std::max(internal,fabsl(ma[i]-F[d][i]));
    }
    R highs[5]={mass_error,maximum,scale,internal,load};MPI_Allreduce(MPI_IN_PLACE,highs,5,MPI_LONG_DOUBLE,MPI_MAX,layout.comm);
    bool pass=fabsl(work)/Uheat<1e-8L&&fabsl(torque)/(2*kinetic)<1e-11L&&highs[0]<1e-11L&&highs[4]>0&&highs[3]/highs[4]<1e-11L;
    if(!layout.rank){FILE*f=fopen((path+"/cell_condensation_actual_solve.json").c_str(),"wx");need(f,"new actual full equation and old-solve comparison");fprintf(f,"{\"actual_force_at_retained_step5\":true,\"time\":%.17g,\"same_selected_rule_generation\":%llu,\"full_original_M_true_residual\":[%.21Lg,%.21Lg],\"full_outer_corrections\":%llu,\"skeleton_inner_iterations_total\":%llu,\"skeleton_corrections_total\":%llu,\"same_force_old_acceleration_M_norm_relative\":%.21Lg,\"same_force_old_acceleration_max_difference\":%.21Lg,\"same_force_old_acceleration_max_scale\":%.21Lg,\"same_force_work_difference_over_U\":%.21Lg,\"same_force_torque_difference_over_2K\":%.21Lg,\"all_cell_interior_reaction_relative\":%.21Lg,\"new_primary_solve_seconds_rank0\":%.9g,\"new_condense_seconds_rank0\":%.9g,\"new_recover_seconds_rank0\":%.9g,\"original_true_apply_seconds_rank0\":%.9g,\"skeleton_inner_seconds_rank0\":%.9g,\"skeleton_halo_seconds_rank0\":%.9g,\"old_full_shadow_seconds_rank0\":%.9g,\"new_setup_seconds_rank0\":%.9g,\"total_hook_seconds_rank0\":%.9g,\"same_actual_RHS_comparison_pass\":%s,\"physical_time_advanced\":false,\"not_a_complete_RK4_cost\":true}\n",accepted_time,(unsigned long long)rule_generation,result.component[0].true_relative_residual,result.component[1].true_relative_residual,(unsigned long long)result.full_corrections,(unsigned long long)result.skeleton_iterations,(unsigned long long)result.skeleton_corrections,highs[0],highs[1],highs[2],fabsl(work)/Uheat,fabsl(torque)/(2*kinetic),highs[3]/highs[4],result.seconds,result.condense_seconds,result.recover_seconds,result.true_apply_seconds,result.inner_seconds,result.halo_seconds,shadow_seconds,selected->setup_seconds,MPI_Wtime()-begin,pass?"true":"false");material_full::durable(f);}
    need(pass,"actual full-space condensation acceleration/work/torque/interior-reaction comparison");return std::move(result.component);
}
int main(int argc,char**argv){
    MPI_Init(&argc,&argv);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        need(argc==8&&ranks==4,"MPI4 usage: actual_rhs original_initial positive_generation original_immutable_manifest step5_prefix partition_run factor_run new_output");double begin=MPI_Wtime();
        int devices=0;cell_native::gpu(cudaGetDeviceCount(&devices),"visible device count");need(devices==1,"one shared visible card for original MPI4");cell_native::gpu(cudaSetDevice(0),"original one-card mapping");cudaDeviceProp prop{};cell_native::gpu(cudaGetDeviceProperties(&prop,0),"actual device properties");
        std::string suffix=std::to_string(rank),output=argv[7];partition_dir=argv[5];factor_dir=argv[6];
        auto*m=new fe::Material;m->saved=argv[1];m->immutable=fe::textfile(argv[3]);m->owned=material_full::load_basis(m->saved+"/full_conditioned_basis_r"+suffix+".bin");
        auto*layout=new md::Ownership(m->owned.global,m->owned.owner,MPI_COMM_WORLD);
        FILE*f=fopen((m->saved+"/full_generating_thermal_r"+suffix+".bin").c_str(),"rb");need(f,"original thermal law header");U th[6];double law[4];material_full::read(f,th,6);material_full::read(f,law,4);need(th[0]==0x4d46544845524d31ULL&&th[1]==1&&th[3]==sizeof(R)&&th[4]==U(rank)&&fclose(f)==0,"retained thermal law identity");m->gamma=law[0];m->P0=law[1];m->eps=law[2];fe::thermal_import(*m,rank);
        fe::load_positive_rule(*m,*layout,argv[2],false);fe::restore_previous(*m,*layout,argv[4]);
        need(m->ledger.steps==5&&m->ledger.stage==0&&m->ledger.time==0.00011263377801805954&&m->numerical.reduced==1&&m->positive.active,"exact retained step5 selected trajectory; no new IC or accepted replay");
        accepted_time=m->ledger.time;next_h=m->ledger.h;rule_generation=m->positive.generation;state_hash=fe::hash(m->owned.f.coefficient.data(),m->owned.f.coefficient.size()*sizeof(STCoefficient));immutable_hash=fe::hash(m->immutable.data(),m->immutable.size());
        f=fopen((m->saved+"/full_fixed_mass_r"+suffix+".bin").c_str(),"rb");need(f,"retained original M");U mh[4];material_full::read(f,mh,4);need(mh[0]==0x4d464d4153533031ULL&&mh[1]==1&&mh[2]==m->owned.f.triangle.size()*ST_MAX*ST_MAX&&mh[3]==sizeof(R),"original mass ABI");md::Vec blocks(mh[2]);material_full::read(f,blocks.data(),blocks.size());need(fgetc(f)==EOF&&fclose(f)==0,"complete retained original M");auto*M=new md::Mass(*layout,m->owned.f.triangle,std::move(blocks));
        GfsEnrichConfig cfg;config_rule(cfg,m->gamma,m->P0,m->eps);auto*gpu=new Device;
        cell_native::gpu(mf_initial_create(m->owned.f.triangle.data(),m->owned.f.triangle.size(),m->generating.data(),m->generating.size(),m->owned.rule.data(),&cfg,&gpu->ptr),"same original persistent physical CUDA context");
        cell_native::gpu(mf_evolved_bind(gpu->ptr,m->generating.data(),m->generating.size(),m->owned.f.cell.size()),"same original generating field binding");fe::upload_positive_rule(*m,*layout,*gpu);fe::use_positive_rule(*m,*layout,*gpu,true);
        m->cached_points=false; // immutable cache is pending; same direct selected physics, not a cost promotion.
        m->quadrature_oracle=true; // retain same-stage source metadata; call rhs directly, no high-rule replay.
        fe::mkdir_collective(output+"/stage1",*layout);m->ledger.stage=1;cell_native::native_rhs_hook=actual_solve;fe::Endpoint endpoint;R dt=1e99L;
        auto k1=fe::rhs(*m,*M,*gpu,cfg,false,output+"/stage1",endpoint,dt,nullptr,false);cell_native::native_rhs_hook=nullptr;
        need(fe::hash(m->owned.f.coefficient.data(),m->owned.f.coefficient.size()*sizeof(STCoefficient))==state_hash&&m->ledger.time==accepted_time&&m->ledger.steps==5,"actual RHS does not mutate accepted physical state");
        Writer k(output+"/prepared_selected_k1_r"+suffix+".bin");U kh[11]={0x4d434e4558544b31ULL,1,U(rank),U(ranks),U(k1.k.size()),sizeof(STCoefficient),sizeof(R),state_hash,immutable_hash,rule_generation,m->ledger.epoch};k.put(kh,11);k.put(&accepted_time,1);k.put(&next_h,1);R rates[4]={k1.work,k1.torque,k1.heat,k1.dt_bound};k.put(rates,4);k.put(k1.k.data(),k1.k.size());k.finish();m->ledger.stage=0;
        if(!rank){f=fopen((output+"/actual_rhs_result.json").c_str(),"wx");need(f,"new actual next RHS completion");fprintf(f,"{\"time\":%.17g,\"accepted_steps\":5,\"new_accepted_steps\":0,\"new_actual_RHS_evaluations\":1,\"same_original_selected_physics\":true,\"selected_points\":%llu,\"same_force_old_resident_full_solve_shadow\":true,\"cell_condensation_full_original_residual_and_physical_cell_gates_pass\":true,\"prepared_selected_k1_bound_to_exact_step5\":true,\"stage_cache_mode\":\"direct selected rule, immutable cache pending\",\"device\":\"%s\",\"MPI_ranks\":4,\"OMP_threads\":%d,\"full_step_cost_pending\":true,\"cache_direct_comparison_pending\":true,\"global_map_guard_pending\":true,\"long_Kepler_validated\":false,\"seconds\":%.9g}\n",accepted_time,(unsigned long long)m->positive.global_points,prop.name,omp_get_max_threads(),MPI_Wtime()-begin);material_full::durable(f);}
        selected.reset();delete gpu;delete M;delete layout;delete m;MPI_Finalize();return 0;
    }catch(const std::exception&e){fprintf(stderr,"rank%d ACTUAL_NEXT_RHS_CONDENSATION_REJECTED %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);return 2;}
}
