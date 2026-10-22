#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
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
#include "material_compact_solver.h"
namespace md=material_distributed;
namespace st=shared_trace;
namespace shared_trace_rhs {using Vec=std::vector<R>;}
#include "material_geometric_pressure_geometry.h"
namespace gp=geometric_pressure;
struct Device {void*ptr=nullptr;~Device(){mf_initial_destroy(ptr);}};
static R value(MFDoubleDouble x){return R(x.hi)+x.lo;}
static void config_rule(GfsEnrichConfig&cfg,double gamma,double pressure,double eps){
    cfg={};cfg.na=64;cfg.nr=28;cfg.field.ic.kind=2;cfg.field.ic.gamma=gamma;cfg.field.ic.pressure=pressure;cfg.field.ic.cx=cfg.field.ic.cy=12;cfg.field.ic.rin=2;cfg.field.ic.rout=8;cfg.field.ic.eps=eps;
    std::vector<double>x,w;shared_trace::gauss(64,x,w);for(int i=0;i<64;i++){cfg.an[i]=x[i];cfg.aw[i]=w[i];}shared_trace::gauss(28,x,w);for(int i=0;i<28;i++){cfg.rn[i]=x[i];cfg.rw[i]=w[i];}
}
#include "material_full_curvature_driver.h"
struct OriginalMass {
    const material_full::Owned&m;const md::Ownership&layout;
    std::vector<R>block;std::vector<U>row;std::vector<uint32_t>triangle;std::vector<uint8_t>basis;
    OriginalMass(const material_full::Owned&state,const md::Ownership&l,const std::string&path):m(state),layout(l),row(m.global.size()+1){
        FILE*f=fopen(path.c_str(),"rb");need(f,"authoritative original M");U h[4];material_full::read(f,h,4);need(h[0]==0x4d464d4153533031ULL&&h[1]==1&&h[2]==m.f.triangle.size()*ST_MAX*ST_MAX&&h[3]==sizeof(R),"original M ABI");block.resize(h[2]);material_full::read(f,block.data(),block.size());need(fgetc(f)==EOF&&fclose(f)==0,"complete original M");
        for(const auto&t:m.f.triangle)for(int a=0;a<t.n;a++)row[t.id[a]+1]++;
        for(size_t i=1;i<row.size();i++)row[i]+=row[i-1];triangle.resize(row.back());basis.resize(row.back());auto at=row;
        for(size_t s=0;s<m.f.triangle.size();s++)for(int a=0;a<m.f.triangle[s].n;a++){U k=at[m.f.triangle[s].id[a]]++;triangle[k]=s;basis[k]=a;}
    }
    md::Vec apply(const md::Vec&x)const{
        md::Vec input=x;layout.pull(input);md::Vec y(input.size());
        #pragma omp parallel for schedule(static)
        for(size_t i=0;i<y.size();i++){R sum=0;for(U k=row[i];k<row[i+1];k++){size_t s=triangle[k];int a=basis[k];const auto&t=m.f.triangle[s];const R*M=block.data()+s*ST_MAX*ST_MAX;for(int b=0;b<t.n;b++)sum+=M[a*ST_MAX+b]*input[t.id[b]];}y[i]=sum;}
        layout.sum(y);return y;
    }
};
struct Writer {
    FILE*f;U hash=14695981039346656037ULL;
    explicit Writer(const std::string&p):f(fopen(p.c_str(),"wbx")){need(f,"new compact solved output");}
    template<class T>void put(const T*p,size_t n){need(fwrite(p,sizeof(T),n,f)==n,"compact solved write");hash=compact_io::hash(p,n*sizeof(T),hash);}
    void finish(){need(fwrite(&hash,sizeof(hash),1,f)==1,"compact solved checksum");material_full::durable(f);f=nullptr;}
    ~Writer(){if(f)fclose(f);}
};
int main(int argc,char**argv){
    MPI_Init(&argc,&argv);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        need(argc==5&&ranks==4,"usage MPI4: initial_dir injection_dir mass_dir new_output_dir");double begin=MPI_Wtime();std::string suffix=std::to_string(rank),initial=argv[1],outdir=argv[4];
        int devices=0;compact_solver::gpu(cudaGetDeviceCount(&devices),"device count");need(devices>0,"native CUDA card required");compact_solver::gpu(cudaSetDevice(0),"shared one-card MPI4 mapping");
        compact_io::Injection t(std::string(argv[2])+"/compact_injection_r"+suffix+".bin",rank,ranks);
        auto m=material_full::load_basis(initial+"/full_conditioned_basis_r"+suffix+".bin");need(t.fine_global==m.global,"same original fine coefficient identities");
        // Ownership pointers are explicitly released on normal completion;
        // exceptional paths abort MPI before collective communicator cleanup.
        auto*fine_layout=new md::Ownership(m.global,m.owner,MPI_COMM_WORLD);
        auto*full_layout=new md::Ownership(t.global,t.owner,MPI_COMM_WORLD);
        std::vector<U>poly_id(t.global.begin(),t.global.end()-2);std::vector<int>poly_owner(t.owner.begin(),t.owner.end()-2);
        auto*poly_layout=new md::Ownership(poly_id,poly_owner,MPI_COMM_WORLD);need(poly_layout->global_owned==t.polynomial,"actual polynomial owned dimension");
        compact_solver::Matrix A(*poly_layout,t,std::string(argv[3])+"/compact_mass_elements_r"+suffix+".bin");
        OriginalMass original(m,*fine_layout,initial+"/full_fixed_mass_r"+suffix+".bin");
        compact_io::Reader rhs(std::string(argv[3])+"/compact_retained_rhs_r"+suffix+".bin");U rh[6];rhs.get(rh,6);need(rh[0]==0x4d43524853303031ULL&&rh[1]==1&&rh[2]==U(rank)&&rh[3]==U(ranks)&&rh[4]==t.global.size()&&rh[5]==sizeof(R),"retained weak RHS identity");
        std::array<md::Vec,2>force;for(auto&f:force){f.resize(t.global.size());rhs.get(f.data(),f.size());for(int d=0;d<2;d++)f[t.border[d]]/=A.scale[d];auto before=f;full_layout->pull(f);need(f==before,"normalized retained RHS owner/ghost match");}rhs.finish();
        auto exact=[&](const md::Vec&normalized){
            auto raw=normalized;for(int d=0;d<2;d++)raw[t.border[d]]/=A.scale[d];full_layout->pull(raw);
            auto fine=t.forward(raw);auto load=original.apply(fine);auto compact=t.restrict_assembled_owned_load(load,fine_layout->owned);full_layout->sum(compact);
            for(int d=0;d<2;d++)compact[t.border[d]]/=A.scale[d];return compact;
        };
        std::array<md::Vec,2>solution;R true_residual[2]={};U outer_count[2]={},iterations=0,polynomial_corrections=0,ghost_bytes=0;R schur[4]={};double inner_seconds=0,halo_seconds=0,schur_seconds=0,solve_seconds=0;
        {
            compact_solver::Resident resident(A);double solve_begin=MPI_Wtime();compact_solver::Bordered border(A,resident);std::copy(border.S,border.S+4,schur);schur_seconds=border.setup_seconds;
            for(int component=0;component<2;component++){
                auto&x=solution[component];x.resize(force[component].size());auto r=force[component];R norm=full_layout->dot(r,r);need(norm>0,"nonzero retained projected physical load");
                for(int correction=0;correction<3;correction++){
                    auto dx=border.solve(r,t);for(size_t i=0;i<x.size();i++)x[i]+=dx[i];auto ax=exact(x);for(size_t i=0;i<r.size();i++)r[i]=force[component][i]-ax[i];outer_count[component]++;
                    true_residual[component]=sqrtl(full_layout->dot(r,r)/norm);if(true_residual[component]<1e-16L)break;
                }
                need(true_residual[component]<1e-16L,"actual T^T original binary80 M T equation residual1e-16");
            }
            solve_seconds=MPI_Wtime()-solve_begin;iterations=resident.iterations;polynomial_corrections=resident.corrections;inner_seconds=resident.inner_seconds;halo_seconds=resident.halo_seconds;ghost_bytes=resident.ghost_bytes;
        }
        Writer solved(outdir+"/compact_initial_acceleration_r"+suffix+".bin");U sh[8]={0x4d43414343454c31ULL,1,U(rank),U(ranks),U(t.global.size()),t.total,sizeof(R),1};solved.put(sh,8);solved.put(A.scale,2);solved.put(t.global.data(),t.global.size());for(auto&x:solution)solved.put(x.data(),x.size());solved.finish();
        FILE*af=fopen((outdir+"/full_initial_acceleration_r"+suffix+".bin").c_str(),"wbx");need(af,"new lifted fine acceleration");U ah[4]={0x4d46414343454c31ULL,1,U(m.global.size()),sizeof(R)};material_full::write(af,ah,4);
        for(auto&x:solution){auto raw=x;for(int d=0;d<2;d++)raw[t.border[d]]/=A.scale[d];full_layout->pull(raw);auto lifted=t.forward(raw);material_full::write(af,lifted.data(),lifted.size());}material_full::durable(af);
        R position_error=0,velocity_error=0;
        for(int d=0;d<4;d++){
            md::Vec c(t.global.size());for(size_t i=0;i<c.size();i++)c[i]=d<2?t.initial[i].x[d]:t.initial[i].v[d-2];full_layout->pull(c);auto field=t.forward(c);
            for(size_t i=0;i<field.size();i++)if(d<2){position_error=std::max(position_error,fabsl(R(double(field[i]))-m.f.coefficient[i].x[d]));m.f.coefficient[i].x[d]=double(field[i]);}else{velocity_error=std::max(velocity_error,fabsl(R(double(field[i]))-m.f.coefficient[i].v[d-2]));m.f.coefficient[i].v[d-2]=double(field[i]);}
        }
        // This is an initial-state representation for the compact branch,
        // persisted in the existing fine reader format for the unchanged
        // curvature kernel. It is not an accepted time-evolution checkpoint.
        material_full::save_basis(m,outdir+"/full_conditioned_basis_r"+suffix+".bin");
        R errors[2]={position_error,velocity_error};MPI_Allreduce(MPI_IN_PLACE,errors,2,MPI_LONG_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        U nonzeros=A.value.size();MPI_Allreduce(MPI_IN_PLACE,&nonzeros,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
        if(rank==0){
            FILE*f=fopen((outdir+"/compact_initial_solve_result.json").c_str(),"wx");need(f,"new actual compact CUDA solve report");
            fprintf(f,"{\"scope\":\"CUDA resident DD assembled polynomial solve plus rank2 Schur and original binary80 T^T M T outer residual; no evolution\",\"global_scalar_dofs\":%llu,\"global_partial_CSR_entries\":%llu,\"true_equation_residual\":[%.21Lg,%.21Lg],\"full_outer_corrections\":[%llu,%llu],\"polynomial_inner_iterations_total\":%llu,\"polynomial_corrections_total\":%llu,\"computed_Schur\":[%.21Lg,%.21Lg,%.21Lg,%.21Lg],\"border_norms\":[%.21Lg,%.21Lg],\"inner_vectors_resident_on_GPU\":true,\"ghost_DD_pairs_only_during_inner_iterations\":true,\"rank0_ghost_pair_bytes\":%llu,\"rank0_inner_seconds\":%.9g,\"rank0_halo_seconds\":%.9g,\"rank0_Schur_setup_seconds\":%.9g,\"rank0_solve_including_original_true_residual_seconds\":%.9g,\"initial_position_coefficient_difference\":%.21Lg,\"initial_velocity_coefficient_difference\":%.21Lg,\"initial_pressure_curvature_pending\":true,\"retained_force_inherits_original_solve_residual\":true,\"seconds_before_curvature\":%.9g,\"accepted_steps\":0,\"long_Kepler_validated\":false}\n",(unsigned long long)t.total,(unsigned long long)nonzeros,true_residual[0],true_residual[1],(unsigned long long)outer_count[0],(unsigned long long)outer_count[1],(unsigned long long)iterations,(unsigned long long)polynomial_corrections,schur[0],schur[1],schur[2],schur[3],A.scale[0],A.scale[1],(unsigned long long)ghost_bytes,inner_seconds,halo_seconds,schur_seconds,solve_seconds,errors[0],errors[1],MPI_Wtime()-begin);material_full::durable(f);
        }
        delete poly_layout;delete full_layout;delete fine_layout;
        // Read actual gamma/P0/eps from the immutable original thermal input.
        FILE*thermal=fopen((initial+"/full_generating_thermal_r"+suffix+".bin").c_str(),"rb");need(thermal,"original generating thermal law");U th[6];double law[4];material_full::read(thermal,th,6);material_full::read(thermal,law,4);need(th[0]==0x4d46544845524d31ULL&&th[1]==1&&th[3]==sizeof(R)&&th[4]==U(rank),"immutable thermal header");fclose(thermal);
        need(chdir(outdir.c_str())==0,"new curvature output working directory");
        int result=curvature_from_saved(outdir,MPI_COMM_WORLD,law[0],law[1],law[2]);
        MPI_Finalize();return result;
    }catch(const std::exception&e){fprintf(stderr,"rank%d COMPACT_INITIAL_REJECTED %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);return 2;}
}
