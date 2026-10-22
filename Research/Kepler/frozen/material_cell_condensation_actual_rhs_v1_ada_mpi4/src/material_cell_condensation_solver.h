#ifndef MATERIAL_CELL_CONDENSATION_SOLVER_H
#define MATERIAL_CELL_CONDENSATION_SOLVER_H
#include <memory>
#include <cuda_runtime.h>
#include "material_compact_gpu.h"
#include "material_cell_condensation_gpu.h"
#include "material_cell_factors_io.h"
#include "material_distributed_mass.h"
namespace cell_native {
using R=long double;using U=uint64_t;using Vec=std::vector<R>;using DD=MFDoubleDouble;
namespace md=material_distributed;namespace cc=cell_condensation;
inline DD encode(R x){double hi=double(x);return {hi,double(x-R(hi))};}
inline DD pair_add(DD a,DD b){double s=a.hi+b.hi,z=s-a.hi,e=(a.hi-(s-z))+(b.hi-z);e+=a.lo+b.lo;double h=s+e;return {h,e-(h-s)};}
inline void gpu(int code,const char*s){if(code){fprintf(stderr,"CELL_CONDENSATION_CUDA_REJECTED %s code%d\n",s,code);MPI_Abort(MPI_COMM_WORLD,2);}}
struct Matrix {
    const md::Ownership&layout;std::vector<U>row;std::vector<uint32_t>column;Vec value,inverse_diagonal;
    Matrix(const md::Ownership&l,const cc::Partition&p,const std::string&path):layout(l){
        cc::SkeletonMass saved(p,path,l.rank,l.ranks);row=std::move(saved.row);column=std::move(saved.column);value=std::move(saved.value);inverse_diagonal=std::move(saved.diagonal);
        for(auto&v:inverse_diagonal){need(v>0,"positive original skeleton diagonal");v=1/v;}
    }
    Vec apply(const Vec&x)const{
        Vec input=x;layout.pull(input);Vec y(x.size());
        #pragma omp parallel for schedule(static)
        for(size_t i=0;i<y.size();i++){R sum=0;for(U j=row[i];j<row[i+1];j++)sum+=value[j]*input[column[j]];y[i]=sum;}
        layout.sum(y);return y;
    }
    U device_bytes()const{
        U ns=0,na=0;for(int p=0;p<layout.ranks;p++){ns+=layout.serve[p].size();na+=layout.ask[p].size();}
        U n=layout.global.size(),no=layout.owned.size();
        return (n+1)*8+column.size()*20+n*16+no*8+(ns+na)*8+5*n*16+std::max<U>(1,std::max(ns,na))*16+((no+255)/256)*16+16;
    }
};
#include "material_cell_skeleton_resident.h"
struct Solution {
    std::array<md::Solution,2>component;U full_corrections=0,skeleton_iterations=0,skeleton_corrections=0;
    double seconds=0,condense_seconds=0,recover_seconds=0,true_apply_seconds=0,halo_seconds=0,inner_seconds=0;
};
struct Solver {
    const md::Mass&original;cc::Partition partition;md::Ownership layout;Matrix matrix;
    void*ptr=nullptr;std::unique_ptr<Resident>resident;
    std::array<std::vector<DD>,2>input,output,reduced;
    U planned_device_bytes=0,actual_local_device_bytes=0;double setup_seconds=0;
    Solver(const md::Mass&M,const std::string&partition_dir,const std::string&factor_dir,const std::string&report_dir):
        original(M),partition(partition_dir+"/cell_partition_r"+std::to_string(M.layout.rank)+".bin",M.layout.rank,M.layout.ranks),
        layout(partition.global_G,partition.owner_G,M.layout.comm),
        matrix(layout,partition,factor_dir+"/skeleton_mass_r"+std::to_string(M.layout.rank)+".bin"){
        double start=MPI_Wtime();need(partition.nf==M.layout.global.size()&&layout.global_owned==partition.total,"same full original/skeleton owner extent");
        for(size_t i=0;i<partition.ng;i++)need(partition.global_G[i]==M.layout.global[partition.old_G[i]]&&partition.owner_G[i]==M.layout.owner[partition.old_G[i]],"unchanged original trace provenance");
        cc::Factors factors(partition,factor_dir+"/cell_factors_r"+std::to_string(layout.rank)+".bin",layout.rank,layout.ranks);
        std::vector<MCCell>cell;std::vector<uint32_t>I,G;std::vector<U>row(partition.ng+1),entry;std::vector<uint8_t>owned(partition.ng);
        for(size_t c=0;c<partition.cell.size();c++){
            const auto&a=partition.cell[c];cell.push_back({U(I.size()),U(G.size()),factors.L_offset[c],factors.IG_offset[c],uint32_t(a.I.size()),uint32_t(a.G.size())});
            I.insert(I.end(),a.I.begin(),a.I.end());G.insert(G.end(),a.G.begin(),a.G.end());for(auto g:a.G)row[g+1]++;
        }
        for(size_t i=1;i<row.size();i++)row[i]+=row[i-1];entry.resize(G.size());auto at=row;
        for(size_t k=0;k<G.size();k++)entry[at[G[k]]++]=k;
        for(auto i:layout.owned)owned[i]=1;
        U local=cell.size()*sizeof(MCCell)+I.size()*4+G.size()*4+partition.ng*4+partition.ng+row.size()*8+entry.size()*8+
            (factors.L.size()+factors.IG.size())*16+4*partition.nf*16+2*G.size()*16+2*partition.ng*16+matrix.device_bytes();
        MPI_Allreduce(&local,&planned_device_bytes,1,MPI_UINT64_T,MPI_SUM,layout.comm);
        size_t free_before=0,total=0;gpu(cudaMemGetInfo(&free_before,&total),"fresh cell solver free memory");int capacity=free_before>=total/5&&planned_device_bytes<=free_before-total/5,all;
        MPI_Allreduce(&capacity,&all,1,MPI_INT,MPI_MIN,layout.comm);need(all,"all-ranks current device memory budget with twenty percent headroom");
        std::vector<DD>L(factors.L.size()),IG(factors.IG.size());
        #pragma omp parallel for schedule(static)
        for(size_t i=0;i<L.size();i++)L[i]=encode(factors.L[i]);
        #pragma omp parallel for schedule(static)
        for(size_t i=0;i<IG.size();i++)IG[i]=encode(factors.IG[i]);
        gpu(mcc_create(partition.nf,partition.ng,cell.data(),cell.size(),I.data(),I.size(),G.data(),G.size(),partition.old_G.data(),owned.data(),row.data(),entry.data(),L.data(),L.size(),IG.data(),IG.size(),&ptr),"all original cell factors and recovery buffers");
        resident=std::make_unique<Resident>(matrix);actual_local_device_bytes=mcc_device_bytes(ptr)+matrix.device_bytes();need(actual_local_device_bytes==local,"exact explicit CUDA allocation accounting");
        MPI_Barrier(layout.comm);size_t free_after=0;gpu(cudaMemGetInfo(&free_after,&total),"actual post-allocation device memory");capacity=free_after>=total/5;MPI_Allreduce(&capacity,&all,1,MPI_INT,MPI_MIN,layout.comm);need(all,"actual full cell/skeleton GPU allocation retains twenty percent total headroom");
        for(int d=0;d<2;d++){input[d].resize(partition.nf);output[d].resize(partition.nf);reduced[d].resize(partition.ng);}
        setup_seconds=MPI_Wtime()-start;
        if(!layout.rank){FILE*f=fopen((report_dir+"/cell_solver_allocation.json").c_str(),"wx");need(f,"new actual cell solver allocation evidence");fprintf(f,"{\"planned_aggregate_new_device_bytes\":%llu,\"device_total_bytes\":%llu,\"device_free_before_new_solver_bytes\":%llu,\"device_free_after_all_ranks_bytes\":%llu,\"physical_context_already_allocated\":true,\"all_original_modes_retained\":true,\"twenty_percent_headroom_pass\":true,\"setup_seconds_rank0\":%.9g}\n",(unsigned long long)planned_device_bytes,(unsigned long long)total,(unsigned long long)free_before,(unsigned long long)free_after,setup_seconds);material_full::durable(f);}
    }
    ~Solver(){resident.reset();mcc_destroy(ptr);}
    Solver(const Solver&)=delete;Solver&operator=(const Solver&)=delete;
    Solution solve(const std::array<Vec,2>&force){
        double start=MPI_Wtime(),halo0=resident->halo_seconds,inner0=resident->inner_seconds;U it0=resident->iterations,correction0=resident->corrections;Solution result;
        std::array<Vec,2>r=force;R norm[2],previous[2]={1,1};bool pass[2]={};
        for(int d=0;d<2;d++){need(force[d].size()==partition.nf,"full force extent");norm[d]=original.layout.dot(force[d],force[d]);need(norm[d]>0&&std::isfinite(norm[d]),"positive actual original force norm");result.component[d].x.resize(partition.nf);}
        for(int correction=0;correction<3;correction++){
            double stamp=MPI_Wtime();
            #pragma omp parallel for schedule(static)
            for(size_t i=0;i<partition.nf;i++)for(int d=0;d<2;d++)input[d][i]=encode(pass[d]?0:r[d][i]);
            gpu(mcc_rhs(ptr,input[0].data(),input[1].data(),reduced[0].data(),reduced[1].data()),"actual original source RHS condensation");
            std::array<Vec,2>b={Vec(partition.ng),Vec(partition.ng)};
            for(int d=0;d<2;d++){
                #pragma omp parallel for schedule(static)
                for(size_t i=0;i<partition.ng;i++)b[d][i]=R(reduced[d][i].hi)+reduced[d][i].lo;
                layout.sum(b[d]);
            }
            result.condense_seconds+=MPI_Wtime()-stamp;
            for(int d=0;d<2;d++){
                auto g=resident->solve(b[d]);
                #pragma omp parallel for schedule(static)
                for(size_t i=0;i<partition.ng;i++)reduced[d][i]=encode(g[i]);
            }
            stamp=MPI_Wtime();gpu(mcc_recover(ptr,reduced[0].data(),reduced[1].data(),output[0].data(),output[1].data()),"all original interior accelerations recovered");
            #pragma omp parallel for schedule(static)
            for(size_t i=0;i<partition.nf;i++)for(int d=0;d<2;d++)result.component[d].x[i]+=R(output[d][i].hi)+output[d][i].lo;
            result.recover_seconds+=MPI_Wtime()-stamp;stamp=MPI_Wtime();
            for(int d=0;d<2;d++){
                auto ax=original.apply(result.component[d].x);
                #pragma omp parallel for schedule(static)
                for(size_t i=0;i<partition.nf;i++)r[d][i]=force[d][i]-ax[i];
                R error=sqrtl(original.layout.dot(r[d],r[d])/norm[d]);need(std::isfinite(error),"finite original full Ma-f residual");result.component[d].true_relative_residual=error;
                pass[d]=error<1e-16L;need(pass[d]||error<previous[d],"strict decrease of original full equation residual");previous[d]=error;
            }
            result.true_apply_seconds+=MPI_Wtime()-stamp;result.full_corrections++;if(pass[0]&&pass[1])break;
        }
        need(pass[0]&&pass[1],"three exact-cell corrections satisfy original full binary80 Ma-f below1e-16");
        result.skeleton_iterations=resident->iterations-it0;result.skeleton_corrections=resident->corrections-correction0;
        result.halo_seconds=resident->halo_seconds-halo0;result.inner_seconds=resident->inner_seconds-inner0;result.seconds=MPI_Wtime()-start;
        for(auto&a:result.component)a.iterations=result.skeleton_iterations;return result;
    }
};
using Hook=std::array<md::Solution,2>(*)(const md::Mass&,const std::array<Vec,2>&,const std::string&,R,R,const std::vector<STCoefficient>&);
inline Hook native_rhs_hook=nullptr;
}
#endif
