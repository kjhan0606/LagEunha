// Two-rank linear-algebra fixture. NO saved source, remap or physical step.
#include "material_compact_mass.h"
#include <cuda_runtime.h>
#include <iostream>
#include <iomanip>
using R=long double;
namespace md=material_distributed;
static void need(bool b,const char*s){md::require(b,s);}
#include "material_full_gpu_mass_host.h"
int main(int argc,char**argv){
    MPI_Init(&argc,&argv);int rank=0,nranks=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&nranks);
    try{need(nranks==2,"fixture requires exactly two MPI ranks");{
        need(cudaSetDevice(0)==cudaSuccess,"Ada fixture GPU mapping");
        std::vector<uint64_t>global;std::vector<int>owner;
        for(int j=0;j<10;j++){global.push_back(10*rank+j);owner.push_back(rank);}
        global.push_back(20);owner.push_back(0);global.push_back(21);owner.push_back(1);
        md::Ownership layout(global,owner,MPI_COMM_WORLD);
        compact_mass::LocalElement e;e.n=12;e.component=0;for(int j=0;j<12;j++)e.id[j]=j;
        for(int j=0;j<10;j++)e.mass[j*12+j]=1;
        e.mass[(10+rank)*12+10+rank]=1;e.mass[10+rank]=e.mass[(10+rank)*12]=.05L;
        std::array<R,ST_MAX*ST_MAX>old_block={};
        for(int a=0;a<12;a++)for(int b=0;b<12;b++)old_block[a*ST_MAX+b]=e.mass[a*12+b];
        int rejected=0;try{md::Factor unsuitable(old_block.data(),12);}catch(const std::exception&){rejected=1;}
        int global_rejected=0;MPI_Allreduce(&rejected,&global_rejected,1,MPI_INT,MPI_SUM,layout.comm);
        need(global_rejected==2,"both local physical blocks are deliberately rank deficient");
        auto prepared=compact_mass::prepare(layout,{e},{{0,{20,21}}});auto&M=*prepared->exact;
        auto unchanged=[&](){bool same=M.block.size()==old_block.size();for(size_t j=0;j<old_block.size();j++)same&=M.block[j]==old_block[j];return same;};
        need(unchanged(),"surrogate does not alter physical M");
        md::Vec expected(layout.global.size());for(size_t j=0;j<expected.size();j++)expected[j]=R(layout.global[j]+1)/17;
        auto rhs=M.apply(expected);full_gpu_mass::Operator op(M);need(op.conversion_error==0,"original DD upload exactly restores native inputs");
        auto gpu_product=op.apply(expected);R product_error=0;
        for(size_t j:layout.owned)product_error=std::max(product_error,fabsl(gpu_product[j]-rhs[j]));
        MPI_Allreduce(MPI_IN_PLACE,&product_error,1,MPI_LONG_DOUBLE,MPI_MAX,layout.comm);
        auto answer=full_gpu_mass::solve(op,rhs);R solution_error=0;
        for(size_t j:layout.owned)solution_error=std::max(solution_error,fabsl(answer.x[j]-expected[j]));
        MPI_Allreduce(MPI_IN_PLACE,&solution_error,1,MPI_LONG_DOUBLE,MPI_MAX,layout.comm);
        need(product_error<1e-16L&&solution_error<1e-14L&&answer.true_relative_residual<1e-16L&&op.resident,"original CUDA resident/MPI solve closure");
        need(unchanged(),"physical M preserved after original resident solve");
        uint64_t bytes=op.resident->ghost_pair_bytes;MPI_Allreduce(MPI_IN_PLACE,&bytes,1,MPI_UINT64_T,MPI_SUM,layout.comm);
        need(bytes>0,"actual sparse MPI ghost exchange exercised");
        if(rank==0)std::cout<<std::setprecision(17)<<"{\"status\":\"COMPACT_SEMIDEFINITE_LOCAL_MASS_ORIGINAL_RESIDENT_MPI2_ADA_PASS_NO_PHYSICAL_RUN\",\"ranks\":2,\"global_owned_dofs\":"<<layout.global_owned<<",\"rank_deficient_local_blocks\":"<<global_rejected<<",\"physical_mass_unchanged\":true,\"binary80_DD_conversion_error\":0,\"GPU_matvec_max_difference\":"<<double(product_error)<<",\"resident_inner_iterations\":"<<answer.inner_iterations<<",\"outer_corrections\":"<<answer.outer_corrections<<",\"true_relative_residual\":"<<double(answer.true_relative_residual)<<",\"solution_max_difference\":"<<double(solution_error)<<",\"sparse_ghost_pair_bytes\":"<<bytes<<",\"actual_source_geometry_proven\":false,\"remap_committed\":false,\"new_physical_step\":false}\n";
    }MPI_Finalize();return 0;
    }catch(const std::exception&e){std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}
}
