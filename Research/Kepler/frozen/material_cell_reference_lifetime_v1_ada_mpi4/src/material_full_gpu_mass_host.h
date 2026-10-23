#ifndef MATERIAL_FULL_GPU_MASS_HOST_H
#define MATERIAL_FULL_GPU_MASS_HOST_H
#include "material_full_gpu_mass.h"
// This first offload keeps CG vectors and global dots in binary80 on the host;
// local block products/inverses execute in DD on the GPU. Halo communication
// also retains binary80. It does not claim device-resident vector/halo PCG.
namespace full_gpu_mass {
static MFDoubleDouble encode(R x){double h=double(x);return {h,double(x-R(h))};}
struct Resident;
struct Operator {
    const md::Mass&exact;void*ptr=nullptr;std::vector<MFDoubleDouble>input,output;
    double upload_seconds=0,kernel_seconds=0,operation_seconds=0;
    uint64_t apply_calls=0,inverse_calls=0;R conversion_error=0;Resident*resident=nullptr;
    explicit Operator(const md::Mass&M):exact(M),input(M.layout.global.size()),output(input.size()){
        double begin=MPI_Wtime();std::vector<MFDoubleDouble>block(M.block.size()),factor(block.size()),partition(M.partition.size());R maximum=0;
        #pragma omp parallel for schedule(static) reduction(max:maximum)
        for(size_t i=0;i<block.size();i++){R a=M.block[i],b=M.factor[i/(ST_MAX*ST_MAX)].L[i%(ST_MAX*ST_MAX)];block[i]=encode(a);factor[i]=encode(b);maximum=std::max({maximum,fabsl((R(block[i].hi)+block[i].lo)-a),fabsl((R(factor[i].hi)+factor[i].lo)-b)});}
        for(size_t i=0;i<partition.size();i++){partition[i]=encode(M.partition[i]);maximum=std::max(maximum,fabsl((R(partition[i].hi)+partition[i].lo)-M.partition[i]));}conversion_error=maximum;need(conversion_error==0,"binary80 fixed M/factor/partition exactly restored by hi+lo conversion");int bad=mf_mass_create(M.triangle.data(),M.triangle.size(),M.layout.global.size(),block.data(),factor.data(),partition.data(),M.row.data(),M.support_triangle.data(),M.support_row.data(),M.support_triangle.size(),&ptr)!=0;int any;MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,M.layout.comm);need(!any,"persistent DD GPU mass block/factor upload");upload_seconds=MPI_Wtime()-begin;
    }
    ~Operator();Operator(const Operator&)=delete;Operator&operator=(const Operator&)=delete;
    md::Vec operation(const md::Vec&x,bool inverse){double begin=MPI_Wtime();md::Vec halo=x,y(x.size());exact.layout.pull(halo);
        #pragma omp parallel for schedule(static)
        for(size_t i=0;i<x.size();i++)input[i]=encode(halo[i]);double seconds=0;int bad=mf_mass_apply(ptr,input.data(),output.data(),inverse?1:0,&seconds)!=0;
        if(!bad){
            #pragma omp parallel for schedule(static) reduction(+:bad)
            for(size_t i=0;i<x.size();i++){y[i]=R(output[i].hi)+output[i].lo;if(!std::isfinite(y[i]))bad++;}}
        int any;MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,exact.layout.comm);need(!any,"collective finite DD GPU mass operation");exact.layout.sum(y);kernel_seconds+=seconds;operation_seconds+=MPI_Wtime()-begin;if(inverse)inverse_calls++;else apply_calls++;return y;
    }
    md::Vec apply(const md::Vec&x){return operation(x,false);}
    md::Vec precondition(const md::Vec&x){return operation(x,true);}
};
struct Solution {md::Vec x;uint64_t inner_iterations=0,outer_corrections=0;R true_relative_residual=0;};
static md::Vec inner_host(Operator&A,const md::Vec&b,uint64_t&iterations){
    const auto&layout=A.exact.layout;md::Vec x(b.size()),r=b;R norm=layout.dot(b,b);need(norm>0&&std::isfinite(norm),"positive finite GPU correction RHS norm");auto z=A.precondition(r),p=z;R rz=layout.dot(r,z);need(rz>0,"positive GPU DD inner preconditioned metric");
    for(int it=0;it<256;it++){auto Ap=A.apply(p);R den=layout.dot(p,Ap);need(den>0&&std::isfinite(den),"positive GPU DD inner operator metric");R alpha=rz/den;for(size_t j=0;j<x.size();j++){x[j]+=alpha*p[j];r[j]-=alpha*Ap[j];}iterations++;if(layout.dot(r,r)<norm*1e-36L)return x;z=A.precondition(r);R next=layout.dot(r,z);need(next>0&&std::isfinite(next),"positive GPU DD inner next metric");R beta=next/rz;for(size_t j=0;j<p.size();j++)p[j]=z[j]+beta*p[j];rz=next;}
    need(false,"bounded 256-iteration DD GPU inner correction failed");return x;
}
#include "material_full_mass_resident_host.h"
static Solution solve_with_inner(Operator&A,const md::Vec&b,bool select_device){
    Solution out;out.x.resize(b.size());const auto&layout=A.exact.layout;R norm=layout.dot(b,b);need(norm>=0&&std::isfinite(norm),"finite GPU outer exact RHS");if(norm==0)return out;md::Vec r=b;R previous=1;
    for(int outer=0;outer<3;outer++){auto delta=select_device?inner(A,r,out.inner_iterations):inner_host(A,r,out.inner_iterations);for(size_t j=0;j<out.x.size();j++)out.x[j]+=delta[j];out.outer_corrections++;auto Mx=A.exact.apply(out.x);for(size_t j=0;j<r.size();j++)r[j]=b[j]-Mx[j];out.true_relative_residual=sqrtl(layout.dot(r,r)/norm);need(std::isfinite(out.true_relative_residual),"finite binary80 authoritative GPU outer true residual");if(out.true_relative_residual<1e-16L)return out;need(out.true_relative_residual<previous,"GPU outer true residual must strictly decrease");previous=out.true_relative_residual;}
    need(false,"three DD GPU corrections did not satisfy original binary80 true equation bound");return out;
}
static Solution solve(Operator&A,const md::Vec&b){return solve_with_inner(A,b,true);}
static Solution solve_host_reference(Operator&A,const md::Vec&b){return solve_with_inner(A,b,false);}
}
#endif
