#ifndef MATERIAL_DISTRIBUTED_MASS_H
#define MATERIAL_DISTRIBUTED_MASS_H
// Owned/ghost material mass operator. Local element blocks only; no global
// replicated CSR/force array and no all-reduce of a whole global vector.
// Native binary80 accumulation/true residual preserve the CONTROL precision.
#include <mpi.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <climits>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>
#include "material_shared_trace_batch.h"
namespace material_distributed {
using R=long double;using Vec=std::vector<R>;
inline void require(bool p,const char*s){if(!p)throw std::runtime_error(s);}
// All peer lengths are64-bit; each actual MPI message is bounded by INT_MAX.
template<class T>inline void exchange(const std::vector<std::vector<T>>&send,
        std::vector<std::vector<T>>&receive,MPI_Datatype type,MPI_Comm comm,int rank,int tag){
    std::vector<MPI_Request>request;
    for(size_t p=0;p<receive.size();p++)if(int(p)!=rank){
        for(size_t at=0;at<receive[p].size();){int n=int(std::min<size_t>(INT_MAX,receive[p].size()-at));MPI_Request r;require(MPI_Irecv(receive[p].data()+at,n,type,p,tag,comm,&r)==MPI_SUCCESS,"material halo receive");request.push_back(r);at+=n;}
    }
    for(size_t p=0;p<send.size();p++)if(int(p)!=rank){
        for(size_t at=0;at<send[p].size();){int n=int(std::min<size_t>(INT_MAX,send[p].size()-at));MPI_Request r;require(MPI_Isend(send[p].data()+at,n,type,p,tag,comm,&r)==MPI_SUCCESS,"material halo send");request.push_back(r);at+=n;}
    }
    require(request.size()<=INT_MAX,"bounded material halo requests");if(!request.empty())require(MPI_Waitall(int(request.size()),request.data(),MPI_STATUSES_IGNORE)==MPI_SUCCESS,"complete material halo exchange");
}
struct Ownership {
    MPI_Comm comm=MPI_COMM_NULL;int rank=0,ranks=0;
    std::vector<uint64_t>global;std::vector<int>owner;std::vector<size_t>owned;
    // ask[p] are local ghosts requested from p; serve[p] are local owned
    // entries requested by p. Peer order is fixed for deterministic sums.
    std::vector<std::vector<size_t>>ask,serve;
    mutable std::vector<std::vector<R>>pull_send,pull_receive,sum_send,sum_receive;
    uint64_t global_owned=0,ghost_count=0;
    Ownership(std::vector<uint64_t>ids,std::vector<int>owners,MPI_Comm input):global(std::move(ids)),owner(std::move(owners)){
        require(MPI_Comm_dup(input,&comm)==MPI_SUCCESS,"material halo communicator");MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
        require(global.size()==owner.size()&&std::is_sorted(global.begin(),global.end())&&std::adjacent_find(global.begin(),global.end())==global.end(),"unique sorted local-to-global DOF IDs");
        ask.resize(ranks);serve.resize(ranks);std::vector<std::vector<uint64_t>>request(ranks),incoming(ranks);
        for(size_t i=0;i<global.size();i++){require(owner[i]>=0&&owner[i]<ranks,"valid material DOF owner");if(owner[i]==rank)owned.push_back(i);else{ask[owner[i]].push_back(i);request[owner[i]].push_back(global[i]);ghost_count++;}}
        std::vector<uint64_t>out(ranks),in(ranks);for(int p=0;p<ranks;p++)out[p]=request[p].size();MPI_Alltoall(out.data(),1,MPI_UINT64_T,in.data(),1,MPI_UINT64_T,comm);
        for(int p=0;p<ranks;p++){require(in[p]<=size_t(std::numeric_limits<size_t>::max()/sizeof(uint64_t)),"bounded halo metadata allocation");incoming[p].resize(in[p]);}
        exchange(request,incoming,MPI_UINT64_T,comm,rank,100);
        for(int p=0;p<ranks;p++)for(uint64_t id:incoming[p]){auto it=std::lower_bound(global.begin(),global.end(),id);require(it!=global.end()&&*it==id&&owner[it-global.begin()]==rank,"every requested ghost has its unique local owner");serve[p].push_back(it-global.begin());}
        pull_send.resize(ranks);pull_receive.resize(ranks);sum_send.resize(ranks);sum_receive.resize(ranks);
        for(int p=0;p<ranks;p++){pull_send[p].resize(serve[p].size());pull_receive[p].resize(ask[p].size());sum_send[p].resize(ask[p].size());sum_receive[p].resize(serve[p].size());}
        uint64_t n=owned.size();MPI_Allreduce(&n,&global_owned,1,MPI_UINT64_T,MPI_SUM,comm);require(global_owned>0,"nonempty global material DOF space");
    }
    ~Ownership(){if(comm!=MPI_COMM_NULL)MPI_Comm_free(&comm);}
    Ownership(const Ownership&)=delete;Ownership&operator=(const Ownership&)=delete;
    void pull(Vec&v)const{
        require(v.size()==global.size(),"material ghost vector extent");
        for(int p=0;p<ranks;p++)for(size_t j=0;j<serve[p].size();j++)pull_send[p][j]=v[serve[p][j]];
        exchange(pull_send,pull_receive,MPI_LONG_DOUBLE,comm,rank,101);
        for(int p=0;p<ranks;p++)for(size_t j=0;j<ask[p].size();j++)v[ask[p][j]]=pull_receive[p][j];
    }
    void sum(Vec&v,bool replicate_ghosts=true)const{
        require(v.size()==global.size(),"local material force contribution extent");
        for(int p=0;p<ranks;p++)for(size_t j=0;j<ask[p].size();j++)sum_send[p][j]=v[ask[p][j]];
        exchange(sum_send,sum_receive,MPI_LONG_DOUBLE,comm,rank,102);
        for(int p=0;p<ranks;p++)for(size_t j=0;j<serve[p].size();j++)v[serve[p][j]]+=sum_receive[p][j];
        if(replicate_ghosts)pull(v);
    }
    R dot(const Vec&a,const Vec&b)const{
        require(a.size()==global.size()&&b.size()==a.size(),"owned global dot extent");R local=0,total=0;
        // Keep one fixed summation order. No redundant ghost contribution.
        for(size_t i:owned)local+=a[i]*b[i];MPI_Allreduce(&local,&total,1,MPI_LONG_DOUBLE,MPI_SUM,comm);return total;
    }
};
struct Factor {
    std::array<R,ST_MAX*ST_MAX>L{};int n=0;
    Factor()=default;
    Factor(const R*M,int count):n(count){
        require(n>0&&n<=ST_MAX,"native material local mass extent");
        for(int i=0;i<n;i++)for(int j=0;j<=i;j++){require(M[i*ST_MAX+j]==M[j*ST_MAX+i],"symmetric fixed material mass block");R a=M[i*ST_MAX+j];for(int k=0;k<j;k++)a-=L[i*ST_MAX+k]*L[j*ST_MAX+k];require(i!=j||(a>0&&std::isfinite(a)),"positive fixed material local Cholesky metric");L[i*ST_MAX+j]=i==j?sqrtl(a):a/L[j*ST_MAX+j];}
    }
    void solve(R*x)const{
        for(int i=0;i<n;i++){for(int j=0;j<i;j++)x[i]-=L[i*ST_MAX+j]*x[j];x[i]/=L[i*ST_MAX+i];}
        for(int i=n-1;i>=0;i--){for(int j=i+1;j<n;j++)x[i]-=L[j*ST_MAX+i]*x[j];x[i]/=L[i*ST_MAX+i];}
    }
};
struct Mass {
    const Ownership&layout;const std::vector<STTriangle>&triangle;
    Vec block,partition;std::vector<Factor>factor;
    // Fixed element scratch is reused across preconditioner calls. Mass is
    // called collectively by the one integrator thread, never reentrantly.
    mutable Vec solved_scratch;
    mutable double apply_seconds=0,precondition_seconds=0;
    mutable uint64_t apply_calls=0,precondition_calls=0;
    double factor_seconds=0,setup_seconds=0;
    std::vector<uint64_t>row;std::vector<uint32_t>support_triangle;
    std::vector<uint8_t>support_row;
    Mass(const Ownership&o,const std::vector<STTriangle>&t,Vec material_blocks):layout(o),triangle(t),block(std::move(material_blocks)),partition(o.global.size()),row(o.global.size()+1){
        double begin=MPI_Wtime();require(t.size()<=UINT32_MAX&&block.size()==t.size()*size_t(ST_MAX)*ST_MAX,"owned matrix-free block extent");factor.resize(t.size());solved_scratch.resize(t.size()*ST_MAX);int failed=0;
        #pragma omp parallel for schedule(static) reduction(+:failed)
        for(size_t s=0;s<t.size();s++)try{factor[s]=Factor(block.data()+s*ST_MAX*ST_MAX,t[s].n);}catch(...){failed++;}
        int any;MPI_Allreduce(&failed,&any,1,MPI_INT,MPI_MAX,layout.comm);require(!any,"all local fixed mass factors positive");factor_seconds=MPI_Wtime()-begin;
        for(size_t s=0;s<t.size();s++){const auto&T=t[s];for(int a=0;a<T.n;a++){require(T.id[a]>=0&&size_t(T.id[a])<partition.size(),"owned triangle local DOF index");row[T.id[a]+1]++;partition[T.id[a]]++;}}
        for(size_t i=1;i<row.size();i++)row[i]+=row[i-1];support_triangle.resize(row.back());support_row.resize(row.back());auto at=row;
        for(size_t s=0;s<t.size();s++)for(int a=0;a<t[s].n;a++){auto j=at[t[s].id[a]]++;support_triangle[j]=s;support_row[j]=a;}
        layout.sum(partition);for(auto&v:partition){require(v>0,"every local owned/ghost DOF has material support");v=1/v;}
        setup_seconds=MPI_Wtime()-begin;
    }
    Vec apply(const Vec&x)const{
        double begin=MPI_Wtime();apply_calls++;
        Vec input=x;layout.pull(input);Vec y(x.size());
        #pragma omp parallel for schedule(static)
        for(size_t i=0;i<y.size();i++){R value=0;for(uint64_t j=row[i];j<row[i+1];j++){size_t s=support_triangle[j];int a=support_row[j];const auto&T=triangle[s];const R*M=block.data()+s*ST_MAX*ST_MAX;for(int b=0;b<T.n;b++)value+=M[a*ST_MAX+b]*input[T.id[b]];}y[i]=value;}
        layout.sum(y);apply_seconds+=MPI_Wtime()-begin;return y;
    }
    Vec precondition(const Vec&r)const{
        double begin=MPI_Wtime();precondition_calls++;Vec input=r;layout.pull(input);auto&solved=solved_scratch;Vec z(r.size());
        #pragma omp parallel for schedule(static)
        for(size_t s=0;s<triangle.size();s++){R*x=solved.data()+s*ST_MAX;for(int a=0;a<triangle[s].n;a++){int i=triangle[s].id[a];x[a]=input[i]*partition[i];}factor[s].solve(x);}
        #pragma omp parallel for schedule(static)
        for(size_t i=0;i<z.size();i++){R value=0;for(uint64_t j=row[i];j<row[i+1];j++)value+=solved[size_t(support_triangle[j])*ST_MAX+support_row[j]];z[i]=partition[i]*value;}
        layout.sum(z);precondition_seconds+=MPI_Wtime()-begin;return z;
    }
};
struct Solution {Vec x;uint64_t iterations=0;R true_relative_residual=0;};
inline Solution solve(const Mass&M,const Vec&b){
    Solution out;out.x.resize(b.size());R norm=M.layout.dot(b,b);require(norm>=0&&std::isfinite(norm),"finite native mass RHS norm");if(norm==0)return out;
    Vec r=b,z=M.precondition(r),p=z;R rz=M.layout.dot(r,z);require(rz>0,"positive distributed preconditioned norm");
    for(uint64_t i=0;i<4*M.layout.global_owned;i++){
        auto Ap=M.apply(p);R den=M.layout.dot(p,Ap);require(den>0&&std::isfinite(den),"positive distributed material mass operator");R alpha=rz/den;
        for(size_t j=0;j<b.size();j++){out.x[j]+=alpha*p[j];r[j]-=alpha*Ap[j];}out.iterations=i+1;
        if(M.layout.dot(r,r)<norm*1e-32L){
            auto exact=M.apply(out.x);for(size_t j=0;j<exact.size();j++)exact[j]=b[j]-exact[j];out.true_relative_residual=sqrtl(M.layout.dot(exact,exact)/norm);
            if(out.true_relative_residual<1e-16L)return out;
            r=std::move(exact);z=M.precondition(r);p=z;rz=M.layout.dot(r,z);require(rz>0,"positive restarted true residual metric");continue;
        }
        z=M.precondition(r);R next=M.layout.dot(r,z);require(next>0&&std::isfinite(next),"positive distributed local mass preconditioner");R beta=next/rz;for(size_t j=0;j<p.size();j++)p[j]=z[j]+beta*p[j];rz=next;
    }
    require(false,"distributed mass solve did not meet original true equation tolerance");return out;
}
}
#endif
