#ifndef MATERIAL_COMPACT_SOLVER_H
#define MATERIAL_COMPACT_SOLVER_H
#include <cuda_runtime.h>
#include "material_compact_gpu.h"
#include "material_compact_injection_io.h"
#include "material_distributed_mass.h"
namespace compact_solver {
using R=long double;using U=uint64_t;using Vec=std::vector<R>;using DD=MFDoubleDouble;
namespace md=material_distributed;
inline DD encode(R x){double hi=double(x);return {hi,double(x-R(hi))};}
inline DD pair_add(DD a,DD b){double s=a.hi+b.hi,z=s-a.hi,e=(a.hi-(s-z))+(b.hi-z);e+=a.lo+b.lo;double h=s+e;return {h,e-(h-s)};}
inline void gpu(int code,const char*s){if(code){fprintf(stderr,"COMPACT_CUDA_REJECTED %s code%d\n",s,code);MPI_Abort(MPI_COMM_WORLD,2);}}
struct Matrix {
    const md::Ownership&layout;
    std::vector<U>row;std::vector<uint32_t>column;Vec value,inverse_diagonal;
    std::array<Vec,2>B;R C[4]={},scale[2]={};U elements=0;
    Matrix(const md::Ownership&l,const compact_io::Injection&t,const std::string&path):layout(l){
        compact_io::Reader in(path);U h[8];in.get(h,8);
        need(h[0]==0x4d434d4153533031ULL&&h[1]==1&&h[2]==U(l.rank)&&h[3]==U(l.ranks)&&h[5]==t.total&&h[6]==sizeof(R)&&h[7]==t.old_total,"compact mass record identity");elements=h[4];
        struct Triplet {uint32_t i,j,sequence;R value;};std::vector<Triplet>entry;entry.reserve(elements*100);
        for(auto&b:B)b.resize(l.global.size());
        for(U s=0;s<elements;s++){
            U d[3];in.get(d,3);size_t np=d[2],n=np+2;need((d[1]==1||d[1]==4)&&np>=10&&np<=13,"compact mass block support");
            std::vector<uint32_t>index(np);std::vector<R>mass(n*n),force(2*n);in.get(index.data(),np);in.get(mass.data(),mass.size());in.get(force.data(),force.size());
            for(size_t a=0;a<np;a++){
                need(index[a]<l.global.size()&&std::isfinite(mass[a*n+a]),"compact polynomial index/diagonal");
                for(size_t b=0;b<np;b++)if(mass[a*n+b]!=0){need(entry.size()<UINT32_MAX,"bounded local sparse insertion order");entry.push_back({index[a],index[b],uint32_t(entry.size()),mass[a*n+b]});}
                for(int b=0;b<2;b++)B[b][index[a]]+=mass[a*n+np+b];
            }
            for(int a=0;a<2;a++)for(int b=0;b<2;b++)C[2*a+b]+=mass[(np+a)*n+np+b];
        }in.finish();
        std::sort(entry.begin(),entry.end(),[](const Triplet&a,const Triplet&b){return a.i!=b.i?a.i<b.i:a.j!=b.j?a.j<b.j:a.sequence<b.sequence;});
        row.resize(l.global.size()+1);
        for(size_t k=0;k<entry.size();){uint32_t i=entry[k].i,j=entry[k].j;R sum=0;do{sum+=entry[k++].value;}while(k<entry.size()&&entry[k].i==i&&entry[k].j==j);column.push_back(j);value.push_back(sum);row[i+1]++;}
        for(size_t i=1;i<row.size();i++)row[i]+=row[i-1];std::vector<Triplet>().swap(entry);
        inverse_diagonal.resize(l.global.size());for(size_t i=0;i<l.global.size();i++)for(U k=row[i];k<row[i+1];k++)if(column[k]==i)inverse_diagonal[i]+=value[k];
        l.sum(inverse_diagonal);for(R&x:inverse_diagonal){need(x>0&&std::isfinite(x),"positive assembled polynomial diagonal");x=1/x;}
        for(auto&b:B)l.sum(b);MPI_Allreduce(MPI_IN_PLACE,C,4,MPI_LONG_DOUBLE,MPI_SUM,l.comm);
        for(int d=0;d<2;d++){need(C[3*d]>0,"positive global raw residual norm");scale[d]=sqrtl(C[3*d]);for(R&x:B[d])x/=scale[d];}
        for(int a=0;a<2;a++)for(int b=0;b<2;b++)C[2*a+b]/=scale[a]*scale[b];
    }
    Vec apply(const Vec&x)const{
        Vec input=x;layout.pull(input);Vec y(x.size());
        #pragma omp parallel for schedule(static)
        for(size_t i=0;i<y.size();i++){R sum=0;for(U k=row[i];k<row[i+1];k++)sum+=value[k]*input[column[k]];y[i]=sum;}
        layout.sum(y);return y;
    }
};
struct Resident {
    const Matrix&A;void*ptr=nullptr;MPI_Datatype pair=MPI_DATATYPE_NULL;
    DD*send=nullptr,*receive=nullptr;std::array<std::vector<U>,2>offset;
    std::vector<DD>input,output,scalars;U iterations=0,corrections=0,ghost_bytes=0;
    double inner_seconds=0,halo_seconds=0;
    explicit Resident(const Matrix&matrix):A(matrix),input(matrix.layout.global.size()),output(input.size()),scalars(matrix.layout.ranks){
        const auto&l=A.layout;std::vector<DD>values(A.value.size()),diagonal(input.size());
        int bad=0;
        #pragma omp parallel for schedule(static) reduction(+:bad)
        for(size_t i=0;i<values.size();i++){values[i]=encode(A.value[i]);if(R(values[i].hi)+values[i].lo!=A.value[i])bad++;}
        need(!bad,"exact binary80 to DD CSR restore");for(size_t i=0;i<diagonal.size();i++)diagonal[i]=encode(A.inverse_diagonal[i]);
        std::array<std::vector<U>,2>index;
        for(int k=0;k<2;k++){offset[k].resize(l.ranks+1);for(int p=0;p<l.ranks;p++){const auto&list=k?l.ask[p]:l.serve[p];for(auto i:list)index[k].push_back(i);offset[k][p+1]=index[k].size();}}
        std::vector<U>owned(l.owned.begin(),l.owned.end());gpu(cm_create(input.size(),A.row.data(),A.column.data(),values.data(),diagonal.data(),owned.data(),owned.size(),index[0].data(),index[0].size(),index[1].data(),index[1].size(),&ptr),"resident compact CSR allocation");
        size_t n=std::max(size_t(1),std::max(index[0].size(),index[1].size()));gpu(cudaMallocHost((void**)&send,n*sizeof(DD)),"pinned compact send");gpu(cudaMallocHost((void**)&receive,n*sizeof(DD)),"pinned compact receive");
        need(MPI_Type_contiguous(2,MPI_DOUBLE,&pair)==MPI_SUCCESS&&MPI_Type_commit(&pair)==MPI_SUCCESS,"compact DD MPI pair");
    }
    ~Resident(){cm_destroy(ptr);cudaFreeHost(send);cudaFreeHost(receive);if(pair!=MPI_DATATYPE_NULL)MPI_Type_free(&pair);}
    void exchange(int vector,bool pull){
        double begin=MPI_Wtime();const auto&l=A.layout;int outgoing=pull?0:1,incoming=1-outgoing;
        gpu(cm_pack(ptr,vector,outgoing,send),"compact ghost pack");std::vector<MPI_Request>request;
        for(int p=0;p<l.ranks;p++)if(p!=l.rank){
            for(U at=offset[incoming][p];at<offset[incoming][p+1];){int n=int(std::min<U>(INT_MAX,offset[incoming][p+1]-at));MPI_Request r;need(MPI_Irecv(receive+at,n,pair,p,pull?610:611,l.comm,&r)==MPI_SUCCESS,"compact ghost receive");request.push_back(r);at+=n;}
        }
        for(int p=0;p<l.ranks;p++)if(p!=l.rank){
            for(U at=offset[outgoing][p];at<offset[outgoing][p+1];){int n=int(std::min<U>(INT_MAX,offset[outgoing][p+1]-at));MPI_Request r;need(MPI_Isend(send+at,n,pair,p,pull?610:611,l.comm,&r)==MPI_SUCCESS,"compact ghost send");request.push_back(r);at+=n;}
        }
        need(request.size()<=INT_MAX,"bounded compact MPI requests");if(!request.empty())need(MPI_Waitall(int(request.size()),request.data(),MPI_STATUSES_IGNORE)==MPI_SUCCESS,"compact halo complete");
        for(int p=0;p<l.ranks;p++){size_t at=offset[incoming][p],n=offset[incoming][p+1]-at;gpu(cm_scatter(ptr,vector,incoming,at,n,pull?0:1,receive+at),"compact deterministic peer scatter");}
        halo_seconds+=MPI_Wtime()-begin;ghost_bytes+=(offset[0].back()+offset[1].back())*sizeof(DD);
    }
    R dot(int a,int b){DD local;gpu(cm_dot(ptr,a,b,&local),"compact owned DD dot");need(MPI_Allgather(&local,1,pair,scalars.data(),1,pair,A.layout.comm)==MPI_SUCCESS,"compact scalar gather");DD sum={0,0};for(auto x:scalars)sum=pair_add(sum,x);return R(sum.hi)+sum.lo;}
    void product(bool inverse){
        exchange(inverse?MF_CG_R:MF_CG_P,true);gpu(cm_product(ptr,inverse),"compact CSR/diagonal product");
        // Jacobi uses an already globally assembled diagonal: do not sum it
        // over peers again. Matrix partial rows do require one owner sum.
        if(!inverse){exchange(MF_CG_AP,false);exchange(MF_CG_AP,true);}
    }
    Vec inner(const Vec&b){
        double begin=MPI_Wtime();for(size_t i=0;i<b.size();i++)input[i]=encode(b[i]);gpu(cm_begin(ptr,input.data()),"compact RHS upload");
        R norm=dot(MF_CG_R,MF_CG_R);if(norm==0)return Vec(b.size());need(norm>0,"positive compact RHS norm");product(true);gpu(cm_update(ptr,0,{0,0}),"compact initial search direction");R rz=dot(MF_CG_R,MF_CG_Z);need(rz>0,"positive compact preconditioned norm");
        for(int k=0;k<256;k++){
            product(false);R denominator=dot(MF_CG_P,MF_CG_AP);need(denominator>0&&std::isfinite(denominator),"positive assembled compact polynomial metric");gpu(cm_update(ptr,1,encode(rz/denominator)),"compact x/r update");iterations++;
            if(dot(MF_CG_R,MF_CG_R)<norm*1e-36L){gpu(cm_finish(ptr,output.data()),"compact solution download");Vec out(output.size());for(size_t i=0;i<out.size();i++)out[i]=R(output[i].hi)+output[i].lo;inner_seconds+=MPI_Wtime()-begin;return out;}
            product(true);R next=dot(MF_CG_R,MF_CG_Z);need(next>0,"positive compact next recurrence");gpu(cm_update(ptr,2,encode(next/rz)),"compact search direction update");rz=next;
        }need(false,"fixed 256-iteration compact DD inner bound");return {};
    }
    Vec solve(const Vec&b){
        Vec x(b.size()),r=b;R norm=A.layout.dot(b,b);if(norm==0)return x;
        for(int outer=0;outer<3;outer++){
            Vec correction=inner(r);corrections++;for(size_t i=0;i<x.size();i++)x[i]+=correction[i];auto ax=A.apply(x);for(size_t i=0;i<r.size();i++)r[i]=b[i]-ax[i];
            if(sqrtl(A.layout.dot(r,r)/norm)<1e-16L)return x;
        }need(false,"compact binary80 polynomial true residual exceeds1e-16");return {};
    }
};
struct Bordered {
    Matrix&A;Resident&gpu_solver;std::array<Vec,2>W;R S[4]={};double setup_seconds=0;
    Bordered(Matrix&a,Resident&resident):A(a),gpu_solver(resident){
        double begin=MPI_Wtime();for(int d=0;d<2;d++)W[d]=gpu_solver.solve(A.B[d]);
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)S[2*i+j]=A.C[2*i+j]-A.layout.dot(A.B[i],W[j]);
        R scale=std::max(fabsl(S[0]),fabsl(S[3]));need(scale>0&&fabsl(S[1]-S[2])<1e-14L*scale,"computed Schur symmetry accuracy");
        need(S[0]>0&&S[3]>0&&S[0]*S[3]-S[1]*S[2]>0,"positive computed global rank2 Schur complement");setup_seconds=MPI_Wtime()-begin;
    }
    Vec solve(const Vec&full,const compact_io::Injection&t){
        Vec poly(full.begin(),full.begin()+A.layout.global.size());Vec y=gpu_solver.solve(poly);R b[2];
        for(int d=0;d<2;d++)b[d]=full[t.border[d]]-A.layout.dot(A.B[d],y);
        R determinant=S[0]*S[3]-S[1]*S[2],z[2]={(S[3]*b[0]-S[1]*b[1])/determinant,(-S[2]*b[0]+S[0]*b[1])/determinant};
        Vec out(full.size());for(size_t i=0;i<y.size();i++)out[i]=y[i]-W[0][i]*z[0]-W[1][i]*z[1];for(int d=0;d<2;d++)out[t.border[d]]=z[d];return out;
    }
};
}
#endif
