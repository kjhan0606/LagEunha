#ifndef MATERIAL_FULL_MASS_RESIDENT_HOST_H
#define MATERIAL_FULL_MASS_RESIDENT_HOST_H
#include "material_full_mass_resident.h"
// Included inside full_gpu_mass after Operator and the reference host inner.
// DD pairs are staged only for sparse ghosts; all five full vectors stay on GPU.
struct Resident {
    Operator&A;void*ptr=nullptr;MFDoubleDouble*send=nullptr,*receive=nullptr;
    std::array<std::vector<uint64_t>,2>offset;
    std::vector<MFDoubleDouble>scalars;
    MPI_Datatype pair=MPI_DATATYPE_NULL;
    double halo_seconds=0,dot_seconds=0,total_seconds=0;
    uint64_t ghost_pair_bytes=0,full_vector_pair_bytes=0;
    void check(int code,const char*message){int bad=code!=0,any;MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,A.exact.layout.comm);need(!any,message);}
    explicit Resident(Operator&op):A(op),scalars(op.exact.layout.ranks){
        const auto&l=A.exact.layout;std::array<std::vector<uint64_t>,2>index;
        for(int kind=0;kind<2;kind++){offset[kind].resize(l.ranks+1);for(int p=0;p<l.ranks;p++){const auto&list=kind?l.ask[p]:l.serve[p];for(auto i:list)index[kind].push_back(i);offset[kind][p+1]=index[kind].size();}}
        std::vector<uint64_t>owned(l.owned.begin(),l.owned.end());check(mf_resident_create(A.ptr,owned.data(),owned.size(),index[0].data(),index[0].size(),index[1].data(),index[1].size(),&ptr),"persistent DD resident vectors/halo indices");
        size_t count=std::max(size_t(1),std::max(index[0].size(),index[1].size()));
        check(cudaMallocHost((void**)&send,count*sizeof(MFDoubleDouble)),"pinned resident ghost send staging");check(cudaMallocHost((void**)&receive,count*sizeof(MFDoubleDouble)),"pinned resident ghost receive staging");
        need(MPI_Type_contiguous(2,MPI_DOUBLE,&pair)==MPI_SUCCESS&&MPI_Type_commit(&pair)==MPI_SUCCESS,"MPI DD pair datatype");
        static_assert(sizeof(MFDoubleDouble)==2*sizeof(double),"DD MPI wire ABI");
    }
    ~Resident(){mf_resident_destroy(ptr);cudaFreeHost(send);cudaFreeHost(receive);if(pair!=MPI_DATATYPE_NULL)MPI_Type_free(&pair);}
    Resident(const Resident&)=delete;Resident&operator=(const Resident&)=delete;
    static MFDoubleDouble add_pair(MFDoubleDouble a,MFDoubleDouble b){double s=a.hi+b.hi,z=s-a.hi,e=(a.hi-(s-z))+(b.hi-z);e+=a.lo+b.lo;double h=s+e;return {h,e-(h-s)};}
    void exchange(int vector,bool pull){
        double start=MPI_Wtime();const auto&l=A.exact.layout;int outgoing=pull?0:1,incoming=1-outgoing;
        check(mf_resident_pack(ptr,vector,outgoing,send),"resident sparse DD ghost pack");std::vector<MPI_Request>request;
        for(int p=0;p<l.ranks;p++)if(p!=l.rank){for(uint64_t at=offset[incoming][p];at<offset[incoming][p+1];){int n=int(std::min<uint64_t>(INT_MAX,offset[incoming][p+1]-at));MPI_Request r;need(MPI_Irecv(receive+at,n,pair,p,pull?410:411,l.comm,&r)==MPI_SUCCESS,"resident DD ghost receive");request.push_back(r);at+=n;}}
        for(int p=0;p<l.ranks;p++)if(p!=l.rank){for(uint64_t at=offset[outgoing][p];at<offset[outgoing][p+1];){int n=int(std::min<uint64_t>(INT_MAX,offset[outgoing][p+1]-at));MPI_Request r;need(MPI_Isend(send+at,n,pair,p,pull?410:411,l.comm,&r)==MPI_SUCCESS,"resident DD ghost send");request.push_back(r);at+=n;}}
        need(request.size()<=INT_MAX,"bounded resident ghost requests");if(!request.empty())need(MPI_Waitall(int(request.size()),request.data(),MPI_STATUSES_IGNORE)==MPI_SUCCESS,"resident DD ghost transfer complete");
        // A given peer requests each ID once. Launch peers in canonical order;
        // corner contributions from distinct peers never race via atomics.
        for(int p=0;p<l.ranks;p++){size_t at=offset[incoming][p],n=offset[incoming][p+1]-at;check(mf_resident_scatter(ptr,vector,incoming,at,n,pull?0:1,receive+at),"deterministic resident DD peer scatter");}
        ghost_pair_bytes+=(offset[0].back()+offset[1].back())*sizeof(MFDoubleDouble);halo_seconds+=MPI_Wtime()-start;
    }
    void product(bool inverse){int input=inverse?MF_CG_R:MF_CG_P,output=inverse?MF_CG_Z:MF_CG_AP;exchange(input,true);check(mf_resident_product(ptr,inverse?1:0),"resident DD local fixed mass/inverse");exchange(output,false);exchange(output,true);}
    R dot(int a,int b){double start=MPI_Wtime();MFDoubleDouble local;check(mf_resident_dot(ptr,a,b,&local),"resident owned-only DD dot");const auto&l=A.exact.layout;need(MPI_Allgather(&local,1,pair,scalars.data(),1,pair,l.comm)==MPI_SUCCESS,"resident DD scalar reduction");MFDoubleDouble sum={0,0};for(auto x:scalars)sum=add_pair(sum,x);R result=R(sum.hi)+sum.lo;need(std::isfinite(result),"finite resident global DD scalar");dot_seconds+=MPI_Wtime()-start;return result;}
    md::Vec inner(const md::Vec&b,uint64_t&iterations){
        double begin=MPI_Wtime();
        #pragma omp parallel for schedule(static)
        for(size_t i=0;i<b.size();i++)A.input[i]=encode(b[i]);check(mf_resident_begin(ptr,A.input.data()),"resident correction RHS upload");full_vector_pair_bytes+=b.size()*sizeof(MFDoubleDouble);
        R norm=dot(MF_CG_R,MF_CG_R);need(norm>0,"positive resident DD correction RHS");product(true);check(mf_resident_update(ptr,0,{0,0}),"resident initial direction");R rz=dot(MF_CG_R,MF_CG_Z);need(rz>0,"positive resident DD initial preconditioned metric");
        for(int it=0;it<256;it++){
            product(false);R den=dot(MF_CG_P,MF_CG_AP);need(den>0,"positive resident DD operator metric");check(mf_resident_update(ptr,1,encode(rz/den)),"resident DD x/r update");iterations++;
            if(dot(MF_CG_R,MF_CG_R)<norm*1e-36L){check(mf_resident_finish(ptr,A.output.data()),"resident final correction download");full_vector_pair_bytes+=b.size()*sizeof(MFDoubleDouble);md::Vec x(b.size());int bad=0;
                #pragma omp parallel for schedule(static) reduction(+:bad)
                for(size_t i=0;i<x.size();i++){x[i]=R(A.output[i].hi)+A.output[i].lo;if(!std::isfinite(x[i]))bad++;}check(bad,"finite resident actual correction");double elapsed=MPI_Wtime()-begin;total_seconds+=elapsed;A.operation_seconds+=elapsed;return x;}
            product(true);R next=dot(MF_CG_R,MF_CG_Z);need(next>0,"positive resident DD next preconditioned metric");check(mf_resident_update(ptr,2,encode(next/rz)),"resident DD direction update");rz=next;
        }
        need(false,"bounded 256-iteration resident DD correction failed");return {};
    }
};
Operator::~Operator(){delete resident;mf_mass_destroy(ptr);}
static md::Vec inner(Operator&A,const md::Vec&b,uint64_t&iterations){
    const char*choice=getenv("GFS_MATERIAL_DEVICE_MASS");if(choice&&atoi(choice)>0){if(!A.resident)A.resident=new Resident(A);return A.resident->inner(b,iterations);}return inner_host(A,b,iterations);
}
#endif
