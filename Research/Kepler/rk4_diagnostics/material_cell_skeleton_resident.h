// DD CSR resident implementation reused from the previously verified backend.
// No compact trial-space construction or residual border is included.
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
        std::vector<U>owned(l.owned.begin(),l.owned.end());gpu(cm_create(input.size(),A.row.data(),A.column.data(),values.data(),diagonal.data(),owned.data(),owned.size(),index[0].data(),index[0].size(),index[1].data(),index[1].size(),&ptr),"resident skeleton CSR allocation");
        size_t n=std::max(size_t(1),std::max(index[0].size(),index[1].size()));gpu(cudaMallocHost((void**)&send,n*sizeof(DD)),"pinned skeleton send");gpu(cudaMallocHost((void**)&receive,n*sizeof(DD)),"pinned skeleton receive");
        need(MPI_Type_contiguous(2,MPI_DOUBLE,&pair)==MPI_SUCCESS&&MPI_Type_commit(&pair)==MPI_SUCCESS,"skeleton DD MPI pair");
    }
    ~Resident(){cm_destroy(ptr);cudaFreeHost(send);cudaFreeHost(receive);if(pair!=MPI_DATATYPE_NULL)MPI_Type_free(&pair);}
    void exchange(int vector,bool pull){
        double begin=MPI_Wtime();const auto&l=A.layout;int outgoing=pull?0:1,incoming=1-outgoing;
        gpu(cm_pack(ptr,vector,outgoing,send),"skeleton ghost pack");std::vector<MPI_Request>request;
        for(int p=0;p<l.ranks;p++)if(p!=l.rank){
            for(U at=offset[incoming][p];at<offset[incoming][p+1];){int n=int(std::min<U>(INT_MAX,offset[incoming][p+1]-at));MPI_Request r;need(MPI_Irecv(receive+at,n,pair,p,pull?610:611,l.comm,&r)==MPI_SUCCESS,"skeleton ghost receive");request.push_back(r);at+=n;}
        }
        for(int p=0;p<l.ranks;p++)if(p!=l.rank){
            for(U at=offset[outgoing][p];at<offset[outgoing][p+1];){int n=int(std::min<U>(INT_MAX,offset[outgoing][p+1]-at));MPI_Request r;need(MPI_Isend(send+at,n,pair,p,pull?610:611,l.comm,&r)==MPI_SUCCESS,"skeleton ghost send");request.push_back(r);at+=n;}
        }
        need(request.size()<=INT_MAX,"bounded skeleton MPI requests");if(!request.empty())need(MPI_Waitall(int(request.size()),request.data(),MPI_STATUSES_IGNORE)==MPI_SUCCESS,"skeleton halo complete");
        for(int p=0;p<l.ranks;p++){size_t at=offset[incoming][p],n=offset[incoming][p+1]-at;gpu(cm_scatter(ptr,vector,incoming,at,n,pull?0:1,receive+at),"skeleton deterministic peer scatter");}
        halo_seconds+=MPI_Wtime()-begin;ghost_bytes+=(offset[0].back()+offset[1].back())*sizeof(DD);
    }
    R dot(int a,int b){DD local;gpu(cm_dot(ptr,a,b,&local),"skeleton owned DD dot");need(MPI_Allgather(&local,1,pair,scalars.data(),1,pair,A.layout.comm)==MPI_SUCCESS,"skeleton scalar gather");DD sum={0,0};for(auto x:scalars)sum=pair_add(sum,x);return R(sum.hi)+sum.lo;}
    void product(bool inverse){
        exchange(inverse?MF_CG_R:MF_CG_P,true);gpu(cm_product(ptr,inverse),"skeleton CSR/diagonal product");
        // Jacobi uses an already globally assembled diagonal: do not sum it
        // over peers again. Matrix partial rows do require one owner sum.
        if(!inverse){exchange(MF_CG_AP,false);exchange(MF_CG_AP,true);}
    }
    Vec inner(const Vec&b){
        double begin=MPI_Wtime();for(size_t i=0;i<b.size();i++)input[i]=encode(b[i]);gpu(cm_begin(ptr,input.data()),"skeleton RHS upload");
        R norm=dot(MF_CG_R,MF_CG_R);if(norm==0)return Vec(b.size());need(norm>0,"positive skeleton RHS norm");product(true);gpu(cm_update(ptr,0,{0,0}),"skeleton initial search direction");R rz=dot(MF_CG_R,MF_CG_Z);need(rz>0,"positive skeleton preconditioned norm");
        for(int k=0;k<256;k++){
            product(false);R denominator=dot(MF_CG_P,MF_CG_AP);need(denominator>0&&std::isfinite(denominator),"positive assembled skeleton skeleton metric");gpu(cm_update(ptr,1,encode(rz/denominator)),"skeleton x/r update");iterations++;
            if(dot(MF_CG_R,MF_CG_R)<norm*1e-36L){gpu(cm_finish(ptr,output.data()),"skeleton solution download");Vec out(output.size());for(size_t i=0;i<out.size();i++)out[i]=R(output[i].hi)+output[i].lo;inner_seconds+=MPI_Wtime()-begin;return out;}
            product(true);R next=dot(MF_CG_R,MF_CG_Z);need(next>0,"positive skeleton next recurrence");gpu(cm_update(ptr,2,encode(next/rz)),"skeleton search direction update");rz=next;
        }need(false,"fixed 256-iteration skeleton DD inner bound");return {};
    }
    Vec solve(const Vec&b){
        Vec x(b.size()),r=b;R norm=A.layout.dot(b,b);if(norm==0)return x;
        for(int outer=0;outer<3;outer++){
            Vec correction=inner(r);corrections++;for(size_t i=0;i<x.size();i++)x[i]+=correction[i];auto ax=A.apply(x);for(size_t i=0;i<r.size();i++)r[i]=b[i]-ax[i];
            if(sqrtl(A.layout.dot(r,r)/norm)<1e-16L)return x;
        }need(false,"skeleton binary80 skeleton true residual exceeds1e-16");return {};
    }
};
