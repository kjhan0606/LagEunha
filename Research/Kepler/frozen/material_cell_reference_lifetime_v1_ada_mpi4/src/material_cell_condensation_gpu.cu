// Reuse the original physical kernels and assembled DD CSR backend. No use of
// the rejected compact injection/extra residual modes or a different field.
#include "material_compact_gpu.cu"
#include "material_cell_condensation_gpu.h"
enum{MCC_MAX_I=463};
struct MCCContext {
    size_t nf=0,ng=0,nc=0,ng_sum=0;uint64_t bytes=0;
    MCCell*cell=nullptr;uint32_t*I=nullptr,*G=nullptr,*old_G=nullptr;
    uint8_t*owned_G=nullptr;uint64_t*row=nullptr,*entry=nullptr;
    DD*L=nullptr,*IG=nullptr,*rhs[2]={},*solution[2]={},*contribution[2]={},*skeleton[2]={};
};
__device__ static size_t mcc_packed(size_t i,size_t j){return i*(i+1)/2+j;}
__device__ static DD mcc_warp_sum(DD x){
    int lane=threadIdx.x;
    for(int d=16;d;d/=2){DD y={__shfl_down_sync(0xffffffff,x.hi,d),__shfl_down_sync(0xffffffff,x.lo,d)};if(lane+d<32)x=add(x,y);}
    return x;
}
// One complete warp per source. No atomic updates or inter-source interior
// coupling. Both vector components use the very same packed factor.
__device__ static void mcc_solve(const MCCell&c,const DD*L,DD*x,DD*y){
    int lane=threadIdx.x;
    for(uint32_t i=0;i<c.ni;i++){
        DD sx=dd(0),sy=dd(0);
        for(uint32_t k=lane;k<i;k+=32){DD a=L[c.L+mcc_packed(i,k)];sx=add(sx,mul(a,x[k]));sy=add(sy,mul(a,y[k]));}
        sx=mcc_warp_sum(sx);sy=mcc_warp_sum(sy);
        if(lane==0){DD a=L[c.L+mcc_packed(i,i)];x[i]=divide(sub(x[i],sx),a);y[i]=divide(sub(y[i],sy),a);}__syncwarp();
    }
    for(int i=int(c.ni)-1;i>=0;i--){
        DD sx=dd(0),sy=dd(0);
        for(uint32_t k=i+1+lane;k<c.ni;k+=32){DD a=L[c.L+mcc_packed(k,i)];sx=add(sx,mul(a,x[k]));sy=add(sy,mul(a,y[k]));}
        sx=mcc_warp_sum(sx);sy=mcc_warp_sum(sy);
        if(lane==0){DD a=L[c.L+mcc_packed(i,i)];x[i]=divide(sub(x[i],sx),a);y[i]=divide(sub(y[i],sy),a);}__syncwarp();
    }
}
__global__ static void mcc_condense_kernel(const MCCell*cell,const uint32_t*I,const DD*L,const DD*IG,
        const DD*fx,const DD*fy,DD*gx,DD*gy){
    const auto&c=cell[blockIdx.x];int lane=threadIdx.x;__shared__ DD x[MCC_MAX_I],y[MCC_MAX_I];
    for(uint32_t i=lane;i<c.ni;i+=32){uint32_t old=I[c.I+i];x[i]=fx[old];y[i]=fy[old];}__syncwarp();mcc_solve(c,L,x,y);
    for(uint32_t g=lane;g<c.ng;g+=32){DD sx=dd(0),sy=dd(0);for(uint32_t i=0;i<c.ni;i++){DD a=IG[c.IG+size_t(i)*c.ng+g];sx=add(sx,mul(a,x[i]));sy=add(sy,mul(a,y[i]));}gx[c.G+g]=neg(sx);gy[c.G+g]=neg(sy);}
}
__global__ static void mcc_gather_kernel(const uint64_t*row,const uint64_t*entry,const uint32_t*old_G,const uint8_t*owned,
        const DD*cx,const DD*cy,const DD*fx,const DD*fy,DD*gx,DD*gy,size_t n){
    size_t g=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(g>=n)return;
    // f_G is already globally assembled. Add it exactly once at its owner;
    // sum only the source-local negative GI MII^-1 fI contributions over MPI.
    DD x=owned[g]?fx[old_G[g]]:dd(0),y=owned[g]?fy[old_G[g]]:dd(0);
    for(uint64_t k=row[g];k<row[g+1];k++){x=add(x,cx[entry[k]]);y=add(y,cy[entry[k]]);}gx[g]=x;gy[g]=y;
}
__global__ static void mcc_copy_trace_kernel(const uint32_t*old_G,const DD*gx,const DD*gy,DD*ax,DD*ay,size_t n){size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n){ax[old_G[i]]=gx[i];ay[old_G[i]]=gy[i];}}
__global__ static void mcc_recover_kernel(const MCCell*cell,const uint32_t*I,const uint32_t*G,const DD*L,const DD*IG,
        const DD*fx,const DD*fy,const DD*gx,const DD*gy,DD*ax,DD*ay){
    const auto&c=cell[blockIdx.x];int lane=threadIdx.x;__shared__ DD x[MCC_MAX_I],y[MCC_MAX_I];
    for(uint32_t i=lane;i<c.ni;i+=32){uint32_t old=I[c.I+i];DD sx=fx[old],sy=fy[old];for(uint32_t g=0;g<c.ng;g++){uint32_t j=G[c.G+g];DD a=IG[c.IG+size_t(i)*c.ng+g];sx=sub(sx,mul(a,gx[j]));sy=sub(sy,mul(a,gy[j]));}x[i]=sx;y[i]=sy;}
    __syncwarp();mcc_solve(c,L,x,y);
    for(uint32_t i=lane;i<c.ni;i+=32){uint32_t old=I[c.I+i];ax[old]=x[i];ay[old]=y[i];}
}
extern "C" void mcc_destroy(void*ptr){auto*s=static_cast<MCCContext*>(ptr);if(!s)return;cudaFree(s->cell);cudaFree(s->I);cudaFree(s->G);cudaFree(s->old_G);cudaFree(s->owned_G);cudaFree(s->row);cudaFree(s->entry);cudaFree(s->L);cudaFree(s->IG);for(int d=0;d<2;d++){cudaFree(s->rhs[d]);cudaFree(s->solution[d]);cudaFree(s->contribution[d]);cudaFree(s->skeleton[d]);}delete s;}
extern "C" uint64_t mcc_device_bytes(void*ptr){auto*s=static_cast<MCCContext*>(ptr);return s?s->bytes:0;}
extern "C" int mcc_create(size_t nf,size_t ng,const MCCell*cell,size_t nc,const uint32_t*I,size_t ni,const uint32_t*G,size_t ngsum,
        const uint32_t*old_G,const uint8_t*owned,const uint64_t*row,const uint64_t*entry,const DD*L,size_t nl,const DD*IG,size_t nig,void**out){
    *out=nullptr;if(!nf||!ng||!nc||nf>UINT32_MAX||ng>UINT32_MAX||nc>UINT32_MAX||ni+ng!=nf||row[ng]!=ngsum)return cudaErrorInvalidValue;
    for(size_t c=0;c<nc;c++)if(!cell[c].ni||cell[c].ni>MCC_MAX_I||!cell[c].ng||cell[c].I+cell[c].ni>ni||cell[c].G+cell[c].ng>ngsum||cell[c].L+size_t(cell[c].ni)*(cell[c].ni+1)/2>nl||cell[c].IG+size_t(cell[c].ni)*cell[c].ng>nig)return cudaErrorInvalidValue;
    auto*s=new MCCContext;s->nf=nf;s->ng=ng;s->nc=nc;s->ng_sum=ngsum;cudaError_t e=cudaSuccess;
    #define MCC_ALLOC(field,count) do{if(e==cudaSuccess){size_t bytes=(count)*sizeof(*s->field);e=cudaMalloc((void**)&s->field,bytes);if(e==cudaSuccess)s->bytes+=bytes;}}while(0)
    #define MCC_COPY(field,input,count) do{MCC_ALLOC(field,count);if(e==cudaSuccess)e=cudaMemcpy(s->field,input,(count)*sizeof(*s->field),cudaMemcpyHostToDevice);}while(0)
    MCC_COPY(cell,cell,nc);MCC_COPY(I,I,ni);MCC_COPY(G,G,ngsum);MCC_COPY(old_G,old_G,ng);MCC_COPY(owned_G,owned,ng);MCC_COPY(row,row,ng+1);MCC_COPY(entry,entry,ngsum);MCC_COPY(L,L,nl);MCC_COPY(IG,IG,nig);
    for(int d=0;d<2;d++){MCC_ALLOC(rhs[d],nf);MCC_ALLOC(solution[d],nf);MCC_ALLOC(contribution[d],ngsum);MCC_ALLOC(skeleton[d],ng);}
    #undef MCC_COPY
    #undef MCC_ALLOC
    if(e!=cudaSuccess){mcc_destroy(s);return e;}*out=s;return 0;
}
extern "C" int mcc_rhs(void*ptr,const DD*fx,const DD*fy,DD*gx,DD*gy){
    auto*s=static_cast<MCCContext*>(ptr);if(!s)return cudaErrorInvalidValue;cudaError_t e;
    #define MCC_TRY(call) do{e=(call);if(e!=cudaSuccess)return e;}while(0)
    MCC_TRY(cudaMemcpy(s->rhs[0],fx,s->nf*sizeof(DD),cudaMemcpyHostToDevice));MCC_TRY(cudaMemcpy(s->rhs[1],fy,s->nf*sizeof(DD),cudaMemcpyHostToDevice));
    mcc_condense_kernel<<<s->nc,32>>>(s->cell,s->I,s->L,s->IG,s->rhs[0],s->rhs[1],s->contribution[0],s->contribution[1]);MCC_TRY(cudaGetLastError());
    mcc_gather_kernel<<<(s->ng+255)/256,256>>>(s->row,s->entry,s->old_G,s->owned_G,s->contribution[0],s->contribution[1],s->rhs[0],s->rhs[1],s->skeleton[0],s->skeleton[1],s->ng);MCC_TRY(cudaGetLastError());
    MCC_TRY(cudaMemcpy(gx,s->skeleton[0],s->ng*sizeof(DD),cudaMemcpyDeviceToHost));MCC_TRY(cudaMemcpy(gy,s->skeleton[1],s->ng*sizeof(DD),cudaMemcpyDeviceToHost));return 0;
    #undef MCC_TRY
}
extern "C" int mcc_recover(void*ptr,const DD*gx,const DD*gy,DD*ax,DD*ay){
    auto*s=static_cast<MCCContext*>(ptr);if(!s)return cudaErrorInvalidValue;cudaError_t e;
    #define MCC_TRY(call) do{e=(call);if(e!=cudaSuccess)return e;}while(0)
    MCC_TRY(cudaMemcpy(s->skeleton[0],gx,s->ng*sizeof(DD),cudaMemcpyHostToDevice));MCC_TRY(cudaMemcpy(s->skeleton[1],gy,s->ng*sizeof(DD),cudaMemcpyHostToDevice));
    mcc_copy_trace_kernel<<<(s->ng+255)/256,256>>>(s->old_G,s->skeleton[0],s->skeleton[1],s->solution[0],s->solution[1],s->ng);MCC_TRY(cudaGetLastError());
    mcc_recover_kernel<<<s->nc,32>>>(s->cell,s->I,s->G,s->L,s->IG,s->rhs[0],s->rhs[1],s->skeleton[0],s->skeleton[1],s->solution[0],s->solution[1]);MCC_TRY(cudaGetLastError());
    MCC_TRY(cudaMemcpy(ax,s->solution[0],s->nf*sizeof(DD),cudaMemcpyDeviceToHost));MCC_TRY(cudaMemcpy(ay,s->solution[1],s->nf*sizeof(DD),cudaMemcpyDeviceToHost));return 0;
    #undef MCC_TRY
}
