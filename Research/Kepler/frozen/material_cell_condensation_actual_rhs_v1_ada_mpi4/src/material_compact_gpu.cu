// Original native physics, DD arithmetic and resident vector kernels preserved.
#include "material_full_initial_cuda.cu"
#include "material_compact_gpu.h"
struct CMContext {
    size_t n=0,no=0,blocks=0,nlist[2]={};
    uint64_t*row=nullptr,*owned=nullptr,*list[2]={};uint32_t*column=nullptr;
    DD*value=nullptr,*inverse_diagonal=nullptr,*v[5]={},*packed=nullptr,*partial=nullptr,*scalar=nullptr;
};
__global__ static void cm_matvec(const uint64_t*row,const uint32_t*col,const DD*value,const DD*x,DD*y,size_t n){
    size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=n)return;DD sum=dd(0);
    for(uint64_t k=row[i];k<row[i+1];k++)sum=add(sum,mul(value[k],x[col[k]]));y[i]=sum;
}
__global__ static void cm_diagonal(const DD*inverse,const DD*r,DD*z,size_t n){size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)z[i]=mul(inverse[i],r[i]);}
extern "C" void cm_destroy(void*ptr){auto*s=static_cast<CMContext*>(ptr);if(!s)return;cudaFree(s->row);cudaFree(s->column);cudaFree(s->value);cudaFree(s->inverse_diagonal);cudaFree(s->owned);for(auto p:s->list)cudaFree(p);for(auto p:s->v)cudaFree(p);cudaFree(s->packed);cudaFree(s->partial);cudaFree(s->scalar);delete s;}
extern "C" int cm_create(size_t n,const uint64_t*row,const uint32_t*column,const DD*value,const DD*inverse,
        const uint64_t*owned,size_t no,const uint64_t*serve,size_t ns,const uint64_t*ask,size_t na,void**out){
    *out=nullptr;if(!n||!no)return cudaErrorInvalidValue;auto*s=new CMContext;s->n=n;s->no=no;s->blocks=(no+255)/256;s->nlist[0]=ns;s->nlist[1]=na;cudaError_t e=cudaSuccess;
    #define TRY(call) do{if(e==cudaSuccess)e=(call);}while(0)
    #define COPY(field,input,count) TRY(cudaMalloc((void**)&s->field,(count)*sizeof(*s->field)));TRY(cudaMemcpy(s->field,input,(count)*sizeof(*s->field),cudaMemcpyHostToDevice))
    COPY(row,row,n+1);COPY(column,column,row[n]);COPY(value,value,row[n]);COPY(inverse_diagonal,inverse,n);COPY(owned,owned,no);
    for(int k=0;k<2;k++)if(s->nlist[k]){TRY(cudaMalloc((void**)&s->list[k],s->nlist[k]*sizeof(uint64_t)));TRY(cudaMemcpy(s->list[k],k?ask:serve,s->nlist[k]*sizeof(uint64_t),cudaMemcpyHostToDevice));}
    for(auto&p:s->v)TRY(cudaMalloc((void**)&p,n*sizeof(DD)));
    TRY(cudaMalloc((void**)&s->packed,std::max(size_t(1),std::max(ns,na))*sizeof(DD)));TRY(cudaMalloc((void**)&s->partial,s->blocks*sizeof(DD)));TRY(cudaMalloc((void**)&s->scalar,sizeof(DD)));
    #undef COPY
    #undef TRY
    if(e!=cudaSuccess){cm_destroy(s);return e;}*out=s;return 0;
}
extern "C" int cm_begin(void*ptr,const DD*rhs){auto*s=static_cast<CMContext*>(ptr);if(!s)return cudaErrorInvalidValue;auto e=cudaMemcpy(s->v[MF_CG_R],rhs,s->n*sizeof(DD),cudaMemcpyHostToDevice);if(e!=cudaSuccess)return e;resident_begin<<<(s->n+255)/256,256>>>(s->v[0],s->v[2],s->v[3],s->v[4],s->n);return cudaGetLastError();}
extern "C" int cm_product(void*ptr,int inverse){auto*s=static_cast<CMContext*>(ptr);if(!s||inverse<0||inverse>1)return cudaErrorInvalidValue;
    if(inverse)cm_diagonal<<<(s->n+255)/256,256>>>(s->inverse_diagonal,s->v[MF_CG_R],s->v[MF_CG_Z],s->n);
    else cm_matvec<<<(s->n+255)/256,256>>>(s->row,s->column,s->value,s->v[MF_CG_P],s->v[MF_CG_AP],s->n);return cudaGetLastError();}
extern "C" int cm_pack(void*ptr,int vector,int list,DD*out){auto*s=static_cast<CMContext*>(ptr);if(!s||vector<0||vector>=5||list<0||list>1)return cudaErrorInvalidValue;size_t n=s->nlist[list];if(!n)return 0;resident_pack<<<(n+255)/256,256>>>(s->v[vector],s->list[list],s->packed,n);auto e=cudaGetLastError();if(e!=cudaSuccess)return e;return cudaMemcpy(out,s->packed,n*sizeof(DD),cudaMemcpyDeviceToHost);}
extern "C" int cm_scatter(void*ptr,int vector,int list,size_t start,size_t n,int plus,const DD*in){auto*s=static_cast<CMContext*>(ptr);if(!s||vector<0||vector>=5||list<0||list>1||start>s->nlist[list]||n>s->nlist[list]-start)return cudaErrorInvalidValue;if(!n)return 0;auto e=cudaMemcpy(s->packed,in,n*sizeof(DD),cudaMemcpyHostToDevice);if(e!=cudaSuccess)return e;resident_scatter<<<(n+255)/256,256>>>(s->v[vector],s->list[list]+start,s->packed,n,plus);return cudaGetLastError();}
extern "C" int cm_update(void*ptr,int mode,DD scale){auto*s=static_cast<CMContext*>(ptr);if(!s||mode<0||mode>2)return cudaErrorInvalidValue;resident_update<<<(s->n+255)/256,256>>>(s->v[0],s->v[1],s->v[2],s->v[3],s->v[4],s->n,mode,scale);return cudaGetLastError();}
extern "C" int cm_dot(void*ptr,int a,int b,DD*out){auto*s=static_cast<CMContext*>(ptr);if(!s||a<0||a>=5||b<0||b>=5)return cudaErrorInvalidValue;resident_dot<<<s->blocks,256>>>(s->v[a],s->v[b],s->owned,s->no,s->partial);auto e=cudaGetLastError();if(e!=cudaSuccess)return e;resident_dot_finish<<<1,256>>>(s->partial,s->blocks,s->scalar);e=cudaGetLastError();if(e!=cudaSuccess)return e;return cudaMemcpy(out,s->scalar,sizeof(DD),cudaMemcpyDeviceToHost);}
extern "C" int cm_finish(void*ptr,DD*out){auto*s=static_cast<CMContext*>(ptr);if(!s)return cudaErrorInvalidValue;return cudaMemcpy(out,s->v[MF_CG_X],s->n*sizeof(DD),cudaMemcpyDeviceToHost);}
