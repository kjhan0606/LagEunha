// DD element algebra and deterministic support gather. No atomic scatter.
// The original binary80 M remains authoritative in the host outer residual.
#include "material_full_gpu_mass.h"
struct MFMassTriangle {int id[ST_MAX],n;};
struct MFMassContext {
    MFMassTriangle*t=nullptr;DD*M=nullptr,*L=nullptr,*partition=nullptr,*input=nullptr,*output=nullptr,*solved=nullptr;
    uint64_t*row=nullptr;uint32_t*support_triangle=nullptr;uint8_t*support_row=nullptr;
    size_t nt=0,nd=0,ns=0;cudaEvent_t begin=nullptr,end=nullptr;
};
__global__ static void mf_mass_matvec(const MFMassTriangle*t,const DD*M,const uint64_t*row,
        const uint32_t*support_triangle,const uint8_t*support_row,const DD*x,DD*y,size_t nd){
    size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=nd)return;DD sum=dd(0);
    for(uint64_t j=row[i];j<row[i+1];j++){size_t s=support_triangle[j];int a=support_row[j];const auto&T=t[s];const DD*block=M+s*ST_MAX*ST_MAX;for(int b=0;b<T.n;b++)sum=add(sum,mul(block[a*ST_MAX+b],x[T.id[b]]));}y[i]=sum;
}
__global__ static void mf_mass_local_inverse(const MFMassTriangle*t,const DD*L,const DD*partition,
        const DD*r,DD*solved,size_t nt){
    size_t s=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(s>=nt)return;const auto&T=t[s];const DD*block=L+s*ST_MAX*ST_MAX;DD x[ST_MAX];for(int a=0;a<T.n;a++)x[a]=mul(partition[T.id[a]],r[T.id[a]]);
    for(int a=0;a<T.n;a++){for(int b=0;b<a;b++)x[a]=sub(x[a],mul(block[a*ST_MAX+b],x[b]));x[a]=divide(x[a],block[a*ST_MAX+a]);}
    for(int a=T.n-1;a>=0;a--){for(int b=a+1;b<T.n;b++)x[a]=sub(x[a],mul(block[b*ST_MAX+a],x[b]));x[a]=divide(x[a],block[a*ST_MAX+a]);}
    for(int a=0;a<T.n;a++)solved[s*ST_MAX+a]=x[a];
}
__global__ static void mf_mass_inverse_gather(const DD*partition,const uint64_t*row,
        const uint32_t*support_triangle,const uint8_t*support_row,const DD*solved,DD*z,size_t nd){
    size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=nd)return;DD sum=dd(0);for(uint64_t j=row[i];j<row[i+1];j++)sum=add(sum,solved[size_t(support_triangle[j])*ST_MAX+support_row[j]]);z[i]=mul(partition[i],sum);
}
extern "C" void mf_mass_destroy(void*ptr){auto*s=static_cast<MFMassContext*>(ptr);if(!s)return;cudaFree(s->t);cudaFree(s->M);cudaFree(s->L);cudaFree(s->partition);cudaFree(s->input);cudaFree(s->output);cudaFree(s->solved);cudaFree(s->row);cudaFree(s->support_triangle);cudaFree(s->support_row);if(s->begin)cudaEventDestroy(s->begin);if(s->end)cudaEventDestroy(s->end);delete s;}
extern "C" int mf_mass_create(const STTriangle*t,size_t nt,size_t nd,const DD*M,const DD*L,
        const DD*partition,const uint64_t*row,const uint32_t*support_triangle,
        const uint8_t*support_row,size_t ns,void**output){
    *output=nullptr;auto*s=new MFMassContext;s->nt=nt;s->nd=nd;s->ns=ns;MFMassTriangle*compact=new MFMassTriangle[nt];for(size_t i=0;i<nt;i++){compact[i]={};compact[i].n=t[i].n;for(int j=0;j<t[i].n;j++)compact[i].id[j]=t[i].id[j];}cudaError_t e=cudaSuccess;
    #define TRY(call) do{if(e==cudaSuccess)e=(call);}while(0)
    #define ALLOC_COPY(field,source,count) TRY(cudaMalloc((void**)&s->field,(count)*sizeof(*s->field)));TRY(cudaMemcpy(s->field,source,(count)*sizeof(*s->field),cudaMemcpyHostToDevice))
    ALLOC_COPY(t,compact,nt);ALLOC_COPY(M,M,nt*ST_MAX*ST_MAX);ALLOC_COPY(L,L,nt*ST_MAX*ST_MAX);ALLOC_COPY(partition,partition,nd);ALLOC_COPY(row,row,nd+1);ALLOC_COPY(support_triangle,support_triangle,ns);ALLOC_COPY(support_row,support_row,ns);
    TRY(cudaMalloc((void**)&s->input,nd*sizeof(DD)));TRY(cudaMalloc((void**)&s->output,nd*sizeof(DD)));TRY(cudaMalloc((void**)&s->solved,nt*ST_MAX*sizeof(DD)));TRY(cudaEventCreate(&s->begin));TRY(cudaEventCreate(&s->end));
    #undef ALLOC_COPY
    #undef TRY
    delete[]compact;if(e!=cudaSuccess){mf_mass_destroy(s);return e;}*output=s;return 0;
}
extern "C" int mf_mass_apply(void*ptr,const DD*input,DD*output,int inverse,double*seconds){
    auto*s=static_cast<MFMassContext*>(ptr);if(!s||inverse<0||inverse>1)return cudaErrorInvalidValue;cudaError_t e;
    #define TRY(call) do{e=(call);if(e!=cudaSuccess)return e;}while(0)
    TRY(cudaMemcpy(s->input,input,s->nd*sizeof(DD),cudaMemcpyHostToDevice));TRY(cudaEventRecord(s->begin));
    if(inverse){mf_mass_local_inverse<<<(s->nt+127)/128,128>>>(s->t,s->L,s->partition,s->input,s->solved,s->nt);TRY(cudaGetLastError());mf_mass_inverse_gather<<<(s->nd+255)/256,256>>>(s->partition,s->row,s->support_triangle,s->support_row,s->solved,s->output,s->nd);}
    else mf_mass_matvec<<<(s->nd+255)/256,256>>>(s->t,s->M,s->row,s->support_triangle,s->support_row,s->input,s->output,s->nd);
    TRY(cudaGetLastError());TRY(cudaEventRecord(s->end));TRY(cudaEventSynchronize(s->end));float ms;TRY(cudaEventElapsedTime(&ms,s->begin,s->end));*seconds=ms/1000.;TRY(cudaMemcpy(output,s->output,s->nd*sizeof(DD),cudaMemcpyDeviceToHost));
    #undef TRY
    return 0;
}
