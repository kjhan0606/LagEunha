// Inner PCG vectors stay on device. Only ghost DD pairs and scalars move
// during iterations; the original binary80 operator governs outer residuals.
#include "material_full_mass_resident.h"
#include <algorithm>
struct MFResident {
    MFMassContext*mass=nullptr;DD*v[5]={};uint64_t*owned=nullptr,*list[2]={};
    DD*packed=nullptr,*partial=nullptr,*scalar=nullptr;size_t no=0,nlist[2]={},blocks=0;
};
__global__ static void resident_begin(DD*x,DD*z,DD*p,DD*ap,size_t n){size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)x[i]=z[i]=p[i]=ap[i]=dd(0);}
__global__ static void resident_pack(const DD*x,const uint64_t*index,DD*out,size_t n){size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)out[i]=x[index[i]];}
__global__ static void resident_scatter(DD*x,const uint64_t*index,const DD*in,size_t n,int plus){size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)x[index[i]]=plus?add(x[index[i]],in[i]):in[i];}
__global__ static void resident_update(DD*x,DD*r,DD*z,DD*p,const DD*ap,size_t n,int mode,DD scale){
    size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=n)return;
    if(mode==0)p[i]=z[i];else if(mode==1){x[i]=add(x[i],mul(scale,p[i]));r[i]=sub(r[i],mul(scale,ap[i]));}else p[i]=add(z[i],mul(scale,p[i]));
}
__global__ static void resident_dot(const DD*a,const DD*b,const uint64_t*owned,size_t n,DD*out){
    __shared__ DD sum[256];size_t j=size_t(blockIdx.x)*256+threadIdx.x;int t=threadIdx.x;sum[t]=j<n?mul(a[owned[j]],b[owned[j]]):dd(0);__syncthreads();
    for(int stride=128;stride;stride>>=1){if(t<stride)sum[t]=add(sum[t],sum[t+stride]);__syncthreads();}if(t==0)out[blockIdx.x]=sum[0];
}
__global__ static void resident_dot_finish(const DD*in,size_t n,DD*out){
    __shared__ DD sum[256];int t=threadIdx.x;DD value=dd(0);for(size_t i=t;i<n;i+=256)value=add(value,in[i]);sum[t]=value;__syncthreads();
    for(int stride=128;stride;stride>>=1){if(t<stride)sum[t]=add(sum[t],sum[t+stride]);__syncthreads();}if(t==0)*out=sum[0];
}
extern "C" void mf_resident_destroy(void*ptr){auto*s=static_cast<MFResident*>(ptr);if(!s)return;for(auto p:s->v)cudaFree(p);cudaFree(s->owned);for(auto p:s->list)cudaFree(p);cudaFree(s->packed);cudaFree(s->partial);cudaFree(s->scalar);delete s;}
extern "C" int mf_resident_create(void*mass,const uint64_t*owned,size_t no,const uint64_t*serve,size_t nserve,const uint64_t*ask,size_t nask,void**out){
    *out=nullptr;auto*M=static_cast<MFMassContext*>(mass);if(!M||!no)return cudaErrorInvalidValue;auto*s=new MFResident;s->mass=M;s->no=no;s->nlist[0]=nserve;s->nlist[1]=nask;s->blocks=(no+255)/256;cudaError_t e=cudaSuccess;
    #define TRY(call) do{if(e==cudaSuccess)e=(call);}while(0)
    for(auto&p:s->v)TRY(cudaMalloc((void**)&p,M->nd*sizeof(DD)));
    TRY(cudaMalloc((void**)&s->owned,no*sizeof(uint64_t)));TRY(cudaMemcpy(s->owned,owned,no*sizeof(uint64_t),cudaMemcpyHostToDevice));
    for(int i=0;i<2;i++)if(s->nlist[i]){TRY(cudaMalloc((void**)&s->list[i],s->nlist[i]*sizeof(uint64_t)));TRY(cudaMemcpy(s->list[i],i?ask:serve,s->nlist[i]*sizeof(uint64_t),cudaMemcpyHostToDevice));}
    TRY(cudaMalloc((void**)&s->packed,std::max(size_t(1),std::max(nserve,nask))*sizeof(DD)));TRY(cudaMalloc((void**)&s->partial,s->blocks*sizeof(DD)));TRY(cudaMalloc((void**)&s->scalar,sizeof(DD)));
    #undef TRY
    if(e!=cudaSuccess){mf_resident_destroy(s);return e;}*out=s;return 0;
}
extern "C" int mf_resident_begin(void*ptr,const DD*rhs){auto*s=static_cast<MFResident*>(ptr);if(!s)return cudaErrorInvalidValue;auto e=cudaMemcpy(s->v[MF_CG_R],rhs,s->mass->nd*sizeof(DD),cudaMemcpyHostToDevice);if(e!=cudaSuccess)return e;resident_begin<<<(s->mass->nd+255)/256,256>>>(s->v[0],s->v[2],s->v[3],s->v[4],s->mass->nd);return cudaGetLastError();}
extern "C" int mf_resident_product(void*ptr,int inverse){
    auto*s=static_cast<MFResident*>(ptr);if(!s||inverse<0||inverse>1)return cudaErrorInvalidValue;auto*M=s->mass;
    if(inverse){mf_mass_local_inverse<<<(M->nt+127)/128,128>>>(M->t,M->L,M->partition,s->v[MF_CG_R],M->solved,M->nt);auto e=cudaGetLastError();if(e!=cudaSuccess)return e;mf_mass_inverse_gather<<<(M->nd+255)/256,256>>>(M->partition,M->row,M->support_triangle,M->support_row,M->solved,s->v[MF_CG_Z],M->nd);}
    else mf_mass_matvec<<<(M->nd+255)/256,256>>>(M->t,M->M,M->row,M->support_triangle,M->support_row,s->v[MF_CG_P],s->v[MF_CG_AP],M->nd);
    return cudaGetLastError();
}
extern "C" int mf_resident_pack(void*ptr,int vector,int list,DD*out){auto*s=static_cast<MFResident*>(ptr);if(!s||vector<0||vector>=5||list<0||list>1)return cudaErrorInvalidValue;size_t n=s->nlist[list];if(!n)return 0;resident_pack<<<(n+255)/256,256>>>(s->v[vector],s->list[list],s->packed,n);auto e=cudaGetLastError();if(e!=cudaSuccess)return e;return cudaMemcpy(out,s->packed,n*sizeof(DD),cudaMemcpyDeviceToHost);}
extern "C" int mf_resident_scatter(void*ptr,int vector,int list,size_t start,size_t n,int plus,const DD*in){auto*s=static_cast<MFResident*>(ptr);if(!s||vector<0||vector>=5||list<0||list>1||start>s->nlist[list]||n>s->nlist[list]-start)return cudaErrorInvalidValue;if(!n)return 0;auto e=cudaMemcpy(s->packed,in,n*sizeof(DD),cudaMemcpyHostToDevice);if(e!=cudaSuccess)return e;resident_scatter<<<(n+255)/256,256>>>(s->v[vector],s->list[list]+start,s->packed,n,plus);return cudaGetLastError();}
extern "C" int mf_resident_update(void*ptr,int mode,DD scale){auto*s=static_cast<MFResident*>(ptr);if(!s||mode<0||mode>2)return cudaErrorInvalidValue;resident_update<<<(s->mass->nd+255)/256,256>>>(s->v[0],s->v[1],s->v[2],s->v[3],s->v[4],s->mass->nd,mode,scale);return cudaGetLastError();}
extern "C" int mf_resident_dot(void*ptr,int a,int b,DD*out){auto*s=static_cast<MFResident*>(ptr);if(!s||a<0||a>=5||b<0||b>=5)return cudaErrorInvalidValue;resident_dot<<<s->blocks,256>>>(s->v[a],s->v[b],s->owned,s->no,s->partial);auto e=cudaGetLastError();if(e!=cudaSuccess)return e;resident_dot_finish<<<1,256>>>(s->partial,s->blocks,s->scalar);e=cudaGetLastError();if(e!=cudaSuccess)return e;return cudaMemcpy(out,s->scalar,sizeof(DD),cudaMemcpyDeviceToHost);}
extern "C" int mf_resident_finish(void*ptr,DD*out){auto*s=static_cast<MFResident*>(ptr);if(!s)return cudaErrorInvalidValue;return cudaMemcpy(out,s->v[MF_CG_X],s->mass->nd*sizeof(DD),cudaMemcpyDeviceToHost);}
