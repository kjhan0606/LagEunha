// Full production moment construction reads the same bounded native cache.
// Matrix rows stream through a fixed output window; no full point atlas is kept.
#include "material_full_cubature_gpu.h"
struct MFCubatureAux {double P,U,g;};
__global__ static void cubature_auxiliary(const STTriangle*t,const MFEvolvedCell*cell,MFCachedPoint*cache,const uint32_t*prefix,MFCubatureAux*aux,size_t start,size_t count){
    size_t at=blockIdx.x;if(at>=count)return;const auto&C=cell[t[start+at].cell];
    for(size_t q=prefix[at]+threadIdx.x;q<prefix[at+1];q+=blockDim.x){auto&p=cache[q];auto&a=aux[q];a={};if(p.bad)continue;double j=p.p.J+C.delta;if(!(j>0)||!(p.j0>0)){p.bad=5;continue;}double gamma=mf_config.field.ic.gamma,P0=mf_config.field.ic.pressure,lr=log1p((j-p.j0)/p.j0),z=p.p.position[0]*p.p.position[0]+p.p.position[1]*p.p.position[1]+mf_config.field.ic.eps*mf_config.field.ic.eps;a.P=P0*exp(-gamma*lr);a.U=P0*p.j0/(gamma-1)*exp((1-gamma)*lr);a.g=pow(z,-1.5);if(!(a.P>0)||!(a.U>0)||!isfinite(a.g))p.bad=6;}
}
__device__ static DD cubature_entry(const MFCachedPoint&p,const MFCubatureAux&a,int n,int row){
    if(p.bad)return {nan(""),0};if(row==0)return dd(1);row--;int mass=n*(n+1)/2;
    if(row<mass){int j=0;while(row>=n-j){row-=n-j;j++;}int k=j+row;return mul(mul(dd(p.rho),dd(p.N[j].v)),dd(p.N[k].v));}row-=mass;int det=n*(n-1)/2;
    if(row<det){int j=0;while(row>=n-j-1){row-=n-j-1;j++;}int k=j+1+row;return sub(mul(dd(p.N[j].x),dd(p.N[k].y)),mul(dd(p.N[j].y),dd(p.N[k].x)));}row-=det;
    if(row==0)return dd(a.U);if(row==1)return dd(a.P);row-=2;
    if(row<2*n){int d=row/n,j=row%n;DD derivative=d==0?sub(mul(dd(p.N[j].x),dd(p.p.F[3])),mul(dd(p.N[j].y),dd(p.p.F[2]))):sub(mul(dd(p.N[j].y),dd(p.p.F[0])),mul(dd(p.N[j].x),dd(p.p.F[1])));return neg(mul(dd(a.P),derivative));}row-=2*n;
    if(row<2*n){int d=row/n,j=row%n;return mul(mul(dd(p.rho),dd(p.N[j].v)),dd(-p.p.position[d]*a.g));}
    double z=p.p.position[0]*p.p.position[0]+p.p.position[1]*p.p.position[1]+mf_config.field.ic.eps*mf_config.field.ic.eps;return dd(-p.rho/sqrt(z));
}
__global__ static void cubature_rows_kernel(const STTriangle*t,const MFCachedPoint*cache,const MFCubatureAux*aux,const uint32_t*prefix,size_t start,size_t local,size_t row_start,size_t rows,DD*out){
    size_t n=prefix[local+1]-prefix[local],q=size_t(blockIdx.x)*blockDim.x+threadIdx.x,r=blockIdx.y;if(q<n&&r<rows){size_t at=prefix[local]+q;out[r*n+q]=cubature_entry(cache[at],aux[at],t[start+local].n,row_start+r);}
}
__global__ static void cubature_weights_kernel(const MFCachedPoint*cache,size_t at,size_t n,double*out){size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)out[i]=cache[at+i].weight;}
extern "C" int mf_cubature_reserve(void*ptr,size_t pairs){auto*s=static_cast<Context*>(ptr);if(!s||!s->cache||s->cubature_aux||!pairs)return cudaErrorInvalidValue;auto e=cudaMalloc((void**)&s->cubature_aux,s->cache_capacity*sizeof(MFCubatureAux));if(e==cudaSuccess)e=cudaMalloc((void**)&s->cubature_output,pairs*sizeof(DD));if(e==cudaSuccess)s->cubature_pairs=pairs;return e;}
extern "C" int mf_cubature_tile(void*ptr){auto*s=static_cast<Context*>(ptr);if(!s||!s->cubature_aux||!s->cache_count)return cudaErrorInvalidValue;cubature_auxiliary<<<s->cache_count,128>>>(s->t,s->cell,s->cache,s->cache_prefix,s->cubature_aux,s->cache_start,s->cache_count);return cudaGetLastError();}
extern "C" int mf_cubature_rows(void*ptr,size_t triangle,size_t first,size_t rows,DD*out,double*seconds){
    auto*s=static_cast<Context*>(ptr);if(!s||!s->cubature_aux||triangle<s->cache_start||triangle>=s->cache_start+s->cache_count||!rows)return cudaErrorInvalidValue;size_t local=triangle-s->cache_start;uint32_t offsets[2];auto e=cudaMemcpy(offsets,s->cache_prefix+local,2*sizeof(uint32_t),cudaMemcpyDeviceToHost);if(e!=cudaSuccess)return e;size_t n=offsets[1]-offsets[0];if(!n||rows*n>s->cubature_pairs)return cudaErrorInvalidValue;
    #define TRY(call) do{e=(call);if(e!=cudaSuccess)return e;}while(0)
    TRY(cudaEventRecord(s->begin));cubature_rows_kernel<<<dim3((n+127)/128,rows),128>>>(s->t,s->cache,s->cubature_aux,s->cache_prefix,s->cache_start,local,first,rows,s->cubature_output);TRY(cudaGetLastError());TRY(cudaEventRecord(s->end));TRY(cudaEventSynchronize(s->end));float ms;TRY(cudaEventElapsedTime(&ms,s->begin,s->end));*seconds=ms/1000.;TRY(cudaMemcpy(out,s->cubature_output,rows*n*sizeof(DD),cudaMemcpyDeviceToHost));
    #undef TRY
    return 0;
}
extern "C" int mf_cubature_weights(void*ptr,size_t triangle,double*out,size_t n){auto*s=static_cast<Context*>(ptr);if(!s||!s->cubature_aux||triangle<s->cache_start||triangle>=s->cache_start+s->cache_count||n*sizeof(double)>s->cubature_pairs*sizeof(DD))return cudaErrorInvalidValue;uint32_t offsets[2];auto e=cudaMemcpy(offsets,s->cache_prefix+triangle-s->cache_start,2*sizeof(uint32_t),cudaMemcpyDeviceToHost);if(e!=cudaSuccess)return e;if(n!=offsets[1]-offsets[0])return cudaErrorInvalidValue;auto*d=reinterpret_cast<double*>(s->cubature_output);cubature_weights_kernel<<<(n+255)/256,256>>>(s->cache,offsets[0],n,d);e=cudaGetLastError();if(e!=cudaSuccess)return e;return cudaMemcpy(out,d,n*sizeof(double),cudaMemcpyDeviceToHost);}
