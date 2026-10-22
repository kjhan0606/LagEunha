#include <cuda_runtime.h>
#include "material_shared_trace_batch.h"
__global__ static void evaluate(const STTriangle*t,const STCoefficient*c,const STPoint*q,size_t n,STPhysics p,STOutput*out){size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)out[i]=st_evaluate(t,c,q[i],p);}
extern "C" int st_cuda_evaluate(const STTriangle*t,size_t nt,const STCoefficient*c,size_t nc,const STPoint*q,size_t nq,const STPhysics*p,STOutput*out,double*seconds){
    STTriangle*dt=nullptr;STCoefficient*dc=nullptr;STPoint*dq=nullptr;STOutput*dr=nullptr;cudaEvent_t a,b;cudaError_t e=cudaSetDevice(0);if(e!=cudaSuccess)return e;
    #define ST_CUDA(call) do{e=(call);if(e!=cudaSuccess)goto cleanup;}while(0)
    ST_CUDA(cudaMalloc((void**)&dt,nt*sizeof(*dt)));ST_CUDA(cudaMalloc((void**)&dc,nc*sizeof(*dc)));ST_CUDA(cudaMalloc((void**)&dq,nq*sizeof(*dq)));ST_CUDA(cudaMalloc((void**)&dr,nq*sizeof(*dr)));
    ST_CUDA(cudaMemcpy(dt,t,nt*sizeof(*dt),cudaMemcpyHostToDevice));ST_CUDA(cudaMemcpy(dc,c,nc*sizeof(*dc),cudaMemcpyHostToDevice));ST_CUDA(cudaMemcpy(dq,q,nq*sizeof(*dq),cudaMemcpyHostToDevice));
    ST_CUDA(cudaEventCreate(&a));ST_CUDA(cudaEventCreate(&b));ST_CUDA(cudaEventRecord(a));evaluate<<<(nq+127)/128,128>>>(dt,dc,dq,nq,*p,dr);ST_CUDA(cudaGetLastError());ST_CUDA(cudaEventRecord(b));ST_CUDA(cudaEventSynchronize(b));
    {float ms;ST_CUDA(cudaEventElapsedTime(&ms,a,b));*seconds=ms/1000.;}ST_CUDA(cudaEventDestroy(a));ST_CUDA(cudaEventDestroy(b));ST_CUDA(cudaMemcpy(out,dr,nq*sizeof(*dr),cudaMemcpyDeviceToHost));
    cleanup:cudaFree(dt);cudaFree(dc);cudaFree(dq);cudaFree(dr);return e;
    #undef ST_CUDA
}
