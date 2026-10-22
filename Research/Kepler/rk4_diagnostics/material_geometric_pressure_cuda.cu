#include <cuda_runtime.h>
#include "material_geometric_pressure_batch.h"
__global__ static void evaluate(const STTriangle*t,const STCoefficient*c,const STPoint*q,size_t n,STPhysics p,const GPInput*input,GPOutput*out){size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)out[i]=gp_evaluate(t,c,q[i],p,input[i]);}
extern "C" int gp_cuda_evaluate(const STTriangle*t,size_t nt,const STCoefficient*c,size_t nc,const STPoint*q,size_t nq,const STPhysics*p,const GPInput*input,GPOutput*out,double*seconds){
    STTriangle*dt=nullptr;STCoefficient*dc=nullptr;STPoint*dq=nullptr;GPInput*di=nullptr;GPOutput*dr=nullptr;cudaEvent_t a=nullptr,b=nullptr;cudaError_t error=cudaSetDevice(0);if(error!=cudaSuccess)return error;
    #define GP_CUDA(call) do{error=(call);if(error!=cudaSuccess)goto cleanup;}while(0)
    GP_CUDA(cudaMalloc((void**)&dt,nt*sizeof(*dt)));GP_CUDA(cudaMalloc((void**)&dc,nc*sizeof(*dc)));GP_CUDA(cudaMalloc((void**)&dq,nq*sizeof(*dq)));GP_CUDA(cudaMalloc((void**)&di,nq*sizeof(*di)));GP_CUDA(cudaMalloc((void**)&dr,nq*sizeof(*dr)));
    GP_CUDA(cudaMemcpy(dt,t,nt*sizeof(*dt),cudaMemcpyHostToDevice));GP_CUDA(cudaMemcpy(dc,c,nc*sizeof(*dc),cudaMemcpyHostToDevice));GP_CUDA(cudaMemcpy(dq,q,nq*sizeof(*dq),cudaMemcpyHostToDevice));GP_CUDA(cudaMemcpy(di,input,nq*sizeof(*di),cudaMemcpyHostToDevice));
    GP_CUDA(cudaEventCreate(&a));GP_CUDA(cudaEventCreate(&b));GP_CUDA(cudaEventRecord(a));evaluate<<<(nq+127)/128,128>>>(dt,dc,dq,nq,*p,di,dr);GP_CUDA(cudaGetLastError());GP_CUDA(cudaEventRecord(b));GP_CUDA(cudaEventSynchronize(b));{float ms;GP_CUDA(cudaEventElapsedTime(&ms,a,b));*seconds=ms/1000.;}GP_CUDA(cudaMemcpy(out,dr,nq*sizeof(*dr),cudaMemcpyDeviceToHost));
    cleanup:if(a)cudaEventDestroy(a);if(b)cudaEventDestroy(b);cudaFree(dt);cudaFree(dc);cudaFree(dq);cudaFree(di);cudaFree(dr);return error;
    #undef GP_CUDA
}
