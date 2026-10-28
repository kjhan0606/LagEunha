#include <cuda_runtime.h>
#include <stddef.h>
#include "moment_gravity_step.h"
__global__ void ms_kernel(size_t n,const double *in,double *out,int *invalid) {
    size_t k=(size_t)blockIdx.x*blockDim.x+threadIdx.x;
    if(k<n && ms_step(in+19*k,out+17*k)) atomicOr(invalid,1);
}
extern "C" int ms_cuda_batch(size_t n,const double *in,double *out) {
    if(n==0)return 0;
    double *a=NULL,*b=NULL;int *invalid=NULL,host_invalid=0;
    cudaError_t error=cudaMalloc(&a,19*n*sizeof(double));
    if(error!=cudaSuccess)return (int)error;
    error=cudaMalloc(&b,17*n*sizeof(double));
    if(error==cudaSuccess)error=cudaMalloc(&invalid,sizeof(int));
    if(error==cudaSuccess)error=cudaMemset(invalid,0,sizeof(int));
    if(error==cudaSuccess)error=cudaMemcpy(a,in,19*n*sizeof(double),cudaMemcpyHostToDevice);
    if(error==cudaSuccess){ms_kernel<<<(unsigned)((n+255)/256),256>>>(n,a,b,invalid);error=cudaGetLastError();}
    if(error==cudaSuccess)error=cudaDeviceSynchronize();
    if(error==cudaSuccess)error=cudaMemcpy(&host_invalid,invalid,sizeof(int),cudaMemcpyDeviceToHost);
    if(error==cudaSuccess && !host_invalid)error=cudaMemcpy(out,b,17*n*sizeof(double),cudaMemcpyDeviceToHost);
    cudaError_t fi=cudaFree(invalid),fb=cudaFree(b),fa=cudaFree(a);
    if(error==cudaSuccess)error=fi;
    if(error==cudaSuccess)error=fb;
    if(error==cudaSuccess)error=fa;
    return error==cudaSuccess ? (host_invalid ? -1 : 0) : (int)error;
}
