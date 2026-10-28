#include <cuda_runtime.h>
#include <stddef.h>
#include "moment_quadrupole_gravity.h"
__global__ void mq_kernel(size_t n,const double *in,double *out,int *invalid) {
    size_t k=(size_t)blockIdx.x*blockDim.x+threadIdx.x;
    if(k<n && mq_source(in+18*k,out+12*k)) atomicOr(invalid,1);
}
extern "C" int mq_cuda_batch(size_t n,const double *in,double *out) {
    if(n==0)return 0;
    double *a=NULL,*b=NULL;int *invalid=NULL,host_invalid=0;
    cudaError_t error=cudaMalloc(&a,18*n*sizeof(double));
    if(error!=cudaSuccess)return (int)error;
    error=cudaMalloc(&b,12*n*sizeof(double));
    if(error==cudaSuccess)error=cudaMalloc(&invalid,sizeof(int));
    if(error==cudaSuccess)error=cudaMemset(invalid,0,sizeof(int));
    if(error==cudaSuccess)error=cudaMemcpy(a,in,18*n*sizeof(double),cudaMemcpyHostToDevice);
    if(error==cudaSuccess){mq_kernel<<<(unsigned)((n+255)/256),256>>>(n,a,b,invalid);error=cudaGetLastError();}
    if(error==cudaSuccess)error=cudaDeviceSynchronize();
    if(error==cudaSuccess)error=cudaMemcpy(&host_invalid,invalid,sizeof(int),cudaMemcpyDeviceToHost);
    if(error==cudaSuccess && !host_invalid)error=cudaMemcpy(out,b,12*n*sizeof(double),cudaMemcpyDeviceToHost);
    cudaError_t fi=cudaFree(invalid),fb=cudaFree(b),fa=cudaFree(a);
    if(error==cudaSuccess)error=fi;
    if(error==cudaSuccess)error=fb;
    if(error==cudaSuccess)error=fa;
    return error==cudaSuccess ? (host_invalid ? -1 : 0) : (int)error;
}
