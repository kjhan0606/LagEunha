#include <cuda_runtime.h>
#include <stddef.h>
#include "weight_remap_moments.h"

__global__ void wm_kernel(size_t n,const double *in,double *out) {
    size_t k=(size_t)blockIdx.x*blockDim.x+threadIdx.x;
    if(k<n) wm_atom(in+7*k,out+17*k);
}
extern "C" int wm_cuda_batch(size_t n,const double *in,double *out) {
    if(n==0) return 0;
    double *a=NULL,*b=NULL;
    cudaError_t error=cudaMalloc(&a,7*n*sizeof(double));
    if(error!=cudaSuccess) return (int)error;
    error=cudaMalloc(&b,17*n*sizeof(double));
    if(error==cudaSuccess) error=cudaMemcpy(a,in,7*n*sizeof(double),cudaMemcpyHostToDevice);
    if(error==cudaSuccess) {
        wm_kernel<<<(unsigned)((n+255)/256),256>>>(n,a,b);
        error=cudaGetLastError();
    }
    if(error==cudaSuccess) error=cudaDeviceSynchronize();
    if(error==cudaSuccess) error=cudaMemcpy(out,b,17*n*sizeof(double),cudaMemcpyDeviceToHost);
    cudaError_t fb=cudaFree(b),fa=cudaFree(a);
    if(error==cudaSuccess) error=fb;
    if(error==cudaSuccess) error=fa;
    return (int)error;
}

__global__ void wm_lift_kernel(size_t n,const double *in,double *out,int *invalid) {
    size_t k=(size_t)blockIdx.x*blockDim.x+threadIdx.x;
    if(k<n && wm_lift_atom(in+16*k,out+17*k)) atomicOr(invalid,1);
}
extern "C" int wm_lift_cuda_batch(size_t n,const double *in,double *out) {
    if(n==0) return 0;
    double *a=NULL,*b=NULL;int *invalid=NULL,host_invalid=0;
    cudaError_t error=cudaMalloc(&a,16*n*sizeof(double));
    if(error!=cudaSuccess) return (int)error;
    error=cudaMalloc(&b,17*n*sizeof(double));
    if(error==cudaSuccess) error=cudaMalloc(&invalid,sizeof(int));
    if(error==cudaSuccess) error=cudaMemset(invalid,0,sizeof(int));
    if(error==cudaSuccess) error=cudaMemcpy(a,in,16*n*sizeof(double),cudaMemcpyHostToDevice);
    if(error==cudaSuccess) {
        wm_lift_kernel<<<(unsigned)((n+255)/256),256>>>(n,a,b,invalid);
        error=cudaGetLastError();
    }
    if(error==cudaSuccess) error=cudaDeviceSynchronize();
    if(error==cudaSuccess) error=cudaMemcpy(&host_invalid,invalid,sizeof(int),cudaMemcpyDeviceToHost);
    if(error==cudaSuccess && !host_invalid) error=cudaMemcpy(out,b,17*n*sizeof(double),cudaMemcpyDeviceToHost);
    cudaError_t fi=cudaFree(invalid),fb=cudaFree(b),fa=cudaFree(a);
    if(error==cudaSuccess) error=fi;
    if(error==cudaSuccess) error=fb;
    if(error==cudaSuccess) error=fa;
    if(error!=cudaSuccess) return (int)error;
    return host_invalid ? -1 : 0;
}
