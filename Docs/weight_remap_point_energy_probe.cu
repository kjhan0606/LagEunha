#include <cuda_runtime.h>
#include <stddef.h>
#include "weight_remap_point_energy.h"
__global__ void wp_kernel(size_t n,const double *in,double *out) {
    size_t k=(size_t)blockIdx.x*blockDim.x+threadIdx.x;
    if(k<n) wp_pair(in+9*k,out+6*k);
}
extern "C" int wp_cuda_batch(size_t n,const double *in,double *out) {
    if(!n) return 0;
    double *a=NULL,*b=NULL;
    cudaError_t e=cudaMalloc(&a,9*n*sizeof(double));
    if(e==cudaSuccess) e=cudaMalloc(&b,6*n*sizeof(double));
    if(e==cudaSuccess) e=cudaMemcpy(a,in,9*n*sizeof(double),cudaMemcpyHostToDevice);
    if(e==cudaSuccess) {
        wp_kernel<<<(unsigned)((n+255)/256),256>>>(n,a,b);
        e=cudaGetLastError();
    }
    if(e==cudaSuccess) e=cudaDeviceSynchronize();
    if(e==cudaSuccess) e=cudaMemcpy(out,b,6*n*sizeof(double),cudaMemcpyDeviceToHost);
    cudaError_t eb=cudaFree(b),ea=cudaFree(a);
    if(e==cudaSuccess) e=eb;
    if(e==cudaSuccess) e=ea;
    return (int)e;
}
