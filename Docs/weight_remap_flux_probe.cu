/* Standalone real-device adapter; no fallback to CPU on CUDA errors. */
#include <cuda_runtime.h>
#include <stddef.h>
#include "weight_remap_flux.h"

__global__ void wr_kernel(int mode, size_t n, const double *in, double *out) {
    size_t k = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (k >= n) return;
    const double *a = in + 7*k;
    WrFlux f = wr_pair_flux(mode,a[0],a[1],a[2],a[3],a[4],a[5],a[6]);
    out[4*k] = f.mass; out[4*k+1] = f.px;
    out[4*k+2] = f.py; out[4*k+3] = f.energy;
}

extern "C" int wr_cuda_batch(int mode, size_t n, const double *in, double *out) {
    if (mode != WR_DONOR && mode != WR_CENTERED && mode != WR_KINETIC_THERMAL_DONOR) return -1;
    if (n == 0) return 0;
    double *a = NULL, *b = NULL;
    cudaError_t error = cudaMalloc(&a,7*n*sizeof(double));
    if (error != cudaSuccess) return (int)error;
    error = cudaMalloc(&b,4*n*sizeof(double));
    if (error == cudaSuccess) error = cudaMemcpy(a,in,7*n*sizeof(double),cudaMemcpyHostToDevice);
    if (error == cudaSuccess) {
        wr_kernel<<<(unsigned)((n+255)/256),256>>>(mode,n,a,b);
        error = cudaGetLastError();
    }
    if (error == cudaSuccess) error = cudaDeviceSynchronize();
    if (error == cudaSuccess) error = cudaMemcpy(out,b,4*n*sizeof(double),cudaMemcpyDeviceToHost);
    cudaError_t free_b = cudaFree(b), free_a = cudaFree(a);
    if (error == cudaSuccess) error = free_b;
    if (error == cudaSuccess) error = free_a;
    return (int)error;
}
