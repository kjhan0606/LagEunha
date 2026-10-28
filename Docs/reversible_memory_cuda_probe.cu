#include <cuda_runtime.h>
#include "reversible_memory.h"

__global__ void probe_kernel(int n,const double *x,const double *y,int owner,
                            double coverage,double *out,double *work,double *meta,int *status){
    if(threadIdx.x==0 && blockIdx.x==0)
        *status=rm_quadratic_row(n,x,y,owner,coverage,12.,out,work,meta);
}
extern "C" int rm_host(int n,const double *x,const double *y,int owner,double coverage,
                        double *out,double *work,double *meta){
    return rm_quadratic_row(n,x,y,owner,coverage,12.,out,work,meta);
}
extern "C" int rm_device(int n,const double *x,const double *y,int owner,double coverage,
                          double *out,double *work,double *meta){
    (void)work;
    if(n<1 || n>4096)return -9999;
    double *buffer=nullptr;int *code=nullptr;int result=0;cudaError_t err;
    #define CUDA_CHECK(expr) do{err=(expr);if(err!=cudaSuccess){result=-1000-(int)err;goto done;}}while(0)
    CUDA_CHECK(cudaMalloc((void**)&buffer,(10*n+2)*sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&code,sizeof(int)));
    CUDA_CHECK(cudaMemcpy(buffer,x,n*sizeof(double),cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(buffer+n,y,n*sizeof(double),cudaMemcpyHostToDevice));
    probe_kernel<<<1,1>>>(n,buffer,buffer+n,owner,coverage,buffer+2*n,buffer+3*n,buffer+10*n,code);
    CUDA_CHECK(cudaGetLastError());CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(&result,code,sizeof(int),cudaMemcpyDeviceToHost));
    if(result==0){
        CUDA_CHECK(cudaMemcpy(out,buffer+2*n,n*sizeof(double),cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(meta,buffer+10*n,2*sizeof(double),cudaMemcpyDeviceToHost));
    }
done:
    if(code)cudaFree(code);
    if(buffer)cudaFree(buffer);
    return result;
    #undef CUDA_CHECK
}
