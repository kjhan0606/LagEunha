#include <cuda_runtime.h>
#include "reversible_memory.h"

__global__ void gradient_kernel(int n,const double*x,const double*y,double coverage,
    double*out,double*work,double*meta,int*status){
    if(threadIdx.x==0 && blockIdx.x==0)
        *status=rm_quadratic_row_gradient(n,x,y,0,coverage,12.,out,work,meta,out+n,out+2*n);
}
extern "C" int gradient_host(int n,const double*x,const double*y,double coverage,
    double*out,double*work,double*meta){
    return rm_quadratic_row_gradient(n,x,y,0,coverage,12.,out,work,meta,out+n,out+2*n);
}
extern "C" int gradient_device(int n,const double*x,const double*y,double coverage,
    double*out,double*work,double*meta){
    (void)work;
    if(n<1 || n>4096)return -9999;
    double *buffer=nullptr; int *code=nullptr;int result=0;cudaError_t err;
    #define CHECK(expr) do{err=(expr);if(err!=cudaSuccess){result=-1000-(int)err;goto done;}}while(0)
    CHECK(cudaMalloc((void**)&buffer,(12*n+2)*sizeof(double)));
    CHECK(cudaMalloc((void**)&code,sizeof(int)));
    CHECK(cudaMemcpy(buffer,x,n*sizeof(double),cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(buffer+n,y,n*sizeof(double),cudaMemcpyHostToDevice));
    gradient_kernel<<<1,1>>>(n,buffer,buffer+n,coverage,buffer+2*n,buffer+5*n,buffer+12*n,code);
    CHECK(cudaGetLastError());CHECK(cudaDeviceSynchronize());
    CHECK(cudaMemcpy(&result,code,sizeof(int),cudaMemcpyDeviceToHost));
    if(result==0){
        CHECK(cudaMemcpy(out,buffer+2*n,3*n*sizeof(double),cudaMemcpyDeviceToHost));
        CHECK(cudaMemcpy(meta,buffer+12*n,2*sizeof(double),cudaMemcpyDeviceToHost));
    }
done:
    if(code)cudaFree(code);if(buffer)cudaFree(buffer);
    return result;
    #undef CHECK
}
