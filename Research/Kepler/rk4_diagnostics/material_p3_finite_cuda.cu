#include <cuda_runtime.h>
#include <chrono>
#include "material_p3_finite_batch.h"
struct Workspace {
    P3FinitePoint *point=nullptr;P3FiniteMap *map=nullptr;P3FiniteOutput *output=nullptr;
    P3FiniteMoment *moment=nullptr;
    size_t points=0,maps=0;cudaEvent_t start,stop;
};
__global__ static void evaluate_kernel(const P3FinitePoint *point,const P3FiniteMap *map,P3FinitePhysics physics,P3FiniteOutput *output,size_t n){
    size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x;if(i<n)output[i]=p3_finite_evaluate(&map[point[i].triangle],&point[i],&physics);
}
__global__ static void moment_kernel(const P3FiniteMap *map,P3FiniteMoment *output,size_t n){
    size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;P3FiniteMoment out={};
    if(map[i].mass[0][0]>0){double constant[12];for(int j=0;j<12;j++)constant[j]=j<10?1:0;out.valid=gfs_material12_moments(&map[i].map,map[i].mass,constant,&out.moment);}
    output[i]=out;
}
extern "C" const char *p3_finite_cuda_error(int code){return cudaGetErrorString((cudaError_t)code);}
extern "C" int p3_finite_cuda_create(const P3FinitePoint *point,size_t n,size_t nm,void **handle,double *seconds){
    auto begin=std::chrono::steady_clock::now();Workspace *w=new Workspace;w->points=n;w->maps=nm;*handle=w;
    cudaError_t e=cudaSetDevice(0);if(e!=cudaSuccess)return e;
    e=cudaMalloc((void**)&w->point,n*sizeof(P3FinitePoint));if(e!=cudaSuccess)return e;
    e=cudaMalloc((void**)&w->map,nm*sizeof(P3FiniteMap));if(e!=cudaSuccess)return e;
    e=cudaMalloc((void**)&w->output,n*sizeof(P3FiniteOutput));if(e!=cudaSuccess)return e;
    e=cudaMalloc((void**)&w->moment,nm*sizeof(P3FiniteMoment));if(e!=cudaSuccess)return e;
    e=cudaEventCreate(&w->start);if(e!=cudaSuccess)return e;e=cudaEventCreate(&w->stop);if(e!=cudaSuccess)return e;
    e=cudaMemcpy(w->point,point,n*sizeof(P3FinitePoint),cudaMemcpyHostToDevice);if(e!=cudaSuccess)return e;
    *seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();return 0;
}
extern "C" int p3_finite_cuda_eval(void *handle,const P3FiniteMap *maps,const P3FinitePhysics *physics,P3FiniteOutput *output,P3FiniteMoment *moment,double *kernel_seconds,double *total_seconds){
    auto begin=std::chrono::steady_clock::now();auto *w=(Workspace*)handle;cudaError_t e=cudaMemcpy(w->map,maps,w->maps*sizeof(P3FiniteMap),cudaMemcpyHostToDevice);if(e!=cudaSuccess)return e;
    e=cudaEventRecord(w->start);if(e!=cudaSuccess)return e;evaluate_kernel<<<(w->points+127)/128,128>>>(w->point,w->map,*physics,w->output,w->points);e=cudaGetLastError();if(e!=cudaSuccess)return e;
    moment_kernel<<<(w->maps+127)/128,128>>>(w->map,w->moment,w->maps);e=cudaGetLastError();if(e!=cudaSuccess)return e;
    e=cudaEventRecord(w->stop);if(e!=cudaSuccess)return e;e=cudaEventSynchronize(w->stop);if(e!=cudaSuccess)return e;float ms;e=cudaEventElapsedTime(&ms,w->start,w->stop);if(e!=cudaSuccess)return e;*kernel_seconds=ms/1000.;
    e=cudaMemcpy(output,w->output,w->points*sizeof(P3FiniteOutput),cudaMemcpyDeviceToHost);if(e!=cudaSuccess)return e;e=cudaMemcpy(moment,w->moment,w->maps*sizeof(P3FiniteMoment),cudaMemcpyDeviceToHost);if(e!=cudaSuccess)return e;*total_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();return 0;
}
extern "C" int p3_finite_cuda_destroy(void *handle){
    auto *w=(Workspace*)handle;if(!w)return 0;cudaError_t first=cudaSuccess,e;
    e=cudaFree(w->point);if(first==cudaSuccess)first=e;e=cudaFree(w->map);if(first==cudaSuccess)first=e;e=cudaFree(w->output);if(first==cudaSuccess)first=e;
    e=cudaFree(w->moment);if(first==cudaSuccess)first=e;
    e=cudaEventDestroy(w->start);if(first==cudaSuccess)first=e;e=cudaEventDestroy(w->stop);if(first==cudaSuccess)first=e;delete w;return first;
}
