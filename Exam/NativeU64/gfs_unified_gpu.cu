#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <vector>
#include <cstdint>
#include <algorithm>
extern "C" {
#include "exam_gpu.h"
}
#define CUDA_CHECK(call) do { \
    cudaError_t e=(call); \
    if(e!=cudaSuccess){ \
        fprintf(stderr,"[UNIFIED_GPU] %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); \
        exit(99); \
    } \
} while(0)
#include "gfs_unified_gpu.inc"
extern "C" double gpu_launch_unified_force_kernel(GPUContext *ctx,int n,int nt,
        const GPUPhysicsParams *params)
{
    return gfs_unified_gpu_launch(ctx,n,nt,params);
}
