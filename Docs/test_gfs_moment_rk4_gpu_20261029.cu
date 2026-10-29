#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "Exam/gfs_moment_rk4.h"
static void ck(cudaError_t e){if(e!=cudaSuccess){fprintf(stderr,"CUDA: %s\n",cudaGetErrorString(e));exit(1);}}
__global__ void run(GfsMomentState *out,int *ok){
    GfsMomentState state;gfs_moment_zero(&state);state.a[GFM_M]=1;
    state.a[GFM_QXX]=state.a[GFM_QYY]=1./12.;state.a[GFM_E]=10;state.a[GFM_PX]=1;
    for(int i=0;i<10;i++){
        GfsMomentRK4 r={};if(!gfs_moment_rk4_begin(&r,&state,.1,1)){*ok=0;return;}
        for(int j=0;j<4;j++){
            GfsMomentState rhs;gfs_moment_zero(&rhs);
            rhs.a[GFM_PX]=-r.stage.a[GFM_PY];rhs.a[GFM_PY]=r.stage.a[GFM_PX];
            if(!gfs_moment_rk4_add_rhs(&r,&rhs)){*ok=0;return;}
        }
        if(!gfs_moment_rk4_candidate(&r,&state)){*ok=0;return;}
    }
    *out=state;*ok=1;
}
int main(){
    GfsMomentState *dstate,state;int *dok,ok;
    ck(cudaMalloc(&dstate,sizeof(state)));ck(cudaMalloc(&dok,sizeof(ok)));
    run<<<1,1>>>(dstate,dok);ck(cudaGetLastError());
    ck(cudaMemcpy(&state,dstate,sizeof(state),cudaMemcpyDeviceToHost));
    ck(cudaMemcpy(&ok,dok,sizeof(ok),cudaMemcpyDeviceToHost));
    if(ok!=1 || std::hypot(state.a[GFM_PX]-std::cos(1.),state.a[GFM_PY]-std::sin(1.))>9e-7 || state.a[GFM_M]!=1)return 2;
    ck(cudaFree(dstate));ck(cudaFree(dok));puts("GPU moment RK4 stages, fixed mass and analytic oscillator pass");
}
