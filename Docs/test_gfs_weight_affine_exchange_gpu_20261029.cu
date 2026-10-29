#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "Exam/gfs_weight_mass_exchange.h"
static void ck(cudaError_t e){if(e!=cudaSuccess){fprintf(stderr,"CUDA: %s\n",cudaGetErrorString(e));exit(1);}}
__global__ void kernel(GfsMomentState *a,GfsMomentState *b,int *ok){
    gfs_moment_zero(a);gfs_moment_zero(b);double G[]={.2,-.4,.6,.1};
    *ok=gfs_weight_exchange_segment_affine(a,b,-.3,.2,1.9,-1.4,-.2,.3,1.7,-.4,.12,-.04,2,1.5,.3,-.7,G,4);
}
int main(){
    GfsMomentState *da,*db,a,b,refa,refb;int *dok,ok;
    ck(cudaMalloc(&da,sizeof(a)));ck(cudaMalloc(&db,sizeof(b)));ck(cudaMalloc(&dok,sizeof(int)));
    kernel<<<1,1>>>(da,db,dok);ck(cudaGetLastError());
    ck(cudaMemcpy(&a,da,sizeof(a),cudaMemcpyDeviceToHost));ck(cudaMemcpy(&b,db,sizeof(b),cudaMemcpyDeviceToHost));ck(cudaMemcpy(&ok,dok,sizeof(ok),cudaMemcpyDeviceToHost));
    double G[]={.2,-.4,.6,.1};gfs_moment_zero(&refa);gfs_moment_zero(&refb);
    if(!ok || !gfs_weight_exchange_segment_affine(&refa,&refb,-.3,.2,1.9,-1.4,-.2,.3,1.7,-.4,.12,-.04,2,1.5,.3,-.7,G,4))return 2;
    for(int k=0;k<GFM_NSTATE;k++)if(!std::isfinite(a.a[k])||!std::isfinite(b.a[k])||
        fabs(a.a[k]-refa.a[k])>1e-12 || fabs(b.a[k]-refb.a[k])>1e-12)return 3;
    ck(cudaFree(da));ck(cudaFree(db));ck(cudaFree(dok));puts("GPU affine weight-exchange moments match host reference");
}
