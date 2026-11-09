#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>
#include "../gfs_material_node2d.h"
__global__ void solve(double *out,int *ok){
 int k=blockIdx.x*blockDim.x+threadIdx.x;if(k>=256)return;
 GfsMaterialHalfEdge e[4]={{{1,0},.5,1,1,{0,0}},{{-1,0},.5,.125,.1,{0,0}},{{0,1},.5,1,1,{0,0}},{{0,-1},.5,1,1,{0,0}}};
 double boost=.01*k;for(int j=0;j<4;j++)e[j].u[0]+=boost;
 double v[2]={boost,0},r;int it;ok[k]=gfs_material_node_solve(4,e,1.4,1,1,v,&it,&r);out[k]=v[0]-boost;
}
int main(void){
 double *d,h[256];int *d_ok,ok[256];
 if(cudaMalloc((void**)&d,sizeof(h))!=cudaSuccess||cudaMalloc((void**)&d_ok,sizeof(ok))!=cudaSuccess)return 1;
 solve<<<2,128>>>(d,d_ok);if(cudaGetLastError()!=cudaSuccess)return 2;
 if(cudaMemcpy(h,d,sizeof(h),cudaMemcpyDeviceToHost)!=cudaSuccess||cudaMemcpy(ok,d_ok,sizeof(ok),cudaMemcpyDeviceToHost)!=cudaSuccess)return 3;
 double err=0;for(int i=0;i<256;i++){if(!ok[i])return 4;err=fmax(err,fabs(h[i]-.92745262004895002));}
 cudaFree(d);cudaFree(d_ok);printf("PASS 256 device exact pressure node solves: maximum velocity difference %.17g\n",err);return err>1e-12;
}
