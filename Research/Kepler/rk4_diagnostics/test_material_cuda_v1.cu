#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <algorithm>
#include "gfs_material_map.h"
static void check(cudaError_t e){if(e!=cudaSuccess){fprintf(stderr,"CUDA %s\n",cudaGetErrorString(e));exit(1);}}
struct Input {GfsMaterialMap map;double N[5],Nx[5],Ny[5],rho,w,sigma[3],gx,gy;};
struct Output {double force[2][5],Udot,Gpower,J,work,torque;int valid;};
__host__ __device__ static Output evaluate(const Input &p){
    Output o={};GfsMaterialPoint point;o.valid=gfs_material_point(&p.map,p.N,p.Nx,p.Ny,&point);
    if(!o.valid)return o;
    o.valid=gfs_material_stress_point(&p.map,p.N,p.Nx,p.Ny,p.rho,p.w,p.sigma,p.gx,p.gy,o.force,&o.Udot,&o.Gpower);
    o.J=point.J;o.work=o.Udot-o.Gpower;
    o.torque=-p.w*p.rho*(point.position[0]*p.gy-point.position[1]*p.gx);
    for(int a=0;a<p.map.n;a++){
        o.work+=p.map.v[0][a]*o.force[0][a]+p.map.v[1][a]*o.force[1][a];
        o.torque+=p.map.x[0][a]*o.force[1][a]-p.map.x[1][a]*o.force[0][a];
    }return o;
}
__global__ static void kernel(const Input *input,Output *out,int n){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)out[i]=evaluate(input[i]);}
int main(){
    const int n=65536;std::vector<Input> inputs(n);std::vector<Output> cpu(n),gpu(n);unsigned long long seed=317;
    auto random=[&](){seed=seed*6364136223846793005ULL+1442695040888963407ULL;return (seed>>11)*0x1.0p-53;};
    for(auto &p:inputs){
        double x=random()-.5,y=random()-.5;p.map.n=5;p.map.x[0][1]=p.map.x[1][2]=1;
        p.map.x[0][3]=.02*(random()-.5);p.map.x[1][4]=.02*(random()-.5);
        for(int k=0;k<2;k++)for(int a=0;a<5;a++)p.map.v[k][a]=random()-.5;
        double N[5]={1,x,y,x*x,x*y},Nx[5]={0,1,0,2*x,y},Ny[5]={0,0,1,0,x};
        for(int a=0;a<5;a++){p.N[a]=N[a];p.Nx[a]=Nx[a];p.Ny[a]=Ny[a];}
        p.rho=.1+2*random();p.w=.001+random();p.sigma[0]=-1.e-6;p.sigma[1]=0;p.sigma[2]=-1.e-6;p.gx=random()-.5;p.gy=random()-.5;
    }
    auto start=std::chrono::steady_clock::now();for(int i=0;i<n;i++)cpu[i]=evaluate(inputs[i]);
    double cpu_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    Input *di;Output *doo;check(cudaMalloc(&di,n*sizeof(Input)));check(cudaMalloc(&doo,n*sizeof(Output)));
    check(cudaMemcpy(di,inputs.data(),n*sizeof(Input),cudaMemcpyHostToDevice));cudaEvent_t a,b;check(cudaEventCreate(&a));check(cudaEventCreate(&b));
    check(cudaEventRecord(a));kernel<<<(n+127)/128,128>>>(di,doo,n);check(cudaGetLastError());check(cudaEventRecord(b));check(cudaEventSynchronize(b));float ms;check(cudaEventElapsedTime(&ms,a,b));
    check(cudaMemcpy(gpu.data(),doo,n*sizeof(Output),cudaMemcpyDeviceToHost));double difference=0,work=0,torque=0,minJ=1;int invalid=0;
    for(int i=0;i<n;i++){
        invalid+=!cpu[i].valid || !gpu[i].valid;minJ=fmin(minJ,gpu[i].J);work=fmax(work,fabs(gpu[i].work));torque=fmax(torque,fabs(gpu[i].torque));
        for(int k=0;k<2;k++)for(int j=0;j<5;j++)difference=fmax(difference,fabs(gpu[i].force[k][j]-cpu[i].force[k][j]));
        difference=fmax(difference,fabs(gpu[i].Udot-cpu[i].Udot));difference=fmax(difference,fabs(gpu[i].Gpower-cpu[i].Gpower));difference=fmax(difference,fabs(gpu[i].J-cpu[i].J));
    }
    cudaDeviceProp prop;check(cudaGetDeviceProperties(&prop,0));check(cudaFree(di));check(cudaFree(doo));
    printf("{\"scope\":\"native CUDA shared material-map/stress helper with65536 manufactured points; kernel timing excludes allocation, transfers, basis/mass assembly and MPI; no Kepler evolution\",\"device\":\"%s\",\"samples\":%d,\"invalid\":%d,\"max_CPU_GPU_difference\":%.17g,\"max_stress_work_identity_error\":%.17g,\"max_internal_stress_torque_error\":%.17g,\"minimum_J\":%.17g,\"CPU_single_thread_seconds\":%.17g,\"GPU_kernel_ms\":%.17g}\n",prop.name,n,invalid,difference,work,torque,minJ,cpu_seconds,(double)ms);
    return invalid || difference>1.e-12 || work>1.e-12 || torque>1.e-12;
}
