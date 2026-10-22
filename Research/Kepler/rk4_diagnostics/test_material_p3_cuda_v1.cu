/* Actual native P3-query device arithmetic; not a CUDA hydro/runtime benchmark. */
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include "gfs_material_p3.h"
#include "gfs_material_p3_radial_input.h"
struct Query {GfsMaterialP3Triangle triangle;double x,y,eps,reference[6];};
struct Output {double N[12],Nx[12],Ny[12];int valid;};
static void need(bool ok,const char *message){if(!ok){fprintf(stderr,"FAIL %s\n",message);exit(1);}}
static void cuda_need(cudaError_t error){if(error!=cudaSuccess){fprintf(stderr,"CUDA %s\n",cudaGetErrorString(error));exit(1);}}
__host__ __device__ static Output evaluate(Query q){
    Output out={};double one[2]={1,1};
    for(int a=0;a<10;a++){double x[2],D[4];gfs_material_p3_node(&q.triangle,a,x);if(!gfs_material_p3_radial_remainder(q.triangle.geometry.x[0],x[0],x[1],q.eps,q.triangle.remainder[a],D))return out;}
    double r[2],Dr[4];if(!gfs_material_p3_radial_remainder(q.triangle.geometry.x[0],q.x,q.y,q.eps,r,Dr))return out;gfs_material_p3_shape_remainder(&q.triangle,q.x,q.y,r,Dr,one,out.N,out.Nx,out.Ny);out.valid=1;return out;
}
__global__ static void kernel(const Query *query,Output *out,int n){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)out[i]=evaluate(query[i]);}
int main(int argc,char **argv){
    need(argc==3,"two frozen actual-query CSVs");std::vector<Query> query;
    for(int file=1;file<argc;file++){
        std::ifstream source(argv[file]);need(source.good(),"actual native P3 query input");std::string line;std::getline(source,line);
        while(std::getline(source,line)){
            if(line.empty())continue;std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);int triangle;Query q={};row>>triangle;for(int a=0;a<3;a++)row>>q.triangle.geometry.x[a][0]>>q.triangle.geometry.x[a][1];row>>q.x>>q.y>>q.eps;for(int k=0;k<6;k++)row>>q.reference[k];need(!row.fail(),"actual query CSV parse");need(gfs_material_triangle_init(&q.triangle.geometry),"positive actual P3 query triangle");query.push_back(q);
        }
    }
    int n=query.size();need(n>0,"nonempty actual query set");std::vector<Output> cpu(n),gpu(n);for(int i=0;i<n;i++)cpu[i]=evaluate(query[i]);Query *device_query;Output *device_output;cuda_need(cudaSetDevice(0));cudaDeviceProp property;cuda_need(cudaGetDeviceProperties(&property,0));cuda_need(cudaMalloc((void**)&device_query,n*sizeof(Query)));cuda_need(cudaMalloc((void**)&device_output,n*sizeof(Output)));cuda_need(cudaMemcpy(device_query,query.data(),n*sizeof(Query),cudaMemcpyHostToDevice));kernel<<<(n+63)/64,64>>>(device_query,device_output,n);cuda_need(cudaGetLastError());cuda_need(cudaDeviceSynchronize());cuda_need(cudaMemcpy(gpu.data(),device_output,n*sizeof(Output),cudaMemcpyDeviceToHost));cuda_need(cudaFree(device_query));cuda_need(cudaFree(device_output));
    double polynomial=0,residual=0,packet=0;int invalid=0;
    for(int i=0;i<n;i++){
        if(!cpu[i].valid || !gpu[i].valid){invalid++;continue;}double h=sqrt(2*query[i].triangle.geometry.area),scale=std::max(fabs(cpu[i].N[10]),fabs(cpu[i].N[11]));for(int j=10;j<12;j++)scale=std::max({scale,h*fabs(cpu[i].Nx[j]),h*fabs(cpu[i].Ny[j])});need(scale>0,"actual residual scale");
        for(int j=0;j<10;j++)polynomial=std::max({polynomial,fabs(cpu[i].N[j]-gpu[i].N[j]),h*fabs(cpu[i].Nx[j]-gpu[i].Nx[j]),h*fabs(cpu[i].Ny[j]-gpu[i].Ny[j])});
        for(int j=10;j<12;j++)residual=std::max({residual,fabs(cpu[i].N[j]-gpu[i].N[j])/scale,h*fabs(cpu[i].Nx[j]-gpu[i].Nx[j])/scale,h*fabs(cpu[i].Ny[j]-gpu[i].Ny[j])/scale});
        double values[6]={cpu[i].N[10],cpu[i].N[11],cpu[i].Nx[10],cpu[i].Ny[10],cpu[i].Nx[11],cpu[i].Ny[11]};for(int k=0;k<6;k++)packet=std::max(packet,(k<2?1:h)*fabs(values[k]-query[i].reference[k])/scale);
    }
    printf("{\"scope\":\"actual native P3 input residual and polynomial helper on CUDA; no device mass/Schur solve, hydro evolution or performance claim\",\"device\":\"%s\",\"queries\":%d,\"invalid\":%d,\"CPU_GPU_polynomial_scaled_difference\":%.17g,\"CPU_GPU_residual_normalized_difference\":%.17g,\"CPU_native_packet_normalized_difference\":%.17g}\n",property.name,n,invalid,polynomial,residual,packet);
    need(invalid==0 && polynomial<1.e-10 && residual<1.e-9 && packet<1.e-9,"actual P3 CPU/CUDA arithmetic gate");return 0;
}
