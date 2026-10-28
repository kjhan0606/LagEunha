#include <cstdio>
#include <cstring>
#include <cmath>
#include <omp.h>
#include <cuda_runtime.h>
#include "tcf_full_face_probe.h"
#define CK(x) do{cudaError_t e=(x);if(e!=cudaSuccess){fprintf(stderr,"%s\n",cudaGetErrorString(e));return 2;}}while(0)
__global__ void kernel(const double*a,double*o,int n){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)tcf_full_face_probe(a+13*i,o+6*i);}
int main(){FILE*f=fopen("inputs.bin","rb");if(!f)return 3;fseek(f,0,SEEK_END);long size=ftell(f);rewind(f);int n=size/(13*sizeof(double));double*a=new double[13*n],*one=new double[6*n],*four=new double[6*n],*gpu=new double[6*n];if(fread(a,sizeof(double),13*n,f)!=(size_t)(13*n))return 4;fclose(f);int teams[2];omp_set_dynamic(0);
for(int p=0;p<2;p++){omp_set_num_threads(p?4:1);double*o=p?four:one;
#pragma omp parallel
{
#pragma omp single
teams[p]=omp_get_num_threads();
#pragma omp for
for(int i=0;i<n;i++)tcf_full_face_probe(a+13*i,o+6*i);}}
double*da,*dd;CK(cudaMalloc(&da,size));CK(cudaMalloc(&dd,6*n*sizeof(double)));CK(cudaMemcpy(da,a,size,cudaMemcpyHostToDevice));kernel<<<(n+127)/128,128>>>(da,dd,n);CK(cudaGetLastError());CK(cudaDeviceSynchronize());CK(cudaMemcpy(gpu,dd,6*n*sizeof(double),cudaMemcpyDeviceToHost));double delta=0,scaled=0;int finite=1;for(int i=0;i<6*n;i++){finite=finite&&std::isfinite(one[i])&&std::isfinite(gpu[i]);delta=fmax(delta,fabs(one[i]-gpu[i]));scaled=fmax(scaled,fabs(one[i]-gpu[i])/fmax(1.,fabs(one[i])));}f=fopen("cpu.bin","wb");fwrite(one,sizeof(double),6*n,f);fclose(f);f=fopen("gpu.bin","wb");fwrite(gpu,sizeof(double),6*n,f);fclose(f);printf("{\"rows\":%d,\"teams\":[%d,%d],\"finite\":%d,\"omp_bitwise\":%s,\"cuda_bitwise\":%s,\"max_abs\":%.17g,\"max_scaled\":%.17g}\n",n,teams[0],teams[1],finite,memcmp(one,four,6*n*sizeof(double))?"false":"true",memcmp(one,gpu,6*n*sizeof(double))?"false":"true",delta,scaled);return !(finite&&teams[0]==1&&teams[1]==4&&scaled<1e-12);}
