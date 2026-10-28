#include <cuda_runtime.h>
#include <stdio.h>
__global__ void geometry_execution_marker(double*out){int i=threadIdx.x;if(i<4)out[i]=42.+i;}
int main(){int runtime=0,driver=0;cudaRuntimeGetVersion(&runtime);cudaDriverGetVersion(&driver);double*device=nullptr;double host[4]={-1,-1,-1,-1};
 int a=cudaMalloc(&device,sizeof(host)),b=cudaMemset(device,0xA5,sizeof(host));geometry_execution_marker<<<1,4>>>(device);int c=cudaGetLastError(),d=cudaDeviceSynchronize(),e=cudaMemcpy(host,device,sizeof(host),cudaMemcpyDeviceToHost);
 printf("runtime=%d driver=%d statuses=%d,%d,%d,%d,%d output=%.17g,%.17g,%.17g,%.17g\n",runtime,driver,a,b,c,d,e,host[0],host[1],host[2],host[3]);if(device)cudaFree(device);return a||b||c||d||e||host[0]!=42.||host[3]!=45.;}
