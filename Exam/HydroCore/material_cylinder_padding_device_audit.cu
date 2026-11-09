#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cuda_runtime.h>
#include "gfs_material_cylinder_boundary2d.h"
struct Input {
 CylinderPaddingOwner map;double x[2],v[2],gv[2],W,q;
};
static __host__ __device__ int apply(const Input&a,GfsMaterialCylinderImage2D*out){
 double center[2]={4,4};return gfs_material_cylinder_image2d(&a.map,1,9,a.x,a.v,
  a.W,a.q,a.gv,center,1,8,8,out);
}
__global__ void check_kernel(const Input*in,GfsMaterialCylinderImage2D*out,int*status,int n){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)status[i]=apply(in[i],out+i);
}
#define CUDA(x) do {cudaError_t e=(x);if(e!=cudaSuccess){fprintf(stderr,"CUDA: %s\n",cudaGetErrorString(e));return 1;}}while(0)
int main(){
 const int n=4096;Input*in=(Input*)calloc(n,sizeof(Input));
 auto*cpu=(GfsMaterialCylinderImage2D*)calloc(n,sizeof(GfsMaterialCylinderImage2D));
 auto*gpu=(GfsMaterialCylinderImage2D*)calloc(n,sizeof(GfsMaterialCylinderImage2D));
 int*ok=(int*)calloc(n,sizeof(int));if(!in||!cpu||!gpu||!ok)return 1;
 for(int i=0;i<n;i++){
  int kind=1+i%4,shift=kind==1?0:kind==2?((i/4)%2?1:-1):i%3-1;
  in[i].map={9,0,kind,shift};double theta=.017*i,r=1+((i%17==0)?1e-12:.1+.001*(i%20));
  in[i].x[0]=kind==1?4+r*cos(theta):kind==4?7.9:.1;
  in[i].x[1]=kind==1?4+r*sin(theta):3.9;
  in[i].v[0]=.3+.001*i;in[i].v[1]=-.2+.0001*i;
  in[i].gv[0]=kind==3?1:in[i].v[0];in[i].gv[1]=kind==3?0:in[i].v[1];
  in[i].W=.01+.0001*i;in[i].q=.02-.00001*i;
  if(!apply(in[i],cpu+i))return 2;
 }
 Input*di;GfsMaterialCylinderImage2D*dg;int*ds;
 CUDA(cudaMalloc(&di,n*sizeof(Input)));CUDA(cudaMalloc(&dg,n*sizeof(*gpu)));CUDA(cudaMalloc(&ds,n*sizeof(int)));
 CUDA(cudaMemcpy(di,in,n*sizeof(Input),cudaMemcpyHostToDevice));
 check_kernel<<<(n+127)/128,128>>>(di,dg,ds,n);CUDA(cudaGetLastError());
 CUDA(cudaMemcpy(gpu,dg,n*sizeof(*gpu),cudaMemcpyDeviceToHost));CUDA(cudaMemcpy(ok,ds,n*sizeof(int),cudaMemcpyDeviceToHost));
 double error=0,normal_error=0;int periodic_positive=0,periodic_negative=0,near_count=0,nonzero_low=0;
 for(int i=0;i<n;i++){
  if(!ok[i]||cpu[i].weight_owner!=gpu[i].weight_owner||cpu[i].wall_owner!=gpu[i].wall_owner||cpu[i].periodic_owner!=gpu[i].periodic_owner)return 3;
  for(int k=0;k<2;k++){
   error=fmax(error,fabs(cpu[i].position[k]-gpu[i].position[k]));error=fmax(error,fabs(cpu[i].low[k]-gpu[i].low[k]));
   error=fmax(error,fabs(cpu[i].mesh_velocity[k]-gpu[i].mesh_velocity[k]));error=fmax(error,fabs(cpu[i].fluid_velocity[k]-gpu[i].fluid_velocity[k]));
  }
  error=fmax(error,fabs(cpu[i].weight-gpu[i].weight));error=fmax(error,fabs(cpu[i].weight_rate-gpu[i].weight_rate));
  if(in[i].map.kind==CYL_IMAGE_PERIODIC){
   periodic_positive+=in[i].map.y_shift>0;periodic_negative+=in[i].map.y_shift<0;
  }
  if(in[i].map.kind==CYL_IMAGE_RADIAL&&i%17==0){
   near_count++;nonzero_low+=gpu[i].low[0]!=0||gpu[i].low[1]!=0;
   double dx=(gpu[i].position[0]-in[i].x[0])+gpu[i].low[0];
   double dy=(gpu[i].position[1]-in[i].x[1])+gpu[i].low[1];
   double d=hypot(dx,dy),r=hypot(in[i].x[0]-4,in[i].x[1]-4);
   if(!(d>0))return 5;
   normal_error=fmax(normal_error,hypot(dx/d+(in[i].x[0]-4)/r,dy/d+(in[i].x[1]-4)/r));
  }
 }
 CUDA(cudaFree(di));CUDA(cudaFree(dg));CUDA(cudaFree(ds));free(in);free(cpu);free(gpu);free(ok);
 printf("{\"cases\":%d,\"max_host_device_error\":%.17g,\"periodic_positive\":%d,\"periodic_negative\":%d,\"near_contact_cases\":%d,\"near_contact_nonzero_low\":%d,\"near_contact_normal_error\":%.17g,\"time_evolution\":false}\n",n,error,periodic_positive,periodic_negative,near_count,nonzero_low,normal_error);
 return error<2e-13&&normal_error<2e-13&&periodic_positive>0&&periodic_negative>0&&near_count>0&&nonzero_low>0?0:4;
}
