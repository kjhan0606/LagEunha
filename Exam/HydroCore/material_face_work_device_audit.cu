#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include "../gfs_face_work.h"
__global__ void evaluate(double *out) {
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=4096)return;
 double a=.001*(i+1),sx=cos(a),sy=sin(a);
 out[2*i]=gfs_face_work2d(1+a,sx,sy,.1,-.2,a,-a);
 out[2*i+1]=gfs_face_work_defect2d(1+a,sx,sy,.1,-.2,a,-a,.5*a,-.3*a);
}
int main(void){
 double *d,h[8192];cudaError_t e=cudaMalloc((void**)&d,sizeof(h));
 if(e!=cudaSuccess){fprintf(stderr,"%s\n",cudaGetErrorString(e));return 1;}
 evaluate<<<16,256>>>(d);e=cudaGetLastError();if(e!=cudaSuccess)return 2;
 e=cudaMemcpy(h,d,sizeof(h),cudaMemcpyDeviceToHost);if(e!=cudaSuccess)return 3;
 double err=0;
 for(int i=0;i<4096;i++){
  double a=.001*(i+1),sx=cos(a),sy=sin(a);
  double host[2]={gfs_face_work2d(1+a,sx,sy,.1,-.2,a,-a),gfs_face_work_defect2d(1+a,sx,sy,.1,-.2,a,-a,.5*a,-.3*a)};
  for(int k=0;k<2;k++)err=fmax(err,fabs(h[2*i+k]-host[k])/fmax(1.,fabs(host[k])));
 }
 cudaFree(d);printf("4096 face C/CUDA work and defect relative error %.17g\n",err);
 return err>1e-12;
}
