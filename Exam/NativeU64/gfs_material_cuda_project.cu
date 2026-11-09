/* Uses the existing host projector verbatim on device; no generated duplicate. */
#include <cuda_runtime.h>
#include <stdint.h>
#include <algorithm>
#include <limits.h>
#include "../gfs_material_cuda_project.h"
#include "../gfs_material_scalar_project2d.h"
__global__ static void project_cells(int n,int degree,const int*off,const GfsTraceDirected2D*edge,const GfsTraceDirected2D*fit,const GfsFvPrimitive*c,const unsigned char*walls,double*A,double*b,double*grad,int*status){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
 int ok=1;for(int scalar=0;scalar<2;scalar++)ok &= gfs_project_one_scalar(i,scalar,off,edge,fit,c,walls,A+(size_t)i*degree*6,b+(size_t)i*degree*3,grad);
 status[i]=ok;
}
extern "C" int gfs_cuda_project_gradients_samples(int n,int ns,int nf,const int*off,const GfsTraceDirected2D*edge,const GfsTraceDirected2D*fit,const GfsFvPrimitive*c,const unsigned char*walls,double*out){
 /* Shared scalar routines index grad with 4*i and bound their active-set
  * iterations by 32*(constraints+1), using int arithmetic. */
 if(n<1||n>INT_MAX/4||ns<n||nf<1||!off||!edge||!fit||!c||!out||off[0]!=0||off[n]!=nf)return 0;
 int degree=0,ok=0;
 for(int i=0;i<n;i++){
  if(off[i+1]<off[i]||off[i+1]>nf)return 0;
  degree=std::max(degree,off[i+1]-off[i]);
  for(int k=off[i];k<off[i+1];k++)
   if(edge[k].i!=i||edge[k].j<0||edge[k].j>=ns)return 0;
 }
 if(!degree||degree>(INT_MAX/32-1)/3||
    (size_t)n>SIZE_MAX/(size_t)degree/6/sizeof(double)||
    (size_t)ns>SIZE_MAX/sizeof(GfsFvPrimitive)||
    (size_t)nf>SIZE_MAX/sizeof(GfsTraceDirected2D))return 0;
 int *o=NULL,*status=NULL;GfsTraceDirected2D *e=NULL,*f=NULL;GfsFvPrimitive *q=NULL;unsigned char*w=NULL;double*A=NULL,*b=NULL,*g=NULL;
 int*host=(int*)malloc((size_t)n*sizeof(int));if(!host)return 0;
#define ALLOC(ptr,count) if(cudaMalloc((void**)&ptr,(size_t)(count)*sizeof(*ptr))!=cudaSuccess)goto done
#define COPY(dst,src,count) if(cudaMemcpy(dst,src,(size_t)(count)*sizeof(*dst),cudaMemcpyHostToDevice)!=cudaSuccess)goto done
 ALLOC(o,(size_t)n+1);ALLOC(e,nf);ALLOC(f,nf);ALLOC(q,ns);ALLOC(g,4*(size_t)n);ALLOC(status,n);ALLOC(A,(size_t)n*degree*6);ALLOC(b,(size_t)n*degree*3);
 COPY(o,off,(size_t)n+1);COPY(e,edge,nf);COPY(f,fit,nf);COPY(q,c,ns);
 if(walls){ALLOC(w,nf);COPY(w,walls,nf);}
 project_cells<<<(n-1)/128+1,128>>>(n,degree,o,e,f,q,w,A,b,g,status);
 if(cudaGetLastError()!=cudaSuccess||cudaMemcpy(host,status,(size_t)n*sizeof(int),cudaMemcpyDeviceToHost)!=cudaSuccess)goto done;
 for(int i=0;i<n;i++)if(host[i]!=1){fprintf(stderr,"[CUDA_PROJECT_FAIL] cell=%d status=%d\n",i,host[i]);goto done;}
 if(cudaMemcpy(out,g,4*(size_t)n*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess)goto done;ok=1;
 done:cudaFree(o);cudaFree(e);cudaFree(f);cudaFree(q);cudaFree(w);cudaFree(A);cudaFree(b);cudaFree(g);cudaFree(status);free(host);return ok;
}
extern "C" int gfs_cuda_project_gradients(int n,int nf,const int*off,const GfsTraceDirected2D*edge,const GfsTraceDirected2D*fit,const GfsFvPrimitive*c,const unsigned char*walls,double*out){
 return gfs_cuda_project_gradients_samples(n,n,nf,off,edge,fit,c,walls,out);
}
