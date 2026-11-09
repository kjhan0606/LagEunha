#include <cuda_runtime.h>
#include <stdio.h>
#include <limits.h>
#include <time.h>
static double weight_wall_seconds(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return (double)t.tv_sec+1e-9*t.tv_nsec;}
#include "../gfs_material_cuda_weight.h"
/* Jacobi-PCG on the same directed Laplacian, with projected preconditioner.
 * Scalars stay on device except the convergence certificate. No tolerance
 * relaxation or solver fallback. Final caller also certifies original L. */
static __device__ double block_sum(double v,double*work){int i=threadIdx.x;work[i]=v;__syncthreads();for(int d=blockDim.x/2;d;d/=2){if(i<d)work[i]+=work[i+d];__syncthreads();}return work[0];}
static __global__ void mv(int n,const int*off,const int*col,const double*a,const double*x,double*y){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n){double z=0;for(int k=off[i];k<off[i+1];k++)z+=a[k]*(x[i]-x[col[k]]);y[i]=z;}}
static __global__ void residual_mean(int n,const double*r,const double*d,double*s,double norm,int iteration){if(iteration>=0&&s[6])return;__shared__ double work[256];double rr=0,z=0;for(int i=threadIdx.x;i<n;i+=blockDim.x){rr+=r[i]*r[i];z+=r[i]/d[i];}double nr=block_sum(rr,work);if(threadIdx.x==0)s[2]=nr;__syncthreads();double mean=block_sum(z,work);if(threadIdx.x==0){s[3]=mean/n;if(iteration<0)s[6]=0;else if(nr<=1e-20*norm){s[6]=1;s[7]=iteration+1;}}}
static __global__ void precondition_scalar(int n,const double*r,const double*d,double*s,int restart){if(s[6])return;__shared__ double work[256];double v=0;for(int i=threadIdx.x;i<n;i+=blockDim.x)v+=r[i]*(r[i]/d[i]-s[3]);double next=block_sum(v,work);if(threadIdx.x==0){s[4]=restart?0:next/s[0];s[0]=next;if(!(next>0)||!isfinite(next)||!isfinite(s[4]))s[5]=0;}}
static __global__ void search_update(int n,const double*r,const double*d,double*p,const double*s,int restart){if(s[6])return;int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)p[i]=r[i]/d[i]-s[3]+(restart?0:s[4]*p[i]);}
static __global__ void curvature(int n,const double*p,const double*ap,double*s){if(s[6])return;__shared__ double work[256];double v=0;for(int i=threadIdx.x;i<n;i+=blockDim.x)v+=p[i]*ap[i];double c=block_sum(v,work);if(threadIdx.x==0){s[1]=s[0]/c;if(!(c>0)||!isfinite(c)||!isfinite(s[1]))s[5]=0;}}
static __global__ void solution_update(int n,double*q,double*r,const double*p,const double*ap,const double*s){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n&&s[5]&&!s[6]){q[i]+=s[1]*p[i];r[i]-=s[1]*ap[i];}}
static __global__ void true_residual(int n,const double*rhs,const double*ap,double*r){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)r[i]=rhs[i]-ap[i];}
#ifdef GFS_MATERIAL_CUDA_PERSISTENT
#include "gfs_material_weight_persistent.cuh"
#endif
/* Experimental numerical initial guesses only; every output still uses the
 * original operator and original relative true-residual acceptance. Per-host
 * thread and device history is discarded on size/device changes. */
struct WeightHistory {
 int n=0,device=-1,count=0,next=0; double *x[2]={NULL,NULL};
 ~WeightHistory(){free(x[0]);free(x[1]);}
};
static thread_local WeightHistory weight_history;
extern "C" int gfs_material_cuda_weight_pcg(int n,int nf,const GfsTraceDirected2D*e,const unsigned char*wall,const double*rhs,const double*diag,double*out,int*iterations,double*relative){
 if(n<1||nf<1||!e||!rhs||!diag||!out||!iterations||!relative)return 0;
 const char *recycle_option=getenv("GFS_MATERIAL_WEIGHT_RECYCLE");
 bool recycle=recycle_option&&!strcmp(recycle_option,"1");
 double *guess=NULL,*initial_r=NULL,*candidate=NULL;double best_residual=0;int selected=-1,device=-1;
 double t_start=weight_wall_seconds(),t_pack=0,t_upload=0,t_plan=0,t_iter=0,t_copy=0;
 int ok=0,*off=(int*)calloc(n+1,sizeof(int)),*cursor=(int*)malloc(n*sizeof(int)),*col=(int*)malloc(nf*sizeof(int)),*dr=NULL,*dc=NULL;double*coef=(double*)malloc(nf*sizeof(double)),*da=NULL,*buf=NULL,*sc=NULL;double norm=0,hs[8]={0,0,0,0,0,1,0,0};int blocks=(n+255)/256;double *q=NULL,*r=NULL,*p=NULL,*ap=NULL,*dd=NULL,*br=NULL;
#ifdef GFS_MATERIAL_CUDA_GRAPH
 cudaStream_t stream=NULL;cudaGraph_t graph=NULL;cudaGraphExec_t executable=NULL;int capturing=0;
#endif
#ifdef GFS_MATERIAL_CUDA_PERSISTENT
 double*partial=NULL;
#endif
 if(!off||!cursor||!col||!coef)goto done;
 for(int i=0;i<n;i++){if(!(diag[i]>0)||!isfinite(diag[i])||!isfinite(rhs[i]))goto done;norm+=rhs[i]*rhs[i];}if(!(norm>0)||!isfinite(norm))goto done;
 for(int k=0;k<nf;k++){if(e[k].i<0||e[k].i>=n||e[k].j<0||e[k].j>=n)goto done;off[e[k].i+1]++;}for(int i=1;i<=n;i++)off[i]+=off[i-1];memcpy(cursor,off,n*sizeof(int));
 for(int k=0;k<nf;k++){double d=hypot(e[k].dx,e[k].dy);if(!(d>0)||!isfinite(d))goto done;int z=cursor[e[k].i]++;col[z]=e[k].j;coef[z]=(wall&&wall[k])?0:e[k].area/(2*d);if(!isfinite(coef[z])||coef[z]<0)goto done;}
 if(recycle){
  if(cudaGetDevice(&device)!=cudaSuccess)goto done;
  if(weight_history.n!=n||weight_history.device!=device){
   free(weight_history.x[0]);free(weight_history.x[1]);
   weight_history.x[0]=(double*)calloc(n,sizeof(double));
   weight_history.x[1]=(double*)calloc(n,sizeof(double));
   weight_history.n=n;weight_history.device=device;weight_history.count=0;weight_history.next=0;
  }
  if(!weight_history.x[0]||!weight_history.x[1])goto done;
  guess=(double*)calloc(n,sizeof(double));initial_r=(double*)malloc(n*sizeof(double));candidate=(double*)malloc(n*sizeof(double));
  if(!guess||!initial_r||!candidate)goto done;
  memcpy(initial_r,rhs,n*sizeof(double));best_residual=norm;
  for(int h=0;h<weight_history.count;h++){
   double dot=0,aa=0;
   for(int i=0;i<n;i++){
    double v=0;for(int k=off[i];k<off[i+1];k++)v+=coef[k]*(weight_history.x[h][i]-weight_history.x[h][col[k]]);
    candidate[i]=v;dot+=rhs[i]*v;aa+=v*v;
   }
   if(!(aa>0)||!isfinite(aa))continue;
   double scale=dot/aa,res=0;
   for(int i=0;i<n;i++){double r0=rhs[i]-scale*candidate[i];res+=r0*r0;}
   /* Avoid a zero-residual PCG entry; cold solve remains valid in that case. */
   if(isfinite(res)&&res<best_residual&&res>1e-20*norm){
    best_residual=res;selected=h;
    for(int i=0;i<n;i++){guess[i]=scale*weight_history.x[h][i];initial_r[i]=rhs[i]-scale*candidate[i];}
   }
  }
 }
 t_pack=weight_wall_seconds();
#define WCALL(call) do{cudaError_t err=(call);if(err!=cudaSuccess){fprintf(stderr,"[MATERIAL_WEIGHT_CUDA_FAIL] %s\n",cudaGetErrorString(err));ok=0;goto done;}}while(0)
 WCALL(cudaMalloc((void**)&dr,(n+1)*sizeof(int)));WCALL(cudaMalloc((void**)&dc,nf*sizeof(int)));WCALL(cudaMalloc((void**)&da,nf*sizeof(double)));WCALL(cudaMalloc((void**)&buf,(size_t)n*6*sizeof(double)));WCALL(cudaMalloc((void**)&sc,8*sizeof(double)));q=buf;r=q+n;p=r+n;ap=p+n;dd=ap+n;br=dd+n;
 WCALL(cudaMemcpy(dr,off,(n+1)*sizeof(int),cudaMemcpyHostToDevice));WCALL(cudaMemcpy(dc,col,nf*sizeof(int),cudaMemcpyHostToDevice));WCALL(cudaMemcpy(da,coef,nf*sizeof(double),cudaMemcpyHostToDevice));if(recycle){WCALL(cudaMemcpy(q,guess,n*sizeof(double),cudaMemcpyHostToDevice));}else{WCALL(cudaMemset(q,0,n*sizeof(double)));}WCALL(cudaMemcpy(r,recycle?initial_r:rhs,n*sizeof(double),cudaMemcpyHostToDevice));WCALL(cudaMemcpy(br,rhs,n*sizeof(double),cudaMemcpyHostToDevice));WCALL(cudaMemcpy(dd,diag,n*sizeof(double),cudaMemcpyHostToDevice));WCALL(cudaMemcpy(sc,hs,sizeof(hs),cudaMemcpyHostToDevice));t_upload=weight_wall_seconds();
#ifdef GFS_MATERIAL_CUDA_PERSISTENT
 {
  int dev,cooperative,sm,per_sm;
  WCALL(cudaGetDevice(&dev));WCALL(cudaDeviceGetAttribute(&cooperative,cudaDevAttrCooperativeLaunch,dev));
  if(!cooperative)goto done;
  WCALL(cudaDeviceGetAttribute(&sm,cudaDevAttrMultiProcessorCount,dev));
  WCALL(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm,weight_persistent<false>,256,0));
  int resident=blocks<sm*per_sm?blocks:sm*per_sm;if(resident<1)goto done;
  WCALL(cudaMalloc((void**)&partial,resident*sizeof(double)));t_plan=weight_wall_seconds();
  int remaining=n>256?n:256;*iterations=0;*relative=1;
  while(remaining>0){
   void*args[]={&n,&dr,&dc,&da,&dd,&q,&r,&p,&ap,&sc,&partial,&norm,&remaining};
   WCALL(cudaLaunchCooperativeKernel((void*)weight_persistent<false>,dim3(resident),dim3(256),args,0,0));
   WCALL(cudaMemcpy(hs,sc,sizeof(hs),cudaMemcpyDeviceToHost));
   int used=(int)hs[7];if(!hs[5]||used<1||used>remaining)goto done;
   *iterations+=used;remaining-=used;
   /* Same true-residual replacement/restart and tolerance as batched PCG. */
   mv<<<blocks,256>>>(n,dr,dc,da,q,ap);true_residual<<<blocks,256>>>(n,br,ap,r);residual_mean<<<1,256>>>(n,r,dd,sc,norm,-1);
   WCALL(cudaGetLastError());WCALL(cudaMemcpy(hs,sc,sizeof(hs),cudaMemcpyDeviceToHost));
   *relative=sqrt(hs[2]/norm);if(!isfinite(*relative))goto done;
   if(*relative<=1e-10){ok=1;break;}
  }
 }
#elif defined(GFS_MATERIAL_CUDA_GRAPH)
 WCALL(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
 residual_mean<<<1,256,0,stream>>>(n,r,dd,sc,norm,-1);precondition_scalar<<<1,256,0,stream>>>(n,r,dd,sc,1);search_update<<<blocks,256,0,stream>>>(n,r,dd,p,sc,1);
 WCALL(cudaStreamSynchronize(stream));
 WCALL(cudaStreamBeginCapture(stream,cudaStreamCaptureModeThreadLocal));capturing=1;
 for(int j=0;j<8;j++){
  mv<<<blocks,256,0,stream>>>(n,dr,dc,da,p,ap);curvature<<<1,256,0,stream>>>(n,p,ap,sc);solution_update<<<blocks,256,0,stream>>>(n,q,r,p,ap,sc);residual_mean<<<1,256,0,stream>>>(n,r,dd,sc,norm,j);precondition_scalar<<<1,256,0,stream>>>(n,r,dd,sc,0);search_update<<<blocks,256,0,stream>>>(n,r,dd,p,sc,0);
 }
 WCALL(cudaStreamEndCapture(stream,&graph));capturing=0;WCALL(cudaGraphInstantiateWithFlags(&executable,graph,0));
 t_plan=weight_wall_seconds();
 *iterations=0;*relative=1;
 for(int it=0;it<(n>256?n:256);){
  int batch=((n>256?n:256)-it>=8)?8:(n>256?n:256)-it;
  if(batch==8){WCALL(cudaGraphLaunch(executable,stream));}
  else for(int j=0;j<batch;j++){
   mv<<<blocks,256,0,stream>>>(n,dr,dc,da,p,ap);curvature<<<1,256,0,stream>>>(n,p,ap,sc);solution_update<<<blocks,256,0,stream>>>(n,q,r,p,ap,sc);residual_mean<<<1,256,0,stream>>>(n,r,dd,sc,norm,j);precondition_scalar<<<1,256,0,stream>>>(n,r,dd,sc,0);search_update<<<blocks,256,0,stream>>>(n,r,dd,p,sc,0);
  }
  WCALL(cudaGetLastError());WCALL(cudaMemcpyAsync(hs,sc,sizeof(hs),cudaMemcpyDeviceToHost,stream));WCALL(cudaStreamSynchronize(stream));*iterations=it+(hs[6]?(int)hs[7]:batch);it+=batch;
  if(!hs[5]||!isfinite(hs[2]))goto done;
  if(hs[6]||sqrt(hs[2]/norm)<=1e-10){
   mv<<<blocks,256,0,stream>>>(n,dr,dc,da,q,ap);true_residual<<<blocks,256,0,stream>>>(n,br,ap,r);residual_mean<<<1,256,0,stream>>>(n,r,dd,sc,norm,-1);
   WCALL(cudaGetLastError());WCALL(cudaMemcpyAsync(hs,sc,sizeof(hs),cudaMemcpyDeviceToHost,stream));WCALL(cudaStreamSynchronize(stream));*relative=sqrt(hs[2]/norm);
   if(isfinite(*relative)&&*relative<=1e-10){ok=1;break;}
   precondition_scalar<<<1,256,0,stream>>>(n,r,dd,sc,1);search_update<<<blocks,256,0,stream>>>(n,r,dd,p,sc,1);
  }
 }
#else
 residual_mean<<<1,256>>>(n,r,dd,sc,norm,-1);precondition_scalar<<<1,256>>>(n,r,dd,sc,1);search_update<<<blocks,256>>>(n,r,dd,p,sc,1);
 t_plan=weight_wall_seconds();
 *iterations=0;*relative=1;
 for(int it=0;it<(n>256?n:256);it++){
  mv<<<blocks,256>>>(n,dr,dc,da,p,ap);curvature<<<1,256>>>(n,p,ap,sc);solution_update<<<blocks,256>>>(n,q,r,p,ap,sc);residual_mean<<<1,256>>>(n,r,dd,sc,norm,it);
  /* Batch host observations only: kernels freeze q/r at the first recursive
   * convergence, then the same true-residual certificate decides acceptance. */
  if((it+1)%8&&it+1<(n>256?n:256)){precondition_scalar<<<1,256>>>(n,r,dd,sc,0);search_update<<<blocks,256>>>(n,r,dd,p,sc,0);continue;}
  WCALL(cudaGetLastError());WCALL(cudaMemcpy(hs,sc,sizeof(hs),cudaMemcpyDeviceToHost));*iterations=hs[6]?(int)hs[7]:it+1;if(!hs[5]||!isfinite(hs[2]))goto done;
  if(hs[6]||sqrt(hs[2]/norm)<=1e-10){
   mv<<<blocks,256>>>(n,dr,dc,da,q,ap);true_residual<<<blocks,256>>>(n,br,ap,r);residual_mean<<<1,256>>>(n,r,dd,sc,norm,-1);WCALL(cudaGetLastError());WCALL(cudaMemcpy(hs,sc,sizeof(hs),cudaMemcpyDeviceToHost));*relative=sqrt(hs[2]/norm);
   if(isfinite(*relative)&&*relative<=1e-10){ok=1;break;}
   precondition_scalar<<<1,256>>>(n,r,dd,sc,1);search_update<<<blocks,256>>>(n,r,dd,p,sc,1);
  }else{precondition_scalar<<<1,256>>>(n,r,dd,sc,0);search_update<<<blocks,256>>>(n,r,dd,p,sc,0);}
 }
#endif
 t_iter=weight_wall_seconds();
 if(ok)WCALL(cudaMemcpy(out,q,n*sizeof(double),cudaMemcpyDeviceToHost));
 t_copy=weight_wall_seconds();
 done:
#ifdef GFS_MATERIAL_CUDA_GRAPH
 if(capturing){cudaGraph_t failed=NULL;cudaStreamEndCapture(stream,&failed);if(failed)cudaGraphDestroy(failed);}
 if(stream)cudaStreamSynchronize(stream);if(executable)cudaGraphExecDestroy(executable);if(graph)cudaGraphDestroy(graph);if(stream)cudaStreamDestroy(stream);
#endif
#ifdef GFS_MATERIAL_CUDA_PERSISTENT
 cudaFree(partial);
#endif
 cudaFree(dr);cudaFree(dc);cudaFree(da);cudaFree(buf);cudaFree(sc);free(off);free(cursor);free(col);free(coef);
 if(ok&&recycle){
  memcpy(weight_history.x[weight_history.next],out,n*sizeof(double));
  weight_history.next=(weight_history.next+1)%2;if(weight_history.count<2)weight_history.count++;
  if(getenv("GFS_MATERIAL_WEIGHT_PROFILE"))fprintf(stderr,"[MATERIAL_WEIGHT_RECYCLE] selected=%d initial_relative=%.17g iterations=%d certified=%.17g\n",selected,sqrt(best_residual/norm),*iterations,*relative);
 }
 free(guess);free(initial_r);free(candidate);
 if(ok&&getenv("GFS_MATERIAL_WEIGHT_PROFILE"))fprintf(stderr,"[MATERIAL_WEIGHT_PROFILE] n=%d iterations=%d host_pack_ms=%.9g upload_ms=%.9g plan_ms=%.9g iterate_ms=%.9g copy_ms=%.9g cleanup_ms=%.9g total_ms=%.9g\n",n,*iterations,1e3*(t_pack-t_start),1e3*(t_upload-t_pack),1e3*(t_plan-t_upload),1e3*(t_iter-t_plan),1e3*(t_copy-t_iter),1e3*(weight_wall_seconds()-t_copy),1e3*(weight_wall_seconds()-t_start));
 return ok;
#undef WCALL
}

/* Anchored SPD graph. Interior edges occur in both row lists; prescribed
 * exterior edges use col=-1, meaning the unknown exterior rate is zero.
 * Known exterior rates are already on rhs. There is NO mean subtraction. */
extern "C" int gfs_material_cuda_exterior_pcg(int n,int nt,int nf,
 const GfsMaterialExteriorEdge*edge,const double*rhs,double*out,int*iterations,double*relative){
#ifndef GFS_MATERIAL_CUDA_PERSISTENT
 (void)n;(void)nt;(void)nf;(void)edge;(void)rhs;(void)out;(void)iterations;(void)relative;
 return 0;
#else
 if(n<1||n>=INT_MAX||nt<=n||nf<1||nf>INT_MAX/2||!edge||!rhs||!out||!iterations||!relative)return 0;
 int ok=0,nnz=0,*off=NULL,*cursor=NULL,*col=NULL,*dr=NULL,*dc=NULL;
 double *diag=NULL,*coef=NULL,*host=NULL,*da=NULL,*buf=NULL,*sc=NULL,*partial=NULL;
 double *q=NULL,*r=NULL,*p=NULL,*ap=NULL,*dd=NULL,norm=0,hs[8]={0};
 int blocks=(n-1)/256+1;double start=weight_wall_seconds(),packed=0,uploaded=0,iterated=0;
 *iterations=0;*relative=INFINITY;
 off=(int*)calloc((size_t)n+1,sizeof(int));cursor=(int*)malloc((size_t)n*sizeof(int));
 col=(int*)malloc((size_t)2*nf*sizeof(int));coef=(double*)malloc((size_t)2*nf*sizeof(double));
 diag=(double*)malloc((size_t)n*sizeof(double));host=(double*)malloc((size_t)2*n*sizeof(double));
 if(!off||!cursor||!col||!coef||!diag||!host||!gfs_material_exterior_diagonal(n,nt,nf,edge,diag))goto done;
 for(int i=0;i<n;i++){if(!isfinite(rhs[i]))goto done;norm+=rhs[i]*rhs[i];}
 if(!isfinite(norm))goto done;
 if(norm==0){memset(out,0,(size_t)n*sizeof(double));*relative=0;ok=1;goto done;}
 for(int k=0;k<nf;k++){off[edge[k].i+1]++;if(edge[k].j<n)off[edge[k].j+1]++;}
 for(int i=1;i<=n;i++)off[i]+=off[i-1];nnz=off[n];memcpy(cursor,off,(size_t)n*sizeof(int));
 for(int k=0;k<nf;k++){
  int i=edge[k].i,j=edge[k].j,z=cursor[i]++;col[z]=j<n?j:-1;coef[z]=edge[k].coefficient;
  if(j<n){z=cursor[j]++;col[z]=i;coef[z]=edge[k].coefficient;}
 }
 packed=weight_wall_seconds();
#define ECALL(call) do{cudaError_t err=(call);if(err!=cudaSuccess){fprintf(stderr,"[MATERIAL_EXTERIOR_CUDA_FAIL] %s\n",cudaGetErrorString(err));ok=0;goto done;}}while(0)
 ECALL(cudaMalloc((void**)&dr,((size_t)n+1)*sizeof(int)));ECALL(cudaMalloc((void**)&dc,(size_t)nnz*sizeof(int)));
 ECALL(cudaMalloc((void**)&da,(size_t)nnz*sizeof(double)));ECALL(cudaMalloc((void**)&buf,(size_t)5*n*sizeof(double)));ECALL(cudaMalloc((void**)&sc,sizeof(hs)));
 q=buf;r=q+n;p=r+n;ap=p+n;dd=ap+n;
 ECALL(cudaMemcpy(dr,off,((size_t)n+1)*sizeof(int),cudaMemcpyHostToDevice));ECALL(cudaMemcpy(dc,col,(size_t)nnz*sizeof(int),cudaMemcpyHostToDevice));
 ECALL(cudaMemcpy(da,coef,(size_t)nnz*sizeof(double),cudaMemcpyHostToDevice));ECALL(cudaMemcpy(dd,diag,(size_t)n*sizeof(double),cudaMemcpyHostToDevice));
 ECALL(cudaMemset(q,0,(size_t)n*sizeof(double)));ECALL(cudaMemcpy(r,rhs,(size_t)n*sizeof(double),cudaMemcpyHostToDevice));
 {
 int device,cooperative,sm,per_sm;
 ECALL(cudaGetDevice(&device));ECALL(cudaDeviceGetAttribute(&cooperative,cudaDevAttrCooperativeLaunch,device));if(!cooperative)goto done;
 ECALL(cudaDeviceGetAttribute(&sm,cudaDevAttrMultiProcessorCount,device));
 ECALL(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm,weight_persistent<true>,256,0));
 int resident=blocks<sm*per_sm?blocks:sm*per_sm;if(resident<1)goto done;
 ECALL(cudaMalloc((void**)&partial,(size_t)resident*sizeof(double)));uploaded=weight_wall_seconds();
 int remaining=n>256?n:256;
 while(remaining>0){
  void*args[]={&n,&dr,&dc,&da,&dd,&q,&r,&p,&ap,&sc,&partial,&norm,&remaining};
  ECALL(cudaLaunchCooperativeKernel((void*)weight_persistent<true>,dim3(resident),dim3(256),args,0,0));
  ECALL(cudaMemcpy(hs,sc,sizeof(hs),cudaMemcpyDeviceToHost));
  int used=(int)hs[7];if(!hs[5]||used<1||used>remaining)goto done;
  *iterations+=used;remaining-=used;
  ECALL(cudaMemcpy(host,q,(size_t)n*sizeof(double),cudaMemcpyDeviceToHost));
  /* The ORIGINAL edge action certifies the GPU result. Roundoff from CSR
   * reduction is never used to loosen the true-residual threshold. */
  gfs_material_exterior_matvec(n,nf,edge,host,host+n);
  double rr=0;for(int i=0;i<n;i++){host[n+i]=rhs[i]-host[n+i];rr+=host[n+i]*host[n+i];}
  *relative=sqrt(rr/norm);if(!isfinite(*relative))goto done;
  if(*relative<=1e-10){memcpy(out,host,(size_t)n*sizeof(double));ok=1;break;}
  ECALL(cudaMemcpy(r,host+n,(size_t)n*sizeof(double),cudaMemcpyHostToDevice));
 }
 }
 iterated=weight_wall_seconds();
 done:
 cudaFree(dr);cudaFree(dc);cudaFree(da);cudaFree(buf);cudaFree(sc);cudaFree(partial);
 free(off);free(cursor);free(col);free(diag);free(coef);free(host);
 if(ok&&norm>0&&getenv("GFS_MATERIAL_WEIGHT_PROFILE"))fprintf(stderr,"[MATERIAL_EXTERIOR_CUDA_PROFILE] n=%d iterations=%d relative=%.17g pack_ms=%.9g upload_ms=%.9g iterate_certificate_ms=%.9g cleanup_ms=%.9g total_ms=%.9g\n",n,*iterations,*relative,1e3*(packed-start),1e3*(uploaded-packed),1e3*(iterated-uploaded),1e3*(weight_wall_seconds()-iterated),1e3*(weight_wall_seconds()-start));
 return ok;
#undef ECALL
#endif
}
