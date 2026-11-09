/* Native single-rank owner chart; included after CUDA_CHECK. */
#include "../gfs_power_chart2d.h"
static double chart_period[2]={0,0};
static int *chart_device_owner=NULL,chart_capacity=0,chart_device=-1;
extern "C" void gpu_upload_power_coordinate_lows(GPUContext*ctx,int nt,
 const double*lowx,const double*lowy){
 if(!ctx||!ctx->initialized||nt<0||((lowx==NULL)!=(lowy==NULL))){fprintf(stderr,"[POWER_COORDINATES] invalid input\n");exit(99);}
 if(!lowx){ctx->power_low_count=0;return;}
 if(nt>ctx->power_low_capacity){
  if(ctx->d_power_lowx)CUDA_CHECK(cudaFree(ctx->d_power_lowx));
  if(ctx->d_power_lowy)CUDA_CHECK(cudaFree(ctx->d_power_lowy));
  CUDA_CHECK(cudaMalloc(&ctx->d_power_lowx,nt*sizeof(double)));
  CUDA_CHECK(cudaMalloc(&ctx->d_power_lowy,nt*sizeof(double)));
  ctx->power_low_capacity=nt;
 }
 for(int i=0;i<nt;i++)if(!isfinite(lowx[i])||!isfinite(lowy[i])){fprintf(stderr,"[POWER_COORDINATES] nonfinite low part\n");exit(99);}
 CUDA_CHECK(cudaMemcpyAsync(ctx->d_power_lowx,lowx,nt*sizeof(double),cudaMemcpyHostToDevice,(cudaStream_t)ctx->stream));
 CUDA_CHECK(cudaMemcpyAsync(ctx->d_power_lowy,lowy,nt*sizeof(double),cudaMemcpyHostToDevice,(cudaStream_t)ctx->stream));
 ctx->power_low_count=nt;
}
extern "C" void gpu_power_chart_domain(double lx,double ly){
 if(lx<0||ly<0||(lx==0&&ly!=0)||!isfinite(lx)||!isfinite(ly)){fprintf(stderr,"[POWER_CHART] invalid domain\n");exit(99);}
 chart_period[0]=lx;chart_period[1]=ly;
}
extern "C" const int *gfs_power_chart_owners(GPUContext*,int,int,double,double);
static void prepare_power_chart(GPUContext *ctx,int n,int nt,cudaStream_t stream){
 if(ctx->power_low_count&&ctx->power_low_count!=nt){fprintf(stderr,"[POWER_COORDINATES] stale state count\n");exit(99);}
 if(chart_period[0]==0)return; /* Explicit nonperiodic physical exterior. */
 if(!(chart_period[0]>0)||n<1||nt<n){fprintf(stderr,"[POWER_CHART] domain/count missing\n");exit(99);}
 const int *alias=gfs_power_chart_owners(ctx,n,nt,chart_period[0],chart_period[1]);
 if(!alias){fprintf(stderr,"[POWER_CHART] certified owner resolver failed\n");exit(99);}
 int dev;CUDA_CHECK(cudaGetDevice(&dev));
 if(chart_device!=-1&&chart_device!=dev){fprintf(stderr,"[POWER_CHART] device changed\n");exit(99);}chart_device=dev;
 if(nt>chart_capacity){if(chart_device_owner)CUDA_CHECK(cudaFree(chart_device_owner));CUDA_CHECK(cudaMalloc(&chart_device_owner,nt*sizeof(int)));chart_capacity=nt;}
 CUDA_CHECK(cudaMemcpyAsync(chart_device_owner,alias,nt*sizeof(int),cudaMemcpyHostToDevice,stream));
}

