/* Two-cell device smoke check for the isolated Jacobi volume-rate kernel. */
#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "exam_gpu.h"

extern "C" int gpu_probe_volume_rate_jacobi(GPUContext *,
    const unsigned char *, int, int, int, double, double *,
    unsigned long long *);

#define CHECK(call) do { cudaError_t e=(call); if(e!=cudaSuccess) { \
    fprintf(stderr,"CUDA %s: %s\n",#call,cudaGetErrorString(e)); return 2; \
} } while(0)
#define MANAGED(ptr,n) CHECK(cudaMallocManaged((void**)&(ptr),(n)*sizeof(*(ptr))))

int main(void)
{
    GPUContext ctx;
    memset(&ctx,0,sizeof(ctx));
    ctx.initialized=1;
    ctx.d_faces.n_particles=2;
    ctx.d_faces.n_faces_total=2;
    MANAGED(ctx.d_faces.face_offset,3);
    MANAGED(ctx.d_faces.neighbor_idx,2);
    MANAGED(ctx.d_faces.is_ghost,2);
    MANAGED(ctx.d_faces.c1x,2); MANAGED(ctx.d_faces.c1y,2);
    MANAGED(ctx.d_faces.c2x,2); MANAGED(ctx.d_faces.c2y,2);
    MANAGED(ctx.d_parts.x,2); MANAGED(ctx.d_parts.y,2);
    MANAGED(ctx.d_parts.vx,2); MANAGED(ctx.d_parts.vy,2);
    ctx.d_faces.face_offset[0]=0;
    ctx.d_faces.face_offset[1]=1;
    ctx.d_faces.face_offset[2]=2;
    ctx.d_faces.neighbor_idx[0]=1;
    ctx.d_faces.neighbor_idx[1]=0;
    ctx.d_faces.is_ghost[0]=ctx.d_faces.is_ghost[1]=0;
    ctx.d_faces.c1x[0]=.5; ctx.d_faces.c1y[0]=-.5;
    ctx.d_faces.c2x[0]=.5; ctx.d_faces.c2y[0]=.5;
    ctx.d_faces.c1x[1]=-.5; ctx.d_faces.c1y[1]=.5;
    ctx.d_faces.c2x[1]=-.5; ctx.d_faces.c2y[1]=-.5;
    ctx.d_parts.x[0]=0.; ctx.d_parts.x[1]=1.;
    ctx.d_parts.y[0]=ctx.d_parts.y[1]=0.;
    ctx.d_parts.vx[0]=0.; ctx.d_parts.vx[1]=1.;
    ctx.d_parts.vy[0]=ctx.d_parts.vy[1]=0.;

    unsigned char active[2]={1,1};
    double rate[2]={0.,0.};
    unsigned long long unsupported=999;
    int status=gpu_probe_volume_rate_jacobi(&ctx,active,2,2,16,.8,rate,&unsupported);
    printf("both_active status=%d unsupported=%llu rates=%.17g %.17g\n",
           status,unsupported,rate[0],rate[1]);
    if(status || unsupported || fabs(rate[0]+.5)>2e-4 ||
       fabs(rate[1]-.5)>2e-4) return 3;
    active[1]=0;
    rate[0]=rate[1]=0.;
    status=gpu_probe_volume_rate_jacobi(&ctx,active,2,2,16,.8,rate,&unsupported);
    printf("one_active status=%d unsupported=%llu rates=%.17g %.17g\n",
           status,unsupported,rate[0],rate[1]);
    if(status || unsupported || fabs(rate[0]+1.)>1e-9 || rate[1]!=0.) return 4;
    ctx.d_faces.is_ghost[0]=1;
    status=gpu_probe_volume_rate_jacobi(&ctx,active,2,2,16,.8,rate,&unsupported);
    printf("unsupported_wall status=%d unsupported=%llu\n",status,unsupported);
    if(status!=-2 || unsupported<1) return 5;
    return 0;
}
