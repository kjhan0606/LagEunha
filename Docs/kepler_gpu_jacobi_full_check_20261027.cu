/* Full accepted t=6.5 FaceCSR comparison against CPU Jacobi16 reference. */
#include <cuda_runtime.h>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "exam_gpu.h"

extern "C" int gpu_probe_volume_rate_jacobi(GPUContext *,
    const unsigned char *, int, int, int, double, double *,
    unsigned long long *);

template<class T> static T *read_managed(FILE *fp, size_t n)
{
    T *p=0;
    cudaError_t e=cudaMallocManaged((void**)&p,n*sizeof(T));
    if(e!=cudaSuccess || fread(p,sizeof(T),n,fp)!=n) {
        fprintf(stderr,"allocation/read failure: %s\n",cudaGetErrorString(e));
        exit(2);
    }
    return p;
}

int main(int argc, char **argv)
{
    if(argc!=2) { fprintf(stderr,"usage: %s exported.bin\n",argv[0]); return 2; }
    FILE *fp=fopen(argv[1],"rb");
    if(!fp) { perror(argv[1]); return 2; }
    int header[4];
    if(fread(header,sizeof(int),4,fp)!=4 || header[0]!=0x4B574452 ||
       header[1]!=55351 || header[2]!=297114 || header[3]!=16) {
        fprintf(stderr,"bad exported CSR header\n"); return 2;
    }
    const int n=header[1], nf=header[2];
    GPUContext ctx;
    memset(&ctx,0,sizeof(ctx));
    ctx.initialized=1;
    ctx.d_faces.n_particles=n;
    ctx.d_faces.n_faces_total=nf;
    ctx.d_faces.face_offset=read_managed<int>(fp,n+1);
    ctx.d_faces.neighbor_idx=read_managed<int>(fp,nf);
    ctx.d_faces.is_ghost=read_managed<int>(fp,nf);
    ctx.d_faces.c1x=read_managed<double>(fp,nf);
    ctx.d_faces.c1y=read_managed<double>(fp,nf);
    ctx.d_faces.c2x=read_managed<double>(fp,nf);
    ctx.d_faces.c2y=read_managed<double>(fp,nf);
    ctx.d_parts.x=read_managed<double>(fp,n);
    ctx.d_parts.y=read_managed<double>(fp,n);
    ctx.d_parts.vx=read_managed<double>(fp,n);
    ctx.d_parts.vy=read_managed<double>(fp,n);
    unsigned char *active=read_managed<unsigned char>(fp,n);
    double *expected=read_managed<double>(fp,n);
    if(fgetc(fp)!=EOF || ctx.d_faces.face_offset[n]!=nf) {
        fprintf(stderr,"export length or offset mismatch\n"); return 2;
    }
    fclose(fp);
    double *rate=(double*)malloc((size_t)n*sizeof(double));
    if(!rate) return 2;
    unsigned long long unsupported=999;
    auto start=std::chrono::steady_clock::now();
    int status=gpu_probe_volume_rate_jacobi(&ctx,active,n,n,16,.8,rate,&unsupported);
    auto end=std::chrono::steady_clock::now();
    double seconds=std::chrono::duration<double>(end-start).count();
    double error2=0.,reference2=0.,max_abs=0.;
    int nonfinite=0, active_count=0;
    for(int i=0;i<n;++i) {
        if(!active[i]) continue;
        ++active_count;
        if(!std::isfinite(rate[i])) { ++nonfinite; continue; }
        double error=rate[i]-expected[i];
        error2+=error*error;
        reference2+=expected[i]*expected[i];
        if(std::fabs(error)>max_abs) max_abs=std::fabs(error);
    }
    double relative_l2=std::sqrt(error2/reference2);
    printf("case=%s active=%d faces=%d status=%d unsupported=%llu nonfinite=%d elapsed_seconds=%.9g max_abs=%.17g relative_l2=%.17g\n",
           argv[1],active_count,nf,status,unsupported,nonfinite,seconds,max_abs,relative_l2);
    if(status || unsupported || nonfinite || active_count!=49038 ||
       max_abs>1e-10 || relative_l2>1e-9) return 3;
    return 0;
}
