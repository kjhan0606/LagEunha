/* Check projected CUDA Wdot and pressure gradient on accepted W0 CSR. */
#include <cuda_runtime.h>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "exam_gpu.h"

extern "C" int gpu_probe_energy_neutral_rate(GPUContext *,
    const unsigned char *,int,int,int,double,double *,double *,double *,
    unsigned long long *);

template<class T> static T *read_managed(FILE *stream,size_t count)
{
    T *p=0;
    cudaError_t e=cudaMallocManaged((void**)&p,count*sizeof(T));
    if(e!=cudaSuccess || fread(p,sizeof(T),count,stream)!=count) {
        fprintf(stderr,"managed allocation/read failed: %s\n",cudaGetErrorString(e));
        exit(2);
    }
    return p;
}

int main(int argc,char **argv)
{
    if(argc!=2) { fprintf(stderr,"usage: %s exported.bin\n",argv[0]); return 2; }
    FILE *stream=fopen(argv[1],"rb");
    if(!stream) { perror(argv[1]); return 2; }
    int header[4];
    if(fread(header,sizeof(int),4,stream)!=4 ||
       header[0]!=0x4B574550 || header[1]!=55351 ||
       header[2]!=297114 || header[3]!=16) {
        fprintf(stderr,"bad energy-neutral CSR header\n"); return 2;
    }
    const int n=header[1],nf=header[2];
    GPUContext ctx;
    memset(&ctx,0,sizeof(ctx));
    ctx.initialized=1;
    ctx.d_faces.n_particles=n;
    ctx.d_faces.n_faces_total=nf;
    ctx.d_faces.face_offset=read_managed<int>(stream,n+1);
    ctx.d_faces.neighbor_idx=read_managed<int>(stream,nf);
    ctx.d_faces.is_ghost=read_managed<int>(stream,nf);
    ctx.d_faces.c1x=read_managed<double>(stream,nf);
    ctx.d_faces.c1y=read_managed<double>(stream,nf);
    ctx.d_faces.c2x=read_managed<double>(stream,nf);
    ctx.d_faces.c2y=read_managed<double>(stream,nf);
    ctx.d_parts.x=read_managed<double>(stream,n);
    ctx.d_parts.y=read_managed<double>(stream,n);
    ctx.d_parts.vx=read_managed<double>(stream,n);
    ctx.d_parts.vy=read_managed<double>(stream,n);
    ctx.d_parts.pressure=read_managed<float>(stream,n);
    unsigned char *active=read_managed<unsigned char>(stream,n);
    double *reference_rate=read_managed<double>(stream,n);
    double *reference_gradient=read_managed<double>(stream,n);
    double *reference_projected=read_managed<double>(stream,n);
    if(fgetc(stream)!=EOF || ctx.d_faces.face_offset[n]!=nf) {
        fprintf(stderr,"bad CSR length/offset\n"); return 2;
    }
    fclose(stream);
    double *projected=(double*)malloc((size_t)n*sizeof(double));
    double *gradient=(double*)malloc((size_t)n*sizeof(double));
    if(!projected || !gradient) return 2;
    unsigned long long unsupported=999;
    double alpha=NAN;
    auto start=std::chrono::steady_clock::now();
    int status=gpu_probe_energy_neutral_rate(&ctx,active,n,n,16,.8,
                                             projected,gradient,&alpha,&unsupported);
    auto end=std::chrono::steady_clock::now();
    double seconds=std::chrono::duration<double>(end-start).count();
    double max_rate=0.,max_gradient=0.,max_projected=0.;
    double dot=0.,norm=0.,reference_dot=0.,reference_norm=0.;
    int nonfinite=0,active_count=0;
    for(int i=0;i<n;++i) {
        if(!active[i]) continue;
        ++active_count;
        if(!std::isfinite(projected[i]) || !std::isfinite(gradient[i])) {
            ++nonfinite; continue;
        }
        max_gradient=fmax(max_gradient,fabs(gradient[i]-reference_gradient[i]));
        max_projected=fmax(max_projected,fabs(projected[i]-reference_projected[i]));
        dot+=gradient[i]*projected[i];
        norm+=gradient[i]*gradient[i];
        reference_dot+=reference_gradient[i]*reference_projected[i];
        reference_norm+=reference_gradient[i]*reference_gradient[i];
        double reconstructed_rate=projected[i]+alpha*gradient[i];
        max_rate=fmax(max_rate,fabs(reconstructed_rate-reference_rate[i]));
    }
    double relative_orthogonality=fabs(dot)/sqrt(norm);
    double reference_relative_orthogonality=fabs(reference_dot)/sqrt(reference_norm);
    printf("case=%s active=%d faces=%d status=%d unsupported=%llu nonfinite=%d elapsed_seconds=%.9g max_rate=%.17g max_gradient=%.17g max_projected=%.17g relative_orthogonality=%.17g reference_relative_orthogonality=%.17g alpha=%.17g\n",
           argv[1],active_count,nf,status,unsupported,nonfinite,seconds,
           max_rate,max_gradient,max_projected,relative_orthogonality,
           reference_relative_orthogonality,alpha);
    if(status || unsupported || nonfinite || active_count!=49038 ||
       max_rate>1e-10 || max_gradient>1e-10 || max_projected>1e-10 ||
       relative_orthogonality>1e-12) return 3;
    return 0;
}
