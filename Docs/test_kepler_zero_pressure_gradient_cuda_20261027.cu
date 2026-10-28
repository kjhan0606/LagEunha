#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "gfs_energy_neutral_probe.cu"

static void check(cudaError_t e,const char *where) {
    if(e!=cudaSuccess) {
        std::fprintf(stderr,"%s: %s\n",where,cudaGetErrorString(e));
        std::exit(1);
    }
}

template<class T> static T *device_copy(const T *host,size_t n) {
    T *device=nullptr;
    check(cudaMalloc((void**)&device,n*sizeof(T)),"cudaMalloc");
    check(cudaMemcpy(device,host,n*sizeof(T),cudaMemcpyHostToDevice),"cudaMemcpy");
    return device;
}

int main(int argc,char **argv) {
    if(argc!=2 || (std::strcmp(argv[1],"old") && std::strcmp(argv[1],"fixed"))) {
        std::fprintf(stderr,"Expected old or fixed\n");return 2;
    }
    const int offsets[3]={0,1,2},neighbors[2]={1,0},ghost[2]={0,0};
    const double x[2]={0.,1.},y[2]={0.,0.},vx[2]={0.,0.},vy[2]={0.,1.};
    const double c1x[2]={.5,-.5},c2x[2]={.5,-.5};
    const double c1y[2]={-.3,-.3},c2y[2]={.7,.7};
    const float pressure[2]={1.e-6f,1.e-6f};
    const unsigned char active[2]={1,1};
    GPUContext ctx{};
    ctx.initialized=1;ctx.stream=nullptr;
    ctx.d_faces.n_particles=2;ctx.d_faces.n_faces_total=2;
    ctx.d_faces.face_offset=device_copy(offsets,3);
    ctx.d_faces.neighbor_idx=device_copy(neighbors,2);
    ctx.d_faces.is_ghost=device_copy(ghost,2);
    ctx.d_faces.c1x=device_copy(c1x,2);ctx.d_faces.c2x=device_copy(c2x,2);
    ctx.d_faces.c1y=device_copy(c1y,2);ctx.d_faces.c2y=device_copy(c2y,2);
    ctx.d_parts.x=device_copy(x,2);ctx.d_parts.y=device_copy(y,2);
    ctx.d_parts.vx=device_copy(vx,2);ctx.d_parts.vy=device_copy(vy,2);
    ctx.d_parts.pressure=device_copy(pressure,2);
    double raw[2]={0.,0.},rate[2]={0.,0.},gradient[2]={1.,1.},alpha=42.;
    unsigned long long unsupported=0;
    int raw_status=gpu_probe_volume_rate_jacobi(&ctx,active,2,2,16,.67,raw,&unsupported);
    int status=gpu_probe_energy_neutral_rate(&ctx,active,2,2,16,.67,
                                             rate,gradient,&alpha,&unsupported);
    bool common=raw_status==0 && unsupported==0 && gradient[0]==0. && gradient[1]==0.
                && std::isfinite(raw[0]) && std::isfinite(raw[1])
                && std::fabs(raw[0])+std::fabs(raw[1])>1e-3;
    bool old=std::strcmp(argv[1],"old")==0;
    bool ok=common && (old ? status==-3 :
        (status==0 && alpha==0. && std::memcmp(raw,rate,sizeof raw)==0));
    std::printf("mode=%s status=%d raw_status=%d gradient=(%.17g,%.17g) "
                "raw=(%.17g,%.17g) rate=(%.17g,%.17g) alpha=%.17g\n",
                argv[1],status,raw_status,gradient[0],gradient[1],raw[0],raw[1],
                rate[0],rate[1],alpha);
    cudaFree(ctx.d_parts.pressure);cudaFree(ctx.d_parts.vy);cudaFree(ctx.d_parts.vx);
    cudaFree(ctx.d_parts.y);cudaFree(ctx.d_parts.x);
    cudaFree(ctx.d_faces.c2y);cudaFree(ctx.d_faces.c1y);
    cudaFree(ctx.d_faces.c2x);cudaFree(ctx.d_faces.c1x);
    cudaFree(ctx.d_faces.is_ghost);cudaFree(ctx.d_faces.neighbor_idx);
    cudaFree(ctx.d_faces.face_offset);
    if(!ok) { std::fprintf(stderr,"FAIL: exact uniform-pressure GPU path\n");return 1; }
    std::printf("PASS: %s exact uniform-pressure GPU path\n",argv[1]);
    return 0;
}
