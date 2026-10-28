#include "gfs_spin_face_exchange_20261027.h"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>

struct result {
    gfs_spin_cell a,b,invalid_a,invalid_b;
    int status,invalid_status;
};

__global__ void test_kernel(result *out) {
    gfs_spin_cell a={0.,0.,2.,.4,.2,10.,.01,.2};
    gfs_spin_cell b={1.,0.,1.5,-.1,.3,8.,-.02,.15};
    const double dm=.02;
    gfs_spin_face_transfer f={.45,.2,dm,.008,-.004,.03,
                              a.inertia*(a.m-dm)/a.m,b.inertia*(b.m+dm)/b.m};
    out->status=gfs_spin_face_exchange(&a,&b,f);
    out->a=a;out->b=b;
    gfs_spin_cell x={0.,0.,2.,.4,.2,10.,.01,.2};
    gfs_spin_cell y={1.,0.,1.5,-.1,.3,8.,-.02,.15};
    f.dm=3.;
    out->invalid_status=gfs_spin_face_exchange(&x,&y,f);
    out->invalid_a=x;out->invalid_b=y;
}

static double angular(const gfs_spin_cell &a,const gfs_spin_cell &b) {
    return a.x*a.py-a.y*a.px+a.spin+b.x*b.py-b.y*b.px+b.spin;
}

static bool close(double a,double b) { return std::fabs(a-b)<1e-12; }

int main() {
    result *device=nullptr,host{};
    cudaError_t err=cudaMalloc((void**)&device,sizeof(result));
    if(err!=cudaSuccess) { std::fprintf(stderr,"cudaMalloc: %s\n",cudaGetErrorString(err));return 1; }
    test_kernel<<<1,1>>>(device);
    err=cudaGetLastError();
    if(err==cudaSuccess)err=cudaDeviceSynchronize();
    if(err==cudaSuccess)err=cudaMemcpy(&host,device,sizeof(result),cudaMemcpyDeviceToHost);
    cudaFree(device);
    if(err!=cudaSuccess) { std::fprintf(stderr,"CUDA test: %s\n",cudaGetErrorString(err));return 1; }
    const gfs_spin_cell a0={0.,0.,2.,.4,.2,10.,.01,.2};
    const gfs_spin_cell b0={1.,0.,1.5,-.1,.3,8.,-.02,.15};
    const gfs_spin_cell &a=host.a,&b=host.b;
    bool ok=host.status==0 && host.invalid_status==-2 &&
        close(a.m+b.m,a0.m+b0.m) && close(a.px+b.px,a0.px+b0.px) &&
        close(a.py+b.py,a0.py+b0.py) && close(a.energy+b.energy,a0.energy+b0.energy) &&
        close(a.m*a.x+b.m*b.x,a0.m*a0.x+b0.m*b0.x) &&
        close(a.m*a.y+b.m*b.y,a0.m*a0.y+b0.m*b0.y) &&
        close(angular(a,b),angular(a0,b0)) &&
        gfs_spin_thermal_energy(a)>0. && gfs_spin_thermal_energy(b)>0. &&
        close(host.invalid_a.m,a0.m) && close(host.invalid_b.m,b0.m) &&
        close(host.invalid_a.spin,a0.spin) && close(host.invalid_b.spin,b0.spin);
    if(!ok) { std::fprintf(stderr,"FAIL: device spin face exchange invariants\n");return 1; }
    std::printf("PASS: CUDA pair mass, first moment, momentum, energy, spin and angular momentum\n");
    return 0;
}
