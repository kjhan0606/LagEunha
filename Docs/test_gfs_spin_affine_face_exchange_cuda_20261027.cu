#include "gfs_spin_affine_face_exchange_20261027.h"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>

__device__ double angular(const gfs_spin_affine_cell &a) {
    return a.mechanical.x*a.mechanical.py-a.mechanical.y*a.mechanical.px+a.mechanical.spin;
}

__global__ void test_kernel(int *result) {
    gfs_spin_affine_cell a{{-1.,0.,1.,.1,.2,2.,.1,1.2},.25};
    gfs_spin_affine_cell b{{1.,0.,1.3,-.1,.05,2.3,-.2,1.6},.3};
    const auto x=a,y=b;
    gfs_spin_affine_transfer f{{.1,.4,.03,.05,.03,.02,1.2,1.6},-.02,.02};
    int status=gfs_spin_affine_face_exchange(&a,&b,f);
    if(!status && (fabs(a.mechanical.m+b.mechanical.m-x.mechanical.m-y.mechanical.m)>1e-12 ||
       fabs(a.mechanical.energy+b.mechanical.energy-x.mechanical.energy-y.mechanical.energy)>1e-12 ||
       fabs(angular(a)+angular(b)-angular(x)-angular(y))>1e-12 ||
       !(gfs_spin_affine_thermal(a)>0.) || !(gfs_spin_affine_thermal(b)>0.)))status=1;
    if(!status) {
        const double before=a.mechanical.energy;f.nonrigid_delta_i=100.;
        if(gfs_spin_affine_face_exchange(&a,&b,f)==0 || a.mechanical.energy!=before)status=2;
    }
    *result=status;
}

int main() {
    int *device=nullptr,status=-1;
    if(cudaMalloc(&device,sizeof(int))!=cudaSuccess)return 10;
    test_kernel<<<1,1>>>(device);
    if(cudaGetLastError()!=cudaSuccess || cudaDeviceSynchronize()!=cudaSuccess ||
       cudaMemcpy(&status,device,sizeof(int),cudaMemcpyDeviceToHost)!=cudaSuccess)return 11;
    cudaFree(device);
    if(status){std::fprintf(stderr,"device test failed: %d\n",status);return 12;}
    std::puts("PASS: GPU nonrigid/spin/thermal pair conservation and rollback");
}
