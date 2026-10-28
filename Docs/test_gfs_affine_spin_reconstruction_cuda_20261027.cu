#include "gfs_affine_spin_reconstruction_20261027.h"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>

__global__ void test_kernel(int *result) {
    gfs_affine_spin_input x{2.,.3,1.,.2,-1.1,.4,.6,0.};
    gfs_affine_spin_output a;
    int status=0;
    if(!gfs_affine_spin_reconstruct(x,&a))status=1;
    if(!status) {
        const double residual=a.ayx*x.qxx-a.axy*x.qyy+(a.ayy-a.axx)*x.qxy;
        if(fabs(residual)>1e-12 || a.spin_energy!=0.)status=2;
    }
    if(!status) {
        x.spin=3.1;
        if(!gfs_affine_spin_reconstruct(x,&a))status=3;
        const double spin=a.ayx*x.qxx-a.axy*x.qyy+(a.ayy-a.axx)*x.qxy;
        const double energy=gfs_affine_energy(a.axx,a.axy,a.ayx,a.ayy,x.qxx,x.qxy,x.qyy);
        if(fabs(spin-x.spin)>1e-12 || fabs(energy-a.nonrigid_energy-a.spin_energy)>1e-12)status=4;
    }
    if(!status) {
        x.spin=0.;
        if(!gfs_affine_spin_reconstruct_excess(x,&a) ||
           a.axx!=x.axx || a.axy!=x.axy || a.ayx!=x.ayx || a.ayy!=x.ayy)status=6;
        x.spin=.7;
        if(!gfs_affine_spin_reconstruct_excess(x,&a))status=7;
        const double spin=a.ayx*x.qxx-a.axy*x.qyy+(a.ayy-a.axx)*x.qxy;
        if(fabs(spin-a.reference_spin-x.spin)>1e-12)status=8;
    }
    if(!status) {
        const double previous=a.axy;
        x.qxy=2.;
        if(gfs_affine_spin_reconstruct(x,&a) || a.axy!=previous)status=5;
    }
    *result=status;
}

int main() {
    int *device=nullptr,status=-1;
    if(cudaMalloc(&device,sizeof(int))!=cudaSuccess)return 10;
    test_kernel<<<1,1>>>(device);
    if(cudaGetLastError()!=cudaSuccess || cudaDeviceSynchronize()!=cudaSuccess ||
       cudaMemcpy(&status,device,sizeof(int),cudaMemcpyDeviceToHost)!=cudaSuccess) return 11;
    cudaFree(device);
    if(status){std::fprintf(stderr,"device test failed: %d\n",status);return status;}
    std::puts("PASS: GPU affine spin projection, energy split and invalid rollback");
}
