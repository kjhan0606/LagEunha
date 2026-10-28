#include "gfs_weight_motion_remap_20261027.h"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>

__device__ bool near(double x,double y){return fabs(x-y)<1e-12*(1.+fabs(y));}
__device__ double angular(const gfs_spin_affine_cell &a){
    const auto &s=a.mechanical;return s.x*s.py-s.y*s.px+s.spin;
}
__global__ void test_kernel(int *result){
    gfs_spin_affine_cell a{{-.5,0.,1.,.4,0.,2.08,0.,1.},0.};
    gfs_spin_affine_cell b{{.5,0.,1.,.4,0.,2.08,0.,1.},0.};
    const auto olda=a,oldb=b;
    gfs_weight_motion_input q{0.,0.,1.,1.,0.,0.,.2,0.,0.,1.,1.,1.,
                              {1.,.4,0.,2.08,0.},{1.,.4,0.,2.08,0.}};
    gfs_spin_affine_transfer f;
    if(!gfs_weight_motion_transfer(q,&f) || !near(f.mechanical.dm,-.1)) {*result=1;return;}
    if(gfs_spin_affine_face_exchange(&a,&b,f)) {*result=2;return;}
    if(!near(a.mechanical.m/1.1,1.) || !near(b.mechanical.m/.9,1.) ||
       !near(gfs_spin_affine_thermal(a)/a.mechanical.m,2.) ||
       !near(gfs_spin_affine_thermal(b)/b.mechanical.m,2.) ||
       !near(a.mechanical.energy+b.mechanical.energy,olda.mechanical.energy+oldb.mechanical.energy) ||
       !near(angular(a)+angular(b),angular(olda)+angular(oldb))) {*result=3;return;}
    q.wdot_i=0.;q.wdot_j=0.;
    if(!gfs_weight_motion_transfer(q,&f) || f.mechanical.dm!=0. ||
       f.mechanical.dpx!=0. || f.mechanical.denergy!=0.) {*result=4;return;}
    q.wdot_i=-.2;q.wdot_j=.2;q.state_i.density=2.;
    if(!gfs_weight_motion_transfer(q,&f) || !near(f.mechanical.dm,.4)) {*result=5;return;}
    const auto oldf=f;q.site_distance=0.;
    if(gfs_weight_motion_transfer(q,&f) || f.mechanical.dm!=oldf.mechanical.dm) {*result=6;return;}
    q.site_distance=1.;q.wdot_i=0.;q.wdot_j=0.;q.weight_i=.2;q.weight_j=0.;q.dt=.1;
    q.normal_separation_rate=.4;q.state_i.density=2.;
    if(!gfs_weight_motion_transfer(q,&f) || !near(f.mechanical.dm,.008)) {*result=7;return;}
    q.weight_i=q.weight_j=0.;
    if(!gfs_weight_motion_transfer(q,&f) || f.mechanical.dm!=0.) {*result=8;return;}
    *result=0;
}
int main(){
    int *device=nullptr,status=-1;
    if(cudaMalloc(&device,sizeof(int))!=cudaSuccess)return 10;
    test_kernel<<<1,1>>>(device);
    if(cudaGetLastError()!=cudaSuccess || cudaDeviceSynchronize()!=cudaSuccess ||
       cudaMemcpy(&status,device,sizeof(int),cudaMemcpyDeviceToHost)!=cudaSuccess)return 11;
    cudaFree(device);
    if(status){std::fprintf(stderr,"device test failed: %d\n",status);return 12;}
    std::puts("PASS: GPU full weight-motion face-speed, mass/energy/angular momentum exchange and W0 reduction");
}
