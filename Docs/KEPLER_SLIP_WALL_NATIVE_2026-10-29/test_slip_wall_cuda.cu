/* Standalone primitive test, not a fluid evolution or RK4 event test. */
#include "Exam/kepler_slip_wall.h"
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>

struct Input { double x,y,vx,vy; };

__global__ void impacts(const Input *in, GfsSlipImpulse *out, int *status, int n) {
    const int i=blockIdx.x*blockDim.x+threadIdx.x;
    if (i<n) status[i]=gfs_slip_circle_impact(in[i].x,in[i].y,3.,3.,2.75,
                                             .017,in[i].vx,in[i].vy,out+i);
}

static void check(cudaError_t e) {
    if(e!=cudaSuccess) {fprintf(stderr,"CUDA test failed: %s\n",cudaGetErrorString(e));exit(2);}
}

static bool near(double a,double b) {
    return fabs(a-b)<=256.*DBL_EPSILON*fmax(1.,fmax(fabs(a),fabs(b)));
}

int main() {
    constexpr int n=4096;
    std::vector<Input> in(n);
    std::vector<GfsSlipImpulse> got(n);
    std::vector<int> status(n);
    for(int i=0;i<n;++i) {
        const double t=(i+.37)*6.2831853071795864769/n,nx=cos(t),ny=sin(t);
        const double vnormal=(i%2 ? -.43 : .43);
        in[i]={3.+2.75*nx,3.+2.75*ny,vnormal*nx-.71*ny,vnormal*ny+.71*nx};
    }
    Input *di; GfsSlipImpulse *dout; int *ds;
    check(cudaMalloc(&di,n*sizeof(Input)));
    check(cudaMalloc(&dout,n*sizeof(GfsSlipImpulse)));
    check(cudaMalloc(&ds,n*sizeof(int)));
    check(cudaMemcpy(di,in.data(),n*sizeof(Input),cudaMemcpyHostToDevice));
    impacts<<<(n+127)/128,128>>>(di,dout,ds,n);
    check(cudaGetLastError());check(cudaDeviceSynchronize());
    check(cudaMemcpy(got.data(),dout,n*sizeof(GfsSlipImpulse),cudaMemcpyDeviceToHost));
    check(cudaMemcpy(status.data(),ds,n*sizeof(int),cudaMemcpyDeviceToHost));
    double worst_dk=0.,worst_torque=0.;
    for(int i=0;i<n;++i) {
        GfsSlipImpulse expected;
        const int s=gfs_slip_circle_impact(in[i].x,in[i].y,3.,3.,2.75,.017,in[i].vx,in[i].vy,&expected);
        const auto &a=got[i];
        const double rx=in[i].x-3.,ry=in[i].y-3.;
        const bool pass=status[i]==s && status[i]>=0 && near(a.vx_after,expected.vx_after)
            && near(a.vy_after,expected.vy_after) && near(a.wall_jx,expected.wall_jx)
            && near(a.wall_jy,expected.wall_jy) && near(a.fluid_kinetic_roundoff,0.)
            && near(a.wall_torque_impulse,0.)
            && near(.017*(a.vx_after-in[i].vx)+a.wall_jx,0.)
            && near(.017*(a.vy_after-in[i].vy)+a.wall_jy,0.)
            && near(rx*a.vy_after-ry*a.vx_after,rx*in[i].vy-ry*in[i].vx);
        if(!pass) {fprintf(stderr,"Primitive parity/conservation failure at %d\n",i);return 1;}
        worst_dk=fmax(worst_dk,fabs(a.fluid_kinetic_roundoff));
        worst_torque=fmax(worst_torque,fabs(a.wall_torque_impulse));
    }
    cudaDeviceProp props;check(cudaGetDeviceProperties(&props,0));
    printf("PASS CUDA primitive: GPU=%s cases=%d max_kinetic_roundoff=%.17g max_wall_torque=%.17g\n",
           props.name,n,worst_dk,worst_torque);
    check(cudaFree(di));check(cudaFree(dout));check(cudaFree(ds));
    return 0;
}
