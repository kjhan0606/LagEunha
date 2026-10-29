#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "Exam/gfs_weight_exchange_gpu.cuh"

static void ck(cudaError_t e){if(e!=cudaSuccess){fprintf(stderr,"CUDA: %s\n",cudaGetErrorString(e));exit(1);}}
static void near(double a,double b){if(!std::isfinite(a)||fabs(a-b)>1e-12*(1+fabs(b))){fprintf(stderr,"%.17g != %.17g\n",a,b);exit(2);}}
int main(){
    const int nc=4,nf=6;
    const double ox[nc]={0.,1.,0.,1.},oy[nc]={0.,0.,1.,1.};
    const int ii[nf]={0,0,0,1,1,2},jj[nf]={1,2,3,2,3,3};
    GfsWeightExchangeFace f[nf]={};
    GfsMomentState expected[nc]={};
    for(int k=0;k<nf;k++){
        f[k].i=ii[k];f[k].j=jj[k];f[k].wdot_i=.02*(1+ii[k]);f[k].wdot_j=.02*(1+jj[k]);
        f[k].distance=std::hypot(ox[ii[k]]-ox[jj[k]],oy[ii[k]]-oy[jj[k]]);
        f[k].rho=1+.1*k;f[k].ux=.2-.03*k;f[k].uy=-.1+.04*k;f[k].specific_energy=2+.2*k;
        f[k].ax=.2+.11*k;f[k].ay=.5-.07*k;f[k].bx=f[k].ax+.2;f[k].by=f[k].ay+.3;
        f[k].ox_i=ox[ii[k]];f[k].oy_i=oy[ii[k]];f[k].ox_j=ox[jj[k]];f[k].oy_j=oy[jj[k]];
        double flux[4];
        if(!gfs_weight_exchange_flux(f[k].wdot_i,f[k].wdot_j,f[k].distance,f[k].rho,f[k].ux,f[k].uy,f[k].specific_energy,flux) ||
           !gfs_weight_exchange_segment(&expected[ii[k]],&expected[jj[k]],f[k].ax,f[k].ay,f[k].bx,f[k].by,
                f[k].ox_i,f[k].oy_i,f[k].ox_j,f[k].oy_j,flux))return 3;
    }
    GfsWeightExchangeFace *df;GfsMomentState *dr;int *invalid;
    ck(cudaMalloc(&df,sizeof(f)));ck(cudaMalloc(&dr,sizeof(expected)));ck(cudaMalloc(&invalid,sizeof(int)));
    ck(cudaMemcpy(df,f,sizeof(f),cudaMemcpyHostToDevice));ck(cudaMemset(dr,0,sizeof(expected)));ck(cudaMemset(invalid,0,sizeof(int)));
    gfs_weight_exchange_faces_kernel<<<1,32>>>(df,nf,nc,dr,invalid);ck(cudaGetLastError());
    GfsMomentState got[nc];int bad;
    ck(cudaMemcpy(got,dr,sizeof(got),cudaMemcpyDeviceToHost));ck(cudaMemcpy(&bad,invalid,sizeof(bad),cudaMemcpyDeviceToHost));
    if(bad)return 4;
    double total_mass=0,total_px=0,total_py=0,total_energy=0,total_L=0;
    for(int k=0;k<nc;k++){
        for(int q=0;q<GFM_NSTATE;q++)near(got[k].a[q],expected[k].a[q]);
        total_mass+=got[k].a[GFM_M];total_px+=got[k].a[GFM_PX];total_py+=got[k].a[GFM_PY];
        total_energy+=got[k].a[GFM_E];total_L+=gfs_moment_Lz(&got[k],ox[k],oy[k]);
    }
    near(total_mass,0);near(total_px,0);near(total_py,0);near(total_energy,0);near(total_L,0);
    for(int k=0;k<nf;k++)f[k].wdot_i=f[k].wdot_j=0;
    ck(cudaMemcpy(df,f,sizeof(f),cudaMemcpyHostToDevice));ck(cudaMemset(dr,0,sizeof(expected)));ck(cudaMemset(invalid,0,sizeof(int)));
    gfs_weight_exchange_faces_kernel<<<1,32>>>(df,nf,nc,dr,invalid);ck(cudaGetLastError());
    ck(cudaMemcpy(got,dr,sizeof(got),cudaMemcpyDeviceToHost));ck(cudaMemcpy(&bad,invalid,sizeof(bad),cudaMemcpyDeviceToHost));
    if(bad)return 5;
    for(int k=0;k<nc;k++)for(int q=0;q<GFM_NSTATE;q++)near(got[k].a[q],0);
    f[1].j=nc;
    ck(cudaMemcpy(df,f,sizeof(f),cudaMemcpyHostToDevice));ck(cudaMemset(dr,0,sizeof(expected)));ck(cudaMemset(invalid,0,sizeof(int)));
    gfs_weight_exchange_faces_kernel<<<1,32>>>(df,nf,nc,dr,invalid);ck(cudaGetLastError());
    ck(cudaMemcpy(&bad,invalid,sizeof(bad),cudaMemcpyDeviceToHost));if(bad!=1)return 6;
    ck(cudaFree(df));ck(cudaFree(dr));ck(cudaFree(invalid));
    puts("GPU/CPU face contributions, conservation, Wdot=0 and invalid-face checks pass");
}
