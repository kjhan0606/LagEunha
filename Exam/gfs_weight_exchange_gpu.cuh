#ifndef GFS_WEIGHT_EXCHANGE_GPU_CUH
#define GFS_WEIGHT_EXCHANGE_GPU_CUH
#include "gfs_weight_mass_exchange.h"

/* Experimental interface for a prevalidated list of unique interior faces.
 * The caller owns construction, orientation, halo identity, MPI face owner,
 * state/weight staging, and accumulation into a disposable RHS buffer.
 * This kernel is not wired to the production solver. */
struct GfsWeightExchangeFace {
    int i,j;
    double wdot_i,wdot_j,distance,rho,ux,uy,specific_energy;
    double ax,ay,bx,by,ox_i,oy_i,ox_j,oy_j;
};

__global__ static void gfs_weight_exchange_faces_kernel(
    const GfsWeightExchangeFace *faces,int nface,int ncell,
    GfsMomentState *rhs,int *invalid){
    int k=blockIdx.x*blockDim.x+threadIdx.x;
    if(k>=nface)return;
    GfsWeightExchangeFace f=faces[k];
    if(f.i<0 || f.j<0 || f.i>=ncell || f.j>=ncell || f.i==f.j){atomicAdd(invalid,1);return;}
    double flux[4];GfsMomentState di,dj;
    gfs_moment_zero(&di);gfs_moment_zero(&dj);
    if(!gfs_weight_exchange_flux(f.wdot_i,f.wdot_j,f.distance,f.rho,f.ux,f.uy,f.specific_energy,flux) ||
       !gfs_weight_exchange_segment(&di,&dj,f.ax,f.ay,f.bx,f.by,f.ox_i,f.oy_i,f.ox_j,f.oy_j,flux)){
        atomicAdd(invalid,1);return;
    }
    for(int q=0;q<GFM_NSTATE;q++)if(!isfinite(di.a[q]) || !isfinite(dj.a[q])){
        atomicAdd(invalid,1);return;
    }
    for(int q=0;q<GFM_NSTATE;q++){
        atomicAdd(&rhs[f.i].a[q],di.a[q]);
        atomicAdd(&rhs[f.j].a[q],dj.a[q]);
    }
}
#endif
