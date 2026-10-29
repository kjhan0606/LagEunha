/* Test bridge: exercise the production candidate's actual scalar flux code. */
#include "Exam/gfs_weight_mass_exchange.h"
int weight_exchange_batch(int n,const double *di,const double *dj,const double *distance,double *flux){
    for(int k=0;k<n;k++){
        double f[4];
        if(!gfs_weight_exchange_flux(di[k],dj[k],distance[k],1,0,0,1,f))return 0;
        flux[k]=f[0];
    }
    return 1;
}
int weight_exchange_segment_batch(int n,const double *di,const double *dj,const double *distance,
    const double *ax,const double *ay,const double *bx,const double *by,double *out_i,double *out_j){
    for(int k=0;k<n;k++){
        double f[4];GfsMomentState i,j;gfs_moment_zero(&i);gfs_moment_zero(&j);
        if(!gfs_weight_exchange_flux(di[k],dj[k],distance[k],1,0,0,1,f) ||
           !gfs_weight_exchange_segment(&i,&j,ax[k],ay[k],bx[k],by[k],0,0,0,0,f))return 0;
        for(int q=0;q<GFM_NSTATE;q++){out_i[k*GFM_NSTATE+q]=i.a[q];out_j[k*GFM_NSTATE+q]=j.a[q];}
    }
    return 1;
}
