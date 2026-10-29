/* Frozen-geometry audit bridge; calls the actual candidate face operator. */
#include "Exam/gfs_weight_mass_exchange.h"
int affine_batch(int n,const double *di,const double *dj,const double *distance,
    const double *ax,const double *ay,const double *bx,const double *by,
    double *out_i,double *out_j){
    const double g[4]={.31,-.72,.45,-.12}, b[2]={.7,-.4};
    for(int k=0;k<n;k++){
        double x=.5*(ax[k]+bx[k]),y=.5*(ay[k]+by[k]);
        GfsMomentState i,j;gfs_moment_zero(&i);gfs_moment_zero(&j);
        if(!gfs_weight_exchange_segment_affine(&i,&j,ax[k],ay[k],bx[k],by[k],
            0,0,0,0,di[k],dj[k],distance[k],1,
            b[0]+g[0]*x+g[1]*y,b[1]+g[2]*x+g[3]*y,g,2))return 0;
        for(int q=0;q<GFM_NSTATE;q++){out_i[k*GFM_NSTATE+q]=i.a[q];out_j[k*GFM_NSTATE+q]=j.a[q];}
    }
    return 1;
}
