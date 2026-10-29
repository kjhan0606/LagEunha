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
