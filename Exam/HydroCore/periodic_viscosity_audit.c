#include <stdio.h>
#include <math.h>
#include "../gfs_viscosity_implicit.h"
/* An interaction across a periodic cut conserves energy and linear
 * momentum. Raw-coordinate angular momentum has a cut flux. */
int main(void){
    int failed=0;
    for(int dim=2;dim<=3;dim++){
        double mass[2]={.3,.7},U[2]={.01,.02},x[6]={0},v[6]={0};
        x[0]=.98;x[1]=.2;x[dim]=.02;x[dim+1]=.25;
        if(dim==3){x[2]=.9;x[5]=.1;}
        v[0]=1;v[dim]=-1;
        double shift[3]={1,0,dim==3?1:0},r[3]={0},length=0;
        for(int a=0;a<dim;a++){r[a]=x[dim+a]+shift[a]-x[a];length+=r[a]*r[a];}
        length=sqrt(length);
        GfsViscosityEdge edge={0,1,{0,0,0},10,{shift[0],shift[1],shift[2]}};
        for(int a=0;a<dim;a++)edge.normal[a]=r[a]/length;
        GfsViscosityWorkspace *w=gfs_viscosity_workspace_create();GfsViscosityStats s;
        int rc=gfs_viscosity_be(w,dim,2,1,&edge,mass,x,v,U,.1,&s);
        gfs_viscosity_workspace_free(w);if(rc)return 2;
        double dp=0,dJ=0,flux=0;
        for(int a=0;a<3;a++){dp=fmax(dp,fabs(s.momentum_defect[a]));dJ=fmax(dJ,fabs(s.angular_defect[a]));flux=fmax(flux,fabs(s.periodic_torque[a]));}
        int pass=fabs(s.energy_defect)<1e-12&&dp<1e-12&&dJ<1e-12&&flux>1e-3&&s.min_heat>0;
        printf("{\"dim\":%d,\"energy_defect\":%.17g,\"momentum_defect\":%.17g,\"angular_defect_with_periodic_flux\":%.17g,\"nonzero_periodic_torque_flux\":%.17g,\"pass\":%s}\n",dim,s.energy_defect,dp,dJ,flux,pass?"true":"false");
        failed|=!pass;
    }
    return failed;
}
