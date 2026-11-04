#include <stdio.h>
#include <math.h>
#include "../gfs_variational_core.h"

static double norm(double v[3]) {
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}
int main(void) {
    int fail=0;
    double worst_shear=0,worst_energy=0,worst_torque=0,worst_frame=0;
    int checks=0;
    for(int dim=1;dim<=3;dim++) {
        double r[3]={.7,dim>1?.4:0,dim>2?-.3:0},c[3];
        double vi[3]={.3,-.8,.5},vj[3]={0},A[3][3]={{-.2,2,.4},{-.8,.1,-.5},{.3,.7,.1}};
        for(int a=0;a<3;a++){
            if(a>=dim)vi[a]=0;
            c[a]=.5*r[a];vj[a]=vi[a];
            for(int b=0;b<dim&&a<dim;b++)vj[a]+=A[a][b]*r[b];
        }
        GfsVariationalFace f;
        if(!gfs_variational_face(dim,r,c,1,&f))return 2;
        double dvface[3]={0};
        /* Both affine reconstructions reach exactly the same point. */
        for(int a=0;a<dim;a++){
            double left=vi[a],right=vj[a];
            for(int b=0;b<dim;b++){left+=A[a][b]*c[b];right+=A[a][b]*(c[b]-r[b]);}
            dvface[a]=right-left;
        }
        double q=gfs_variational_pair_pressure(&f,vi,vj,dvface,1,1,.125,sqrt(8),0,1,1);
        worst_shear=fmax(worst_shear,q);checks++;
        if(q>1e-14)fail=1;

        /* Unresolved approach: zero limited slopes, exact opposite impulses. */
        for(int a=0;a<dim;a++){vi[a]=f.normal[a];vj[a]=-f.normal[a];dvface[a]=vj[a]-vi[a];}
        q=gfs_variational_pair_pressure(&f,vi,vj,dvface,1,1,.125,sqrt(8),0,1,1);
        double F[3]={0},cross[3],power=0,closing=0;
        for(int a=0;a<dim;a++){
            F[a]=-q*f.area*f.normal[a];
            power+=F[a]*(vi[a]-vj[a]);
            closing+=(vi[a]-vj[a])*f.normal[a];
        }
        double heat=q*f.area*closing;
        cross[0]=r[1]*F[2]-r[2]*F[1];
        cross[1]=r[2]*F[0]-r[0]*F[2];
        cross[2]=r[0]*F[1]-r[1]*F[0];
        worst_energy=fmax(worst_energy,fabs(power+heat));
        worst_torque=fmax(worst_torque,norm(cross));
        checks+=3;
        if(!(q>0)||!(heat>0)||fabs(power+heat)>1e-14||norm(cross)>1e-14)fail=1;
        for(int a=0;a<dim;a++){vi[a]+=123.456;vj[a]+=123.456;}
        double qboost=gfs_variational_pair_pressure(&f,vi,vj,dvface,1,1,.125,sqrt(8),0,1,1);
        worst_frame=fmax(worst_frame,fabs(q-qboost));checks++;
        if(fabs(q-qboost)>64*2.220446049250313e-16*fmax(1.0,q))fail=1;
        for(int a=0;a<dim;a++){vi[a]=-f.normal[a];vj[a]=f.normal[a];}
        q=gfs_variational_pair_pressure(&f,vi,vj,dvface,1,1,.125,sqrt(8),0,1,1);
        checks++;if(q!=0)fail=1;
    }
    printf("{\"checks\":%d,\"max_affine_pair_pressure\":%.17g,\"max_energy_rhs_residual\":%.17g,\"max_pair_torque\":%.17g,\"max_galilean_pressure_difference\":%.17g,\"scope\":\"algebra only; slope fitting, limiters, shocks and Kepler not qualified\",\"pass\":%s}\n",
        checks,worst_shear,worst_energy,worst_torque,worst_frame,fail?"false":"true");
    return fail?1:0;
}
