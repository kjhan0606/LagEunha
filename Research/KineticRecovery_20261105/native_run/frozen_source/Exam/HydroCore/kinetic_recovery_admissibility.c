/* Exact-state admissibility gate for the specified unlimited-LSQ recovery.
 * No PDE, limiter, floor, energy reset, global solve or evolved reservoir.
 * The fit is the existing dimension-common native face-weighted helper. */
#include <stdio.h>
#include <math.h>
#include "../gfs_velocity_fit.h"

int main(void)
{
    const double gamma=(double)(float)(5./3.),rho=1,P=1,Mach=8;
    const double amplitude=Mach*sqrt(gamma*P/rho);
    double jump_ratio[2];int emitted=0;
    printf("{\"scope\":\"exact-state pressure admissibility of unlimited centroid-LSQ kinetic recovery; zero PDE steps\",\"gamma\":%.17g,\"Mach\":%.17g,\"cases\":[",gamma,Mach);
    for(int level=0;level<2;level++){
        int n=level?32:16;double h=1./n,Q=h*h*h/12,m=rho*h;
        for(int mode=0;mode<3;mode++){
            const char *name=mode==0?"constant_velocity":mode==1?"single_velocity_jump":"affine_integral_control";
            double pmin=INFINITY,pmax=-INFINITY,gmax=0,initial_min=INFINITY,formula_error=0;
            int negative=0;double adjacent_ratio=1;
            for(int i=0;i<n;i++){
                double x=(i+.5)*h,u,ul,ur,true_variance=0;
                if(mode==0)u=ul=ur=amplitude;
                else if(mode==1){
                    u=i<n/2?amplitude:-amplitude;
                    ul=(i-1)<n/2?amplitude:-amplitude;
                    ur=(i+1)<n/2?amplitude:-amplitude;
                }else{
                    u=.2+1.3*x;ul=u-1.3*h;ur=u+1.3*h;
                    true_variance=.5*rho*1.3*1.3*Q;
                }
                GfsVelocityFit fit;gfs_velocity_fit_init(&fit,1);
                double left[3]={-h,0,0},right[3]={h,0,0};
                double vl[3]={ul-u,0,0},vr[3]={ur-u,0,0};
                if(!gfs_velocity_fit_add_face(&fit,left,vl,1)||!gfs_velocity_fit_add_face(&fit,right,vr,1))return 2;
                double G[3][3];if(!gfs_velocity_fit_solve(&fit,G,0))return 2;
                double thermal=P*h/(gamma-1),Kmean=.5*m*u*u;
                double E=thermal+Kmean+true_variance;
                double Kfit=.5*rho*G[0][0]*G[0][0]*Q;
                double Prec=(gamma-1)*((E-Kmean)-Kfit)/h;
                double Pphysical=(gamma-1)*((E-Kmean)-true_variance)/h;
                double expected=P;
                if(mode==1&&(i==n/2-1||i==n/2))expected=P-(gamma-1)*rho*amplitude*amplitude/24;
                formula_error=fmax(formula_error,fabs(Prec-expected));
                initial_min=fmin(initial_min,Pphysical);pmin=fmin(pmin,Prec);pmax=fmax(pmax,Prec);
                gmax=fmax(gmax,fabs(G[0][0]));negative+=Prec<0;
                if(mode==1&&i==n/2-1)adjacent_ratio=Prec/P;
                if(fabs(Pphysical-P)>1e-11)return 2;
            }
            if(formula_error>1e-11)return 2;
            if(mode==1){if(negative!=2||!(pmin<0))return 2;jump_ratio[level]=adjacent_ratio;}
            else if(negative||fabs(pmin-P)>1e-11||fabs(pmax-P)>1e-11)return 2;
            printf("%s{\"cells\":%d,\"state\":\"%s\",\"h\":%.17g,\"Q\":%.17g,\"maximum_abs_native_gradient\":%.17g,\"physical_pressure_min\":%.17g,\"recovered_pressure_min\":%.17g,\"recovered_pressure_max\":%.17g,\"negative_recovery_cells\":%d,\"jump_adjacent_pressure_ratio\":%.17g,\"independent_formula_max_error\":%.17g}",emitted?",":"",n,name,h,Q,gmax,initial_min,pmin,pmax,negative,adjacent_ratio,formula_error);emitted++;
        }
    }
    if(fabs(jump_ratio[1]-jump_ratio[0])>1e-11)return 2;
    printf("],\"predicted_jump_pressure_ratio\":%.17g,\"jump_refinement_ratio_difference\":%.17g,\"decision\":\"REJECT_UNLIMITED_LSQ_KINETIC_RECOVERY\",\"accepted_steps\":0,\"energy_or_momentum_changed\":false,\"new_limiter_or_floor\":false,\"scope_of_rejection\":\"Specified unlimited-gradient algebraic recovery only; no claim against all admissible reconstructions or finite-volume models.\"}\n",1-gamma*(gamma-1)*Mach*Mach/24,jump_ratio[1]-jump_ratio[0]);
    return 0;
}
