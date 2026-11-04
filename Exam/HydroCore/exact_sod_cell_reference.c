/* Exact Riemann reference integrated on requested actual cell intervals.
 * Eight-point Gauss is applied separately between all physical waves.
 * This is an analytic reference, not a PDE or production state update. */
#define main archived_hydro1d_main
#include "../../Hydro1DExam/laguerre_sod.c"
#undef main
int main(int argc,char **argv){
    if(argc!=3)return 2;double t=atof(argv[1]),gamma=atof(argv[2]);if(!(t>0)||!(gamma>1))return 2;
    setup_problem(PROB_SOD);GAM=gamma;GM1=gamma-1;GP1=gamma+1;exact_riemann(&ps_exact,&vs_exact);
    double cl=eos_c(rhoL_ic,pL_ic),cr=eos_c(rhoR_ic,pR_ic),cs=cl*pow(ps_exact/pL_ic,GM1/(2*GAM));
    double waves[4]={xmid1-cl*t,xmid1+(vs_exact-cs)*t,xmid1+vs_exact*t,
        xmid1+cr*sqrt(GP1/(2*GAM)*ps_exact/pR_ic+GM1/(2*GAM))*t};
    const double node[8]={-.9602898564975363,-.7966664774136267,-.5255324099163290,-.1834346424956498,.1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363};
    const double weight[8]={.1012285362903763,.2223810344533745,.3137066458778873,.3626837833783620,.3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763};
    double a,b;
    while(scanf("%lf %lf",&a,&b)==2){
        if(!(b>a))return 2;double mass=0,momentum=0,pressure=0,kinetic=0,left=a;
        for(int segment=0;segment<5;segment++){
            double right=segment<4?fmin(b,waves[segment]):b;if(right<=left)continue;
            for(int q=0;q<8;q++){
                double x=.5*(left+right)+.5*(right-left)*node[q],rho,v,p;eval_solution(x,t,&rho,&v,&p);
                double w=.5*(right-left)*weight[q];mass+=w*rho;momentum+=w*rho*v;pressure+=w*p;kinetic+=.5*w*rho*v*v;
            }
            left=right;if(left>=b)break;
        }
        printf("%.17g,%.17g,%.17g,%.17g\n",mass/(b-a),momentum/mass,pressure/(b-a),kinetic-.5*momentum*momentum/mass);
    }
    return 0;
}
