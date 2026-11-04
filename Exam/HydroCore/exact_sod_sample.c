#define main archived_hydro1d_main
#include "../../Hydro1DExam/laguerre_sod.c"
#undef main
/* Sample the existing analytic Riemann solver at requested coordinates. */
int main(int argc,char **argv){
    if(argc!=3)return 2;
    double t=atof(argv[1]),gamma=atof(argv[2]);if(!(t>0)||!(gamma>1))return 2;
    setup_problem(PROB_SOD);GAM=gamma;GM1=gamma-1;GP1=gamma+1;
    exact_riemann(&ps_exact,&vs_exact);
    double x;
    while(scanf("%lf",&x)==1){double rho,v,p;eval_solution(x,t,&rho,&v,&p);printf("%.17g,%.17g,%.17g,%.17g\n",x,rho,v,p);}
    return 0;
}
