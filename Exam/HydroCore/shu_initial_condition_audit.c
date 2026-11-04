/* Independent quadrature of the original sin(5*pi*x) Shu--Osher IC.
 * Includes the real native initializer; no evolution or state correction. */
#define main archived_hydro_suite_main
#include "../../Hydro1DExam/laguerre_sod.c"
#undef main

int main(void)
{
    const double xg[4]={.1834346424956498,.5255324099163290,
                       .7966664774136267,.9602898564975363};
    const double wg[4]={.3626837833783620,.3137066458778873,
                       .2223810344533745,.1012285362903763};
    double density_L1=0,max_density=0,max_pressure=0,max_velocity=0;
    setup_problem(PROB_SHUOSHER);unified1d=1;gfs_unified1d_init(200);
    for(int i=0;i<N;i++){
        double a=face[i],b=face[i+1],integral=0;
        for(int k=0;k<4;k++)for(int sign=-1;sign<=1;sign+=2){
            double x=.5*(a+b)+.5*(b-a)*sign*xg[k];
            integral+=.5*(b-a)*wg[k]*(x < -4 ? 3.857143 : 1+.2*sin(5*M_PI*x));
        }
        double exact=integral/(b-a),error=fabs(P[i].rho-exact);
        max_density=fmax(max_density,error);density_L1+=(b-a)*error;
        max_pressure=fmax(max_pressure,fabs(P[i].p-(P[i].x < -4 ? 10.33333 : 1)));
        max_velocity=fmax(max_velocity,fabs(P[i].v-(P[i].x < -4 ? 2.629369 : 0)));
    }
    printf("{\"scope\":\"native initial projection vs independent eight-point Gauss integration\","
           "\"cells\":%d,\"density_wavenumber\":\"5*pi\",\"density_L1\":%.17g,"
           "\"density_max_error\":%.17g,\"pressure_max_error\":%.17g,"
           "\"velocity_max_error\":%.17g,\"PDE_steps\":0,\"pass\":%s}\n",
           N,density_L1,max_density,max_pressure,max_velocity,
           max_density<1e-12&&max_pressure<1e-12&&max_velocity<1e-12?"true":"false");
    return max_density<1e-12&&max_pressure<1e-12&&max_velocity<1e-12?0:1;
}
