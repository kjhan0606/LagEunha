/* Compare fixed Newton-correction counts on one shock/wall suite. */
#include <stdio.h>
#include "../gfs_kep_es_ale.h"
#include "../gfs_hllc_pressure_ale.h"

static void normal(int d,double n[3]){n[0]=1;n[1]=n[2]=0;if(d==2){n[0]=.6;n[1]=.8;}if(d==3){n[0]=1./3;n[1]=2./3;n[2]=2./3;}}
static double wall_curve(double p,double rho,double p0,double g){double A=2/((g+1)*rho),B=(g-1)/(g+1)*p0;return (p-p0)*sqrt(A/(p+B));}
static double wall_pressure(double rho,double p0,double u,double g){double lo=p0,hi=2*p0;while(wall_curve(hi,rho,p0,g)<u)hi*=2;for(int k=0;k<100;k++){double m=.5*(lo+hi);if(wall_curve(m,rho,p0,g)<u)lo=m;else hi=m;}return .5*(lo+hi);}
static int flux(const GfsFvFan *fan,double wn,double f[5]){GfsFvSample s=gfs_fv_sample(fan,wn);if(!s.valid)return 0;for(int c=0;c<5;c++)f[c]=s.flux[c];return 1;}
static void phys(int d,GfsFvPrimitive p,const double n[3],double g,double wn,double f[5]){double U[5],F[5];gfs_fv_state_flux(d,p,n,g,U,F);for(int c=0;c<5;c++)f[c]=F[c]-wn*U[c];}
static double ferr(const double a[5],const double b[5]){double e=0;for(int c=0;c<5;c++)e=fmax(e,fabs(a[c]-b[c])/(1+fabs(a[c])+fabs(b[c])));return e;}
int main(void){
    const double g=1.4,machs[3]={1.05,1.5,3},boosts[3]={0,.5,2};
    printf("{\"scope\":\"face-only correction-count comparison; no evolution\",\"correction_counts\":[");
    for(int steps=0;steps<=2;steps++){
        double es=0,ew=0,em=0,ee=0;int fallback=0;int first=steps==0;
        for(int d=1;d<=3;d++){
            double n[3];normal(d,n);
            for(int im=0;im<3;im++)for(int ib=0;ib<3;ib++){
                double M=machs[im],boost=boosts[ib],c=sqrt(g),rr=(g+1)*M*M/((g-1)*M*M+2),pr=1+2*g/(g+1)*(M*M-1);
                GfsFvPrimitive l={1,1,{0}},r={rr,pr,{0}};for(int a=0;a<d;a++){l.u[a]=(boost+M*c)*n[a];r.u[a]=(boost+M*c/rr)*n[a];}
                GfsFvFan fan=gfs_fv_fan_pressure_corrected_steps(d,l,r,n,g,steps);double f[5],ref[5];if(!fan.valid||!flux(&fan,boost,f))return 2;
                fallback+=fan.hll;phys(d,l,n,g,boost,ref);es=fmax(es,ferr(f,ref));
            }
            for(int im=0;im<3;im++){
                double M=machs[im],u=M*sqrt(g),p=wall_pressure(1,1,u,g);GfsFvPrimitive l={1,1,{0}},r={1,1,{0}};
                for(int a=0;a<d;a++){l.u[a]=u*n[a];r.u[a]=-u*n[a];}
                GfsFvFan fan=gfs_fv_fan_pressure_corrected_steps(d,l,r,n,g,steps);double f[5];if(!fan.valid||!flux(&fan,0,f))return 3;
                double pn=0;for(int a=0;a<d;a++)pn+=f[a+1]*n[a];ew=fmax(ew,fabs(pn-p)/p);em=fmax(em,fabs(f[0]));ee=fmax(ee,fabs(f[4]));
            }
        }
        if(!first)printf(",");
        printf("{\"pressure_newton_corrections\":%d,\"max_stationary_shock_flux_error\":%.17g,\"max_wall_pressure_relative_error\":%.17g,\"max_wall_mass_flux\":%.17g,\"max_wall_energy_flux\":%.17g,\"HLL_fallbacks\":%d}",steps,es,ew,em,ee,fallback);
    }
    printf("]}\n");return 0;
}
