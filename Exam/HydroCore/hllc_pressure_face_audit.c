/* Frozen-face audit of a single pressure-corrected HLLC ALE fan. */
#include <stdio.h>
#include "../gfs_kep_es_ale.h"
#include "../gfs_hllc_pressure_ale.h"
#include "../gfs_exact_riemann_ale.h"

static void normal(int d,double n[3]){
    n[0]=1;n[1]=n[2]=0;
    if(d==2){n[0]=.6;n[1]=.8;}
    if(d==3){n[0]=1./3;n[1]=2./3;n[2]=2./3;}
}
static double dot(const double *a,const double *b,int d){double q=0;for(int k=0;k<d;k++)q+=a[k]*b[k];return q;}
static int sample(const GfsFvFan *fan,double wn,double f[5]){
    GfsFvSample s=gfs_fv_sample(fan,wn);if(!s.valid)return 0;
    for(int c=0;c<5;c++)f[c]=s.flux[c];
    return 1;
}
static void physical_ale(int d,GfsFvPrimitive q,const double n[3],double g,double wn,double f[5]){
    double U[5],F[5];gfs_fv_state_flux(d,q,n,g,U,F);for(int c=0;c<5;c++)f[c]=F[c]-wn*U[c];
}
static double relative_error(const double a[5],const double b[5]){
    double e=0;for(int c=0;c<5;c++)e=fmax(e,fabs(a[c]-b[c])/(1+fabs(a[c])+fabs(b[c])));return e;
}
static double entropy_defect(int d,GfsFvPrimitive l,GfsFvPrimitive r,const double n[3],
        double gamma,double wn,const double f[5]){
    double vl[5],vr[5],q=0,unl=0,unr=0;
    gfs_kep_entropy_variables(d,l,gamma,vl);gfs_kep_entropy_variables(d,r,gamma,vr);
    for(int c=0;c<5;c++)q+=(vr[c]-vl[c])*f[c];
    for(int a=0;a<d;a++){unl+=l.u[a]*n[a];unr+=r.u[a]*n[a];}
    return q-(r.rho*(unr-wn)-l.rho*(unl-wn));
}
static double wall_curve(double p,double rho,double p0,double gamma){
    double A=2/((gamma+1)*rho),B=(gamma-1)/(gamma+1)*p0;
    return (p-p0)*sqrt(A/(p+B));
}
static double wall_pressure(double rho,double p0,double u,double gamma){
    double lo=p0,hi=2*p0;while(wall_curve(hi,rho,p0,gamma)<u)hi*=2;
    for(int k=0;k<100;k++){double m=.5*(lo+hi);if(wall_curve(m,rho,p0,gamma)<u)lo=m;else hi=m;}
    return .5*(lo+hi);
}

int main(void){
    const double gamma=1.4,machs[3]={1.05,1.5,3.0},boosts[3]={0,.5,2};
    int bad=0;double max_shock_pc=0,max_shock_hllc=0,max_entropy_positive=0;
    double max_wall_pc=0,max_wall_hllc=0,max_wall_mass=0,max_wall_energy=0;
    double max_contact_error=0;
    printf("{\"scope\":\"frozen face audit only; no PDE evolution; pressure-corrected HLLC uses two safeguarded pressure Newton corrections\",\"dimensions\":[");
    for(int d=1;d<=3;d++){
        double n[3];normal(d,n);if(d>1)printf(",");printf("{\"dim\":%d,\"normal_shocks\":[",d);
        int first=1;
        for(int im=0;im<3;im++)for(int ib=0;ib<3;ib++){
            double M=machs[im],boost=boosts[ib],c=sqrt(gamma);
            double ratio=((gamma+1)*M*M)/((gamma-1)*M*M+2);
            double p_ratio=1+2*gamma/(gamma+1)*(M*M-1),un_l=boost+M*c,un_r=boost+M*c/ratio;
            GfsFvPrimitive l={1,1,{0}},r={ratio,p_ratio,{0}};
            for(int a=0;a<d;a++){l.u[a]=un_l*n[a];r.u[a]=un_r*n[a];}
            double fc[5],fh[5],fe[5],fphys[5],wn=boost;
            GfsFvFan pc=gfs_fv_fan_pressure_corrected(d,l,r,n,gamma),h=gfs_fv_fan(d,l,r,n,gamma);
            if(!sample(&pc,wn,fc)||!sample(&h,wn,fh)||!gfs_exact_riemann_ale_flux(d,l,r,n,gamma,wn,fe,NULL)){bad=1;continue;}
            physical_ale(d,l,n,gamma,wn,fphys);
            double ep=relative_error(fc,fphys),eh=relative_error(fh,fphys),ed=entropy_defect(d,l,r,n,gamma,wn,fc);
            max_shock_pc=fmax(max_shock_pc,ep);max_shock_hllc=fmax(max_shock_hllc,eh);
            max_entropy_positive=fmax(max_entropy_positive,ed/(1+fabs(fphys[4])));
            if(!first)printf(",");
            first=0;
            printf("{\"Mach\":%.17g,\"boost_over_sound\":%.17g,\"pressure_corrected_flux_error\":%.17g,\"HLLC_flux_error\":%.17g,\"pressure_corrected_entropy_defect\":%.17g}",M,boost,ep,eh,ed);
        }
        printf("],\"reflecting_walls\":[");first=1;
        for(int im=0;im<3;im++){
            double M=machs[im],u=M*sqrt(gamma),pstar=wall_pressure(1,1,u,gamma);
            GfsFvPrimitive l={1,1,{0}},r={1,1,{0}};
            for(int a=0;a<d;a++){l.u[a]=u*n[a];r.u[a]=-u*n[a];}
            GfsFvFan pc=gfs_fv_fan_pressure_corrected(d,l,r,n,gamma),h=gfs_fv_fan(d,l,r,n,gamma);
            double fc[5],fh[5];if(!sample(&pc,0,fc)||!sample(&h,0,fh)){bad=1;continue;}
            double pp=dot(fc+1,n,d),ph=dot(fh+1,n,d),ep=fabs(pp-pstar)/pstar,eh=fabs(ph-pstar)/pstar;
            max_wall_pc=fmax(max_wall_pc,ep);max_wall_hllc=fmax(max_wall_hllc,eh);
            max_wall_mass=fmax(max_wall_mass,fabs(fc[0]));max_wall_energy=fmax(max_wall_energy,fabs(fc[4]));
            if(fabs(fc[0])>1e-12||fabs(fc[4])>1e-12)bad=1;
            if(!first)printf(",");
            first=0;
            printf("{\"Mach\":%.17g,\"exact_pressure\":%.17g,\"pressure_corrected_pressure\":%.17g,\"pressure_corrected_relative_error\":%.17g,\"HLLC_pressure\":%.17g,\"HLLC_relative_error\":%.17g}",M,pstar,pp,ep,ph,eh);
        }
        printf("],\"contact_shear_control\":");
        GfsFvPrimitive l={4,1,{0}},r={1,1,{0}};GfsFvFan fan=gfs_fv_fan_pressure_corrected(d,l,r,n,gamma);double f[5];
        if(!sample(&fan,0,f)){bad=1;printf("false");continue;}
        double e=fmax(fabs(f[0]),fabs(f[4]));for(int a=0;a<d;a++)e=fmax(e,fabs(f[a+1]-n[a]));
        if(d>=2){
            double t[3]={-n[1],n[0],0};if(d==3){t[0]=n[1];t[1]=-n[0];t[2]=0;}
            for(int a=0;a<d;a++)r.u[a]=.7*t[a];
            fan=gfs_fv_fan_pressure_corrected(d,l,r,n,gamma);
            if(!sample(&fan,0,f)){bad=1;e=HUGE_VAL;}
            double es=fmax(fabs(f[0]),fabs(f[4]));for(int a=0;a<d;a++)es=fmax(es,fabs(f[a+1]-n[a]));e=fmax(e,es);
        }
        max_contact_error=fmax(max_contact_error,e);if(e>2e-12)bad=1;printf("%.17g",e);
        printf("}");
    }
    printf("],\"max_pressure_corrected_shock_flux_error\":%.17g,\"max_HLLC_shock_flux_error\":%.17g,\"max_pressure_corrected_entropy_positive_part\":%.17g,\"max_pressure_corrected_wall_pressure_relative_error\":%.17g,\"max_HLLC_wall_pressure_relative_error\":%.17g,\"max_pressure_corrected_wall_mass_flux\":%.17g,\"max_pressure_corrected_wall_energy_flux\":%.17g,\"max_pressure_corrected_contact_shear_error\":%.17g,\"algebraic_contacts_and_wall_fluxes_pass\":%s,\"no_evolution_acceptance_threshold_registered\":true,\"pde_evolution_verified\":false}\n",max_shock_pc,max_shock_hllc,max_entropy_positive,max_wall_pc,max_wall_hllc,max_wall_mass,max_wall_energy,max_contact_error,bad?"false":"true");
    return bad?1:0;
}
