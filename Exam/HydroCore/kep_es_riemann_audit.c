/* Cheap 1D/2D/3D exact normal-shock and reflecting-wall Riemann probes. */
#include <stdio.h>
#include "../gfs_kep_es_ale.h"
#include "../gfs_exact_riemann_ale.h"

static double dot(const double *a,const double *b,int d){double x=0;for(int k=0;k<d;k++)x+=a[k]*b[k];return x;}
static void normal(int d,double n[3]){
    n[0]=1;n[1]=n[2]=0;
    if(d==2){n[0]=.6;n[1]=.8;}
    if(d==3){n[0]=1./3;n[1]=2./3;n[2]=2./3;}
}
static void entropy(int d,GfsFvPrimitive p,double g,double v[5]){
    gfs_kep_entropy_variables(d,p,g,v);
}
static double entropy_defect(int d,GfsFvPrimitive l,GfsFvPrimitive r,
        const double n[3],double g,double wn,const double f[5]){
    double vl[5],vr[5],du=0,unl=0,unr=0;
    entropy(d,l,g,vl);entropy(d,r,g,vr);
    for(int c=0;c<5;c++)du+=(vr[c]-vl[c])*f[c];
    for(int a=0;a<d;a++){unl+=l.u[a]*n[a];unr+=r.u[a]*n[a];}
    return du-(r.rho*(unr-wn)-l.rho*(unl-wn));
}
static double wall_curve(double p,double rho,double p0,double gamma){
    double A=2/((gamma+1)*rho),B=(gamma-1)/(gamma+1)*p0;
    return (p-p0)*sqrt(A/(p+B));
}
static double wall_star_pressure(double rho,double p0,double speed,double gamma){
    double lo=p0,hi=2*p0;
    while(wall_curve(hi,rho,p0,gamma)<speed)hi*=2;
    for(int k=0;k<100;k++){
        double mid=.5*(lo+hi);
        if(wall_curve(mid,rho,p0,gamma)<speed)lo=mid;else hi=mid;
    }
    return .5*(lo+hi);
}
static void physical_ale_flux(int d,GfsFvPrimitive p,const double n[3],double g,double wn,double f[5]){
    double U[5],F[5];gfs_fv_state_flux(d,p,n,g,U,F);for(int c=0;c<5;c++)f[c]=F[c]-wn*U[c];
}
static int hllc_ale_flux(int d,GfsFvPrimitive l,GfsFvPrimitive r,const double n[3],
        double g,double wn,double f[5]){
    GfsFvFan fan=gfs_fv_fan(d,l,r,n,g);
    GfsFvSample sample=gfs_fv_sample(&fan,wn);
    if(!sample.valid)return 0;
    for(int c=0;c<5;c++)f[c]=sample.flux[c];
    return 1;
}
int main(void){
    const double gamma=1.4;int bad=0;double max_rh=0,max_shock_flux_error=0,max_entropy_positive=0;
    double max_hllc_shock_flux_error=0,max_hllc_entropy_positive=0;
    double max_exact_shock_flux_error=0,max_exact_entropy_positive=0;int max_exact_iterations=0;
    double max_wall_mass=0,max_wall_energy=0,max_wall_pressure_rel=0,max_hllc_wall_pressure_rel=0;
    double max_hllc_wall_mass=0,max_hllc_wall_energy=0;
    double max_exact_wall_pressure_rel=0,max_exact_star_pressure_error=0;
    double max_exact_contact_shear_error=0;
    const double machs[3]={1.05,1.5,3.0};
    const double boosts[3]={0,.5,2.0};
    printf("{\"scope\":\"frozen face probes comparing entropy-variable, existing HLLC/HLL, and exact Godunov ALE flux; not PDE shock evolution\",\"gamma\":%.17g,\"normal_shocks\":[",gamma);
    int first=1;
    for(int d=1;d<=3;d++){
        double n[3];normal(d,n);
        for(int k=0;k<3;k++){
          for(int b=0;b<3;b++){
            double M=machs[k],c=sqrt(gamma),wn=boosts[b]*c,u1=M*c;
            double ratio=(gamma+1)*M*M/((gamma-1)*M*M+2);
            double p2=(2*gamma*M*M-(gamma-1))/(gamma+1),u2=u1/ratio;
            double M2=u2/sqrt(gamma*p2/ratio);
            GfsFvPrimitive l={1,1,{0}},r={ratio,p2,{0}};
            for(int a=0;a<d;a++){l.u[a]=(u1+wn)*n[a];r.u[a]=(u2+wn)*n[a];}
            double fn[5],fh[5],fe[5],fl[5],fr[5],unused=0;GfsExactStar se;
            if(!gfs_kep_es_flux(d,l,r,n,gamma,wn,fn,&unused,NULL)){bad=1;continue;}
            if(!hllc_ale_flux(d,l,r,n,gamma,wn,fh)){bad=1;continue;}
            if(!gfs_exact_riemann_ale_flux(d,l,r,n,gamma,wn,fe,&se)){bad=1;continue;}
            if(se.iterations>max_exact_iterations)max_exact_iterations=se.iterations;
            physical_ale_flux(d,l,n,gamma,wn,fl);physical_ale_flux(d,r,n,gamma,wn,fr);
            double scale=1,err=0,shock_flux_error=0,hllc_flux_error=0,exact_flux_error=0;for(int cidx=0;cidx<5;cidx++){
                scale+=fabs(fl[cidx])+fabs(fr[cidx]);err=fmax(err,fabs(fl[cidx]-fr[cidx]));
                shock_flux_error=fmax(shock_flux_error,fabs(fn[cidx]-fl[cidx])/(1+fabs(fn[cidx])+fabs(fl[cidx])));
                hllc_flux_error=fmax(hllc_flux_error,fabs(fh[cidx]-fl[cidx])/(1+fabs(fh[cidx])+fabs(fl[cidx])));
                exact_flux_error=fmax(exact_flux_error,fabs(fe[cidx]-fl[cidx])/(1+fabs(fe[cidx])+fabs(fl[cidx])));
            }
            double rh=err/scale;max_rh=fmax(max_rh,rh);
            max_shock_flux_error=fmax(max_shock_flux_error,shock_flux_error);
            max_hllc_shock_flux_error=fmax(max_hllc_shock_flux_error,hllc_flux_error);
            max_exact_shock_flux_error=fmax(max_exact_shock_flux_error,exact_flux_error);
            double ed=entropy_defect(d,l,r,n,gamma,wn,fn);
            double edh=entropy_defect(d,l,r,n,gamma,wn,fh);
            double ede=entropy_defect(d,l,r,n,gamma,wn,fe);
            max_entropy_positive=fmax(max_entropy_positive,ed/(1+fabs(fl[4])+fabs(fr[4])));
            max_hllc_entropy_positive=fmax(max_hllc_entropy_positive,edh/(1+fabs(fl[4])+fabs(fr[4])));
            max_exact_entropy_positive=fmax(max_exact_entropy_positive,ede/(1+fabs(fl[4])+fabs(fr[4])));
            if(rh>2e-13||ed>2e-11*(1+fabs(fl[4])+fabs(fr[4])))bad=1;
            if(!first)printf(",");
            first=0;
            printf("{\"dim\":%d,\"upstream_Mach_in_shock_frame\":%.17g,\"galilean_boost_over_sound\":%.17g,\"downstream_Mach_in_shock_frame\":%.17g,\"rho_ratio\":%.17g,\"p_ratio\":%.17g,\"RH_flux_relative_defect\":%.17g,\"candidate_ALE_flux_difference_vs_exact_RH\":%.17g,\"HLLC_ALE_flux_difference_vs_exact_RH\":%.17g,\"exact_Godunov_flux_difference_vs_exact_RH\":%.17g,\"candidate_entropy_defect\":%.17g,\"HLLC_entropy_defect\":%.17g,\"exact_Godunov_entropy_defect\":%.17g,\"exact_pressure_iterations\":%d}",d,M,boosts[b],M2,ratio,p2,rh,shock_flux_error,hllc_flux_error,exact_flux_error,ed,edh,ede,se.iterations);
          }
        }
    }
    printf("],\"reflecting_walls\":[");first=1;
    const double wall_machs[3]={.5,1.5,3.0};
    for(int d=1;d<=3;d++){
        double n[3];normal(d,n);
        for(int k=0;k<3;k++){
            double M=wall_machs[k],c=sqrt(gamma),u=M*c,pstar=wall_star_pressure(1,1,u,gamma);
            GfsFvPrimitive l={1,1,{0}},r={1,1,{0}};
            for(int a=0;a<d;a++){l.u[a]=u*n[a];r.u[a]=-u*n[a];}
            double f[5],unused=0;
            if(!gfs_kep_es_flux(d,l,r,n,gamma,0,f,&unused,NULL)){bad=1;continue;}
            double pnum=dot(f+1,n,d),mass=fabs(f[0]),energy=fabs(f[4]);
            max_wall_mass=fmax(max_wall_mass,mass);max_wall_energy=fmax(max_wall_energy,energy);
            double pressure_rel=fabs(pnum-pstar)/pstar;
            max_wall_pressure_rel=fmax(max_wall_pressure_rel,pressure_rel);
            double fh[5],fe[5];GfsExactStar se;
            if(!hllc_ale_flux(d,l,r,n,gamma,0,fh)||!gfs_exact_riemann_ale_flux(d,l,r,n,gamma,0,fe,&se)){bad=1;continue;}
            double ph=dot(fh+1,n,d),pressure_hllc_rel=fabs(ph-pstar)/pstar;
            max_hllc_wall_pressure_rel=fmax(max_hllc_wall_pressure_rel,pressure_hllc_rel);
            max_hllc_wall_mass=fmax(max_hllc_wall_mass,fabs(fh[0]));
            max_hllc_wall_energy=fmax(max_hllc_wall_energy,fabs(fh[4]));
            double pe=dot(fe+1,n,d),pressure_exact_rel=fabs(pe-pstar)/pstar;
            max_exact_wall_pressure_rel=fmax(max_exact_wall_pressure_rel,pressure_exact_rel);
            max_exact_star_pressure_error=fmax(max_exact_star_pressure_error,fabs(se.p-pstar)/pstar);
            if(se.iterations>max_exact_iterations)max_exact_iterations=se.iterations;
            if(mass>2e-12||energy>2e-12||fabs(fh[0])>2e-12||fabs(fh[4])>2e-12||fabs(fe[0])>2e-12||fabs(fe[4])>2e-12)bad=1;
            if(!first)printf(",");
            first=0;
            printf("{\"dim\":%d,\"incident_Mach\":%.17g,\"exact_star_pressure\":%.17g,\"candidate_wall_pressure\":%.17g,\"candidate_pressure_relative_error\":%.17g,\"HLLC_wall_pressure\":%.17g,\"HLLC_pressure_relative_error\":%.17g,\"exact_Godunov_wall_pressure\":%.17g,\"exact_Godunov_pressure_relative_error\":%.17g,\"candidate_mass_flux\":%.17g,\"candidate_energy_flux\":%.17g,\"HLLC_mass_flux\":%.17g,\"HLLC_energy_flux\":%.17g,\"exact_Godunov_mass_flux\":%.17g,\"exact_Godunov_energy_flux\":%.17g,\"exact_pressure_iterations\":%d}",d,M,pstar,pnum,pressure_rel,ph,pressure_hllc_rel,pe,pressure_exact_rel,f[0],f[4],fh[0],fh[4],fe[0],fe[4],se.iterations);
        }
    }
    for(int d=1;d<=3;d++){
        double n[3];normal(d,n);double f[5]={0};GfsExactStar s;
        GfsFvPrimitive l={4,1,{0}},r={1,1,{0}};
        if(!gfs_exact_riemann_ale_flux(d,l,r,n,gamma,0,f,&s))bad=1;
        max_exact_contact_shear_error=fmax(max_exact_contact_shear_error,fabs(f[0]));
        for(int a=0;a<d;a++)max_exact_contact_shear_error=fmax(max_exact_contact_shear_error,fabs(f[a+1]-n[a]));
        max_exact_contact_shear_error=fmax(max_exact_contact_shear_error,fabs(f[4]));

        for(int a=0;a<d;a++)l.u[a]=r.u[a]=.23*(a+1);
        double wn=dot(l.u,n,d);
        if(!gfs_exact_riemann_ale_flux(d,l,r,n,gamma,wn,f,&s))bad=1;
        max_exact_contact_shear_error=fmax(max_exact_contact_shear_error,fabs(f[0]));
        for(int a=0;a<d;a++)max_exact_contact_shear_error=fmax(max_exact_contact_shear_error,fabs(f[a+1]-n[a]));
        max_exact_contact_shear_error=fmax(max_exact_contact_shear_error,fabs(f[4]-wn));

        if(d>=2){
            double t[3]={n[1],-n[0],0};double tn=sqrt(dot(t,t,d));
            if(d==3&&tn<1e-12){t[0]=1;t[1]=0;t[2]=0;tn=1;}
            for(int a=0;a<d;a++)t[a]/=tn;
            l=(GfsFvPrimitive){1,1,{0}};r=(GfsFvPrimitive){1,1,{0}};
            for(int a=0;a<d;a++){l.u[a]=.4*n[a]-.4*t[a];r.u[a]=.4*n[a]+.4*t[a];}
            if(!gfs_exact_riemann_ale_flux(d,l,r,n,gamma,.4,f,&s))bad=1;
            max_exact_contact_shear_error=fmax(max_exact_contact_shear_error,fabs(f[0]));
            for(int a=0;a<d;a++)max_exact_contact_shear_error=fmax(max_exact_contact_shear_error,fabs(f[a+1]-n[a]));
            max_exact_contact_shear_error=fmax(max_exact_contact_shear_error,fabs(f[4]-.4));
        }
    }
    if(max_exact_contact_shear_error>2e-12)bad=1;
    printf("],\"max_RH_input_defect\":%.17g,\"max_candidate_exact_ALE_shock_flux_difference\":%.17g,\"max_HLLC_exact_ALE_shock_flux_difference\":%.17g,\"max_exact_Godunov_ALE_shock_flux_difference\":%.17g,\"max_candidate_entropy_positive_part\":%.17g,\"max_HLLC_entropy_positive_part\":%.17g,\"max_exact_Godunov_entropy_positive_part\":%.17g,\"max_candidate_wall_pressure_relative_error\":%.17g,\"max_HLLC_wall_pressure_relative_error\":%.17g,\"max_exact_Godunov_wall_pressure_relative_error\":%.17g,\"max_exact_star_pressure_relative_error\":%.17g,\"max_exact_contact_shear_flux_error\":%.17g,\"max_exact_pressure_iterations\":%d,\"max_candidate_wall_mass_flux\":%.17g,\"max_candidate_wall_energy_flux\":%.17g,\"max_HLLC_wall_mass_flux\":%.17g,\"max_HLLC_wall_energy_flux\":%.17g,\"face_algebra_invariants_pass\":%s,\"shock_and_wall_flux_differences_are_diagnostic_only\":true,\"no_physical_acceptance_threshold_registered\":true,\"pde_evolution_verified\":false}\n",max_rh,max_shock_flux_error,max_hllc_shock_flux_error,max_exact_shock_flux_error,max_entropy_positive,max_hllc_entropy_positive,max_exact_entropy_positive,max_wall_pressure_rel,max_hllc_wall_pressure_rel,max_exact_wall_pressure_rel,max_exact_star_pressure_error,max_exact_contact_shear_error,max_exact_iterations,max_wall_mass,max_wall_energy,max_hllc_wall_mass,max_hllc_wall_energy,bad?"false":"true");
    return bad?1:0;
}
