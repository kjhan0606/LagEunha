/* Broad frozen-state audit of the one-correction pressure HLLC candidate. */
#include <stdio.h>
#include <stdint.h>
#include "../gfs_kep_es_ale.h"
#include "../gfs_hllc_pressure_ale.h"
#include "../gfs_exact_riemann_ale.h"

enum {N=6000};
static uint64_t rng=UINT64_C(0xd1b54a32d192ed03);
static double uniform(void){rng^=rng>>12;rng^=rng<<25;rng^=rng>>27;return (rng*UINT64_C(2685821657736338717)>>11)*(1.0/9007199254740992.0);}
static void normal(int d,double n[3]){n[0]=1;n[1]=n[2]=0;if(d==2){n[0]=.6;n[1]=.8;}if(d==3){n[0]=1./3;n[1]=2./3;n[2]=2./3;}}
static double entropy_defect(int d,GfsFvPrimitive l,GfsFvPrimitive r,const double n[3],double g,double wn,const double f[5]){
    double vl[5],vr[5],q=0,ul=0,ur=0;gfs_kep_entropy_variables(d,l,g,vl);gfs_kep_entropy_variables(d,r,g,vr);
    for(int c=0;c<5;c++)q+=(vr[c]-vl[c])*f[c];
    for(int a=0;a<d;a++){ul+=l.u[a]*n[a];ur+=r.u[a]*n[a];}
    return q-(r.rho*(ur-wn)-l.rho*(ul-wn));
}
static double flux_error(const double a[5],const double b[5]){double e=0;for(int c=0;c<5;c++)e=fmax(e,fabs(a[c]-b[c])/(1+fabs(a[c])+fabs(b[c])));return e;}
int main(void){
    const double gamma=1.4;int invalid_pc=0,invalid_exact=0,exact_vacuum=0,entropy_positive=0;
    double max_flux_error=0,max_hllc_error=0,max_entropy_positive=0;long long exact_iterations=0;
    printf("{\"scope\":\"frozen random Riemann states only; rho,p in [1e-3,1e3], component velocities in [-10,10], not PDE qualification\",\"states_per_dim\":%d,\"dimensions\":[",N);
    for(int d=1;d<=3;d++){
        double n[3];normal(d,n);if(d>1)printf(",");
        double dim_max=0,dim_hllc=0,dim_ent=0,worst_wn=0,worst_star_p=0;
        GfsFvPrimitive worst_l={0},worst_r={0};int worst_hll=0,dim_invalid=0,dim_fallback=0,dim_positive=0;
        for(int i=0;i<N;i++){
            GfsFvPrimitive l={exp(log(1e-3)+uniform()*log(1e6)),exp(log(1e-3)+uniform()*log(1e6)),{0}};
            GfsFvPrimitive r={exp(log(1e-3)+uniform()*log(1e6)),exp(log(1e-3)+uniform()*log(1e6)),{0}};
            for(int a=0;a<d;a++){l.u[a]=20*(uniform()-.5);r.u[a]=20*(uniform()-.5);}
            double wn=10*(uniform()-.5),fp[5],fh[5],fe[5];
            GfsFvFan pc=gfs_fv_fan_pressure_corrected_steps(d,l,r,n,gamma,2),h=gfs_fv_fan(d,l,r,n,gamma);
            GfsExactStar star;int ok_exact=gfs_exact_riemann_ale_flux(d,l,r,n,gamma,wn,fe,&star);
            if(!ok_exact){invalid_exact++;continue;}if(star.vacuum){exact_vacuum++;continue;}
            exact_iterations+=star.iterations;
            if(!pc.valid){invalid_pc++;dim_invalid++;continue;}if(pc.hll)dim_fallback++;
            GfsFvSample sp=gfs_fv_sample(&pc,wn),sh=gfs_fv_sample(&h,wn);if(!sp.valid||!sh.valid){invalid_pc++;dim_invalid++;continue;}
            for(int c=0;c<5;c++){fp[c]=sp.flux[c];fh[c]=sh.flux[c];}
            double e=flux_error(fp,fe),eh=flux_error(fh,fe),ed=entropy_defect(d,l,r,n,gamma,wn,fp)/(1+fabs(fe[4]));
            if(e>dim_max){dim_max=e;worst_l=l;worst_r=r;worst_wn=wn;worst_star_p=star.p;worst_hll=pc.hll;}
            dim_hllc=fmax(dim_hllc,eh);dim_ent=fmax(dim_ent,ed);
            if(ed>1e-12)dim_positive++;
        }
        printf("{\"dim\":%d,\"max_pressure_corrected_flux_error_vs_exact\":%.17g,\"max_HLLC_flux_error_vs_exact\":%.17g,\"max_entropy_defect_positive_part\":%.17g,\"entropy_positive_count\":%d,\"invalid_candidate_count\":%d,\"candidate_HLL_fallback_count\":%d,\"worst_state\":{\"left_rho\":%.17g,\"left_p\":%.17g,\"left_u\":[%.17g,%.17g,%.17g],\"right_rho\":%.17g,\"right_p\":%.17g,\"right_u\":[%.17g,%.17g,%.17g],\"wn\":%.17g,\"exact_star_p\":%.17g,\"candidate_HLL_fallback\":%d}}",d,dim_max,dim_hllc,fmax(0,dim_ent),dim_positive,dim_invalid,dim_fallback,worst_l.rho,worst_l.p,worst_l.u[0],worst_l.u[1],worst_l.u[2],worst_r.rho,worst_r.p,worst_r.u[0],worst_r.u[1],worst_r.u[2],worst_wn,worst_star_p,worst_hll);
        max_flux_error=fmax(max_flux_error,dim_max);max_hllc_error=fmax(max_hllc_error,dim_hllc);max_entropy_positive=fmax(max_entropy_positive,fmax(0,dim_ent));
        entropy_positive+=dim_positive;
    }
    printf("],\"invalid_exact_count\":%d,\"exact_vacuum_count\":%d,\"mean_exact_pressure_iterations_nonvacuum\":%.17g,\"max_candidate_flux_error\":%.17g,\"max_HLLC_flux_error\":%.17g,\"aggregate_entropy_positive_count\":%d,\"max_entropy_positive_part\":%.17g}\n",invalid_exact,exact_vacuum,(double)exact_iterations/(3*N-exact_vacuum-invalid_exact),max_flux_error,max_hllc_error,entropy_positive,max_entropy_positive);
    return invalid_pc?1:0;
}
