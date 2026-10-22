#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "gfs_ale_fan.h"
static void require(bool b,const char *s){if(!b){fprintf(stderr,"FAIL %s\n",s);exit(1);}}
struct Vec {double x,y;};
static Vec add(Vec a,Vec b){return {a.x+b.x,a.y+b.y};}
static Vec sub(Vec a,Vec b){return {a.x-b.x,a.y-b.y};}
static Vec times(Vec a,double t){return {a.x*t,a.y*t};}
static double dot(Vec a,Vec b){return a.x*b.x+a.y*b.y;}
static Vec velocity(const double B[4],Vec x,Vec b){return {B[0]*x.x+B[1]*x.y+b.x,B[2]*x.x+B[3]*x.y+b.y};}
struct Result {double formula_error=0,minimum_RMS=0,expected_RMS=0,massflux=0,first_moment=0,expected_first=0,euler_flux_error=0,weight_shift=0;};
static Result check(const double B[4],double angle,double d,double length,Vec midpoint,Vec translation){
    Vec n={cos(angle),sin(angle)},t={-sin(angle),cos(angle)};
    Vec qi=sub(midpoint,times(n,.5*d)),qj=add(midpoint,times(n,.5*d));
    Vec vi=velocity(B,qi,translation),vj=velocity(B,qj,translation);
    Vec Dn={B[0]*n.x+.5*(B[1]+B[2])*n.y,.5*(B[1]+B[2])*n.x+B[3]*n.y};
    double slope=2*dot(t,Dn),rho=1.3,P=1.e-6,gamma=5./3;
    const double x[4]={-.8611363115940525752,-.3399810435848562648,.3399810435848562648,.8611363115940525752};
    const double w[4]={.3478548451374538574,.6521451548625461426,.6521451548625461426,.3478548451374538574};
    double mismatch[4],wn[4],s[4];Vec point[4],u[4];Result out;double mean=0;
    for(int k=0;k<4;k++){
        s[k]=.5*length*x[k];point[k]=add(midpoint,times(t,s[k]));u[k]=velocity(B,point[k],translation);
        wn[k]=dot(vj,n)-dot(sub(vj,vi),sub(point[k],qi))/d;
        mismatch[k]=dot(u[k],n)-wn[k];mean+=.5*w[k]*mismatch[k];
        out.formula_error=std::max(out.formula_error,fabs(mismatch[k]-slope*s[k]));
        /* Independently sample the production full Euler lab fan at the
         * actual geometric face speed. Matching states are not zero ALE flux
         * when the power face does not follow this affine fluid field. */
        GfsAlePrimitive p={rho,P,u[k].x,u[k].y};auto fan=gfs_ale_fan_build(p,p,n.x,n.y,gamma);
        auto sample=gfs_ale_fan_sample(&fan,wn[k]);require(sample.valid,"native matching-state ALE fan");
        double E=P/(gamma-1)+.5*rho*dot(u[k],u[k]),un=dot(u[k],n);
        double exact[4]={rho*mismatch[k],rho*u[k].x*mismatch[k]+P*n.x,
            rho*u[k].y*mismatch[k]+P*n.y,E*mismatch[k]+P*un};
        for(int c=0;c<4;c++)out.euler_flux_error=std::max(out.euler_flux_error,fabs(sample.flux[c]-exact[c])/(1+fabs(exact[c])));
        out.massflux+=.5*length*w[k]*sample.flux[0];out.first_moment+=.5*length*w[k]*s[k]*sample.flux[0];
    }
    out.weight_shift=2*d*mean;double rms=0;
    for(int k=0;k<4;k++)rms+=.5*w[k]*pow(mismatch[k]-mean,2);
    out.minimum_RMS=sqrt(rms);out.expected_RMS=fabs(slope)*length/sqrt(12.);
    out.expected_first=rho*slope*length*length*length/12.;
    require(out.formula_error<1.e-12 && out.euler_flux_error<1.e-12,"affine geometry and full Euler ALE identity");
    require(fabs(out.minimum_RMS-out.expected_RMS)<1.e-12,"best scalar weight-rate cannot alter tangential slope");
    require(fabs(out.massflux)<1.e-12 && fabs(out.first_moment-out.expected_first)<1.e-12,"zero net mass does not imply zero spatial transport");
    return out;
}
int main(){
    const char *names[5]={"rigid_rotation","isotropic_dilation","simple_shear","symmetric_strain","rotation_plus_strain"};
    const double B[5][4]={{0,-.7,.7,0},{.3,0,0,.3},{0,1.2,0,0},{0,1.2,1.2,0},{.1,-.4,.9,-.1}};
    printf("{\"scope\":\"native affine-face kinematic counterexample and matching-state full Euler flux; no integration in time or IC modification\",\"cases\":[");
    for(int i=0;i<5;i++){
        auto r=check(B[i],0,.7,.4,{.13,-.09},{.2,-.1});if(i)printf(",");
        printf("{\"name\":\"%s\",\"minimum_RMS_normal_speed_mismatch_after_scalar_weight_rate\":%.17g,\"expected_RMS\":%.17g,\"optimal_Wdot_difference\":%.17g,\"net_mass_flux\":%.17g,\"tangential_first_mass_flux_moment\":%.17g,\"expected_first_moment\":%.17g,\"Euler_flux_error\":%.17g}",
            names[i],r.minimum_RMS,r.expected_RMS,r.weight_shift,r.massflux,r.first_moment,r.expected_first,r.euler_flux_error);
    }
    double formula=0,flux=0,rms=0,mass=0,moment=0;unsigned long long seed=42;
    auto random=[&](){seed=seed*6364136223846793005ULL+1442695040888963407ULL;return (seed>>11)*0x1.0p-53;};
    for(int k=0;k<2000;k++){
        double A[4];for(double &v:A)v=2*random()-1;
        auto r=check(A,6.28*random(),.4+random(),.2+random(),{random()-.5,random()-.5},{random()-.5,random()-.5});
        formula=std::max(formula,r.formula_error);flux=std::max(flux,r.euler_flux_error);
        rms=std::max(rms,fabs(r.minimum_RMS-r.expected_RMS));mass=std::max(mass,fabs(r.massflux));moment=std::max(moment,fabs(r.first_moment-r.expected_first));
    }
    printf("],\"random_cases\":2000,\"max_formula_error\":%.17g,\"max_native_Euler_flux_error\":%.17g,\"max_RMS_formula_error\":%.17g,\"max_net_mass_flux\":%.17g,\"max_first_moment_error\":%.17g}\n",formula,flux,rms,mass,moment);return 0;
}
