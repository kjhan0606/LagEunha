#ifndef GFS_HLLC_PRESSURE_ALE_H
#define GFS_HLLC_PRESSURE_ALE_H
/* Experimental pressure-corrected HLLC fan. The star pressure is initialized
 * with closed-form two-shock/two-rarefaction estimates and receives at most
 * two safeguarded Newton corrections to the ideal-gas matching equation;
 * that one estimate supplies the outer HLLC wave speeds.
 * The conservative HLLC Rankine--Hugoniot states and ALE sampling are shared
 * with the existing fan representation. This is a single approximate
 * Riemann solver, not a test-specific stress or sensor. */
#include "gfs_exact_riemann_ale.h"

static inline int gfs_hllc_pressure_guess_steps(GfsFvPrimitive l,GfsFvPrimitive r,
        double ul,double ur,double cl,double cr,double gamma,int corrections,double *pstar){
    if(corrections<0||corrections>2)return 0;
    double gm1=gamma-1,a=gm1/(2*gamma);
    double vacuum_limit=2*(cl+cr)/gm1;
    if(ur-ul>=vacuum_limit)return 0;
    double pvrs=.5*(l.p+r.p)-.125*(ur-ul)*(l.rho+r.rho)*(cl+cr);
    double p=pvrs;
    if(ul>ur){
        double AL=2/((gamma+1)*l.rho),AR=2/((gamma+1)*r.rho);
        double BL=(gamma-1)/(gamma+1)*l.p,BR=(gamma-1)/(gamma+1)*r.p;
        double wl=sqrt(AL),wr=sqrt(AR),sum=wl+wr;
        double Bbar=(wl*BL+wr*BR)/sum;
        double T=wl*(BL+l.p)+wr*(BR+r.p);
        double y=(ul-ur+sqrt((ul-ur)*(ul-ur)+4*sum*T))/(2*sum);
        double two_shock=y*y-Bbar;
        if(two_shock>fmax(l.p,r.p)&&isfinite(two_shock))p=two_shock;
        else if(!(pvrs>fmin(l.p,r.p))){
            double numerator=cl+cr-.5*gm1*(ur-ul);
            double denominator=cl*pow(l.p,-a)+cr*pow(r.p,-a);
            if(!(numerator>0)||!(denominator>0))return 0;
            p=pow(numerator/denominator,1/a);
        }
    }else if(!(pvrs>fmin(l.p,r.p))){
        double numerator=cl+cr-.5*gm1*(ur-ul);
        double denominator=cl*pow(l.p,-a)+cr*pow(r.p,-a);
        if(!(numerator>0)||!(denominator>0))return 0;
        p=pow(numerator/denominator,1/a);
    }
    if(!(p>0)||!isfinite(p))return 0;
    for(int k=0;k<corrections;k++){
        double f=gfs_exact_wave(p,l,cl,gamma)+gfs_exact_wave(p,r,cr,gamma)+ur-ul;
        if(!isfinite(f))return 0;
        double scale=cl+cr+fabs(ul)+fabs(ur)+sqrt(gamma*l.p/l.rho)+sqrt(gamma*r.p/r.rho);
        if(fabs(f)<=2e-14*scale)break;
        double fp=gfs_exact_wave_derivative(p,l,cl,gamma)+gfs_exact_wave_derivative(p,r,cr,gamma);
        if(!(fp>0)||!isfinite(fp))return 0;
        double next=p-f/fp;
        if(!(next>0)||!isfinite(next))next=f<0?2*p:.5*p;
        if(!(next>0)||!isfinite(next))return 0;
        p=next;
    }
    if(!(p>0)||!isfinite(p))return 0;
    *pstar=p;return 1;
}

static inline GfsFvFan gfs_fv_fan_pressure_corrected_steps(int dim,GfsFvPrimitive left,
        GfsFvPrimitive right,const double normal[3],double gamma,int corrections){
    GfsFvFan fan={0};fan.dim=dim;
    if(!gfs_fv_input(dim,left)||!gfs_fv_input(dim,right)||!(gamma>1)||!isfinite(gamma))return fan;
    double norm=0,ul=0,ur=0;
    for(int a=0;a<dim;a++){
        if(!isfinite(normal[a]))return fan;
        norm+=normal[a]*normal[a];ul+=left.u[a]*normal[a];ur+=right.u[a]*normal[a];
    }
    if(fabs(norm-1)>1e-12)return fan;
    gfs_fv_state_flux(dim,left,normal,gamma,fan.state[0],fan.flux[0]);
    gfs_fv_state_flux(dim,right,normal,gamma,fan.state[3],fan.flux[3]);
    double cl=sqrt(gamma*left.p/left.rho),cr=sqrt(gamma*right.p/right.rho),pstar;
    if(!(cl>0)||!(cr>0)||!isfinite(cl)||!isfinite(cr))return fan;
    if(!gfs_hllc_pressure_guess_steps(left,right,ul,ur,cl,cr,gamma,corrections,&pstar))
        return gfs_fv_fan(dim,left,right,normal,gamma);
    double ql=1,qr=1;
    if(pstar>left.p)ql=sqrt(1+(gamma+1)/(2*gamma)*(pstar/left.p-1));
    if(pstar>right.p)qr=sqrt(1+(gamma+1)/(2*gamma)*(pstar/right.p-1));
    double sl=ul-cl*ql,sr=ur+cr*qr;
    if(!(sr>sl)||!isfinite(sl)||!isfinite(sr))return gfs_fv_fan(dim,left,right,normal,gamma);
    double al=left.rho*(sl-ul),ar=right.rho*(sr-ur);
    double denominator=al-ar;
    if(denominator==0||!isfinite(denominator))return gfs_fv_fan(dim,left,right,normal,gamma);
    double sm=(right.p-left.p+al*ul-ar*ur)/denominator;
    double ps=left.p+al*(sm-ul);
    int good=isfinite(sm)&&isfinite(ps)&&ps>0&&sm>=sl&&sm<=sr;
    if(good)for(int side=0;side<2;side++){
        int src=side?3:0,dst=side?2:1;
        const double *state=fan.state[src];double *qs=fan.state[dst];
        double sk=side?sr:sl,uk=side?ur:ul,pk=side?right.p:left.p;
        double factor=state[0]*(sk-uk)/(sk-sm);
        qs[0]=factor;
        for(int a=0;a<dim;a++)qs[a+1]=factor*(state[a+1]/state[0]+(sm-uk)*normal[a]);
        qs[4]=factor*(state[4]/state[0]+(sm-uk)*(sm+pk/(state[0]*(sk-uk))));
        double thermal=gfs_fv_internal(dim,qs);
        if(!(thermal>0)||!isfinite(thermal)||!isfinite(qs[4]))good=0;
        for(int c=0;c<5;c++)fan.flux[dst][c]=fan.flux[src][c]+sk*(qs[c]-state[c]);
    }
    if(!good)return gfs_fv_fan(dim,left,right,normal,gamma);
    fan.speed[0]=sl;fan.speed[1]=sm;fan.speed[2]=sr;fan.valid=1;return fan;
}

static inline GfsFvFan gfs_fv_fan_pressure_corrected(int dim,GfsFvPrimitive left,
        GfsFvPrimitive right,const double normal[3],double gamma){
    return gfs_fv_fan_pressure_corrected_steps(dim,left,right,normal,gamma,2);
}

#endif
