/* Native algebraic audits of the new cell-pair core. No geometry/PDE. */
#include <stdio.h>
#include "../gfs_kep_ale.h"
static double error(double a,double b){return fabs(a-b)/(1+fabs(a)+fabs(b));}
static void entropy_variables(int d,GfsFvPrimitive p,double gamma,double v[5]){
    double v2=0,beta=.5*p.rho/p.p;for(int c=0;c<5;c++)v[c]=0;
    for(int a=0;a<d;a++){v2+=p.u[a]*p.u[a];v[a+1]=2*beta*p.u[a];}
    v[0]=(gamma-log(p.p)+gamma*log(p.rho))/(gamma-1)-beta*v2;v[4]=-2*beta;
}
int main(void){
    int bad=0,cases=0;double reversal=0,boost=0,entropy=0,kinetic=0,derivative=0,uniform=0,contact=0,stress=0,logarithmic=0;
    const double gamma=1.4;
    for(int d=1;d<=3;d++)for(int test=0;test<80;test++){
        GfsFvPrimitive l={1,1,{0}},r={.125,.1,{0}};double n[3]={0},neg[3]={0},shift[3]={0},norm=0,bn=0;
        l.rho=exp(2*sin(.71*test));r.rho=exp(2*cos(.49*test));l.p=exp(2*cos(.39*test));r.p=exp(2*sin(.27*test));
        if(test==0){l.rho=r.rho=1;l.p=r.p=1e-6;}
        if(test==1){l.rho=1;r.rho=1+1e-14;l.p=1;r.p=1+2e-14;}
        for(int a=0;a<d;a++){
            n[a]=a+1;norm+=n[a]*n[a];l.u[a]=1.3*sin(.7*(test+1)*(a+1));r.u[a]=1.1*cos(.9*(test+1)*(a+1));shift[a]=.8*(a+1);
        }
        for(int a=0;a<d;a++){n[a]/=sqrt(norm);neg[a]=-n[a];bn+=shift[a]*n[a];}
        double wn=.6*sin(.2*test),f[5],rev[5],fb[5],rho,unused;GfsFvPrimitive lb=l,rb=r;
        for(int a=0;a<d;a++){lb.u[a]+=shift[a];rb.u[a]+=shift[a];}
        if(!gfs_kep_flux(d,l,r,n,gamma,wn,0,f,&rho)||!gfs_kep_flux(d,r,l,neg,gamma,-wn,0,rev,&unused)||!gfs_kep_flux(d,lb,rb,n,gamma,wn+bn,0,fb,&unused)){bad=1;continue;}
        cases++;
        double energy=f[4];for(int a=0;a<d;a++)energy+=shift[a]*f[a+1]+.5*shift[a]*shift[a]*f[0];
        boost=fmax(boost,error(energy,fb[4]));boost=fmax(boost,error(f[0],fb[0]));
        for(int a=0;a<d;a++)boost=fmax(boost,error(f[a+1]+shift[a]*f[0],fb[a+1]));
        for(int c=0;c<5;c++)reversal=fmax(reversal,error(f[c],-rev[c]));
        double vl[5],vr[5],jump=0,scale=1,ln=0,rn=0,du=0,v2diff=0,kinrate=0;
        entropy_variables(d,l,gamma,vl);entropy_variables(d,r,gamma,vr);
        for(int c=0;c<5;c++){jump+=(vr[c]-vl[c])*f[c];scale+=fabs((vr[c]-vl[c])*f[c]);}
        for(int a=0;a<d;a++){ln+=l.u[a]*n[a];rn+=r.u[a]*n[a];du+=(r.u[a]-l.u[a])*n[a];v2diff+=r.u[a]*r.u[a]-l.u[a]*l.u[a];kinrate+=(r.u[a]-l.u[a])*f[a+1];}
        double potential=r.rho*(rn-wn)-l.rho*(ln-wn);
        entropy=fmax(entropy,fabs(jump-potential)/(scale+fabs(potential)));
        double pressure=(.5*l.rho+.5*r.rho)/(l.rho/(2*l.p)+r.rho/(2*r.p));
        kinetic=fmax(kinetic,error(kinrate-.5*v2diff*f[0],pressure*du));
        double plus[5],minus[5],h=1e-5;gfs_kep_flux(d,l,r,n,gamma,wn+h,0,plus,&unused);gfs_kep_flux(d,l,r,n,gamma,wn-h,0,minus,&unused);
        derivative=fmax(derivative,error((plus[0]-minus[0])/(2*h),-rho));
        double pi=du<0?.3:0,sf[5];gfs_kep_flux(d,l,r,n,gamma,wn,pi,sf,&unused);
        double production=0;for(int c=0;c<5;c++)production+=(vr[c]-vl[c])*(sf[c]-f[c]);
        double expected=pi*(l.rho/l.p+r.rho/r.p)*.5*du;
        /* 2*beta_bar = (rho_L/P_L+rho_R/P_R)/2. */
        stress=fmax(stress,error(production,expected));if(production>1e-9*(1+fabs(expected)))bad=1;
        GfsFvPrimitive same=l;double uf[5],state[5],physical[5];gfs_kep_flux(d,l,same,n,gamma,wn,0,uf,&unused);gfs_fv_state_flux(d,l,n,gamma,state,physical);
        for(int c=0;c<5;c++)uniform=fmax(uniform,error(uf[c],physical[c]-wn*state[c]));
        same.rho=r.rho;gfs_kep_flux(d,l,same,n,gamma,ln,0,uf,&unused);contact=fmax(contact,fabs(uf[0]));
        for(int a=0;a<d;a++)contact=fmax(contact,error(uf[a+1],l.p*n[a]));contact=fmax(contact,error(uf[4],l.p*ln));
        long double a=l.rho,b=r.rho,exact=a==b?a:(b-a)/(logl(b)-logl(a));logarithmic=fmax(logarithmic,error(rho,(double)exact));
        for(int a=d;a<3;a++)if(f[a+1]!=0)bad=1;
    }
    bad|=cases!=240||reversal>1e-12||boost>1e-12||entropy>1e-11||kinetic>1e-12||derivative>1e-8||uniform>1e-12||contact>1e-12||stress>1e-9||logarithmic>1e-12;
    printf("{\"scope\":\"native D1/D2/D3 cell-pair KEP/EC ALE algebra only; no full timestep or PDE\",\"cases\":%d,\"reversal\":%.17g,\"Galilean\":%.17g,\"ALE_entropy\":%.17g,\"kinetic_chain\":%.17g,\"mass_derivative\":%.17g,\"uniform\":%.17g,\"translating_contact\":%.17g,\"compression_entropy\":%.17g,\"logmean_long_double\":%.17g,\"pass\":%s}\n",cases,reversal,boost,entropy,kinetic,derivative,uniform,contact,stress,logarithmic,bad?"false":"true");return bad?1:0;
}
