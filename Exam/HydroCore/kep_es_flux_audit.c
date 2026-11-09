/* Frozen algebra audit for the complete KEP/entropy-stable ALE face flux. */
#include <stdio.h>
#include <float.h>
#include "../gfs_kep_es_ale.h"

static double err(double a,double b){return fabs(a-b)/(1+fabs(a)+fabs(b));}
static double dot(const double *a,const double *b,int d){double z=0;for(int k=0;k<d;k++)z+=a[k]*b[k];return z;}
static void unit(double *a,int d){double z=sqrt(dot(a,a,d));for(int k=0;k<d;k++)a[k]/=z;}

int main(void){
    const double gamma=1.4;int bad=0,cases=0;
    double reversal=0,boost=0,entropy_positive=0,kinetic_positive=0;
    double mass_derivative=0,uniform=0,contact=0,shear=0,wall=0;
    for(int d=1;d<=3;d++)for(int t=0;t<64;t++){
        GfsFvPrimitive l={1,1,{0}},r={1,1,{0}},lb,rb;
        double n[3]={0},neg[3]={0},b[3]={0},bn=0;
        for(int a=0;a<d;a++){
            n[a]=sin(.71*(t+1)*(a+1))+.37*cos(.29*(t+2)*(a+1));
            l.u[a]=1.4*sin(.63*(t+1)*(a+1));r.u[a]=1.1*cos(.47*(t+2)*(a+1));
            b[a]=.31*(a+1)-.17*t;
        }
        unit(n,d);for(int a=0;a<d;a++){neg[a]=-n[a];bn+=b[a]*n[a];}
        l.rho=exp(1.8*sin(.31*t)+.2);r.rho=exp(1.7*cos(.23*t)-.1);
        l.p=exp(1.6*cos(.37*t)-.3);r.p=exp(1.5*sin(.19*t)-.2);
        if(t==0){l.rho=1;r.rho=.125;l.p=r.p=1e-6;}
        double ul=dot(l.u,n,d),ur=dot(r.u,n,d),wn=.6*sin(.2*t);
        double f[5],fr[5],fb[5],unused=0;
        if(!gfs_kep_es_flux(d,l,r,n,gamma,wn,f,&unused,NULL)||
           !gfs_kep_es_flux(d,r,l,neg,gamma,-wn,fr,&unused,NULL)){bad=1;continue;}
        lb=l;rb=r;for(int a=0;a<d;a++){lb.u[a]+=b[a];rb.u[a]+=b[a];}
        if(!gfs_kep_es_flux(d,lb,rb,n,gamma,wn+bn,fb,&unused,NULL)){bad=1;continue;}
        cases++;
        for(int c=0;c<5;c++)reversal=fmax(reversal,err(f[c],-fr[c]));
        boost=fmax(boost,err(f[0],fb[0]));
        for(int a=0;a<d;a++)boost=fmax(boost,err(f[a+1]+b[a]*f[0],fb[a+1]));
        double eb=f[4];for(int a=0;a<d;a++)eb+=b[a]*f[a+1]+.5*b[a]*b[a]*f[0];
        boost=fmax(boost,err(eb,fb[4]));

        double vl[5],vr[5],dvf=0,scale=1;
        gfs_kep_entropy_variables(d,l,gamma,vl);gfs_kep_entropy_variables(d,r,gamma,vr);
        for(int c=0;c<5;c++){dvf+=(vr[c]-vl[c])*f[c];scale+=fabs((vr[c]-vl[c])*f[c]);}
        double potential=r.rho*(ur-wn)-l.rho*(ul-wn);
        double entropy_defect=dvf-potential;
        entropy_positive=fmax(entropy_positive,entropy_defect/(scale+fabs(potential)));
        if(entropy_defect>2e-11*(1+scale+fabs(potential)))bad=1;

        double du=ur-ul,dk=0,kin=0;
        for(int a=0;a<d;a++){double dvel=r.u[a]-l.u[a];dk+=r.u[a]*r.u[a]-l.u[a]*l.u[a];kin+=dvel*f[a+1];}
        kin-=.5*dk*f[0];
        double phat=(.5*l.rho+.5*r.rho)/(l.rho/(2*l.p)+r.rho/(2*r.p));
        double kd=kin-phat*du;
        kinetic_positive=fmax(kinetic_positive,kd/(1+fabs(kin)+fabs(phat*du)));
        if(kd>2e-11*(1+fabs(kin)+fabs(phat*du)))bad=1;

        double fp[5],fm[5],h=1e-3*(1+fabs(wn));
        if(!gfs_kep_es_flux(d,l,r,n,gamma,wn+h,fp,&unused,NULL)||
           !gfs_kep_es_flux(d,l,r,n,gamma,wn-h,fm,&unused,NULL)){bad=1;continue;}
        double rholog=gfs_kep_logmean(l.rho,r.rho);
        mass_derivative=fmax(mass_derivative,err((fp[0]-fm[0])/(2*h),-rholog));
        GfsFvPrimitive same=l;double fs[5],state[5],physical[5];
        if(!gfs_kep_es_flux(d,l,same,n,gamma,wn,fs,&unused,NULL))bad=1;
        gfs_fv_state_flux(d,l,n,gamma,state,physical);
        for(int c=0;c<5;c++)uniform=fmax(uniform,err(fs[c],physical[c]-wn*state[c]));
    }

    for(int d=1;d<=3;d++){
        double n[3]={.31,-.72,.44};unit(n,d);
        GfsFvPrimitive l={4,1,{0}},r={1,1,{0}};double f[5],rho;
        if(!gfs_kep_es_flux(d,l,r,n,gamma,0,f,&rho,NULL))bad=1;
        for(int a=0;a<d;a++)contact=fmax(contact,err(f[a+1],n[a]));
        contact=fmax(contact,fabs(f[0]));contact=fmax(contact,fabs(f[4]));
        for(int a=0;a<d;a++)l.u[a]=r.u[a]=.23*(a+1);
        double wn=dot(l.u,n,d);
        if(!gfs_kep_es_flux(d,l,r,n,gamma,wn,f,&rho,NULL))bad=1;
        contact=fmax(contact,fabs(f[0]));
        for(int a=0;a<d;a++)contact=fmax(contact,err(f[a+1],n[a]));
        contact=fmax(contact,err(f[4],wn));
    }

    for(int d=2;d<=3;d++){
        double n[3]={.31,-.72,.44};unit(n,d);
        double t[3]={n[1],-n[0],0};if(d==3&&dot(t,t,d)<1e-12){t[0]=1;t[1]=0;t[2]=0;}unit(t,d);
        GfsFvPrimitive l={1,1,{0}},r={1,1,{0}};double f[5],rho;
        for(int a=0;a<d;a++){l.u[a]=.4*n[a]-.4*t[a];r.u[a]=.4*n[a]+.4*t[a];}
        if(!gfs_kep_es_flux(d,l,r,n,gamma,.4,f,&rho,NULL))bad=1;
        shear=fmax(shear,fabs(f[0]));for(int a=0;a<d;a++)shear=fmax(shear,err(f[a+1],n[a]));
        shear=fmax(shear,fabs(f[4]-.4));
        for(int a=0;a<d;a++){l.u[a]=.7*n[a];r.u[a]=-.7*n[a];}
        if(!gfs_kep_es_flux(d,l,r,n,gamma,0,f,&rho,NULL))bad=1;
        wall=fmax(wall,fabs(f[0]));wall=fmax(wall,fabs(f[4]));
    }
    bad|=cases!=192||reversal>2e-12||boost>2e-11||entropy_positive>2e-11||
         kinetic_positive>2e-10||mass_derivative>2e-8||uniform>2e-12||
         contact>2e-12||shear>2e-12||wall>2e-12;
    printf("{\"scope\":\"D1-D3 KEP/EC ALE plus rank-two entropy-stable acoustic dissipation; frozen face algebra only\",\"pairs\":%d,\"reversal\":%.17g,\"Galilean\":%.17g,\"entropy_positive_part\":%.17g,\"kinetic_positive_part\":%.17g,\"mass_derivative\":%.17g,\"uniform\":%.17g,\"contact\":%.17g,\"tangential_shear\":%.17g,\"wall_mass_energy\":%.17g,\"pass\":%s}\n",cases,reversal,boost,entropy_positive,kinetic_positive,mass_derivative,uniform,contact,shear,wall,bad?"false":"true");
    return bad?1:0;
}
