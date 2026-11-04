/* Dimensional native FV fan/cone audit; no geometry or fluid evolution. */
#include <stdio.h>
#include <math.h>
#include "../gfs_fv_ale.h"
#ifdef __cplusplus
#include "../gfs_ale_fan.h"
#endif

static double discrepancy(double a,double b){return fabs(a-b)/(1+fabs(a)+fabs(b));}
static void boost_state(int dim,double q[5],const double b[3]){
    double E=q[4];
    for(int a=0;a<dim;a++){E+=b[a]*q[a+1]+.5*b[a]*b[a]*q[0];q[a+1]+=b[a]*q[0];}
    q[4]=E;
}
int main(void)
{
    int bad=0,fallback=0,cases=0;
    double reversal=0,galilean=0,uniform=0,derivative=0,embedding=0,cone=0,old2d=0;
    const double gamma=1.4;
    for(int dim=1;dim<=3;dim++)for(int test=0;test<40;test++){
        double n[3]={0},neg[3]={0},b[3]={0};double nn=0,bn=0;
        GfsFvPrimitive l={0},r={0};
        l.rho=1+.2*sin(.7*test);r.rho=.5+.15*cos(.4*test);l.p=1;r.p=.3;
        for(int a=0;a<dim;a++){
            n[a]=1+a;nn+=n[a]*n[a];
            l.u[a]=.8*sin(.7*(test+1)*(a+1));r.u[a]=.9*cos(.9*(test+1)*(a+1));
            b[a]=.7*(a+1);
        }
        for(int a=0;a<dim;a++){n[a]/=sqrt(nn);neg[a]=-n[a];bn+=b[a]*n[a];}
        if(test==0){for(int a=0;a<dim;a++){l.u[a]=-10*n[a];r.u[a]=10*n[a];}r.rho=l.rho;r.p=l.p;}
        if(test==1)r=l;
        GfsFvFan f=gfs_fv_fan(dim,l,r,n,gamma),rev=gfs_fv_fan(dim,r,l,neg,gamma);
        GfsFvPrimitive lb=l,rb=r;
        for(int a=0;a<dim;a++){lb.u[a]+=b[a];rb.u[a]+=b[a];}
        GfsFvFan boosted=gfs_fv_fan(dim,lb,rb,n,gamma);
        if(!f.valid||!rev.valid||!boosted.valid){bad=1;continue;}
        if(test==0&&!f.hll)bad=1;
        fallback+=f.hll;cases++;
        for(int s=0;s<7;s++){
            double wn=f.speed[0]+(s-.5)*(f.speed[2]-f.speed[0])/5;
            GfsFvSample q=gfs_fv_sample(&f,wn),qr=gfs_fv_sample(&rev,-wn),qb=gfs_fv_sample(&boosted,wn+bn);
            if(!q.valid||!qr.valid||!qb.valid){bad=1;continue;}
            for(int c=0;c<5;c++)reversal=fmax(reversal,discrepancy(q.flux[c],-qr.flux[c]));
            double transformed[5]={0};for(int c=0;c<5;c++)transformed[c]=q.flux[c];
            boost_state(dim,transformed,b);
            for(int c=0;c<5;c++)galilean=fmax(galilean,discrepancy(transformed[c],qb.flux[c]));
            double delta=1e-6*(1+fabs(wn));GfsFvSample plus=gfs_fv_sample(&f,wn+delta);
            if(plus.region==q.region)derivative=fmax(derivative,discrepancy((plus.flux[0]-q.flux[0])/delta,-q.density));
            if(test==1){
                double state[5]={0},flux[5]={0},un=0,v2=0;
                for(int a=0;a<dim;a++){un+=l.u[a]*n[a];v2+=l.u[a]*l.u[a];}
                state[0]=l.rho;flux[0]=l.rho*(un-wn);
                for(int a=0;a<dim;a++){state[a+1]=l.rho*l.u[a];flux[a+1]=flux[0]*l.u[a]+l.p*n[a];}
                state[4]=l.p/(gamma-1)+.5*l.rho*v2;
                flux[4]=(state[4]+l.p)*un-wn*state[4];
                for(int c=0;c<5;c++)uniform=fmax(uniform,discrepancy(flux[c],q.flux[c]));
            }
#ifdef __cplusplus
            if(dim==2&&!f.hll){
                GfsAlePrimitive a={l.rho,l.p,l.u[0],l.u[1]},z={r.rho,r.p,r.u[0],r.u[1]};
                GfsAleFan old=gfs_ale_fan_build(a,z,n[0],n[1],gamma);
                GfsAleFanSample result=gfs_ale_fan_sample(&old,wn);
                if(!result.valid||old.hll)bad=1;
                else for(int c=0;c<4;c++)old2d=fmax(old2d,discrepancy(result.flux[c],q.flux[c==3?4:c]));
            }
#endif
        }
    }
    /* Embedding the same planar Riemann data must leave active components
     * unchanged when the unused direction is introduced. */
    for(int dim=1;dim<=3;dim++){
        GfsFvPrimitive l={1,1,{.2,0,0}},r={.125,.1,{-.1,0,0}};
        double n[3]={1,0,0};GfsFvFan f=gfs_fv_fan(dim,l,r,n,gamma),one=gfs_fv_fan(1,l,r,n,gamma);
        for(int k=0;k<4;k++){
            GfsFvSample a=gfs_fv_sample(&f,-2+1.1*k),b=gfs_fv_sample(&one,-2+1.1*k);
            embedding=fmax(embedding,discrepancy(a.flux[0],b.flux[0]));
            embedding=fmax(embedding,discrepancy(a.flux[1],b.flux[1]));
            embedding=fmax(embedding,discrepancy(a.flux[4],b.flux[4]));
            for(int d=1;d<dim;d++)embedding=fmax(embedding,fabs(a.flux[d+1]));
        }
        for(int test=0;test<5;test++){
            double q[5]={0},v[5]={0},expected=0;q[0]=1;
            if(test==0){q[4]=2.5;v[0]=-1;v[4]=-2.5;expected=1;}
            if(test==1){q[4]=1;v[1]=1;expected=sqrt(2.);}
            if(test==2){q[1]=3;q[4]=5;v[4]=-1;expected=.5;}
            if(test==3){q[1]=1;q[4]=1;v[0]=1;v[1]=-1;expected=2+sqrt(5.);}
            if(test==4){q[4]=1;v[0]=-1;v[4]=-2;expected=.5;}
            double bound=gfs_fv_admissible_interval(dim,q,v),inside[5]={0},outside[5]={0};
            cone=fmax(cone,discrepancy(bound,expected));
            for(int c=0;c<5;c++){
                inside[c]=fma(.9*bound,v[c],q[c]);outside[c]=fma(1.01*bound,v[c],q[c]);
            }
            if(!(bound>0)||!(gfs_fv_internal(dim,inside)>0)||gfs_fv_internal(dim,outside)>0)bad=1;
            double boost[3]={.7,-.4,.2};boost_state(dim,q,boost);boost_state(dim,v,boost);
            cone=fmax(cone,discrepancy(bound,gfs_fv_admissible_interval(dim,q,v)));
        }
        double q[5]={0},zero[5]={0};q[0]=q[4]=1;
        if(!isinf(gfs_fv_admissible_interval(dim,q,zero)))bad=1;
        q[4]=-1;if(gfs_fv_admissible_interval(dim,q,zero)>=0)bad=1;
        if(dim<3){
            q[4]=1;zero[dim+1]=1;
            if(gfs_fv_admissible_interval(dim,q,zero)>=0)bad=1;
        }
    }
    const int old2d_checked=
#ifdef __cplusplus
        1;
#else
        0;
#endif
    int pass=!bad&&reversal<1e-12&&galilean<1e-12&&uniform<1e-12&&derivative<2e-8&&embedding<1e-12&&cone<1e-12&&old2d<1e-12;
    printf("{\"scope\":\"common 1D/2D/3D FV fan and forward-ray admissibility only; no PDE\","
           "\"fan_cases\":%d,\"HLL_fallback_cases\":%d,\"reversal_error\":%.17g,"
           "\"Galilean_flux_error\":%.17g,\"uniform_physical_flux_error\":%.17g,"
           "\"mass_derivative_error\":%.17g,\"dimensional_embedding_error\":%.17g,"
           "\"analytic_cone_boundary_error\":%.17g,\"original_2D_HLLC_error\":%.17g,\"original_2D_checked\":%s,"
           "\"invalid_count\":%d,\"PDE_steps\":0,\"pass\":%s}\n",cases,fallback,reversal,galilean,uniform,derivative,embedding,cone,old2d,old2d_checked?"true":"false",bad,pass?"true":"false");
    return pass?0:1;
}
