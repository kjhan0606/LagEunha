#ifndef GFS_FV_ALE_H
#define GFS_FV_ALE_H
/* Experimental dimensional Euler FV primitives. Production branches are
 * unchanged. A fixed payload [rho, rho*ux, rho*uy, rho*uz, E] embeds d+2
 * independent variables; unused momentum slots are exactly zero. It also
 * avoids runtime component remapping inside CUDA local aggregate arrays.
 * E-Kmean defines the declared FV thermal approximation;
 * no reconstructed variance, independent thermal state, or energy repair.
 * One frozen lab-frame fan supplies every ALE flux component F-w*U.
 */
#include <math.h>
#include <float.h>
#ifdef __CUDACC__
#define GFS_FV_HD static __host__ __device__ inline
#else
#define GFS_FV_HD static inline
#endif

typedef struct {double rho,p,u[3];} GfsFvPrimitive;
typedef struct {
    double speed[3],state[4][5],flux[4][5];
    int dim,valid,hll;
} GfsFvFan;
typedef struct {double flux[5],density;int valid,region;} GfsFvSample;

GFS_FV_HD double gfs_fv_internal(int dim,const double q[5])
{
    if(dim<1||dim>3||!(q[0]>0)||!isfinite(q[0]))return -1;
    if((dim<3&&q[3]!=0)||(dim<2&&q[2]!=0))return -1;
    double thermal=q[4];
    thermal=fma(-.5*q[1],q[1]/q[0],thermal);
    if(dim>=2)thermal=fma(-.5*q[2],q[2]/q[0],thermal);
    if(dim==3)thermal=fma(-.5*q[3],q[3]/q[0],thermal);
    return thermal;
}

GFS_FV_HD int gfs_fv_primitive(int dim,const double q[5],double volume,
        double gamma,GfsFvPrimitive *p)
{
    if(dim<1||dim>3||!(volume>0)||!isfinite(volume)||!(gamma>1)||!isfinite(gamma))return 0;
    double thermal=gfs_fv_internal(dim,q);
    if(!(thermal>0)||!isfinite(thermal))return 0;
    p->rho=q[0]/volume;p->p=(gamma-1)*thermal/volume;
    for(int a=0;a<3;a++)p->u[a]=a<dim?q[1+a]/q[0]:0;
    if(!(p->rho>0)||!(p->p>0)||!isfinite(p->rho)||!isfinite(p->p))return 0;
    for(int a=0;a<dim;a++)if(!isfinite(p->u[a]))return 0;
    return 1;
}

GFS_FV_HD int gfs_fv_input(int dim,GfsFvPrimitive p)
{
    if(dim<1||dim>3||!(p.rho>0)||!(p.p>0)||!isfinite(p.rho)||!isfinite(p.p))return 0;
    for(int a=0;a<dim;a++)if(!isfinite(p.u[a]))return 0;
    return 1;
}

GFS_FV_HD void gfs_fv_state_flux(int dim,GfsFvPrimitive p,const double normal[3],
        double gamma,double state[5],double flux[5])
{
    double un=0,kinetic=0;
    for(int a=0;a<dim;a++){un+=p.u[a]*normal[a];kinetic+=p.u[a]*p.u[a];}
    for(int a=0;a<5;a++)state[a]=flux[a]=0;
    state[0]=p.rho;flux[0]=p.rho*un;
    for(int a=0;a<dim;a++){
        state[a+1]=p.rho*p.u[a];flux[a+1]=state[a+1]*un+p.p*normal[a];
    }
    state[4]=p.p/(gamma-1)+.5*p.rho*kinetic;
    flux[4]=(state[4]+p.p)*un;
}

GFS_FV_HD GfsFvFan gfs_fv_fan(int dim,GfsFvPrimitive left,GfsFvPrimitive right,
        const double normal[3],double gamma)
{
    GfsFvFan fan={0};fan.dim=dim;
    if(!gfs_fv_input(dim,left)||!gfs_fv_input(dim,right)||!(gamma>1)||!isfinite(gamma))return fan;
    double unit=0,ul=0,ur=0;
    for(int a=0;a<dim;a++){
        if(!isfinite(normal[a]))return fan;
        unit+=normal[a]*normal[a];ul+=left.u[a]*normal[a];ur+=right.u[a]*normal[a];
    }
    if(fabs(unit-1)>1e-12)return fan;
    gfs_fv_state_flux(dim,left,normal,gamma,fan.state[0],fan.flux[0]);
    gfs_fv_state_flux(dim,right,normal,gamma,fan.state[3],fan.flux[3]);
    double cl=sqrt(gamma*left.p/left.rho),cr=sqrt(gamma*right.p/right.rho);
    double sl=fmin(ul-cl,ur-cr),sr=fmax(ul+cl,ur+cr);
    if(!(sr>sl)||!isfinite(sl)||!isfinite(sr))return fan;
    double al=left.rho*(sl-ul),ar=right.rho*(sr-ur);
    double sm=(right.p-left.p+al*ul-ar*ur)/(al-ar);
    double ps=left.p+al*(sm-ul);
    int good=isfinite(sm)&&isfinite(ps)&&ps>0&&sm>=sl&&sm<=sr;
    if(good)for(int side=0;side<2;side++){
        int src=side?3:0,dst=side?2:1;
        const double *q=fan.state[src];double *qs=fan.state[dst];
        double sk=side?sr:sl,uk=side?ur:ul,pk=side?right.p:left.p;
        double factor=q[0]*(sk-uk)/(sk-sm);
        qs[0]=factor;
        for(int a=0;a<dim;a++)qs[a+1]=factor*(q[a+1]/q[0]+(sm-uk)*normal[a]);
        qs[4]=factor*(q[4]/q[0]+(sm-uk)*(sm+pk/(q[0]*(sk-uk))));
        double thermal=gfs_fv_internal(dim,qs);
        if(!(thermal>0)||!isfinite(thermal)||!isfinite(qs[4]))good=0;
        for(int c=0;c<5;c++)fan.flux[dst][c]=fan.flux[src][c]+sk*(qs[c]-q[c]);
    }
    if(!good){
        /* Local Lax--Friedrichs outer bounds in the normal mean frame.
         * Keep the COMPLETE HLL numerical flux, including its residual. */
        double mean=.5*(ul+ur),a=fmax(fabs(ul-mean)+cl,fabs(ur-mean)+cr);
        sl=mean-a;sr=mean+a;sm=0;fan.hll=1;
        if(!(sr>sl)||!isfinite(sr)||!isfinite(sl))return fan;
        for(int c=0;c<5;c++){
            fan.state[1][c]=(sr*fan.state[3][c]-sl*fan.state[0][c]
                            -(fan.flux[3][c]-fan.flux[0][c]))/(sr-sl);
            fan.flux[1][c]=(sr*fan.flux[0][c]-sl*fan.flux[3][c]
                           +sl*sr*(fan.state[3][c]-fan.state[0][c]))/(sr-sl);
        }
        double thermal=gfs_fv_internal(dim,fan.state[1]);
        if(!(thermal>0)||!isfinite(thermal))return fan;
    }
    fan.speed[0]=sl;fan.speed[1]=sm;fan.speed[2]=sr;
    fan.valid=1;return fan;
}

GFS_FV_HD GfsFvSample gfs_fv_sample(const GfsFvFan *fan,double wn)
{
    GfsFvSample sample={0};if(!fan->valid||!isfinite(wn))return sample;
    int k=wn<=fan->speed[0]?0:wn>=fan->speed[2]?3:
          fan->hll?1:wn<=fan->speed[1]?1:2;
    for(int c=0;c<5;c++){
        sample.flux[c]=fma(-wn,fan->state[k][c],fan->flux[k][c]);
        if(!isfinite(sample.flux[c]))return sample;
    }
    sample.region=k;sample.density=fan->state[k][0];
    sample.valid=sample.density>0&&isfinite(sample.density);return sample;
}

/* First boundary of m(s)>0 and 2*m(s)*E(s)-|p(s)|^2>0 on the
 * forward-Euler ray Q+s*rate. A negative return denotes invalid input;
 * infinity denotes no finite forward boundary. Callers apply a fixed safety
 * factor and validate the actual stage and geometric V, never repair Q.
 * This is a floating-point calculation, not an interval certificate or GCL.
 */
GFS_FV_HD double gfs_fv_admissible_interval(int dim,const double q[5],const double rate[5])
{
    if(!(gfs_fv_internal(dim,q)>0))return -1;
    if((dim<3&&rate[3]!=0)||(dim<2&&rate[2]!=0))return -2;
    double qs=0,rs=0;
    for(int c=0;c<5;c++){
        if(!isfinite(q[c])||!isfinite(rate[c]))return -2;
        qs=fmax(qs,fabs(q[c]));rs=fmax(rs,fabs(rate[c]));
    }
    if(rs==0)return HUGE_VAL;
    double x[5]={0},v[5]={0};
    for(int c=0;c<5;c++){x[c]=q[c]/qs;v[c]=rate[c]/rs;}
    double C=2*x[0]*gfs_fv_internal(dim,x);
    double B=2*fma(x[0],v[4],v[0]*x[4]);
    double A=2*v[0]*v[4];
    for(int a=0;a<dim;a++){
        A=fma(-v[1+a],v[1+a],A);B=fma(-2*x[1+a],v[1+a],B);
    }
    if(!(C>0)||!isfinite(A)||!isfinite(B)||!isfinite(C))return -3;
    double boundary=v[0]<0?-x[0]/v[0]:HUGE_VAL;
    if(A==0){if(B<0)boundary=fmin(boundary,-C/B);}
    else if(A<0||B<0){
        double disc=fma(B,B,-4*A*C);
        if(disc>=0){
            double z=-.5*(B+copysign(sqrt(disc),B));
            if(z!=0){
                double r0=z/A,r1=C/z;
                if(r0>0)boundary=fmin(boundary,r0);
                if(r1>0)boundary=fmin(boundary,r1);
            }
        }
    }
    if(isinf(boundary))return HUGE_VAL;
    return boundary*(qs/rs);
}

#undef GFS_FV_HD
#endif
