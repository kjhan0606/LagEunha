#ifndef GFS_ALE_FAN_H
#define GFS_ALE_FAN_H
#include "gfs_ale_flux.h"
#ifdef __CUDACC__
#define GFS_FAN_HD __host__ __device__
#else
#define GFS_FAN_HD
#endif
/* One fixed LAB-FRAME HLLC/HLL Riemann fan, sampled at normal mesh speed w.
 * States/reconstruction stay fixed inside a weight-rate solve. All four
 * components use F(w)=F_k-w U_k. Consequently d F_mass/dw=-rho_k<0.
 * Re-solving a translated Riemann problem at every Newton iteration can
 * obscure this property with frame-dependent branch roundoff.
 */
typedef struct {
    double speed[3],u[4][4],f[4][4]; /* L,Lstar,Rstar,R; HLL uses slot1 */
    int valid,hll;
} GfsAleFan;
typedef struct {
    double flux[4],density;
    int valid,mode,region;
} GfsAleFanSample;

GFS_FAN_HD static inline double gfs_fan_internal(const double u[4]){
    if(!(u[0]>0))return -1;
    return fma(-.5*u[2],u[2]/u[0],fma(-.5*u[1],u[1]/u[0],u[3]));
}

GFS_FAN_HD static inline GfsAleFan gfs_ale_fan_build(GfsAlePrimitive l,
    GfsAlePrimitive r,double nx,double ny,double gamma){
    GfsAleFan fan={};
    if(!gfs_ale_primitive_valid(l) || !gfs_ale_primitive_valid(r) || !(gamma>1) ||
       fabs(nx*nx+ny*ny-1)>1.e-12)return fan;
    gfs_ale_state(l,0,0,nx,ny,gamma,fan.u[0],fan.f[0]);
    gfs_ale_state(r,0,0,nx,ny,gamma,fan.u[3],fan.f[3]);
    double vl=l.vx*nx+l.vy*ny,vr=r.vx*nx+r.vy*ny;
    double cl=sqrt(gamma*l.p/l.rho),cr=sqrt(gamma*r.p/r.rho);
    double sl=fmin(vl-cl,vr-cr),sr=fmax(vl+cl,vr+cr);
    if(!(sr>sl) || !isfinite(sl) || !isfinite(sr))return fan;
    fan.speed[0]=sl;fan.speed[2]=sr;
    double al=l.rho*(sl-vl),ar=r.rho*(sr-vr),den=al-ar;
    double sm=(r.p-l.p+al*vl-ar*vr)/den,ps=l.p+al*(sm-vl);
    int good=isfinite(sm) && ps>0 && isfinite(ps) && sm>=sl && sm<=sr;
    fan.speed[1]=sm;
    if(good)for(int side=0;side<2;side++){
        int src=side?3:0,dst=side?2:1;
        double sk=side?sr:sl,vk=side?vr:vl,pk=side?r.p:l.p;
        const double *q=fan.u[src];double *qs=fan.u[dst];
        double factor=q[0]*(sk-vk)/(sk-sm);
        qs[0]=factor;qs[1]=factor*(q[1]/q[0]+(sm-vk)*nx);
        qs[2]=factor*(q[2]/q[0]+(sm-vk)*ny);
        qs[3]=factor*(q[3]/q[0]+(sm-vk)*(sm+pk/(q[0]*(sk-vk))));
        if(!(gfs_fan_internal(qs)>0) || !isfinite(qs[3]))good=0;
        for(int c=0;c<4;c++)fan.f[dst][c]=fan.f[src][c]+sk*(qs[c]-q[c]);
    }
    if(!good){
        fan.hll=1;fan.speed[1]=0;
        for(int c=0;c<4;c++){
            fan.u[1][c]=(sr*fan.u[3][c]-sl*fan.u[0][c]-(fan.f[3][c]-fan.f[0][c]))/(sr-sl);
            fan.f[1][c]=(sr*fan.f[0][c]-sl*fan.f[3][c]+sl*sr*(fan.u[3][c]-fan.u[0][c]))/(sr-sl);
        }
        if(!(gfs_fan_internal(fan.u[1])>0) || !isfinite(fan.u[1][3]))return fan;
    }
    fan.valid=1;return fan;
}

GFS_FAN_HD static inline GfsAleFanSample gfs_ale_fan_sample(const GfsAleFan *fan,double wn){
    GfsAleFanSample q={};if(!fan->valid || !isfinite(wn))return q;
    int k,mode;
    if(wn<=fan->speed[0]){k=0;mode=1;}
    else if(wn>=fan->speed[2]){k=3;mode=2;}
    else if(fan->hll){k=1;mode=5;}
    else if(wn<=fan->speed[1]){k=1;mode=3;}
    else{k=2;mode=4;}
    q.density=fan->u[k][0];q.region=k;q.mode=mode;
    for(int c=0;c<4;c++){
        q.flux[c]=fma(-wn,fan->u[k][c],fan->f[k][c]);
        if(!isfinite(q.flux[c]))return q;
    }
    q.valid=q.density>0;return q;
}

/* A canonical face quadrature point, outward i->j. w0 is the actual
 * geometric normal speed with Wdot=0; inv2d=1/(2|xj-xi|). The same fan/flux
 * supplies every Euler moment and both signs of the mass ledger. */
typedef struct {
    GfsAleFan fan;
    double weight,w0,inv2d,xi,yi,xj,yj;
    int i,j;
} GfsWeightPoint;

GFS_FAN_HD static inline GfsAleFanSample gfs_weight_flux(const GfsWeightPoint *q,
    double zi,double zj){
    return gfs_ale_fan_sample(&q->fan,q->w0+q->inv2d*(zi-zj));
}

/* Residual R_i = mu_i z_i - sum_face A Fmass(w0+(zi-zj)/(2d)).
 * Pair Jacobian is a[[1,-1],[-1,1]], a=A*rho/(2d)>0.
 * With every mu_i>0 the complete Jacobian is symmetric positive definite,
 * including HLL fallback. No hard fixed-mass DAE or mass-flux deletion. */
GFS_FAN_HD static inline double gfs_weight_jacobian(const GfsWeightPoint *q,
    double zi,double zj){
    GfsAleFanSample f=gfs_weight_flux(q,zi,zj);
    return f.valid?q->weight*q->inv2d*f.density:-1;
}
#undef GFS_FAN_HD
#endif
