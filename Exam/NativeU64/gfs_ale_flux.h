#ifndef GFS_ALE_FLUX_H
#define GFS_ALE_FLUX_H
/* Diagnostic conservative Euler flux in an arbitrary translating face frame.
 * HLLC with Davis bounds; admissibility failure falls back to the complete
 * HLL flux, never to a pressure-only average. No state is evolved here.
 * All four conserved components use the same Riemann problem.
 */
#include <math.h>
#ifdef __CUDACC__
#define GFS_ALE_HD __host__ __device__
#else
#define GFS_ALE_HD
#endif
typedef struct { double rho,p,vx,vy; } GfsAlePrimitive;
typedef struct {
    double flux[4],pressure[4],advection[4],residual[4];
    int valid,mode;
} GfsAleFlux;

GFS_ALE_HD static inline int gfs_ale_primitive_valid(GfsAlePrimitive q){
    return q.rho>0 && q.p>0 && isfinite(q.rho) && isfinite(q.p) &&
        isfinite(q.vx) && isfinite(q.vy);
}

GFS_ALE_HD static inline void gfs_ale_state(GfsAlePrimitive q,double wx,double wy,
    double nx,double ny,double gamma,double u[4],double f[4]){
    double vx=q.vx-wx,vy=q.vy-wy,vn=vx*nx+vy*ny;
    u[0]=q.rho;u[1]=q.rho*vx;u[2]=q.rho*vy;
    u[3]=q.p/(gamma-1)+.5*q.rho*(vx*vx+vy*vy);
    f[0]=q.rho*vn;f[1]=u[1]*vn+q.p*nx;
    f[2]=u[2]*vn+q.p*ny;f[3]=(u[3]+q.p)*vn;
}

GFS_ALE_HD static inline GfsAleFlux gfs_ale_flux(GfsAlePrimitive l,GfsAlePrimitive r,
    double nx,double ny,double wx,double wy,double gamma){
    GfsAleFlux out={};
    if(!gfs_ale_primitive_valid(l) || !gfs_ale_primitive_valid(r) || !(gamma>1)) return out;
    double ql[4],qr[4],fl[4],fr[4],qs[4],fs[4];
    gfs_ale_state(l,wx,wy,nx,ny,gamma,ql,fl);
    gfs_ale_state(r,wx,wy,nx,ny,gamma,qr,fr);
    double ul=(l.vx-wx)*nx+(l.vy-wy)*ny,ur=(r.vx-wx)*nx+(r.vy-wy)*ny;
    double cl=sqrt(gamma*l.p/l.rho),cr=sqrt(gamma*r.p/r.rho);
    double sl=fmin(ul-cl,ur-cr),sr=fmax(ul+cl,ur+cr),ps=0;
    if(sl>=0){
        for(int k=0;k<4;k++){qs[k]=ql[k];fs[k]=fl[k];}
        ps=l.p;out.mode=1;
    }else if(sr<=0){
        for(int k=0;k<4;k++){qs[k]=qr[k];fs[k]=fr[k];}
        ps=r.p;out.mode=2;
    }else{
        double al=l.rho*(sl-ul),ar=r.rho*(sr-ur),den=al-ar;
        double sm=(r.p-l.p+al*ul-ar*ur)/den;
        ps=l.p+al*(sm-ul);
        int good=isfinite(sm) && ps>0 && isfinite(ps) && sm>=sl && sm<=sr;
        if(good){
            const double *q=sm>=0?ql:qr,*f=sm>=0?fl:fr;
            double sk=sm>=0?sl:sr,uk=sm>=0?ul:ur,pk=sm>=0?l.p:r.p;
            double rho=q[0],factor=rho*(sk-uk)/(sk-sm);
            qs[0]=factor;
            qs[1]=factor*(q[1]/rho+(sm-uk)*nx);
            qs[2]=factor*(q[2]/rho+(sm-uk)*ny);
            qs[3]=factor*(q[3]/rho+(sm-uk)*(sm+pk/(rho*(sk-uk))));
            double eint=qs[3]-.5*(qs[1]*qs[1]+qs[2]*qs[2])/qs[0];
            good=qs[0]>0 && eint>0 && isfinite(qs[0]) && isfinite(eint);
            if(good){
                for(int k=0;k<4;k++) fs[k]=f[k]+sk*(qs[k]-q[k]);
                out.mode=sm>=0?3:4;
            }
        }
        if(!good){
            for(int k=0;k<4;k++){
                qs[k]=(sr*qr[k]-sl*ql[k]-(fr[k]-fl[k]))/(sr-sl);
                fs[k]=(sr*fl[k]-sl*fr[k]+sl*sr*(qr[k]-ql[k]))/(sr-sl);
            }
            double eint=qs[3]-.5*(qs[1]*qs[1]+qs[2]*qs[2])/qs[0];
            ps=(gamma-1)*eint;
            if(!(qs[0]>0) || !(ps>0) || !isfinite(ps)) return out;
            out.mode=5;
        }
    }
    /* Transform the *flux*, including mass, from the moving frame. */
    out.flux[0]=fs[0];out.flux[1]=fs[1]+wx*fs[0];out.flux[2]=fs[2]+wy*fs[0];
    out.flux[3]=fs[3]+wx*fs[1]+wy*fs[2]+.5*(wx*wx+wy*wy)*fs[0];
    double vx=qs[1]/qs[0]+wx,vy=qs[2]/qs[0]+wy;
    double elab=qs[3]+wx*qs[1]+wy*qs[2]+.5*(wx*wx+wy*wy)*qs[0];
    double wn=wx*nx+wy*ny;
    out.pressure[0]=0;out.pressure[1]=ps*nx;out.pressure[2]=ps*ny;out.pressure[3]=ps*wn;
    out.advection[0]=out.flux[0];out.advection[1]=out.flux[0]*vx;
    out.advection[2]=out.flux[0]*vy;out.advection[3]=out.flux[0]*(elab+ps)/qs[0];
    for(int k=0;k<4;k++){
        out.residual[k]=out.flux[k]-out.pressure[k]-out.advection[k];
        if(!isfinite(out.flux[k]) || !isfinite(out.residual[k])) return GfsAleFlux{};
    }
    /* HLLC residual is roundoff. HLL fallback can have a nonzero difference
     * between its numerical flux and F(its intermediate averaged state).
     * Keep that difference explicit; do not call it physical advection. */
    out.valid=1;return out;
}
#undef GFS_ALE_HD
#endif
