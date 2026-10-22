/* Algebra audit of Astra's proposed unresolved-kinetic tensor closure.
 * This is not a simulation and is not installed in the hydro integrator. */
#include <cmath>
#include <cstdio>
#include <cstdlib>
struct Q {double rho,p,vx,vy,xx,xy,yy;};
static void flux(Q q,double nx,double ny,double wn,double gamma,double u[4],double f[4]){
    double vn=q.vx*nx+q.vy*ny,a=vn-wn;
    double tx=q.xx*nx+q.xy*ny,ty=q.xy*nx+q.yy*ny;
    u[1]=q.rho*(q.vx*q.vx+q.xx);u[2]=q.rho*(q.vx*q.vy+q.xy);
    u[3]=q.rho*(q.vy*q.vy+q.yy);u[0]=.5*(u[1]+u[3])+q.p/(gamma-1);
    f[0]=u[0]*a+q.p*vn+q.rho*(q.vx*tx+q.vy*ty);
    f[1]=u[1]*a+2*q.rho*q.vx*tx+2*q.p*q.vx*nx;
    f[2]=u[2]*a+q.rho*(q.vx*ty+q.vy*tx)+q.p*(q.vx*ny+q.vy*nx);
    f[3]=u[3]*a+2*q.rho*q.vy*ty+2*q.p*q.vy*ny;
}
int main(){
    const double gamma=5./3.;double max_physical=0,max_hll=0,max_source=0;
    for(int k=0;k<100;k++){
        Q l={1.2,1.e-6,.2+.01*k,-.4,.07,.02,.03};
        Q r={.7,2.e-6,-.3,.1-.007*k,.02,-.01,.04};
        double angle=.07*k,nx=cos(angle),ny=sin(angle),wn=.03*k-.7;
        double ul[4],ur[4],fl[4],fr[4];flux(l,nx,ny,wn,gamma,ul,fl);flux(r,nx,ny,wn,gamma,ur,fr);
        double a=l.vx*nx+l.vy*ny-wn,expected=l.p/(gamma-1)*a;
        max_physical=fmax(max_physical,fabs(fl[0]-.5*(fl[1]+fl[3])-expected));
        const double sl=-5,sr=5;double h[4];
        for(int j=0;j<4;j++)h[j]=(sr*fl[j]-sl*fr[j]+sl*sr*(ur[j]-ul[j]))/(sr-sl);
        double b=r.vx*nx+r.vy*ny-wn;
        double hu=(sr*l.p*a-sl*r.p*b+sl*sr*(r.p-l.p))/((gamma-1)*(sr-sl));
        max_hll=fmax(max_hll,fabs(h[0]-.5*(h[1]+h[3])-hu));
        double div=.1-.003*k,gx=.03,gy=-.07;
        double se=l.rho*(l.vx*gx+l.vy*gy);
        double spi_trace=2*l.p*div+2*se;
        max_source=fmax(max_source,fabs(se-.5*spi_trace+l.p*div));
    }
    if(max_physical>1.e-14 || max_hll>1.e-14 || max_source>1.e-14)return 1;
    printf("{\"status\":\"identity_confirmed\",\"cases\":100,\"physical_flux_residual\":%.17g,\"HLL_flux_residual\":%.17g,\"source_residual\":%.17g,\"implication\":\"U=E-half_trace_Pi has only thermal advection and pressure-divergence work; unresolved kinetic dissipation has no R-to-U thermalization term\",\"scope\":\"algebra only, not a shock simulation\"}\n",max_physical,max_hll,max_source);
}
