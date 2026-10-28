#ifndef LAGFORCE_MOMENT_QUADRUPOLE_GRAVITY_H
#define LAGFORCE_MOMENT_QUADRUPOLE_GRAVITY_H
#include <math.h>
#ifdef __CUDACC__
#define MQ_INLINE static __host__ __device__ __forceinline__
#else
#define MQ_INLINE static inline
#endif
/* Experimental second-order external-potential Hamiltonian, GM=1.
 * Psi=m*phi(c)+H_phi(c):S/2, S=Q-m*c*c^T.
 * F=-m*grad(phi)-grad(H):S/2; affine acceleration A=-H.
 * J=c*F^T+S*A^T, work=F.vbar+trace(A*C), C=D-q*p^T/m.
 * The third radial derivative in F is essential for torque and work.
 * Input17 raw moments, rc. Output F[2],J[4],work,Psi,A[4].
 * One radial evaluation per cell; no spatial quadrature or matrix inverse.
 * This is a truncated Hamiltonian, NOT exact integrated gravity or a proof
 * of equivalence to the old generator-point force. Production disabled.
 */
MQ_INLINE int mq_source(const double *a,double *b) {
    double m=a[0],rc=a[17];if(!(m>0 && rc>0))return 1;
    double x=a[1]/m,y=a[2]/m,r=sqrt(x*x+y*y),u=r/rc;
    double phi,k,l,t;
    if(u>=1) {
        double ir=1/r,ir2=ir*ir;
        phi=-ir;k=ir*ir2;l=-3*k*ir2;t=-5*l*ir2;
    } else {
        double irc=1/rc,irc2=irc*irc,irc3=irc*irc2,irc5=irc3*irc2,irc7=irc5*irc2;
        if(u<=0.5) {
            phi=(16.0/3*u*u-48.0/5*u*u*u*u+32.0/5*u*u*u*u*u-14.0/5)*irc;
            k=(32.0/3-192.0/5*u*u+32*u*u*u)*irc3;
            l=(-384.0/5+96*u)*irc5;t=u>0 ? 96/u*irc7 : 0;
        } else {
            double iu=1/u,iu2=iu*iu,iu3=iu2*iu,iu5=iu3*iu2,iu7=iu5*iu2;
            phi=(iu/15+32.0/3*u*u-16*u*u*u+48.0/5*u*u*u*u-32.0/15*u*u*u*u*u-16.0/5)*irc;
            k=(-iu3/15+64.0/3-48*u+192.0/5*u*u-32.0/3*u*u*u)*irc3;
            l=(iu5/5-48*iu+384.0/5-32*u)*irc5;
            t=(-iu7+48*iu3-32*iu)*irc7;
        }
    }
    double sx=a[3]-m*x*x,sxy=a[4]-m*x*y,sy=a[5]-m*y*y;
    double scx=sx*x+sxy*y,scy=sxy*x+sy*y,trace=sx+sy,radial=x*scx+y*scy;
    double factor=-m*k-0.5*(l*trace+t*radial);
    double fx=factor*x-l*scx,fy=factor*y-l*scy;
    double A00=-k-l*x*x,A01=-l*x*y,A11=-k-l*y*y;
    double vx=a[6]/m,vy=a[7]/m;
    double C00=a[8]-a[1]*vx,C01=a[9]-a[1]*vy,C10=a[10]-a[2]*vx,C11=a[11]-a[2]*vy;
    b[0]=fx;b[1]=fy;
    b[2]=x*fx+sx*A00+sxy*A01;b[3]=x*fy+sx*A01+sxy*A11;
    b[4]=y*fx+sxy*A00+sy*A01;b[5]=y*fy+sxy*A01+sy*A11;
    b[6]=fx*vx+fy*vy+A00*C00+A01*(C01+C10)+A11*C11;
    b[7]=m*phi+0.5*(k*trace+l*radial);
    b[8]=A00;b[9]=A01;b[10]=A01;b[11]=A11;
    return 0;
}
#undef MQ_INLINE
#endif
