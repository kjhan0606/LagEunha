#ifndef MATERIAL_NATIVE_FORWARD_H
#define MATERIAL_NATIVE_FORWARD_H
#include "material_forward_ball.h"
namespace material_forward {
/* Mirrors the actual frozen helper's operation order, branch/series stop,
 * shape factors, sums, inverse and origin addition. Native .n values must
 * equal curved::native_point at every saved point. */
static Point<double> native(const curved::Record &r,double px,double py,double eps,bool coordinate_cast){
    using T=double;using S=B<T>;S X(px,coordinate_cast?cast_radius(px):0),Y(py,coordinate_cast?cast_radius(py):0),ax(r.triangle.geometry.x[0][0]),ay(r.triangle.geometry.x[0][1]);
    S dx=X-ax,dy=Y-ay,z=ax*ax+ay*ay+S(eps)*S(eps),L=S(2)*(ax*dx+ay*dy)/z,Q=(dx*dx+dy*dy)/z,s=L+Q;
    T abs_s=fabs(s.n);need(abs_s<.75,"same native remainder branch");S tail,dtail,term;int last=4;
    if(s.n!=0){
        term=S(.56396484375)*s*s*s*s;tail=term;dtail=S(4)*term/s;bool converged=false;
        for(int n=4;n<256;n++){
            T value_bound=fabs(term.n)*abs_s/(1-abs_s),derivative_bound=fabs(term.n)*((n+1)/(1-abs_s)+abs_s/((1-abs_s)*(1-abs_s)));
            if(value_bound<=DBL_EPSILON*fabs(tail.n) && derivative_bound<=DBL_EPSILON*fabs(dtail.n)){converged=true;last=n;break;}
            term=term*((S(-.75-n)/S(n+1))*s);tail=tail+term;dtail=dtail+S(n+1)*term/s;
        }need(converged,"same native remainder convergence");
        /* Exact infinite-series remainder after the actual stop, covering
         * both the represented value and the analytic derivative. */
        T sh=FP<T>::up(FP<T>::abs(s.n)+s.e),th=FP<T>::up(FP<T>::abs(term.n)+term.e);need(sh<.75,"whole forward series domain");
        T den=FP<T>::down(T(1)-sh),vt=FP<T>::up(times_up(th,sh)/den),dt=times_up(th,plus_up(FP<T>::up(T(last+1)/den),FP<T>::up(sh/FP<T>::down(den*den))));tail.e=plus_up(tail.e,vt);dtail.e=plus_up(dtail.e,dt);
    }else{
        // A rounded zero may represent a nonzero exact s. Cover the
        // analytic tail about zero while preserving the actual no-loop n.
        T sh=s.e;need(sh<.75,"zero-centered forward series domain");
        T den=FP<T>::down(T(1)-sh),cube=times_up(times_up(sh,sh),sh),fourth=times_up(cube,sh);
        tail.e=FP<T>::up(times_up(T(.56396484375),fourth)/den);
        dtail.e=times_up(times_up(T(.56396484375),cube),plus_up(FP<T>::up(T(4)/den),FP<T>::up(sh/FP<T>::down(den*den))));
    }
    S A=S(.65625)*Q*Q+S(-.6015625)*(S(3)*L*L*Q+S(3)*L*Q*Q+Q*Q*Q)+tail;
    S Bv=S(.65625)*(S(2)*L*Q+Q*Q)+S(-.6015625)*s*s*s+tail;
    S AL=S(-.6015625)*(S(6)*L*Q+S(3)*Q*Q)+dtail,AQ=S(2)*S(.65625)*Q+S(-.6015625)*(S(3)*L*L+S(6)*L*Q+S(3)*Q*Q)+dtail;
    S BL=S(2)*S(.65625)*Q+S(3)*S(-.6015625)*s*s+dtail,BQ=S(2)*S(.65625)*(L+Q)+S(3)*S(-.6015625)*s*s+dtail;
    S omega=power(z,T(-.75)),ja[2]={-ay,ax},jd[2]={-dy,dx},q[2],Dq[4];
    for(int k=0;k<2;k++){q[k]=omega*(ja[k]*A+jd[k]*Bv);for(int d=0;d<2;d++){S ld=S(2)*(d?ay:ax)/z,qd=S(2)*(d?dy:dx)/z,dA=AL*ld+AQ*qd,dB=BL*ld+BQ*qd,dJ((k==0 && d==1)?-1:(k==1 && d==0)?1:0);Dq[2*k+d]=omega*(ja[k]*dA+jd[k]*dB+dJ*Bv);}}
    S grad[3][2];for(int k=0;k<3;k++)for(int d=0;d<2;d++)grad[k][d]=S(r.triangle.geometry.grad[k][d]);
    for(int d=0;d<2;d++){auto mismatch=curved_bounds::I(r.triangle.geometry.grad[0][d])+curved_bounds::I(r.triangle.geometry.grad[1][d])+curved_bounds::I(r.triangle.geometry.grad[2][d]);grad[0][d].e=FP<T>::up(T(curved_bounds::absmax(mismatch)));}
    S l[3],N[12],Nx[12],Ny[12];l[1]=grad[1][0]*dx+grad[1][1]*dy;l[2]=grad[2][0]*dx+grad[2][1]*dy;l[0]=S(1)-l[1]-l[2];
    for(int a=0;a<10;a++){int id[3];gfs_material_p3_index(a,id);S p[3],d[3];for(int k=0;k<3;k++){int n=id[k];if(n==0){p[k]=S(1);d[k]=S(0);}else if(n==1){p[k]=S(3)*l[k];d[k]=S(3);}else if(n==2){p[k]=S(1.5)*l[k]*(S(3)*l[k]-S(1));d[k]=S(9)*l[k]-S(1.5);}else{p[k]=S(.5)*l[k]*(S(3)*l[k]-S(1))*(S(3)*l[k]-S(2));d[k]=(S(13.5)*l[k]-S(9))*l[k]+S(1);}}N[a]=p[0]*p[1]*p[2];for(int k=0;k<3;k++){S v=d[k]*p[(k+1)%3]*p[(k+2)%3];Nx[a]=Nx[a]+v*grad[k][0];Ny[a]=Ny[a]+v*grad[k][1];}}
    for(int k=0;k<2;k++){S rem=q[k]-S(r.triangle.remainder[0][k]),rx=Dq[2*k],ry=Dq[2*k+1];for(int a=1;a<10;a++){S dq=S(r.triangle.remainder[a][k])-S(r.triangle.remainder[0][k]);rem=rem-N[a]*dq;rx=rx-Nx[a]*dq;ry=ry-Ny[a]*dq;}N[10+k]=rem/S(r.scale[k]);Nx[10+k]=rx/S(r.scale[k]);Ny[10+k]=ry/S(r.scale[k]);}
    Point<T> p={};for(int a=0;a<r.map.n;a++)for(int k=0;k<2;k++){p.position[k]=p.position[k]+S(r.map.x[k][a])*N[a];p.velocity[k]=p.velocity[k]+S(r.map.v[k][a])*N[a];p.F[2*k]=p.F[2*k]+S(r.map.x[k][a])*Nx[a];p.F[2*k+1]=p.F[2*k+1]+S(r.map.x[k][a])*Ny[a];p.V[2*k]=p.V[2*k]+S(r.map.v[k][a])*Nx[a];p.V[2*k+1]=p.V[2*k+1]+S(r.map.v[k][a])*Ny[a];}finish(p);for(int k=0;k<2;k++)p.position[k]=p.position[k]+S(r.origin[k]);return p;
}
}
#endif
