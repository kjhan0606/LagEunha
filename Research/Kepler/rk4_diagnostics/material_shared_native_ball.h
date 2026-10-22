#ifndef MATERIAL_SHARED_NATIVE_BALL_H
#define MATERIAL_SHARED_NATIVE_BALL_H
#include "material_shared_function.h"
// Native radial operation graph retained from material_native_forward.h;
// same stop and analytic infinite-series remainder, no new libm assumption.
namespace shared_native_ball {
namespace sf=shared_function;namespace mf=material_forward;
using B=mf::B<double>;using D=mf::D<double>;
static std::array<D,2> radial(const double*anchor,D Xj,D Yj,double eps){
    using T=double;using S=B;
    S X=Xj.v,Y=Yj.v,ax(anchor[0]),ay(anchor[1]);
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
        T sh=mf::FP<T>::up(mf::FP<T>::abs(s.n)+s.e),th=mf::FP<T>::up(mf::FP<T>::abs(term.n)+term.e);need(sh<.75,"whole forward series domain");
        T den=mf::FP<T>::down(T(1)-sh),vt=mf::FP<T>::up(mf::times_up(th,sh)/den),dt=mf::times_up(th,mf::plus_up(mf::FP<T>::up(T(last+1)/den),mf::FP<T>::up(sh/mf::FP<T>::down(den*den))));tail.e=mf::plus_up(tail.e,vt);dtail.e=mf::plus_up(dtail.e,dt);
    }else{
        // A rounded zero may represent a nonzero exact s. Cover the
        // analytic tail about zero while preserving the actual no-loop n.
        T sh=s.e;need(sh<.75,"zero-centered forward series domain");
        T den=mf::FP<T>::down(T(1)-sh),cube=mf::times_up(mf::times_up(sh,sh),sh),fourth=mf::times_up(cube,sh);
        tail.e=mf::FP<T>::up(mf::times_up(T(.56396484375),fourth)/den);
        dtail.e=mf::times_up(mf::times_up(T(.56396484375),cube),mf::plus_up(mf::FP<T>::up(T(4)/den),mf::FP<T>::up(sh/mf::FP<T>::down(den*den))));
    }
    S A=S(.65625)*Q*Q+S(-.6015625)*(S(3)*L*L*Q+S(3)*L*Q*Q+Q*Q*Q)+tail;
    S Bv=S(.65625)*(S(2)*L*Q+Q*Q)+S(-.6015625)*s*s*s+tail;
    S AL=S(-.6015625)*(S(6)*L*Q+S(3)*Q*Q)+dtail,AQ=S(2)*S(.65625)*Q+S(-.6015625)*(S(3)*L*L+S(6)*L*Q+S(3)*Q*Q)+dtail;
    S BL=S(2)*S(.65625)*Q+S(3)*S(-.6015625)*s*s+dtail,BQ=S(2)*S(.65625)*(L+Q)+S(3)*S(-.6015625)*s*s+dtail;
    S omega=mf::power(z,T(-.75)),ja[2]={-ay,ax},jd[2]={-dy,dx},q[2],Dq[4];
    for(int k=0;k<2;k++){q[k]=omega*(ja[k]*A+jd[k]*Bv);for(int d=0;d<2;d++){S ld=S(2)*(d?ay:ax)/z,qd=S(2)*(d?dy:dx)/z,dA=AL*ld+AQ*qd,dB=BL*ld+BQ*qd,dJ((k==0 && d==1)?-1:(k==1 && d==0)?1:0);Dq[2*k+d]=omega*(ja[k]*dA+jd[k]*dB+dJ*Bv);}}
    std::array<D,2>out;for(int k=0;k<2;k++)out[k]=D(q[k],Dq[2*k]*Xj.x+Dq[2*k+1]*Yj.x,Dq[2*k]*Xj.y+Dq[2*k+1]*Yj.y);return out;
}
static D ipow(D a,int n){D b(B(1));for(int j=0;j<n;j++)b=b*a;return b;}
static D lagrange(int node,D t){D p(B(1));for(int j=0;j<4;j++)if(j!=node)p=p*(D(B(3))*t-D(B(j)))/D(B(node-j));return p;}
static std::array<D,2> edge_raw(const sf::Edge&e,D t){
    D X=D(B(e.native.A[0]))+t*D(e.Dcache[0]),Y=D(B(e.native.A[1]))+t*D(e.Dcache[1]);auto out=radial(e.native.A,X,Y,e.native.eps);for(int j=0;j<4;j++){auto L=lagrange(j,t);for(int k=0;k<2;k++)out[k]=out[k]-L*D(e.qcache[j][k]);}if((t.v.n==0||t.v.n==1)&&t.v.e==0)for(int k=0;k<2;k++)out[k].v=B(0);return out;
}
static mf::Point<double> evaluate(const sf::Map&m,const STPoint&q,sf::S X,sf::S Y){
    const auto&t=m.state->triangle[m.source];need(size_t(q.triangle)==m.source,"native ball source identity");auto original=sf::barycentric(m,X,Y);D l[3],N[ST_MAX];for(int j=0;j<3;j++)l[j]=D(sf::cache(q.l[j],original[j].v),m.gradcache[j][0],m.gradcache[j][1]);
    for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);N[a]=D(B(1));for(int k=0;k<3;k++)for(int j=0;j<index[k];j++)N[a]=N[a]*(D(B(3))*l[k]-D(B(j)))/D(B(index[k]-j));}
    auto raw=radial(t.vertex[0],D(sf::cache(q.X,X),B(1),B(0)),D(sf::cache(q.Y,Y),B(0),B(1)),t.eps);for(int a=0;a<10;a++)for(int k=0;k<2;k++)raw[k]=raw[k]-N[a]*D(m.qcache[a][k]);
    for(int e=0;e<3;e++){D s=l[t.endpoint[e][0]]+l[t.endpoint[e][1]],extension[2];if(s.v.n!=0){auto r=edge_raw(m.edge[e],l[t.endpoint[e][1]]/s);for(int k=0;k<2;k++)extension[k]=ipow(s,4)*r[k];}else need(s.v.e==0,"native opposite-vertex branch needs canonical vertex chart");for(int k=0;k<2;k++)raw[k]=raw[k]-extension[k];for(int j=0;j<t.edge[e].rank;j++){N[t.offset[e]+j]=D();for(int k=0;k<2;k++)N[t.offset[e]+j]=N[t.offset[e]+j]+D(B(t.edge[e].encode[j][k]))*extension[k];}}
    if(q.boundary>=0)for(int k=0;k<2;k++)raw[k].v=B(0);for(int k=0;k<2;k++)raw[k]=raw[k]-D(B(t.bubble_shift[k]))*N[9];for(int j=0;j<2;j++){N[t.bubble_offset+j]=D();for(int k=0;k<2;k++)N[t.bubble_offset+j]=N[t.bubble_offset+j]+D(B(t.bubble_encode[j][k]))*raw[k];}
    mf::Point<double>p;const auto&origin=m.state->coefficient[t.id[0]];for(int a=0;a<t.n;a++)for(int k=0;k<2;k++){const auto&c=m.state->coefficient[t.id[a]];B x=B(c.x[k])-B(a<10?origin.x[k]:0),v=B(c.v[k])-B(a<10?origin.v[k]:0);p.position[k]=p.position[k]+x*N[a].v;p.velocity[k]=p.velocity[k]+v*N[a].v;p.F[2*k]=p.F[2*k]+x*N[a].x;p.F[2*k+1]=p.F[2*k+1]+x*N[a].y;p.V[2*k]=p.V[2*k]+v*N[a].x;p.V[2*k+1]=p.V[2*k+1]+v*N[a].y;}for(int k=0;k<2;k++){p.position[k]=p.position[k]+B(origin.x[k]);p.velocity[k]=p.velocity[k]+B(origin.v[k]);}mf::finish(p);return p;
}
}
#endif
