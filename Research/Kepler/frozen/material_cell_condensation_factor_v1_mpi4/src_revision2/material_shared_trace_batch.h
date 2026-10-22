#ifndef MATERIAL_SHARED_TRACE_BATCH_H
#define MATERIAL_SHARED_TRACE_BATCH_H
#include "gfs_material_p3.h"
#include "gfs_material_p3_radial_input.h"
#include "gfs_material_map12.h"
#include <stddef.h>
#ifdef __CUDACC__
#define ST_HD __host__ __device__
#else
#define ST_HD
#endif

// Fixed reference functions; no force/acceleration/future orbit is a basis.
// Cache entries approximate analytic IC evaluations, not kinematic state.
// All dynamic coefficients live ONCE in the shared global coefficient table.
enum {ST_MAX=18};
struct STEdge {
    double A[2],D[2],qnode[4][2],eps;
    double encode[2][2],decode[2][2];
    int rank,id[2],poly[4];
};
struct STTriangle {
    double vertex[3][2],grad[3][2],qnode[10][2],eps;
    double bubble_shift[2],bubble_encode[2][2],bubble_decode[2][2];
    STEdge edge[3];int endpoint[3][2],id[ST_MAX],offset[3],bubble_offset,n,cell;
};
struct STCoefficient {double x[2],v[2];};
struct STPoint {double X,Y,l[3],weight,rho;int triangle,boundary;};
struct STOutput {
    GfsMaterialPoint12 point;
    double force[2][ST_MAX],U,Udot,gpower,gtorque,pressure;
    double identity_energy,identity_angular,scale_energy,scale_angular;
    int valid;
};
struct STPhysics {double gamma,P0,eps;int gravity;};
struct STJet {double v,x,y;ST_HD STJet(double a=0,double b=0,double c=0):v(a),x(b),y(c){}};
ST_HD static inline STJet operator+(STJet a,STJet b){return {a.v+b.v,a.x+b.x,a.y+b.y};}
ST_HD static inline STJet operator-(STJet a){return {-a.v,-a.x,-a.y};}
ST_HD static inline STJet operator-(STJet a,STJet b){return a+(-b);}
ST_HD static inline STJet operator*(STJet a,STJet b){return {a.v*b.v,a.x*b.v+a.v*b.x,a.y*b.v+a.v*b.y};}
ST_HD static inline STJet operator/(STJet a,STJet b){double v=a.v/b.v;return {v,(a.x-v*b.x)/b.v,(a.y-v*b.y)/b.v};}
ST_HD static inline STJet st_ipow(STJet a,int n){STJet b(1);for(int i=0;i<n;i++)b=b*a;return b;}
ST_HD static inline STJet st_lagrange(int node,STJet t){STJet p(1);for(int j=0;j<4;j++)if(j!=node)p=p*(STJet(3)*t-STJet(j))/STJet(node-j);return p;}
ST_HD static inline int st_initial_remainder(const double*anchor,STJet X,STJet Y,double eps,STJet*q){
    double value[2],gradient[4];if(!gfs_material_p3_radial_remainder(anchor,X.v,Y.v,eps,value,gradient))return 0;
    for(int k=0;k<2;k++)q[k]={value[k],gradient[2*k]*X.x+gradient[2*k+1]*Y.x,gradient[2*k]*X.y+gradient[2*k+1]*Y.y};return 1;
}
ST_HD static inline int st_edge_raw(const STEdge&e,STJet t,STJet*r){
    STJet X=STJet(e.A[0])+t*STJet(e.D[0]),Y=STJet(e.A[1])+t*STJet(e.D[1]);if(!st_initial_remainder(e.A,X,Y,e.eps,r))return 0;
    for(int j=0;j<4;j++){STJet L=st_lagrange(j,t);for(int k=0;k<2;k++)r[k]=r[k]-L*STJet(e.qnode[j][k]);}
    // Exact endpoint identity of v0-I3(v0). This is evaluation of a known
    // zero trace; no state coefficient or physical energy is changed.
    if(t.v==0||t.v==1)for(int k=0;k<2;k++)r[k].v=0;
    return 1;
}
ST_HD static inline int st_shape(const STTriangle&t,const STPoint&p,STJet*N){
    STJet l[3];for(int j=0;j<3;j++)l[j]={p.l[j],t.grad[j][0],t.grad[j][1]};
    for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);N[a]=STJet(1);for(int k=0;k<3;k++)for(int j=0;j<index[k];j++)N[a]=N[a]*(STJet(3)*l[k]-STJet(j))/STJet(index[k]-j);}
    STJet raw[2];if(!st_initial_remainder(t.vertex[0],STJet(p.X,1,0),STJet(p.Y,0,1),t.eps,raw))return 0;
    for(int a=0;a<10;a++)for(int k=0;k<2;k++)raw[k]=raw[k]-N[a]*STJet(t.qnode[a][k]);
    for(int e=0;e<3;e++){
        STJet s=l[t.endpoint[e][0]]+l[t.endpoint[e][1]],extension[2];
        if(s.v!=0){STJet r[2];if(!st_edge_raw(t.edge[e],l[t.endpoint[e][1]]/s,r))return 0;for(int k=0;k<2;k++)extension[k]=st_ipow(s,4)*r[k];}
        for(int k=0;k<2;k++)raw[k]=raw[k]-extension[k];
        for(int k=0;k<t.edge[e].rank;k++){N[t.offset[e]+k]=STJet();for(int c=0;c<2;c++)N[t.offset[e]+k]=N[t.offset[e]+k]+STJet(t.edge[e].encode[k][c])*extension[c];}
    }
    // Analytic b_T=v0-I3(v0)-sum E_e(r_e) is zero on EVERY edge.
    // Boundary topology is explicit, rather than inferred from a tolerance.
    if(p.boundary>=0)for(int k=0;k<2;k++)raw[k].v=0;
    for(int k=0;k<2;k++)raw[k]=raw[k]-STJet(t.bubble_shift[k])*N[9];
    for(int k=0;k<2;k++){N[t.bubble_offset+k]=STJet();for(int c=0;c<2;c++)N[t.bubble_offset+k]=N[t.bubble_offset+k]+STJet(t.bubble_encode[k][c])*raw[c];}
    return 1;
}
ST_HD static inline int st_point(const STTriangle&t,const STCoefficient*c,const STJet*N,GfsMaterialPoint12*p){
    *p={};double V[4]={0};const auto&origin=c[t.id[0]];for(int a=0;a<t.n;a++)for(int k=0;k<2;k++){
        const auto&q=c[t.id[a]];double x=q.x[k]-(a<10?origin.x[k]:0),v=q.v[k]-(a<10?origin.v[k]:0);
        p->position[k]+=x*N[a].v;p->velocity[k]+=v*N[a].v;
        p->F[2*k]+=x*N[a].x;p->F[2*k+1]+=x*N[a].y;V[2*k]+=v*N[a].x;V[2*k+1]+=v*N[a].y;
    }
    for(int k=0;k<2;k++){p->position[k]+=origin.x[k];p->velocity[k]+=origin.v[k];}
    p->J=p->F[0]*p->F[3]-p->F[1]*p->F[2];if(!(p->J>0)||!isfinite(p->J))return 0;
    double A[4]={p->F[3]/p->J,-p->F[1]/p->J,-p->F[2]/p->J,p->F[0]/p->J};for(int k=0;k<2;k++){p->grad_v[2*k]=V[2*k]*A[0]+V[2*k+1]*A[2];p->grad_v[2*k+1]=V[2*k]*A[1]+V[2*k+1]*A[3];}return 1;
}
ST_HD static inline int st_trace(const STEdge&e,const STCoefficient*c,double t,double*x,double*v){
    STJet raw[2];if(!st_edge_raw(e,STJet(t),raw))return 0;const auto&origin=c[e.poly[0]];
    for(int k=0;k<2;k++){x[k]=v[k]=0;for(int j=0;j<4;j++){double N=st_lagrange(j,STJet(t)).v;x[k]+=(c[e.poly[j]].x[k]-origin.x[k])*N;v[k]+=(c[e.poly[j]].v[k]-origin.v[k])*N;}for(int j=0;j<e.rank;j++){double N=0;for(int a=0;a<2;a++)N+=e.encode[j][a]*raw[a].v;x[k]+=c[e.id[j]].x[k]*N;v[k]+=c[e.id[j]].v[k]*N;}x[k]+=origin.x[k];v[k]+=origin.v[k];}return 1;
}
ST_HD static inline STOutput st_evaluate(const STTriangle*t,const STCoefficient*c,const STPoint&q,const STPhysics&physics){
    STOutput out={};const auto&triangle=t[q.triangle];STJet N[ST_MAX];if(!st_shape(triangle,q,N)||!st_point(triangle,c,N,&out.point))return out;
    if(q.boundary>=0&&!st_trace(triangle.edge[q.boundary],c,q.l[triangle.endpoint[q.boundary][1]],out.point.position,out.point.velocity))return out;
    out.pressure=physics.P0*pow(out.point.J,-physics.gamma);out.U=q.weight*physics.P0/(physics.gamma-1)*pow(out.point.J,1-physics.gamma);if(!(out.pressure>0)||!(out.U>=0))return out;
    double x=out.point.position[0],y=out.point.position[1],gx=0,gy=0;if(physics.gravity){double g=pow(x*x+y*y+physics.eps*physics.eps,-1.5);gx=-x*g;gy=-y*g;}
    double sigma[3]={-out.pressure,0,-out.pressure};
    // Reuse the user's native symmetric stress-work operator in blocks of
    // at most12 basis functions. Udot/gpower count ONCE, not per block.
    for(int begin=0;begin<triangle.n;begin+=12){int n=triangle.n-begin;if(n>12)n=12;double A[12],Ax[12],Ay[12],f[2][12]={},u=0,g=0;
        for(int j=0;j<n;j++){A[j]=N[begin+j].v;Ax[j]=N[begin+j].x;Ay[j]=N[begin+j].y;}
        if(!gfs_material12_stress_at_point(n,&out.point,A,Ax,Ay,q.rho,q.weight,sigma,gx,gy,f,&u,&g))return out;
        if(begin==0){out.Udot=u;out.gpower=g;}for(int j=0;j<n;j++)for(int k=0;k<2;k++)out.force[k][begin+j]=f[k][j];
    }
    out.gtorque=q.weight*q.rho*(x*gy-y*gx);out.identity_energy=out.Udot-out.gpower;out.identity_angular=-out.gtorque;out.scale_energy=fabs(out.Udot)+fabs(out.gpower);out.scale_angular=fabs(out.gtorque);
    for(int a=0;a<triangle.n;a++){auto z=c[triangle.id[a]];double w=z.v[0]*out.force[0][a]+z.v[1]*out.force[1][a],T=z.x[0]*out.force[1][a]-z.x[1]*out.force[0][a];out.identity_energy+=w;out.identity_angular+=T;out.scale_energy+=fabs(w);out.scale_angular+=fabs(T);}out.valid=1;return out;
}
extern "C" int st_cuda_evaluate(const STTriangle*,size_t,const STCoefficient*,size_t,const STPoint*,size_t,const STPhysics*,STOutput*,double*);
#undef ST_HD
#endif
