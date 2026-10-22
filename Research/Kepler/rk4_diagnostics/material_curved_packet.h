#ifndef MATERIAL_CURVED_PACKET_H
#define MATERIAL_CURVED_PACKET_H
#include <array>
#include "gfs_material_p3.h"
#include "gfs_material_p3_radial_input.h"
#include "gfs_material_map12.h"
namespace curved {
using R=long double;
struct Record {uint64_t cell;int component;GfsMaterialP3Triangle triangle;GfsMaterialMap12 map;double origin[2],scale[2];};
struct Packet {double parameters[12];uint64_t cells,dofs;std::vector<Record> record;};
static Packet read(const char *path){
    FILE *f=fopen(path,"rb");need(f,"native curved map control packet");uint64_t h[8];need(fread(h,sizeof(h),1,f)==1 && h[0]==0x4746503346494e31ULL && h[1]==1 && h[5]==12 && h[6]==10 && h[7]==864,"curved map packet format");Packet out={};out.cells=h[3];out.dofs=h[4];need(fread(out.parameters,sizeof(out.parameters),1,f)==1,"curved control physical parameters");out.record.resize(h[2]);
    for(auto &r:out.record){
        uint64_t tag[3];double geometry[13];need(fread(tag,sizeof(tag),1,f)==1 && tag[2]==12 && fread(geometry,sizeof(geometry),1,f)==1,"curved control record");r.cell=tag[0];r.component=tag[1];r.map.n=12;
        for(int a=0;a<3;a++)for(int k=0;k<2;k++){r.triangle.geometry.x[a][k]=geometry[2*a+k];r.triangle.geometry.grad[a][k]=geometry[6+2*a+k];}r.triangle.geometry.area=geometry[12];
        need(fread(r.triangle.v0,sizeof(r.triangle.v0),1,f)==1 && fread(r.triangle.remainder,sizeof(r.triangle.remainder),1,f)==1 && fread(r.map.x,sizeof(r.map.x),1,f)==1 && fread(r.map.v,sizeof(r.map.v),1,f)==1 && fread(r.origin,sizeof(r.origin),1,f)==1 && fread(r.scale,sizeof(r.scale),1,f)==1,"curved control fixed basis and actual state");
    }need(fgetc(f)==EOF,"curved map packet exact end");fclose(f);return out;
}
static GfsMaterialPoint12 native_point(const Record &r,double X,double Y,double eps){
    double q[2],Dq[4],N[12],Nx[12],Ny[12];need(gfs_material_p3_radial_remainder(r.triangle.geometry.x[0],X,Y,eps,q,Dq),"curved initial basis domain");gfs_material_p3_shape_remainder(&r.triangle,X,Y,q,Dq,r.scale,N,Nx,Ny);GfsMaterialPoint12 p;need(gfs_material12_point(&r.map,N,Nx,Ny,&p),"curved actual positive map");p.position[0]+=r.origin[0];p.position[1]+=r.origin[1];return p;
}
struct Dual {
    R v,d;Dual(R value=0,R derivative=0):v(value),d(derivative){}
};
static Dual operator+(Dual a,Dual b){return {a.v+b.v,a.d+b.d};}
static Dual operator-(Dual a,Dual b){return {a.v-b.v,a.d-b.d};}
static Dual operator*(Dual a,Dual b){return {a.v*b.v,a.d*b.v+a.v*b.d};}
struct Poly {
    int degree;std::array<Dual,7> c;explicit Poly(int d=0):degree(d),c{}{}
};
static Poly constant(R a){Poly p;p.c[0]=Dual(a);return p;}
static Poly add(const Poly &a,const Poly &b){Poly p(std::max(a.degree,b.degree));for(int i=0;i<=p.degree;i++)p.c[i]=a.c[i]+b.c[i];return p;}
static Poly scale(const Poly &a,R b){Poly p(a.degree);for(int i=0;i<=a.degree;i++)p.c[i]=a.c[i]*Dual(b);return p;}
static Poly scale_dual(const Poly &a,Dual b){Poly p(a.degree);for(int i=0;i<=a.degree;i++)p.c[i]=a.c[i]*b;return p;}
static Poly sub(const Poly &a,const Poly &b){return add(a,scale(b,-1));}
static Poly mul(const Poly &a,const Poly &b){need(a.degree+b.degree<=6,"exact ray polynomial degree cap, no truncation");Poly p(a.degree+b.degree);for(int i=0;i<=a.degree;i++)for(int j=0;j<=b.degree;j++)p.c[i+j]=p.c[i+j]+a.c[i]*b.c[j];return p;}
static Dual value(const Poly &p,R u){Dual v=p.c[p.degree];for(int i=p.degree-1;i>=0;i--)v=v*Dual(u)+p.c[i];return v;}
static R derivative(const Poly &p,R u){R v=0;for(int i=p.degree;i>=1;i--)v=v*u+i*p.c[i].v;return v;}
struct Ray {R angle;Poly denominator;Poly position[2],velocity[2];};
/* At FIXED reference radius, ANY actual state in this P3+initial-residual
 * space has cubic trigonometric x/v, not an exact prescribed Kepler orbit.
 * u=tan((theta-angle)/2) makes every power-plane preimage degree<=6 over
 * the positive denominator (1+u^2)^3. Dual coefficients give exact dr.
 */
static Ray ray(const Record &r,R radius,R angle,R eps){
    Ray out;out.angle=angle;Poly D(2);D.c[0]=Dual(1);D.c[2]=Dual(1);out.denominator=D;Poly D2=mul(D,D),D3=mul(D2,D);Dual rad(radius,1);R cs=cosl(angle),sn=sinl(angle),ax=r.triangle.geometry.x[0][0],ay=r.triangle.geometry.x[0][1],za=ax*ax+ay*ay+eps*eps;
    Poly X(2),Y(2);X.c[0]=rad*Dual(cs);X.c[1]=rad*Dual(-2*sn);X.c[2]=rad*Dual(-cs);Y.c[0]=rad*Dual(sn);Y.c[1]=rad*Dual(2*cs);Y.c[2]=rad*Dual(-sn);
    Poly dx=sub(X,scale(D,ax)),dy=sub(Y,scale(D,ay));Poly l1=add(scale(dx,r.triangle.geometry.grad[1][0]),scale(dy,r.triangle.geometry.grad[1][1])),l2=add(scale(dx,r.triangle.geometry.grad[2][0]),scale(dy,r.triangle.geometry.grad[2][1])),l0=sub(sub(D,l1),l2);Poly bary[3]={l0,l1,l2},shape[10];
    for(int a=0;a<10;a++){shape[a]=constant(1);int index[3];gfs_material_p3_index(a,index);for(int k=0;k<3;k++)for(int j=0;j<index[k];j++)shape[a]=mul(shape[a],scale(sub(scale(bary[k],3),scale(D,j)),1.L/(index[k]-j)));}
    Poly L=scale(add(scale(dx,2*ax),scale(dy,2*ay)),1/za),Q=scale(add(mul(dx,dx),mul(dy,dy)),1/za),L2=mul(L,L),L3=mul(L2,L);
    const R c1=-.75L,c2=.65625L,c3=-.6015625L;
    Poly A3=add(add(add(D3,scale(mul(L,D2),c1)),scale(mul(Q,D),c1)),add(add(scale(mul(L2,D),c2),scale(mul(L,Q),2*c2)),scale(L3,c3)));
    Poly A2=add(add(D2,scale(mul(L,D),c1)),add(scale(Q,c1),scale(L2,c2)));
    Poly T3[2]={scale(add(scale(A3,-ay),mul(scale(dy,-1),A2)),powl(za,-.75L)),scale(add(scale(A3,ax),mul(dx,A2)),powl(za,-.75L))};
    R z=radius*radius+eps*eps,omega=powl(z,-.75L);Dual o(omega,-1.5L*radius*omega/z);Poly v0[2]={scale_dual(mul(scale(Y,-1),D2),o),scale_dual(mul(X,D2),o)},residual[2]={sub(v0[0],T3[0]),sub(v0[1],T3[1])};
    for(int k=0;k<2;k++)for(int a=0;a<10;a++)residual[k]=sub(residual[k],scale(shape[a],r.triangle.remainder[a][k]));
    for(int k=0;k<2;k++){
        out.position[k]=scale(D3,r.origin[k]);out.velocity[k]=constant(0);
        for(int a=0;a<10;a++){out.position[k]=add(out.position[k],scale(shape[a],r.map.x[k][a]));out.velocity[k]=add(out.velocity[k],scale(shape[a],r.map.v[k][a]));}
        for(int j=0;j<2;j++){out.position[k]=add(out.position[k],scale(residual[j],r.map.x[k][10+j]/(R)r.scale[j]));out.velocity[k]=add(out.velocity[k],scale(residual[j],r.map.v[k][10+j]/(R)r.scale[j]));}
    }return out;
}
static std::array<R,3> evaluate(const Poly &p,R u){
    R D=1+u*u,D2=D*D,D3=D2*D;Dual v=value(p,u);return {v.v/D3,v.d/D3,derivative(p,u)/(2*D2)-3*u*v.v/D3};
}
static Poly power_plane(const Ray &ray,R ax,R ay,R b){return sub(add(scale(ray.position[0],ax),scale(ray.position[1],ay)),scale(mul(mul(ray.denominator,ray.denominator),ray.denominator),b));}
}
#endif
