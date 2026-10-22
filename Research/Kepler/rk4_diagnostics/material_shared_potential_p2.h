#ifndef MATERIAL_SHARED_POTENTIAL_P2_H
#define MATERIAL_SHARED_POTENTIAL_P2_H
#include "material_shared_clip.h"

// A proof of the SAME Gauss16 potential quadrature error. It neither
// changes a stored quadrature value nor projects/refits a physical field.
namespace shared_target { namespace potential_p2 {
using J=sd::Jet;

static J power(const J&x,Q exponent){
    J out(holomorphic_power(x.v,exponent));
    S first=S(exponent)*holomorphic_power(x.v,exponent-Q(1));
    S second=S(exponent*(exponent-Q(1)))*holomorphic_power(x.v,exponent-Q(2));
    for(int i=0;i<2;i++){
        out.g[i]=first*x.g[i];
        for(int j=0;j<2;j++)out.H[i][j]=first*x.H[i][j]+second*x.g[i]*x.g[j];
    }
    return out;
}
static J root(const J&x){
    J out(complex_sqrt(x.v));
    S first=S(Q(.5L))*holomorphic_power(x.v,Q(-.5L));
    S second=S(Q(-.25L))*holomorphic_power(x.v,Q(-1.5L));
    for(int i=0;i<2;i++){
        out.g[i]=first*x.g[i];
        for(int j=0;j<2;j++)out.H[i][j]=first*x.H[i][j]+second*x.g[i]*x.g[j];
    }
    return out;
}
static J exponential(const J&x){
    J out(holomorphic_exp(x.v));
    for(int i=0;i<2;i++){
        out.g[i]=out.v*x.g[i];
        for(int j=0;j<2;j++)out.H[i][j]=out.v*(x.H[i][j]+x.g[i]*x.g[j]);
    }
    return out;
}
static J density(const curved::Packet&p,const J&radius){
    J inner=exponential(J(S(2))*(radius-J(sf::exact(p.parameters[4])))/J(sf::exact(p.parameters[8])));
    J outer=exponential(J(S(2))*(radius-J(sf::exact(p.parameters[5])))/J(sf::exact(p.parameters[8])));
    return J(sf::exact(p.parameters[6]))+J(sf::exact(p.parameters[9]))/(radius+J(sf::exact(p.parameters[7])))*inner/(J(S(1))+inner)/(J(S(1))+outer);
}

// Derivatives of the SAME certified implicit branch with respect to the
// normalized free coordinate u. No derivative of an interval radius is
// taken, and no root is selected from a different event slab.
static J graph(const Face&f,const sd::Chart&chart,const Axis&a,int index,S free,S du,bool complex,bool regular_q){
    S value,first,second;
    if(complex&&regular_q&&a.dependent==1){
        S center=shared_target::graph(f,chart,a,index,S(free.n));
        auto branch=graph_q_disk(f,chart,a,index,free,center);
        auto phi=face_q(f,chart,free,branch.coordinate);
        S q=branch.certificate.value;
        S dq=-phi.g[0]/phi.g[1];
        S ddq=-(phi.H[0][0]+S(2)*phi.H[0][1]*dq+phi.H[1][1]*dq*dq)/phi.g[1];
        value=q/free;first=(dq-value)/free;second=(ddq-S(2)*first)/free;
    }else{
        if(complex){
            S center=shared_target::graph(f,chart,a,index,S(free.n));
            value=graph_disk(f,chart,a,index,free,center).value;
        }else value=shared_target::graph(f,chart,a,index,free);
        auto phi=in_axis(f,chart,a.dependent,free,value);
        int d=a.dependent,u=1-d;
        first=-phi.g[u]/phi.g[d];
        second=-(phi.H[u][u]+S(2)*phi.H[u][d]*first+phi.H[d][d]*first*first)/phi.g[d];
    }
    J out(value,first*du);out.H[0][0]=second*du*du;return out;
}

struct Geometry {sd::Point point;J s,t,X[2],measure;};
static Geometry geometry(const sd::Compiled&m,const std::vector<Face>&faces,const Axis&a,const ClipPanel&panel,S u,S v,bool complex,bool regular_q){
    Geometry out;J U(u,S(1)),V(v,S(0),S(1));
    const auto&chart=panel.reference;
    if(panel.whole){
        S ds=panel.ordinary.shi-panel.ordinary.slo,dt=panel.ordinary.thi-panel.ordinary.tlo;
        out.s=J(panel.ordinary.slo)+J(ds)*U;
        out.t=J(panel.ordinary.tlo)+J(dt)*V;
        out.point=sd::evaluate(m,chart,out.s.v,out.t.v);
        out.measure=J(chart.det*ds*dt)*out.s;
    }else{
        S du=panel.strip.hi-panel.strip.lo,de=panel.eta_hi-panel.eta_lo;
        J free=J(panel.strip.lo)+J(du)*U,eta=J(panel.eta_lo)+J(de)*V;
        J lower(S(0)),upper(S(1));
        if(panel.strip.lower>=0)lower=graph(faces[panel.strip.lower],chart,a,panel.strip.lower,free.v,du,complex,regular_q);
        if(panel.strip.upper>=0)upper=graph(faces[panel.strip.upper],chart,a,panel.strip.upper,free.v,du,complex,regular_q);
        J width=upper-lower,dependent=lower+eta*width;
        out.s=a.dependent==1?free:dependent;out.t=a.dependent==1?dependent:free;
        auto original=strip_point(m,chart,faces,a,panel.strip,free.v,eta.v,complex,regular_q);
        out.point=original.value;
        out.measure=J(chart.det*du*de)*out.s*width;
        // The existing correlated construction encloses the same value;
        // retain its tighter value disk. Derivatives above enclose the
        // same geometric function on the identical branch.
        out.measure.v=original.jacobian*de;
    }
    auto r=frame(m,chart,S(0));
    for(int k=0;k<2;k++){
        out.X[k]=J(r.A[k])+out.s*(J(r.D[k])+out.t*J(r.E[k]));
        out.X[k].v=out.point.X[k];
    }
    return out;
}

static J integrand(const curved::Packet&p,const sd::Compiled&m,const std::vector<Face>&faces,const Axis&a,const ClipPanel&panel,S u,S v,bool complex,bool regular_q,const Thermal&context){
    auto g=geometry(m,faces,a,panel,u,v,complex,regular_q);
    J position[2];
    for(int k=0;k<2;k++){
        position[k].v=g.point.field.position[k];
        for(int i=0;i<2;i++){
            for(int x=0;x<2;x++)position[k].g[i]=position[k].g[i]+g.point.field.F[2*k+x]*g.X[x].g[i];
            for(int j=0;j<2;j++)for(int x=0;x<2;x++){
                position[k].H[i][j]=position[k].H[i][j]+g.point.field.F[2*k+x]*g.X[x].H[i][j];
                for(int y=0;y<2;y++)position[k].H[i][j]=position[k].H[i][j]+g.point.Hx[k][x][y]*g.X[x].g[i]*g.X[y].g[j];
            }
        }
    }
    bool origin=true;for(int k=0;k<2;k++)origin=origin&&m.base.vertex[panel.reference.apex][k].n==0&&m.base.vertex[panel.reference.apex][k].e==0;
    J radius;
    if(origin){
        auto r=frame(m,panel.reference,S(0));
        J D[2];for(int k=0;k<2;k++)D[k]=J(r.D[k])+g.t*J(r.E[k]);
        radius=g.s*root(D[0]*D[0]+D[1]*D[1]);
    }else radius=root(g.X[0]*g.X[0]+g.X[1]*g.X[1]);
    radius.v=reference_radius(m,panel.reference,g.s.v,g.t.v,g.point);
    J soft(sf::exact(p.parameters[3]));
    J potential=power(position[0]*position[0]+position[1]*position[1]+soft*soft,Q(-.5L));
    // rho0/J * J is cancelled analytically before ANY enclosures.
    J out=-density(p,radius)*potential*g.measure;
    auto original=reference_moments(p,g.point,radius.v,context);
    S expected=original[POTENTIAL]*g.measure.v;
    need(fabsq(out.v.n-expected.n)<=sf::add(out.v.e,expected.e),"normalized potential AD2 value agrees with unchanged reference integrand");
    return out;
}

static Q residual_bound(const J&center,const J&domain,const Q A[2][2],const Q*radius,Q*difference_out=nullptr){
    Q bound=center.v.e;
    for(int i=0;i<2;i++)bound=sf::add(bound,sf::mul(center.g[i].e,radius[i]));
    for(int i=0;i<2;i++)for(int j=0;j<2;j++){
        Q difference=sf::mag(domain.H[i][j]-S(A[i][j]));
        if(difference_out)difference_out[2*i+j]=difference;
        bound=sf::add(bound,sf::mul(Q(.5L),sf::mul(difference,sf::mul(radius[i],radius[j]))));
    }
    return bound;
}
static J polynomial_fixture(S u,S v){
    J U(u-S(Q(.5L)),S(1)),V(v-S(Q(.5L)),S(0),S(1));
    return J(S(1))+J(S(2))*U+J(S(3))*V+J(S(4))*U*U+J(S(5))*U*V+J(S(6))*V*V;
}
static bool unit_contract(){
    auto center=polynomial_fixture(S(Q(.5L)),S(Q(.5L)));
    Q A[2][2]={{center.H[0][0].n,center.H[0][1].n},{center.H[0][1].n,center.H[1][1].n}};
    for(int axis=0;axis<2;axis++){
        Q radius[2]={axis?Q(.5L):Q(1.25L),axis?Q(1.25L):Q(.5L)};
        auto domain=polynomial_fixture(S(Q(.5L),radius[0]),S(Q(.5L),radius[1]));
        need(residual_bound(center,domain,A,radius)<1e-25Q,"P2 subtraction removes exact normalized quadratic up to its forward arithmetic enclosure");
    }
    auto rule=gauss16();S sum;
    for(int i=0;i<16;i++)for(int j=0;j<16;j++)sum=sum+rule.w[i]*rule.w[j]*polynomial_fixture(rule.x[i],rule.x[j]).v;
    S exact=S(1)+S(10)/S(12);
    need(fabsq(sum.n-exact.n)<=sf::add(sum.e,exact.e),"finite Gauss16 weights/nodes enclose exact normalized P2 integral");
    return true;
}

struct Certificate {
    bool pass=false,center_complete=false,axis_complete[2]={false,false};uint64_t evaluations=0;
    Q remainder=0,coefficient[6]={},center_error[3]={},axis_majorant[2]={};
    Q hessian_difference[2][4]={};char reason[192]={};
};
static Certificate certify(const curved::Packet&p,const sd::Compiled&m,const std::vector<Face>&faces,const Axis&a,const ClipPanel&panel,bool regular_q,const Thermal&context){
    Certificate out;
    try{
        out.evaluations++;auto center=integrand(p,m,faces,a,panel,S(Q(.5L)),S(Q(.5L)),false,regular_q,context);
        out.coefficient[0]=center.v.n;out.coefficient[1]=center.g[0].n;out.coefficient[2]=center.g[1].n;
        out.coefficient[3]=center.H[0][0].n;out.coefficient[4]=center.H[0][1].n;out.coefficient[5]=center.H[1][1].n;
        out.center_error[0]=center.v.e;out.center_error[1]=center.g[0].e;out.center_error[2]=center.g[1].e;
        out.center_complete=true;
        Q A[2][2]={{out.coefficient[3],out.coefficient[4]},{out.coefficient[4],out.coefficient[5]}};
        for(int axis=0;axis<2;axis++){
            Q radius[2]={axis?Q(.5L):Q(1.25L),axis?Q(1.25L):Q(.5L)};
            S u(Q(.5L),radius[0]),v(Q(.5L),radius[1]);
            out.evaluations++;auto domain=integrand(p,m,faces,a,panel,u,v,true,regular_q,context);
            // The same product of disks is convex about c, so it encloses
            // c+t*(z-c) for the entire complex Taylor integral, 0<=t<=1.
            out.axis_majorant[axis]=residual_bound(center,domain,A,radius,out.hessian_difference[axis]);
            out.axis_complete[axis]=true;
        }
        out.remainder=sf::mul(sf::mag(cauchy_coefficient()),std::max(out.axis_majorant[0],out.axis_majorant[1]));
        need(finiteq(out.remainder)&&out.remainder>=0,"finite outward full-integrand P2 remainder");
        out.pass=true;snprintf(out.reason,sizeof(out.reason),"full_normalized_potential_P2_complex_Hessian_integral_remainder");
    }catch(const std::exception&e){snprintf(out.reason,sizeof(out.reason),"%s",e.what());}
    return out;
}
} }
#endif
