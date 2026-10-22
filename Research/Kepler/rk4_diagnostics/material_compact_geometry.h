#ifndef MATERIAL_COMPACT_GEOMETRY_H
#define MATERIAL_COMPACT_GEOMETRY_H
// Geometry-completed P3, per ASTRA_COMPACT_CURVED_DOMAIN_ADDENDUM_20261022.md.
// This is basis/kinematics/work algebra, NOT a domain-coverage certificate,
// positive curved quadrature, EOS, remap, or a committed representation.
#include "material_compact_p3.h"
#include "material_compact_arithmetic.h"
#include <cfloat>
#include <type_traits>
#ifdef __CUDACC__
#define COMPACT_GEOMETRY_HD __host__ __device__
#else
#define COMPACT_GEOMETRY_HD
#endif

namespace compact_geometry {
using compact_p3::R;
using compact_p3::require;
constexpr int local_size=12;
// Both fields are densities with respect to the STRAIGHT PARAMETER triangle
// measure dX. Initial det(DG) belongs in their projection, not in mass twice.
enum class DensityMeasure : uint32_t { parameter_reference_volume=1 };
constexpr DensityMeasure density_measure=DensityMeasure::parameter_reference_volume;
struct FrozenTrace {
    STEdge function={};
    STCoefficient coefficient[6]={};
    uint64_t edge=UINT64_MAX,source_hash=0,source_epoch=0;
};
struct ComponentModes {
    uint64_t component=UINT64_MAX,id[2]={UINT64_MAX,UINT64_MAX};
    double scale[2]={1,1};
    // phi_k = h_k/scale[k]; x-mode coefficients initially scale[k]*e_k.
    // Neither tiny scales nor singular modes are removed by a tolerance.
    // Whole-component rank/positive mass must be certified by the transaction.
};
struct Patch {
    compact_p3::Triangle polynomial;
    ComponentModes modes;
    FrozenTrace trace;
    double parameter[2]={0,1};
    double interpolation[10][2]={};
    uint32_t curved=0;
};
static_assert(std::is_trivially_copyable<Patch>::value,"immutable GPU geometry descriptor");
struct Basis {
    STJet N[local_size];
    uint64_t id[local_size]={};
    int n=10;
};
COMPACT_GEOMETRY_HD static inline bool finite_jet(const STJet&q){
    return isfinite(q.v)&&isfinite(q.x)&&isfinite(q.y);
}
// Evaluate the SAME accepted source trace, with derivatives by the existing
// jet algebra. No IC velocity, orbit, or new evolution is inserted as a field.
// The old analytic functions occur only in this frozen exterior geometry.
COMPACT_GEOMETRY_HD static inline int trace(const FrozenTrace&f,STJet u,STJet*x){
    if(f.function.rank<0||f.function.rank>2||!finite_jet(u)||u.v<0||u.v>1)return 0;
    STJet raw[2];if(!st_edge_raw(f.function,u,raw))return 0;
    for(int k=0;k<2;k++){
        double origin=f.coefficient[0].x[k];x[k]=STJet(origin);
        for(int a=0;a<4;a++)x[k]=x[k]+st_lagrange(a,u)*STJet(f.coefficient[a].x[k]-origin);
        for(int a=0;a<f.function.rank;a++){
            STJet N;for(int d=0;d<2;d++)N=N+STJet(f.function.encode[a][d])*raw[d];
            x[k]=x[k]+STJet(f.coefficient[4+a].x[k])*N;
        }
        if(!finite_jet(x[k]))return 0;
    }
    return 1;
}
// Evaluate Gamma(lo+t*(hi-lo))-((1-t)*B+t*C) WITHOUT forming two absolute
// positions and subtracting them. A polynomial trace is decomposed into its
// original endpoint line plus its two P3 curvature coefficients. Coefficient
// cancellation is in the original CUDA two-double arithmetic. Analytic scalar
// edge functions/jet evaluations remain the original functions, not a fit.
// Endpoint residuals are retained; no snapping, zero clamp or domain change.
COMPACT_GEOMETRY_HD static inline int relative_trace(const FrozenTrace&f,STJet t,
    double lo,double hi,const double*B,const double*C,STJet*out){
    using namespace compact_arithmetic;
    if(f.function.rank<0||f.function.rank>2||!finite_jet(t)||t.v<0||t.v>1||!(lo>=0&&lo<=1&&hi>=0&&hi<=1))return 0;
    STJet u=STJet(lo)+t*STJet(hi-lo),L1=st_lagrange(1,u),L2=st_lagrange(2,u),raw[2];
    if(f.function.rank&&!st_edge_raw(f.function,u,raw))return 0;
    Jet T(Pair(t.v),Pair(t.x),Pair(t.y)),J1(Pair(L1.v),Pair(L1.x),Pair(L1.y)),J2(Pair(L2.v),Pair(L2.x),Pair(L2.y));
    for(int k=0;k<2;k++){
        Pair a(f.coefficient[0].x[k]),d=Pair(f.coefficient[3].x[k])-a;
        Pair curve1=(Pair(f.coefficient[1].x[k])-a)-d/Pair(3);
        Pair curve2=(Pair(f.coefficient[2].x[k])-a)-Pair(2)*d/Pair(3);
        Pair origin=a+Pair(lo)*d-Pair(B[k]);
        Pair slope=(Pair(hi)-Pair(lo))*d-(Pair(C[k])-Pair(B[k]));
        Jet value=Jet(origin)+T*slope+J1*curve1+J2*curve2;
        for(int j=0;j<f.function.rank;j++){
            Jet N;for(int z=0;z<2;z++)N=N+Jet(Pair(raw[z].v),Pair(raw[z].x),Pair(raw[z].y))*Pair(f.function.encode[j][z]);
            value=value+N*Pair(f.coefficient[4+j].x[k]);
        }
        out[k]={rounded(value.v),rounded(value.x),rounded(value.y)};if(!finite_jet(out[k]))return 0;
    }
    return 1;
}
// Reference barycentrics l=(1-s,s(1-t),st). G=X+E, where
// E=s^4*(Gamma(t)-((1-t)*B+t*C)). A,B,C are actual canonical patch vertices.
// s=0 has zero value AND first derivative; no arbitrary epsilon apex branch.
COMPACT_GEOMETRY_HD static inline int extension(const Patch&p,const double*l,STJet*E){
    E[0]=E[1]=STJet();if(!p.curved)return 1;
    STJet L[3];for(int a=0;a<3;a++)L[a]={l[a],p.polynomial.grad[a][0],p.polynomial.grad[a][1]};
    STJet s=L[1]+L[2];if(s.v==0)return 1;
    STJet t=L[2]/s,relative[2];
    if(!relative_trace(p.trace,t,p.parameter[0],p.parameter[1],p.polynomial.vertex[1],p.polynomial.vertex[2],relative))return 0;
    for(int k=0;k<2;k++){
        E[k]=st_ipow(s,4)*relative[k];if(!finite_jet(E[k]))return 0;
    }
    return 1;
}
COMPACT_GEOMETRY_HD static inline int basis(const Patch&p,const double*l,Basis*out){
    *out={};compact_p3::shape(l,p.polynomial.grad,out->N);
    for(int a=0;a<10;a++)out->id[a]=p.polynomial.id[a];
    if(!p.curved)return 1;
    STJet E[2];if(!extension(p,l,E))return 0;
    out->n=local_size;
    for(int k=0;k<2;k++){
        if(!(p.modes.scale[k]>0)||!isfinite(p.modes.scale[k]))return 0;
        // G_P3 is the nodal P3 interpolant of G. The remainder is not
        // silently frozen outside the physical trial/test/velocity space.
        STJet h=E[k];for(int a=0;a<10;a++)h=h-out->N[a]*STJet(p.interpolation[a][k]);
        out->N[10+k]=h/STJet(p.modes.scale[k]);out->id[10+k]=p.modes.id[k];
        if(!finite_jet(out->N[10+k]))return 0;
    }
    return 1;
}
COMPACT_GEOMETRY_HD static inline int point(const Basis&b,const STCoefficient*c,GfsMaterialPoint12*out){
    *out={};double V[4]={};const auto&origin=c[b.id[0]];
    for(int a=0;a<b.n;a++)for(int k=0;k<2;k++){
        const auto&q=c[b.id[a]];double x=q.x[k]-(a<10?origin.x[k]:0),v=q.v[k]-(a<10?origin.v[k]:0);
        out->position[k]+=b.N[a].v*x;out->velocity[k]+=b.N[a].v*v;
        out->F[2*k]+=b.N[a].x*x;out->F[2*k+1]+=b.N[a].y*x;
        V[2*k]+=b.N[a].x*v;V[2*k+1]+=b.N[a].y*v;
    }
    for(int k=0;k<2;k++){out->position[k]+=origin.x[k];out->velocity[k]+=origin.v[k];}
    out->J=out->F[0]*out->F[3]-out->F[1]*out->F[2];if(!(out->J>0)||!isfinite(out->J))return 0;
    double A[4]={out->F[3]/out->J,-out->F[1]/out->J,-out->F[2]/out->J,out->F[0]/out->J};
    for(int k=0;k<2;k++){
        out->grad_v[2*k]=V[2*k]*A[0]+V[2*k+1]*A[2];
        out->grad_v[2*k+1]=V[2*k]*A[1]+V[2*k+1]*A[3];
        if(!isfinite(out->position[k])||!isfinite(out->velocity[k]))return 0;
    }
    for(int k=0;k<4;k++)if(!isfinite(out->F[k])||!isfinite(out->grad_v[k]))return 0;
    return 1;
}
struct StressValue {
    double force[2][local_size]={},Udot=0,gravity_power=0;
    double energy_residual=0,torque_residual=0,energy_scale=0,torque_scale=0;
    int valid=0;
};
// Local symmetric-stress work algebra only. The geometric-pressure volume
// derivative/nonlocal terms and the entropy law must ALSO use this basis in
// the complete EOS operator; this helper does not implement those terms.
COMPACT_GEOMETRY_HD static inline StressValue stress(const Basis&b,const STCoefficient*c,
    const GfsMaterialPoint12&p,double rho_reference,double weight,
    const double*sigma,double gx,double gy){
    StressValue out;double N[12]={},Nx[12]={},Ny[12]={};
    if(!isfinite(weight)||!isfinite(rho_reference)||!isfinite(gx)||!isfinite(gy))return out;
    for(int k=0;k<3;k++)if(!isfinite(sigma[k]))return out;
    for(int a=0;a<b.n;a++){N[a]=b.N[a].v;Nx[a]=b.N[a].x;Ny[a]=b.N[a].y;}
    if(!gfs_material12_stress_at_point(b.n,&p,N,Nx,Ny,rho_reference,weight,sigma,gx,gy,out.force,&out.Udot,&out.gravity_power))return out;
    double torque=weight*rho_reference*(p.position[0]*gy-p.position[1]*gx);
    out.energy_residual=out.Udot-out.gravity_power;out.torque_residual=-torque;
    out.energy_scale=fabs(out.Udot)+fabs(out.gravity_power);out.torque_scale=fabs(torque);
    for(int a=0;a<b.n;a++){
        const auto&q=c[b.id[a]];
        double work=q.v[0]*out.force[0][a]+q.v[1]*out.force[1][a];
        double moment=q.x[0]*out.force[1][a]-q.x[1]*out.force[0][a];
        out.energy_residual+=work;out.torque_residual+=moment;out.energy_scale+=fabs(work);out.torque_scale+=fabs(moment);
        if(!isfinite(work)||!isfinite(moment))return StressValue{};
    }
    out.valid=isfinite(out.energy_residual)&&isfinite(out.torque_residual);return out;
}
struct BoundaryValue {
    double force[2][local_size]={},work=0,torque=0,work_residual=0,torque_residual=0;
    int valid=0;
};
COMPACT_GEOMETRY_HD static inline BoundaryValue exterior(const Basis&b,const STCoefficient*c,
    const GfsMaterialPoint12&p,double nx,double ny,double weight,double pressure){
    BoundaryValue out;if(!(weight>=0)||!(pressure>0)||!isfinite(weight)||!isfinite(pressure)||!isfinite(nx)||!isfinite(ny))return out;
    double sigma[3]={-pressure,0,-pressure},traction[2],ratio;
    if(!gfs_material12_piola_traction(&p,sigma,nx,ny,traction,&ratio))return out;
    out.work=weight*(p.velocity[0]*traction[0]+p.velocity[1]*traction[1]);
    out.torque=weight*(p.position[0]*traction[1]-p.position[1]*traction[0]);
    out.work_residual=-out.work;out.torque_residual=-out.torque;
    for(int a=0;a<b.n;a++){
        const auto&q=c[b.id[a]];
        for(int k=0;k<2;k++)out.force[k][a]=weight*b.N[a].v*traction[k];
        out.work_residual+=q.v[0]*out.force[0][a]+q.v[1]*out.force[1][a];
        out.torque_residual+=q.x[0]*out.force[1][a]-q.x[1]*out.force[0][a];
    }
    out.valid=isfinite(out.work)&&isfinite(out.torque)&&isfinite(out.work_residual)&&isfinite(out.torque_residual);return out;
}
inline FrozenTrace freeze_trace(const STEdge&e,const std::vector<STCoefficient>&accepted,
    uint64_t edge,uint64_t hash,uint64_t epoch){
    require(e.rank>=0&&e.rank<=2,"complete accepted exterior trace rank");
    FrozenTrace out;out.function=e;out.edge=edge;out.source_hash=hash;out.source_epoch=epoch;
    for(int a=0;a<4+e.rank;a++){
        int id=a<4?e.poly[a]:e.id[a-4];require(id>=0&&size_t(id)<accepted.size(),"frozen exterior uses actual accepted coefficients");
        out.coefficient[a]=accepted[id];for(int k=0;k<2;k++)require(std::isfinite(out.coefficient[a].x[k])&&std::isfinite(out.coefficient[a].v[k]),"finite frozen exterior function");
    }
    for(int a=0;a<4;a++)out.function.poly[a]=a;
    for(int a=0;a<e.rank;a++)out.function.id[a]=4+a;
    return out;
}
inline void check_modes(const Patch&p,uint64_t scalar_extent){
    if(!p.curved)return;
    require(p.curved==1&&p.modes.component!=UINT64_MAX,"actual curved component identity");
    require(p.modes.id[0]<scalar_extent&&p.modes.id[1]<scalar_extent&&p.modes.id[0]!=p.modes.id[1],"two shared component mode owners");
    for(int k=0;k<2;k++){
        require(p.modes.scale[k]>0&&std::isfinite(p.modes.scale[k]),"positive finite mode normalization, no singular-value truncation");
        for(auto id:p.polynomial.id)require(id<scalar_extent&&id!=p.modes.id[k],"disjoint polynomial and geometry DOFs");
    }
}
inline void interpolate_chart(Patch&p){
    require(p.curved==1&&p.parameter[0]>=0&&p.parameter[0]<=1&&p.parameter[1]>=0&&p.parameter[1]<=1&&p.parameter[0]!=p.parameter[1],"oriented canonical exterior subinterval");
    STJet ends[2][2];
    for(int a=0;a<2;a++){
        require(trace(p.trace,STJet(p.parameter[a]),ends[a]),"actual accepted exterior endpoint");
        for(int k=0;k<2;k++)require(ends[a][k].v==p.polynomial.vertex[1+a][k],"same evaluated canonical curved endpoint, no tolerance merge");
    }
    for(int a=0;a<10;a++){
        int index[3];gfs_material_p3_index(a,index);double l[3];for(int k=0;k<3;k++)l[k]=index[k]/3.0;
        STJet E[2];require(extension(p,l,E),"finite canonical collar nodal interpolation");
        for(int k=0;k<2;k++)p.interpolation[a][k]=E[k].v;
    }
    // C0/event enclosure and J/coverage/rank are separate prerequisites. This
    // construction and equality of rounded endpoints do NOT prove them.
}
inline STCoefficient initial_mode(const ComponentModes&m,int k){
    require(k==0||k==1,"geometry mode component");STCoefficient out={};out.x[k]=m.scale[k];return out;
}
inline std::array<double,2> initial_polynomial_position(const Patch&p,int a){
    require(a>=0&&a<10,"compact polynomial node");int index[3];gfs_material_p3_index(a,index);
    std::array<double,2>out={};for(int k=0;k<2;k++){
        R x=0;for(int j=0;j<3;j++)x+=R(index[j])*p.polynomial.vertex[j][k]/3;
        out[k]=double(x)+p.interpolation[a][k];
    }return out;
}
// A LOCAL block of the supplied positive functional. For curved patches the
// caller must supply its actual nonpolynomial quadrature with error accounting;
// the straight degree9 rule is not a curved-collar accuracy certificate.
inline std::array<R,144> mass(const Patch&p,const compact_p3::ThermalDensity&f,
    const std::vector<compact_p3::QuadraturePoint>&quadrature){
    require(compact_p3::nonnegative_controls(f.mass),"positive reference-measure density candidate");
    require(!quadrature.empty(),"explicit positive integration functional");std::array<R,144>M={};
    for(const auto&q:quadrature){
        double total=0;for(double l:q.l){require(l>=0&&std::isfinite(l),"inside mass quadrature point");total+=l;}
        require(fabs(total-1)<=32*DBL_EPSILON&&q.weight>0&&std::isfinite(q.weight),"positive parameter-measure quadrature");
        Basis b;require(basis(p,q.l,&b),"same completed mass basis");
        double B[10];compact_p3::bernstein(q.l,B);R rho=0;for(int a=0;a<10;a++)rho+=R(f.mass[a])*B[a];require(rho>0,"positive represented reference density");
        for(int a=0;a<b.n;a++)for(int j=0;j<=a;j++)M[a*12+j]+=R(q.weight)*rho*R(b.N[a].v)*b.N[j].v;
    }
    for(int a=0;a<12;a++)for(int j=0;j<a;j++)M[j*12+a]=M[a*12+j];return M;
}
} // namespace compact_geometry
#undef COMPACT_GEOMETRY_HD
#endif
