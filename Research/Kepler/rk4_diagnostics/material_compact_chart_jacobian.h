#ifndef MATERIAL_COMPACT_CHART_JACOBIAN_H
#define MATERIAL_COMPACT_CHART_JACOBIAN_H
// Continuous whole-reference Jacobian bounds for the selected curved-collar
// chart G(s,t)=A+s[(1-t)B+tC]+s^4[Gamma(t)-((1-t)B+tC)].
#include "material_compact_chart_certificate.h"

namespace compact_chart {
namespace sf=shared_function;
using Q=sf::Q;using D=sf::D;
struct JacobianBudget {size_t maximum_boxes_per_patch=4096;unsigned maximum_depth=32;};
struct JacobianCertificate {
    uint64_t source_registry_hash=0,weight_hash=0,source_epoch=0,source_steps=0,chart_hash=0;
    double source_time=0;
    std::vector<Q>patch_lower;
    size_t boxes=0,subdivisions=0,curved_patches=0;
    Q minimum_lower=0;
    bool all_reference_J_positive=false;
};
inline S cross2(const std::array<S,2>&a,const std::array<S,2>&b){return a[0]*b[1]-a[1]*b[0];}
inline bool intersects(const S&a,const S&b){return std::max(sf::lower(a),sf::lower(b))<=std::min(sf::upper(a),sf::upper(b));}

inline std::array<D,2> curve_position(const cp::Curve&curve,D parameter){
    const auto&c=curve;std::array<D,2>raw;
    if(c.scalar.native.rank)raw=sf::edge_raw(c.scalar,parameter);
    std::array<D,2>out;
    for(int k=0;k<2;k++){
        S origin=sf::exact(c.coefficient[0].x[k]),delta=sf::exact(c.coefficient[3].x[k])-origin;
        S bend1=(sf::exact(c.coefficient[1].x[k])-origin)-delta/S(3);
        S bend2=(sf::exact(c.coefficient[2].x[k])-origin)-S(2)*delta/S(3);
        D x=D(origin)+parameter*D(delta)+sf::lagrange(1,parameter)*D(bend1)+sf::lagrange(2,parameter)*D(bend2);
        for(int j=0;j<c.scalar.native.rank;j++){
            D n;for(int z=0;z<2;z++)n=n+D(sf::exact(c.scalar.native.encode[j][z]))*raw[z];
            x=x+D(sf::exact(c.coefficient[4+j].x[k]))*n;
        }
        out[k]=x;
    }
    return out;
}

// Evaluate Gamma(u0+t*(u1-u0)) - (B+t*(C-B)) before forming either large
// absolute position. The same small canonical correction is used in the
// chart map, its Jacobian, and the Green-boundary area identity.
inline std::array<D,2> curve_relative(const Patch&patch,const cp::Curve&curve,D t){
    const auto&edge=patch.edge[1];S u0=edge.parameter[0],du=edge.parameter[1]-u0;
    D parameter=D(u0)+t*D(du);std::array<D,2>raw={},raw0={},raw1={};
    if(curve.scalar.native.rank){
        raw=sf::edge_raw(curve.scalar,parameter);
        raw0=sf::edge_raw(curve.scalar,D(u0));
        raw1=sf::edge_raw(curve.scalar,D(u0+du));
    }
    auto gamma0=curve_position(curve,D(u0)),gamma1=curve_position(curve,D(u0+du));
    D factor=D(du*du)*t*(t-D(S(1)));
    D residual1=factor*(D(-S(45)/S(2)+S(81)*u0/S(2))+D(S(27)*du/S(2))*(t+D(S(1))));
    D residual2=factor*(D(S(18)-S(81)*u0/S(2))-D(S(27)*du/S(2))*(t+D(S(1))));
    std::array<D,2>out;
    for(int k=0;k<2;k++){
        S origin=sf::exact(curve.coefficient[0].x[k]),delta=sf::exact(curve.coefficient[3].x[k])-origin;
        S bend1=(sf::exact(curve.coefficient[1].x[k])-origin)-delta/S(3);
        S bend2=(sf::exact(curve.coefficient[2].x[k])-origin)-S(2)*delta/S(3);
        S end0=gamma0[k].v-patch.vertex[1].x[k],end1=gamma1[k].v-patch.vertex[2].x[k];
        D value=(D(S(1))-t)*D(end0)+t*D(end1)+residual1*D(bend1)+residual2*D(bend2);
        for(int j=0;j<curve.scalar.native.rank;j++){
            D current;S first,last;
            for(int z=0;z<2;z++){
                S coefficient=sf::exact(curve.scalar.native.encode[j][z]);
                current=current+D(coefficient)*raw[z];
                first=first+coefficient*raw0[z].v;last=last+coefficient*raw1[z].v;
            }
            value=value+D(sf::exact(curve.coefficient[4+j].x[k]))*(current-(D(S(1))-t)*D(first)-t*D(last));
        }
        out[k]=value;
    }
    return out;
}

// Second-order interval jet in the collar parameter. It provides a rigorous
// centered form for thin curved ears, where direct interval evaluation of
// Gamma(t)-chord(t) loses correlation despite a positive pointwise Jacobian.
struct J2 {S v,d,dd;J2(S value=S(),S first=S(),S second=S()):v(value),d(first),dd(second){}};
inline J2 operator+(J2 a,J2 b){return {a.v+b.v,a.d+b.d,a.dd+b.dd};}
inline J2 operator-(J2 a){return {-a.v,-a.d,-a.dd};}
inline J2 operator-(J2 a,J2 b){return a+(-b);}
inline J2 operator*(J2 a,J2 b){return {a.v*b.v,a.d*b.v+a.v*b.d,a.dd*b.v+S(2)*a.d*b.d+a.v*b.dd};}
inline J2 operator/(J2 a,S b){return {a.v/b,a.d/b,a.dd/b};}
inline J2 jet_power(J2 a,Q exponent){
    S value=sf::mf::power(a.v,exponent),first=S(exponent)*value/a.v;
    S second=S(exponent*(exponent-Q(1)))*value/(a.v*a.v);
    return {value,first*a.d,first*a.dd+second*a.d*a.d};
}
inline J2 jet_lagrange(int j,J2 t){
    J2 out(S(1));for(int k=0;k<4;k++)if(j!=k)
        out=out*((J2(S(3))*t-J2(S(k)))/S(j-k));
    return out;
}
inline std::array<J2,2> jet_edge_raw(const cp::Curve&curve,J2 u){
    const auto&e=curve.scalar;
    J2 X=J2(e.A[0])+u*J2(e.Dref[0]);
    J2 Y=J2(e.A[1])+u*J2(e.Dref[1]);
    S eps=sf::exact(e.native.eps);
    J2 w=jet_power(X*X+Y*Y+J2(eps*eps),Q(-.75L));
    std::array<J2,2>out={-Y*w,X*w};
    for(int j=0;j<4;j++){
        J2 L=jet_lagrange(j,u);
        for(int k=0;k<2;k++)out[k]=out[k]-L*J2(e.node[j][k]);
    }
    return out;
}
inline std::array<J2,2> curve_relative_second(
        const Patch&patch,const cp::Curve&curve,S parameter){
    const auto&edge=patch.edge[1];S u0=edge.parameter[0],du=edge.parameter[1]-u0;
    J2 t(parameter,S(1));J2 u=J2(u0)+t*J2(du);
    std::array<J2,2>raw={};std::array<D,2>raw0={},raw1={};
    if(curve.scalar.native.rank){
        raw=jet_edge_raw(curve,u);
        raw0=sf::edge_raw(curve.scalar,D(u0));
        raw1=sf::edge_raw(curve.scalar,D(u0+du));
    }
    auto gamma0=curve_position(curve,D(u0));
    auto gamma1=curve_position(curve,D(u0+du));
    J2 factor=J2(du*du)*t*(t-J2(S(1)));
    J2 residual1=factor*(J2(-S(45)/S(2)+S(81)*u0/S(2))+
        J2(S(27)*du/S(2))*(t+J2(S(1))));
    J2 residual2=factor*(J2(S(18)-S(81)*u0/S(2))-
        J2(S(27)*du/S(2))*(t+J2(S(1))));
    std::array<J2,2>out;
    for(int k=0;k<2;k++){
        S origin=sf::exact(curve.coefficient[0].x[k]);
        S delta=sf::exact(curve.coefficient[3].x[k])-origin;
        S bend1=(sf::exact(curve.coefficient[1].x[k])-origin)-delta/S(3);
        S bend2=(sf::exact(curve.coefficient[2].x[k])-origin)-S(2)*delta/S(3);
        S end0=gamma0[k].v-patch.vertex[1].x[k];
        S end1=gamma1[k].v-patch.vertex[2].x[k];
        J2 value=(J2(S(1))-t)*J2(end0)+t*J2(end1)+
            residual1*J2(bend1)+residual2*J2(bend2);
        for(int j=0;j<curve.scalar.native.rank;j++){
            J2 current;S first,last;
            for(int z=0;z<2;z++){
                S coefficient=sf::exact(curve.scalar.native.encode[j][z]);
                current=current+J2(coefficient)*raw[z];
                first=first+coefficient*raw0[z].v;
                last=last+coefficient*raw1[z].v;
            }
            value=value+J2(sf::exact(curve.coefficient[4+j].x[k]))*
                (current-(J2(S(1))-t)*J2(first)-t*J2(last));
        }
        out[k]=value;
    }
    return out;
}

struct RelativeBounds {std::array<S,2> position,derivative;};
inline RelativeBounds curve_relative_centered_bounds(
        const Patch&patch,const cp::Curve&curve,Q t0,Q t1){
    Q tm=(t0+t1)/2;S T=sf::range(t0,t1),h=T-S(tm);
    auto full=curve_relative_second(patch,curve,T);
    auto middle=curve_relative_second(patch,curve,S(tm));
    auto direct=curve_relative(patch,curve,D(T,S(1),S(0)));
    RelativeBounds out;
    for(int k=0;k<2;k++){
        S centered=middle[k].v+middle[k].d*h+full[k].dd*h*h/S(2);
        S derivative=middle[k].d+full[k].dd*h;
        Q lo=std::max(sf::lower(centered),sf::lower(direct[k].v));
        Q hi=std::min(sf::upper(centered),sf::upper(direct[k].v));
        require(lo<=hi,"centered and direct canonical collar positions intersect");
        out.position[k]=sf::range(lo,hi);
        lo=std::max(sf::lower(derivative),sf::lower(direct[k].x));
        hi=std::min(sf::upper(derivative),sf::upper(direct[k].x));
        require(lo<=hi,"centered and direct canonical collar derivatives intersect");
        out.derivative[k]=sf::range(lo,hi);
    }
    return out;
}
inline std::array<S,2> chart_J_centered_coeff(const Patch&patch,const cp::Curve&curve,
        size_t curve_index,Q t0,Q t1){
    require(patch.edge[1].curved&&patch.edge[1].boundary&&
        patch.edge[1].source==curve_index,"same curved collar for centered Jacobian");
    S det=orientation(patch.vertex[0],patch.vertex[1],patch.vertex[2]);
    require(sf::lower(det)>0,"positive thin collar chord orientation");
    S T=sf::range(t0,t1);
    auto relative=curve_relative_centered_bounds(patch,curve,t0,t1);
    const auto&R=relative.position,&Rprime=relative.derivative;
    std::array<S,2>A=patch.vertex[0].x,B=patch.vertex[1].x,C=patch.vertex[2].x;
    std::array<S,2>Dch={C[0]-B[0],C[1]-B[1]},v;
    for(int k=0;k<2;k++)v[k]=B[k]-A[k]+T*Dch[k];
    S correction=cross2(v,Rprime)+S(4)*cross2(R,Dch);
    return {correction/det,S(4)*cross2(R,Rprime)/det};
}
inline S chart_J_centered_interval(const Patch&patch,const cp::Curve&curve,
        size_t curve_index,Q s0,Q s1,Q t0,Q t1){
    auto coefficient=chart_J_centered_coeff(patch,curve,curve_index,t0,t1);
    S s=sf::range(s0,s1),s3=s*s*s,s6=s3*s3;
    S jacobian=S(1)+s3*coefficient[0]+s6*coefficient[1];
    require(__builtin_isfinite(jacobian.n)&&__builtin_isfinite(jacobian.e),
        "finite centered whole-chart Jacobian interval");
    return jacobian;
}
inline Q chart_J_column_lower(const Patch&patch,const cp::Curve&curve,
        size_t curve_index,Q t0,Q t1){
    auto coefficient=chart_J_centered_coeff(patch,curve,curve_index,t0,t1);
    Q a=sf::lower(coefficient[0]),b=sf::lower(coefficient[1]);
    // z=s^3 spans [0,1]; J=1+a z+b z^2. The following cases minimize a
    // downward coefficient bound over the entire radial coordinate.
    if(a>=0&&b>=0)return Q(1);
    if(b<=0)return std::min(Q(1),sf::lower(S(1)+S(a)+S(b)));
    if(a>=0)return Q(1);
    if(sf::lower(S(-a)-S(2)*S(b))>=0)
        return sf::lower(S(1)+S(a)+S(b));
    return sf::lower(S(1)-S(a)*S(a)/(S(4)*S(b)));
}

inline S chart_J_interval(const Patch&patch,const cp::Curve&curve,size_t curve_index,Q s0,Q s1,Q t0,Q t1){
    require(patch.edge[1].curved&&patch.edge[1].boundary&&patch.edge[1].source==curve_index,"one selected source curve on the collar's opposite edge");
    const auto&edge=patch.edge[1];S det=orientation(patch.vertex[0],patch.vertex[1],patch.vertex[2]);
    require(sf::lower(det)>0,"positive whole-chord collar reference orientation");
    S u0=edge.parameter[0],u1=edge.parameter[1],du=u1-u0;
    require(cp::strict_sign(du)!=0,"nonzero certified source-curve parameter span");
    D t(sf::range(t0,t1),S(1),S(0));
    auto relative=curve_relative(patch,curve,t);auto gamma0=curve_position(curve,D(u0));auto gamma1=curve_position(curve,D(u1));
    for(int k=0;k<2;k++){
        require(intersects(gamma0[k].v,patch.vertex[1].x[k])&&intersects(gamma1[k].v,patch.vertex[2].x[k]),"implicit curve endpoint functions agree with shared chart node enclosures");
    }
    std::array<S,2>A=patch.vertex[0].x,B=patch.vertex[1].x,C=patch.vertex[2].x,Dch={C[0]-B[0],C[1]-B[1]};
    S T=sf::range(t0,t1),s=sf::range(s0,s1),s3=s*s*s,s6=s3*s3;
    std::array<S,2>R,Rprime,v;
    for(int k=0;k<2;k++){
        R[k]=relative[k].v;Rprime[k]=relative[k].x;
        v[k]=B[k]-A[k]+T*Dch[k];
    }
    // Divide the determinant by the reference barycentric factor s
    // symbolically. This avoids an artificial 0/0 at the chart apex s=0.
    S correction=cross2(v,Rprime)+S(4)*cross2(R,Dch);
    S jacobian=S(1)+s3*correction/det+S(4)*s6*cross2(R,Rprime)/det;
    require(__builtin_isfinite(jacobian.n)&&__builtin_isfinite(jacobian.e),"finite outward whole-chart Jacobian interval");
    return jacobian;
}

inline JacobianCertificate certify_chart_jacobians(const compact_target::Registry&registry,
        const Candidate&candidate,const JacobianBudget&budget){
    require(candidate.source_registry_hash==registry.wire_hash&&candidate.source_epoch==registry.source_epoch&&
        candidate.source_steps==registry.source_steps&&candidate.source_time==registry.source_time&&
        budget.maximum_boxes_per_patch>0&&budget.maximum_boxes_per_patch<=65536&&budget.maximum_depth>0&&budget.maximum_depth<=64,
        "same canonical source chart and finite continuous-J proof budget");
    JacobianCertificate out;out.patch_lower.reserve(candidate.patch.size());out.minimum_lower=HUGE_VALQ;
    struct Box{Q s0,s1,t0,t1;unsigned depth;};
    for(size_t patch_index=0;patch_index<candidate.patch.size();patch_index++){
        const auto&patch=candidate.patch[patch_index];
        size_t curved=0,curve_edge=SIZE_MAX;for(size_t k=0;k<3;k++)if(patch.edge[k].curved){curved++;curve_edge=k;}
        require(curved<=1&&(!curved||curve_edge==1),"one curved opposite edge on a collar chart");
        if(!curved){Q lower=Q(1);out.patch_lower.push_back(lower);out.minimum_lower=std::min(out.minimum_lower,lower);continue;}
        const auto&descriptor=patch.edge[curve_edge];require(descriptor.source<registry.exterior.size(),"curve chart belongs to one actual source exterior edge");
        auto curve=cp::curve(registry.exterior[descriptor.source]);require(curve.edge==registry.exterior[descriptor.source].edge,"same canonical source curve identity");
        std::vector<Box>stack={{0,1,0,1,0}};Q patch_lower=HUGE_VALQ;size_t patch_boxes=0;
        while(!stack.empty()){
            Box b=stack.back();stack.pop_back();
            if(patch_boxes>=budget.maximum_boxes_per_patch){
                S local=chart_J_interval(patch,curve,descriptor.source,b.s0,b.s1,b.t0,b.t1);
                Q sm=(b.s0+b.s1)/2,tm=(b.t0+b.t1)/2;
                S center=chart_J_interval(patch,curve,descriptor.source,sm,sm,tm,tm);
                S area=orientation(patch.vertex[0],patch.vertex[1],patch.vertex[2]);
                char message[480];std::snprintf(message,sizeof(message),
                    "finite interval-box budget for collar cell=%llu patch=%zu edge=%llu s=[%.9Lg,%.9Lg] t=[%.9Lg,%.9Lg] depth=%u local_J=[%.9Lg,%.9Lg] center_J=[%.9Lg,%.9Lg] area=[%.9Lg,%.9Lg] u=[%.9Lg,%.9Lg]",
                    static_cast<unsigned long long>(candidate.cell),patch_index,
                    static_cast<unsigned long long>(curve.edge),
                    static_cast<long double>(b.s0),static_cast<long double>(b.s1),
                    static_cast<long double>(b.t0),static_cast<long double>(b.t1),b.depth,
                    static_cast<long double>(sf::lower(local)),static_cast<long double>(sf::upper(local)),
                    static_cast<long double>(sf::lower(center)),static_cast<long double>(sf::upper(center)),
                    static_cast<long double>(sf::lower(area)),static_cast<long double>(sf::upper(area)),
                    static_cast<long double>(patch.edge[1].parameter[0].n),
                    static_cast<long double>(patch.edge[1].parameter[1].n));
                throw std::runtime_error(message);
            }
            patch_boxes++;out.boxes++;
            Q lower=chart_J_column_lower(patch,curve,descriptor.source,b.t0,b.t1);
            if(lower>0){patch_lower=std::min(patch_lower,lower);continue;}
            require(b.depth<budget.maximum_depth,"continuous curved-chart Jacobian remains unresolved at the finite proof depth");
            Box left=b,right=b;left.depth=right.depth=b.depth+1;
            Q mid=(b.t0+b.t1)/2;left.t1=mid;right.t0=mid;
            stack.push_back(right);stack.push_back(left);out.subdivisions++;
        }
        require(patch_lower>0,"whole-reference positive curved-chart Jacobian lower bound");out.curved_patches++;out.patch_lower.push_back(patch_lower);out.minimum_lower=std::min(out.minimum_lower,patch_lower);
    }
    require(!out.patch_lower.empty()&&out.minimum_lower>0,"positive continuous Jacobian certificate for every straight and curved chart patch");
    out.all_reference_J_positive=true;return out;
}
inline JacobianCertificate certify_all_chart_jacobians(const compact_target::Registry&registry,
        const std::vector<Candidate>&candidate,const JacobianBudget&budget){
    require(candidate.size()==registry.site.size()&&!candidate.empty(),"one whole-chart Jacobian certificate set for every target site");
    JacobianCertificate out;out.minimum_lower=HUGE_VALQ;out.source_registry_hash=registry.wire_hash;
    out.weight_hash=candidate.front().weight_hash;out.source_epoch=registry.source_epoch;out.source_steps=registry.source_steps;out.source_time=registry.source_time;
    for(size_t j=0;j<candidate.size();j++){
        require(candidate[j].cell==registry.site[j].cell&&candidate[j].source_registry_hash==registry.wire_hash&&
            candidate[j].weight_hash==out.weight_hash&&candidate[j].source_epoch==registry.source_epoch&&
            candidate[j].source_steps==registry.source_steps&&candidate[j].source_time==registry.source_time,
            "sorted full chart collection for one immutable source and weight epoch");
        auto part=certify_chart_jacobians(registry,candidate[j],budget);
        out.boxes+=part.boxes;out.subdivisions+=part.subdivisions;out.curved_patches+=part.curved_patches;
        out.minimum_lower=std::min(out.minimum_lower,part.minimum_lower);out.patch_lower.insert(out.patch_lower.end(),part.patch_lower.begin(),part.patch_lower.end());
    }
    require(!out.patch_lower.empty()&&out.minimum_lower>0,"all target chart Jacobians are positive over their full reference triangles");
    out.chart_hash=chart_collection_hash(candidate);out.all_reference_J_positive=true;return out;
}
}
#endif
