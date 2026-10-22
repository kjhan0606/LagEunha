#ifndef MATERIAL_COMPACT_CHART_AREA_H
#define MATERIAL_COMPACT_CHART_AREA_H
// Positive physical-area nodes for the selected chart and an independent
// Green-boundary area enclosure. This is target geometry integration, not a
// source/target overlap, physical-field quadrature, or coverage certificate.
#include "material_compact_chart_jacobian.h"
#include "material_shared_trace_factory.h"

namespace compact_chart {
struct AreaBudget {
    size_t maximum_interior_boxes_per_patch=65536,maximum_boundary_boxes_per_patch=65536;
    unsigned maximum_depth=32;
    Q maximum_absolute_area_error_per_patch=1e-3Q;
};
struct AreaNode {
    Q s=0,t=0,weight=0,error=0;
    S physical_x[2];
};
struct AreaRule {
    uint64_t chart_hash=0,source_registry_hash=0,weight_hash=0,source_epoch=0,source_steps=0,law_hash=0;
    double source_time=0;
    size_t flat_patch=0,interior_boxes=0,interior_subdivisions=0,boundary_boxes=0,boundary_subdivisions=0;
    std::vector<AreaNode>node;
    S interior_area,boundary_area;
    Q common_lower=0,common_upper=0,minimum_node_weight_lower=0;
    bool positive_nodes=false,independent_area_agreement=false,continuous_physical_coverage=false;
};
struct GaussNode {
    double l[3]={};
    Q s=0,t=0;
    S reference_weight,physical_weight,physical_x[2];
};
struct GaussRule {
    uint64_t chart_hash=0,source_epoch=0,weight_hash=0;
    size_t flat_patch=0;
    int order=0;
    std::vector<GaussNode>node;
    S reference_area,physical_area;
    bool positive_reference_and_physical_weights=false;
    // Order comparison is an accuracy diagnostic, not a certified bound on
    // arbitrary source-field integrals or a continuous overlap certificate.
    bool field_integration_error_certified=false,continuous_physical_coverage=false;
};

inline S area_cross(const std::array<S,2>&a,const std::array<S,2>&b){return a[0]*b[1]-a[1]*b[0];}
inline std::array<S,2> chart_position_box(const Patch&patch,const cp::Curve*curve,S ss,S tt){
    std::array<S,2>out;S s2=ss*ss,s4=s2*s2;
    const auto&A=patch.vertex[0].x,&B=patch.vertex[1].x,&C=patch.vertex[2].x;
    for(int k=0;k<2;k++)out[k]=A[k]+ss*(B[k]-A[k]+tt*(C[k]-B[k]));
    if(curve){
        auto relative=curve_relative(patch,*curve,D(tt));
        for(int k=0;k<2;k++)out[k]=out[k]+s4*relative[k].v;
    }
    return out;
}
inline std::array<S,2> chart_position(const Patch&patch,const cp::Curve*curve,Q s,Q t){
    return chart_position_box(patch,curve,S(s),S(t));
}
inline S area_box(const Patch&patch,const cp::Curve*curve,size_t curve_index,Q s0,Q s1,Q t0,Q t1){
    require(s0>=0&&s0<s1&&s1<=1&&t0>=0&&t0<t1&&t1<=1,"finite chart-area box lies inside the whole reference triangle");
    S determinant=orientation(patch.vertex[0],patch.vertex[1],patch.vertex[2]);
    require(sf::lower(determinant)>0,"positive reference-chord determinant for physical area");
    S relative_J=curve?chart_J_interval(patch,*curve,curve_index,s0,s1,t0,t1):S(1);
    S integral_s=(S(s1)*S(s1)-S(s0)*S(s0))/S(2);
    return determinant*relative_J*integral_s*S(t1-t0);
}
inline S boundary_integrand(const Patch&patch,const cp::Curve&curve,Q t0,Q t1){
    auto relative=curve_relative_centered_bounds(patch,curve,t0,t1);
    S dx=patch.vertex[2].x[0]-patch.vertex[1].x[0],dy=patch.vertex[2].x[1]-patch.vertex[1].x[1];
    return relative.position[0]*dy-relative.position[1]*dx+
        (relative.position[0]*relative.derivative[1]-
         relative.position[1]*relative.derivative[0])/S(2);
}
inline S boundary_area(const Patch&patch,const cp::Curve*curve,const AreaBudget&budget,
        size_t&boxes,size_t&subdivisions){
    const auto&A=patch.vertex[0].x,&B=patch.vertex[1].x,&C=patch.vertex[2].x;
    S correction;
    if(curve){
        struct Interval{Q lo,hi;unsigned depth;};std::vector<Interval>stack={{0,1,0}};
        while(!stack.empty()){
            auto b=stack.back();stack.pop_back();
            require(boxes<budget.maximum_boundary_boxes_per_patch,"finite independent curved-boundary area box budget");boxes++;
            S integral=boundary_integrand(patch,*curve,b.lo,b.hi)*S(b.hi-b.lo);
            Q width=sf::mf::FP<Q>::up(sf::upper(integral)-sf::lower(integral));
            Q allowed=sf::mf::FP<Q>::down(budget.maximum_absolute_area_error_per_patch*(b.hi-b.lo));
            if(width<=allowed){correction=correction+integral;continue;}
            require(b.depth<budget.maximum_depth,"independent curved-boundary area unresolved at finite proof depth");
            Q midpoint=(b.lo+b.hi)/2;require(midpoint>b.lo&&midpoint<b.hi,"representable boundary area subdivision");
            stack.push_back({midpoint,b.hi,b.depth+1});stack.push_back({b.lo,midpoint,b.depth+1});subdivisions++;
        }
    }
    return (area_cross(A,B)+area_cross(B,C)+area_cross(C,A))/S(2)+correction;
}

class AreaBuilder {
    const compact_target::Registry&registry_;
    const std::vector<Candidate>&candidate_;
    const Certificate&topology_;
    const JacobianCertificate&jacobian_;
    const AreaBudget budget_;
public:
    AreaBuilder(const compact_target::Registry&registry,const std::vector<Candidate>&candidate,
            const Certificate&topology,const JacobianCertificate&jacobian,const AreaBudget&budget):
        registry_(registry),candidate_(candidate),topology_(topology),jacobian_(jacobian),budget_(budget){
        require_same_chart_complex(registry_,candidate_,topology_);
        require(topology_.disk_topology&&topology_.exact_C0_connectivity&&topology_.oriented_source_boundary_complete&&
            jacobian_.all_reference_J_positive&&jacobian_.chart_hash==topology_.chart_hash&&
            jacobian_.source_registry_hash==topology_.source_registry_hash&&jacobian_.weight_hash==topology_.weight_hash&&
            jacobian_.source_epoch==topology_.source_epoch&&jacobian_.source_steps==topology_.source_steps&&
            jacobian_.source_time==topology_.source_time&&jacobian_.patch_lower.size()==topology_.patch.size(),
            "same sealed positive-J disk chart before physical-area quadrature");
        require(budget_.maximum_interior_boxes_per_patch>0&&budget_.maximum_interior_boxes_per_patch<=1048576&&
            budget_.maximum_boundary_boxes_per_patch>0&&budget_.maximum_boundary_boxes_per_patch<=1048576&&
            budget_.maximum_depth>0&&budget_.maximum_depth<=64&&
            budget_.maximum_absolute_area_error_per_patch>0&&__builtin_isfinite(budget_.maximum_absolute_area_error_per_patch),
            "explicit finite positive physical-area integration budgets");
    }
    AreaRule patch_rule(size_t flat_patch)const{
        require(flat_patch<topology_.patch.size()&&jacobian_.patch_lower[flat_patch]>0,"positive-J chart patch in the sealed global order");
        const auto&ref=topology_.patch[flat_patch];const auto&patch=candidate_[ref.candidate].patch[ref.patch];
        size_t curved_count=0;for(const auto&e:patch.edge)curved_count+=e.curved;
        require(curved_count<=1&&(!curved_count||patch.edge[1].curved),"one curved opposite edge in physical-area chart");
        cp::Curve descriptor;const cp::Curve*curve=nullptr;size_t curve_index=SIZE_MAX;
        if(curved_count){curve_index=patch.edge[1].source;require(curve_index<registry_.exterior.size(),"selected original exterior curve for physical-area rule");
            descriptor=cp::curve(registry_.exterior[curve_index]);curve=&descriptor;}
        AreaRule out;out.chart_hash=topology_.chart_hash;out.source_registry_hash=topology_.source_registry_hash;
        out.weight_hash=topology_.weight_hash;out.source_epoch=topology_.source_epoch;out.source_steps=topology_.source_steps;
        out.source_time=topology_.source_time;out.law_hash=topology_.law_hash;out.flat_patch=flat_patch;
        struct Box{Q s0,s1,t0,t1;unsigned depth;};std::vector<Box>stack={{0,1,0,1,0}};
        out.minimum_node_weight_lower=HUGE_VALQ;
        while(!stack.empty()){
            auto b=stack.back();stack.pop_back();
            require(out.interior_boxes<budget_.maximum_interior_boxes_per_patch,"finite interior physical-area box budget");out.interior_boxes++;
            S area=area_box(patch,curve,curve_index,b.s0,b.s1,b.t0,b.t1);
            Q lo=sf::lower(area),hi=sf::upper(area),width=sf::mf::FP<Q>::up(hi-lo);
            Q fraction=(b.s1-b.s0)*(b.t1-b.t0);
            Q allowed=sf::mf::FP<Q>::down(budget_.maximum_absolute_area_error_per_patch*fraction);
            if(lo>0&&width<=allowed){
                AreaNode node;node.s=(b.s0+b.s1)/2;node.t=(b.t0+b.t1)/2;
                auto position=chart_position(patch,curve,node.s,node.t);for(int k=0;k<2;k++)node.physical_x[k]=position[k];
                node.weight=(lo+hi)/2;node.error=std::max(sf::mf::FP<Q>::up(node.weight-lo),sf::mf::FP<Q>::up(hi-node.weight));
                require(node.weight>node.error&&node.error>=0,"strictly positive physical-area quadrature node and enclosing error");
                out.minimum_node_weight_lower=std::min(out.minimum_node_weight_lower,lo);
                out.interior_area=out.interior_area+area;out.node.push_back(node);continue;
            }
            require(b.depth<budget_.maximum_depth,"interior physical-area quadrature unresolved at finite proof depth");
            Box left=b,right=b;left.depth=right.depth=b.depth+1;
            if(b.s1-b.s0>b.t1-b.t0){Q middle=(b.s0+b.s1)/2;require(middle>b.s0&&middle<b.s1,"representable s subdivision");left.s1=middle;right.s0=middle;}
            else{Q middle=(b.t0+b.t1)/2;require(middle>b.t0&&middle<b.t1,"representable t subdivision");left.t1=middle;right.t0=middle;}
            stack.push_back(right);stack.push_back(left);out.interior_subdivisions++;
        }
        require(!out.node.empty()&&sf::lower(out.interior_area)>0,"positive interior area on the complete reference patch");
        out.boundary_area=boundary_area(patch,curve,budget_,out.boundary_boxes,out.boundary_subdivisions);
        out.common_lower=std::max(sf::lower(out.interior_area),sf::lower(out.boundary_area));
        out.common_upper=std::min(sf::upper(out.interior_area),sf::upper(out.boundary_area));
        require(out.common_lower>0&&out.common_lower<=out.common_upper,
            "independent boundary and interior physical-area enclosures agree");
        out.positive_nodes=true;out.independent_area_agreement=true;
        out.continuous_physical_coverage=false;return out;
    }
    GaussRule positive_gauss_rule(size_t flat_patch,int order)const{
        require(flat_patch<topology_.patch.size()&&jacobian_.patch_lower[flat_patch]>0&&order>=2&&order<=64,
            "finite positive Gauss order on a sealed whole-chart positive-J patch");
        const auto&ref=topology_.patch[flat_patch];const auto&patch=candidate_[ref.candidate].patch[ref.patch];
        size_t curved_count=0;for(const auto&e:patch.edge)curved_count+=e.curved;
        require(curved_count<=1&&(!curved_count||patch.edge[1].curved),"one curved opposite edge for positive Gauss rule");
        cp::Curve descriptor;const cp::Curve*curve=nullptr;size_t curve_index=SIZE_MAX;
        if(curved_count){curve_index=patch.edge[1].source;require(curve_index<registry_.exterior.size(),"same accepted exterior curve for Gauss rule");
            descriptor=cp::curve(registry_.exterior[curve_index]);curve=&descriptor;}
        std::vector<double>x,w;shared_trace::gauss(order,x,w);
        require(x.size()==size_t(order)&&w.size()==x.size(),"complete original positive Gauss rule");
        GaussRule out;out.chart_hash=topology_.chart_hash;out.source_epoch=topology_.source_epoch;
        out.weight_hash=topology_.weight_hash;out.flat_patch=flat_patch;out.order=order;
        out.node.reserve(size_t(order)*order);
        S determinant=orientation(patch.vertex[0],patch.vertex[1],patch.vertex[2]);
        require(sf::lower(determinant)>0,"positive canonical reference triangle for Gauss rule");
        for(int i=0;i<order;i++)for(int j=0;j<order;j++){
            require(x[i]>0&&x[i]<1&&x[j]>0&&x[j]<1&&w[i]>0&&w[j]>0&&
                std::isfinite(x[i])&&std::isfinite(x[j])&&std::isfinite(w[i])&&std::isfinite(w[j]),
                "interior positive original Gauss nodes and weights");
            GaussNode node;node.s=Q(x[i]);node.t=Q(x[j]);
            node.l[0]=1-x[i];node.l[1]=x[i]*(1-x[j]);node.l[2]=x[i]*x[j];
            node.reference_weight=determinant*S(node.s)*sf::exact(w[i])*sf::exact(w[j]);
            S relative_J=curve?chart_J_interval(patch,*curve,curve_index,node.s,node.s,node.t,node.t):S(1);
            node.physical_weight=node.reference_weight*relative_J;
            require(sf::lower(node.reference_weight)>0&&sf::lower(node.physical_weight)>0,
                "positive reference and physical measure at every curved Gauss node");
            auto point=chart_position(patch,curve,node.s,node.t);
            for(int k=0;k<2;k++)node.physical_x[k]=point[k];
            out.reference_area=out.reference_area+node.reference_weight;
            out.physical_area=out.physical_area+node.physical_weight;
            out.node.push_back(node);
        }
        out.positive_reference_and_physical_weights=true;return out;
    }
};
}
#endif
