#ifndef MATERIAL_COMPACT_MOMENTS_H
#define MATERIAL_COMPACT_MOMENTS_H
// Algebraic source/target accounting on one already validated overlap stream.
// Physical sample construction and quadrature/field error certificates belong
// to the future actual overlap and source-field evaluators.
#include "material_compact_overlap.h"
#include <utility>

namespace compact_overlap {
enum MomentIndex : size_t {M=0,QX,QY,QXX,QXY,QYY,PX,PY,L,K,U,PHI,NMOMENTS};
using Moment=std::array<R,NMOMENTS>;

struct PhysicalSample {
    uint64_t point_id=UINT64_MAX;
    R rho=0,internal_energy_density=0,vx=0,vy=0,phi=0;
};
struct FieldStream {
    Seal seal;
    std::vector<PhysicalSample>sample;
};
struct Ledger {
    std::vector<Moment>outgoing_source_triangle,incoming_target_cell;
    Moment packet_total{};
    Marginals area;
    R maximum_roundoff_gap=0;
    bool source_target_algebraic_balance=false;
    bool physical_quadrature_error_certified=false;
};

class StreamAccumulator {
    Validator validator_;
    Ledger out_;
    bool finished_=false;
public:
    StreamAccumulator(const Seal&expected_seal,const Seal&packet_seal,const Seal&field_seal,
        const std::vector<uint64_t>&source_triangle_ids,const std::vector<uint64_t>&target_cell_ids,
        const std::vector<uint64_t>&target_patch_cell_ids,
        const std::vector<R>&expected_source_area,const std::vector<R>&expected_source_area_error,
        const std::vector<R>&expected_target_area,const std::vector<R>&expected_target_area_error,
        const Budget&budget):validator_(expected_seal,packet_seal,source_triangle_ids,target_cell_ids,
            target_patch_cell_ids,expected_source_area,expected_source_area_error,
            expected_target_area,expected_target_area_error,budget){
        require_same_seal(packet_seal,field_seal);
        out_.outgoing_source_triangle.resize(source_triangle_ids.size());
        out_.incoming_target_cell.resize(target_cell_ids.size());
    }
    void append(const Point&p,const PhysicalSample&f){
        require(!finished_,"physical moment stream remains open");
        require(f.point_id==p.id,"field sample and common packet point IDs agree in order");
        require(f.rho>0&&std::isfinite(f.rho)&&f.internal_energy_density>0&&
            std::isfinite(f.internal_energy_density)&&std::isfinite(f.vx)&&
            std::isfinite(f.vy)&&std::isfinite(f.phi),"positive finite source physical density and heat; finite velocity and potential");
        // The source coordinate is the single integration location. The
        // target coordinate is only an enclosure check, never a second RHS.
        const R x=p.source_x[0],y=p.source_x[1],w=p.physical_area_weight;
        const R m=w*f.rho,px=m*f.vx,py=m*f.vy;
        const Moment v={m,m*x,m*y,m*x*x,m*x*y,m*y*y,px,py,x*py-y*px,
            R(.5)*m*(f.vx*f.vx+f.vy*f.vy),w*f.internal_energy_density,m*f.phi};
        for(R a:v)require(std::isfinite(a),"finite physical moment contribution");
        const auto index=validator_.append(p);
        for(size_t j=0;j<NMOMENTS;j++){
            out_.outgoing_source_triangle[index.source][j]+=v[j];
            out_.incoming_target_cell[index.target][j]+=v[j];
            out_.packet_total[j]+=v[j];
        }
    }
    Ledger finish(){
        require(!finished_,"physical moment stream finalized exactly once");finished_=true;
        out_.area=validator_.finish();
        for(size_t j=0;j<NMOMENTS;j++){
            R outgoing=0,incoming=0;
            for(const Moment&m:out_.outgoing_source_triangle)outgoing+=m[j];
            for(const Moment&m:out_.incoming_target_cell)incoming+=m[j];
            const R gap=fabsl(outgoing-incoming);
            out_.maximum_roundoff_gap=std::max(out_.maximum_roundoff_gap,gap);
            const R scale=std::max({fabsl(outgoing),fabsl(incoming),fabsl(out_.packet_total[j]),R(1)});
            require(gap<=128*std::numeric_limits<R>::epsilon()*scale,
                "common packet source/target moment sums agree to bounded floating-point roundoff");
        }
        out_.source_target_algebraic_balance=true;
        // The shared stream proves an algebraic identity. It does not prove
        // true continuous integrals or projected-field accuracy.
        out_.physical_quadrature_error_certified=false;
        return std::move(out_);
    }
};

inline Ledger accumulate(const Seal&expected_seal,const Seal&packet_seal,
        const std::vector<Point>&points,const FieldStream&fields,
        const std::vector<uint64_t>&source_triangle_ids,const std::vector<uint64_t>&target_cell_ids,
        const std::vector<uint64_t>&target_patch_cell_ids,
        const std::vector<R>&expected_source_area,const std::vector<R>&expected_source_area_error,
        const std::vector<R>&expected_target_area,const std::vector<R>&expected_target_area_error,
        const Budget&budget){
    require(fields.sample.size()==points.size(),"one source-field sample per common overlap packet");
    StreamAccumulator stream(expected_seal,packet_seal,fields.seal,source_triangle_ids,
        target_cell_ids,target_patch_cell_ids,expected_source_area,expected_source_area_error,
        expected_target_area,expected_target_area_error,budget);
    for(size_t i=0;i<points.size();i++)stream.append(points[i],fields.sample[i]);
    return stream.finish();
}
}
#endif
