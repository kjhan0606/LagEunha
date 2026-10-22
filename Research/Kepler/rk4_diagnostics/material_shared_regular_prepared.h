#ifndef MATERIAL_SHARED_REGULAR_PREPARED_H
#define MATERIAL_SHARED_REGULAR_PREPARED_H
#include "material_shared_weighted_angular.h"
#include "material_shared_origin_quadratic.h"
namespace prepared_regular {
struct Boundary {uint64_t slot,source,target,chart,panel,axis;S free;target::RegularGraph lo,hi;target::Implicit radial_lo,radial_hi;};
static Boundary prepare(const sd::Compiled&m,const std::vector<target::Face>&faces,const target::ClippedChart&chart,const TypedPanel&b,int axis,S free){Boundary z={};z.slot=b.slot;z.source=b.source;z.target=b.target;z.chart=b.geometry.strip.chart;z.panel=b.index;z.axis=axis;z.free=free;const auto&strip=b.geometry.strip;
 if(strip.dependent==1){need(sf::lower(free)>0,"prepared regular angular graph radial pole excluded");z.hi.coordinate.slope=S(1);auto graph=[&](int i){auto center=target::graph(faces[i],chart.reference,chart.coordinate,i,S(free.n));return target::graph_q_disk(faces[i],chart.reference,chart.coordinate,i,free,center);};if(strip.lower>=0)z.lo=graph(strip.lower);if(strip.upper>=0)z.hi=graph(strip.upper);}
 else{z.radial_hi.value=S(1);auto graph=[&](int i){auto center=target::graph(faces[i],chart.reference,chart.coordinate,i,S(free.n));return target::graph_disk(faces[i],chart.reference,chart.coordinate,i,free,center);};if(strip.lower>=0)z.radial_lo=graph(strip.lower);if(strip.upper>=0)z.radial_hi=graph(strip.upper);}
 return z;
}
static target::StripPoint point(const sd::Compiled&m,const target::ClippedChart&chart,const target::Strip&strip,S free,S eta,const Boundary&b){target::StripPoint out;auto certificate=[&](const target::Implicit&p){out.contraction=std::max(out.contraction,p.contraction);out.invariance=std::max(out.invariance,p.invariance);out.forcing=std::max(out.forcing,p.forcing);out.correction_radius=std::max(out.correction_radius,p.correction_radius);out.prior_radius=std::max(out.prior_radius,p.prior_radius);out.real_root_radius=std::max(out.real_root_radius,p.real_root_radius);};
 if(strip.dependent==1){certificate(b.lo.certificate);certificate(b.hi.certificate);sd::AffineAngular mixed=b.lo.coordinate;mixed.interpolated=true;mixed.eta=eta;mixed.width_slope=b.hi.coordinate.slope-b.lo.coordinate.slope;mixed.width_intercept=b.hi.coordinate.intercept-b.lo.coordinate.intercept;mixed.width_correction=b.hi.coordinate.correction-b.lo.coordinate.correction;
  S width=free*mixed.width_slope+mixed.width_intercept+mixed.width_correction;out.s=free;out.t=target::weighted_angular(b.lo.coordinate,b.hi.coordinate,free,eta);out.value=[&]{sd::TracePhase phase(54);return sd::evaluate(m,chart.reference,free,out.t,&mixed);}();out.jacobian=chart.reference.det*width*(strip.hi-strip.lo);
 }else{certificate(b.radial_lo);certificate(b.radial_hi);S width=b.radial_hi.value-b.radial_lo.value;out.s=b.radial_lo.value+eta*width;out.t=free;out.value=[&]{sd::TracePhase phase(44);return sd::evaluate(m,chart.reference,out.s,out.t);}();out.jacobian=out.value.reference_jacobian*width*(strip.hi-strip.lo);}
 return out;
}
}
#endif
