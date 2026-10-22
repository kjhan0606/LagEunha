#ifndef MATERIAL_SHARED_REAL_SLICE_H
#define MATERIAL_SHARED_REAL_SLICE_H
#include "material_shared_clip.h"
namespace shared_target {
// A real boundary solve depends on the free Gauss coordinate, not eta.
// Reuse the exact same input-family enclosure for every eta on that slice.
// Runtime-only cache. Never serialize its pointers or reuse it across maps,
// chart/strip objects, or a reconstructed RK stage.
struct RealSlice {
    S free,lower,upper,width;
    const sd::Compiled*map=nullptr;
    const sd::Chart*chart=nullptr;
    const Axis*axis=nullptr;
    const Strip*strip=nullptr;
};
static RealSlice real_slice(const sd::Chart&chart,const std::vector<Face>&faces,const Axis&a,const Strip&strip,S free){
    need(!faces.empty()&&faces.front().map&&strip.dependent==a.dependent,"real slice binds one actual map and strip axis");
    RealSlice out;out.free=free;out.lower=S(0);out.upper=S(1);
    out.map=faces.front().map;out.chart=&chart;out.axis=&a;out.strip=&strip;
    for(int index:{strip.lower,strip.upper})if(index>=0)need(size_t(index)<faces.size()&&faces[index].map==out.map,"cached real branches belong to one map");
    if(strip.lower>=0)out.lower=graph(faces[strip.lower],chart,a,strip.lower,free);
    if(strip.upper>=0)out.upper=graph(faces[strip.upper],chart,a,strip.upper,free);
    out.width=out.upper-out.lower;need(sf::lower(out.width)>0,"positive retained real dependent-coordinate slice");return out;
}
static StripPoint evaluate_real_slice(const sd::Compiled&m,const sd::Chart&chart,const Axis&a,const Strip&strip,const RealSlice&slice,S eta){
    need(slice.map==&m&&slice.chart==&chart&&slice.axis==&a&&slice.strip==&strip,"real slice must be consumed within its original immutable map/chart/strip lifetime");
    S dependent=slice.lower+eta*slice.width;StripPoint out;
    out.s=a.dependent==1?slice.free:dependent;out.t=a.dependent==1?dependent:slice.free;
    out.value=[&]{sd::TracePhase phase(44);return sd::evaluate(m,chart,out.s,out.t);}();
    out.jacobian=out.value.reference_jacobian*slice.width*(strip.hi-strip.lo);
    need(sf::lower(out.s)>0&&sf::lower(out.t)>0&&sf::upper(out.s)<1&&sf::upper(out.t)<1&&sf::lower(out.jacobian)>0,"positive interior cached real slice with original orientation");return out;
}
}
#endif
