#ifndef MATERIAL_COMPACT_SOURCE_FIELDS_H
#define MATERIAL_COMPACT_SOURCE_FIELDS_H
// Pointwise accepted-source field reconstruction for a known original
// triangle. The caller must construct the geometric overlap and supply a
// reference point with a justified enclosure; this does not certify any
// quadrature error or settle source-edge ties.
#include "material_shared_function.h"
#include <cmath>
#include <stdexcept>

namespace compact_source_fields {
namespace sf=shared_function;
using Q=sf::Q;using S=sf::S;
inline void require(bool ok,const char*why){if(!ok)throw std::runtime_error(why);}

struct Law {
    double gamma=0,P0=0,eps=0;
    double current_delta=0,generating_delta=0;
    // Same accepted Kepler IC, with exact binary64 values used in its native
    // gfs_enrich_reference implementation.
    double rin=2,rout=8,density_floor=.001,radial_soft=1e-6,taper_width=.25;
};
struct BoundedField {
    S physical_x[2],velocity[2],rho_reference,rho_physical,heat_physical,phi;
    S current_j,generating_j;
    bool strictly_inside_source=false;
    // The analytical ball is not a certificate for the native double exp/
    // hypot implementation or for integrated physical-field quadrature.
    bool native_field_error_certified=false,quadrature_error_certified=false;
};
struct Containment {
    S reference_X,reference_Y;
    Q inverse_radius=0;
    bool entire_physical_box_inside=false;
};

inline S nonnegative_root(S x){
    require(sf::upper(x)>=0,"nonnegative radial square enclosure");
    if(sf::lower(x)>0)return sf::mf::root(x);
    Q upper=sf::mf::FP<Q>::root(sf::upper(x));
    upper=sf::mf::FP<Q>::up(upper+sf::mf::libm_error(upper));
    Q middle=upper/2;
    return S(middle,sf::mf::FP<Q>::up(upper-middle));
}
inline S taper(S z){return (S(1)+sf::mf::tanh_(z))/S(2);}
inline S reference_density(S X,S Y,const Law&law){
    require(law.rin>0&&law.rout>law.rin&&law.density_floor>0&&law.radial_soft>0&&law.taper_width>0,
        "positive original radial Kepler density parameters");
    S r=nonnegative_root(X*X+Y*Y);
    S inner=taper((r-sf::exact(law.rin))/sf::exact(law.taper_width));
    S outer=taper((sf::exact(law.rout)-r)/sf::exact(law.taper_width));
    S rho=sf::exact(law.density_floor)+sf::exact(law.rin)/(r+sf::exact(law.radial_soft))*inner*outer;
    require(sf::lower(rho)>0,"positive interval original reference density");
    return rho;
}
inline BoundedField evaluate(const sf::Map&current,const sf::Map&generating,S X,S Y,const Law&law){
    require(current.state==generating.state&&current.source==generating.source&&
        law.gamma>1&&law.P0>0&&law.eps>0&&std::isfinite(law.gamma)&&
        std::isfinite(law.P0)&&std::isfinite(law.eps)&&
        std::isfinite(law.current_delta)&&std::isfinite(law.generating_delta),
        "same original triangle and finite accepted thermal/gravity law");
    auto bary=sf::barycentric(current,X,Y);
    for(int k=0;k<3;k++)require(sf::lower(bary[k].v)>0,
        "source quadrature point strictly inside its known original triangle");
    auto now=sf::evaluate(current,X,Y),initial=sf::evaluate(generating,X,Y);
    BoundedField out;out.current_j=now.J+sf::exact(law.current_delta);
    out.generating_j=initial.J+sf::exact(law.generating_delta);
    require(sf::lower(out.current_j)>0&&sf::lower(out.generating_j)>0,
        "positive accepted current and generating corrected volume factors");
    out.rho_reference=reference_density(X,Y,law);
    out.rho_physical=out.rho_reference/out.current_j;
    // The accepted CUDA RHS stores e0_ref=P0*j0/(gamma-1), then evolves it
    // with (j/j0)^(1-gamma); division by j gives physical heat density.
    out.heat_physical=sf::exact(law.P0)/sf::exact(law.gamma-1)*
        sf::mf::power(out.generating_j/out.current_j,Q(law.gamma));
    for(int k=0;k<2;k++){out.physical_x[k]=now.position[k];out.velocity[k]=now.velocity[k];}
    S r2=out.physical_x[0]*out.physical_x[0]+out.physical_x[1]*out.physical_x[1]+
        sf::exact(law.eps)*sf::exact(law.eps);
    out.phi=-S(1)/sf::mf::root(r2);
    require(sf::lower(out.rho_physical)>0&&sf::lower(out.heat_physical)>0,
        "positive reconstructed physical density and heat");
    out.strictly_inside_source=true;
    return out;
}
inline BoundedField evaluate_at_physical(const sf::Map&current,const sf::Map&generating,
        const sf::RealBounds&geometry,S x,S y,const Law&law){
    require(geometry.L<1,"whole-source contraction bound for a known original triangle");
    auto inverse=sf::inverse_inside(current,geometry,x,y);
    auto out=evaluate(current,generating,inverse.X,inverse.Y,law);
    require(sf::lower(out.physical_x[0])<=sf::upper(x)&&sf::lower(x)<=sf::upper(out.physical_x[0])&&
        sf::lower(out.physical_x[1])<=sf::upper(y)&&sf::lower(y)<=sf::upper(out.physical_x[1]),
        "source inverse and forward maps enclose the same physical query point");
    return out;
}
inline Containment enclose_physical_box(const sf::Map&current,
        const sf::RealBounds&geometry,S x,S y){
    require(current.state&&geometry.L<1&&sf::lower(x)<=sf::upper(x)&&sf::lower(y)<=sf::upper(y),
        "known complete source map and finite physical query box");
    Containment out;
    try{
        auto inverse=sf::inverse_inside(current,geometry,x,y);
        out.reference_X=inverse.X;out.reference_Y=inverse.Y;
        out.inverse_radius=inverse.radius;out.entire_physical_box_inside=inverse.inside;
    }catch(const std::exception&){
        // A failed sufficient enclosure is unresolved, never a proof that
        // the source triangle and target physical box are disjoint.
    }
    return out;
}
}
#endif
