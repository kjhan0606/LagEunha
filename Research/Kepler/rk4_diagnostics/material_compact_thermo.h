#ifndef MATERIAL_COMPACT_THERMO_H
#define MATERIAL_COMPACT_THERMO_H
#include "material_compact_geometry.h"
#include "material_full_initial_batch.h"
#ifdef __CUDACC__
#define COMPACT_THERMO_HD __host__ __device__
#else
#define COMPACT_THERMO_HD
#endif
namespace compact_thermo {
struct Config {double gamma=0,background_pressure=0,eps=0;};
// Samples are grouped by patch with 64-bit offsets, reference measure dX.
// Explicit positive functional: NOT automatically certified by its point count.
struct QuadraturePoint {double l[3],weight;};
struct Thermal {
    double j=0,j0=0,energy_reference=0,energy_dot_reference=0;
    double pressure=0,delta_pressure=0,sound_squared=0;
    int valid=0;
};
COMPACT_THERMO_HD static inline Thermal evaluate(const Config&cfg,double e0_reference,
    double rho_reference,double j0,double j,double jdot){
    Thermal out;
    if(!(cfg.gamma>1)||!(cfg.background_pressure>0)||!(e0_reference>0)||!(rho_reference>0)||!(j0>0)||!(j>0))return out;
    if(!isfinite(cfg.gamma)||!isfinite(cfg.background_pressure)||!isfinite(e0_reference)||!isfinite(rho_reference)||!isfinite(j0)||!isfinite(j)||!isfinite(jdot))return out;
    double lr=log1p((j-j0)/j0);
    // e0 is the immutable PROJECTED physical internal energy per reference
    // volume, never reset to background_pressure*j0/(gamma-1).
    double generating_pressure=(cfg.gamma-1)*e0_reference/j0;
    out.j=j;out.j0=j0;
    out.delta_pressure=(generating_pressure-cfg.background_pressure)+generating_pressure*expm1(-cfg.gamma*lr);
    out.pressure=cfg.background_pressure+out.delta_pressure;
    out.energy_reference=e0_reference*exp((1-cfg.gamma)*lr);
    out.energy_dot_reference=-out.pressure*jdot;
    out.sound_squared=cfg.gamma*out.pressure*j/rho_reference;
    out.valid=isfinite(lr)&&out.pressure>0&&isfinite(out.pressure)&&out.energy_reference>0&&isfinite(out.energy_reference)&&isfinite(out.energy_dot_reference)&&out.sound_squared>0&&isfinite(out.sound_squared);
    return out;
}
} // namespace compact_thermo
#undef COMPACT_THERMO_HD
#endif
