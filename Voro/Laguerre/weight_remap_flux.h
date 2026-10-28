#ifndef LAGFORCE_WEIGHT_REMAP_FLUX_H
#define LAGFORCE_WEIGHT_REMAP_FLUX_H

/* Experimental pure weight-remap flux, NOT enabled in production hydro.
 * Positive rate transfers mass from j to i. Apply the returned tuple to i
 * and its negative to j. No pressure force or gravity term is included.
 *
 * CENTERED: dp = dm*(vi+vj)/2 implies zero pair kinetic-energy derivative
 * when K=sum |p|^2/(2m). Choosing dE=dm*(mean(e)+vi.vj/2) also transports
 * internal energy without donor mixing at the semidiscrete level.
 * KINETIC_THERMAL_DONOR uses the same centered momentum and cross kinetic
 * flux, but upwinds specific thermal energy. When rate is geometric swept
 * volume times donor density, it preserves uniform pressure under pure
 * weight remapping at the semidiscrete level, including density jumps.
 * Caller must select one of the modes below and supply finite states.
 * This is NOT a positivity, entropy, angular-momentum, finite-step, or
 * geometric-consistency guarantee. Do not replace the hydro solver with it.
 */
#ifdef __CUDACC__
#define WR_INLINE static __host__ __device__ __forceinline__
#else
#define WR_INLINE static inline
#endif

enum { WR_DONOR = 0, WR_CENTERED = 1, WR_KINETIC_THERMAL_DONOR = 2 };
typedef struct { double mass, px, py, energy; } WrFlux;

WR_INLINE WrFlux wr_pair_flux(int mode, double rate,
        double vix, double viy, double vjx, double vjy,
        double ei, double ej) {
    WrFlux f;
    f.mass = rate;
    if (mode == WR_CENTERED || mode == WR_KINETIC_THERMAL_DONOR) {
        double thermal = mode == WR_CENTERED ? 0.5 * (ei + ej)
                         : (rate >= 0.0 ? ej : ei);
        f.px = rate * (0.5 * (vix + vjx));
        f.py = rate * (0.5 * (viy + vjy));
        f.energy = rate * (thermal + 0.5 * (vix * vjx + viy * vjy));
    } else {
        double vx = rate >= 0.0 ? vjx : vix;
        double vy = rate >= 0.0 ? vjy : viy;
        double e = rate >= 0.0 ? ej : ei;
        f.px = rate * vx;
        f.py = rate * vy;
        f.energy = rate * (e + 0.5 * (vx * vx + vy * vy));
    }
    return f;
}

#undef WR_INLINE
#endif
