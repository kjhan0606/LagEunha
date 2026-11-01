#ifndef GFS_WEIGHT_STATE_H
#define GFS_WEIGHT_STATE_H
#include <stdint.h>
#include <math.h>

/* W has units L^2; inertia has units mass/L^2. Inertia is persistent and
 * fixed after explicit initialization, never inferred again on restart. */
#define GFS_WEIGHT_STATE_TAG UINT64_C(0x574754494e455231)
typedef struct {
    double W, pi, inertia;
    uint64_t tag;
} GfsWeightState;

/* Particle-owned stage increments survive the native particle permutation
 * and whole-payload rollback. No ID-order sidecar or legacy scratch alias. */
typedef struct {
    double k1W,k1pi,k2W,k2pi,k3W,k3pi,k4W,k4pi;
} GfsWeightRK4;

static inline double gfs_weight_energy(const GfsWeightState *s){
    return 0.5*s->pi*(s->pi/s->inertia);
}

static inline int gfs_weight_state_valid(const GfsWeightState *s){
    return s->tag==GFS_WEIGHT_STATE_TAG && isfinite(s->W) && isfinite(s->pi) &&
        isfinite(s->inertia) && s->inertia>0 && isfinite(gfs_weight_energy(s));
}
#endif
