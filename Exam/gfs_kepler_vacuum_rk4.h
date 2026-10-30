#ifndef GFS_KEPLER_VACUUM_RK4_H
#define GFS_KEPLER_VACUUM_RK4_H

#include "gfs_kepler_vacuum_stencil.h"
#include "gfs_vacuum_riemann.h"

typedef struct {
    double rho, pressure, vx, vy;
} GfsKeplerVacuumFaceTrace;

/* RK4 trial coordinates for the geometry-only vacuum generators. The array
 * is deliberately separate from VORO_TBP/TBPP: those buffers are rebuilt at
 * every stage and include MPI halos. This helper manages kinematics only; it
 * does not determine boundary velocities, migrate ownership, or checkpoint
 * the state. The caller must provide one consistent velocity per stable ID. */
typedef struct {
    int stable_id;
    double x, y;
    double kx[4], ky[4];
} GfsKeplerVacuumRK4Site;

typedef struct {
    GfsKeplerVacuumRK4Site *sites;
    size_t count;
    int active;
    int next_stage;
} GfsKeplerVacuumRK4;

/* Compute the kinematic velocity of each shared mirror generator from the
 * exact planar gas-to-vacuum edge speed. `face_trace` is laid out as four
 * entries per stencil site, matching stencil->incidences; only entries below
 * incidence_count are read. Traces must be reconstructed on the gas side of
 * each face. On failure, output arrays are scratch and must be discarded. */
static inline int gfs_kepler_vacuum_rk4_boundary_velocities(
    const GfsKeplerVacuumStencil *stencil,
    const GfsKeplerVacuumFaceTrace *face_trace, double gamma,
    double *site_vx, double *site_vy, double *edge_speed)
{
    if (!stencil || !stencil->sites || !stencil->incidences ||
        stencil->n_sites == 0 || !face_trace || !site_vx || !site_vy ||
        !edge_speed || !isfinite(gamma) || !(gamma > 1.0))
        return 0;

    for (size_t s = 0; s < stencil->n_sites; ++s) {
        const GfsKeplerVacuumSite *site = &stencil->sites[s];
        const GfsKeplerVacuumIncidence *inc = &stencil->incidences[4 * s];
        if (site->incidence_count == 1) {
            const GfsKeplerVacuumFaceTrace *trace = &face_trace[4 * s];
            if (!gfs_vacuum_mirror_kinematics(
                    trace->rho, trace->pressure, trace->vx, trace->vy,
                    gamma, (double)inc[0].nx, (double)inc[0].ny,
                    &edge_speed[4 * s], &site_vx[s], &site_vy[s]))
                return 0;
        } else if (site->incidence_count == 2) {
            const GfsKeplerVacuumFaceTrace *t0 = &face_trace[4 * s];
            const GfsKeplerVacuumFaceTrace *t1 = &face_trace[4 * s + 1];
            if (!isfinite(t0->rho) || !isfinite(t0->pressure) ||
                !isfinite(t0->vx) || !isfinite(t0->vy) ||
                !isfinite(t1->rho) || !isfinite(t1->pressure) ||
                !isfinite(t1->vx) || !isfinite(t1->vy) ||
                !(t0->rho > 0.0) || !(t0->pressure > 0.0) ||
                !(t1->rho > 0.0) || !(t1->pressure > 0.0))
                return 0;
            double n0 = hypot((double)inc[0].nx, (double)inc[0].ny);
            double n1 = hypot((double)inc[1].nx, (double)inc[1].ny);
            if (!(n0 > 0.0) || !(n1 > 0.0)) return 0;
            double un0 = (t0->vx * (double)inc[0].nx +
                          t0->vy * (double)inc[0].ny) / n0;
            double un1 = (t1->vx * (double)inc[1].nx +
                          t1->vy * (double)inc[1].ny) / n1;
            double c0 = sqrt(gamma * t0->pressure / t0->rho);
            double c1 = sqrt(gamma * t1->pressure / t1->rho);
            if (!gfs_vacuum_mirror_kinematics_pair(
                    un0, c0, (double)inc[0].nx, (double)inc[0].ny,
                    un1, c1, (double)inc[1].nx, (double)inc[1].ny,
                    gamma, &edge_speed[4 * s], &edge_speed[4 * s + 1],
                    &site_vx[s], &site_vy[s]))
                return 0;
        } else {
            return 0;
        }
    }
    return 1;
}

static inline void gfs_kepler_vacuum_rk4_destroy(
    GfsKeplerVacuumRK4 *state)
{
    if (!state) return;
    free(state->sites);
    *state = (GfsKeplerVacuumRK4){0};
}

/* Copy the audited stencil into persistent stage state. `out` must be zero
 * initialized (or previously destroyed) before its first use. */
static inline int gfs_kepler_vacuum_rk4_init(
    const GfsKeplerVacuumStencil *stencil, GfsKeplerVacuumRK4 *out)
{
    GfsKeplerVacuumRK4 result = {0};
    if (!stencil || !out || out->sites || out->count || out->active ||
        out->next_stage || !stencil->sites || stencil->n_sites == 0 ||
        stencil->n_sites > SIZE_MAX / sizeof(*result.sites))
        return 0;

    result.sites = (GfsKeplerVacuumRK4Site *)calloc(
        stencil->n_sites, sizeof(*result.sites));
    if (!result.sites) return 0;
    result.count = stencil->n_sites;
    for (size_t s = 0; s < result.count; ++s) {
        const GfsKeplerVacuumSite *source = &stencil->sites[s];
        if (source->stable_id < 0 || !isfinite(source->x) ||
            !isfinite(source->y)) {
            gfs_kepler_vacuum_rk4_destroy(&result);
            return 0;
        }
        result.sites[s].stable_id = source->stable_id;
        result.sites[s].x = source->x;
        result.sites[s].y = source->y;
    }
    *out = result;
    return 1;
}

/* Begin a trial. Rejected trials never alter accepted coordinates, so
 * rollback only drops the stage increments; a retry can use a smaller dt. */
static inline int gfs_kepler_vacuum_rk4_begin(GfsKeplerVacuumRK4 *state)
{
    if (!state || !state->sites || state->count == 0 || state->active)
        return 0;
    for (size_t s = 0; s < state->count; ++s) {
        if (!isfinite(state->sites[s].x) || !isfinite(state->sites[s].y))
            return 0;
    }
    for (size_t s = 0; s < state->count; ++s) {
        for (int k = 0; k < 4; ++k) {
            state->sites[s].kx[k] = 0.0;
            state->sites[s].ky[k] = 0.0;
        }
    }
    state->active = 1;
    state->next_stage = 0;
    return 1;
}

/* Return the generator positions at the current RK stage (0..3):
 * x_n, x_n+k1/2, x_n+k2/2, x_n+k3. */
static inline int gfs_kepler_vacuum_rk4_stage_positions(
    const GfsKeplerVacuumRK4 *state, int stage, double *x, double *y)
{
    if (!state || !state->active || stage != state->next_stage ||
        stage < 0 || stage > 3 || !x || !y)
        return 0;
    for (size_t s = 0; s < state->count; ++s) {
        const GfsKeplerVacuumRK4Site *site = &state->sites[s];
        double dx = 0.0, dy = 0.0;
        if (stage == 1) {
            dx = 0.5 * site->kx[0];
            dy = 0.5 * site->ky[0];
        } else if (stage == 2) {
            dx = 0.5 * site->kx[1];
            dy = 0.5 * site->ky[1];
        } else if (stage == 3) {
            dx = site->kx[2];
            dy = site->ky[2];
        }
        x[s] = site->x + dx;
        y[s] = site->y + dy;
        if (!isfinite(x[s]) || !isfinite(y[s])) return 0;
    }
    return 1;
}

/* Record dx/dt,dy/dt for the stage as increments over `dt`. All entries are
 * validated before any are written, preserving the state on bad input. */
static inline int gfs_kepler_vacuum_rk4_record_rates(
    GfsKeplerVacuumRK4 *state, int stage, double dt,
    const double *vx, const double *vy)
{
    if (!state || !state->active || stage != state->next_stage ||
        stage < 0 || stage > 3 || !(dt > 0.0) || !isfinite(dt) ||
        !vx || !vy)
        return 0;
    for (size_t s = 0; s < state->count; ++s) {
        if (!isfinite(vx[s]) || !isfinite(vy[s]) ||
            !isfinite(vx[s] * dt) || !isfinite(vy[s] * dt))
            return 0;
    }
    for (size_t s = 0; s < state->count; ++s) {
        state->sites[s].kx[stage] = vx[s] * dt;
        state->sites[s].ky[stage] = vy[s] * dt;
    }
    ++state->next_stage;
    return 1;
}

/* Commit only after all four derivatives have been recorded. */
static inline int gfs_kepler_vacuum_rk4_accept(GfsKeplerVacuumRK4 *state)
{
    if (!state || !state->active || state->next_stage != 4) return 0;
    for (size_t s = 0; s < state->count; ++s) {
        const GfsKeplerVacuumRK4Site *site = &state->sites[s];
        double dx = (site->kx[0] + 2.0 * site->kx[1] +
                     2.0 * site->kx[2] + site->kx[3]) / 6.0;
        double dy = (site->ky[0] + 2.0 * site->ky[1] +
                     2.0 * site->ky[2] + site->ky[3]) / 6.0;
        if (!isfinite(site->x + dx) || !isfinite(site->y + dy)) return 0;
    }
    for (size_t s = 0; s < state->count; ++s) {
        GfsKeplerVacuumRK4Site *site = &state->sites[s];
        site->x += (site->kx[0] + 2.0 * site->kx[1] +
                    2.0 * site->kx[2] + site->kx[3]) / 6.0;
        site->y += (site->ky[0] + 2.0 * site->ky[1] +
                    2.0 * site->ky[2] + site->ky[3]) / 6.0;
    }
    state->active = 0;
    state->next_stage = 0;
    return 1;
}

static inline int gfs_kepler_vacuum_rk4_rollback(
    GfsKeplerVacuumRK4 *state)
{
    if (!state || !state->active) return 0;
    state->active = 0;
    state->next_stage = 0;
    return 1;
}

#endif
