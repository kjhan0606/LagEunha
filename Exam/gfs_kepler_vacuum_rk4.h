#ifndef GFS_KEPLER_VACUUM_RK4_H
#define GFS_KEPLER_VACUUM_RK4_H

#include "gfs_kepler_vacuum_stencil.h"
#include "gfs_vacuum_riemann.h"

typedef struct {
    double rho, pressure, vx, vy;
} GfsKeplerVacuumFaceTrace;

/* RK4 trial coordinates and face constraints for zero-mass vacuum generators. The array
 * is deliberately separate from VORO_TBP/TBPP: those buffers are rebuilt at
 * every stage and include MPI halos. */
typedef struct {
    int stable_id;
    double x, y;
    double stage_x, stage_y;
    double kx[4], ky[4];
    double row[2][2], rhs[2], distance[2], edge_speed[2], vi[2][2];
    size_t nface, expected_faces;
    unsigned char face_seen[2];
    double vx, vy;
    int rhs_valid, invalid;
} GfsKeplerVacuumRK4Site;

typedef struct {
    GfsKeplerVacuumRK4Site *sites;
    size_t count;
    size_t n_single, n_pair;
    double max_face_speed_residual;
    double max_face_speed_relative_residual;
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
            !isfinite(source->y) || source->incidence_count == 0 ||
            source->incidence_count > 2) {
            gfs_kepler_vacuum_rk4_destroy(&result);
            return 0;
        }
        result.sites[s].stable_id = source->stable_id;
        result.sites[s].x = source->x;
        result.sites[s].y = source->y;
        result.sites[s].stage_x = source->x;
        result.sites[s].stage_y = source->y;
        result.sites[s].expected_faces = source->incidence_count;
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
    for (size_t s = 0; s < state->count; ++s) {
        GfsKeplerVacuumRK4Site *site = state->sites + s;
        site->stage_x = site->x;
        site->stage_y = site->y;
    }
    return 1;
}

static inline void gfs_kepler_vacuum_rk4_begin_force(
    GfsKeplerVacuumRK4 *state)
{
    if (!state) return;
    state->n_single = state->n_pair = 0;
    state->max_face_speed_residual = 0.0;
    state->max_face_speed_relative_residual = 0.0;
    for (size_t s = 0; s < state->count; ++s) {
        state->sites[s].nface = 0;
        state->sites[s].face_seen[0] = 0;
        state->sites[s].face_seen[1] = 0;
        state->sites[s].rhs_valid = 0;
        state->sites[s].invalid = 0;
    }
}

/* Add the exact zero-weight face-speed equation at the face centroid `c`:
 *   ((c-x_i)·v_i - (c-x_j)·v_j)/|x_j-x_i| = w_edge.
 * This includes the off-center term, so shear/face rotation is represented.
 * Each site accepts one face (minimum-change velocity) or two independent
 * faces (exact 2x2 solve); missing, singular, or over-subscribed topology fails. */
static inline int gfs_kepler_vacuum_rk4_add_face(
    GfsKeplerVacuumRK4 *state, size_t site_id, size_t face_slot,
    const double xi[2], const double xj[2], const double c[2],
    const double vi[2], double w_edge)
{
    if (!state || !state->active || site_id >= state->count ||
        !xi || !xj || !c || !vi || !isfinite(w_edge)) return 0;
    GfsKeplerVacuumRK4Site *site = state->sites + site_id;
    const double dx = xj[0]-xi[0], dy = xj[1]-xi[1];
    const double dist = hypot(dx,dy);
    const double qx = c[0]-xj[0], qy = c[1]-xj[1];
    const double ax = c[0]-xi[0], ay = c[1]-xi[1];
    if (!(dist > 0.0) || !(qx*qx+qy*qy > 0.0) ||
        !isfinite(dist) || !isfinite(qx) || !isfinite(qy) ||
        !isfinite(ax) || !isfinite(ay) ||
        !isfinite(vi[0]) || !isfinite(vi[1]) ||
        face_slot >= site->expected_faces || face_slot >= 2 ||
        site->face_seen[face_slot]) {
        site->invalid = 1;
        return 0;
    }
    const size_t f = face_slot;
    site->row[f][0] = qx;
    site->row[f][1] = qy;
    site->rhs[f] = ax*vi[0] + ay*vi[1] - dist*w_edge;
    site->distance[f] = dist;
    site->edge_speed[f] = w_edge;
    site->vi[f][0] = vi[0];
    site->vi[f][1] = vi[1];
    site->face_seen[f] = 1;
    ++site->nface;
    if (!isfinite(site->rhs[f])) {
        site->invalid = 1;
        return 0;
    }
    return 1;
}

static inline int gfs_kepler_vacuum_rk4_solve_faces(
    GfsKeplerVacuumRK4 *state)
{
    if (!state || !state->active) return 0;
    state->n_single = state->n_pair = 0;
    state->max_face_speed_residual = 0.0;
    state->max_face_speed_relative_residual = 0.0;
    for (size_t s = 0; s < state->count; ++s) {
        GfsKeplerVacuumRK4Site *site = state->sites+s;
        if (site->invalid || site->expected_faces == 0 ||
            site->nface != site->expected_faces) return 0;
        for (size_t f = 0; f < site->expected_faces; ++f)
            if (!site->face_seen[f]) return 0;
        if (site->nface == 1) {
            state->n_single++;
            const double qx=site->row[0][0], qy=site->row[0][1];
            const double q2=qx*qx+qy*qy;
            const double residual=site->rhs[0]-qx*site->vi[0][0]-qy*site->vi[0][1];
            if (!(q2>0.0) || !isfinite(q2)) return 0;
            site->vx=site->vi[0][0]+residual*qx/q2;
            site->vy=site->vi[0][1]+residual*qy/q2;
        } else {
            state->n_pair++;
            const double a=site->row[0][0], b=site->row[0][1];
            const double c=site->row[1][0], d=site->row[1][1];
            const double det=a*d-b*c;
            const double scale=hypot(a,b)*hypot(c,d);
            if (!(scale>0.0) || fabs(det)<=1.0e-12*scale || !isfinite(det))
                return 0;
            site->vx=(site->rhs[0]*d-b*site->rhs[1])/det;
            site->vy=(a*site->rhs[1]-site->rhs[0]*c)/det;
        }
        if (!isfinite(site->vx) || !isfinite(site->vy)) return 0;
        for (size_t f=0; f<site->nface; ++f) {
            const double residual=fabs(
                (site->row[f][0]*site->vx+site->row[f][1]*site->vy-site->rhs[f]) /
                site->distance[f]);
            const double relative=residual/fmax(fabs(site->edge_speed[f]),1.0e-12);
            if (!isfinite(residual) || !isfinite(relative)) return 0;
            if (residual>state->max_face_speed_residual)
                state->max_face_speed_residual=residual;
            if (relative>state->max_face_speed_relative_residual)
                state->max_face_speed_relative_residual=relative;
        }
        site->rhs_valid=1;
    }
    return 1;
}

/* Return the generator positions at the current RK stage (0..3):
 * x_n, x_n+k1/2, x_n+k2/2, x_n+k3. After K4 is recorded, stage 0 restores
 * the coordinates to x_n before the fluid solver's final weighted blend;
 * next_stage remains 4 so the complete trial can still be accepted. */
static inline int gfs_kepler_vacuum_rk4_stage_positions(
    GfsKeplerVacuumRK4 *state, int stage, double *x, double *y)
{
    const int k4_coordinate_rewind =
        stage == 0 && state && state->active && state->next_stage == 4;
    if (!state || !state->active ||
        (stage != state->next_stage && !k4_coordinate_rewind) ||
        stage < 0 || stage > 3 || ((x == NULL) != (y == NULL)))
        return 0;
    for (size_t s = 0; s < state->count; ++s) {
        GfsKeplerVacuumRK4Site *site = &state->sites[s];
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
        site->stage_x = site->x + dx;
        site->stage_y = site->y + dy;
        if (!isfinite(site->stage_x) || !isfinite(site->stage_y)) return 0;
        if (x) { x[s] = site->stage_x; y[s] = site->stage_y; }
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

static inline int gfs_kepler_vacuum_rk4_record_solved_rates(
    GfsKeplerVacuumRK4 *state, int stage, double dt)
{
    if (!state || !state->active || stage != state->next_stage ||
        stage < 0 || stage > 3 || !(dt > 0.0) || !isfinite(dt)) return 0;
    for (size_t s=0; s<state->count; ++s) {
        GfsKeplerVacuumRK4Site *site=state->sites+s;
        if (!site->rhs_valid || !isfinite(site->vx*dt) || !isfinite(site->vy*dt))
            return 0;
    }
    for (size_t s=0; s<state->count; ++s) {
        GfsKeplerVacuumRK4Site *site=state->sites+s;
        site->kx[stage]=site->vx*dt;
        site->ky[stage]=site->vy*dt;
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
        site->stage_x=site->x;
        site->stage_y=site->y;
    }
    state->active = 0;
    state->next_stage = 0;
    return 1;
}

static inline int gfs_kepler_vacuum_rk4_rollback(
    GfsKeplerVacuumRK4 *state)
{
    if (!state || !state->active) return 0;
    for (size_t s=0; s<state->count; ++s) {
        state->sites[s].stage_x=state->sites[s].x;
        state->sites[s].stage_y=state->sites[s].y;
    }
    state->active = 0;
    state->next_stage = 0;
    return 1;
}

#endif
