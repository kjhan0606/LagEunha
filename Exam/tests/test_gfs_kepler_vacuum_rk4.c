#include "../gfs_kepler_vacuum_rk4.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

int main(void)
{
    GfsKeplerVacuumStencil stencil = {0};
    GfsKeplerVacuumRK4 state = {0};
    GfsKeplerVacuumFaceTrace *face_trace = NULL;
    double *x = NULL, *y = NULL, *x0 = NULL, *y0 = NULL;
    double *vx = NULL, *vy = NULL, *edge_speed = NULL;
    const double dt = 0.125;
    const double gamma = 5.0 / 3.0;

    assert(gfs_kepler_vacuum_stencil_build_periodic(
        320, 5.0, 0.5, 2.0, 1, 2.5, 2.5, &stencil));
    assert(gfs_kepler_vacuum_rk4_init(&stencil, &state));
    assert(state.count == 908);

    x = (double *)malloc(state.count * sizeof(*x));
    y = (double *)malloc(state.count * sizeof(*y));
    x0 = (double *)malloc(state.count * sizeof(*x0));
    y0 = (double *)malloc(state.count * sizeof(*y0));
    vx = (double *)malloc(state.count * sizeof(*vx));
    vy = (double *)malloc(state.count * sizeof(*vy));
    face_trace = (GfsKeplerVacuumFaceTrace *)malloc(
        4 * state.count * sizeof(*face_trace));
    edge_speed = (double *)malloc(4 * state.count * sizeof(*edge_speed));
    assert(x && y && x0 && y0 && vx && vy && face_trace && edge_speed);
    for (size_t s = 0; s < state.count; ++s) {
        x0[s] = state.sites[s].x;
        y0[s] = state.sites[s].y;
        for (size_t q = 0; q < stencil.sites[s].incidence_count; ++q) {
            GfsKeplerVacuumFaceTrace *t = &face_trace[4 * s + q];
            t->rho = 1.0 + 0.01 * (double)s + 0.03 * (double)q;
            t->pressure = 0.2 + 0.001 * (double)s + 0.01 * (double)q;
            t->vx = -0.4 + 0.0007 * (double)s + 0.05 * (double)q;
            t->vy = 0.3 - 0.0003 * (double)s - 0.02 * (double)q;
        }
    }

    /* Face-local states produce one exact fan-edge speed per incidence. A
     * shared two-face site must satisfy both independent normal constraints. */
    assert(gfs_kepler_vacuum_rk4_boundary_velocities(
        &stencil, face_trace, gamma, vx, vy, edge_speed));
    for (size_t s = 0; s < state.count; ++s) {
        const GfsKeplerVacuumSite *site = &stencil.sites[s];
        const GfsKeplerVacuumIncidence *inc = &stencil.incidences[4 * s];
        for (size_t q = 0; q < site->incidence_count; ++q) {
            const GfsKeplerVacuumFaceTrace *t = &face_trace[4 * s + q];
            double nmag = hypot((double)inc[q].nx, (double)inc[q].ny);
            double nx = (double)inc[q].nx / nmag;
            double ny = (double)inc[q].ny / nmag;
            double un = t->vx * nx + t->vy * ny;
            double c = sqrt(gamma * t->pressure / t->rho);
            double expected_edge = un + 2.0 * c / (gamma - 1.0);
            double mirror_normal = vx[s] * nx + vy[s] * ny;
            assert(fabs(edge_speed[4 * s + q] - expected_edge) < 1e-12);
            assert(fabs(mirror_normal - (2.0 * expected_edge - un)) < 1e-12);
        }
    }
    for (size_t s = 0; s < state.count; ++s) {
        if (stencil.sites[s].incidence_count == 2) {
            double saved_pressure = face_trace[4 * s + 1].pressure;
            face_trace[4 * s + 1].pressure = -1.0;
            assert(!gfs_kepler_vacuum_rk4_boundary_velocities(
                &stencil, face_trace, gamma, vx, vy, edge_speed));
            face_trace[4 * s + 1].pressure = saved_pressure;
            break;
        }
    }

    /* Exercise the off-center geometric face-speed equations used by the
     * solver. Each synthetic target velocity determines the exact required
     * face speed; one-face sites choose a minimum-change velocity, while
     * two-face sites recover both independent constraints. */
    assert(gfs_kepler_vacuum_rk4_begin(&state));
    int checked_incomplete_topology = 0;
    for (int stage = 0; stage < 4; ++stage) {
        assert(gfs_kepler_vacuum_rk4_stage_positions(&state, stage, x, y));
        gfs_kepler_vacuum_rk4_begin_force(&state);
        for (size_t s = 0; s < state.count; ++s) {
            const GfsKeplerVacuumSite *site = &stencil.sites[s];
            const GfsKeplerVacuumIncidence *inc =
                &stencil.incidences[4 * s];
            const double target_vx = 0.37 + 0.0001 * (double)s;
            const double target_vy = -0.29 + 0.0002 * (double)s;
            for (size_t q = 0; q < site->incidence_count; ++q) {
                const size_t guarded_id = (size_t)inc[q].fluid_id;
                const size_t ci = guarded_id % (size_t)stencil.computational_n;
                const size_t cj = guarded_id / (size_t)stencil.computational_n;
                const size_t fluid_id =
                    (cj - (size_t)stencil.guard) *
                        (size_t)stencil.target_n +
                    ci - (size_t)stencil.guard;
                size_t mapped_slot = SIZE_MAX;
                assert(gfs_kepler_vacuum_stencil_incidence_slot_for_particle(
                    &stencil, s, fluid_id, &mapped_slot));
                assert(mapped_slot == q);
                const double nx = (double)inc[q].nx;
                const double ny = (double)inc[q].ny;
                const double xj[2] = {x[s], y[s]};
                const double xi[2] = {x[s] - nx * stencil.dx,
                                      y[s] - ny * stencil.dx};
                const double tangent = (q == 0 ? 0.17 : -0.13) * stencil.dx;
                const double face[2] = {
                    0.5 * (xi[0] + xj[0]) - ny * tangent,
                    0.5 * (xi[1] + xj[1]) + nx * tangent
                };
                const double vi[2] = {
                    -0.21 + 0.00003 * (double)s + 0.02 * (double)q,
                    0.16 - 0.00002 * (double)s - 0.01 * (double)q
                };
                const double distance = hypot(xj[0] - xi[0],
                                              xj[1] - xi[1]);
                const double w = ((face[0] - xi[0]) * vi[0] +
                                  (face[1] - xi[1]) * vi[1] -
                                  (face[0] - xj[0]) * target_vx -
                                  (face[1] - xj[1]) * target_vy) / distance;
                edge_speed[4 * s + q] = w;
                assert(gfs_kepler_vacuum_rk4_add_face(
                    &state, s, q, xi, xj, face, vi, w));
            }
        }
        if (!checked_incomplete_topology) {
            for (size_t s = 0; s < state.count; ++s) {
                if (state.sites[s].expected_faces == 2) {
                    --state.sites[s].nface;
                    assert(!gfs_kepler_vacuum_rk4_solve_faces(&state));
                    ++state.sites[s].nface;
                    checked_incomplete_topology = 1;
                    break;
                }
            }
            assert(checked_incomplete_topology);
        }
        assert(gfs_kepler_vacuum_rk4_solve_faces(&state));
        for (size_t s = 0; s < state.count; ++s) {
            const GfsKeplerVacuumSite *site = &stencil.sites[s];
            const GfsKeplerVacuumIncidence *inc =
                &stencil.incidences[4 * s];
            const double target_vx = 0.37 + 0.0001 * (double)s;
            const double target_vy = -0.29 + 0.0002 * (double)s;
            for (size_t q = 0; q < site->incidence_count; ++q) {
                const double nx = (double)inc[q].nx;
                const double ny = (double)inc[q].ny;
                const double xj[2] = {x[s], y[s]};
                const double xi[2] = {x[s] - nx * stencil.dx,
                                      y[s] - ny * stencil.dx};
                const double tangent = (q == 0 ? 0.17 : -0.13) * stencil.dx;
                const double face[2] = {
                    0.5 * (xi[0] + xj[0]) - ny * tangent,
                    0.5 * (xi[1] + xj[1]) + nx * tangent
                };
                const double vi[2] = {
                    -0.21 + 0.00003 * (double)s + 0.02 * (double)q,
                    0.16 - 0.00002 * (double)s - 0.01 * (double)q
                };
                const double distance = hypot(xj[0] - xi[0],
                                              xj[1] - xi[1]);
                const double calculated =
                    ((face[0] - xi[0]) * vi[0] +
                     (face[1] - xi[1]) * vi[1] -
                     (face[0] - xj[0]) * state.sites[s].vx -
                     (face[1] - xj[1]) * state.sites[s].vy) / distance;
                assert(fabs(calculated - edge_speed[4 * s + q]) < 1e-11);
            }
            if (site->incidence_count == 2) {
                assert(fabs(state.sites[s].vx - target_vx) < 1e-11);
                assert(fabs(state.sites[s].vy - target_vy) < 1e-11);
            }
        }
        assert(gfs_kepler_vacuum_rk4_record_solved_rates(
            &state, stage, dt));
    }
    assert(gfs_kepler_vacuum_rk4_rollback(&state));

    /* A repeated incidence slot is rejected even if the site would still
     * receive the expected number of faces. */
    {
        GfsKeplerVacuumRK4 duplicate = {0};
        assert(gfs_kepler_vacuum_rk4_init(&stencil, &duplicate));
        assert(gfs_kepler_vacuum_rk4_begin(&duplicate));
        gfs_kepler_vacuum_rk4_begin_force(&duplicate);
        int checked_duplicate = 0;
        for (size_t s = 0; s < stencil.n_sites; ++s) {
            if (stencil.sites[s].incidence_count != 2) continue;
            const GfsKeplerVacuumIncidence *inc =
                &stencil.incidences[4 * s];
            const double xj[2] = {
                duplicate.sites[s].x, duplicate.sites[s].y};
            const double xi[2] = {
                xj[0] - (double)inc[0].nx * stencil.dx,
                xj[1] - (double)inc[0].ny * stencil.dx};
            const double face[2] = {
                0.5 * (xi[0] + xj[0]), 0.5 * (xi[1] + xj[1])};
            const double vi[2] = {0.0, 0.0};
            assert(gfs_kepler_vacuum_rk4_add_face(
                &duplicate, s, 0, xi, xj, face, vi, 0.0));
            assert(!gfs_kepler_vacuum_rk4_add_face(
                &duplicate, s, 0, xi, xj, face, vi, 0.0));
            checked_duplicate = 1;
            break;
        }
        assert(checked_duplicate);
        gfs_kepler_vacuum_rk4_destroy(&duplicate);
    }
    {
        size_t slot = 0;
        assert(!gfs_kepler_vacuum_stencil_incidence_slot_for_particle(
            &stencil, 0, SIZE_MAX, &slot));
    }
    gfs_kepler_vacuum_stencil_destroy(&stencil);

    /* Integrate dx/dt=x, dy/dt=3 with the same stage abscissae as the fluid
     * RK4 loop. The RK4 polynomial is exact through the fourth-order term. */
    assert(gfs_kepler_vacuum_rk4_begin(&state));
    for (int stage = 0; stage < 4; ++stage) {
        assert(gfs_kepler_vacuum_rk4_stage_positions(&state, stage, x, y));
        for (size_t s = 0; s < state.count; ++s) {
            vx[s] = x[s];
            vy[s] = 3.0;
        }
        assert(gfs_kepler_vacuum_rk4_record_rates(
            &state, stage, dt, vx, vy));
    }
    /* The fluid RK4 path undoes the K4 input shift before its final blend.
     * Mirror sites must return to x_n for that gap without losing K1..K4. */
    assert(state.next_stage == 4 && state.active);
    assert(gfs_kepler_vacuum_rk4_stage_positions(&state, 0, x, y));
    assert(state.next_stage == 4 && state.active);
    for (size_t s = 0; s < state.count; ++s) {
        assert(fabs(x[s] - x0[s]) < 1e-14);
        assert(fabs(y[s] - y0[s]) < 1e-14);
    }
    assert(gfs_kepler_vacuum_rk4_accept(&state));
    for (size_t s = 0; s < state.count; ++s) {
        double growth = 1.0 + dt + 0.5 * dt * dt +
                        dt * dt * dt / 6.0 + dt * dt * dt * dt / 24.0;
        assert(fabs(state.sites[s].x - x0[s] * growth) < 1e-13);
        assert(fabs(state.sites[s].y - (y0[s] + 3.0 * dt)) < 1e-13);
        assert(state.sites[s].stable_id >= 0);
    }

    /* A rejected step leaves accepted positions unchanged and can be retried. */
    assert(gfs_kepler_vacuum_rk4_begin(&state));
    assert(gfs_kepler_vacuum_rk4_stage_positions(&state, 0, x, y));
    for (size_t s = 0; s < state.count; ++s) {
        vx[s] = 100.0;
        vy[s] = -100.0;
    }
    vx[17] = NAN;
    assert(!gfs_kepler_vacuum_rk4_record_rates(&state, 0, dt, vx, vy));
    assert(state.next_stage == 0);
    assert(gfs_kepler_vacuum_rk4_rollback(&state));
    assert(!state.active && state.next_stage == 0);
    for (size_t s = 0; s < state.count; ++s) {
        assert(fabs(state.sites[s].x - x[s]) < 1e-14);
        assert(fabs(state.sites[s].y - y[s]) < 1e-14);
    }

    assert(gfs_kepler_vacuum_rk4_begin(&state));
    for (int stage = 0; stage < 4; ++stage) {
        assert(gfs_kepler_vacuum_rk4_stage_positions(&state, stage, x, y));
        for (size_t s = 0; s < state.count; ++s) {
            vx[s] = 0.0;
            vy[s] = 1.0;
        }
        assert(gfs_kepler_vacuum_rk4_record_rates(
            &state, stage, 0.5 * dt, vx, vy));
    }
    assert(gfs_kepler_vacuum_rk4_accept(&state));
    for (size_t s = 0; s < state.count; ++s)
        assert(fabs(state.sites[s].y - (y0[s] + 3.0 * dt + 0.5 * dt)) < 1e-13);

    free(x);
    free(y);
    free(x0);
    free(y0);
    free(vx);
    free(vy);
    free(face_trace);
    free(edge_speed);
    gfs_kepler_vacuum_rk4_destroy(&state);
    assert(state.sites == NULL && state.count == 0);
    return 0;
}
