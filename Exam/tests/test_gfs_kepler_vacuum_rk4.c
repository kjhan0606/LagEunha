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
