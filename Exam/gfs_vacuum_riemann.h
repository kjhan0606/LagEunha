#ifndef GFS_VACUUM_RIEMANN_H
#define GFS_VACUUM_RIEMANN_H

#include <float.h>
#include <math.h>

#if defined(__CUDACC__)
#define GFS_VAC_HD __host__ __device__ __forceinline__
#else
#define GFS_VAC_HD static inline
#endif

/* State in the right-going expansion of a uniform ideal gas into vacuum.
 * xi is the self-similar coordinate (x - x0) / t, or the normal speed of a
 * moving face when sampling the exact Riemann solution. */
typedef struct {
    double rho;
    double pressure;
    double velocity;
    double sound_speed;
    int vacuum;
} GfsVacuumFanState;

/* CUDA-safe finite check: NaNs fail self-equality and infinities exceed the
 * largest finite double. */
GFS_VAC_HD int gfs_vacuum_finite(double value)
{
    return value == value && value <= DBL_MAX && value >= -DBL_MAX;
}

/* Solve the isentropic rarefaction fan for gas on the left and vacuum on the
 * right.  The returned state is the exact Euler solution at similarity speed
 * xi.  A return value of zero means the input state is not physical. */
GFS_VAC_HD int gfs_vacuum_fan_state(
    double rho0, double pressure0, double velocity0, double gamma,
    double xi, GfsVacuumFanState *state)
{
    if (!state || !gfs_vacuum_finite(rho0) || !gfs_vacuum_finite(pressure0)
        || !gfs_vacuum_finite(velocity0) || !gfs_vacuum_finite(gamma)
        || !gfs_vacuum_finite(xi) || !(rho0 > 0.0)
        || !(pressure0 > 0.0) || !(gamma > 1.0))
        return 0;

    const double c0 = sqrt(gamma * pressure0 / rho0);
    if (!(c0 > 0.0) || !gfs_vacuum_finite(c0))
        return 0;

    const double head = velocity0 - c0;
    const double tail = velocity0 + 2.0 * c0 / (gamma - 1.0);
    if (!gfs_vacuum_finite(head) || !gfs_vacuum_finite(tail))
        return 0;
    if (xi <= head) {
        state->rho = rho0;
        state->pressure = pressure0;
        state->velocity = velocity0;
        state->sound_speed = c0;
        state->vacuum = 0;
        return 1;
    }
    if (xi >= tail) {
        state->rho = 0.0;
        state->pressure = 0.0;
        state->velocity = tail;
        state->sound_speed = 0.0;
        state->vacuum = 1;
        return 1;
    }

    const double c = ((gamma - 1.0) / (gamma + 1.0)) * (tail - xi);
    const double ratio = c / c0;
    state->sound_speed = c;
    state->velocity = xi + c;
    state->rho = rho0 * pow(ratio, 2.0 / (gamma - 1.0));
    state->pressure = pressure0 * pow(ratio, 2.0 * gamma / (gamma - 1.0));
    state->vacuum = 0;
    return 1;
}

/* Kinematics for an equal-weight Voronoi mirror site across a planar
 * gas-vacuum edge. The primitive state is the gas-side trace at that face;
 * using cell-center values is only a first-order approximation. n is oriented
 * from gas into vacuum and may be unnormalized.
 * The physical vacuum edge moves at u_n + 2 c/(gamma-1). Choosing the mirror
 * site's normal velocity as 2*w_edge-u_n makes the Voronoi face speed
 * (u_n+u_ghost,n)/2 exactly w_edge. Tangential velocity is copied from gas;
 * curvature and Laguerre weight-rate terms are outside this planar W=0 rule. */
GFS_VAC_HD int gfs_vacuum_mirror_kinematics(
    double rho, double pressure, double vx, double vy, double gamma,
    double nx, double ny, double *edge_speed,
    double *ghost_vx, double *ghost_vy)
{
    if (!edge_speed || !ghost_vx || !ghost_vy
        || !gfs_vacuum_finite(rho) || !gfs_vacuum_finite(pressure)
        || !gfs_vacuum_finite(vx) || !gfs_vacuum_finite(vy)
        || !gfs_vacuum_finite(gamma) || !gfs_vacuum_finite(nx)
        || !gfs_vacuum_finite(ny) || !(rho > 0.0)
        || !(pressure > 0.0) || !(gamma > 1.0))
        return 0;

    const double nmag = sqrt(nx * nx + ny * ny);
    if (!(nmag > 0.0) || !gfs_vacuum_finite(nmag))
        return 0;

    const double nxhat = nx / nmag;
    const double nyhat = ny / nmag;
    const double un = vx * nxhat + vy * nyhat;
    const double c = sqrt(gamma * pressure / rho);
    if (!(c > 0.0) || !gfs_vacuum_finite(un) || !gfs_vacuum_finite(c))
        return 0;

    const double w = un + 2.0 * c / (gamma - 1.0);
    const double delta_vn = 2.0 * (w - un);
    const double vgx = vx + delta_vn * nxhat;
    const double vgy = vy + delta_vn * nyhat;
    if (!gfs_vacuum_finite(w) || !gfs_vacuum_finite(vgx)
        || !gfs_vacuum_finite(vgy))
        return 0;
    *edge_speed = w;
    *ghost_vx = vgx;
    *ghost_vy = vgy;
    return 1;
}

/* Euler flux through a face moving at face_speed, in the face-normal frame.
 * tangential_velocity is retained in total energy for a multidimensional
 * application.  At the physical vacuum edge all three fluxes are zero. */
GFS_VAC_HD int gfs_vacuum_fan_ale_flux(
    double rho0, double pressure0, double velocity0, double gamma,
    double face_speed, double tangential_velocity,
    double *mass_flux, double *momentum_flux, double *energy_flux)
{
    if (!mass_flux || !momentum_flux || !energy_flux
        || !gfs_vacuum_finite(face_speed)
        || !gfs_vacuum_finite(tangential_velocity))
        return 0;

    GfsVacuumFanState state;
    if (!gfs_vacuum_fan_state(rho0, pressure0, velocity0, gamma,
                              face_speed, &state))
        return 0;
    if (state.vacuum) {
        *mass_flux = 0.0;
        *momentum_flux = 0.0;
        *energy_flux = 0.0;
        return 1;
    }

    const double rel_velocity = state.velocity - face_speed;
    const double total_energy = state.pressure / (gamma - 1.0)
        + 0.5 * state.rho * (state.velocity * state.velocity
                            + tangential_velocity * tangential_velocity);
    *mass_flux = state.rho * rel_velocity;
    *momentum_flux = state.rho * state.velocity * rel_velocity
                   + state.pressure;
    *energy_flux = total_energy * rel_velocity
                 + state.pressure * state.velocity;
    return 1;
}

#undef GFS_VAC_HD

#endif
