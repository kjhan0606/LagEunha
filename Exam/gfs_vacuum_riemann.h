#ifndef GFS_VACUUM_RIEMANN_H
#define GFS_VACUUM_RIEMANN_H

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

/* Solve the isentropic rarefaction fan for gas on the left and vacuum on the
 * right.  The returned state is the exact Euler solution at similarity speed
 * xi.  A return value of zero means the input state is not physical. */
GFS_VAC_HD int gfs_vacuum_fan_state(
    double rho0, double pressure0, double velocity0, double gamma,
    double xi, GfsVacuumFanState *state)
{
    if (!state || !(rho0 > 0.0) || !(pressure0 > 0.0) || !(gamma > 1.0))
        return 0;

    const double c0 = sqrt(gamma * pressure0 / rho0);
    if (!(c0 > 0.0))
        return 0;

    const double head = velocity0 - c0;
    const double tail = velocity0 + 2.0 * c0 / (gamma - 1.0);
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

/* Euler flux through a face moving at face_speed, in the face-normal frame.
 * tangential_velocity is retained in total energy for a multidimensional
 * application.  At the physical vacuum edge all three fluxes are zero. */
GFS_VAC_HD int gfs_vacuum_fan_ale_flux(
    double rho0, double pressure0, double velocity0, double gamma,
    double face_speed, double tangential_velocity,
    double *mass_flux, double *momentum_flux, double *energy_flux)
{
    if (!mass_flux || !momentum_flux || !energy_flux)
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
