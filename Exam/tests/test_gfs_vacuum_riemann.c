#include "../gfs_vacuum_riemann.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int close_rel(double a, double b, double tol)
{
    return fabs(a - b) <= tol * fmax(1.0, fmax(fabs(a), fabs(b)));
}

static double simpson_fan_density(double rho0, double p0, double u0,
                                  double gamma, double left, double right,
                                  int intervals)
{
    assert(intervals > 0 && intervals % 2 == 0);
    const double step = (right - left) / intervals;
    double sum = 0.0;
    GfsVacuumFanState state;
    assert(gfs_vacuum_fan_state(rho0, p0, u0, gamma, left, &state));
    sum += state.rho;
    assert(gfs_vacuum_fan_state(rho0, p0, u0, gamma, right, &state));
    sum += state.rho;
    for (int i = 1; i < intervals; ++i) {
        const double xi = left + i * step;
        assert(gfs_vacuum_fan_state(rho0, p0, u0, gamma, xi, &state));
        sum += (i % 2 ? 4.0 : 2.0) * state.rho;
    }
    return sum * step / 3.0;
}

static void test_fan_edges_and_entropy(void)
{
    const double rho0 = 1.0;
    const double p0 = 1.0 / (5.0 / 3.0);
    const double gamma = 5.0 / 3.0;
    const double u0 = 0.0;
    const double c0 = 1.0;
    const double head = u0 - c0;
    const double tail = u0 + 2.0 * c0 / (gamma - 1.0);
    GfsVacuumFanState s;

    assert(gfs_vacuum_fan_state(rho0, p0, u0, gamma, head, &s));
    assert(!s.vacuum);
    assert(close_rel(s.rho, rho0, 1e-14));
    assert(close_rel(s.pressure, p0, 1e-14));
    assert(close_rel(s.velocity, u0, 1e-14));
    assert(close_rel(s.sound_speed, c0, 1e-14));

    /* For gamma=5/3 and xi=1, c/c0=1/2 exactly. */
    assert(gfs_vacuum_fan_state(rho0, p0, u0, gamma, 1.0, &s));
    assert(!s.vacuum);
    assert(close_rel(s.sound_speed, 0.5, 1e-14));
    assert(close_rel(s.velocity, 1.5, 1e-14));
    assert(close_rel(s.rho, 0.125, 1e-14));
    assert(close_rel(s.pressure, 0.01875, 1e-14));
    assert(close_rel(s.pressure / pow(s.rho, gamma),
                     p0 / pow(rho0, gamma), 1e-13));

    assert(gfs_vacuum_fan_state(rho0, p0, u0, gamma, tail, &s));
    assert(s.vacuum && s.rho == 0.0 && s.pressure == 0.0);
    assert(gfs_vacuum_fan_state(rho0, p0, u0, gamma, tail + 1.0, &s));
    assert(s.vacuum);
}

static void test_fan_mass_integral(void)
{
    const double rho0 = 1.0;
    const double gamma = 5.0 / 3.0;
    const double p0 = 1.0 / gamma;
    const double u0 = 0.0;
    const double c0 = 1.0;
    const double head = u0 - c0;
    const double tail = u0 + 2.0 * c0 / (gamma - 1.0);
    const double fan_mass_rate = simpson_fan_density(
        rho0, p0, u0, gamma, head, tail, 10000);

    /* The fan replaces the advected initial gas interval [head,u0]. */
    assert(close_rel(fan_mass_rate, rho0 * (u0 - head), 2e-13));
}

static void test_nonzero_ale_flux_inside_fan(void)
{
    const double gamma = 5.0 / 3.0;
    const double rho0 = 1.0;
    const double c0 = 1.0;
    const double p0 = rho0 * c0 * c0 / gamma;
    double mass_flux, momentum_flux, energy_flux;

    assert(gfs_vacuum_fan_ale_flux(rho0, p0, 0.0, gamma, 0.0, 0.2,
                                   &mass_flux, &momentum_flux, &energy_flux));
    assert(close_rel(mass_flux, 0.31640625, 1e-14));
    assert(close_rel(momentum_flux, 0.3796875, 1e-14));
    assert(close_rel(energy_flux, 0.36228515625, 1e-14));
}

static void test_vacuum_mirror_kinematics(void)
{
    const double gamma = 5.0 / 3.0;
    const double rho = 1.0;
    const double pressure = 1.0 / gamma;
    const double vx = 0.4;
    const double vy = -0.2;
    double edge_speed, ghost_vx, ghost_vy;

    /* Outer edge: n points along +x and c=1. */
    assert(gfs_vacuum_mirror_kinematics(rho, pressure, vx, vy, gamma,
                                       1.0, 0.0, &edge_speed,
                                       &ghost_vx, &ghost_vy));
    assert(close_rel(edge_speed, 3.4, 1e-14));
    assert(close_rel(0.5 * (vx + ghost_vx), edge_speed, 1e-14));
    assert(close_rel(ghost_vy, vy, 1e-14));

    /* Inner edge: orient n into the hole, along -x. */
    assert(gfs_vacuum_mirror_kinematics(rho, pressure, vx, vy, gamma,
                                       -2.0, 0.0, &edge_speed,
                                       &ghost_vx, &ghost_vy));
    assert(close_rel(edge_speed, 2.6, 1e-14));
    assert(close_rel(0.5 * (vx + ghost_vx) * -1.0,
                     edge_speed, 1e-14));
    assert(close_rel(ghost_vy, vy, 1e-14));

    /* An oblique, unnormalized normal checks projection and normalization. */
    assert(gfs_vacuum_mirror_kinematics(rho, pressure, vx, vy, gamma,
                                       3.0, 4.0, &edge_speed,
                                       &ghost_vx, &ghost_vy));
    assert(close_rel(edge_speed, 3.08, 1e-14));
    assert(close_rel(0.5 * ((vx + ghost_vx) * 0.6
                          + (vy + ghost_vy) * 0.8),
                     edge_speed, 1e-14));
    assert(close_rel(-0.8 * ghost_vx + 0.6 * ghost_vy,
                     -0.8 * vx + 0.6 * vy, 1e-14));

    assert(!gfs_vacuum_mirror_kinematics(rho, pressure, vx, vy, gamma,
                                        0.0, 0.0, &edge_speed,
                                        &ghost_vx, &ghost_vy));
    assert(!gfs_vacuum_mirror_kinematics(rho, pressure, NAN, vy, gamma,
                                        1.0, 0.0, &edge_speed,
                                        &ghost_vx, &ghost_vy));
    assert(!gfs_vacuum_mirror_kinematics(rho, pressure, vx, vy, gamma,
                                        INFINITY, 0.0, &edge_speed,
                                        &ghost_vx, &ghost_vy));
}

static void test_shifted_state_and_ale_flux(void)
{
    const double gamma = 1.4;
    const double rho0 = 2.0;
    const double c0 = 0.8;
    const double p0 = rho0 * c0 * c0 / gamma;
    const double u0 = -0.3;
    const double tail = u0 + 2.0 * c0 / (gamma - 1.0);
    GfsVacuumFanState s;
    double fm, fp, fe;

    assert(gfs_vacuum_fan_state(rho0, p0, u0, gamma, u0, &s));
    assert(!s.vacuum);
    assert(close_rel(s.sound_speed, 2.0 * c0 / (gamma + 1.0), 1e-14));
    assert(close_rel(s.velocity, u0 + s.sound_speed, 1e-14));
    assert(close_rel(s.pressure / pow(s.rho, gamma),
                     p0 / pow(rho0, gamma), 1e-13));

    assert(gfs_vacuum_fan_ale_flux(rho0, p0, u0, gamma, tail, 0.25,
                                   &fm, &fp, &fe));
    assert(fm == 0.0 && fp == 0.0 && fe == 0.0);
}

static void test_reject_invalid_input(void)
{
    GfsVacuumFanState s;
    assert(!gfs_vacuum_fan_state(0.0, 1.0, 0.0, 5.0 / 3.0, 0.0, &s));
    assert(!gfs_vacuum_fan_state(1.0, 1.0, 0.0, 1.0, 0.0, &s));
    assert(!gfs_vacuum_fan_state(1.0, 1.0, 0.0, 5.0 / 3.0, NAN, &s));
    assert(!gfs_vacuum_fan_state(1.0, 1.0, INFINITY, 5.0 / 3.0, 0.0, &s));
    assert(!gfs_vacuum_fan_state(INFINITY, 1.0, 0.0, 5.0 / 3.0, 0.0, &s));

    double fm, fp, fe;
    assert(!gfs_vacuum_fan_ale_flux(1.0, 1.0, 0.0, 5.0 / 3.0,
                                    NAN, 0.0, &fm, &fp, &fe));
    assert(!gfs_vacuum_fan_ale_flux(1.0, 1.0, 0.0, 5.0 / 3.0,
                                    0.0, INFINITY, &fm, &fp, &fe));
}

int main(void)
{
    test_fan_edges_and_entropy();
    test_fan_mass_integral();
    test_nonzero_ale_flux_inside_fan();
    test_vacuum_mirror_kinematics();
    test_shifted_state_and_ale_flux();
    test_reject_invalid_input();
    puts("gfs_vacuum_riemann: PASS");
    return 0;
}
