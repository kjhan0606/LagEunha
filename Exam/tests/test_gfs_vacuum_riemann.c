#include "../gfs_vacuum_riemann.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static int close_rel(double a, double b, double tol)
{
    return fabs(a - b) <= tol * fmax(1.0, fmax(fabs(a), fabs(b)));
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
}

int main(void)
{
    test_fan_edges_and_entropy();
    test_shifted_state_and_ale_flux();
    test_reject_invalid_input();
    puts("gfs_vacuum_riemann: PASS");
    return 0;
}
