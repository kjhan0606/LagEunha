#include "gfs_kepler_kdk.h"
#include <stdio.h>

#define CHECK(c) do { if (!(c)) { \
    fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); \
    return 1; \
} } while (0)

static int circular_orbit_metrics(double eta, double *position_error,
                                  double *energy_error, double *angular_error,
                                  unsigned *max_level_used)
{
    const double r = 1.0, eps = 0.01175;
    const double omega = pow(r*r + eps*eps, -0.75);
    const double period = 2.0 * acos(-1.0) / omega;
    const double macro_dt = period / 32.0;
    double x = r, y = 0.0, vx = 0.0, vy = r * omega;
    const double E0 = 0.5 * vy*vy - 1.0/sqrt(r*r + eps*eps);
    const double L0 = x*vy - y*vx;
    *max_level_used = 0;

    /* 100 circular periods at common macro endpoints; gravity is subcycled
     * independently with dyadic KDK intervals inside each macro interval. */
    for (int orbit = 0; orbit < 100; ++orbit) {
        for (int step = 0; step < 32; ++step) {
            unsigned level = 0;
            int rc=gfs_kepler_kdk_advance(&x, &y, &vx, &vy, 0.0, 0.0,
                    eps, eta, macro_dt, 12, &level);
            if (rc != GFS_TIMEBINS_OK) return 0;
            if (level > *max_level_used) *max_level_used = level;
        }
    }

    double E = 0.5 * (vx*vx + vy*vy) - 1.0/sqrt(x*x + y*y + eps*eps);
    double L = x*vy - y*vx;
    *position_error = hypot(x-r, y) / r;
    *energy_error = fabs((E-E0)/E0);
    *angular_error = fabs((L-L0)/L0);
    return 1;
}

int main(void)
{
    double x_coarse, e_coarse, l_coarse, x_fine, e_fine, l_fine;
    unsigned b_coarse=0, b_fine=0;
    CHECK(circular_orbit_metrics(0.01,&x_coarse,&e_coarse,&l_coarse,&b_coarse));
    CHECK(circular_orbit_metrics(0.005,&x_fine,&e_fine,&l_fine,&b_fine));
    CHECK(b_coarse == 5 && b_fine == 6);
    CHECK(x_fine < x_coarse);
    CHECK(x_coarse/x_fine > 3.5 && x_coarse/x_fine < 4.5);
    CHECK(x_fine < 0.005);
    CHECK(e_coarse < 2e-5 && e_fine < 2e-5);
    CHECK(l_coarse < 1e-11 && l_fine < 1e-11);

    unsigned unchanged = 77;
    CHECK(gfs_timebins_level_for_limit(0.001, 1.0, 4, &unchanged)
          == GFS_TIMEBINS_TOO_FINE);
    CHECK(unchanged == 77);
    printf("PASS local KDK: 100 Plummer circular periods, eta=.01/.005 rel_x=%.3e/%.3e ratio=%.3f rel_E=%.3e rel_L=%.3e bins=%u/%u\n",
           x_coarse, x_fine, x_coarse/x_fine, e_fine, l_fine, b_coarse, b_fine);
    return 0;
}
