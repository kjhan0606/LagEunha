#ifndef GFS_KEPLER_KDK_H
#define GFS_KEPLER_KDK_H

/* Per-particle KDK subcycling for the GM=1 softened central Kepler force.
 * Hydrodynamic kicks remain a separate, globally synchronized operator. */
#include <math.h>
#include "gfs_timebins.h"

static inline int gfs_kepler_kdk_accel(double x, double y,
                                       double cx, double cy, double eps,
                                       double *ax, double *ay)
{
    if (!ax || !ay || !isfinite(x) || !isfinite(y) ||
        !isfinite(cx) || !isfinite(cy) || !isfinite(eps) || eps < 0.0)
        return 0;
    double dx = x - cx, dy = y - cy;
    double r2 = dx*dx + dy*dy;
    if (eps > 0.0) {
        double R2 = r2 + eps*eps;
        double f = -1.0 / (R2 * sqrt(R2));
        *ax = dx * f; *ay = dy * f;
    } else if (r2 > 1e-16) {
        double invr = 1.0 / sqrt(r2);
        double f = -invr * invr * invr;
        *ax = dx * f; *ay = dy * f;
    } else {
        *ax = 0.0; *ay = 0.0;
    }
    return isfinite(*ax) && isfinite(*ay);
}

/* One common macro interval, with a dyadic local KDK step selected from the
 * softened orbital time eta*(r^2+eps^2)^(3/4) at the interval start. */
static inline int gfs_kepler_kdk_advance(double *x, double *y,
                                         double *vx, double *vy,
                                         double cx, double cy, double eps,
                                         double eta, double macro_dt,
                                         unsigned max_level,
                                         unsigned *level_used)
{
    if (!x || !y || !vx || !vy || !level_used ||
        !isfinite(*x) || !isfinite(*y) || !isfinite(*vx) || !isfinite(*vy) ||
        !isfinite(cx) || !isfinite(cy) || !isfinite(eps) || eps < 0.0 ||
        !isfinite(eta) || !(eta > 0.0) ||
        !isfinite(macro_dt) || !(macro_dt > 0.0))
        return GFS_TIMEBINS_INVALID;
    double dx = *x - cx, dy = *y - cy;
    double r2soft = dx*dx + dy*dy + eps*eps;
    if (!isfinite(r2soft) || !(r2soft > 0.0)) return GFS_TIMEBINS_INVALID;
    double dt_limit = eta * pow(r2soft, 0.75);
    unsigned level = 0;
    int rc = gfs_timebins_level_for_limit(dt_limit, macro_dt,
                                           max_level, &level);
    if (rc != GFS_TIMEBINS_OK) return rc;

    unsigned nsub = 1u << level;
    double h = macro_dt / (double)nsub;
    for (unsigned k = 0; k < nsub; ++k) {
        double ax, ay;
        if (!gfs_kepler_kdk_accel(*x, *y, cx, cy, eps, &ax, &ay))
            return GFS_TIMEBINS_INVALID;
        *vx += 0.5 * h * ax;
        *vy += 0.5 * h * ay;
        *x += h * *vx;
        *y += h * *vy;
        if (!gfs_kepler_kdk_accel(*x, *y, cx, cy, eps, &ax, &ay))
            return GFS_TIMEBINS_INVALID;
        *vx += 0.5 * h * ax;
        *vy += 0.5 * h * ay;
    }
    *level_used = level;
    return GFS_TIMEBINS_OK;
}

#endif
