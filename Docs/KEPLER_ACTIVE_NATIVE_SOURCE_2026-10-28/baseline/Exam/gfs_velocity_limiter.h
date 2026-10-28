#ifndef GFS_VELOCITY_LIMITER_H
#define GFS_VELOCITY_LIMITER_H

#include <float.h>
#include <math.h>

#ifdef __CUDACC__
#define GFS_VEL_HD __host__ __device__
#else
#define GFS_VEL_HD
#endif

/* Limit the velocity tensor using generator-midpoint increments. The
 * admissible velocity region is the intersection of sample projection
 * strips, including their perpendiculars. One common alpha preserves
 * coordinate rotation and boost invariance. For an affine field,
 * B*d_j/2 = (v_j-v_i)/2 lies in every strip, hence alpha=1.
 *
 * This is a reconstruction limiter, not a proof of an Euler maximum
 * principle at a skew face centroid. Density/P positivity limiting and
 * shock/conservation gates remain necessary. Arrays include every real
 * face neighbor; no fixed-capacity truncation is permitted. */
GFS_VEL_HD static inline double gfs_velocity_matrix_condition(
    double xx, double xy, double yy)
{
    double scale = fmax(fabs(xx), fmax(fabs(xy), fabs(yy)));
    if (!(scale > 0.0) || !isfinite(scale)) return INFINITY;
    xx /= scale; xy /= scale; yy /= scale;
    double det = xx*yy - xy*xy;
    if (!(xx > 0.0) || !(yy > 0.0) || !(det > 0.0)) return INFINITY;
    double largest = 0.5*(xx+yy+hypot(xx-yy,2.0*xy));
    return largest*largest/det;
}

/* A negative return reports an invalid/insufficiently resolved fit; callers
 * must report/stop, never disguise that condition as a zero gradient. */
GFS_VEL_HD static inline double gfs_velocity_midpoint_alpha(
    int n, const double *dx, const double *dy,
    const double *dvx, const double *dvy,
    double bxx, double bxy, double byx, double byy,
    double condition, double input_velocity_scale)
{
    if (!isfinite(bxx) || !isfinite(bxy) || !isfinite(byx) || !isfinite(byy) ||
        !isfinite(condition) || condition*DBL_EPSILON > 1.e-6 || n <= 0)
        return -1.0;
    double scale2 = 0.0, term_scale = 0.0, prediction2 = 0.0, alpha = 1.0;
    for (int j = 0; j < n; j++) {
        double s2 = dvx[j]*dvx[j] + dvy[j]*dvy[j];
        if (!isfinite(s2) || !isfinite(dx[j]) || !isfinite(dy[j])) return -1.0;
        if (s2 > scale2) scale2 = s2;
        double tx = 0.5*(fabs(bxx*dx[j])+fabs(bxy*dy[j]));
        double ty = 0.5*(fabs(byx*dx[j])+fabs(byy*dy[j]));
        double ts = hypot(tx,ty);
        if (ts > term_scale) term_scale = ts;
        double px = 0.5*(bxx*dx[j]+bxy*dy[j]);
        double py = 0.5*(byx*dx[j]+byy*dy[j]);
        double ps = px*px+py*py;
        if (!isfinite(ps)) return -1.0;
        if (ps > prediction2) prediction2 = ps;
    }
    if (!(scale2 > 0.0)) return prediction2 == 0.0 ? 1.0 : -1.0;

    /* Dot products with collinear velocity data can have zero-width bounds.
     * Bound errors using input subtraction, the actual LSF matrix condition,
     * and absolute products in B*d, rather than the interval width. */
    double vector_error = DBL_EPSILON *
        (32.0*(1.0+condition)*sqrt(scale2) +
         8.0*fabs(input_velocity_scale) + 8.0*term_scale);

    for (int k = 0; k < n; k++) {
        double nk2 = dvx[k]*dvx[k] + dvy[k]*dvy[k];
        if (!(nk2 > 0.0)) continue;
        double tol = vector_error * sqrt(nk2);
        for (int perpendicular = 0; perpendicular < 2; perpendicular++) {
            double nx = perpendicular ? -dvy[k] : dvx[k];
            double ny = perpendicular ?  dvx[k] : dvy[k];
            double lo = 0.0, hi = 0.0, trial_lo = 0.0, trial_hi = 0.0;
            for (int j = 0; j < n; j++) {
                double actual = nx*dvx[j] + ny*dvy[j];
                double px = 0.5 * (bxx*dx[j] + bxy*dy[j]);
                double py = 0.5 * (byx*dx[j] + byy*dy[j]);
                double trial = nx*px + ny*py;
                if (!isfinite(trial)) return -1.0;
                if (actual < lo) lo = actual;
                if (actual > hi) hi = actual;
                if (trial < trial_lo) trial_lo = trial;
                if (trial > trial_hi) trial_hi = trial;
            }
            if (trial_hi > hi + tol) {
                double a = (hi + tol) / trial_hi;
                if (a < alpha) alpha = a;
            }
            if (trial_lo < lo - tol) {
                double a = (lo - tol) / trial_lo;
                if (a < alpha) alpha = a;
            }
        }
    }
    return fmax(0.0, fmin(1.0, alpha));
}

#undef GFS_VEL_HD
#endif
