#ifndef GFS_SEED_BIN_CERTIFICATE_20261027_H
#define GFS_SEED_BIN_CERTIFICATE_20261027_H

#include <float.h>
#include <math.h>

#if defined(__CUDACC__)
#define GFS_CERT_HD __host__ __device__
#else
#define GFS_CERT_HD
#endif

/* True only when no point stored in the bin can pass the existing clipped
 * power-bisector/circumradius rejection. Boundary bins are treated as
 * unbounded along that coordinate, because linked-cell construction may
 * clamp points outside the nominal grid into edge bins. An uncertainty
 * margin makes floating-point ambiguity fall back to normal enumeration. */
GFS_CERT_HD static inline bool gfs_seed_bin_cannot_cut(
    double xi, double yi, int cx, int cy, int mx, int my,
    double xmin, double ymin, double cellsize, double owner_w2,
    double global_max_w2, double maxdist2)
{
    if (!(cellsize > 0.) || cx < 0 || cx >= mx || cy < 0 || cy >= my ||
        !isfinite(xi) || !isfinite(yi) || !isfinite(xmin) || !isfinite(ymin) ||
        !isfinite(cellsize) || !isfinite(owner_w2) || !isfinite(global_max_w2) ||
        !isfinite(maxdist2) || maxdist2 < 0.) return false;
    const double xlo = xmin + cx*cellsize, xhi = xlo + cellsize;
    const double ylo = ymin + cy*cellsize, yhi = ylo + cellsize;
    const double uncertainty = 64.*DBL_EPSILON*(1. + fabs(xi) + fabs(yi) +
                                fabs(xlo) + fabs(xhi) + fabs(ylo) + fabs(yhi));
    double dx = 0., dy = 0.;
    if (cx != 0 && cx != mx-1)
        dx = fmax(0., fmax(xlo-xi, xi-xhi)-uncertainty);
    if (cy != 0 && cy != my-1)
        dy = fmax(0., fmax(ylo-yi, yi-yhi)-uncertainty);
    const double d2 = dx*dx + dy*dy;
    if (!(d2 > 0.) || !isfinite(d2)) return false;
    /* For every candidate, wi-wj >= -delta. The clamped offset function
     * max(.05d,min(.95d,.5(d-delta/d))) is monotone increasing in d. */
    const double delta = fmax(0., global_max_w2-owner_w2);
    double fraction = .5-.5*delta/d2;
    fraction = fmax(.05, fmin(.95, fraction));
    const double lower2 = fraction*fraction*d2;
    if (!isfinite(lower2)) return false;
    return lower2 > maxdist2*(1. + 128.*DBL_EPSILON) + DBL_MIN;
}

#undef GFS_CERT_HD
#endif
