#ifndef GFS_PRESSURE_NEUTRAL_PROJECTION_20261027_H
#define GFS_PRESSURE_NEUTRAL_PROJECTION_20261027_H

#include <math.h>

/* Host projection onto h·rate=0. An exactly zero gradient imposes no
 * constraint and must leave every rate unchanged. For a nonzero gradient,
 * retain the original accumulation and update order. */
static inline int gfs_pressure_neutral_projection(
    double *rate, const double *gradient, const unsigned char *active,
    int n, double *alpha_out)
{
    if (!rate || !gradient || !active || !alpha_out || n <= 0) return -1;
    long double dot=0., norm=0.;
    int count=0;
    for (int i=0; i<n; ++i) {
        if (!active[i]) continue;
        ++count;
        if (!isfinite(gradient[i]) || !isfinite(rate[i])) return -2;
        dot+=(long double)gradient[i]*rate[i];
        norm+=(long double)gradient[i]*gradient[i];
    }
    if (!count) return -3;
    if (!(norm>0.)) {
        *alpha_out=0.;
        return 0;
    }
    double alpha=(double)(dot/norm);
    if (!isfinite(alpha)) return -4;
    for (int i=0; i<n; ++i) {
        if (active[i]) rate[i]-=alpha*gradient[i];
    }
    *alpha_out=alpha;
    return 0;
}

#endif
