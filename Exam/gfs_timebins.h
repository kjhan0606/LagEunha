#ifndef GFS_TIMEBINS_H
#define GFS_TIMEBINS_H

/*
 * Deterministic power-of-two local-time-bin planner.
 *
 * This header only assigns compatible time bins. It does not advance fluid
 * state or synchronize finite-volume face fluxes; callers must do both before
 * using the bins in an integrator. The CSR graph must contain both directions
 * of every cell-neighbor relation, including partition/halo neighbors.
 */
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#define GFS_TIMEBINS_MAX_LEVEL 30u

enum {
    GFS_TIMEBINS_OK = 0,
    GFS_TIMEBINS_INVALID = -1,
    GFS_TIMEBINS_NOMEM = -2,
    GFS_TIMEBINS_TOO_FINE = -3
};

/* Return the coarsest allowed level whose interval is <= dt_limit. */
static inline int gfs_timebins_level_for_limit(double dt_limit,
                                               double macro_dt,
                                               unsigned max_level,
                                               unsigned *level_out)
{
    if (!level_out || !isfinite(dt_limit) || !(dt_limit > 0.0) ||
        !isfinite(macro_dt) || !(macro_dt > 0.0) ||
        max_level > GFS_TIMEBINS_MAX_LEVEL) return GFS_TIMEBINS_INVALID;
    double interval = macro_dt;
    unsigned level = 0;
    while (interval > dt_limit && level < max_level) {
        interval *= 0.5;
        ++level;
    }
    if (interval > dt_limit) return GFS_TIMEBINS_TOO_FINE;
    *level_out = level;
    return GFS_TIMEBINS_OK;
}

/*
 * Assign the coarsest power-of-two interval not exceeding each local limit,
 * then refine bins until every adjacent pair differs by at most one level.
 *
 * level 0 has interval macro_dt; level l has interval macro_dt / 2^l.
 * For MPI closure, n/CSR must include halo vertices and dt limits exchanged
 * from their owners; only the levels for owned vertices are authoritative.
 * Output is transactional: level_out is unchanged on any error.
 */
static inline int gfs_timebins_assign(
    size_t n,
    const double *dt_local,
    double macro_dt,
    unsigned max_level,
    const size_t *neighbor_offset,
    const uint32_t *neighbor_index,
    size_t neighbor_count,
    unsigned *level_out)
{
    if (!isfinite(macro_dt) || !(macro_dt > 0.0) ||
        max_level > GFS_TIMEBINS_MAX_LEVEL) return GFS_TIMEBINS_INVALID;
    if (n == 0) return (neighbor_count == 0) ? GFS_TIMEBINS_OK : GFS_TIMEBINS_INVALID;
    if (!dt_local || !level_out || !neighbor_offset ||
        (neighbor_count && !neighbor_index)) return GFS_TIMEBINS_INVALID;
    if (n > SIZE_MAX / sizeof(unsigned) || n > SIZE_MAX / sizeof(size_t))
        return GFS_TIMEBINS_NOMEM;
    if (neighbor_offset[0] != 0 || neighbor_offset[n] != neighbor_count)
        return GFS_TIMEBINS_INVALID;
    for (size_t i = 0; i < n; ++i) {
        if (neighbor_offset[i] > neighbor_offset[i + 1] ||
            neighbor_offset[i + 1] > neighbor_count) return GFS_TIMEBINS_INVALID;
        for (size_t q = neighbor_offset[i]; q < neighbor_offset[i + 1]; ++q)
            if ((size_t)neighbor_index[q] >= n) return GFS_TIMEBINS_INVALID;
    }

    unsigned *work = (unsigned *)malloc(n * sizeof(*work));
    size_t *queue = (size_t *)malloc(n * sizeof(*queue));
    unsigned char *queued = (unsigned char *)calloc(n, sizeof(*queued));
    if (!work || !queue || !queued) {
        free(work); free(queue); free(queued);
        return GFS_TIMEBINS_NOMEM;
    }

    for (size_t i = 0; i < n; ++i) {
        if (!isfinite(dt_local[i]) || !(dt_local[i] > 0.0)) {
            free(work); free(queue); free(queued);
            return GFS_TIMEBINS_INVALID;
        }
        unsigned level = 0;
        int rc = gfs_timebins_level_for_limit(dt_local[i], macro_dt,
                                               max_level, &level);
        if (rc != GFS_TIMEBINS_OK) {
            free(work); free(queue); free(queued);
            return rc;
        }
        work[i] = level;
    }

    /* Propagate fine-step requirements outward over the symmetric graph. */
    size_t head = 0, tail = 0, pending = 0;
    for (size_t i = 0; i < n; ++i) {
        queue[tail] = i;
        tail = (tail + 1) % n;
        queued[i] = 1;
        ++pending;
    }
    while (pending) {
        size_t i = queue[head];
        head = (head + 1) % n;
        --pending;
        queued[i] = 0;
        for (size_t q = neighbor_offset[i]; q < neighbor_offset[i + 1]; ++q) {
            size_t j = (size_t)neighbor_index[q];
            if (work[j] + 1u < work[i]) {
                work[j] = work[i] - 1u;
                if (!queued[j]) {
                    queue[tail] = j;
                    tail = (tail + 1) % n;
                    queued[j] = 1;
                    ++pending;
                }
            }
        }
    }

    for (size_t i = 0; i < n; ++i) level_out[i] = work[i];
    free(work); free(queue); free(queued);
    return GFS_TIMEBINS_OK;
}

static inline double gfs_timebins_interval(double macro_dt, unsigned level)
{
    return ldexp(macro_dt, -(int)level);
}

/* tick ranges from 0 through 2^max_level, inclusive of the macro endpoint. */
static inline int gfs_timebins_active(unsigned level, unsigned max_level,
                                     uint64_t tick)
{
    if (max_level > GFS_TIMEBINS_MAX_LEVEL || level > max_level) return 0;
    uint64_t macro_ticks = UINT64_C(1) << max_level;
    if (tick > macro_ticks) return 0;
    uint64_t interval_ticks = UINT64_C(1) << (max_level - level);
    return tick % interval_ticks == 0;
}

#endif
