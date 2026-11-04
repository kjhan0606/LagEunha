#ifndef GFS_KDK_LTS_H
#define GFS_KDK_LTS_H

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#define GFS_KDK_LTS_MAX_BIN 30

typedef struct {
    unsigned long long count[GFS_KDK_LTS_MAX_BIN + 1];
    unsigned long long ncell;
    double active_fraction;
    unsigned max_bin;
} GfsKdkLtsHistogram;

/* Largest power-of-two level whose step does not exceed candidate_dt. */
static inline int gfs_kdk_lts_bin_from_dt(double candidate_dt,
                                         double dt_base,
                                         unsigned max_bin)
{
    if (!isfinite(candidate_dt) || !isfinite(dt_base) ||
        candidate_dt <= 0.0 || dt_base <= 0.0 ||
        max_bin > GFS_KDK_LTS_MAX_BIN)
        return -1;

    double ratio = candidate_dt / dt_base;
    if (!isfinite(ratio) || ratio < 1.0)
        return 0;

    unsigned bin = 0;
    while (bin < max_bin && ratio >= 2.0) {
        ratio *= 0.5;
        ++bin;
    }
    return (int)bin;
}

/*
 * Limit neighboring bins so adjacent Voronoi cells differ by at most one.
 * edge_i/edge_j are local indices into bin[]. Cross-rank edges are added by
 * MPI halo exchange in the solver integration, not by this pure helper.
 */
static inline int gfs_kdk_lts_limit_neighbor_bins(unsigned char *bin,
                                                   size_t nedge,
                                                   const int *edge_i,
                                                   const int *edge_j,
                                                   unsigned max_bin)
{
    if (!bin || (nedge && (!edge_i || !edge_j)) ||
        max_bin > GFS_KDK_LTS_MAX_BIN)
        return 0;

    for (size_t e = 0; e < nedge; ++e)
        if (edge_i[e] < 0 || edge_j[e] < 0 ||
            bin[edge_i[e]] > max_bin || bin[edge_j[e]] > max_bin)
            return 0;

    for (unsigned pass = 0; pass <= max_bin; ++pass) {
        int changed = 0;
        for (size_t e = 0; e < nedge; ++e) {
            unsigned char *a = &bin[edge_i[e]];
            unsigned char *b = &bin[edge_j[e]];
            if ((unsigned)*a > (unsigned)*b + 1u) {
                *a = (unsigned char)(*b + 1u);
                changed = 1;
            }
            if ((unsigned)*b > (unsigned)*a + 1u) {
                *b = (unsigned char)(*a + 1u);
                changed = 1;
            }
        }
        if (!changed)
            return 1;
    }

    for (size_t e = 0; e < nedge; ++e)
        if ((unsigned)bin[edge_i[e]] + 1u < (unsigned)bin[edge_j[e]] ||
            (unsigned)bin[edge_j[e]] + 1u < (unsigned)bin[edge_i[e]])
            return 0;
    return 1;
}

/* Tick zero is a synchronization point, so every cell is active there. */
static inline int gfs_kdk_lts_active_at_tick(unsigned bin, uint64_t tick)
{
    if (bin > GFS_KDK_LTS_MAX_BIN)
        return 0;
    const uint64_t period = UINT64_C(1) << bin;
    return (tick % period) == 0;
}

static inline void gfs_kdk_lts_histogram_init(GfsKdkLtsHistogram *h)
{
    if (!h) return;
    for (unsigned i = 0; i <= GFS_KDK_LTS_MAX_BIN; ++i)
        h->count[i] = 0;
    h->ncell = 0;
    h->active_fraction = 0.0;
    h->max_bin = 0;
}

static inline int gfs_kdk_lts_histogram_add(GfsKdkLtsHistogram *h,
                                            unsigned bin)
{
    if (!h || bin > GFS_KDK_LTS_MAX_BIN)
        return 0;
    ++h->count[bin];
    ++h->ncell;
    if (bin > h->max_bin) h->max_bin = bin;
    h->active_fraction += ldexp(1.0, -(int)bin);
    return 1;
}

#endif
