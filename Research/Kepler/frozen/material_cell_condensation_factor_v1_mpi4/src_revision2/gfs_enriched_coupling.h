#ifndef GFS_ENRICHED_COUPLING_H
#define GFS_ENRICHED_COUPLING_H
#include "gfs_enriched_euler.h"
typedef void (*GfsEnrichSum)(void *user,double *values,int count);
enum {
    EC_POINTS,EC_PAIRS,EC_NEWTON,EC_CG,EC_BACKTRACKS,EC_RESIDUAL,
    EC_MASS_SUM,EC_GAUGE,EC_SYMMETRY,EC_RAYLEIGH,EC_SECANT,
    EC_HLL_POINTS,EC_MAX_RHS_MINUS_GEOMETRY,EC_NREPORT
};
/* Full initial moving-cell Euler RHS and screened weight-rate solve.
 * Allreduce is supplied by the native MPI caller. This is an RHS audit,
 * not a dynamic-weight geometry or time-stepping implementation. */
int gpu_enriched_coupled_rhs(GPUContext *ctx,int np,int n,int rank,
    const unsigned long long *ids,const GfsEnrichAudit *global,
    const GfsEnrichConfig *cfg,const double perturb[2],
    GfsEnrichSum sum,void *user,double *wdot,double *rhs,double *geometry,
    double report[EC_NREPORT]);
#endif
