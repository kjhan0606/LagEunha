#ifndef EXAM_FULL_MOMENT_COUPLED_CONTROL_H
#define EXAM_FULL_MOMENT_COUPLED_CONTROL_H
#include "full_moment_state_io.h"
#include "full_moment_coupled_state.h"
/* Callbacks execute collectively when used by MPI. evaluate must perform
 * exactly ONE fresh midpoint provider/map and return the actual physical map
 * residual, not an iterate-change surrogate. finalize performs agreement,
 * journal and admissibility checks with ZERO additional physical providers.
 * No callback may publish an accepted checkpoint. The caller commits only
 * after this controller returns success on every rank. */
typedef struct {
 void *context;
 int (*evaluate)(void*,const FullMomentState*,const FullMomentState*,int,
                 FullMomentState*,FullMomentCoupledResidual*,
                 long double*,size_t);
 int (*finalize)(void*,const FullMomentState*,const FullMomentState*,
                 const FullMomentState*,const FullMomentCoupledResidual*,int);
 int (*type2_fit)(void*,size_t,int,const long double*,
                  const long double *const[3],double[3],int*,int*,
                  long double*,long double*);
 int (*collective_status)(void*,int);
 int (*provider_count)(void*);
 double (*elapsed_seconds)(void*);
 int maximum_step_providers,maximum_total_providers;
 double wall_limit_seconds;
} FullMomentCoupledCallbacks;
typedef struct {
 int iterations,providers_before,providers_after,converged,final_checks_passed;
 int callback_status;
 FullMomentCoupledResidual last_residual;
 FullMomentCoupledResidual iteration_residual[10];
 int anderson_depth[10],anderson_qr_rank[10],anderson_fallback[10];
 int prediction_miss_streak[10];
 double anderson_gamma[10][3],anderson_lambda[10],anderson_condition[10];
 long double anderson_predicted_norm2[10];
} FullMomentCoupledControlReport;
enum {
 FM_COUPLED_OK=0,FM_COUPLED_INVALID=1,FM_COUPLED_EVALUATION_FAILED=2,
 FM_COUPLED_PROVIDER_BUDGET=3,FM_COUPLED_NOT_CONVERGED=4,
 FM_COUPLED_COUNTER_VIOLATION=5,FM_COUPLED_FINAL_CHECK_FAILED=6,
 FM_COUPLED_WALL_BUDGET=7
};
enum {
 FM_AA_FALLBACK_NONE=0,FM_AA_FALLBACK_NO_HISTORY=1,
 FM_AA_FALLBACK_QR_OR_COEFFICIENT=2,FM_AA_FALLBACK_POSITIVITY_BACKTRACK=3,
 FM_AA_FALLBACK_PURE_MAP=4,FM_AA_FALLBACK_METRIC_RESET=5,
 FM_AA_FALLBACK_PREDICTION_RESET=6
};
/* Initial endpoint guess is the unchanged base state. Each Picard iteration
 * evaluates M_h(base,F((base+trial)/2)). Success transfers the VERIFIED map
 * output, never the trial, into out. On failure out/base are unchanged. */
int full_moment_coupled_control(const FullMomentState *base,
 const FullMomentCoupledCallbacks *callbacks,FullMomentState *out,
 FullMomentCoupledControlReport *report);
/* Safeguarded Type-II Anderson least-squares fit. history[0] is the most
 * recent prior residual and up to three prior residuals may be supplied.
 * Rank loss/large coefficients reduce depth by dropping the oldest column.
 * Returns 1 for a usable fit, 0 for a safe pure-map fallback. */
int full_moment_coupled_type2_fit(size_t residual_count,int maximum_depth,
 const long double *current,const long double *const history[3],
 double gamma[3],int *used_depth,int *qr_rank,
 long double *predicted_norm2,long double *condition);
#endif
