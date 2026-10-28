#ifndef EXAM_FULL_MOMENT_COUPLED_STATE_H
#define EXAM_FULL_MOMENT_COUPLED_STATE_H
#include <stddef.h>
#include <stdint.h>
/* Physical centered moments of one fixed-mass cell, decoded from raw state.
 * No use of stored L, no positivity floor, no projection to a reference. */
typedef struct {
 long double center[2],velocity[2],S[4],C[4],R[4],thermal;
} FullMomentCenteredState;
int full_moment_centered_decode(const double state[20],double mass,
 FullMomentCenteredState *out);
/* Trial arithmetic midpoint of x/W/raw/ledger/Psi. Cholesky is reconstructed
 * from the midpoint raw moments, not averaged from the endpoint factors.
 * Arrays have 20*n+4, n, 4*n elements. Outputs remain unchanged on failure.
 * This function performs no provider call and accepts no time step. */
int full_moment_coupled_midpoint(int n,const double *base,const double *trial,
 const double *mass,const long double *base_Psi,const long double *trial_Psi,
 double *midpoint,long double *midpoint_Psi,double *midpoint_L,int *failed_cell);
typedef struct {
 long double origin[2],scale,reference[6],projection[6],sigma,rank_error;
 int dependent;
} FullMomentPotentialNorm;
enum {
 FM_RES_X,FM_RES_W,FM_RES_CENTER,FM_RES_VELOCITY,FM_RES_S,FM_RES_C,FM_RES_R,
 FM_RES_THERMAL,FM_RES_PSI_REMAINDER,FM_RES_DEPENDENT_PSI,FM_RES_COUNT
};
typedef struct {
 long double maximum[FM_RES_COUNT];
 long long cell[FM_RES_COUNT];
 int dependent_cells;
 /* Solver coordinates use the timestep-fixed provider metric. Acceptance
  * maxima above always use the actual provider fit from this iteration. */
 long double vector_norm2,psi_remainder_numerator_max;
 long double actual_sigma_min,actual_sigma_max,actual_rank_margin_min;
 uint64_t actual_metric_fingerprint,anderson_metric_fingerprint;
 int metric_epoch,metric_reset_cells;
} FullMomentCoupledResidual;
/* Fixed-mass map residual, not an iterate-difference proxy. Areas and velocity/
 * thermal scales come from the base. Potential norm parameters must come from
 * the actual provider fit and retain its declared polynomial/rank contract.
 * Dependent Psi constraints use the existing binary80 evaluation bound.
 * No scale is increased using a numerical floor. */
int full_moment_coupled_residual(int n,double gamma,const double *base,
 const double *mass,const double *base_area,const double *trial,const double *mapped,
 const long double *trial_Psi,const long double *mapped_Psi,
 const FullMomentPotentialNorm *potential,FullMomentCoupledResidual *out,int *failed_cell);
/* Signed 21-component-per-cell residual vector. Physical scales and the
 * provider potential fit are supplied by the caller. The controller uses a
 * fit frozen at the accepted timestep's iteration-0 provider for all
 * Anderson coordinates; actual acceptance remains in the function above. */
int full_moment_coupled_residual_vector(int n,double gamma,const double *base,
 const double *mass,const double *base_area,const double *trial,const double *mapped,
 const long double *trial_Psi,const long double *mapped_Psi,
 const FullMomentPotentialNorm *potential,long double *out,size_t out_count,
 int *failed_cell);
int full_moment_coupled_converged(const FullMomentCoupledResidual *residual);
#endif
