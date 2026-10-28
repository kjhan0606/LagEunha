#ifndef EXAM_FULL_MOMENT_UPDATE_H
#define EXAM_FULL_MOMENT_UPDATE_H
#include <stdint.h>
/* Construct tentative owned candidates. MPI collective acceptance and full QR
 * conditioning audits are the caller's responsibility. Psi stays long double.
 * All outputs remain unchanged if any input/cell check fails.
 */
int full_moment_build_candidate(int n,double h,const double *base_state_20,
 const double *fixed_mass,const double *stage_94,const double *generator_2,
 const double *weight_rate,const long double *base_Psi,const double *Psi_rate,
 const double *factor_result_63,double *candidate_state_20,
 long double *candidate_Psi,double *candidate_factors_4);
#endif
