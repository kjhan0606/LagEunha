#ifndef EXAM_FULL_MOMENT_NATIVE_FACTOR_H
#define EXAM_FULL_MOMENT_NATIVE_FACTOR_H
#include <mpi.h>
#include "full_moment_provider.h"
#include "full_moment_state_io.h"
/* Tentative owned midpoint factor/state update only. The controller must run
 * geometry, constrained-density and conservation gates before acceptance. */
int full_moment_native_factor(MPI_Comm,FullMomentProvider*,FullMomentProviderRequest*,
 const FullMomentState *base,const FullMomentState *linear,const char *stage_root,
 const char *owned_root,const char *rates_root,const char *candidate_root,
 uint64_t candidate_epoch,double h,double candidate_time,int threads,int shadow_audit);
#endif
