#ifndef EXAM_FULL_MOMENT_PROVIDER_MPI_H
#define EXAM_FULL_MOMENT_PROVIDER_MPI_H
#include <mpi.h>
#include "full_moment_provider.h"
/* All ranks ACK validated local release before the rank0 shared deletion. */
int full_moment_collective_release(MPI_Comm,FullMomentProvider*,FullMomentProviderRequest*,const char *geometry);
#endif
