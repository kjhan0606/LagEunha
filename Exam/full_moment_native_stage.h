#ifndef EXAM_FULL_MOMENT_NATIVE_STAGE_H
#define EXAM_FULL_MOMENT_NATIVE_STAGE_H
#include <mpi.h>
#include <stdint.h>
int full_moment_native_select_device(MPI_Comm);
/* Collective native GPU stage on one host; root global face data, rank-owned
 * real gravity buffers, canonical order and physical boundary roles intact. */
int full_moment_native_stage(MPI_Comm,const char *global_root,const char *owned_root,
 uint64_t epoch,uint64_t run_token,int threads);
#endif
