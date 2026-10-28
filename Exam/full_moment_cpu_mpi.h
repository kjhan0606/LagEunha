#ifndef EXAM_FULL_MOMENT_CPU_MPI_H
#define EXAM_FULL_MOMENT_CPU_MPI_H
#include <mpi.h>
int full_moment_cpu_mpi_validation(const char*,MPI_Comm);
int full_moment_join_rows(MPI_Comm,int rank,int ranks,const char *local,const char *global);
#endif
