#ifndef EXAM_FULL_MOMENT_KEPLER_H
#define EXAM_FULL_MOMENT_KEPLER_H
#include <mpi.h>
int full_moment_kepler_is_config(const char *path);
int full_moment_kepler_run(const char *path,MPI_Comm communicator);
#endif
