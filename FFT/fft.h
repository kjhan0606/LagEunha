#ifndef FFTW3_H
#include "fftw3.h"
#include "fftw3-mpi.h"
#define FFTW3_H
#endif
void mpi_fftw_initialize(int , char **);
void set_mpi_fftw_3d_info( SimParameters *);
void set_mpi_fftw_3d_plans( SimParameters *);
void do_fftw_forward(FFTWGridInfo *, FFTW_INFO *);
void do_fftw_backward(FFTWGridInfo *, FFTW_INFO *);
void malloc_fftw_den(FFTWGridInfo *, FFTW_INFO *);
void free_fftw_den(FFTWGridInfo *, FFTW_INFO *);
/*
void mpi_fftw_finalize(FFTW_INFO *);
*/
void mpi_fftw_finalize(void);
void set_mpi_fftw_3d_all_info(FFTWGridInfo *, FFTW_INFO *, MPI_Comm );
