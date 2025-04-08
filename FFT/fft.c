/* For historical reason in developing the GOTPM and Eunha codes that are based on the Fortran language,
 * the array for the fftw is in column-major order. 
 * Therefore A(x,y,z) = A(x+mx*(y+ny*z)) when expressed in the one-dimensional array.
 * To sum up, the order of dimensions in the code is reversed actually compared to the C convention,
 * while we decided to use the Fortran convention.*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include<complex.h>
#include<omp.h>
#include "eunha.h"
/*
#include "fftw3.h"
#include "fftw3-mpi.h"
*/





static int fftwinitflag = 1;
int threads_ok;

void mpi_fftw_initialize(int argc, char **argv){
	if(fftwinitflag){
#ifdef _OPENMP
		int provided;
		if(!MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,&provided) ){
			DEBUGPRINT0("Something wrong in the MPI_Init_thread\n");
		}
		DEBUGPRINT("Provided thread num %d & set nthreads= %d\n",provided, omp_get_num_threads());
		threads_ok = provided >= MPI_THREAD_FUNNELED;
		if(threads_ok) threads_ok = fftwf_init_threads();
		fftwf_mpi_init();
		if(threads_ok) fftwf_plan_with_nthreads(omp_get_num_threads());
#else
		MPI_Init(&argc, &argv);
		fftwf_mpi_init();
#endif
		fftwinitflag = 0;
	}
}
#ifdef OBSOLTE_NOT_USE_THIS_ANYMORE
void set_mpi_fftw_3d_info( FFTWGridInfo *fftwgrid,FFTW_INFO *fftwinfo, MPI_Comm com){
	fftwinfo->com = com;
	fftwinfo->local_grid_size = 
		fftwf_mpi_local_size_3d( fftwgrid->nz,fftwgrid->ny,fftwgrid->nx/2+1,
			fftwinfo->com,
			&(fftwinfo->local_nz),
			&(fftwinfo->local_z_start) 
			);
	fftwinfo->local_grid_size_after_transpose = 
		fftwf_mpi_local_size_3d_transposed(
			fftwgrid->nz,fftwgrid->ny,fftwgrid->nx/2+1, 
			fftwinfo->com,
			&(fftwinfo->local_nz), &(fftwinfo->local_z_start),
			&(fftwinfo->local_ny_after_transpose), 
			&(fftwinfo->local_y_start_after_transpose) 
			);
}
void set_mpi_fftw_3d_plans(FFTWGridInfo *fftwgrid,FFTW_INFO *fftwinfo, MPI_Comm com){
	fftwinfo->p = fftwf_mpi_plan_dft_r2c_3d(
				fftwgrid->nz, fftwgrid->ny, fftwgrid->nx, 
				fftwgrid->den.rden, 
				fftwgrid->den.cden, 
				fftwinfo->com,  
				FFTW_MPI_TRANSPOSED_OUT
				|
				FFTW_ESTIMATE
				);
	fftwinfo->ip = fftwf_mpi_plan_dft_c2r_3d(
			fftwgrid->nz, fftwgrid->ny, fftwgrid->nx, 
			fftwgrid->den.cden,
			fftwgrid->den.rden, 
			fftwinfo->com,  
			FFTW_MPI_TRANSPOSED_IN
			|
			FFTW_ESTIMATE
			);
	fftwinfo->np = fftwf_mpi_plan_dft_r2c_3d(
				fftwgrid->nz, fftwgrid->ny, fftwgrid->nx, 
				fftwgrid->den.rden, 
				fftwgrid->den.cden, 
				fftwinfo->com,  
				FFTW_ESTIMATE
				);
	fftwinfo->inp = fftwf_mpi_plan_dft_c2r_3d(
			fftwgrid->nz, fftwgrid->ny, fftwgrid->nx, 
			fftwgrid->den.cden,
			fftwgrid->den.rden, 
			fftwinfo->com,  
			FFTW_ESTIMATE
			);
}
#else
void set_mpi_fftw_3d_info( SimParameters *simpar){
	FFTWGridInfo *fftwgrid;
	FFTW_INFO *fftwinfo;
	FFTW_COMM(simpar) = MPI_COMM(simpar);
	LOCAL_GRID_SIZE(simpar) = 
		fftwf_mpi_local_size_3d( NZ(simpar),NY(simpar),NX(simpar)/2+1, 
				FFTW_COMM(simpar),
				&LOCAL_NZ(simpar),
				&LOCAL_Z_START(simpar)
				);
#ifdef DEBUG
	DEBUGPRINT("-P%d has fftw domain nx/ny/nz= %ld %ld %ld with local_z %ld %ld\n",MYID(simpar),NX(simpar),NY(simpar),NZ(simpar),LOCAL_NZ(simpar), LOCAL_Z_START(simpar));
#endif
	LOCAL_GRID_SIZE_AFTER_TRANSPOSE(simpar) = 
		fftwf_mpi_local_size_3d_transposed(
			NZ(simpar),NY(simpar),NX(simpar)/2+1, 
			FFTW_COMM(simpar),
			&LOCAL_NZ(simpar),
			&LOCAL_Z_START(simpar), 
			&LOCAL_NY_AFTER_TRANSPOSE(simpar),
			&LOCAL_Y_START_AFTER_TRANSPOSE(simpar)
			);
#ifdef DEBUG
	DEBUGPRINT("+P%d has fftw domain nx/ny/nz= %ld %ld %ld with local_z %ld %ld\n",MYID(simpar),NX(simpar),NY(simpar),NZ(simpar),LOCAL_NY_AFTER_TRANSPOSE(simpar), LOCAL_Y_START_AFTER_TRANSPOSE(simpar));
#endif
}
void set_mpi_fftw_3d_plans(SimParameters *simpar){
	FFTW_COMM(simpar) = MPI_COMM(simpar);
	FFTW_RDEN(simpar)  = (DenType*)malloc(sizeof(fftwf_complex)*LOCAL_GRID_SIZE(simpar));

#ifdef DEBUG
	if(MYID(simpar)==0) DEBUGPRINT("P%d initilize the memory for fftw %ld\n",MYID(simpar), LOCAL_GRID_SIZE(simpar));
#endif
	FFTW_F_PLAN(simpar) = fftwf_mpi_plan_dft_r2c_3d(
			NZ(simpar), NY(simpar), NX(simpar),
			FFTW_RDEN(simpar),
			FFTW_CDEN(simpar),
			FFTW_COMM(simpar), 
			FFTW_MPI_TRANSPOSED_OUT | FFTW_ESTIMATE | FFTW_UNALIGNED 
			);
#ifdef DEBUG
	if(MYID(simpar)==0) DEBUGPRINT("P%d passed forward plan  %ld\n",MYID(simpar), LOCAL_GRID_SIZE(simpar));
#endif
	FFTW_B_PLAN(simpar) = fftwf_mpi_plan_dft_c2r_3d(
			NZ(simpar), NY(simpar), NX(simpar),
			FFTW_CDEN(simpar),
			FFTW_RDEN(simpar),
			FFTW_COMM(simpar), 
			FFTW_MPI_TRANSPOSED_IN | FFTW_ESTIMATE | FFTW_UNALIGNED 
			);
#ifdef DEBUG
	if(MYID(simpar)==0) DEBUGPRINT("P%d passed backward plan  %ld\n",MYID(simpar), LOCAL_GRID_SIZE(simpar));
#endif
	FFTW_NF_PLAN(simpar) = fftwf_mpi_plan_dft_r2c_3d(
			NZ(simpar), NY(simpar), NX(simpar),
			FFTW_RDEN(simpar),
			FFTW_CDEN(simpar),
			FFTW_COMM(simpar), 
			FFTW_ESTIMATE | FFTW_UNALIGNED
			);
#ifdef DEBUG
	if(MYID(simpar)==0) DEBUGPRINT("P%d passed forward-out-of-place plan  %ld\n",MYID(simpar), LOCAL_GRID_SIZE(simpar));
#endif
	FFTW_NB_PLAN(simpar) = fftwf_mpi_plan_dft_c2r_3d(
			NZ(simpar), NY(simpar), NX(simpar),
			FFTW_CDEN(simpar),
			FFTW_RDEN(simpar),
			FFTW_COMM(simpar), 
			FFTW_ESTIMATE | FFTW_UNALIGNED
			);
	free(FFTW_RDEN(simpar));
#ifdef DEBUG
	if(MYID(simpar)==0) DEBUGPRINT("P%d passed backward-out-of-place plan  %ld\n",MYID(simpar), LOCAL_GRID_SIZE(simpar));
#endif
}
#endif

void do_fftw_forward(FFTWGridInfo *fftwgrid, FFTW_INFO *fftwinfo){
	fftwf_mpi_execute_dft_r2c(
			fftwinfo->np, 
			(fftwgrid->den).rden, 
			(fftwgrid->den).cden
			);
}

void do_fftw_backward(FFTWGridInfo *fftwgrid, FFTW_INFO *fftwinfo){
	fftwf_mpi_execute_dft_c2r(
			fftwinfo->inp, 
			(fftwgrid->den).cden, 
			(fftwgrid->den).rden
			);
}
void malloc_fftw_den(FFTWGridInfo *fftwgrid, FFTW_INFO *fftwinfo){
	fftwgrid->den.rden = (DenType*) fftwf_malloc(sizeof(DenType)*fftwinfo->local_grid_size*2);
}
void free_fftw_den(FFTWGridInfo *fftwgrid, FFTW_INFO *fftwinfo){
	fftwf_free(fftwgrid->den.rden);
}

void mpi_fftw_finalize(void){
	if(fftwinitflag == 0){
		/*
		fftwf_mpi_cleanup();
		fftwf_destroy_plan(fftwinfo->p);
		fftwf_destroy_plan(fftwinfo->ip);
		*/
		fftwinitflag = 1;
		MPI_Finalize();
	}
}

