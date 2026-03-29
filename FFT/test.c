#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include<unistd.h>
#ifdef _OPENMP
#include<omp.h>
#endif
#include "eunha.h"
#include "fft.h"

#include "mpirms.h"


#define NXYZ 256



int threads_ok = 1;


SimParameters par,*simpar;
float ran2(long *);

long iseed = -92;

int threads_ok;

int main(int argc, char **argv){
	int myid,nid;
	ptrdiff_t tsize;

	/*
	mpi_fftw_initialize(argc, argv);
	*/

	MPI_Init(&argc, &argv);
	if(threads_ok) threads_ok = fftwf_init_threads();
	fftwf_mpi_init();
	fftwf_plan_with_nthreads(4);


	simpar = &par;

	FFTW_COMM(simpar) = MPI_COMM_WORLD;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	iseed = iseed*(myid+1);

	if(1){ 
#ifdef DEBUG
		char hostname[190]; 
		hostname[189] = '\0'; 
		gethostname(hostname,189); 
		printf("P%d host name is %s with pid= %d\n",myid,hostname, getpid());
		if(myid==-1){
			int kkk = 1;
			while(kkk) {
				kkk = 1;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	}


	MPI_COMM(simpar) = MPI_COMM_WORLD;
	/* note the integer type of nx/ny/nz */
	int nx,ny,nz;
	NX(simpar) = NY(simpar)= NZ(simpar) = 1024;
	int mx = 2*(NX(simpar)/2+1);

	/* note the integer type of local variables in fftw */
	ptrdiff_t local_grid_size,local_z_start, local_nz;
	ptrdiff_t local_grid_size_after_transpose,local_y_start_after_transpose, local_ny_after_transpose;

	LOCAL_GRID_SIZE(simpar) = fftwf_mpi_local_size_3d(NZ(simpar),NY(simpar),NX(simpar)/2+1, MPI_COMM_WORLD, &local_nz,
			&local_z_start);
	tsize = NX(simpar)*NY(simpar)*NZ(simpar);
	LOCAL_GRID_SIZE_AFTER_TRANSPOSE(simpar) = fftwf_mpi_local_size_3d_transposed(
			NZ(simpar),NY(simpar),NX(simpar)/2+1, FFTW_COMM(simpar),
			&LOCAL_NZ(simpar),
			&LOCAL_Z_START(simpar), 
			&LOCAL_NY_AFTER_TRANSPOSE(simpar), 
			&LOCAL_Y_START_AFTER_TRANSPOSE(simpar));

	/*
	float *rden = (float *)malloc(sizeof(float)*LOCAL_GRID_SIZE(simpar)*2);
	fftwf_complex *cden = (fftwf_complex*)rden;
	*/
	FFTW_RDEN(simpar) = (DenType*)malloc(sizeof(fftwf_complex)*LOCAL_GRID_SIZE(simpar));

	fftwf_plan p,ip;
	FFTW_NF_PLAN(simpar) = 
		fftwf_mpi_plan_dft_r2c_3d( NZ(simpar),NY(simpar),NX(simpar), FFTW_RDEN(simpar), FFTW_CDEN(simpar), MPI_COMM_WORLD,
                FFTW_ESTIMATE
				|
				FFTW_UNALIGNED
                );
	FFTW_NB_PLAN(simpar) = fftwf_mpi_plan_dft_c2r_3d( NZ(simpar),NY(simpar),NX(simpar), FFTW_CDEN(simpar), FFTW_RDEN(simpar), MPI_COMM_WORLD,
                FFTW_ESTIMATE
				|
				FFTW_UNALIGNED
                );
	/*
	FFTW_F_PLAN(simpar) = 
		fftwf_mpi_plan_dft_r2c_3d( NZ(simpar),NY(simpar),NX(simpar), rden, cden, MPI_COMM_WORLD,
                FFTW_MPI_TRANSPOSED_OUT
                |
                FFTW_ESTIMATE
				|
				FFTW_UNALIGNED
                );
	FFTW_B_PLAN(simpar) = fftwf_mpi_plan_dft_c2r_3d( NZ(simpar),NY(simpar),NX(simpar), cden, rden, MPI_COMM_WORLD,
                FFTW_MPI_TRANSPOSED_IN
                |
                FFTW_ESTIMATE
				|
				FFTW_UNALIGNED
                );
				*/
	float *rrden = (float *)malloc(sizeof(float)*LOCAL_GRID_SIZE(simpar)*2);
	fftwf_complex *ccden = (fftwf_complex*)rrden;

	ptrdiff_t i,j,k;
	for(i=0;i<LOCAL_GRID_SIZE(simpar)*2;i++){
		rrden[i] = ran2(&iseed);
	}
	printf("-P%d has %g %g %g %g\n",myid,rrden[0],rrden[100],rrden[42],rrden[34]);fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	fftwf_mpi_execute_dft_r2c(FFTW_NF_PLAN(simpar), rrden,ccden);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("0P%d has %g %g %g %g\n",myid,rrden[0],rrden[100],rrden[42],rrden[34]);fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),ccden, rrden);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("+P%d has %g %g %g %g\n",myid,rrden[0]/tsize,rrden[100]/tsize,rrden[42]/tsize,rrden[34]/tsize);fflush(stdout);


	MPI_Finalize();

}
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
