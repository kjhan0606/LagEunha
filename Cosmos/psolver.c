#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
//#include "mpirms.h"
#include "pmmain.h"
#include "timerutil.h"


void psolver(SimParameters *simpar, GridInfo *fftwgrid){
	ptrdiff_t i,j,k;
	ptrdiff_t local_ny_after_transpose = LOCAL_NY_AFTER_TRANSPOSE(simpar);
	ptrdiff_t local_y_start_after_transpose = LOCAL_Y_START_AFTER_TRANSPOSE(simpar);
	ptrdiff_t nx,ny,nz,mx;
	double rnx,rny,rnz;
	float cputime0[2];
	float cputime1[2];

	rnx = nx = NX(simpar);
	rny = ny = NY(simpar);
	rnz = nz = NZ(simpar);
	mx = 2*(nx/2+1);
	float xyy = rnx/rny;
	float xyz = rnx/rnz;
	double fact = -1.L/M_PI;
	double sfact = -pow(M_PI*GRV_EFOLD(simpar)/rnx, 2);
	double rngc = rnx*rny*rnz;
	double irngc = 1.L/rngc;
	double rngm6 = 1.L/rngc/rngc;
	double poten,poten1, poten2;
	poten1 = poten2 = 0;
#ifdef DEBUG
	DEBUGPRINT("P%d has %g %g\n",MYID(simpar), GRV_EPSILON(simpar), GRV_EFOLD(simpar));
#endif


	/*
	TIMER_START(1);
	*/


	DenType *den = (DenType*)(fftwgrid+1);
		


	for(j=0;j<local_ny_after_transpose;j++){
		float wy = xyy*(j+local_y_start_after_transpose);
		if(wy >rny/2) wy = wy - rny;
#ifdef _OPENMP
#pragma omp parallel for private(i,k) reduction (+:poten1, poten2)
#endif
		for(k=0;k<nz;k++){
			float wz = xyz*k;
			if(wz >rnz/2) wz = wz-rnz;
			for(i=0;i<=nx/2;i++) {
				double wave = i*i + wy*wy + wz*wz;
				if(wave == 0.L) {
					den[2*i   + mx*(k+nz*j)] = 0;
					den[2*i+1 + mx*(k+nz*j)] = 0;
				}
				else {
					double factr = fact/wave;
					double factor = factr*exp(wave*sfact);
					double delksq = (den[2*i + mx*(k+nz*j)]* den[2*i + mx*(k+nz*j)] + 
						den[2*i+1 + mx*(k+nz*j)]* den[2*i+1 + mx*(k+nz*j)])*rngm6;
					if(i==0) poten1 += factor*delksq;
					else poten2 += factor*delksq;
					den[2*i   + mx*(k+nz*j)] = factor* den[2*i   + mx*(k+nz*j)]*irngc;
					den[2*i+1 + mx*(k+nz*j)] = factor* den[2*i+1 + mx*(k+nz*j)]*irngc;
				}
			}
		}
	}
	poten = 0.5*(poten1+ 2*poten2);
	double tpoten;
	MPI_Reduce(&poten, &tpoten, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM(simpar));
	poten = tpoten;
	MPI_Bcast(&poten, 1, MPI_DOUBLE, 0, MPI_COMM(simpar));
	/*
	TIMER_STOP(1);
	if(MYID(simpar)==0) printf("%d CPU(Poisson) = %g\n",STEPCOUNT(simpar), ELAPSED_TIME(1));
	*/
}
