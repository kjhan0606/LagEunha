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

/*
float cputime0[1];
float cputime1[1];
*/

#define MAX(a,b) ( (a)>(b)? (a):(b) )

#define ROneStepEvolv(simpar, type, pfact) for(i=0;i<type##_NP(simpar);i++){\
		PosType xd,yd,zd;\
		xd = (type##_BP(simpar)+i)->vx * pfact;\
		yd = (type##_BP(simpar)+i)->vy * pfact;\
		zd = (type##_BP(simpar)+i)->vz * pfact;\
		PosType rd = xd*xd + yd*yd + zd*zd;\
		ddmax = MAX(ddmax,rd);\
		(type##_BP(simpar)+i)->x += xd;\
		(type##_BP(simpar)+i)->y += yd;\
		(type##_BP(simpar)+i)->z += zd;\
	}

#define Periodicity(simpar, type)\
	for(i=0;i<type##_NP(simpar);i++){\
		(type##_BP(simpar)+i)->x = fmod((type##_BP(simpar)+i)->x+rnx, rnx);\
		(type##_BP(simpar)+i)->y = fmod((type##_BP(simpar)+i)->y+rny, rny);\
		(type##_BP(simpar)+i)->z = fmod((type##_BP(simpar)+i)->z+rnz, rnz);\
	}

void onestepforwardposition(SimParameters *simpar){
	ptrdiff_t i,j,k;
	/*
	ptrdiff_t local_ny_after_transpose = LOCAL_NY_AFTER_TRANSPOSE(simpar);
	ptrdiff_t local_y_start_after_transpose = LOCAL_Y_START_AFTER_TRANSPOSE(simpar);
	*/
	ptrdiff_t nx,ny,nz,mx;
	PosType rnx,rny,rnz;
	rnx = nx = NX(simpar);
	rny = ny = NY(simpar);
	rnz = nz = NZ(simpar);
	float pfact = (Evol_PFACT(simpar) = ASTEP(simpar)*NX(simpar));
	float tdmax = 0;
#ifdef DEBUG
	DEBUGPRINT("P%d has pfact= %g\n",MYID(simpar), pfact);
#endif
	/*
#ifdef _OPENMP
#pragma omp parallel 
#endif
*/
	{
		float ddmax = -9.e20;
#ifdef _OPENMP
#pragma omp parallel for private(i) reduction(max:ddmax)
#endif
		ROneStepEvolv(simpar, DM, pfact);
#ifdef _OPENMP
#pragma omp parallel for private(i)reduction(max:ddmax)
#endif
		ROneStepEvolv(simpar, SPH, pfact);
#ifdef _OPENMP
#pragma omp parallel for private(i) reduction(max:ddmax)
#endif
		ROneStepEvolv(simpar, STAR, pfact);
#ifdef _OPENMP
#pragma omp parallel for private(i)reduction(max:ddmax)
#endif
		ROneStepEvolv(simpar, AGN, pfact);
		{
			/*
#ifdef _OPENMP
#pragma omp critical
#endif
*/
			if(ddmax > tdmax) tdmax = ddmax;
		}
	}

	float dmax = sqrt(tdmax);
	MPI_Reduce(&dmax, &tdmax, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM(simpar));
	if(IndT_FLAG(simpar)=='Y'){
		if(MYID(simpar)==0) printf("%d.%d Maximum Displacement = %g\n",STEPCOUNT(simpar),IndT_NOWTSUBDIV(simpar), tdmax);
	}
	else {
		if(MYID(simpar)==0) printf("Step %d Maximum Displacement = %g\n",STEPCOUNT(simpar), tdmax);
	}
#ifndef XYZDBL
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	Periodicity(simpar, DM);
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	Periodicity(simpar, SPH);
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	Periodicity(simpar, STAR);
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	Periodicity(simpar, AGN);
#endif
}

#define TreeROneStepEvolv(simpar, type, pfact) for(i=0;i<type##_NP(simpar);i++){\
	PosType xd,yd,zd;\
	xd = (type##_TBP(simpar)+i)->vx * pfact;\
	yd = (type##_TBP(simpar)+i)->vy * pfact;\
	zd = (type##_TBP(simpar)+i)->vz * pfact;\
	PosType rd = xd*xd + yd*yd + zd*zd;\
	ddmax = MAX(ddmax,rd);\
	(type##_TBP(simpar)+i)->x += xd;\
	(type##_TBP(simpar)+i)->y += yd;\
	(type##_TBP(simpar)+i)->z += zd;\
}
#define TreePeriodicity(simpar, type) for(i=0;i<type##_NP(simpar);i++){\
	(type##_TBP(simpar)+i)->x = fmod((type##_TBP(simpar)+i)->x+rnx, rnx);\
	(type##_TBP(simpar)+i)->y = fmod((type##_TBP(simpar)+i)->y+rny, rny);\
	(type##_TBP(simpar)+i)->z = fmod((type##_TBP(simpar)+i)->z+rnz, rnz);\
}

void treeonestepforwardposition(SimParameters *simpar){
	ptrdiff_t i,j,k;
	/*
	ptrdiff_t local_ny_after_transpose = LOCAL_NY_AFTER_TRANSPOSE(simpar);
	ptrdiff_t local_y_start_after_transpose = LOCAL_Y_START_AFTER_TRANSPOSE(simpar);
	*/
	ptrdiff_t nx,ny,nz,mx;
	PosType rnx,rny,rnz;
	rnx = nx = NX(simpar);
	rny = ny = NY(simpar);
	rnz = nz = NZ(simpar);
	float pfact = Evol_PFACT(simpar);
	float tdmax = 0;
	/*
#ifdef _OPENMP
#pragma omp parallel 
#endif
*/
	{
		float ddmax = -9.e20;
#ifdef _OPENMP
#pragma omp parallel for private(i) reduction(max:ddmax)
#endif
		TreeROneStepEvolv(simpar, DM, pfact);
#ifdef _OPENMP
#pragma omp parallel for private(i) reduction(max:ddmax)
#endif
		TreeROneStepEvolv(simpar, SPH, pfact);
#ifdef _OPENMP
#pragma omp parallel for private(i) reduction(max:ddmax)
#endif
		TreeROneStepEvolv(simpar, STAR, pfact);
#ifdef _OPENMP
#pragma omp parallel for private(i) reduction(max:ddmax)
#endif
		TreeROneStepEvolv(simpar, AGN, pfact);
		{
			/*
#ifdef _OPENMP
#pragma omp critical
#endif
*/
			if(ddmax > tdmax) tdmax = ddmax;
		}
	}

	float dmax = sqrt(tdmax);
	MPI_Reduce(&dmax, &tdmax, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM(simpar));
	if(IndT_FLAG(simpar)=='Y'){
		if(MYID(simpar)==0) printf("%d.%d Maximum Displacement = %g with pfact= %g\n",STEPCOUNT(simpar),
				IndT_NOWTSUBDIV(simpar), tdmax, pfact);
	}
	else {
		if(MYID(simpar)==0) printf("Step %d Maximum Displacement = %g with pfact= %g\n",STEPCOUNT(simpar), tdmax, pfact);
	}
#ifndef XYZDBL
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	TreePeriodicity(simpar, DM);
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	TreePeriodicity(simpar, SPH);
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	TreePeriodicity(simpar, STAR);
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	TreePeriodicity(simpar, AGN);
#endif
}
