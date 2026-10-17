#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
//#include "mpirms.h"
#include "fft.h"





void Zeldovich(SimParameters *simpar, DenType *den, DenType *fx, DenType *fy, DenType *fz){
	float fNL = FNL(simpar);
	float gNL = GNL(simpar);
	ptrdiff_t i,j,k;
	ptrdiff_t nx = NX(simpar);
	ptrdiff_t ny = NY(simpar);
	ptrdiff_t nz = NZ(simpar);
	ptrdiff_t local_nz = LOCAL_NZ(simpar);
	ptrdiff_t local_z_start = LOCAL_Z_START(simpar);

	ptrdiff_t mx = 2*(nx/2+1);
	float ky, xyy, kz, xyz;


	xyy = (float)nx/(float)ny;
	xyz = (float)nx/(float)nz;

	double rincube = 1./(double)nx/(double)ny/(double)nz; 



	fftwf_mpi_execute_dft_c2r(FFTWINFO(simpar).inp,(fftwf_complex*)den, den);


	for(i=0;i<mx*ny*local_nz;i++) den[i] = den[i]*rincube;


	double phistd, phimean, tphistd, tphimean;

	phistd = phimean = 0;

	for(k=0;k<local_nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++){
		phimean += den[i+mx*(j+ny*k)];
		phistd  += den[i+mx*(j+ny*k)] * den[i+mx*(j+ny*k)];
	}

	MPI_Comm comm = MPI_COMM(simpar);

	MPI_Reduce(&phistd, &tphistd, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	phistd = tphistd;
	MPI_Bcast(&phistd, 1, MPI_DOUBLE, 0, comm);

	MPI_Reduce(&phimean, &tphimean, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	phimean = tphimean;
	phimean = phimean*rincube;
	MPI_Bcast(&phimean, 1, MPI_DOUBLE, 0, comm);
	phistd = sqrt(phistd*rincube - phimean*phimean);
	for(k=0;k<local_nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++){
		ptrdiff_t ipos = i+mx*(j+ny*k);
		den[ipos] += fNL*(den[ipos]*den[ipos] -phistd*phistd)
			+ gNL*(den[ipos]*den[ipos]*den[ipos] - 3*phistd*phistd*den[ipos]);
	}

	fftwf_mpi_execute_dft_r2c(FFTWINFO(simpar).np, den, (fftwf_complex*)den);
	float twopi = 2*M_PI;

#ifdef _OPENMP
#pragma omp parallel for private(ky, kz)
#endif
	for(i=0;i<nx/2;i++){
		for(j=0;j<ny;j++){
			ky = xyy*j;
			if(ky > ny/2) ky = ky -ny;
			for(k=0;k<local_nz;k++){
				kz = xyz*( k + local_z_start);
				if(kz > nz/2) kz = kz - nz;
				ptrdiff_t ipos = 2*i + mx*(j+ny*k);
				fx[ipos]   = -twopi* i*den[ipos+1];
				fx[ipos+1] =  twopi* i*den[ipos];
				fy[ipos]   = -twopi*ky*den[ipos+1];
				fy[ipos+1] =  twopi*ky*den[ipos];
				fz[ipos]   = -twopi*kz*den[ipos+1];
				fz[ipos+1] =  twopi*kz*den[ipos];
			}
		}
	}
	fftwf_mpi_execute_dft_c2r(FFTWINFO(simpar).inp, (fftwf_complex*)fx, fx);
	fftwf_mpi_execute_dft_c2r(FFTWINFO(simpar).inp, (fftwf_complex*)fy, fy);
	fftwf_mpi_execute_dft_c2r(FFTWINFO(simpar).inp, (fftwf_complex*)fz, fz);


	for(i=0;i<mx*ny*local_nz;i++){
		fx[i] = fx[i]*rincube;
		fy[i] = fy[i]*rincube;
		fz[i] = fz[i]*rincube;
	}
}
