#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<complex.h>
#include<mpi.h>
#ifdef _OPENMP
#include<omp.h>
#endif
#include "eunha.h"
#include "fft.h"
//#include "mpirms.h"


/* This function shrink data from phi(mx,ny,local_nz) to phi(nx,ny,local_nz) */
void selfdiffphi(SimParameters *simpar, DenType *phi, enum dimension ixyz){
	float twopi = 2*M_PI;
	float fourpi2 = twopi*twopi;

	ptrdiff_t i,j,k;
	ptrdiff_t nx,ny,nz;
	ptrdiff_t local_nz = LOCAL_NZ(simpar);
	ptrdiff_t local_z_start = LOCAL_Z_START(simpar);

	PosType rnx,rny,rnz, xyz,xyy;

	rnx = nx = NX(simpar);
	rny = ny = NY(simpar);
	rnz = nz = NZ(simpar);


	double rincube = 1.L/rnx/rny/rnz;

	DEBUGPRINT("P%d has rincube %g\n", MYID(simpar), rincube);

	ptrdiff_t mx = 2*(nx/2+1);

	xyz = rnx/rnz;
	xyy = rnx/rny;
#ifdef DEBUG
	DEBUGPRINT("-P%d has data %g %g %g %g\n", MYID(simpar), phi[0],phi[1],phi[2345],phi[1235]);
#endif

	if(MYID(simpar)==-1){
		int kkk = 1;
		while(kkk){
			kkk = 2;
		}
	}


	if(ixyz == X) {
		for(k=0; k<local_nz;k++){
			ptrdiff_t kz = xyz*(k + local_z_start);
			if(kz > nz/2) kz = kz - nz;
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
			for(j=0;j<ny;j++){
				ptrdiff_t ky = xyy*j;
				if(ky > ny/2) ky = ky - ny;
				for(i=0;i<=nx/2;i++){
					PosType tmp1 = phi[2*i     + mx*(j + ny*k)];
					PosType tmp2 = phi[2*i + 1 + mx*(j + ny*k)];
					phi[2*i    + mx*(j+ny*k)] = -twopi*i*tmp2;
					phi[2*i +1 + mx*(j+ny*k)] =  twopi*i*tmp1;
				}

			}
		}
	}
	else if(ixyz ==Y) {
		for(k=0;k< local_nz;k++){
			ptrdiff_t kz = xyz*(k + local_z_start);
			if(kz > nz/2) kz = kz - nz;
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
			for(j=0;j<ny;j++){
				ptrdiff_t ky = xyy*j;
				if(ky > ny/2) ky = ky - ny;
				for(i=0;i<=nx/2;i++){
					PosType tmp1 = phi[2*i     + mx*(j + ny*k)];
					PosType tmp2 = phi[2*i + 1 + mx*(j + ny*k)];
					phi[2*i    + mx*(j+ny*k)] = -twopi*ky*tmp2;
					phi[2*i +1 + mx*(j+ny*k)] =  twopi*ky*tmp1;
				}

			}
		}
	}
	else if(ixyz ==Z) {
		for(k=0;k< local_nz;k++){
			ptrdiff_t kz = xyz*(k + local_z_start);
			if(kz > nz/2) kz = kz - nz;
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
			for(j=0;j<ny;j++){
				ptrdiff_t ky = xyy*j;
				if(ky > ny/2) ky = ky - ny;
				for(i=0;i<=nx/2;i++){
					PosType tmp1 = phi[2*i     + mx*(j + ny*k)];
					PosType tmp2 = phi[2*i + 1 + mx*(j + ny*k)];
					phi[2*i    + mx*(j+ny*k)] = -twopi*kz*tmp2;
					phi[2*i +1 + mx*(j+ny*k)] =  twopi*kz*tmp1;
				}

			}
		}
	}
	else {
		DEBUGPRINT("Please check input ixyz flag 1/2/3 :: %d\n", ixyz);
		/*
		mpi_fftw_finalize(&FFTWINFO(simpar));
		*/
		mpi_fftw_finalize();
		exit(99);
	}
	long iseed = -1*(MYID(simpar)+10);


	if(0){
		void dumpphi(SimParameters *, float *);
		dumpphi(simpar,phi);
	}

	fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),(fftwf_complex*)phi, phi);


	if(0){
		void dumpden(SimParameters *, float *, ptrdiff_t);
		dumpden(simpar,phi, mx);
	}


	for(k=0;k<local_nz;k++) for(j=0;j<ny;j++) for(i=0;i<mx;i++) 
			phi[i+nx*(j+ny*k)] = phi[i+mx*(j+ny*k)] * rincube;

#ifdef DEBUG
	/*
	if(0){
		double minval,maxval;
		minval = 1.e20L;
		maxval = -1e20L;
		for(k=0;k<local_nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++){
			minval = MIN(minval, phi[i+mx*(j+ny*k)]);
			maxval = MAX(maxval, phi[i+mx*(j+ny*k)]);
		}
		DEBUGPRINT("P%d has min max value %g %g\n",MYID(simpar),minval,maxval);
	}
	MPI_Finalize();
	exit(99);
	*/
#endif
}

void diffphi(SimParameters *simpar, DenType *phi, DenType *dphidx,enum dimension ixyz){
	float twopi = 2*M_PI;
	float fourpi2 = twopi*twopi;
	PosType tmp1, tmp2;

	ptrdiff_t i,j,k;
	ptrdiff_t nx,ny,nz;
	ptrdiff_t local_nz = LOCAL_NZ(simpar);
	ptrdiff_t local_z_start = LOCAL_Z_START(simpar);

	PosType rnx,rny,rnz, xyz,xyy;

	rnx = nx = NX(simpar);
	rny = ny = NY(simpar);
	rnz = nz = NZ(simpar);
	double rincube = 1.L/rnx/rny/rnz;
	int mx = 2*(nx/2+1);

	xyz = rnx/rnz;
	xyy = rnx/rny;


	if(ixyz == X) {
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
		for(k=0; k<local_nz;k++){
			ptrdiff_t kz = xyz*(k + local_z_start);
			if(kz > nz/2) kz = kz - nz;
			for(j=0;j<ny;j++){
				ptrdiff_t ky = xyy*j;
				if(ky > ny/2) ky = ky - ny;
				for(i=0;i<=nx/2;i++){
					tmp1 = phi[2*i     + mx*(j + ny*k)];
					tmp2 = phi[2*i + 1 + mx*(j + ny*k)];
					dphidx[2*i    + mx*(j+ny*k)] = -twopi*i*tmp2;
					dphidx[2*i +1 + mx*(j+ny*k)] =  twopi*i*tmp1;
				}

			}
		}
	}
	else if(ixyz ==Y) {
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
		for(k=0; k<local_nz;k++){
			ptrdiff_t kz = xyz*(k + local_z_start);
			if(kz > nz/2) kz = kz - nz;
			for(j=0;j<ny;j++){
				ptrdiff_t ky = xyy*j;
				if(ky > ny/2) ky = ky - ny;
				for(i=0;i<=nx/2;i++){
					tmp1 = phi[2*i     + mx*(j + ny*k)];
					tmp2 = phi[2*i + 1 + mx*(j + ny*k)];
					dphidx[2*i    + mx*(j+ny*k)] = -twopi*ky*tmp2;
					dphidx[2*i +1 + mx*(j+ny*k)] =  twopi*ky*tmp1;
				}

			}
		}
	}
	else if(ixyz ==Z) {
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
		for(k=0;k< local_nz;k++){
			ptrdiff_t kz = xyz*(k + local_z_start);
			if(kz > nz/2) kz = kz - nz;
			for(j=0;j<ny;j++){
				ptrdiff_t ky = xyy*j;
				if(ky > ny/2) ky = ky - ny;
				for(i=0;i<=nx/2;i++){
					tmp1 = phi[2*i     + mx*(j + ny*k)];
					tmp2 = phi[2*i + 1 + mx*(j + ny*k)];
					dphidx[2*i    + mx*(j+ny*k)] = -twopi*kz*tmp2;
					dphidx[2*i +1 + mx*(j+ny*k)] =  twopi*kz*tmp1;
				}

			}
		}
	}
	else {
		DEBUGPRINT("Please check input ixyz flag 1/2/3 :: %d\n", ixyz);
		/*
		mpi_fftw_finalize(&FFTWINFO(simpar));
		*/
		mpi_fftw_finalize();
		exit(99);
	}
	fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),(fftwf_complex*)dphidx,dphidx);
	for(k=0;k<local_nz;k++) 
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
		for(j=0;j<ny;j++) for(i=0;i<nx;i++) dphidx[i+mx*(j+ny*k)] *= rincube;
}



void TwoLPT(SimParameters *simpar,  DenType *phi1, DenType *work1, DenType *work2, DenType *phi2){
	double Cvel = 2.99792458E5L;
	float twopi = 2*M_PI;
	float fourpi2 = twopi*twopi;
	float growthfactor = COS_GROWTH(simpar);




	ptrdiff_t i,j,k;
	ptrdiff_t nx,ny,nz;
	ptrdiff_t local_nz = LOCAL_NZ(simpar);
	ptrdiff_t local_z_start = LOCAL_Z_START(simpar);
	double VelPot2RealPot, iVelPot2RealPot;
	float fNL = FNL(simpar);
	float gNL = GNL(simpar);

	PosType rnx,rny,rnz, xyz,xyy;

	rnx = nx = NX(simpar);
	rny = ny = NY(simpar);
	rnz = nz = NZ(simpar);
	ptrdiff_t mx = 2*(nx/2+1);

	xyz = rnx/rnz;
	xyy = rnx/rny;

	double rincube = 1.L/(double)nx/(double)ny/(double)nz;
	double sqrrincube = rincube * rincube;
	double phistd, phimean, tphistd, tphimean;
	double fNLzcorr, gNLzcorr;

	if(MYID(simpar)==0) DEBUGPRINT("P%d: now I am here with fNL/gNL= %g %g\n", MYID(simpar), fNL, gNL);

	if(fNL != 0 || gNL !=0){
		VelPot2RealPot = 3./2.*pow(100*COSBOXSIZE(simpar)/Cvel, 2)*OMEP(simpar)*AMAX(simpar);
		iVelPot2RealPot = 1.L/VelPot2RealPot;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
		for(k=0;k<local_nz;k++) for(j=0;j<ny;j++)for(i=0;i<nx;i++)
			phi1[i+mx*(j+ny*k)] *= VelPot2RealPot;

		fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),(fftwf_complex*)phi1,phi1);
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
		for(k=0;k<local_nz;k++) for(j=0;j<ny;j++)for(i=0;i<nx;i++)
			phi1[i+mx*(j+ny*k)] *= rincube;
		phistd = 0;
		phimean = 0;
		for(k=0;k<local_nz;k++) for(j=0;j<ny;j++)for(i=0;i<nx;i++){
			ptrdiff_t ipos = i+mx*(j+ny*k);
			phistd += phi1[ipos]*phi1[ipos];
			phimean += phi1[ipos];
		}
		MPI_Reduce(&phistd, &tphistd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM(simpar));
		MPI_Bcast(&tphistd, 1, MPI_DOUBLE, 0, MPI_COMM(simpar));
		phistd = tphistd;
		MPI_Reduce(&phimean, &tphimean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM(simpar));
		MPI_Bcast(&tphimean, 1, MPI_DOUBLE, 0, MPI_COMM(simpar));
		phimean = tphimean*rincube;

		phistd = sqrt(phistd*rincube - phimean*phimean);
		fNLzcorr = 1./growthfactor/AMAX(simpar);
		gNLzcorr = pow(1./growthfactor/AMAX(simpar), 2);
		for(k=0;k<local_nz;k++) for(j=0;j<ny;j++)for(i=0;i<nx;i++){
			ptrdiff_t ipos = i+mx*(j+ny*k);
			phi1[ipos] += fNLzcorr*fNL*(pow(phi1[ipos],2)-phistd*phistd)
				+gNLzcorr*gNL*(pow(phi1[ipos],3)-3*phistd*phistd*phi1[ipos]);
		}
		fftwf_mpi_execute_dft_r2c(FFTW_NF_PLAN(simpar),phi1,(fftwf_complex*)phi1);
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
		for(k=0;k<local_nz;k++) for(j=0;j<ny;j++)for(i=0;i<nx;i++)
			phi1[i+mx*(j+ny*k)] *= iVelPot2RealPot;
	}
	/*
	for(i=0;i<mx*ny*local_nz;i++){
		if(isnan(phi1[i]) || isinf(phi1[i])){
			DEBUGPRINT(stderr,"P%d Error in put value %g :: %ld\n",MYID(simpar),phi1[i],i);
			exit(99);
		}
	}
	*/
	/* phi, 11 * phi,22 component */
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
	for(k=0;k<local_nz;k++){
		ptrdiff_t kz = xyz*(k+local_z_start);
		if(kz > nz/2) kz -= nz;
		for(j=0;j<ny;j++){
			ptrdiff_t ky = xyy*j;
			if(ky > ny/2) ky -= ny; 
			ptrdiff_t offset = mx*(j+ny*k);
			for(i=0;i<=nx/2;i++){
				ptrdiff_t ipos = 2*i + offset;
				work1[ipos  ] = -fourpi2*i*i*phi1[ipos];     /* phi,11 */
				work1[ipos+1] = -fourpi2*i*i*phi1[ipos+1];   /* phi,11 */
				work2[ipos  ] = -fourpi2*ky*ky*phi1[ipos];   /* phi,22 */
				work2[ipos+1] = -fourpi2*ky*ky*phi1[ipos+1]; /* phi,22 */
			}
		}
	}
	fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),(fftwf_complex*)work1, work1);
	fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),(fftwf_complex*)work2, work2);
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
	for(k=0;k<local_nz;k++) 
		for(j=0;j<ny;j++)for(i=0;i<nx;i++){
			ptrdiff_t ipos = i + mx*(j+ny*k);
			phi2[ipos] = work1[ipos]*work2[ipos]*sqrrincube;
		}
	/* phi, 11 * phi,33 component */
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
	for(k=0;k<local_nz;k++){
		ptrdiff_t kz = xyz*(k+local_z_start);
		if(kz > nz/2) kz -= nz;
		for(j=0;j<ny;j++){
			ptrdiff_t ky = xyy*j;
			if(ky > ny/2) ky -= ny; 
			ptrdiff_t offset = mx*(j+ny*k);
			for(i=0;i<=nx/2;i++){
				ptrdiff_t ipos = 2*i + offset;
				work2[ipos]   = -fourpi2*kz*kz*phi1[ipos];   /* phi,33 */
				work2[ipos+1] = -fourpi2*kz*kz*phi1[ipos+1]; /* phi,33 */
			}
		}
	}
	fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),(fftwf_complex*)work2, work2);
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
	for(k=0;k<local_nz;k++) 
		for(j=0;j<ny;j++)for(i=0;i<nx;i++){
			ptrdiff_t ipos = i + mx*(j+ny*k);
			phi2[ipos] += work1[ipos]*work2[ipos]*sqrrincube;
		}
	/* phi, 22 * phi,33 component */
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
	for(k=0;k<local_nz;k++){
		ptrdiff_t kz = xyz*(k+local_z_start);
		if(kz > nz/2) kz -= nz;
		for(j=0;j<ny;j++){
			ptrdiff_t ky = xyy*j;
			if(ky > ny/2) ky -= ny; 
			for(i=0;i<=nx/2;i++){
				ptrdiff_t ipos = 2*i + mx*(j+ny*k);
				work1[ipos]   = -fourpi2*ky*ky*phi1[ipos];   /* phi,22 */
				work1[ipos+1] = -fourpi2*ky*ky*phi1[ipos+1]; /* phi,22 */
			}
		}
	}
	fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),(fftwf_complex*)work1, work1);
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
	for(k=0;k<local_nz;k++) 
		for(j=0;j<ny;j++)for(i=0;i<nx;i++){
		ptrdiff_t ipos = i + mx*(j+ny*k);
		phi2[ipos] += work1[ipos]*work2[ipos]*sqrrincube;
	}
	/* phi, 21 * phi,21 component */
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
	for(k=0;k<local_nz;k++){
		ptrdiff_t kz = xyz*(k+local_z_start);
		if(kz > nz/2) kz -= nz;
		for(j=0;j<ny;j++){
			ptrdiff_t ky = xyy*j;
			if(ky > ny/2) ky -= ny; 
			for(i=0;i<=nx/2;i++){
				ptrdiff_t ipos = 2*i + mx*(j+ny*k);
				work1[ipos]   = -fourpi2*ky*i*phi1[ipos];   /* phi,21 */
				work1[ipos+1] = -fourpi2*ky*i*phi1[ipos+1]; /* phi,21 */
			}
		}
	}
	fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),(fftwf_complex*)work1, work1);
	for(k=0;k<local_nz;k++) 
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
		for(j=0;j<ny;j++)for(i=0;i<nx;i++){
			ptrdiff_t ipos = i + mx*(j+ny*k);
			phi2[ipos] -= work1[ipos]*work1[ipos]*sqrrincube;
		}
	/* phi, 31 * phi,31 component */
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
	for(k=0;k<local_nz;k++){
		ptrdiff_t kz = xyz*(k+local_z_start);
		if(kz > nz/2) kz -= nz;
		for(j=0;j<ny;j++){
			ptrdiff_t ky = xyy*j;
			if(ky > ny/2) ky -= ny; 
			for(i=0;i<=nx/2;i++){
				ptrdiff_t ipos = 2*i + mx*(j+ny*k);
				work1[ipos]   = -fourpi2*kz*i*phi1[ipos];   /* phi,31 */
				work1[ipos+1] = -fourpi2*kz*i*phi1[ipos+1]; /* phi,31 */
			}
		}
	}
	fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),(fftwf_complex*)work1, work1);
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
	for(k=0;k<local_nz;k++) 
		for(j=0;j<ny;j++)for(i=0;i<nx;i++){
			ptrdiff_t ipos = i + mx*(j+ny*k);
			phi2[ipos] -= work1[ipos]*work1[ipos]*sqrrincube;
		}
	/* phi, 32 * phi,32 component */
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
	for(k=0;k<local_nz;k++){
		ptrdiff_t kz = xyz*(k+local_z_start);
		if(kz > nz/2) kz -= nz;
		for(j=0;j<ny;j++){
			ptrdiff_t ky = xyy*j;
			if(ky > ny/2) ky -= ny; 
			for(i=0;i<=nx/2;i++){
				ptrdiff_t ipos = 2*i + mx*(j+ny*k);
				work1[ipos]   = -fourpi2*kz*ky*phi1[ipos];   /* phi,32 */
				work1[ipos+1] = -fourpi2*kz*ky*phi1[ipos+1]; /* phi,32 */
			}
		}
	}
	fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),(fftwf_complex*)work1, work1);
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
	for(k=0;k<local_nz;k++) 
		for(j=0;j<ny;j++)for(i=0;i<nx;i++){
			ptrdiff_t ipos = i + mx*(j+ny*k);
			phi2[ipos] -= work1[ipos]*work1[ipos]*sqrrincube;
		}

	/* go back to the Fourier space */
	fftwf_mpi_execute_dft_r2c(FFTW_NF_PLAN(simpar),phi2, (fftwf_complex*)phi2);

#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
	for(k=0;k<local_nz;k++){
		ptrdiff_t kz = xyz*(k+local_z_start);
		if(kz > nz/2) kz -= nz;
		for(j=0;j<ny;j++){
			ptrdiff_t ky = xyy*j;
			if(ky > ny/2) ky -= ny; 
			ptrdiff_t ipos = mx*(j+ny*k);
			for(i=0;i<=nx/2;i++){
				double wave2 = -(i*i + ky*ky + kz*kz)*fourpi2;
				double invwave2 = 1.L/wave2;
				if(wave2 < 0.L) {
					phi2[2*i   + ipos] *= invwave2;
					phi2[2*i+1 + ipos] *= invwave2;
				}
				else {
					phi2[2*i   + ipos] = 0;
					phi2[2*i+1 + ipos] = 0;
				}
			}
		}
	}

}
