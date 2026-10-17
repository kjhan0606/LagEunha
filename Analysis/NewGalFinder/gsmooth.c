#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <omp.h>
#ifdef FFTW2
#include<srfftw.h>
#include<sfftw.h>
#else
#include <fftw3.h>
#endif
#include "ramses.h"


void gaussian_Smoothing(float *denGrid,int nx,int ny,int nz, 
		double cellsize, float RG){
	long i,j,k,mx,mxh;
	mx = 2*(nx/2+1);
	mxh = (nx/2+1);

	long ncells =  mx*ny*nz;

	double scale;
	double NORM;
	int nthreads=1;

	scale = 1.L/((double)nx*(double)ny*(double)nz);

#ifdef FFTW2
	fftw_real *a = (fftw_real*)fftw_malloc(sizeof(fftw_real)*ncells);
	fftw_complex *A = (fftw_complex*)a;
	rfftwnd_plan p, pinv;
	p = rfftw3d_create_plan(nz,ny,nx,FFTW_REAL_TO_COMPLEX,
			FFTW_ESTIMATE | FFTW_IN_PLACE);
	pinv = rfftw3d_create_plan(nz, ny, nx, FFTW_COMPLEX_TO_REAL,
                                FFTW_ESTIMATE | FFTW_IN_PLACE);
#else
#pragma omp parallel
	{
#pragma omp critical
		{
			nthreads = omp_get_num_threads();
		}

	}
	fftwf_plan_with_nthreads(nthreads);

	float *a = (float *)fftwf_malloc(sizeof(float )*ncells);
	fftwf_complex *A = (fftwf_complex*)a;
	fftwf_plan p, pinv;
	p = fftwf_plan_dft_r2c_3d(nz, ny, nx, a, A,  FFTW_ESTIMATE);
	pinv = fftwf_plan_dft_c2r_3d(nz, ny, nx, A, a,  FFTW_ESTIMATE);
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<ncells;i++) a[i] = denGrid[i];

	LOGPRINT("just before forward fftw for nx= %d %d %d with nthreads= %d\n", nx,ny,nz, nthreads);
#ifdef FFTWV2
	rfftwnd_one_real_to_complex(p,a,NULL);
#else
	fftwf_execute(p);
#endif
	LOGPRINT("just after forward fftw for nx= %d %d %d\n", nx,ny,nz);

	RG = RG/cellsize;

	double expfac = -(M_PI*RG)*(M_PI*RG);
	int nzh = nz/2;
	int nyh = ny/2;
	int nxh = nx/2;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
	for(k=0;k<nz;k++){
		long k1=k;
		if(k > nzh) k1 = k-nz;
		double zfact = exp(k1*k1*expfac/(double)(nz*nz));
		for(j=0;j<ny;j++){
			long j1 = j;
			if(j > nyh) j1 = j-ny;
			double yfact = exp(j1*j1*expfac/(double)(ny*ny));
			for(i=0;i<=nx/2;i++){
				long i2 = i*2;
				double xfact = exp(i*i*expfac/(double)(nx*nx));
				double factor = xfact*yfact*zfact;
				a[i2+mx*(j+ny*k)] *= factor;
				a[i2+1+mx*(j+ny*k)] *= factor;
			}
		}
	}

#ifdef FFTWV2
	rfftwnd_one_complex_to_real(pinv,A, NULL);
#else
	fftwf_execute(pinv);
#endif


#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<ncells;i++) denGrid[i] = a[i]*scale;
#ifdef FFTWV2
	fftw_free(a);
	rfftwnd_destroy_plan(p);
	rfftwnd_destroy_plan(pinv);
#else
	fftwf_free(a);
	fftwf_destroy_plan(p);
	fftwf_destroy_plan(pinv);
#endif
	DEBUGPRINT("Now exit the Gaussian Smoothing Job\n");
}
