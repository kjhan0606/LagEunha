#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
#include "fft.h"
/*
#include "mpirms.h"
*/
#include "pmmain.h"
#include "cosmology.h"
#include "CosmosEvolFactor.h"
#include "denrw.h"
#include "timerutil.h"


#define MIN(a,b) ( (a)<(b)?(a):(b))
#define MAX(a,b) ( (a)>(b)?(a):(b))



void SetEntireSimBox(SimParameters *simpar){
	COS_SIMBOX(simpar).x.min = 0;
	COS_SIMBOX(simpar).y.min = 0;
	COS_SIMBOX(simpar).z.min = 0;
	COS_SIMBOX(simpar).x.max = NX(simpar);
	COS_SIMBOX(simpar).y.max = NY(simpar);
	COS_SIMBOX(simpar).z.max = NZ(simpar);
}
/*
void UpdateDdinfoBox(SimParameters *simpar){
	int nddinfo = NDDINFO(simpar);
	int i = nddinfo -1;
	for(i=nddinfo-1;i>0;i--){
		if(i==nddinfo-1){
			DM_DDINFO(simpar)[i].lgroup.xyz.xmin = COS_SIMBOX(simpar).x.min;
			DM_DDINFO(simpar)[i].lgroup.xyz.ymin = COS_SIMBOX(simpar).y.min;
			DM_DDINFO(simpar)[i].lgroup.xyz.zmin = COS_SIMBOX(simpar).z.min;
			DM_DDINFO(simpar)[i].lgroup.xyz.xmax = COS_SIMBOX(simpar).x.max;
			DM_DDINFO(simpar)[i].lgroup.xyz.ymax = COS_SIMBOX(simpar).y.max;
			DM_DDINFO(simpar)[i].lgroup.xyz.zmax = COS_SIMBOX(simpar).z.max;
		}
		{
			PosType xmin,ymin,zmin,xmax,ymax,zmax;
			MPI_Comm com = DM_DDINFO(simpar)[i].com;
			MPI_Reduce(&xmin, &(DM_DDINFO(simpar)[i].lgroup.xyz.xmin), 1, MPI_POSTYPE, MPI_MIN, 0, com);
			MPI_Reduce(&ymin, &(DM_DDINFO(simpar)[i].lgroup.xyz.ymin), 1, MPI_POSTYPE, MPI_MIN, 0, com);
			MPI_Reduce(&zmin, &(DM_DDINFO(simpar)[i].lgroup.xyz.zmin), 1, MPI_POSTYPE, MPI_MIN, 0, com);
			MPI_Reduce(&xmax, &(DM_DDINFO(simpar)[i].lgroup.xyz.xmax), 1, MPI_POSTYPE, MPI_MAX, 0, com);
			MPI_Reduce(&ymax, &(DM_DDINFO(simpar)[i].lgroup.xyz.ymax), 1, MPI_POSTYPE, MPI_MAX, 0, com);
			MPI_Reduce(&zmax, &(DM_DDINFO(simpar)[i].lgroup.xyz.zmax), 1, MPI_POSTYPE, MPI_MAX, 0, com);
			MPI_Bcast(&xmin, 1, MPI_POSTYPE, 0, com);
			MPI_Bcast(&ymin, 1, MPI_POSTYPE, 0, com);
			MPI_Bcast(&zmin, 1, MPI_POSTYPE, 0, com);
			MPI_Bcast(&xmax, 1, MPI_POSTYPE, 0, com);
			MPI_Bcast(&ymax, 1, MPI_POSTYPE, 0, com);
			MPI_Bcast(&zmax, 1, MPI_POSTYPE, 0, com);
			DM_DDINFO(simpar)[i-1].lgroup.xyz.xmin = xmin;
			DM_DDINFO(simpar)[i-1].lgroup.xyz.ymin = ymin;
			DM_DDINFO(simpar)[i-1].lgroup.xyz.zmin = zmin;
			DM_DDINFO(simpar)[i-1].lgroup.xyz.xmax = xmax;
			DM_DDINFO(simpar)[i-1].lgroup.xyz.ymax = ymax;
			DM_DDINFO(simpar)[i-1].lgroup.xyz.zmax = zmax;
		}
		STAR_DDINFO(simpar)[i].lgroup = AGN_DDINFO(simpar)[i].lgroup = SPH_DDINFO(simpar)[i].lgroup = DM_DDINFO(simpar)[i].lgroup;
		TSTAR_DDINFO(simpar)[i].lgroup = TAGN_DDINFO(simpar)[i].lgroup = TSPH_DDINFO(simpar)[i].lgroup = TDM_DDINFO(simpar)[i].lgroup
			= DM_DDINFO(simpar)[i].lgroup;
	}
}
*/

void *MakeDenGrid4TSC(SimParameters *simpar, GridInfo *dengrid, ptrdiff_t bbuffer, ptrdiff_t ubuffer){
	dengrid->nx = NX(simpar);
	dengrid->ny = NY(simpar);
	dengrid->nz = NZ(simpar);
	dengrid->ix = rint(SIM_LXMIN(simpar,dm))   -bbuffer;
	dengrid->iy = rint(SIM_LYMIN(simpar,dm))   -bbuffer;
	dengrid->iz = rint(SIM_LZMIN(simpar,dm))   -bbuffer;
	dengrid->jx = rint(SIM_LXMAX(simpar,dm))-1 +ubuffer;
	dengrid->jy = rint(SIM_LYMAX(simpar,dm))-1 +ubuffer;
	dengrid->jz = rint(SIM_LZMAX(simpar,dm))-1 +ubuffer;

#ifdef DEBUG
	if(1){
		PosType xmin,ymin,zmin,xmax,ymax,zmax;
		xmin = ymin = zmin = 1.E20;
		xmax = ymax = zmax = -1.E20;
		ptrdiff_t i,j,k;
		for(i=0;i<DM_NP(simpar);i++){
			xmin = MIN(xmin, XofP(simpar, DM_BP(simpar)+i));
			xmax = MAX(xmax, XofP(simpar, DM_BP(simpar)+i));
			ymin = MIN(ymin, YofP(simpar, DM_BP(simpar)+i));
			ymax = MAX(ymax, YofP(simpar, DM_BP(simpar)+i));
			zmin = MIN(zmin, ZofP(simpar, DM_BP(simpar)+i));
			zmax = MAX(zmax, ZofP(simpar, DM_BP(simpar)+i));
		}
		DEBUGPRINT("P%d has pregion %g %g : %g %g : %g %g\n",MYID(simpar),xmin,xmax, ymin,ymax,zmin,zmax);
		DEBUGPRINT("P%d has region %g %g : %g %g : %g %g\n",MYID(simpar), 
				SIM_LXMIN(simpar,dm), SIM_LXMAX(simpar,dm),
				SIM_LYMIN(simpar,dm), SIM_LYMAX(simpar,dm),
				SIM_LZMIN(simpar,dm), SIM_LZMAX(simpar,dm));
	}
#endif 
	dengrid->npix = (dengrid->jx-dengrid->ix+1)*(dengrid->jy-dengrid->iy+1)*(dengrid->jz-dengrid->iz+1);
	dengrid = (GridInfo*)realloc(dengrid,sizeof(GridInfo)+sizeof(DenType)*dengrid->npix);
	return dengrid;
}
#ifdef _OPENMP
void TSC(SimParameters *simpar, GridInfo *dengrid){
	DenType *den = (DenType*)(dengrid + 1);
	ptrdiff_t nx,ny,nz,ngrids,j;
	PosType xstart,ystart, zstart;
	nx = dengrid->jx - dengrid->ix +1;
	ny = dengrid->jy - dengrid->iy +1;
	nz = dengrid->jz - dengrid->iz +1;
	ngrids = dengrid->npix;

	xstart = dengrid->ix;
	ystart = dengrid->iy;
	zstart = dengrid->iz;

	for(j=0;j<ngrids;j++) den[j] = 0;

	DenType *threadDen;
	int nthreads;
	PosType zdomainsize = SIM_LZMAX(simpar,dm) - SIM_LZMIN(simpar,dm);
	/*
	if(MYID(simpar)==0) {
		int kkk = 1;
		while(kkk) {
			kkk = 1;
		}
	}
	MPI_Barrier(MPI_COMM(simpar));
	*/


#pragma omp parallel 
	{
#pragma omp master
		nthreads = omp_get_num_threads();
	}

#pragma omp parallel 
	{
		int ithread = omp_get_thread_num();
		PosType zslicewidth = zdomainsize/nthreads;
		PosType zstart_thread, zend_thread;
		zstart_thread = SIM_LZMIN(simpar,dm) + ithread*zslicewidth;
		zend_thread = SIM_LZMIN(simpar,dm) + (ithread+1)*zslicewidth;
		if(ithread == nthreads-1) zend_thread = SIM_LZMAX(simpar,dm);
		int iiz,ijz;
		iiz = rint(zstart_thread) - 1;
		ijz = rint(zend_thread)   -1  +2;
#ifdef DEBUG
		DEBUGPRINT("P%.2d.%.2d has variables iiz/ijz= %d %d\n", MYID(simpar),ithread, iiz,ijz);
#endif
		DenType *denthread = (DenType*)calloc(nx *ny *(ijz-iiz+1), sizeof(DenType));
		if(denthread == NULL) {
			DEBUGPRINT("Error in allocating denthread: %p\n", denthread);
		}
		PosType zzstart = iiz;
		if(DM_NP(simpar)>0) OMP_Get_TSC_Den(simpar, DM, denthread, xstart,ystart,zzstart, zstart_thread,zend_thread);
		if(SPH_NP(simpar)>0) OMP_Get_TSC_Den(simpar, SPH, denthread, xstart,ystart,zzstart,zstart_thread,zend_thread);
		if(VORO_NP(simpar)>0) OMP_Get_TSC_Den(simpar, VORO, denthread, xstart,ystart,zzstart,zstart_thread,zend_thread);
		if(STAR_NP(simpar)>0) OMP_Get_TSC_Den(simpar, STAR, denthread, xstart,ystart,zzstart,zstart_thread,zend_thread);
		if(AGN_NP(simpar)>0) OMP_Get_TSC_Den(simpar, AGN,  denthread, xstart,ystart,zzstart,zstart_thread,zend_thread);
#ifdef DEBUG
		DEBUGPRINT("P%.2d.%.2d exit with variables iiz/ijz= %d %d\n", MYID(simpar),ithread, iiz,ijz);
#endif
#pragma omp critical
		{
			size_t i,j,k;
			size_t zoffset = iiz-dengrid->iz;
			for(k=0;k<(ijz-iiz+1);k++) for(i=0;i<nx*ny;i++){
				den[i+nx*ny*(k+zoffset)] += denthread[i+nx*ny*k];
			}

		}
		free(denthread);
	}
}
#else
/* This routine assumes all species have the same domain */
void TSC(SimParameters *simpar, GridInfo *dengrid){
	DenType *den = (DenType*)(dengrid + 1);
	ptrdiff_t nx,ny,nz,ngrids,j;
	PosType xstart,ystart, zstart;
	nx = dengrid->jx - dengrid->ix +1;
	ny = dengrid->jy - dengrid->iy +1;
	nz = dengrid->jz - dengrid->iz +1;
	ngrids = dengrid->npix;

	xstart = dengrid->ix;
	ystart = dengrid->iy;
	zstart = dengrid->iz;
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));
	if(1){
		PosType xmin,ymin,zmin,xmax,ymax,zmax;
		xmin = ymin = zmin = 1.E20;
		xmax = ymax = zmax = -1.E20;
		ptrdiff_t i,j,k;
		for(i=0;i<DM_NP(simpar);i++){
			/*
			xmin = MIN(xmin, XofP(simpar, DM_BP(simpar)+i));
			xmax = MAX(xmax, XofP(simpar, DM_BP(simpar)+i));
			*/
			xmin = MIN(xmin, (DM_BP(simpar)+i)->x);
			xmax = MAX(xmax, (DM_BP(simpar)+i)->x);
			ymin = MIN(ymin, YofP(simpar, DM_BP(simpar)+i));
			ymax = MAX(ymax, YofP(simpar, DM_BP(simpar)+i));
			zmin = MIN(zmin, ZofP(simpar, DM_BP(simpar)+i));
			zmax = MAX(zmax, ZofP(simpar, DM_BP(simpar)+i));
		}
		DEBUGPRINT("P%d has pregion %g %g : %g %g : %g %g for %ld %ld %ld\n",MYID(simpar),xmin,xmax, ymin,ymax,zmin,zmax, NX(simpar),NY(simpar),NZ(simpar));
		DEBUGPRINT("P%d has region %g %g : %g %g : %g %g\n",MYID(simpar), 
				SIM_LXMIN(simpar,dm), SIM_LXMAX(simpar,dm),
				SIM_LYMIN(simpar,dm), SIM_LYMAX(simpar,dm),
				SIM_LZMIN(simpar,dm), SIM_LZMAX(simpar,dm));
	}
#endif 

#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));
	DEBUGPRINT("P%d is now inside TSC: %g %g %g ::: nx/ny/nz= %ld %ld %ld\n", MYID(simpar), xstart,ystart,zstart, nx,ny,nz);
#endif

#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));
	DEBUGPRINT("P%d has dmnp = %ld sphnp= %ld: xyzstart: %g %g %g\n",MYID(simpar), DM_NP(simpar), SPH_NP(simpar), xstart,ystart,zstart);
	{
		PosType xmin,ymin,zmin,xmax,ymax,zmax;
		xmin = ymin = zmin = 1.e20;
		xmax = ymax = zmax =-1.e20;
		size_t i;
		for(i=0;i<DM_NP(simpar);i++){
			xmin = MIN(xmin, XofP(simpar,DM_BP(simpar)+i));
			xmax = MAX(xmax, XofP(simpar,DM_BP(simpar)+i));
			ymin = MIN(ymin, YofP(simpar,DM_BP(simpar)+i));
			ymax = MAX(ymax, YofP(simpar,DM_BP(simpar)+i));
			zmin = MIN(zmin, ZofP(simpar,DM_BP(simpar)+i));
			zmax = MAX(zmax, ZofP(simpar,DM_BP(simpar)+i));
		}
		DEBUGPRINT("P%d has %g %g %g ::: %g %g %g\n",MYID(simpar),xmin,ymin,zmin,xmax,ymax,zmax);

	}
#endif
	for(j=0;j<ngrids;j++) den[j] = 0;

	if(DM_NP(simpar)>0) Get_TSC_Den(simpar, DM, den, xstart,ystart,zstart);
	if(SPH_NP(simpar)>0) Get_TSC_Den(simpar, SPH, den, xstart,ystart,zstart);
	if(VORO_NP(simpar)>0) Get_TSC_Den(simpar, VORO, den, xstart,ystart,zstart);
	if(STAR_NP(simpar)>0) Get_TSC_Den(simpar, STAR, den, xstart,ystart,zstart);
	if(AGN_NP(simpar)>0) Get_TSC_Den(simpar, AGN,  den, xstart,ystart,zstart);
}

#endif

GridInfo *ReallocDen4FFTW_Format(SimParameters *simpar, GridInfo *fftwgrid){
	ptrdiff_t nx,ny,nz;
	nx = NX(simpar);
	ny = NY(simpar);
	nz = LOCAL_NZ(simpar);
	ptrdiff_t mx = 2*(nx/2+1);

	fftwgrid = (GridInfo*)realloc(fftwgrid, sizeof(GridInfo)+sizeof(fftwf_complex)*LOCAL_GRID_SIZE(simpar));
	DenType *den = (DenType*) ( (fftwgrid+1));
	ptrdiff_t i,j,k;
	for(k=nz-1;k>=0;k--) for(j=ny-1;j>=0;j--){
		DenType tmp[nx];
		memmove(tmp, den+nx*(j+ny*k), sizeof(DenType)*nx);
		memmove(den+mx*(j+ny*k), tmp, sizeof(DenType)*nx);
	}

	return fftwgrid;
}


GridInfo *MakeGrids4FFTWDD(SimParameters *simpar, GridInfo *dengrid){
	ptrdiff_t nx,ny,nz,mx;
	GridInfo *fftwgrid = (GridInfo*)malloc(sizeof(GridInfo)+sizeof(fftwf_complex)*LOCAL_GRID_SIZE(simpar));
	fftwgrid->nx = (nx = NX(simpar));
	fftwgrid->ny = (ny = NY(simpar));
	fftwgrid->nz = (nz = NZ(simpar));
	fftwgrid->nxny = (nx *ny);
	fftwgrid->ix = fftwgrid->iy = 0;
	fftwgrid->jx = NX(simpar)-1;
	fftwgrid->jy = NY(simpar)-1;
	fftwgrid->npix = 2*LOCAL_GRID_SIZE(simpar);
	if(LOCAL_NZ(simpar) >0) {
		fftwgrid->iz = LOCAL_Z_START(simpar);
		fftwgrid->jz = LOCAL_Z_START(simpar) + LOCAL_NZ(simpar)-1;
	}
	else {
		fftwgrid->iz = 10*nz;
		fftwgrid->jz = -10*nz;
	}
	DenType *den = (DenType*)(fftwgrid + 1);
	ptrdiff_t i,j,k;
	for(i=0;i<2*LOCAL_GRID_SIZE(simpar);i++) den[i] = -1;
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));
	DEBUGPRINT("P%d is now before gmigrate %d %d %d :: %d %d %d :: %ld\n",
			MYID(simpar), fftwgrid->ix, fftwgrid->iy, fftwgrid->iz,
			fftwgrid->jx, fftwgrid->jy, fftwgrid->jz, fftwgrid->npix);
#endif
	gmigrate(dengrid, fftwgrid,DM_DDINFO(simpar), NDDINFO(simpar));
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));
	DEBUGPRINT("+P%d is now after gmigrate %d %d %d :: %d %d %d :: %ld\n",
			MYID(simpar), fftwgrid->ix, fftwgrid->iy, fftwgrid->iz,
			fftwgrid->jx, fftwgrid->jy, fftwgrid->jz, fftwgrid->npix);
#endif
	fftwgrid->npix = 2*LOCAL_GRID_SIZE(simpar);
	den = (DenType*)(fftwgrid + 1);
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));
	DEBUGPRINT("P%d is now passing throught the den mesh arrangement.\n",MYID(simpar));
#endif
			

	if( (PMSTATUS(simpar) == PUSH && STEPCOUNT(simpar) == 1) || PMSTATUS(simpar) == PULL) volumeout(simpar, fftwgrid, xzslice);

	mx = 2*(nx/2+1);
	for(k=LOCAL_NZ(simpar)-1;k>=0;k--) for(j=ny-1;j>=0;j--){
		DenType tmp[nx];
		memmove(tmp, den+nx*(j+ny*k), sizeof(DenType)*nx);
		memmove( den+mx*(j+ny*k), tmp,sizeof(DenType)*nx);

	}
	return fftwgrid;
}
GridInfo *MakeDen4FDA(SimParameters *simpar, GridInfo *fftwgrid){
	ptrdiff_t i,j,k;
	ptrdiff_t nx,ny,nz,mx;
	nx = NX(simpar);
	ny = NY(simpar);
	nz = NZ(simpar);
	mx = 2*(nx/2+1);
	GridInfo *ingrid,  *fda4grid;

	DenType *den = (DenType*)(fftwgrid+1);
	for(k=0;k<LOCAL_NZ(simpar);k++) for(j=0;j<NY(simpar);j++){
		DenType tmp[nx];
		memmove(tmp, den+mx*(j+ny*k), sizeof(DenType)*nx);
		memmove(den+nx*(j+ny*k), tmp, sizeof(DenType)*nx);
	}


	den = (DenType*)(fftwgrid+1);
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));DEBUGPRINT("P%d has %d %d : %d %d : %d %d : %g %g %g %g\n",MYID(simpar),
			fftwgrid->ix, fftwgrid->jx, fftwgrid->iy, fftwgrid->jy, fftwgrid->iz, fftwgrid->jz,
			den[fftwgrid->npix/4],den[fftwgrid->npix/3],den[fftwgrid->npix/2], den[fftwgrid->npix-1]);
#endif
	
	fda4grid = malloc(sizeof(GridInfo));
	fda4grid->nx = NX(simpar);
	fda4grid->ny = NY(simpar);
	fda4grid->nz = NZ(simpar);
	fda4grid->ix = rint(SIM_LXMIN(simpar,dm))  -3;
	fda4grid->iy = rint(SIM_LYMIN(simpar,dm))  -3;
	fda4grid->iz = rint(SIM_LZMIN(simpar,dm))  -3;
	fda4grid->jx = rint(SIM_LXMAX(simpar,dm))-1+4;
	fda4grid->jy = rint(SIM_LYMAX(simpar,dm))-1+4;
	fda4grid->jz = rint(SIM_LZMAX(simpar,dm))-1+4;
	fda4grid->npix = (ptrdiff_t)(fda4grid->jx-fda4grid->ix+1)*(ptrdiff_t)(fda4grid->jy-fda4grid->iy+1)*(ptrdiff_t)(fda4grid->jz-fda4grid->iz+1);
	fda4grid = realloc(fda4grid, sizeof(GridInfo)+sizeof(DenType)* fda4grid->npix);
	den = (DenType*)(fda4grid+1);
	for(i=0;i<fda4grid->npix;i++) den[i] = 0;

	den = (DenType*)(fftwgrid+1);
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));DEBUGPRINT("P%d has %d %d : %d %d : %d %d : %g %g %g %g\n",MYID(simpar),
			fftwgrid->ix, fftwgrid->jx, fftwgrid->iy, fftwgrid->jy, fftwgrid->iz, fftwgrid->jz,
			den[fftwgrid->npix/4],den[fftwgrid->npix/3],den[fftwgrid->npix/2], den[fftwgrid->npix-1]);
#endif
	gmigrate(fftwgrid, fda4grid, DM_DDINFO(simpar), NDDINFO(simpar));
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));DEBUGPRINT("P%d is now after gmigrate \n",MYID(simpar));
#endif

	return fda4grid;
}

void pmmain(SimParameters *simpar, EvolFact evolfact){
	DoDeInfo *ddinfo;
	int nddinfo;
	ptrdiff_t i;
	DenType *den;

	/*
	float fact1, fact2;
	if(PMSTATUS(simpar) == PUSH){
		fact1 = evolfact.fact1_push;
		fact2 = evolfact.fact2_push;
		float a1, a2;
		a1 = ANOW(simpar); a2 = a1 + ASTEP(simpar)/2;
		fact2 = fact2/GetFact1(a1,a1);
	}
	else if (PMSTATUS(simpar) == PULL){
		fact1 = evolfact.fact1_pull;
		fact2 = evolfact.fact2_pull;
	}
	else {
		fact1 = evolfact.fact1;
		fact2 = evolfact.fact2;
	}
	*/

	/*
	GridInfo *dengrid, *fdagrid;
	ptrdiff_t ngrids;
	*/
	float cputime0[1];
	float cputime1[1];

#ifdef DEBUG
	if(0){
		void dumppart(SimParameters *);
		dumppart(simpar);
	}
#endif


	/*
	UpdateSimBox(simpar);
	*/
	GridInfo *dengrid = (GridInfo*)malloc(sizeof(GridInfo));

	TIMER_START(0);
	dengrid = (GridInfo *)MakeDenGrid4TSC(simpar, dengrid, 1, 2);
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));
	DEBUGPRINT("P%d is now before TSC\n",MYID(simpar)); 
#endif
	TSC(simpar, dengrid);
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));
	DEBUGPRINT("P%d is now passing through TSC\n",MYID(simpar)); 
#endif
	GridInfo *fftwgrid = MakeGrids4FFTWDD(simpar, dengrid);
	TIMER_STOP(0);
	if(MYID(simpar) ==0) printf("--- Step %d CPU(TSC)= %g\n",STEPCOUNT(simpar), ELAPSED_TIME(0));


	den = (DenType*)(fftwgrid+1);

	if(PMSTATUS(simpar) == PULL || PMSTATUS(simpar) == KICK){
		void savexzslice(SimParameters *, DenType *);
		savexzslice(simpar,den);
	}

	/*
	void dumpden(SimParameters *, float *); 
	dumpden(simpar,den);
	MPI_Finalize();
	exit(99);
	*/



#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));DEBUGPRINT("P%d after FFT C2R %g %g %g %g\n",
			MYID(simpar),den[fftwgrid->npix/2],den[fftwgrid->npix/3],den[fftwgrid->npix/4],den[fftwgrid->npix-3]);
#endif
	/*
	fftwf_mpi_execute_dft_r2c(FFTW_F_PLAN(simpar), (DenType*)(fftwgrid+1), (fftwf_complex*)(fftwgrid+1));
	*/
	TIMER_START(0);
	fftwf_mpi_execute_dft_r2c(FFTW_F_PLAN(simpar), den, (fftwf_complex*)den);
	TIMER_STOP(0);
	if(MYID(simpar) ==0) printf("--- Step %d CPU(F-FFT)= %g\n",STEPCOUNT(simpar), ELAPSED_TIME(0));


	void CorrMeasure(SimParameters *, GridInfo *); 
	if( (PMSTATUS(simpar) == PUSH && STEPCOUNT(simpar) == 1) || PMSTATUS(simpar) == PULL || PMSTATUS(simpar)== KICK) CorrMeasure(simpar, fftwgrid);


	void psolver(SimParameters *, GridInfo*);
	TIMER_START(0);
	psolver(simpar, fftwgrid);
	TIMER_STOP(0);
	if(MYID(simpar) ==0) printf("--- Step %d CPU(psolv)= %g\n",STEPCOUNT(simpar), ELAPSED_TIME(0));

	TIMER_START(0);

	fftwf_mpi_execute_dft_c2r(FFTW_B_PLAN(simpar), (fftwf_complex*)den,den);

	TIMER_STOP(0);
	if(MYID(simpar) ==0) printf("--- Step %d CPU(B-FFT)= %g\n",STEPCOUNT(simpar), ELAPSED_TIME(0));
	/*
	den = (DenType*)(fftwgrid+1); 
	double rincube = 1.L/(NX(simpar)*NY(simpar)*NZ(simpar));
	for(i=0;i<fftwgrid->npix;i++){
		den[i] = den[i]*rincube;
	}
	den = (DenType*)(fftwgrid+1); 
	*/
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));DEBUGPRINT("P%d after FFT C2R %g %g %g %g\n", MYID(simpar),den[0],den[1],den[2],den[3]);
#endif

	GridInfo *fdagrid = MakeDen4FDA(simpar, fftwgrid);
	den = (DenType*)(fdagrid+1);
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));DEBUGPRINT("P%d after FFT C2R %g %g %g %g\n",
			MYID(simpar),den[fdagrid->npix/2],den[fdagrid->npix/3],den[fdagrid->npix/4],den[fdagrid->npix-1]);
#endif

	/*
	void fda4(SimParameters *, GridInfo *, float , float );
	*/
	TIMER_START(0);
	void fda4(SimParameters *, GridInfo *, EvolFact);
	fda4(simpar, fdagrid, evolfact);
	TIMER_STOP(0);
	if(MYID(simpar) ==0) printf("--- Step %d CPU(4FDA)= %g\n",STEPCOUNT(simpar), ELAPSED_TIME(0));
#ifdef DEBUG
	DEBUGPRINT("P%d is now at pmmain\n",MYID(simpar));
#endif
	free(fdagrid);
}
