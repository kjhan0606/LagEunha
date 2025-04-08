#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "mpirms.h"
#include "pmmain.h"
#include "cosmology.h"
#include "CosmosEvolFactor.h"


#define MIN(a,b) ( (a)<(b)?(a):(b))
#define MAX(a,b) ( (a)>(b)?(a):(b))



void UpdateSimBox(SimParameters *simpar){
	int nddinfo = NDDINFO(simpar);
	DoDeInfo *ddinfo = DM_DDINFO(simpar) + (nddinfo-1);
	COS_SIMBOX(simpar).x.min = ddinfo->lgroup.xyz.xmin;
	COS_SIMBOX(simpar).y.min = ddinfo->lgroup.xyz.ymin;
	COS_SIMBOX(simpar).z.min = ddinfo->lgroup.xyz.zmin;
	COS_SIMBOX(simpar).x.max = ddinfo->lgroup.xyz.xmax;
	COS_SIMBOX(simpar).y.max = ddinfo->lgroup.xyz.ymax;
	COS_SIMBOX(simpar).z.max = ddinfo->lgroup.xyz.zmax;
}

void *MakeDenGrid4TSC(SimParameters *simpar, GridInfo *dengrid, ptrdiff_t bbuffer, ptrdiff_t ubuffer){
	dengrid->ix = rint(COS_SIMBOX(simpar).x.min)-bbuffer;
	dengrid->iy = rint(COS_SIMBOX(simpar).y.min)-bbuffer;
	dengrid->iz = rint(COS_SIMBOX(simpar).z.min)-bbuffer;
	dengrid->jx = rint(COS_SIMBOX(simpar).x.max)+ubuffer;
	dengrid->jy = rint(COS_SIMBOX(simpar).y.max)+ubuffer;
	dengrid->jz = rint(COS_SIMBOX(simpar).z.max)+ubuffer;

	if(0){
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
				COS_SIMBOX(simpar).x.min, COS_SIMBOX(simpar).x.max, 
				COS_SIMBOX(simpar).y.min, COS_SIMBOX(simpar).y.max, 
				COS_SIMBOX(simpar).z.min, COS_SIMBOX(simpar).z.max);
	}
	dengrid->npix = (dengrid->jx-dengrid->ix+1)*(dengrid->jy-dengrid->iy+1)*(dengrid->jz-dengrid->iz+1);
	dengrid = (GridInfo*)realloc(dengrid,sizeof(GridInfo)+sizeof(DenType)*dengrid->npix);
	return dengrid;
}


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

	DEBUGPRINT("P%d has %d %d : %d %d : %d %d\n",MYID(simpar), dengrid->ix, dengrid->jx,
			dengrid->iy, dengrid->jy, dengrid->iz, dengrid->jz);

	for(j=0;j<ngrids;j++) den[j] = 0;
	if(DM_NP(simpar)>0) Get_TSC_Den(simpar, DM, den, xstart,ystart,zstart);
	DEBUGPRINT("P%d has dmnp = %ld sphnp= %ld\n",MYID(simpar), DM_NP(simpar), SPH_NP(simpar));
	if(SPH_NP(simpar)>0) Get_TSC_Den(simpar, SPH, den, xstart,ystart,zstart);
	if(STAR_NP(simpar)>0) Get_TSC_Den(simpar, STAR, den, xstart,ystart,zstart);
	if(AGN_NP(simpar)>0) Get_TSC_Den(simpar, AGN,  den, xstart,ystart,zstart);

}

FFTWGridInfo GridInfo2FFTWGridInfo4fftw(SimParameters *simpar, GridInfo *dengrid){
	FFTWGridInfo fftwgrid;
	fftwgrid.nx = NX(simpar);
	fftwgrid.ny = NY(simpar);
	fftwgrid.nz = NZ(simpar);
	fftwgrid.nspace = NSPACE(simpar);
	fftwgrid.fftw_info = FFTWINFO(simpar);
	fftwgrid.gridinfo = *dengrid;
	fftwgrid.den.rden = (DenType*)(dengrid + 1);

	return fftwgrid;
}
FFTWGridInfo ReallocDen4FFTW_Format(SimParameters *simpar, GridInfo *fftwgrid){
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

	FFTWINFO(simpar).com = MPI_COMM(simpar);
	FFTWGridInfo res = GridInfo2FFTWGridInfo4fftw(simpar,fftwgrid);
	return res;
}



FFTWGridInfo MakeGrids4FFTWDD(SimParameters *simpar, GridInfo *dengrid){
	GridInfo *fftwgrid = malloc(sizeof(GridInfo)+sizeof(fftwf_complex)*LOCAL_GRID_SIZE(simpar));
	fftwgrid->ix = fftwgrid->iy = 0;
	fftwgrid->jx = NX(simpar)-1;
	fftwgrid->jy = NY(simpar)-1;
	fftwgrid->iz = LOCAL_Z_START(simpar);
	fftwgrid->jz = LOCAL_Z_START(simpar) + LOCAL_NZ(simpar)-1;
	DenType *den = (DenType*)(fftwgrid + 1);
	ptrdiff_t i;
	for(i=0;i<2*LOCAL_GRID_SIZE(simpar);i++) den[i] = 0;
	gmigrate(dengrid, fftwgrid,DM_DDINFO(simpar), NDDINFO(simpar));
	FFTWGridInfo res = ReallocDen4FFTW_Format(simpar,fftwgrid);
	return res;
}
GridInfo *MakeDen4FDA(SimParameters *simpar, FFTWGridInfo *fftwgrid){
	ptrdiff_t i,j,k;
	ptrdiff_t nx,ny,nz,mx;
	nx = NX(simpar);
	ny = NY(simpar);
	nz = NZ(simpar);
	mx = 2*(nx/2+1);
	DEBUGPRINT("P%d here \n",MYID(simpar));
	GridInfo *ingrid, *outgrid, *fda4grid;
	ingrid = malloc(sizeof(GridInfo)+sizeof(fftwf_complex)*LOCAL_GRID_SIZE(simpar));
	DEBUGPRINT("P%d here \n",MYID(simpar));
	ingrid->ix = ingrid->iy = 0;
	ingrid->jx = nx-1; ingrid->jy = ny-1;
	ingrid->iz = LOCAL_Z_START(simpar);
	ingrid->jz = LOCAL_Z_START(simpar) + LOCAL_NZ(simpar)-1;

	DenType *den = fftwgrid->den.rden;
	DenType *den2 = (DenType*)(ingrid+1);
	DEBUGPRINT("P%d here \n",MYID(simpar));
	for(k=0;k<LOCAL_Z_START(simpar);k++) for(j=0;j<NY(simpar);j++){
		memmove(den2+nx*(j+ny*k), den + mx*(j+ny*k), sizeof(DenType)*nx);
	}
	DEBUGPRINT("P%d is freeing the fftw data \n",MYID(simpar));
	free((char*)FFTW_RDEN(simpar) - sizeof(GridInfo));
	DEBUGPRINT("P%d is freeing the fftw data \n",MYID(simpar));
	
	outgrid = malloc(sizeof(GridInfo));
	outgrid->ix = rint(COS_SIMBOX(simpar).x.min)-2;
	outgrid->iy = rint(COS_SIMBOX(simpar).y.min)-2;
	outgrid->iz = rint(COS_SIMBOX(simpar).z.min)-2;
	outgrid->jx = rint(COS_SIMBOX(simpar).x.max)+2;
	outgrid->jy = rint(COS_SIMBOX(simpar).y.max)+2;
	outgrid->jz = rint(COS_SIMBOX(simpar).z.max)+2;

	fda4grid = getoutgridnpix(outgrid, 0.);
	gmigrate(ingrid, fda4grid, DM_DDINFO(simpar), NDDINFO(simpar));

	return fda4grid;
}

void pmmain(SimParameters *simpar, EvolFact evolfact){
	DoDeInfo *ddinfo;
	float fact1, fact2;
	int nddinfo;

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

	GridInfo *dengrid, *fdagrid;
	ptrdiff_t ngrids;

	UpdateSimBox(simpar);
	dengrid = (GridInfo*)malloc(sizeof(GridInfo));


	dengrid = (GridInfo *)MakeDenGrid4TSC(simpar, dengrid, 1, 2);

	TSC(simpar, dengrid);

	FFTWGridInfo fftwgrid = MakeGrids4FFTWDD(simpar, dengrid);

	fftwf_mpi_execute_dft_r2c(FFTW_B_PLAN(simpar), fftwgrid.den.rden, (fftwf_complex*)(fftwgrid.den.rden));

	psolver(simpar, &fftwgrid);

	MPI_Barrier(MPI_COMM(simpar)); DEBUGPRINT("P%d is now at pmmain %ld %g\n",MYID(simpar), DM_NP(simpar), XofP(simpar, DM_BP(simpar)));

	fftwf_mpi_execute_dft_c2r(FFTW_F_PLAN(simpar), (fftwf_complex*)(fftwgrid.den.cden),fftwgrid.den.rden);

	MPI_Barrier(MPI_COMM(simpar)); DEBUGPRINT("P%d is now at pmmain %ld %g\n",MYID(simpar), DM_NP(simpar), XofP(simpar, DM_BP(simpar)));

	fdagrid = MakeDen4FDA(simpar, &fftwgrid);

	DEBUGPRINT("P%d is now at pmmain\n",MYID(simpar));

	fda4(simpar, fdagrid, fact1, fact2);
	DEBUGPRINT("P%d is now at pmmain\n",MYID(simpar));
}
