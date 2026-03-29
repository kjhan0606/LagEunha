#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
/*
#include "mpirms.h"
*/
#include "cosmosic.h"
#include "flrw.h"
#include "parallelIO.h"



#define MAX(a,b) ( (a)>(b) ? (a): (b) )
#define MIN(a,b) ( (a)>(b)?(b):(a))

/*
static SimParameters lsimpar;

void getpowerparams_(int *nspace, float *boxsize, float *amax,
		float *npow, float *hubble, float *bias8,
		long *iseed, float *omepb, int *powreadflag, char *powfilename,
		char *asciipkfile){
	SimParameters *simpar = &lsimpar;
	*nspace = NSPACE(simpar);
	*boxsize = COSBOXSIZE(simpar);
	*amax = AMAX(simpar);
	*npow = NPOW(simpar);
	*hubble = HUBBLE(simpar);
	*bias8 = BIAS8(simpar);
	*iseed = ISEED(simpar);
	*omepb = OMEPB(simpar);
	*powreadflag = POWREADFLAG(simpar);
	sprintf(powfilename, "%s", POWFILENAME(simpar));
	sprintf(asciipkfile, "%s", INPAPKFILENAME(simpar));
}
*/
#define DDBBUGG do{\
		DEBUGPRINT("P%d has amplitudes %g %g %g %g for type= %d\n", \
				MYID(simpar),COS_DAMP1(simpar,itype),COS_DAMP2(simpar,itype), \
				COS_VAMP1(simpar,itype), COS_VAMP2(simpar,itype),itype);\
		ptrdiff_t mx = 2*(NX(simpar)/2+1); \
		if(0){\
			DEBUGPRINT("P%d is now dumping IDen.dat for investigating\n", MYID(simpar));\
			ptrdiff_t mx = 2*(NX(simpar)/2+1);\
			double rincube = 1.L/NX(simpar)/NY(simpar)/NZ(simpar);\
			fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),(fftwf_complex*)d2, d2);\
			for(k=0;k<local_nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++){\
   	         d2[i+mx*(j+ny*k)] *= rincube;\
   	     }\
			FILE *wp;\
   	     if(MYID(simpar)==0){\
   	         wp = fopen("IDen.dat","w");\
   	         int nnx,nny,nnz;\
   	         nnx = NX(simpar);\
   	         nny = NY(simpar);\
   	         nnz = NZ(simpar);\
   	         fwrite(&nnx,sizeof(int),1,wp);\
   	         fwrite(&nny,sizeof(int),1,wp);\
   	         fwrite(&nnz,sizeof(int),1,wp);\
   	         fclose(wp);\
   	     }\
   	     for(i=0;i<NID(simpar);i++){\
   	         if(MYID(simpar)==i){\
   	             wp = fopen("IDen.dat","a");\
   	             for(k=0;k<local_nz;k++) for(j=0;j<ny;j++){\
   	                 fwrite(d2+mx*(j+ny*k),sizeof(float), nx,wp);\
   	             }\
   	             fclose(wp);\
   	         }\
   	         MPI_Barrier(MPI_COMM(simpar));\
   	     }\
   	     MPI_Finalize();\
   	     exit(99);\
		}\
		double valmax = 0;\
		for(k=0;k<LOCAL_NZ(simpar);k++) for(j=0;j<NY(simpar);j++) for(i=0;i<NX(simpar);i++){\
			valmax = MAX(valmax,d1[i+mx*(j+ny*k)]);\
		}\
		DEBUGPRINT("P%d has max d1 %g\n",MYID(simpar), valmax);\
		valmax = 0;\
		for(k=0;k<LOCAL_NZ(simpar);k++) for(j=0;j<NY(simpar);j++) for(i=0;i<NX(simpar);i++){\
			valmax = MAX(valmax,d2[i+mx*(j+ny*k)]);\
		}\
		DEBUGPRINT("P%d has max d2 %g\n",MYID(simpar), valmax);\
} while(0)

void expectDMandSphNumParticle(SimParameters *simpar){
	DM_NP(simpar) = NX(simpar)*NY(simpar)*LOCAL_NZ(simpar);
	if(GAS_TYPE(simpar) == 'S') {
		SPH_NP(simpar) = NX(simpar)*NY(simpar)*LOCAL_NZ(simpar);
		SPH_INIT_MASS(simpar) = 1 - OMEPB(simpar)/OMEP(simpar);
	}
	else if(GAS_TYPE(simpar) == 'V') {
		VORO_NP(simpar) = NX(simpar)*NY(simpar)*LOCAL_NZ(simpar);
		VORO_INIT_MASS(simpar) = 1 - OMEPB(simpar)/OMEP(simpar);
	}
	else {
		SPH_NP(simpar) = 0;
		VORO_NP(simpar) = 0;
		DM_MASS(simpar,0) = 1;
	}
	MPI_Allreduce(&DM_NP(simpar), &DM_TNP(simpar), 1, MPI_LONG, MPI_SUM, MPI_COMM(simpar));
#ifdef DEBUG
	if(MYID(simpar)==0) DEBUGPRINT("Total DM np = %ld\n",DM_TNP(simpar));
#endif
	MPI_Allreduce(&SPH_NP(simpar), &SPH_TNP(simpar), 1, MPI_LONG, MPI_SUM, MPI_COMM(simpar));
	MPI_Allreduce(&VORO_NP(simpar), &VORO_TNP(simpar), 1, MPI_LONG, MPI_SUM, MPI_COMM(simpar));
#ifdef DEBUG
	if(MYID(simpar)==0) DEBUGPRINT("Total SPH np = %ld & voro np = %ld\n",
			SPH_TNP(simpar), VORO_TNP(simpar));
#endif
}

/* This is moved to ../CosmosRW/readGraficFile.c.
GridInfo *makeFftDomainVolume(SimParameters *simpar, size_t densize){
	GridInfo *a = (GridInfo*)malloc(sizeof(GridInfo)+sizeof(DenType)*densize);
	a->nx = NX(simpar); a->ny = NY(simpar); a->nz = NZ(simpar);
	a->ix = a->iy = 0;
	a->jx = NX(simpar)-1;
	a->jy = NY(simpar)-1;
	a->iz = LOCAL_Z_START(simpar);
	a->jz = a->iz + LOCAL_NZ(simpar)-1;
	a->npix = (ptrdiff_t)(a->jx-a->ix+1) * (ptrdiff_t)(a->jy-a->iy+1) * (ptrdiff_t)(a->jz-a->iz+1);
	return a;
}
GridInfo *makeRmsDomainVolume(SimParameters *simpar, enum mtype itype, int nbuff){
	GridInfo *a = (GridInfo*)malloc(sizeof(GridInfo));
	a->nx = NX(simpar); a->ny = NY(simpar); a->nz = NZ(simpar);

	if(itype <=cdm) {
		a->ix = rint(SIM_LXMIN(simpar,dm)) - nbuff;
		a->iy = rint(SIM_LYMIN(simpar,dm)) - nbuff;
		a->iz = rint(SIM_LZMIN(simpar,dm)) - nbuff;
		a->jx = rint(SIM_LXMAX(simpar,dm)) - 1 + nbuff;
		a->jy = rint(SIM_LYMAX(simpar,dm)) - 1 + nbuff;
		a->jz = rint(SIM_LZMAX(simpar,dm)) - 1 + nbuff;
	}
	else{
		a->ix = rint(SIM_LXMIN(simpar,sph)) - nbuff;
		a->iy = rint(SIM_LYMIN(simpar,sph)) - nbuff;
		a->iz = rint(SIM_LZMIN(simpar,sph)) - nbuff;
		a->jx = rint(SIM_LXMAX(simpar,sph)) - 1 + nbuff;
		a->jy = rint(SIM_LYMAX(simpar,sph)) - 1 + nbuff;
		a->jz = rint(SIM_LZMAX(simpar,sph)) - 1 + nbuff;
	}

	a->npix = (ptrdiff_t)(a->jx-a->ix+1) *(ptrdiff_t)(a->jy-a->iy+1) *(ptrdiff_t)(a->jz-a->iz+1);
	a = realloc(a, sizeof(GridInfo)+sizeof(DenType)*a->npix);
	return a;
}
*/

void TwoLPTmain(SimParameters *simpar){
	size_t i,j,k;
	enum mtype itype;


	expectDMandSphNumParticle(simpar);

	for(itype=cdm;itype<star;itype++)
	{
		if(itype ==cdm && DM_TNP(simpar) ==0) continue;
		else if(itype ==sph && SPH_TNP(simpar) ==0) continue;
		else if(itype ==voro && VORO_TNP(simpar) ==0) continue;

		size_t densize = 2*LOCAL_GRID_SIZE(simpar);
		densize = MAX(10, densize);

		DenType *d1 = (DenType*)malloc(2*sizeof(DenType)*LOCAL_GRID_SIZE(simpar));
		DenType *d2 = (DenType*)malloc(2*sizeof(DenType)*LOCAL_GRID_SIZE(simpar));
		DenType *denpad = (DenType*)malloc(sizeof(DenType)*2*NY(simpar)*NZ(simpar));

		int nx = NX(simpar);
		int ny = NY(simpar);
		int local_nz = LOCAL_NZ(simpar);
		int local_z_start = LOCAL_Z_START(simpar);

		DenType *fx = (DenType*)malloc(2*sizeof(DenType)*LOCAL_GRID_SIZE(simpar));
		DenType *fy = (DenType*)malloc(2*sizeof(DenType)*LOCAL_GRID_SIZE(simpar));
		DenType *fz = (DenType*)malloc(2*sizeof(DenType)*LOCAL_GRID_SIZE(simpar));
//void PMSeedForce(SimParameters *, DenType *, DenType *, DenType *, DenType *, DenType *, DenType *, enum mtype);
		PMSeedForce(simpar, d1,fx,fy,fz, d2, denpad, itype);
		free(fx);free(fy);free(fz);
		free(denpad);
#ifdef DEBUG
		DDBBUGG;
#endif
		float damp1, vamp1, damp2, vamp2;
		damp1 = COS_DAMP1(simpar, itype);
		damp2 = COS_DAMP2(simpar, itype);
		vamp1 = COS_VAMP1(simpar, itype);
		vamp2 = COS_VAMP2(simpar, itype);


		GridInfo *ffx, *ffy, *ffz;
		if(itype==cdm){
			TwoLPT_GetVR(simpar,ffx, d1, d2, vamp1, vamp2, VX,dm);
			TwoLPT_GetVR(simpar,ffy, d1, d2, vamp1, vamp2, VY,dm);
			TwoLPT_GetVR(simpar,ffz, d1, d2, vamp1, vamp2, VZ,dm);
		}
		else if(itype ==sph){
			TwoLPT_GetVR(simpar,ffx, d1, d2, vamp1, vamp2, VX,sph);
			TwoLPT_GetVR(simpar,ffy, d1, d2, vamp1, vamp2, VY,sph);
			TwoLPT_GetVR(simpar,ffz, d1, d2, vamp1, vamp2, VZ,sph);
		}
		Set_Particle_2LPT_Vel(simpar, ffx,ffy,ffz,itype);
		float amp1, amp2;
		amp1 = AI(simpar)*NX(simpar)*damp1;
		amp2 = AI(simpar)*NX(simpar)*damp2;
#ifdef DEBUG
		DEBUGPRINT("P%d has amplitudes : %g %g %g %g\n",MYID(simpar), amp1,amp2,vamp1,vamp2);
#endif
		if(itype ==cdm){
			TwoLPT_GetVR(simpar,ffx, d1, d2, amp1, amp2, X,dm);
			TwoLPT_GetVR(simpar,ffy, d1, d2, amp1, amp2, Y,dm);
			TwoLPT_GetVR(simpar,ffz, d1, d2, amp1, amp2, Z,dm);
		}
		else if(itype ==sph){
			TwoLPT_GetVR(simpar,ffx, d1, d2, amp1, amp2, X,sph);
			TwoLPT_GetVR(simpar,ffy, d1, d2, amp1, amp2, Y,sph);
			TwoLPT_GetVR(simpar,ffz, d1, d2, amp1, amp2, Z,sph);
		}
#ifdef DEBUG
		if(0){
			void dumpphi(SimParameters *, float *);
			dumpphi(simpar, d1);
		}
#endif

		Set_Particle_2LPT_Pos(simpar, ffx,ffy,ffz,itype);
		free(ffz);free(ffy);free(ffx);
		free(d2);free(d1);

		if(itype ==sph){
			for(i=0;i<SPH_NP(simpar);i++){
				SPH_MASS(simpar,i) = OMEPB(simpar)/OMEP(simpar);
			}
		}
#ifdef DEBUG
		{
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
		}
#endif
		GRIDINFO(simpar).nx = NX(simpar);
		GRIDINFO(simpar).ny = NY(simpar);
		GRIDINFO(simpar).nz = NZ(simpar);
		GRIDINFO(simpar).nxny = NX(simpar)*NY(simpar);

		if(itype==cdm) pmigrate((void**)(&DM_BP(simpar)), &DM_NP(simpar), DM_DDINFO(simpar), &GRIDINFO(simpar));
		else if(itype==sph) pmigrate((void**)(&SPH_BP(simpar)), &SPH_NP(simpar), SPH_DDINFO(simpar), &GRIDINFO(simpar));

#ifdef DEBUG
		{
			MPI_Barrier(MPI_COMM(simpar));
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
					SIM_LXMIN(simpar,dm), SIM_LXMAX(simpar, dm),
					SIM_LYMIN(simpar,dm), SIM_LYMAX(simpar, dm),
					SIM_LZMIN(simpar,dm), SIM_LZMAX(simpar, dm));
			MPI_Barrier(MPI_COMM(simpar));
		}
#endif
	}
#ifdef DEBUG
	fflush(stdout); fflush(stderr);
	MPI_Barrier(MPI_COMM(simpar));
	if(MYID(simpar)==0) DEBUGPRINT("Exiting 2LPT main part %ld & %ld\n", DM_NP(simpar), SPH_NP(simpar));
#endif
}

void dumpphi(SimParameters *simpar, float *den){ 
	ptrdiff_t i,j,k;
	ptrdiff_t local_nz = LOCAL_NZ(simpar);
	ptrdiff_t nx,ny,nz;
	nx = NX(simpar);
	ny = NY(simpar);
	nz = NZ(simpar);
	DEBUGPRINT("P%d is now dumping IDen.dat for investigating\n", MYID(simpar)); 
	ptrdiff_t mx = 2*(NX(simpar)/2+1); 
	double rincube = 1.L/NX(simpar)/NY(simpar)/NZ(simpar); 
	fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),(fftwf_complex*)den, den); 
	for(k=0;k<local_nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++){ 
		den[i+mx*(j+ny*k)] *= rincube; 
	} 
	FILE *wp; 
	if(MYID(simpar)==0){ 
		wp = fopen("IDen.dat","w"); 
		int nnx,nny,nnz; 
		nnx = NX(simpar); 
		nny = NY(simpar); 
		nnz = NZ(simpar); 
		fwrite(&nnx,sizeof(int),1,wp); 
		fwrite(&nny,sizeof(int),1,wp); 
		fwrite(&nnz,sizeof(int),1,wp); 
		fclose(wp); 
	} 
	for(i=0;i<NID(simpar);i++){ 
		if(MYID(simpar)==i){ 
			wp = fopen("IDen.dat","a"); 
			for(k=0;k<local_nz;k++) for(j=0;j<ny;j++){ 
				fwrite(den+mx*(j+ny*k),sizeof(float), nx,wp); 
			} 
			fclose(wp); 
		} 
		MPI_Barrier(MPI_COMM(simpar)); 
	} 
	MPI_Finalize(); 
	exit(99);
}
void dumppart(SimParameters *simpar){
	ptrdiff_t i,j,k;
	FILE *wp;

	for(i=0;i<NID(simpar);i++){
		if(MYID(simpar)==i){
			if(MYID(simpar) ==0) wp = fopen("ParticleCheck.dat","w");
			else wp = fopen("ParticleCheck.dat","a");
			for(j=0;j<DM_NP(simpar);j++){
				if(ZofP(simpar, DM_BP(simpar)+j) < 2){
					fprintf(wp,"%g %g %g %g %g %g :: %ld\n",
							XofP(simpar,DM_BP(simpar)+j),
							YofP(simpar,DM_BP(simpar)+j),
							ZofP(simpar,DM_BP(simpar)+j),
							(DM_BP(simpar)+j)->vx,
							(DM_BP(simpar)+j)->vy,
							(DM_BP(simpar)+j)->vz,
							PINDX(DM_BP(simpar)+j));
				}
			}
			fclose(wp);
		}
		MPI_Barrier(MPI_COMM(simpar));
	}
	MPI_Finalize();
	exit(9);
}
void dumpden(SimParameters *simpar, float *den, ptrdiff_t mx){ 
	ptrdiff_t i,j,k;
	ptrdiff_t local_nz = LOCAL_NZ(simpar);
	ptrdiff_t nx,ny,nz;
	nx = NX(simpar);
	ny = NY(simpar);
	nz = NZ(simpar);
	DEBUGPRINT("P%d is now dumping IDen.dat for investigating\n", MYID(simpar)); 
	FILE *wp; 
	if(MYID(simpar)==0){ 
		wp = fopen("IDen.dat","w"); 
		int nnx,nny,nnz; 
		nnx = NX(simpar); 
		nny = NY(simpar); 
		nnz = NZ(simpar); 
		fwrite(&nnx,sizeof(int),1,wp); 
		fwrite(&nny,sizeof(int),1,wp); 
		fwrite(&nnz,sizeof(int),1,wp); 
		fclose(wp); 
	} 
	for(i=0;i<NID(simpar);i++){ 
		if(MYID(simpar)==i){ 
			wp = fopen("IDen.dat","a"); 
			for(k=0;k<local_nz;k++) for(j=0;j<ny;j++){ 
				fwrite(den+mx*(j+ny*k),sizeof(float), nx,wp); 
			} 
			fclose(wp); 
		} 
		MPI_Barrier(MPI_COMM(simpar)); 
	} 
	MPI_Finalize(); 
	exit(99);
}
