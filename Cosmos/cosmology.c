#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
//#include "mpirms.h"
#include "cosmology.h"
#include "cosmosic.h"
#include "indT.h"
#include "eunha.sub.h"
#include "kjhrw.h"
#include "fft.h"
#include "parallelIO.h"

#define MIN(a,b) ( (a)<(b)?(a):(b))
#define MAX(a,b) ( (a)>(b)?(a):(b))

void CosmologicalIC(SimParameters *simpar){
	if(PORDER(simpar) ==1) {
		DEBUGPRINT("Still not implemented and terminated : porder = %d\n", PORDER(simpar));
		mpi_fftw_finalize();
	}
	else if(PORDER(simpar) ==2){
		if(FLAG_IC(simpar) ==0) TwoLPTmain(simpar);
		else if(FLAG_IC(simpar) ==1) readGraficFile(simpar);
	}
	else {
		DEBUGPRINT("P%d has error in PORDER %d\n", MYID(simpar), PORDER(simpar));
		mpi_fftw_finalize();
	}
}
#define MakeAllDDFuncDDinfo(type)  \
EunhaParticleFuncs(type##particletype);\
DDRule3D(type##particletype, a, nmem, com);\
InSideOrNot(type##particletype, a, lbox, width, simbox, mflag, pfflag, padptl);\
determineEdgePtl(type##particletype, a, simbox, width);\
Dmuladd(type##particletype, a, b, fact, xyz, iflag);\
EunhaParticleFuncs(tree##type##particletype);\
DDRule3D(tree##type##particletype, a, nmem, com);\
InSideOrNot(tree##type##particletype, a, lbox, width, simbox, mflag, pfflag, padptl);\
determineEdgePtl(tree##type##particletype, a, simbox, width);\
Dmuladd(tree##type##particletype, a, b, fact, xyz, iflag);

MakeAllDDFuncDDinfo(dm);
MakeAllDDFuncDDinfo(sph);
MakeAllDDFuncDDinfo(voro);
MakeAllDDFuncDDinfo(star);
MakeAllDDFuncDDinfo(agn);
	


/*
EunhaParticleFuncs(dmparticletype);
EunhaParticleFuncs(sphparticletype);
EunhaParticleFuncs(starparticletype);
EunhaParticleFuncs(agnparticletype);

EunhaParticleFuncs(treedmparticletype);
EunhaParticleFuncs(treesphparticletype);
EunhaParticleFuncs(treestarparticletype);
EunhaParticleFuncs(treeagnparticletype);


DDRule3D(dmparticletype, a, nmem, com);
DDRule3D(sphparticletype, a, nmem, com);
DDRule3D(starparticletype, a, nmem, com);
DDRule3D(agnparticletype, a, nmem, com);

DDRule3D(treedmparticletype, a, nmem, com);
DDRule3D(treesphparticletype, a, nmem, com);
DDRule3D(treestarparticletype, a, nmem, com);
DDRule3D(treeagnparticletype, a, nmem, com);


InSideOrNot(dmparticletype, a, lbox, width, simbox, mflag, pflag);
InSideOrNot(sphparticletype, a, lbox, width, simbox, mflag, pflag);
InSideOrNot(starparticletype, a, lbox, width, simbox, mflag, pflag);
InSideOrNot(agnparticletype, a, lbox, width, simbox, mflag, pflag);

InSideOrNot(treedmparticletype, a, lbox, width, simbox, mflag, pflag);
InSideOrNot(treesphparticletype, a, lbox, width, simbox, mflag, pflag);
InSideOrNot(treestarparticletype, a, lbox, width, simbox, mflag, pflag);
InSideOrNot(treeagnparticletype, a, lbox, width, simbox, mflag, pflag);


determineEdgePtl(dmparticletype, a, simbox, width);
determineEdgePtl(starparticletype, a, simbox, width);
determineEdgePtl(sphparticletype, a, simbox, width);
determineEdgePtl(agnparticletype, a, simbox, width);

determineEdgePtl(treedmparticletype, a, simbox, width);
determineEdgePtl(treestarparticletype, a, simbox, width);
determineEdgePtl(treesphparticletype, a, simbox, width);
determineEdgePtl(treeagnparticletype, a, simbox, width);


Dmuladd(dmparticletype, a, b, fact, xyz, iflag);
Dmuladd(sphparticletype, a, b, fact, xyz, iflag);
Dmuladd(starparticletype, a, b, fact, xyz, iflag);
Dmuladd(agnparticletype, a, b, fact, xyz, iflag);

Dmuladd(treedmparticletype, a, b, fact, xyz, iflag);
Dmuladd(treesphparticletype, a, b, fact, xyz, iflag);
Dmuladd(treestarparticletype, a, b, fact, xyz, iflag);
Dmuladd(treeagnparticletype, a, b, fact, xyz, iflag);
*/


void BuildCosmosDDInfo(SimParameters *simpar){
	int nprime;
	PrimeNumber prime[100];
	nprime = getprimenumber(NID(simpar), prime);
	
	MakeDoDeInfo3D(NID(simpar), prime, dmparticletype, x,y,z, DM_DDINFO(simpar), NDDINFO(simpar));
	MakeDoDeInfo3D(NID(simpar), prime, sphparticletype, x,y,z, SPH_DDINFO(simpar), NDDINFO(simpar));
	MakeDoDeInfo3D(NID(simpar), prime, voroparticletype, x,y,z, VORO_DDINFO(simpar), NDDINFO(simpar));
	MakeDoDeInfo3D(NID(simpar), prime, starparticletype, x,y,z, STAR_DDINFO(simpar), NDDINFO(simpar));
	MakeDoDeInfo3D(NID(simpar), prime, agnparticletype, x,y,z, AGN_DDINFO(simpar), NDDINFO(simpar));


	MakeDoDeInfo3D(NID(simpar), prime, treedmparticletype, x,y,z, TDM_DDINFO(simpar), NDDINFO(simpar));
	MakeDoDeInfo3D(NID(simpar), prime, treesphparticletype, x,y,z, TSPH_DDINFO(simpar), NDDINFO(simpar));
	MakeDoDeInfo3D(NID(simpar), prime, treevoroparticletype, x,y,z, TVORO_DDINFO(simpar), NDDINFO(simpar));
	MakeDoDeInfo3D(NID(simpar), prime, treestarparticletype, x,y,z, TSTAR_DDINFO(simpar), NDDINFO(simpar));
	MakeDoDeInfo3D(NID(simpar), prime, treeagnparticletype, x,y,z, TAGN_DDINFO(simpar), NDDINFO(simpar));

	MakeDoDeFunc3D(dmparticletype, DM_DDFUNC(simpar), xcompare,ycompare, zcompare, xpinner, ypinner, zpinner, 
			muladd, DDRule3D, InSideBox, determineEdgePtl, shift2pos, DM_DDINFO(simpar));
	MakeDoDeFunc3D(sphparticletype, SPH_DDFUNC(simpar), xcompare,ycompare, zcompare, xpinner, ypinner, zpinner, 
			muladd, DDRule3D, InSideBox, determineEdgePtl, shift2pos, SPH_DDINFO(simpar));
	MakeDoDeFunc3D(voroparticletype, VORO_DDFUNC(simpar), xcompare,ycompare, zcompare, xpinner, ypinner, zpinner, 
			muladd, DDRule3D, InSideBox, determineEdgePtl, shift2pos, VORO_DDINFO(simpar));
	MakeDoDeFunc3D(starparticletype, STAR_DDFUNC(simpar), xcompare,ycompare, zcompare, xpinner, ypinner, zpinner, 
			muladd, DDRule3D, InSideBox, determineEdgePtl, shift2pos, STAR_DDINFO(simpar));
	MakeDoDeFunc3D(agnparticletype, AGN_DDFUNC(simpar), xcompare,ycompare, zcompare, xpinner, ypinner, zpinner, 
			muladd, DDRule3D, InSideBox, determineEdgePtl, shift2pos, AGN_DDINFO(simpar));

	MakeDoDeFunc3D(treedmparticletype, TDM_DDFUNC(simpar), xcompare,ycompare, zcompare, xpinner, ypinner, zpinner, 
			muladd, DDRule3D, InSideBox, determineEdgePtl, shift2pos, TDM_DDINFO(simpar));
	MakeDoDeFunc3D(treesphparticletype, TSPH_DDFUNC(simpar), xcompare,ycompare, zcompare, xpinner, ypinner, zpinner, 
			muladd, DDRule3D, InSideBox, determineEdgePtl, shift2pos, TSPH_DDINFO(simpar));
	MakeDoDeFunc3D(treevoroparticletype, TVORO_DDFUNC(simpar), xcompare,ycompare, zcompare, xpinner, ypinner, zpinner, 
			muladd, DDRule3D, InSideBox, determineEdgePtl, shift2pos, TVORO_DDINFO(simpar));
	MakeDoDeFunc3D(treestarparticletype, TSTAR_DDFUNC(simpar), xcompare,ycompare, zcompare, xpinner, ypinner, zpinner, 
			muladd, DDRule3D, InSideBox, determineEdgePtl, shift2pos, TSTAR_DDINFO(simpar));
	MakeDoDeFunc3D(treeagnparticletype, TAGN_DDFUNC(simpar), xcompare,ycompare, zcompare, xpinner, ypinner, zpinner, 
			muladd, DDRule3D, InSideBox, determineEdgePtl, shift2pos, TAGN_DDINFO(simpar));
}
/*
void NewDDs(SimParameters *simpar){
	mpirms( (void**)(&DM_BP(simpar)),&DM_NP(simpar), sizeof(dmparticletype), &DM_DDFUNC(simpar), DM_DDINFO(simpar),COM(simpar),
			&GRIDINFO(simpar));
	ExtractLocalDomainVolume(DM_DDINFO(simpar), NDDINFO(simpar),COS_SIMBOX(simpar));
	if(SPH_TNP(simpar)>0){
		pmigrate((void**)(&SPH_BP(simpar)), &SPH_NP(simpar), SPH_DDINFO(simpar), &GRIDINFO(simpar));
		MPI_Allreduce(&SPH_NP(simpar), &SPH_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
	}
	if(STAR_TNP(simpar)>0){ 
		pmigrate((void**)(&STAR_BP(simpar)), &STAR_NP(simpar), STAR_DDINFO(simpar), &GRIDINFO(simpar));
		MPI_Allreduce(&STAR_NP(simpar), &STAR_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
	}
	if(AGN_TNP(simpar)>0){ 
		pmigrate((void**)(&AGN_BP(simpar)), &AGN_NP(simpar), AGN_DDINFO(simpar), &GRIDINFO(simpar));
		MPI_Allreduce(&AGN_NP(simpar), &AGN_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
	}
}
*/

void AllParticleMigrate(SimParameters *simpar){
	MPI_Allreduce(&DM_NP(simpar), &DM_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
	MPI_Allreduce(&SPH_NP(simpar), &SPH_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
	MPI_Allreduce(&VORO_NP(simpar), &VORO_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
	MPI_Allreduce(&STAR_NP(simpar), &STAR_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
	MPI_Allreduce(&AGN_NP(simpar), &AGN_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
	if(DM_TNP(simpar)>0) {
		pmigrate((void**)(&DM_BP(simpar)), &DM_NP(simpar), DM_DDINFO(simpar), &GRIDINFO(simpar));
	}
	if(SPH_TNP(simpar)>0){
		pmigrate((void**)(&SPH_BP(simpar)), &SPH_NP(simpar), SPH_DDINFO(simpar), &GRIDINFO(simpar));
	}
	if(VORO_TNP(simpar)>0){
		pmigrate((void**)(&VORO_BP(simpar)), &VORO_NP(simpar), VORO_DDINFO(simpar), &GRIDINFO(simpar));
	}
	if(STAR_TNP(simpar)>0){ 
		pmigrate((void**)(&STAR_BP(simpar)), &STAR_NP(simpar), STAR_DDINFO(simpar), &GRIDINFO(simpar));
	}
	if(AGN_TNP(simpar)>0){ 
		pmigrate((void**)(&AGN_BP(simpar)), &AGN_NP(simpar), AGN_DDINFO(simpar), &GRIDINFO(simpar));
	}
	NPSUM(simpar) = DM_NP(simpar) + SPH_NP(simpar) + STAR_NP(simpar) + AGN_NP(simpar) + VORO_NP(simpar);
}

void TreeAllParticleMigrate(SimParameters *simpar){
	MPI_Allreduce(&DM_NP(simpar), &DM_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
	MPI_Allreduce(&SPH_NP(simpar), &SPH_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
	MPI_Allreduce(&VORO_NP(simpar), &VORO_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
	MPI_Allreduce(&STAR_NP(simpar), &STAR_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
	MPI_Allreduce(&AGN_NP(simpar), &AGN_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
#ifdef DEBUG
	DEBUGPRINT("P%d has %ld : %ld : %ld : %ld : %ld\n", MYID(simpar), DM_NP(simpar), SPH_NP(simpar), 
			VORO_NP(simpar), AGN_NP(simpar), STAR_NP(simpar));
#endif
//	HAMB(simpar);
	if(DM_TNP(simpar)>0) {
		pmigrate((void**)(&DM_TBP(simpar)), &DM_NP(simpar), TDM_DDINFO(simpar), &GRIDINFO(simpar));
	}
	if(SPH_TNP(simpar)>0){
		pmigrate((void**)(&SPH_TBP(simpar)), &SPH_NP(simpar), TSPH_DDINFO(simpar), &GRIDINFO(simpar));
	}
	if(VORO_TNP(simpar)>0){
		pmigrate((void**)(&VORO_TBP(simpar)), &VORO_NP(simpar), TVORO_DDINFO(simpar), &GRIDINFO(simpar));
	}

	if(STAR_TNP(simpar)>0){ 
		pmigrate((void**)(&STAR_TBP(simpar)), &STAR_NP(simpar), TSTAR_DDINFO(simpar), &GRIDINFO(simpar));
	}
	if(AGN_TNP(simpar)>0){ 
		pmigrate((void**)(&AGN_TBP(simpar)), &AGN_NP(simpar), TAGN_DDINFO(simpar), &GRIDINFO(simpar));
	}
	NPSUM(simpar) = DM_NP(simpar) + SPH_NP(simpar) + STAR_NP(simpar) + AGN_NP(simpar) + VORO_NP(simpar);
#ifdef DEBUG
	DEBUGPRINT("P%d has %ld : %ld : %ld : %ld : %ld\n", MYID(simpar), DM_NP(simpar), SPH_NP(simpar), 
			VORO_NP(simpar),
			AGN_NP(simpar), STAR_NP(simpar));
#endif
}


void AllParticlePadding(SimParameters *simpar, float width){
	ppadding(DM_BP(simpar), DM_NP(simpar), (void**)(&DM_BPP(simpar)), &DM_NPAD(simpar), 
			DM_DDINFO(simpar), NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar), 3);
	ppadding(SPH_BP(simpar), SPH_NP(simpar), (void**)(&SPH_BPP(simpar)), &SPH_NPAD(simpar), 
			SPH_DDINFO(simpar), NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar), 3); 
	ppadding(VORO_BP(simpar), VORO_NP(simpar), (void**)(&VORO_BPP(simpar)), &VORO_NPAD(simpar), 
			VORO_DDINFO(simpar), NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar), 3); 
	ppadding(STAR_BP(simpar), STAR_NP(simpar), (void**)(&STAR_BPP(simpar)), &STAR_NPAD(simpar), 
			STAR_DDINFO(simpar), NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar), 3);
	ppadding(AGN_BP(simpar), AGN_NP(simpar), (void**)(&AGN_BPP(simpar)), &AGN_NPAD(simpar), 
			AGN_DDINFO(simpar), NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar), 3);
}

void TreeAllParticlePadding(SimParameters *simpar, float width){
	ppadding(DM_TBP(simpar), DM_NP(simpar), (void**)(&DM_TBPP(simpar)), &DM_NPAD(simpar), 
			TDM_DDINFO(simpar), NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar), 3);
	ppadding(SPH_TBP(simpar), SPH_NP(simpar), (void**)(&SPH_TBPP(simpar)), &SPH_NPAD(simpar), 
			TSPH_DDINFO(simpar), NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar), 3); 
	ppadding(VORO_TBP(simpar), VORO_NP(simpar), (void**)(&VORO_TBPP(simpar)), &VORO_NPAD(simpar), 
			TVORO_DDINFO(simpar), NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar), 3); 
	ppadding(STAR_TBP(simpar), STAR_NP(simpar), (void**)(&STAR_TBPP(simpar)), &STAR_NPAD(simpar), 
			TSTAR_DDINFO(simpar), NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar), 3);
	ppadding(AGN_TBP(simpar), AGN_NP(simpar), (void**)(&AGN_TBPP(simpar)), &AGN_NPAD(simpar), 
			TAGN_DDINFO(simpar), NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar), 3);
}


void Cosmos_IndT(SimParameters *simpar, int);
void GOTPM_Cosmos_GlobalT(SimParameters *simpar, int);


void StartCosmosRkSDD(SimParameters *);
void SetEntireSimBox(SimParameters *);
/*
void UpdateDdinfoBox(SimParameters *);
*/

void InitializeReadStart(SimParameters *simpar){
	DM_MASS(simpar,0) = 1;
	/*
	UpdateDdinfoBox(simpar);
	*/
}

int RunCosmos(SimParameters *simpar, int icont)
{
#ifdef GOTPM
	if(MYID(simpar)==0){
		printf("GOTPM Simulation is Starting\n");
	}
#endif
	set_mpi_fftw_3d_info(simpar);
	set_mpi_fftw_3d_plans(simpar);

	SetEnvGridInfo(simpar);

	StartCosmosRkSDD(simpar); /* This makes all the ddinfo (pivot) */

	if(icont == 0) {
		CosmologicalIC(simpar);
		SetEntireSimBox(simpar);
	}
	else {
		InitializeReadStart(simpar);
		jread(simpar);



	}


#ifdef DEBUG
	DEBUGPRINT("P%d: now %g %g %g %g %g %g\n", MYID(simpar),
			COS_SIMBOX(simpar).x.min,
			COS_SIMBOX(simpar).y.min,
			COS_SIMBOX(simpar).z.min,
			COS_SIMBOX(simpar).x.max,
			COS_SIMBOX(simpar).y.max,
			COS_SIMBOX(simpar).z.max
			);
#endif

	if(IndT_FLAG(simpar) == 'Y') Cosmos_IndT(simpar, icont);
	else GOTPM_Cosmos_GlobalT(simpar,icont);

	return 1;
}
void SetEnvGridInfo(SimParameters *simpar){
	GridInfo *grid = &GRIDINFO(simpar);
	PosNX(grid) = NX(simpar);
	PosNY(grid) = NY(simpar);
	PosNZ(grid) = NZ(simpar);
	PosNXNY(grid) = NX(simpar) * NY(simpar);
}



/*
void RunCosSim( SimParameters *simpar, int icont){

	if(IndT_FLAG(simpar)) Cosmos_IndT(simpar, icont);
	else Cosmos_GlobalT(simpar, icont);

	for(STEPCOUNT(simpar)= STEPNUM(simpar); STEPCOUNT(simpar) <= NSTEP(simpar);STEPCOUNT(simpar)++){
	}
}
*/
