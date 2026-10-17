#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<string.h>
#include<mpi.h>
#include "eunha.h"
//#include "mpirms.h"
#include "cosmology.h"
#include "globalT.h"
#include "indT.h"
#include "CosmosEvolFactor.h"
#include "domaindecomp.h"

#include "cosmosic.h"
#include "timerutil.h"
#include "flags.h"
#include "PreFoF.h"
#include "kjhrw.h"

/*
*/

/*
typedef struct Half_Step{
	int first, second;
}Half_Step;
*/

#define MIN(a,b) ( (a)<(b)?(a):(b))
#define MAX(a,b) ( (a)>(b)?(a):(b))



#include "aux.h"


void TreeAndPreFoFmain(SimParameters *simpar){
	float cputime0[5], cputime1[5];
	TIMER_START(4);
	if(STEPCOUNT(simpar) <10) THETA(simpar) = 0.3;
	else if(STEPCOUNT(simpar) <20) THETA(simpar) = 0.35;
	else if(STEPCOUNT(simpar) <60) THETA(simpar) = 0.4;
	else THETA(simpar) = 0.45;

#ifdef DEBUG
	{
		PosType xxmin,xxmax,yymin,yymax,zzmin,zzmax; 
		xxmin = yymin = zzmin = 1.e20; 
		xxmax = yymax = zzmax = -1.e20; 
		ptrdiff_t i; 
		for(i=0;i<DM_NP(simpar);i++){ 
			xxmin = MIN(xxmin, XofP(simpar, DM_BP(simpar) + i)); 
			yymin = MIN(yymin, YofP(simpar, DM_BP(simpar) + i)); 
			zzmin = MIN(zzmin, ZofP(simpar, DM_BP(simpar) + i)); 
			xxmax = MAX(xxmax, XofP(simpar, DM_BP(simpar) + i)); 
			yymax = MAX(yymax, YofP(simpar, DM_BP(simpar) + i)); 
			zzmax = MAX(zzmax, ZofP(simpar, DM_BP(simpar) + i)); 
		} 
		MPI_Barrier(MPI_COMM(simpar)); 
		DEBUGPRINT("P%d %g %g : %g %g : %g %g\n",MYID(simpar),xxmin,xxmax,yymin,yymax,zzmin,zzmax);

	}
#endif
	PM2TreeConversion(simpar);
#ifdef DEBUG
	{
		PosType xxmin,xxmax,yymin,yymax,zzmin,zzmax; 
		xxmin = yymin = zzmin = 1.e20; 
		xxmax = yymax = zzmax = -1.e20; 
		ptrdiff_t i; 
		for(i=0;i<DM_NP(simpar);i++){ 
			xxmin = MIN(xxmin, XofP(simpar, DM_TBP(simpar) + i)); 
			yymin = MIN(yymin, YofP(simpar, DM_TBP(simpar) + i)); 
			zzmin = MIN(zzmin, ZofP(simpar, DM_TBP(simpar) + i)); 
			xxmax = MAX(xxmax, XofP(simpar, DM_TBP(simpar) + i)); 
			yymax = MAX(yymax, YofP(simpar, DM_TBP(simpar) + i)); 
			zzmax = MAX(zzmax, ZofP(simpar, DM_TBP(simpar) + i)); 
		} 
		MPI_Barrier(MPI_COMM(simpar)); 
		DEBUGPRINT("P%d %g %g : %g %g : %g %g\n",MYID(simpar),xxmin,xxmax,yymin,yymax,zzmin,zzmax);
	}
#endif

	if(PMSTATUS(simpar) == PUSH) TREE_HALF(simpar, ANOW(simpar), ASTEP(simpar), KICK);
	else if(PMSTATUS(simpar) == PULL) TREE_HALF(simpar, ANOW(simpar), ASTEP(simpar), PULL);
	else TREE_HALF(simpar, ANOW(simpar), ASTEP(simpar), KICK);


	if(PM_FIRST_HALF(simpar) == 'Y' && PM_SECOND_HALF(simpar) == 'N' && CONT_FLAGPREFOF(simpar) =='Y'){
		if(MYID(simpar) ==0) printf("Entering into the PreFoF\n");
		TIMER_START(1);
		PreFoF(simpar);
		TIMER_STOP(1);
		if(MYID(simpar)==0) printf("%d CPU(PreFoF) = %f\n", STEPCOUNT(simpar), ELAPSED_TIME(1));
	}
	Tree2PMConversion(simpar);
	TIMER_STOP(4);
	if(MYID(simpar)==0)printf("--- Step %d CPU(TREE) = %f\n",STEPCOUNT(simpar), ELAPSED_TIME(4));
}

void HubbleDragDM(SimParameters *simpar, float fact1){
	ptrdiff_t i;
	dmparticletype *bp = DM_BP(simpar);
	for(i=0;i<DM_NP(simpar);i++){
		(bp+i)->vx *= fact1;
		(bp+i)->vy *= fact1;
		(bp+i)->vz *= fact1;
	}
}


void PM_PULL_AND_PUSH(SimParameters *simpar){ 
	float cputime0[5], cputime1[5];
	float anow = ANOW(simpar); float astep = ASTEP(simpar);
	TIMER_START(4);
	EvolFact evolfact; 
	if(BGEXPAND(simpar) == 'Y') evolfact = GetEvolFactor(simpar);
	else evolfact = StaticGetEvolFactor(simpar);
	if(PM_FIRST_HALF(simpar) == 'Y' && PM_SECOND_HALF(simpar)== 'N') FLAG4SAVEXZSLICE(simpar) = 'Y';
	else FLAG4SAVEXZSLICE(simpar) = 'N';
	float vfact1;
	if(PMSTATUS(simpar)==PUSH) vfact1 = evolfact.fact1_push;
	else if(PMSTATUS(simpar)==PULL) vfact1 = evolfact.fact1_pull;
	else vfact1 = evolfact.fact1;
#ifdef DEBUG
	DEBUGPRINT("P%d has velocity factor vfact1 = %g\n",MYID(simpar),vfact1);
#endif
	HubbleDragDM(simpar, vfact1);
	void pmmain(SimParameters *, EvolFact);
	if(GRAVITY(simpar) == 'Y') pmmain(simpar,evolfact);
	TIMER_STOP(4);
	if(MYID(simpar)==0)printf("--- Step %d CPU(PM) = %f\n",STEPCOUNT(simpar), ELAPSED_TIME(4));
}

void Cosmos_GlobalT(SimParameters *simpar, int icont){
	float cputime0[5], cputime1[5];
	REDSHIFT(simpar) = AMAX(simpar)/ANOW(simpar)  -1;


	for(STEPCOUNT(simpar)= STEPNUM(simpar); STEPCOUNT(simpar) <= NSTEP(simpar);STEPCOUNT(simpar)++){


		if(icont !=1) GSAVING_ALL_DATA(simpar);
		icont = 0;
		if(STEPCOUNT(simpar)==1) dumpmain(simpar);


		TIMER_START(3);
/*
		if(ANIM_FLAG(simpar)) mkanimate(simpar);
*/
		/*
		if(OBS_FLAG(simpar)) observerextmain(simpar);
		*/

		PM_FIRST_HALF(simpar) = (PM_SECOND_HALF(simpar) = 'N');
		flagPreFoF(simpar);
		flagsyncpdata(simpar);

		if(CONT_FLAGSYNCP(simpar) == 'Y' || CONT_FLAGPREFOF(simpar) == 'Y'){
			PMSTATUS(simpar) = PULL;
			PM_FIRST_HALF(simpar) = 'Y';
			PM_SECOND_HALF(simpar) = 'N';
		}
		if(CONT_FLAGSYNCP(simpar) == 'N' && CONT_FLAGPREFOF(simpar) == 'N'){
			PM_FIRST_HALF(simpar) = PM_SECOND_HALF(simpar) = 'N';
			PMSTATUS(simpar) = KICK;
		}

halfevolution:
		PM_PULL_AND_PUSH(simpar);
#ifdef DEBUG
		DEBUGPRINT("P%d is at %d %d %d\n",MYID(simpar), STEPCOUNT(simpar), STEPNUM(simpar), NSTEP(simpar));
#endif


		TreeAndPreFoFmain(simpar);

#ifdef DEBUG
		DEBUGPRINT("P%d is at %d %d %d\n",MYID(simpar), STEPCOUNT(simpar), STEPNUM(simpar), NSTEP(simpar));
#endif

#ifdef DEBUG
		DEBUGPRINT("P%d is at %d %d %d\n",MYID(simpar), STEPCOUNT(simpar), STEPNUM(simpar), NSTEP(simpar));
#endif

		if(PM_FIRST_HALF(simpar) == 'Y' && PM_SECOND_HALF(simpar) == 'N'){

			dumpmain(simpar);

			PM_SECOND_HALF(simpar) = 'Y';

			PMSTATUS(simpar) = PUSH;
			goto halfevolution;
		}
		/*
		if(OBS_FLAG(simpar)) observersavemain(simpar);
		*/

		void onestepforwardposition(SimParameters *);
		onestepforwardposition(simpar);

		AllParticleMigrate(simpar);


		ANOW(simpar) += ASTEP(simpar);
		MPI_Barrier(MPI_COMM(simpar));
		DomainDecomp(simpar,0);
#ifdef DEBUG
		DEBUGPRINT("P%d is at %d %d %d\n",MYID(simpar), STEPCOUNT(simpar), STEPNUM(simpar), NSTEP(simpar));
#endif
		TIMER_STOP(3);
		if(MYID(simpar)==0){
			printf("Step %d CPU = %f\n", STEPCOUNT(simpar),ELAPSED_TIME(3)); 
			printf("###########################################\n\n");
		}
	}
}

void GOTPM_Cosmos_GlobalT(SimParameters *simpar, int icont){
	float cputime0[4], cputime1[4];
	REDSHIFT(simpar) = AMAX(simpar)/ANOW(simpar)  -1;


	for(STEPCOUNT(simpar)= STEPNUM(simpar); STEPCOUNT(simpar) <= NSTEP(simpar);STEPCOUNT(simpar)++){

		TIMER_START(3);

		if(icont !=1) GSAVING_ALL_DATA(simpar);
		icont = 0;


		/* This is only for the first step */
		if(STEPCOUNT(simpar)==1) {
			PM_FIRST_HALF(simpar) = 'Y';
			PM_SECOND_HALF(simpar) = 'N';
			dumpmain(simpar);
		}

		{ 
			PM_FIRST_HALF(simpar) = (PM_SECOND_HALF(simpar) = 'N'); 
			flagPreFoF(simpar); 
			flagsyncpdata(simpar);

			if(CONT_FLAGSYNCP(simpar) == 'Y' || CONT_FLAGPREFOF(simpar) == 'Y'){
				PMSTATUS(simpar) = PULL;
				PM_FIRST_HALF(simpar) = 'Y';
				PM_SECOND_HALF(simpar) = 'N';
			}
			if(CONT_FLAGSYNCP(simpar) == 'N' && CONT_FLAGPREFOF(simpar) == 'N'){
				PM_FIRST_HALF(simpar) = PM_SECOND_HALF(simpar) = 'N';
				PMSTATUS(simpar) = KICK;
			}
		}

		PM_PULL_AND_PUSH(simpar);

#ifndef PMonly
		TreeAndPreFoFmain(simpar);
#else
		if(CONT_FLAGPREFOF(simpar) =='Y'){
			if(MYID(simpar) ==0) printf("Entering into the PreFoF\n"); 
			TIMER_START(1); 
			PM2TreeConversion(simpar);
			PreFoF(simpar); 
			Tree2PMConversion(simpar);
			TIMER_STOP(1); 
			if(MYID(simpar)==0) printf("%d CPU(PreFoF) = %f\n", STEPCOUNT(simpar), ELAPSED_TIME(1));
		}
#endif

		dumpmain(simpar);

		if(PM_FIRST_HALF(simpar) == 'Y' && PM_SECOND_HALF(simpar) == 'N'){
			PM_SECOND_HALF(simpar) = 'Y';
			PMSTATUS(simpar) = PUSH;
			PM_PULL_AND_PUSH(simpar);
#ifndef PMonly
			TreeAndPreFoFmain(simpar);
#endif
		}

		/*
		if(OBS_FLAG(simpar)) observersavemain(simpar);
		*/

		void onestepforwardposition(SimParameters *);
		onestepforwardposition(simpar);

		AllParticleMigrate(simpar);


		ANOW(simpar) += ASTEP(simpar);
		MPI_Barrier(MPI_COMM(simpar));
		DomainDecomp(simpar, 0);
#ifdef DEBUG
		DEBUGPRINT("P%d is at %d %d %d\n",MYID(simpar), STEPCOUNT(simpar), STEPNUM(simpar), NSTEP(simpar));
#endif
		TIMER_STOP(3);
		if(MYID(simpar)==0){
			printf("Step %d CPU = %f\n", STEPCOUNT(simpar),ELAPSED_TIME(3)); 
			printf("###########################################\n\n");
		}
	}
}
