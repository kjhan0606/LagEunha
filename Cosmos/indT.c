#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
//#include "mpirms.h"
#include "cosmology.h"
#include "indT.h"
#include "cosmosic.h"
#include "CosmosEvolFactor.h"

#include "aux.h"
#include "flags.h"
#include "fft.h"
#include "kjhrw.h"

void DetermineIndT(SimParameters *, float );

void Cosmos_IndT(SimParameters *simpar, int icont){
	IndT_NSUBSTEPCOUNT(simpar) = 0;
	REDSHIFT(simpar) = AMAX(simpar)/ANOW(simpar)  -1;
	/* This is the outer loop for global time step */
	for(STEPCOUNT(simpar)= STEPNUM(simpar); STEPCOUNT(simpar) <= NSTEP(simpar);STEPCOUNT(simpar)++){
		if(STEPCOUNT(simpar) == 1) dumpmain(simpar);
		int startTsubdiv;
		if(icont != 1){
			PM_HALF(simpar, PUSH);
			PM2TreeConversion(simpar);
			startTsubdiv = 0;
			float astep = ASTEP(simpar);
			DetermineIndT(simpar, astep);
		}
		else {
			startTsubdiv = IndT_NOWTSUBDIV(simpar);
			if(MYID(simpar) ==0) 
				printf("read maxTsubpower/nowTsubdiv= %d/%d\n",
						IndT_MAXTSUBPOWER(simpar),IndT_NOWTSUBDIV(simpar));
		}
		if(MYID(simpar)==0){
			printf("The maximum subpower of indT = %d at a= %g : %d in nstep %d : %d\n",
					IndT_MAXTSUBPOWER(simpar),ANOW(simpar), STEPCOUNT(simpar), 
					NSTEP(simpar),IndT_IFLAGSYNCPDATA(simpar));
		}
	
		IndT_DAMIN(simpar) = ASTEP(simpar)/(1<< IndT_MAXTSUBPOWER(simpar));
		IndT_NSUBSTEP(simpar) = IndT_NSUBSTEPCOUNT(simpar) = (1 << IndT_MAXTSUBPOWER(simpar));
	
		/* This is the inner loop for the individual time step */
		for(IndT_NOWTSUBDIV(simpar)= startTsubdiv; IndT_NOWTSUBDIV(simpar)<=(1<<IndT_MAXTSUBPOWER(simpar)); IndT_NOWTSUBDIV(simpar)++){
			REDSHIFT(simpar) = AMAX(simpar)/ANOW(simpar) - 1;
			PrintStatus(simpar);
			IndT_ISUBSTEP(simpar) = IndT_NOWTSUBDIV(simpar);

			IndT_IFLAGPREFOF(simpar)    = flagPreFoF(simpar);
			IndT_IFLAGSYNCPDATA(simpar) = flagsyncpdata(simpar);

			if(icont !=1) SAVING_ALL_DATA(simpar);
			icont = 0;
#ifdef DEBUG
			DEBUGPRINT("P%d has np = %ld\n", MYID(simpar), DM_NP(simpar));
#endif

			if(IndT_NOWTSUBDIV(simpar) ==0) {
				TREE_HALF(simpar, ANOW(simpar), IndT_DAMIN(simpar), PUSH);
			}
			else if(IndT_NOWTSUBDIV(simpar) == (1<<IndT_MAXTSUBPOWER(simpar))){
				TREE_HALF(simpar, ANOW(simpar), IndT_DAMIN(simpar), PULL);
			}
			else {
				if(IndT_IFLAGPREFOF(simpar) || IndT_IFLAGSYNCPDATA(simpar)){
					TREE_HALF(simpar, ANOW(simpar), IndT_DAMIN(simpar), PULL);
					DumpMainAndPreFoF(simpar);
					TREE_HALF(simpar, ANOW(simpar), IndT_DAMIN(simpar), PUSH);
				}
				else TREE_HALF(simpar, ANOW(simpar), IndT_DAMIN(simpar), KICK);

			}
			/*
			if(0) OBSERVER_S;
			*/


			if(IndT_NOWTSUBDIV(simpar) != (1<<IndT_MAXTSUBPOWER(simpar)) ){
				Evol_PFACT(simpar) = IndT_DAMIN(simpar) * NX(simpar);
				void treeonestepforwardposition(SimParameters *);
				treeonestepforwardposition(simpar);
				TreeAllParticleMigrate(simpar);
				ANOW(simpar) += IndT_DAMIN(simpar);
			}
			IndT_TNUMCOUNT(simpar) ++;
		}
		CleanIndTMemChip(simpar);
		Tree2PMConversion(simpar);
		PM_HALF(simpar, PULL);
		DumpMainAndPreFoF(simpar);
		/*
		if(SIMMODEL(simpar) != Cosmos) DomainDecomp;
		*/
	}
/*
	mpi_fftw_finalize(FFTWINFO(simpar));
*/
	mpi_fftw_finalize();
	exit(0);

}
