#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include "eunha.h"
#include "cosmosic.h"

void CosmologicalIC(SimParameters *simpar){


	if(PORDER(simpar)==1){
		Zeldovichmain(simpar);
	}
	else if(PORDER(simpar)==2){
		TwoLPTmain(simpar);
	}
	if(DM_NP(simpar)){
		ptrdiff_t i;
		for(i=0;i<DM_NP(simpar);i++){
			CLEAR_FLAG(DM_BP(simpar)+i);
			SET_FLAG(DM_BP(simpar)+i, DMflag);
			UNSET_FLAG(DM_BP(simpar)+i, BoundaryGhostflag);
		}
	}
	if(SPH_NP(simpar)){
		ptrdiff_t i;
		for(i=0;i<SPH_NP(simpar);i++){
			CLEAR_FLAG(SPH_BP(simpar)+i);
			SET_FLAG(SPH_BP(simpar)+i, SPHflag);
			UNSET_FLAG(SPH_BP(simpar)+i, BoundaryGhostflag);
		}
	}
	NPSUM(simpar) = DM_NP(simpar) + SPH_NP(simpar);
}
