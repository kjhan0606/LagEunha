#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>

#include "eunha.h"
#include "kjhrw.h"
void dumpmain(SimParameters *simpar){
	if(( PM_FIRST_HALF(simpar) == 'Y' && PM_SECOND_HALF(simpar) == 'N') || STEPCOUNT(simpar) ==1) {
	/* Now positions and velocities of all simulation particles are 
	 * synchronized */
		if(CONT_FLAGSYNCP(simpar) == 'Y' || STEPCOUNT(simpar) == 1){
			char saveprefix[80];
			strcpy(saveprefix,RV_FILEPREFIX(simpar));
			sprintf(RV_FILE(simpar),"Sync%s",RV_FILEPREFIX(simpar));
			jwrite(simpar);
			strcpy(RV_FILEPREFIX(simpar),saveprefix);
		}
	}
}
