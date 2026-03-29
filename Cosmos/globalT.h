#define GSAVING_ALL_DATA(simpar) do {\
	    int suddenstop;\
	    int flagsuddenstopGOTPM(SimParameters *);\
	    flagsuddenstopGOTPM(simpar);\
	    if(CONT_FLAGCONTINUE(simpar) == 'E') { \
			MPI_Finalize();exit(0); \
		} \
		else if(CONT_FLAGCONTINUE(simpar) =='D' || CONT_FLAGCONTINUE(simpar) == 'S'){\
			char saveprefix[80];\
			strcpy(saveprefix,RV_FILEPREFIX(simpar));\
			sprintf(RV_FILE(simpar),"%s",RV_FILEPREFIX(simpar));\
			jwrite(simpar); \
			if(CONT_FLAGCONTINUE(simpar) == 'S') { \
				if(MYID(simpar) ==0) DEBUGPRINT0("Now Closing Due to Stop Sign\n");\
				MPI_Finalize();exit(0);\
			} \
		} \
}while(0)

/*
#define APM_PULL_AND_PUSH(simpar) do {\
	float anow = ANOW(simpar); float astep = ASTEP(simpar);\
	EvolFact evolfact; void pmmain(SimParameters *, EvolFact);\
	if(BGEXPAND(simpar) == 'Y') evolfact = GetEvolFactor(simpar);\
	else evolfact = StaticGetEvolFactor(simpar);\
	PMSTATUS(simpar) = PULL;\
	if(PM_FIRST_HALF(simpar) == 'Y' && PM_SECOND_HALF(simpar)== 'N') FLAG4SAVEXZSLICE(simpar) = 'Y';\
	else FLAG4SAVEXZSLICE(simpar) = 'N';\
	debugparticles(simpar); DEBUGPRINT("P%d has %c\n",MYID(simpar), GRAVITY(simpar));\
	if(GRAVITY(simpar) == 'Y') pmmain(simpar,evolfact);\
} while(0)
*/
