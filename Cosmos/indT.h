void Cosmos_IndT(SimParameters *, int );

/* isub % 2^(npow-ii) == 0? 1 : 0 */
/* modulus(isub, 2^(npow-ii)) ==0 ? 1 : 0 */
#define IsNowSphStep(isub,jj,npow) (( isub & (( 1<< (npow-jj))-1)) ==0 ? 1 : 0)
#define isnowdmstep(isub,ii,npow) (( isub & (( 1<< (npow-ii))-1)) ==0 ? 1 : 0)
#define IsNowOneOfTwoSteps(isub,ii,jj,npow) (( isub & (( 1<< (npow-ii))-1)) ==0 ? 1 : IsNowSphStep(isub,jj,npow))
#define isnowstep(isub,ii,npow) (( isub & (( 1<< (npow-ii))-1)) ==0 ? 1 : 0)
#define IsNowOnlyFixedStep(isub,ii,ifix,npow) (( isub & (( 1<< (npow-ii))-1)) ==0 ? 0 : isnowstep(isub,ifix,npow))

#define Teta 0.9
#define GetNDiv(simpar,dist) ((dist)/(GRV_EPSILON(simpar)*Teta))




#ifndef GOTPM
/* This is to use the index memory chip space for the individual time step */
#define GetTsubPower(p) ((p)->indt[0])
#define SetTsubPower(p,v) do {(p)->indt[0] = v; } while(0)
#define GetSphTsubPower(p) ((p)->indt[1])
#define GetSphNextTsubPower(p) ((p)->indt[2])
#define GetFixedPower(p) ((p)->indt[3])
#define SetFixedPower(p,v) do {(p)->indt[3] = v; } while(0)


#define ShiftSphIndT(p,ndiv) ((p)->indt[1] = (((p)->indt[1])+ndiv))
#define SetSphTsubPower(p,v) do {\
	if(v>20) DEBUGPRINT("Strange Value, sphtsbpower=%d\n", v);\
	(p)->indt[1] = v; \
} while(0)
#define SetSphNextTsubPower(p,v) do {\
	if(v>20) DEBUGPRINT("Strange Value, sphtsbpower=%d\n", v);\
	(p)->indt[2] = v; \
} while(0)

#define GPU_SetSphNextTsubPower(p,v) do {\
    (p)->indt[2] = v; \
} while(0)

#define IsStarSinkCandidate(p) ((p)->indt[2])
#define StarSinkCandidate(p) do {\
	    (p)->indt[2] = 1; \
} while(0)
#define NotStarSinkCandidate(p) do {\
	    (p)->indt[2] = 0; \
} while(0)

#define SetNbodyIndT(simpar,bp) {\
	    float pfact = ASTEP(simpar)*NX(simpar);\
	    float vx,vy,vz,vr,dist;\
	    vx = bp->vx;\
	    vy = bp->vy;\
	    vz = bp->vz;\
	    vr = sqrtf(vx*vx+vy*vy+vz*vz);\
	    dist = vr* pfact;\
	    int ndiv = GetNDiv(simpar, dist);\
	    int isubpow = MAX(0,ceil(logf(ndiv)/logf(2.f)));\
	    SetTsubPower(bp,isubpow);\
} while(0)
#define TREE_HALF(simpar,anow, astep, PUSHPULLKICK) do{\
	void TreeMain(SimParameters *, DeterminedEvolFact *);\
	if(MYID(simpar)==0) printf("Now entering into the Tree "#PUSHPULLKICK" of indT Style.\n");\
	if(BGEXPAND(simpar) == 'Y') GetEvolFactorArr(simpar, anow,astep);\
	else StaticGetEvolFactorArr(simpar, anow, astep);\
	CONT_HALFSTEP(simpar) = PUSHPULLKICK;\
	TreeMain(simpar, EvolFact##PUSHPULLKICK);\
}while(0)
#else
/* NULL DEFINITION */
/* This is to use the index memory chip space for the individual time step */
#define GetTsubPower(p) (0)
#define SetTsubPower(p,v) do {} while(0)
#define GetSphTsubPower(p) (0)
#define GetSphNextTsubPower(p) (0)
#define GetFixedPower(p) (0)
#define SetFixedPower(p,v) do {} while(0)


#define ShiftSphIndT(p,ndiv) (0)
#define SetSphTsubPower(p,v) do {\
} while(0)
#define SetSphNextTsubPower(p,v) do {\
} while(0)

#define GPU_SetSphNextTsubPower(p,v) do {\
} while(0)

#define IsStarSinkCandidate(p) (0)
#define StarSinkCandidate(p) do {\
} while(0)
#define NotStarSinkCandidate(p) do {\
} while(0)

#define SetNbodyIndT(simpar,bp) {\
} while(0)

#define TREE_HALF(simpar,anow, astep, PUSHPULLKICK) do{\
	void TreeMain(SimParameters *, DeterminedEvolFact *);\
	/*\
	if(BGEXPAND(simpar) == 'Y') GetEvolFactorArr(simpar, anow,astep);\
	else StaticGetEvolFactorArr(simpar, anow, astep);\
	*/\
	if(BGEXPAND(simpar) == 'Y') {\
    	EvolFact evolfactor = GetEvolFactor(simpar);\
		EvolFact##PUSHPULLKICK[0].fact1 = evolfactor.fact1;\
		EvolFact##PUSHPULLKICK[0].fact2 = evolfactor.fact2;\
		if(MYID(simpar)==0) printf("Now entering into the Tree "#PUSHPULLKICK" of GOTPM Style in Expanding BG.: %g %g\n", evolfactor.fact1, evolfactor.fact2);\
	}\
	else {\
		EvolFact evolfactor = StaticGetEvolFactor(simpar);\
		EvolFact##PUSHPULLKICK[0].fact1 = evolfactor.fact1;\
		EvolFact##PUSHPULLKICK[0].fact2 = evolfactor.fact2;\
		if(MYID(simpar)==0) printf("Now entering into the Tree "#PUSHPULLKICK" of GOTPM Style in nonExpanding BG.\n");\
	}\
	PMSTATUS(simpar) = PUSHPULLKICK;\
	TreeMain(simpar, EvolFact##PUSHPULLKICK);\
}while(0)

/* end of GOTPM */
#endif







#define PM_HALF_PUSH(simpar) do {\
	float anow = ANOW(simpar); float astep = ASTEP(simpar);\
	EvolFact evolfact; void pmmain(SimParameters *, EvolFact);\
	if(BGEXPAND(simpar) == 'Y') evolfact = GetEvolFactor(simpar);\
	else evolfact = StaticGetEvolFactor(simpar);\
	PMSTATUS(simpar) = PUSH;\
	FLAG4SAVEXZSLICE(simpar) = 'Y';\
	if(GRAVITY(simpar) == 'Y') pmmain(simpar, evolfact);\
} while(0)

#define PM_HALF_PULL(simpar) do {\
	float anow = ANOW(simpar); float astep = ASTEP(simpar);\
	EvolFact evolfact; void pmmain(SimParameters *, EvolFact);\
	if(BGEXPAND(simpar) == 'Y') evolfact = GetEvolFactor(simpar);\
	else evolfact = StaticGetEvolFactor(simpar);\
	PMSTATUS(simpar) = PULL;\
	FLAG4SAVEXZSLICE(simpar) = 'N';\
	if(GRAVITY(simpar) == 'Y') pmmain(simpar,evolfact);\
} while(0)


#define PM_HALF(simpar,PUSHPULL) do{\
	EvolFact evolfact; void pmmain(SimParameters *, EvolFact);\
	if(BGEXPAND(simpar) == 'Y') evolfact = GetEvolFactor(simpar);\
	else evolfact = StaticGetEvolFactor(simpar);\
	PMSTATUS(simpar) = PUSHPULL;\
	if(GRAVITY(simpar) == 'Y') pmmain(simpar, evolfact);\
} while(0)


/*
#define PRE_FOF(simpar) do {\
	if(IFLAGPREFOF(simpar)) {\
		PreFoF(simpar);\
		if(MYID(simpar)==0) printf("Now entering into the PreFoF\n");\
		TIMER_START(1);\
		PreFoF(simpar);\
		TIMER_STOP(1);\
		if(MYID(simpar)==0) printf("PreFoF CPU= %f\n",ELAPSED_TIME(1));\
	}\
}while(0)
*/






#define PrintStatus(simpar) do {\
	if(MYID(simpar)==0){\
		FILE *wp;\
		if(wp=fopen("SimulationStatus.output","a")){\
			fprintf(wp,"Stepcount= %4d nsubstep= %3d nowTsubdiv= %3d maxTsubpower= %3d : ", STEPCOUNT(simpar), \
					IndT_NSUBSTEP(simpar), IndT_NOWTSUBDIV(simpar),IndT_MAXTSUBPOWER(simpar));\
			fprintf(wp,"Now at a= %12.7g and z= %12.7g\n",ANOW(simpar),(AMAX(simpar)/ANOW(simpar)-1));\
			fclose(wp);\
		}\
	}\
}while(0)


/*
#define RetrieveIndTData(simpar) do{\
	damin = DAMIN(simpar);\
	iflagPreFoF = IFLAGPREFOF(simpar);\
	iflagsyncpdata = IFLAGSYNCPDATA(simpar);\
	maxTsubpower = IndT_MAXTSUBPOWER(simpar);\
	startTsubdiv = IndT_STARTTSUBDIV(simpar);\
	nsubstepcount = IndT_NSUBSTEPCOUNT(simpar);\
}while(0)
*/



#ifdef UPDATEGLASS

#define SAVING_ALL_DATA(simpar) do {\
	    int suddenstop;\
	    int flagsuddenstop(int,int);\
	    suddenstop = flagsuddenstop(STEPCOUNT(simpar),IndT_NOWTSUBDIV(simpar));\
	    if(suddenstop == 3) { \
			        MPI_Finalize();exit(0); \
			    } \
	    else if( suddenstop) { \
			        void update_glacial_data(SimParameters *); update_glacial_data(&simpar);\
			        jwrite(simpar);\
			        PTYPE(simpar) = TREETYPE;\
			        if(suddenstop == 1) { \
						            MPI_Finalize();exit(0);\
						        } \
			    } \
}while(0)
#else

#define SAVING_ALL_DATA(simpar) do {\
	    int flagsuddenstopInd(SimParameters *);\
	    flagsuddenstopInd(simpar);\
	    if(CONT_FLAGCONTINUE(simpar) == 'E') { \
			MPI_Finalize();exit(0); \
		} \
		else if(CONT_FLAGCONTINUE(simpar) =='D' || CONT_FLAGCONTINUE(simpar) == 'S'){\
			char saveprefix[89];\
			strcpy(saveprefix,RV_FILEPREFIX(simpar));\
			sprintf(RV_FILE(simpar),"%s.%.5d.%.5d",RV_FILEPREFIX(simpar), STEPCOUNT(simpar),\
					IndT_NOWTSUBDIV(simpar));\
			jwrite(simpar); \
			if(CONT_FLAGCONTINUE(simpar) == 'N') { \
				MPI_Finalize();exit(0);\
			} \
		} \
}while(0)
#endif


