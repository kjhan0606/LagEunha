#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
#define DEFINE_SIM_INDT
#include "CosmosEvolFactor.h"
#undef DEFINE_SIM_INDT
#include "indT.h"





#define MAX(a,b) ((a)>(b) ? (a):(b))
#define MIN(a,b) ((a)<(b) ? (a):(b))


float simparqsimp(float (*func)(SimParameters *, float), float , float , SimParameters *);
float ap(SimParameters *, float );
float HofEz(SimParameters *, float );

float a1ap(SimParameters *simpar, float a){
    return 1./(a*ap(simpar,a));
}
float GetFact2BHF(SimParameters *simpar, float a1, float a2){
    float res1,res2; 
    res1 = 1./(a1*a1*ap(simpar,a1));
    res2 = simparqsimp(a1ap,a1,a2, simpar);
    return res1*res2*NX(simpar);
}
float GetFact2(SimParameters *simpar, float a1, float a2){
    float res1,res2; 
    res1 = 1./(a2*a2*ap(simpar, a2));
    res2 = simparqsimp(a1ap,a1,a2, simpar);
    return res1*res2*NX(simpar);
}
float GetFact2Saitoh(SimParameters *simpar, float a1, float a2, float a3){
	float res1,res2;
	res1 = 1./(a2*a2*ap(simpar, a2));
	res2 = simparqsimp(a1ap,a1,a3, simpar);
	return res1*res2*NX(simpar);
}
float GetFact1(SimParameters *simpar, float a1,float a2){
    float Fact1,omepk;
    float H1,H2;
	/*
    a1 = a1/AMAX(simpar);
    a2 = a2/AMAX(simpar);
	float red1 = 1./a1 - 1;
	float red2 = 1./a2 - 1;
	*/
	/*
    omepk = 1.-OMEP(simpar)-OMEPLAM(simpar);
    H1 = sqrt(OMEP(simpar)/(a1*a1*a1)+OMEPLAM(simpar)+omepk/(a1*a1));
    H2 = sqrt(OMEP(simpar)/(a2*a2*a2)+OMEPLAM(simpar)+omepk/(a2*a2));
	*/
    H1 = HofEz(simpar, a1);
    H2 = HofEz(simpar, a2);
    Fact1 = (a1*a1*a1*H1)/(a2*a2*a2*H2);
    return Fact1;
}
float Hzzp1(SimParameters *simpar, float anow){
	/*
	float red = AMAX(simpar)/anow - 1;
	float result = 1./((1+red)*HofEz(simpar,anow));
	*/
	float result = anow/AMAX(simpar)/HofEz(simpar,anow);
	return result;
}
double GetTime(SimParameters *simpar, float a1,float a2){
	/*
	float z1,z2;
	z2 = AMAX(simpar)/a1 -1;
	z1 = AMAX(simpar)/a2 -1;
	*/
    float omepk = 1.-OMEP(simpar)-OMEPLAM(simpar);
	float res = simparqsimp(Hzzp1,a1,a2, simpar);

	double factor1 = Mpccgs/(H0*HUBBLE(simpar));
	double result = res* factor1;

	return result;
}
double StaticGetRealTime(SimParameters *simpar, float a1,float a2){
	float res = a2-a1;
	res = res * FFTIME(simpar);
	return res;
}
float StaticGetFact1(SimParameters *simpar, float a1,float a2){
    return 1.;
}
float StaticGetFact2(SimParameters *simpar, float a1, float a2){
	return (a2-a1)*NX(simpar);
}
float StaticGetFact2BHF(SimParameters *simpar, float a1, float a2){
    float res1,res2; 
    return (a2-a1)*NX(simpar);
}

void GetEvolFactorArr(SimParameters *simpar, float anow,float astep){
    double afact,bfact,fact1,fact2;
    float da;
    int i,j,k;
	int nowTsubdiv = IndT_NOWTSUBDIV(simpar);
    
#ifdef DEBUG
    DEBUGPRINT("anow/astep %g %g\n",anow,astep);
#endif
    if(IndT_MAXTSUBPOWER(simpar) > (MaxEvolArrSize-1) ){
        DEBUGPRINT("Error in the MaxEvolArrSize %d: %d\n",MaxEvolArrSize, IndT_MAXTSUBPOWER(simpar));
        MPI_Finalize();
        exit(99);
    }
    for(i=0;i<MaxEvolArrSize;i++){
		/*
        da = astep* (1<<(IndT_MAXTSUBPOWER(simpar)-i));
		*/
        da = ASTEP(simpar)/(1<<i);
        EvolFactKICK[i].fact1 = GetFact1(simpar,anow-0.5*da,anow+0.5*da);
        EvolFactPULL[i].fact1 = GetFact1(simpar,anow-0.5*da,anow       );
        EvolFactPUSH[i].fact1 = GetFact1(simpar,anow       ,anow+0.5*da);
        EvolFactKICK[i].fact2 = GetFact2(simpar,anow-0.5*da,anow+0.5*da);
        EvolFactPULL[i].fact2 = GetFact2(simpar,anow-0.5*da,anow       );
        EvolFactPUSH[i].fact2 = GetFact2(simpar,anow       ,anow+0.5*da);
        dlnaKICK[i]           = da/anow;
        dlnaPULL[i]           = 0.5*da/anow;
        dlnaPUSH[i]           = 0.5*da/anow;
        dTKICK[i]           = GetTime(simpar,anow-0.5*da,anow+0.5*da);
        dTPULL[i]           = GetTime(simpar,anow-0.5*da,anow       );
        dTPUSH[i]           = GetTime(simpar,anow       ,anow+0.5*da);
    } 
	if(CONT_HALFSTEP(simpar) == KICK){
		for(i=0;i<MaxEvolArrSize;i++){
			for(j=i;j<MaxEvolArrSize;j++){
				if(IsNowOneOfTwoSteps(nowTsubdiv,i,j,IndT_MAXTSUBPOWER(simpar))){
					for(k=j-1;k<MaxEvolArrSize;k++){
						double a1 = ANOW(simpar) - ASTEP(simpar)*(1./(1<<(j+1)));
						double a2 = ANOW(simpar) + ASTEP(simpar)*(1./(1<<(k+1)));
						/* The number of substeps to the just-passed nbody i-power step */
						int ijdiff = nowTsubdiv & (( 1 << (IndT_MAXTSUBPOWER(simpar) -i))-1);
						/* The number of substeps to the nearest nbody i-power step */
						if(ijdiff ==0){
							EvolFact2ChangeSubStep[i][j][k].fact2 = GetFact2BHF(simpar,a1,ANOW(simpar)) + 
								GetFact2(simpar,ANOW(simpar),a2);
						}
						else if(ijdiff > (1<<(IndT_MAXTSUBPOWER(simpar)-i))/2) {
							EvolFact2ChangeSubStep[i][j][k].fact2 = GetFact2BHF(simpar,a1,a2);
						}
						else if(ijdiff == (1<<(IndT_MAXTSUBPOWER(simpar)-i))/2) {
							EvolFact2ChangeSubStep[i][j][k].fact2 = GetFact2(simpar,a1,ANOW(simpar)) + 
								GetFact2BHF(simpar,ANOW(simpar),a2);
						}
						else {
							EvolFact2ChangeSubStep[i][j][k].fact2 = GetFact2(simpar,a1,a2);
						}
						dlnaKICKChangeSubStep[i][j][k] = (a2-a1)/ANOW(simpar);
						dTKICKChangeSubStep[i][j][k] = GetTime(simpar,a1,a2);
					}
				}
			}
		}
	}
}
void StaticGetEvolFactorArr(SimParameters *simpar, float anow,float astep){
    float da;
    int i,j,k;
    
    if(IndT_MAXTSUBPOWER(simpar) > (MaxEvolArrSize-1) ){
        DEBUGPRINT("Error in the MaxEvolArrSize %d: %d\n",MaxEvolArrSize, IndT_MAXTSUBPOWER(simpar));
        MPI_Finalize();
        exit(99);
    }
    for(i=0;i<MaxEvolArrSize;i++){
		/*
        da = astep* (1<<(IndT_MAXTSUBPOWER(simpar)-i));
		*/
        da = ASTEP(simpar)/(1<<i);
        EvolFactKICK[i].fact1 = StaticGetFact1(simpar, anow-0.5*da,anow+0.5*da);
        EvolFactPULL[i].fact1 = StaticGetFact1(simpar, anow-0.5*da,anow       );
        EvolFactPUSH[i].fact1 = StaticGetFact1(simpar, anow       ,anow+0.5*da);
        EvolFactKICK[i].fact2 = StaticGetFact2(simpar, anow-0.5*da,anow+0.5*da);
        EvolFactPULL[i].fact2 = StaticGetFact2(simpar, anow-0.5*da,anow       );
        EvolFactPUSH[i].fact2 = StaticGetFact2(simpar, anow       ,anow+0.5*da);
        dlnaKICK[i]           = da; /* This is time in unit of freefall time*/
        dlnaPULL[i]           = 0.5*da; /* This is time in unit of freefall time*/
        dlnaPUSH[i]           = 0.5*da; /* This is time in unit of freefall time*/
        dTKICK[i]           = StaticGetRealTime(simpar, anow-0.5*da,anow+0.5*da);
        dTPULL[i]           = StaticGetRealTime(simpar,anow-0.5*da,anow       );
        dTPUSH[i]           = StaticGetRealTime(simpar,anow       ,anow+0.5*da);
    }
	if(CONT_HALFSTEP(simpar) == KICK){
		for(i=0;i<MaxEvolArrSize;i++){
			for(j=i;j<MaxEvolArrSize;j++){
				if(IsNowOneOfTwoSteps(nowTsubdiv,i,j,IndT_MAXTSUBPOWER(simpar))){
					for(k=j-1;k<MaxEvolArrSize;k++){
						double a1 = ANOW(simpar) - ASTEP(simpar)*(1./(1<<(j+1)));
						double a2 = ANOW(simpar) + ASTEP(simpar)*(1./(1<<(k+1)));
						EvolFact2ChangeSubStep[i][j][k].fact2 = StaticGetFact2(simpar, a1,a2);
						dlnaKICKChangeSubStep[i][j][k] = (a2-a1);
						dTKICKChangeSubStep[i][j][k] = StaticGetRealTime(simpar,a1,a2);
					}
				}
			}
		}
	}
}
EvolFact StaticGetEvolFactor(SimParameters *simpar){
	float anow, astep;
	anow = ANOW(simpar); astep = ASTEP(simpar);
	static float fact1,fact2;
    EvolFact evolfactor;
    fact1 = GetFact1(simpar,anow-0.5*astep,anow+0.5*astep);
    fact2 = GetFact2(simpar,anow-0.5*astep,anow+0.5*astep);

    evolfactor.fact1 = fact1;
    evolfactor.fact2 = fact2;
    evolfactor.pfact = astep * NX(simpar); /* simpar.nx = simpar.ny = simpar.nz */
    evolfactor.fact1_push = StaticGetFact1(simpar, anow          ,anow+0.5*astep);
    evolfactor.fact1_pull = StaticGetFact1(simpar, anow-0.5*astep,anow          );
    evolfactor.fact2_push = StaticGetFact2(simpar, anow          ,anow+0.5*astep);
    evolfactor.fact2_pull = StaticGetFact2(simpar, anow-0.5*astep,anow          );

    return evolfactor;
}
/* These factors should be used for the PM force update only. */
EvolFact GetEvolFactor(SimParameters *simpar){
    EvolFact evolfactor;
	float anow = ANOW(simpar);
	float astep = ASTEP(simpar);

	/* For the cosmological expansion, the PM force should be modified from the original factors
	 * because the hubble flow effect is considered in the tree correction */
	/*################################################################# */
	/*################################################################# */
	/*################################################################# */
	/* 1/a2^2/a2p \int_a1^a2 should be changed to 1/a1^2/a1p\int_a1^a2 */
	/*################################################################# */
	/*################################################################# */
	/*################################################################# */
	/*################################################################# */

    evolfactor.fact1 = GetFact1(simpar,anow-0.5*astep,anow+0.5*astep);;
    evolfactor.fact2 = GetFact2(simpar,anow-0.5*astep,anow+0.5*astep);;
    evolfactor.pfact = astep * NX(simpar); /* simpar.nx = simpar.ny = simpar.nz */
    evolfactor.fact1_push = GetFact1(simpar,anow          ,anow+0.5*astep);
    evolfactor.fact1_pull = GetFact1(simpar,anow-0.5*astep,anow          );
    evolfactor.fact2_push = GetFact2(simpar,anow          ,anow+0.5*astep);
    evolfactor.fact2_pull = GetFact2(simpar,anow-0.5*astep,anow          );

    return evolfactor;
}


#define MeasureIndT(simpar,type,pfact,TBP){\
	ptrdiff_t _i;\
	PosType mxdist = 0,tmaxdist,tmindist;\
	for(_i=0;_i<type##_NP(simpar);_i++){\
		float vx,vy,vz;\
		vx = (type##_##TBP(simpar)+_i)->vx;\
		vy = (type##_##TBP(simpar)+_i)->vy;\
		vz = (type##_##TBP(simpar)+_i)->vz;\
		vr = sqrtf(vx*vx+vy*vy+vz*vz);\
		dist = vr*pfact;\
		mxdist = MAX(mxdist, dist);\
		ndiv = GetNDiv(simpar,dist);\
		isubpow = MAX(0,ceil(logf(ndiv)/logf(2.f)));\
		SetTsubPower(type##_##TBP(simpar)+_i,isubpow);\
		IndT_MAXTSUBPOWER(simpar) = MAX(IndT_MAXTSUBPOWER(simpar),isubpow);\
	}\
	MPI_Reduce(&mxdist,&tmaxdist, 1, MPI_POSTYPE, MPI_MAX, 0, MPI_COMM(simpar));\
	MPI_Reduce(&mxdist,&tmindist, 1, MPI_POSTYPE, MPI_MIN, 0, MPI_COMM(simpar));\
	if(MYID(simpar)==0) printf("IndT max/min dist %g %gof "#type" from np= %ld \n",tmaxdist,tmindist, type##_NP(simpar));\
}
#define HydroIndT(simpar,type,TBP){\
	ptrdiff_t _i;\
	for(_i=0;_i<type##_NP(simpar);_i++){\
		IndT_MAXTSUBPOWER(simpar) = MAX(IndT_MAXTSUBPOWER(simpar),GetSphTsubPower(type##_##TBP(simpar)+_i));\
	}\
}
void DetermineIndT(SimParameters *simpar, float astep){
	long long i;
	int isubpow;
	float ndiv;
	float pfact;
	double vr,dist;
	int tmaxTsubpower;
	IndT_MAXTSUBPOWER(simpar) = 0;
	pfact = astep*NX(simpar);

	if(PTYPE(simpar) == TREETYPE){
		MeasureIndT(simpar, STAR,pfact,TBP);
		MeasureIndT(simpar,SPH,pfact,TBP);
		MeasureIndT(simpar, DM,pfact,TBP);
		MeasureIndT(simpar, AGN,pfact,TBP);
		if(STEPCOUNT(simpar) !=1) HydroIndT(simpar, SPH,TBP);
	}
	else {
		MeasureIndT(simpar,STAR,pfact,BP);
		MeasureIndT(simpar,SPH,pfact,BP);
		MeasureIndT(simpar,DM,pfact,BP);
		MeasureIndT(simpar,AGN,pfact,BP);
		if(STEPCOUNT(simpar) != 1) HydroIndT(simpar, SPH,BP);
	}

	MPI_Allreduce(&IndT_MAXTSUBPOWER(simpar),&tmaxTsubpower,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
	IndT_MAXTSUBPOWER(simpar) = tmaxTsubpower;
	if(MYID(simpar) ==0) {
		IndT_MAXTSUBPOWER(simpar) = tmaxTsubpower;
		DEBUGPRINT("The maximum TsubPower %d with pfact= %g\n",IndT_MAXTSUBPOWER(simpar), pfact);
	}
	/*
	MPI_Bcast(&IndT_MAXTSUBPOWER(simpar),1,MPI_INT,0,MPI_COMM_WORLD);
	*/
}

void CleanIndTMemChip(SimParameters *simpar){
	long i;
	for(i=0;i<DM_NP(simpar);i++){
		SetTsubPower(DM_TBP(simpar)+i,0);
	}
	for(i=0;i<STAR_NP(simpar);i++){
		SetTsubPower(STAR_TBP(simpar)+i,0);
	}
	for(i=0;i<SPH_NP(simpar);i++){
		SetTsubPower(SPH_TBP(simpar)+i,0);
	}
	for(i=0;i<AGN_NP(simpar);i++){
		SetTsubPower(AGN_TBP(simpar)+i,0);
	}
}
void CleanHydroIndT(SimParameters *simpar){
	long i;
	for(i=0;i<STAR_NP(simpar);i++) SetTsubPower(STAR_TBP(simpar)+i,0);
	for(i=0;i<SPH_NP(simpar);i++) SetTsubPower(SPH_TBP(simpar)+i,0);
	for(i=0;i<AGN_NP(simpar);i++) SetTsubPower(AGN_TBP(simpar)+i,0);
}
void CleanHydroNextIndT(SimParameters *simpar){
	long i;
	for(i=0;i<STAR_NP(simpar);i++) SetSphTsubPower(STAR_TBP(simpar)+i,0);
	for(i=0;i<SPH_NP(simpar);i++) SetSphTsubPower(SPH_TBP(simpar)+i,0);
	for(i=0;i<AGN_NP(simpar);i++) SetSphTsubPower(AGN_TBP(simpar)+i,0);
}
void CopyHydroNext2NowIndT(SimParameters *simpar){
	long i;
	for(i=0;i<STAR_NP(simpar);i++) SetTsubPower(STAR_TBP(simpar)+i,GetSphTsubPower(STAR_TBP(simpar)+i));
	for(i=0;i<SPH_NP(simpar);i++) SetTsubPower(SPH_TBP(simpar)+i,GetSphTsubPower(SPH_TBP(simpar)+i));
	for(i=0;i<AGN_NP(simpar);i++) SetTsubPower(AGN_TBP(simpar)+i,GetSphTsubPower(AGN_TBP(simpar)+i));
}
