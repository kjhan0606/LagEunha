#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>

#include "eunha.h"
#include "indT.h"
#include "CosmosEvolFactor.h"
#include "Treemain.h"
#include "timerutil.h"

/*
float cputime0[10];
float cputime1[10];
*/

#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)? (a):(b)

#define InitSphTsub(simpar,TYPE,type){\
	ptrdiff_t i;\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	for(i=0;i<TYPE##_NP(simpar);i++){\
		SetSphTsubPower(bp+i,0);\
		SetSphNextTsubPower(bp+i,0);\
	}\
}

#define updateindT(simpar, TYPE,subpower) do {\
	ptrdiff_t i;\
	int newMaxTsubPower=0;\
	for(i = 0;i<TYPE##_NP(simpar);i++){\
		int ii = GetTsubPower(TYPE##_TBP(simpar)+i);\
		int jj = GetSphTsubPower(TYPE##_TBP(simpar)+i);\
		int kk = GetSphNextTsubPower(TYPE##_TBP(simpar)+i);\
		if(kk>=jj ){\
			SetSphTsubPower(TYPE##_TBP(simpar)+i,kk);\
		}\
		else if(IsNowSphStep(IndT_NOWTSUBDIV(simpar),(jj-1),IndT_MAXTSUBPOWER(simpar))){\
			SetSphTsubPower(TYPE##_TBP(simpar)+i,(jj-1));\
		}\
		jj = GetSphTsubPower(TYPE##_TBP(simpar)+i);\
		newMaxTsubPower = MAX(newMaxTsubPower,ii);\
		newMaxTsubPower = MAX(newMaxTsubPower,jj);\
		SetSphNextTsubPower(TYPE##_TBP(simpar)+i,0);\
		subpower = newMaxTsubPower;\
	}\
}while(0)

void UpdateIndT(SimParameters *simpar){
	int newmaxTsubpower= 0;
	updateindT(simpar,SPH,newmaxTsubpower);
	if(newmaxTsubpower != IndT_MAXTSUBPOWER(simpar)) {
		if(newmaxTsubpower>IndT_MAXTSUBPOWER(simpar)) 
			IndT_NOWTSUBDIV(simpar) = IndT_NOWTSUBDIV(simpar) << (newmaxTsubpower-IndT_MAXTSUBPOWER(simpar));
		else IndT_NOWTSUBDIV(simpar) = IndT_NOWTSUBDIV(simpar) >> (IndT_MAXTSUBPOWER(simpar)-newmaxTsubpower);
		IndT_MAXTSUBPOWER(simpar) = newmaxTsubpower;
	}
	{ 
		int tmaxTsubpower,tnowTsubdiv;
		MPI_Reduce(&IndT_MAXTSUBPOWER(simpar),&tmaxTsubpower,1, MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
		MPI_Reduce(&IndT_NOWTSUBDIV(simpar),&tnowTsubdiv,1, MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
		if(MYID(simpar)==0) {
			IndT_MAXTSUBPOWER(simpar) = tmaxTsubpower;
			IndT_NOWTSUBDIV(simpar) = tnowTsubdiv;
		}
		MPI_Bcast(&IndT_MAXTSUBPOWER(simpar),1, MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&IndT_NOWTSUBDIV(simpar),1, MPI_INT,0,MPI_COMM_WORLD);
		IndT_DAMIN(simpar) = ASTEP(simpar)/(1<<IndT_MAXTSUBPOWER(simpar));
	}
}
#define HubbleFlowForNbodyStep_InitAccel(simpar,TYPE,type,EvolFactArr) do{\
	ptrdiff_t i;\
	int imax;\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	for(i = 0;i<TYPE##_NP(simpar);i++){\
		int ii = GetTsubPower(bp+i);\
		/*\
		int jj = GetSphTsubPower(simpar.type.u.tbp+i);\
		*/\
		(bp+i)->ax = 0;\
		(bp+i)->ay = 0;\
		(bp+i)->az = 0;\
		SetSphNextTsubPower(bp+i,0);\
		/*\
		if(IsNowOneOfTwoSteps(nowTsubdiv,ii,jj,maxTsubpower)){\
		*/\
		if(isnowdmstep(IndT_NOWTSUBDIV(simpar),ii,IndT_MAXTSUBPOWER(simpar))){\
			float fact1 = EvolFactArr[ii].fact1;\
			(bp+i)->vx = fact1*(bp+i)->vx;\
			(bp+i)->vy = fact1*(bp+i)->vy;\
			(bp+i)->vz = fact1*(bp+i)->vz;\
		}\
	}\
} while(0)

#define HubbleFlowForNbodyStep(simpar,TYPE,type,EvolFactArr) do{\
	ptrdiff_t i;\
	int imax;\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	for(i = 0;i<TYPE##_NP(simpar);i++){\
		int ii = GetTsubPower(bp+i);\
		if(isnowdmstep(IndT_NOWTSUBDIV(simpar),ii,IndT_MAXTSUBPOWER(simpar))){\
			float fact1 = EvolFactArr[ii].fact1;\
			(bp+i)->vx = fact1*(bp+i)->vx;\
			(bp+i)->vy = fact1*(bp+i)->vy;\
			(bp+i)->vz = fact1*(bp+i)->vz;\
		}\
	}\
} while(0)

#define HubbleFlowForSmallerStep(simpar, TYPE, type,EvolFactArr) do{\
	ptrdiff_t i;\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	for(i = 0;i<TYPE##_NP(simpar);i++){\
		int ii = GetTsubPower(bp+i);\
		int jj = GetSphTsubPower(bp+i);\
		if(IsNowOneOfTwoSteps(IndT_NOWTSUBDIV(simpar),ii,jj,IndT_MAXTSUBPOWER(simpar))){\
			int ismaller;\
			if(ii>jj) ismaller=ii;\
			else ismaller = jj;\
			float fact1 = EvolFactArr[ismaller].fact1;\
			(bp+i)->vx = fact1*(bp+i)->vx;\
			(bp+i)->vy = fact1*(bp+i)->vy;\
			(bp+i)->vz = fact1*(bp+i)->vz;\
		}\
	}\
} while(0)

#define KickUpdate(simpar, TYPE,type,EvolFactArr) do{\
	int ncount = 0;\
	ptrdiff_t i;\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	for(i = 0;i<TYPE##_NP(simpar);i++){\
		int ii = GetTsubPower(bp+i);\
		int jj = GetSphTsubPower(bp+i);\
		if(ii>jj) jj = ii;\
		if(IsNowSphStep(IndT_NOWTSUBDIV(simpar),jj,IndT_MAXTSUBPOWER(simpar))){\
			int kk = GetSphNextTsubPower(bp+i);\
			if(kk>=jj){\
				float vfact2 = EvolFact2ChangeSubStep[ii][jj][kk].fact2;\
				(bp+i)->vx += (bp+i)->ax*vfact2;\
				(bp+i)->vy += (bp+i)->ay*vfact2;\
				(bp+i)->vz += (bp+i)->az*vfact2;\
				ncount ++;\
			}\
			else if(IsNowSphStep(IndT_NOWTSUBDIV(simpar),(jj-1),IndT_MAXTSUBPOWER(simpar))){\
				float vfact2 = EvolFact2ChangeSubStep[ii][jj][jj-1].fact2;\
				(bp+i)->vx += (bp+i)->ax*vfact2;\
				(bp+i)->vy += (bp+i)->ay*vfact2;\
				(bp+i)->vz += (bp+i)->az*vfact2;\
				ncount ++;\
			}\
			else {\
				float vfact2 = EvolFact2ChangeSubStep[ii][jj][jj].fact2;\
				(bp+i)->vx += (bp+i)->ax*vfact2;\
				(bp+i)->vy += (bp+i)->ay*vfact2;\
				(bp+i)->vz += (bp+i)->az*vfact2;\
			}\
		}\
	}\
	if(ncount>0) printf("P%d has changed ind. time step of %d in %ld %s particles\n",\
			MYID(simpar),ncount,TYPE##_NP(simpar),#type);\
} while(0)

#define PushUpdate(simpar, TYPE,type,EvolFactArr) do{\
	ptrdiff_t i;\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	for(i = 0;i<TYPE##_NP(simpar);i++){\
		int ii = GetTsubPower(bp+i);\
		int jj = GetSphTsubPower(bp+i);\
		int kk = GetSphNextTsubPower(bp+i);\
		if(ii>jj) jj = ii;\
		if(kk>jj) jj = kk;\
		if(IsNowSphStep(IndT_NOWTSUBDIV(simpar),jj,IndT_MAXTSUBPOWER(simpar))){\
			float vfact2 = EvolFactArr[jj].fact2;\
			(bp+i)->vx += (bp+i)->ax*vfact2;\
			(bp+i)->vy += (bp+i)->ay*vfact2;\
			(bp+i)->vz += (bp+i)->az*vfact2;\
		}\
	}\
} while(0)

#define PullUpdate(simpar, TYPE, type, EvolFactArr) do {\
	ptrdiff_t i;\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	for(i = 0;i<TYPE##_NP(simpar);i++){\
		int ii = GetTsubPower(bp+i);\
		int jj = GetSphTsubPower(bp+i);\
		if(ii>jj) jj = ii;\
		if(IsNowSphStep(IndT_NOWTSUBDIV(simpar),jj,IndT_MAXTSUBPOWER(simpar))){\
			float vfact2 = EvolFactArr[jj].fact2;\
			(bp+i)->vx += (bp+i)->ax*vfact2;\
			(bp+i)->vy += (bp+i)->ay*vfact2;\
			(bp+i)->vz += (bp+i)->az*vfact2;\
		}\
	}\
}while(0)

void EntropyPullUpdate(SimParameters *simpar) {
	ptrdiff_t i;
	treesphparticletype *sphp = SPH_TBP(simpar);
	float gammam1=GAS_GAMMA(simpar)-1;
	for(i = 0;i<SPH_NP(simpar);i++){
		int ii = GetTsubPower(sphp+i);
		int jj = GetSphTsubPower(sphp+i);
		if(ii>jj) jj = ii;
		if(IsNowSphStep(IndT_NOWTSUBDIV(simpar),jj,IndT_MAXTSUBPOWER(simpar))){
			float fact1 = dlnaPULL[jj];
			float dEntropy = (sphp+i)->dAs*fact1;
			float dEoverE = dEntropy/(sphp+i)->Entropy;
			float dToverT = dEoverE;
			float gasden = GAS_RHOS2RHOR(simpar)*((sphp+i)->rho);
			(sphp+i)->Entropy += dEntropy; 
//jhshin1
			(sphp+i)->Entropy += (sphp+i)->delta_e;
			float temp = gammam1*(sphp+i)->Entropy*pow((sphp+i)->rho,GAS_GAMMA(simpar)-1) 
				*mHoverkB*(sphp+i)->mu; 
			if(temp < GAS_MINTEMP(simpar)){
				(sphp+i)->Entropy= kBovermH/gammam1*GAS_MINTEMP(simpar)*
					pow((sphp+i)->rho,1-GAS_GAMMA(simpar))/(sphp+i)->mu; 
			}
		}
//jhshin2
	}
}

void EntropyPushUpdate(SimParameters *simpar){
	ptrdiff_t i;
	treesphparticletype *sphp = SPH_TBP(simpar);
	float gammam1=GAS_GAMMA(simpar)-1;
	for(i = 0;i<SPH_NP(simpar);i++){
		int ii = GetTsubPower(sphp+i);
		int jj = GetSphTsubPower(sphp+i);
		int kk = GetSphNextTsubPower(sphp+i);
		if(ii>jj) jj = ii;
		if(kk>jj) jj = kk;
		if(isnowstep(IndT_NOWTSUBDIV(simpar),jj,IndT_MAXTSUBPOWER(simpar))){
			float fact1 = dlnaPUSH[jj];
			float dEntropy = (sphp+i)->dAs*fact1;
			float dEoverE = dEntropy/(sphp+i)->Entropy;
			float dToverT = dEoverE;
			float gasden = GAS_RHOS2RHOR(simpar)*((sphp+i)->rho);
			(sphp+i)->Entropy += dEntropy;
//jhshin1
			(sphp+i)->Entropy += (sphp+i)->delta_e; 
			float temp = gammam1*(sphp+i)->Entropy*pow((sphp+i)->rho,GAS_GAMMA(simpar)-1) 
				*mHoverkB*(sphp+i)->mu;
			if(temp < GAS_MINTEMP(simpar)){ 
				(sphp+i)->Entropy= kBovermH/gammam1*GAS_MINTEMP(simpar)*
					pow((sphp+i)->rho,1-GAS_GAMMA(simpar)) /(sphp+i)->mu; 
			}
//jhshin2
		}
	}
}

void EntropyKickUpdate(SimParameters *simpar) {
	ptrdiff_t i;
	treesphparticletype *sphp = SPH_TBP(simpar);
	float gammam1=GAS_GAMMA(simpar)-1;
	for(i = 0;i<SPH_NP(simpar);i++){
		int ii = GetTsubPower(sphp+i);
		int jj = GetSphTsubPower(sphp+i);
		if(ii>jj) jj = ii;
		float fact1,vfact2,minentropy;
		if(IsNowSphStep(IndT_NOWTSUBDIV(simpar),jj,IndT_MAXTSUBPOWER(simpar))){
			int kk = GetSphNextTsubPower(sphp+i);
			float fact1;
			if(kk>=jj) {
				fact1 = dlnaKICKChangeSubStep[ii][jj][kk];
			}
			else if(IsNowSphStep(IndT_NOWTSUBDIV(simpar),(jj-1),IndT_MAXTSUBPOWER(simpar))){
				fact1 = dlnaKICKChangeSubStep[ii][jj][jj-1];
			}
			else {
				fact1 = dlnaKICKChangeSubStep[ii][jj][jj];
			}


			float dEntropy = (sphp+i)->dAs*fact1;
			float dEoverE = dEntropy/(sphp+i)->Entropy;
			float dToverT = dEoverE;
			float gasden = GAS_RHOS2RHOR(simpar)*((sphp+i)->rho);
			(sphp+i)->Entropy += dEntropy;
//jhshin1 
			(sphp+i)->Entropy += (sphp+i)->delta_e; 
			float temp = gammam1*(sphp+i)->Entropy*pow((sphp+i)->rho,GAS_GAMMA(simpar)-1) 
				*mHoverkB*(sphp+i)->mu; 
			if(temp < GAS_MINTEMP(simpar)){ 
				(sphp+i)->Entropy= kBovermH/gammam1*GAS_MINTEMP(simpar)*
					pow((sphp+i)->rho,1-GAS_GAMMA(simpar)) /(sphp+i)->mu; 
			}
//jhshin2
		}
	}
}

void  SphSubStepVelUpdate(SimParameters *simpar,DeterminedEvolFact *EvolFactArr){
	if(CONT_HALFSTEP(simpar) == KICK){
		KickUpdate(simpar,SPH,sph,EvolFactArr);
		if(GAS_TYPE(simpar) == 'S') EntropyKickUpdate(simpar);
	}
	else if(CONT_HALFSTEP(simpar) == HALFPUSH){
		PushUpdate(simpar,SPH,sph,EvolFactArr);
		if(GAS_TYPE(simpar) == 'S') EntropyPushUpdate(simpar);
	}
	else{
		PullUpdate(simpar,SPH,sph,EvolFactArr);
		if(GAS_TYPE(simpar) == 'S') EntropyPullUpdate(simpar);
	}
}
/*
*/

#define Init_Accel(simpar, TYPE,type) do{\
	ptrdiff_t i;\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	for(i = 0;i<TYPE##_NP(simpar);i++){\
		(bp+i)->ax = 0;\
		(bp+i)->ay = 0;\
		(bp+i)->az = 0;\
		SetSphNextTsubPower(bp+i,0);\
	}\
} while(0)
#define CheckAccel(simpar, TYPE,type) do {\
	float ampmax = 0;\
	int ii;\
	ptrdiff_t i;\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	for(i = 0;i<TYPE##_NP(simpar);i++){\
		float amp = (bp+i)->ax*(bp+i)->ax;\
		amp +=  (bp+i)->ay*(bp+i)->ay;\
		amp +=  (bp+i)->az*(bp+i)->az;\
		if(amp > ampmax) {\
			ampmax = amp;\
			ii = i;\
		}\
	}\
	ampmax = sqrt(ampmax);\
	if(TYPE##_NP(simpar)>0)printf("P%d has maximum accel amplitude %g of %s at %ld\n",\
			MYID(simpar),ampmax,#type,ii);\
} while(0)


int FindMaxTsubPowrofDM(SimParameters *simpar){
	int maxdmpower=0;
	int tmaxdmpower;
	long long i;
	for(i=0;i<DM_NP(simpar);i++){
		int isub = GetTsubPower(DM_TBP(simpar)+i);
		if(isub>maxdmpower) maxdmpower = isub;
	}
	for(i=0;i<SPH_NP(simpar);i++){
		int isub = GetTsubPower(SPH_TBP(simpar)+i);
		if(isub>maxdmpower) maxdmpower = isub;
	}
	for(i=0;i<STAR_NP(simpar);i++){
		int isub = GetTsubPower(STAR_TBP(simpar)+i);
		if(isub>maxdmpower) maxdmpower = isub;
	}
	MPI_Reduce(&maxdmpower,&tmaxdmpower,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
	MPI_Bcast(&tmaxdmpower,1,MPI_INT,0,MPI_COMM_WORLD);
	return tmaxdmpower;
}

#define SaitohNbodyMultiStepKick(simpar, TYPE,type,evolfactor) do {\
	ptrdiff_t i;\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	for(i=0;i<TYPE##_NP(simpar);i++){\
		int ii = GetTsubPower(bp+i);\
		int jj = GetSphTsubPower(bp+i);\
		if(IsNowOneOfTwoSteps(IndT_NOWTSUBDIV(simpar),ii,jj,IndT_MAXTSUBPOWER(simpar))){\
			/*\
			if(jj>=ii) {\
				float fact2 = EvolFact2SaitohScheme[ii][jj].fact2;\
				(simpar.type.u.tbp+i)->vx += fact2*(simpar.type.u.tbp+i)->ax;\
				(simpar.type.u.tbp+i)->vy += fact2*(simpar.type.u.tbp+i)->ay;\
				(simpar.type.u.tbp+i)->vz += fact2*(simpar.type.u.tbp+i)->az;\
			}\
			else {\
			}\
			*/\
			float fact2;\
			if(jj>=ii) fact2 = evolfactor[jj];\
			else fact2 = evolfactor[ii];\
			(bp+i)->vx += fact2*(bp+i)->ax;\
			(bp+i)->vy += fact2*(bp+i)->ay;\
			(bp+i)->vz += fact2*(bp+i)->az;\
		}\
	}\
} while(0)
#define NbodyStepVelUpdate(simpar, TYPE,type,evolfactor) do{\
	ptrdiff_t i;\
	double mean,std,amp;\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	for(i=0;i<TYPE##_NP(simpar);i++){\
		int ii = GetTsubPower(bp+i);\
		if(isnowdmstep(IndT_NOWTSUBDIV(simpar),ii,IndT_MAXTSUBPOWER(simpar)))\
		{\
			float fact2 = evolfactor[ii].fact2;\
			(bp+i)->vx += fact2*(bp+i)->ax;\
			(bp+i)->vy += fact2*(bp+i)->ay;\
			(bp+i)->vz += fact2*(bp+i)->az;\
		}\
	}\
} while(0)

/*
void NbodyVelUpdate(SimParameters *simpar,DeterminedEvolFact *evolfactor){
	NbodyStepVelUpdate(simpar, DM,dm,evolfactor);
	if(GAS_SPHFLAG(simpar) == 'Y' && SPH_NP(simpar) >0) NbodyStepVelUpdate(simpar, SPH,sph,evolfactor);
	NbodyStepVelUpdate(simpar, STAR,star,evolfactor);
	NbodyStepVelUpdate(simpar, AGN,agn,evolfactor);
}
*/

#ifndef GOTPM

void TreeMain(SimParameters *simpar,
		/*
		DeterminedEvolFact *evolfactkick,
		DeterminedEvolFact *evolfactpush,
		DeterminedEvolFact *evolfactpull
		*/
		DeterminedEvolFact *evolfactor
		){
	long nplong;
	int myid;
	long i;
	void treecorrection(SimParameters *,  DeterminedEvolFact *);
	float cputime0[10];
	float cputime1[10];

	/*
	if(CONT_HALFSTEP(simpar) == HALFPUSH) evolfactor = evolfactpush;
	else if(CONT_HALFSTEP(simpar) == HALFPULL) evolfactor = evolfactpull;
	else if(CONT_HALFSTEP(simpar) == KICK) evolfactor = evolfactkick;
	*/

	if(STEPCOUNT(simpar) < 10) THETA(simpar) = 0.3;
	else if(STEPCOUNT(simpar) < 20) THETA(simpar) = 0.35;
	else if(STEPCOUNT(simpar) < 60) THETA(simpar) = 0.4;
	else THETA(simpar) = 0.45;


	if(MYID(simpar)==0){
		printf("Tree fact2 = %g %g %g at a= %g\n",
				evolfactor[0].fact2, evolfactor[1].fact2, evolfactor[2].fact2,ANOW(simpar));
		fflush(stdout);
	}
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
		DEBUGPRINT("P%d %g %g : %g %g : %g %g ::: with np= %ld\n",
				MYID(simpar),xxmin,xxmax,yymin,yymax,zzmin,zzmax, DM_NP(simpar));

	}
#endif

	Init_Accel(simpar,DM,dm);
	Init_Accel(simpar,SPH,sph);
	Init_Accel(simpar,STAR,star);
	Init_Accel(simpar,AGN,agn);

#ifdef DEBUG
	{
		DoDeInfo *dd = DM_DDINFO(simpar);
		ptrdiff_t nddinfo = NDDINFO(simpar);
		for(i=0;i<nddinfo;i++){
			DEBUGPRINT("P%d has %g %g : %g %g : %g %g\n",MYID(simpar),dd->lgroup.xyz.xmin, dd->lgroup.xyz.xmax,
					dd->lgroup.xyz.ymin, dd->lgroup.xyz.ymax,
					dd->lgroup.xyz.zmin, dd->lgroup.xyz.zmax);
			dd ++;
		}
	}
#endif
	HubbleFlowForNbody(simpar, evolfactor);

	int maxTsubPowerofDM =FindMaxTsubPowrofDM(simpar);
	if(isnowdmstep(IndT_NOWTSUBDIV(simpar),maxTsubPowerofDM,IndT_MAXTSUBPOWER(simpar))){
		if(GRAVITY(simpar) == 'Y'){
			TIMER_START(2);
			if(MYID(simpar)==0){
				printf("now tree nbody force update\n");
			}
			treecorrection(simpar, evolfactor);
			TIMER_STOP(2);
#ifdef DEBUG
			if(MYID(simpar)==0) DEBUGPRINT("CPU(Tree Correction) = %f\n",ELAPSED_TIME(2));
#endif
		}


		/* This is the external force part */
		if(EXTERNALFORCE(simpar)!=0){
			if(SIMMODEL(simpar)==Static){
				void External_Force(float,float,float);
				float X00,Y00,Z00;
				if(MYID(simpar)==0){
					printf("now external force update\n");
				}
				X00 = NX(simpar)/2.;
				Y00 = NY(simpar)/2.;
				Z00 = NZ(simpar)/2.;
				External_Force(X00,Y00,Z00);
			}
			else if(SIMMODEL(simpar)==RT){
				void RTExternForce(SimParameters *);
				RTExternForce(simpar);
			}
		}

	}
	/* SPH Force measurement */
	if(GAS_TYPE(simpar) == 'S'){

		TIMER_START(1);
		if(CONT_HALFSTEP(simpar) != KICK) {
			InitSphTsub(simpar,SPH,sph);
		}

		if(0){
			float velmax,velmin;
			velmax = -1E20; velmin = 1e20;
			treesphparticletype *bp = SPH_TBP(simpar);
			for(i=0;i<SPH_NP(simpar);i++){
				velmax = MAX(velmax,(bp+i)->vx);
				velmin = MIN(velmin,(bp+i)->vx);
			}
			printf("P%d vel min/max = %g %g \n",MYID(simpar),velmin,velmax);
			MPI_Finalize();exit(99);
		}
		void sphmain(); sphmain();



		TIMER_STOP(1);
		if(MYID(simpar)==0) fprintf(stdout,"SPH CPU= %f \n", ELAPSED_TIME(1));
	}

	{
		NbodyStepVelUpdate(simpar, DM,dm,evolfactor);
		if(GAS_TYPE(simpar) == 'S' && SPH_NP(simpar) >0) {
			NbodyStepVelUpdate(simpar, SPH,sph,evolfactor);
			NbodyStepVelUpdate(simpar, STAR,star,evolfactor);
			NbodyStepVelUpdate(simpar, AGN,agn,evolfactor);
		}
	}

	if(GAS_TYPE(simpar) == 'S') SphSubStepVelUpdate(simpar,evolfactor);
	if(GAS_TYPE(simpar) == 'S') UpdateIndT(simpar);


	/*
	if(simpar.star.tnp >0) {
		void mergingstar(SimParameters);
		mergingstar(simpar);
	}
	*/
}

#else
void TreeMain(SimParameters *simpar,
		/*
		DeterminedEvolFact *evolfactkick,
		DeterminedEvolFact *evolfactpush,
		DeterminedEvolFact *evolfactpull
		*/
		DeterminedEvolFact *evolfactor
		){
	int myid;
	ptrdiff_t i;
	void treecorrection(SimParameters *,  DeterminedEvolFact *);
	float cputime0[10];
	float cputime1[10];

	/*
	if(CONT_HALFSTEP(simpar) == HALFPUSH) evolfactor = evolfactpush;
	else if(CONT_HALFSTEP(simpar) == HALFPULL) evolfactor = evolfactpull;
	else if(CONT_HALFSTEP(simpar) == KICK) evolfactor = evolfactkick;
	*/

	if(STEPCOUNT(simpar) < 10) THETA(simpar) = 0.3;
	else if(STEPCOUNT(simpar) < 20) THETA(simpar) = 0.35;
	else if(STEPCOUNT(simpar) < 60) THETA(simpar) = 0.4;
	else THETA(simpar) = 0.45;


	if(MYID(simpar)==0){
		printf("Tree fact1/2 = %g %g at a= %g\n",
				evolfactor[0].fact1, evolfactor[0].fact2, ANOW(simpar));
		fflush(stdout);
	}
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
	{
		DoDeInfo *dd = DM_DDINFO(simpar);
		ptrdiff_t nddinfo = NDDINFO(simpar);
		for(i=0;i<nddinfo;i++){
			DEBUGPRINT("P%d has %g %g : %g %g : %g %g\n",MYID(simpar),dd->lgroup.xyz.xmin, dd->lgroup.xyz.xmax,
					dd->lgroup.xyz.ymin, dd->lgroup.xyz.ymax,
					dd->lgroup.xyz.zmin, dd->lgroup.xyz.zmax);
			dd ++;
		}
	}
#endif

#ifndef GOTPM
	{
		/* Hubble Drag for N-body step */
		ptrdiff_t i; int imax; 
		treedmparticletype *bp = DM_TBP(simpar); 
		float fact1 = evolfactor[0].fact1;
		for(i = 0;i<DM_NP(simpar);i++){ 
			(bp+i)->vx = fact1*(bp+i)->vx; 
			(bp+i)->vy = fact1*(bp+i)->vy; 
			(bp+i)->vz = fact1*(bp+i)->vz;
		}
	}
#endif

	TIMER_START(2);
	if(MYID(simpar)==0){
		printf("now GOTPM-type tree nbody force update\n");
	}


	treecorrection(simpar, evolfactor);

	TIMER_STOP(2);
	if(MYID(simpar)==0) fprintf(stdout,"Treecorrectingtime CPU= %f \n", ELAPSED_TIME(2));

}
#endif
