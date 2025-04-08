typedef struct EvolFact{
	float afact,bfact,pfact;
	float fact1,fact2;
	float fact1_push,fact2_push;
	float fact1_pull,fact2_pull;
} EvolFact;
typedef struct DeterminedEvolFact{
	float fact1,fact2;
} DeterminedEvolFact;

float GetFact1(SimParameters *, float ,float );
void CleanIndTMemChip(SimParameters *);

EvolFact GetEvolFactor(SimParameters *);
EvolFact StaticGetEvolFactor(SimParameters *);
void GetEvolFactorArr(SimParameters *, float, float);
void StaticGetEvolFactorArr(SimParameters *, float, float);

#define MaxEvolArrSize 20
#ifdef DEFINE_SIM_INDT
#define EXTERN
#else
#define EXTERN extern
#endif




EXTERN int  maxTsubpower;
EXTERN int  nowTsubdiv;
EXTERN int  nsubstepcount;

EXTERN double  damin;


EXTERN DeterminedEvolFact EvolFactKICK[MaxEvolArrSize];
EXTERN DeterminedEvolFact EvolFactPULL[MaxEvolArrSize];
EXTERN DeterminedEvolFact EvolFactPUSH[MaxEvolArrSize];
EXTERN double dlnaKICK[MaxEvolArrSize];
EXTERN double dlnaPULL[MaxEvolArrSize];
EXTERN double dlnaPUSH[MaxEvolArrSize];
EXTERN double dTKICK[MaxEvolArrSize];
EXTERN double dTPULL[MaxEvolArrSize];
EXTERN double dTPUSH[MaxEvolArrSize];
EXTERN DeterminedEvolFact EvolFact2ChangeSubStep[MaxEvolArrSize][MaxEvolArrSize][MaxEvolArrSize];
EXTERN double dlnaKICKChangeSubStep[MaxEvolArrSize][MaxEvolArrSize][MaxEvolArrSize];
EXTERN double dTKICKChangeSubStep[MaxEvolArrSize][MaxEvolArrSize][MaxEvolArrSize];

EXTERN DeterminedEvolFact EvolFact2SaitohScheme[MaxEvolArrSize][MaxEvolArrSize];
#undef EXTERN














#define DumpMainAndPreFoF(simpar) do{\
	/*\
	void dumpmain(int);\
	float PreFoF(treepmparticletype *, long , int);\
	if(IFLAGSYNCPDATA(simpar)) dumpmain(simpar);\
    if(IFLAGPREFOF(simpar) && DM_NP(simpar) > 0){\
		if(PTYPE(simpar) == PMTYPE) {\
			pmparticles2treeparticles;\
			PreFoF(simpar.dm.u.tbp,simpar.dm.np,1);\
			treeparticles2pmparticles;\
		}\
		else\
		PreFoF(simpar.dm.u.tbp,simpar.dm.np,1);\
	}\
	*/\
}while(0)

