#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<mpi.h>
#include<omp.h>


#include "eunha.h"
#include "aux.h"


/*
#define pm2tree(simpar, TYPE,type) do{\
	ptrdiff_t i;\
	TYPE##_TBP(simpar) = (tree##type##particletype*)realloc(TYPE##_TBP(simpar), sizeof(tree##type##particletype)*TYPE##_NP(simpar));\
	for(i=TYPE##_NP(simpar)-1;i>=0;i--){\
		type##particletype tmp = TYPE##_BP(simpar)[i];\
		type##particletype *a = (type##particletype*) ( &(TYPE##_TBP(simpar)[i].u4if) );\
		*a = tmp;\
		SET_FLAG(a,TYPE##flag);\
	}\
}while(0)

#define tree2pm(simpar, TYPE,type) do{\
	ptrdiff_t i;\
	for(i=0;i<TYPE##_NP(simpar);i++){\
		type##particletype tmp = *(type##particletype*)(  &(TYPE##_TBP(simpar)[i].u4if) );\
		type##particletype *a = TYPE##_BP(simpar) + i;\
		*a = tmp;\
	}\
	TYPE##_BP(simpar) = (type##particletype*)realloc(TYPE##_BP(simpar), sizeof(type##particletype)*TYPE##_NP(simpar));\
}while(0)
*/
#define pm2tree(simpar, TYPE,type) do{\
	ptrdiff_t i,nsize,offset;\
	nsize = sizeof(type##particletype);\
	offset = offsetof(tree##type##particletype,u4if);\
	TYPE##_TBP(simpar) = (tree##type##particletype*)realloc(TYPE##_TBP(simpar), sizeof(tree##type##particletype)*TYPE##_NP(simpar));\
	for(i=TYPE##_NP(simpar)-1;i>=0;i--){\
		type##particletype *tmp = TYPE##_BP(simpar) + i;\
		type##particletype *a = (type##particletype*) ( (char*)(TYPE##_TBP(simpar)+i) + offset);\
		memmove(a, tmp, nsize);\
	}\
}while(0)

#define tree2pm(simpar, TYPE,type) do{\
	ptrdiff_t i,nsize, offset;\
	nsize = sizeof(type##particletype);\
	offset = offsetof(tree##type##particletype,u4if);\
	for(i=0;i<TYPE##_NP(simpar);i++){\
		type##particletype *tmp = (type##particletype*) ( (char*)(TYPE##_TBP(simpar)+i) + offset );\
		type##particletype *a = TYPE##_BP(simpar) + i;\
		memmove(a, tmp, nsize);\
	}\
	TYPE##_BP(simpar) = (type##particletype*)realloc(TYPE##_BP(simpar), sizeof(type##particletype)*TYPE##_NP(simpar));\
}while(0)



void PM2TreeConversion(SimParameters *simpar){
	PTYPE(simpar) = TREETYPE;
	if(DM_NP(simpar) >0) pm2tree(simpar, DM,dm);
	if(SPH_NP(simpar) >0) pm2tree(simpar, SPH,sph);
	if(STAR_NP(simpar) >0) pm2tree(simpar, STAR,star);
	if(AGN_NP(simpar) >0) pm2tree(simpar, AGN,agn);
	/*
	ptrdiff_t i;
	DM_TBP(simpar) = (treedmparticletype*)realloc(DM_TBP(simpar), sizeof(treedmparticletype)*DM_NP(simpar));
	DEBUGPRINT("P%d has %ld : %p %p \n",MYID(simpar), DM_NP(simpar), DM_BP(simpar), DM_TBP(simpar));
	for(i=DM_NP(simpar)-1;i>=0;i--){
		dmparticletype tmp = DM_BP(simpar)[i];
		dmparticletype *a = (dmparticletype*) (&(DM_TBP(simpar)[i].u4if));
		*a = tmp;
	}
	*/
}
void Tree2PMConversion(SimParameters *simpar){
	PTYPE(simpar) = PMTYPE;
	if(DM_NP(simpar) >0) tree2pm(simpar, DM,dm);
	if(SPH_NP(simpar) >0) tree2pm(simpar, SPH,sph);
	if(STAR_NP(simpar) >0) tree2pm(simpar, STAR,star);
	if(AGN_NP(simpar) >0) tree2pm(simpar, AGN,agn);
}
