#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include "eunha.h"
#include "sph.h"
#include "indT.h"
#include "hydroBasicCell.h"
#include "aux.h"


void StaticDetermineSphFactor();

void RTEquilibrium(SimParameters *simpar){
	size_t SPHBuildLinkedList(SimParameters *);
	void DestroyLinkedList(SimParameters *);
	void FindSphDensity(SimParameters *);

	void PM2TreeConversion(SimParameters *);
	PM2TreeConversion(simpar);

	int Numnear = SPH_NUMNEAR(simpar);
	StaticDetermineSphFactor();
	SPHBuildLinkedList(simpar);
	FindSphDensity(simpar);
	DestroyHydroLinkedCell(simpar);
	int i;
	for(i=0;i<SPH_NP(simpar);i++){
		(SPH_TBP(simpar)+i)->temp = RT_INITTEMP(simpar)/((SPH_TBP(simpar)+i)->rho);
		SPH_TBP(simpar)[i].Entropy = kBovermH*SPH_TBP(simpar)[i].temp/SPH_TBP(simpar)[i].mu *
					pow(SPH_TBP(simpar)[i].rho, 1 - GAS_GAMMA(simpar));
	}
	Tree2PMConversion(simpar);
}


/*
void  mkGlacialRTInit(SimParameters *simpar){
	float A,lambda,x,y;
	float vfact=0.1;
	A=0.025 * vfact;
	lambda = 1./3.;
	float zstart,zfinal;
	zstart = simpar->zmin;
	zfinal = simpar->zmax;
		void read_glacial_and_inflate(char *, float , float);
		if(simpar->myid==0) printf("P%d just before reading Glacial particle data\n",simpar->myid);
		read_glacial_and_inflate(simpar->glacialheaderfilename,zstart,zfinal);
	}

	{
		long long npsph = DM_NP(simpar) +SPH_NP(simpar);
		SPH_BP(simpar) = (sphparticletype *)Realloc(SPH_BP(simpar),sizeof(sphparticletype)*npsph);
		long long i;
		for(i=0;i<DM_NP(simpar);i++){
			(SPH_BP(simpar)+ SPH_NP(simpar) +i)->x = (DM_BP(simpar)+i)->x;
			(SPH_BP(simpar)+ SPH_NP(simpar) +i)->y = (DM_BP(simpar)+i)->y;
			(SPH_BP(simpar)+ SPH_NP(simpar) +i)->z = (DM_BP(simpar)+i)->z;
		}
		Free(DM_BP(simpar));
		STAR_NP(simpar) = DM_NP(simpar) =0;
		SPH_NP(simpar) =npsph;
		for(i=0;i<SPH_NP(simpar);i++){
			(SPH_BP(simpar)+ i)->vx = (SPH_BP(simpar)+ i)->vy = (SPH_BP(simpar)+ i)->vz = 0.;
			(SPH_BP(simpar)+ i)->mu = RT_MU(simpar);
			y = (SPH_BP(simpar)+i)->y /(float) (NY(simpar)) - 0.5;
			if(fabs(y)<0.5){
				(SPH_BP(simpar)+i)->mass = 1;
				(SPH_BP(simpar)+i)->rho = 1.;
				SPH_BP(simpar)[i].temp = RT_INITTEMP(simpar);
			}
			else {
				(SPH_BP(simpar)+i)->mass = 2;
				(SPH_BP(simpar)+i)->rho = 2;
				SPH_BP(simpar)[i].temp = RT_INITTEMP(simpar)/2;
			}
			SPH_BP(simpar)[i].Entropy = kBovermH*SPH_BP(simpar)[i].temp/SPH_BP(simpar)[i].mu *
				pow(SPH_BP(simpar)[i].rho, 1 - simpar->sph.gamma);
		}
		SPH_TBP(simpar) = (treesphparticletype*)(SPH_BP(simpar));
	}
	printf("P%d has particle %g %g %g %g %g %g\n",simpar->myid,SPH_BP(simpar)->x,SPH_BP(simpar)->y,SPH_BP(simpar)->z,SPH_BP(simpar)->vx,SPH_BP(simpar)->vy,SPH_BP(simpar)->vz);
	RTEquilibrium(simpar);
}
*/


void  mkRTInit(SimParameters *simpar){

	long long i,j,k,np = 0;
	float A,lambda,x,y,z;


	DM_NP(simpar) = STAR_NP(simpar) = 0;

	SPH_NP(simpar) = (long long)NX(simpar)*(long long)NY(simpar)*(long long)LOCAL_NZ(simpar);

	SPH_BP(simpar) = (sphparticletype*)malloc(sizeof(sphparticletype)*SPH_NP(simpar));

	for(k=0;k<LOCAL_NZ(simpar);k++){
		for(j=0;j<NY(simpar);j++){
			for(i=0;i<NX(simpar);i++){
				x = i;
				y = j  +  sin(2*M_PI*x/NX(simpar)) * exp(-pow((j-0.5*NY(simpar))/2.,2.L));
				SPH_BP(simpar)[np].x = i;
				SPH_BP(simpar)[np].y = y;
				SPH_BP(simpar)[np].z = k+ LOCAL_Z_START(simpar);
				(SPH_BP(simpar)+ np)->vx = (SPH_BP(simpar)+ np)->vy = (SPH_BP(simpar)+ np)->vz = 0.;
				SPH_BP(simpar)[np].mu = RT_MU(simpar);
				SPH_BP(simpar)[np].temp = RT_INITTEMP(simpar);
				if(j<NY(simpar)/2) {
					SPH_BP(simpar)[np].mass = 2;
					SPH_BP(simpar)[np].temp = RT_INITTEMP(simpar)/2.;
				}
				else {
					SPH_BP(simpar)[np].mass = 1;
					SPH_BP(simpar)[np].temp = RT_INITTEMP(simpar);
				}
				CHANGEINDX(SPH_BP(simpar)+np, 1);
				SPH_BP(simpar)[np].Entropy = kBovermH*SPH_BP(simpar)[np].temp/SPH_BP(simpar)[np].mu *
					pow(SPH_BP(simpar)[np].rho, 1 - GAS_GAMMA(simpar));
				np++;
			}
		}
	}
	SPH_NP(simpar) = np;
	SPH_BP(simpar) = (sphparticletype*)realloc(SPH_BP(simpar),sizeof(sphparticletype)*SPH_NP(simpar));
	SPH_TBP(simpar) = (treesphparticletype*)(SPH_BP(simpar));
	RTEquilibrium(simpar);
}


void RTExternForce(SimParameters *simpar){
	long long i;
	for(i=0;i<SPH_NP(simpar);i++){
		SPH_TBP(simpar)[i].ay += -RT_FORCE(simpar)/NX(simpar);
	}
}
