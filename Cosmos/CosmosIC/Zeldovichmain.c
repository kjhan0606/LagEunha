#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "cosmosic.h"
//#include "mpirms.h"
#include "zeld.h"

#define MAX(a,b) ( (a)>(b)?(a):(b) )

	



void AllParticleMigrate(SimParameters *);
/*
void AllParticleMigrate(SimParameters *simpar){
	pmigrate((void**)&DM_BP(simpar),  &DM_NP(simpar), DM_DDINFO(simpar));
	pmigrate((void**)&SPH_BP(simpar),  &SPH_NP(simpar), SPH_DDINFO(simpar));
	pmigrate((void**)&STAR_BP(simpar),  &STAR_NP(simpar), STAR_DDINFO(simpar));
	pmigrate((void**)&AGN_BP(simpar),  &AGN_NP(simpar), AGN_DDINFO(simpar));
}
*/

void Zeldovichmain(SimParameters *simpar){
	size_t i,j,k;
	enum mtype itype;


	for(itype=cdm;itype<star;itype++){
		if(itype ==cdm && DM_TNP(simpar) ==0) continue;
		else if(itype ==sph && SPH_TNP(simpar) ==0) continue;

		size_t densize = 2*LOCAL_GRID_SIZE(simpar);
		size_t fftwsize = 2*(NX(simpar)/2+1) * NY(simpar) * (LOCAL_NZ(simpar) + 5);
		densize = MAX(10, densize);
		DenType *d1 = (DenType*)malloc(sizeof(DenType)*densize);
		DenType *fx = (DenType*)malloc(sizeof(DenType)*densize);
		DenType *fy = (DenType*)malloc(sizeof(DenType)*densize);
		DenType *fz = (DenType*)malloc(sizeof(DenType)*densize);
		DenType *denpad = (DenType*)malloc(sizeof(DenType)*2*NX(simpar)*NY(simpar));

		enum mtype iitype = itype;

		if(GAS_TYPE(simpar) == 'N') iitype = cdm;

		int nx = NX(simpar);
		int ny = NY(simpar);
		int local_nz = LOCAL_NZ(simpar);
		int local_z_start = LOCAL_Z_START(simpar);

//		void PMSeedForce(SimParameters *, DenType *, DenType *, DenType *, DenType *, DenType *, DenType *, int );
		PMSeedForce(simpar, d1,fx, fy, fz, NULL, denpad, itype);
		free(denpad);
		free(d1);

		float pamp, vamp;
		pamp = COS_PAMP(simpar, itype);
		vamp = COS_VAMP1(simpar, itype);

		GridInfo *ffx,*ffy, *ffz;


		if(itype ==cdm){
			Zeld_GetVR(simpar,ffx, fx, vamp,dm); free(fx);
			Zeld_GetVR(simpar,ffy, fy, vamp,dm); free(fy);
			Zeld_GetVR(simpar,ffz, fz, vamp,dm); free(fz);
		}
		else if(itype ==sph){
			Zeld_GetVR(simpar,ffx, fx, vamp,sph); free(fx);
			Zeld_GetVR(simpar,ffy, fy, vamp,sph); free(fy);
			Zeld_GetVR(simpar,ffz, fz, vamp,sph); free(fz);
		}

		free(fx);free(fy);free(fz);
		Set_Particle_Zeld(simpar, ffx,ffy,ffz, itype);
		if(itype==cdm) pmigrate((void**)(&DM_BP(simpar)), &DM_NP(simpar), DM_DDINFO(simpar), &GRIDINFO(simpar));
		else if(itype==sph) pmigrate((void**)(&SPH_BP(simpar)), &SPH_NP(simpar), SPH_DDINFO(simpar), &GRIDINFO(simpar));
	}
}

