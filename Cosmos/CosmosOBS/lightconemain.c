#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>

#include "pmheader.h"
#include "pm_common.h"

#define nobs 32

double xobs[nobs],yobs[nobs],zobs[nobs];
float maxd;


void Obs0ESLightConeData(treepmparticletype *,long,int,float,
		float,float,float,int,int,double *,double *,double *);
void BossESLightConeData(treepmparticletype *,long,int,float,
		float,float,float,int,int,double *,double *,double *, int, float );
void DiskESLightConeData(treepmparticletype *,long,int,float,
		float,float,float,int);
void Obs0SavingLightConeData(treepmparticletype *,long ,int ,float, 
		float , float ,float ,int , double * , double  *, double *);
void BossSavingLightConeData(treepmparticletype *,long ,int ,float, 
		float , float ,float ,int , double * , double  *, double *, int , float);
void DiskSavingLightConeData(treepmparticletype *,long ,int ,float, 
		float , float ,float);

void LightConeES(treepmparticletype *treeparticles, long np,
		int tnstep,float anow, float preastep, float nowastep,
		int flagsync){
	long i,j,k;
	int iobs;
	float xshift,yshift,zshift,shift;
	float redi;

	maxd = simpar.nspace * 0.1;


	iobs = 0;
	shift = simpar.nx/4.;
	for(k=0;k<4;k++){ 
		for(j=0;j<4;j++){ 
			xshift = ((j+k)%2)*shift;
			for(i=0;i<2;i++){
				xobs[iobs] = i*simpar.nx/2. + xshift;
				yobs[iobs] = j*simpar.ny/4.;
				zobs[iobs] = k*simpar.nz/4.;
				iobs ++;
			}
		}
	}
	/*
	iobs = 1;
	xobs[iobs] = simpar.nx/2;
	yobs[iobs] = simpar.ny/2;
	zobs[iobs] = simpar.nz/2;
	*/
	redi = 0.8;

	TIMER_START(62);
	BossESLightConeData(treeparticles,np,tnstep,anow,
			preastep,nowastep,maxd,flagsync,0,xobs,yobs, zobs,iobs,redi);
	TIMER_STOP(62);
	if(simpar.myid==0) fprintf(stdout,"BOSS ES Survey CPU= %f \n", ELAPSED_TIME(62));
	TIMER_START(62);
	DiskESLightConeData(treeparticles,np,tnstep,anow,
			preastep,nowastep,maxd,flagsync);
	TIMER_STOP(62);
	if(simpar.myid==0) fprintf(stdout,"DISK ES Survey CPU= %f \n", ELAPSED_TIME(62));
}
void LightConeS(treepmparticletype *treeparticles,long np, int tnstep,
		float amax, float anow,float preastep, float nowastep){
	float redi;
	int iobs = nobs;
	redi = 0.8;
	TIMER_START(62);
	BossSavingLightConeData(treeparticles,np,tnstep,amax,anow,
			preastep,nowastep,0,xobs,yobs, zobs,iobs,redi);
	TIMER_STOP(62);
	if(simpar.myid==0) fprintf(stdout,"BOSS S Survey CPU= %f \n", ELAPSED_TIME(62));
	TIMER_START(62);
	DiskSavingLightConeData(treeparticles,np,tnstep,amax,anow,
			preastep,nowastep);
	TIMER_STOP(62);
	if(simpar.myid==0) fprintf(stdout,"DISK S Survey CPU= %f \n", ELAPSED_TIME(62));

}
