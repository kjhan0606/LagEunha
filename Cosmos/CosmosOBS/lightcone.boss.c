/*
 * 이 코드는 임의의 X0,Y0,Z0 에 위치한 관측자가 시뮬레이션 박스를 기반으로
 * 시뮬레이션 입자들을 관측했을때  입자의 분포를 나타낸다.
 * 최총 결과는 comoving space 에서의 시뮬레이션에서 사용한 거리단위로 나온다.
 * 현재는 theta, phi  의 값이 1사분면에서만 계산이 가능하다.
 * 만약 2,3,4 사분면에서도 가능하게 할려면, min,max 의 값들을 적절히
 * 바꾸어야 한다. 그리고 마지막 부분에 해당되 입자를 찾을 때 phi 값을 
 * 계산하는데,
 * 여기에서는 그냥 1,4사분면이라고 가정을 했다. 
 * wflag 의 값은 2진법으로 표현해야 한다.
 * 09/08/2002 김주한 
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include "Memory.h"
#include "pmheader.h"
#include "lightcone.boss.h"
#define MIN(A,B) ((A) <(B) ? (A):(B))
#define PI 3.1415926535L
double Theta,Phi,DTheta,DPhi;
MPI_Status status;
char surveyfilename[100];
int *freework;
typedef struct Spherical{
	double r,theta,phi;
} Spherical;
double r1, r2;
double theta,phi,theta1,phi1,theta2;
float omega0,H,lambda0;
float qsimp(float (*)(float),float,float);
float trapzd(float (*)(float), float , float , int );
int ssorttype(const void *, const void *);
int bosssorttype(const void *a,const void *b){
	bossparticletype *aa,*bb;
	aa = (bossparticletype *)a;
	bb = (bossparticletype *)b;
	if(aa->iobs < bb->iobs) return -1;
	else if(aa->iobs > bb->iobs) return +1;
	else {
		if(aa->type < bb->type) return -1;
		else if(aa->type > bb->type) return +1;
		return 0;
	}
}
#define NPSTEP 400000
int nnp[3];
/* Particle data type for 108 BOSS survey */
static double rmin2,rmax2;
static bossparticletype *bossparticles;

static int snlcp,maxlcp,tsnlcp,mbox;
/* if dist >= r1-maxd && dist < r1 --> located in inner buffer zone 
   if dist >= r1 && dist < r2 --> located in the main zone
   if dist >= r2 && dist < r2+maxd --> located in the outter buffer zone
   */
#define twoPI (6.283185306L)
#define deg2rad (0.017453292516L)


double SX0,SY0,SZ0;
void a2comovingpixel(double *, double *, double *,int , float ,float , float ,float , 
		float ,float ,float , float );
void BossESLightConeData(treepmparticletype *treeparticles,long np,int tnstep,float a,
		float preastep, float astep,float maxd,int flagsync,int obsidstart,
		double *sx0, double *sy0, double *sz0,int nobs, float redi){
	int nx,ny,nz, nspace;
	float amax;
	float boxsize,hboxsize;
	float omep,omeplam,hubble,omepb;
	int myid,nid;
	treepmparticletype *bp;
	float redshift;
	double r,dr1,dr2,distsq,dist,rot;
	double xp,yp,zp;
	double xpp,ypp,zpp;
	double rminsq,rmaxsq;
	int nbox,ibox;
	long i,j,k;
	long ix,iy,iz;
	lightconetype *lightcones;

	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;
	nspace = simpar.nspace;
	boxsize = simpar.boxsize;
	hboxsize = boxsize*0.5;
	omep = simpar.omep;
	omeplam = simpar.omeplam;
	omepb = simpar.omepb;
	hubble = simpar.hubble;
	myid = simpar.myid;
	nid = simpar.nid;
	amax = simpar.amax;

	if(1){
		float redshift;
		redshift = amax/a-1.;
		if(redshift > redi) return;
	}
	/*
	nbox = 3*3*2;
	box = (InBox *)Malloc(sizeof(InBox)*nbox,PPTR(box));
	theta1 = 48.L*deg2rad;
	theta2 = 132.L*deg2rad;
	phi1 = 65.L*deg2rad;
	rot = 0.L;
	*/
	a2comovingpixel(&r,&dr1,&dr2,nx,amax,a,preastep,astep,
			omep,omeplam,boxsize,hubble);
	r1 = (r-dr1);
	if(r1 <0.) r1 = 0;
	r2 = (r+dr2);
	rmin2 = r1*r1;
	rmax2 = r2*r2;
	rminsq = (r1-maxd)*(r1-maxd);
	if(r1<maxd) rminsq = 0.;
	rmaxsq = (r2+maxd)*(r2+maxd);
	if(myid==0) printf("Spherical Boss Survey %g %g %g\n",r,dr1,dr2);
	snlcp = 0;
	maxlcp = NPSTEP;
	bossparticles = (bossparticletype *)Malloc(sizeof(bossparticletype)*maxlcp, PPTR(bossparticles));
	{
		bp = treeparticles;
		for(i=0;i<np;i++){
			for(j=0;j<nobs;j++){
				SX0 = sx0[j]; SY0 = sy0[j]; SZ0 = sz0[j];
				xp = XofP(bp) - SX0;
				yp = YofP(bp) - SY0;
				zp = ZofP(bp) - SZ0;
				/* This should be dealt with much care. */
				BOSSBCondition(xp,hboxsize,boxsize);
				BOSSBCondition(yp,hboxsize,boxsize);
				BOSSBCondition(zp,hboxsize,boxsize);

				distsq = xp*xp+yp*yp+zp*zp;
				if(distsq >= rminsq && distsq < rmaxsq){
					dist = sqrt(distsq);
					bossparticles[snlcp].vx = bp->vx;
					bossparticles[snlcp].vy = bp->vy;
					bossparticles[snlcp].vz = bp->vz;
					bossparticles[snlcp].bp = bp;
					if(dist>=r1&& dist < r2) bossparticles[snlcp].type = 0;
					else if(dist >=r2 && dist < r2+maxd) bossparticles[snlcp].type = 1;
					else bossparticles[snlcp].type = 2;
					bossparticles[snlcp].iobs = j;
					snlcp ++; 
					if(snlcp >= maxlcp){
						maxlcp += NPSTEP;
						bossparticles = (bossparticletype*)Realloc(bossparticles,
								sizeof(bossparticletype)*maxlcp);
					}
				}
			}
			bp++;
		}
	}
	if(snlcp >0) bossparticles = (bossparticletype*)Realloc(bossparticles, sizeof(bossparticletype)*snlcp);
	else 
		Free(bossparticles);

	if(flagsync){
		int ioffset=0;
		if(snlcp>0) qsort(bossparticles,snlcp,sizeof(bossparticletype),bosssorttype);
		for(j=0;j<nobs;j++){
			int npiobs=0;

			SX0 = sx0[j]; SY0 = sy0[j]; SZ0 = sz0[j];

			nnp[0] = nnp[1] = nnp[2] = 0;
			for(i=ioffset;i<snlcp;i++){
				if(bossparticles[i].iobs == j) npiobs ++;
			}
			for(i=ioffset;i<ioffset+npiobs;i++){
				nnp[bossparticles[i].type]++;
			}
			if(npiobs > 0) lightcones = (lightconetype *)
				Malloc(sizeof(lightconetype)*npiobs,PPTR(lightcones));
			k=ioffset;
			for(i=0;i<npiobs;i++,k++){

				bp = bossparticles[k].bp;
				xp = XofP(bp) - SX0;
				yp = YofP(bp) - SY0;
				zp = ZofP(bp) - SZ0;
				/* This should be dealt with much care. */
				BOSSBCondition(xp,hboxsize,boxsize);
				BOSSBCondition(yp,hboxsize,boxsize);
				BOSSBCondition(zp,hboxsize,boxsize);

				lightcones[i].x = xp;
				lightcones[i].y = yp;
				lightcones[i].z = zp;
				lightcones[i].vx = bossparticles[k].vx;
				lightcones[i].vy = bossparticles[k].vy;
				lightcones[i].vz = bossparticles[k].vz;
#ifdef INDEX
				lightcones[i].indx = bp->indx;
#endif
			}
			{
				lightconetype *pp;
		
				sprintf(surveyfilename,"Bosslightcone.%.3d.main.%.5d",(int)j+obsidstart,tnstep);
				pp = lightcones;
				writedownlightconedata(surveyfilename,pp,nnp[0],sx0[j],sy0[j],sz0[j]);
		
				sprintf(surveyfilename,"Bosslightcone.%.3d.outer.%.5d",(int)j+obsidstart,tnstep);
				pp += nnp[0];
				writedownlightconedata(surveyfilename,pp,nnp[1],sx0[j],sy0[j],sz0[j]);
		
				sprintf(surveyfilename,"Bosslightcone.%.3d.inner.%.5d",(int)j+obsidstart,tnstep);
				pp += nnp[1];
				writedownlightconedata(surveyfilename,pp,nnp[2],sx0[j],sy0[j],sz0[j]);
			}
			if(npiobs > 0 ) Free(lightcones);
			ioffset += npiobs;
		}
	}
	else {
		OBSTEMWRITE(bossparticles,snlcp,myid,nid,obsidstart);
	}
	if(snlcp >0) Free(bossparticles);
	if(myid==0) printf("Boss SURVEY %d detected\n",snlcp);
}
#include "indT.h"
void BossSavingLightConeData(treepmparticletype *treeparticles,long np,int tnstep, 
		float amax, float a, float preastep,float astep,
		int obsidstart, double *sx0, double *sy0, double *sz0,int nobs,float redi){
	treepmparticletype *bp;
	int nx,ny,nz;
	float redshift;
	double r,dr1,dr2;
#ifdef XYZDBL
	double xp,yp,zp,distsq,dist;
#else
	float xp,yp,zp,distsq,dist;
#endif
	float rminsq,rmaxsq;
	long i,j,k;
	float ax,ay,az;
	bossparticletype *tp;
	lightconetype *lightcones;
	float ainv,app,apsq,afact,bfact,vfact1h,vfact2h,rng;
	float asteph;
	float vfact1,vfact2;
	float gv1,gv2;
	float omei,omep,omepb,omeplam,hubble;
	int myid,nid;
	char nowtsubpower;
	float boxsize,hboxsize;

	boxsize = simpar.boxsize;
	hboxsize = simpar.boxsize*0.5;
	omei = simpar.omei;
	omep = simpar.omep;
	omepb = simpar.omepb;
	omeplam = simpar.omeplam;
	hubble = simpar.hubble;
	myid = simpar.myid;
	nid = simpar.nid;
	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;

	if(1){
		float redshift;
		redshift = amax/a-1.;
		if(redshift > redi) return;
	}

	OBSTEMREAD(bossparticles,snlcp,myid,nid,obsidstart);

	MPI_Reduce(&snlcp,&tsnlcp,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&tsnlcp,1,MPI_INT,0,MPI_COMM_WORLD);
	if(tsnlcp==0) return;



	/* recover the synchronized velocity */
	tp = bossparticles;
	for(i=0;i<snlcp;i++){
		bp = tp->bp;
		nowtsubpower = GetTsubPower(bp);
		if(isnowdmstep(nowTsubdiv,nowtsubpower,maxTsubpower)){
			float vfact,vfact2,vfact3;
			vfact2 = EvolFactKick[nowtsubpower].fact2;
			vfact3 = EvolFactPull[nowtsubpower].fact2;
			vfact = vfact3/vfact2;
			ax = (bp->vx-tp->vx)*vfact;
			ay = (bp->vy-tp->vy)*vfact;
			az = (bp->vz-tp->vz)*vfact;
			tp->vx += ax;
			tp->vy += ay;
			tp->vz += az;
		}
		tp++;
	}
	{
		int ioffset=0;
		if(snlcp>0) qsort(bossparticles,snlcp,sizeof(bossparticletype),bosssorttype);
		for(j=0;j<nobs;j++){
			int npiobs = 0;
			nnp[0] = nnp[1] = nnp[2] = 0;

			SX0 = sx0[j]; SY0 = sy0[j]; SZ0 = sz0[j];

			for(i=ioffset;i<snlcp;i++){
				if(bossparticles[i].iobs == j) npiobs ++;
			}
			for(i=ioffset;i<ioffset+npiobs;i++){
				nnp[bossparticles[i].type]++;
			}
			if(npiobs > 0) lightcones = (lightconetype *) 
				Malloc(sizeof(lightconetype)*npiobs,PPTR(lightcones));
			k=ioffset;
			for(i=0;i<npiobs;i++,k++){

				bp = bossparticles[k].bp;
				xp = XofP(bp) - SX0;
				yp = YofP(bp) - SY0;
				zp = ZofP(bp) - SZ0;
				/* This should be dealt with much care. */
				BOSSBCondition(xp,hboxsize,boxsize);
				BOSSBCondition(yp,hboxsize,boxsize);
				BOSSBCondition(zp,hboxsize,boxsize);

				lightcones[i].x = xp;
				lightcones[i].y = yp;
				lightcones[i].z = zp;
				lightcones[i].vx = bossparticles[k].vx;
				lightcones[i].vy = bossparticles[k].vy;
				lightcones[i].vz = bossparticles[k].vz;
#ifdef INDEX
				lightcones[i].indx = bp->indx;
#endif
			}
			{
				lightconetype *pp;
		
				sprintf(surveyfilename,"Bosslightcone.%.3d.main.%.5d",(int)j+obsidstart,tnstep);
				pp = lightcones;
				writedownlightconedata(surveyfilename,pp,nnp[0],sx0[j],sy0[j],sz0[j]);
		
				sprintf(surveyfilename,"Bosslightcone.%.3d.outer.%.5d",(int)j+obsidstart,tnstep);
				pp += nnp[0];
				writedownlightconedata(surveyfilename,pp,nnp[1],sx0[j],sy0[j],sz0[j]);
		
				sprintf(surveyfilename,"Bosslightcone.%.3d.inner.%.5d",(int)j+obsidstart,tnstep);
				pp += nnp[1];
				writedownlightconedata(surveyfilename,pp,nnp[2],sx0[j],sy0[j],sz0[j]);
			}
			if(npiobs > 0 ) Free(lightcones);
			ioffset += npiobs;
		}
		Free(bossparticles);
	}
}
#undef npstep
