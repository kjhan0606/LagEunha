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
#include <math.h>
#include <mpi.h>
#include "Memory.h"
#include "pmheader.h"
#include "indT.h"
#include "lightcone.boss.h"
#define MIN(A,B) ((A) <(B) ? (A):(B))
#define PI (3.14159265358979323846L)
MPI_Status status;
int nnp[3];
double Theta,Phi,DTheta,DPhi;
char surveyfilename[100];
int *freework;
typedef struct Spherical{
	double r,theta,phi;
} Spherical;
/*
typedef struct InBox{
	int nx,ny,nz;
} InBox;
*/
InBox *box;
double r1, r2;
double theta,phi,theta1,phi1,theta2,phi2;
float omega0,H,lambda0;
float qsimp(float (*)(float),float,float);
float trapzd(float (*)(float), float , float , int );
/* ov: velocity data at half-step before the current step */
/* nv: velocity data at half-step after the current step */

int pslicesorttype(const void *a,const void *b){
	slcparticletype *aa,*bb;
	aa = (slcparticletype *)a;
	bb = (slcparticletype *)b;
	if(aa->type < bb->type) return -1;
	else if(aa->type > bb-> type) return +1;
	else return 0;
}

#define NPSTEP 10000

static double rmin2,rmax2;
static slcparticletype *slcparticles;
static int snlcp,maxlcp,tsnlcp,mbox;
/* if dist >= r1-maxd && dist < r1 --> located in inner buffer zone 
   if dist >= r1 && dist < r2 --> located in the main zone
   if dist >= r2 && dist < r2+maxd --> located in the outter buffer zone
   */
#define ZWIDTH (20.)
float SX0,SY0,SZ0;
#define twoPI (6.28318530717958647692L)
#define deg2rad (.01745329251994329576L)
void a2comovingpixel(double *, double *, double *,int , float ,float , float ,float ,
		float ,float ,float , float );
void DiskESLightConeData(treepmparticletype *treeparticles,long np,
		int tnstep, float anow,float preastep, float astep,float maxd, int flagsync){
	int nx,ny,nz, nspace;
	float  amax, a,  boxsize;
	float omep,omeplam,hubble,omepb;
	int myid,nid;
	treepmparticletype *bp;
	float redshift;
	double r,dr1,dr2,distsq,dist,rot;
	double xp,yp,zp;
	double xpp,ypp,zpp;
	double rminsq,rmaxsq;
	double xrot,yrot,zrot,phimin,zwidth,mzwidth,xppp,yppp,zppp;
	int nbox,ibox;
	long i,j,k;
	long ix,iy,iz,nxy;

	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;
	nspace = simpar.nspace;
	amax = simpar.amax;
	a = anow;
	boxsize = simpar.boxsize;
	omep = simpar.omep;
	omeplam = simpar.omeplam;
	omepb = simpar.omeplam;
	hubble = simpar.hubble;
	myid = simpar.myid;
	nid = simpar.nid;


	SX0 = nx/2;
	SY0 = ny/2+40.4;
	SZ0 = nz/2;

	nbox = 9;
	box = (InBox *)Malloc(sizeof(InBox)*nbox,PPTR(box));
	/*
	theta1 = 75.L*deg2rad;
	theta2 = 105.L*deg2rad;
	phi1 = 30.L*deg2rad;
	phi2 = 330.L*deg2rad;
	*/

	/*
	zrot = 45.L*deg2rad;
	yrot = acos(sqrt(2.)/sqrt(3.));
	*/
	xrot = atan(0.1L);
	yrot = atan(0.23333L);
	phimin = 33.L*deg2rad;
	/*
	nbox = 0;
	for(k=-1;k<=0;k++){
		box[ibox].nx = 0.;
		box[ibox].ny = 0.;
		box[ibox].nz = k*nz-SZ0;
		nbox ++;
	}
	*/
	for(j=0;j<3;j++) for(i=0;i<3;i++){
		box[i].nx = nx*(i-1)-SX0;
		box[i].ny = ny*(j-1)-SY0;
		box[i].nz =         -SZ0;
	}

	zwidth = ZWIDTH;
	mzwidth = -1.*zwidth;

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
	if(myid==0) printf("Disk SURVEY %g %g %g at observer=%g %g %g ....",r,dr1,dr2,
			SX0,SY0,SZ0);
	snlcp = 0;
	maxlcp = NPSTEP;
	slcparticles = (slcparticletype *)Malloc(sizeof(slcparticletype)*maxlcp,
			PPTR(slcparticles));
	for(j=0;j<nbox;j++){
		bp = treeparticles;
		for(i=0;i<np;i++){
			xp = (double)XofP(bp) + (double)box[j].nx;
			yp = (double)YofP(bp) + (double)box[j].ny;
			zp = (double)ZofP(bp) + (double)box[j].nz;
			distsq = xp*xp + yp*yp + zp*zp;
			if(distsq >= rminsq && distsq < rmaxsq){
				xpp =  xp;
				ypp =  cos(xrot)*yp + sin(xrot)*zp;
				zpp = -sin(xrot)*yp + cos(xrot)*zp; 
	
	
				yppp = ypp;
				zppp =  cos(yrot)*zpp + sin(yrot)*xpp;
				xppp = -sin(yrot)*zpp + cos(yrot)*xpp;
	
				if(zppp >=mzwidth && zppp < zwidth){
					/*
					double distrho;
					distrho = sqrt(xppp*xppp+yppp*yppp);
					phi = acos(xppp/distrho);
					if(phi <= phimin)
					*/
					{
						dist = sqrt(distsq);
						slcparticles[snlcp].x = xp;
						slcparticles[snlcp].y = yp;
						slcparticles[snlcp].z = zp;
						/* Write down the velocity data at half-step before 
						 * identifying particles */
						slcparticles[snlcp].vx = bp->vx;
						slcparticles[snlcp].vy = bp->vy;
						slcparticles[snlcp].vz = bp->vz;
						slcparticles[snlcp].bp = bp;
						if(dist>=r1&& dist < r2) slcparticles[snlcp].type = 0;
						else if(dist >=r2 && dist < r2+maxd) slcparticles[snlcp].type = 1;
						else slcparticles[snlcp].type = 2;
						snlcp ++;
						if(snlcp >= maxlcp){
							maxlcp += NPSTEP;
							slcparticles = (slcparticletype*)Realloc(slcparticles,
									sizeof(slcparticletype)*maxlcp);
	
						}
					}
				}
			}
			bp++;
		}
	}
	if(snlcp>0) slcparticles = (slcparticletype*)Realloc(slcparticles,
		sizeof(slcparticletype)*snlcp);
	else
		Free(slcparticles);
	if(flagsync){
		lightconetype *lightcones;
		int ioffset = 0;
		if(snlcp >0) qsort(slcparticles,snlcp,sizeof(slcparticletype),pslicesorttype);
		if(snlcp>0) lightcones = (lightconetype *)Malloc(sizeof(lightconetype)*snlcp,PPTR(lightcones));
		nnp[0] = nnp[1] = nnp[2] = 0;
		for(i=0;i<snlcp;i++) nnp[slcparticles[i].type]++;
		for(i=0;i<snlcp;i++){
			bp = slcparticles[i].bp;
			lightcones[i].x = slcparticles[i].x;
			lightcones[i].y = slcparticles[i].y;
			lightcones[i].z = slcparticles[i].z;
			lightcones[i].vx = slcparticles[i].vx;
			lightcones[i].vy = slcparticles[i].vy;
			lightcones[i].vz = slcparticles[i].vz;
#ifdef INDEX
			lightcones[i].indx = bp->indx;
#endif
		}
		{
			lightconetype *pp;

			sprintf(surveyfilename,"Disklightcone.main.%.5d",tnstep);
			pp = lightcones;
			writedownlightconedata(surveyfilename,pp,nnp[0],SX0,SY0,SZ0);
	
			sprintf(surveyfilename,"Disklightcone.outer.%.5d",tnstep);
			pp += nnp[0];
			writedownlightconedata(surveyfilename,pp,nnp[1],SX0,SY0,SZ0);
	
			sprintf(surveyfilename,"Disklightcone.inner.%.5d",tnstep);
			pp += nnp[1];
			writedownlightconedata(surveyfilename,pp,nnp[2],SX0,SY0,SZ0);
		}
		if(snlcp>0) Free(lightcones);
	}
	else { 
		SLCP_OBSTEMWRITE(slcparticles,snlcp,myid,nid,1000);
	}
	if(snlcp >0) Free(slcparticles);
	Free(box);
	if(myid==0) printf("Disk SURVEY %d detected\n",snlcp);
}


void DiskSavingLightConeData(treepmparticletype *treeparticles,long np,
		int tnstep, float amax, float a, float astep,float nextastep){
	float boxsize, hboxsize;
	float omei;
	float omep,omepb,omeplam,hubble;
	int myid,nid;
	int nx, ny, nz;
	treepmparticletype *bp;
	float redshift;
	double r,dr1,dr2;
	double xp,yp,zp,distsq,dist;
	float rminsq,rmaxsq;
	int nbox;
	long i,j,k;
	float ax,ay,az;
	slcparticletype *tp;
	lightconetype *lightcones;
	float ainv,app,apsq,afact,bfact,vfact1h,vfact2h,rng;
	float asteph;
	float vfact1,vfact2;
	float gv1,gv2;
	char nowtsubpower;
	MPI_Status status;

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


	SLCP_OBSTEMREAD(slcparticles,snlcp,myid,nid,1000);
	MPI_Reduce(&snlcp,&tsnlcp,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&tsnlcp,1,MPI_INT,0,MPI_COMM_WORLD);
	if(tsnlcp==0) return;


	tp = slcparticles;
	for(i=0;i<snlcp;i++){
		bp = tp->bp;
		nowtsubpower = GetTsubPower(bp);
		if(isnowdmstep(nowTsubdiv,nowtsubpower,maxTsubpower)){
			float vfact,vfact2,vfact3;
			float ax,ay,az;
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
		tp ++;
	}



	nnp[0] = nnp[1] = nnp[2] =0;
	for(i=0;i<snlcp;i++){
		nnp[slcparticles[i].type]++;
	}
	if(snlcp>0) qsort(slcparticles,snlcp,sizeof(slcparticletype),pslicesorttype);
	if(snlcp>0) lightcones = (lightconetype *)Malloc(sizeof(lightconetype)*snlcp,
			PPTR(lightcones));
	for(i=0;i<snlcp;i++){
		bp = slcparticles[i].bp;
		lightcones[i].x = slcparticles[i].x;
		lightcones[i].y = slcparticles[i].y;
		lightcones[i].z = slcparticles[i].z;
		lightcones[i].vx = slcparticles[i].vx;
		lightcones[i].vy = slcparticles[i].vy;
		lightcones[i].vz = slcparticles[i].vz;
#ifdef INDEX
		lightcones[i].indx = bp->indx;
#endif
	}
	{
		lightconetype *pp;

		sprintf(surveyfilename,"Disklightcone.main.%.5d",tnstep);
		pp = lightcones;
		writedownlightconedata(surveyfilename,pp,nnp[0],SX0,SY0,SZ0);

		sprintf(surveyfilename,"Disklightcone.outer.%.5d",tnstep);
		pp += nnp[0];
		writedownlightconedata(surveyfilename,pp,nnp[1],SX0,SY0,SZ0);

		sprintf(surveyfilename,"Disklightcone.inner.%.5d",tnstep);
		pp += nnp[1];
		writedownlightconedata(surveyfilename,pp,nnp[2],SX0,SY0,SZ0);
	}
	if(snlcp>0) Free(lightcones);
	Free(slcparticles);
}

#undef npstep
