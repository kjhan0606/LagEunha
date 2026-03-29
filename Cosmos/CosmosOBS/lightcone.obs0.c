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
#include "lightcone.h"
#define MIN(A,B) ((A) <(B) ? (A):(B))
#define PI 3.1415926535L
double Theta,Phi,DTheta,DPhi;
MPI_Status status;
char surveyfilename[100];
int *freework;
typedef struct Spherical{
	double r,theta,phi;
} Spherical;
InBox *box;
double r1, r2;
double theta,phi,theta1,phi1,theta2;
float omega0,H,lambda0;
float qsimp(float (*)(float),float,float);
float trapzd(float (*)(float), float , float , int );
int ssorttype(const void *, const void *);
#define NPSTEP 10000
int nnp[3];

static double rmin2,rmax2;
static slcparticletype *slcparticles;

static int snlcp,maxlcp,tsnlcp,mbox;
/* if dist >= r1-maxd && dist < r1 --> located in inner buffer zone 
   if dist >= r1 && dist < r2 --> located in the main zone
   if dist >= r2 && dist < r2+maxd --> located in the outter buffer zone
   */
/*
#define SX0 (nx/4.)
#define SY0 (ny/4.)
#define SZ0 (nz/4.)
*/
#define twoPI (6.283185306L)
#define deg2rad (0.017453292516L)
void a2comovingpixel(double *, double *, double *,int , float ,float , float ,float , 
		float ,float ,float , float );
void Obs0ESLightConeData(treepmparticletype *treeparticles,int np,int mstep,float a,
		float beforeastep, float astep,float maxd,int saveflag,int obsid,
		double SX0, double SY0, double SZ0){
	int nx,ny,nz, nspace;
	int amax, boxsize;
	float omep,omeplam,hubble,omepb;
	int myid,nid;
	treepmparticletype *bp;
	float redshift;
	double r,dr1,dr2,distsq,dist,rot;
	double xp,yp,zp;
	double xpp,ypp,zpp;
	double rminsq,rmaxsq;
	int nbox,ibox;
	int i,j,k;
	long ix,iy,iz,nxy;
	lightconetype *lightcones;

	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;
	nspace = simpar.nspace;
	boxsize = simpar.boxsize;
	omep = simpar.omep;
	omeplam = simpar.omeplam;
	omepb = simpar.omepb;
	hubble = simpar.hubble;
	myid = simpar.myid;
	nid = simpar.nid;
	amax = simpar.amax;

	{
		float redshift;
		redshift = amax/a-1.;
		if(redshift > 1.0) return;
	}
	nxy = nx*ny;
	/*
	nbox = 3*3*2;
	box = (InBox *)Malloc(sizeof(InBox)*nbox,PPTR(box));
	theta1 = 48.L*deg2rad;
	theta2 = 132.L*deg2rad;
	phi1 = 65.L*deg2rad;
	rot = 0.L;
	*/
	a2comovingpixel(&r,&dr1,&dr2,nx,amax,a,beforeastep,astep,
			omep,omeplam,boxsize,hubble);
	r1 = (r-dr1);
	if(r1 <0.) r1 = 0;
	r2 = (r+dr2);
	rmin2 = r1*r1;
	rmax2 = r2*r2;
	rminsq = (r1-maxd)*(r1-maxd);
	if(r1<maxd) rminsq = 0.;
	rmaxsq = (r2+maxd)*(r2+maxd);
	if(myid==0) printf("Slice SURVEY %g %g %g at observer=%g %g %g ....",r,dr1,dr2,SX0,SY0,SZ0);
	snlcp = 0;
	maxlcp = NPSTEP;
	slcparticles = (slcparticletype *)Malloc(sizeof(slcparticletype)*maxlcp, PPTR(slcparticles));
	for(j=0;j<1;j++){
		bp = treeparticles;
		for(i=0;i<np;i++){
			xp = XofP(bp) - SX0;
			yp = YofP(bp) - SY0;
			zp = ZofP(bp) - SZ0;
			distsq = xp*xp+yp*yp+zp*zp;
			if(distsq >= rminsq && distsq < rmaxsq){
				dist = sqrt(distsq);
				slcparticles[snlcp].x = xp;
				slcparticles[snlcp].y = yp;
				slcparticles[snlcp].z = zp;
				slcparticles[snlcp].vx = bp->vx;
				slcparticles[snlcp].vy = bp->vy;
				slcparticles[snlcp].vz = bp->vz;
				slcparticles[snlcp].bp = bp;
				if(dist>=r1&& dist < r2) slcparticles[snlcp].type = 0;
				else if(dist >=r2 && dist < r2+maxd) slcparticles[snlcp].type = 1;
				else slcparticles[snlcp].type = 2;
#ifdef INDEX
				slcparticles[snlcp].indx = bp->indx;
#endif
				snlcp ++;
				if(snlcp >= maxlcp){
					maxlcp += NPSTEP;
					slcparticles = (slcparticletype*)Realloc(slcparticles,
							sizeof(slcparticletype)*maxlcp);
				}
			}
			bp++;
		}
	}
	slcparticles = (slcparticletype*)Realloc(slcparticles, sizeof(slcparticletype)*snlcp);
	if(saveflag){
		nnp[0] = nnp[1] = nnp[2] = 0;
		for(i=0;i<snlcp;i++){
			nnp[slcparticles[i].type]++;
		}
		if(snlcp>0) qsort(slcparticles,snlcp,sizeof(slcparticletype),ssorttype);
		if(snlcp>0) lightcones = (lightconetype *)Malloc(sizeof(lightconetype)*snlcp, PPTR(lightcones));
		for(i=0;i<snlcp;i++){
			lightcones[i].x = slcparticles[i].x;
			lightcones[i].y = slcparticles[i].y;
			lightcones[i].z = slcparticles[i].z;
			lightcones[i].vx = slcparticles[i].vx;
			lightcones[i].vy = slcparticles[i].vy;
			lightcones[i].vz = slcparticles[i].vz;
#ifdef INDEX
			lightcones[i].indx = slcparticles[i].indx;
#endif
		}
		{
			lightconetype *pp;
	
			sprintf(surveyfilename,"Olightcone0.main.%.5d",mstep);
			pp = lightcones;
			writedownlightconedata(surveyfilename,pp,nnp[0]);
	
			sprintf(surveyfilename,"Olightcone0.outer.%.5d",mstep);
			pp += nnp[0];
			writedownlightconedata(surveyfilename,pp,nnp[1]);
	
			sprintf(surveyfilename,"Olightcone0.inner.%.5d",mstep);
			pp += nnp[1];
			writedownlightconedata(surveyfilename,pp,nnp[2]);
		}
	
		if(snlcp > 0) Free(lightcones);
		Free(slcparticles);
	}
	else {
		OBSTEMWRITE(slcparticles,snlcp,myid,nid,obsid);
		Free(slcparticles);
	}

	if(myid==0) printf("OSlice SURVEY %d detected\n",snlcp);
}
#include "indT.h"
void Obs0SavingLightConeData(treepmparticletype *treeparticles,int np,int mstep, 
		float amax, float a, float beforeastep,float astep,
		int obsid, double SX0, double SY0, double SZ0){
	treepmparticletype *bp;
	int nx,ny,nz;
	float redshift;
	double r,dr1,dr2;
	float xp,yp,zp,distsq,dist;
	float rminsq,rmaxsq;
	int i,j,k;
	float ax,ay,az;
	slcparticletype *tp;
	lightconetype *lightcones;
	float ainv,app,apsq,afact,bfact,vfact1h,vfact2h,rng;
	float asteph;
	float vfact1,vfact2;
	float gv1,gv2;
	float boxsize,omei,omep,omepb,omeplam,hubble;
	int myid,nid;
	boxsize = simpar.boxsize;
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

	{
		float redshift;
		redshift = amax/a-1.;
		if(redshift > 1.0) return;
	}

	OBSTEMREAD(slcparticles,snlcp,myid,nid,obsid);

	MPI_Reduce(&snlcp,&tsnlcp,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(&tsnlcp,1,MPI_INT,0,MPI_COMM_WORLD);
	if(tsnlcp==0) return;


	/*
	asteph = astep*0.5;
	rng = nx;
	ainv = 1./a;
	app  = -4.*PI/3.*(ainv*ainv-2*omeplam/omep*a/amax/amax/amax);
	apsq = 8.*PI/3.*(ainv+1./omei-1+omeplam/omep*(a*a-1)/amax/amax/amax);
	afact = (2.+app*a/apsq)/2.*ainv;
	bfact = ainv*ainv*ainv/apsq;
	vfact1 = (1.-afact*astep)/(1.+afact*astep);
	vfact2 = bfact*astep/(1.+afact*astep)*rng;
	vfact1h = (1.-afact*asteph)/(1.+afact*asteph);
	vfact2h = bfact*asteph/(1.+afact*asteph)*rng;
	gv1 = vfact1/vfact1h;
	gv2 = (vfact2-vfact1*vfact2h)/vfact2;
	*/


	tp = slcparticles;
	for(i=0;i<snlcp;i++){
		bp = tp->bp;
		if(isnowdmstep(nowTsubdiv,GetTsubPower(bp),maxTsubpower)){
			float vfact2,vfact3;
			vfact2 = EvolFactKick[GetTsubPower(bp)].fact2;
			ax = (tp->nvx-tp->ovx)/vfact2;
			ay = (tp->nvy-tp->ovy)/vfact2;
			az = (tp->nvz-tp->ovz)/vfact2;
			vfact3 = EvolFactPull[GetTsubPower(bp)].fact2;
			tp->vx = tp->ovx + ax * vfact3;
			tp->vy = tp->ovy + ay * vfact3;
			tp->vz = tp->ovz + az * vfact3;
		}
		else {
			tp->vx = tp->ovx;
			tp->vy = tp->ovy;
			tp->vz = tp->ovz;
		}
		tp++;
	}
	nnp[0] = nnp[1] = nnp[2] =0;
	for(i=0;i<snlcp;i++){
		nnp[slcparticles[i].type]++;
	}
	if(snlcp>0) qsort(slcparticles,snlcp,sizeof(slcparticletype),ssorttype);
	if(snlcp>0) lightcones = (lightconetype *)Malloc(sizeof(lightconetype)*snlcp,
			PPTR(lightcones));
	for(i=0;i<snlcp;i++){
		lightcones[i].x = slcparticles[i].x;
		lightcones[i].y = slcparticles[i].y;
		lightcones[i].z = slcparticles[i].z;
		lightcones[i].vx = slcparticles[i].vx;
		lightcones[i].vy = slcparticles[i].vy;
		lightcones[i].vz = slcparticles[i].vz;
#ifdef INDEX
		lightcones[i].indx = slcparticles[i].indx;
#endif
	}
	{
		lightconetype *pp;

		sprintf(surveyfilename,"Olightcone0.main.%.5d",mstep);
		pp = lightcones;
		writedownlightconedata(surveyfilename,pp,nnp[0]);

		sprintf(surveyfilename,"Olightcone0.outer.%.5d",mstep);
		pp += nnp[0];
		writedownlightconedata(surveyfilename,pp,nnp[1]);

		sprintf(surveyfilename,"Olightcone0.inner.%.5d",mstep);
		pp += nnp[1];
		writedownlightconedata(surveyfilename,pp,nnp[2]);
	}
	if(snlcp>0) Free(lightcones);
	Free(slcparticles);
}
#undef npstep
