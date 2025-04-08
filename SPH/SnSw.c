#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "nnost.h"
#include "sph.h"
#include "indT.h"
#include "hydroBasicCell.h"
#include "CosmosEvolFactor.h"

#define MAX(a,b) (a)>(b)?(a):(b)

#define epsilon2 0.01

#define invLog2 1.4426950408


float ran3( long *);
long isseed=-987;
int seedflag=1;
void calcsnfb_(float *, float*, double*, float *, float*,float *,float*);
void calcswfb_(float *, float *,float *, float *,float *);
void calcsnprob_( float *, float *, float *, float *);

/* This updates mass, metallicity and dAs. */

int Complete_SN_SW(SimParameters *simpar, particle *threadpos, TPtlStruct **bpneighbor,int Numnear,double *dT){
	int i,j,k;
	int nexplosion=0;

	treestarparticletype *bp;
	bp = LinkedTreeStarType(threadpos);
	float tdelay = bp->sftime - LOOKBACKTIME(simpar); /* in Myr */
	float metal = bp->metallicity;
	int ii = GetTsubPower(bp);
	float tstep = dT[ii];
	float prob;
#ifndef DDDDDDDDD	
	if(tdelay > 10  && bp->SNmark == NEVER){
		float xxx=0;
		calcsnprob_(&metal,&tdelay,&tstep,&prob);
		xxx = ran3(&isseed);
		if(xxx < prob){
			double energyE;
			float massE,xE,yE,zE;
			float massI;
			float gamma = GAS_GAMMA(simpar);
			massI = bp->mass * SPH_MASSINSOLARMASS(simpar);
			// metal = Z/(x+Y+Z) 
			// massI = in Msun
			// energyE = in erg
			// massE = in Msun
			// xE = in Msun 
			// yE = in Msun 
			// zE = in Msun 
			if(MYID(simpar)==0){
                               	printf("SN feedback(ssp) %g %g %g \n",bp->mass,SPH_MASSINSOLARMASS(simpar),massI);
			}
			calcsnfb_(&metal, &massI, &energyE, &massE, &xE, &yE, &zE);
			// massE : ejected mass in simulation mass 
			massE = massE/SPH_MASSINSOLARMASS(simpar);
			bp->mass -= massE; // loosing stell mass due to the SN exposion 
			// zE : ejected metals in simulation mass 
			zE = zE/SPH_MASSINSOLARMASS(simpar);
			float xi,yi,zi;
			xi = threadpos->x;
			yi = threadpos->y;
			zi = threadpos->z;
			float sumsolid=0;
			for(i=0;i<Numnear;i++){
				float dist2;
				dist2  = (bpneighbor[i]->x - xi)*(bpneighbor[i]->x - xi);
				dist2 += (bpneighbor[i]->y - yi)*(bpneighbor[i]->y - yi);
				dist2 += (bpneighbor[i]->z - zi)*(bpneighbor[i]->z - zi);
				dist2 += epsilon2;
				sumsolid += pow(bpneighbor[i]->rho,-0.666666666)/dist2;
			}
			float gammam1 = GAS_GAMMA(simpar)-1;
			for(i=0;i<Numnear;i++){
				float myfrac;
				float dist2;
				dist2  = (bpneighbor[i]->x - xi)*(bpneighbor[i]->x - xi);
				dist2 += (bpneighbor[i]->y - yi)*(bpneighbor[i]->y - yi);
				dist2 += (bpneighbor[i]->z - zi)*(bpneighbor[i]->z - zi);
				dist2 += epsilon2;
				treesphparticletype *tbp = LinkedTreeSphType(bpneighbor[i]);
				double pmassincgs = tbp->mass*SPH_MASSINSOLARMASS(simpar)*onesolarmass;
				myfrac = pow(bpneighbor[i]->rho,-0.666666666)/dist2 / sumsolid;
				float Mz = tbp->mass * tbp->metallicity;
				tbp->mass += myfrac * massE;
				tbp->metallicity = (Mz+myfrac*zE)/tbp->mass;
				tbp->Entropy +=myfrac*gammam1*energyE/
					(pmassincgs*pow(tbp->rho,gammam1));
				/*
				if(tbp->indx == 47695661821465){
					DEBUGPRINT("SN effect %ld :  %g %g %g %g %g %g\n",tbp->indx,myfrac,Mz,zE,energyE,pmassincgs,tbp->rho);
				}
				*/
			}

			/* This is for Time Step Limiter (Saitoh) */
			for(i=0;i<Numnear;i++){
				treesphparticletype *pj =  LinkedTreeSphType(bpneighbor[i]);
				double Csi = sqrt(pj->Entropy*GAS_GAMMA(simpar)*pow(pj->rho, GAS_GAMMA(simpar)-1));
				float dai = GAS_COURANT(simpar)*pj->hsml/(Csi/GAS_G2(simpar))/NX(simpar);
				float ndiv = ASTEP(simpar)/dai;
				int inew = max(0,ceil(logf(ndiv)*invLog2));
				short int iold = GetSphTsubPower(pj);
				short int iold2 = GetSphNextTsubPower(pj);
				iold = (iold>iold2? iold: iold2);
				if(inew > iold && inew >0 ) SetSphNextTsubPower(pj,inew);
				/*
				if(pj->indx == 47695661821465){
					DEBUGPRINT("SN effect %ld :  %d %d %d : %g\n",pj->indx,inew,iold,iold2,Csi);
				}
				*/
			}


			bp->SNmark=ONCE;
			bp->swlasttime = (tdelay>60 ? 60:tdelay);
			nexplosion = 1;
		}
	}else if(bp->SNmark == ONCE){
		float massI,massE,zE; 
		massI = GAS_MSTAR(simpar) * SPH_INITMASS(simpar);
		calcswfb_(&metal,&massI,&tdelay, &(bp->swlasttime),&massE);
		if(massE>0){
			massE = massE/SPH_MASSINSOLARMASS(simpar);
			bp->mass -= massE;
			zE = bp->metallicity * massE;
			float xi,yi,zi;
			xi = threadpos->x;
			yi = threadpos->y;
			zi = threadpos->z;
			float sumsolid=0;
			for(i=0;i<Numnear;i++){
				float dist2;
				dist2  = (bpneighbor[i]->x - xi)*(bpneighbor[i]->x - xi);
				dist2 += (bpneighbor[i]->y - yi)*(bpneighbor[i]->y - yi);
				dist2 += (bpneighbor[i]->z - zi)*(bpneighbor[i]->z - zi);
				dist2 += epsilon2;
				sumsolid += pow(bpneighbor[i]->rho,-0.666666666)/dist2;
			}
			for(i=0;i<Numnear;i++){
				float myfrac;
				float dist2;
				dist2  = (bpneighbor[i]->x - xi)*(bpneighbor[i]->x - xi);
				dist2 += (bpneighbor[i]->y - yi)*(bpneighbor[i]->y - yi);
				dist2 += (bpneighbor[i]->z - zi)*(bpneighbor[i]->z - zi);
				dist2 += epsilon2;
				treesphparticletype *tbp = LinkedTreeSphType(bpneighbor[i]);
				myfrac = pow(bpneighbor[i]->rho,-0.666666666)/dist2 / sumsolid;
				float Mz = tbp->mass * tbp->metallicity;
				tbp->mass += myfrac * massE;
				tbp->metallicity = (Mz+myfrac*zE)/tbp->mass;
			}
		}
	}
#endif
	
	return nexplosion;
}

#define NSBUFF 100000
void MakeStar(SimParameters *simpar){
	long long nsbuff=NSBUFF,nssize,i;
	double lookbacktime; /* in real time */
//jhshin1
	float sf_epsilon=1e-8;  /* offset distance gas-> sf */	
//jhshin2
        treesphparticletype *p = SPH_TBP(simpar);
	treestarparticletype *sbp;
	int invMstar = rint(1./GAS_MSTAR(simpar));
	long newstar=0;

	if(seedflag){
		isseed *= (MYID(simpar)+1);
		seedflag = 0;
	}
	if(STAR_NP(simpar)==0){
		nssize = nsbuff;
		STAR_TBP(simpar)  = (treestarparticletype*)malloc(sizeof(treestarparticletype)*nssize);
		sbp = STAR_TBP(simpar);
	}
	else {
		nssize = STAR_NP(simpar)+nsbuff;
		STAR_TBP(simpar)  = (treestarparticletype*)realloc(STAR_TBP(simpar),sizeof(treestarparticletype)*nssize);
		sbp = STAR_TBP(simpar) + STAR_NP(simpar);
	}
	if(MYID(simpar)==0){
		printf("######################################\n");
		printf("Entering into star formation \n");
		printf("######################################\n");
	}
	{
		float red;
		float Hzzp1(float);
		float qsimp(float (*func)(float),float,float);
		double H0incgs = H0/Mpccgs*HUBBLE(simpar);
		red = AMAX(simpar)/ANOW(simpar)-1;
		if(BGEXPAND(simpar) == 'Y') lookbacktime = qsimp(Hzzp1,0,red)/(H0incgs*Myr);
		else lookbacktime = (AMAX(simpar)-ANOW(simpar))*FFTIME(simpar)/Myr;/* rescale to the real time */
	}
   	p = SPH_TBP(simpar);
   	for(i = 0;i<SPH_NP(simpar);i++){
		int ii = GetTsubPower(p);
		int jj = GetSphTsubPower(p);
		int kk = GetSphNextTsubPower(p);
		if(ii>jj) jj = ii;
		if(p->divVel <0 && IsNowSphStep(nowTsubdiv,jj,IndT_MAXTSUBPOWER(simpar))){ 
			if(p->rho > GAS_SFRHO577(simpar) || BGEXPAND(simpar) == 'N'){ 
				float Temp = GetT(simpar, p); 
				if(Temp < GAS_SFTEMP(simpar)){ 
					double gasdenincgs = p->rho*GAS_RHOS2RHOR(simpar); 
					double totalhydrogenmass = gasdenincgs*(1-GAS_YP(simpar)-p->metallicity); 
					double nH = totalhydrogenmass/mH; 
					if(nH > GAS_SFGASDEN(simpar)){ 
						double tdynamic = 1./sqrt(P_G*gasdenincgs);
						double Dt=0;
						if(CONT_HALFSTEP(simpar)==KICK){
							if(kk>=jj){
								Dt = dTKICKChangeSubStep[ii][jj][kk];
							}
							else if(IsNowSphStep(nowTsubdiv,(jj-1),IndT_MAXTSUBPOWER(simpar))){
								Dt = dTKICKChangeSubStep[ii][jj][jj-1];
							}
							else{
								Dt = dTKICKChangeSubStep[ii][jj][jj];
							}
						}
						else if(CONT_HALFSTEP(simpar)==HALFPUSH){
							if(jj>kk) Dt = dTPUSH[jj];
							else Dt = dTPUSH[kk];
						}
						else {
							Dt = dTPULL[jj];
						}
						float prob = (1-exp(-GAS_CSTAR(simpar)/GAS_MSTAR(simpar)*Dt/tdynamic));
						float xxx = ran3(&isseed);
						if(xxx < prob){
							float gasmass = p->mass;
							float orggasmass = p->mass /(1-GAS_MSTAR(simpar)*p->nstar);
							p->mass -= orggasmass * GAS_MSTAR(simpar);
							p->nstar ++;
							sbp->mass = orggasmass*GAS_MSTAR(simpar);
							sbp->metallicity = p->metallicity;
							sbp->sftime = lookbacktime;
							sbp->x = p->x;
							sbp->y = p->y;
							sbp->z = p->z;
							sbp->vx = p->vx;
							sbp->vy = p->vy;
							sbp->vz = p->vz;
							sbp->ax = p->ax;
							sbp->ay = p->ay;
							sbp->az = p->az;
							size_t newindx = PINDX(p) * invMstar + p->nstar;
							CHANGEINDX(p, newindx);
							SetNbodyIndT(simpar,sbp);
							/*
							sbp->indt[0] =  p->indt[0];
							sbp->indt[1] =  p->indt[1];
							*/
							sbp->SNmark = NEVER;
							newstar ++;
							sbp++;
							STAR_NP(simpar)++;
							if(STAR_NP(simpar) >= nssize){
								nssize += NSBUFF;
								STAR_TBP(simpar) = (treestarparticletype*)realloc(STAR_TBP(simpar),
									sizeof(treestarparticletype)*nssize);
								sbp = STAR_TBP(simpar) + STAR_NP(simpar);
							}
						} 
								}
                		}
            		}
        	}
        	p++;
    	}
   	p = SPH_TBP(simpar);
	for(i = 0;i<SPH_NP(simpar);i++){
		if(SPH_TBP(simpar)[i].nstar < invMstar){
			*p = *(SPH_TBP(simpar)+i);
			p++;
		}
	}
	{
		long tnewstar;
		MPI_Reduce(&newstar,&tnewstar,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Bcast(&tnewstar,1,MPI_LONG,0,MPI_COMM_WORLD);
		if(MYID(simpar) ==0 && tnewstar>0) {
			float time;
			if(BGEXPAND(simpar) == 'Y') time = AMAX(simpar)/ANOW(simpar)-1;
			else time = ANOW(simpar);
			printf("###################################################\n");
			printf("np= %ld stars are made at z(t)= %g\n",tnewstar,time);
			printf("###################################################\n\n");
		}
	}
	SPH_NP(simpar) = (p-SPH_TBP(simpar));
	SPH_TBP(simpar) = (treesphparticletype*)realloc(SPH_TBP(simpar), sizeof(treesphparticletype)*SPH_NP(simpar));

	if(STAR_NP(simpar)>0) STAR_TBP(simpar) = (treestarparticletype*)realloc(STAR_TBP(simpar),
				sizeof(treestarparticletype)*STAR_NP(simpar));
	else free(STAR_TBP(simpar));
	{
		long long tsnp;
		MPI_Reduce(&STAR_NP(simpar),&STAR_TNP(simpar),1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Bcast(&STAR_TNP(simpar),1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	}
}   
