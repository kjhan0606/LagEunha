#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "nnost.h"
#include "sph.h"
#include "fixedneighbor.h"
#include "indT.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))


#define invLog2 1.4426950408


long long NsizeBytes;

NeighborList *neighborlist = NULL;
long long Nneighborlist,maxNneighborlist = 1000000;
int iup,idown,iwidth;
treesphparticletype **sphPadd,**tbp;
long long *nrequest,*request;


float rho_w(float, float);
float drhodr_w(float, float);



static int sortrequest(const void *a, const void *b){
	long long *aa,*bb;
	aa = (long long*)a;
	bb = (long long*)b;
	if(*aa > *bb) return 1;
	else if (*aa < *bb) return -1;
	else return 0;
}
static int sortNeighborListpi(const void *a, const void *b){
	NeighborList *aa,*bb;
	aa = (NeighborList*)a;
	bb = (NeighborList*)b;
	if(aa->pi>bb->pi) return 1;
	else if(aa->pi < bb->pi) return -1;
	else return 0;
}
static int sortNeighborListpidnj(const void *a, const void *b){
	NeighborList *aa,*bb;
	aa = (NeighborList*)a;
	bb = (NeighborList*)b;
	if(aa->pid>bb->pid)  return 1;
	else if(aa->pid < bb->pid) return -1;
	else {
		if(aa->nj>bb->nj) return 1;
		else if(aa->nj<bb->nj) return -1;
		else return 0;
	}
}
treesphparticletype *localPstart;
treesphparticletype *localPfinal;

void InitFixedNeighborList(SimParameters *simpar){
	localPstart = SPH_TBP(simpar);
	localPfinal = localPstart + SPH_NP(simpar);
	NsizeBytes = sizeof(NeighborList)*SPH_NUMNEAR(simpar);
	Nneighborlist = 0;
	neighborlist = (NeighborList *)malloc(maxNneighborlist *NsizeBytes);
	maxNneighborlist = 1000000;
}
void FreeFixedNeighborList(SimParameters *simpar){
	if(neighborlist != NULL){
		Nneighborlist = 0;
		free(neighborlist);
	}
	if(request != NULL) free(request);
	if(nrequest != NULL) free(nrequest);
}

void GetFixedSubStepInfo(SimParameters *simpar){
	int ii,jj,kk,i;
	int nbodypow,sphpow,fixedpow;
	treesphparticletype *bp = SPH_TBP(simpar);
		
	ii = jj = kk = 0;
	for(i=0;i<SPH_NP(simpar);i++){
		int ia = GetSphTsubPower(bp);
		int ib = GetFixedPower(bp);
		int ic = GetTsubPower(bp);
		if(ia > ii) ii = ia;
		if(ib > jj) jj = ib;
		if(ic > kk) kk = ic;
	}
	MPI_Reduce(&ii,&sphpow,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
	MPI_Reduce(&jj,&fixedpow,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
	MPI_Reduce(&kk,&nbodypow,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
	MPI_Bcast(&sphpow,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&fixedpow,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nbodypow,1,MPI_INT,0,MPI_COMM_WORLD);
	IndT_NSUBGAS(simpar) = sphpow;
	IndT_NSUBFIXED(simpar) = fixedpow;
	IndT_NSUBNBODY(simpar) = nbodypow;
}

void InsertFixedNeighborList(SimParameters *simpar, treesphparticletype *p,int Num_near, TPtlStruct **bpneighbor){
	long long i;
	{
		for(i=0;i<Num_near;i++){
			treesphparticletype *pj = LinkedTreeSphType(bpneighbor[i]);
			if(pj>=localPstart && pj<localPfinal) {
				insert_local_neighbor(simpar,p,Nneighborlist,pj);
			}
			else {
				fprintf(stderr,"There is something wrong here\n");
				exit(99);
			}
			Nneighborlist ++;
		}
	}
	if(Nneighborlist >= maxNneighborlist-SPH_NUMNEAR(simpar)){
		maxNneighborlist += 1000000;
		neighborlist = (NeighborList *)realloc(neighborlist,maxNneighborlist*NsizeBytes);
	}
}

void GetSphParticleAddress(SimParameters *simpar){
	int min,max,tmin,tmax;
	min = NID(simpar);
	max = -NID(simpar);
	long long i,j,ineighborlist=0;
	/*
	for(i=0;i<Nneighborlist;i++, ineighborlist++){
		if(neighborlist[i].pid >=0) neighborlist[ineighborlist] = neighborlist[i];
	}
	Nneighborlist = ineighborlist;
	*/
	for(i=0;i<Nneighborlist;i++){
		NeighborList *p = neighborlist+i;
		int diff = p->pid-MYID(simpar);
		min = MIN(min,diff);
		max = MAX(max,diff);
	}
	MPI_Reduce(&min,&tmin,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
	MPI_Reduce(&max,&tmax,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
	MPI_Bcast(&tmin,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&tmax,1,MPI_INT,0,MPI_COMM_WORLD);

	iup = tmax;
	idown = -tmin;
	iwidth = MAX(idown,iup);
	if(MYID(simpar)==0) 
		printf("P%d has Nneighborlist %ld : iup/idown iwidth = %d %d %d \n",MYID(simpar),Nneighborlist,iup,idown,iwidth);

	
	sphPadd = (treesphparticletype**)malloc(sizeof(treesphparticletype*)*(2*iwidth+1));
	tbp = sphPadd + iwidth;

	for(i=-iwidth;i<=iwidth;i++){
		MPI_Status status;
		if(i==0) continue;
		int src = (MYID(simpar)+i + NID(simpar))%NID(simpar);
		int tgt = (MYID(simpar)-i + NID(simpar))%NID(simpar);
		MPI_Sendrecv(SPH_TBP(simpar),sizeof(treesphparticletype*),MPI_BYTE,tgt,0,
				&tbp[i],sizeof(treesphparticletype*),MPI_BYTE,src,0,MPI_COMM_WORLD,&status);
	}
	tbp[0] = SPH_TBP(simpar);

	for(i=0;i<Nneighborlist;i++){
		NeighborList *p = neighborlist+i;
		int diff = p->pid-MYID(simpar);
		p->nj = p->pj- tbp[diff];
		if(MYID(simpar) == p->pid) p->njnow = p->nj;
	}
	free(sphPadd);
	if(MYID(simpar)==0) printf("Now exiting GetSphParticleAddress\n");

}
void GetInfoOfNeighborWorkList(SimParameters *simpar){
	long long i,j;
	long long nrequestrunning=0;

	qsort(neighborlist,Nneighborlist,sizeof(NeighborList),sortNeighborListpidnj);
	nrequest = (long long *)malloc(sizeof(long long)*(2*iwidth+1));
	request = (long long *)malloc(1);

	for(i=-iwidth;i<=iwidth;i++){
		MPI_Status status;
		if(i==0) continue;
		int src = (MYID(simpar)+i + NID(simpar))%NID(simpar);
		int tgt = (MYID(simpar)-i + NID(simpar))%NID(simpar);
		long long nsend=0,*send;
		nrequest[i+iwidth] = 0;
		for(j=0;j<Nneighborlist;j++){
			if(neighborlist[j].pid == tgt){
				nsend ++;
			}
		}
		if(nsend>0) send = (long long*)malloc(sizeof(long long)*nsend);
		nsend = 0;
		for(j=0;j<Nneighborlist;j++){
			if(neighborlist[j].pid == tgt){
				send[nsend++] = neighborlist[j].nj;
			}
		}
		/* erase duplicated array in send */
		long long prev=0,curr=1;
		while(curr< nsend){
			for(; send[prev] == send[curr] && curr < nsend ;curr++);
			if(curr!=nsend) send[++prev] = send[curr];
		} 
		nsend = prev;
		for(j=0;j<nsend;j++){
			if(send[j]<0)printf("P%d erorr %ld : %ld\n",MYID(simpar),send[j],j);
		}

		MPI_Sendrecv(&nsend,1,MPI_LONG_LONG,src,0,nrequest+i+iwidth,1,MPI_LONG_LONG,tgt,0,MPI_COMM_WORLD,&status);
		if(nrequest[i+iwidth]>0) request = (long long*)realloc(request,sizeof(long long)*(nrequest[i+iwidth]+nrequestrunning));
		MPI_Sendrecv(send,nsend,MPI_LONG_LONG,src,0,request+nrequestrunning,nrequest[i+iwidth],
				MPI_LONG_LONG,tgt,0,MPI_COMM_WORLD,&status);
		if(nsend > 0) free(send);

		nrequestrunning += nrequest[i+iwidth];
	}
}

void DoSphDenFixedNeighbor(SimParameters *simpar){
	int i,j;
	for(i=0;i<Nneighborlist;){
		treesphparticletype *pi = neighborlist[i].pi;
		if(!IsNowOneOfTwoSteps(IndT_NOWTSUBDIV(simpar),GetTsubPower(pi),GetSphTsubPower(pi),IndT_MAXTSUBPOWER(simpar))){
			i++;
			continue;
		}
		float tmpx,tmpy,tmpz,r,hsml;
		float tmpvx,tmpvy,tmpvz,vfact,divVel,curlVx,curlVy,curlVz,Cs; 
		float massj,massi,density=0,density09=0,hsml09; 
		divVel = curlVx = curlVy = curlVz = 0;
		hsml = pi->hsml;
		hsml09 = hsml*0.99;
		for(j=i;neighborlist[j].pi== pi && j < Nneighborlist;j++){
			long long nj = neighborlist[j].njnow;
			float x,y,z,vx,vy,vz,mass;
			treesphparticletype *pj;
			if(neighborlist[j].pid == MYID(simpar)) pj = (treesphparticletype *)neighborlist[j].pj;
			else fprintf(stderr,"There is something wrong here\n");
			tmpx = pj->x - pi->x;
			tmpy = pj->y - pi->y;
			tmpz = pj->z - pi->z;
			if(tmpx < -NX(simpar)*0.5) tmpx += NX(simpar);
			if(tmpy < -NY(simpar)*0.5) tmpy += NY(simpar);
			if(tmpz < -NZ(simpar)*0.5) tmpz += NZ(simpar);
			if(tmpx < NX(simpar)*0.5) tmpx -= NX(simpar);
			if(tmpy < NY(simpar)*0.5) tmpy -= NY(simpar);
			if(tmpz < NZ(simpar)*0.5) tmpz -= NZ(simpar);
			r = tmpx*tmpx; r+= tmpz*tmpz; r+= tmpy*tmpy;
			r = sqrtf(r);
			density += rho_w(r,hsml)*pj->mass;
			density09 += rho_w(r,hsml09)*pj->mass;

			tmpvx = GAS_G12RATIO(simpar)*(pj->vx-pi->vx) + tmpx*GAS_G1EXPANSION(simpar);
			tmpvy = GAS_G12RATIO(simpar)*(pj->vy-pi->vy) + tmpy*GAS_G1EXPANSION(simpar);
			tmpvz = GAS_G12RATIO(simpar)*(pj->vz-pi->vz) + tmpz*GAS_G1EXPANSION(simpar);
			if(r>0){
				vfact = pj->mass * (-drhodr_w(r,hsml))/r;
				divVel += vfact * (tmpx*tmpvx + tmpy*tmpvy + tmpz*tmpvz);
				curlVx += vfact * ( tmpvy*tmpz - tmpvz * tmpy);
				curlVy += vfact * ( tmpvz*tmpx - tmpvx * tmpz);
				curlVz += vfact * ( tmpvx*tmpy - tmpvy * tmpx);
			}
		}

		pi->rho = density;
		pi->fi = 1./(1+1./3.*(log(density)-log(density09))/(log(hsml)-log(hsml09)));
		pi->divVel = divVel;
		float curlVel = sqrtf(curlVx*curlVx + curlVy*curlVy + curlVz*curlVz);
		Cs = sqrt(pi->Entropy*GAS_GAMMA(simpar)*pow(pi->rho,GAS_GAMMA(simpar)-1));
		divVel = fabsf(divVel)/density;
		curlVel = fabsf(curlVel)/density;
		pi->F1 = divVel/(divVel + curlVel + 0.0001*Cs/hsml/GAS_G1(simpar));

		i = j;
	}
}

void DoSphComFixedNeighbor(SimParameters *simpar){
	int i,j;
	for(i=0;i<Nneighborlist;){
		treesphparticletype *pi = neighborlist[i].pi;
		if(!IsNowOneOfTwoSteps(IndT_NOWTSUBDIV(simpar),GetTsubPower(pi),GetSphTsubPower(pi),IndT_MAXTSUBPOWER(simpar))){
			i++;
			continue;
		}
		float tmpx,tmpy,tmpz,rij,hsmli,hsmlj,invrij; 
		float dx,dy,dz; 
		float tmpvx,tmpvy,tmpvz,vfact,divVel,curlVx,curlVy,curlVz,Cs; 
		float Tempiovermui,massi,massj,Tempjovermuj,Viscfactor; 
		float fi,rhoi,mui,Csi,invCsi,gamma; 
		float Ei,Ej; 
		float fj,rhoj,muj,Csj; 
		float MySphForce,G1overCsi; 
		float dAs,alpha,F1; 
		float dlna,Tmui,Tmuj; 
		float G2overG1 = GAS_G12RATIO(simpar); 
		gamma = GAS_GAMMA(simpar);
		alpha = GAS_ALPHA(simpar);

		hsmli = pi->hsml;
		F1 = pi->F1;

		rhoi = pi->rho;
		Ei = pi->Entropy;
		Tempiovermui = Ei*mHoverkB*pow(rhoi,gamma-1);
		Csi = sqrtf(gamma*Tempiovermui*kBovermH);
		Viscfactor = GAS_VISFORCEFACTOR(simpar) * Tempiovermui;
		invCsi = 1./Csi;
		massi = pi->mass;
		MySphForce = (pi->fi)*GAS_SPHFORCEFACTOR(simpar)*Tempiovermui/rhoi;
		G1overCsi = GAS_G1(simpar)/Csi;
		dAs = 0;
		float maxSignalVel = 0;
		double ax,ay,az;
		ax = ay = az = 0;
		for(j=i;neighborlist[j].pi == pi && j <Nneighborlist;j++){
			long long nj = neighborlist[j].njnow;
			float F12, F2;
			treesphparticletype *pj;
			if(neighborlist[j].pid == MYID(simpar)) pj = (treesphparticletype*)neighborlist[j].pj;
			else {
				DEBUGPRINT2("Something wrong here\n");
			}
			massj = pj->mass;
			rhoj = pj->rho;
			Ej = pj->Entropy;
			F2 = pj->F1;
			hsmlj = pj->hsml;
			dx = pj->x - pi->x;
			dy = pj->y - pi->y;
			dz = pj->z - pi->z;
			if(dx < -0.5*NX(simpar)) dx += NX(simpar);
			if(dy < -0.5*NY(simpar)) dy += NY(simpar);
			if(dz < -0.5*NZ(simpar)) dz += NZ(simpar);
			if(dx >  0.5*NX(simpar)) dx -= NX(simpar);
			if(dy >  0.5*NY(simpar)) dy -= NY(simpar);
			if(dz >  0.5*NZ(simpar)) dz -= NZ(simpar);
			rij = dx*dx + dy*dy + dz*dz;
			if(rij>0) {
				rij = sqrtf(rij);
				invrij = 1./rij;
				tmpx = invrij*dx;
				tmpy = invrij*dy;
				tmpz = invrij*dz;

				Tempjovermuj = Ej*mHoverkB*pow(rhoj,gamma-1);

				float drhodri = -drhodr_w(rij,hsmli);
				float drhodrj = -drhodr_w(rij,hsmlj);
				vfact = MySphForce*drhodri;

				float fact1 = vfact*massj;
				float fact2 = (pj->fi)*GAS_SPHFORCEFACTOR(simpar)*massj*Tempjovermuj/rhoj*drhodrj;
				float fact12 = fact1 + fact2;
				ax -= tmpx*fact12;
				ay -= tmpy*fact12;
				az -= tmpz*fact12;

				float v21x  = dx*GAS_G1EXPANSION(simpar)+(pj->vx-pi->vx)*G2overG1;
				float v21y  = dy*GAS_G1EXPANSION(simpar)+(pj->vy-pi->vy)*G2overG1;
				float v21z  = dz*GAS_G1EXPANSION(simpar)+(pj->vz-pi->vz)*G2overG1;


				float wsij = tmpx*v21x + tmpy*v21y + tmpz*v21z;
				Csj = sqrtf(gamma*Tempjovermuj*kBovermH);

				float wrijoverCsi = G1overCsi*wsij;
				float SigVel = (1+Csj*invCsi-3*wrijoverCsi);
				float SigVelCsi = SigVel*Csi;
				if(SigVelCsi > maxSignalVel) maxSignalVel = SigVelCsi;
				if(wsij<0){
					float invrhoij = 2./(rhoi+rhoj); 
					float drhodrij = 0.5*(drhodri+drhodrj); 
					float xGradW = tmpx*drhodrij; 
					float yGradW = tmpy*drhodrij; 
					float zGradW = tmpz*drhodrij; 
					F12 = (F1+F2)*0.5; 
					float ViscForceJ = massj*SigVel*wrijoverCsi*invrhoij*F12; 
					float Visforcefact = Viscfactor*ViscForceJ; 
					ax += Visforcefact*xGradW; 
					ay += Visforcefact*yGradW; 
					az += Visforcefact*zGradW; 
					dAs += ViscForceJ*wsij*drhodrij;
				}
			}
		}
		pi->ax += ax;
		pi->ay += ay;
		pi->az += az;
		dAs = GAS_ENTROPYEVOLFACT(simpar) * dAs;
		pi->dAs += Ei*(GAS_G1EXPANSION(simpar)*3*(1.-gamma)+dAs);

		if(maxSignalVel>0){
			float dai = GAS_COURANT(simpar)*hsmli/(maxSignalVel/GAS_G2(simpar))/NX(simpar);
			float ndiv = ASTEP(simpar)/dai;
			int isubpow = MAX(0,ceil(logf(ndiv)*invLog2));
			int imax = MAX(GetSphTsubPower(pi),isubpow);
			if(imax>32) {
				DEBUGPRINT("Error in setting sph next step sub power:  %g %g : %d %d : %d\n", 
						dai,ndiv,isubpow,imax,GetSphTsubPower(pi));
				DEBUGPRINT("Error :  %g %g : %g %g \n",GAS_COURANT(simpar),hsmli,maxSignalVel,GAS_G2(simpar));
				DEBUGPRINT("Error :  Csi = %g T/mu= %g \n",Csi,Tempiovermui);
				exit(999);
			}
			else if(imax>=SPH_CONSTNEIGHPOW(simpar)){
			}
			SetSphNextTsubPower(pi,imax);
			if(imax>SPH_TIMESTEPLIMITER(simpar)) Fixed_TimeStep_limiter(simpar, imax,i,pi);
		}

	}
}


int determine_fixed_timestep(SimParameters *simpar){
	int maxfixedpower=0,maxsphsubpower=0;
	long long i;
	float fact1 = (NX(simpar)*ASTEP(simpar))/GAS_DURANT(simpar);
	treesphparticletype *bp = SPH_TBP(simpar);
	for(i=0;i<SPH_NP(simpar);i++){
		/*
		float ndiv = fact1*fabsf(bp->divVel) *bp->hsml;
		int isubpow = MAX(0,ceilf(logf(ndiv)*invLog2));
		int sphisubpow = GetSphTsubPower(bp);
		isubpow = MAX(isubpow,sphisubpow);
		SetFixedPower(bp,isubpow);
		*/
		int sphisubpow = GetSphTsubPower(bp);
		int isubpow = GetFixedPower(bp);

		maxfixedpower = MAX(maxfixedpower,isubpow);
		maxsphsubpower = MAX(maxsphsubpower,sphisubpow);
		bp++;
	}
	int tmaxfixedpower,tmaxsphsubpower;
	MPI_Reduce(&maxfixedpower,&tmaxfixedpower,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
	MPI_Bcast(&tmaxfixedpower,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Reduce(&maxsphsubpower,&tmaxsphsubpower,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
	MPI_Bcast(&tmaxsphsubpower,1,MPI_INT,0,MPI_COMM_WORLD);
	int maxdiff = tmaxsphsubpower-tmaxfixedpower;
	if(MYID(simpar)==0) {
		printf("Fixed Neight Test: fixed= %d sph=%d : diff = %d\n",tmaxfixedpower,tmaxsphsubpower,maxdiff);
	}
	return maxdiff;
}

void fixedsphmain(SimParameters *simpar){
	int Numnear = SPH_NUMNEAR(simpar);
	float astep = ASTEP(simpar)/pow(2.,IndT_MAXTSUBPOWER(simpar));

	if(BGEXPAND(simpar) == 'Y'){
		void GetEvolFactorArr(SimParameters *,float, float );
		void DetermineSphFactor(SimParameters *);
		DetermineSphFactor(simpar);
		GetEvolFactorArr(simpar, ANOW(simpar),astep);
	}
	else {
		void StaticGetEvolFactorArr(SimParameters *,float, float);
		void StaticDetermineSphFactor(SimParameters *);
		StaticDetermineSphFactor(simpar);
		StaticGetEvolFactorArr(simpar, ANOW(simpar),astep);
	}
	Init_dAs(simpar);

	DoSphDenFixedNeighbor(simpar);

	/* Only necessary at the starting of the simulation.
	if(init_flag == 1){
		double redshift = simpar->amax/simpar->anow-1;
		void InitializeMu(double); InitializeMu(redshift);
		if(simpar->BGCosmology=='Y'){
			simpar->sph.initmass = simpar->sph.massinsolarmass;
			treesphparticletype *sp = simpar->sph.u.tbp;
			float temp;
			long ii;
			for(ii=0;ii<simpar->sph.np;ii++){
				temp = GetT(sp);
				sp->mu = getmu(temp,sp->rho);
				sp++;
			}
		}
		else if(simpar->SimModel==Galaxy){
			if(simpar->stepcount ==1){
				long long ii;
				if(MYID(simpar)==0) printf("jhshin's Entropy initializing.\n");
				treesphparticletype *sp = simpar->sph.u.tbp;
				for(ii=0;ii<simpar->sph.np;ii++){
					sp->mu = 1;
					sp->Entropy = kBovermH/(sp->mu)*1.E4*pow(sp->rho,1-simpar->sph.gamma);
					sp->temp = GetT(sp);
					sp++;
				}
			}
			else {
				long long ii;
				treesphparticletype *sp = simpar->sph.u.tbp;
				for(ii=0;ii<simpar->sph.np;ii++) sp->mu = 1;
			}

		}
		init_flag = 0;
	}
	*/
	void Calc_mutempnden(SimParameters *); Calc_mutempnden(simpar);
	if(GAS_BGHEATFLAG(simpar)=='Y'|| GAS_COOLFLAG(simpar)=='Y') {
		void CoolingAndHeating(SimParameters *);
		CoolingAndHeating(simpar);
	}
	DoSphComFixedNeighbor(simpar);
	if(GAS_SFFLAG(simpar) == 'Y') {
		void MakeStar(SimParameters *);
		MakeStar(simpar);
	}

}
