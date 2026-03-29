/* A part of the SPH code to find density and update several
 * important information */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<omp.h>
#include<mpi.h>
#include "eunha.h"
#include "hydroBasicCell.h"
#include "nnost.h"
#include "indT.h"
#include "sph.h"
#include "CosmosEvolFactor.h"
#include "timerutil.h"
#ifdef USE_GPU
#include "nvgpu.h"
void gputreeforce(particle *,int, TStruct *,int , TPtlStruct *,int, float, int, int,int);
void Direct_Nbody(particle *,int, TPtlStruct *,int , int, int,int);
#endif


float cputime0[5], cputime1[10];


#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define invLog2 1.4426950408


#ifndef _OPENMP
#define omp_get_thread_num() nullfct0()
#define omp_get_num_threads() nullfct1()
#endif

void i_force_spline(long,float);

size_t STARBuildLinkedList(SimParameters *simpar){
	size_t mx,my,mz;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	mz = BASICCELL_MZ(simpar);
	PosType xmin,ymin,zmin;
	xmin = SIM_LXMIN(simpar,sph);
	ymin = SIM_LYMIN(simpar,sph);
	zmin = SIM_LZMIN(simpar,sph);
	PosType CellWidth = BASICCELL_CELLWIDTH(simpar);

	HydroTreeLinkedCell *STAR_BasicCell = STAR_BASICCELL(simpar);

	ptrdiff_t i;
	for(i=0;i<mx*my*mz;i++){
		STAR_BasicCell[i].link = NULL;
		STAR_BasicCell[i].nmem = 0;
		STAR_BasicCell[i].calflag = 0;
	}
	LinkParticles(simpar, STAR, star, STAR_BasicCell, mx,my,xmin,ymin,zmin, CellWidth);
	size_t npmaxcell = -1;
	for(i=0;i<mx*my*mz;i++){
		npmaxcell = MAX(npmaxcell, STAR_BasicCell[i].nmem);
	}
	return npmaxcell;
}


void CalculateAllSph(SimParameters *simpar, particle *pos, TPtlStruct **bpneighbor,int Num_neighbor){
	size_t i;
	PosType tmpx,tmpy,tmpz,rij,hsmli,hsmlj,invrij;
	PosType dx,dy,dz;
	float tmpvx,tmpvy,tmpvz,vfact,divVel,curlVx,curlVy,curlVz,Cs;
	float Tempiovermui,massi,massj,Tempjovermuj,Viscfactor;
	float fi,rhoi,mui,Csi,invCsi,gamma;
	float Ei,Ej;
	float fj,rhoj,muj,Csj;
	float MySphForce,G1overCsi;
	float dAs,alpha,F1;
	float dlna,Tmui,Tmuj;
	float G2overG1 = GAS_G12RATIO(simpar);
	treesphparticletype *p;

	p = LinkedTreeSphType(pos);
	gamma = GAS_GAMMA(simpar);
	alpha = GAS_ALPHA(simpar);

	hsmli= p->hsml;
	F1 = p->F1;

	rhoi = p->rho;
	Ei = p->Entropy;
	Tempiovermui = Ei*mHoverkB*pow(rhoi,gamma-1);
	Csi = sqrtf(gamma*Tempiovermui*kBovermH); 
	/* GAS_VISFORCEFACTOR(simpar) = 7.8e-9 */
	/* Viscfactor = 1.5E-6 */
	Viscfactor = GAS_VISFORCEFACTOR(simpar) * Tempiovermui;
	invCsi = 1./Csi;
	massi = p->mass;

	MySphForce = (p->fi)*GAS_SPHFORCEFACTOR(simpar)*Tempiovermui/rhoi;
	G1overCsi = GAS_G1(simpar)/Csi;
	dAs = 0;
	float maxSignalVel=0;
	float maxwsij = 0;
	double ax,ay,az;
	ax = ay = az =0;
	for(i=0;i<Num_neighbor;i++){
		float F12,F2;
		TPtlStruct *pj;
		pj = bpneighbor[i];
		massj = pj->mass;
		rhoj = pj->rho;
		Ej = pj->Entropy;
		F2 = pj->F1;
		hsmlj = pj->hsml;
		dx = pj->x - pos->x;/* This should not be "p", because it is rescaled. */
		dy = pj->y - pos->y;
		dz = pj->z - pos->z;
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
			float fact12 = fact1+fact2;
			ax -= tmpx*fact12;
			ay -= tmpy*fact12;
			az -= tmpz*fact12;
	
			float v21x = dx*GAS_G1EXPANSION(simpar)+(pj->vx-p->vx)*G2overG1;
			float v21y = dy*GAS_G1EXPANSION(simpar)+(pj->vy-p->vy)*G2overG1;
			float v21z = dz*GAS_G1EXPANSION(simpar)+(pj->vz-p->vz)*G2overG1;
			/* wsij ~ 100 */
			float wsij = tmpx*v21x + tmpy*v21y+ tmpz*v21z;
			Csj = sqrtf(gamma*Tempjovermuj*kBovermH);
	
			float wrijoverCsi = G1overCsi*wsij;
			float SigVel = (1+Csj*invCsi-3*wrijoverCsi);
			float SigVelCsi = SigVel*Csi;
			if(SigVelCsi > maxSignalVel) maxSignalVel = SigVelCsi;
			if(wsij<0) {
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
	p->ax += ax;
	p->ay += ay;
	p->az += az;
	dAs = GAS_ENTROPYEVOLFACT(simpar)*dAs;
	p->dAs += Ei*(GAS_G1EXPANSION(simpar)*3*(1.-gamma)+dAs);

	if(maxSignalVel>0){
		float dai = GAS_COURANT(simpar)*hsmli/(maxSignalVel/GAS_G2(simpar))/NX(simpar);
		float ndiv = ASTEP(simpar)/dai;
		int isubpow = max(0,ceil(logf(ndiv)*invLog2));
		int imax = max(GetSphTsubPower(p),isubpow);
		if(imax>32) {
			DEBUGPRINT("Error in setting sph next step sub power:  %g %g : %d %d : %d\n",
				dai,ndiv,isubpow,imax,GetSphTsubPower(p));
			DEBUGPRINT("Error :  %g %g : %g %g \n",GAS_COURANT(simpar),hsmli,maxSignalVel,GAS_G2(simpar));
			DEBUGPRINT("Error :  Csi = %g T/mu= %g \n",Csi,Tempiovermui);
			exit(999);
		}
		else if(imax>16) {
			DEBUGPRINT("Error in setting sph next step sub power:  %ld : %g %g : %d %d : %d\n",
					PINDX(p), dai,ndiv,isubpow,imax,GetSphTsubPower(p));
			DEBUGPRINT("Error :  %g %g : %g %g \n",GAS_COURANT(simpar),hsmli,maxSignalVel,GAS_G2(simpar));
			DEBUGPRINT("Error :  Csi = %g T/mu= %g \n",Csi,Tempiovermui);
			exit(999);
		}
		else if(imax >= SPH_CONSTNEIGHPOW(simpar)){
		}

		SetSphNextTsubPower(p,imax);
		/*
		if(imax > GetTsubPower(p)) SetSphTsubPower(p,imax);
		*/
		if(imax > SPH_TIMESTEPLIMITER(simpar)) TimeStep_limiter(imax,bpneighbor,Num_neighbor);
	}
	if(GAS_FNFLAG(simpar) == 'Y'){ 
		/*
		maxwsij = 3*maxwsij*GAS_G1(simpar)/GAS_G2(simpar);
		float dai = simpar.sph.Durant*hsmli/(maxwsij*NX(simpar));
		float ndiv = ASTEP(simpar)/dai;
		int isubpow = max(0,ceilf(logf(ndiv)*invLog2));
		SetFixedPower(p,isubpow);
		if(isubpow > 20) 
			printf("PP%d has %g %g %g, %g %g %d %g %g %d\n",MYID(simpar),maxwsij,GAS_G1(simpar),GAS_G2(simpar),
				simpar.sph.Durant,hsmli,NX(simpar),ASTEP(simpar),dai,isubpow);
		*/
		void InsertFixedNeighborList(SimParameters *, treesphparticletype *, int, TPtlStruct **);
		if(GetSphTsubPower(p) > IndT_NSUBFIXED(simpar) ) InsertFixedNeighborList(simpar,p,Num_neighbor, bpneighbor);
	}


}

void Complete_GetDenAndOthers(SimParameters *simpar, particle *pos, TPtlStruct **bpneighbor,int Num_neighbor){
	int i;
	float tmpx,tmpy,tmpz,r,hsml;
	float tmpvx,tmpvy,tmpvz,vfact,divVel,curlVx,curlVy,curlVz,Cs;
	float density=0,density09=0,hsml09;
	divVel = curlVx = curlVy = curlVz = 0;
	hsml= LinkedTreeSphType(pos)->hsml;
	hsml09= hsml*0.99;
	for(i=0;i<Num_neighbor;i++){
		tmpx = bpneighbor[i]->x - pos->x;
		tmpy = bpneighbor[i]->y - pos->y;
		tmpz = bpneighbor[i]->z - pos->z;
		r = tmpx*tmpx; r += tmpz*tmpz; r += tmpy*tmpy;
		r = sqrtf(r);
		density += rho_w(r,hsml);
		density09 += rho_w(r,hsml09);

		tmpvx = ANOW(simpar)*(LinkedTreeSphType(bpneighbor[i])->vx - LinkedTreeSphType(pos)->vx) + tmpx*GAS_G1EXPANSION(simpar);
		tmpvy = ANOW(simpar)*(LinkedTreeSphType(bpneighbor[i])->vy - LinkedTreeSphType(pos)->vy) + tmpy*GAS_G1EXPANSION(simpar);
		tmpvz = ANOW(simpar)*(LinkedTreeSphType(bpneighbor[i])->vz - LinkedTreeSphType(pos)->vz) + tmpz*GAS_G1EXPANSION(simpar);
		if(r>0){
			vfact = bpneighbor[i]->mass *(-drhodr_w(r,hsml))/r;
			divVel += vfact * (tmpx*tmpvx + tmpy*tmpvy + tmpz*tmpvz);
			curlVx += vfact * (tmpvy * tmpz - tmpvz * tmpy);
			curlVy += vfact * (tmpvz * tmpx - tmpvx * tmpz);
			curlVz += vfact * (tmpvx * tmpy - tmpvy * tmpx);
		}
	}
	LinkedTreeSphType(pos)->rho = density;
	LinkedTreeSphType(pos)->fi =  1./(1+1./3.*(log(density)-log(density09))/(log(hsml)-log(hsml09)));

	float curlVel = sqrtf(curlVx*curlVx + curlVy*curlVy + curlVz*curlVz);
	/*
	Cs = sqrt(LinkedTreeSphType(pos)->Temp*GAS_GAMMA(simpar)/GetMu(LinkedTreeSphType(pos))*kBovermH);
	*/
	Cs = sqrt(LinkedTreeSphType(pos)->Entropy*GAS_GAMMA(simpar)*pow(density,GAS_GAMMA(simpar)-1));
	divVel = fabsf(divVel)/density;
	curlVel = fabsf(curlVel)/density;
	LinkedTreeSphType(pos)->F1 = divVel/(divVel + curlVel + 0.0001*Cs/hsml/GAS_G1(simpar));
}
void Complete_AllDetermineSph(SimParameters *simpar,
		particle *pos,int ngrv, TPtlStruct *grv, int Num_neighbor, TPtlStruct **bpneighbor){
	int i;
	float hsml=-2.e10;
	for(i=0;i<ngrv;i++){
		float tmpx, tmpy,tmpz,r;
		tmpx = grv[i].x - pos->x;
		tmpy = grv[i].y - pos->y;
		tmpz = grv[i].z - pos->z;
		r = sqrtf(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);
		hsml = max(hsml,r);
	}
	if(hsml ==0) hsml = BASICCELL_CELLWIDTH(simpar);
	LinkedTreeSphType(pos)->hsml = hsml;
	for(i=0;i<ngrv;i++) bpneighbor[i] = grv + i;
	CalculateAllSph(simpar, pos,bpneighbor,ngrv);
}
int Complete_SN_SW(SimParameters *simpar, particle *, TPtlStruct **,int ,double *);
void Complete_All_SN_SW(SimParameters *simpar, particle *pos,int ngrv, TPtlStruct *grv, TPtlStruct **bpneighbor, 
		double *dT){
	int i;
	/*
	float hsml=-2.e10;
	for(i=0;i<ngrv;i++){
		float tmpx, tmpy,tmpz,r;
		tmpx = grv[i].r[0] - pos->x;
		tmpy = grv[i].r[1] - pos->y;
		tmpz = grv[i].r[2] - pos->z;
		r = sqrtf(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);
		hsml = max(hsml,r);
	}
	LinkedTreeSphType(pos)->hsml = hsml;
	*/
	for(i=0;i<ngrv;i++) bpneighbor[i] = grv + i;
	Complete_SN_SW(simpar, pos,bpneighbor,ngrv,dT);
}

void Complete_DetermineSph(SimParameters *simpar, particle *pos, PosType neardist, TPtlStruct **bpneighbor, int Num_neighbor){
	LinkedTreeSphType(pos)->hsml = neardist;
	CalculateAllSph(simpar, pos,bpneighbor,Num_neighbor);
}

long SPHDump2PeriodicPosExceptZ(long, long, long,long,long,long,long,long,long,
		int,int,int, int, particle *,double,double,double,linkedlisttype**, long *,
		short int *);

long complete_SPHDump2Grv(SimParameters *simpar, size_t ix0, size_t iy0, size_t iz0,
		TPtlStruct *grv){
	size_t ix,iy,iz;
	PosType cellwidth = BASICCELL_CELLWIDTH(simpar);
	HydroTreeLinkedCell *BasicCell = SPH_BASICCELL(simpar);
	PosType xtran, ytran, ztran;
	xtran = -(ix0+0.5)*cellwidth;
	ytran = -(iy0+0.5)*cellwidth;
	ztran = -(iz0+0.5)*cellwidth;
	size_t ngrv = 0;

	linkedlisttype *p;
	int imin,jmin,kmin,imax,jmax,kmax;
	imin = MAX(ix0-1, 0);
	jmin = MAX(iy0-1, 0);
	kmin = MAX(iz0-1, 0);
	imax = MIN(ix0+2, BASICCELL_MX(simpar));
	jmax = MIN(iy0+2, BASICCELL_MY(simpar));
	kmax = MIN(iz0+2, BASICCELL_MZ(simpar));
	for(iz=kmin;iz<kmax;iz++){
        for(iy=jmin;iy<jmax;iy++){ 
            for(ix=imin;ix<imax;ix++){ 
				size_t mcell = ix*BASICCELL_MX(simpar)*(iy+BASICCELL_MY(simpar)*iz);
				p=BasicCell[mcell].link;
                while(p){ 
					treesphparticletype *tp = (treesphparticletype*)p;
                    grv[ngrv].type = TYPE_PTL; 
					if(IS_FLAG(p,BoundaryGhostflag)){
	                    grv[ngrv].x = p->x+xtran; 
						grv[ngrv].y = p->y+ytran;
						grv[ngrv].z = p->z+ztran; 
					}
					else {
	                    grv[ngrv].x = XofP(simpar, p)+xtran; 
						grv[ngrv].y = YofP(simpar, p)+ytran;
						grv[ngrv].z = ZofP(simpar, p)+ztran; 
					}
					grv[ngrv].mass = tp->mass;
					grv[ngrv].Entropy = tp->Entropy;
					grv[ngrv].rho = tp->rho;
					grv[ngrv].F1 = tp->F1;
					grv[ngrv].fi = tp->fi;
					grv[ngrv].hsml = tp->hsml;
					grv[ngrv].metallicity = tp->metallicity;
					grv[ngrv].vx = tp->vx;
					grv[ngrv].vy = tp->vy;
					grv[ngrv].vz = tp->vz;
					grv[ngrv].bp = p;
                    p = p->next; 
                    ngrv++; 
                } 
            } 
        } 
    } 
	return ngrv;
}
size_t STARDump2Pos(SimParameters *simpar, size_t ix0,size_t iy0,size_t iz0,
		particle *pos){
	size_t mcell = ix0+BASICCELL_MX(simpar)*(iy0+BASICCELL_MY(simpar)*iz0);
	PosType cellwidth = BASICCELL_CELLWIDTH(simpar);
	HydroTreeLinkedCell *BasicCell = STAR_BASICCELL(simpar);
	PosType xtran,ytran,ztran;
	xtran = -(ix0+0.5)*cellwidth;
	ytran = -(iy0+0.5)*cellwidth;
	ztran = -(iz0+0.5)*cellwidth;
	size_t npos = 0;
	linkedlisttype *p=BasicCell[mcell].link; 
	while(p){
#ifdef INDTIME
		if(isnowdmstep(IndT_NOWTSUBDIV(simpar),GetTsubPower(p),IndT_MAXTSUBPOWER(simpar))&&
				!IS_FLAG(p,BoundaryGhostflag))
#endif
		{
			pos[npos].x = XofP(simpar,p) + xtran;
			pos[npos].y = YofP(simpar,p) + ytran;
			pos[npos].z = ZofP(simpar,p) + ztran;
			pos[npos].bp = p;
			npos ++;
		}
		p = p->next;
	}
	return npos;
}
void UpdateTimeLimiter4Ghosts(SimParameters *simpar){
	size_t i,j;
	treesphparticletype *ghost = (treesphparticletype*)malloc(sizeof(treesphparticletype) *
			SPH_NPP(simpar));
	for(i=0;i<SPH_NPP(simpar);i++){
		ghost[i] = *(SPH_TBP(simpar)+i);
	}
	ptrdiff_t np = SPH_NPP(simpar);
	for(i=0;i<np;i++){
		(ghost+i)->x = fmod(XofP(simpar,ghost+i)+COSBOXSIZE(simpar), COSBOXSIZE(simpar));
		(ghost+i)->y = fmod(YofP(simpar,ghost+i)+COSBOXSIZE(simpar), COSBOXSIZE(simpar));
		(ghost+i)->z = fmod(ZofP(simpar,ghost+i)+COSBOXSIZE(simpar), COSBOXSIZE(simpar));
	}
	pmigrate( (void*)(&ghost), &np, SPH_DDINFO(simpar), &GRIDINFO(simpar) );
	for(i=0;i<np;i++){
		treesphparticletype *bp = ghost[i].bp;
		bp->mass = ghost[i].mass;
		bp->metallicity = ghost[i].metallicity;
		bp->dAs= ghost[i].dAs;
		for(j=0;j<4;j++) bp->indt[j] = ghost[i].indt[j];
	}
	free(ghost);
}
size_t com_SPHDump2Grv(SimParameters *simpar, size_t ix0,size_t iy0,size_t iz0,
		        TPtlStruct *grv) { 
	size_t ix,iy,iz,mcell; 
	PosType cellwidth = BASICCELL_CELLWIDTH(simpar); 
	HydroTreeLinkedCell *BasicCell = SPH_BASICCELL(simpar); 
	PosType xtran,ytran,ztran; 
	xtran = -(ix0+0.5)*cellwidth; 
	ytran = -(iy0+0.5)*cellwidth; 
	ztran = -(iz0+0.5)*cellwidth; 
	size_t ngrv = 0; 
	for(iz=iz0-1;iz<=iz0+1;iz++){ 
		for(iy=iy0-1;iy<=iy0+1;iy++){ 
			for(ix=ix0-1;ix<=ix0+1;ix++){ 
				mcell = ix+BASICCELL_MX(simpar)*(iy+BASICCELL_MY(simpar)*iz); 
				linkedlisttype *p=BasicCell[mcell].link; 
				while(p){ 
					grv[ngrv].type = TYPE_PTL; 
					grv[ngrv].x = XofP(simpar,p)+xtran; 
					grv[ngrv].y = YofP(simpar,p)+ytran; 
					grv[ngrv].z = ZofP(simpar,p)+ztran; 
					grv[ngrv].mass = p->mass; 
					grv[ngrv].bp = p; 

					treesphparticletype *tp = (treesphparticletype*)p;

					grv[ngrv].Entropy = tp->Entropy;
					grv[ngrv].rho = tp->rho;
					grv[ngrv].F1 = tp->F1;
					grv[ngrv].fi = tp->fi;
					grv[ngrv].hsml = tp->hsml;
					grv[ngrv].metallicity = tp->metallicity;
					grv[ngrv].vx = tp->vx;
					grv[ngrv].vy = tp->vy;
					grv[ngrv].vz = tp->vz;
					grv[ngrv].bp = p;
                    p = p->next; 
                    ngrv++; 
				} 
			} 
		} 
	} 
	return ngrv;
}

float CompleteSph(SimParameters *simpar,int Numnear, linkedlisttype **BasicCell, 
		linkedlisttype **StarBasicCell){
	size_t npmaxcell;
	float vfact2;
	float mass,zstart,zheight;
	ptrdiff_t nx,ny,nz,nspace;
	size_t i,j,k;
	size_t ix0,iy0,iz0;
	size_t ix2,iy2,iz2;
	size_t ncell,mcell,mp;
	size_t nsphpos,nstarpos,ngrv;
	size_t mpos,mgrv;
	PosType xtran,ytran,ztran;
	/* Time duration in terms of Myr */
	double dT[MaxEvolArrSize];


	TIMER_START(2);

	TStruct *TREECELL;
	linkedlisttype *p;
	size_t rnmax;
	if(MYID(simpar)==0){
		printf("###########################################\n");
		printf("SPH: measuring force\n");
		printf("###########################################\n\n\n");
	}
	mass = NSPACE(simpar)*NSPACE(simpar)*NSPACE(simpar);
	nx = NX(simpar);
	ny = NY(simpar);
	nz = NZ(simpar);
	nspace = NSPACE(simpar);

	PosType rspheresq = rsphere*rsphere;

	PosType cellwidth = BASICCELL_CELLWIDTH(simpar);
	PosType xmin,ymin,zmin,xmax,ymax, zmax;
	xmin = SIM_LXMIN(simpar,sph) - cellwidth; 
	ymin = SIM_LYMIN(simpar,sph) - cellwidth; 
	zmin = SIM_LZMIN(simpar,sph) - cellwidth; 
	xmax = SIM_LXMAX(simpar,sph) + cellwidth; 
	ymax = SIM_LYMAX(simpar,sph) + cellwidth; 
	zmax = SIM_LZMAX(simpar,sph) + cellwidth;
	size_t mx = (BASICCELL_MX(simpar) = ceil( (xmax-xmin)/cellwidth)); 
	size_t my = (BASICCELL_MY(simpar) = ceil( (ymax-ymin)/cellwidth)); 
	size_t mz = (BASICCELL_MZ(simpar) = ceil( (zmax-zmin)/cellwidth));


	if(CONT_HALFSTEP(simpar) == KICK) {
		for(i=0;i<MaxEvolArrSize;i++)  dT[i] = dTKICK[i] /Myr;
	}
	else if(CONT_HALFSTEP(simpar) == HALFPUSH) {
		for(i=0;i<MaxEvolArrSize;i++) dT[i] = dTPUSH[i] /Myr;
	}
	else {
		for(i=0;i<MaxEvolArrSize;i++) dT[i] = dTPULL[i] /Myr;
	}


	/* This is to remember the original particle before padding */
	for(i=0;i<SPH_NP(simpar);i++) (SPH_TBP(simpar)+i)->bp = (void*)(SPH_TBP(simpar)+i);
	ComPaddingHydroParticles(simpar);
	flagboundaryghostparticles(simpar);

	SPH_BASICCELL(simpar) = (HydroTreeLinkedCell*)malloc(sizeof(HydroTreeLinkedCell)*mx*my*mz);
	npmaxcell = SPHBuildLinkedList(simpar);

	int nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel
	{
		if(omp_get_thread_num()==0) nthreads = omp_get_num_threads();
	}
#endif
	particle *pos = (particle*)malloc(sizeof(particle)*npmaxcell*nthreads);
	TStruct *TREE  = (TStruct*)malloc(sizeof(TStruct)*npmaxcell*nthreads*27/3);
	TPtlStruct *grv  = (TPtlStruct*)malloc(sizeof(TPtlStruct)*npmaxcell*nthreads*27);
	TPtlStruct **bpneighbor = (TPtlStruct **)malloc(sizeof(TPtlStruct*)*Numnear*npmaxcell*omp_get_num_threads());



	size_t ii;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 2)
#endif
	for(ii=0;ii<mx*my*mz;ii++){
		size_t iz = ii/(mx*my); if(iz==mz-1 || iz==0) continue;
		size_t iy = (ii%(mx*ny))/my; if(iy==my-1 || iy==0) continue;
		size_t ix = (ii%mx); if(ix==mx-1 || ix==0) continue;
		int idthread = omp_get_thread_num();
		particle *threadpos = pos + npmaxcell*idthread;
		TStruct *TREECELL = TREE + npmaxcell*idthread * (27/3);
		TPtlStruct *threadgrv = grv + npmaxcell*idthread * 27;
		TPtlStruct **threadbpneighbor =  bpneighbor + npmaxcell *Numnear*idthread;
		size_t nnode = npmaxcell*27;
		size_t npos = SPHDump2Pos(simpar, ix,iy,iz,threadpos);
		if(npos>0) {
			size_t ngrv = com_SPHDump2Grv(simpar, ix,iy,iz,threadgrv);
			if(ngrv>0) {
				if(ngrv >300){
					Make_NN_Tree(TREECELL, nnode, threadgrv, ngrv, RECURSIVE);
					for(i=0;i<npos;i++){
						int numnear;
						PosType neardist;
						numnear = Find_Near(threadpos+i, Numnear, TREECELL, threadgrv,
								&neardist, threadbpneighbor);
						Complete_DetermineSph(simpar, threadpos+i, neardist, threadbpneighbor, numnear);
					}
				}
				else if(ngrv > Numnear){
					for(i=0;i<npos;i++){
						PosType neardist;
						int numnear = Direct_Find_Near(threadpos+i, Numnear, threadgrv, ngrv, &neardist,
								threadbpneighbor);
						Complete_DetermineSph(simpar, threadpos+i, neardist, threadbpneighbor, numnear);
					}
				}
				else {
					if(npos>0){
						for(i=0;i<npos;i++) Complete_AllDetermineSph(simpar,
								threadpos+i, ngrv, threadgrv, Numnear, threadbpneighbor);
					}
				}
			}
		}

	}

	UpdateTimeLimiter4Ghosts(simpar);
	free(bpneighbor); 
	free(TREE); 
	free(grv); 
	free(pos);


	STAR_BASICCELL(simpar) = (HydroTreeLinkedCell*)malloc(sizeof(HydroTreeLinkedCell)*mx*my*mz);
	size_t starnpmaxcell = STARBuildLinkedList(simpar);
	npmaxcell = MAX(npmaxcell, starnpmaxcell);

	pos = (particle*)malloc(sizeof(particle)*npmaxcell*nthreads);
	TREE  = (TStruct*)malloc(sizeof(TStruct)*npmaxcell*nthreads*27/3);
	grv  = (TPtlStruct*)malloc(sizeof(TPtlStruct)*npmaxcell*nthreads*27);


	ptrdiff_t nsn=0;
	bpneighbor = (TPtlStruct **)malloc(
			sizeof(TPtlStruct*)*Numnear*npmaxcell*omp_get_num_threads());
	if(GAS_SNFBFLAG(simpar) == 'Y'){
		if(MYID(simpar)==0){
			printf("Now entering into SN and SW routine \n");
		}
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 2)
#endif
		for(ii=0;ii<mx*my*mz;ii++){
			size_t iz = ii/(mx*my); if(iz==mz-1 || iz==0) continue;
			size_t iy = (ii%(mx*ny))/my; if(iy==my-1 || iy==0) continue;
			size_t ix = (ii%mx); if(ix==mx-1 || ix==0) continue;
			int idthread = omp_get_thread_num();
			particle *threadpos = pos + npmaxcell*idthread;
			TStruct *TREECELL = TREE + npmaxcell*idthread * (27/3);
			TPtlStruct *threadgrv = grv + npmaxcell*idthread * 27;
			TPtlStruct **threadbpneighbor = bpneighbor + npmaxcell *Numnear*idthread;
			size_t nnode = npmaxcell*27;
			size_t npos = STARDump2Pos(simpar, ix,iy,iz,threadpos);
			if(npos>0) {
				size_t ngrv = com_SPHDump2Grv(simpar, ix,iy,iz,threadgrv);
				if(ngrv>0){
					if(ngrv>Numnear){
						int numnear;
						PosType neardist;
						Make_NN_Tree(TREECELL, nnode, threadgrv, ngrv, RECURSIVE);
						for(i=0;i<npos;i++){
							numnear = Find_Near(threadpos+i, Numnear, TREECELL, threadgrv,
									&neardist, threadbpneighbor);
							nsn += Complete_SN_SW(simpar, threadpos+i, threadbpneighbor,numnear,dT);
						}
					}
					else {
						for(i=0;i<npos;i++) 
							Complete_All_SN_SW(simpar, threadpos+i,ngrv,threadgrv,threadbpneighbor,dT); 
					}
				}
			}
		}
		UpdateTimeLimiter4Ghosts(simpar);
	}
	free(grv);free(bpneighbor);free(TREE);free(pos);
	ptrdiff_t tnsn;
	MPI_Reduce(&nsn,&tnsn,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	if(MYID(simpar)==0) {
		if(tnsn >0){
			float time;
			if(BGEXPAND(simpar)=='Y') time = AMAX(simpar)/ANOW(simpar) -1;
			else time = ANOW(simpar);
			printf("Total Number of SuperNovae explosion %ld at z/t = %g\n",tnsn,time);
		}
	}
	TIMER_STOP(2);
	float maintreetime = ELAPSED_TIME(2);
	return maintreetime;
}
