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
#include "indT.h"
#include "hydroBasicCell.h"
#include "nnost.h"
#include "sph.h"

#ifdef USE_GPU
#include "nvgpu.h"
void gputreeforce(particle *,int, TStruct *,int , TPtlStruct *,int, float, int, int,int);
void Direct_Nbody(particle *,int, TPtlStruct *,int , int, int,int);
#endif

#define MAX(a,b) ( (a)>(b) ? (a):(b) )
#define MIN(a,b) ( (a)<(b) ? (a):(b) )


#define invLog2 1.4426950408






#ifdef NO_TREE_WALK
#define flag_TW (0)
#else
#define flag_TW (1)
#endif

#ifndef _OPENMP
#define omp_get_thread_num() nullfct0()
#define omp_get_num_threads() nullfct1()
#endif

#define p8overpi (2.546479089)
#define p16overpi (5.0929581789)
#define p48overpi (15.2788745)
#define m48overpi (-15.2788745)
#define m24overpi (-7.639437268)


float rho_w(float r, float h){
    float hinv = 1./h;
    float rhinv = r*hinv;
    float result;
    if(rhinv < 0.5) {
        float rhinv2 = rhinv*rhinv;
        float frac = 1.+6*rhinv2*(rhinv-1.);
        float hinv3 = hinv*hinv*hinv;
        result= p8overpi*hinv3*frac;
    }
    else if(rhinv<1.) {
        float frac = (1.-rhinv);
        float frac3 = frac*frac*frac;
        float hinv3 = hinv*hinv*hinv;
        result= p16overpi*hinv3*frac3;
    }
    else {
        result = 0;
    }
    return result;
}
float drhodr_w(float r,float h){
    float hinv = 1./h;
    float rhinv = r*hinv;
    float result;
    if(rhinv < 0.5) {
        float hinv2 = hinv*hinv;
        float hinv4 = hinv2*hinv2;
        result = p8overpi*hinv4*rhinv*(18.*rhinv-12.);
    }
    else if(rhinv<1.) {
        float hinv2 = hinv*hinv;
        float hinv4 = hinv2*hinv2;
        float fact = (rhinv-1.);
        result= m48overpi*hinv4*fact*fact;
    }
    else result = 0;
    return result;
}
/* This function calculates the density (rho), and fi */
void GetDenAndOthers(SimParameters *simpar, particle *pos, TPtlStruct **bpneighbor,int Num_neighbor){
	int i;
	float tmpx,tmpy,tmpz,r,hsml;
	float tmpvx,tmpvy,tmpvz,vfact,divVel,curlVx,curlVy,curlVz,Cs;
	float massj,massi,density=0,density09=0,hsml09;
	divVel = curlVx = curlVy = curlVz = 0;
	hsml= LinkedTreeSphType(pos)->hsml;
	hsml09 = hsml*0.99;
	treesphparticletype *p = LinkedTreeSphType(pos);
	float maxwsij = 0;
	for(i=0;i<Num_neighbor;i++){
		tmpx = bpneighbor[i]->x - pos->x;
		tmpy = bpneighbor[i]->y - pos->y;
		tmpz = bpneighbor[i]->z - pos->z;
		r = tmpx*tmpx; r += tmpz*tmpz; r += tmpy*tmpy;
		r = sqrtf(r);
		density += rho_w(r,hsml)*bpneighbor[i]->mass;
		density09 += rho_w(r,hsml09)*bpneighbor[i]->mass;


		float dvx = (LinkedTreeSphType(bpneighbor[i])->vx - p->vx);
		float dvy = (LinkedTreeSphType(bpneighbor[i])->vy - p->vy);
		float dvz = (LinkedTreeSphType(bpneighbor[i])->vz - p->vz);

		tmpvx = GAS_G12RATIO(simpar)*dvx + tmpx*GAS_G1EXPANSION(simpar);
		tmpvy = GAS_G12RATIO(simpar)*dvy + tmpy*GAS_G1EXPANSION(simpar);
		tmpvz = GAS_G12RATIO(simpar)*dvz + tmpz*GAS_G1EXPANSION(simpar);


		if(r>0){
			/* These two lines are inserted for the Fixed Neighbor scheme. */
			float wsij = fabs((tmpx*dvx + tmpy*dvy + tmpz*dvz)/r);
			if(wsij > maxwsij) maxwsij = wsij;

			vfact = bpneighbor[i]->mass *(-drhodr_w(r,hsml))/r;
			divVel += vfact * (tmpx*tmpvx + tmpy*tmpvy + tmpz*tmpvz);
			curlVx += vfact * (tmpvy * tmpz - tmpvz * tmpy);
			curlVy += vfact * (tmpvz * tmpx - tmpvx * tmpz);
			curlVz += vfact * (tmpvx * tmpy - tmpvy * tmpx);
		}
	}
	if(GAS_FNFLAG(simpar) == 'Y'){
		float dai = GAS_DURANT(simpar)*hsml/(maxwsij*NX(simpar));
		float ndiv = ASTEP(simpar)/dai;
		int isubpow = max(0,ceilf(logf(ndiv)*invLog2));
		SetFixedPower(p,isubpow);
	}

	LinkedTreeSphType(pos)->rho = density;
	LinkedTreeSphType(pos)->fi = 1./(1+1./3.*(log(density)-log(density09))/(log(hsml)-log(hsml09)));
	LinkedTreeSphType(pos)->divVel = divVel;

	float curlVel = sqrtf(curlVx*curlVx + curlVy*curlVy + curlVz*curlVz);
	Cs = sqrt(LinkedTreeSphType(pos)->Entropy*GAS_GAMMA(simpar)*pow(LinkedTreeSphType(pos)->rho,GAS_GAMMA(simpar)-1));
	divVel = fabsf(divVel)/density;
	curlVel = fabsf(curlVel)/density;
	LinkedTreeSphType(pos)->F1 = divVel/(divVel + curlVel + 0.0001*Cs/hsml/GAS_G1(simpar));
}

void AllDetermineSphDen(SimParameters *simpar,particle *pos,int ngrv, TPtlStruct *grv, int Num_neighbor, TPtlStruct **bpneighbor){
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
	if(hsml ==0) hsml = 1;
	LinkedTreeSphType(pos)->hsml = hsml;
	for(i=0;i<ngrv;i++) bpneighbor[i] = grv + i;
	GetDenAndOthers(simpar, pos,bpneighbor,ngrv);
}
void DetermineSphDen(SimParameters *simpar,
		particle *pos, float maxdistofneigh, TPtlStruct **bpneighbor, int Num_neighbor){
	LinkedTreeSphType(pos)->hsml = maxdistofneigh;
	GetDenAndOthers(simpar, pos,bpneighbor,Num_neighbor);
}
void GetDenAndOthers1(SimParameters *simpar, particle *pos, TPtlStruct **bpneighbor,int Num_neighbor){
	int i;
	float tmpx,tmpy,tmpz,r,hsml;
	float tmpvx,tmpvy,tmpvz,vfact,divVel,curlVx,curlVy,curlVz,Cs;
	float density=0,density09=0,hsml09;
	divVel = curlVx = curlVy = curlVz = 0;
	hsml= LinkedTreeSphType(pos)->hsml;
	hsml09 = hsml*0.99;
	for(i=0;i<Num_neighbor;i++){
		tmpx = bpneighbor[i]->x - pos->x;
		tmpy = bpneighbor[i]->y - pos->y;
		tmpz = bpneighbor[i]->z - pos->z;
		printf("s%d tmpx/y/z = %g %g %g\n",i,tmpx,tmpy,tmpz);
		r = tmpx*tmpx; r += tmpz*tmpz; r += tmpy*tmpy;
		r = sqrtf(r);
		printf("s%d r = %g \n",i,r);
		density += rho_w(r,hsml);
		density09 += rho_w(r,hsml09);
		printf("s%d density = %g \n",i,density);

		tmpvx = GAS_G12RATIO(simpar)*(LinkedTreeSphType(bpneighbor[i])->vx - LinkedTreeSphType(pos)->vx) + tmpx*GAS_G1EXPANSION(simpar);
		tmpvy = GAS_G12RATIO(simpar)*(LinkedTreeSphType(bpneighbor[i])->vy - LinkedTreeSphType(pos)->vy) + tmpy*GAS_G1EXPANSION(simpar);
		tmpvz = GAS_G12RATIO(simpar)*(LinkedTreeSphType(bpneighbor[i])->vz - LinkedTreeSphType(pos)->vz) + tmpz*GAS_G1EXPANSION(simpar);
		printf("s%d tmpvx/vy/vz = %g %g %g \n",i,tmpvx,tmpvy,tmpvz);
		if(r>0){
			vfact = bpneighbor[i]->mass *(-drhodr_w(r,hsml))/r;
			printf("s%d vfact/hsml = %g %g \n",i,vfact,hsml);
			divVel += vfact * (tmpx*tmpvx + tmpy*tmpvy + tmpz*tmpvz);
			curlVx += vfact * (tmpvy * tmpz - tmpvz * tmpy);
			curlVy += vfact * (tmpvz * tmpx - tmpvx * tmpz);
			curlVz += vfact * (tmpvx * tmpy - tmpvy * tmpx);
		}
	}
	LinkedTreeSphType(pos)->rho = density;
	LinkedTreeSphType(pos)->fi = 1./(1+1./3.*(log(density)-log(density09))/(log(hsml)-log(hsml09)));

	float curlVel = sqrtf(curlVx*curlVx + curlVy*curlVy + curlVz*curlVz);
	/*
	Cs = sqrt(LinkedTreeSphType(pos)->Temp/GetMu(LinkedTreeSphType(pos))*kBovermH);
	*/
	Cs = sqrt(LinkedTreeSphType(pos)->Entropy*GAS_GAMMA(simpar)*pow(LinkedTreeSphType(pos)->rho,GAS_GAMMA(simpar)-1));
	printf("Cs= %g\n",Cs);
	divVel = fabsf(divVel)/density;
	curlVel = fabsf(curlVel)/density;
	printf("divVel/curlVel= %g %g\n",divVel,curlVel);
	LinkedTreeSphType(pos)->F1 = divVel/(divVel + curlVel + 0.0001*Cs/hsml/GAS_G1(simpar));
	/*
	LinkedTreeSphType(pos)->fi = 1.;
	*/
	/*
	printf("density %g Cs= %g divVel/curlVel = %g %g fi= %g ::::: gamma/g1/kBovermHi %g %g %g: Entropy= %g \n",
			density,Cs,divVel,curlVel,LinkedTreeSphType(pos)->fi, GAS_GAMMA(simpar),GAS_G1(simpar),kBovermH,
			LinkedTreeSphType(pos)->Entropy);
	MPI_Finalize();exit(999);
	*/
}
void DetermineSphDen1(SimParameters *simpar,particle *pos, float maxdistofneigh, TPtlStruct **bpneighbor, int Num_neighbor){
	LinkedTreeSphType(pos)->hsml = maxdistofneigh;
	GetDenAndOthers1(simpar, pos,bpneighbor,Num_neighbor);
}

size_t SPHBuildLinkedList(SimParameters *simpar){
	size_t mx,my,mz;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	mz = BASICCELL_MZ(simpar);
	PosType xmin,ymin,zmin;
	xmin = SIM_LXMIN(simpar, sph);
	ymin = SIM_LYMIN(simpar, sph);
	zmin = SIM_LZMIN(simpar, sph);
	PosType CellWidth = BASICCELL_CELLWIDTH(simpar);

	HydroTreeLinkedCell *SPH_BasicCell = SPH_BASICCELL(simpar);
#ifdef NEEDLESS
	HydroTreeLinkedCell *OTHERS_BasicCell = OTHERS_BASICCELL(simpar);
#endif

	ptrdiff_t i;
	for(i=0;i<mx*my*mz;i++){
		SPH_BasicCell[i].link = NULL;
		SPH_BasicCell[i].nmem = 0;
		SPH_BasicCell[i].calflag = 0;
#ifdef NEEDLESS
		OTHERS_BasicCell[i].link = NULL;
		OTHERS_BasicCell[i].nmem = 0;
		OTHERS_BasicCell[i].calflag = 0;
#endif
	}
	LinkParticles(simpar, SPH, sph, SPH_BasicCell, mx,my,xmin,ymin,zmin, CellWidth);
#ifdef NEEDLESS
	LinkParticles(simpar, STAR, star, OTHERS_BasicCell, mx,my,xmin,ymin,zmin, CellWidth);
	LinkParticles(simpar, AGN, agn, OTHERS_BasicCell, mx,my,xmin,ymin,zmin, CellWidth);
#endif
	size_t npmaxcell = -1;
	for(i=0;i<mx*my*mz;i++){
		npmaxcell = MAX(npmaxcell, SPH_BasicCell[i].nmem);
	}
	return npmaxcell;
}
size_t SPHDump2Pos(SimParameters *simpar, size_t ix0,size_t iy0,size_t iz0,
		particle *pos){
	size_t mcell = ix0+BASICCELL_MX(simpar)*(iy0+BASICCELL_MY(simpar)*iz0);
	PosType cellwidth = BASICCELL_CELLWIDTH(simpar);
	HydroTreeLinkedCell *BasicCell = SPH_BASICCELL(simpar);
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

size_t SPHDump2Grv(SimParameters *simpar, size_t ix0,size_t iy0,size_t iz0,
		TPtlStruct *grv) {
	size_t ix,iy,iz,mcell;
	PosType cellwidth = BASICCELL_CELLWIDTH(simpar);
	HydroTreeLinkedCell *BasicCell = SPH_BASICCELL(simpar);
	PosType xtran,ytran,ztran;
	xtran = -(ix0+0.5)*cellwidth;
	ytran = -(iy0+0.5)*cellwidth;
	ztran = -(iz0+0.5)*cellwidth;

	size_t ngrv = 0;
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
				mcell = ix+BASICCELL_MX(simpar)*(iy+BASICCELL_MY(simpar)*iz);
				linkedlisttype *p=BasicCell[mcell].link; 
                while(p){ 
					if(IS_FLAG(p,BoundaryGhostflag)){
						grv[ngrv].x = p->x+xtran; 
						grv[ngrv].y = p->y+ytran; 
						grv[ngrv].z = p->z+ztran; 
					}
					else { 
						grv[ngrv].x = XofP(simpar,p)+xtran; 
						grv[ngrv].y = YofP(simpar,p)+ytran; 
						grv[ngrv].z = ZofP(simpar,p)+ztran; 
					}
                    grv[ngrv].type = TYPE_PTL; 
					grv[ngrv].mass = p->mass;
					grv[ngrv].bp = p;
                    p = p->next; 
                    ngrv++; 
                } 
            } 
        } 
    } 
	return ngrv;
}
#define FlagBoundaryGhost(simpar, TYPE, type) do{\
	tree##type##particletype *p = TYPE##_TBP(simpar);\
	for(i=0;i<TYPE##_NP(simpar);i++){\
		SET_FLAG(p, BoundaryGhostflag);\
		p++;\
	}\
}while(0)
void flagboundaryghostparticles(SimParameters *simpar){
	size_t i;
	FlagBoundaryGhost(simpar, SPH, sph);
	FlagBoundaryGhost(simpar, AGN, agn);
	FlagBoundaryGhost(simpar, STAR, star);
}
void FindSphDensity(SimParameters *simpar){
	int Numnear = SPH_NUMNEAR(simpar);
	size_t npmaxcell;
	HydroTreeLinkedCell *BasicCell = SPH_BASICCELL(simpar);
#ifdef NEEDLESS
	HydroTreeLinkedCell *StarBasicCell = OTHERS_BASICCELL(simpar);
#endif
	float vfact2;
	float mass,zstart,zheight;
	size_t nx,ny,nz,nspace;
	ptrdiff_t i,j,k;
	ptrdiff_t ix0,iy0,iz0;
	ptrdiff_t ix2,iy2,iz2;
	ptrdiff_t ncell,mcell,mp;
	ptrdiff_t npos,ngrv;
	size_t mpos,mgrv;
	double xtran,ytran,ztran;
#if defined (NO_TREE_WALK) || (TREE_WALK_BUT_CONTACT_LIST)
	long mbuff;
	float *buff;
#endif
	if(MYID(simpar)==0) {
		printf("###########################################\n");
		printf("SPH: measuring local density\n");
		printf("###########################################\n\n\n");
	}
	nx = NX(simpar);
	ny = NY(simpar);
	nz = NZ(simpar);
	nspace = NSPACE(simpar);

	PosType rspheresq = rsphere*rsphere;
	/*
	determine_mpi_long();
	*/
	int nid = NID(simpar);
	int myid = MYID(simpar);

	PosType cellwidth = BASICCELL_CELLWIDTH(simpar);


	PosType xmin,ymin,zmin,xmax,ymax,zmax;
	xmin = SIM_LXMIN(simpar,sph) - cellwidth;
	ymin = SIM_LYMIN(simpar,sph) - cellwidth;
	zmin = SIM_LZMIN(simpar,sph) - cellwidth;
	xmax = SIM_LXMAX(simpar,sph) + cellwidth;
	ymax = SIM_LYMAX(simpar,sph) + cellwidth;
	zmax = SIM_LZMAX(simpar,sph) + cellwidth;

	size_t mx = (BASICCELL_MX(simpar) = ceil( (xmax-xmin)/cellwidth));
	size_t my = (BASICCELL_MY(simpar) = ceil( (ymax-ymin)/cellwidth));
	size_t mz = (BASICCELL_MZ(simpar) = ceil( (zmax-zmin)/cellwidth));


	PaddingHydroParticles(simpar);

	flagboundaryghostparticles(simpar);

	SPH_BASICCELL(simpar) = (HydroTreeLinkedCell*)malloc(sizeof(HydroTreeLinkedCell)*mx*my*mz);
#ifdef NEEDLESS
	OTHERS_BASICCELL(simpar) = (HydroTreeLinkedCell*)malloc(sizeof(HydroTreeLinkedCell)*mx*my*mz);
#endif
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

	size_t ii;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 2)
#endif
	for(ii=0;ii<mx*my*mz;ii++){
		TPtlStruct *bpneighbor[MAX_NUM_NEAR];
		size_t iz = ii/(mx*my); if(iz==mz-1 || iz==0) continue;
		size_t iy = (ii%(mx*my))/my; if(iy==my-1 || iy==0) continue;
		size_t ix = ii%mx; if(ix == mx-1 || ix ==0) continue;
		int idthread = omp_get_thread_num();
		particle *threadpos = pos + npmaxcell*idthread;
		TStruct *TREECELL = TREE + npmaxcell*idthread * (27/3);
		TPtlStruct *threadgrv = grv + npmaxcell*idthread * 27;
		size_t nnode = npmaxcell*27;
		size_t npos = SPHDump2Pos(simpar, ix,iy,iz,threadpos);
		if(npos>0){
			size_t ngrv = SPHDump2Grv(simpar, ix,iy,iz,threadgrv);
			if(ngrv>0){
#ifdef USE_GPU
				if(npos>600) {
					Make_NN_Tree(TREECELL,nnode, threadgrv,ngrv, RECURSIVE);
					void GpuDeterminSphDen(particle *, long  , int , TStruct *, int ,TPtlStruct *, 
							int ,SimParameters , int , int );
					GpuDeterminSphDen(threadpos,npos,Numnear,TREECELL,
							threadgrv,ngrv, simpar, GPUSPERNODE,0);// 0 is density measurement flag 
				}
				else
#endif
				if(npos>=23)
				{
					Make_NN_Tree(TREECELL,nnode, threadgrv,ngrv, RECURSIVE);

					for(i=0;i<npos;i++){
						int numnear;
						PosType maxdistofneigh;
						numnear = Find_Near(threadpos+i,Numnear,TREECELL,threadgrv,
							&maxdistofneigh,bpneighbor);
						DetermineSphDen(simpar, threadpos+i,maxdistofneigh, bpneighbor,numnear);
					}
				}
				else if(ngrv> Numnear)
				{
					for(i=0;i<npos;i++){
						int numnear;
						PosType maxdistofneigh;
						numnear = Direct_Find_Near(threadpos+i,Numnear, threadgrv,ngrv,&maxdistofneigh,
								bpneighbor);
						DetermineSphDen(simpar, threadpos+i,maxdistofneigh, bpneighbor,numnear);
					}
				}
				else {
					for(i=0;i<npos;i++) 
						AllDetermineSphDen(simpar, threadpos+i,ngrv,threadgrv,Numnear,bpneighbor);

				}
			}
		}

	}
	free(grv);free(TREE);free(pos);

}
