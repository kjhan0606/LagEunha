/* It adopts the tidal radius of the elliptical orbits. May 2008 */
/* It implements more improved version of peak identification
 * including the merging of underpopulated peaks. May, 23, 2008 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <mpi.h>
#include <omp.h>
#include "Memory.h"
#include "ramses.h"
#include "tree.h"
#include "defs.h"
#include "params.h"
/*
#define free(a) Free(a)
*/

/*
#define PEAKTHRESHOLD 5.E2L
#define BOUNDITER 4
#define COREDENRESOLUTION (1.e1L)
#define MINCORENMEM 30
#define NSHELLDIVIDE 30
#define MAXNUMCORE 200000
*/


static	int maxnumcore = MAXNUMCORE;


#define BOUND 1
#define UNBOUND 0 
#define MAX_TIDAL_R 2.e5
extern int myid,nid;
/* ????Àº ???????? ??À» ??Á¤?Ñ´?. */
extern double onesolarmass;
extern double com2real,real2com,potentfact;
extern double pntmass;
extern double r1kineticfact,r2kineticfact;
extern float amax,a,rng,size,hubble;
extern int ng,nspace;
extern int nx,ny,nz;
extern float omep,omeplam,acoeff[8];
extern float epsilon;

double Xmin,Ymin,Zmin;


/*
#define COREIDFROMESCORE(i) ((int)(score[i].core-core))
*/
#define KP2BPID(kp,i) ((int)(kp[i].bp-bp))
#define SCOREID2COREID(i) ((int)(score[i].core-core))


int sorttype(const void *a,const void *b){
	float *aa,*bb;
	aa = (float*)a;
	bb = (float*)b;
	if(aa<bb) return -1;
	else if(aa>bb) return 1;
	else return 0;
}

int numstack;
#include "header.h"
int Motherrank = motherrank;
extern float m_tidal[NUM_MASS],r_tidal[NUM_MASS];
typedef struct Vector3d{
	double x,y,z;
} Vector3d;
typedef struct Kptype{
	double x,y,z;
	float vx,vy,vz,mass;
	SimpleBasicParticleType *bp;
}Kptype;

void mklocalize(SimpleBasicParticleType *bp,lint np,float *xinit,float *yinit,
		float *zinit,float *xmax,float *ymax,float *zmax){
	float min,max;
	SimpleBasicParticleType *tmpr;
	float *tmp,*tmparr;
	lint i,j;
	int *idx;
	float maxdist,tmpdist;
	int upbound;
	*xinit = size; *xmax = -size;
	*yinit = size; *ymax = -size;
	*zinit = size; *zmax = -size;
	for(i=0;i<np;i++){
		tmpr = bp+i;
		*xinit = MIN(*xinit,tmpr->x);
		*xmax  = MAX(*xmax ,tmpr->x);
		*yinit = MIN(*yinit,tmpr->y);
		*ymax  = MAX(*ymax ,tmpr->y);
		*zinit = MIN(*zinit,tmpr->z);
		*zmax  = MAX(*zmax ,tmpr->z);
	}
	if(floor(*xinit) <= 0+nblur && ceil(*xmax) >= size-nblur){
		tmparr = (float *) Malloc(sizeof(float)*np,PPTR(tmparr));
		tmp = tmparr;tmpr = bp;
		for(i=0;i<np;i++) *(tmp++) = (tmpr++)->x;
		qsort(tmparr,np,sizeof(float),sorttype);
		maxdist = -10.;
		for(i=0;i<np-1;i++){
			tmpdist = tmparr[i+1] - tmparr[i];
			if(tmpdist > maxdist) {
				upbound = i;
				maxdist = tmpdist;
			}
		}
		for(i=0;i<np;i++){
			if(bp[i].x <= tmparr[upbound]) bp[i].x += size;
		}
		/*
		for(i=0;i<=upbound;i++){
			bp[i].x += nx;
		}
		*/
		*xinit = bp[upbound+1].x;
		*xmax = bp[upbound].x;
		Free(tmparr);
	}
	if(floor(*yinit) <= 0+nblur && ceil(*ymax) >= size-nblur){
		tmparr = (float *) Malloc(sizeof(float)*np,PPTR(tmparr));
		tmp = tmparr;
		tmpr = bp;
		for(i=0;i<np;i++) *(tmp++) = tmpr++->y;
		qsort(tmparr,np,sizeof(float),sorttype);
		maxdist = -10.;
		for(i=0;i<np-1;i++){
			tmpdist = tmparr[i+1] - tmparr[i];
			if(tmpdist > maxdist) {
				upbound = i;
				maxdist = tmpdist;
			}
		}
		for(i=0;i<np;i++){
			if(bp[i].y <= tmparr[upbound]) bp[i].y += size;
		}
		/*
		for(i=0;i<=upbound;i++){
			bp[i].y += ny;
		}
		*/
		*yinit = bp[upbound+1].y;
		*ymax = bp[upbound].y;
		Free(tmparr);
	}
	if(floor(*zinit) <= 0+nblur && ceil(*zmax) >= size-nblur){
		tmparr = (float *) Malloc(sizeof(float)*np,PPTR(tmparr));
		tmp = tmparr;
		tmpr = bp;
		for(i=0;i<np;i++) *(tmp++) = tmpr++->z;
		qsort(tmparr,np,sizeof(float),sorttype);
		maxdist = -10.;
		for(i=0;i<np-1;i++){
			tmpdist = tmparr[i+1] - tmparr[i];
			if(tmpdist > maxdist) {
				upbound = i;
				maxdist = tmpdist;
			}
		}
		for(i=0;i<np;i++){
			if(bp[i].z <= tmparr[upbound]) bp[i].z += size;
		}
		/*
		for(i=0;i<=upbound;i++){
			bp[i].z += nz;
		}
		*/
		*zinit = bp[upbound+1].z;
		*zmax = bp[upbound].z;
		Free(tmparr);
	}
}

#define NOT 0 /* This is used for the flag of WorkingParticle */

#define WP_RESET 0x00
#define WP_PEAK  0x01
#define WP_VISITED 0x02
#define WP_CORE  0x04
#define WP_SHELL 0x08
#define WP_BOUND 0x10
#define WP_REMAINING 0x20


#define WP_VISITEDT 0x01
/*
#define WP_BOUND 0x40
#define WP_REST  0x80
*/

#define RESET_WHOLE_FLAGS_WP(x) (wp[x].flag  = WP_RESET)

#define SET_VISITED(x) (wp[x].flag |= WP_VISITED)
#define SET_PEAK(x) (wp[x].flag |= WP_PEAK)
#define SET_CORE(x) (wp[x].flag |= WP_CORE)
#define SET_SHELL(x) (wp[x].flag |= WP_SHELL)
#define SET_BOUND(x) (wp[x].flag |= WP_BOUND)
#define SET_REMAINING(x) (wp[x].flag |= WP_REMAINING)

#define UNSET_VISITED(x) (wp[x].flag &= (~WP_VISITED))
#define UNSET_PEAK(x) (wp[x].flag &= (~WP_PEAK))
#define UNSET_CORE(x) (wp[x].flag &= (~WP_CORE))
#define UNSET_SHELL(x) (wp[x].flag &= (~WP_SHELL))
#define UNSET_BOUND(x) (wp[x].flag &= (~WP_BOUND))
#define UNSET_REMAINING(x) (wp[x].flag &= (~WP_REMAINING))

#define IS_PEAK(x) (wp[x].flag & WP_PEAK)
#define IS_VISITED(x) (wp[x].flag & WP_VISITED)
#define IS_CORE(x) (wp[x].flag & WP_CORE)
#define IS_SHELL(x) (wp[x].flag & WP_SHELL)
#define IS_BOUND(x) (wp[x].flag & WP_BOUND)
#define IS_REMAINING(x) (wp[x].flag & WP_REMAINING)

#define TOGGLE_PEAK(x) (wp[x].flag ^= WP_PEAK)
#define TOGGLE_VISITED(x) (wp[x].flag ^= WP_VISITED)
#define TOGGLE_CORE(x) (wp[x].flag ^= WP_CORE)
#define TOGGLE_SHELL(x) (wp[x].flag ^= WP_SHELL)
#define TOGGLE_BOUND(x) (wp[x].flag ^= WP_BOUND)
#define TOGGLE_REMAINING(x) (wp[x].flag ^= WP_REMAINING)

////////
#define RESET_WHOLE_FLAGS_WPT(x) (int _i;for(_i=0;_i<MAXTHREADS;_i++) wp[x].flag2[_i]  = WP_RESET)
#define SET_VISITEDT(x,it) (wp[x].flag2[it] |= WP_VISITEDT)
#define UNSET_VISITEDT(x,it) (wp[x].flag2[it] &= (~WP_VISITEDT))
#define IS_VISITEDT(x,it) (wp[x].flag2[it] & WP_VISITEDT)
#define TOGGLE_VISITEDT(x,it) (wp[x].flag2[it] ^= WP_VISITEDT)

#define SET_MEMBER_ID(a,b) (wp[a].haloid=b)
////////
//

/*
unsigned char setbits(unsigned char x,int p) {
    unsigned char b;
	b = x|(~((~0)<<1) << (p-1));
	return b;
}
*/


typedef struct WorkingParticle{
	float den;
	int haloid;
	unsigned char flag;
	unsigned char flag2[MAXTHREADS];
} WorkingParticle;
WorkingParticle *wp;


#define CORE_RESET 0x00
#define CORE_ENCLOSED 0x01
#define CORE_CONFIRMED 0x02
#define REAL_CORE 0x04

#define RESET_CORE_ENCLOSED(x) ( core[x].flag = CORE_RESET)
#define SET_CORE_ENCLOSED(x) ( core[x].flag |= CORE_ENCLOSED)
#define UNSET_CORE_ENCLOSED(x) ( core[x].flag &= (~CORE_ENCLOSED))
#define IS_CORE_ENCLOSED(x) ( core[x].flag & CORE_ENCLOSED)
#define TOGGLE_CORE_ENCLOSED(x) ( core[x].flag ^= CORE_ENCLOSED)

#define RESET_REAL_CORE(x) ( core[x].flag = CORE_RESET)
#define SET_REAL_CORE(x) ( core[x].flag |= REAL_CORE)
#define UNSET_REAL_CORE(x) ( core[x].flag &= (~REAL_CORE))
#define IS_REAL_CORE(x) ( core[x].flag & REAL_CORE)
#define TOGGLE_REAL_CORE(x) ( core[x].flag ^= REAL_CORE)
/* SCORE means the shell cores */
#define RESET_SCORE_ENCLOSED(x) ( score[x].flag = CORE_RESET)
#define SET_SCORE_ENCLOSED(x) ( score[x].flag |= CORE_ENCLOSED)
#define UNSET_SCORE_ENCLOSED(x) ( score[x].flag &= (~CORE_ENCLOSED))
#define IS_SCORE_ENCLOSED(x) ( score[x].flag & CORE_ENCLOSED)
#define TOGGLE_SCORE_ENCLOSED(x) ( score[x].flag ^= CORE_ENCLOSED)

#define RESET_SCORE_CONFIRMED(x) ( score[x].flag = CORE_RESET)
#define SET_SCORE_CONFIRMED(x) ( score[x].flag |= CORE_CONFIRMED)
#define UNSET_SCORE_CONFIRMED(x) ( score[x].flag &= (~CORE_CONFIRMED))
#define IS_SCORE_CONFIRMED(x) ( score[x].flag & CORE_CONFIRMED)
#define TOGGLE_SCORE_CONFIRMED(x) ( score[x].flag ^= CORE_CONFIRMED)


/*
typedef struct Coretype{
	int peak;
	int nummem, numstar;
	float starmass;
	float coredensity;
	float cx,cy,cz;
	float cvx,cvy,cvz;
	float Rtidal,density;
	unsigned char flag;
}Coretype;
*/

typedef struct Coresorttype{
	int nmem;
	Coretype *core;
}Coresorttype;
typedef struct Coresortdentype{
	int nmem;
	float density;
	Coretype *core;
	unsigned char flag;
}Coresortdentype;


#define NSHELL2P(x) (shell[x].nlinkedparticle) // number of particles enclosed by the shell
#define SHELL2P(x) (shell[x].linkedparticle) // shell to particle linked list
#define NSHELL2C(x) (shell[x].nlinkedcore) // number of cores enclosed by the shell
#define SHELL2C(x) (shell[x].linkedcore) // shell to the core's id number
typedef struct Shelltype{
	int nlinkedcore;
	int *linkedcore;
	int nlinkedparticle;
	int *linkedparticle;
}Shelltype;
Shelltype shell[1000000];
int nshell,numcore;

int sortcoreden(const void *a, const void *b){
	Coretype *aa,*bb;
	aa = (Coretype *)a;
	bb = (Coretype *)b;
	if(aa->density < bb->density) return 1;
	else if(aa->density > bb->density) return -1;
	else return 0;
}

int sortcorenumstar(const void *a, const void *b){
	Coretype *aa,*bb;
	aa = (Coretype *)a;
	bb = (Coretype *)b;
	if(aa->numstar < bb->numstar) return 1;
	else if(aa->numstar > bb->numstar) return -1;
	else return 0;
}

int  MergingPeak(SimpleBasicParticleType *bp,int np,Coretype *core,int numcore, int iflag){
	int i,j,k;
	float fof_link = MERGINGPEAKLENGTH;
	FoFTPtlStruct *ptl = (FoFTPtlStruct *) Malloc(sizeof(FoFTPtlStruct)*numcore,PPTR(ptl));
	particle *linked = (particle *)Malloc(sizeof(particle)*numcore,PPTR(linked));
	size_t nnode = MAX(65*10000,numcore);
	FoFTStruct *TREE = (FoFTStruct *)Malloc(sizeof(FoFTStruct)*nnode,PPTR(TREE));
	if(iflag ==0) { 
		qsort(core,numcore,sizeof(Coretype),sortcoreden);
	}
	else if(iflag == 1){
		qsort(core,numcore,sizeof(Coretype),sortcorenumstar);
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<numcore;i++){
		ptl[i].type = TYPE_PTL;
		ptl[i].x = bp[core[i].peak].x;
		ptl[i].y = bp[core[i].peak].y;
		ptl[i].z = bp[core[i].peak].z;
		ptl[i].link02 = fof_link;
		ptl[i].indx = 1;
		ptl[i].sibling = &ptl[i+1];
		ptl[i].included = NO;
	}
	ptl[numcore-1].sibling = NULL;

	int recursiveflag;
	if(nnode > 65*10000) recursiveflag= PTHREAD;
	else recursiveflag = RECURSIVE;
	recursiveflag = SERIALIZED;
	DEBUGPRINT("Before core merging %d\n", numcore);
	FoF_Make_Tree(TREE,nnode, ptl,numcore, recursiveflag);

	for(i=0;i<numcore;i++){
		if(ptl[i].included == NO){
			particle p;
			p.x = ptl[i].x;
			p.y = ptl[i].y;
			p.z = ptl[i].z;
			p.link02 = fof_link;
			int nlink = new_fof_link(&p,fof_link,TREE,ptl,linked);
			ptl[i].included = NO;
		}
	}
	j = 0;
	if(iflag==0){
		for(i=0;i<numcore;i++){
			if(ptl[i].included == NO){
				core[j] = core[i];
				j++;
			}
		}
	}
	else if(iflag==1){
		for(i=0;i<numcore;i++){
			if(ptl[i].included == NO){
				core[j] = core[i];
				SET_PEAK(core[i].peak);
				j++;
			}
			else {
				UNSET_PEAK(core[i].peak);
			}
		}
	}

	numcore = j;

	Free(TREE); Free(linked); Free(ptl);

	return numcore;
}
int finddenpeak(float *den,int numneigh,long long *neighbor,int np,
		Coretype **Core, int jflag, SimpleBasicParticleType *bp){
	Coretype *core = *Core;
	long long i,j,k;
	float *me,*you;
	int iflag;
	numcore = 0;
	DEBUGPRINT("Now before finddenpeak with %d particles with flag= %d\n", np,jflag);
	if(jflag == 1){
		for(i=0;i<np;i++){
			if(den[i] > PEAKTHRESHOLD && bp[i].type == TYPE_STAR)
			{
				iflag = 1;
				k = i*numneigh;
				for(j=0;j<numneigh;j++){
					if(den[neighbor[k+j]] > den[i]) {
						iflag = 0;
						break;
					}
				}
				if(iflag == 1){
					core[numcore].peak = i;
					core[numcore].cx = bp[i].x;
					core[numcore].cy = bp[i].y;
					core[numcore].cz = bp[i].z;
					core[numcore].density = den[i];
					if(numcore >= maxnumcore-10) {
						maxnumcore += MAXNUMCORE;
						*Core = Realloc(*Core, sizeof(Coretype)*maxnumcore);
						core = *Core;
					}
					/*
						printf("Warning.. insufficient cores: %d with max ncores= %d for i= %d np= %d\n", numcore, MAXNUMCORE, i, np);
					*/
					numcore ++;
				}
			}
		}
	}
	else {
		for(i=0;i<np;i++){
			if(den[i] > PEAKTHRESHOLD)
			{
				iflag = 1;
				k = i*numneigh;
				for(j=0;j<numneigh;j++){
					if(den[neighbor[k+j]] > den[i]) {
						iflag = 0;
						break;
					}
				}
				if(iflag == 1){
					if(numcore >= MAXNUMCORE){
						fprintf(stderr,"Error exceeding the number of cores: %d :: %lld   %d   \n", numcore, i, np);
						exit(999);
					}
					core[numcore].peak = i;
					core[numcore].cx = bp[i].x;
					core[numcore].cy = bp[i].y;
					core[numcore].cz = bp[i].z;
					core[numcore].density = den[i];

					if(numcore >= maxnumcore-10) {
						maxnumcore += MAXNUMCORE;
						core = Realloc(core, sizeof(Coretype)*maxnumcore);
						*Core = core;
					}
					numcore ++;
				}
			}
		}
	}
	DEBUGPRINT("Now before merging peak. numcore= %d\n", numcore);
//	if(numcore > 10) numcore = MergingPeak(bp,np,core,numcore,0);
	DEBUGPRINT("Now after merging peak\n");
	return numcore;
}
void Dump2WorkingParticle(float *density,int np,Coretype *core,int numcore){
	int i,now;
	for(i=0;i<np;i++) {
		wp[i].den = density[i];
		RESET_WHOLE_FLAGS_WP(i);
	}
	for(i=0;i<numcore;i++) SET_PEAK(core[i].peak);
}
int CountMemInRadius2(SimpleBasicParticleType *bp,int np,int *plist,
		float cx,float cy,float cz,float radius2){
	int i,j,k;
	int ncount;
	float tmpx,tmpy,tmpz,dist2;
	ncount = 0;
	for(i=0;i<np;i++){
		tmpx = bp[plist[i]].x-cx; 
		tmpy = bp[plist[i]].y-cy; 
		tmpz = bp[plist[i]].z-cz;
		dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
		if(dist2 < radius2) ncount ++;
	}
	return ncount;
}
int GetShellParticleFromRestParticleLIST(SimpleBasicParticleType  *bp,int np,int ishell, int *plist,
		Coretype *icore){
	int i,j,k,bid;
	int nmem,omem;
	float cx,cy,cz,tmpx,tmpy,tmpz,tidalr2,dist2;

	nmem = 0;
	/* Include shell particles that are not bound to any halos. */
	for(i=0;i<NSHELL2P(ishell);i++)
		if(IS_BOUND((bid=SHELL2P(ishell)[i]))==NOT) plist[nmem++] = bid;
	for(i=0;i<np;i++) if(IS_REMAINING(i) != NOT) plist[nmem++] = i;

	cx = icore->cx; cy = icore->cy; cz = icore->cz;
	tidalr2 = icore->Rtidal;
	tidalr2 = tidalr2*tidalr2;
	omem = nmem;
	nmem = 0;
	for(i=0;i<omem;i++){
		tmpx = bp[plist[i]].x -cx;
		tmpy = bp[plist[i]].y -cy;
		tmpz = bp[plist[i]].z -cz;
		dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
		if(dist2 < tidalr2) {
			plist[nmem++] = plist[i];
		}
	}
	return nmem;
}
int GetMemCandidateList(int np,int haloid, int *plist,int ishell){
	int i,j,k;
	int nmem,nlcore;

	nmem  = 0;
	for(i=0;i<np;i++) if(wp[i].haloid == haloid) plist[nmem++] = i;
	for(i=ishell+1;i<nshell;i++){
		nlcore = NSHELL2C(i);
		for(j = 0;j<nlcore;j++)
			if(SHELL2C(i)[j] == haloid)
				for(k=0;k<NSHELL2P(i);k++)
					plist[nmem++] = SHELL2P(i)[k];

	}
	return nmem;
}
int GET_NONENCLOSED_ID(Coretype  *core, int numcore,float threshold){
	int i;
	for(i=0;i<numcore;i++) {
		if(IS_CORE_ENCLOSED(i) == NOT && core[i].coredensity > threshold) {
			return i;
		}
	}
	return -1;
}
void SaveRemainingParticles2LastShell(int np,Coretype  *core,int numcore){
	int i,j,k;
	int nrest,ishell;

	nrest = 0;
	for(i=0;i<np;i++) if(IS_CORE(i) == NOT && IS_SHELL(i) ==NOT) nrest++;
	if(nrest ==0) return;

	if(nshell >0) {
		ishell = nshell -1;
		if(SHELL2P(ishell)) 
			SHELL2P(ishell) = Realloc(SHELL2P(ishell),sizeof(int)*(NSHELL2P(ishell)+nrest));
		else SHELL2P(ishell) = Malloc(sizeof(int)*nrest,PPTR(SHELL2P(ishell)));

		for(i=0;i<np;i++)
			if(IS_CORE(i) == NOT && IS_SHELL(i) ==NOT) {
				SET_SHELL(i);
				UNSET_REMAINING(i);
				SHELL2P(ishell)[(NSHELL2P(ishell)++)] = i;
				SET_MEMBER_ID(i,NOT_HALO_MEMBER);
			}
	}
	else if(nshell ==0){
		SHELL2P(nshell) = Malloc(sizeof(int)*nrest,PPTR(SHELL2P(nshell)));
		NSHELL2P(nshell) = 0;
		for(i=0;i<np;i++)
			if(IS_CORE(i) == NOT && IS_SHELL(i) ==NOT) {
				SET_SHELL(i);
				UNSET_REMAINING(i);
				SET_MEMBER_ID(i,NOT_HALO_MEMBER);
				SHELL2P(nshell)[(NSHELL2P(nshell)++)] = i;
			}
		NSHELL2C(nshell) = numcore;
		SHELL2C(nshell) = Malloc(sizeof(int)*numcore,PPTR(SHELL2C(nshell)));
		for(i=0;i<numcore;i++) SHELL2C(nshell)[i] = i;
		nshell++;
	}
	else {
		fprintf(stderr,"Strange value of nshell in SaveRemainingParticles2LastShell\n");
		fprintf(stderr,"Ok.. Now exiting with nshell =%d\n",nshell);
		exit(999);
	}
}
void GetShellParticlesPeaks(int np,long long *neighbor,int NumNeighbor,
		Coretype *core,int numcore, float dthreshold){
	long long i,j,k,jj,kk;
	int icore,*contactlist,nlist=0;
	int numenclosedpeaks;
	int now,new;
	int nowcore;

	contactlist = (int *)Malloc(sizeof(int)*np,PPTR(contactlist));
	for(i=0;i<numcore;i++) RESET_CORE_ENCLOSED(i);

	for(j=0;j<np;j++) UNSET_VISITED(j);

	while((nowcore = GET_NONENCLOSED_ID(core,numcore,dthreshold))>=0){
		for(j=0;j<nlist;j++) {
	 		int jj = contactlist[j];
			UNSET_VISITED(jj);// unmarking particles (initialization)
		}
		now = nlist = 0;
		SET_VISITED((contactlist[nlist++]= core[nowcore].peak));
		while(now<nlist){
			kk = (long long)contactlist[now]*(long long)NumNeighbor;
			for(k=0;k<NumNeighbor;k++){
				new = neighbor[kk];
				if(wp[new].den > dthreshold && IS_VISITED(new) == NOT){
					SET_VISITED(new);
					contactlist[nlist++] = new;
				}
				kk++;
			}
			now++;
		}
		if(nlist==0) continue;
		/* erase shell or core particles from this shell-particle list
		 * and turn on the shell flag */
		SHELL2P(nshell) = Malloc(sizeof(int)*nlist,PPTR(SHELL2P(nshell)));
		now = 0;
		for(j=0;j<nlist;j++){
			if(IS_SHELL(contactlist[j])==NOT && IS_CORE(contactlist[j])==NOT){
				SHELL2P(nshell)[now++] = contactlist[j];
				SET_SHELL(contactlist[j]);
				UNSET_REMAINING(contactlist[j]);
				SET_MEMBER_ID(contactlist[j],NOT_HALO_MEMBER);
			}
		}
		NSHELL2P(nshell) = now;
		SHELL2P(nshell) = Realloc(SHELL2P(nshell),sizeof(int)*NSHELL2P(nshell));

		NSHELL2C(nshell) = 0;
		SHELL2C(nshell) = (int *)Malloc(sizeof(int)*numcore,PPTR(SHELL2C(nshell)));
		for(j=0;j<nlist;j++){
			now = contactlist[j];
			if(IS_PEAK(now)!=NOT){
				for(icore=0;icore<numcore;icore++){
					if(core[icore].peak == now){
						SET_CORE_ENCLOSED(icore);
						SHELL2C(nshell)[NSHELL2C(nshell)++] = icore;
						break;
					}
				}
			}
		}
		SHELL2C(nshell) = Realloc(SHELL2C(nshell),sizeof(int)*NSHELL2C(nshell));
		nshell ++;
	}
	Free(contactlist);
}
#ifdef SCORE_NMEM
int coresortden(const void *a,const void *b){
	Coresortdentype *aa,*bb;
	aa = (Coresortdentype *)a;
	bb = (Coresortdentype *)b;
	if(aa->nmem < bb->nmem) return 1;
	else if(aa->nmem > bb->nmem) return -1;
	else return 0;
}
#else
int coresortden(const void *a,const void *b){
	Coresortdentype *aa,*bb;
	aa = (Coresortdentype *)a;
	bb = (Coresortdentype *)b;
	if(aa->density < bb->density) return 1;
	else if(aa->density > bb->density) return -1;
	else return 0;
}
#endif
/* This function is trying to merge the underpopulated cores to
 * fulfill the core criterion (> MINCOREMEM) */
int DMSmartFinding(SimpleBasicParticleType *bp,int np,Coretype *core,int numcore,
		long long *neighbor,int NumNeighbor){
	int i,j,k;
	int num;
	int *contactlist;
	Coresortdentype *score;
	float minden;
	num = 0;
	for(i=0;i<numcore;i++) if(core[i].nummem < MINCORENMEM) {
#ifdef DEBUG
		printf("DMSmartingFinding: %d'th core is being erased for %d < mincorenemme= %d\n", i,core[i].nummem,MINCORENMEM);
#endif
		num++;
	}
#ifdef DEBUG
	printf("DMSmartFinding: The number of erased core is %d from %d for mincore = %d\n",num,numcore, MINCORENMEM);
#endif
	if(num <2) {
		num = 0;
		for(i=0;i<numcore;i++) 
			if(core[i].nummem >= MINCORENMEM) core[num++]=core[i];
		return num;
	}
	score = (Coresortdentype *)Malloc(sizeof(Coresortdentype)*num,PPTR(score));
	num = 0;
	/* Dump the core data to the sorted core array */
	for(i=0;i<numcore;i++) if(core[i].nummem < MINCORENMEM){
		score[num].nmem = core[i].nummem;
		score[num].core = core + i;
		score[num].density = wp[core[i].peak].den; /* peak density */
		num++; /* num is the total number of cores that is underpopulated
				  while numcore is the total number of cores. */
	}
	minden = 2.e23;
	for(i=0;i<np;i++) minden = MIN(minden,wp[i].den);
	/* Sort score in reverse order of peak density */
	qsort(score,num,sizeof(Coresortdentype),coresortden);
	/* Unset the peak flag for the underpopulated cores.
	   Erase underpopulated peaks from peak (core) list and only consider the sorted cores. */
	for(i=0;i<numcore;i++) 
		if(core[i].nummem < MINCORENMEM) UNSET_PEAK(core[i].peak);
	for(i=0;i<num;i++) {
		RESET_SCORE_ENCLOSED(i);
		RESET_SCORE_CONFIRMED(i);
	}
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
	for(j=0;j<np;j++) {
		UNSET_VISITED(j);
	}

	/* Now rock'n roll */
	contactlist = (int *)Malloc(sizeof(int)*np,PPTR(contactlist));
	for(i=0;i<num;i++){/* From the highest density peak that are underpopulated */
		if(IS_SCORE_ENCLOSED(i) == NOT){
			float upden = (score[i].core)->coredensity;
			float downden = minden;
			float denthr;
			int ncontact=0, now,new;
			int mcontact;
			float corestarmass = 0;
			int iter = 0;
			while(iter==0 || fabs(upden-downden)/downden>COREDENRESOLUTION){
 				for(j=0;j<ncontact;j++) {
					int jj = contactlist[j];
					UNSET_VISITED(jj);/* Set all the particle not to be visited */
				}
				int breakflag = ncontact = now = 0;
				mcontact = 0;

				corestarmass = 0;

				denthr = 0.25*(3.*upden+downden);
				SET_VISITED((contactlist[ncontact++] = (score[i].core)->peak));
				while(now<ncontact && breakflag ==0){
					long long kk = (long long)contactlist[now]*(long long)NumNeighbor;
					for(k=0;k<NumNeighbor;k++){
						new = neighbor[kk];
						if(wp[new].den>denthr){
							if(IS_VISITED(new) ==NOT && IS_PEAK(new) == NOT){
								SET_VISITED(new);
								contactlist[ncontact++] = new;
								if(bp[new].type == TYPE_STAR) mcontact ++;
								corestarmass += bp[new].mass;
							}
							else if(IS_VISITED(new) == NOT && IS_PEAK(new) !=NOT){
								downden = denthr;
								breakflag = 1;
								break;
							}
						}
						kk++;
					}
					now++;
				}
				if(breakflag==0) {
					upden = denthr;
					/* This is newly inserted because we don't have to
					 * find the exact value of the threashold. Rather than that
					 * it is sufficient to satisfy the number of core particles
					 * should be larger than and equal to "MINCORENMEM"
					 * */
					if(mcontact >= MINCORENMEM && corestarmass >= MINCORESTARMASS) break;
				}
				iter ++;
			}
			(score[i].core)->coredensity = (denthr = upden);

			/* Now scoop up core particles */
//			for(j=0;j<np;j++) UNSET_VISITED(j);
			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				UNSET_VISITED(jj);
			}
			ncontact = now = 0;
			mcontact = 0;

			corestarmass = bp[contactlist[ncontact]].mass;

			SET_VISITED((contactlist[ncontact++] = (score[i].core)->peak));
			while(now<ncontact){
				long long kk = (long long)contactlist[now]*(long long)NumNeighbor;
				for(k=0;k<NumNeighbor;k++){
					new = neighbor[kk];
					if(wp[new].den > denthr && IS_VISITED(new) == NOT){
						SET_VISITED(new);
						if(bp[new].type == TYPE_STAR) mcontact ++;
						contactlist[ncontact++] = new;
						corestarmass += bp[new].mass;
					}
					kk++;
				}
				now++;
			}

			if(mcontact >= MINCORENMEM && corestarmass >= MINCORESTARMASS) {/* restore this as a meaningful peak */
#ifdef DEBUG
				printf("New merging core is detected with np= %d in %d cores\n",ncontact,i );
#endif
				SET_PEAK((score[i].core)->peak);
				SET_SCORE_CONFIRMED(i);
			}
			else {
				UNSET_PEAK((score[i].core)->peak);
				UNSET_SCORE_CONFIRMED(i);
			}
			for(j=i+1;j<num;j++){
				k = (score[j].core)->peak;
				if(IS_VISITED(k) != NOT ) {
					/* delete this peak forever since there is no possibility 
					 * for this peak to get a meaning. */
					SET_SCORE_ENCLOSED(j);
					UNSET_SCORE_CONFIRMED(j);
				}
			}
			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				UNSET_VISITED(jj);
			}
		}
	}
	Free(contactlist);

	for(i=0;i<numcore;i++) SET_REAL_CORE(i);
	for(i=0;i<num;i++) 
		if(IS_SCORE_CONFIRMED(i) == NOT) UNSET_REAL_CORE(SCOREID2COREID(i));
	num = 0;
	for(i=0;i<numcore;i++){
		if(IS_REAL_CORE(i) != NOT){/* If this is real core */
			core[num++] = core[i];
		}
	}
#ifdef DEBUG
	printf("Total number of survived cores is %d from %d\n",num,numcore);
#endif
	numcore = num;
	return numcore;
}
/* This function is trying to merge the underpopulated cores to
 * fulfill the core criterion (> MINCOREMEM) */
int SmartFinding(SimpleBasicParticleType *bp,int np,Coretype *core,int numcore,
		long long *neighbor,int NumNeighbor){
	long long i,j,k;
	int num;
	int *contactlist;
	Coresortdentype *score;
	float minden;
	num = 0;
	for(i=0;i<numcore;i++) if(core[i].numstar < MINCORENMEM || core[i].starmass < MINCORESTARMASS) num++;
	DEBUGPRINT("SmartFinding: The number of erased core is %d from %d for mincore = %d\n",num,numcore, MINCORENMEM);
	if(num <2) {
		num = 0;
		for(i=0;i<numcore;i++) 
			if(core[i].numstar >= MINCORENMEM && core[i].starmass >= MINCORESTARMASS) core[num++]=core[i];
		return num;
	}
	score = (Coresortdentype *)Malloc(sizeof(Coresortdentype)*num,PPTR(score));
	num = 0;
	/* Dump the core data to the sorted core array */
	for(i=0;i<numcore;i++) if(core[i].numstar < MINCORENMEM || core[i].starmass < MINCORESTARMASS){
		score[num].nmem = core[i].nummem;
		score[num].core = core + i;
		score[num].density = wp[core[i].peak].den; /* peak density */
		DEBUGPRINT("C%lld has num= %d peakden= %g coreden= %g\n", i, core[i].nummem,
				wp[core[i].peak].den, core[i].coredensity);
		num++; /* num is the total number of cores that is underpopulated
				  while numcore is the total number of cores. */
	}
	minden = 2.e23;
	for(i=0;i<np;i++) minden = MIN(minden,wp[i].den);
	/* Sort score in reverse order of peak density */
	qsort(score,num,sizeof(Coresortdentype),coresortden);
	/* Unset the peak flag for the underpopulated cores.
	   Erase underpopulated peaks from peak (core) list and only consider the sorted cores. */
	for(i=0;i<numcore;i++) 
		if(core[i].numstar < MINCORENMEM || core[i].starmass < MINCORESTARMASS) UNSET_PEAK(core[i].peak);
	for(i=0;i<num;i++) {
		RESET_SCORE_ENCLOSED(i);
		RESET_SCORE_CONFIRMED(i);
	}
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
	for(j=0;j<np;j++) {
		UNSET_VISITED(j);
	}

	// Now rock'n roll 
	contactlist = (int *)Malloc(sizeof(int)*np,PPTR(contactlist));
	for(i=0;i<num;i++){/* From the highest density peak that are underpopulated */
		if(IS_SCORE_ENCLOSED(i) == NOT){
			float upden = (score[i].core)->coredensity;
			float downden = minden;
			float denthr;
			int ncontact=0, now,new;
			int mcontact;
			float corestarmass = 0;
			int iter = 0;
//			while(iter==0 || fabs(upden-downden)/denthr>COREDENRESOLUTION)
			denthr = 0.5*(upden+downden);
			DEBUGPRINT("SC%lld has upden= %g downden= %g denthr= %g res= %g\n", 
						i, upden,downden, denthr, COREDENRESOLUTION);
			do {
 				for(j=0;j<ncontact;j++) {
					int jj = contactlist[j];
					UNSET_VISITED(jj);/* Set all the particle not to be visited */
				}
				int breakflag = ncontact = now = 0;
				mcontact = 0;

				corestarmass = 0;

				denthr = 0.5*(upden+downden);
				SET_VISITED((contactlist[ncontact++] = (score[i].core)->peak));
				while(now<ncontact && breakflag ==0){
					long long kk = (long long)contactlist[now]*(long long)NumNeighbor;
					for(k=0;k<NumNeighbor;k++){
						new = neighbor[kk];
						if(wp[new].den>denthr){
							if(IS_VISITED(new) ==NOT && IS_PEAK(new) == NOT){
								SET_VISITED(new);
								contactlist[ncontact++] = new;
								if(bp[new].type == TYPE_STAR) {
									mcontact ++; corestarmass += bp[new].mass;
								}
							}
							else if(IS_VISITED(new) == NOT && IS_PEAK(new) !=NOT){
								downden = denthr;
								breakflag = 1;
								break;
							}
						}
						kk++;
					}
					now++;
				}
				if(breakflag==0) {
					upden = denthr;
					/* This is newly inserted because we don't have to
					 * find the exact value of the threashold. Rather than that
					 * it is sufficient to satisfy the number of core particles
					 * should be larger than and equal to "MINCORENMEM"
					 * */
					if(mcontact >= MINCORENMEM && corestarmass >= MINCORESTARMASS) break;
				}
				iter ++;
			} while(fabs(upden-downden)/denthr>COREDENRESOLUTION 
					&& fabs(minden-(score[i].core)->coredensity)>1.);
			(score[i].core)->coredensity = (denthr = upden);

			// Now scoop up core particles
//			for(j=0;j<np;j++) UNSET_VISITED(j);
			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				UNSET_VISITED(jj);
			}
			ncontact = now = 0;
			mcontact = 0;

			corestarmass = bp[contactlist[ncontact]].mass;

			SET_VISITED((contactlist[ncontact++] = (score[i].core)->peak));
			while(now<ncontact){
				long long kk = (long long)contactlist[now]*(long long)NumNeighbor;
				for(k=0;k<NumNeighbor;k++){
					new = neighbor[kk];
					if(wp[new].den > denthr && IS_VISITED(new) == NOT){
						SET_VISITED(new);
						if(bp[new].type == TYPE_STAR) mcontact ++;
						contactlist[ncontact++] = new;
						corestarmass += bp[new].mass;
					}
					kk++;
				}
				now++;
			}

			if(mcontact >= MINCORENMEM && corestarmass >= MINCORESTARMASS) {
				// restore this as a meaningful peak 
				DEBUGPRINT("New merging core is detected with np= %d in %lld cores\n",ncontact,i );
				SET_PEAK((score[i].core)->peak);
				SET_SCORE_CONFIRMED(i);
			}
			else {
				UNSET_PEAK((score[i].core)->peak);
				UNSET_SCORE_CONFIRMED(i);
			}
			for(j=i+1;j<num;j++){
				k = (score[j].core)->peak;
				if(IS_VISITED(k) != NOT ) {
					/* delete this peak forever since there is no possibility 
					 * for this peak to get a meaning. */
					SET_SCORE_ENCLOSED(j);
					UNSET_SCORE_CONFIRMED(j);
				}
			}
			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				UNSET_VISITED(jj);
			}
		}
	}
	Free(contactlist);

	for(i=0;i<numcore;i++) SET_REAL_CORE(i);
	for(i=0;i<num;i++) 
		if(IS_SCORE_CONFIRMED(i) == NOT) UNSET_REAL_CORE(SCOREID2COREID(i));
	num = 0;
	for(i=0;i<numcore;i++){
		if(IS_REAL_CORE(i) != NOT){/* If this is real core */
			core[num++] = core[i];
		}
	}
	DEBUGPRINT("Total number of survived cores is %d from %d\n",num,numcore);
	numcore = num;
	return numcore;
}
/* Comerford MNRAS 379 190 */
double get_concentration(double mass){
	double mstar= 1.3E13L;
	double redshift;
	double res;
	redshift = amax/a - 1.L;
	res = 14.5 *pow(mass/mstar,-0.15)/(1.L+redshift);
	return res;
}
double get_alpha(double c,double s){
	double res;
	{
		res = 2.L-(c*s/(1.L+c*s))*(c*s/(1.L+c*s));
		res = res/(log(1.L+c*s)-c*s/(1.L+c*s));
		return res;
	}
}
double Get_Vc2(double mass, double Rv){
	double res;
	res = G*mass/Rv;
	return res;
}
/* Definition of v for flat cosmologies */
#define RHOC 2.7755E11L
#define vv (18.L*M_PI*M_PI+82.L*(omep-1.L)-39.L*pow(omep-1.L,2.L))
double gnfw(double c){
	double res;
	res =  1./(log(1.L+c)-c/(1.L+c));
	return res;
}
double GetRvir(double Mvir){
	double res;
	res = pow(3.*Mvir/(4.L*M_PI*vv*RHOC),0.3333333333333L);
	return res;
}
double GetBGpotent(double Mvir){
	double massincgs = Mvir*Msun;
	double c = get_concentration(Mvir);
	double gc = gnfw(c);
	double Rvir = GetRvir(Mvir);
	double Vcir2 = Get_Vc2(Mvir, Rvir)*potentfact;
	double halobgpotent = Vcir2*(1-gc*log(1+c));
	return halobgpotent;
}
double get_tidal_ellipse(int np, int myid, int mp, int uid, int mpall, Coretype *core, float Mvir,
		double mratio, double bmass){
	double alpha,ellipse_contribution;
	double w2,wx,wy,wz,dx,dy,dz,dist2,dvx,dvy,dvz;
	double dist,pot,res,s,Rvir,v2,nfwc;
	double kpdratio;

	dx = core[myid].cx - core[uid].cx;
	dy = core[myid].cy - core[uid].cy;
	dz = core[myid].cz - core[uid].cz;
	dist = sqrt(dx*dx+dy*dy+dz*dz);
	dvx = core[myid].cvx - core[uid].cvx;
	dvy = core[myid].cvy - core[uid].cvy;
	dvz = core[myid].cvz - core[uid].cvz;
	v2 = dvx*dvx+dvy*dvy+dvz*dvz;
	wx = dy*dvz-dz*dvy;
	wy = dz*dvx-dx*dvz;
	wz = dx*dvy-dy*dvx;
	w2 = (wx*wx+wy*wy+wz*wz)/(dist*dist);
	pot = Mvir/dist;
	kpdratio = r2kineticfact*r2kineticfact/potentfact;
	ellipse_contribution = w2/pot * kpdratio;
	/* This is for the case when the minor body is not bound to the major body. */
	/* This is the case when the pair do not bound to each other. 
	 * In this case the contribution of the orbital motion to the
	 * tidal radius is assumed to be negligible. */
/*
#error This should be modified further.
*/

	Rvir = pow(3.*Mvir/(4.L*M_PI*vv*RHOC),0.3333333333333L);
	/* Change to the simulation length scale */
	/*
	Rvir = Rvir/size*rng;
	*/
	s = dist/Rvir;

	nfwc = get_concentration(bmass);
#ifdef OLD_BETA
	if(ellipse_contribution > 2.L){
		ellipse_contribution = 0.L;
	}
#else
	if(s>1) {
		double fact;
		fact = v2/(bmass/Rvir)*kpdratio;
		if(fact > 2.L) ellipse_contribution = 0.L;
	}
	else {
		double fact;
		fact = v2*s/(gnfw(nfwc)*log(1.L+nfwc*s)*bmass/Rvir )*kpdratio;
		if(fact > 2.L) ellipse_contribution = 0.L;
	}
#endif
//	alpha = get_alpha(get_concentration(bmass),s);
	alpha = get_alpha(nfwc,s);
	res = pow(mratio/(ellipse_contribution+alpha),0.33333333333333L);
#ifdef DEBUG
	/*
	printf(" ellipse_c= %g alpha= %g c= %g s= %g\n",ellipse_contribution,alpha,get_concentration(mpall),s);
	*/
#endif
	return res;
}

int  FindCoreDensity(SimpleBasicParticleType *bp,int np,long long *neighbor,
		int NumNeighbor,Coretype *core,int numcore){
	int i,j,k;
//	float denthr,upden,downden;
	int *Tcontactlist;
	int newnmem,iflag;
	float minden;

	iflag = 0;
recycling:
	for(i=0;i<np;i++) SET_MEMBER_ID(i,NOT_HALO_MEMBER);
	minden = 1.e23;
	for(i=0;i<np;i++) minden = MIN(minden,wp[i].den);
	minden = MAX(DENFLOOR, minden);
	DEBUGPRINT("the minimum density %g \n", minden);

	int nthreads=1;
#ifdef _OPENMP
#pragma omp parallel
	{
		int it = omp_get_thread_num();
		if(it ==0) nthreads = omp_get_num_threads();
	}
	nthreads = MIN(nthreads, MAXTHREADS);
#endif



	size_t numlinkingwatershedding= MIN(np, MAXNUMWATERSHEDDING);

	Tcontactlist = (int *)Malloc(sizeof(int)*numlinkingwatershedding*nthreads,
			PPTR(Tcontactlist));

#ifdef _OPENMP
#pragma omp parallel for private(j,i)
#endif
	for(j=0;j<np;j++) {
		for(i=0;i<nthreads;i++) UNSET_VISITEDT(j,i);
	}

	for(i=0;i<np;i++) RESET_WHOLE_FLAGS_WP(i);
	for(i=0;i<numcore;i++) SET_PEAK(core[i].peak);


	/*
#ifdef _OPENMP
#pragma omp parallel private(i,j,k) num_threads(nthreads)
#endif
	{
		int it = 0;
#ifdef _OPENMP
		it = omp_get_thread_num();
#endif
		int *contactlist = Tcontactlist + (long)it*numlinkingwatershedding;
		for(i=it;i<numcore;i+=nthreads){
		*/
	{
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k) num_threads(nthreads) schedule(dynamic)
#endif
		for(i=0;i<numcore;i++){
			int it = omp_get_thread_num();
			int *contactlist = Tcontactlist + (long long)it*numlinkingwatershedding;
			// Now find core density threshold 
			float upden = wp[core[i].peak].den;
			float downden = minden;
			float denthr;
			int mcontact;
			int ncontact=0, now;
			do{
				// initialization before a search for the core density 
 				for(j=0;j<ncontact;j++) {
					int jj = contactlist[j];
					UNSET_VISITEDT(jj,it);// Set all the particle not to be visited 
				}
				int breakflag = 0; 
				ncontact = now = 0;
				// Trial value 
				denthr = 0.5*(upden+downden); 
				// Now peak particle is included. 
				contactlist[ncontact++] = core[i].peak;
				SET_VISITEDT(core[i].peak,it); 
				while(now < ncontact && breakflag ==0){
					long long kk = (long long)contactlist[now]*(long long)NumNeighbor;
					for(k=0;k<NumNeighbor;k++){
						int new = neighbor[kk];
						if(wp[new].den>denthr) {
							if(IS_VISITEDT(new,it) == NOT && IS_PEAK(new) == NOT){
								SET_VISITEDT(new,it);
								contactlist[ncontact++] = new;
							}
							else if(IS_VISITEDT(new,it) == NOT && IS_PEAK(new) != NOT){
								breakflag = 1;
								break;
							}
						}
						kk ++;
					}
					now++;
				}
				if(breakflag==0) upden = denthr;
				else if(breakflag==1) downden = denthr;
			}while(fabs((upden-downden)/denthr)>COREDENRESOLUTION);
			core[i].coredensity = (denthr = upden);
			/* Now scoop up core particles */
			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				UNSET_VISITEDT(jj,it);
			}
			{
				int breakflag = 0; 
				ncontact = now = 0;
				// Trial value 
				// Now peak particle is included. 
				contactlist[ncontact++] = core[i].peak;
				SET_VISITEDT(core[i].peak,it); 
				while(now < ncontact){
					long long kk = (long long)contactlist[now]*(long long)NumNeighbor;
					for(k=0;k<NumNeighbor;k++){
						int new = neighbor[kk];
						if(wp[new].den>denthr) {
							if(IS_VISITEDT(new,it) == NOT && IS_PEAK(new) == NOT){
								SET_VISITEDT(new,it);
								contactlist[ncontact++] = new;
							}
							else if(IS_VISITEDT(new,it) == NOT && IS_PEAK(new) != NOT){
								breakflag = 1;
								break;
							}
						}
						kk ++;
					}
					now++;
				}
				if(breakflag ==1) { //no coredensity is found. discard this peak.
					ncontact = 0;
				}
			}
			/* Now scoop up core particles */
			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				UNSET_VISITEDT(jj,it);
			}

			/*
			ncontact = now = 0;
			float corestarmass = 0;
			mcontact = 0;
			if(bp[core[i].peak].type == TYPE_STAR) {
				mcontact ++;
				corestarmass += bp[core[i].peak].mass;
			}
			SET_VISITEDT((contactlist[ncontact++] = core[i].peak),it);
			while(now<ncontact){
				long long kk = (long long)contactlist[now]*(long long)NumNeighbor;
				for(k=0;k<NumNeighbor;k++){
					int new = neighbor[kk];
					if(wp[new].den > denthr && IS_VISITEDT(new,it) == NOT){

						if(bp[new].type== TYPE_STAR) {
							mcontact ++;
							corestarmass += bp[new].mass;
						}

						SET_VISITEDT(new,it);
						contactlist[ncontact++] = new;
					}
					kk++;
				}
				now++;
			}
			*/

			mcontact = 0;
			float corestarmass = 0;
			for(j=0;j<ncontact;j++){
				int jj = contactlist[j];
				if(bp[jj].type == TYPE_STAR){
					mcontact ++;
					corestarmass += bp[jj].mass;
				}
			}

			for(j=0;j<ncontact;j++) {// make sure these lines are thread-safe
				SET_CORE(contactlist[j]);
				SET_BOUND(contactlist[j]);
				UNSET_REMAINING(contactlist[j]);
				SET_MEMBER_ID(contactlist[j],i);
			}

			for(j=0;j<ncontact;j++) {
				int jj = contactlist[j];
				UNSET_VISITEDT(jj,it);
			}

			/* Now core[i].center is temporarily considered as the center for tidal radius calculation.*/
			/* But in case of kinetic energy measurement the center of mass is instead used. */
			if(iflag != 0)
			{
				double cx,cy,cz;
				double cvx,cvy,cvz;
				dptype tmass=0,amass=0;
				cx = cy = cz = cvx = cvy = cvz = 0;
				if(MINCORESTARMASS>0){
					for(j=0;j<ncontact;j++){
						if(bp[contactlist[j]].type == TYPE_STAR){
							amass = bp[contactlist[j]].mass;
							tmass += amass;
							cx += bp[contactlist[j]].x*amass;
							cy += bp[contactlist[j]].y*amass;
							cz += bp[contactlist[j]].z*amass;
							cvx += bp[contactlist[j]].vx*amass;
							cvy += bp[contactlist[j]].vy*amass;
							cvz += bp[contactlist[j]].vz*amass;
						}
					}
				}
				else {
					for(j=0;j<ncontact;j++){
						amass = bp[contactlist[j]].mass;
						tmass += amass;
						cx += bp[contactlist[j]].x*amass;
						cy += bp[contactlist[j]].y*amass;
						cz += bp[contactlist[j]].z*amass;
						cvx += bp[contactlist[j]].vx*amass;
						cvy += bp[contactlist[j]].vy*amass;
						cvz += bp[contactlist[j]].vz*amass;
					}
				}
				cx = cx/tmass;
				cy = cy/tmass;
				cz = cz/tmass;
				cvx = cvx/tmass;
				cvy = cvy/tmass;
				cvz = cvz/tmass;
				core[i].cx = cx; core[i].cy = cy; core[i].cz = cz;
				core[i].cvx = cvx; core[i].cvy = cvy; core[i].cvz = cvz;
			}
			/*
			if(i==1){ 
				FILE *wwp = fopen("C1.out","w");
				for(j=0;j<ncontact;j++) fprintf(wwp,"%g %g %g %g %d\n",
						bp[contactlist[j]].x/TSC_CELL_SIZE+NCELLBUFF/2.,
						bp[contactlist[j]].y/TSC_CELL_SIZE+NCELLBUFF/2.,
						bp[contactlist[j]].z/TSC_CELL_SIZE+NCELLBUFF/2.,
						wp[contactlist[j]].den,
						contactlist[j]
						);
				fclose(wwp);
			}
			if(i==697){ 
				FILE *wwp = fopen("C697.out","w");
				for(j=0;j<ncontact;j++) fprintf(wwp,"%g %g %g %g %d\n",
						bp[contactlist[j]].x/TSC_CELL_SIZE+NCELLBUFF/2.,
						bp[contactlist[j]].y/TSC_CELL_SIZE+NCELLBUFF/2.,
						bp[contactlist[j]].z/TSC_CELL_SIZE+NCELLBUFF/2.,
						wp[contactlist[j]].den,
						contactlist[j]
						);
				fclose(wwp);
			}
			*/
			core[i].nummem = ncontact;
			core[i].numstar = mcontact;
			core[i].starmass = corestarmass;
			float xpix,ypix,zpix;
			xpix = core[i].cx/TSC_CELL_SIZE+NCELLBUFF/2.;
			ypix = core[i].cy/TSC_CELL_SIZE+NCELLBUFF/2.;
			zpix = core[i].cz/TSC_CELL_SIZE+NCELLBUFF/2.;
			DEBUGPRINT("C%d has number %d in %d with mstar= %g den= %g/%g at %g %g %g :: %g %g %g\n",i,
					core[i].numstar, core[i].nummem, core[i].starmass,
					core[i].density, core[i].coredensity,
					core[i].cx, core[i].cy, core[i].cz, xpix,ypix,zpix);
		}
	}
	Free(Tcontactlist);
	DEBUGPRINT("Now Found the core densities for %d cores\n",numcore);
	for(i=0;i<numcore;i++){
		DEBUGPRINT("C%d has nstar= %d npall= %d xyz= %g %g %g\n", 
				i, core[i].numstar,core[i].nummem, core[i].cx, core[i].cy, core[i].cz);
	}
	if(iflag ==0){
		// Deleting peaks of core particles less than minimum number 
		iflag = 1;
#ifdef OLD
		newnmem  =0;
		for(i=0;i<numcore;i++){
			if(core[i].numstar >= MINCORENMEM && core[i].starmass >= MINCORESTARMASS) core[newnmem++] = core[i];
		}
		numcore = newnmem;
		for(i=0;i<np;i++) RESET_WHOLE_FLAGS_WP(i);
		for(i=0;i<numcore;i++) SET_PEAK(core[i].peak);
#else
//		if(MINSTELLARMASS >0) numcore = SmartFinding(bp,np,core,numcore,neighbor,NumNeighbor);
//		else  numcore = DMSmartFinding(bp,np,core,numcore,neighbor,NumNeighbor);
//		if(MINSTELLARMASS<=0) numcore = DMSmartFinding(bp,np,core,numcore,neighbor,NumNeighbor);

//		DEBUGPRINT("Now before merging peak\n");
//		if(numcore > 10) numcore = MergingPeak(bp,np,core,numcore,1);
//		DEBUGPRINT("Now after merging peak\n");

#endif

//		goto recycling; // go and restart again .
	}


	return numcore;
}

void Self_Halo_Potent(int np,Vector3d *r,float *mass,float *penergy){
    TStruct *TREE;
    TPtlStruct *ptl;
    particle p;
    float theta2=0.5;
	dptype Mvir = 0;
    int i,j,k;
	if(np < 100){
		float epsilon2 = EPSILON * EPSILON;
		for(i=0;i<np;i++){
			Mvir += mass[i];
		}
		double halobgpotent= GetBGpotent(Mvir);
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
		for(i=0;i<np;i++){
			float potent = 0;
			for(j=0;j<np;j++){
				float ptlmass = mass[j]; 
				float tmpx = r[i].x -  r[j].x;
				float tmpy = r[i].y -  r[j].y;
				float tmpz = r[i].z -  r[j].z;
				float dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;  
				if(dist2 >0) potent += -ptlmass/sqrt(dist2+epsilon2);
			}
			penergy[i] = potent*potentfact + halobgpotent;
		}
	}
	else { 
		size_t nnode = MAX(65*10000,np*0.5);
		TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE)); 
		ptl = (TPtlStruct *) Malloc(sizeof(TPtlStruct)*np,PPTR(ptl)); 
		for(i=0;i<np;i++){ 
			ptl[i].type = TYPE_PTL; 
			ptl[i].x = r[i].x; ptl[i].y = r[i].y; ptl[i].z = r[i].z; 
			ptl[i].mass = mass[i]; 
			ptl[i].sibling = &ptl[i+1]; 
			Mvir += mass[i]; 
		}
		double halobgpotent= GetBGpotent(Mvir);
		ptl[np-1].sibling = NULL;
	
		int recursiveflag;
		if(nnode > 65*10000) recursiveflag = PTHREAD;
		else recursiveflag = RECURSIVE;

		recursiveflag = SERIALIZED;
		Make_Tree(TREE,nnode,ptl,np,theta2,recursiveflag);
#ifdef _OPENMP
#pragma omp parallel for private(i,p) schedule(guided)
#endif
		for(i=0;i<np;i++){ 
			p.x = ptl[i].x; p.y = ptl[i].y; p.z = ptl[i].z; 
			penergy[i] = treeplumpotential(&p,theta2,TREE,ptl)*potentfact + halobgpotent; 
		} 
		Free(ptl); 
		Free(TREE);
	}
}
void External_Halo_Potent(int nend,Vector3d *r,float *mass, float *penergy, int snp, Vector3d *sr, float *smass){
    TStruct *TREE;
    TPtlStruct *ptl;
    particle p;
    float theta2=0.5;
    int i,j,k;
    static int np;

	if(nend < 1000 ){
		float epsilon2 = EPSILON * EPSILON;
		dptype Mvir = 0;
		for(i=0;i<snp;i++){
			Mvir += smass[i];
		}
		double halobgpotent= GetBGpotent(Mvir);
		if(nend >= MAXTHREADS){
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
			for(i=0;i<nend;i++){
				float potent = 0;
				for(j=0;j<snp;j++){
					float ptlmass = smass[j]; 
					float tmpx = r[i].x -  sr[j].x;
					float tmpy = r[i].y -  sr[j].y;
					float tmpz = r[i].z -  sr[j].z;
					float dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;  
					if(dist2 >0) potent += -ptlmass/sqrt(dist2+epsilon2);
				}
				penergy[i] = potent*potentfact + halobgpotent;
			}
		}
		else {
			for(i=0;i<nend;i++) penergy[i] = 0;
#ifdef _OPENMP
#pragma omp parallel for private(i,j) reduction(+:penergy[:nend])
#endif
			for(j=0;j<snp;j++){
				float ptlmass = smass[j];
				for(i=0;i<nend;i++){
					float potent;
					float tmpx = r[i].x -  sr[j].x;
					float tmpy = r[i].y -  sr[j].y;
					float tmpz = r[i].z -  sr[j].z;
					float dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;  
					if(dist2 >0) potent = -ptlmass/sqrt(dist2+epsilon2);
					penergy[i] += potent*potentfact + halobgpotent;
				}
			}
		}
	}
	else {
		size_t nnode = MAX(65*10000,snp*0.5);
	    TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE)); 
		ptl = (TPtlStruct *) Malloc(sizeof(TPtlStruct)*snp,PPTR(ptl)); 
		dptype Mvir = 0; 
#ifdef _OPENMP
#pragma omp parallel for private(i) reduction(+:Mvir)
#endif
		for(i=0;i<snp;i++){ 
			ptl[i].type = TYPE_PTL; 
			ptl[i].x = sr[i].x; ptl[i].y = sr[i].y; ptl[i].z = sr[i].z; 
			ptl[i].mass = smass[i]; 
			ptl[i].sibling = &ptl[i+1]; 
			Mvir += smass[i]; 
		} 
		ptl[snp-1].sibling = NULL; 
		double halobgpotent =  GetBGpotent(Mvir); 

		int recursiveflag;
		if(nnode > 65*10000) recursiveflag = PTHREAD;
		else recursiveflag = RECURSIVE;
		recursiveflag = SERIALIZED;
		Make_Tree(TREE,nnode,ptl,snp,theta2,recursiveflag);
#ifdef _OPENMP
#pragma omp parallel for private(i,p) schedule(dynamic)
#endif
		for(i=0;i<nend;i++){ 
			p.x = r[i].x; p.y = r[i].y; p.z = r[i].z; 
			penergy[i] = treeplumpotential(&p,theta2,TREE,ptl)*potentfact + halobgpotent; 
		} 
		Free(ptl); 
		Free(TREE);
	}
}
int CheckExternalTE(Kptype *tkp,int tnp,unsigned char *bndflag,
		Kptype *skp,int snp, Coretype *icore){
	int i,j,k,nbound;
	float *kenergy,*penergy;
	double cx,cy,cz,cvx,cvy,cvz;
	float distx,disty,distz;
	float tmpvx,tmpvy,tmpvz;
	Vector3d *tr,*tvr;
	Vector3d *sr;
	float *mass;

	tr = (Vector3d*)Malloc(sizeof(Vector3d)*tnp,PPTR(tr));
	tvr = (Vector3d*)Malloc(sizeof(Vector3d)*tnp,PPTR(tvr));
	mass = (float*)Malloc(sizeof(float)*tnp,PPTR(mass));
	kenergy = (float*)Malloc(sizeof(float)*tnp,PPTR(kenergy));
	penergy = (float*)Malloc(sizeof(float)*tnp,PPTR(penergy));

	/*
	cx = cy = cz = cvx = cvy = cvz = 0.L;
	for(i=0;i<tnp;i++){
		cx += (tr[i].x = tkp[i].x);
		cy += (tr[i].y = tkp[i].y);
		cz += (tr[i].z = tkp[i].z);
		cvx += (tvr[i].x = tkp[i].vx);
		cvy += (tvr[i].y = tkp[i].vy);
		cvz += (tvr[i].z = tkp[i].vz);
	}
	cx = cx/(double)tnp; cy = cy/(double)tnp; cz = cz/(double)tnp;
	cvx = cvx/(double)tnp; cvy = cvy/(double)tnp; cvz = cvz/(double)tnp;
	*/
	for(i=0;i<tnp;i++){
		tr[i].x = tkp[i].x;
		tr[i].y = tkp[i].y;
		tr[i].z = tkp[i].z;
		tvr[i].x = tkp[i].vx;
		tvr[i].y = tkp[i].vy;
		tvr[i].z = tkp[i].vz;
		mass[i] = tkp[i].mass;
	}
	cx = icore->cx; cy = icore->cy; cz = icore->cz;
	cvx = icore->cvx; cvy = icore->cvy; cvz = icore->cvz;

#ifdef _OPENMP
#pragma omp parallel for private(i,distx,disty, distz,tmpvx,tmpvy,tmpvz)
#endif
	for(i=0;i<tnp;i++){
		distx = (tr[i].x-cx)*r1kineticfact;
		disty = (tr[i].y-cy)*r1kineticfact;
		distz = (tr[i].z-cz)*r1kineticfact;
		tmpvx = distx + (tvr[i].x-cvx)*r2kineticfact;
		tmpvy = disty + (tvr[i].y-cvy)*r2kineticfact;
		tmpvz = distz + (tvr[i].z-cvz)*r2kineticfact;
		kenergy[i] = 0.5*(tmpvx*tmpvx+tmpvy*tmpvy+tmpvz*tmpvz);
	}

	sr = (Vector3d*)Malloc(sizeof(Vector3d)*snp,PPTR(sr));
	float *smass = (float*)Malloc(sizeof(float)*snp,PPTR(smass));
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
	for(i=0;i<snp;i++){
		sr[i].x = skp[i].x; sr[i].y = skp[i].y; sr[i].z = skp[i].z;
		smass[i] = skp[i].mass;
	}
	if(tnp>0) External_Halo_Potent(tnp,tr,mass,penergy,snp,sr, smass);
	Free(smass);
	Free(sr);

	nbound = 0;
	for(i=0;i<tnp;i++){
		if(kenergy[i] > -penergy[i]) {
			bndflag[i] = UNBOUND;
		}
		else {
			bndflag[i] = BOUND;
			nbound ++;
		}
	}
	Free(penergy); Free(kenergy);
	Free(mass); 
	Free(tvr);Free(tr);
	return nbound;
}
int MULTIHALOSHELL(Kptype *tkp,int tnp, Kptype *skp, int snp,
		SimpleBasicParticleType *bp,int np,Coretype *icore){
	int i,j,k;
	unsigned char *bndflag;
	int nmem,omem;
	bndflag = (unsigned char *)Malloc(sizeof(unsigned char)*tnp,PPTR(bndflag));
	for(i=0;i<tnp;i++) bndflag[i] = BOUND;
	omem = CheckExternalTE(tkp,tnp,bndflag,skp,snp,icore);
	nmem = 0;
	for(i=0;i<tnp;i++) {
		if(bndflag[i] == BOUND) {
			SET_BOUND(KP2BPID(tkp,i));
			tkp[nmem++] = tkp[i];
		}
		else {
			SET_REMAINING(KP2BPID(tkp,i));
			UNSET_BOUND(KP2BPID(tkp,i));
			/* */
			SET_MEMBER_ID(KP2BPID(tkp,i),NOT_HALO_MEMBER);
		}
	}
	Free(bndflag);
	return nmem;
}

int CheckSelfTE(Kptype *kp,int np,unsigned char *bndflag, Coretype *icore){
	int i,j,k,nbound;
	float *kenergy,*penergy;
	double cx,cy,cz,cvx,cvy,cvz;
	float distx,disty,distz;
	float tmpvx,tmpvy,tmpvz;
	Vector3d *r,*vr;
	float *mass;

	r = (Vector3d*)Malloc(sizeof(Vector3d)*np,PPTR(r));
	vr = (Vector3d*)Malloc(sizeof(Vector3d)*np,PPTR(vr));
	mass = (float*)Malloc(sizeof(float)*np,PPTR(mass));
	kenergy = (float*)Malloc(sizeof(float)*np,PPTR(kenergy));
	penergy = (float*)Malloc(sizeof(float)*np,PPTR(penergy));
	/*
	if(icore !=NULL){
		for(i=0;i<np;i++){
			r[i].x = kp[i].x;
			r[i].y = kp[i].y;
			r[i].z = kp[i].z;
			vr[i].x = kp[i].vx;
			vr[i].y = kp[i].vy;
			vr[i].z = kp[i].vz;
			mass[i] = kp[i].mass;
		}
		cx = icore->cx; cy = icore->cy; cz = icore->cz;
		cvx = icore->cvx; cvy = icore->cvy; cvz = icore->cvz;
	}
	else 
	*/
	{
		for(i=0;i<np;i++){
			r[i].x = kp[i].x;
			r[i].y = kp[i].y;
			r[i].z = kp[i].z;
			vr[i].x = kp[i].vx;
			vr[i].y = kp[i].vy;
			vr[i].z = kp[i].vz;
			mass[i] = kp[i].mass;
		}
		dptype tmass = 0;
		cx = cy = cz = cvx = cvy = cvz = 0.L;
		for(i=0;i<np;i++){
			if(kp[i].bp->type == TYPE_STAR){
				cx +=  r[i].x*kp[i].mass;
				cy +=  r[i].y*kp[i].mass;
				cz +=  r[i].z*kp[i].mass;
				cvx += vr[i].x*kp[i].mass;
				cvy += vr[i].y*kp[i].mass;
				cvz += vr[i].z*kp[i].mass;
				tmass += kp[i].mass;
			}
		}
		cx = cx/tmass; cy = cy/tmass; cz = cz/tmass;
		cvx = cvx/tmass; cvy = cvy/tmass; cvz = cvz/tmass;
	}

	for(i=0;i<np;i++){
		distx = (r[i].x-cx)*r1kineticfact;
		disty = (r[i].y-cy)*r1kineticfact;
		distz = (r[i].z-cz)*r1kineticfact;
		tmpvx = distx + (vr[i].x-cvx)*r2kineticfact;
		tmpvy = disty + (vr[i].y-cvy)*r2kineticfact;
		tmpvz = distz + (vr[i].z-cvz)*r2kineticfact;
		kenergy[i] = 0.5*(tmpvx*tmpvx+tmpvy*tmpvy+tmpvz*tmpvz);
	}
	Self_Halo_Potent(np,r,mass,penergy);
	nbound = 0;
	for(i=0;i<np;i++){
		if(kenergy[i] > -penergy[i]) {
			bndflag[i] = UNBOUND;
		}
		else {
			bndflag[i] = BOUND;
			nbound ++;
		}
	}
	/*
#ifdef DEBUG
	for(j=0;j<np;j++){
		printf("P%d of np = %d is ke %g and pe %g bnd? = %d\n",j,np,
				kenergy[j],penergy[j],bndflag[j]);
	}
#endif
*/
	Free(penergy); Free(kenergy);
	Free(mass);
	Free(vr);Free(r);
	DEBUGPRINT(" Total bound : %d  among %d\n", nbound,np);
	return nbound;
}
int ALONEHALO(Kptype *kp,int np,SimpleBasicParticleType *bp,int haloid){
	int i,j,k;
	unsigned char *bndflag;
	int nmem,onmem;
	bndflag = (unsigned char *)Malloc(sizeof(unsigned char)*np,PPTR(bndflag));
	for(i=0;i<np;i++) bndflag[i] = BOUND;
	onmem = np;
	for(i=0;i<BOUNDITER;i++){
		nmem = 0;
		for(j=0;j<onmem;j++){
			if(bndflag[j] == BOUND){
				kp[nmem++] = kp[j];
			}
		}
		if(nmem==0) {
			Free(bndflag);
			return nmem;
		}
		nmem = CheckSelfTE(kp,(k=nmem),bndflag,NULL);
#ifdef DEBUG
		printf("iterating to find bound particles %d:    %d from %d\n",i,nmem,onmem);
#endif
		/*
		nmem = CheckSelfTE(kp,(k=nmem),bndflag);
		*/
		if(nmem == onmem  || nmem ==0) {
			break;
		}
		onmem = nmem;
	}
	Free(bndflag);
	for(i=0;i<nmem;i++) SET_MEMBER_ID(KP2BPID(kp,i),haloid);
//	for(i=0;i<np;i++) SET_MEMBER_ID(i,haloid);
#ifdef DEBUG
	printf("alone halo loses member particles from %d to %d\n",np,nmem);
#endif
	return nmem;
}
int SINGLEHALO(int nkp,Kptype *kp,int np,SimpleBasicParticleType *bp,int haloid, Coretype *core){
	int i,j,k;
	unsigned char *bndflag;
	int nmem;
	bndflag = (unsigned char*)Malloc(sizeof(unsigned char)*np,PPTR(bndflag));
	for(i=0;i<nkp;i++) bndflag[i] = BOUND;

	nmem = CheckSelfTE(kp,nkp,bndflag,core);

	for(i=0;i<nkp;i++) {
		if(bndflag[i]==BOUND) {
			SET_MEMBER_ID(KP2BPID(kp,i),haloid);
			UNSET_REMAINING(KP2BPID(kp,i));
		}
		else if(wp[KP2BPID(kp,i)].haloid == haloid){
			SET_MEMBER_ID(KP2BPID(kp,i),NOT_HALO_MEMBER);
			SET_REMAINING(KP2BPID(kp,i));
		}
	}
	Free(bndflag);
	return nmem;
}
int coresort(const void *a,const void *b){
	Coresorttype *aa,*bb;
	aa = (Coresorttype *)a;
	bb = (Coresorttype *)b;
	if(aa->nmem < bb->nmem) return -1;
	else if(aa->nmem > bb->nmem) return +1;
	else return 0;
}
typedef struct Memtype{
	double cx,cy,cz;
	Vector3d *List;
	float *mass;
	int num; 
	float tmass;
} Memtype;
/*
Memtype member[40000];
*/
#define Ndist 256
#define Nc 32
#define NMratio 256
float tidalR[NMratio][Ndist][Nc];
int jtflag = 1;
void init_nfw_tidal_table(){
	if(jtflag == 0) return;
	
}
int iflag = 1;

void AdGetTidalRCenterCore(Coretype *core,int numcore,
		Coresorttype *score,int ncore,SimpleBasicParticleType *bp,int np){
	int i,j,k;
#ifdef OLD_TIDAL
	if(iflag ==1) {
		void mkRtidal(void);
		(void)mkRtidal();
		iflag = 0;
	}
#endif

	Memtype *core2member;
	core2member = (Memtype *)Malloc(sizeof(Memtype)*numcore,PPTR(core2member));
	for(i=0;i<numcore;i++) {
		core2member[i].num = 0;
		core2member[i].tmass = 0;
	}
	for(i=0;i<np;i++) if(wp[i].haloid>=0) (core2member[wp[i].haloid].num)++;
//	DEBUGPRINT(" C2 core has %d members\n", core2member[2].num);
	for(i=0;i<numcore;i++){
		if(core2member[i].num ==0) {
			core2member[i].List = NULL;
			core2member[i].mass = NULL;
		}
		else {
			core2member[i].List = (Vector3d*) Malloc(sizeof(Vector3d)*core2member[i].num,PPTR(core2member[i].List));
			core2member[i].mass = (float*) Malloc(sizeof(float)*core2member[i].num,PPTR(core2member[i].mass));
			core2member[i].num = 0;
			core2member[i].tmass = 0;
		}
	}
	for(i=0;i<np;i++)
		if(wp[i].haloid>=0) {
			k = wp[i].haloid;
			((core2member[k].List)+core2member[k].num)->x = bp[i].x;
			((core2member[k].List)+core2member[k].num)->y = bp[i].y;
			((core2member[k].List)+core2member[k].num)->z = bp[i].z;
			*((core2member[k].mass)+core2member[k].num) = bp[i].mass;
			(core2member[k].num)++;
			core2member[k].tmass += bp[i].mass;
		}
#ifdef _OPENMP
#pragma omp parallel private(i,k)
#endif
	{
#ifdef _OPENMP
		int it = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		for(i=it;i<ncore;i+=nthreads)
#else
		for(i=0;i<ncore;i++)
#endif
		{
			int coreID = SCOREID2COREID(i);
			int npCore = core2member[coreID].num;
			core[coreID].Rtidal = MAX_TIDAL_R;
	
			double coreX = core[coreID].cx; double coreY = core[coreID].cy; double coreZ = core[coreID].cz;
			core[coreID].nummem = npCore;
			if(npCore>0){
				int jdcore;
				for(jdcore=i+1;jdcore<ncore;jdcore++){
					int countCoreID = SCOREID2COREID(jdcore);
					int bnpCore = core2member[countCoreID].num;
					double bmass = core2member[countCoreID].tmass;
					double bcoreX = core[countCoreID].cx; double bcoreY = core[countCoreID].cy; double bcoreZ = core[countCoreID].cz;
					double r2 = (bcoreX-coreX)*(bcoreX-coreX)+
						(bcoreY-coreY)*(bcoreY-coreY)+
						(bcoreZ-coreZ)*(bcoreZ-coreZ);
					if(bnpCore>0){
#ifndef OLD_TIDAL
						float Mvir=0;
						int inum = 0;
						for(k=0;k<bnpCore;k++){
							double tmpx = ((core2member[countCoreID].List+k)->x) - bcoreX;
							double tmpy = ((core2member[countCoreID].List+k)->y) - bcoreY;
							double tmpz = ((core2member[countCoreID].List+k)->z) - bcoreZ;
							double dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
							if(dist2< r2) {
								Mvir += core2member[countCoreID].mass[k];
								inum++;
							}
						}
						if(inum >0){
							int mtmp = core2member[coreID].tmass /Mvir * NUM_MASS;
							if(mtmp >= NUM_MASS) {
								core[coreID].Rtidal = sqrt(r2)*0.5; //  middle 
							}
							else {
								double tidal_radius;
								double mratio = core2member[coreID].tmass/Mvir;
								tidal_radius  = get_tidal_ellipse(npCore,coreID,inum,countCoreID,bnpCore,core, Mvir, mratio,bmass);
								tidal_radius = MIN(tidal_radius,0.5);
								core[coreID].Rtidal = MIN(core[coreID].Rtidal, sqrt(r2)*tidal_radius);
							}
						}
						else{
							core[coreID].Rtidal = MIN(core[coreID].Rtidal, MAX_TIDAL_R);
						}
#else
						float mratio = core2member[coreID].tmass/bmass;
						float Rvir = pow(3.*bmass/(4.L*M_PI*vv*RHOC),0.3333333333333L);
						float dist_over_Rv = sqrt(r2)/Rvir;
						float nfw_c = get_concentration(core2member[countCoreID].tmass);
						float nfw_rtidal(float, float, float);
						float tidal_radius = nfw_rtidal(mratio, dist_over_Rv, nfw_c);
						tidal_radius = MIN(tidal_radius,0.5);
						core[coreID].Rtidal = MIN(core[coreID].Rtidal, sqrt(r2)*tidal_radius);
	
#endif
					}
				}
			}
			else {
				core[coreID].Rtidal = 0.;
			}
			DEBUGPRINT("C%d has tidal radius %g for mass = %g\n",coreID,core[coreID].Rtidal,core2member[coreID].tmass);
		}
	}

	for(i=numcore-1;i>=0;i--) {
		Free(core2member[i].List);
		Free(core2member[i].mass);
	}
	Free(core2member);
}
int GetMemList(int np,int haloid, int *plist){
	int i,j,k;
	int nmem,nlcore;
	nmem  = 0;
	for(i=0;i<np;i++) if(wp[i].haloid == haloid) plist[nmem++] = i;
	return nmem;
}

#define getcenter(list,nlist,xx,yy,zz,c,id){\
	int kk;\
	xx=yy=zz=0.L;\
	for(kk=0;kk<nlist;kk++){\
		xx+=bp[list[kk]].x; yy+=bp[list[kk]].y; zz+=bp[list[kk]].z;\
	}\
	xx = xx/nlist;yy=yy/nlist;zz=zz/nlist;\
	c[id].cx = xx; c[id].cy = yy; c[id].cz = zz;\
}
/*
void GetTidalRCenterCore(Coretype *core,int numcore,
		Coresorttype *score,int ncore,SimpleBasicParticleType *bp,int np){
	int i,j,k;
	int sid,bid,inum;
	int *slist,*blist,snp,bnp,mtmp;
	float dist,tmpx,tmpy,tmpz,dist2,r2,r1;
	double scx,scy,scz,bcx,bcy,bcz;
	for(i=0;i<ncore;i++){
		sid = SCOREID2COREID(i);
		slist = (int*)Malloc(sizeof(int)*np,PPTR(slist));
		snp = GetMemList(np,sid, slist);
		slist = Realloc(slist,sizeof(int)*snp);
		getcenter(slist,snp,scx,scy,scz,core,sid);
		core[sid].Rtidal = MAX_TIDAL_R;
		if(snp>0){
			for(j=i+1;j<ncore;j++){
				bid = SCOREID2COREID(j);
				blist = (int*)Malloc(sizeof(int)*np,PPTR(blist));
				bnp = GetMemList(np,bid, blist);
				blist = Realloc(blist,sizeof(int)*bnp);
				getcenter(blist,bnp,bcx,bcy,bcz,core,bid);
				r2 = (bcx-scx)*(bcx-scx)+(bcy-scy)*(bcy-scy)+(bcz-scz)*(bcz-scz);
				if(bnp>0){
					inum = 0;
					float emass = 0;
					for(k=0;k<bnp;k++){
						tmpx = bp[blist[k]].x - bcx;
						tmpy = bp[blist[k]].y - bcy;
						tmpz = bp[blist[k]].z - bcz;
						dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
						if(dist2< r2) {
							inum++;
							emass += bp[blist[k]].mass;
						}
					}
					if(inum >0){
						mtmp = (int)((float) snp / (float) inum * (float) NUM_MASS);
						if(mtmp >= NUM_MASS) {
						}
						else 
							core[sid].Rtidal = MIN(core[sid].Rtidal, sqrt(r2)*r_tidal[mtmp]);
					}
					else{
						core[sid].Rtidal = MIN(core[sid].Rtidal, MAX_TIDAL_R);
					}

				}
				Free(blist);
			}
		}
#ifdef DEBUG
		printf("C%d has tidal radius %g in r1=%g\n",sid,core[sid].Rtidal,sqrt(r2));
#endif
		Free(slist);
	}
}
*/
void UnboundShellP2Rest(int ishell,SimpleBasicParticleType *bp){
	int i,j,k;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<NSHELL2P(ishell);i++) 
		if(IS_BOUND(SHELL2P(ishell)[i]) == NOT) {
			SET_MEMBER_ID(SHELL2P(ishell)[i],NOT_HALO_MEMBER);
			SET_REMAINING(SHELL2P(ishell)[i]);
		}
}

#define DumpKPfromBP(kp,list,nlist){\
	int ii;\
	for(ii=0;ii<nlist;ii++){\
		kp[ii].x = bp[list[ii]].x;\
		kp[ii].y = bp[list[ii]].y;\
		kp[ii].z = bp[list[ii]].z;\
		kp[ii].vx = bp[list[ii]].vx;\
		kp[ii].vy = bp[list[ii]].vy;\
		kp[ii].vz = bp[list[ii]].vz;\
		kp[ii].mass = bp[list[ii]].mass;\
		kp[ii].bp = bp+list[ii];\
	}\
}

void PDumpKPfromBP(Kptype *kp, int *list, int nlist,SimpleBasicParticleType *bp){
	int ii;
#ifdef _OPENMP
#pragma omp parallel for private(ii) schedule(guided)
#endif
	for(ii=0;ii<nlist;ii++){
		kp[ii].x = bp[list[ii]].x;
		kp[ii].y = bp[list[ii]].y;
		kp[ii].z = bp[list[ii]].z;
		kp[ii].vx = bp[list[ii]].vx;
		kp[ii].vy = bp[list[ii]].vy;
		kp[ii].vz = bp[list[ii]].vz;
		kp[ii].mass = bp[list[ii]].mass;
		kp[ii].bp = bp+list[ii];
	}
}
int findUpperLink(ompFoFParticleType *ptl,int i){
	if(ptl[i].imother != i) ptl[i].imother = findUpperLink(ptl,ptl[i].imother);
	return ptl[i].imother;
}
int withinFoFRange(ompFoFParticleType *ptl, int i,int j){
	POSTYPE tmpx = ptl[i].x-ptl[j].x;
	POSTYPE tmpy = ptl[i].y-ptl[j].y;
	POSTYPE tmpz = ptl[i].z-ptl[j].z;
	POSTYPE dist = sqrt(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);
	if(dist < 0.5*( ptl[i].link02 + ptl[j].link02)){
		return 1;
	}
	else {
		return 0;
	}
	return 0;
}
void uniteFoF(ompFoFParticleType *ptl, int i,int j){
	int root_i = findUpperLink(ptl, i);
	int root_j = findUpperLink(ptl, j);
	if(root_i != root_j) {
		if(root_i < root_j) ptl[root_j].imother = root_i;
		else ptl[root_i].imother = root_j;
	}
}


typedef struct ompLinkedList{
	int bpid;
	particle p;
} ompLinkedList;


void ompEnabledMemberFoF(SimpleBasicParticleType *bp, int np, int numcore, 
		Coretype *core, int ptype){
	FoFTPtlStruct *ptl;
	int i,j,k;
	int icore;
	float fof_link = FOFLINK4MEMBERSHIP;
	ptl = (FoFTPtlStruct*) Malloc(sizeof(FoFTPtlStruct)*np, PPTR(ptl));
	for(icore=0;icore<numcore;icore++){
		int num = 0;
		if(ptype == TYPE_STAR){
			for(j=0;j<np;j++){
				if(wp[j].haloid == icore && bp[j].type == TYPE_STAR){
					ptl[num].x = bp[j].x;
					ptl[num].y = bp[j].y;
					ptl[num].z = bp[j].z;
					ptl[num].link02 = bp[j].link02;
					ptl[num].indx = j;
					ptl[num].haloindx = -1;
					num++;
				}
			}
		}
		else {
			for(j=0;j<np;j++){
				if(wp[j].haloid == icore){
					ptl[num].x = bp[j].x;
					ptl[num].y = bp[j].y;
					ptl[num].z = bp[j].z;
					ptl[num].link02 = bp[j].link02;
					ptl[num].indx = j;
					ptl[num].haloindx = -1;
					num++;
				}
			}
		}
		if(num == 0) continue;
        else if(num == 1) {
            for(j=0;j<num;j++){
                SET_MEMBER_ID(ptl[j].indx,NOT_HALO_MEMBER);
                UNSET_BOUND(ptl[j].indx);
                SET_REMAINING(ptl[j].indx);
			}
			continue;
		}
		size_t nnode = MAX(65*10000,num/2);
		FoFTStruct *TREE;
		particle *linked;
		DEBUGPRINT("Making FoF Tree in the parallel mode for C%d with np= %d & nnode= %zu currentMemStack= %lld\n", 
				icore, num, nnode, CurMemStack());
		TREE = (FoFTStruct *)Malloc(sizeof(FoFTStruct)*nnode,PPTR(TREE));
		linked = (particle *)Malloc(sizeof(particle)*num,PPTR(linked));

		int nlink, ilink;
		nlink = 0;
		int haloindx = 0;
#ifdef _OPENMP
#pragma omp parallel private(i,j,k)
#endif
		{
			int pid = omp_get_thread_num();
			int nthreads = omp_get_num_threads();
			int nstep = (num+nthreads-1)/nthreads;
			int istart = nstep * pid;
			int iend = nstep *(pid+1);
			iend = MIN(iend, num);
			int mynp = (iend-istart);
			if(mynp>0){
				FoFTPtlStruct *myptl = ptl + istart;
				/*
				for(i=0;i<mynp;i++){
					myptl[i].sibling =  myptl+i+1;
					myptl[i].included =  NO;
				}
				myptl[mynp-1].sibling = NULL;
				*/
				for(i=pid;i<num;i+=nthreads){
					ptl[i].sibling = ptl+i+nthreads;
					ptl[i].included = NO;
				}
				ptl[i-nthreads].sibling = NULL;

				int nodestep = (nnode+nthreads-1)/nthreads;
				FoFTStruct *mytree = TREE + nodestep*pid;
				int recursiveflag = SERIALIZED;
				/*
				FoF_Make_Tree(mytree,nodestep, myptl,mynp,recursiveflag);
				*/
				FoF_Make_Tree(mytree,nodestep, ptl+pid,mynp,recursiveflag);
			}
		}

		int iloop=0;
		int jjob=0;
		int nmaxlink = 0;
		do {
			ilink = 0;
			if(iloop ==0){
					linked[0].x = bp[core[icore].peak].x;
					linked[0].y = bp[core[icore].peak].y;
					linked[0].z = bp[core[icore].peak].z;
					linked[0].link02 = bp[core[icore].peak].link02;
					nlink = 0; ilink = 1;
			}
			else {
				for(j=0;j<num;j++){
					if(ptl[j].included == NO){
						linked[0].x = ptl[j].x;
						linked[0].y = ptl[j].y;
						linked[0].z = ptl[j].z;
						linked[0].link02 = ptl[j].link02;
						ptl[j].included = YES;
						ptl[j].haloindx = haloindx;
						nlink = 0; ilink = 1;
						jjob = j;
						break;
					}
				}
			}
			iloop =1;
			if(ilink == 0) break;
			while(ilink){
#ifdef _OPENMP
#pragma omp parallel  private(j)
#endif
				{
					int pid = omp_get_thread_num();
		            int nthreads = omp_get_num_threads();
					int nstep = (num+nthreads-1)/nthreads;
					int nodestep = (nnode+nthreads-1)/nthreads;
					FoFTStruct *mytree = TREE + nodestep*pid;
					int istart = nstep * pid;
					int iend = nstep *(pid+1);
					iend =MIN(iend,num);
					int mynp = (iend-istart);
					if(mynp>0){
						FoFTPtlStruct *myptl = ptl + istart;
						for(j=nlink;j<nlink+ilink;j++){
							particle p;
							p.x = linked[j].x;
							p.y = linked[j].y;
							p.z = linked[j].z;
							p.link02 = linked[j].link02;
							destroy_omp_fof_link(&p,fof_link,mytree,myptl);
						}
					}
				}

				nlink += ilink;
				ilink = 0;
				for(j=0;j<num;j++){
					if(ptl[j].included == NEW){
						linked[nlink+ilink].x = ptl[j].x;
						linked[nlink+ilink].y = ptl[j].y;
						linked[nlink+ilink].z = ptl[j].z;
						linked[nlink+ilink].link02 = ptl[j].link02;
						ptl[j].included = YES;
						ptl[j].haloindx = haloindx;
						ilink ++;
					}
				}
				DEBUGPRINT("C%d has %d group of ilink= %d nlink= %d for num= %d pxyz= %g %g %g\n", icore, haloindx, ilink,nlink, num,
						bp[core[icore].peak].x,
						bp[core[icore].peak].y,
						bp[core[icore].peak].z
						);
			}
			haloindx ++;
			nmaxlink = MAX(nmaxlink, nlink);
			if(nmaxlink>= (num-jjob)/2 ) {
				nlink += ilink;
				break;
			}
		} while(1);



		int jmax = 0;
		int maxcount=0;
		for(j = 0;j<haloindx;j++) {
			int count = 0;
			int k;
			for(i=0;i<num;i++){
				if(ptl[i].haloindx == j) count ++;
			}
			if(count > maxcount){
				jmax = j;
				maxcount = count;
			}
			if(count > num/2) {
				break;
			}
		}
		for(j=0;j<num;j++){
			if(ptl[j].haloindx != jmax) {
				SET_MEMBER_ID(ptl[j].indx,NOT_HALO_MEMBER);
                UNSET_BOUND(ptl[j].indx);
                SET_REMAINING(ptl[j].indx);
			}
		}
		if(ptype == TYPE_STAR) 
			DEBUGPRINT("C%d 's # of star members changes from %d to %d\n",icore,num,maxcount);
		else
			DEBUGPRINT("C%d 's # of AllType members changes from %d to %d\n",icore,num,maxcount);
		Free(linked);
		Free(TREE);
	}
	Free(ptl);
}

void MemberStarFoF(SimpleBasicParticleType *bp,int np, int numcore, Coretype *core){
	float fof_link;
	int i,j,k;
	int num;
	FoFTPtlStruct *ptl;
	FoFTStruct *TREE;
	particle *linked,p;

	fof_link = FOFLINK4MEMBERSHIP;
	ptl = (FoFTPtlStruct *) Malloc(sizeof(FoFTPtlStruct)*np,PPTR(ptl));
	linked = (particle *)Malloc(sizeof(particle)*np,PPTR(linked));
	size_t nnode = MAX(65*10000,np);
	TREE = (FoFTStruct *)Malloc(sizeof(FoFTStruct)*nnode,PPTR(TREE));
	for(i=0;i<numcore;i++){
		num = 0;
		for(j=0;j<np;j++){
			if(wp[j].haloid==i && bp[j].type == TYPE_STAR){
				ptl[num].type = TYPE_PTL;
				ptl[num].x = bp[j].x;
				ptl[num].y = bp[j].y;
				ptl[num].z = bp[j].z;
				ptl[num].link02 = bp[j].link02;
				ptl[num].indx = j;
				ptl[num].sibling = &ptl[num+1];
				ptl[num].included = NO;
				num++;
			}
		}
		if(num == 0) continue;
		else if(num == 1) {
			for(j=0;j<num;j++){
				SET_MEMBER_ID(ptl[j].indx,NOT_HALO_MEMBER);
				UNSET_BOUND(ptl[j].indx);
				SET_REMAINING(ptl[j].indx);
			}
			continue;
		}
		ptl[num-1].sibling = NULL;
		int recursiveflag;
		/*
		if(nnode>65*10000) {
			recursiveflag = PTHREAD;
		}
		else {
			recursiveflag = RECURSIVE;
		}
		*/
		DEBUGPRINT("C%d is now doing the stellar" 
				" FoF with num= %d with nnode= %zu\n", 
				i,num,nnode);
		recursiveflag = SERIALIZED;
		FoF_Make_Tree(TREE,nnode, ptl,num,recursiveflag);

		p.x = bp[core[i].peak].x;
		p.y = bp[core[i].peak].y;
		p.z = bp[core[i].peak].z;
		p.link02 = bp[core[i].peak].link02;
		int nlink = destroy_new_fof_link(&p,fof_link,TREE,ptl,linked);
		DEBUGPRINT("C%d has nlink %d in nstar= %d : %g %g %g link02=%g\n", 
				i, nlink, num, p.x,p.y,p.z,p.link02);
		if(nlink < num*0.5) {
			int imax=0,mlink=0,ilink;
			/* This should be checked */
			for(k=0;k<num;k++) {
				ptl[k].sibling = &ptl[k+1];
				ptl[k].included = NO;
			}
			ptl[num-1].sibling = NULL;
		    FoF_Make_Tree(TREE,nnode, ptl,num,recursiveflag);
			/* This should be checked */
            int run_nlink = 0;
			int residual;



			for(j=0;j<num;j++){

				if(ptl[j].included == YES) continue;

				p.x = ptl[j].x;
				p.y = ptl[j].y;
				p.z = ptl[j].z;
				p.link02 = ptl[j].link02;
				/*
				for(k=0;k<num;k++) ptl[k].included = NO;
				*/

				ilink = destroy_new_fof_link(&p,fof_link,TREE,ptl,linked);
#ifndef OLD
                run_nlink += ilink;
				residual = num-run_nlink;
               
				/* update the maximum linked group */
				if(ilink >= mlink) {
					mlink = ilink;
					imax = j;
				}
				if(residual < mlink){
					break;
				}
#else
				if(ilink >= num*0.5) {
					mlink = ilink;
					imax = j;
					break;
				}
				else if(ilink > mlink){
					mlink = ilink;
					imax = j;
				}
#endif
			}
			if(1){
				p.x = ptl[imax].x;
				p.y = ptl[imax].y;
				p.z = ptl[imax].z;
				p.link02 = ptl[imax].link02;
				for(k=0;k<num;k++) {
					ptl[k].sibling = &ptl[k+1];
					ptl[k].included = NO;
				}
				ptl[num-1].sibling = NULL;
				FoF_Make_Tree(TREE,nnode, ptl,num,recursiveflag);

				ilink = destroy_new_fof_link(&p,fof_link,TREE,ptl,linked);
			}
		}

		k = num;
		for(j=0;j<num;j++){
			if(ptl[j].included == NO){
				SET_MEMBER_ID(ptl[j].indx,NOT_HALO_MEMBER);
				UNSET_BOUND(ptl[j].indx);
				SET_REMAINING(ptl[j].indx);
				k--;
			}
		}
		DEBUGPRINT("C%d 's # of star members changes from %d to %d\n",i,num,k);
	}
	Free(TREE);Free(linked);Free(ptl);
}
void MemberFoF(SimpleBasicParticleType *bp,int np, int numcore, Coretype *core){
	float fof_link;
	int i,j,k;
	int num;
	FoFTPtlStruct *ptl;
	FoFTStruct *TREE;
	particle *linked,p;

	fof_link = FOFLINK4MEMBERSHIP;
	ptl = (FoFTPtlStruct *) Malloc(sizeof(FoFTPtlStruct)*np,PPTR(ptl));
	linked = (particle *)Malloc(sizeof(particle)*np,PPTR(linked));
	size_t nnode = MAX(65*10000,np);
	TREE = (FoFTStruct *)Malloc(sizeof(FoFTStruct)*nnode,PPTR(TREE));
	for(i=0;i<numcore;i++){
		num = 0;
		for(j=0;j<np;j++){
			if(wp[j].haloid==i){
				ptl[num].type = TYPE_PTL;
				ptl[num].x = bp[j].x;
				ptl[num].y = bp[j].y;
				ptl[num].z = bp[j].z;
				ptl[num].link02 = bp[j].link02;
				ptl[num].indx = j;
				ptl[num].sibling = &ptl[num+1];
				ptl[num].included = NO;
				num++;
			}
		}
		if(num == 0) continue;
		else if(num == 1) {
			for(j=0;j<num;j++){
				SET_MEMBER_ID(ptl[j].indx,NOT_HALO_MEMBER);
				UNSET_BOUND(ptl[j].indx);
				SET_REMAINING(ptl[j].indx);
			}
			continue;
		}
		ptl[num-1].sibling = NULL;

		int recursiveflag;
		/*
		if(nnode>65*10000) {
			recursiveflag = PTHREAD;
		}
		else {
			recursiveflag = RECURSIVE;
		}
		*/
		recursiveflag = SERIALIZED;
		FoF_Make_Tree(TREE,nnode,ptl,num,recursiveflag);

		p.x = bp[core[i].peak].x;
		p.y = bp[core[i].peak].y;
		p.z = bp[core[i].peak].z;
		p.link02 = bp[core[i].peak].link02;
		int nlink = destroy_new_fof_link(&p,fof_link,TREE,ptl,linked);
		if(nlink < num*0.5) {
			int imax=0,mlink=0,ilink;
			/* This should be checked */
			for(k=0;k<num;k++) {
				ptl[k].sibling = &ptl[k+1];
				ptl[k].included = NO;
			}
			ptl[num-1].sibling = NULL;
			FoF_Make_Tree(TREE,nnode,ptl,num,recursiveflag);
			/* This should be checked */


            int run_nlink = 0;
			int residual;

			for(j=0;j<num;j++){

				if(ptl[j].included == YES) continue;

				p.x = ptl[j].x;
				p.y = ptl[j].y;
				p.z = ptl[j].z;
				p.link02 = ptl[j].link02;
				/*
				for(k=0;k<num;k++) ptl[k].included = NO;
				*/

				ilink = destroy_new_fof_link(&p,fof_link,TREE,ptl,linked);
#ifdef DEBUG
#endif

#ifndef OLD
                run_nlink += ilink;
				residual = num-run_nlink;
               
				/* update the maximum linked group */
				if(ilink >= mlink) {
					mlink = ilink;
					imax = j;
				}
				if(residual < mlink){
					break;
				}
#else
				if(ilink >= num*0.5) {
					mlink = ilink;
					imax = j;
					break;
				}
				else if(ilink > mlink){
					mlink = ilink;
					imax = j;
				}
#endif
			}
			if(mlink <= 5){
				for(j=0;j<num;j++) ptl[j].included = NO;
			}
			else {
				p.x = ptl[imax].x;
				p.y = ptl[imax].y;
				p.z = ptl[imax].z;
				p.link02 = ptl[imax].link02;
				for(k=0;k<num;k++) {
					ptl[k].sibling = &ptl[k+1];
					ptl[k].included = NO;
				}
				ptl[num-1].sibling = NULL;
				FoF_Make_Tree(TREE,nnode,ptl,num,recursiveflag);
				ilink = destroy_new_fof_link(&p,fof_link,TREE,ptl,linked);
			}
		}

		k = num;
		for(j=0;j<num;j++){
			if(ptl[j].included == NO){
				SET_MEMBER_ID(ptl[j].indx,NOT_HALO_MEMBER);
				UNSET_BOUND(ptl[j].indx);
				SET_REMAINING(ptl[j].indx);
				k--;
			}
		}
		DEBUGPRINT("C%d 's # of total members changes from %d to %d\n",i,num,k);
	}
	Free(TREE);Free(linked);Free(ptl);
}


int	AdGetMemberRestParticleList(SimpleBasicParticleType *bp,int np,int hid,int *list,
		Coretype *core, Coresorttype *score,int now,int numcore){
	int i,j,k,num,onum,inum;
	int id;
	float cx,cy,cz;
	float rt2;
	float tmpx,tmpy,tmpz,dist2;
	cx = core[hid].cx; cy = core[hid].cy; cz = core[hid].cz;
	rt2 = (core[hid].Rtidal)*(core[hid].Rtidal);
	
	num = 0;
	for(i=0;i<np;i++){
		if(wp[i].haloid==hid){
			list[num++] = i;
		}
	}
	onum = num;
	for(i=0;i<np;i++){
		if(wp[i].haloid != hid){
			tmpx =  bp[i].x-cx; tmpy =  bp[i].y-cy; tmpz =  bp[i].z-cz;
			dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
			if(dist2 < rt2) list[num++] = i;
		}
	}
    inum = onum;
	for(i=onum;i<num;i++){
		id = wp[list[i]].haloid;
		if(id==NOT_HALO_MEMBER) {
			list[inum++] = list[i];
		}
		else {
			for(j=now+1;j<numcore;j++){
				if(id == SCOREID2COREID(j)){
					list[inum++] = list[i];
					break;
				}
			}
		}
	}
	num = inum;
	return num;
}
int	GetMemberRestParticleList(SimpleBasicParticleType *bp,int np,int hid,int *list,
		Coretype *core, Coresorttype *score,int now,int numcore){
	int i,j,k,num,onum,inum;
	int id;
	float cx,cy,cz;
	float rt2;
	float tmpx,tmpy,tmpz,dist2;
	cx = core[hid].cx; cy = core[hid].cy; cz = core[hid].cz;
	rt2 = (core[hid].Rtidal)*(core[hid].Rtidal);
	
	num = 0;
	for(i=0;i<np;i++){
		tmpx =  bp[i].x-cx; tmpy =  bp[i].y-cy; tmpz =  bp[i].z-cz;
		dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
		/* Include all the particles within the tidal radius irrepective of
		 * whether they are core particles to other halo candidates. */
		if(dist2 < rt2) {
			list[num++] = i;
		}
		/* Erase from member list halo particles that are beyond the tidal radius. */
		else if(wp[i].haloid == hid){
			SET_MEMBER_ID(i,NOT_HALO_MEMBER);
			SET_REMAINING(i);
		}
	}
	inum = 0;
	for(i=0;i<num;i++){
		id = wp[list[i]].haloid;
		/* If it is the member particle or free, simply include it to the list. */
		if(id==NOT_HALO_MEMBER || id == hid) {
			list[inum++] = list[i];
		}
		/* If it is not a core particle of another halo and it is 
		 * a member particle of another larger system, then include it. */
		else if(IS_CORE(list[i]) == NOT) {
			for(j=now+1;j<numcore;j++){
				if(id == SCOREID2COREID(j)){
					list[inum++] = list[i];
					break;
				}
			}
		}
	}
	num = inum;
	return num;
}
/*
int	GetMemberRestParticleList(SimpleBasicParticleType *bp,int np,int hid,int *list,
		Coretype *core, Coresorttype *score,int now,int numcore){
	int i,j,k,num,onum,inum;
	int id;
	float cx,cy,cz,rt2;
	float tmpx,tmpy,tmpz,dist2;
	num = 0;
	cx = core[hid].cx; cy = core[hid].cy; cz = core[hid].cz;
	rt2 = core[hid].Rtidal*core[hid].Rtidal;
	for(i=0;i<np;i++){
		if(wp[i].haloid==hid || IS_REMAINING(i) != NOT){
			list[num++] = i;
		}
	}
	onum = num;
	for(i=0;i<np;i++){
		if(wp[i].haloid!=hid && IS_REMAINING(i) == NOT){
			tmpx =  bp[i].x-cx; tmpy =  bp[i].y-cy; tmpz =  bp[i].z-cz;
			dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
			if(dist2 <= rt2) list[num++] = i;
		}
	}
	inum = onum;
	for(i=onum;i<num;i++){
		id = wp[list[i]].haloid;
		if(id==NOT_HALO_MEMBER) continue;
		for(j=now+1;j<numcore;j++){
			if(id == SCOREID2COREID(j)){
				list[inum++] = list[i];
				break;
			}
		}
	}
	num = inum;
	return num;
}
*/


#define min(a,b) (a)<(b)? (a):(b)

double Dump2BP(FoFTPtlStruct *rbp,SimpleBasicParticleType *bp,lint np){
	lint i;
	double xmin,ymin,zmin;
	double tmass = 0;
	xmin = ymin = zmin =  1.E30L;
	for(i=0;i<np;i++){
		xmin = MIN(xmin,rbp[i].x);
		ymin = MIN(ymin,rbp[i].y);
		zmin = MIN(zmin,rbp[i].z);
	}
	Xmin = xmin;
	Ymin = ymin;
	Zmin = zmin;
	for(i=0;i<np;i++){
		bp[i].type = rbp[i].type;
		bp[i].mass = rbp[i].mass;
		bp[i].x = rbp[i].x - xmin;
		bp[i].y = rbp[i].y - ymin;
		bp[i].z = rbp[i].z - zmin;
		bp[i].vx = rbp[i].vx;
		bp[i].vy = rbp[i].vy;
		bp[i].vz = rbp[i].vz;
		bp[i].indx = rbp[i].indx;
		bp[i].link02 = rbp[i].link02;
		tmass += rbp[i].mass;
	}
	return tmass;
}
float findstarmass(SimpleBasicParticleType *bp, int np){
    float starmass=0;
    int i,j,k;
    for(i=0;i<np;i++){
        if(bp[i].type == TYPE_STAR)
            starmass += bp[i].mass;
    }
    return starmass;

}
int findstarnum(SimpleBasicParticleType *bp, int np){
    int starmass=0;
    int i,j,k;
    for(i=0;i<np;i++){
        if(bp[i].type == TYPE_STAR)
            starmass += 1;
    }
    return starmass;

}



int subhalo_den(FoFTPtlStruct *rbp, lint np,lint *p2halo){
	int i,j,k,bid;
	float xinit,yinit,zinit;
	float xmax,ymax,zmax;
	long long *neighbor;
	int NumNeighbor,numcore;
	float dthreshold;
	float *density;
	Kptype *kp,*skp,*tkp;
	Coretype *core;
	SimpleBasicParticleType *bp;

	bp = (SimpleBasicParticleType*)Malloc(sizeof(SimpleBasicParticleType)*np,PPTR(bp));
	double Mvir = Dump2BP(rbp,bp,np);
	double halobgpotent = GetBGpotent(Mvir); 

	for(i=0;i<np;i++){
		p2halo[i] = NOT_HALO_MEMBER;
	}
//	NumNeighbor = MIN(MINCORENMEM,np);
	NumNeighbor = MIN(NUMNEIGHBOR,np);
	if(1){
		mklocalize(bp,np,&xinit,&yinit, &zinit,&xmax,&ymax,&zmax);
	}
	{
		if(findstarnum(bp,np)<= NUMNEIGHBOR || findstarmass(bp,np)<MINSTELLARMASS){
			/*
			neighbor = (int*)Malloc(sizeof(int)*np*NumNeighbor,PPTR(neighbor));
			density = (float*)Malloc(sizeof(float)*np,PPTR(density));
			core = (Coretype*)Malloc(sizeof(Coretype)*maxnumcore,PPTR(core));

			void lagFindStellarCore(SimpleBasicParticleType *, int, int, float *, 
					Coretype **, int *, int, int **, int);
			lagFindStellarCore(bp,np,NumNeighbor,density, &core, &numcore, maxnumcore,
					&neighbor, TYPE_ALL);
					*/
			neighbor = (long long*)Malloc(sizeof(long long)*np*NumNeighbor,PPTR(neighbor));
			density = (float*)Malloc(sizeof(float)*np,PPTR(density));
			void findsphdensity(SimpleBasicParticleType *,int ,long long *, int , float *);
			findsphdensity(bp,np,neighbor,NumNeighbor,density);
			core = (Coretype*)Malloc(sizeof(Coretype)*maxnumcore,PPTR(core));
			numcore = finddenpeak(density,NumNeighbor,neighbor,np,&core,0,bp);
		}
		else {
			neighbor = (long long*)Malloc(sizeof(long long)*np*NumNeighbor,PPTR(neighbor));
			density = (float*)Malloc(sizeof(float)*np,PPTR(density));
			core = (Coretype*)Malloc(sizeof(Coretype)*maxnumcore,PPTR(core));
#ifdef ADV
			void lagFindStellarCore(SimpleBasicParticleType *, int, int, float *, 
					Coretype **, int *, int, long long *, int);
			lagFindStellarCore(bp,np,NumNeighbor,density, &core, &numcore, maxnumcore,
					neighbor, TYPE_STAR);
#else
			neighbor = (long long*)Malloc(sizeof(long long)*np*(long)NumNeighbor,PPTR(neighbor));
			void starfindsphdensity(SimpleBasicParticleType *,int ,long long *, int , float *);
			starfindsphdensity(bp,np,neighbor,NumNeighbor,density);
			numcore = finddenpeak(density,NumNeighbor,neighbor,np,&core,1,bp);
#endif
		}
		DEBUGPRINT("density calculates\n");
	}
	{
		DEBUGPRINT("%d numcore detected\n",numcore);
		for(i=0;i<numcore;i++){
			DEBUGPRINT("C%d has ipeak= %d\n", i, core[i].peak);
		}
		wp = (WorkingParticle *)Malloc(sizeof(WorkingParticle)*np,PPTR(wp));
	}
	if(numcore == 0) {
		/*
#ifdef NOBACKGROUND
		for(i=0;i<np;i++){
			p2halo[i] = 0;
		}
#endif
		Free(wp);Free(density);Free(neighbor);
		Free(bp);
		return 0;
		*/
		for(i=0;i<np;i++) p2halo[i] = 0;
		return 1;
	}
renumcore :
	if(numcore ==1) {
		int nmem;
		/*
		Free(neighbor);

		kp = (Kptype*)Malloc(sizeof(Kptype)*np,PPTR(kp));
		for(i=0;i<np;i++){
			kp[i].x = bp[i].x; kp[i].y = bp[i].y; kp[i].z = bp[i].z;
			kp[i].vx = bp[i].vx; kp[i].vy = bp[i].vy; kp[i].vz = bp[i].vz;
			kp[i].mass = bp[i].mass;
			kp[i].bp = bp+i;
		}
		for(i=0;i<np;i++) SET_MEMBER_ID(i,NOT_HALO_MEMBER);
		nmem = ALONEHALO(kp,np,bp,0);
		Free(kp);

		MemberFoF(bp,np,numcore,core);
		for(i=0;i<np;i++){
#ifdef NOBACKGROUND
			p2halo[i] = 0;
#else
			if(wp[i].haloid==0) p2halo[i] = 0;
			else p2halo[i] = NOT_HALO_MEMBER;
#endif
		}

		Free(wp);
		*/
		for(i=0;i<np;i++) p2halo[i] = 0;
		return 1;
	}
	else if(numcore > 1) {
		float minshellden,maxshellden;
		int ishell,nshelldivide;

		Dump2WorkingParticle(density,np,core,numcore);

		numcore = FindCoreDensity(bp,np,neighbor,NumNeighbor,core,numcore);
		Free(density);
		if(numcore ==1) {
			goto renumcore;
		}
		else if(numcore ==0) {
			for(i=0;i<np;i++) p2halo[i] = 0;
			return 1;
			/*
			Free(wp);Free(density);Free(neighbor);
			Free(bp);
			return numcore;
			*/
		}
		DEBUGPRINT("total number of cores changes to %d\n",numcore);fflush(stdout);

		minshellden = 1.e27;
		maxshellden = -1.e27;
		for(i=0;i<numcore;i++){ /* set the maximum value with the maximum core density */
			maxshellden = MAX(maxshellden,core[i].coredensity);
		}
		for(i=0;i<np;i++){ /* set the minium value with the minimum particle density */
			minshellden = MIN(minshellden,wp[i].den);
		}
		minshellden = MAX(DENFLOOR, minshellden);
		minshellden = log10(minshellden+1.);
		maxshellden = log10(maxshellden+1.);

		nshell = 0;
		if(np>1000) {
			nshelldivide = MAX(NSHELLDIVIDE,NSHELLDIVIDE*log10((double)np*1.5));
			nshelldivide = MIN(15, nshelldivide);
//			nshelldivide = MAX(NSHELLDIVIDE,NSHELLDIVIDE*log10((double)np*2));
		}
		else nshelldivide = NSHELLDIVIDE;

//		numstack = CurMemStack();

		for(i=0;i<nshelldivide;i++){
			dthreshold = minshellden + (nshelldivide-1-i)*
				(maxshellden-minshellden)/(float)nshelldivide;
			dthreshold = pow(10.,dthreshold)-1.;
			GetShellParticlesPeaks(np,neighbor,NumNeighbor,core,numcore,dthreshold);
		}
		SaveRemainingParticles2LastShell(np,core,numcore);
#ifdef DEBUG
		for(i=0;i<nshell;i++)
			DEBUGPRINT("S%d has %d particles & %d cores\n",i, NSHELL2P(i),NSHELL2C(i));
#endif
		/*
		goto gogo;
		*/
		for(ishell=0;ishell<nshell;ishell++){
//			DEBUGPRINT("C0 has posxyz= %g %g %g\n", bp[core[0].peak].x, bp[core[0].peak].y, bp[core[0].peak].z);
			int icore;
			int *halonmem;
			float tradius;
			Coresorttype *score;
			score = (Coresorttype*)Malloc(sizeof(Coresorttype)*NSHELL2C(ishell),PPTR(score));

			halonmem = (int *) Malloc(sizeof(int)*numcore,PPTR(halonmem));
			for(j=0;j<numcore;j++) halonmem[j] = 0;
			for(j=0;j<np;j++)
				if(wp[j].haloid>=0) halonmem[wp[j].haloid]++;
			for(j=0;j<NSHELL2C(ishell);j++) {
				score[j].nmem =  halonmem[SHELL2C(ishell)[j]];
				score[j].core = core+SHELL2C(ishell)[j];
			}
			Free(halonmem);
			qsort(score,NSHELL2C(ishell),sizeof(Coresorttype),coresort);
#ifdef DEBUG
			for(j=0;j<NSHELL2C(ishell);j++) {
				k = SCOREID2COREID(j);
				DEBUGPRINT("C%d has nmem=%d\n",k,core[k].nummem);
			}
#endif
			AdGetTidalRCenterCore(core,numcore, score,NSHELL2C(ishell),bp,np);

			for(j=0;j<NSHELL2C(ishell);j++) {
				int *tlist,ntarget;
				int *slist,nsource,nmem;
				icore = SCOREID2COREID(j);
				tlist = (int*)Malloc(sizeof(int)*np,PPTR(tlist));
				slist = (int*)Malloc(sizeof(int)*np,PPTR(slist));
				ntarget = GetShellParticleFromRestParticleLIST(bp,np,ishell, tlist,score[j].core);
				if(ntarget ==0) {
					Free(slist);Free(tlist);
					continue;
				}
				nsource = GetMemCandidateList(np,icore, slist,ishell);
				skp = (Kptype*)Malloc(sizeof(Kptype)*nsource,PPTR(skp));
				tkp = (Kptype*)Malloc(sizeof(Kptype)*ntarget,PPTR(tkp));
				if(nsource < 10000) {
					DumpKPfromBP(skp,slist,nsource);
				}
				else {
					PDumpKPfromBP(skp,slist,nsource,bp);
				}
				if(ntarget < 10000) {
					DumpKPfromBP(tkp,tlist,ntarget);
				}
				else {
					PDumpKPfromBP(tkp,tlist,ntarget,bp);
				}

				nmem = MULTIHALOSHELL(tkp,ntarget,skp,nsource,bp,np,core+icore);

				for(i=0;i<nmem;i++){
					bid = KP2BPID(tkp,i);
					SET_BOUND(bid);
					UNSET_REMAINING(bid);
					SET_MEMBER_ID(bid,icore);
				}
				Free(tkp);Free(skp);
				Free(slist); Free(tlist);
				DEBUGPRINT("S%d/S%d: C%d ntarget=%d nsource=%d & finally get nmem = %d\n",
						ishell,nshell,icore,ntarget,nsource,nmem);
			}
			UnboundShellP2Rest(ishell,bp);/* Turning on the rest flag for unbound shell particles */
			Free(score);
		}

//		InitialOldMemStack(numstack);
  		for(i=nshell-1;i>=0;i--){
  			Free(SHELL2C(i));
  			Free(SHELL2P(i));
  		}

		if(np <NOMPFoF){
			if(MINSTELLARMASS>=0) MemberStarFoF(bp,np,numcore,core);
			MemberFoF(bp,np,numcore,core);
		}
		else {
			if(MINSTELLARMASS>=0) ompEnabledMemberFoF(bp,np,numcore,core, TYPE_STAR);
			ompEnabledMemberFoF(bp,np,numcore,core, TYPE_ALL);
		}

		for(i=0;i<BOUNDITER;i++) {
			int *tlist,ntarget,jcore,*halonmem;
			int kk,nbound;
#ifdef DEBUG
			DEBUGPRINT("\n\n\n");
			DEBUGPRINT("#### BOUNDITER %d",i);
			DEBUGPRINT("\n\n\n");
#endif
			Coresorttype *score;
			score = (Coresorttype*)Malloc(sizeof(Coresorttype)*numcore,PPTR(score));
			halonmem = (int *) Malloc(sizeof(int)*numcore,PPTR(halonmem));
			for(j=0;j<numcore;j++) halonmem[j] = 0;
			for(j=0;j<np;j++)
				if(wp[j].haloid>=0) halonmem[wp[j].haloid]++;
			for(j=0;j<numcore;j++) {
				score[j].nmem =  halonmem[j];
				score[j].core = core+j;
			}
			Free(halonmem);

			qsort(score,numcore,sizeof(Coresorttype),coresort);
			AdGetTidalRCenterCore(core,numcore, score,numcore,bp,np);
#ifdef DEBUG
			for(j=0;j<numcore;j++){
				jcore = SCOREID2COREID(j);
				DEBUGPRINT("C%d has nmem= %d & Rtidal= %g c %g %g %g\n",
						jcore,score[j].nmem,score[j].core->Rtidal,
						core[jcore].cx,core[jcore].cy,core[jcore].cz);
			}
#endif
			tlist = (int *)Malloc(sizeof(int)*np,PPTR(tlist));
			for(j=0;j<numcore;j++){
				Kptype *kp;
				jcore = SCOREID2COREID(j);
				if(score[j].nmem==0) continue;
				/*
				ntarget = AdGetMemberRestParticleList(bp,np,jcore,tlist,core,score,j,numcore);
				*/
				ntarget = GetMemberRestParticleList(bp,np,jcore,tlist,core,score,j,numcore);
				if(ntarget ==0) continue;
				kp = (Kptype*)Malloc(sizeof(Kptype)*ntarget,PPTR(kp));
				if(ntarget < 100000) {
					DumpKPfromBP(kp,tlist,ntarget);
				}
				else {
					PDumpKPfromBP(kp,tlist,ntarget,bp);
				}
				nbound = SINGLEHALO(ntarget,kp,np,bp,jcore, core+jcore);
				DEBUGPRINT("C%d has bound particles of np= %d from mp= %d within rtidal= %g c %g %g %g\n",
						jcore,nbound, ntarget,core[jcore].Rtidal,
						core[jcore].cx+Xmin,core[jcore].cy+Ymin,core[jcore].cz+Zmin);
				Free(kp);
			}
			Free(tlist); Free(score);
		}


		if(np <NOMPFoF){
			if(MINSTELLARMASS>=0) MemberStarFoF(bp,np,numcore,core);
			MemberFoF(bp,np,numcore,core);
		}
		else {
			if(MINSTELLARMASS>=0) ompEnabledMemberFoF(bp,np,numcore,core, TYPE_STAR);
			ompEnabledMemberFoF(bp,np,numcore,core, TYPE_ALL);
		}

#ifdef DEBUG
		/*
		MPI_Finalize();exit(99);
		*/
#endif
gogo:
		for(i=0;i<np;i++){
			p2halo[i] = wp[i].haloid;
		}

		int *hid_num = (int*)Malloc(sizeof(int)*numcore, PPTR(hid_num));
		int *hid_new = (int*)Malloc(sizeof(int)*numcore, PPTR(hid_new));
		for(i=0;i<numcore;i++){
			hid_num[i] = 0;
			hid_new[i] = -1;
		}
		for(i=0;i<np;i++){
			if(wp[i].haloid >=0) hid_num[wp[i].haloid]++;
		}
		int jpeak = 0;
		int numpeak = 0;
		j= 0;
		for(i=0;i<numcore;i++){
			if(hid_num[i] > numpeak){
				numpeak = hid_num[i];
				jpeak = i;
			}
			if(hid_num[i] >0) {
				hid_new[i] = j;
				j++;
			}
		}
		numcore = j;
		for(i=0;i<np;i++){
			if(p2halo[i]!= NOT_HALO_MEMBER) p2halo[i] = hid_new[p2halo[i]];
#ifdef NOBACKGROUND
			else p2halo[i] = hid_new[jpeak];
#endif
		}


		Free(hid_new);Free(hid_num);



		Free(wp); Free(neighbor);
	}
	/*
	*/

	Free(core);
	Free(bp);
	return numcore;
}
#undef vv
#undef RHOC
#undef Msun
