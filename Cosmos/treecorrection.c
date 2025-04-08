#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
//#include "mpirms.h"
#include "indT.h"
#include "ost.h"
#include "Model4TreeForce.h"
#include "timerutil.h"
#include "CosmosEvolFactor.h"
#include "BasicCell.h"

#define DIRECTSUM 16
/* This value should be determined from the gputest.c */
#ifdef USE_GPU
void gputreeforce();
void Direct_Nbody();
#define GPUCPUDIRECTSUM (25)
int GPUDIRECTSUM;
#endif

#define MIN(a, b) ( (a)>(b)?(b): (a) )
#define MAX(a, b) ( (a)<(b)?(b): (a) )

#ifndef _OPENMP
#define omp_get_thread_num() (0)
#define omp_get_num_threads() (1)
#endif





/*
void BuildLinkedList(SimParameters *simpar, ptrdiff_t mx, ptrdiff_t my, ptrdiff_t mz,
		PosType xmin, PosType ymin, PosType zmin){
	TreeLinkedCell *BasicCell = BASICCELL(simpar);
	ptrdiff_t i,j,k,ix,iy,iz;
	for(i=0;i<mx*my*mz;i++) {
		BasicCell[i].link = NULL;
		BasicCell[i].nmem = 0;
	}
	LinkParticles(simpar, DM,dm, BasicCell,mx,my,xmin,ymin,zmin);
	LinkParticles(simpar, SPH,sph, BasicCell,mx,my,xmin,ymin,zmin);
	LinkParticles(simpar, STAR,star, BasicCell,mx,my,xmin,ymin,zmin);
	LinkParticles(simpar, AGN,agn, BasicCell,mx,my,xmin,ymin,zmin);
}
*/

size_t Dump2Pos(SimParameters *simpar, 
		ptrdiff_t ix,ptrdiff_t iy,ptrdiff_t iz, ptrdiff_t mx, ptrdiff_t my,
		ptrdiff_t mz, Position *pos){
	size_t npos = 0;
#ifndef GOTPM
	int maxTsubpower = IndT_MAXTSUBPOWER(simpar);
	int nowTsubdiv = IndT_NOWTSUBDIV(simpar);
#endif
	TreeLinkedCell *BasicCell = BASICCELL(simpar);
	PosType cellxpos,cellypos, cellzpos;
	cellxpos = CellWidth*(ix - 0.5L) + SIM_LXMIN(simpar,dm);
	cellypos = CellWidth*(iy - 0.5L) + SIM_LYMIN(simpar,dm);
	cellzpos = CellWidth*(iz - 0.5L) + SIM_LZMIN(simpar,dm);
	npos = BasicCell[ix+mx*(iy+my*iz)].nmem;
	if(npos ==0) return npos;
	Position *p = pos;
	npos = 0;
	linkedlisttype *tmp = BasicCell[ix+mx*(iy+my*iz)].link;
	while(tmp) {
		if(!IS_FLAG_ON(tmp,BoundaryGhostflag)){
#ifndef GOTPM
			if(isnowdmstep(nowTsubdiv, GetTsubPower(tmp), maxTsubpower))
#endif
			{
				p[npos].x = XofP(simpar,tmp) - cellxpos;
				p[npos].y = YofP(simpar,tmp) - cellypos;
				p[npos].z = ZofP(simpar,tmp) - cellzpos;
				if( IS_FLAG_ON(tmp, DMflag) ) p[npos].mass = TDM_MASS(simpar,0);
				else p[npos].mass = LINKEDTYPE_MASS(simpar,tmp);
#ifndef GOTPM
				if( IS_FLAG_ON(tmp, DMflag)) p[npos].axyz = &(( (treedmparticletype*)tmp)->ax);
				else if(IS_FLAG_ON(tmp, SPHflag)) p[npos].axyz = &(( (treesphparticletype*)tmp)->ax);
				else if(IS_FLAG_ON(tmp, VOROflag)) p[npos].axyz = &(( (treevoroparticletype*)tmp)->ax);
				else if(IS_FLAG_ON(tmp, STARflag)) p[npos].axyz = &(( (treestarparticletype*)tmp)->ax);
				else if(IS_FLAG_ON(tmp, AGNflag)) p[npos].axyz = &(( (treeagnparticletype*)tmp)->ax);
				else if(IS_FLAG_ON(tmp, VOROflag)) p[npos].axyz = &(( (treevoroparticletype*)tmp)->ax);
#else
				p[npos].axyz = &(((treedmparticletype*)tmp)->vx);
#endif
				npos++;
			}
		}
		tmp = tmp->next;
	}
	return npos;
}
size_t  Dump2GrvP(SimParameters *simpar, ptrdiff_t ix, ptrdiff_t iy,ptrdiff_t iz,
		ptrdiff_t mx,ptrdiff_t my,ptrdiff_t mz,TPtlStruct *grvp){
	TreeLinkedCell *BasicCell = BASICCELL(simpar);
	size_t ngrv;
	/*
	linkedlisttype *sdm, *fdm;
	sdm = (linkedlisttype*)DM_TBP(simpar);
	fdm = (linkedlisttype*)(DM_TBP(simpar) + DM_NP(simpar));
	*/
	PosType cellxpos,cellypos, cellzpos;
	cellxpos = CellWidth*(ix - 0.5L) + SIM_LXMIN(simpar,dm);
	cellypos = CellWidth*(iy - 0.5L) + SIM_LYMIN(simpar,dm);
	cellzpos = CellWidth*(iz - 0.5L) + SIM_LZMIN(simpar,dm);
	ptrdiff_t imin,imax,jmin,jmax,kmin,kmax;
	imin = MAX(0,ix-1);
	jmin = MAX(0,iy-1);
	kmin = MAX(0,iz-1);
	imax = MIN(mx,ix+2);
	jmax = MIN(my,iy+2);
	kmax = MIN(mz,iz+2);
	ptrdiff_t i,j,k;
	ngrv = 0;
	for(k=kmin;k<kmax;k++)for(j=jmin;j<jmax;j++) for(i=imin;i<imax;i++){
		ngrv += BasicCell[i+mx*(j+my*k)].nmem;
	}
	TPtlStruct *g = grvp;
	ngrv = 0;
	for(k=kmin;k<kmax;k++)for(j=jmin;j<jmax;j++) for(i=imin;i<imax;i++){
		linkedlisttype *tmp = BasicCell[i+mx*(j+my*k)].link;
		while(tmp){
			if(IS_FLAG_ON(tmp,BoundaryGhostflag)){
				g[ngrv].x = tmp->x - cellxpos;
				g[ngrv].y = tmp->y - cellypos;
				g[ngrv].z = tmp->z - cellzpos;
			}
			else {
				g[ngrv].x = XofP(simpar, tmp) - cellxpos;
				g[ngrv].y = YofP(simpar, tmp) - cellypos;
				g[ngrv].z = ZofP(simpar, tmp) - cellzpos;
			}
			if( IS_FLAG_ON(tmp, DMflag)) g[ngrv].mass = TDM_MASS(simpar,0);
			else g[ngrv].mass = LINKEDTYPE_MASS(simpar, tmp);
			g[ngrv].indx = PINDX(tmp);
			ngrv ++;
			tmp = tmp->next;
		}
	}
	return ngrv;

}
void free_padding_particles(SimParameters *simpar){
	if(DM_NPP(simpar)) free(DM_TBPP(simpar));
	if(SPH_NPP(simpar)) free(SPH_TBPP(simpar));
	if(VORO_NPP(simpar)) free(VORO_TBPP(simpar));
	if(STAR_NPP(simpar)) free(STAR_TBPP(simpar));
	if(AGN_NPP(simpar)) free(AGN_TBPP(simpar));
	DM_NPP(simpar) = SPH_NPP(simpar) = VORO_NPP(simpar) = STAR_NPP(simpar) = AGN_NPP(simpar) = 0;
}
/*
void PadParticle2XYZDBL(SimParameters *simpar){ 
	PosType px,py,pz; 
	DEBUGPRINT0("###################################################\n"); 
	DEBUGPRINT0("You need to change the format of pad particles into XYZDBL format if needed\n"); 
	DEBUGPRINT0("###################################################\n"); 
	{ 
		treedmparticletype *bp = DM_TBPP(simpar); 
		ptrdiff_t i,np = DM_NPAD(simpar); 
		for(i=0;i<np;i++){ 
			px = XofP(simpar, bp); 
			py = YofP(simpar, bp); 
			pz = ZofP(simpar, bp); 
			bp ++; 
		} 
	}
}
*/


void DirectSummation(SimParameters *simpar, size_t ngrv, TPtlStruct *grvp, size_t npos, Position *pos){
	size_t i,j;
	double rspheresq = RSPHERE(simpar)*RSPHERE(simpar);
#ifdef GOTPM
	float vfact2 = Evol_FACT2(simpar);
#endif

	for(i=0;i<npos;i++){
		float ax,ay,az;
		ax = ay = az = 0;
		PosType px = pos[i].x;
		PosType py = pos[i].y;
		PosType pz = pos[i].z;
		for(j=0;j<ngrv;j++){
			PosType tmpx = px - grvp[j].x;
			PosType tmpy = py - grvp[j].y;
			PosType tmpz = pz - grvp[j].z;
			double dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
			if(dist2 <= rspheresq){
				float dist = sqrt(dist2);
				/*
				int ntmp = dist*invran2nran;
				PosType fplmf = forcecorrectdiff(ntmp,0) + forcecorrectslope(ntmp,0) * (dist-ntmp*ran2nran);
				*/

				float  fplmf = GetPointGrav(dist);
				fplmf = fplmf* grvp[j].mass;
				ax += tmpx*fplmf;
				ay += tmpy*fplmf;
				az += tmpz*fplmf;
			}
		}
#ifndef GOTPM
		pos->axyz[0] += ax;
		pos->axyz[1] += ay;
		pos->axyz[2] += az;
#else
		pos->axyz[0] += vfact2*ax;
		pos->axyz[1] += vfact2*ay;
		pos->axyz[2] += vfact2*az;
#endif
	}
}

int initflag=0;



void treecorrection(SimParameters *simpar, DeterminedEvolFact *EvolFactArr){
#ifdef GOTPM
	Evol_FACT2(simpar) = EvolFactArr[0].fact2;
#ifdef DEBUG
	DEBUGPRINT("P%d is now has vfact2 = %g\n",MYID(simpar), Evol_FACT2(simpar));
#endif
#endif
	PosType zmin, zmax, ymin, ymax, xmin, xmax;
	float theta,theta2;
	ptrdiff_t i,j,k;
	ptrdiff_t mx,my,mz;
	float rsphere, rspheresq;
	rsphere = RSPHERE(simpar);
	rspheresq = rsphere*rsphere;

	if(offsetof(treesphparticletype,ax) != offsetof(treestarparticletype,ax) ||
			offsetof(treesphparticletype, ax) != offsetof(treeagnparticletype, ax)){
		DEBUGPRINT("Error in the offset of ax %ld : %ld : %ld\n",
				offsetof(treesphparticletype,ax),offsetof(treestarparticletype,ax),
				offsetof(treeagnparticletype, ax));
		MPI_Finalize();
		exit(999);
	}


	theta = THETA(simpar);
	theta2 = theta*theta;



	/*
	xmin = COS_XMIN(simpar) - CellWidth;
	ymin = COS_YMIN(simpar) - CellWidth;
	zmin = COS_ZMIN(simpar) - CellWidth;
	xmax = COS_XMAX(simpar) + CellWidth;
	ymax = COS_YMAX(simpar) + CellWidth;
	zmax = COS_ZMAX(simpar) + CellWidth;
	*/
	xmin = SIM_LXMIN(simpar,dm) - CellWidth;
	ymin = SIM_LYMIN(simpar,dm) - CellWidth;
	zmin = SIM_LZMIN(simpar,dm) - CellWidth;
	xmax = SIM_LXMAX(simpar,dm) + CellWidth;
	ymax = SIM_LYMAX(simpar,dm) + CellWidth;
	zmax = SIM_LZMAX(simpar,dm) + CellWidth;


	mx = ceil((xmax-xmin)/CellWidth);
	my = ceil((ymax-ymin)/CellWidth);
	mz = ceil((zmax-zmin)/CellWidth);
	if(initflag==0){
		if(MYID(simpar)==0){
			printf("#################################################\n");
			printf("P%d Now initializing force array\n",MYID(simpar));
			printf("#################################################\n");
		}
		void i_force_spline(SimParameters *);
		i_force_spline(simpar);
		initflag = 1;
	}


#ifdef DEBUG
	{
		DoDeInfo *dd = TDM_DDINFO(simpar);
		ptrdiff_t nddinfo = NDDINFO(simpar);
		for(i=0;i<nddinfo;i++){
			DEBUGPRINT("P%d has %g %g : %g %g : %g %g\n",MYID(simpar),dd->lgroup.xyz.xmin, dd->lgroup.xyz.xmax,
					dd->lgroup.xyz.ymin, dd->lgroup.xyz.ymax,
					dd->lgroup.xyz.zmin, dd->lgroup.xyz.zmax);
			dd ++;
		}
	}
#endif
	{
		ptrdiff_t i;
		treedmparticletype *bp = DM_TBP(simpar);
		for(i=0;i<DM_NP(simpar);i++){
			CLEAR_FLAG(bp);
			SET_FLAG(bp,DMflag);
			bp++;
		}
	}


	void PaddingParticles(SimParameters *);
	PaddingParticles(simpar);


#ifdef DDEBUG
	{
		MPI_Barrier(MPI_COMM(simpar)); DEBUGPRINT("P%d %g %g : %g %g : %g %g ::: %ld %p\n",MYID(simpar),xmin,xmax,ymin,ymax,zmin,zmax, DM_NPAD(simpar), DM_TBPP(simpar));
		PosType xxmin,xxmax,yymin,yymax,zzmin,zzmax;
		xxmin = yymin = zzmin = 1.e20;
		xxmax = yymax = zzmax = -1.e20;
		treedmparticletype *bp;
		ptrdiff_t i;

		bp = DM_TBP(simpar);
		for(i=0;i<DM_NP(simpar);i++){
			CLEAR_FLAG(bp);
			SET_FLAG(bp, DMflag);
			if(!IS_FLAG_ON(bp,DMflag)){
				DEBUGPRINT("P%d found error in the particle definition: %d\n",MYID(simpar), (bp)->u4if.Flag[ENDIAN_OFFSET]);
				exit(999);
			}
			bp++;
		}

		bp = DM_TBPP(simpar);
		for(i=0;i<DM_NPAD(simpar);i++){
			xxmin = MIN(xxmin, bp->x);
			yymin = MIN(yymin, bp->y);
			zzmin = MIN(zzmin, bp->z);
			xxmax = MAX(xxmax, bp->x);
			yymax = MAX(yymax, bp->y);
			zzmax = MAX(zzmax, bp->z);
			if(!IS_FLAG_ON(bp,DMflag)){
				DEBUGPRINT("P%d found error in the particle definition: %d\n",MYID(simpar), IS_FLAG_ON(bp,DMflag));
				exit(999);
			}
			bp++;
		}
		MPI_Barrier(MPI_COMM(simpar)); DEBUGPRINT("P%d %g %g : %g %g : %g %g\n",MYID(simpar),xxmin,xxmax,yymin,yymax,zzmin,zzmax);
	}
#endif



	BASICCELL(simpar) = (TreeLinkedCell*)malloc(sizeof(TreeLinkedCell)*mx*my*mz);

	/*
	void HAMB(SimParameters *);
	HAMB(simpar);
	*/

	BuildLinkedList(simpar,mx,my,mz, xmin,ymin,zmin, zmax);

#ifdef DEBUG
	{
		for(i=0;i<mx*my*mz;i++){
			struct linkedlisttype *next = BASICCELL(simpar)[i].link;
			ptrdiff_t np = 0;
			while(next) {
				next= next->next;
				np ++;
			}
			if(np != BASICCELL(simpar)[i].nmem){
				fprintf(stderr,"P%d Error in counting %ld : %ld\n",MYID(simpar), np, BASICCELL(simpar)[i].nmem);
				exit(999);
			}
		}
		MPI_Barrier(MPI_COMM(simpar)); DEBUGPRINT("P%d here with mx/my/mz= %ld %ld %ld\n",MYID(simpar), mx,my,mz);
	}
#endif




#ifdef USE_GPU
	GPUDIRECTSUM = 2000 - 5000*(theta-0.2);
	GPUDIRECTSUM = min(GPUDIRECTSUM, 3000);
	GPUDIRECTSUM = max(GPUDIRECTSUM,  800);
#endif
	size_t mpos = 0;
	for(i=0;i<mx*my*mz;i++) mpos = MAX(mpos, BASICCELL(simpar)[i].nmem);

	/*
	*/

	Position *T_pos;
	TPtlStruct *T_grvp;
	TStruct *T_TREECELL;
	float *T_work;

#ifdef _OPENMP
#pragma omp parallel 
#endif
	{
		size_t linesize = (mpos/8+1)*8;
		if(omp_get_thread_num() ==0){
			T_pos = (Position*)malloc(sizeof(Position)*mpos*omp_get_num_threads());
			T_grvp = (TPtlStruct*)malloc(sizeof(TPtlStruct)*mpos*27*omp_get_num_threads());
			T_TREECELL = (TStruct*)malloc(sizeof(TStruct)*mpos*27*omp_get_num_threads());
			T_work = (float*)malloc(sizeof(float)*mpos*27*15*linesize*omp_get_num_threads());
		}
#ifdef _OPENMP
#pragma omp barrier
#endif
		Position *pos = T_pos + mpos * omp_get_thread_num();
		TPtlStruct *grvp = T_grvp + mpos* 27 * omp_get_thread_num();
		TStruct  *TREECELL = T_TREECELL + mpos* 27 * omp_get_thread_num();
		float  *work = T_work + mpos* 27 * 15 * linesize * omp_get_thread_num();
		ptrdiff_t ii;
#ifdef _OPENMP
#pragma omp for  schedule(dynamic)
#endif
		for(ii=0;ii<mx*my*mz;ii++){
			ptrdiff_t iz = ii/(mx*my); if(iz==mz-1 || iz ==0) continue;
			ptrdiff_t iy = (ii%(mx*my))/mx; if(iy==my-1 || iy ==0) continue;
			ptrdiff_t ix = ii%mx; if(ix==mx-1 || ix ==0) continue;
			size_t npos = Dump2Pos(simpar, ix,iy,iz,mx,my,mz, pos);
			if(npos>0){
				size_t ngrv = Dump2GrvP(simpar, ix,iy,iz,mx,my,mz,grvp);
				if(ngrv>0)
#ifndef USE_GPU
				{
					if(ngrv >=DIRECTSUM) {
						size_t linesize = (ngrv/8+1)*8;
						Make_Tree3(TREECELL,ngrv, grvp, ngrv, theta2,RECURSIVE);
						ptrdiff_t jj;
						for(jj=0;jj<npos;jj++){
							treeforce3(simpar, pos+jj, TREECELL, linesize, rspheresq, work);
						}
					}
					else {
						DirectSummation(simpar, ngrv, grvp, npos, pos);
					}
				}
#else
				{
					if(ngrv<=GPUOCPOUDIRECTSUM) DirectSummation(ngrv, grvp, npos, pos);
					else if(ngrv <= GPUDIRECTSUM) Direct_GPU_NBODY(pos, npos, grvp, ngrvp, GPUSPERNODE, 
							MYID(simpar), NX(simpar));
					else {
						ntreenodes = Make_Tree3(TREECELL, grv, ngrv, box, theta, flag_TW);
						gputreeforce();
					}
				}
#endif
	
			}
		}
	}
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		if(omp_get_thread_num()==0){
			free(T_pos);
			free(T_grvp);
			free(T_TREECELL);
		}
	}

#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar)); DEBUGPRINT("P%d is now after the tree correction \n",MYID(simpar));
#endif

	free(BASICCELL(simpar));
	free_padding_particles(simpar);
}
