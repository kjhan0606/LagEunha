#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>
#include<string.h>
#include "eunha.h"
//#include "mpirms.h"
#include "cosmology.h"
#include "mpiaux.h"




void PaddingParticles(SimParameters *simpar) { 
	SimBoxRange box; 
	box.x.min= box.y.min= box.z.min= 0; 
	box.x.max= NX(simpar); box.y.max= NY(simpar); box.z.max= NZ(simpar);
	PosType width = RSPHERE(simpar);
	NXNY(simpar) = NX(simpar)*NY(simpar); 
	if(MYID(simpar)==0) {
		DEBUGPRINT("P%d has width = %g\n",MYID(simpar), width);
	} 
	DM_TBPP(simpar) = NULL; SPH_TBPP(simpar) = NULL; VORO_TBPP(simpar) = NULL;  STAR_TBPP(simpar) = NULL;  AGN_TBPP(simpar) = NULL; 
	DM_NPAD(simpar) = SPH_NPAD(simpar) = VORO_NPAD(simpar) =  STAR_NPAD(simpar) = AGN_NPAD(simpar) = 0; 
	ppadding(DM_TBP(simpar), DM_NP(simpar), (void**)(&DM_TBPP(simpar)), 
			&DM_NPAD(simpar), TDM_DDINFO(simpar), NDDINFO(simpar), 
			box, width,&GRIDINFO(simpar), 3);
	ppadding(SPH_TBP(simpar), SPH_NP(simpar), (void**)(&SPH_TBPP(simpar)), 
			&SPH_NPAD(simpar), TSPH_DDINFO(simpar), NDDINFO(simpar), 
			box, width, &GRIDINFO(simpar), 3); 
	ppadding(VORO_TBP(simpar), VORO_NP(simpar), (void**)(&VORO_TBPP(simpar)), 
			&VORO_NPAD(simpar), TVORO_DDINFO(simpar), NDDINFO(simpar), 
			box, width, &GRIDINFO(simpar), 3); 
	ppadding(STAR_TBP(simpar), STAR_NP(simpar), (void**)(&STAR_TBPP(simpar)), &STAR_NPAD(simpar), TSTAR_DDINFO(simpar), NDDINFO(simpar), 
			box, width, &GRIDINFO(simpar), 3); 
	ppadding(AGN_TBP(simpar), AGN_NP(simpar), (void**)(&AGN_TBPP(simpar)), &AGN_NPAD(simpar), TAGN_DDINFO(simpar), NDDINFO(simpar), 
			box, width, &GRIDINFO(simpar), 3);
	ptrdiff_t i; 
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<DM_NPAD(simpar);i++) { 
		CLEAR_FLAG(DM_TBPP(simpar)+i);
		SET_FLAG(DM_TBPP(simpar)+i, BoundaryGhostflag); 
		SET_FLAG(DM_TBPP(simpar)+i, DMflag); 
	} 
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<SPH_NPAD(simpar);i++) { 
		CLEAR_FLAG(SPH_TBPP(simpar)+i); 
		SET_FLAG(SPH_TBPP(simpar)+i, BoundaryGhostflag); 
		SET_FLAG(SPH_TBPP(simpar)+i, SPHflag);
	} 
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NPAD(simpar);i++) { 
		CLEAR_FLAG(VORO_TBPP(simpar)+i); 
		SET_FLAG(VORO_TBPP(simpar)+i, BoundaryGhostflag); 
		SET_FLAG(VORO_TBPP(simpar)+i, VOROflag);
	} 
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<STAR_NPAD(simpar);i++) {
		CLEAR_FLAG(STAR_TBPP(simpar)+i); 
		SET_FLAG(STAR_TBPP(simpar)+i, STARflag); 
		SET_FLAG(STAR_TBPP(simpar)+i, BoundaryGhostflag); 
	} 
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<AGN_NPAD(simpar);i++){ 
		CLEAR_FLAG(AGN_TBPP(simpar)+i); 
		SET_FLAG(AGN_TBPP(simpar)+i, AGNflag); 
		SET_FLAG(AGN_TBPP(simpar)+i, BoundaryGhostflag);
	}
}



/* Declaration of prototypes of functions */
DefineProtoTypeFunctions(dm);
DefineProtoTypeFunctions(sph);
DefineProtoTypeFunctions(voro);
DefineProtoTypeFunctions(star);
DefineProtoTypeFunctions(agn);
DefineProtoTypeFunctions(treedm);
DefineProtoTypeFunctions(treesph);
DefineProtoTypeFunctions(treevoro);
DefineProtoTypeFunctions(treestar);
DefineProtoTypeFunctions(treeagn);



void BroadCastDDFuncAndDDInfo(SimParameters *simpar){
	CopyDDFuncDDInfoFromDM(simpar,TDM,treedm);
	CopyDDFuncDDInfoFromDM(simpar,SPH,sph);
	CopyDDFuncDDInfoFromDM(simpar,VORO,voro);
	CopyDDFuncDDInfoFromDM(simpar,STAR,star);
	CopyDDFuncDDInfoFromDM(simpar,AGN,agn);
	CopyDDFuncDDInfoFromDM(simpar,TSPH,treesph);
	CopyDDFuncDDInfoFromDM(simpar,TVORO,treevoro);
	CopyDDFuncDDInfoFromDM(simpar,TSTAR,treestar);
	CopyDDFuncDDInfoFromDM(simpar,TAGN,treeagn);
	/*
	int i,j;
	DoDeInfo *sph,*agn, *star, *dm;
	for(i=0;i<NDDINFO(simpar);i++){
		dm = (DM_DDINFO(simpar)) + i;
		star = (STAR_DDINFO(simpar)) + i;
		agn = (AGN_DDINFO(simpar)) + i;
		sph = (SPH_DDINFO(simpar)) + i;
		sph->npivot = agn->npivot = star->npivot = dm->npivot;
		sph->lgroup = agn->lgroup = star->lgroup = dm->lgroup;
		sph->pivot = (char*)malloc(sizeof(char)*sph->n_size*sph->npivot);
		star->pivot = (char*)malloc(sizeof(char)*sph->n_size*sph->npivot);
		agn->pivot = (char*)malloc(sizeof(char)*sph->n_size*sph->npivot);

		for(j=0;j<NPIVOT(simpar,i);j++){
			dmparticletype *dmp;
			sphparticletype *gasp;
			starparticletype *starp;
			agnparticletype *agnp;
			dmp = (dmparticletype*)(dm->pivot)+j;
			gasp = (sphparticletype*)(sph->pivot)+j;
			starp = (starparticletype*)(star->pivot)+j;
			agnp = (agnparticletype*)(agn->pivot)+j;
			gasp->x = starp->x = agnp->x = dmp->x; 
			gasp->y = starp->y = agnp->y = dmp->y; 
			gasp->z = starp->z = agnp->z = dmp->z; 
		}
	}
	*/
}





void BuildSimpleRMS(SimParameters *simpar, SimBoxRange box,  size_t n_size, DoDeFunc *ddfunc, DoDeInfo *ddinfo, MPI_Comm Comm){
	MPI_Status status;
	int src,dest;
	ptrdiff_t nrecv,nsend;
	ptrdiff_t i,j,k;
	int myid,nid;
	int igroup,subgroupsize,subgroupid;
	ptrdiff_t nowmem=0;


	SimBoxRange sbox = box;

	MPI_Comm_rank(Comm,&myid);
	MPI_Comm_size(Comm,&nid);

	int nsubgroup = getNextPrimeNumber(nid);
	if(nsubgroup == 1) {
		return;
	}
	NDDINFO(simpar) ++;
	int npivot = nsubgroup-1;
	dmparticletype pivot[npivot];

	int idir = getdirection(box);
	if(idir==0) ddfunc->xyzchip = 'x';
	else if(idir==1) ddfunc->xyzchip = 'y';
	else if(idir==2) ddfunc->xyzchip = 'z';
	else  ddfunc->xyzchip = 'w';

#ifdef DEBUG
	DEBUGPRINT("P%d divided along with %c\n",myid, ddfunc->xyzchip);
#endif
	subgroupsize = nid/nsubgroup;
	subgroupid = myid/subgroupsize;
	char xyzchip = ddfunc->xyzchip;


	/*
	void HAMB(SimParameters *);
	HAMB(simpar);
	*/

	for(i=0;i<npivot;i++){
		pivot[i].x = pivot[i].y = pivot[i].z = 0;

		if(idir ==0) {
			pivot[i].x = (box.x.max-box.x.min)/(npivot+1.L) *(i+1) + box.x.min;
  			DEBUGPRINT("P%d divided along %c with %g\n",myid, ddfunc->xyzchip, pivot[i].x);
		}
		else if(idir ==1)  {
			pivot[i].y = (box.y.max-box.y.min)/(npivot+1.L) *(i+1) + box.y.min;
			DEBUGPRINT("P%d divided along %c with %g\n",myid, ddfunc->xyzchip, pivot[i].y);
		}
		else if(idir ==2)  {
			pivot[i].z = (box.z.max-box.z.min)/(npivot+1.L) *(i+1) + box.z.min;
			DEBUGPRINT("P%d divided along %c with %g\n",myid, ddfunc->xyzchip, pivot[i].z);
		}
#ifdef XYZDBL
		CHANGEINDX( (pivot+i), 0L);
#endif
	}

	if(idir==0) {
		sbox.x.min = (box.x.max-box.x.min)/(npivot+1.L) * subgroupid + box.x.min;
		sbox.x.max = (box.x.max-box.x.min)/(npivot+1.L) * (subgroupid+1) + box.x.min;
	}
	else if(idir==1) {
		sbox.y.min = (box.y.max-box.y.min)/(npivot+1.L) * subgroupid + box.y.min;
		sbox.y.max = (box.y.max-box.y.min)/(npivot+1.L) * (subgroupid+1) + box.y.min;
	}
	else if(idir==2) {
		sbox.z.min = (box.z.max-box.z.min)/(npivot+1.L) * subgroupid + box.z.min;
		sbox.z.max = (box.z.max-box.z.min)/(npivot+1.L) * (subgroupid+1) + box.z.min;
	}
	{ 
		ddinfo->com = Comm;
		ddinfo->myid = myid;
		ddinfo->nid = nid;
		ddinfo->nsubgroup = nsubgroup;
		ddinfo->npivot = npivot;
		ddinfo->idirection = idir;
		ddinfo->subgroupsize = subgroupsize;
		ddinfo->subgroupid = subgroupid;
		ddinfo->xyzchip = ddfunc->xyzchip;
		ddinfo->xcompare = ddfunc->xcompare;
		ddinfo->ycompare = ddfunc->ycompare;
		ddinfo->zcompare = ddfunc->zcompare;
		ddinfo->wcompare = ddfunc->wcompare;

		ddinfo->xpinner = ddfunc->xpinner;
		ddinfo->ypinner = ddfunc->ypinner;
		ddinfo->zpinner = ddfunc->zpinner;
		ddinfo->muladd = ddfunc->muladd;
		ddinfo->wdist = ddfunc->wdist;

		ddinfo->insidebox = ddfunc->insidebox;
		ddinfo->edgeptl = ddfunc->edgeptl;
		ddinfo->n_size = n_size;
		memmove(ddinfo->pivot, pivot, npivot*n_size);
	}
	int (*compare)(GridInfo *, const void *, const void *);
	if(xyzchip=='x') compare = ddfunc->xcompare;
	else if(xyzchip=='y') compare = ddfunc->ycompare;
	else if(xyzchip=='z') compare = ddfunc->zcompare;
	else compare = ddfunc->wcompare;
	{
		MPI_Comm newcom;
		int key = myid % subgroupsize;
		MPI_Comm_split(Comm, subgroupid, key, &newcom);
		BuildSimpleRMS(simpar, sbox,  n_size, ddfunc, ddinfo+1,newcom);
		if(subgroupsize==1) MPI_Comm_free(&newcom);
	}
	return;
} 

void refine_kth_smallest(void *base, ptrdiff_t nmem, ptrdiff_t n_size, 
		DoDeFunc *ddfunc,
		int kth, int npivot, void *result, int myid, int nid,
		MPI_Comm com, GridInfo *gridinfo){
	ptrdiff_t ileft,iright,left,right,tleft;
	char swaptmp[n_size];
	char trypivotval[n_size];
	left = 0; right = nmem-1;
	ptrdiff_t tnmem;
	MPI_Reduce(&nmem, &tnmem, 1, MPI_INT64_T, MPI_SUM, 0,com);
	MPI_Bcast(&tnmem, 1, MPI_INT64_T,0,com);
	ptrdiff_t ntarget = (kth+1)*( (double)tnmem/(double)(npivot+1) );
	ptrdiff_t Nresolution = tnmem * DomainError;
	char xyzchip = ddfunc->xyzchip;
	int (*compare)(GridInfo *,  const void *, const void *);
	if(xyzchip == 'x') compare = (ddfunc->xcompare);
	else if(xyzchip == 'y') compare = (ddfunc->ycompare);
	else if(xyzchip == 'z') compare = (ddfunc->zcompare);
	else if(xyzchip == 'w') compare = (ddfunc->wcompare);
	int iter = 0;
	while(1){
		ileft = left;
		iright = right;
		memmove(trypivotval,result, n_size);
		/*
		MPI_Guess_kth_value(base,nmem,n_size, left,right, ddfunc, trypivotval,myid,nid,com);
		*/
		/*
#ifdef DEBUG
		if(myid==0) printf("Determined pivot val = %g ::: %ld %ld \n",(((AA*)trypivotval)->x), left,right);
#endif
*/
		do {
			while( ileft<nmem-1 && compare( gridinfo, (char*)base+ileft*n_size, trypivotval) <=0 ) ileft++;
			while( iright > 0 && compare(gridinfo,  trypivotval, (char*)base + iright*n_size)<=0 ) iright--;
			if(ileft<=iright){
				SWAP( (char*)base+ileft*n_size, (char*)base+iright*n_size, swaptmp);
				ileft++; iright--;
			}
		} while(ileft<iright);

		MPI_Reduce(&ileft, &tleft, 1, MPI_INT64_T, MPI_SUM,0, com);
		MPI_Bcast(&tleft, 1, MPI_INT64_T, 0, com);

		if(tleft > ntarget && llabs(tleft-ntarget) > Nresolution) right = iright;
		else if(tleft < ntarget && llabs(tleft-ntarget) > Nresolution) left = ileft;
		else break;
#ifdef DEBUG
		printf("P%d has tleft/ntarget %ld %ld  :: left/right %ld %ld\n",myid, tleft, ntarget, left,right);
#endif
		iter ++;
		if(iter > 32) break;
	}
	memmove(result, trypivotval, n_size);
}

void refine_pivotvals(void *base, ptrdiff_t nmem, ptrdiff_t n_size, DoDeFunc *ddfunc,
		int npivot, void *pivot,int myid, int nid, MPI_Comm com, GridInfo *gridinfo){
	int kth;
	for(kth=0;kth<npivot;kth++){
		char *pos =  (char*)pivot + kth*n_size;
		refine_kth_smallest(base, nmem, n_size, ddfunc, kth, npivot,pos, myid, nid, com, gridinfo);
	}
	return;
}



#ifdef DO_NOT_USE_IT_BECAUSE_IT_IS_UNDER_DEVELOPMENT

void refinempirks(void **ibase, ptrdiff_t *mmem, ptrdiff_t n_size, DoDeFunc *ddfunc, DoDeInfo *ddinfo,
		MPI_Comm Comm, GridInfo *gridinfo){
	MPI_Status status;
	int src,dest;
	ptrdiff_t nrecv,nsend;
	ptrdiff_t i,j,k;
	int myid,nid;
	ptrdiff_t nmem = *mmem;
	void *base = *ibase;
	char swaptmp[n_size];
	int igroup,subgroupsize,subgroupid;
	ptrdiff_t nowmem=0;
	void *rbase = NULL;
	char tlocalmin[n_size];
	char tlocalmax[n_size];

	MPI_Comm_rank(Comm,&myid);
	MPI_Comm_size(Comm,&nid);

	int nsubgroup = getNextPrimeNumber(nid);
	if(nsubgroup == 1) {
		return;
	}
	int npivot = nsubgroup-1;
	char aaa[n_size*npivot];
	void *pivot = (void *) aaa;

	/*
	int idir = ddfunc->divdir(base, nmem,Comm);
	if(idir==0) ddfunc->xyzchip = 'x';
	else if(idir==1) ddfunc->xyzchip = 'y';
	else if(idir==2) ddfunc->xyzchip = 'z';
	else  ddfunc->xyzchip = 'w';
	*/
	int idir;
	if(ddfunc->xyzchip == 'x') idir = 0;
	else if(ddfunc->xyzchip == 'y') idir = 1;
	else if(ddfunc->xyzchip == 'z') idir = 2;
	else idir = 3;


#ifdef DEBUG
	printf("P%d divided along with %c\n",myid, ddfunc->xyzchip);
#endif

	refine_pivotvals(base, nmem, n_size, ddfunc,npivot,pivot,myid, nid, Comm, gridinfo);

	subgroupsize = nid/nsubgroup;
	subgroupid = myid/subgroupsize;
	char xyzchip = ddfunc->xyzchip;

	{ 
		ddinfo->com = Comm;
		ddinfo->myid = myid;
		ddinfo->nid = nid;
		ddinfo->nsubgroup = nsubgroup;
		ddinfo->npivot = npivot;
		ddinfo->idirection = idir;
		ddinfo->subgroupsize = subgroupsize;
		ddinfo->subgroupid = subgroupid;
		ddinfo->xyzchip = ddfunc->xyzchip;
		ddinfo->xcompare = ddfunc->xcompare;
		ddinfo->ycompare = ddfunc->ycompare;
		ddinfo->zcompare = ddfunc->zcompare;
		ddinfo->wcompare = ddfunc->wcompare;

		ddinfo->xpinner = ddfunc->xpinner;
		ddinfo->ypinner = ddfunc->ypinner;
		ddinfo->zpinner = ddfunc->zpinner;
		ddinfo->muladd = ddfunc->muladd;
		ddinfo->wdist = ddfunc->wdist;

		ddinfo->insidebox = ddfunc->insidebox;
		ddinfo->edgeptl = ddfunc->edgeptl;
		ddinfo->n_size = n_size;
		memmove(ddinfo->pivot, pivot, npivot*n_size);
	}
	int (*compare)(GridInfo *, const void *, const void *);
	if(xyzchip=='x') compare = ddfunc->xcompare;
	else if(xyzchip=='y') compare = ddfunc->ycompare;
	else if(xyzchip=='z') compare = ddfunc->zcompare;
	else compare = ddfunc->wcompare;

	for(igroup=1;igroup<nsubgroup;igroup++){
		
		dest = (myid + subgroupsize*igroup + nid)%nid;
		src =  (myid - subgroupsize*igroup + nid)%nid;
		char localmin[n_size], localmax[n_size];
		int targetsubgid = (subgroupid + igroup+nsubgroup)%nsubgroup;

		char *left,*right;
		left = (char*)base;
		right = (char*)base + nmem*n_size;

		if(targetsubgid==0){
			memmove(tlocalmax, (char*)pivot+targetsubgid*n_size,n_size);
			for(;left<right;){
				if(compare(gridinfo, left, tlocalmax)<0){
					right -= n_size;
					SWAP(left,right, swaptmp);
				}
				else left+=n_size;
			}
		}
		else if(targetsubgid == nsubgroup-1){
			memmove(tlocalmin, (char*)pivot+(targetsubgid-1)*n_size,n_size);
			for(;left<right;){
				if(compare(gridinfo, left, tlocalmin) >= 0){
					right -= n_size;
					SWAP(left,right, swaptmp);
				}
				else left+=n_size;
			}
		}
		else {
			memmove(tlocalmin, (char*)pivot+(targetsubgid-1)*n_size,n_size);
			memmove(tlocalmax, (char*)pivot+(targetsubgid+0)*n_size,n_size);
			{
				for(;left<right;){
					if(compare(gridinfo, left,tlocalmin)>=0 && compare(gridinfo, left,tlocalmax)<0){
						right -= n_size;
						SWAP(left,right, swaptmp);
					}
					else left += n_size;
				}
			}
		}
		ptrdiff_t  nsend = nmem-((char*)right-(char *)base)/n_size;
		MPI_Sendrecv(&nsend,1,MPI_INT64_T,dest,0,&nrecv,1,MPI_INT64_T,src,0,Comm,&status);
		if(rbase == NULL) rbase = (void *)malloc(sizeof(char)*nrecv*n_size);
		else rbase = (void *)realloc(rbase,sizeof(char)*(nrecv + nowmem) *n_size);
		my_MPI_Sendrecv(right,nsend*sizeof(char)*n_size,MPI_BYTE,dest,0,
				(char*)rbase+nowmem*n_size,nrecv*sizeof(char)*n_size,MPI_BYTE,src,0,Comm,&status);
		nmem = nmem-nsend;
		base = (void*)realloc(base, nmem * n_size);
		nowmem += nrecv;
		if(myid==0) printf("+: %d : (%ld) <-> (%ld) %ld\n",igroup,nrecv, nsend,(ptrdiff_t)nowmem); fflush(stdout);
	}
    rbase = (void *)realloc(rbase, sizeof(char)*(nowmem+nmem)*n_size);
    memmove((char*)rbase + nowmem*n_size, base, nmem*n_size);
    nowmem += nmem;
	free(base);

	{
		MPI_Comm newcom;
		int key = myid % subgroupsize;
		MPI_Comm_split(Comm, subgroupid, key, &newcom);
		refinempirks(&rbase, &nowmem, n_size, ddfunc, ddinfo+1,newcom, gridinfo);
		if(subgroupsize==1) MPI_Comm_free(&newcom);
		/*
		MPI_Comm_free(&newcom);
		*/
	}
	*mmem = nowmem;
	*ibase =  rbase;
	return;
} 


#endif
