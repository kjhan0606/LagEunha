/* alpha version of parallel domain deocomposition of the Recursive MultiSection (RMS).
   It now allows only the one-dimensional data.
   Further implementation is needed to allow for arbitrary dimensional data. */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<limits.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "mpirms.h"
#include "mpiaux.h"
#define MAX(a,b) ((a)>(b)?(a):(b))

#ifdef TEST
#include "test.h"
#else
typedef dmparticletype AA;
#endif

#ifdef XYZDBL
static int xyzdblflag = 1;
#else
static int xyzdblflag = 0;
#endif

int getNextPrimeNumber(int number){
	int i;
	int result = 1;
	for(i=2;i<=number;){
		while((number % i) ==0){
			result = MAX(i,result);
			number = number/i;
		}
		if(i==2) i = 3;
		else i += 2;
	}
	return result;
}

int getprimenumber(int number,PrimeNumber *prime){
	int i,j,k,denom,nprime;
	int res;


	nprime = 0;
	prime[nprime].factor = 0;
	prime[nprime].prime = 0;

	do{
		int pr=getNextPrimeNumber(number);
		if(prime[nprime].factor ==0){
			prime[nprime].prime = pr;
			prime[nprime].factor = 1;
		}
		else if(pr==prime[nprime].prime) {
			prime[nprime].factor ++;
		}
		else {
			nprime ++;
			prime[nprime].prime = pr;
			prime[nprime].factor = 1;
		}
		number = number/pr;
	} while (number>1);
	nprime ++;

	return nprime;
}

void ExtractLocalDomainVolume(DoDeInfo *ddinfo, int nddinfo, SimBoxRange box){
	int i;
	ddinfo->lgroup.xyz.xmin = box.x.min; 
	ddinfo->lgroup.xyz.ymin = box.y.min; 
	ddinfo->lgroup.xyz.zmin = box.z.min; 
	ddinfo->lgroup.xyz.wmin = box.w.min;

	ddinfo->lgroup.xyz.xmax = box.x.max; 
	ddinfo->lgroup.xyz.ymax = box.y.max; 
	ddinfo->lgroup.xyz.zmax = box.z.max; 
	ddinfo->lgroup.xyz.wmax = box.w.max;
	for(i=0;i<nddinfo;i++){
		if(i!=0) ddinfo[i].lgroup.r = ddinfo[i-1].lgroup.r;
		int xyzoffset;
		if(ddinfo[i].xyzchip == 'x') xyzoffset = ddinfo[i].memoffset.xyz.x;
		else if(ddinfo[i].xyzchip == 'y') xyzoffset = ddinfo[i].memoffset.xyz.y;
		else if(ddinfo[i].xyzchip == 'z') xyzoffset = ddinfo[i].memoffset.xyz.z;
		else xyzoffset = ddinfo[i].memoffset.xyz.w;
		if(ddinfo[i].subgroupid==0) {
			ddinfo[i].lgroup.r.rmax[ddinfo[i].idirection] = *(float*)((char*)
					(&ddinfo[i].pivot[ddinfo[i].subgroupid*ddinfo[i].n_size])+xyzoffset);
#ifdef DEBUG
			DEBUGPRINT("Now p%d rmax %g %g %g\n",i,
					*(float*)(ddinfo[i].pivot+ddinfo[i].subgroupid*ddinfo[i].n_size+ddinfo[i].memoffset.xyz.x),
					*(float*)(ddinfo[i].pivot+ddinfo[i].subgroupid*ddinfo[i].n_size+ddinfo[i].memoffset.xyz.y),
					*(float*)(ddinfo[i].pivot+ddinfo[i].subgroupid*ddinfo[i].n_size+ddinfo[i].memoffset.xyz.z)
					);
			DEBUGPRINT("+Now p%d rmax %g %g %g :: %ld : %ld %ld %ld\n",i,
					( (AA*)(ddinfo[i].pivot+ddinfo[i].subgroupid * ddinfo[i].n_size))->x,
					( (AA*)(ddinfo[i].pivot+ddinfo[i].subgroupid * ddinfo[i].n_size))->y,
					( (AA*)(ddinfo[i].pivot+ddinfo[i].subgroupid * ddinfo[i].n_size))->z,
					ddinfo[i].n_size, ddinfo[i].memoffset.xyz.x, ddinfo[i].memoffset.xyz.y, ddinfo[i].memoffset.xyz.z
					);
#endif
		}
		else if(ddinfo[i].subgroupid != ddinfo[i].nsubgroup-1){
			ddinfo[i].lgroup.r.rmax[ddinfo[i].idirection] = *(float*)((char*)
					(&ddinfo[i].pivot[ddinfo[i].subgroupid*ddinfo[i].n_size])+xyzoffset);
			ddinfo[i].lgroup.r.rmin[ddinfo[i].idirection] = *(float*)((char*)
					(&ddinfo[i].pivot[(ddinfo[i].subgroupid-1)*ddinfo[i].n_size])+xyzoffset);
		}
		else {
			ddinfo[i].lgroup.r.rmin[ddinfo[i].idirection] = *(float*)((char*)
					(&ddinfo[i].pivot[(ddinfo[i].subgroupid-1)*ddinfo[i].n_size])+xyzoffset);
#ifdef DEBUG
			DEBUGPRINT("Now p%d rmin %g %g %g\n",i,
					*(float*)((char*)(&ddinfo[i].pivot[(ddinfo[i].subgroupid-1)*ddinfo[i].n_size])+ddinfo[i].memoffset.xyz.x),
					*(float*)((char*)(&ddinfo[i].pivot[(ddinfo[i].subgroupid-1)*ddinfo[i].n_size])+ddinfo[i].memoffset.xyz.y),
					*(float*)((char*)(&ddinfo[i].pivot[(ddinfo[i].subgroupid-1)*ddinfo[i].n_size])+ddinfo[i].memoffset.xyz.z)
				  );
#endif
		}
	}
}

void Old_MPI_Guess_kth_value(void *a, ptrdiff_t nmem, ptrdiff_t n_size,  ptrdiff_t left, ptrdiff_t right,
		DoDeFunc *ddfunc,
		void *result, int myid, int nid, MPI_Comm com, GridInfo *gridinfo){
	int isum=0;
	int (*muladd)(GridInfo *, const void *, const void *, float, char, int);
	muladd = ddfunc->muladd;
	char xyzchip = ddfunc->xyzchip;
	char guess[nid*n_size];
	char myguess[n_size];
	char *aa = (char*)result;
	int i;
	ptrdiff_t trypivot = (left+right)/2;
	for(i=0;i<n_size;i++) aa[i] = '\0'; /* Initialization is needed */
	if(trypivot <0 || left >= nmem-1 || right <1) { /* if it is negative member or pivot is beyond the range 
												   of the local domain*/
		float aa,fact = 0./0.;
		muladd(gridinfo, myguess, myguess, fact,xyzchip, 1);
	}
	else {
		memmove(myguess, (char*)a+trypivot*n_size,n_size); 
	}
	MPI_Gather(myguess, n_size, MPI_BYTE, guess, n_size, MPI_BYTE,0, com);

	if(myid==0) {
		int i;
		for(i=0;i<nid;i++) if(muladd(gridinfo,  result, guess+i*n_size, 1,xyzchip, 1)) isum ++;
		float norm = 1./isum;
		norm -= 1.;
		muladd(gridinfo, result, result, norm, xyzchip, -1);
#ifdef DEBUG
		DEBUGPRINT("P%d has pivot value %g %g %g trypivot %ld ::: %c ::: %g %d\n",myid, (((AA*)result))->x, (((AA*)result))->y,
			(((AA*)result))->z,trypivot, xyzchip, norm, isum);
#endif
	}
	MPI_Bcast(result, n_size, MPI_BYTE, 0, com);
}
void Old_Find_kth_smallest(void *base, ptrdiff_t nmem, ptrdiff_t n_size, DoDeFunc *ddfunc,
		int kth, int npivot, void *result, int myid, int nid, MPI_Comm com, GridInfo *gridinfo){
	ptrdiff_t ileft,iright,left,right,tleft, tright;
	char swaptmp[n_size];
	char trypivotval[n_size];
	left = 0; right = nmem-1;
	ptrdiff_t tnmem;
	MPI_Reduce(&nmem, &tnmem, 1, MPI_INT64_T, MPI_SUM, 0,com);
	MPI_Bcast(&tnmem, 1, MPI_INT64_T,0,com);
	ptrdiff_t ntarget = (kth+1)*( (double)tnmem/(double)(npivot+1) );
	ptrdiff_t Nresolution = tnmem * DomainError/nid;
	char xyzchip = ddfunc->xyzchip;
	int (*compare)(GridInfo *, const void *, const void *);
	if(xyzchip == 'x') compare = (ddfunc->xcompare);
	else if(xyzchip == 'y') compare = (ddfunc->ycompare);
	else if(xyzchip == 'z') compare = (ddfunc->zcompare);
	else if(xyzchip == 'w') compare = (ddfunc->wcompare);
	int iter = 0;
	while(1){
		ileft = left;
		iright = right;
		Old_MPI_Guess_kth_value(base,nmem,n_size, left,right, ddfunc, trypivotval,myid,nid,com,gridinfo);
#ifdef DEBUG
		if(myid==0) DEBUGPRINT("Determined pivot val = %g %g %g::: %ld %ld \n",(((AA*)trypivotval)->x),
				(((AA*)trypivotval)->y),(((AA*)trypivotval)->z), left,right);
#endif
		do {
			while( ileft<nmem-1 && compare(gridinfo, (char*)base+ileft*n_size, trypivotval) <0 ) ileft++;
			while( iright>0 &&  compare(gridinfo, trypivotval, (char*)base + iright*n_size)<=0 ) iright--;
			if(ileft<iright){
				SWAP( (char*)base+ileft*n_size, (char*)base+iright*n_size, swaptmp);
				ileft++; iright--;
			}
		} while(ileft<iright);

		MPI_Reduce(&ileft, &tleft, 1, MPI_INT64_T, MPI_SUM,0, com);
		MPI_Bcast(&tleft, 1, MPI_INT64_T, 0, com);

		if(tleft > ntarget && llabs(tleft-ntarget) > Nresolution) right = iright;
		else if(tleft < ntarget && llabs(tleft-ntarget) > Nresolution) left = ileft;
		else break;
		iter ++;
		if(iter >=MAX_ITER) break;
	}
#ifdef DEBUG
	{
		float aa = 100*(float) tleft/(float) ntarget;
		int nyid;
		MPI_Comm_rank(MPI_COMM_WORLD,&nyid);
		DEBUGPRINT("P%d has : %g percent ::: tleft/ntarget %ld %ld  :: left/right %ld %ld from %ld in %d with nid= %d found at iter= %d\n",
				nyid, aa,
				tleft, ntarget, left,right, tnmem, kth, nid,iter);
	}
#endif
	memmove(result, trypivotval, n_size);
}

typedef struct WithinRanges{
	ptrdiff_t nleft;
	char pivot[2048];
} WithinRanges;

void MPI_Guess_kth_value(void *a, ptrdiff_t nmem, ptrdiff_t n_size,  
		int nid, MPI_Comm com, GridInfo *gridinfo, int (*compare)(GridInfo *, const void *, const void *), 
		ptrdiff_t ntarget, WithinRanges *lowerupper){
	char guess[nid*n_size];
	char myguess[n_size];
	int i;
	if(n_size > 2048){
		DEBUGPRINT0("Error: increase the max_n_size\n");
		MPI_Finalize();exit(999);
	}
	lowerupper[0].nleft = 0;
	lowerupper[1].nleft = LONG_MAX;
	ptrdiff_t trypivot = nmem/2;
	do{
		memmove(myguess, (char*)a+trypivot*n_size,n_size); 
		MPI_Gather(myguess, n_size, MPI_BYTE, guess, n_size, MPI_BYTE,0, com);
		MPI_Bcast(guess, nid*n_size, MPI_BYTE, 0, com);
		for(i=0;i<nid;i++){
			ptrdiff_t j,nleft = 0, tnleft;
	
			for(j=0;j<nmem;j++) if(compare(gridinfo, (char*)a+j*n_size, (char*)guess + i*n_size) < 0) nleft ++;
			MPI_Reduce(&nleft, &tnleft, 1, MPI_INT64_T, MPI_SUM, 0, com);
			MPI_Bcast(&tnleft, 1, MPI_INT64_T,  0, com);
			if(tnleft > lowerupper[0].nleft && tnleft < ntarget){
				lowerupper[0].nleft = tnleft;
				memcpy(lowerupper[0].pivot, (char*) guess + i*n_size, n_size);
			}
			if(tnleft < lowerupper[1].nleft && tnleft >= ntarget){
				lowerupper[1].nleft = tnleft;
				memcpy(lowerupper[1].pivot, (char*) guess + i*n_size, n_size);
			}
	
		}
		trypivot = (trypivot*7+ 17 + nmem)%nmem;
	} while(lowerupper[0].nleft == 0 || lowerupper[1].nleft == LONG_MAX);
/*
#ifdef DEBUG
	DEBUGPRINT("ntarget = %ld ::; between %ld and %ld\n",ntarget, lowerupper[0].nleft,
			lowerupper[1].nleft);
#endif
*/
	return;
}
ptrdiff_t CoarseSort(void *base, ptrdiff_t nmem, ptrdiff_t n_size,  
		int nid, MPI_Comm com, GridInfo *gridinfo, int (*compare)(GridInfo *, const void *, const void *), 
		void *pivot){
	ptrdiff_t ilower, iupper;
	ilower = 0; iupper = nmem - 1;
	do {
		char swaptmp[n_size];
		while( ilower<nmem-1 && compare(gridinfo, (char*)base+ilower*n_size, pivot) <0 ) ilower++;
		while( iupper>0 &&  compare(gridinfo, pivot, (char*)base + iupper*n_size)<=0 ) iupper--;
		if(ilower<iupper){
			SWAP( (char*)base+ilower*n_size, (char*)base+iupper*n_size, swaptmp);
			ilower++; iupper--;
		}
	} while(ilower<iupper);
	return ilower;
}
void GetPivotVal(void *base, ptrdiff_t iupper, ptrdiff_t ilower, int n_size, int myid, int nid, MPI_Comm com, 
		void *trypivotval){

	ptrdiff_t iwidth = (iupper - ilower);
	ptrdiff_t tiwidth[nid];
	ptrdiff_t trypivot;
	int id=-1;
	MPI_Gather(&iwidth, 1, MPI_INT64_T, tiwidth, 1, MPI_INT64_T, 0, com);
	if(myid==0){
		int i;
		ptrdiff_t imax = 0;
		for(i=0;i<nid;i++){
			if(tiwidth[i] > imax) {
				imax = tiwidth[i]; id = i;
			}
		}
	}
	MPI_Bcast(&id, 1, MPI_INT, 0, com);
	if(myid == id){
		trypivot = (iupper + ilower) /2;
		memcpy(trypivotval, (char*)base + trypivot * n_size, n_size);
#ifdef DEBUG
		DEBUGPRINT("P%d : Determined pivot val = %g %g %g::: %ld %ld \n",id,(((AA*)trypivotval)->x),
				(((AA*)trypivotval)->y),(((AA*)trypivotval)->z), ilower,iupper);
#endif
	}
	MPI_Bcast(trypivotval, n_size, MPI_BYTE, id, com);
}

void Find_kth_smallest(void *base, ptrdiff_t nmem, ptrdiff_t n_size, 
		DoDeFunc *ddfunc,
		int kth, int npivot, void *result, int myid, int nid,
		MPI_Comm com, GridInfo *gridinfo){
	ptrdiff_t ilower,iupper,left,right,tleft, tright, ipivot;
	char swaptmp[n_size];
	char trypivotval[n_size];
	left = 0; right = nmem-1;
	ptrdiff_t tnmem;
	MPI_Reduce(&nmem, &tnmem, 1, MPI_INT64_T, MPI_SUM, 0,com);
	MPI_Bcast(&tnmem, 1, MPI_INT64_T,0,com);
	ptrdiff_t ntarget = (kth+1)*( (double)tnmem/(double)(npivot+1) );
	ptrdiff_t Nresolution = tnmem * DomainError/nid;
	char xyzchip = ddfunc->xyzchip;
	int (*compare)(GridInfo *, const void *, const void *);
	if(xyzchip == 'x') compare = (ddfunc->xcompare);
	else if(xyzchip == 'y') compare = (ddfunc->ycompare);
	else if(xyzchip == 'z') compare = (ddfunc->zcompare);
	else if(xyzchip == 'w') compare = (ddfunc->wcompare);
	WithinRanges lowerupper[2];
	MPI_Guess_kth_value(base,nmem,n_size, nid, com, gridinfo, compare, ntarget, lowerupper);
	ilower = CoarseSort(base,nmem,n_size, nid, com, gridinfo, compare, lowerupper[0].pivot);
	iupper = CoarseSort(base,nmem,n_size, nid, com, gridinfo, compare, lowerupper[1].pivot);
	int iter = 0;
	do{
		GetPivotVal(base, iupper, ilower, n_size, myid,nid, com, trypivotval);
		ipivot = CoarseSort(base, nmem, n_size, nid, com, gridinfo, compare, trypivotval);

		MPI_Reduce(&ipivot, &tleft, 1, MPI_INT64_T, MPI_SUM,0, com);
		MPI_Bcast(&tleft, 1, MPI_INT64_T, 0, com);

		if(tleft > ntarget && llabs(tleft-ntarget) > Nresolution) {
			iupper = ipivot;
			memcpy(lowerupper[1].pivot, trypivotval, n_size);
			lowerupper[1].nleft = tleft;
		}
		else if(tleft < ntarget && llabs(tleft-ntarget) > Nresolution) {
			ilower = ipivot;
			memcpy(lowerupper[0].pivot, trypivotval, n_size);
			lowerupper[0].nleft = tleft;
		}
		else break;
		iter ++;
		if(iter >=MAX_ITER) break;
#ifdef DEBUG
		if(myid==0) DEBUGPRINT(" lower / upper = %ld / %ld :: Nresolution = %ld\n",
				lowerupper[0].nleft, lowerupper[1].nleft, Nresolution);
#endif
	} while( (lowerupper[1].nleft - lowerupper[0].nleft) > Nresolution);
#ifdef DEBUG
	{
		float aa = 100*(float) tleft/(float) ntarget;
		int nyid;
		MPI_Comm_rank(MPI_COMM_WORLD,&nyid);
		DEBUGPRINT("P%d has : %g percent ::: tleft/ntarget %ld %ld  :: left/right %ld %ld from %ld in %d with nid= %d found at iter= %d\n",
				nyid, aa,
				tleft, ntarget, left,right, tnmem, kth, nid,iter);
	}
#endif
	memmove(result, trypivotval, n_size);
}

void Find_pivotvals(void *base, ptrdiff_t nmem, ptrdiff_t n_size, DoDeFunc *ddfunc,
		int npivot, void *pivot,int myid, int nid, MPI_Comm com, GridInfo *gridinfo){
	int kth;
	for(kth=0;kth<npivot;kth++){
		char *pos =  (char*)pivot + kth*n_size;
		Find_kth_smallest(base, nmem, n_size, ddfunc, kth, npivot,pos, myid, nid, com, gridinfo);
	}
	return;
}



void mpirms(void **ibase, ptrdiff_t *mmem, ptrdiff_t n_size, DoDeFunc *ddfunc, DoDeInfo *ddinfo,
		MPI_Comm Comm, GridInfo *gridinfo , enum DDdirection dddirection, enum RMSComBuild newrmsbuild){
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

	int idir ;

	if(dddirection == newway){
		idir = ddfunc->divdir(gridinfo, base, nmem,Comm);
		if(idir==0) ddfunc->xyzchip = 'x';
		else if(idir==1) ddfunc->xyzchip = 'y';
		else if(idir==2) ddfunc->xyzchip = 'z';
		else  ddfunc->xyzchip = 'w';
	}
	else if(dddirection == oldway){
		ddfunc->xyzchip = ddinfo->xyzchip;
		if(ddinfo->xyzchip == 'x') idir = 0;
		else if(ddinfo->xyzchip == 'y') idir = 1;
		else if(ddinfo->xyzchip == 'z') idir = 2;
		else if(ddinfo->xyzchip == 'w') idir = 3;
	}


#ifdef DEBUG
	if(myid==0) DEBUGPRINT("P%d divided along with %c\n",myid, ddfunc->xyzchip);
#endif

	Find_pivotvals(base, nmem, n_size, ddfunc,npivot,pivot,myid, nid, Comm, gridinfo);


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
#ifdef DEBUG
		for(i=0;i<npivot;i++){
			DEBUGPRINT0("##############################\n");
			DEBUGPRINT("P%d has %ld pivot results    %g %g %g in direction of %c\n",
					myid, i, ((AA*)(ddinfo->pivot) + i)->x, ((AA*)(ddinfo->pivot) + i)->y, 
					((AA*)(ddinfo->pivot) + i)->z, ddfunc->xyzchip);
			DEBUGPRINT0("##############################\n");
		}
#endif
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

		ptrdiff_t  nsend=0;
		char *left,*right;
		if(nmem > 0){
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
			nsend = nmem-((char*)right-(char *)base)/n_size;
		}

		MPI_Sendrecv(&nsend,1,MPI_INT64_T,dest,0,&nrecv,1,MPI_INT64_T,src,0,Comm,&status);
		if(rbase == NULL) rbase = (void *)malloc(sizeof(char)*nrecv*n_size);
		else rbase = (void *)realloc(rbase,sizeof(char)*(nrecv + nowmem) *n_size);

		my_MPI_Sendrecv(right,nsend*sizeof(char)*n_size,MPI_BYTE,dest,0,
				(char*)rbase+nowmem*n_size,nrecv*sizeof(char)*n_size,MPI_BYTE,src,0,Comm,&status);
		nmem = nmem-nsend;
		base = (void*)realloc(base, nmem * n_size);


		nowmem += nrecv;
#ifdef DEBUG
		DEBUGPRINT("P%d +: %d : s(P%d;%ld) <-> r(P%d;%ld) %ld :::: %ld particles afterward.\n", myid,igroup,
				dest,nsend, src,
				nrecv, nowmem, nmem); fflush(stdout);
#endif
	}
    base = (void *)realloc(base, sizeof(char)*(nowmem+nmem)*n_size);
    memcpy((char*)base + nmem*n_size, rbase, nowmem*n_size);
    nowmem += nmem;
	free(rbase);

	{
		MPI_Comm newcom;
		if(newrmsbuild == OldCom && subgroupsize != 1){
			newcom = (ddinfo+1)->com;
		}
		else {
			int key = myid % subgroupsize;
			MPI_Comm_split(Comm, subgroupid, key, &newcom);
		}
		mpirms(&base, &nowmem, n_size, ddfunc, ddinfo+1,newcom, gridinfo, dddirection, newrmsbuild);
		if(subgroupsize ==1) MPI_Comm_free(&newcom);
	}
	*mmem = nowmem;
	*ibase =  base;
	/*
    rbase = (void *)realloc(rbase, sizeof(char)*(nowmem+nmem)*n_size);
    memmove((char*)rbase + nowmem*n_size, base, nmem*n_size);
    nowmem += nmem;
	free(base);

	{
		MPI_Comm newcom;
		if(newrmsbuild == OldCom && subgroupsize != 1){
			newcom = (ddinfo+1)->com;
		}
		else {
			int key = myid % subgroupsize;
			MPI_Comm_split(Comm, subgroupid, key, &newcom);
		}
		mpirms(&rbase, &nowmem, n_size, ddfunc, ddinfo+1,newcom, gridinfo, dddirection, newrmsbuild);
		if(subgroupsize ==1) MPI_Comm_free(&newcom);
	}
	*mmem = nowmem;
	*ibase =  rbase;
	*/
	return;
} 


void pmigrate(void **ibase,ptrdiff_t *mmem, DoDeInfo *ddinfo, GridInfo *gridinfo){
	MPI_Status status;
	MPI_Comm com;
	int src,dest;
	ptrdiff_t nrecv,nsend;
	ptrdiff_t i,j,k;
	int myid,nid;
	ptrdiff_t nmem = *mmem;
	void *base = *ibase;
	ptrdiff_t n_size = ddinfo->n_size;
	char swaptmp[n_size];
	int igroup,subgroupsize,subgroupid;
	ptrdiff_t nowmem=0;
	void *rbase = NULL;
	char tlocalmin[n_size];
	char tlocalmax[n_size];




	myid = ddinfo->myid;
	nid = ddinfo->nid;
	com = ddinfo->com;

	int nsubgroup = ddinfo->nsubgroup;

	/*
	 * */

	int npivot = ddinfo->npivot;
	void *pivot = ddinfo->pivot;


	int idir = ddinfo->idirection;

	subgroupsize = ddinfo->subgroupsize;
	subgroupid = ddinfo->subgroupid;
	char xyzchip = ddinfo->xyzchip;
	int (*compare)(GridInfo *, const void *, const void *);
	if(xyzchip == 'x') compare = ddinfo->xcompare;
	else if(xyzchip == 'y') compare = ddinfo->ycompare;
	else if(xyzchip == 'z') compare = ddinfo->zcompare;
	else compare = ddinfo->wcompare;

	for(igroup=1;igroup<nsubgroup;igroup++){
		dest = (myid + subgroupsize*igroup + nid)%nid;
		src =  (myid - subgroupsize*igroup + nid)%nid;
		char localmin[n_size], localmax[n_size];
		int targetsubgid = (subgroupid + igroup+nsubgroup)%nsubgroup;

		char *left,*right;
		ptrdiff_t nsend;
		if(nmem >0) {
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
			nsend = nmem-((char*)right-(char *)base)/n_size;
		}
		else{
			nsend = 0;
		}

		MPI_Sendrecv(&nsend,1,MPI_INT64_T,dest,0,&nrecv,1,MPI_INT64_T,src,0,com,&status);
		if(rbase == NULL) rbase = (void *)malloc(sizeof(char)*nrecv*n_size);
		else rbase = (void *)realloc(rbase,sizeof(char)*(nrecv + nowmem) *n_size);
		my_MPI_Sendrecv(right,nsend*sizeof(char)*n_size,MPI_BYTE,dest,0,
				(char*)rbase+nowmem*n_size,nrecv*sizeof(char)*n_size,MPI_BYTE,src,0,com,&status);
		nmem = nmem-nsend;
		base = (void*)realloc(base, nmem * n_size);
		nowmem += nrecv;
#ifdef DEBUG
		if(myid==0) DEBUGPRINT("+: %d : (%ld) <-> (%ld) %ld\n",igroup,nrecv, nsend,nowmem); fflush(stdout);
#endif
	}
    base = (void *)realloc(base, sizeof(char)*(nowmem+nmem)*n_size);
    memcpy((char*)base + nmem*n_size, rbase, nowmem*n_size);
    nowmem += nmem;
	free(rbase);
	if(subgroupsize > 1) pmigrate(&base, &nowmem, ddinfo+1,gridinfo);
	*mmem = nowmem;
	*ibase =  base;
	/*
    rbase = (void *)realloc(rbase, sizeof(char)*(nowmem+nmem)*n_size);
    memmove((char*)rbase + nowmem*n_size, base, nmem*n_size);
    nowmem += nmem;
	free(base);
	if(subgroupsize > 1) pmigrate(&rbase, &nowmem, ddinfo+1,gridinfo);
	*mmem = nowmem;
	*ibase =  rbase;
	*/
	return;
}

void *ExtractPtls2Send(void *base, ptrdiff_t nmem, void *padbase, ptrdiff_t npad, DoDeInfo *ddinfo, 
		BoxMinMax targetbox, PosType width, SimBoxRange simbox,
		ptrdiff_t *sizeinbyte2send, ptrdiff_t npad_in_current_level, GridInfo *gridinfo){
	char *run;
	ptrdiff_t i;
	ptrdiff_t isend = 0;
	int n_size = ddinfo->n_size;
	int (*insidebox)(GridInfo *, const void *, BoxMinMax *, PosType *,  SimBoxRange *, int , int) = ddinfo->insidebox;
	int (*edgeptl)(GridInfo *, const void *, SimBoxRange *, PosType *) = ddinfo->edgeptl;

	ptrdiff_t maxsize = MAX(nmem*0.1+npad_in_current_level,10000);
	ptrdiff_t stepsendsize = maxsize*0.2;
	void *sendbuff = (void *)malloc(maxsize*n_size);
	run = (char *)padbase;
	isend = 0;
	/* second last argument value of insidebox(gridinfo, run, &targetbox, &width,&simbox,0, 13) 
	 * has multiple meanging... If it is 0, then it deals with a pad particle, xPos does not applied and
	 * the second argument should be modified .*/
	for(i=0;i<npad_in_current_level;i++){
		if(insidebox(gridinfo, run, &targetbox, &width,&simbox,0, 13)) {
			memmove((char*)sendbuff+isend*n_size, run, n_size);
			isend ++;
		}
		run += n_size;
	}

	run = (char*)base;
	for(i=0;i<nmem;i++){
		int pmin,pmax;
		if(edgeptl(gridinfo, run,&simbox,&width)) {
			pmin = 0;
			pmax = 26;
		}
		else {
			pmin = pmax = 13;
		}
		char buff[n_size];
		memmove(buff, run, n_size);
		int pflag;
		for(pflag=pmin;pflag<=pmax;pflag++) {
			if(insidebox(gridinfo, buff, &targetbox, &width,&simbox,1, pflag)) 
			{
				if(isend >=maxsize) {
					maxsize += stepsendsize;
					sendbuff = (void*)realloc(sendbuff, maxsize*n_size);
				}
				memmove((char*)sendbuff+isend*n_size, buff, n_size);
				isend ++;
				memmove(buff, run, n_size);
			}
		}
		run += n_size;
	}
	sendbuff = (void*)realloc(sendbuff, isend*n_size);
	*sizeinbyte2send = isend*n_size;
	return sendbuff;
}

void DeleteOutofLocalBox(void **padbase, ptrdiff_t *npad, ptrdiff_t n_size, PosType width, BoxMinMax *lbox, 
		SimBoxRange *simbox,
		int (*insidebox)(GridInfo *, const void *, BoxMinMax *, PosType *, SimBoxRange *,int , int),
			GridInfo *gridinfo){
	ptrdiff_t i;
	char *left, *right;
	left = right = (char*)(*padbase);
	for(i=0;i<*npad;i++){
		if(insidebox(gridinfo, right, lbox, &width, simbox, 0, 13)) {
			memmove(left,right,n_size);
			left += n_size;
		}
		right += n_size;
	}
	*npad = (left  - (char*)(*padbase))/n_size;
	*padbase = (void *) realloc( *padbase, *npad * n_size);
}

void justppadding(void *base, ptrdiff_t nmem, void **padbase, ptrdiff_t *npad, 
		DoDeInfo *ddinfo, SimBoxRange box, PosType width, GridInfo *gridinfo){
	MPI_Status status;
	MPI_Comm com;
	int src, dest;
	ptrdiff_t nrec, nsend;
	ptrdiff_t i,j,k;
	int myid, nid;
	ptrdiff_t n_size = ddinfo->n_size;
	int igroup, subgroupsize, subgroupid;
	BoxMinMax mylocalbox,yourlocalbox;
	mylocalbox= ddinfo->lgroup.xyz;
	myid = ddinfo->myid;
	nid = ddinfo->nid;
	com = ddinfo->com;
	int nsubgroup = ddinfo->nsubgroup;
	subgroupsize = ddinfo->subgroupsize;
	subgroupid = ddinfo->subgroupid;

	ptrdiff_t npad_in_current_level = *npad;
	for(igroup = 1; igroup < nsubgroup; igroup ++){
		dest = (myid + subgroupsize * igroup + nid)%nid;
		src  = (myid - subgroupsize * igroup + nid)%nid;
		int targetsubgid = (subgroupid + igroup+nsubgroup)%nsubgroup;
		MPI_Sendrecv(&mylocalbox, sizeof(BoxMinMax), MPI_BYTE, src,0,
				&yourlocalbox, sizeof(BoxMinMax), MPI_BYTE, dest, 0,
				com, &status);
		ptrdiff_t sizeinbyte2send, sizeinbyte2recv;
#ifdef DEBUG
		DEBUGPRINT("P%d has mylocal %g %g : %g %g : %g %g\n", myid, mylocalbox.xmin, mylocalbox.xmax,
				mylocalbox.ymin, mylocalbox.ymax, mylocalbox.zmin, mylocalbox.zmax);
		DEBUGPRINT("P%d has targetlocalbox %g %g : %g %g : %g %g ::: P%d\n", myid, yourlocalbox.xmin, yourlocalbox.xmax,
				yourlocalbox.ymin, yourlocalbox.ymax, yourlocalbox.zmin, yourlocalbox.zmax, dest);
#endif
		void *buff2send = ExtractPtls2Send(base, nmem, *padbase,*npad, ddinfo,yourlocalbox, 
				width, box, &sizeinbyte2send, npad_in_current_level, gridinfo);
		MPI_Sendrecv(&sizeinbyte2send,1,MPI_INT64_T,dest,0, 
				&sizeinbyte2recv,1,MPI_INT64_T,src,0, 
				com,&status);
		*padbase = (void *) realloc( (*padbase),(*npad*n_size + sizeinbyte2recv));
		void *buff2get =  (void*)( (char*)(*padbase) + n_size * (*npad));
		my_MPI_Sendrecv(buff2send,sizeinbyte2send,MPI_BYTE, dest, 0, 
				buff2get, sizeinbyte2recv, MPI_BYTE, src, 0, 
				com, &status);
		*npad += sizeinbyte2recv/n_size;
		if(buff2send) free(buff2send);
	}
	DeleteOutofLocalBox(padbase, npad, n_size, width, &mylocalbox, &box, ddinfo->insidebox , gridinfo);
	if(subgroupsize >1) justppadding(base, nmem, padbase, npad, ddinfo+1, box,width, gridinfo);
}
void copyselfperiodic(void *base, ptrdiff_t nmem, void **padbase, ptrdiff_t *npad, DoDeInfo *ddinfo, SimBoxRange simbox,PosType width,
		GridInfo *gridinfo){
	int n_size = ddinfo->n_size;
	int (*insidebox)(GridInfo *, const void *, BoxMinMax *, PosType *,  SimBoxRange *, int , int) = ddinfo->insidebox;
	int (*edgeptl)(GridInfo *, const void *, SimBoxRange *, PosType *) = ddinfo->edgeptl;
	ptrdiff_t i,maxsize;
	char *run;
	BoxMinMax mylocalbox = ddinfo->lgroup.xyz;
#ifdef DEBUG
	DEBUGPRINT("P%d has mylocalbox %g %g : %g %g : %g %g\n", ddinfo->myid, mylocalbox.xmin, mylocalbox.xmax,
			mylocalbox.ymin, mylocalbox.ymax, mylocalbox.zmin, mylocalbox.zmax);
#endif

	maxsize = *npad + nmem * 0.1;
	ptrdiff_t stepsize = maxsize * 0.1;
	*padbase = (void*)realloc(*padbase, maxsize*n_size);

	run = (char*)base;
	for(i=0;i<nmem;i++){
		if(edgeptl(gridinfo, run,&simbox,&width)) {
			char buff[n_size];
			memmove(buff, run, n_size);
			int pflag;
			for(pflag=0;pflag<27;pflag++) {
				if(pflag==13) continue;
				else if(insidebox(gridinfo, buff, &mylocalbox, &width,&simbox,1, pflag)){
					if(*npad >=maxsize) {
						maxsize += stepsize;
						*padbase = (void*)realloc(*padbase, maxsize*n_size);
					}
					memmove((char*)(*padbase)+(*npad)*n_size, buff, n_size);
					(*npad) ++;
					memmove(buff, run, n_size);
				} 
			}
		}
		run += n_size;
	}
	*padbase = (void *)realloc(*padbase, *npad*n_size);
}
void old_ppadding(void *base,ptrdiff_t nmem, void **padbase, ptrdiff_t *npad, 
		DoDeInfo *ddinfo, int nddinfo, SimBoxRange box, PosType width, GridInfo *gridinfo){
#ifdef DEBUG
	if(ddinfo->myid==0 ){
		fprintf(stdout,"################################################\n");
		fprintf(stdout,"\n\n\nNow particle padding with width %g\n\n\n\n",width);
		fprintf(stdout,"################################################\n");
	}
#endif
	justppadding(base, nmem, padbase, npad, ddinfo, box, width, gridinfo);
	copyselfperiodic(base,nmem,padbase, npad, ddinfo + nddinfo-1, box,width, gridinfo);
}
void *SimpleExtractPtls2Send(void *padbase, ptrdiff_t npad, DoDeInfo *ddinfo, 
		BoxMinMax targetbox, PosType width, SimBoxRange simbox,
		ptrdiff_t *sizeinbyte2send, ptrdiff_t npad_in_current_level, GridInfo *gridinfo){
	char *run;
	ptrdiff_t i;
	ptrdiff_t isend = 0;
	int n_size = ddinfo->n_size;
	int (*insidebox)(GridInfo *, const void *, BoxMinMax *, PosType *,  SimBoxRange *, int , int) = ddinfo->insidebox;
	int (*edgeptl)(GridInfo *, const void *, SimBoxRange *, PosType *) = ddinfo->edgeptl;

	ptrdiff_t maxsize = MAX(npad_in_current_level,10000);
	ptrdiff_t stepsendsize = maxsize*0.2;
	void *sendbuff = (void *)malloc(maxsize*n_size);
	run = (char *)padbase;
	isend = 0;
	/* second last argument value of insidebox(gridinfo, run, &targetbox, &width,&simbox,0, 13) 
	 * has multiple meanging... If it is 0, then it deals with a pad particle, xPos does not applied and
	 * the second argument should be modified .*/
	for(i=0;i<npad_in_current_level;i++){
		if(insidebox(gridinfo, run, &targetbox, &width,&simbox,0, 13)) {
			memcpy((char*)sendbuff+isend*n_size, run, n_size);
			isend ++;
			if(isend >= maxsize) {
				maxsize += stepsendsize;
				sendbuff = (void *)realloc(sendbuff,maxsize*n_size);
			}
		}
		run += n_size;
	}
	sendbuff = (void*)realloc(sendbuff, isend*n_size);
	*sizeinbyte2send = isend*n_size;
	return sendbuff;
}
void adv_justppadding(void **padbase, ptrdiff_t *npad, 
		DoDeInfo *ddinfo, SimBoxRange box, PosType width, GridInfo *gridinfo){
	MPI_Status status;
	MPI_Comm com;
	int src, dest;
	ptrdiff_t nrec, nsend;
	ptrdiff_t i,j,k;
	int myid, nid;
	ptrdiff_t n_size = ddinfo->n_size;
	int igroup, subgroupsize, subgroupid;
	BoxMinMax mylocalbox,yourlocalbox;
	mylocalbox= ddinfo->lgroup.xyz;
	myid = ddinfo->myid;
	nid = ddinfo->nid;
	com = ddinfo->com;
	int nsubgroup = ddinfo->nsubgroup;
	subgroupsize = ddinfo->subgroupsize;
	subgroupid = ddinfo->subgroupid;

	ptrdiff_t npad_in_current_level = *npad;
	for(igroup = 1; igroup < nsubgroup; igroup ++){
		dest = (myid + subgroupsize * igroup + nid)%nid;
		src  = (myid - subgroupsize * igroup + nid)%nid;
		int targetsubgid = (subgroupid + igroup+nsubgroup)%nsubgroup;
		MPI_Sendrecv(&mylocalbox, sizeof(BoxMinMax), MPI_BYTE, src,0,
				&yourlocalbox, sizeof(BoxMinMax), MPI_BYTE, dest, 0,
				com, &status);
		ptrdiff_t sizeinbyte2send, sizeinbyte2recv;
		void *buff2send = SimpleExtractPtls2Send(*padbase,*npad, ddinfo,yourlocalbox, 
				width, box, &sizeinbyte2send, npad_in_current_level, gridinfo);
		MPI_Sendrecv(&sizeinbyte2send,1,MPI_INT64_T,dest,0, 
				&sizeinbyte2recv,1,MPI_INT64_T,src,0, 
				com,&status);
		*padbase = (void *) realloc( (*padbase),(*npad*n_size + sizeinbyte2recv));
		void *buff2get =  (void*)( (char*)(*padbase) + n_size * (*npad));
		my_MPI_Sendrecv(buff2send,sizeinbyte2send,MPI_BYTE, dest, 0, 
				buff2get, sizeinbyte2recv, MPI_BYTE, src, 0, 
				com, &status);
		*npad += sizeinbyte2recv/n_size;
		if(buff2send) free(buff2send);
	}
	DeleteOutofLocalBox(padbase, npad, n_size, width, &mylocalbox, &box, ddinfo->insidebox , gridinfo);
	if(subgroupsize >1) adv_justppadding(padbase, npad, ddinfo+1, box,width, gridinfo);
}
void DeleteInsideLocalBox(void **padbase, ptrdiff_t *npad, ptrdiff_t n_size, BoxMinMax *lbox, 
		SimBoxRange *simbox,
		int (*insidebox)(GridInfo *, const void *, BoxMinMax *, PosType *, SimBoxRange *,int , int),
			GridInfo *gridinfo){
	ptrdiff_t i;
	char *left, *right;
	PosType ww = 0;
	left = right = (char*)(*padbase);
	for(i=0;i<*npad;i++){
		if(!insidebox(gridinfo, right, lbox, &ww, simbox, 0, 13)) {
			memmove(left,right,n_size);
			left += n_size;
		}
		right += n_size;
	}
	*npad = (left  - (char*)(*padbase))/n_size;
	*padbase = (void *) realloc( *padbase, *npad * n_size);
}

void ppadding(void *base,ptrdiff_t nmem, void **padbase, ptrdiff_t *npad, 
		DoDeInfo *ddinfo, int nddinfo, SimBoxRange simbox, PosType width, GridInfo *gridinfo){
	BoxMinMax mybox = ddinfo[nddinfo-1].lgroup.xyz;
	BoxMinMax mycorebox = mybox;
	mycorebox.xmin += width;
	mycorebox.ymin += width;
	mycorebox.zmin += width;
	mycorebox.xmax -= width;
	mycorebox.ymax -= width;
	mycorebox.zmax -= width;
	BoxMinMax Bbox,Ebox;
	Bbox.xmin = -width;
	Bbox.ymin = -width;
	Bbox.zmin = -width;
	Bbox.xmax = simbox.x.max +width;
	Bbox.ymax = simbox.y.max +width;
	Bbox.zmax = simbox.z.max +width;
	/*
	Ebox.xmin = 0;
	Ebox.ymin = 0;
	Ebox.zmin = 0;
	Ebox.xmax = simbox.x.max;
	Ebox.ymax = simbox.y.max;
	Ebox.zmax = simbox.z.max;
	*/
	int n_size = ddinfo->n_size;
	int (*insidebox)(GridInfo *, const void *, BoxMinMax *, PosType *,  SimBoxRange *, int , int) = ddinfo->insidebox;
	int (*edgeptl)(GridInfo *, const void *, SimBoxRange *, PosType *) = ddinfo->edgeptl;
	ptrdiff_t i,maxsize;
	maxsize = 100+ nmem * 0.1;
	ptrdiff_t stepsize = maxsize * 0.1;
	*padbase = (void*)malloc(maxsize*n_size);
	PosType ww = 0;
	int iedge;

	char *run = (char*)base;
	for(i=0;i<nmem;i++){
		if(insidebox(gridinfo, run, &mybox, &ww, &simbox, 0, 13) && 
				!insidebox(gridinfo, run, &mycorebox, &ww, &simbox, 0, 13)){
			if(*npad >= maxsize - 26){
				maxsize += stepsize;
				*padbase = (void*)realloc(*padbase, maxsize*n_size);
			}
			memcpy((char*)(*padbase)+(*npad)*n_size, run, n_size);
			(*npad)++;
			int pflag;
			char buff[n_size];
			for(pflag =0; pflag < 27; pflag ++){
				if(pflag == 13) continue;
				memcpy(buff, run, n_size);
				if(
						/*
						insidebox(gridinfo, buff, &Bbox, &ww, &simbox, 0, pflag) 
						&&
						!insidebox(gridinfo,buff, &Ebox, &ww, &simbox,0,pflag)
						){
					insidebox(gridinfo, buff, &Bbox, &ww, &simbox, 1, pflag);
					*/
						insidebox(gridinfo, buff, &Bbox, &ww, &simbox, 1, pflag)
				){

					memcpy((char*)(*padbase)+(*npad)*n_size, buff, n_size);
					(*npad)++;
				}
			}
		}
		run += n_size;
	}
	*padbase = (void *)realloc(*padbase, *npad*n_size);
	adv_justppadding(padbase, npad,ddinfo, simbox, width, gridinfo);
	DeleteInsideLocalBox(padbase, npad, n_size, &mybox, &simbox, ddinfo->insidebox, gridinfo);
}



void buildrmscom(DoDeFunc *ddfunc, DoDeInfo *ddinfo, MPI_Comm Comm){
	ptrdiff_t i,j,k;
	int myid,nid;
	int subgroupsize,subgroupid;

	MPI_Comm_rank(Comm,&myid);
	MPI_Comm_size(Comm,&nid);

	int nsubgroup = getNextPrimeNumber(nid);
	if(nsubgroup == 1) {
		return;
	}
	subgroupsize = nid/nsubgroup;
	subgroupid = myid/subgroupsize;
	{
		MPI_Comm newcom;
		int key = myid % subgroupsize;
		MPI_Comm_split(Comm, subgroupid, key, &newcom);
		buildrmscom(ddfunc, ddinfo+1,newcom);
		if(subgroupsize == 1) MPI_Comm_free(&newcom);
	}
	return;
} 

