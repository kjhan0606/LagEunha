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
#include "mpirms.h"
#include "../MpiAux/mpiaux.h"


#define SWAP(a,b,stmp) do{\
	memcpy(stmp,a,n_size);\
	memcpy(a,b,n_size);\
	memcpy(b,stmp,n_size);\
}while(0)


#define MAX(a,b) ((a)>(b)?(a):(b))

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

	/*
	for(i=2;i<=number;){
		while((number % i) ==0){
			prime[nprime].prime = i;
			prime[nprime].factor ++;
			number = number/i;
		}
		if(prime[nprime].factor > 0) nprime ++;
		prime[nprime].factor = 0;
		if(i==2) i = 3;
		else i += 2;
	}
	*/

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

void ExtractLocalDomainVolume(DoDeInfo *ddinfo, int nddinfo, SimBox box){
	int i;
	ddinfo->lgroup.xyz.xmin = box.x.min; ddinfo->lgroup.xyz.ymin = box.y.min; ddinfo->lgroup.xyz.zmin = box.z.min; ddinfo->lgroup.xyz.wmin = box.w.min;
	ddinfo->lgroup.xyz.xmax = box.x.max; ddinfo->lgroup.xyz.ymax = box.y.max; ddinfo->lgroup.xyz.zmax = box.z.max; ddinfo->lgroup.xyz.wmax = box.w.max;
	for(i=0;i<nddinfo;i++){
		if(i!=0) ddinfo[i].lgroup.r = ddinfo[i-1].lgroup.r;
		int xyzoffset;
		if(ddinfo[i].xyzchip == 'x') xyzoffset = ddinfo[i].memoffset.xyz.x;
		else if(ddinfo[i].xyzchip == 'y') xyzoffset = ddinfo[i].memoffset.xyz.y;
		else if(ddinfo[i].xyzchip == 'z') xyzoffset = ddinfo[i].memoffset.xyz.z;
		else xyzoffset = ddinfo[i].memoffset.xyz.w;
		if(ddinfo[i].subgroupid==0) {
			ddinfo[i].lgroup.r.rmax[ddinfo[i].idirection] = *(PosType*)((char*)
					(&ddinfo[i].pivot[ddinfo[i].subgroupid*ddinfo[i].n_size])+xyzoffset);
		}
		else if(ddinfo[i].subgroupid != ddinfo[i].nsubgroup-1){
			ddinfo[i].lgroup.r.rmax[ddinfo[i].idirection] = *(PosType*)((char*)
					(&ddinfo[i].pivot[ddinfo[i].subgroupid*ddinfo[i].n_size])+xyzoffset);
			ddinfo[i].lgroup.r.rmin[ddinfo[i].idirection] = *(PosType*)((char*)
					(&ddinfo[i].pivot[(ddinfo[i].subgroupid-1)*ddinfo[i].n_size])+xyzoffset);
		}
		else {
			ddinfo[i].lgroup.r.rmin[ddinfo[i].idirection] = *(PosType*)((char*)
					(&ddinfo[i].pivot[(ddinfo[i].subgroupid-1)*ddinfo[i].n_size])+xyzoffset);
		}
	}
}


int example_muladd(const void *a, const void *b, float factor, int chip){
	float *aa = (float*)a;
	float *bb = (float*)b;
	if(isnan(*bb) || isnan(*aa)) {
		return 0;
	}
	else{
		(*aa) += (*bb)*factor;
		return 1;
	}
}

void MPI_Guess_kth_value(void *a, long nmem, long n_size,  long left, long right,
		DoDeFunc *ddfunc,
		void *result, int myid, int nid, MPI_Comm com){
	int isum=0;

	int (*muladd)(const void *, const void *, float, char);
	muladd = ddfunc->muladd;
	char xyzchip = ddfunc->xyzchip;

	char guess[nid*n_size];
	char myguess[n_size];
	long trypivot = (left+right)/2;

	char *aa = (char*)result;
	int i;
	for(i=0;i<n_size;i++) aa[i] = '\0'; /* Initialization is needed */
	if(trypivot <0 || left >= nmem-1 || right <1) { /* if it is negative member or pivot is beyond the range 
												   of the local domain*/
		float fact = NAN;
		muladd(myguess, myguess, fact,xyzchip);
	}
	else {
		memmove(myguess, (char*)a+trypivot*n_size,n_size); 
	}
	MPI_Gather(myguess, n_size, MPI_BYTE, guess, n_size, MPI_BYTE,0, com);

	if(myid==0) {
		int i;
		for(i=0;i<nid;i++) if(muladd( result, guess+i*n_size, 1,xyzchip)) isum ++;
		float norm = 1./isum;
		norm -= 1.;
		muladd( result, result, norm,xyzchip);
	}
	MPI_Bcast(result, n_size, MPI_BYTE, 0, com);
#ifdef DEBUG
	printf("P%d has pivot value %g trypivot %ld ::: %c \n",myid, (((AA*)result))->x, trypivot, xyzchip);
#endif
}

void Find_kth_smallest(void *base, long nmem, long n_size, 
		DoDeFunc *ddfunc,
		int kth, int npivot, void *result, int myid, int nid,
		MPI_Comm com){
	long ileft,iright,left,right,tleft;
	char swaptmp[n_size];
	char trypivotval[n_size];
	left = 0; right = nmem-1;
	long tnmem;
	MPI_Reduce(&nmem, &tnmem, 1, MPI_LONG, MPI_SUM, 0,com);
	MPI_Bcast(&tnmem, 1, MPI_LONG,0,com);
	long ntarget = (kth+1)*( (double)tnmem/(double)(npivot+1) );
	long Nresolution = tnmem * DomainError;
	char xyzchip = ddfunc->xyzchip;
	int (*compare)( const void *, const void *);
	if(xyzchip == 'x') compare = (ddfunc->xcompare);
	else if(xyzchip == 'y') compare = (ddfunc->ycompare);
	else if(xyzchip == 'z') compare = (ddfunc->zcompare);
	if(xyzchip == 'w') compare = (ddfunc->wcompare);
	while(1){
		ileft = left;
		iright = right;
		MPI_Guess_kth_value(base,nmem,n_size, left,right, ddfunc, trypivotval,myid,nid,com);
#ifdef DEBUG
		if(myid==0) printf("Determined pivot val = %g ::: %ld %ld \n",(((AA*)trypivotval)->x), left,right);
#endif
		do {
			while( compare( (char*)base+ileft*n_size, trypivotval) <=0 && ileft < nmem-1) ileft++;
			while( compare( trypivotval, (char*)base + iright*n_size)<=0 && iright >0) iright--;
			if(ileft<=iright){
				SWAP( (char*)base+ileft*n_size, (char*)base+iright*n_size, swaptmp);
				ileft++; iright--;
			}
		} while(ileft<iright);

		MPI_Reduce(&ileft, &tleft, 1, MPI_LONG, MPI_SUM,0, com);
		MPI_Bcast(&tleft, 1, MPI_LONG, 0, com);

		if(tleft > ntarget && llabs(tleft-ntarget) > Nresolution) right = iright;
		else if(tleft < ntarget && llabs(tleft-ntarget) > Nresolution) left = ileft;
		else break;
#ifdef DEBUG
		printf("P%d has tleft/ntarget %ld %ld  :: left/right %ld %ld\n",myid, tleft, ntarget, left,right);
#endif
	}
	memmove(result, trypivotval, n_size);
}

void Find_pivotvals(void *base, long nmem, long n_size, DoDeFunc *ddfunc,
		int npivot, void *pivot,int myid, int nid, MPI_Comm com){
	int kth;
	for(kth=0;kth<npivot;kth++){
		char *pos =  (char*)pivot + kth*n_size;
		Find_kth_smallest(base, nmem, n_size, ddfunc, kth, npivot,pos, myid, nid, com);
#ifdef DEBUG
		printf("P%d has results    %g\n",myid, ((AA*)pos)->x);
#endif
	}
	return;
}



void mpirms(void **ibase,long *mmem, long n_size, DoDeFunc *ddfunc, DoDeInfo *ddinfo,
		MPI_Comm Comm){
	MPI_Status status;
	int src,dest;
	long nrecv,nsend;
	long i,j,k;
	int myid,nid;
	long nmem = *mmem;
	void *base = *ibase;
	char swaptmp[n_size];
	int igroup,subgroupsize,subgroupid;
	long nowmem=0;
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

	int idir = ddfunc->divdir(base, nmem,Comm);
	if(idir==0) ddfunc->xyzchip = 'x';
	else if(idir==1) ddfunc->xyzchip = 'y';
	else if(idir==2) ddfunc->xyzchip = 'z';
	else  ddfunc->xyzchip = 'w';


#ifdef DEBUG
	printf("P%d divided along with %c\n",myid, ddfunc->xyzchip);
#endif



	Find_pivotvals(base, nmem, n_size, ddfunc,npivot,pivot,myid, nid, Comm);

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

		ddinfo->xdist = ddfunc->xdist;
		ddinfo->ydist = ddfunc->ydist;
		ddinfo->zdist = ddfunc->zdist;
		ddinfo->wdist = ddfunc->wdist;

		ddinfo->n_size = n_size;
		memmove(ddinfo->pivot, pivot, npivot*n_size);
	}
	int (*compare)(const void *, const void *);
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
				if(compare(left, tlocalmax)<0){
					right -= n_size;
					SWAP(left,right, swaptmp);
				}
				else left+=n_size;
			}
		}
		else if(targetsubgid == nsubgroup-1){
			memmove(tlocalmin, (char*)pivot+(targetsubgid-1)*n_size,n_size);
			for(;left<right;){
				if(compare(left, tlocalmin) >= 0){
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
					if(compare(left,tlocalmin)>=0 && compare(left,tlocalmax)<0){
						right -= n_size;
						SWAP(left,right, swaptmp);
					}
					else left += n_size;
				}
			}
		}
		long  nsend = nmem-((char*)right-(char *)base)/n_size;
		MPI_Sendrecv(&nsend,1,MPI_LONG,dest,0,&nrecv,1,MPI_LONG,src,0,Comm,&status);
		if(rbase == NULL) rbase = (void *)malloc(sizeof(char)*nrecv*n_size);
		else rbase = (void *)realloc(rbase,sizeof(char)*(nrecv + nowmem) *n_size);
		my_MPI_Sendrecv(right,nsend*sizeof(char)*n_size,MPI_BYTE,dest,0,
				(char*)rbase+nowmem*n_size,nrecv*sizeof(char)*n_size,MPI_BYTE,src,0,Comm,&status);
		nmem = nmem-nsend;
		base = (void*)realloc(base, nmem * n_size);
		nowmem += nrecv;
		if(myid==0) printf("+: %d : (%ld) <-> (%ld) %ld\n",igroup,nrecv, nsend,(long)nowmem); fflush(stdout);
	}
    rbase = (void *)realloc(rbase, sizeof(char)*(nowmem+nmem)*n_size);
    memmove((char*)rbase + nowmem*n_size, base, nmem*n_size);
    nowmem += nmem;
	free(base);

	{
		MPI_Comm newcom;
		int key = myid % subgroupsize;
		MPI_Comm_split(Comm, subgroupid, key, &newcom);
		mpirms(&rbase, &nowmem, n_size, ddfunc, ddinfo+1,newcom);
		/*
		MPI_Comm_free(&newcom);
		*/
	}
	*mmem = nowmem;
	*ibase =  rbase;
	return;
} 


void pmigrate(void **ibase,long *mmem, DoDeInfo *ddinfo){
	MPI_Status status;
	MPI_Comm com;
	int src,dest;
	long nrecv,nsend;
	long i,j,k;
	int myid,nid;
	long nmem = *mmem;
	void *base = *ibase;
	long n_size = ddinfo->n_size;
	char swaptmp[n_size];
	int igroup,subgroupsize,subgroupid;
	long nowmem=0;
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
	int (*compare)(const void *, const void *);
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
		left = (char*)base;
		right = (char*)base + nmem*n_size;

		if(targetsubgid==0){
			memmove(tlocalmax, (char*)pivot+targetsubgid*n_size,n_size);
			for(;left<right;){
				if(compare(left, tlocalmax)<0){
					right -= n_size;
					SWAP(left,right, swaptmp);
				}
				else left+=n_size;
			}
		}
		else if(targetsubgid == nsubgroup-1){
			memmove(tlocalmin, (char*)pivot+(targetsubgid-1)*n_size,n_size);
			for(;left<right;){
				if(compare(left, tlocalmin) >= 0){
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
					if(compare(left,tlocalmin)>=0 && compare(left,tlocalmax)<0){
						right -= n_size;
						SWAP(left,right, swaptmp);
					}
					else left += n_size;
				}
			}
		}
		long  nsend = nmem-((char*)right-(char *)base)/n_size;
		MPI_Sendrecv(&nsend,1,MPI_LONG,dest,0,&nrecv,1,MPI_LONG,src,0,com,&status);
		if(rbase == NULL) rbase = (void *)malloc(sizeof(char)*nrecv*n_size);
		else rbase = (void *)realloc(rbase,sizeof(char)*(nrecv + nowmem) *n_size);
		my_MPI_Sendrecv(right,nsend*sizeof(char)*n_size,MPI_BYTE,dest,0,
				(char*)rbase+nowmem*n_size,nrecv*sizeof(char)*n_size,MPI_BYTE,src,0,com,&status);
		nmem = nmem-nsend;
		base = (void*)realloc(base, nmem * n_size);
		nowmem += nrecv;
		if(myid==0) printf("+: %d : (%ld) <-> (%ld) %ld\n",igroup,nrecv, nsend,(long)nowmem); fflush(stdout);
	}
    rbase = (void *)realloc(rbase, sizeof(char)*(nowmem+nmem)*n_size);
    memmove((char*)rbase + nowmem*n_size, base, nmem*n_size);
    nowmem += nmem;
	free(base);

	if(subgroupsize > 1) pmigrate(&rbase, &nowmem, ddinfo+1);

	*mmem = nowmem;
	*ibase =  rbase;
	return;
} 


void justppadding(void *base,long nmem, void **padbase, long *npad, DoDeInfo *ddinfo, int nddinfo, SimBox box, PosType width){
	MPI_Status status;
	MPI_Comm com;
	int src,dest;
	long nrecv,nsend;
	long i,j,k;
	int myid,nid;
	long n_size = ddinfo->n_size;
	char swaptmp[n_size];
	int igroup,subgroupsize,subgroupid;
	char tsubgrpmin[n_size];
	char tsubgrpmax[n_size];
	char mylocalmin[n_size];
	char mylocalmax[n_size];

	myid = ddinfo->myid;
	nid = ddinfo->nid;
	com = ddinfo->com;

	int nsubgroup = ddinfo->nsubgroup;

	int npivot = ddinfo->npivot;
	void *pivot = ddinfo->pivot;


	int idir = ddinfo->idirection;

	subgroupsize = ddinfo->subgroupsize;
	subgroupid = ddinfo->subgroupid;
	char xyzchip = ddinfo->xyzchip;
	int (*distmeasure)(const void *, const void *, PosType,  Range *);
	int ioff;
	Range *simRange;
	if(xyzchip == 'x') {
		distmeasure = ddinfo->xdist;
		simRange = &(box.x);
		ioff = 0;
	}
	else if(xyzchip == 'y') {
		distmeasure = ddinfo->ydist;
		simRange = &(box.y);
		ioff = 1;
	}
	else if(xyzchip == 'z') {
		distmeasure = ddinfo->zdist;
		simRange = &(box.z);
		ioff = 2;
	}
	else {
		distmeasure = ddinfo->wdist;
		simRange = &(box.w);
		ioff = 3;
	}

	*(PosType*)(mylocalmin + ddinfo->memoffset.r[ioff]) = ddinfo->lgroup.r.rmin[ioff];
	*(PosType*)(mylocalmax + ddinfo->memoffset.r[ioff]) = ddinfo->lgroup.r.rmax[ioff];
	printf("P%d has min/max %g %g\n",myid,ddinfo->lgroup.r.rmin[ioff], ddinfo->lgroup.r.rmax[ioff]);

	for(igroup=1;igroup<nsubgroup;igroup++){
		dest = (myid + subgroupsize*igroup + nid)%nid;
		src =  (myid - subgroupsize*igroup + nid)%nid;
		char localmin[n_size], localmax[n_size];
		int targetsubgid = (subgroupid + igroup+nsubgroup)%nsubgroup;
		MPI_Sendrecv(mylocalmin,n_size,MPI_BYTE,src,1,tsubgrpmin,n_size,MPI_BYTE,dest, 1, com, &status);
		MPI_Sendrecv(mylocalmax,n_size,MPI_BYTE,src,1,tsubgrpmax,n_size,MPI_BYTE,dest, 1, com, &status);



		{
			char *left,*right;
			long nright;
	
			left = (char*)(*padbase);
			right = (char*)(*padbase) + (*npad)*n_size;


			for(;left<right;){
				if(distmeasure(left,tsubgrpmin,width,simRange)>=0 || distmeasure(left,tsubgrpmax,width,simRange)>0){
					right -= n_size;
					SWAP(left,right, swaptmp);
				}
				else left += n_size;
			}
			long nsend = (*npad) -((char*)right-(char*)(*padbase))/n_size;
			nright = ( right -(char*)(*padbase))/n_size;

			MPI_Sendrecv(&nsend,1,MPI_LONG, dest, 0, &nrecv, 1, MPI_LONG, src, 0, com, &status);
			printf("-P%d has min/max %g %g of P%d :::: nsend/nrecv %ld %ld :: range %g %g\n",
					myid,*(PosType*)(tsubgrpmin + ddinfo->memoffset.r[ioff]),
					*(PosType*)(tsubgrpmax + ddinfo->memoffset.r[ioff]), dest, nsend, nrecv, simRange->min, simRange->max);
			if((*padbase) == NULL) (*padbase) = (void *)malloc(sizeof(char)*nrecv*n_size);
			else (*padbase) = (void *)realloc((*padbase), sizeof(char)*(nrecv +(*npad))*n_size);
			right = ((char*) (*padbase) + nright*n_size);

			my_MPI_Sendrecv(right,nsend*sizeof(char)*n_size,MPI_BYTE,dest,0,
					(char*)(*padbase)+(*npad)*n_size,nrecv*sizeof(char)*n_size,MPI_BYTE,src,0,com,&status);
			(*npad) += nrecv;
		}
		{
			char *left,*right;
			left = (char*)base;
			right = (char*)base + nmem*n_size;

			for(;left<right;){
				if(distmeasure(left,tsubgrpmin,width, simRange)>=0 || distmeasure(left,tsubgrpmax,width ,simRange)>0)
				{
					right -= n_size;
					SWAP(left,right, swaptmp);
				}
				else left += n_size;
			}
			long  nsend = nmem-((char*)right-(char *)base)/n_size;
			MPI_Sendrecv(&nsend,1,MPI_LONG,dest,0,&nrecv,1,MPI_LONG,src,0,com,&status);

			printf("0P%d has min/max %g %g of P%d :::: nsend/nrecv %ld %ld :: range %g %g\n",
					myid,*(PosType*)(tsubgrpmin + ddinfo->memoffset.r[ioff]),
					*(PosType*)(tsubgrpmax + ddinfo->memoffset.r[ioff]), dest, nsend, nrecv, simRange->min, simRange->max);
	
			if( (*padbase) == NULL) (*padbase) = (void *)malloc(sizeof(char)*nrecv*n_size);
			else (*padbase) = (void *)realloc((*padbase),sizeof(char)*(nrecv + (*npad)) *n_size);
			my_MPI_Sendrecv(right,nsend*sizeof(char)*n_size,MPI_BYTE,dest,0,
					(char*)(*padbase)+(*npad)*n_size,nrecv*sizeof(char)*n_size,MPI_BYTE,src,0,com,&status);
			(*npad) += nrecv;
		}



		if(myid==0) printf("+: %d : (%ld) <-> (%ld) %ld\n",igroup,nrecv, nsend,(long)(*npad)); fflush(stdout);
	}

	if(subgroupsize > 1) justppadding(base, nmem, padbase, (void*)npad,ddinfo+1, nddinfo, box, width);
	return;
} 
void DeleteOutOfBound(void **padbase, long *npad, DoDeInfo *ddinfo, int nddinfo, SimBox box, PosType width){
	int (*xdistmeasure)(const void *, const void *, PosType,  Range *);
	int (*ydistmeasure)(const void *, const void *, PosType,  Range *);
	int (*zdistmeasure)(const void *, const void *, PosType,  Range *);
	int ioff;
	Range *xsimRange;
	Range *ysimRange;
	Range *zsimRange;
	xdistmeasure = ddinfo->xdist;
	xsimRange = &(box.x);
	ydistmeasure = ddinfo->ydist;
	ysimRange = &(box.y);
	zdistmeasure = ddinfo->zdist;
	zsimRange = &(box.z);
	char *left,*ok, *right;
	long n_size = ddinfo->n_size;

	char pmax[n_size], pmin[n_size];


	*(PosType*)(pmax+ddinfo->memoffset.xyz.x) = ddinfo[nddinfo-1].lgroup.r.rmax[0];
	*(PosType*)(pmax+ddinfo->memoffset.xyz.y) = ddinfo[nddinfo-1].lgroup.r.rmax[1];
	*(PosType*)(pmax+ddinfo->memoffset.xyz.z) = ddinfo[nddinfo-1].lgroup.r.rmax[2];
	*(PosType*)(pmin+ddinfo->memoffset.xyz.x) = ddinfo[nddinfo-1].lgroup.r.rmin[0];
	*(PosType*)(pmin+ddinfo->memoffset.xyz.y) = ddinfo[nddinfo-1].lgroup.r.rmin[1];
	*(PosType*)(pmin+ddinfo->memoffset.xyz.z) = ddinfo[nddinfo-1].lgroup.r.rmin[2];

	left = (char*)(*padbase);
	right = (char*)(*padbase) + (*npad)*n_size;
	ok = left;



	printf("-P%d has padding particles %ld\n",ddinfo->myid, *npad);

	for(;left<right;){
		if(xdistmeasure(left,pmin,width,xsimRange)>=0 || xdistmeasure(left,pmax,width,xsimRange)>0)
			if(ydistmeasure(left,pmin,width,ysimRange)>=0 || ydistmeasure(left,pmax,width,ysimRange)>0)
				if(zdistmeasure(left,pmin,width,zsimRange)>=0 || zdistmeasure(left,pmax,width,zsimRange)>0){
					memmove(ok, left, n_size);
					ok += n_size;
				}
		left += n_size;
	}
	*npad = ( ok -  (char*) (*padbase))/n_size;
	printf("+P%d has padding particles %ld\n",ddinfo->myid, *npad);
	(*padbase) = (void *)realloc( (*padbase), sizeof(char)*(*npad) * n_size);
}

void ppadding(void *base,long nmem, void **padbase, long *npad, DoDeInfo *ddinfo, int nddinfo, SimBox box, PosType width){
	justppadding(base,nmem, padbase, npad, ddinfo, nddinfo, box, width);
	DeleteOutOfBound(padbase, npad, ddinfo, nddinfo, box,width);
}
