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

	int (*muladd)(const void *, const void *, float, int);
	muladd = ddfunc->muladd;
	int xyzchip = ddfunc->xyzchip;

	char guess[nid*n_size];
	char myguess[n_size];
	long trypivot = (left+right)/2;

	char *aa = (char*)result;
	int i;
	for(i=0;i<n_size;i++) aa[i] = '\0'; /* Initialization is needed */
#ifdef DEBUG
	printf("P%d has pivot value %g trypivot %ld ::: %d : %d %d %d\n",myid, (((AA*)a)+trypivot)->x, trypivot, xyzchip,
			ddfunc->xchip, ddfunc->ychip, ddfunc->zchip);
#endif

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
	int (*compare)( const void *, const void *, int);
	int xyzchip = ddfunc->xyzchip;
	compare = (ddfunc->compare);
	while(1){
		ileft = left;
		iright = right;
		MPI_Guess_kth_value(base,nmem,n_size, left,right, ddfunc, trypivotval,myid,nid,com);
#ifdef DEBUG
		if(myid==0) printf("Determined pivot val = %g ::: %ld %ld \n",*(float*)((char*)trypivotval+xyzchip), left,right);
#endif
		do {
			while( compare( (char*)base+ileft*n_size, trypivotval,xyzchip ) <=0 && ileft < nmem-1) ileft++;
			while( compare( trypivotval, (char*)base + iright*n_size,xyzchip )<=0 && iright >0) iright--;
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
	int (*compare)(const void *, const void *, int);
	compare = ddfunc->compare;

	int idir = ddfunc->divdir(base, nmem,Comm);
	if(idir==0) ddfunc->xyzchip = ddfunc->xchip;
	else if(idir==1) ddfunc->xyzchip = ddfunc->ychip;
	else  ddfunc->xyzchip = ddfunc->zchip;



	Find_pivotvals(base, nmem, n_size, ddfunc,npivot,pivot,myid, nid, Comm);

	subgroupsize = nid/nsubgroup;
	subgroupid = myid/subgroupsize;
	int xyzchip = ddfunc->xyzchip;

	{ 
		ddinfo->com = Comm;
		ddinfo->myid = myid;
		ddinfo->nid = nid;
		ddinfo->nsubgroup = nsubgroup;
		ddinfo->npivot = npivot;
		ddinfo->idirection = idir;
		ddinfo->subgroupsize = subgroupsize;
		ddinfo->subgroupid = subgroupid;
		ddinfo->xyzchip = xyzchip;
		ddinfo->compare = compare;
		memmove(ddinfo->pivot, pivot, npivot*n_size);
	}

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
				if(compare(left, tlocalmax,xyzchip)<0){
					right -= n_size;
					SWAP(left,right, swaptmp);
				}
				else left+=n_size;
			}
		}
		else if(targetsubgid == nsubgroup-1){
			memmove(tlocalmin, (char*)pivot+(targetsubgid-1)*n_size,n_size);
			for(;left<right;){
				if(compare(left, tlocalmin,xyzchip) >= 0){
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
					if(compare(left,tlocalmin,xyzchip)>=0 && compare(left,tlocalmax,xyzchip)<0){
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
		MPI_Comm_free(&newcom);
	}
	*mmem = nowmem;
	*ibase =  rbase;
	return;
} 

void migrate(void **ibase,long *mmem, long n_size, DoDeInfo *ddinfo, MPI_Comm Comm){
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

	myid = ddinfo->myid;
	nid = ddinfo->nid;

	int nsubgroup = ddinfo->nsubgroup;
	if(nsubgroup == 1) {
		return;
	}

	int npivot = ddinfo->npivot;
	void *pivot = ddinfo->pivot;

	int (*compare)(const void *, const void *, int);
	compare = ddinfo->compare;

	int idir = ddinfo->idirection;

	subgroupsize = ddinfo->subgroupsize;
	subgroupid = ddinfo->subgroupid;
	int xyzchip = ddinfo->xyzchip;



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
				if(compare(left, tlocalmax,xyzchip)<0){
					right -= n_size;
					SWAP(left,right, swaptmp);
				}
				else left+=n_size;
			}
		}
		else if(targetsubgid == nsubgroup-1){
			memmove(tlocalmin, (char*)pivot+(targetsubgid-1)*n_size,n_size);
			for(;left<right;){
				if(compare(left, tlocalmin,xyzchip) >= 0){
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
					if(compare(left,tlocalmin,xyzchip)>=0 && compare(left,tlocalmax,xyzchip)<0){
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

	MPI_Comm newcom = ddinfo->com;
	migrate(&rbase, &nowmem, n_size, ddinfo+1,newcom);

	*mmem = nowmem;
	*ibase =  rbase;
	return;
} 

