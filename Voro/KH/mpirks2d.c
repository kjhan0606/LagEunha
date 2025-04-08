typedef struct PrimeNumber{
	int prime, factor;
} PrimeNumber;


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

void buildrks2dcom(DDinfo2D *ddinfo, MPI_Comm com){
	ptrdiff_t i,j,k;
	int myid, nid;

	MPI_Comm_rank(com, &myid);
	MPI_Comm_size(com, &nid);


	int nSubGroup = getNextPrimeNumber(nid);
	if(nSubGroup ==1) return;
	int subGroupSize = nid/nSubGroup;
	int subGroupId = myid/subGroupSize;

	MPI_Comm newcom;
	int key = myid % subGroupSize;
	MPI_Comm_split(com, subGroupId, key, &newcom);
	buildrks2dcom(ddinfo+1, newcom);
	if(subGroupSize == 1) MPI_Comm_free(&newcom);
	return;
}

void mpirks2D(void **base, ptrdiff_t *mmem, ptrdiff_t n_size, 
		DDFunc2D *ddfunc, DDInfo2D *ddinfo,
		MPI_Comm com, GridInfo *gridinfo, 
		enum DDdirection dd_direction, enum RKSComBuild newrksbuild){
	MPI_Status status;
	int src,dest;
	ptrdiff_t nrecv, nsend;
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

	int idir;
	if(dd_direction == newway){
        idir = ddfunc->divdir(gridinfo, base, nmem,Comm);
        if(idir==0) ddfunc->xyzchip = 'x';
        else if(idir==1) ddfunc->xyzchip = 'y';
        else if(idir==2) ddfunc->xyzchip = 'z';
        else  ddfunc->xyzchip = 'w';
    }
    else if(dd_direction == oldway){
        ddfunc->xyzchip = ddinfo->xyzchip;
        if(ddinfo->xyzchip == 'x') idir = 0;
        else if(ddinfo->xyzchip == 'y') idir = 1;
        else if(ddinfo->xyzchip == 'z') idir = 2;
        else if(ddinfo->xyzchip == 'w') idir = 3;
    }
	Find_pivotvals(base, nmem, n_size, ddfunc,npivot,pivot,myid, nid, Comm, gridinfo);
}
