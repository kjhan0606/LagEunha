#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<limits.h>
#include<sys/types.h>
#include<unistd.h>
#include<mpi.h>
#include "eunha.h"
#include "mpirks.h"
#include "test.h"






float ran2(long *);





#define MIN(a,b) ((a)<(b)? (a):(b))
#define MAX(a,b) ((a)>(b)? (a):(b))


AA *bp;


/* Global value for the RMS */
DoDeInfo *ddinfo;
int nddinfo;
DoDeFunc ddfunc;

/* a macro to generate a function with name "xcompare" for 
 * a type of AA and its member "x", whose type is float */
GenerateMemberCompare(xcompare,AA,x,PosType);
GenerateMemberCompare(ycompare,AA,y,PosType);
GenerateMemberCompare(zcompare,AA,z,PosType);

GenerateInsideRange(xpinner,AA,x,PosType);
GenerateInsideRange(ypinner,AA,y,PosType);
GenerateInsideRange(zpinner,AA,z,PosType);

int AAdivdir3d(GridInfo *gridinfo, const void *a, ptrdiff_t nmem, MPI_Comm com){
	ptrdiff_t i;
	double xmean,xstd;
	double ymean,ystd;
	double zmean,zstd;
	double totalval;
	ptrdiff_t totalnum;
	AA *aa = (AA*)a;
	xmean = xstd = 0;
	ymean = ystd = 0;
	zmean = zstd = 0;
	for(i=0;i<nmem;i++){
		xmean += xPos(gridinfo, aa+i);
		ymean += yPos(gridinfo, aa+i);
		zmean += zPos(gridinfo, aa+i);
		xstd += xPos(gridinfo, aa+i)*xPos(gridinfo, aa+i);
		ystd += yPos(gridinfo, aa+i)*yPos(gridinfo, aa+i);
		zstd += zPos(gridinfo, aa+i)*zPos(gridinfo, aa+i);
	}
	MPI_Reduce(&xmean,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);
	xmean = totalval;

	MPI_Reduce(&ymean,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);
	ymean = totalval;

	MPI_Reduce(&zmean,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);
	zmean = totalval;

	MPI_Reduce(&xstd,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);
	xstd = totalval;

	MPI_Reduce(&ystd,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);
	ystd = totalval;

	MPI_Reduce(&zstd,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);
	zstd = totalval;

	MPI_Reduce(&nmem,&totalnum,1, MPI_INT64_T, MPI_SUM, 0, com);
	MPI_Bcast(&totalnum, 1, MPI_INT64_T, 0, com);
	xmean = xmean/totalnum;
	ymean = ymean/totalnum;
	zmean = zmean/totalnum;
	xstd = xstd/totalnum - xmean*xmean;
	ystd = ystd/totalnum - ymean*ymean;
	zstd = zstd/totalnum - zmean*zmean;
	if(xstd >= ystd && xstd >= zstd) return 0;
	else if(zstd >= ystd && zstd >= xstd) return 2;
	else return 1;
}

int AAdetermineEdgePtl(GridInfo *gridinfo, const void *a, SimBoxRange *simbox, PosType *width){
	AA *aa = ( (AA*) a) ;
	if(xPos(gridinfo, aa)<=*width) return 1;
	else if ( (simbox->x).max - xPos(gridinfo,aa) < *width) return 1;
	if(yPos(gridinfo, aa)<=*width) return 1;
	else if ( (simbox->y).max - yPos(gridinfo, aa) < *width) return 1;
	if(zPos(gridinfo,aa)<=*width) return 1;
	else if ( (simbox->z).max - zPos(gridinfo,aa) < *width) return 1;
	return 0;
}

int AAInSideBox(GridInfo *gridinfo, const void *a, BoxMinMax *lbox, PosType *width, SimBoxRange *simbox, int mflag , int pflag){
	AA *aa = ((AA*)a);
	int i,j,k;
	if(pflag == 13) {
		i = j = k = 0;
	}
	else {
		i = pflag%3-1;
		j = (pflag%9)/3 -1;
		k = (pflag/9)-1;
	}
	PosType ax,ay,az,pbxmin,pbymin,pbzmin,pbxmax,pbymax,pbzmax;
	/*
	if(mflag)
	*/
	{
		ax = xPos(gridinfo,aa);
		ay = yPos(gridinfo,aa);
		az = zPos(gridinfo,aa);
	}
	/*
	else {
		ax = aa->x;
		ay = aa->y;
		az = aa->z;
	}
	*/
	ax +=  i*((simbox->x).max);
	ay +=  j*((simbox->y).max);
	az +=  k*((simbox->z).max);
	if( (ax >= lbox->xmin -(*width) && ax < lbox->xmax +(*width)) && 
			(ay >= lbox->ymin -(*width) && ay <lbox->ymax +(*width)) &&
			(az >= lbox->zmin -(*width) && az <lbox->zmax +(*width))) {
		if(mflag) { 
			aa->x += i*(simbox->x.max); 
			aa->y += j*(simbox->y.max); 
			aa->z += k*(simbox->z.max); 
		} 
		return 1; 
	} 
	else return 0;
}
int AAmuladd(GridInfo *gridinfo, const void *a, const void *b, float fact, char xyz, int iflag){
	PosType aa;
	float *av;
	PosType bb;
	switch (xyz) {
		case 'x': 
			av = &(((AA*)a)->x);
			if(iflag<0) bb = ( (AA*)b)->x;
			else bb = xPos(gridinfo, (AA*)b);
			break;
		case 'y': 
			av = &(((AA*)a)->y);
			if(iflag<0) bb = ( (AA*)b )->y;
			else bb = yPos(gridinfo, (AA*)b);
			break;
		default: 
			av = &(((AA*)a)->z);
			if(iflag<0) bb = ( (AA*)b )->z;
			else bb = zPos(gridinfo, (AA*)b);
	}
	aa = (PosType)(*av);
   	if(iflag){
		if(isnan(aa) || isnan(bb)) return 0;
		else {
#ifdef DEBUG
			printf("###### %g %g %g :::; %p %p %ld\n",*av, bb, fact, a, b, PINDX( (AA*)a));
#endif
			(*av) += bb*fact;
			return 1;
		}
	}
	else {
		Range *bbb = (Range*)b;
		(*av) += bbb->max*fact;
	}
	return 1;
}
int WRONGAAmuladd(GridInfo *gridinfo, const void *a, const void *b, float fact, char xyz, int iflag){
	PosType aa;
	float *av;
	PosType bb;
	switch (xyz) {
		case 'x': 
			av = &(((AA*)a)->x);
			aa = xPos(gridinfo, (AA*)a);
			bb = xPos(gridinfo, (AA*)b);
			break;
		case 'y': 
			av = &(((AA*)a)->y);
			aa = yPos(gridinfo, (AA*)a);
			bb = yPos(gridinfo, (AA*)b);
			break;
		default: 
			av = &(((AA*)a)->z);
			aa = zPos(gridinfo, (AA*)a);
			bb = zPos(gridinfo, (AA*)b);
	}
    if(iflag){
		if(isnan(aa) || isnan(bb)) return 0;
		else {
#ifdef DEBUG
			printf("###### %g %g %g :::; %p %p %ld\n",*av, bb, fact, a, b, PINDX( (AA*)a));
#endif
			(*av) += bb*fact;
			return 1;
		}
	}
	else {
		Range *bbb = (Range*)b;
		(*av) += bbb->max*fact;
	}
	return 1;
}
int oldAAmuladd(GridInfo *gridinfo, const void *a, const void *b, float fact, char xyz, int iflag){
	register float *aa, *bb;
	switch (xyz) {
		case 'x': 
			aa = &(((AA*)a)->x); 
			bb = &(((AA*)b)->x); 
			break;
		case 'y': 
			aa = &(((AA*)a)->y); 
			bb = &(((AA*)b)->y);
			break;
		default: 
			aa = &(((AA*)a)->z); 
			bb = &(((AA*)b)->z);
	}
	if(iflag){
		if(isnan(*aa) || isnan(*bb)) return 0;
		else {
			(*aa) += (*bb)*fact;
			return 1;
		}
	}
	else {
		Range *bbb = (Range*)b;
		(*aa) += bbb->max*fact;
	}
	return 1;
}



int myid,nid;


void check(AA *bp, size_t np, int myid){
	size_t i;
	float min, max;
	min = 1.E20;
	max = -1.E20;
	for(i=0;i<np;i++){
		min = MIN(min,bp->x);
		max = MAX(max,bp->x);
		bp ++;
	}
	printf("P%d has min/max= %g %g with np = %ld\n",myid, min,max,(ptrdiff_t)np);
}
void checkgrid(GridInfo *outgrid){
	DenType min, max;
	min = 1.e20;
	max = -1.e20;
	ptrdiff_t i;
	ptrdiff_t ic,jc;
	ic = jc = 0;
	DenType *oden = (DenType*)(outgrid + 1);
	for(i=0;i<outgrid->npix;i++){
		min = MIN(min, oden[i]);
		max = MAX(max, oden[i]);
		if(oden[i] != 2) {
			/*
			printf("P%d Strange Grid value at %ld\n",myid,i);
			*/
			jc ++;
		}
		else ic ++;
	}
	printf("P%d grid value: min/max = %g %g :: correct/incorrect count %ld/%ld\n",myid,min,max, ic,jc);
}


int main(int argc, char **argv){
	int i,j,k;
	long iseed;


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	{ 
#ifdef DEBUG
		char hostname[190]; 
		hostname[189] = '\0'; 
		gethostname(hostname,189); 
		printf("P%d host name is %s with pid= %d\n",myid,hostname, getpid());
		if(myid==-1){
			int kkk = 1;
			while(kkk) {
				kkk = 1;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	}



	iseed = -(myid+1);

	int nx,ny,nz;
	nx = ny = nz = 128;
	ptrdiff_t np=nx*ny*nz;

	bp = (AA*)malloc(sizeof(AA)*np);

	for(i=0;i<np;i++){
		bp[i].imem = (int) (10000000.0 * (rand() / (RAND_MAX + 1.0)));
		bp[i].x = nx*ran2(&iseed);
		bp[i].y = ny*ran2(&iseed);
		bp[i].z = nz*ran2(&iseed);
		CHANGEINDX(bp+i,0);
	}


	int nprime;
	PrimeNumber prime[100];
	{
		nprime = getprimenumber(nid,prime);
		for(i=0;i<nprime;i++){
			if(myid==0) printf("prime number : %d^%d\n",prime[i].prime,prime[i].factor);
		}
	}



	double stime,etime;
	stime = MPI_Wtime();



	/* macro making space to save ddinfo data */
	SimBoxRange box;
	box.x.min= box.y.min= box.z.min= 0;
	box.x.max= nx; box.y.max= ny; box.z.max= nz;
	MakeDoDeInfo3D(nid,prime, AA, x,y,z,ddinfo, nddinfo);
	/* bp: data, np: number of member, AA: type, (x,y,z): position member, 
	 * ddfunc: input argument with functions,
	 * muladd: multiply and add function for the type,
	 * compare: comparison of the type,
	 * divdir3d: a rule to define how to subdivide,
	 * MPI_COMM_WORLD: top-most communicator */
	MakeDoDeFunc3D(AA,ddfunc,xcompare,ycompare,zcompare,xpinner,ypinner,zpinner,muladd,divdir3d, 
			InSideBox,determineEdgePtl,ddinfo);

	GridInfo gridinfo;
	gridinfo.nx = nx; gridinfo.ny = ny; gridinfo.nz = nz; gridinfo.nxny = nx*ny;

	/*
	BuildSimpleRMS( box, sizeof(AA), &ddfunc,ddinfo,MPI_COMM_WORLD);
	pmigrate((void **)(&bp),&np, ddinfo, &gridinfo);
	*/
	mpirks( (void**)(&bp), &np, sizeof(AA), &ddfunc, ddinfo,MPI_COMM_WORLD, &gridinfo, newway, NewBuildCom);
	etime = MPI_Wtime();
	if(myid==0) printf("\n\nwall clocktime for mpi RMS DD: %g second\n\n\n",etime-stime);

	ExtractLocalDomainVolume(ddinfo, nddinfo,box);

	printf("**P%d has %g <= x < %g : %g <= y < %g : %g <=z < %g with np= %ld\n",myid,
			ddinfo[nddinfo-1].lgroup.xyz.xmin, ddinfo[nddinfo-1].lgroup.xyz.xmax,
			ddinfo[nddinfo-1].lgroup.xyz.ymin, ddinfo[nddinfo-1].lgroup.xyz.ymax,
			ddinfo[nddinfo-1].lgroup.xyz.zmin, ddinfo[nddinfo-1].lgroup.xyz.zmax, np);

	stime = MPI_Wtime();
	pmigrate((void **)(&bp),&np, ddinfo, &gridinfo);
	etime = MPI_Wtime();
	if(myid==0) printf("\n\nwall clocktime for mpi Particle migration: %g second\n\n\n",etime-stime);


	AA *ppad = NULL;
	ptrdiff_t npad = 0;
	PosType width = 4.;
	stime = MPI_Wtime();

	ppadding(bp, np, (void**)(&ppad), &npad, ddinfo, nddinfo, box, width, &gridinfo);

	etime = MPI_Wtime();
	if(myid==0) printf("\n\nwall clocktime for mpi particle padding: %g second\n\n\n",etime-stime);
	stime = MPI_Wtime();
	{
		PosType ymin,ymax,xmin,xmax,zmin,zmax;
		xmin = ymin =zmin = 1.e20; 
		xmax = ymax =zmax = -1e20;
		for(i=0;i<npad;i++){
			xmin = MIN(xmin,ppad[i].x);
			xmax = MAX(xmax,ppad[i].x);
			ymin = MIN(ymin,ppad[i].y);
			ymax = MAX(ymax,ppad[i].y);
			zmin = MIN(zmin,ppad[i].z);
			zmax = MAX(zmax,ppad[i].z);
		}
		printf("##$$ P%d has padding p: rmin/max = %g %g . %g %g . %g %g ::: %ld\n",myid, xmin,xmax,ymin,ymax,zmin,zmax,npad);
	}

	etime = MPI_Wtime();
	if(myid==0) printf("\n\nwall clocktime for mpi multisection sorting: %g second\n\n\n",etime-stime);
	check(bp,np, myid);
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);




	stime = MPI_Wtime();
	GridInfo *ingrid, *outgrid;
	ingrid = (GridInfo*)malloc(sizeof(GridInfo));
	outgrid = (GridInfo*)malloc(sizeof(GridInfo));
	ingrid->nx = nx;
	ingrid->ny = ny;
	ingrid->nz = nz;
	outgrid->nx = nx;
	outgrid->ny = ny;
	outgrid->nz = nz;
	/*
	ingrid->ix = ingrid->iy = 0;
	ingrid->jx = ingrid->jy = nx-1;
	ingrid->iz = nx*(ddinfo->myid/(float)(ddinfo->nid));
	ingrid->jz = nx*((ddinfo->myid+1)/(float)(ddinfo->nid))-1;

	outgrid->ix = outgrid->iz = 0;
	outgrid->jx = outgrid->jz = nx-1;
	outgrid->iy = nx*(ddinfo->myid/(float)(ddinfo->nid));
	outgrid->jy = nx*((ddinfo->myid+1)/(float)(ddinfo->nid))-1;
	*/
	ingrid->ix = ddinfo[nddinfo-1].lgroup.xyz.xmin;
	ingrid->iy = ddinfo[nddinfo-1].lgroup.xyz.ymin;
	ingrid->iz = ddinfo[nddinfo-1].lgroup.xyz.zmin;
	ingrid->jx = ddinfo[nddinfo-1].lgroup.xyz.xmax;
	ingrid->jy = ddinfo[nddinfo-1].lgroup.xyz.ymax;
	ingrid->jz = ddinfo[nddinfo-1].lgroup.xyz.zmax;
	ingrid->jx --;
	ingrid->jy --;
	ingrid->jz --;

	int nbuff=4;
	outgrid->ix = outgrid->iy = 0;
	outgrid->jx = outgrid->jy = nz-1;
	outgrid->iz = nx*(ddinfo->myid/(float)(ddinfo->nid));
	outgrid->jz = nx*((ddinfo->myid+1)/(float)(ddinfo->nid))-1;

	outgrid->ix -= nbuff;
	outgrid->iy -= nbuff;
	outgrid->iz -= nbuff;
	outgrid->jx += nbuff;
	outgrid->jy += nbuff;
	outgrid->jz += nbuff;



	printf("GI:: %d has grid %d %d : %d %d : %d %d\n",myid, ingrid->ix,ingrid->jx, ingrid->iy,ingrid->jy, ingrid->iz, ingrid->jz);

	getingridnpix(ingrid);
	ingrid = realloc(ingrid, sizeof(GridInfo)+ ingrid->npix * sizeof(DenType));
	DenType *iden = (DenType*)(ingrid + 1);

	for(i=0;i<ingrid->npix;i++) iden[i] = 2.;

	outgrid = getoutgridnpix(outgrid, 0.);

	gmigrate(ingrid,outgrid, ddinfo, nddinfo);
	printf("GO:: %d has grid %d %d : %d %d : %d %d\n",myid, outgrid->ix,outgrid->jx, outgrid->iy,outgrid->jy, outgrid->iz, outgrid->jz);
	fflush(stdout);
	etime = MPI_Wtime();
	if(myid==0) printf("\n\nwall clocktime for mpi multisection grid sending: %g second\n\n\n",etime-stime);
	checkgrid(outgrid);
	MPI_Finalize(); 
}
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
