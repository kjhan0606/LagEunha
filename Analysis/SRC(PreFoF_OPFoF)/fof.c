#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "fof.h"
#include "Memory.h"
#include "Time.h"
/*
#define N 116
*/
float ran2(long *);
void check_tree();
particle *p;
void force_spline();
extern double a[8];
float epsilon2=0.04;
float fof_link = 0.2;
double f(double);
int main(int argc,char *argv[]){
	double std,mean;
	int ntmp;
	int num,ii;
	float tmpx,tmpy,tmpz,dist2;
	float fplmf,ptlmass;
	int i,j,k;
	int N,M;
	FoFTPtlStruct *ptl;
	FoFTPtlStruct *ptr;
	particle *linked;
	int nfof;
	FoFBeginEndTree beginend;
	Box box;
	FoFTStruct *TREE;
	float theta = 1.;
	float wtime;
	long iseed=-9;
	FILE *fp;
	(void) ran2(&iseed);
	/*
	i_potent_spline();
	*/
	for(k=50000;k<10000000;k+=10000){
		N = k;
		M = N;
		ptl = (FoFTPtlStruct *) malloc(sizeof(FoFTPtlStruct)*N);
		p = (particle *) malloc(sizeof(particle)*M);
		linked = (particle *) malloc(sizeof(FoFTPtlStruct)*M);
		TREE = (FoFTStruct *) malloc(sizeof(FoFTStruct)*N);
		wtime = WALLCLOCK();
		for(i=0;i<N;i++){
			ptl[i].type = TYPE_PTL;
			p[i].x = ptl[i].r[0] = ran2(&iseed)*8;
			p[i].y = ptl[i].r[1] = ran2(&iseed)*8;
			p[i].z = ptl[i].r[2] = ran2(&iseed)*8;
			ptl[i].sibling = &ptl[i+1];
			ptl[i].indx = i;
		}
		for(i=N;i<M;i++){
			p[i].x = ran2(&iseed)*7;
			p[i].y = ran2(&iseed)*7;
			p[i].z = ran2(&iseed)*7;
		}
		ptl[N-1].sibling = NULL;
		/* */
		(void)WALLCLOCK();
		/* */
		box.x = box.y = box.z = 0.;
		box.width = 8.;
		FoF_Make_Tree(TREE,ptl,N,box);
		printf("%g  ellapsed ..",WALLCLOCK());
		printf("%d :",N);
		/*
		check_tree(TREE,ptl);
		*/
		num=new_fof_link(p,fof_link,TREE,ptl,linked);
		/*
		num = 0;
		for(ii=0;ii<N;ii++) if(ptl[ii].included == YES) num ++;
		*/
		printf("%d linked ",num);
		/*
		for(i=0;i<M;i++){
			potent[i] = treeplumpotential(p+i,theta,TREE,ptl,epsilon);
		}
		*/
		printf("%g  ellapsed \n",WALLCLOCK());
		free(TREE);free(ptl);free(p);free(linked);
	}
	fprintf(stdout,"Terminated with success\n");fflush(stdout);
	return 1;
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
