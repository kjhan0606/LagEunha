#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "ramses.h"
#include "tree.h"
#include "near.h"
#include "params.h"
#include "Memory.h"
#define MIN(a,b) ((a)<(b)? (a): (b))
#define MAX(a,b) ((a)>(b)? (a): (b))
double f(double);

/*
float x,y,z;
*/

/* inverse of pi */
#define AMPW41 (-0.2387324)
#define AMPW42 (0.07957747)


float W4(float dist,float h){
	float twomr,r;
	r = dist/h;
	if(r<1.){
		/*
		return 1/PI*(1.-1.5*dist*dist+0.75*dist*dist*dist);
		*/
		return AMPW41*((2.-r)*r*r-1.333333)/(h*h*h);
	}
	else if(r < 2.){
		twomr = (2.-r)/h;
		/*
		return 1/PI*0.25*twomr*twomr*twomr;
		*/
		return AMPW42*twomr*twomr*twomr;
	}
	else {
		return 0.;
	}
}
#undef AMPW41
#undef AMPW42
void findsphdensity(SimpleBasicParticleType *bp,int np,int *nearindex, int Numnear, 
		float *densph){
	float *h;
	double std,mean;
	int ntmp;
	float tmpx,tmpy,tmpz;
	float *dist2, *mass;
	float fplmf,ptlmass;
	size_t i,j,k;
	int N,M;
	TPtlStruct *ptl;
	BeginEndTree beginend;
	Box box;
	TStruct *TREE;
	float theta = 1.;
	float wtime;
	long iseed=-9;
	float x0,y0,z0,pscale;
	float Rmin,Rmax,size;
	int index;
	particle *p;


	p = (particle *) Malloc(sizeof(particle)*np,PPTR(p));
	for(i=0;i<np;i++){
		p[i].x = bp[i].x;
		p[i].y = bp[i].y;
		p[i].z = bp[i].z;
	}
	ptl = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
	for(i=0;i<np;i++){
		ptl[i].type = TYPE_PTL;
		ptl[i].x = bp[i].x;
		ptl[i].y = bp[i].y;
		ptl[i].z = bp[i].z;
		ptl[i].sibling = ptl+(i+1);
		ptl[i].mass = bp[i].mass;
	}
	ptl[np-1].sibling = NULL;
	TREE = (TStruct *) Malloc(sizeof(TStruct)*MAX(np/2,10000),PPTR(TREE));
	h = (float*) Malloc(sizeof(float)*np,PPTR(h));
	dist2 = (float*) Malloc(sizeof(float)*np*Numnear,PPTR(dist2));
	mass = (float*) Malloc(sizeof(float)*np*Numnear,PPTR(mass));

	/*
	(void) WALLCLOCK();
	*/
	box = findsphbox(ptl,np);
#ifdef DEBUG
	printf("after findphsbox rmin =%g %g %g with width=%g for np = %d\n",box.x,box.y,box.z,box.width, np);
#endif
	Make_Tree_Near(TREE,ptl,np,box);
#ifdef DEBUG
	printf("after Make_Tree_near\n");
#endif
	/*
	printf("CPU(Make_Tree_test)=%g\n",WALLCLOCK());
	*/
	float mindist;
	mindist = 2.e25;
#pragma omp parallel for private(i,j,k) schedule(guided)
	for(i=0;i<np;i++){
		int tmpindx[Numnear];
		int res;
		float neardist;
		k = i*Numnear;
		res = Find_Near(p+i,Numnear,TREE,ptl,&neardist,tmpindx,dist2+k, mass+k);
		h[i] = neardist;
		densph[i] = 0;
		for(j=0;j<res;j++){
			nearindex[k+j] = tmpindx[j];
			dist2[k+j] = sqrt(dist2[k+j]);
			densph[i] += W4(dist2[k+j],h[i]/2.);
		}
	}
	Free(mass);Free(dist2); Free(h); 
	Free(TREE);Free(ptl);Free(p);
	/*
	printf("ins=", tmpindx);
	fprintf(stdout,"Terminated with success\n");fflush(stdout);
	*/
#ifdef DEBUG
	printf("exiting findsphdensity\n");
#endif
}
void starfindsphdensity(SimpleBasicParticleType *bp,int np,int *nearindex, int Numnear, 
		float *densph){
	float *h;
	double std,mean;
	int ntmp;
	float tmpx,tmpy,tmpz;
	float *dist2, *mass;
	float fplmf,ptlmass;
	size_t i,j,k;
	int N,M;
	TPtlStruct *ptl;
	BeginEndTree beginend;
	Box box;
	TStruct *TREE;
	float theta = 1.;
	float wtime;
	long iseed=-9;
	float x0,y0,z0,pscale;
	float Rmin,Rmax,size;
	int index;
	int mp;
	particle *p;


	p = (particle *) Malloc(sizeof(particle)*np,PPTR(p));
	for(i=0;i<np;i++){
		p[i].x = bp[i].x;
		p[i].y = bp[i].y;
		p[i].z = bp[i].z;
	}

	ptl = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
	mp = 0;
	for(i=0;i<np;i++){
		if(bp[i].type == TYPE_STAR){
			ptl[mp].type = TYPE_PTL;
			ptl[mp].x = bp[i].x;
			ptl[mp].y = bp[i].y;
			ptl[mp].z = bp[i].z;
			ptl[mp].sibling = ptl+(mp+1);
			ptl[mp].mass = bp[i].mass;
			mp++;
		}
	}
	ptl[mp-1].sibling = NULL;
	TREE = (TStruct *) Malloc(sizeof(TStruct)*MAX(np/2,10000),PPTR(TREE));
	h = (float*) Malloc(sizeof(float)*np,PPTR(h));
	dist2 = (float*) Malloc(sizeof(float)*np*Numnear,PPTR(dist2));
	mass = (float*) Malloc(sizeof(float)*np*Numnear,PPTR(mass));

	box = findsphbox(ptl,mp);
#ifdef DEBUG
	printf("after starfindphsbox rmin =%g %g %g with width=%g for np = %d\n",box.x,box.y,box.z,box.width, mp);
#endif
	Make_Tree_Near(TREE,ptl,mp,box);
#ifdef DEBUG
	printf("after Star_Make_Tree_near\n");
#endif
	/*
	printf("CPU(Make_Tree_test)=%g\n",WALLCLOCK());
	*/
	float mindist;
	mindist = 2.e25;
	int NumNearDen = NUMNEARDEN;
#pragma omp parallel for private(i,j,k) schedule(guided)
	for(i=0;i<np;i++){
		int res;
		int tmpindx[NumNearDen];
		float tmpd2[NumNearDen];
		float tmpmass[NumNearDen];
		float neardist;
		res = Find_Near(p+i,NumNearDen,TREE,ptl,&neardist,tmpindx,tmpd2, tmpmass);
		densph[i] = 0;
		for(j=0;j<res;j++){
			tmpd2[j] = sqrt(tmpd2[j]);
			densph[i] += W4(tmpd2[j],neardist/2.) * ptl[tmpindx[j]].mass;
		}
	}

	for(i=0;i<np;i++){
		ptl[i].type = TYPE_PTL;
		ptl[i].x = bp[i].x;
		ptl[i].y = bp[i].y;
		ptl[i].z = bp[i].z;
		ptl[i].sibling = ptl+i+1;
		ptl[i].mass = bp[i].mass;
	}
	ptl[np-1].sibling = NULL;
	box = findsphbox(ptl,np);
	Make_Tree_Near(TREE,ptl,np,box);
#ifdef DEBUG
	printf("after Make_Tree_near2\n");
#endif
#pragma omp parallel for private(i,j,k) schedule(guided)
	for(i=0;i<np;i++){
		int res;
		int tmpindx[Numnear];
		float neardist;
		k = i*Numnear;
		res = Find_Near(p+i,Numnear,TREE,ptl,&neardist,tmpindx,dist2+k, mass+k);
		if(res != Numnear){
			printf("error occurred %d %d : %d\n", res, Numnear, np);
			exit(99);
		}
		for(j=0;j<res;j++){
			nearindex[k+j] = tmpindx[j];
		}
	}

	Free(mass);Free(dist2); Free(h); 
	Free(TREE);Free(ptl);Free(p);
}
