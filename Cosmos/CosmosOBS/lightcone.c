/*
 * �� �ڵ�� ������ X0,Y0,Z0 �� ��ġ�� �����ڰ� �ùķ��̼� �ڽ��� �������
 * �ùķ��̼� ���ڵ��� ����������  ������ ������ ��Ÿ����.
 * ���� ����� comoving space ������ �ùķ��̼ǿ��� ����� �Ÿ������� ���´�.
 * ����� theta, phi  �� ���� 1��и鿡���� ����� �����ϴ�.
 * ���� 2,3,4 ��и鿡���� �����ϰ� �ҷ���, min,max �� ������ ������
 * �ٲپ�� �Ѵ�. �׸��� ������ �κп� �ش�� ���ڸ� ã�� �� phi ���� 
 * ����ϴµ�,
 * ���⿡���� �׳� 1,4��и��̶�� ������ �ߴ�. 
 * wflag �� ���� 2�������� ǥ���ؾ� �Ѵ�.
 * 09/08/2002 ������ 
 */
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <mpi.h>
#include "eunha.h"
#include "lightcone.h"
#define MIN(A,B) ((A) <(B) ? (A):(B))
#define PI 3.14159265354979L
double Theta,Phi,DTheta,DPhi;
char surveyfilename[100];
int *freework;
typedef struct Spherical{
	double r,theta,phi;
} Spherical;
unsigned int setbits(unsigned int x,int p) {
    unsigned int  b;
	b = x|(~((~0)<<1) << (p-1));
	return b;
}
InBox *box;
double r1, r2;
double theta1,phi1,theta2,phi2;
float omega0,H,lambda0;
float func1(float red){
    return 2997.92458L/sqrt((1.+red)*(1.+red)*(1.+red)*omega0
             +(1.+red)*(1.+red)*(1.-omega0-lambda0) +lambda0);
}
float qsimp(float (*)(float),float,float);
float trapzd(float (*)(float), float , float , int );
void a2comovingpixel(double *R, double *dR1, double *dR2,int nx, 
		float amax,float a, float astep1,float astep2,
		float omep,float omeplam,float size, float h){
	float initz = 0;
	float shellwidth1,shellwidth2;
	float redshift = amax/a - 1.L;
	shellwidth1 = amax/(a+astep2*0.5L)-1.L;
	shellwidth2 = amax/(a-astep1*0.5L)-1.L;
	if(omep == 1.0) {
		*R = 2997.92458L*2.L * (1.L-1.L/sqrt(1.L+redshift));
		*dR1 = 2997.92458L*2.L*(1.L-1.L/sqrt(1.L+shellwidth1));
		*dR2 = 2997.92458L*2.L*(1.L-1.L/sqrt(1.L+shellwidth2));
	}
	else {
		omega0 = omep;
		lambda0 = omeplam;
		H = h;
		*R = qsimp(func1,initz,redshift);
		*dR1 = qsimp(func1,initz,shellwidth1);
		*dR2 = qsimp(func1,initz,shellwidth2);
	}
	*R = *R/size*nx;
	*dR1 = *dR1/size*nx;
	*dR2 = *dR2/size*nx;
	*dR1 = *R - *dR1;
	*dR2 = *dR2 - *R;
}

int ssorttype(const void *a,const void *b){
	slcparticletype *aa,*bb;
	aa = (slcparticletype *)a;
	bb = (slcparticletype *)b;
	if(aa->type < bb->type) return -1;
	else if(aa->type > bb->type) return +1;
	else return 0;
}

#include <math.h>
#define EPS 1.0e-6
#define JMAX 20

float qsimp(float (*func)(float), float a, float b)
{
	float trapzd(float (*func)(float), float a, float b, int n);
	int j;
	float s,st,ost,os;

	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=trapzd(func,a,b,j);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		os=s;
		ost=st;
	}
	return 0.0;
}
#undef EPS
#undef JMAX
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
#define FUNC(x) ((*func)(x))

float trapzd(float (*func)(float), float a, float b, int n)
{
	float x,tnm,sum,del;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
