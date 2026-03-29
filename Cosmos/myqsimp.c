#include <math.h>
#include "eunha.h"
#define EPS 1.0e-6
#define JMAX 20

float simpartrapzd(float (*)(SimParameters *, float), float , float , int , SimParameters *);

float simparqsimp(float (*func)(SimParameters *, float), float a, float b, SimParameters *simpar)
{
	float trapzd(float (*func)(float), float a, float b, int n, SimParameters);
	void nrerror(char error_text[]);
	int j;
	float s,st,ost,os;

	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=simpartrapzd(func,a,b,j, simpar);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		os=s;
		ost=st;
	}
	nrerror("Too many steps in routine qsimp");
	return 0.0;
}
#undef EPS
#undef JMAX
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
#define FUNC(simpar,x) ((*func)(simpar,x))

float simpartrapzd(float (*func)(SimParameters *, float), float a, float b, int n, SimParameters *simpar)
{
	float x,tnm,sum,del;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(simpar,a)+FUNC(simpar,b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(simpar,x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
