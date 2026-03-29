#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include "eunha.h"
#include "flrw.h"
#define EPS 1.0e-6
#define JMAX 20

/*
float simpartrapzd(float (*)(SimParameters *, float), float , float , int , SimParameters *);
float simparqsimp(float (*func)(SimParameters *, float), float a, float b, SimParameters *simpar)
{
	float trapzd(float (*func)(float), float a, float b, int n, SimParameters);
	void nrerror(char error_text[]);
	int j;
	float s,st,ost,os;

	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		printf("simp %g %g\n", a,b);
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

void nrerror(char error_text[]);
*/


float ap(SimParameters *, float);
float apold(SimParameters *, float);
float HofEz(SimParameters *, float);
float growthall(SimParameters *, float);
float gzDE(SimParameters *, float);
float CwOveraH3(SimParameters *, float);
float Omega_lambda(SimParameters *, float);
float Dmslope(SimParameters *, float);
float DplusAllQuint0(SimParameters *, float );

int main(int argc, char **argv){
	SimParameters simpar, *s;
	float amax;

	s = &simpar;
	AMAX(s) = amax = 101;
	OMEP(s) = 0.26;
	float zi = AMAX(s) - 1;
	OMEPLAM(s) = 0.74;
	WLAM0(s) = -0.5;
	WLAM1(s) = 0.;
	OMEI(s) = OMEP(s) * pow(1+zi,3)/HofEz(s,zi)/HofEz(s,zi);
	float anow = 100;

	anow = 1;
	float a = DplusAllQuint0(s, anow);
	float b = DplusAllQuint0(s, amax);
	printf("growthgen %g  : %g %g : %g\n",growthgen(s,anow)/growthgen(s,amax), a,b,a/b);

	/*
	int i;
	for(i=0;i<amax-1;i++){
		anow = i+1;
		float red = AMAX(s)/anow - 1;
		printf(" %g %g : %g :%g %g\n", anow, red, Dmslope(s, anow), growthgen(s,anow)/growthgen(s,amax), 
				DplusAllQuint0(s, anow)/DplusAllQuint0(s,amax));
	}
	*/
}
