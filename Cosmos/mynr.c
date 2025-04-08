#include <math.h>
#include "eunha.h"
#include "Complex.h"
#include "nrutil.h"
/*
*/

float simpartrapzd(float (*)(SimParameters *, float), float , float , int , SimParameters *);

/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
#include <math.h>
#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5
float simparqromb(float (*func)(SimParameters *, float), float a, float b, SimParameters *simpar)
{
	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	void nrerror(char error_text[]);
	float ss,dss;
	float s[JMAXP+1],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=simpartrapzd(func,a,b,j, simpar);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */

#define EPS 1.0e-4
#define JMAX 20
float simparqsimp(float (*func)(SimParameters *, float), float a, float b, SimParameters *simpar)
{
	void nrerror(char error_text[]);
	int j;
	float s,st,ost,os;

	if(a==b) return 0.0;

	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=simpartrapzd(func,a,b,j, simpar);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) <= EPS*fabs(os)) return s;
		os=s;
		ost=st;
	}
	return 0.0;
}
#undef JMAX
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
#define FUNC(simpar,x) ((*func)(simpar,x))

float simpartrapzd(float (*func)(SimParameters *, float), float a, float b, int n, SimParameters *simpar)
{
	float x,tnm,sum,del;
	static float s;
	int it,j;

	if(a==b) return 0.0;
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
#include <math.h>
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5

float simparqromo(float (*func)(SimParameters *, float), float a, float b,
	float (*choose)(float(*)(SimParameters *, float), SimParameters *, float, float, int), SimParameters *simpar)
{
	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	void nrerror(char error_text[]);
	int j;
	float ss,dss,h[JMAXP+1],s[JMAXP+1];
	if(a==b) return 0.0;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,simpar, a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=h[j]/9.0;
	}
	nrerror("Too many steps in routing qromo");
	return 0.0;
}
#undef JMAX
#undef JMAXP
#undef K
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
#include <math.h>
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
#define FUNC(a,x) ((*funk)(a, 1.0/(x))/((x)*(x)))

float midinf(float (*funk)(SimParameters *, float), SimParameters *simpar, float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,b,a;
	static float s;
	int it,j;

	b=1.0/aa;
	a=1.0/bb;
	if (n == 1) {
		return (s=(b-a)*FUNC(simpar, 0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(simpar, x);
			x += ddel;
			sum += FUNC(simpar, x);
			x += del;
		}
		return (s=(s+(b-a)*sum/tnm)/3.0);
	}
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
#include <math.h>
#define FUNC(a,x) ((*funk)(a,-log(x))/(x))

float midexp(float (*funk)(SimParameters *, float), SimParameters *simpar, float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,a,b;
	static float s;
	int it,j;

	b=exp(-aa);
	a=0.0;
	if (n == 1) {
		return (s=(b-a)*FUNC(simpar, 0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(simpar,x);
			x += ddel;
			sum += FUNC(simpar,x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
#undef EPS
#include <math.h>
#define FUNC(a, x) (2.0*(x)*(*funk)(a, aa+(x)*(x)))

float midsql(float (*funk)(SimParameters *,float), SimParameters *simpar, float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,a,b;
	static float s;
	int it,j;

	b=sqrt(bb-aa);
	a=0.0;
	if (n == 1) {
		return (s=(b-a)*FUNC(simpar, 0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(simpar, x);
			x += ddel;
			sum += FUNC(simpar, x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
#include <math.h>
#define FUNC(a,x) (2.0*(x)*(*funk)(a,bb-(x)*(x)))

float midsqu(float (*funk)(SimParameters *, float), SimParameters *simpar, float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,a,b;
	static float s;
	int it,j;

	b=sqrt(bb-aa);
	a=0.0;
	if (n == 1) {
		return (s=(b-a)*FUNC(simpar, 0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(simpar, x);
			x += ddel;
			sum += FUNC(simpar, x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */


/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file complex.c.  Do not confuse this file with the same-named
   file complex.c that is supplied in the same subdirectory or archive
   as the header file complex.h.  *That* file contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#include <math.h>

fcomplex Cadd(fcomplex a, fcomplex b)
{
	fcomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

fcomplex Csub(fcomplex a, fcomplex b)
{
	fcomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}


fcomplex Cmul(fcomplex a, fcomplex b)
{
	fcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

fcomplex Complex(float re, float im)
{
	fcomplex c;
	c.r=re;
	c.i=im;
	return c;
}

fcomplex Conjg(fcomplex z)
{
	fcomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}

fcomplex Cdiv(fcomplex a, fcomplex b)
{
	fcomplex c;
	float r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

float Cabs(fcomplex z)
{
	float x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

fcomplex Csqrt(fcomplex z)
{
	fcomplex c;
	float x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

fcomplex RCmul(float x, fcomplex a)
{
	fcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}
