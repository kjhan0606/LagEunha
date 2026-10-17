#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include"Complex.h"
#include"nrutil.h"
#include"nr.h"





float omep, omeplam, omepb, omepk;
float x;


float coscon(float y){
	return 1./pow(sqrt(1+x*x*x*pow(y,1.2)),3);
}
float growth(float a, float amax, float Omep){
	float x0 = pow(1./Omep-1, 1./3.);
	float red = amax/a -1;
	float gfac;
	x = pow(1./Omep-1, 1./3.)/(1+red);
	gfac = qsimp(coscon, 0., 1.);
	float growth = (x/x0)*sqrt(x*x*x+1)*gfac;
	return growth;
}


float coscon2(float y){
	omepk = 1 - omep - omeplam;
	float coscon2 = pow(y, 1.5)/pow(omeplam*y*y*y+omepk*y+omep, 1.5);
	return coscon2;
}

float growth2(float a, float amax, float om, float ol){
	float gfac;
	omep = om;
	omeplam = ol;

	omepk = 1-omep-omeplam;

	float aa = a/amax;

	float growth2 = sqrt(omeplam*aa*aa*aa + omepk*aa + omep)/pow(aa, 1.5);

	gfac = qsimp(coscon2, 0., aa);
	growth2 = growth2 * gfac;
	return growth2;
}


float growthgen(float a, float amax, float Omep, float Omeplam, float wlam){
	fcomplex a1, b1, c1, z, a2, b2, c2, hypgeo(fcomplex, fcomplex, fcomplex, fcomplex);
	float xx = (1-Omep)/Omep *pow(a/amax, -3*wlam);

	a1 = Complex(1.5, 0);
	b1 = Complex(-5./6./wlam, 0.);
	c1 = Complex(1.-5./6./wlam, 0.);
	z = Complex(-xx, 0.);
	a2 = Complex(1.5, 0.);
	b2 = Complex(1.-5./6./wlam, 0.);
	c2 = Complex(2.-5./6./wlam, 0.);

	float res = (a/amax)*sqrt(1+xx)* hypgeo(a1,b1,c1,z).r
		+ xx*5*(1+wlam)/(5-6*wlam)*hypgeo(a2,b2,c2,z).i;
	return res;
}
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
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
