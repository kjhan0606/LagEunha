#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include "eunha.h"
#include "Complex.h"
#include "nrutil.h"
#include "nr.h"
#include "flrw.h"



float midinf(float (*)(SimParameters *, float), SimParameters *, float , float , int );
float midsql(float (*)(SimParameters *, float), SimParameters *, float , float , int );
float midsqu(float (*)(SimParameters *, float), SimParameters *, float , float , int );
float midexp(float (*)(SimParameters *, float), SimParameters *, float , float , int );

/* Quintessence model is based on Sefusatti & Vernizzi 2011, arXiv:1101.1026v3 */
/* Note: read just below equation(29) . The Growthfactor for quintessence model is give in 
 * integral form (Eq. 29)*/

float simparqsimp(float (*func)(SimParameters *, float), float, float, SimParameters *);
float simparqromb(float (*func)(SimParameters *, float), float, float, SimParameters *);
float simparqromo(float (*)(SimParameters *, float), float , float ,
		    float (*)(float(*)(SimParameters *, float), SimParameters *, float, float, int), SimParameters *);

float wlambda(SimParameters *simpar, float red){
	float result = WLAM0(simpar) + WLAM1(simpar) * red/(1+red);
	return result;
}
float DEfunc(SimParameters *simpar, float tau){
	float red;
	if(tau==0) red = 1.e20;
	else red = 1/tau - 1;
	float result = 3*(1+wlambda(simpar, red)) *(1+red);
	return result;
}

float gzDE(SimParameters *simpar, float anow){
	float red,integral;
	if(WLAM0(simpar)==-1 && WLAM1(simpar)==0) return 1;
	else if(WLAM1(simpar)==0) {
		return  pow(AMAX(simpar)/anow, 3.*(1+WLAM0(simpar)));
	}
	else if(anow == AMAX(simpar)) {
		integral = 0;
	}
	else if(anow==0) {
		integral = simparqsimp(DEfunc, 0, 1, simpar);
	}
	else {
		float tau = anow/AMAX(simpar);
		if(tau <0.8) integral = simparqsimp(DEfunc, tau, 1, simpar);
		else integral = simparqromb(DEfunc, tau, 1, simpar);
	}
	integral = exp(integral);
	return integral;
}

float HofEz(SimParameters *simpar, float anow){
	float red = AMAX(simpar)/anow -1;
	float zp1 = red + 1;
	return sqrt(zp1*zp1*zp1*OMEP(simpar)+zp1*zp1*(1-OMEP(simpar)-OMEPLAM(simpar))+OMEPLAM(simpar) * gzDE(simpar, anow) );
}

float ap(SimParameters *simpar, float anow){
	float red = AMAX(simpar)/anow-1;
	float hofez = HofEz(simpar, anow);
    float result =  sqrt(8.*M_PI/3./(AMAX(simpar)*OMEP(simpar))/((red+1)*(red+1)) * hofez*hofez);
	return result;
}
float apold(SimParameters *simpar, float anow){
	float red = AMAX(simpar)/anow-1;
    return sqrt(8.*M_PI/3.*(1/anow+1/OMEI(simpar)-1+OMEPLAM(simpar)/OMEP(simpar)*(pow(anow, -1-3*WLAM0(simpar))-1)*pow(AMAX(simpar),3*WLAM0(simpar))));
}


float Omega_matter(SimParameters *simpar, float anow){
	if(anow ==0){
		return 1;
	}
	else {
		float red = AMAX(simpar)/anow-1;
		float zp1 = red + 1;
		float hofez = HofEz(simpar, anow);
		return zp1*zp1*zp1*OMEP(simpar)/(hofez*hofez);
	}
}
float Omega_lambda(SimParameters *simpar, float anow){
	if(anow==0){
		return 0;
	}
	float red = AMAX(simpar)/anow-1;
	float hofez = HofEz(simpar, anow);
	return gzDE(simpar,anow)*OMEPLAM(simpar)/(hofez*hofez);
}

float Cw(SimParameters *simpar, float anow){
	float red = AMAX(simpar)/anow-1;
	return 1. + (1+wlambda(simpar, red))*Omega_lambda(simpar, anow)/Omega_matter(simpar, anow);
}
float CwOveraH3(SimParameters *simpar, float anow){
	if(anow ==0) return 0;
	else {
		return Cw(simpar,anow)/pow(anow*HofEz(simpar,anow),3);
	}
}
float Dall(SimParameters *simpar, float anow){
	float amax = AMAX(simpar);
	float red = amax/anow-1;
	float a1 ;
	if(anow==0) return 0;
	else {
		a1 = simparqromo(CwOveraH3, 0, anow, midsqu,simpar); 
		return 5./2.*OMEP(simpar)*HofEz(simpar, anow)*a1*amax*amax;
	}
}
float Dmslope(SimParameters *simpar, float anow){
	float amax = AMAX(simpar);
	float a1,a2,a3;
	a1 = Dall(simpar, anow);
	float red = amax/anow-1;
	a2 = Omega_matter(simpar, anow);
	return (2.5-1.5*a1*amax/anow) *a2;
}
float growthall(SimParameters *simpar, float anow){
	printf("growthall with anow = %g\n",anow);
	/*
	if(anow < 40)  return simparqromo(Dmslope, 0, anow, midsqu,simpar);
	else return simparqsimp(Dmslope, 0, anow, simpar);
	*/
	return simparqromo(Dmslope, 0, anow, midsqu,simpar)/AMAX(simpar);
}

float coscon(SimParameters *simpar, float y){
	float x = COSCONX(simpar);
	return 1./pow(sqrt(1+x*x*x*pow(y,1.2)),3);
}
float coscon2(SimParameters *simpar, float y){
	float omeplam = OMEPLAM(simpar);
	float omep = OMEP(simpar);
	float omepk = 1 - omeplam - omep;
	float coscon2 = pow(y, 1.5)/pow(omeplam*y*y*y+omepk*y+omep, 1.5);
	return coscon2;
}

float growth(SimParameters *simpar, float anow){
	float amax = AMAX(simpar);
	float Omep = OMEP(simpar);
	float x0 = pow(1./Omep-1, 1./3.); 
	float red = amax/anow -1; 
	float gfac; 
	float x = COSCONX(simpar) = pow(1./Omep-1, 1./3.)/(1+red); 
	gfac = simparqsimp(coscon, 0., 1.,simpar); 
	float growth = (x/x0)*sqrt(x*x*x+1)*gfac; 
	return growth;

}

float growth2(SimParameters *simpar, float anow){
	float gfac; 
	float omep = OMEP(simpar); 
	float omeplam = OMEPLAM(simpar); 
	float omepk = 1-omep-omeplam; 
	float aa = anow/AMAX(simpar); 
	float growth2 = sqrt(omeplam*aa*aa*aa + omepk*aa + omep)/pow(aa, 1.5); 
	gfac = simparqsimp(coscon2, 0., aa,simpar); 
	growth2 = growth2 * gfac; 
	return growth2;
}


float DplusAllQuint0(SimParameters *simpar, float anow){ 
		float Omep = OMEP(simpar);
		float Omeplam = OMEPLAM(simpar);
		float amax = AMAX(simpar);
		float wlam = WLAM0(simpar);
	    fcomplex a1, b1, c1, z, a2, b2, c2, hypgeo(fcomplex, fcomplex, fcomplex, fcomplex); 
		float xx = (1-Omep)/Omep *pow(anow/amax, -3*wlam); 
		a1 = Complex(1.5, 0.);
		b1 = Complex(-5./6./wlam, 0.);
		c1 = Complex(1.-5./6./wlam, 0.);
		z = Complex(-xx, 0.);
		a2 = Complex(1.5, 0.);
		b2 = Complex(1.-5./6./wlam, 0.);
		c2 = Complex(2.-5./6./wlam, 0.); 
		float geo1 = hypgeo(a1,b1,c1,z).r;
		float geo2 = hypgeo(a2,b2,c2,z).r;
		float res = (anow/amax)*sqrt(1+xx)* (geo1 + xx*5*(1+wlam)/(5-6*wlam)*geo2);
		return res;
}

float XgrowthM(SimParameters *simpar, float anow){
	return (2.5 - 1.5*DplusAllQuint0(simpar, anow)*AMAX(simpar)/anow)*Omega_matter(simpar,anow)/AMAX(simpar);
}
float growthgen(SimParameters *simpar, float anow){ 
	return simparqromo(XgrowthM, 0, anow, midsql,simpar);
}
