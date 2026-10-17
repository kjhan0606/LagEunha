#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<omp.h>
#include "eunha.h"
#define Model4TreeForce
#include "Model4TreeForce.h"
#undef Model4TreeForce


#define RANGE 6
#define MIN(a,b) (a)>(b) ? (b): (a)

double LL;
double pmcf(SimParameters *, double);
double PMFITPARA[8] = {
	1.L/256.L/256.L,
	0.78224E+00L,
	0.37971E-06L,
	-0.60338E+00L,
	-0.34419E-07L,
	-0.39741E+01L,
	-0.10607E+01L,
	-0.38145E+00L
};

void i_force_spline(SimParameters *simpar){
	double r;
	int i,j,k;
	double e,dr;
	double rp2,rp1,rm1,rm2;
	double  fp, fpp,fr;
	ptrdiff_t nx = NX(simpar);
	float rsphere = RSPHERE(simpar);

	LL = nx*nx;
	ran2nran=(double)RANGE/(double)NSPLINE;
	invran2nran = 1./ran2nran;

	dr = ran2nran*0.1;
#ifdef _OPENMP
#pragma omp parallel for private(i,r,rm1,rm2,rp1,rp2,fp,fpp,fr)
#endif
	for(i=1;i<NSPLINE;i++){
		r = (double)(i) * ran2nran;
		rm1 = r - dr; rm2 = r - 2*dr; rp1 = r + dr; rp2 = r + 2*dr;
		fp = (pmcf(simpar, rp1) - pmcf(simpar,rm1))/(dr+dr);
		fpp = (pmcf(simpar,rp2) - 2*pmcf(simpar,r ) + pmcf(simpar, rm2))/(4*dr*dr);
		fr = pmcf(simpar, r);
		forcecorrectdiff(i,0) = -(double)(fr/r);
		forcecorrectdiff(i,1) = -(double)(0.5L*(fp*r-fr)/r/r/r);
		forcecorrectdiff(i,2) = -(double)(0.5L*(fpp*r*r - 3*fp*r + 3*fr)/r/r/r/r/r);
	}
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	for(i=0;i<NSPLINE;i++){
		forcecorrectslope(i,0) = (forcecorrectdiff(i+1,0)-forcecorrectdiff(i,0))/ran2nran;
		forcecorrectslope(i,1) = (forcecorrectdiff(i+1,1)-forcecorrectdiff(i,1))/ran2nran;
		forcecorrectslope(i,2) = (forcecorrectdiff(i+1,2)-forcecorrectdiff(i,2))/ran2nran;
	}
	for(i=0;i<3;i++){
		forcecorrectslope(0,i) = forcecorrectslope(1,i);
		forcecorrectdiff(0,i) = forcecorrectdiff(1,i);
	}
	for(i=0;i<3;i++){
		forcecorrectdiff(NSPLINE-1,i)  =0 ;
		forcecorrectslope(NSPLINE-1,i) =0 ;
	}
}
double pmcf(SimParameters *simpar,double r){
	double fr;
	double tmp1;
	double e;
	/*
	if(r == 0.L) {
		fr = 0.;
	}
	else 
	*/
	{
		tmp1 = cosh(PMFITPARA[1]*r);
		fr = r/pow(r*r+(GRV_EPSILON(simpar))*(GRV_EPSILON(simpar)),1.5L)
			- (
					(1L/(r*r)*tanh(PMFITPARA[1]*r)-PMFITPARA[1]/(tmp1*tmp1*r)
					-2L*PMFITPARA[2]/PMFITPARA[0]*r*exp(PMFITPARA[3]*r*r)
					-2L*PMFITPARA[2]*PMFITPARA[3]/PMFITPARA[0]*r*r*r*exp(PMFITPARA[3]*r*r)
					+PMFITPARA[4]/PMFITPARA[0]*(1L+PMFITPARA[5]*r*r+PMFITPARA[6]*r*r*r*r)*exp(PMFITPARA[7]*r*r))
				);
		fr = fr/LL; 
		/* scaling relation for the PM force which is proportional to 1/L**2 */
		/* One is multiplied to pfact and the other one is to fact2 */
	}
	return fr;
}
void i_potent_spline(SimParameters *simpar){
	double x;
	double xstep;
	int i,j,k;
	double e,dx;
	double g(SimParameters *, double);
	double xp2,xp1,xm1,xm2;
	double  gp, gpp;
	ran2nran=(double)RANGE/(double)NSPLINE;
	invran2nran=1./ran2nran;

	for(i=0;i<3;i++){
		forcecorrectslope(0,i) = forcecorrectdiff(0,i) = 0.;
	}
	dx = ran2nran;

#ifdef _OPENMP
#pragma omp parallel for private(i,x,xm1,xm2,xp1,xp2,gp,gpp)
#endif
	for(i=1;i<NSPLINE;i++){
		x = (double)(i) * ran2nran;
		xm1 = x - dx; xm2 = x - 2*dx; xp1 = x + dx; xp2 = x + 2*dx;
		gp = (g(simpar, xp1) - g(simpar, xm1))/(2*dx);
		gpp = (g(simpar, xp2) - 2*g(simpar, x) + g(simpar, xm2))/(4.*dx*dx);
		forcecorrectdiff(i,0) = g(simpar, x);
		forcecorrectdiff(i,1) = 0.5L*gp/x;
		forcecorrectdiff(i,2) = 0.5L*(gpp/x/x - gp/x/x/x);
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=1;i<NSPLINE;i++){
		forcecorrectslope(i,0) = (forcecorrectdiff(i+1,0)-forcecorrectdiff(i,0))/dx;
		forcecorrectslope(i,1) = (forcecorrectdiff(i+1,1)-forcecorrectdiff(i,1))/dx;
		forcecorrectslope(i,2) = (forcecorrectdiff(i+1,2)-forcecorrectdiff(i,2))/dx;
	}
	return;
}
double g(SimParameters *simpar,double x){
	double e,g;
	g = -1.L/sqrt(x*x+(GRV_EPSILON(simpar))*(GRV_EPSILON(simpar)));
	return g;
}
