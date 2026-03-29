#include<stdio.h>
#include<math.h>
#include "ramses.h"
#include "tree.h"
#include "force_spline.h"
/*
#define NSPLINE 100000
#define RANGE2 16
*/
/*
#define RANGE2 1000
*/
#define RANGE 20
float diff[NSPLINE][3],slope[NSPLINE][3];
float ran2nran;
double acoef[8] = {
	1./256./256.,
	0.78224E+00L,
	0.37971E-06L,
	-0.60338E+00L,
	-0.34419E-07L,
	-0.39741E+01L,
	-0.10607E+01L,
	-0.38145E+00L
};
void i_force_spline(){
	double x;
	double xstep;
	int i,j,k;
	double e;
	ran2nran=(float)RANGE/(float)NSPLINE;

	for(i=0;i<3;i++){
		slope[0][i] = diff[0][i] = 0.;
	}

	for(i=1;i<NSPLINE;i++){
		double xp2,xp1,x,xm1,xm2;
		double  dx, fp, fpp;
		/*
		x = (double)(i) * ran2nran;
		x = sqrt(x);
		dx = sqrt((i+1)*ran2nran) - sqrt(i*ran2nran);
		*/
		x = (double)(i) * ran2nran;
		dx = ran2nran;
		xm1 = x - dx; xm2 = x - 2*dx; xp1 = x + dx; xp2 = x + 2*dx;
		fp = (f(xp1) - f(xm1))/(2.*dx);
		fpp = (f(xp2) - 2.*f(x) + f(xm2))/(4.*dx*dx);
		diff[i][0] = f(x)/x;
		diff[i][1] = 0.5L*(fp/x-f(x)/x/x)/x;
		diff[i][2] = 0.5L*(fpp/x - fp/x/x + 3.L*f(x)/x/x/x
				-2.L*fp/x/x)/x/x;
	}
	for(i=1;i<NSPLINE;i++){
		/*
		x = (double)(i) * (double)RANGE2 / (double)NSPLINE;
		xstep = (double)(i+1)*(double)RANGE2/(double)NSPLINE - x;
		*/
		xstep = ran2nran;
		slope[i][0] = (diff[i+1][0]-diff[i][0])/xstep;
		slope[i][1] = (diff[i+1][1]-diff[i][1])/xstep;
		slope[i][2] = (diff[i+1][2]-diff[i][2])/xstep;
	}
	return;
}
/* corrected force for pm force */
double f(double x){
	double fx;
	double tanha2x,expdxx;
	double e;
	e = EPSILON;
	tanha2x = tanh(acoef[1]*x);
	expdxx = exp(acoef[3]*x*x);
	fx = -acoef[0]*( 1.L/pow(x*x+e*e,1.5)  
		- 1.L/x/x*(acoef[1]*x*(1.L-tanha2x*tanha2x)-tanha2x)/x
		- 2.L*acoef[2]/acoef[0]*expdxx*(1.L+acoef[3]*x*x)
		- 1.L/x/acoef[0]*(acoef[4]+acoef[5]*x*x+acoef[6]*x*x*x*x)*exp(acoef[7]*x*x));
	return fx;
}
void i_potent_spline(){
	double x;
	double xstep;
	int i,j,k;
	double e;
	ran2nran=(float)RANGE/(float)NSPLINE;

	for(i=0;i<3;i++){
		slope[0][i] = diff[0][i] = 0.;
	}

	for(i=1;i<NSPLINE;i++){
		double xp2,xp1,x,xm1,xm2;
		double  dx, gp, gpp;
		/*
		x = (double)(i) * ran2nran;
		x = sqrt(x);
		dx = sqrt((i+1)*ran2nran) - sqrt(i*ran2nran);
		*/
		x = (double)(i)*ran2nran;
		dx = ran2nran;
		xm1 = x - dx; xm2 = x - 2*dx; xp1 = x + dx; xp2 = x + 2*dx;
		gp = (g(xp1) - g(xm1))/(2.*dx);
		gpp = (g(xp2) - 2.*g(x) + g(xm2))/(4.*dx*dx);
		diff[i][0] = g(x);
		diff[i][1] = 0.5L*gp/x;
		diff[i][2] = 0.5L*(gpp/x/x - gp/x/x/x);
	}
	for(i=1;i<NSPLINE;i++){
		/*
		x = (double)(i) * (double)RANGE2 / (double)NSPLINE;
		xstep = (double)(i+1)*(double)RANGE2/(double)NSPLINE - x;
		*/
		xstep = ran2nran;
		slope[i][0] = (diff[i+1][0]-diff[i][0])/xstep;
		slope[i][1] = (diff[i+1][1]-diff[i][1])/xstep;
		slope[i][2] = (diff[i+1][2]-diff[i][2])/xstep;
	}
	for(i=0;i<3;i++){
		slope[0][i] = slope[1][i];
		diff[0][i] = diff[1][i];
		slope[NSPLINE-1][i] = slope[NSPLINE-2][i];
		diff[NSPLINE-1][i] = diff[NSPLINE-2][i];
	}
	/*
	diff[NSPLINE-1][0] = 0;
	diff[NSPLINE-1][1] = 0;
	diff[NSPLINE-1][2] = 0;
	slope[NSPLINE-1][0] = 0;
	slope[NSPLINE-1][1] = 0;
	slope[NSPLINE-1][2] = 0;
	*/
	return;
}
double g(double x){
	double e,g;
	g = -1.L/sqrt(x*x+EPSILON*EPSILON);
	return g;
}
