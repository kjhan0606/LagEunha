#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include "eunha.h"
#include "hydroBasicCell.h"
#include "nnost.h"
#include "sph.h"

//jhshin1
#define _MAIN_MU
#include "mu.h"
#undef _MAIN_MU
//jhshin2

//jhshin1
//#define nDelta 128
//#define nT 128
//#define nTMU 128
//jhshin2

#define tempmin 1000.
#define tempmax 10000.
#define log10rhomin -3.
#define log10rhomax 7.

const double C=2.99792458E10L;
const double k_B=1.380658E-16L;
const double h_P=6.6260755E-27L;
const double m_e=9.1093897E-28L;
const double m_H=1.673575E-24L;
const double NOT4=3.9715E0L; /* mass ratio (He/H) atom */
const double sigma=6.6524616E-29L;
const double a=7.565914E-16L;
const double G=6.6742E-8L;
const double eV=1.6021765314E-12L;
const int g1=1, g2 = 4, g3 = 1;
double Tstep,log10tempmin,log10tmumin;
static double Tmumin,Tmumax,Tmustep;

//jhshin1
//double aMU[nT][nDelta];
//double aTMU[nT][nDelta];
//double aT[nTMU][nDelta];
//jhshin2

const double bigH = 100.E3L/(1.E6L*3.0856775807E16L);


double four(double a0,double a1, double a2,double a3,double a4, double x){
	double res;
	res = a4*x*x*x*x+a3*x*x*x+a2*x*x+a1*x+a0;
	return res;
}

#define M 1000
double finderatio(double a0,double a1,double a2,double a3,double a4) {
	double a,y;
	double x = 2.L;
	double b = four(a0,a1,a2,a3,a4,x);
	double utrial = log10(2.L);
	double btrial = log10(1.e-20L);
	double xstep = (utrial-btrial)/(double)M;
	int i;
	for(i=M-1;i>=0;i--){
		x = pow(10.L,xstep*i+btrial);
		a =  four(a0,a1,a2,a3,a4,x);
		if(a*b <= 0.L) break;
	}
	if(i ==0) x=0.L;
	utrial = (xstep*(i+1)+btrial);
	btrial = (xstep*(i)+btrial);

	while(fabs(utrial-btrial) > 1.e-10L){
		x = pow(10.L,0.5L*(utrial+btrial));
		a = four(a0,a1,a2,a3,a4,x);
		y = a*b;
		if(y > 0.L) utrial = log10(x);
		else if ( y < 0.) btrial = log10(x);
		else return x;
	}
	x = pow(10.L,0.5L*(utrial+btrial));
	return x;
}

void InitializeMu(SimParameters *simpar, double red){
	double tol = 1.E-5L;
	double HO = HUBBLE(simpar) * bigH;
	double Yp = GAS_YP(simpar);
	double mu_H = 1.L/(1.L-Yp);
	double mu_T = NOT4/(NOT4*(1.L-Yp));
	double fHe = Yp/(NOT4*(1.L-Yp));
	double nNow;
	if(BGEXPAND(simpar) == 'Y') {
		nNow = 3*HO*HO*OMEP(simpar)/(8*M_PI*G*mu_H*m_H);
	}
	else {
		nNow = TOTMASS(simpar)*onesolarmass/pow(BOXSIZE(simpar)*Mpccgs,3.L);
	}
	double zp1 = red+1;
	double n = nNow*zp1*zp1*zp1;
	double fnu = 21.L/8.L*pow(4.L/11.L,4.L/3.L);
	Tstep = (log10(tempmax)-log10(tempmin))/nT;
	log10tempmin = log10(tempmin);

		Tmumin = tempmin/mumax;
		Tmumax = tempmax/mumin;
		Tmustep = (log10(Tmumax)-log10(Tmumin))/nTMU;
		log10tmumin = log10(Tmumin);

	double log_drhostep = (log10rhomax-log10rhomin)/nDelta;
	double chi1 = 13.54*eV;
	double chi2 = 24.48*eV;
	double chi3 = 54.17*eV;
	double CR = 2*M_PI*(m_e/h_P)*(k_B/h_P);
	int i,j,k;
	for(j = 0;j<nDelta;j++){
		double density = log10rhomin + log_drhostep*j;
		density = pow(10.L,density);
		double n = nNow*zp1*zp1*zp1*density;
		double nH = n/(1.L+fHe);
		double nHe = n-nH;
		for(i=0;i<nT;i++){
			double Temp = log10(tempmin) + Tstep*i;
			Temp = pow(10.L,Temp);
			double f1 = pow(CR*Temp,1.5L)*g1*exp(-chi1/k_B/Temp);
			double f2 = pow(CR*Temp,1.5L)*g2*exp(-chi2/k_B/Temp);
			double f3 = pow(CR*Temp,1.5L)*g3*exp(-chi3/k_B/Temp);
			double p1 = f1/(nH+nHe);
			double p2 = f2/(nH+nHe);
			double p3 = f3/(nH+nHe);
			double pnH = nH/(nH+nHe);
			double pnHe = nHe/(nH+nHe);
			double a4 = 1.L;
			double a3 = (p1+p2);
			double a2 = (p1*p2+p2*p3-p2*pnHe-p1*pnH);
			double a1 = (p1*p2*p3-p1*p2*(pnH+pnHe)-2.L*p2*p3*pnHe);
			double a0 =  -(p1*p2*p3*(pnH+pnHe+pnHe));
			double x;
			x =  finderatio(a0,a1,a2,a3,a4);
			double ne = x*(nH+nHe);
			double nH0 = nH/(1.L+f1/ne);
			double nHp = nH0*f1/ne;
			double nHe0 = nHe/(1.L+f2/ne+f2*f3/ne/ne);
			double nHep = nHe0/ne*f2;
			double nHepp = nHep*f3/ne;
			double beta = (nHp+nHep+nHepp)/(nH+nHe);
			double x_H = nHp/nH;
			double x_He = (nHepp+nHep)/nHe;
			double x_He2 = nHepp/nHe;
			double meanwmol = (NOT4*fHe+1.L)/(1.L+x_H+fHe*(1.L+x_He+x_He2));
			double alpha = 1.14E-19L * 4.309L*pow(Temp/1.L,0.6166L)/
				(1.L+0.6703*pow(Temp/1.E4L,0.53L)) * 1.E6L;

			aMU(i,j) = meanwmol;
			aTMU(i,j) = Temp/meanwmol;
		}
		for(i=0;i<nTMU;i++){
			double Tmu = log10(Tmumin) + Tmustep*i;
			Tmu = pow(10.L,Tmu);
			if(Tmu<aTMU(0,j)) aT(i,j) = Tmu*mumax;
			else if(Tmu>=aTMU(nT-1,j)) aT(i,j) = Tmu*mumin;
			else {
				for(k = 0;k<nT-1;k++){
					if(aTMU(k,j)<=Tmu && aTMU(k+1,j)>Tmu){
						aT(i,j) = pow(10.L,log10(tempmin) + Tstep*k) + 
							(pow(10.L,log10(tempmin) + Tstep*(k+1)) -pow(10.L,log10(tempmin) + Tstep*k))/
							(aTMU(k+1,j) -aTMU(k,j)) * (Tmu-aTMU(k,j));
					}
				}
			}
		}
	}
}

float getmu(SimParameters *simpar, float temp, float rho){
	float result;
	float Yp = GAS_YP(simpar);
	float fHe = Yp/(NOT4*(1.L-Yp));
	
//jhshin1
//	return 1.0;
//jhshin2
 
	if(temp < tempmin) {
		result = (NOT4*fHe+1.L)/(1.L+fHe);
		return result;
	}
	else if( temp > tempmax ) {
		result = (NOT4*fHe+1.L)/(2.L+fHe*3);
		return result;
	}
	else {
		int i,j;
		j = (log10(rho) - log10rhomin)/(log10rhomax-log10rhomin)*nDelta;
		j = max(0,min(j,nDelta-1));
		i = (log10(temp) - log10tempmin)/Tstep;
		i = max(0,min(i,nT-1));
		return aMU(i,j);
	}
	return 0;
}

/*
float GetMu(treesphparticletype *sph){
	float Temp;
	return getmu(sph->Temp,sph->rho);
}
*/
float oldgetT(float Tmu, float rho){
	float result;
	if(Tmu < Tmumin) {
		result = Tmu*mumax;
		return result;
	}
	else if( Tmu > Tmumax ) {
		result = Tmu*mumin;
		return result;
	}
	else {
		int i,j;
		j = (log10(rho) - log10rhomin)/(log10rhomax-log10rhomin)*nDelta;
		j = max(0,min(j,nDelta-1));
		i = (log10(Tmu) - log10tmumin)/Tmustep;
		if(i<0) return Tmu*mumax;
		else if(i>=nTMU) return Tmu*mumin;
		else return aT(i,j);
	}
	return 0;
}
float getT(SimParameters *simpar, treesphparticletype *sph){
	float Tmu = (GAS_GAMMA(simpar)-1)*sph->Entropy*pow(sph->rho,GAS_GAMMA(simpar)-1)*mHoverkB;
	float rho = sph->rho;
	float result;
	if(Tmu < Tmumin) {
		result = Tmu*mumax;
		return result;
	}
	else if( Tmu > Tmumax ) {
		result = Tmu*mumin;
		return result;
	}
	else {
		int i,j;
		j = (log10(rho) - log10rhomin)/(log10rhomax-log10rhomin)*nDelta;
		j = max(0,min(j,nDelta-1));
		i = (log10(Tmu) - log10tmumin)/Tmustep;
		if(i<0) return Tmu*mumax;
		else if(i>=nTMU) return Tmu*mumin;
		else return aT(i,j);
	}
	return 0;
}



float GetT(SimParameters *simpar, treesphparticletype *sph){
	float Tmu = (GAS_GAMMA(simpar)-1)*sph->Entropy*pow(sph->rho,GAS_GAMMA(simpar)-1)*mHoverkB;
//jhshin1
//	return Tmu;
//jhshin2
	float temp = oldgetT(Tmu,sph->rho);
	return temp;
}
float GetTempGivenMu(SimParameters *simpar, treesphparticletype *sph){
	return (GAS_GAMMA(simpar)-1)*sph->Entropy*pow(sph->rho,GAS_GAMMA(simpar)-1)*mHoverkB*sph->mu;
	/*
	return sph->temp;
	*/

}
