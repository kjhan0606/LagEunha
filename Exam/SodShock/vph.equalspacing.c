#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include "vph.h"

#define MAX(a,b) ( (a)>(b)?(a):(b) )
#define MIN(a,b) ( (a)<(b)?(a):(b) )

#define NTIME 1000


#define Courant 0.2L

double Gamma = 1.4;
double dtime = 0.2/NTIME;
double targetT = 0.2;
double ntime=0;
typedef struct gasparticle{
	int id;
	double x,den,mass,pressure,vx,ax,energy,ke,ie,Cs;
	double die, dte, dke, dt;
	double ovx;
	double volume,weight;
	double poverrhogam;
	double viscosity;
} gasparticle;

int sortbp(const void *a, const void *b){
	gasparticle *aa, *bb;
	aa = (gasparticle*)a;
	bb = (gasparticle*)b;
	if(aa->x < bb->x) return -1;
	else if(aa->x > bb->x) return 1;
	else return 0;
}


double vph(gasparticle *bp, int np){
	int i,j;
	double area = 1.L;
	double Dtime = 1.E20;
	int isave = -1;

	for(i=1;i<np-1;i++){
		bp[i].den = bp[i].mass/(bp[i+1].x - bp[i-1].x)*2;
		bp[i].ke = 0.5*bp[i].vx*bp[i].vx;
		bp[i].ie = bp[i].pressure/bp[i].den/(Gamma-1);
		bp[i].dt = Dtime;
	}
	for(i=1;i<np-1;i++){
		double ui, ua, ub,dramp, rvel,muij;
		bp[i].ax = 0;
		double die,dte,dke,fx;
		double pij, dS,dr;
		die = dte = dke = fx = 0;
		for(j=-1;j<2;j+=2){ /* j=-1: from left, j=1: from right */
			double vsig,dt;
			double denij = 0.5*(bp[i+j].den +  bp[i].den);

			double Cs = 0.5*(bp[i+j].Cs + bp[i].Cs);

			pij = (bp[i+j].pressure+bp[i].pressure)/2.L;
			dS = j;
			dr = bp[i+j].x - bp[i].x;
			dramp = dr;
			if(dr <0) dramp = -dr;
			ui = 0.5*(bp[i+j].vx - bp[i].vx);
			ub = bp[i].vx;
			ua = 0.5*(bp[i+j].vx + bp[i].vx);

			rvel = ui * dr;
			muij = rvel*eta/(dr*dr + epsilon *epsilon * eta*eta);
			if(muij<0) pij = pij + (-alpha*Cs*muij + beta * muij*muij)*denij;
			die +=  -pij * ui*dS; /* for the internal energy */
			dte +=  -pij * ua*dS; /* for the total energy */
			dke +=  -pij * ub*dS;
			fx += -pij *dS;
			
			double VdotR = (bp[i+j].vx - bp[i].vx)*dr/dramp;
			vsig = (bp[i+j].Cs + bp[i].Cs - MIN(0, VdotR));
			dt = 2*Courant*dramp/vsig;
			bp[i].dt = MIN(bp[i].dt, dt);
			if(dt <  Dtime){
				Dtime = dt;
				isave = i;
			}
		
		}
		bp[i].die = die;
		bp[i].dte = dte;
		bp[i].dke = dke;
		bp[i].ax = fx/bp[i].mass;
	}
	for(i=1;i<np-1;i++){
		bp[i].x += bp[i].vx*Dtime + 0.5*bp[i].ax*Dtime*Dtime;
		bp[i].ovx = bp[i].vx;
		bp[i].vx += bp[i].ax * Dtime;
	}
	qsort(bp,np,sizeof(gasparticle), sortbp);
	for(i=1;i<np-1;i++){
		bp[i].volume = 0.5*(bp[i+1].x - bp[i-1].x);
		bp[i].den = bp[i].mass/bp[i].volume;
        bp[i].ie += bp[i].die/bp[i].volume/bp[i].den*Dtime; 
        bp[i].ke += bp[i].dke/bp[i].volume/bp[i].den*Dtime;
        bp[i].pressure = bp[i].ie*bp[i].den*(Gamma-1);
		bp[i].den = bp[i].mass/(bp[i+1].x - bp[i-1].x)*2;
        bp[i].poverrhogam = bp[i].pressure/pow(bp[i].den,Gamma);
		bp[i].Cs = sqrt(bp[i].pressure/bp[i].den*Gamma);
	}


	if(0)for(i=1;i<np-1;i++){
		double ui, ua, ub,dramp, rvel,muij;
		bp[i].ax = 0;
		double die,dte,dke,fx;
		double pij, dS,dr;
		die = dte = dke = fx = 0;
		for(j=-1;j<2;j+=2){ // j=-1: from left, j=1: from right */
			double vsig,dt;
			double denij = 0.5*(bp[i+j].den +  bp[i].den);

			double Cs = 0.5*(bp[i+j].Cs + bp[i].Cs);

			pij = (bp[i+j].pressure+bp[i].pressure)/2.L;
			dS = j;
			dr = bp[i+j].x - bp[i].x;
			dramp = dr;
			if(dr <0) dramp = -dr;
			ui = 0.5*(bp[i+j].vx - bp[i].vx);
			ub = bp[i].vx;
			ua = 0.5*(bp[i+j].vx + bp[i].vx);

			rvel = ui * dr;
			muij = rvel*eta/(dr*dr + epsilon *epsilon* eta*eta);
			if(muij<0) pij = pij + (-alpha*Cs*muij + beta * muij*muij)*denij;
			die +=  -pij * ui*dS; // for the internal energy */
			dte +=  -pij * ua*dS; // for the total energy */
			dke +=  -pij * ub*dS;
			fx += -pij *dS;
		}
		bp[i].die = die;
		bp[i].dte = dte;
		bp[i].dke = dke;
		bp[i].ax = fx/bp[i].mass;
	}
	/*
	for(i=1;i<np-1;i++){
        bp[i].ie += bp[i].die/bp[i].volume/bp[i].den*Dtime; 
        bp[i].ke += bp[i].dke/bp[i].volume/bp[i].den*Dtime;
        bp[i].pressure = bp[i].ie*bp[i].den*(Gamma-1);
		bp[i].den = bp[i].mass/(bp[i+1].x - bp[i-1].x)*2;
        bp[i].poverrhogam = bp[i].pressure/pow(bp[i].den,Gamma);
		bp[i].Cs = sqrt(bp[i].pressure/bp[i].den*Gamma);
    }
	*/

	return Dtime;
}


int main(int argc, char **argv){

	int np=100;
	gasparticle *bp;

	targetT = atof(argv[1]);

	bp = (gasparticle*)malloc(sizeof(gasparticle)*np);
	double xsep = 1./np;
	int i;
	for(i=0;i<np;i++){
		bp[i].x = xsep*(i+0.5);
		if(bp[i].x<=0.5){
			bp[i].mass = xsep;
			bp[i].volume = xsep;
			bp[i].pressure = 1.L;
			bp[i].den = 1.L;
			bp[i].Cs = sqrt(bp[i].pressure/bp[i].den*Gamma);
		}
		else {
			bp[i].mass = 0.125L*xsep;
			bp[i].volume = 0.125L*xsep;
			bp[i].pressure = 0.1L;
			bp[i].den = 0.125L;
			bp[i].Cs = sqrt(bp[i].pressure/bp[i].den*Gamma);
		}
		bp[i].ie = bp[i].pressure/ ( (Gamma-1)*bp[i].den );
		bp[i].vx = 0;
		bp[i].id = i;
	}
	for(i=1;i<np-1;i++){
		bp[i].volume = 0.5*(bp[i+1].x - bp[i-1].x);
	}
	time_t time1, time2;
	time1 = clock();
	double time=0;
	int jcount = 0;
	do {
		time += vph(bp,np);
		fprintf(stderr,"Time= %g @stepi= %d\n", time, jcount);
		jcount ++;
	} while (time < targetT);
	time2 = clock();
	fprintf(stderr,"time step count = %d :: %15.7lg\n", jcount, (double)(time2-time1)/CLOCKS_PER_SEC);
	for(i=0;i<np;i++){
		printf("%12.7g %12.7g %12.7g %12.7g %5d %12.7g %12.7g %12.7g %12.7g\n",bp[i].x,bp[i].vx,bp[i].den,bp[i].pressure,bp[i].id,bp[i].mass,bp[i].viscosity,bp[i].ax, bp[i].poverrhogam);
	}
}



