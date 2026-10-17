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


#define Courant 0.4L

double xmin,xmax;

double Gamma = 1.4L;
double dtime = 0.2/NTIME;
double targetT = 0.2;
double ntime=0;
double meand;
// RK4 coefficients
typedef struct RK4{
    double k1x,k1vx;
    double k2x,k2vx;
    double k3x,k3vx;
    double k4x,k4vx;
} RK4;
typedef struct gasparticle{
	int id;
	double x,den,mass,pressure,vx,ax,energy,ke,ie,Cs;
	double oden,ow2, odrad1,odrad2;
	double die, dte, dke, dt;
	double ovx;
	double volume,weight;
	double poverrhogam;
	double viscosity;
	double w2,  w2ceil, drad1, drad2;
	RK4 rk4;
} gasparticle;

int sortbp(const void *a, const void *b){
	gasparticle *aa, *bb;
	aa = (gasparticle*)a;
	bb = (gasparticle*)b;
	if(aa->x < bb->x) return -1;
	else if(aa->x > bb->x) return 1;
	else return 0;
}

double p0 = 0.2L;
double rho0 = 0.2L;
double w2Measure(double meand, double pressure,double rho){
	double res = 0.5L*meand*pow(pressure/p0*rho0/rho,0.5L);
	res = res*res;
	return res;
}
void getW2Ceil(gasparticle *bp, int np){
	int i;
	for(i=0;i<np;i++){

		double dleft, dright;
		if(i==0) {
			dleft = 2*(bp[i].x-xmin);
		}
		else { 
			dleft = bp[i].x - bp[i-1].x;
		}
		if(i == np-1){
			dright = 2*(xmax-bp[i].x);
		}
		else {
			dright = bp[i+1].x - bp[i].x;
		}

		bp[i].w2ceil = MIN(dleft, dright);
		bp[i].w2ceil = bp[i].w2ceil*bp[i].w2ceil; // square value
		bp[i].w2 = w2Measure(meand, bp[i].pressure, bp[i].den);
		bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
	}
}

double getw2(gasparticle *bp, int np, int inow){
	double w2;
	if(inow == -1){
		w2 = bp[0].w2;
	}
	else if(inow == np){
		w2 = bp[np-1].w2;
	}
	else w2 = bp[inow].w2;
	return w2;
}
double getpos(gasparticle *bp, int np, int inow){
	double xpos;
	if(inow == -1){
		xpos = xmin -(bp[0].x-xmin);
	}
	else if(inow == np){
		xpos = xmax+(xmax-bp[np-1].x);
	}
	else xpos = bp[inow].x;
	return xpos;
}

void getDen(gasparticle *bp, int np){
	int i;
	for(i=0;i<np;i++){
		double dpq,w2p, w2q;
		w2p = bp[i].w2;

		dpq = fabs(bp[i].x - getpos(bp, np, i-1));
		w2q = getw2(bp,np,i-1);
		double drad1 = 0.5*dpq*(1+(w2p - w2q)/(dpq*dpq));

		dpq = fabs(getpos(bp, np, i+1)-bp[i].x);
		w2q = getw2(bp,np,i+1);
		double drad2 = 0.5*dpq*(1+(w2p - w2q)/(dpq*dpq));

		bp[i].volume = drad1+drad2;
		bp[i].den = bp[i].mass/bp[i].volume;

		bp[i].drad1 = drad1;
		bp[i].drad2 = drad2;

		/*
        bp[i].pressure = bp[i].ie*bp[i].den*(Gamma-1);
        bp[i].poverrhogam = bp[i].pressure/pow(bp[i].den,Gamma);
		bp[i].Cs = sqrt(bp[i].pressure/bp[i].den*Gamma);
		*/
	}
}

double oDtime = 0.0845;
double Dtime;

double getAccVoro1D(gasparticle *bp, int np){
	int i,j;
	double area = 1.L;
	Dtime = 1.E20;
	int isave = -1;
	getW2Ceil(bp,np);
	getDen(bp,np);
	for(i=1;i<np-1;i++){
        double ui, upq, up,dpq, rvel,muij;
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
            dpq = dr;
            if(dr <0) dpq = -dr;
            up = bp[i].vx;
            upq = 0.5*(bp[i+j].vx + bp[i].vx); // velocity at the Lag. face
            double w2p = bp[i].w2;
            double w2q = bp[i+j].w2;
            ui = upq - bp[i].vx;


            double vol = bp[i].volume;
            double rscale = vol/2;

            rvel = ui * dr;
            muij = rvel*eta/(dr*dr + epsilon *epsilon * eta*eta);
            double wp = sqrt(w2p);
            double wq = sqrt(w2q);
            if(muij<0 && wp+wq > 0.5*dpq) {
                pij += pow(dpq/(wp+wq),0.5)*(-alpha*Cs*muij+beta*muij*muij) * denij;
            }
            die +=  -pij * ui*dS; /* for the internal energy */
            dte +=  -pij * upq*dS; /* for the total energy */
            dke +=  -pij * up*dS;
            fx += -pij *dS;

            double VdotR = (bp[i+j].vx - bp[i].vx)*dr/dpq;
            vsig = (0.5*(bp[i+j].Cs + bp[i].Cs) - MIN(0, VdotR));
            dt = Courant*dpq/vsig;
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
    return Dtime;
}
double vph_rk4(gasparticle *bp, int np){
	int i;
	double Dtime, dt;

	Dtime = getAccVoro1D(bp,np);
	for(i=1;i<np-1;i++){
		bp[i].rk4.k1x = bp[i].vx * Dtime;
		bp[i].rk4.k1vx = bp[i].ax * Dtime;
		bp[i].oden = bp[i].den;
		bp[i].ow2 = bp[i].w2;
		bp[i].odrad1 = bp[i].drad1;
		bp[i].odrad2 = bp[i].drad2;
	}
	for(i=1;i<np-1;i++){ // x => x + k1x/2
		bp[i].x += bp[i].rk4.k1x/2;
	}
	qsort(bp,np,sizeof(gasparticle),sortbp);
	dt = getAccVoro1D(bp,np);
	for(i=1;i<np-1;i++){
		bp[i].rk4.k2x = (bp[i].vx + bp[i].rk4.k1vx/2)*Dtime;
		bp[i].rk4.k2vx = bp[i].ax * Dtime;
	}

	for(i=1;i<np-1;i++){ // x => x + k2x/2
		bp[i].x += (bp[i].rk4.k2x - bp[i].rk4.k1x)/2;
	}
	qsort(bp,np,sizeof(gasparticle),sortbp);
	dt = getAccVoro1D(bp,np);
	for(i=1;i<np-1;i++){
		bp[i].rk4.k3x = (bp[i].vx + bp[i].rk4.k2vx/2)*Dtime;
		bp[i].rk4.k3vx = bp[i].ax * Dtime;
	}


	for(i=1;i<np-1;i++){  // x => x+k3x
		bp[i].x += bp[i].rk4.k3x -bp[i].rk4.k2x/2;
	}
	qsort(bp,np,sizeof(gasparticle),sortbp);
	dt = getAccVoro1D(bp,np);
	for(i=1;i<np-1;i++){
		bp[i].rk4.k4x = (bp[i].vx + bp[i].rk4.k3vx)*Dtime;
		bp[i].rk4.k4vx = bp[i].ax * Dtime;
	}
	
	for(i=1;i<np-1;i++){ // back to the original position 
		bp[i].x -= (bp[i].rk4.k3x);
		bp[i].w2 = bp[i].ow2; // recover the original density due to Lag. characters.
	}
	qsort(bp,np,sizeof(gasparticle),sortbp);
//	dt = getAccVoro1D(bp,np);
	for(i=1;i<np-1;i++){
		bp[i].x +=  (bp[i].rk4.k1x +2*bp[i].rk4.k2x +2*bp[i].rk4.k3x +bp[i].rk4.k4x )/6.;
		bp[i].vx += (bp[i].rk4.k1vx+2*bp[i].rk4.k2vx+2*bp[i].rk4.k3vx+bp[i].rk4.k4vx)/6.;
	}
	getDen(bp,np);
	for(i=1;i<np-1;i++){
		bp[i].ke += bp[i].dke*Dtime;
		bp[i].ie += bp[i].die *Dtime;
		bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1);
		bp[i].poverrhogam = bp[i].pressure/pow(bp[i].den,Gamma);
		bp[i].Cs = sqrt(Gamma*bp[i].pressure/bp[i].den);
	}

	qsort(bp,np,sizeof(gasparticle),sortbp);

	return Dtime;
}

double vph(gasparticle *bp, int np){
	int i,j;
	double area = 1.L;
	Dtime = 1.E20;
	int isave = -1;

	getW2Ceil(bp, np);
	getDen(bp,np);


	for(i=1;i<np-1;i++){
		double ui, upq, up,dpq, rvel,muij;
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
			dpq = dr;
			if(dr <0) dpq = -dr;
			up = bp[i].vx;
			upq = 0.5*(bp[i+j].vx + bp[i].vx); // velocity at the Lag. face
			double w2p = bp[i].w2;
			double w2q = bp[i+j].w2;
			ui = upq - bp[i].vx;


			double vol = bp[i].volume;
			double rscale = vol/2;

			rvel = ui * dr;
			muij = rvel*eta/(dr*dr + epsilon *epsilon * eta*eta);
			double wp = sqrt(w2p);
			double wq = sqrt(w2q);
			if(muij<0 && wp+wq > 0.5*dpq) {
				pij += pow(dpq/(wp+wq),0.5)*(-alpha*Cs*muij+beta*muij*muij) * denij;
			}
			die +=  -pij * ui*dS; /* for the internal energy */
			dte +=  -pij * upq*dS; /* for the total energy */
			dke +=  -pij * up*dS;
			fx += -pij *dS;

			double VdotR = (bp[i+j].vx - bp[i].vx)*dr/dpq;
			vsig = (0.5*(bp[i+j].Cs + bp[i].Cs) - MIN(0, VdotR));
			dt = Courant*dpq/vsig;
			bp[i].dt = MIN(bp[i].dt, dt);
			if(dt <  Dtime){
				oDtime = Dtime;
				Dtime = dt;
				isave = i;
			}
		
		}
		if(i==80){
			fprintf(stderr,"p80 has %g %g %g : den= %g %g %g\n", bp[i].x, bp[i].vx, fx,
					bp[i].den, bp[i].pressure, bp[i].ie);
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
	getDen(bp,np);

	for(i=1;i<np-1;i++){
        bp[i].ie += bp[i].die*Dtime; 
        bp[i].ke += bp[i].dke*Dtime;
        bp[i].pressure = bp[i].ie/bp[i].volume*(Gamma-1);
        bp[i].poverrhogam = bp[i].pressure/pow(bp[i].den,Gamma);
		bp[i].Cs = sqrt(bp[i].pressure/bp[i].den*Gamma);
	}


	return Dtime;
}


int main(int argc, char **argv){

	int np=100;
	gasparticle *bp;

	targetT = atof(argv[1]);

	bp = (gasparticle*)malloc(sizeof(gasparticle)*np);
	double xsep = 1./np;

	int i;
	if(0){
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
			bp[i].ie = bp[i].pressure/ ( (Gamma-1)*bp[i].den )*bp[i].volume*bp[i].den;
			bp[i].vx = 0;
			bp[i].id = i;
		}
	}
	else { /* Gizmo setting */
		double Lden,Rden,LP,RP;
        Lden = 1;
        Rden = 0.25;
        LP = 1;
        RP = 0.1795;
        xsep = 10.L/(np*4/5);
        xmin = -10;
        xmax = 10;
		meand = (xmax-xmin)/np;
        for(i=0;i<np*4/5;i++){
            bp[i].x = xmin + xsep*(i+0.5);
            bp[i].pressure = LP;
            bp[i].den = Lden;
            bp[i].volume = xsep;
            bp[i].mass = xsep*Lden;
        }
        int istart = np*4/5;
        xsep = xmax/(np/5);
        for(i=0;i<np/5;i++){
            bp[i+istart].x = xsep*(i+0.5);
            bp[i+istart].pressure = RP;
            bp[i+istart].den = Rden;
            bp[i+istart].volume = xsep;
            bp[i+istart].mass = xsep*Rden;
        }
        for(i=0;i<np;i++){
			bp[i].Cs = sqrt(bp[i].pressure/bp[i].den*Gamma);
            bp[i].ie = bp[i].pressure/ ( (Gamma-1)*bp[i].den )*bp[i].volume*bp[i].den;
            bp[i].vx = 0;
            bp[i].id = i;
			bp[i].w2 = w2Measure(meand, bp[i].pressure, bp[i].den);
			double dleft = bp[i].x - getpos(bp,np,(i-1));
			double dright = getpos(bp,np,(i+1)) - bp[i].x;
			bp[i].w2ceil = MIN(dleft, dright);
			bp[i].w2ceil = bp[i].w2ceil*bp[i].w2ceil; // square value
			bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
        }

		/*
		for(i=0;i<np;i++){
			double dpq,w2p, w2q;
			w2p = bp[i].w2;

			dpq = fabs(bp[i].x - getpos(bp, np, i-1));
			w2q = getw2(bp,np,i-1);
			double drad1 = 0.5*dpq*(1+(w2p - w2q)/(dpq*dpq));

			dpq = fabs(getpos(bp, np, i+1)-bp[i].x);
			w2q = getw2(bp,np,i+1);
			double drad2 = 0.5*dpq*(1+(w2p - w2q)/(dpq*dpq));

			bp[i].volume = drad1+drad2;
			bp[i].mass = bp[i].den*bp[i].volume;

		}
		*/
    }

	for(i=1;i<np-1;i++){
		bp[i].volume = 0.5*(bp[i+1].x - bp[i-1].x);
		bp[i].mass = bp[i].den * bp[i].volume;
	}
	time_t time1, time2;
	time1 = clock();
	double time=0;
	int jcount = 0;
	do {
		time += vph_rk4(bp,np);
		fprintf(stderr,"Time= %g @stepi= %d\n", time, jcount);
		jcount ++;
	} while (time < targetT);
	time2 = clock();
	fprintf(stderr,"gizmo: time step count = %d :: %15.7lg\n", jcount, (double)(time2-time1)/CLOCKS_PER_SEC);
	for(i=0;i<np;i++){
		printf("%12.7g %12.7g %12.7g %12.7g %5d %12.7g %12.7g %12.7g %12.7g\n",bp[i].x,bp[i].vx,bp[i].den,bp[i].pressure,bp[i].id,bp[i].mass,bp[i].viscosity,bp[i].ax, bp[i].poverrhogam);
	}
}



