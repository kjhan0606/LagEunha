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

#include "init.h"

//double xmin,xmax;

//double Gamma = 1.4;
double dtime = 0.2/NTIME;
double targetT = 0.2;
double ntime=0;
// ie is the total internal energy of the cell.
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

double vph(gasparticle *bp, int np){
	int i,j;
	double area = 1.L;
	double Dtime = 1.E20;
	int isave = -1;

	/*
	for(i=1;i<np-1;i++){
		bp[i].dt = Dtime;
	}
	*/
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

			double vol1, vol2;
			vol1 = 0.5*(bp[i+1].x -bp[i-1].x);
			vol2 = 0.5*(bp[i+j+1].x -bp[i+j-1].x);

			if(i+j+1==np) {
				vol2 = xmax-bp[i+j].x + 0.5*(bp[i+j].x-bp[i+j-1].x);
			}
			else if(i+j-1==-1) {
				vol2 = bp[i+j].x-xmin + 0.5*(bp[i+j+1].x-bp[i+j].x);
			}
			else {
				vol2 = 0.5*(bp[i+j+1].x -bp[i+j-1].x);
			}
			double vol = 0.5*(vol1+vol2);
			double rscale = vol/2;
			rvel = ui * dr;
			muij = rvel*eta/(dr*dr + epsilon *epsilon * eta*eta);
			if(muij<0) pij = pij + (-alpha*Cs*muij + beta * muij*muij)*denij;
			die +=  -pij * ui*dS; /* for the internal energy */
			dte +=  -pij * ua*dS; /* for the total energy */
			dke +=  -pij * ub*dS;
			fx += -pij *dS;
			
			double VdotR = (bp[i+j].vx - bp[i].vx)*dr/dramp;
			vsig = ((bp[i+j].Cs + bp[i].Cs)*0.5 - MIN(0, VdotR));
			dt = Courant*dramp/vsig;
			bp[i].dt = MIN(bp[i].dt, dt);
			if(dt <  Dtime){
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
	for(i=1;i<np-1;i++){
		bp[i].volume = 0.5*(bp[i+1].x - bp[i-1].x);
		bp[i].den = bp[i].mass/bp[i].volume;
		/*
        bp[i].ie += bp[i].die/bp[i].volume/bp[i].den*Dtime; 
        bp[i].ke += bp[i].dke/bp[i].volume/bp[i].den*Dtime;
        bp[i].pressure = bp[i].ie*bp[i].den*(pgamma-1);
		*/
        bp[i].ie += bp[i].die*Dtime; 
        bp[i].ke += bp[i].dke*Dtime;
        bp[i].pressure = bp[i].ie/bp[i].volume*(pgamma-1);

        bp[i].poverrhogam = bp[i].pressure/pow(bp[i].den,pgamma);
		bp[i].Cs = sqrt(bp[i].pressure/bp[i].den*pgamma);
	}


	return Dtime;
}


int main(int argc, char **argv){

	int np=200;
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
				bp[i].Cs = sqrt(bp[i].pressure/bp[i].den*pgamma);
			}
			else {
				bp[i].mass = 0.125L*xsep;
				bp[i].volume = 0.125L*xsep;
				bp[i].pressure = 0.1L;
				bp[i].den = 0.125L;
				bp[i].Cs = sqrt(bp[i].pressure/bp[i].den*pgamma);
			}
			bp[i].ie = bp[i].pressure/ ( (pgamma-1)*bp[i].den ) * bp[i].volume *bp[i].den;
			bp[i].vx = 0;
			bp[i].id = i;
		}
	}
	else { /* Gizmo setting */
//		double Lden,Rden,LP,RP;
		/*
        Lden = 1;
        Rden = 0.25;
        LP = 1;
        RP = 0.1795;
        xmin = -10;
        xmax = 10;
		*/


        xsep = 10.L/(np*4/5);
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
			bp[i].Cs = sqrt(bp[i].pressure/bp[i].den*pgamma);
            bp[i].ie = bp[i].pressure/ ( (pgamma-1)*bp[i].den ) * bp[i].volume *bp[i].den;
            bp[i].vx = 0;
            bp[i].id = i;
        }
    }

	for(i=1;i<np-1;i++){
		bp[i].volume = 0.5*(bp[i+1].x - bp[i-1].x);
		bp[i].mass = bp[i].den*bp[i].volume;
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
	fprintf(stderr,"gizmo: time step count = %d :: %15.7lg\n", jcount, (double)(time2-time1)/CLOCKS_PER_SEC);
	for(i=0;i<np;i++){
		printf("%12.7g %12.7g %12.7g %12.7g %5d %12.7g %12.7g %12.7g %12.7g\n",bp[i].x,bp[i].vx,bp[i].den,bp[i].pressure,bp[i].id,bp[i].mass,bp[i].viscosity,bp[i].ax, bp[i].poverrhogam);
	}
}



