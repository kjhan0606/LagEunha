#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<time.h>

#define MAX(a,b) ( (a)>(b)?(a):(b) )
#define MIN(a,b) ( (a)<(b)?(a):(b) )

#define NTIME 1000
double Gamma = 1.4;
double dtime = 0.2/NTIME;
double targetT = 0.2;
double ntime=0;
typedef struct gasparticle{
	int id;
	double x,den,mass,pressure,vx,ax,energy,ke,ie;
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
	int i;
	double area = 1.L;
	for(i=1;i<np-1;i++){
		bp[i].den = bp[i].mass/(bp[i+1].x - bp[i-1].x)*2;
		bp[i].ke = 0.5*bp[i].vx*bp[i].vx;
		bp[i].ie = bp[i].pressure/bp[i].den/(Gamma-1);
	}
	for(i=1;i<np-1;i++){
		double d = bp[i+1].x-bp[i-1].x;
		double d1 = bp[i].x-bp[i-1].x;
		double d2 = bp[i+1].x-bp[i].x;

		bp[i].ax = 0;

		if(1) {
			bp[i].ax += -(bp[i+1].pressure+bp[i].pressure)/2.L/bp[i].mass * area;
			bp[i].ax += (bp[i-1].pressure+bp[i].pressure)/2.L/bp[i].mass * area;
			double vis1, vis2;
			vis1 = vis2 = 0;
			if(bp[i+1].vx - bp[i].vx <0) vis2 = pow(bp[i+1].vx - bp[i].vx,1.L)/(2*d2);
			if(bp[i-1].vx - bp[i].vx >0) vis1 = pow(bp[i-1].vx - bp[i].vx,1.L)/(2*d1);
//			double viscosity = ( (bp[i+1].vx-bp[i].vx)/(2*d2) + (bp[i-1].vx-bp[i].vx)/(2*d1))*0.5*area;
			double viscosity = (vis1+vis2)*2*area;
			bp[i].ax += viscosity;
			bp[i].viscosity = viscosity;
		}
		else {
			bp[i].ax += -(bp[i+1].pressure)/bp[i].mass * area;
			bp[i].ax += (bp[i-1].pressure)/bp[i].mass * area;
			double viscosity = ( (bp[i+1].vx-bp[i].vx)/(2*d2) + (bp[i-1].vx-bp[i].vx)/(2*d1))*0.5*area;
			bp[i].ax += viscosity;
			bp[i].viscosity = viscosity;
		}
	}
	double Dtime = 1.E20;

	for(i=1;i<np-1;i++){
		double dt,dl,dv,a;

		dl = bp[i+1].x - bp[i].x;
		dl = dl * 0.05;
		a = bp[i].ax - bp[i+1].ax;
		dv = bp[i].vx - bp[i+1].vx;
		if(a ==0 ){
			if(dv>0) {
				dt = dl/dv;
			}
			else 
				continue;
		}
		else {
			double tmp = dv*dv/a/a + 2/a*dl;
			if(tmp>0){
				double dt1,dt2;
				dt1 = -dv/a + sqrt(tmp);
				dt2 = -dv/a - sqrt(tmp);
				if(dt2 >0) dt = dt2;
				else if(dt1>0) dt = dt1;
				else continue;
			}
			else 
				continue;
		}
		if(dt < Dtime) {
			Dtime = dt;
		}
	}
	for(i=1;i<np-1;i++){
		bp[i].x += bp[i].vx*Dtime + 0.5*bp[i].ax*Dtime*Dtime;
		bp[i].ovx = bp[i].vx;
		bp[i].vx += bp[i].ax * Dtime;
	}
	qsort(bp,np,sizeof(gasparticle), sortbp);
	for(i=1;i<np-1;i++){
		bp[i].volume = 0.5*(bp[i+1].x - bp[i-1].x);
	}
	bp[0].volume = bp[1].x - bp[0].x;
	bp[np-1].volume = bp[np-1].x - bp[np-2].x;
	for(i=1;i<np-1;i++){
		double volume;
//		volume = (bp[i+1].x-bp[i-1].x)*0.5;
		volume = bp[i].volume;
		double den = bp[i].mass/volume;
		/* to solve the energy conservation */
		if(0){
			double ui = 0.5*(bp[i+1].vx + bp[i].vx);
			double uj = 0.5*(bp[i-1].vx + bp[i].vx);
			double pi = 0.5*(bp[i].pressure+bp[i+1].pressure);
			double pj = 0.5*(bp[i].pressure+bp[i-1].pressure);
			double de = (-pi*ui+pj*uj)*Dtime; /* Increase of total energy for time interval Dtime; 1.33 */

			double dke= (0.5*bp[i].vx*bp[i].vx - bp[i].ke)*bp[i].mass;
			double die = de-dke; /* increase of internal energy */
			bp[i].ie += die/volume/den; /* increase of internal energy per unit density & and volume*/
			bp[i].volume = volume;
			bp[i].pressure = bp[i].ie*den*(Gamma-1);
			bp[i].poverrhogam = bp[i].pressure/pow(den,Gamma);
		}
		else {
			double ui = 0.5*(bp[i+1].vx + bp[i].vx);
			double uj = 0.5*(bp[i-1].vx + bp[i].vx);
			double pi = 0.5*(bp[i].pressure+bp[i+1].pressure);
			double pj = 0.5*(bp[i].pressure+bp[i-1].pressure);
			double de = (-pi*ui+pj*uj)*Dtime; /* Increase of total energy for time interval Dtime; 1.33 */
			double dke= (0.5*bp[i].vx*bp[i].vx - bp[i].ke)*bp[i].mass;
			double die = de-dke; /* increase of internal energy */
			bp[i].ie += die/volume/den; /* increase of internal energy per unit density & and volume*/
			bp[i].pressure = bp[i].ie*den*(Gamma-1);
			bp[i].volume = volume;
			bp[i].poverrhogam = bp[i].pressure/pow(den,Gamma);
		}
		bp[i].den = den;
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
	for(i=0;i<np;i++){
		bp[i].x = xsep*(i+0.5);
		if(bp[i].x<=0.5){
			bp[i].mass = xsep;
			bp[i].volume = xsep;
			bp[i].pressure = 1.L;
			bp[i].den = 1.L;
		}
		else {
			bp[i].mass = 0.125L*xsep;
			bp[i].volume = 0.125L*xsep;
			bp[i].pressure = 0.1L;
			bp[i].den = 0.125L;
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
		jcount ++;
	} while (time < targetT);
	time2 = clock();
	fprintf(stderr,"time step count = %d :: %g\n", jcount, (float)(time2-time1)/CLOCKS_PER_SEC);
	for(i=0;i<np;i++){
		printf("%g %g %g %g %d %g %g %g %g\n",bp[i].x,bp[i].vx,bp[i].den,bp[i].pressure,bp[i].id,bp[i].mass,bp[i].viscosity,bp[i].ax, bp[i].poverrhogam);
	}
}



