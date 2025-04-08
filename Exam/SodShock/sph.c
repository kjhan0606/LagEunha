/*
icc -o sph.exe sph.c -lm
*/
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<time.h>

/*
#define MAX(a,b) ( (a)>(b)?(a):(b) )
#define MIN(a,b) ( (a)<(b)?(a):(b) )
*/

//#define NTIME 1000
//double dtime = 0.2/NTIME;
//double ntime=0;
//
double GAMMAA = 1.4;
double targetT = 0.2;
int nneigh = 5;
typedef struct gasparticle{
	int id;
	double x,den,mass,pressure,vx,ax,energy,ke,ie,hsml,du;
	double ovx;
	double volume,weight;
	double poverrhogam;
	double viscosity,csound;
} gasparticle;

int sortdist(const void *a, const void *b){
	float *aa, *bb;
	aa = (float*)a;
	bb = (float*)b;
	if(*aa < *bb) return -1;
	else if(*aa > *bb) return 1;
	else return 0;
}
int sortbp(const void *a, const void *b){
	gasparticle *aa, *bb;
	aa = (gasparticle*)a;
	bb = (gasparticle*)b;
	if(aa->x < bb->x) return -1;
	else if(aa->x > bb->x) return 1;
	else return 0;
}
float dW4(float r, float h){
    float q = r/h;
    if(q<1) {
        return 2./3./h/h*(-3*q +2.25*q*q);
    }
    else if (q<2){
        return -2./3./h/h*0.75*(2.-q)*(2.-q);
    }
    else {
        return 0;
    }
}

float W4(float r, float h){
	float q = r/h;
	if(q<1) {
		return 2./3./h*(1-1.5*q*q +0.75*q*q*q);
	}
	else if (q<2){
		return 2./3./h*0.25*(2.-q)*(2.-q)*(2.-q);
	}
	else {
		return 0;
	}
}

void gethsml(gasparticle *bp, int np, int inow, float *rdist){
	if(nneigh%2 ==0){
		fprintf(stderr,"Please input odd number for the neighbors: %d\n",nneigh);
		exit(99);
	}
	rdist[0] =0;
	int ndist = 1;
	int i;
	for(i=1;i<=nneigh;i++){
		int ii,ij;
		ii = inow-i;
		ij = inow+i;
		if(ii<0) {
			ii = -ii-1;
			rdist[ndist++] = (bp[inow].x +  bp[ii].x);
		}
		else if(ii<np) {
			rdist[ndist++] = fabs(bp[inow].x-bp[ii].x);
		}
		else if(ii>=np){
			ii = np - (ii-np) -1;
			rdist[ndist++] = fabs(1-bp[inow].x+1-bp[ii].x);
		}
		if(ij<0) {
			ij = -ij-1;
			rdist[ndist++] = (bp[inow].x + bp[ij].x);
		}
		else if(ij<np) {
			rdist[ndist++] = fabs(bp[inow].x-bp[ij].x);
		}
		else if(ij>=np){
			ij = np - (ij-np) -1;
			rdist[ndist++] = fabs(1-bp[inow].x+1-bp[ij].x);
		}
	}
	qsort(rdist,ndist,sizeof(float), sortdist);
	bp[inow].hsml = rdist[nneigh-1]*0.5;

	
}
float getax(gasparticle *bp, int np, int inow){
	float ax = 0;
	float alpha = 1; /* Monaghan 1992 ARAA 30, 543 */
	float beta = 2;
	float AV = 0;
	float dist;
	double eta = 0.1*bp[inow].hsml;
	int i,ii;
	for(i=1;i<nneigh;i++){
		for(ii=inow-i;ii<=inow+i; ii+=2*i){
			if(ii<0) {
				int ij = -ii-1;
				float dist = bp[inow].x + bp[ij].x;
				float hav = 0.5*(bp[inow].hsml+bp[ij].hsml);
				int ifactor = 1;
				if(dist < hav*2) { 
					double cab,denab,PIab;
//					double muab = bp[inow].hsml * (bp[inow].vx+bp[ij].vx)*(bp[inow].x+bp[ij].x)/(dist*dist+eta*eta); 
					eta = 0.1*hav;
					double muab = hav * (bp[inow].vx+bp[ij].vx)*(bp[inow].x+bp[ij].x)/(dist*dist+eta*eta); 
					PIab = 0; 
					if(muab <0) { 
						cab = 0.5*(bp[inow].csound + bp[ij].csound); 
						denab = 0.5*(bp[inow].den + bp[ij].den); 
						PIab = (-alpha*cab*muab + beta*muab*muab)/denab; 
					} 
					ax -= bp[ij].mass*(bp[ij].pressure/bp[ij].den/bp[ij].den+ 
							bp[inow].pressure/bp[inow].den/bp[inow].den + PIab)*dW4(dist,hav)*ifactor; 
					AV += PIab;
				}
			}
			else if(ii<np){
				float dist = fabs(bp[ii].x - bp[inow].x);
				int ifactor;
				if(bp[inow].x > bp[ii].x) ifactor = 1;
				else ifactor =-1;
				float hav = 0.5*(bp[inow].hsml+bp[ii].hsml);
				if(fabs(dist) < hav*2) { 
					double cab,denab,PIab; 
//					double muab = bp[inow].hsml * (bp[inow].vx-bp[ii].vx)*(bp[inow].x-bp[ii].x)/(dist*dist+eta*eta); 
					eta = 0.1*hav;
					double muab = hav * (bp[inow].vx-bp[ii].vx)*(bp[inow].x-bp[ii].x)/(dist*dist+eta*eta); 
					PIab = 0; 
					if(muab <0) { 
						cab = 0.5*(bp[inow].csound + bp[ii].csound); 
						denab = 0.5*(bp[inow].den + bp[ii].den); 
						PIab = (-alpha*cab*muab + beta*muab*muab)/denab; 
					} 
					ax -= bp[ii].mass*(bp[ii].pressure/bp[ii].den/bp[ii].den+ 
							bp[inow].pressure/bp[inow].den/bp[inow].den + PIab)*dW4(dist,hav)*ifactor; 
					AV += PIab;
				}
			}
			else if(ii>=np) {
				int	ij = np - (ii-np) -1;
				float dist = fabs(1 - bp[inow].x + 1 -bp[ij].x);
                float hav = 0.5*(bp[inow].hsml+bp[ij].hsml);
				int ifactor = -1;
                if(dist < hav*2) {
                    double cab,denab,PIab;
//                    double muab = bp[inow].hsml * (bp[inow].vx+bp[ij].vx)*(-dist)/(dist*dist+eta*eta);
					eta = 0.1*hav;
                    double muab = hav * (bp[inow].vx+bp[ij].vx)*(-dist)/(dist*dist+eta*eta);
                    PIab = 0;
                    if(muab <0) {
                        cab = 0.5*(bp[inow].csound + bp[ij].csound);
                        denab = 0.5*(bp[inow].den + bp[ij].den);
                        PIab = (-alpha*cab*muab + beta*muab*muab)/denab;
                    }
                    ax -= bp[ij].mass*(bp[ij].pressure/bp[ij].den/bp[ij].den+
                            bp[inow].pressure/bp[inow].den/bp[inow].den + PIab)*dW4(dist,hav)*ifactor;
                    AV += PIab;
                }
            }
		}
	}

	bp[inow].viscosity = AV;

	return ax;
}

float getdu(gasparticle *bp, int np, int inow){
    float du = 0;
    float alpha = 1; /* Monaghan 1992 ARAA 30, 543 */
    float beta = 2;
    float AV = 0;
    float dist;
    double eta = 0.1*bp[inow].hsml;
    int i,ii;
    for(i=1;i<nneigh;i++){
        for(ii=inow-i;ii<=inow+i; ii+=2*i){
            if(ii<0) {
                int ij = -ii-1;
                float dist = bp[inow].x + bp[ij].x;
                float hav = 0.5*(bp[inow].hsml+bp[ij].hsml);
				float dvel = bp[inow].vx + bp[ij].vx;
				int ifactor = 1;
                if(dist < hav*2) {
                    double cab,denab,PIab;
//                    double muab = bp[inow].hsml * (bp[inow].vx+bp[ij].vx)*(bp[inow].x+bp[ij].x)/(dist*dist+eta*eta);
					eta = 0.1*hav;
                    double muab = hav * (bp[inow].vx+bp[ij].vx)*(bp[inow].x+bp[ij].x)/(dist*dist+eta*eta);
                    PIab = 0;
                    if(muab <0) {
                        cab = 0.5*(bp[inow].csound + bp[ij].csound);
                        denab = 0.5*(bp[inow].den + bp[ij].den);
                        PIab = (-alpha*cab*muab + beta*muab*muab)/denab;
                    }
                    du += 0.5*bp[ij].mass*(bp[ij].pressure/bp[ij].den/bp[ij].den+
                            bp[inow].pressure/bp[inow].den/bp[inow].den + PIab)*dW4(dist,hav)*dvel*ifactor;
                    AV += PIab;
                }
            }
            else if(ii<np){
                float dist = fabs(bp[ii].x - bp[inow].x);
                int ifactor;
                if(bp[inow].x > bp[ii].x) ifactor = 1;
                else ifactor = -1;
                float hav = 0.5*(bp[inow].hsml+bp[ii].hsml);
				float dvel = bp[inow].vx - bp[ii].vx;
                if(fabs(dist) < hav*2) { 
					double cab,denab,PIab; 
//					double muab = bp[inow].hsml * (bp[inow].vx-bp[ii].vx)*(bp[inow].x-bp[ii].x)/(dist*dist+eta*eta); 
					eta = 0.1*hav;
					double muab = hav * (bp[inow].vx-bp[ii].vx)*(bp[inow].x-bp[ii].x)/(dist*dist+eta*eta); 
					PIab = 0; 
					if(muab <0) { 
						cab = 0.5*(bp[inow].csound + bp[ii].csound); 
						denab = 0.5*(bp[inow].den + bp[ii].den); 
						PIab = (-alpha*cab*muab + beta*muab*muab)/denab; 
					} 
					du += 0.5*bp[ii].mass*(bp[ii].pressure/bp[ii].den/bp[ii].den+ 
							bp[inow].pressure/bp[inow].den/bp[inow].den + PIab)*dW4(dist,hav)*ifactor*dvel; 
					AV += PIab;
   	             }
            }
			else if(ii>=np){
				int ij = np - (ii-np)-1;
                float dist = (1- bp[ij].x) + 1 -bp[inow].x;
                float hav = 0.5*(bp[inow].hsml+bp[ij].hsml);
                float dvel = bp[inow].vx + bp[ij].vx;
				int ifactor = -1;
                if(dist < hav*2) { 
					double cab,denab,PIab;
//                    double muab = bp[inow].hsml * (bp[inow].vx+bp[ij].vx)*(-dist)/(dist*dist+eta*eta);
					eta = 0.1*hav;
                    double muab = hav * (bp[inow].vx+bp[ij].vx)*(-dist)/(dist*dist+eta*eta);
                    PIab = 0;
                    if(muab <0) { 
						cab = 0.5*(bp[inow].csound + bp[ij].csound);
                        denab = 0.5*(bp[inow].den + bp[ij].den);
                        PIab = (-alpha*cab*muab + beta*muab*muab)/denab;
                    }
                    du += 0.5*bp[ij].mass*(bp[ij].pressure/bp[ij].den/bp[ij].den+
                            bp[inow].pressure/bp[inow].den/bp[inow].den + PIab)*dW4(dist,hav)*dvel*ifactor;
                    AV += PIab;
                }
			}
        }
    }


    return du;
}

void  getden(gasparticle *bp, int np, int inow){
    float den = 0;
	float hav = bp[inow].hsml;
	den = bp[inow].mass*W4(0.,hav);
	int i,ii;
	for(i=1;i<nneigh;i++){
		for(ii=inow-i;ii<=inow+i;ii+=2*i){
			if(ii<0) {
				int ij = -ii-1;
				float dist = fabs(bp[inow].x + bp[ij].x);
				hav = 0.5*(bp[inow].hsml+bp[ij].hsml);
				if(dist < hav*2) den += bp[ij].mass*W4(dist,hav);
			}
			else if(ii<np){
				float dist = fabs(bp[ii].x - bp[inow].x);
				hav = 0.5*(bp[inow].hsml+bp[ii].hsml);
				if(dist < hav*2) den += bp[ii].mass*W4(dist,hav);
			}
			else if(ii>=np){
				int	ij = np - (ii-np) -1;
				float dist = fabs(1 - bp[inow].x + 1 -bp[ij].x);
				hav = 0.5*(bp[inow].hsml+bp[ij].hsml);
				if(dist < hav*2) den += bp[ij].mass*W4(dist,hav);
			}
		}
	}
	bp[inow].den = den;
}



double sph(gasparticle *bp, int np){
	int i;
	double area = 1.L;
	float *rdist = (float*)malloc(sizeof(float)*np);
	for(i=0;i<np;i++) gethsml(bp,np, i,rdist);
	for(i=0;i<np;i++){
		getden(bp,np,i);
		bp[i].csound = sqrt(GAMMAA*bp[i].pressure/bp[i].den);
		bp[i].ke = 0.5*bp[i].vx*bp[i].vx;
		bp[i].pressure = bp[i].ie*bp[i].den*(GAMMAA-1);
	}
	for(i=0;i<np;i++){
		bp[i].ax = 1*getax(bp,np,i);
		bp[i].du = 1.*getdu(bp,np,i);
	}
	double Dtime = 1.E20;

	for(i=1;i<np-1;i++){
		double dt,dl,dv,a;
		dt = 1.e20L;

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
		bp[i].ie += bp[i].du*Dtime; /* increase of internal energy per unit density & and volume*/
		bp[i].pressure = bp[i].ie*bp[i].den*(GAMMAA-1);
		bp[i].poverrhogam = bp[i].pressure/pow(bp[i].den,GAMMAA);
	}

	return Dtime;
}


int main(int argc, char **argv){

	int np=100;
	gasparticle *bp;

	if(argc !=2){
		fprintf(stderr,"Please input as: ./sph.exe [time]\n");
		exit(99);
	}
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
		bp[i].ie = bp[i].pressure/ ( (GAMMAA-1)*bp[i].den );
		bp[i].vx = 0;
		bp[i].id = i;
	}
	for(i=1;i<np-1;i++){
		bp[i].volume = 0.5*(bp[i+1].x - bp[i-1].x);
	}
	double time=0;
	int jcount = 0;
	time_t time1, time2;
	time1 = clock();
	while(time < targetT) {
		time += sph(bp,np);
		jcount ++;
	} 
	time2 = clock();
	fprintf(stderr,"time step count = %d:   %lg\n", jcount, (double)(time2-time1)/CLOCKS_PER_SEC);
	for(i=0;i<np;i++){
		printf("%12.7g %12.7g %12.7g %12.7g %5d %12.7g %12.7g %12.7g %12.7g\n",bp[i].x,bp[i].vx,bp[i].den,bp[i].pressure,bp[i].id,bp[i].mass,bp[i].viscosity,bp[i].ax, bp[i].poverrhogam);
	}
}



