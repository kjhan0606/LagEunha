#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>


#define P_G 6.673E-8L
#define Mpccgs 3.085677488036136E24L
#define kpcincgs 3.085677488036136E21L
#define onesolarmass 1.989E33L

#include "pmheader.h"


#define NSIZE 128

double accel[NSIZE];

#define GetExtAccel(type,X0,Y0,Z0) {\
	double x,y,z;\
	long i;\
	double rscale,ellinkpc,ellinMpc;\
	ellinMpc = simpar.boxsize/simpar.nx;\
	ellinkpc = ellinMpc * 1000.;\
	rscale = 1./ellinkpc;\
	for(i=0;i<simpar.##type##.np;i++){\
		x = (simpar.##type##.u.tbp+i)->x-X0;\
		y = (simpar.##type##.u.tbp+i)->y-Y0;\
		z = (simpar.##type##.u.tbp+i)->z-Z0;\
		double r = sqrt(x*x+y*y+z*z);\
		double lnr = 22*(log10(r/rscale) + 3);\
		int ilnr = lnr;\
		float acc;\
		if(ilnr >=0){\
			acc = accel[ilnr] + (accel[ilnr+1]-accel[ilnr]) *(lnr-ilnr);\
		}\
		acc = acc / r;\
		float ax,ay,az;\
		ax = x * acc;\
		ay = y * acc;\
		az = z * acc;\
		(simpar.##type##.u.tbp+i)->ax += ax;\
		(simpar.##type##.u.tbp+i)->ay += ay;\
		(simpar.##type##.u.tbp+i)->az += az;\
	}\
}


int iflag = 1;

void readextforce_(double *,int *);

void External_Force(float X0,float Y0,float Z0){
	int i,nsize=NSIZE;
	long np;


	if(iflag){
		double fscale,pixmass,ellinkpc,mscale,ellMpc,fs,fr,pmassinMsun,meancellmass;
		pixmass = simpar.totalmassinboxinsolarmass/simpar.nx/simpar.ny/simpar.nz;
		mscale = 1./pixmass;
		ellMpc = simpar.boxsize/simpar.nx;
		ellinkpc = simpar.boxsize/simpar.nx * 1000.L;
		/* in M_sun */
		pmassinMsun = 61569.2; 
		meancellmass = pmassinMsun / pow(simpar.nx/100.L,3.L);
		fs  = P_G*meancellmass * onesolarmass/pow(ellinkpc*kpcincgs, 2.L);
		fr =  1.E10L / kpcincgs; /* coefficient in (cm/sec)^2/cm unit */
		fscale = fr/fs;
	  	if(simpar.myid == 0) readextforce_(accel,&nsize);
		MPI_Bcast(accel,nsize,MPI_DOUBLE,0,MPI_COMM_WORLD);
		for(i=0;i<NSIZE;i++) accel[i] = accel[i]*fscale/simpar.nx/simpar.nx;
#ifdef DEBUG
		if(simpar.myid==0){
			float r;
			double rscale,ellinkpc,ellinMpc;
			ellinMpc = simpar.boxsize/simpar.nx;
			ellinkpc = ellinMpc * 1000.;
			rscale = 1./ellinkpc;
			for(i=0;i<128;i++){
				r = i+0.0001;
				float lnr = 22*(log10(r/rscale) + 3);
				int ilnr = lnr;
				float acc;
				acc = accel[ilnr] + (accel[ilnr+1]-accel[ilnr]) *(lnr-ilnr);
				printf("r= %g force= %g : fscale= %g\n",r,acc,fscale);
			}
		}
		MPI_Finalize();exit(0);
#endif
		iflag = 0;
	}
	/*
	if(simpar.myid==0){
		double rscale,ellinkpc,ellinMpc;
		ellinMpc = simpar.boxsize/simpar.nx;
		ellinkpc = ellinMpc * 1000.;
		rscale = 1./ellinkpc;
		for(i=0;i<256;i++){
			float r=i;
			float lnr = 22*(log10(r)/rscale + 3);
			int ilnr = lnr;
			float acc;
			if(ilnr >=0){
				acc = accel[ilnr] + (accel[ilnr+1]-accel[ilnr]) *(lnr-ilnr);
			}
			printf("%g %g\n",r,acc);
		}
	}
	MPI_Finalize();
	exit(99);
	*/

	GetExtAccel(dm,X0,Y0,Z0);
	GetExtAccel(sph,X0,Y0,Z0);
	GetExtAccel(star,X0,Y0,Z0);
}


