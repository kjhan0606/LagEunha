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

#include "eunha.h"


#define NSIZE 128

static double accel[NSIZE];

#define GetExtAccel(TYPE,x00,y00,z00) do {\
	long i;\
	double rscale,ellinkpc,ellinMpc;\
	ellinMpc = BOXSIZE(simpar)/NX(simpar);\
	ellinkpc = ellinMpc * 1000.;\
	rscale = 1./ellinkpc;\
	for(i=0;i<TYPE##_NP(simpar);i++){\
		double x,y,z;\
		x = (TYPE##_TBP(simpar)+i)->x-x00;\
		y = (TYPE##_TBP(simpar)+i)->y-y00;\
		z = (TYPE##_TBP(simpar)+i)->z-z00;\
		double r = sqrt(x*x+y*y+z*z);\
		double lnr = 22*(log10(r/rscale) + 3);\
		int ilnr = lnr;\
		float acc;\
		if(ilnr >=0){\
			acc = accel[ilnr] + (accel[ilnr+1]-accel[ilnr]) *(lnr-ilnr);\
		}\
		else {\
			acc = 0;\
		}\
		acc = acc / r;\
		float ax,ay,az;\
		ax = x * acc;\
		ay = y * acc;\
		az = z * acc;\
		(TYPE##_TBP(simpar)+i)->ax += ax;\
		(TYPE##_TBP(simpar)+i)->ay += ay;\
		(TYPE##_TBP(simpar)+i)->az += az;\
	}\
}while(0)


int iflag = 1;

void readextforce_(double *,int *);

void External_Force(SimParameters *simpar, float x00,float y00,float z00){
	int i,nsize=NSIZE;
	long np;


	if(iflag){
		double fscale,pixmass,ellinkpc,mscale,ellMpc,fs,fr,pmassinMsun,meancellmass;
		pixmass = TOTMASS(simpar)/NX(simpar)/NY(simpar)/NZ(simpar);
		mscale = 1./pixmass;
		ellMpc = BOXSIZE(simpar)/NX(simpar);
		ellinkpc = BOXSIZE(simpar)/NX(simpar) * 1000.L;
		/* in M_sun */
//jhshin1
                /*extreme resolution */
                //pmassinMsun = 5245.21;
                //meancellmass = pmassinMsun / (pow(NX(simpar),3.L)/7864320.0);
		/*high resolution */
		//pmassinMsun = 41961.668; 	
		//meancellmass = pmassinMsun / (pow(NX(simpar),3.L)/983040.0);
		/* law  resolution */
		pmassinMsun = 335693.34;      
		meancellmass = pmassinMsun / (pow(NX(simpar),3.L)/122880.0);
		SPH_INITMASS(simpar) = pmassinMsun;
//jhshin2
		fs  = P_G*meancellmass * onesolarmass/pow(ellinkpc*kpcincgs, 2.L);
		fr =  1.E10L / kpcincgs; /* coefficient in (cm/sec)^2/cm unit */
		fscale = fr/fs;
	  	if(MYID(simpar) == 0) readextforce_(accel,&nsize);
		MPI_Bcast(accel,nsize,MPI_DOUBLE,0,MPI_COMM_WORLD);
		for(i=0;i<NSIZE;i++) accel[i] = accel[i]*fscale/NX(simpar)/NX(simpar);
#ifdef DEBUG
		if(MYID(simpar)==0){
			float r;
			double rscale,ellinkpc,ellinMpc;
			ellinMpc = BOXSIZE(simpar)/NX(simpar);
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
		ellinMpc = BOXSIZE(simpar)/NX(simpar);
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

#ifndef GOTPM
	GetExtAccel(DM,x00,y00,z00);
	GetExtAccel(SPH,x00,y00,z00);
	GetExtAccel(STAR,x00,y00,z00);
	GetExtAccel(AGN,x00,y00,z00);
#endif
}


