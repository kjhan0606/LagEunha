#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>



#define Mpccgs 3.085677488036136E24L
#define onesolarmass 1.989E33L

#include "mpi.h"

#include "eunha.h"


void readstarparticle_(double *, double *, double *, double *, double *, double *, double*, long *);
void readsphparticle_(double *, double *, double *, double *, double *, double *, double*, double *, long *);


double mscale,rscale,vscale,fscale,escale;
double ellinMpc,ellinkpc,pixmass,shift,pixsize;


void jhsread(SimParameters *simpar){
	double *mass,*x,*y,*z,*vx,*vy,*vz,*entropy;
	long i=0,np;
//jhshin1
        /* extreme resolution */
        //DM_TNP(simpar) = 7077888;
        //SPH_TNP(simpar) = 786432;
	/* high resolution */
	//DM_TNP(simpar) = 884736;
        //SPH_TNP(simpar) = 98304;
	/* law  resolution */
	DM_TNP(simpar) = 110592;
	SPH_TNP(simpar) = 12288;
//jhshin2
	DM_BP(simpar) = (dmparticletype *)malloc(sizeof(dmparticletype)*DM_TNP(simpar));
	SPH_BP(simpar) = (sphparticletype *)malloc(sizeof(sphparticletype)*SPH_TNP(simpar));

	ellinMpc = BOXSIZE(simpar)/NX(simpar);
	ellinkpc = ellinMpc*1000;
	pixmass = TOTMASS(simpar)/NX(simpar)/NY(simpar)/NZ(simpar);
	mscale = 1./pixmass;
	rscale = 1./ellinkpc;
	vscale = 1./(ellinMpc * Mpccgs/FFTIME(simpar)*NX(simpar) /1.e5);

	pixsize = ellinMpc * Mpccgs;
	escale = 1./pow(pixmass*onesolarmass/pixsize/pixsize/pixsize,1-GAS_GAMMA(simpar));
	fscale = 1./(pixmass/ellinkpc/ellinkpc);

	shift = NX(simpar)/2;

//jhshin1
        float rhos2rhor = TOTMASS(simpar)*1.989E33L/pow(BOXSIZE(simpar)*3.085677488036136E24L,3.L);
        if(MYID(simpar)==0){
                printf("------scaling---------------\n");
                printf("mscale    = %g \n",mscale);
                printf("rscale    = %g \n",rscale);
                printf("vscale    = %g \n",vscale);
                printf("rhos2rhor = %g \n",rhos2rhor);
                printf("----------------------------\n");
        }
//jhshin2

	mass = (double *)malloc(sizeof(double)*DM_TNP(simpar));
	x = (double *)malloc(sizeof(double)*DM_TNP(simpar));
	y = (double *)malloc(sizeof(double)*DM_TNP(simpar));
	z = (double *)malloc(sizeof(double)*DM_TNP(simpar));
	vx = (double *)malloc(sizeof(double)*DM_TNP(simpar));
	vy = (double *)malloc(sizeof(double)*DM_TNP(simpar));
	vz = (double *)malloc(sizeof(double)*DM_TNP(simpar));
	entropy = (double *)malloc(sizeof(double)*DM_TNP(simpar));
	if(MYID(simpar)==0) readstarparticle_(mass,x,y,z,vx,vy,vz, &(DM_TNP(simpar)));
	MPI_Bcast(mass,DM_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(x,DM_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(y,DM_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(z,DM_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(vx,DM_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(vy,DM_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(vz,DM_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	/*
	if(MYID(simpar)==0){
		for(i=0;i<DM_TNP(simpar);i++){
			float r= sqrt(x[i]*x[i]+y[i]*y[i] + z[i]*z[i]);
			if(r > 30){
				printf("%g : %g %g %g : %g %g %g\n",r,x[i],y[i],z[i],vx[i],vy[i],vz[i]);
			}
		}
	}
	MPI_Finalize();exit(99);
	*/
	np = 0;
	float nzwidth = (float) NZ(simpar)/(float) NID(simpar);
	for(i=0;i<DM_TNP(simpar);i++){
		y[i] = y[i]*rscale +shift;
		int ii = y[i]/nzwidth;
		if(ii == MYID(simpar)){
			mass[i] = mass[i]*mscale;
			x[i] = x[i]*rscale +shift;
			z[i] = z[i]*rscale +shift;
			vx[i] = vx[i]*vscale;
			vy[i] = vy[i]*vscale;
			vz[i] = vz[i]*vscale;

			TDM_MASS(simpar,np) = mass[i];
			(DM_BP(simpar)+np)->x = x[i];
			(DM_BP(simpar)+np)->y = z[i];
			(DM_BP(simpar)+np)->z = y[i];
			(DM_BP(simpar)+np)->vx = vx[i];
			(DM_BP(simpar)+np)->vy = vz[i];
			(DM_BP(simpar)+np)->vz = vy[i];
			CHANGEINDX( DM_BP(simpar) + np, i);
			np++;
		}
	}
	DM_NP(simpar) = np;
	DM_BP(simpar) = (dmparticletype*)realloc(DM_BP(simpar),sizeof(dmparticletype)*DM_NP(simpar));

	if(np>0){
		dmparticletype *bp = DM_BP(simpar);
		printf("P%d has 1st disk particle %g %g %g %g %g %g %g\n",
				MYID(simpar),DM_MASS(simpar,0),bp->x,bp->y,bp->z,bp->vx,bp->vy,bp->vz);
	}


	if(MYID(simpar) == 0) readsphparticle_(mass,x,y,z,vx,vy,vz,entropy,&(SPH_TNP(simpar)));
	MPI_Bcast(mass,SPH_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(x,SPH_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(y,SPH_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(z,SPH_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(vx,SPH_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(vy,SPH_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(vz,SPH_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(entropy,SPH_TNP(simpar),MPI_DOUBLE,0,MPI_COMM_WORLD);
	np = 0;
	nzwidth = (float) NZ(simpar)/(float) NID(simpar);
//jhshin1
        GAS_INITMETAL(simpar) = 0.02;
//jhshin2
	for(i=0;i<SPH_TNP(simpar);i++){
		y[i] = y[i]*rscale +shift;
		int ii = y[i]/nzwidth;
		if(ii == MYID(simpar)){
			mass[i] = mass[i]*mscale;
			x[i] = x[i]*rscale +shift;
			z[i] = z[i]*rscale +shift;
			vx[i] = vx[i]*vscale;
			vy[i] = vy[i]*vscale;
			vz[i] = vz[i]*vscale;

			(SPH_BP(simpar)+np)->mass = mass[i];
			(SPH_BP(simpar)+np)->x = x[i];
			(SPH_BP(simpar)+np)->y = z[i];
			(SPH_BP(simpar)+np)->z = y[i];
			(SPH_BP(simpar)+np)->vx = vx[i];
			(SPH_BP(simpar)+np)->vy = vz[i];
			(SPH_BP(simpar)+np)->vz = vy[i];
			(SPH_BP(simpar)+np)->Entropy = entropy[i] * escale;
			(SPH_BP(simpar)+np)->nstar = 0;
			(SPH_BP(simpar)+np)->metallicity = GAS_INITMETAL(simpar);
			CHANGEINDX( SPH_BP(simpar) + np, i);
			np++;
		}
	}
	if(np>0){
		sphparticletype *bp = SPH_BP(simpar);
		printf("P%d has 1st sph particle %g %g %g %g %g %g %g : %g escale = %g\n",
				MYID(simpar),bp->mass,bp->x,bp->y,bp->z,bp->vx,bp->vy,bp->vz,bp->Entropy,escale);
	}
	SPH_NP(simpar) = np;
	SPH_BP(simpar) = (sphparticletype*)realloc(SPH_BP(simpar),sizeof(sphparticletype)*SPH_NP(simpar));
	free(entropy); free(vz);free(vy);free(vx);free(z);free(y);free(x);free(mass);
}
