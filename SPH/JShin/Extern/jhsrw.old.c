#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>

#include "mpi.h"

#include "pmheader.h"
#include "Memory.h"


void readstarparticle_(double *, double *, double *, double *, double *, double *, double*, long long *);
void readsphparticle_(double *, double *, double *, double *, double *, double *, double*, double *, long long *);


void jhsread(){
	double *mass,*x,*y,*z,*vx,*vy,*vz,*entropy;
	long i=0,np;

	simpar.dm.tnp = 909091;
	simpar.dm.u.tbp = (treepmparticletype *)Malloc(sizeof(treepmparticletype)*simpar.dm.tnp,PPTR(simpar.dm.u.tbp));
	simpar.sph.tnp = 90909;
	simpar.sph.u.tbp = (treesphparticletype *)Malloc(sizeof(treesphparticletype)*simpar.sph.tnp,PPTR(simpar.sph.u.tbp));


	mass = (double *)Malloc(sizeof(double)*simpar.dm.tnp,PPTR(mass));
	x = (double *)Malloc(sizeof(double)*simpar.dm.tnp,PPTR(x));
	y = (double *)Malloc(sizeof(double)*simpar.dm.tnp,PPTR(y));
	z = (double *)Malloc(sizeof(double)*simpar.dm.tnp,PPTR(z));
	vx = (double *)Malloc(sizeof(double)*simpar.dm.tnp,PPTR(vx));
	vy = (double *)Malloc(sizeof(double)*simpar.dm.tnp,PPTR(vy));
	vz = (double *)Malloc(sizeof(double)*simpar.dm.tnp,PPTR(vz));
	entropy = (double *)Malloc(sizeof(double)*simpar.dm.tnp,PPTR(entropy));
	if(simpar.myid==0) readstarparticle_(mass,x,y,z,vx,vy,vz, &(simpar.dm.tnp));
	MPI_Bcast(mass,simpar.dm.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(x,simpar.dm.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(y,simpar.dm.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(z,simpar.dm.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(vx,simpar.dm.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(vy,simpar.dm.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(vz,simpar.dm.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	np = 0;
	float nzwidth = (float) simpar.nz/(float) simpar.nid;
	for(i=0;i<simpar.dm.tnp;i++){
		int ii = y[i]/nzwidth;
		if(ii == simpar.myid){
			(simpar.dm.u.bp+np)->mass = mass[i];
			(simpar.dm.u.bp+np)->x = x[i];
			(simpar.dm.u.bp+np)->y = z[i];
			(simpar.dm.u.bp+np)->z = y[i];
			(simpar.dm.u.bp+np)->vx = vx[i];
			(simpar.dm.u.bp+np)->vy = vz[i];
			(simpar.dm.u.bp+np)->vz = vy[i];
			np++;
		}
	}
	simpar.dm.np = np;
	simpar.dm.u.tbp = (treepmparticletype*)Realloc(simpar.dm.u.tbp,sizeof(treepmparticletype)*simpar.dm.np);


	if(simpar.myid == 0) readsphparticle_(mass,x,y,z,vx,vy,vz,entropy,&(simpar.sph.tnp));
	MPI_Bcast(mass,simpar.sph.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(x,simpar.sph.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(y,simpar.sph.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(z,simpar.sph.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(vx,simpar.sph.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(vy,simpar.sph.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(vz,simpar.sph.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(entropy,simpar.sph.tnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
	np = 0;
	nzwidth = (float) simpar.nz/(float) simpar.nid;
	for(i=0;i<simpar.sph.tnp;i++){
		int ii = y[i]/nzwidth;
		if(ii == simpar.myid){
			(simpar.sph.u.bp+np)->mass = mass[i];
			(simpar.sph.u.bp+np)->x = x[i];
			(simpar.sph.u.bp+np)->y = z[i];
			(simpar.sph.u.bp+np)->z = y[i];
			(simpar.sph.u.bp+np)->vx = vx[i];
			(simpar.sph.u.bp+np)->vy = vz[i];
			(simpar.sph.u.bp+np)->vz = vy[i];
			(simpar.sph.u.bp+np)->Entropy = entropy[i];
			(simpar.sph.u.bp+np)->nstar = 0;
			(simpar.sph.u.bp+np)->metallicity = 0;
			np++;
		}
	}
	simpar.sph.np = np;
	simpar.sph.u.tbp = (treesphparticletype*)Realloc(simpar.sph.u.tbp,sizeof(treesphparticletype)*simpar.sph.np);


	Free(entropy); Free(vz);Free(vy);Free(vx);Free(z);Free(y);Free(x);Free(mass);
}
