#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
#include "voro.h"
#include "kh.h"
#include "DD2d.h"

#define Nneigh 50
Voro2D_point neigh[Nneigh];


void outdata(SimParameters *simpar, int nstep, postype t, postype dt){
	treevoroparticletype *bp = VORO_TBP(simpar);
	int np = VORO_NP(simpar);
	int myid, nid;
	myid = MYID(simpar);
	nid = NID(simpar);
	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);

	char outfile[190];
	sprintf(outfile,"khout.%.6d.dat",nstep);
	int i;
	int snp=0;
	FILE *wp;
	if(myid ==0){
		wp = fopen(outfile,"w");
		int tnp = VORO_TNP(simpar);
		fwrite(&tnp, sizeof(int),1,wp);
		fwrite(bp, sizeof(treevoroparticletype), np, wp);
		long offset = sizeof(int)+sizeof(treevoroparticletype)*VORO_TNP(simpar);
		fseek(wp, offset, SEEK_SET);
		fwrite(&t, sizeof(postype), 1, wp);
		fwrite(&dt, sizeof(postype), 1, wp);
		fwrite(&nx, sizeof(int), 1, wp);
		fwrite(&ny, sizeof(int), 1, wp);
		fwrite(&(SIMBOX(simpar).x.max), sizeof(postype), 1, wp);
		fwrite(&(SIMBOX(simpar).y.max), sizeof(postype), 1, wp);
		fwrite(&GAS_AlphaVis(simpar), sizeof(float), 1, wp);
		fwrite(&GAS_BetaVis(simpar), sizeof(float), 1, wp);
		fclose(wp);
		int rnp;
		rnp = np;
		MPI_Send(&rnp, 1, MPI_INT, myid+1, 0, MPI_COMM_WORLD);
	}
	else{
		int rnp;
		MPI_Status status;
		MPI_Recv(&snp, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD,&status);
		long offset = sizeof(int)+sizeof(treevoroparticletype)*snp;
		wp = fopen(outfile,"r+");
		fseek(wp, offset, SEEK_SET);
		fwrite(bp, sizeof(treevoroparticletype), np, wp);
		fclose(wp);
		rnp = snp + np;
		if(myid != nid-1){
			MPI_Send(&rnp, 1, MPI_INT, myid+1, 0, MPI_COMM_WORLD);
		}
	}
}


void readdata(SimParameters *simpar, postype *t, postype *dt,int nstep){
	char infile[190];
	sprintf(infile,"khout.%.6d.dat", nstep);
	int myid = MYID(simpar);
	if(myid ==0){
		FILE *fp = fopen(infile,"r");
		int np;
		fread(&np, sizeof(int), 1, fp);
		VORO_TNP(simpar) = VORO_NP(simpar) = np;
		VORO_TBP(simpar) = (treevoroparticletype*)malloc(sizeof(treevoroparticletype)*np);
		fread(VORO_TBP(simpar),sizeof(treevoroparticletype), np,fp);
		fread(t,sizeof(postype), 1,fp);
		fread(dt,sizeof(postype), 1,fp);
		fclose(fp);
	}
	else {
		VORO_NP(simpar) = 0;
		VORO_TBP(simpar) = (treevoroparticletype*)malloc(sizeof(treevoroparticletype)*100);;
	}
	KH_TreeAllParticleMigrate(simpar);
	MPI_Bcast(&VORO_TNP(simpar), 1, MPI_INT64_T,0, MPI_COMM(simpar));
	MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
}


int RunKH(SimParameters *simpar, int icont){
	postype t,dt;
	int np;
	treevoroparticletype *bp;
	int icount = 0;;
	int iflag,jflag;


	{
		GridInfo *grid = &(GRIDINFO(simpar));
		PosNX(grid) = NX(simpar);
		PosNY(grid) = NY(simpar);
		PosNZ(grid) = NZ(simpar);
		PosNXNY(grid) = NX(simpar)* NY(simpar);
	}
	void StartKHRkSDD(SimParameters *);
	StartKHRkSDD(simpar); /* This makes all the ddinfo including pivot */


	printf("P%d has volume: rmin= %g %g rmax= %g %g\n", MYID(simpar),
			SIM_LXMIN(simpar,voro), SIM_LYMIN(simpar,voro),
			SIM_LXMAX(simpar,voro), SIM_LYMAX(simpar,voro));
	KH_XMIN(simpar) = SIM_LXMIN(simpar,voro);
	KH_YMIN(simpar) = SIM_LYMIN(simpar,voro);
	KH_XMAX(simpar) = SIM_LXMAX(simpar,voro);
	KH_YMAX(simpar) = SIM_LYMAX(simpar,voro);

	VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*100);



	if(icont ==0) {
		t = 0;
		bp = mkinitial(simpar, &np);
		printf("made initial conditions\n");
	}
	else{
		postype dt;
		int nstep = icont;
		readdata(simpar, &t,&dt, nstep);
		icount = nstep+1;
	}

	do {
//		MkLinkedList(bp,np);
		postype dt = vph2D(simpar);
		t += dt;
		if(MYID(simpar)==0){
			printf("Time is %g with icount = %d\n",t, icount);
			fflush(stdout);
		}


		iflag = t * 10.;
		jflag = (t-dt) * 10.;


//		if(icount%50 ==0) 
		if(iflag != jflag)		
		{
			outdata(simpar,  icount, t, dt);
			makemap(simpar, icount);
		}


		icount ++;

	} while(t<10.);
	return 1;
}
