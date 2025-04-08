#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
#include "voro.h"
#include "nnost.h"
#include "exam.h"
#include "exam2d.h"
#include "gl2d.h"
//#include "DD2d.h"

#define Nneigh 50
Voro2D_point neigh[Nneigh];


void gl2d_outdata(SimParameters *simpar, int nstep, postype t, postype dt){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int myid, nid;
	myid = MYID(simpar);
	nid = NID(simpar);
	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);


	char outfile[190];
	sprintf(outfile,"glout.%.6d.dat",nstep);
	int i;
	int snp=0;
	FILE *wp;
	int tnp = VORO_TNP(simpar);
	for(i=0;i<nid;i++){
		if(myid ==i){
			if(myid==0){
				wp = fopen(outfile,"w"); 
				fwrite(&tnp, sizeof(int),1,wp);
			}
			else wp = fopen(outfile,"a");
			fwrite(bp, sizeof(treevorork4particletype), np, wp);
			fclose(wp);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if(myid==0){
		wp = fopen(outfile,"a");
		fwrite(&t, sizeof(postype), 1, wp);
		fwrite(&nx, sizeof(int), 1, wp);
		fwrite(&ny, sizeof(int), 1, wp);
		fwrite(&(SIMBOX(simpar).x.max), sizeof(postype), 1, wp);
		fwrite(&(SIMBOX(simpar).y.max), sizeof(postype), 1, wp);
		fwrite(&GAS_AlphaVis(simpar), sizeof(float), 1, wp);
		fwrite(&GAS_BetaVis(simpar), sizeof(float), 1, wp);
		fclose(wp);
	}
	/*
	MPI_Barrier(MPI_COMM_WORLD);
	DEBUGPRINT("P%d is now exiting the dump\n", myid);
	*/
}


void gl2d_readdata(SimParameters *simpar, postype *t, int nstep){
	char infile[190];
	sprintf(infile,"glout.%.6d.dat", nstep);
	int myid = MYID(simpar);
	if(myid ==0){
		printf("P0: is now reading starting file: %s\n", infile);
		FILE *fp = fopen(infile,"r");
		int np;
		fread(&np, sizeof(int), 1, fp);
		VORO_TNP(simpar) = VORO_NP(simpar) = np;
		VORORK4_TBP(simpar) = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*np);
		fread(VORORK4_TBP(simpar),sizeof(treevorork4particletype), np,fp);
		fread(t,sizeof(postype), 1,fp);
		fclose(fp);
	}
	else {
		VORO_NP(simpar) = 0;
		VORORK4_TBP(simpar) = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*100);;

	}
	DEBUGPRINT("P%d: after reading iput data\n", myid);
	migrateTreeVoroParticles(simpar);
	MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Barrier(MPI_COMM(simpar));
	BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
}


int Make2DGlass(SimParameters *simpar, int icont){
	postype t,dt;
	int np;
	treevorork4particletype *bp;
	int icount = 0;;
	int iflag,jflag;


	{
		GridInfo *grid = &(GRIDINFO(simpar));
		PosNX(grid) = NX(simpar);
		PosNY(grid) = NY(simpar);
		PosNZ(grid) = NZ(simpar);
		PosNXNY(grid) = NX(simpar)* NY(simpar);
	}
	/*
	void StartGL2DRkSDD(SimParameters *);
	StartGL2DRkSDD(simpar); 
	*/
	startRkSDD2D(simpar,GL2D);


	DEBUGPRINT("P%d has volume: rmin= %g %g rmax= %g %g\n", MYID(simpar),
			SIM_LXMIN(simpar,vorork4), SIM_LYMIN(simpar,vorork4),
			SIM_LXMAX(simpar,vorork4), SIM_LYMAX(simpar,vorork4));
	GL2D_XMIN(simpar) = SIM_LXMIN(simpar,vorork4);
	GL2D_YMIN(simpar) = SIM_LYMIN(simpar,vorork4);
	GL2D_XMAX(simpar) = SIM_LXMAX(simpar,vorork4);
	GL2D_YMAX(simpar) = SIM_LYMAX(simpar,vorork4);


	VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*100);




	if(icont ==0) {
		t = 0;
		bp = gl2d_mkinitial(simpar, &np);
		printf("made initial conditions\n");
	}
	else{
		postype dt;
		int nstep = icont;
		gl2d_readdata(simpar, &t, nstep);
		icount = nstep+1;
	}

	dt=0;
	do {
//		MkLinkedList(bp,np);
//
		iflag = t * 1.;
		jflag = (t-dt) * 1;
		if(iflag != jflag || icount%50 == 0)
		{
			gl2d_outdata(simpar,  icount, t, dt);
//			gl2d_makemap(simpar, icount);
			void gl2d_maketscmap(SimParameters *, int);
			gl2d_maketscmap(simpar, icount);
		}

//		dt = gl2d_vph2D_rk4(simpar);
		dt =  exam2d_vph_rk4(simpar, paddingTreeVoroParticles,gl2d_w2Measure,searchCellNeighbors2D, findCellBP2D);

		t += dt;
		if(MYID(simpar)==0){
			printf("Time is %g with icount = %d\n",t, icount);
			fflush(stdout);
		}
		DEBUGPRINT("P%d is now at t = %g with dt= %g\n", MYID(simpar), t,dt);



		icount ++;

	} while(icount<2000);
	return 1;
}
