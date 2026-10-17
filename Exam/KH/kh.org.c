#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
#include "voro.h"
#include "nnost.h"
#include "gnnost.h"
#include "exam.h"
#include "exam2d.h"
#include "kh.h"
//#include "DD2d.h"

#define Nneigh 50
Voro2D_point neigh[Nneigh];


void kh_outdata(SimParameters *simpar, int nstep, postype t, postype dt){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
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
		fwrite(&dt, sizeof(postype), 1, wp);
		fwrite(&GAS_dtold(simpar), sizeof(float), 1, wp);
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


void kh_readdata(SimParameters *simpar, postype *t, postype *dt, int nstep){
	char infile[190];
	sprintf(infile,"khout.%.6d.dat", nstep);
	int myid = MYID(simpar);
	float dtold;
	if(myid ==0){
		printf("P0: is now reading starting file: %s\n", infile);
		FILE *fp = fopen(infile,"r");
		int np;
		fread(&np, sizeof(int), 1, fp);
		VORO_TNP(simpar) = VORO_NP(simpar) = np;
		VORORK4_TBP(simpar) = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*np);

		if(0){
			treevoroparticletype *bb = (treevoroparticletype *)malloc(sizeof(treevoroparticletype)*np);
			treevorork4particletype *aa = VORORK4_TBP(simpar);
			fread(bb, sizeof(treevoroparticletype), np, fp);
			int i;
			for(i=0;i<np;i++){
				aa[i].u4if = bb[i].u4if;
				aa[i].x = bb[i].x;
				aa[i].y = bb[i].y;
				aa[i].vx = bb[i].vx;
				aa[i].vy = bb[i].vy;
				aa[i].mass = bb[i].mass;
				aa[i].den = bb[i].den;
				aa[i].pressure = bb[i].pressure;
				aa[i].ie = bb[i].ie;
				aa[i].csound = bb[i].csound;
				aa[i].volume = bb[i].volume;
				aa[i].w2 = bb[i].w2;
				aa[i].w2old = bb[i].w2old;
				aa[i].w2ceil = bb[i].w2ceil;
			}
		}
		else {
			fread(VORORK4_TBP(simpar),sizeof(treevorork4particletype), np,fp);
		}
		fread(t,sizeof(postype), 1,fp);
		fread(dt,sizeof(postype), 1,fp);
		fread(&dtold,sizeof(float), 1,fp);

		fclose(fp);
	}
	else {
		VORO_NP(simpar) = 0;
		VORORK4_TBP(simpar) = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*100);;

	}

	DEBUGPRINT("P%d: after reading iput data\n", myid);
	migrateTreeVorork4Particles(simpar);
	MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(dt, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(&dtold, 1, MPI_FLOAT,0, MPI_COMM(simpar));
	GAS_dtold(simpar) = dtold;
	MPI_Barrier(MPI_COMM(simpar));
	BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
}


int RunKH(SimParameters *simpar, int icont){
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
	startRkSDD2D(simpar,KH);


	DEBUGPRINT("P%d has volume: rmin= %g %g rmax= %g %g\n", MYID(simpar),
			SIM_LXMIN(simpar,vorork4), SIM_LYMIN(simpar,vorork4),
			SIM_LXMAX(simpar,vorork4), SIM_LYMAX(simpar,vorork4));
	KH_XMIN(simpar) = SIM_LXMIN(simpar,vorork4);
	KH_YMIN(simpar) = SIM_LYMIN(simpar,vorork4);
	KH_XMAX(simpar) = SIM_LXMAX(simpar,vorork4);
	KH_YMAX(simpar) = SIM_LYMAX(simpar,vorork4);


	VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*100);




	if(icont ==0) {
		t = 0;
		dt=1.e-7;
		bp = kh_mkinitial(simpar, &np);
		printf("made initial conditions\n");
		GAS_dtold(simpar) = dt;
	}
	else{
		int nstep = icont;
		kh_readdata(simpar, &t, &dt, nstep);
		icount = nstep;
	}

	do {

//		if(GAS_EVOLMETHOD(simpar) == 1) dt = kh_vph2D_rk4(simpar);
//		else dt = kh_vph2D(simpar);
		if(GAS_EVOLMETHOD(simpar) == 1) dt = exam2d_vph_rk4(simpar,
				paddingTreeVorork4Particles, kh_w2Measure2D,
				searchCellRk4Neighbors2D, findCellRk4BP2D);
		else dt = exam2d_vph(simpar, kh_w2Measure2D);
		t += dt;


		GAS_dtold(simpar) = dt;


		icount ++;

		iflag = t * 10.;
		jflag = (t-dt) * 10.;

//		if(iflag != jflag || icount%100 ==0 || icount == 1 || icount == 2)
		if(iflag != jflag )
		{
			kh_outdata(simpar,  icount, t, dt);
			kh_makemap(simpar, icount);
		}

		if(MYID(simpar)==0){
			DEBUGPRINT("P%d is now at Time= %g & icount= %d with dTime= %g\n", 
					MYID(simpar), t,icount, dt);
			fflush(stdout);
		}




	} while(t<10.);
	return 1;
}
