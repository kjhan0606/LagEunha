#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
#include "voro.h"
#include "nnost.h"
//#include "gnnost.h"
#include "exam.h"
#include "exam2d.h"
#include "rt.h"
//#include "DD2d.h"

#define Nneigh 50
static Voro2D_point neigh[Nneigh];

extern int  malloc_calls;

void rt_outdata(SimParameters *simpar, int nstep, postype t, postype dt){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int myid, nid;
	myid = MYID(simpar);
	nid = NID(simpar);
	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);

	char outfile[190];
	sprintf(outfile,"rtout.%.6d.dat",nstep);
	int i;
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
			fflush(wp);
			fclose(wp);
		}
		MPI_Barrier(MPI_COMM(simpar));
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
		fflush(wp);
		fclose(wp);
	}
}


void rt_readdata(SimParameters *simpar, postype *t, postype *dt, int nstep){
	float dtold;
	char infile[190];
	sprintf(infile,"rtout.%.6d.dat", nstep);
	int myid = MYID(simpar);
	if(myid ==0){
		FILE *fp = fopen(infile,"r");
		int np;
		fread(&np, sizeof(int), 1, fp);
		VORO_TNP(simpar) = VORO_NP(simpar) = np;
		VORORK4_TBP(simpar) = (treevorork4particletype*)my_malloc(sizeof(treevorork4particletype)*np);
		fread(VORORK4_TBP(simpar),sizeof(treevorork4particletype), np,fp);
		fread(t,sizeof(postype), 1,fp);
		fread(dt,sizeof(postype), 1,fp);
		fread(&GAS_dtold(simpar),sizeof(float), 1,fp);
		fclose(fp);
	}
	else {
		VORO_NP(simpar) = 0;
		VORORK4_TBP(simpar) = (treevorork4particletype*)my_malloc(sizeof(treevorork4particletype)*100);;

	}
	migrateTreeVorork4Particles(simpar);
	MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(dt, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(&GAS_dtold(simpar), 1, MPI_FLOAT,0, MPI_COMM(simpar));
	MPI_Barrier(MPI_COMM(simpar));
	BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
}


int rt_evolBP(treevorork4particletype *bp,postype Lx, postype Ly){
	/*
	if(bp->y>0.05*Ly && bp->y < 0.95*Ly) return 1;
	else return 0;
	*/
	if(IS_FLAG_ON(bp, Wallflag)) return 0;
	else return 1;
}

int RunRT(SimParameters *simpar, int icont){
	postype t;
	postype dt;
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
	startRkSDD2D(simpar,RT); /* This makes all the ddinfo including pivot */
	KH_XMIN(simpar) = SIM_LXMIN(simpar,vorork4);
	KH_YMIN(simpar) = SIM_LYMIN(simpar,vorork4);
	KH_XMAX(simpar) = SIM_LXMAX(simpar,vorork4);
	KH_YMAX(simpar) = SIM_LYMAX(simpar,vorork4);


	VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*100);


	if(icont ==0) {
		t = 0;
		dt = 1.e-5;
		bp = rt_mkinitial(simpar, &np);
		printf("made initial conditions\n");
		GAS_dtold(simpar) = dt;
	}
	else{
		int nstep = icont;
		rt_readdata(simpar, &t, &dt,nstep);
		icount = nstep;
		postype dmean = SIMBOX(simpar).x.max/NX(simpar);
		GAS_dMean(simpar) = dmean;
		postype rho1 = RT_DEN1(simpar), rho2 = RT_DEN2(simpar);
		postype deltay = RT_Deltay(simpar);
		float Ly = SIMBOX(simpar).y.max;
		postype ycen = 0.5*Ly;
		postype getRT_Den(postype , postype , postype , postype ,postype );
//		GAS_invw2Scale(simpar) = 1.L/(RT_Phalf(simpar)/getRT_Den(rho1,rho2,deltay,ycen,ycen));
		GAS_invw2Scale(simpar) = 1.L/(RT_Phalf(simpar)*dmean*dmean);
	}

	do {

		if(GAS_EVOLMETHOD(simpar) == 1) dt = exam2d_vph_rk4_int_rt(simpar,
				paddingTreeVorork4Particles,rt_w2Measure2D,
				searchCellRk4Neighbors2D,findCellRk4BP2D,
				rt_evolBP,
				mkLinkedList2D_rt
				);
		else dt = exam2d_vph_int_rt(simpar, 
				paddingTreeVorork4Particles,rt_w2Measure2D,
				searchCellRk4Neighbors2D,findCellRk4BP2D,
				rt_evolBP,
				mkLinkedList2D_rt
				);
		t += dt;
		GAS_dtold(simpar) = dt;

		icount ++;

		iflag = t * 10.;
		jflag = (t-dt) * 10.;


		if(iflag != jflag || icount == 21045)
		{
			rt_outdata(simpar,  icount, t, dt);
			MPI_Barrier(MPI_COMM_WORLD);
			if(MYID(simpar) ==0) {DEBUGPRINT0("passing saving data\n");fflush(stdout);}
			rt_makemap(simpar, icount);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(MYID(simpar)==0){
            DEBUGPRINT("Time= %g & icount= %d with dTime= %g malloc count= %d\n", 
					t,icount, dt, get_malloc_call_count()); fflush(stdout);
        }
	} while(t<10.);
	return 1;
}
