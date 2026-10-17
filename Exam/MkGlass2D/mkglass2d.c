#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<mpi.h>
#include<omp.h>
#include<math.h>
#include "eunha.h"
#include "voro.h"
#include "voro_eunha.h"
#include "nnost.h"
#include "exam.h"
#include "exam2d.h"
#include "gl2d.h"
//#include "DD2d.h"

/*
#define Nneigh 50
Voro2D_point neigh[Nneigh];
*/

int mkglass2d_evol_BP(treevorork4particletype *bp,postype Lx, postype Ly){
    if(IS_FLAG_ON(bp, Wallflag)) return 0;
    else return 1;
}

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
	migrateTreeVorork4Particles(simpar);
	MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Barrier(MPI_COMM(simpar));
	BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
}


int Make2DGlass(SimParameters *simpar, int icont){
	postype t,dt;
	int np;
	treevorork4particletype *bp;
	int iflag,jflag;
	int icount = 0;;



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

	double meanvol = NX(simpar)*NY(simpar)/(GL2D_XMAX(simpar)-GL2D_XMIN(simpar))/
		(GL2D_YMAX(simpar)-GL2D_YMIN(simpar));

	DEBUGPRINT("P%d has volume: rmin= %g %g rmax= %g %g\n", MYID(simpar),
			SIM_LXMIN(simpar,vorork4), SIM_LYMIN(simpar,vorork4),
			SIM_LXMAX(simpar,vorork4), SIM_LYMAX(simpar,vorork4));
	GL2D_XMIN(simpar) = SIM_LXMIN(simpar,vorork4);
	GL2D_YMIN(simpar) = SIM_LYMIN(simpar,vorork4);
	GL2D_XMAX(simpar) = SIM_LXMAX(simpar,vorork4);
	GL2D_YMAX(simpar) = SIM_LYMAX(simpar,vorork4);
	VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*100);// initialize the mem space



	if(icont ==0) {
		t = 0;
		dt = 1.e-5;
		bp = gl2d_mkinitial(simpar, &np);
		printf("made initial conditions\n");
		GAS_dtold(simpar) = dt;
	}
	else{
		int nstep = icont;
		readRk4Data2D(simpar, &t, &dt, nstep);
		icount = nstep;
		postype dmean = SIMBOX(simpar).x.max/NX(simpar);
        GAS_dMean(simpar) = dmean;
		GAS_invw2Scale(simpar) = 1.L;
	}

	do {
		if(GAS_EVOLMETHOD(simpar) == 1) dt = exam2d_vph_rk4_int(simpar,
                paddingTreeVorork4Particles,rt_w2Measure2D, // w2Measure2D is obsolete.
                searchCellRk4Neighbors2D,findCellRk4BP2D,
                mkglass2d_evol_BP,
                mkLinkedList2D_oExam
                );
        else dt = exam2d_vph_int_rt(simpar,
                paddingTreeVorork4Particles,rt_w2Measure2D,
                searchCellRk4Neighbors2D,findCellRk4BP2D,
                mkglass2d_evol_BP,
                mkLinkedList2D_oExam
                );
		// inserted for mkglass2d
		{
	   		treevorork4particletype *bp = VORORK4_TBP(simpar);
			int i;
			for(i=0;i<VORO_NP(simpar);i++){
				bp[i].vx = bp[i].vy = 0;
				/*
				bp[i].pressure = pow(bp[i].den, Gamma);
				bp[i].ie  = bp[i].pressure *bp[i].volume/(Gamma-1);
				bp[i].csound  = sqrt(Gamma*bp[i].pressure/bp[i].den);
				*/
			}
		}
        t += dt;
        GAS_dtold(simpar) = dt;

        icount ++;
		iflag = t * 10.;
        jflag = (t-dt) * 10.;


        MPI_Barrier(MPI_COMM_WORLD);
        if(iflag != jflag || (icount %50) ==0) {
            dumpRk4Data2D(simpar,  icount, t, dt);
            MPI_Barrier(MPI_COMM_WORLD);
            if(MYID(simpar) ==0) {DEBUGPRINT0("passing saving data\n");fflush(stdout);}
			void gl2d_maketscmap(SimParameters *, int);
            gl2d_maketscmap(simpar, icount);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if(MYID(simpar)==0){
            DEBUGPRINT("Time= %g & icount= %d with dTime= %g malloc count= %d\n",
                    t,icount, dt, get_malloc_call_count()); fflush(stdout);
        }
    } while(t<10.);
    return 1;
}

/*

		iflag = t * 1.;
		jflag = (t-dt) * 1;
		if(iflag != jflag || icount%50 ==0)
		{
			dumpRk4Data2D(simpar,  icount, t, dt);
			void gl2d_maketscmap(SimParameters *, int);
			gl2d_maketscmap(simpar, icount);
		}
		if(dt !=0) GAS_dtold(simpar) = dt;
		else GAS_dtold(simpar) = 1.e-6;

		if(GAS_EVOLMETHOD(simpar) == 1) {
			dt = gl2d_vph2D_rk4(simpar);
		}
		else {
			dt = gl2d_vph2D(simpar);
		}

		t += dt;
		if(MYID(simpar)==0){
			printf("Time is %g with icount = %d\n",t, icount);
			fflush(stdout);
		}
		DEBUGPRINT("P%d is now at t = %g with dt= %g\n", MYID(simpar), t,dt);

		icount ++;

	} while(icount<3000);
	return 1;
}
*/
