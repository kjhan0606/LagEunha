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
#include "kp.h"

#define Nneigh 50
Voro2D_point neigh[Nneigh];

DefineProtoTypeFunctions2D(vorork4);
DefineProtoTypeFunctions2D(treevorork4);

void kp_outdata(SimParameters *simpar, int nstep, postype t, postype dt){
    treevorork4particletype *bp = VORORK4_TBP(simpar);
    int np = VORO_NP(simpar);
    int myid, nid;
    myid = MYID(simpar);
    nid = NID(simpar);
    int nx,ny;
    nx = NX(simpar);
    ny = NY(simpar);

    char outfile[190];
    sprintf(outfile,"kpout.%.6d.dat",nstep);
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
}

void kp_readdata(SimParameters *simpar, postype *t, postype *dt, int nstep){
    float dtold;
    char infile[190];
    sprintf(infile,"kpout.%.6d.dat", nstep);
    int myid = MYID(simpar);
    if(myid ==0){
        FILE *fp = fopen(infile,"r");
        int np;
        fread(&np, sizeof(int), 1, fp);
        VORO_TNP(simpar) = VORO_NP(simpar) = np;
        VORORK4_TBP(simpar) = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*np);
        fread(VORORK4_TBP(simpar),sizeof(treevorork4particletype), np,fp);
        fread(t,sizeof(postype), 1,fp);
        fread(dt,sizeof(postype), 1,fp);
        fread(&dtold,sizeof(float), 1,fp);
        fclose(fp);
    }
    else {
        VORO_NP(simpar) = 0;
        VORORK4_TBP(simpar) = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*100);;

    }
    migrateTreeVorork4Particles(simpar);
    MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
    MPI_Bcast(dt, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
    MPI_Bcast(&dtold, 1, MPI_FLOAT,0, MPI_COMM(simpar));
    GAS_dtold(simpar) = dtold;
    MPI_Barrier(MPI_COMM(simpar));
    BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
}

int kp_evolBP(treevorork4particletype *bp,postype Lx, postype Ly){
	return 1;
}

int RunKepler(SimParameters *simpar, int icont){
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
	startRkSDD2D(simpar,KP); /* This makes all the ddinfo including pivot */


	printf("P%d has volume: rmin= %g %g rmax= %g %g\n", MYID(simpar),
			SIM_LXMIN(simpar,vorork4), SIM_LYMIN(simpar,vorork4),
			SIM_LXMAX(simpar,vorork4), SIM_LYMAX(simpar,vorork4));
	KP_XMIN(simpar) = SIM_LXMIN(simpar,vorork4);
	KP_YMIN(simpar) = SIM_LYMIN(simpar,vorork4);
	KP_XMAX(simpar) = SIM_LXMAX(simpar,vorork4);
	KP_YMAX(simpar) = SIM_LYMAX(simpar,vorork4);

	VORORK4_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*100);

	if(icont ==0) {
        t = 0;
        dt = 1.e-5;
        bp = kp_mkinitial(simpar, &np);
        printf("made initial conditions\n");
        GAS_dtold(simpar) = dt;
    }
    else{
        int nstep = icont;
        kp_readdata(simpar, &t, &dt,nstep);
        icount = nstep;
    }

    do {
        if(GAS_EVOLMETHOD(simpar) == 1) dt = exam2d_vph_rk4_int(simpar,
                paddingTreeVorork4Particles,kp_w2Measure2D,
                searchCellRk4Neighbors2D,findCellRk4BP2D,
                kp_evolBP,
				mkLinkedList2D_oExam
                );
        else dt = exam2d_vph_int(simpar,
                paddingTreeVorork4Particles,kp_w2Measure2D,
                searchCellRk4Neighbors2D,findCellRk4BP2D,
                kp_evolBP,
				mkLinkedList2D_oExam
                );
        t += dt;
        GAS_dtold(simpar) = dt;

        icount ++;

        iflag = t * 10.;
        jflag = (t-dt) * 10.;


        if(iflag != jflag || icount == 146)
        {
            kp_outdata(simpar,  icount, t, dt);
            kp_makemap(simpar, icount);
        }

        if(MYID(simpar)==0){
            DEBUGPRINT("Time= %g & icount= %d with dTime= %g\n", t,icount, dt);
            fflush(stdout);
        }
    } while(t<121.);
    return 1;
}
