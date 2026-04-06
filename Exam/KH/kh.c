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
static Voro2D_point neigh[Nneigh];


void kh_outdata(SimParameters *simpar, int nstep, postype t, postype dt){
	int np = VORO_NP(simpar);
	int myid, nid;
	myid = MYID(simpar);
	nid = NID(simpar);
	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);

	char outfile[190];
	sprintf(outfile,"khout.%.6d.dat",nstep);
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
			int j;
			for(j=0;j<np;j++)
				fwrite(bp_raw + j*p_size, sizeof(treevorork4particletype), 1, wp);
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
	int av_mode = GAS_AVMODE(simpar);
	size_t p_size = (av_mode >= 1) ? sizeof(treevorostressrk4particletype) : sizeof(treevorork4particletype);
	if(myid ==0){
		printf("P0: is now reading starting file: %s\n", infile);
		FILE *fp = fopen(infile,"r");
		int np;
		fread(&np, sizeof(int), 1, fp);
		VORO_TNP(simpar) = VORO_NP(simpar) = np;
		if(av_mode >= 1){
			size_t common_size = offsetof(treevorork4particletype, bp);
			treevorork4particletype *tmp_buf = (treevorork4particletype*)my_malloc(sizeof(treevorork4particletype)*np);
			fread(tmp_buf, sizeof(treevorork4particletype), np, fp);
			VORORK4_TBP(simpar) = (treevorork4particletype*)my_malloc(p_size*np);
			treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
			int ii;
			for(ii=0;ii<np;ii++){
				memset(&sbp[ii], 0, sizeof(treevorostressrk4particletype));
				memcpy(&sbp[ii], &tmp_buf[ii], common_size);
			}
			my_free(tmp_buf);
		} else {
			VORORK4_TBP(simpar) = (treevorork4particletype*)my_malloc(p_size*np);
			fread(VORORK4_TBP(simpar), sizeof(treevorork4particletype), np, fp);
		}
		fread(t,sizeof(postype), 1,fp);
		fread(dt,sizeof(postype), 1,fp);
		fread(&GAS_dtold(simpar),sizeof(float), 1,fp);
		fclose(fp);
	}
	else {
		VORO_NP(simpar) = 0;
		VORORK4_TBP(simpar) = (treevorork4particletype*)my_malloc(p_size*100);
	}

	DEBUGPRINT("P%d: after reading iput data\n", myid);
	migrateTreeVorork4Particles(simpar);
	MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(dt, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(&GAS_dtold(simpar), 1, MPI_FLOAT,0, MPI_COMM(simpar));
	MPI_Barrier(MPI_COMM(simpar));
	BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
}

int kh_evolBP(treevorork4particletype *bp,postype Lx, postype Ly){
	return 1;
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

	int av_mode = GAS_AVMODE(simpar);
	if(av_mode >= 1){
		int ndd = NDDINFO(simpar);
		int k;
		size_t new_nsize = sizeof(treevorostressrk4particletype);
		for(k=0;k<ndd;k++){
			DoDeInfo *di = &TVORORK4_DDINFO(simpar)[k];
			size_t old_nsize = di->n_size;
			di->n_size = new_nsize;
			if(di->npivot > 0 && new_nsize > old_nsize){
				char *new_pivot = (char*)my_malloc(di->npivot * new_nsize);
				memset(new_pivot, 0, di->npivot * new_nsize);
				char *old_pivot = (char*)di->pivot;
				int j;
				for(j=0; j<di->npivot; j++)
					memcpy(new_pivot + j*new_nsize, old_pivot + j*old_nsize, old_nsize);
				my_free(old_pivot);
				di->pivot = new_pivot;
			}
		}
		if(MYID(simpar)==0)
			DEBUGPRINT("KH: using treevorostressrk4particletype for av_mode=%d\n", av_mode);
	}

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
		postype Gamma = GAS_GAMMA(simpar);
		GAS_invw2Scale(simpar) = 1.L / KH_Pressure(simpar);
		GAS_w2Power(simpar) = (Gamma-1)/Gamma;
	}

	do {
		/*
		if(GAS_EVOLMETHOD(simpar) == 1) dt = exam2d_vph_rk4(simpar,
				paddingTreeVorork4Particles, kh_w2Measure2D,
				searchCellRk4Neighbors2D, findCellRk4BP2D);
		else dt = exam2d_vph(simpar, kh_w2Measure2D);
		*/
		if(av_mode >= 1 && GAS_EVOLMETHOD(simpar) == 1)
			dt = exam2d_vph_rk4_int_blend(simpar,
				paddingTreeVorork4Particles,kh_w2Measure2D,
				searchCellRk4Neighbors2D,findCellRk4BP2D,
				kh_evolBP,
				mkLinkedList2D_oExam,
				periodic_postStage_blend
				);
		else if(GAS_EVOLMETHOD(simpar) == 2)
			dt = exam2d_vph_rk4_int_kNN(simpar,
				paddingTreeVorork4Particles,kh_w2Measure2D,
				kh_evolBP,
				NULL, periodic_postStage_kNN
				);
		else if(GAS_EVOLMETHOD(simpar) == 1)
			dt = exam2d_vph_rk4_int(simpar,
				paddingTreeVorork4Particles,kh_w2Measure2D,
				searchCellRk4Neighbors2D,findCellRk4BP2D,
				kh_evolBP,
				mkLinkedList2D_oExam
				);
        else dt = exam2d_vph_int(simpar,
                paddingTreeVorork4Particles,kh_w2Measure2D,
                searchCellRk4Neighbors2D,findCellRk4BP2D,
                kh_evolBP,
				mkLinkedList2D_oExam
                );
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
			printf("Time= %g  step= %d  dt= %g\n", t, icount, dt);
			fflush(stdout);
		}

	} while(t<10.);
	return 1;
}
