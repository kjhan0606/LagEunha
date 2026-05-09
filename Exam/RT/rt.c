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
	int np = VORO_NP(simpar);
	int myid, nid;
	myid = MYID(simpar);
	nid = NID(simpar);
	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);
	/* Actual byte stride per particle (stress-sized when av_mode >= 1) */
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);

	char outfile[190];
	sprintf(outfile,"rtout.%.6d.dat",nstep);
	int i;
	FILE *wp;
	int tnp = VORO_TNP(simpar);
	for(i=0;i<nid;i++){
		if(myid ==i){
			if(myid==0){
				wp = fopen(outfile,"w");
				if(!wp) { fprintf(stderr,"P%d: cannot open %s\n",myid,outfile); MPI_Abort(MPI_COMM(simpar),1); }
				fwrite(&tnp, sizeof(int),1,wp);
			}
			else { wp = fopen(outfile,"a"); if(!wp) { fprintf(stderr,"P%d: cannot open %s\n",myid,outfile); MPI_Abort(MPI_COMM(simpar),1); } }
			/* Write common-prefix per particle using correct stride.
			   Always writes sizeof(treevorork4particletype) bytes per particle
			   for backward-compatible checkpoint format. */
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
		if(!wp) { fprintf(stderr,"P0: cannot open %s for append\n",outfile); MPI_Abort(MPI_COMM(simpar),1); }
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
	int av_mode = GAS_AVMODE(simpar);
	size_t p_size = (av_mode >= 1) ? sizeof(treevorostressrk4particletype) : sizeof(treevorork4particletype);
	if(myid ==0){
		FILE *fp = fopen(infile,"r");
		if(!fp) { fprintf(stderr,"P0: cannot open %s for read\n",infile); MPI_Abort(MPI_COMM(simpar),1); }
		int np;
		fread(&np, sizeof(int), 1, fp);
		VORO_TNP(simpar) = VORO_NP(simpar) = np;
		if(av_mode >= 1){
			/* Checkpoint uses treevorork4particletype layout; expand to stress-sized.
			   Copy only the common prefix (up to the bp pointer field) to avoid
			   overwriting the stress region with stale pointer bytes. */
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
	int av_mode = GAS_AVMODE(simpar);


	{
		GridInfo *grid = &(GRIDINFO(simpar));
		PosNX(grid) = NX(simpar);
		PosNY(grid) = NY(simpar);
		PosNZ(grid) = NZ(simpar);
		PosNXNY(grid) = NX(simpar)* NY(simpar);
	}
	startRkSDD2D(simpar,RT); /* This makes all the ddinfo including pivot */

	/* Override n_size for stress particles when using NS-stress/blend */
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
			DEBUGPRINT("RT: using treevorostressrk4particletype for av_mode=%d\n", av_mode);
	}

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
		GAS_invw2Scale(simpar) = 1.L / RT_Phalf(simpar);
		GAS_w2Power(simpar) = (GAS_GAMMA(simpar)-1)/GAS_GAMMA(simpar);
	}

	do {

		if(av_mode >= 1 && GAS_EVOLMETHOD(simpar) == 1)
			dt = exam2d_vph_rk4_int_blend(simpar,
				paddingTreeVorork4Particles,rt_w2Measure2D,
				searchCellRk4Neighbors2D,findCellRk4BP2D,
				rt_evolBP,
				mkLinkedList2D_rt,
				periodic_postStage_blend
				);
		else if(av_mode >= 1 && GAS_EVOLMETHOD(simpar) == 3)
			dt = exam2d_vph_kdk_int_blend(simpar,
				paddingTreeVorork4Particles,rt_w2Measure2D,
				searchCellRk4Neighbors2D,findCellRk4BP2D,
				rt_evolBP,
				mkLinkedList2D_rt,
				periodic_postStage_blend
				);
		else if(GAS_EVOLMETHOD(simpar) == 1)
			dt = exam2d_vph_rk4_int_rt(simpar,
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
