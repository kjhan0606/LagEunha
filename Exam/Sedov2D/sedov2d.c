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
#include "sedov2d.h"

extern int malloc_calls;

void sedov2d_outdata(SimParameters *simpar, int nstep, postype t, postype dt){
	int np = VORO_NP(simpar);
	int myid = MYID(simpar);
	int nid = NID(simpar);
	int nx = NX(simpar);
	int ny = NY(simpar);
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);

	char outfile[190];
	sprintf(outfile,"sedov2dout.%.6d.dat",nstep);
	int i;
	FILE *wp;
	int tnp = VORO_TNP(simpar);
	for(i=0;i<nid;i++){
		if(myid == i){
			if(myid==0){
				wp = fopen(outfile,"w");
				if(!wp){ fprintf(stderr,"P%d: cannot open %s\n",myid,outfile); MPI_Abort(MPI_COMM(simpar),1); }
				fwrite(&tnp, sizeof(int),1,wp);
			} else {
				wp = fopen(outfile,"a");
				if(!wp){ fprintf(stderr,"P%d: cannot open %s\n",myid,outfile); MPI_Abort(MPI_COMM(simpar),1); }
			}
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
		fflush(wp);
		fclose(wp);
	}
}

void sedov2d_readdata(SimParameters *simpar, postype *t, postype *dt, int nstep){
	char infile[190];
	sprintf(infile,"sedov2dout.%.6d.dat", nstep);
	int myid = MYID(simpar);
	int av_mode = GAS_AVMODE(simpar);
	size_t p_size = (av_mode >= 1) ? sizeof(treevorostressrk4particletype) : sizeof(treevorork4particletype);
	if(myid == 0){
		FILE *fp = fopen(infile,"r");
		if(!fp){ fprintf(stderr,"P0: cannot open %s\n",infile); MPI_Abort(MPI_COMM(simpar),1); }
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
		fread(t,sizeof(postype),1,fp);
		fread(dt,sizeof(postype),1,fp);
		fread(&GAS_dtold(simpar),sizeof(float),1,fp);
		fclose(fp);
	} else {
		VORO_NP(simpar) = 0;
		VORORK4_TBP(simpar) = (treevorork4particletype*)my_malloc(p_size*100);
	}
	migrateTreeVorork4Particles(simpar);
	MPI_Bcast(t, 1, MPI_POSTYPE, 0, MPI_COMM(simpar));
	MPI_Bcast(dt, 1, MPI_POSTYPE, 0, MPI_COMM(simpar));
	MPI_Bcast(&GAS_dtold(simpar), 1, MPI_FLOAT, 0, MPI_COMM(simpar));
	MPI_Barrier(MPI_COMM(simpar));
	BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
}

int sedov2d_evolBP(treevorork4particletype *bp, postype Lx, postype Ly){
	if(IS_FLAG_ON(bp, Wallflag)) return 0;
	else return 1;
}

postype rt_w2Measure2D(SimParameters *, postype, postype, postype);

int RunSedov2D(SimParameters *simpar, int icont){
	postype t, dt;
	int np;
	treevorork4particletype *bp;
	int icount = 0;
	int iflag, jflag;
	int av_mode = GAS_AVMODE(simpar);

	fprintf(stderr,"[Sedov2D-DBG] P%d enter RunSedov2D av_mode=%d Kappa=%g\n",
		MYID(simpar), av_mode, GAS_Kappa(simpar)); fflush(stderr);

	{
		GridInfo *grid = &(GRIDINFO(simpar));
		PosNX(grid) = NX(simpar);
		PosNY(grid) = NY(simpar);
		PosNZ(grid) = NZ(simpar);
		PosNXNY(grid) = NX(simpar)*NY(simpar);
	}
	fprintf(stderr,"[Sedov2D-DBG] P%d before startRkSDD2D xmin=%g xmax=%g ymin=%g ymax=%g\n",
		MYID(simpar), SEDOV2D_XMIN(simpar), SEDOV2D_XMAX(simpar),
		SEDOV2D_YMIN(simpar), SEDOV2D_YMAX(simpar)); fflush(stderr);
	/* Inlined startRkSDD2D macro for debugging */
	buildDDInfo2D(simpar);
	fprintf(stderr,"[Sedov2D-DBG] P%d after buildDDInfo2D\n", MYID(simpar)); fflush(stderr);
	NDDINFO(simpar) = 0;
	SIMBOX(simpar).x.min = SEDOV2D_XMIN(simpar);
	SIMBOX(simpar).x.max = SEDOV2D_XMAX(simpar);
	SIMBOX(simpar).y.min = SEDOV2D_YMIN(simpar);
	SIMBOX(simpar).y.max = SEDOV2D_YMAX(simpar);
	SIMBOX(simpar).z.min = 0;
	SIMBOX(simpar).z.max = 0;
	fprintf(stderr,"[Sedov2D-DBG] P%d simbox set, before buildSimpleRMS\n", MYID(simpar)); fflush(stderr);
	buildSimpleRMS(simpar, SIMBOX(simpar), sizeof(vorork4particletype),
		&VORORK4_DDFUNC(simpar),
		VORORK4_DDINFO(simpar), MPI_COMM(simpar));
	fprintf(stderr,"[Sedov2D-DBG] P%d after buildSimpleRMS NDDINFO=%d\n",
		MYID(simpar), NDDINFO(simpar)); fflush(stderr);
	extractLocalDomainVolume(simpar);
	fprintf(stderr,"[Sedov2D-DBG] P%d after extractLocalDomainVolume\n", MYID(simpar)); fflush(stderr);
	copyDDFuncDDInfoFromTYPE1(simpar, VORORK4, vorork4, TVORORK4, treevorork4);
	fprintf(stderr,"[Sedov2D-DBG] P%d after copyDDFunc DONE startRkSDD2D simbox.x=(%g,%g)\n",
		MYID(simpar), SIMBOX(simpar).x.min, SIMBOX(simpar).x.max); fflush(stderr);

	if(av_mode >= 1){
		int ndd = NDDINFO(simpar);
		size_t new_nsize = sizeof(treevorostressrk4particletype);
		for(int k=0;k<ndd;k++){
			DoDeInfo *di = &TVORORK4_DDINFO(simpar)[k];
			size_t old_nsize = di->n_size;
			di->n_size = new_nsize;
			if(di->npivot > 0 && new_nsize > old_nsize){
				char *new_pivot = (char*)my_malloc(di->npivot * new_nsize);
				memset(new_pivot, 0, di->npivot * new_nsize);
				char *old_pivot = (char*)di->pivot;
				for(int j=0;j<di->npivot;j++)
					memcpy(new_pivot + j*new_nsize, old_pivot + j*old_nsize, old_nsize);
				my_free(old_pivot);
				di->pivot = new_pivot;
			}
		}
		if(MYID(simpar)==0)
			DEBUGPRINT("Sedov2D: using treevorostressrk4particletype for av_mode=%d\n", av_mode);
	}

	VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*100);

	if(icont == 0){
		t = 0;
		dt = 1.e-6;
		bp = sedov2d_mkinitial(simpar, &np);
		if(MYID(simpar)==0) printf("Sedov2D: initial conditions made (np=%d)\n", np);
		GAS_dtold(simpar) = dt;
	} else {
		int nstep = icont;
		sedov2d_readdata(simpar, &t, &dt, nstep);
		icount = nstep;
		postype dmean = SIMBOX(simpar).x.max/NX(simpar);
		GAS_dMean(simpar) = dmean;
		GAS_invw2Scale(simpar) = 1.L / SEDOV2D_PAMB(simpar);
		GAS_w2Power(simpar) = (GAS_GAMMA(simpar)-1)/GAS_GAMMA(simpar);
	}

	do {
		if(MYID(simpar)==0){ fprintf(stderr,"[ITER_DBG] A enter icount=%d t=%g\n", icount, t); fflush(stderr); }
		if(av_mode == 4 && GAS_EVOLMETHOD(simpar) == 1)
			dt = exam2d_vph_rk4_int_lagmfm(simpar,
				paddingTreeVorork4Particles, rt_w2Measure2D,
				sedov2d_evolBP,
				mkLinkedList2D_sedov2d,
				walls_xy_postStage_blend);
		else if(av_mode >= 1 && GAS_EVOLMETHOD(simpar) == 1)
			dt = exam2d_vph_rk4_int_blend(simpar,
				paddingTreeVorork4Particles, rt_w2Measure2D,
				searchCellRk4Neighbors2D, findCellRk4BP2D,
				sedov2d_evolBP,
				mkLinkedList2D_sedov2d,
				walls_xy_postStage_blend);
		else if(av_mode >= 1 && GAS_EVOLMETHOD(simpar) == 3)
			dt = exam2d_vph_kdk_int_blend(simpar,
				paddingTreeVorork4Particles, rt_w2Measure2D,
				searchCellRk4Neighbors2D, findCellRk4BP2D,
				sedov2d_evolBP,
				mkLinkedList2D_sedov2d,
				walls_xy_postStage_blend);
		else
			dt = exam2d_vph_rk4_int(simpar,
				paddingTreeVorork4Particles, rt_w2Measure2D,
				searchCellRk4Neighbors2D, findCellRk4BP2D,
				sedov2d_evolBP,
				mkLinkedList2D_sedov2d);

		if(MYID(simpar)==0){ fprintf(stderr,"[ITER_DBG] B integrator returned dt=%g\n", dt); fflush(stderr); }
		t += dt;
		GAS_dtold(simpar) = dt;
		icount++;

		/* Sedov evolves on ~0.1 timescale; snapshot every 0.01 OR every 50 steps */
		iflag = (int)(t * 100.0);
		jflag = (int)((t-dt) * 100.0);

		if(MYID(simpar)==0){ fprintf(stderr,"[ITER_DBG] C check icount=%d t=%g iflag=%d jflag=%d trigger=%d\n", icount, t, iflag, jflag, (iflag!=jflag||(icount%50)==0||icount==1)); fflush(stderr); }
		if(iflag != jflag || (icount % 50) == 0 || icount == 1){
			if(MYID(simpar)==0){ fprintf(stderr,"[ITER_DBG] D before outdata icount=%d\n", icount); fflush(stderr); }
			sedov2d_outdata(simpar, icount, t, dt);
			if(MYID(simpar)==0){ fprintf(stderr,"[ITER_DBG] E after outdata icount=%d\n", icount); fflush(stderr); }
			MPI_Barrier(MPI_COMM_WORLD);
			sedov2d_makemap(simpar, icount);
			if(MYID(simpar)==0){ fprintf(stderr,"[ITER_DBG] F after makemap icount=%d\n", icount); fflush(stderr); }
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(MYID(simpar)==0){
			fprintf(stderr,"[ITER_DBG] G iter end icount=%d t=%g dt=%g\n", icount, t, dt); fflush(stderr);
			DEBUGPRINT("Sedov2D Time= %g  step= %d  dt= %g  malloc= %d\n",
				t, icount, dt, get_malloc_call_count());
			fflush(stdout);
		}
	} while(t < 0.5);

	return 1;
}
