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


/* Snapshot format magic: -1 in the first int signals new-format snapshot
 * with full per-particle stride (so stress fields survive restart). Old
 * snapshots start with tnp >= 0 in the first int, so detection is safe. */
#define KH_SNAP_MAGIC (-1)

void kh_outdata(SimParameters *simpar, int nstep, postype t, postype dt){
	int np = VORO_NP(simpar);
	int myid, nid;
	myid = MYID(simpar);
	nid = NID(simpar);
	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	int has_stress = (p_size > sizeof(treevorork4particletype)) ? 1 : 0;
	int particle_size = (int)p_size;
	int magic = KH_SNAP_MAGIC;
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
				if(!wp) { fprintf(stderr,"P%d: cannot open %s\n",myid,outfile); MPI_Abort(MPI_COMM(simpar),1); }
				fwrite(&magic, sizeof(int),1,wp);
				fwrite(&tnp, sizeof(int),1,wp);
				fwrite(&particle_size, sizeof(int),1,wp);
				fwrite(&has_stress, sizeof(int),1,wp);
			}
			else { wp = fopen(outfile,"a"); if(!wp) { fprintf(stderr,"P%d: cannot open %s\n",myid,outfile); MPI_Abort(MPI_COMM(simpar),1); } }
			int j;
			for(j=0;j<np;j++)
				fwrite(bp_raw + j*p_size, p_size, 1, wp);
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
#ifdef USE_CUDA
	int use_stress = 1;  /* GPU blend path needs stress-type for all av_modes */
#else
	int use_stress = (av_mode >= 1 || GAS_VISCOSITY(simpar) > 0);
#endif
	size_t p_size = (use_stress) ? sizeof(treevorostressrk4particletype) : sizeof(treevorork4particletype);
	if(myid ==0){
		printf("P0: is now reading starting file: %s\n", infile);
		FILE *fp = fopen(infile,"r");
		int np;
		if(!fp) { fprintf(stderr,"P0: cannot open %s for read\n",infile); MPI_Abort(MPI_COMM(simpar),1); }
		/* Detect format: new snapshots start with magic == -1, old ones with tnp >= 0. */
		int first;
		fread(&first, sizeof(int), 1, fp);
		int file_has_stress;
		size_t file_p_size;
		if(first == KH_SNAP_MAGIC){
			int ps_int, hs_int;
			fread(&np, sizeof(int), 1, fp);
			fread(&ps_int, sizeof(int), 1, fp);
			fread(&hs_int, sizeof(int), 1, fp);
			file_p_size = (size_t)ps_int;
			file_has_stress = hs_int;
			printf("P0: new-format snapshot (np=%d p_size=%zu has_stress=%d)\n",
				np, file_p_size, file_has_stress);
		} else {
			np = first;
			file_p_size = sizeof(treevorork4particletype);
			file_has_stress = 0;
			printf("P0: legacy snapshot (np=%d, no stress fields)\n", np);
		}
		VORO_TNP(simpar) = VORO_NP(simpar) = np;
		VORORK4_TBP(simpar) = (treevorork4particletype*)my_malloc(p_size*np);
		if(use_stress){
			treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
			if(file_has_stress && file_p_size == sizeof(treevorostressrk4particletype)){
				/* Direct read: file and memory layouts match */
				fread(sbp, file_p_size, np, fp);
			} else {
				/* Legacy or smaller: read common prefix only, zero-init the rest,
				 * and seed alpha_cd = cd_amax so CD10 has a non-zero starting
				 * value (decays toward local equilibrium over ~hundreds of steps). */
				size_t common_size = offsetof(treevorork4particletype, bp);
				treevorork4particletype *tmp_buf = (treevorork4particletype*)my_malloc(file_p_size*np);
				fread(tmp_buf, file_p_size, np, fp);
				float cd_amax = GAS_CDAMAX(simpar);
				int ii;
				char *tbuf_raw = (char*)tmp_buf;
				for(ii=0;ii<np;ii++){
					memset(&sbp[ii], 0, sizeof(treevorostressrk4particletype));
					memcpy(&sbp[ii], tbuf_raw + ii*file_p_size, common_size);
					sbp[ii].stress.alpha_cd = cd_amax;
				}
				my_free(tmp_buf);
			}
		} else {
			if(file_p_size == p_size){
				fread(VORORK4_TBP(simpar), p_size, np, fp);
			} else {
				/* File has more (stress) than we need: read each particle, take common */
				size_t common_size = offsetof(treevorork4particletype, bp);
				char *scratch = (char*)my_malloc(file_p_size);
				char *dst = (char*)VORORK4_TBP(simpar);
				int ii;
				for(ii=0;ii<np;ii++){
					fread(scratch, file_p_size, 1, fp);
					memcpy(dst + ii*p_size, scratch, common_size);
				}
				my_free(scratch);
			}
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
	postype nu_phys = GAS_VISCOSITY(simpar);
#ifdef USE_CUDA
	int use_stress = 1;  /* GPU blend path needs stress-type for all av_modes */
#else
	int use_stress = (av_mode >= 1 || nu_phys > 0);
#endif
	if(use_stress){
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
			DEBUGPRINT("KH: using treevorostressrk4particletype (av_mode=%d, nu_phys=%g)\n", av_mode, nu_phys);
	}

	if(MYID(simpar)==0)
		printf("DEBUG: gpu_enabled=%d av_mode=%d use_stress=%d w2_mode=%d\n",
			GAS_GPU_ENABLED(simpar), av_mode, use_stress, GAS_W2MODE(simpar));
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
#if defined(IC_SOD2D)
#  define KH_MKLL      mkLinkedList2D_sod
#  define KH_POSTSTAGE wallx_postStage_blend
#else
#  define KH_MKLL      mkLinkedList2D_oExam
#  define KH_POSTSTAGE periodic_postStage_blend
#endif
		if(av_mode == 4 && GAS_EVOLMETHOD(simpar) == 1)
			dt = exam2d_vph_rk4_int_lagmfm(simpar,
				paddingTreeVorork4Particles,kh_w2Measure2D,
				kh_evolBP,
				KH_MKLL,
				KH_POSTSTAGE
				);
		else
		if(use_stress && GAS_EVOLMETHOD(simpar) == 1)
			dt = exam2d_vph_rk4_int_blend(simpar,
				paddingTreeVorork4Particles,kh_w2Measure2D,
				searchCellRk4Neighbors2D,findCellRk4BP2D,
				kh_evolBP,
				KH_MKLL,
				KH_POSTSTAGE
				);
		else if(use_stress && GAS_EVOLMETHOD(simpar) == 3)
			dt = exam2d_vph_kdk_int_blend(simpar,
				paddingTreeVorork4Particles,kh_w2Measure2D,
				searchCellRk4Neighbors2D,findCellRk4BP2D,
				kh_evolBP,
				KH_MKLL,
				KH_POSTSTAGE
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

		/* snapshot cadence: Sedov needs fine time sampling (short-duration blast). */
#if defined(IC_SEDOV2D) || defined(IC_NOH2D)
		iflag = t * 100.;
		jflag = (t-dt) * 100.;
#else
		iflag = t * 10.;
		jflag = (t-dt) * 10.;
#endif

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
