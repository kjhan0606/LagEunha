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
#include "../gfs_moment.h"
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
		/* The tagged moment payload aliases unused legacy RK scratch. Refuse
		 * every legacy integrator before it can overwrite that state. */
		const char *moment_init=getenv("GFS_MOMENT_INIT_ONLY");
		if(p_size==sizeof(treevorostressrk4particletype) && sizeof(PosType)==sizeof(double) &&
		   offsetof(RK4,w2backup)==24*sizeof(double)){
			treevorostressrk4particletype *check=(void*)VORORK4_TBP(simpar);
			for(int ii=0;ii<np;ii++){
				GfsMomentPayload moment;memcpy(&moment,&check[ii].rk4,sizeof(moment));
				if(moment.tag==GFS_MOMENT_TAG && (!moment_init || atoi(moment_init)<=0)){
					fprintf(stderr,"[MOMENT_STATE] legacy restart rejected: moment evolution is not implemented\n");
					MPI_Abort(MPI_COMM(simpar),99);
				}
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


static postype kh_tmax(void){
	const char *ts = getenv("HYDRO_TSTOP");
	if(ts && ts[0]) return atof(ts);
	return (postype)10;
}

int RunKH(SimParameters *simpar, int icont){
	postype t,dt;
	int np;
	treevorork4particletype *bp;
	int icount = 0;;
	int iflag,jflag;
	const char *vacuum_geometry_env = getenv("EUNHA_KEPLER_VACUUM_GEOMETRY_ONLY");
	int vacuum_geometry_only = vacuum_geometry_env != NULL;
	const char *vacuum_ale_env = getenv("EUNHA_KEPLER_VACUUM_ALE");
	int vacuum_ale = vacuum_ale_env != NULL;
	if(vacuum_geometry_only && strcmp(vacuum_geometry_env, "1") != 0){
		fprintf(stderr, "[KEPLER_VACUUM_GEOMETRY] EUNHA_KEPLER_VACUUM_GEOMETRY_ONLY must equal 1\n");
		MPI_Abort(MPI_COMM(simpar), 96);
		return 0;
	}
	if(vacuum_ale && strcmp(vacuum_ale_env, "1") != 0){
		fprintf(stderr, "[KEPLER_VACUUM_ALE] EUNHA_KEPLER_VACUUM_ALE must equal 1\n");
		MPI_Abort(MPI_COMM(simpar), 97);
		return 0;
	}
	if(vacuum_geometry_only && vacuum_ale){
		fprintf(stderr, "[KEPLER_VACUUM] choose geometry-only or ALE mode, not both\n");
		MPI_Abort(MPI_COMM(simpar), 97);
		return 0;
	}
	if((vacuum_geometry_only || vacuum_ale) && icont != 0){
		fprintf(stderr, "[KEPLER_VACUUM] first boundary implementation supports fresh IC only (icont=0)\n");
		MPI_Abort(MPI_COMM(simpar), 96);
		return 0;
	}


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
		if(vacuum_geometry_only) return 0;
		printf("made initial conditions\n");
		GAS_dtold(simpar) = dt;
		/* Save an exact IC receipt when requested, before any RK4 stage. This
		 * permits mass/volume/pressure audits for non-Cartesian generator sets. */
		const char *dump_initial = getenv("GFS_DUMP_INITIAL");
		if(dump_initial && strcmp(dump_initial, "1") == 0)
			kh_outdata(simpar, 0, t, dt);
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
		/* Optional exact endpoint for controlled comparisons/restart replay.
		 * The native RK4 still selects its physical CFL; only the last step
		 * is shortened. Ordinary runs retain their previous stepping. */
		{
			const char *exact = getenv("GFS_EXACT_STOP");
			double ceiling = exact && atoi(exact) > 0 ? (double)kh_tmax()-(double)t : 0;
			if(exact && atoi(exact) > 0 && !(ceiling > 0)) break;
			const char *cap = getenv("GFS_DIAG_DT_MAX");
			if(cap && cap[0]){
				double requested=atof(cap);
				if(!(requested>0) || !isfinite(requested)){
					fprintf(stderr,"[RK4_DT_CEILING] invalid GFS_DIAG_DT_MAX=%s\n",cap);
					MPI_Abort(MPI_COMM(simpar),97);
				}
				if(!(ceiling>0) || requested<ceiling) ceiling=requested;
			}
			gfs_rk4_set_step_context(icount, (double)t, ceiling);
		}
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
		const char *fixed_poly_env=getenv("GFS_KEPLER_FIXED_POLYGON");
		int fixed_poly=fixed_poly_env && strcmp(fixed_poly_env,"1")==0;
		void (*kh_poststage)(SimParameters *)=fixed_poly
			? fixed_polygon_postStage_blend : KH_POSTSTAGE;
		if(fixed_poly){
			int ranks=0;
			MPI_Comm_size(MPI_COMM(simpar),&ranks);
			if(ranks!=1 || GAS_Kappa(simpar)!=0 || GAS_EVOLMETHOD(simpar)!=1 || !use_stress){
				fprintf(stderr,"[FIXED_POLYGON] first native gate requires one rank, W=0, stress RK4\n");
				MPI_Abort(MPI_COMM(simpar),98);
			}
		}
		if(av_mode == 4 && GAS_EVOLMETHOD(simpar) == 1)
			dt = exam2d_vph_rk4_int_lagmfm(simpar,
				paddingTreeVorork4Particles,kh_w2Measure2D,
				kh_evolBP,
				KH_MKLL,
				kh_poststage
				);
		else
		if(use_stress && GAS_EVOLMETHOD(simpar) == 1)
			dt = exam2d_vph_rk4_int_blend(simpar,
				paddingTreeVorork4Particles,kh_w2Measure2D,
				searchCellRk4Neighbors2D,findCellRk4BP2D,
				kh_evolBP,
				KH_MKLL,
				kh_poststage
				);
		else if(use_stress && GAS_EVOLMETHOD(simpar) == 3)
			dt = exam2d_vph_kdk_int_blend(simpar,
				paddingTreeVorork4Particles,kh_w2Measure2D,
				searchCellRk4Neighbors2D,findCellRk4BP2D,
				kh_evolBP,
				KH_MKLL,
				kh_poststage
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

		/* snapshot cadence: Sedov needs fine time sampling (short-duration blast).
		 * EUNHA_DUMP_DT overrides the default, in time units per snapshot. */
		{
		const char *dump_s = getenv("EUNHA_DUMP_DT");
		double dump_dt = (dump_s && dump_s[0]) ? atof(dump_s) : 0.0;
		if(dump_dt > 0){
			iflag = t / dump_dt;
			jflag = (t-dt) / dump_dt;
		} else {
#if defined(IC_SEDOV2D) || defined(IC_NOH2D)
		iflag = t * 100.;
		jflag = (t-dt) * 100.;
#else
		iflag = t * 10.;
		jflag = (t-dt) * 10.;
#endif
		}
		}

//		if(iflag != jflag || icount%100 ==0 || icount == 1 || icount == 2)
		/* One explicitly requested diagnostic checkpoint, in addition to
		 * the existing time cadence. No extra dump unless configured. */
		const char *checkpoint_s = getenv("GFS_DIAG_CHECKPOINT_STEP");
		int checkpoint = checkpoint_s && atoi(checkpoint_s) == icount;
		/* A requested exact endpoint must be restartable even if division
		 * used by the cadence rounds the final quotient below an integer. */
		const char *exact_output_s = getenv("GFS_EXACT_STOP");
		int exact_output = exact_output_s && atoi(exact_output_s)>0 && t>=kh_tmax();
		if(iflag != jflag || checkpoint || exact_output)
		{
			kh_outdata(simpar,  icount, t, dt);
			kh_makemap(simpar, icount);
		}

		if(MYID(simpar)==0){
			printf("Time= %g  step= %d  dt= %g\n", t, icount, dt);
			fflush(stdout);
		}

	} while(t<kh_tmax());
	return 1;
}
