#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
#include "voro.h"
#include "nnost.h"
#include "gnnost.h"
#include "exam.h"
#include "exam2d.h"
#include "cylinder.h"

/* Cylinder geometry cached for cyl_evolBP (which has no simpar access) */
static postype _cyl_cx, _cyl_cy, _cyl_R;

void cyl_outdata(SimParameters *simpar, int nstep, postype t, postype dt){
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
	sprintf(outfile,"cylout.%.6d.dat",nstep);
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
}


void cyl_readdata(SimParameters *simpar, postype *t, postype *dt, int nstep){
	char infile[190];
	sprintf(infile,"cylout.%.6d.dat", nstep);
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

	DEBUGPRINT("P%d: after reading input data\n", myid);
	migrateTreeVorork4Particles(simpar);
	MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(dt, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(&GAS_dtold(simpar), 1, MPI_FLOAT,0, MPI_COMM(simpar));
	MPI_Barrier(MPI_COMM(simpar));
	BASICCELL_CELLWIDTH(simpar) = CYL_GridSize(simpar);
}

int cyl_evolBP(treevorork4particletype *bp, postype Lx, postype Ly){
	/* x: flag out-of-bounds particles (don't update in final RK4 step) */
	if(bp->x >= Lx || bp->x < 0) return 0;

	/* y: periodic wrapping */
	bp->y = fmod(bp->y + Ly, Ly);

	/* Cylinder boundary: push particles inside cylinder to surface */
	postype dx = bp->x - _cyl_cx;
	postype dy = bp->y - _cyl_cy;
	postype r = sqrt(dx*dx + dy*dy);
	if(r < _cyl_R){
		postype nx = dx/r;
		postype ny = dy/r;
		bp->x = _cyl_cx + nx * _cyl_R;
		bp->y = _cyl_cy + ny * _cyl_R;
		/* Reflect normal velocity component (free-slip) */
		postype vn = bp->vx*nx + bp->vy*ny;
		if(vn < 0){
			bp->vx -= 2*vn*nx;
			bp->vy -= 2*vn*ny;
		}
	}

	return 1;
}


int RunCylinder(SimParameters *simpar, int icont){
	postype t,dt;
	int np;
	treevorork4particletype *bp;
	int icount = 0;
	int iflag,jflag;

	{
		GridInfo *grid = &(GRIDINFO(simpar));
		PosNX(grid) = NX(simpar);
		PosNY(grid) = NY(simpar);
		PosNZ(grid) = NZ(simpar);
		PosNXNY(grid) = NX(simpar)* NY(simpar);
	}
	startRkSDD2D(simpar,CYL);
	DEBUGPRINT("P%d after startRkSDD2D: NDDINFO=%d TVORORK4[0].n_size=%zu VORORK4[0].n_size=%zu\n",
		MYID(simpar), NDDINFO(simpar),
		TVORORK4_DDINFO(simpar)[0].n_size, VORORK4_DDINFO(simpar)[0].n_size);

	int av_mode = GAS_AVMODE(simpar);
	if(av_mode >= 1){
		int ndd = NDDINFO(simpar);
		int k;
		for(k=0;k<=ndd;k++){
			TVORORK4_DDINFO(simpar)[k].n_size = sizeof(treevorostressrk4particletype);
		}
		if(MYID(simpar)==0)
			DEBUGPRINT("CYL: using treevorostressrk4particletype for av_mode=%d\n", av_mode);
	}

	/* For np=1, NDDINFO stays 0 so SIM_LXMIN would access ddinfo[-1].
	   Use SIMBOX directly — startRkSDD2D already set it from CYL_XMIN etc. */
	if(NDDINFO(simpar) > 0){
		CYL_XMIN(simpar) = SIM_LXMIN(simpar,vorork4);
		CYL_YMIN(simpar) = SIM_LYMIN(simpar,vorork4);
		CYL_XMAX(simpar) = SIM_LXMAX(simpar,vorork4);
		CYL_YMAX(simpar) = SIM_LYMAX(simpar,vorork4);
	}
	DEBUGPRINT("P%d local domain: %g %g : %g %g\n", MYID(simpar),
			CYL_XMIN(simpar), CYL_YMIN(simpar),
			CYL_XMAX(simpar), CYL_YMAX(simpar));

	VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*100);

	/* Cache cylinder geometry for cyl_evolBP */
	_cyl_cx = CYL_CX(simpar);
	_cyl_cy = CYL_CY(simpar);
	_cyl_R  = CYL_R(simpar);

	if(icont ==0) {
		t = 0;
		dt=1.e-7;
		bp = cyl_mkinitial(simpar, &np);
		printf("made initial conditions for Cylinder\n");
		GAS_dtold(simpar) = dt;
	}
	else{
		int nstep = icont;
		cyl_readdata(simpar, &t, &dt, nstep);
		icount = nstep;
	}

	postype t_max = AMAX(simpar);
	postype astep = ASTEP(simpar);
	postype out_interval = astep > 0 ? astep : 0.1;

	do {
		if(av_mode >= 1 && GAS_EVOLMETHOD(simpar) == 1)
			dt = exam2d_vph_rk4_int_blend(simpar,
				paddingCylinderParticles, cyl_w2Measure2D,
				searchCellRk4Neighbors2D, findCellRk4BP2D,
				cyl_evolBP,
				mkLinkedList2D_oExam
				);
		else if(GAS_EVOLMETHOD(simpar) == 2)
			dt = exam2d_vph_rk4_int_kNN(simpar,
				paddingCylinderParticles, cyl_w2Measure2D,
				cyl_evolBP
				);
		else if(GAS_EVOLMETHOD(simpar) == 1)
			dt = exam2d_vph_rk4_int(simpar,
				paddingCylinderParticles, cyl_w2Measure2D,
				searchCellRk4Neighbors2D, findCellRk4BP2D,
				cyl_evolBP,
				mkLinkedList2D_oExam
				);
		else dt = exam2d_vph_int(simpar,
				paddingCylinderParticles, cyl_w2Measure2D,
				searchCellRk4Neighbors2D, findCellRk4BP2D,
				cyl_evolBP,
				mkLinkedList2D_oExam
				);
		/* --- Outflow deletion: remove particles with x >= Lx or x < 0 --- */
		{
			size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
			char *bp_raw = (char*)VORORK4_TBP(simpar);
			postype Lx = SIMBOX(simpar).x.max;
			int old_np = VORO_NP(simpar);
			int new_np = 0, ii;
			for(ii = 0; ii < old_np; ii++){
				treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + ii*p_size);
				if(bpi->x >= 0 && bpi->x < Lx){
					if(new_np < ii)
						memmove(bp_raw + new_np*p_size, bp_raw + ii*p_size, p_size);
					new_np++;
				}
			}
			if(new_np < old_np){
				VORO_NP(simpar) = new_np;
				DEBUGPRINT("P%d outflow: deleted %d particles (%d -> %d)\n",
					MYID(simpar), old_np - new_np, old_np, new_np);
			}
		}

		/* --- Inflow injection: constant flux from left boundary --- */
		{
			static postype inject_accum = 0;
			postype Lx = SIMBOX(simpar).x.max;
			postype Ly = SIMBOX(simpar).y.max;
			postype dmean = GAS_dMean(simpar);
			/* Only inject on the process owning the left boundary */
			if(CYL_XMIN(simpar) <= 0){
				inject_accum += CYL_UINF(simpar) * dt;
				if(inject_accum >= dmean){
					inject_accum -= dmean;
					postype rho_inf = CYL_RHO(simpar);
					postype u_inf = CYL_UINF(simpar);
					postype p_inf = CYL_P(simpar);
					postype Gamma = GAS_GAMMA(simpar);
					postype cx = CYL_CX(simpar);
					postype cy = CYL_CY(simpar);
					postype R  = CYL_R(simpar);
					postype cs_inf = sqrt(Gamma*p_inf/rho_inf);
					postype meanvol = Lx*Ly/(NX(simpar)*NY(simpar));
					postype x_inj = dmean*0.5;
					int ny_inj = (int)(Ly/dmean);
					size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;

					/* Count how many can be injected (skip cylinder interior) */
					int ninj = 0, jj;
					for(jj = 0; jj < ny_inj; jj++){
						postype y_inj = (jj + 0.5)*Ly/ny_inj;
						postype dxc = x_inj - cx;
						postype dyc = y_inj - cy;
						if(dxc*dxc + dyc*dyc > (R + 0.5*dmean)*(R + 0.5*dmean))
							ninj++;
					}
					if(ninj > 0){
						int old_np = VORO_NP(simpar);
						int new_np = old_np + ninj;
						VORORK4_TBP(simpar) = (treevorork4particletype*)realloc(
							VORORK4_TBP(simpar), p_size * (new_np + 100));
						char *bp_raw = (char*)VORORK4_TBP(simpar);
						int idx = old_np;
						for(jj = 0; jj < ny_inj; jj++){
							postype y_inj = (jj + 0.5)*Ly/ny_inj;
							postype dxc = x_inj - cx;
							postype dyc = y_inj - cy;
							if(dxc*dxc + dyc*dyc <= (R + 0.5*dmean)*(R + 0.5*dmean))
								continue;
							treevorork4particletype *newp = (treevorork4particletype*)(bp_raw + idx*p_size);
							memset(newp, 0, p_size);
							newp->x = x_inj;
							newp->y = y_inj;
							newp->vx = u_inf;
							newp->vy = 0;
							newp->mass = rho_inf * meanvol;
							newp->den = rho_inf;
							newp->pressure = p_inf;
							newp->ie = p_inf * meanvol / (Gamma - 1);
							newp->csound = cs_inf;
							newp->volume = meanvol;
							newp->u4if.indx = 800000000L + icount*1000 + jj;
							newp->u4if.Flag[ENDIAN_OFFSET] = 0;
							if(GAS_Kappa(simpar) < 0){
								newp->w2 = -GAS_Kappa(simpar);
								newp->w2old = -GAS_Kappa(simpar);
							} else {
								newp->w2 = (dmean*GAS_Kappa(simpar))*(dmean*GAS_Kappa(simpar));
								newp->w2old = newp->w2;
							}
							newp->w2ceil = newp->w2;
							idx++;
						}
						VORO_NP(simpar) = new_np;
						VORORK4_BP(simpar) = (vorork4particletype*)VORORK4_TBP(simpar);
						DEBUGPRINT("P%d inflow: injected %d particles (%d -> %d)\n",
							MYID(simpar), ninj, old_np, new_np);
					}
				}
			}
		}

		/* Update global particle count */
		{
			int local_np = VORO_NP(simpar);
			int total_np;
			MPI_Allreduce(&local_np, &total_np, 1, MPI_INT, MPI_SUM, MPI_COMM(simpar));
			VORO_TNP(simpar) = total_np;
		}

		t += dt;

		GAS_dtold(simpar) = dt;

		icount ++;

		iflag = (int)(t / out_interval);
		jflag = (int)((t-dt) / out_interval);

		if(iflag != jflag)
		{
			cyl_outdata(simpar, icount, t, dt);
			cyl_makemap(simpar, icount);
		}

		if(MYID(simpar)==0){
			printf("Time= %g  step= %d  dt= %g\n", t, icount, dt);
			fflush(stdout);
		}

	} while(t < t_max);
	return 1;
}
