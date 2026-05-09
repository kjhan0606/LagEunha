#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "voro.h"
#include "sedov2d.h"
#include "nnost.h"
#include "exam.h"
#include "exam2d.h"
#include "color.h"

int sedov2d_makemap(SimParameters *simpar, int icount){
	postype cellsize = SEDOV2D_GridSize(simpar);
	postype Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	postype Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	int nximg = NX(simpar);
	int nyimg = NY(simpar);
	float *map = (float*)my_malloc(sizeof(float)*nximg*nyimg);
	float *img = (float*)my_malloc(sizeof(float)*nximg*nyimg);
	int i,j,ii,jj;
	for(i=0;i<nximg*nyimg;i++) map[i] = 0;
	postype pixsize = Lx/nximg;
	postype xmin = SEDOV2D_XMIN(simpar);
	postype ymin = SEDOV2D_YMIN(simpar);
	postype xmax = SEDOV2D_XMAX(simpar);
	postype ymax = SEDOV2D_YMAX(simpar);
	int mx = ceil((xmax-xmin)/cellsize);
	int my = ceil((ymax-ymin)/cellsize);
	CellType *cells = (CellType*)my_malloc(sizeof(CellType)*mx*my);
	for(i=0;i<mx*my;i++){ cells[i].link = NULL; cells[i].nmem = 0; }

	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		int ix = (int)((bpi->x-xmin)/cellsize);
		int iy = (int)((bpi->y-ymin)/cellsize);
		if(ix < 0 || ix >= mx || iy < 0 || iy >= my) continue;
		int index = ix + mx*iy;
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)bpi;
		cells[index].nmem++;
		bpi->next = tmp;
	}

	for(j=0;j<nyimg;j++){
		postype yp = (j+0.5)*pixsize;
		if(yp < ymin || yp >= ymax) continue;
		for(i=0;i<nximg;i++){
			postype xp = (i+0.5)*pixsize;
			if(xp < xmin || xp >= xmax) continue;
			postype nearden = 0;
			postype idist = 1.e20;
			int ix = (int)((xp-xmin)/cellsize);
			int jy = (int)((yp-ymin)/cellsize);
			for(jj=jy-1;jj<=jy+1;jj++){
				if(jj<0 || jj>=my) continue;
				for(ii=ix-1;ii<=ix+1;ii++){
					if(ii<0 || ii>=mx) continue;
					size_t ipixel = ii + mx*jj;
					struct linkedlisttype *tmp = cells[ipixel].link;
					while(tmp){
						treevorork4particletype *tt = (treevorork4particletype*)tmp;
						if(tt->x >= xmin && tt->x < xmax && tt->y >= ymin && tt->y < ymax){
							postype distx = fabs(tt->x - xp);
							postype disty = fabs(tt->y - yp);
							postype dist2 = distx*distx + disty*disty;
							if(dist2 < idist){ idist = dist2; nearden = tt->den; }
						}
						tmp = tmp->next;
					}
				}
			}
			map[i+nximg*(nyimg-j-1)] = nearden;
		}
	}

	MPI_Reduce(map, img, nximg*nyimg, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	if(MYID(simpar)==0){
		char outfile[189];
		char outsao[180];
		sprintf(outfile,"sedov2dmap.%.6d.ppm", icount);
		strcpy(outsao,"rt.sao");
		colorizeit(5, img, nximg, nyimg, outsao, outfile);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(cells) my_free(cells);
	if(img) my_free(img);
	if(map) my_free(map);
	return 0;
}

treevorork4particletype *sedov2d_mkinitial(SimParameters *simpar, int *mp){
	int i,j;
	treevorork4particletype *res;
	postype xmin = SEDOV2D_XMIN(simpar);
	postype ymin = SEDOV2D_YMIN(simpar);
	postype xmax = SEDOV2D_XMAX(simpar);
	postype ymax = SEDOV2D_YMAX(simpar);
	int nx = NX(simpar);
	int ny = NY(simpar);
	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	postype dmean = Lx/nx;
	postype rho_amb = SEDOV2D_RHOAMB(simpar);
	postype P_amb = SEDOV2D_PAMB(simpar);
	postype E_blast = SEDOV2D_E(simpar);
	postype r_blast = SEDOV2D_RBLAST(simpar);
	postype cx = SEDOV2D_CX(simpar);
	postype cy = SEDOV2D_CY(simpar);
	int av_mode = GAS_AVMODE(simpar);
	postype Gamma = GAS_GAMMA(simpar);

	GAS_dMean(simpar) = dmean;

	fprintf(stderr,"[Sedov2D-DBG] P%d enter mkinitial: nx/y=%d %d L=%g %g rho_amb=%g P_amb=%g E_blast=%g r_blast=%g\n",
		MYID(simpar), nx, ny, Lx, Ly, rho_amb, P_amb, E_blast, r_blast); fflush(stderr);

	/* Concentric rings IC overshoots before box clipping; allocate 2x for safety */
	size_t buf_capacity = (size_t)(2 * nx * ny);
	if(av_mode >= 1)
		res = (treevorork4particletype*)my_malloc(sizeof(treevorostressrk4particletype)*buf_capacity);
	else
		res = (treevorork4particletype*)my_malloc(sizeof(treevorork4particletype)*buf_capacity);

	postype meanvol = Lx*Ly/(postype)nx/(postype)ny;

	GAS_invw2Scale(simpar) = 1.L / P_amb;
	GAS_w2Power(simpar) = (Gamma-1)/Gamma;

	/* Energy injection: E_blast deposited as thermal energy in r < r_blast.
	   Total area A = π r_blast². Pressure inside: P_blast/(γ-1) · A = E_blast
	   → P_blast = E_blast (γ-1) / A. */
	postype A_blast = M_PI * r_blast * r_blast;
	postype P_blast = E_blast * (Gamma - 1) / A_blast;

	/* IC generated only on rank 0; migrate distributes after */
	int myid = MYID(simpar);
	int rank0_only = (myid == 0);

	/* IC mode selector via env var SEDOV2D_IC=brick to use legacy nx*ny grid;
	   default is concentric rings. */
	const char *ic_env = getenv("SEDOV2D_IC");
	int use_brick_ic = (ic_env && strcmp(ic_env, "brick") == 0);

	/* Concentric rings IC: place particles on rings r_k = k*dr around (cx,cy),
	   with N_k = round(2π·r_k/dr) particles per ring (uniform area density).
	   Alternating angular phase between rings breaks radial spoke alignment.
	   Particles outside the box are clipped. The central anchor is placed last. */
	size_t np = 0;
	if(rank0_only && !use_brick_ic){
		postype dx_ic = Lx/(postype)nx;
		postype dr = dx_ic;
		/* k_max covers the box corner farthest from center */
		postype r_corner = sqrt((cx-xmin)*(cx-xmin) + (cy-ymin)*(cy-ymin));
		postype r2;
		r2 = sqrt((xmax-cx)*(xmax-cx) + (cy-ymin)*(cy-ymin)); if(r2 > r_corner) r_corner = r2;
		r2 = sqrt((cx-xmin)*(cx-xmin) + (ymax-cy)*(ymax-cy)); if(r2 > r_corner) r_corner = r2;
		r2 = sqrt((xmax-cx)*(xmax-cx) + (ymax-cy)*(ymax-cy)); if(r2 > r_corner) r_corner = r2;
		int k_max = (int)(r_corner/dr) + 2;

		int k;
		for(k=1; k<=k_max; k++){
			postype rk = (postype)k * dr;
			int Nk = (int)(2.0*M_PI*rk/dr + 0.5);
			if(Nk < 6) Nk = 6;
			postype dtheta = 2.0*M_PI/(postype)Nk;
			postype theta_offset = (k%2) ? (dtheta*0.5) : 0.0;
			int m;
			for(m=0; m<Nk; m++){
				postype theta = (postype)m * dtheta + theta_offset;
				postype x = cx + rk*cos(theta);
				postype y = cy + rk*sin(theta);
				if(x < xmin || x >= xmax || y < ymin || y >= ymax) continue;

				postype rho = rho_amb;
				postype P = P_amb;

				UNSET_FLAG(res+np, Wallflag);
				res[np].u4if.indx = (size_t)nx*(size_t)ny + 1 + np;  /* ring indices start above anchor's nx*ny */
				res[np].x = x;
				res[np].y = y;
				res[np].vx = 0; res[np].vy = 0;
				res[np].z = 0; res[np].vz = 0;
				res[np].ax = 0; res[np].ay = 0; res[np].az = 0;

				res[np].mass = rho*meanvol;
				res[np].den = rho;
				res[np].pressure = P;
				res[np].ie = P*meanvol/(Gamma-1);
				res[np].csound = sqrt(Gamma*P/rho);
				res[np].w2 = (Lx/nx*GAS_Kappa(simpar))*(Lx/nx*GAS_Kappa(simpar));
				res[np].w2old = res[np].w2;
				res[np].die = 0;
				if(GAS_Kappa(simpar) < 0){
					res[np].w2 = -GAS_Kappa(simpar);
					res[np].w2old = -GAS_Kappa(simpar);
				} else {
					res[np].avgNeighboringPressure = res[np].pressure;
					res[np].w2 = getw2forHydroParticle(simpar, (res+np), 1);
					res[np].w2old = res[np].w2;
				}
				res[np].w2ceil = res[np].w2;
				np++;
			}
		}
		/* Anchor: marker index = nx*ny. Ring indices start at nx*ny+1, so no collision. */
		size_t anchor_idx = (size_t)nx * (size_t)ny;
		UNSET_FLAG(res+np, Wallflag);
		res[np].u4if.indx = anchor_idx;
		res[np].x = cx; res[np].y = cy;
		res[np].vx = 0; res[np].vy = 0;
		res[np].z = 0; res[np].vz = 0;
		res[np].ax = 0; res[np].ay = 0; res[np].az = 0;
		res[np].mass = rho_amb*meanvol;
		res[np].den = rho_amb;
		res[np].pressure = P_amb;
		res[np].ie = P_amb*meanvol/(Gamma-1);
		res[np].csound = sqrt(Gamma*P_amb/rho_amb);
		res[np].die = 0;
		if(GAS_Kappa(simpar) < 0){
			res[np].w2 = -GAS_Kappa(simpar);
			res[np].w2old = -GAS_Kappa(simpar);
		} else {
			res[np].avgNeighboringPressure = res[np].pressure;
			res[np].w2 = getw2forHydroParticle(simpar, (res+np), 1);
			res[np].w2old = res[np].w2;
		}
		res[np].w2ceil = res[np].w2;
		np++;

		fprintf(stderr,"[Sedov2D-RINGS] k_max=%d np=%zu (target nx*ny+1=%d, ratio=%.4f)\n",
			k_max, np, nx*ny+1, (double)np/(double)(nx*ny+1)); fflush(stderr);
	}

	/* Brick-stagger IC: regular nx*ny grid with half-cell offset per row.
	   Followed by Lloyd glass relaxation downstream (n_glass>0). */
	if(rank0_only && use_brick_ic){
		postype dx_ic = Lx/(postype)nx;
		postype dy_ic = Ly/(postype)ny;
		int ix, iy;
		for(iy=0; iy<ny; iy++){
			postype y = ymin + (iy + 0.5)*dy_ic;
			postype xshift = (iy & 1) ? 0.5*dx_ic : 0.0;
			for(ix=0; ix<nx; ix++){
				postype x = xmin + (ix + 0.5)*dx_ic + xshift;
				if(x >= xmax) x -= Lx;
				if(x < xmin || x >= xmax || y < ymin || y >= ymax) continue;

				postype rho = rho_amb;
				postype P = P_amb;

				UNSET_FLAG(res+np, Wallflag);
				res[np].u4if.indx = (size_t)nx*(size_t)ny + 1 + np;
				res[np].x = x;
				res[np].y = y;
				res[np].vx = 0; res[np].vy = 0;
				res[np].z = 0; res[np].vz = 0;
				res[np].ax = 0; res[np].ay = 0; res[np].az = 0;

				res[np].mass = rho*meanvol;
				res[np].den = rho;
				res[np].pressure = P;
				res[np].ie = P*meanvol/(Gamma-1);
				res[np].csound = sqrt(Gamma*P/rho);
				res[np].w2 = (Lx/nx*GAS_Kappa(simpar))*(Lx/nx*GAS_Kappa(simpar));
				res[np].w2old = res[np].w2;
				res[np].die = 0;
				if(GAS_Kappa(simpar) < 0){
					res[np].w2 = -GAS_Kappa(simpar);
					res[np].w2old = -GAS_Kappa(simpar);
				} else {
					res[np].avgNeighboringPressure = res[np].pressure;
					res[np].w2 = getw2forHydroParticle(simpar, (res+np), 1);
					res[np].w2old = res[np].w2;
				}
				res[np].w2ceil = res[np].w2;
				np++;
			}
		}
		/* Anchor */
		size_t anchor_idx = (size_t)nx * (size_t)ny;
		UNSET_FLAG(res+np, Wallflag);
		res[np].u4if.indx = anchor_idx;
		res[np].x = cx; res[np].y = cy;
		res[np].vx = 0; res[np].vy = 0;
		res[np].z = 0; res[np].vz = 0;
		res[np].ax = 0; res[np].ay = 0; res[np].az = 0;
		res[np].mass = rho_amb*meanvol;
		res[np].den = rho_amb;
		res[np].pressure = P_amb;
		res[np].ie = P_amb*meanvol/(Gamma-1);
		res[np].csound = sqrt(Gamma*P_amb/rho_amb);
		res[np].die = 0;
		if(GAS_Kappa(simpar) < 0){
			res[np].w2 = -GAS_Kappa(simpar);
			res[np].w2old = -GAS_Kappa(simpar);
		} else {
			res[np].avgNeighboringPressure = res[np].pressure;
			res[np].w2 = getw2forHydroParticle(simpar, (res+np), 1);
			res[np].w2old = res[np].w2;
		}
		res[np].w2ceil = res[np].w2;
		np++;

		fprintf(stderr,"[Sedov2D-BRICK] np=%zu (target nx*ny+1=%d)\n",
			np, nx*ny+1); fflush(stderr);
	}

	DEBUGPRINT("P%d Sedov2D: np=%ld meanvol=%g P_blast=%g (rho_amb=%g P_amb=%g)\n",
		MYID(simpar), np, meanvol, P_blast, rho_amb, P_amb);

	if(av_mode >= 1){
		res = (treevorork4particletype*)realloc(res, sizeof(treevorostressrk4particletype)*np);
		treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)res;
		char *old_base = (char*)res;
		size_t old_size = sizeof(treevorork4particletype);
		int ii;
		for(ii=np-1; ii>=0; ii--){
			memmove(&sbp[ii], old_base + ii*old_size, old_size);
			memset(&sbp[ii].stress, 0, sizeof(Stress));
			sbp[ii].bp = NULL;
		}
	} else {
		res = (treevorork4particletype*)realloc(res, sizeof(treevorork4particletype)*np);
	}

	fprintf(stderr,"[Sedov2D-DBG] P%d after IC loop: np=%ld P_blast=%g\n",
		MYID(simpar), np, P_blast); fflush(stderr);

	int nbp = np;
	*mp = nbp;
	VORORK4_TBP(simpar) = res;
	VORORK4_BP(simpar) = (vorork4particletype*)res;
	VORO_NP(simpar) = nbp;

	{
		size_t p_size = (av_mode >= 1) ? sizeof(treevorostressrk4particletype)
		                               : sizeof(treevorork4particletype);
		char *bp_raw = (char*)res;
		for(i=0;i<np;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			bpi->u4if.Flag[ENDIAN_OFFSET] &= (~BoundaryGhostflag);
		}
	}

	fprintf(stderr,"[Sedov2D-DBG] P%d before migrate: np=%d\n", MYID(simpar), VORO_NP(simpar)); fflush(stderr);
	migrateTreeVorork4Particles(simpar);
	fprintf(stderr,"[Sedov2D-DBG] P%d after migrate: np=%d\n", MYID(simpar), VORO_NP(simpar)); fflush(stderr);

	exam2dUpdateVol(simpar, paddingTreeVorork4Particles,
		searchCellRk4Neighbors2D, findCellRk4BP2D,
		mkLinkedList2D_sedov2d);

	/* Glass relaxation: Lloyd-style centroid shifts to homogenize Voronoi cells.
	   The central anchor (indx == nx*ny) is restored to (cx,cy) after every shift
	   so the blast center stays exact.
	   With concentric rings IC, particles are already at uniform density and
	   radially isotropic; Lloyd would distort rings into hex patches and break
	   isotropy. Set n_glass=0 to preserve ring layout. */
	{
		int n_glass = use_brick_ic ? 100 : 0;
		float fshift_save = GAS_FCENTROID(simpar);
		GAS_FCENTROID(simpar) = 0.5f;
		size_t anchor_idx = (size_t)nx * (size_t)ny;
		int it;
		for(it=0; it<n_glass; it++){
			exam2d_centroidShift(simpar, paddingTreeVorork4Particles,
				searchCellRk4Neighbors2D, findCellRk4BP2D,
				mkLinkedList2D_sedov2d);
			/* Pin anchor exactly at (cx,cy) */
			{
				size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
				char *bp_raw = (char*)VORORK4_TBP(simpar);
				int kk, npp = VORO_NP(simpar);
				for(kk=0;kk<npp;kk++){
					treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + kk*p_size);
					if(bpi->u4if.indx == anchor_idx){
						bpi->x = cx; bpi->y = cy;
						bpi->vx = 0; bpi->vy = 0;
					}
				}
			}
			exam2dUpdateVol(simpar, paddingTreeVorork4Particles,
				searchCellRk4Neighbors2D, findCellRk4BP2D,
				mkLinkedList2D_sedov2d);
			if((it+1) % 20 == 0){
				size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
				char *bp_raw = (char*)VORORK4_TBP(simpar);
				int kk, npp = VORO_NP(simpar);
				postype Vmin=1e20, Vmax=0;
				for(kk=0;kk<npp;kk++){
					treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + kk*p_size);
					if(bpi->u4if.indx == anchor_idx) continue;
					if(bpi->volume < Vmin) Vmin = bpi->volume;
					if(bpi->volume > Vmax) Vmax = bpi->volume;
				}
				postype gVmin, gVmax;
				MPI_Allreduce(&Vmin, &gVmin, 1, MPI_POSTYPE, MPI_MIN, MPI_COMM_WORLD);
				MPI_Allreduce(&Vmax, &gVmax, 1, MPI_POSTYPE, MPI_MAX, MPI_COMM_WORLD);
				if(MYID(simpar)==0){
					fprintf(stderr,"[Sedov2D-GLASS] it=%d Vmin=%g Vmax=%g ratio=%g\n",
						it+1, gVmin, gVmax, gVmax/gVmin); fflush(stderr);
				}
			}
		}
		GAS_FCENTROID(simpar) = fshift_save;
	}

	res = VORORK4_TBP(simpar);
	nbp = VORO_NP(simpar);

	{
		size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
		char *bp_raw = (char*)res;
		int entropy_mode = (av_mode >= 1) ? GAS_ENTROPY_MODE(simpar) : 0;
		double Vmin=1e30, Vmax=0, denmin=1e30, denmax=0, Kmax=0;
		int idx_Vmax=-1, idx_denmin=-1, idx_Kmax=-1;
		for(i=0;i<nbp;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			bpi->mass = bpi->den * bpi->volume;
			bpi->ie = bpi->pressure * bpi->volume / (Gamma-1);
			bpi->csound = sqrt(Gamma * bpi->pressure / bpi->den);
			if(entropy_mode == 1){
				treevorostressrk4particletype *sbi = (treevorostressrk4particletype*)bpi;
				sbi->stress.K  = bpi->pressure / pow((double)bpi->den, (double)Gamma);
				sbi->stress.dK = 0;
				if(sbi->stress.K > Kmax){ Kmax = sbi->stress.K; idx_Kmax = i; }
			}
			if(bpi->volume < Vmin) Vmin = bpi->volume;
			if(bpi->volume > Vmax){ Vmax = bpi->volume; idx_Vmax = i; }
			if(bpi->den < denmin){ denmin = bpi->den; idx_denmin = i; }
			if(bpi->den > denmax) denmax = bpi->den;
		}
		fprintf(stderr,"[Sedov2D-IC_DIAG] P%d Vmin=%g Vmax=%g denmin=%g denmax=%g Kmax=%g\n",
			MYID(simpar), Vmin, Vmax, denmin, denmax, Kmax); fflush(stderr);
		if(idx_Vmax >= 0){
			treevorork4particletype *b = (treevorork4particletype*)(bp_raw + idx_Vmax*p_size);
			fprintf(stderr,"[Sedov2D-IC_VMAX] P%d idx=%d x=%g y=%g V=%g den=%g P=%g mass=%g\n",
				MYID(simpar), idx_Vmax, (double)b->x, (double)b->y, b->volume, b->den, b->pressure, b->mass); fflush(stderr);
		}
		if(idx_denmin >= 0){
			treevorork4particletype *b = (treevorork4particletype*)(bp_raw + idx_denmin*p_size);
			fprintf(stderr,"[Sedov2D-IC_DENMIN] P%d idx=%d x=%g y=%g V=%g den=%g P=%g mass=%g\n",
				MYID(simpar), idx_denmin, (double)b->x, (double)b->y, b->volume, b->den, b->pressure, b->mass); fflush(stderr);
		}
		if(idx_Kmax >= 0 && entropy_mode==1){
			treevorork4particletype *b = (treevorork4particletype*)(bp_raw + idx_Kmax*p_size);
			treevorostressrk4particletype *sb = (treevorostressrk4particletype*)b;
			fprintf(stderr,"[Sedov2D-IC_KMAX] P%d idx=%d x=%g y=%g V=%g den=%g P=%g K=%g\n",
				MYID(simpar), idx_Kmax, (double)b->x, (double)b->y, b->volume, b->den, b->pressure, sb->stress.K); fflush(stderr);
		}
	}
	VORO_NPAD(simpar) = 0;

	fprintf(stderr,"[Sedov2D-DBG] P%d after exam2dUpdateVol: np=%d\n", MYID(simpar), VORO_NP(simpar)); fflush(stderr);
	if(GAS_Kappa(simpar) > 0) det2d_dpqRK4(simpar, paddingTreeVorork4Particles);

	/* Blast energy injection on the relaxed mesh.  Compute total volume of
	   particles inside r<r_blast, then set P_blast such that integrated
	   thermal energy = E_blast.  Anchor particle's mass is boosted last. */
	{
		size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
		char *bp_raw = (char*)VORORK4_TBP(simpar);
		int kk, npp = VORO_NP(simpar);
		size_t anchor_idx = (size_t)nx * (size_t)ny;
		postype Vloc = 0;
		int n_loc = 0;
		for(kk=0;kk<npp;kk++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + kk*p_size);
			postype dx = bpi->x - cx;
			postype dy = bpi->y - cy;
			if(dx*dx + dy*dy < r_blast*r_blast){
				Vloc += bpi->volume;
				n_loc++;
			}
		}
		postype Vtot;
		int n_tot;
		MPI_Allreduce(&Vloc, &Vtot, 1, MPI_POSTYPE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&n_loc, &n_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		postype P_blast_inj = E_blast * (Gamma - 1) / Vtot;
		if(MYID(simpar)==0){
			fprintf(stderr,"[Sedov2D-BLAST] n_inj=%d Vtot=%g P_blast=%g\n",
				n_tot, Vtot, P_blast_inj); fflush(stderr);
		}
		/* For LagMFM (av_mode=4), keep anchor mass = 1x to avoid contaminating
		   neighbors' kernel-sum density. Anchor immobility is enforced by
		   position pinning in walls_xy_postStage_blend regardless of mass. */
		int av_mode_now = GAS_AVMODE(simpar);
		postype M_anchor = (av_mode_now == 4) ? rho_amb*meanvol
		                                       : 1.0e3 * rho_amb * meanvol;
		for(kk=0;kk<npp;kk++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + kk*p_size);
			postype dx = bpi->x - cx;
			postype dy = bpi->y - cy;
			int in_blast = (dx*dx + dy*dy < r_blast*r_blast);
			int is_anchor = (bpi->u4if.indx == anchor_idx);
			if(is_anchor){
				bpi->x = cx; bpi->y = cy;
				bpi->vx = 0; bpi->vy = 0;
				bpi->ax = 0; bpi->ay = 0;
				bpi->mass = M_anchor;
				bpi->den = bpi->mass / bpi->volume;
			}
			if(in_blast || is_anchor){
				bpi->pressure = P_blast_inj;
				bpi->ie = P_blast_inj * bpi->volume / (Gamma-1);
				bpi->csound = sqrt(Gamma * P_blast_inj / bpi->den);
				bpi->avgNeighboringPressure = P_blast_inj;
				if(av_mode_now >= 1 && GAS_ENTROPY_MODE(simpar) == 1){
					treevorostressrk4particletype *sbi = (treevorostressrk4particletype*)bpi;
					sbi->stress.K  = P_blast_inj / pow((double)bpi->den, (double)Gamma);
					sbi->stress.dK = 0;
				}
				if(GAS_Kappa(simpar) > 0){
					bpi->w2 = getw2forHydroParticle(simpar, bpi, 1);
					bpi->w2old = bpi->w2;
				}
			}
		}
		if(GAS_Kappa(simpar) > 0) det2d_dpqRK4(simpar, paddingTreeVorork4Particles);
	}

	fprintf(stderr,"[Sedov2D-DBG] P%d mkinitial DONE\n", MYID(simpar)); fflush(stderr);

	return res;
}
