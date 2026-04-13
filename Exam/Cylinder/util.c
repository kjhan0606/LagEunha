#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "voro.h"
#include "nnost.h"
#include "gnnost.h"
#include "exam.h"
#include "exam2d.h"
#include "cylinder.h"
#include "color.h"


static CellType *cells = NULL;

int cyl_makemap(SimParameters *simpar, int icount){
	postype cellsize = CYL_GridSize(simpar);
	postype Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	postype Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	int nximg = NX(simpar);
	int nyimg = NY(simpar);
	float *map = (float*)malloc(sizeof(float)*nximg*nyimg);
	float *img = (float*)malloc(sizeof(float)*nximg*nyimg);
	postype meanvol = Lx*Ly/nximg/nyimg;
	int i,j;
	int ii,jj;
	for(i=0;i<nximg*nyimg;i++) map[i] = 0;
	postype pixsize = Lx/nximg;
	postype xmin,ymin,xmax,ymax;
	xmin = CYL_XMIN(simpar);
	ymin = CYL_YMIN(simpar);
	xmax = CYL_XMAX(simpar);
	ymax = CYL_YMAX(simpar);
	int mx,my;
	mx = ceil((xmax-xmin)/cellsize);
	my = ceil((ymax-ymin)/cellsize);
	cells = (CellType*)malloc(sizeof(CellType)*mx*my);
	for(i=0;i<mx*my;i++){
		cells[i].link = NULL;
		cells[i].nmem = 0;
	}
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		int ix,iy;
		ix = ((bpi->x-xmin)/cellsize);
		iy = ((bpi->y-ymin)/cellsize);
		if(ix<0||ix>=mx||iy<0||iy>=my) continue;
		int index = ix+mx*iy;
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)bpi;
		cells[index].nmem ++;
		bpi->next = tmp;
	}

	postype cx = CYL_CX(simpar);
	postype cy = CYL_CY(simpar);
	postype R  = CYL_R(simpar);

	for(j=0;j<nyimg;j++){
		postype yp = (j+0.5)*pixsize;
		if(yp < ymin || yp >= ymax) continue;
		for(i=0;i<nximg;i++){
			postype xp = (i+0.5)*pixsize;
			if(xp < xmin || xp >= xmax) continue;
			/* Mask cylinder interior — skip nearest-neighbor search */
			postype dxc = xp - cx;
			postype dyc = yp - cy;
			if(dxc*dxc + dyc*dyc <= R*R) continue;
			postype nearden = 0;
			postype idist = 1.e20;
			int jy,ix;
			ix = (xp-xmin)/cellsize;
			jy = (yp-ymin)/cellsize;
			for(jj=jy-1;jj<=jy+1;jj++){
				if(jj<0 || jj>=my) continue;
				for(ii=ix-1;ii<=ix+1;ii++){
					if(ii<0 || ii>=mx) continue;
					size_t ipixel = ii + mx*jj;
					struct linkedlisttype *tmp = cells[ipixel].link;
					while(tmp){
						treevorork4particletype *tt = (treevorork4particletype*)tmp;
						if(tt->x >= xmin && tt->x < xmax && tt->y >=ymin && tt->y < ymax){
							postype distx = fabs(tt->x - xp);
							postype disty = fabs(tt->y - yp);
							postype dist2 = distx*distx + disty*disty;
							if(dist2 < idist) {
								idist = dist2;
								nearden = tt->den;
							}
						}
						tmp = tmp->next;
					}
				}
			}
			map[i+nximg*(nyimg-j-1)] = nearden;
		}
	}
	MPI_Reduce(map, img, nximg*nyimg, MPI_FLOAT, MPI_MAX, 0, MPI_COMM(simpar));

    if(MYID(simpar)==0){
        char outfile[189];
        sprintf(outfile,"cylmap.%.6d.ppm", icount);
        colorizeit(7, img, nximg, nyimg, "cyl.sao", outfile, 0.8, 1.2);
    }
    free(cells);
    free(map); free(img);
	return 0;
}


treevorork4particletype *cyl_mkinitial(SimParameters *simpar, int *mp){
	int i,j;
	treevorork4particletype *res;
	postype rho_inf = CYL_RHO(simpar);
	postype u_inf = CYL_UINF(simpar);
	postype p_inf = CYL_P(simpar);
	postype cx = CYL_CX(simpar);
	postype cy = CYL_CY(simpar);
	postype R  = CYL_R(simpar);
	postype xmin,ymin,xmax,ymax;
	xmin = CYL_XMIN(simpar);
	ymin = CYL_YMIN(simpar);
	xmax = CYL_XMAX(simpar);
	ymax = CYL_YMAX(simpar);
	int nx = NX(simpar);
	int ny = NY(simpar);
	postype Lx = CYL_SIMBOX(simpar).x.max - CYL_SIMBOX(simpar).x.min;
	postype Ly = CYL_SIMBOX(simpar).y.max - CYL_SIMBOX(simpar).y.min;
	postype dmean = Lx/nx;
	GAS_dMean(simpar) = dmean;
	int av_mode = GAS_AVMODE(simpar);

	postype Gamma = GAS_GAMMA(simpar);

	if(av_mode >= 1)
		res = (treevorork4particletype*)my_malloc(sizeof(treevorostressrk4particletype)*nx*ny);
	else
		res = (treevorork4particletype*)my_malloc(sizeof(treevorork4particletype)*nx*ny);
	postype meanvol = Lx*Ly/nx/ny;
	size_t np = 0;
	for(j=0;j<ny;j++){
		postype y = (postype)(j+0.5)*Ly/(postype)ny;
		for(i=0;i<nx;i++){
			postype x = (postype)(i+0.5)*Lx/(postype)nx;
			/* skip particles inside or too close to cylinder surface.
			   Gap of 0.5*dmean ensures ghost-real distance >= dmean,
			   giving reasonable CFL timestep near the cylinder. */
			postype dx = x - cx;
			postype dy = y - cy;
			postype r = sqrt(dx*dx + dy*dy);
			if(r <= R + 0.5*dmean) continue;

			size_t indx = i+nx*j;
			if(x>=xmin && x < xmax && y>=ymin && y < ymax){
				res[np].u4if.indx  = indx;
				res[np].u4if.Flag[ENDIAN_OFFSET]  = 0;
				res[np].x  = x;
				res[np].y  = y;
				res[np].vx  = u_inf;
				res[np].vy  = 0;
				res[np].mass  = rho_inf*meanvol;
				res[np].den  = rho_inf;
				res[np].pressure  = p_inf;
				res[np].ie  = p_inf*meanvol/(Gamma-1);
				res[np].csound  = sqrt(Gamma*p_inf/rho_inf);
				res[np].die  = 0;
				if(GAS_Kappa(simpar) < 0) {
					res[np].w2 = -GAS_Kappa(simpar);
					res[np].w2old = -GAS_Kappa(simpar);
				}
				else {
					res[np].w2  = (dmean*GAS_Kappa(simpar))*(dmean*GAS_Kappa(simpar));
					res[np].w2old  = res[np].w2;
					res[np].avgNeighboringPressure = p_inf;
					res[np].w2 = getw2forHydroParticle(simpar, res+np, 0);
					res[np].w2old = res[np].w2;
				}
				res[np].w2ceil = res[np].w2;
				np++;
			}
		}
	}
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

	GAS_invw2Scale(simpar) = 1.L / p_inf;
	GAS_w2Power(simpar) = (Gamma-1)/Gamma;

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
	DEBUGPRINT("P%d has np = %ld in box %g %g : %g %g\n", MYID(simpar), np, xmin,ymin,xmax,ymax);

	migrateTreeVorork4Particles(simpar);

	void cyl_findVol(SimParameters *);
	cyl_findVol(simpar);

	res = VORORK4_TBP(simpar);
	nbp = VORO_NP(simpar);

	{
		size_t p_size = (av_mode >= 1) ? sizeof(treevorostressrk4particletype)
		                               : sizeof(treevorork4particletype);
		char *bp_raw = (char*)res;
		for(i=0;i<nbp;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			bpi->ie = bpi->pressure*bpi->volume/(Gamma-1);
			bpi->csound = sqrt(Gamma*bpi->pressure/bpi->den);
		}
	}
	VORO_NPAD(simpar) = 0;

	/* Debug: check for negative ie after IC */
	{
		size_t p_size2 = (av_mode >= 1) ? sizeof(treevorostressrk4particletype)
		                               : sizeof(treevorork4particletype);
		char *bp_raw2 = (char*)res;
		int neg_ie_count = 0;
		for(i=0;i<nbp;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw2 + i*p_size2);
			if(bpi->ie < 0 || bpi->volume <= 0){
				if(neg_ie_count < 5)
					DEBUGPRINT("P%d IC bad: i=%d PINDX=%ld ie=%g vol=%g den=%g P=%g pos=(%g,%g)\n",
						MYID(simpar), i, PINDX(bpi), bpi->ie, bpi->volume, bpi->den,
						bpi->pressure, bpi->x, bpi->y);
				neg_ie_count++;
			}
		}
		if(neg_ie_count > 0) DEBUGPRINT("P%d IC: %d particles have bad ie or volume\n", MYID(simpar), neg_ie_count);
		else DEBUGPRINT("P%d IC: all %d particles have positive ie and volume\n", MYID(simpar), nbp);
	}

	DEBUGPRINT("P%d has w2 values P= %g w2= %g w2ceil = %g den= %g meanV= %g vol= %g Cs= %g\n",
			MYID(simpar), res[0].pressure, res[0].w2, res[0].w2ceil, res[0].den,
			meanvol, res[0].volume, res[0].csound);
	if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar, paddingCylinderParticles);
	return res;
}


void paddingCylinderParticles(SimParameters *simpar, postype width){
	/* 1. Standard MPI padding (skip for single-process: NDDINFO=0) */
	if(NDDINFO(simpar) > 0)
		paddingTreeVorork4Particles(simpar, width);
	else {
		VORORK4_TBPP(simpar) = NULL;
		VORO_NPAD(simpar) = 0;
	}

	/* 1a. Filter out periodically-wrapped MPI ghosts in x.
	   ppadding() assumes periodic boundaries, but cylinder x is non-periodic.
	   Remove ghosts whose x falls outside [0,Lx]. y is periodic so no y filter. */
	{
		size_t p_sz = TVORORK4_DDINFO(simpar)[0].n_size;
		postype Lx_f = SIMBOX(simpar).x.max;
		char *pp_raw = (char*)VORORK4_TBPP(simpar);
		int npad_orig = VORO_NPAD(simpar);
		int kept = 0, ii;
		for(ii = 0; ii < npad_orig; ii++){
			treevorork4particletype *gpi = (treevorork4particletype*)(pp_raw + ii*p_sz);
			if(gpi->x >= 0 && gpi->x <= Lx_f){
				if(kept < ii)
					memmove(pp_raw + kept*p_sz, pp_raw + ii*p_sz, p_sz);
				kept++;
			}
		}
		VORO_NPAD(simpar) = kept;
/*		DEBUGPRINT("P%d paddingCyl: MPI ghosts %d -> %d after x-periodic filter\n",
			MYID(simpar), npad_orig, kept); */
	}

	int np = VORO_NP(simpar);
	int npad_mpi = VORO_NPAD(simpar);
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;

	postype cx = CYL_CX(simpar);
	postype cy = CYL_CY(simpar);
	postype R  = CYL_R(simpar);
	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	postype u_inf = CYL_UINF(simpar);
	postype rho_inf = CYL_RHO(simpar);
	postype p_inf = CYL_P(simpar);
	postype Gamma = GAS_GAMMA(simpar);
	postype cs_inf = sqrt(Gamma*p_inf/rho_inf);

/*	DEBUGPRINT("P%d paddingCyl: width=%g Lx=%g Ly=%g cx=%g cy=%g R=%g p_size=%zu sizeof(tree)=%zu\n",
		MYID(simpar), width, Lx, Ly, cx, cy, R, p_size, sizeof(treevorork4particletype)); */

	/* Source particles: local + MPI padding (need wall ghosts for both) */
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	char *bpp_raw = (char*)VORORK4_TBPP(simpar);
	int nsrc = np + npad_mpi;
	int i;

	postype xmin_local = CYL_XMIN(simpar);
	postype ymin_local = CYL_YMIN(simpar);
	postype xmax_local = CYL_XMAX(simpar);
	postype ymax_local = CYL_YMAX(simpar);

	/* Helper to get source particle by index (local or MPI padding) */
	#define SRC_PTL(idx) ((treevorork4particletype*)( \
		(idx) < np ? (bp_raw + (idx)*p_size) : (bpp_raw + ((idx)-np)*p_size) ))

	/* 2. Count ghost particles needed from all source particles */
	int ncyl_ghost = 0, ntop_ghost = 0, nbot_ghost = 0;
	int nleft_ghost = 0, nright_ghost = 0;
	int ncorner_ghost = 0;
	for(i=0;i<nsrc;i++){
		treevorork4particletype *bpi = SRC_PTL(i);
		postype dx = bpi->x - cx;
		postype dy = bpi->y - cy;
		postype r = sqrt(dx*dx + dy*dy);
		if(r > R && r < R + 2*width) ncyl_ghost++;
		int near_bot = (bpi->y < width && bpi->y >= 0);
		int near_top = (bpi->y > Ly - width && bpi->y <= Ly);
		int near_left = (bpi->x < width && bpi->x >= 0);
		int near_right = (bpi->x > Lx - width && bpi->x <= Lx);
		nbot_ghost += near_bot;
		ntop_ghost += near_top;
		nleft_ghost += near_left;
		nright_ghost += near_right;
		/* Corner ghosts: particles near two boundaries */
		if(near_bot && near_left) ncorner_ghost++;
		if(near_bot && near_right) ncorner_ghost++;
		if(near_top && near_left) ncorner_ghost++;
		if(near_top && near_right) ncorner_ghost++;
	}

	int nghost_total = ncyl_ghost + ntop_ghost + nbot_ghost + nleft_ghost + nright_ghost + ncorner_ghost;

	/* 3. Reallocate padding buffer (use my_malloc to avoid heap corruption with -O2) */
	{
		size_t new_total = npad_mpi + nghost_total + 100;
		char *new_raw = (char*)my_malloc(p_size * new_total);
		if(npad_mpi > 0 && bpp_raw)
			memcpy(new_raw, bpp_raw, p_size * npad_mpi);
		if(bpp_raw) my_free(bpp_raw);
		bpp_raw = new_raw;
		VORORK4_TBPP(simpar) = (treevorork4particletype*)bpp_raw;
	}
	bp_raw = (char*)VORORK4_TBP(simpar);

	int gidx = npad_mpi;
	/* All boundary ghosts must have PINDX == MAX_INDEX so that
	   the force loop skips AV and applies CFL floor for them */

	/* 4. Create cylinder reflection ghosts (free-slip) */
	for(i=0;i<nsrc;i++){
		treevorork4particletype *bpi = SRC_PTL(i);
		postype dx = bpi->x - cx;
		postype dy = bpi->y - cy;
		postype r = sqrt(dx*dx + dy*dy);
		if(r > R && r < R + 2*width){
			treevorork4particletype *ghost = (treevorork4particletype*)(bpp_raw + gidx*p_size);
			memcpy(ghost, bpi, p_size);
			postype nx = dx/r;
			postype ny = dy/r;
			postype r_mirror = 2*R - r;
			ghost->x = cx + nx*r_mirror;
			ghost->y = cy + ny*r_mirror;
			postype vn = bpi->vx*nx + bpi->vy*ny;
			ghost->vx = bpi->vx - 2*vn*nx;
			ghost->vy = bpi->vy - 2*vn*ny;
			ghost->u4if.indx = MAX_INDEX;
			ghost->u4if.Flag[ENDIAN_OFFSET] |= BoundaryGhostflag;
			gidx++;
		}
	}

	/* 5. Bottom periodic ghost (y < width, y >= 0) — copy shifted up by Ly */
	for(i=0;i<nsrc;i++){
		treevorork4particletype *bpi = SRC_PTL(i);
		if(bpi->y < width && bpi->y >= 0){
			treevorork4particletype *ghost = (treevorork4particletype*)(bpp_raw + gidx*p_size);
			memcpy(ghost, bpi, p_size);
			ghost->y = bpi->y + Ly;
			/* Periodic ghost: keep original PINDX so force loop
			   treats it as a normal neighbor (AV on, no CFL floor) */
			gidx++;
		}
	}

	/* 6. Top periodic ghost (y > Ly - width, y <= Ly) — copy shifted down by Ly */
	for(i=0;i<nsrc;i++){
		treevorork4particletype *bpi = SRC_PTL(i);
		if(bpi->y > Ly - width && bpi->y <= Ly){
			treevorork4particletype *ghost = (treevorork4particletype*)(bpp_raw + gidx*p_size);
			memcpy(ghost, bpi, p_size);
			ghost->y = bpi->y - Ly;
			/* Periodic ghost: keep original PINDX */
			gidx++;
		}
	}

	/* 7. Left inflow ghost (x < width, x >= 0) — fixed freestream state */
	{
		postype ie_per_vol = p_inf / (Gamma-1);  /* P/(γ-1): ie per unit volume */
		for(i=0;i<nsrc;i++){
			treevorork4particletype *bpi = SRC_PTL(i);
			if(bpi->x < width && bpi->x >= 0){
				treevorork4particletype *ghost = (treevorork4particletype*)(bpp_raw + gidx*p_size);
				memcpy(ghost, bpi, p_size);
				ghost->x = -bpi->x;
				ghost->vx = u_inf;
				ghost->vy = 0;
				ghost->den = rho_inf;
				ghost->pressure = p_inf;
				ghost->csound = cs_inf;
				ghost->ie = ie_per_vol * ghost->volume;
				ghost->u4if.indx = MAX_INDEX;
				ghost->u4if.Flag[ENDIAN_OFFSET] |= BoundaryGhostflag;
				gidx++;
			}
		}
	}

	/* 8. Right outflow ghost (x > Lx - width, x <= Lx) — extrapolation (copy nearest) */
	for(i=0;i<nsrc;i++){
		treevorork4particletype *bpi = SRC_PTL(i);
		if(bpi->x > Lx - width && bpi->x <= Lx){
			treevorork4particletype *ghost = (treevorork4particletype*)(bpp_raw + gidx*p_size);
			memcpy(ghost, bpi, p_size);
			ghost->x = 2*Lx - bpi->x;
			ghost->u4if.indx = MAX_INDEX;
			ghost->u4if.Flag[ENDIAN_OFFSET] |= BoundaryGhostflag;
			gidx++;
		}
	}

	/* 9. Corner ghosts (double reflections for particles near two boundaries) */
	{
		postype ie_per_vol = p_inf / (Gamma-1);
		for(i=0;i<nsrc;i++){
			treevorork4particletype *bpi = SRC_PTL(i);
			int near_bot = (bpi->y < width && bpi->y >= 0);
			int near_top = (bpi->y > Ly - width && bpi->y <= Ly);
			int near_left = (bpi->x < width && bpi->x >= 0);
			int near_right = (bpi->x > Lx - width && bpi->x <= Lx);
			if(near_bot && near_left){
				treevorork4particletype *ghost = (treevorork4particletype*)(bpp_raw + gidx*p_size);
				memcpy(ghost, bpi, p_size);
				ghost->x = -bpi->x;
				ghost->y = bpi->y + Ly;
				ghost->vx = u_inf;
				ghost->vy = bpi->vy;
				ghost->den = rho_inf;
				ghost->pressure = p_inf;
				ghost->csound = cs_inf;
				ghost->ie = ie_per_vol * ghost->volume;
				ghost->u4if.indx = MAX_INDEX;
				ghost->u4if.Flag[ENDIAN_OFFSET] |= BoundaryGhostflag;
				gidx++;
			}
			if(near_bot && near_right){
				treevorork4particletype *ghost = (treevorork4particletype*)(bpp_raw + gidx*p_size);
				memcpy(ghost, bpi, p_size);
				ghost->x = 2*Lx - bpi->x;
				ghost->y = bpi->y + Ly;
				ghost->u4if.indx = MAX_INDEX;
				ghost->u4if.Flag[ENDIAN_OFFSET] |= BoundaryGhostflag;
				gidx++;
			}
			if(near_top && near_left){
				treevorork4particletype *ghost = (treevorork4particletype*)(bpp_raw + gidx*p_size);
				memcpy(ghost, bpi, p_size);
				ghost->x = -bpi->x;
				ghost->y = bpi->y - Ly;
				ghost->vx = u_inf;
				ghost->vy = bpi->vy;
				ghost->den = rho_inf;
				ghost->pressure = p_inf;
				ghost->csound = cs_inf;
				ghost->ie = ie_per_vol * ghost->volume;
				ghost->u4if.indx = MAX_INDEX;
				ghost->u4if.Flag[ENDIAN_OFFSET] |= BoundaryGhostflag;
				gidx++;
			}
			if(near_top && near_right){
				treevorork4particletype *ghost = (treevorork4particletype*)(bpp_raw + gidx*p_size);
				memcpy(ghost, bpi, p_size);
				ghost->x = 2*Lx - bpi->x;
				ghost->y = bpi->y - Ly;
				ghost->u4if.indx = MAX_INDEX;
				ghost->u4if.Flag[ENDIAN_OFFSET] |= BoundaryGhostflag;
				gidx++;
			}
		}
	}

	#undef SRC_PTL

	/* 10. Update total padding count */
	VORO_NPAD(simpar) = gidx;
/*	DEBUGPRINT("P%d paddingCyl: npad_mpi=%d ncyl=%d nbot=%d ntop=%d nleft=%d nright=%d ncorner=%d total=%d\n",
		MYID(simpar), npad_mpi, ncyl_ghost, nbot_ghost, ntop_ghost, nleft_ghost, nright_ghost, ncorner_ghost, gidx); */
}


void cyl_findVol(SimParameters *simpar){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	det2d_dpqRK4(simpar, paddingCylinderParticles);

	postype xmin,ymin,xmax,ymax, cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar) = CYL_GridSize(simpar);
	xmin = CYL_XMIN(simpar)-cellsize;
	ymin = CYL_YMIN(simpar)-cellsize;
	xmax = CYL_XMAX(simpar)+cellsize;
	ymax = CYL_YMAX(simpar)+cellsize;

	int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
	int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);
	CellType *cells = VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*mx*my);
	mkLinkedList2D_oExam(simpar, cellsize,xmin,ymin,xmax,ymax, paddingCylinderParticles);

	int iy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(iy=0;iy<my;iy++){
		int mp=2000;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		int ix;
		for(ix=0;ix<mx;ix++){
			int np;
			treevorork4particletype *p = findCellRk4BP2D(simpar,ix,iy,&np);
			if(!p || np == 0) continue;
			int nneigh;
			Voro2D_point *neighbors = searchCellRk4Neighbors2D(simpar,ix,iy,&nneigh);
			Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			int i;
			for(i=0;i<np;i++){
				Voro2D_point center;
				center.x = p[i].x;
				center.y = p[i].y;
				center.indx = PINDX(p+i);
				center.csound = p[i].csound;
				center.w2 = p[i].w2;

				treevorork4particletype *ibp = p[i].bp;
				ibp->dt = 1.e10;
				int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
				ibp->volume = Area2DPolygon(vorocorner,mp);
				if(ibp->volume <= 0) ibp->volume = 1e-30;
				ibp->den = ibp->mass / ibp->volume;
			}
			free(neighwork);
			free(neighbors);
			free(p);
		}
		free(vorocorner);
	}
	free(VORO_BASICCELL(simpar));

	VORO_NPAD(simpar) = 0;
	free(VORORK4_TBPP(simpar));
	VORORK4_TBPP(simpar) = NULL;
}
