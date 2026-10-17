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

	/* Rings fit in 2*nx*ny. Fixed-angle rays use one N for every shell, so the
	   count is N * k_max with N ~ π nx. */
	const char *ic_env_cap = getenv("SEDOV2D_IC");
	int rays_cap = (ic_env_cap && strcmp(ic_env_cap, "rays") == 0);
	size_t buf_capacity = (size_t)((rays_cap ? 6 : 2) * nx * ny);
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

	/* IC mode selector via env var SEDOV2D_IC.
	     "brick"     — legacy nx*ny Cartesian (mode-4 lattice imprint)
	     "isoshell"  — equal-Δr concentric shells with random per-shell phase.
	     "sunflower" — Vogel spiral r_i=c√i, θ_i=i·φ_g (golden angle).  No
	                   radial discreteness; SO(2) symmetric in continuum limit.
	     "glass"     — random uniform + Lloyd relaxation (amorphous packing).
	                   Controls: SEDOV2D_GLASS_SEED (default 42),
	                   SEDOV2D_GLASS_NLLOYD (default 100).
	     "rays"      — one angular set for every radius. θ_j=(j+1/2)2π/N,
	                   mass ∝ r so the area density stays uniform. No ray on
	                   a coordinate axis when that angle is not in the set.
	     default     — concentric rings with alternating phase. */
	const char *ic_env = getenv("SEDOV2D_IC");
	int use_brick_ic     = (ic_env && strcmp(ic_env, "brick") == 0);
	int use_isoshell_ic  = (ic_env && strcmp(ic_env, "isoshell") == 0);
	int use_sunflower_ic = (ic_env && strcmp(ic_env, "sunflower") == 0);
	int use_glass_ic     = (ic_env && strcmp(ic_env, "glass") == 0);
	int use_rays_ic      = (ic_env && strcmp(ic_env, "rays") == 0);
	int use_shell_ic     = !use_brick_ic && !use_isoshell_ic
	                       && !use_sunflower_ic && !use_glass_ic
	                       && !use_rays_ic;

	/* Concentric rings IC: place particles on rings r_k = k*dr around (cx,cy),
	   with N_k = round(2π·r_k/dr) particles per ring (uniform area density).
	   Alternating angular phase between rings breaks radial spoke alignment.
	   Particles outside the box are clipped. The central anchor is placed last. */
	size_t np = 0;
	if(rank0_only && use_shell_ic){
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

	/* Fixed-angle rays: every shell uses the same θ_j. Mass ∝ r. */
	if(rank0_only && use_rays_ic){
		postype dr = Lx/(postype)nx;
		int N = (int)(M_PI * (postype)nx + 0.5);
		if(N < 8) N = 8;
		postype dtheta = 2.0*M_PI/(postype)N;
		postype r_corner = sqrt((cx-xmin)*(cx-xmin) + (cy-ymin)*(cy-ymin));
		postype r2;
		r2 = sqrt((xmax-cx)*(xmax-cx) + (cy-ymin)*(cy-ymin)); if(r2 > r_corner) r_corner = r2;
		r2 = sqrt((cx-xmin)*(cx-xmin) + (ymax-cy)*(ymax-cy)); if(r2 > r_corner) r_corner = r2;
		r2 = sqrt((xmax-cx)*(xmax-cx) + (ymax-cy)*(ymax-cy)); if(r2 > r_corner) r_corner = r2;
		int k_max = (int)(r_corner/dr) + 2;
		int k, j;
		for(k=1; k<=k_max; k++){
			postype rk = (postype)k * dr;
			postype mass_k = rho_amb * rk * dr * dtheta;
			for(j=0; j<N; j++){
				postype theta = ((postype)j + 0.5) * dtheta;
				postype x = cx + rk*cos(theta);
				postype y = cy + rk*sin(theta);
				if(x < xmin || x >= xmax || y < ymin || y >= ymax) continue;
				if(np + 1 >= buf_capacity){
					fprintf(stderr,"[Sedov2D-RAYS] buffer full np=%zu\n", np);
					break;
				}
				UNSET_FLAG(res+np, Wallflag);
				res[np].u4if.indx = (size_t)nx*(size_t)ny + 1 + np;
				res[np].x = x;
				res[np].y = y;
				res[np].vx = 0; res[np].vy = 0;
				res[np].z = 0; res[np].vz = 0;
				res[np].ax = 0; res[np].ay = 0; res[np].az = 0;
				res[np].mass = mass_k;
				res[np].den = rho_amb;
				res[np].pressure = P_amb;
				res[np].ie = P_amb * (mass_k/rho_amb) / (Gamma-1);
				res[np].csound = sqrt(Gamma*P_amb/rho_amb);
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
		fprintf(stderr,"[Sedov2D-RAYS] N=%d k_max=%d np=%zu dtheta_deg=%.4f\n",
			N, k_max, np, (double)(dtheta*180.0/M_PI)); fflush(stderr);
	}

	/* Isoshell IC: stratified concentric shells r_k = k·Δr with N_k particles
	   on each shell at angular positions θ_{k,j} = ψ_k + (j+0.5)·Δθ_k, where
	   ψ_k ∼ U[0, Δθ_k) is an independent random phase per shell.  Continuum
	   limit: each shell has C_∞ symmetry on average; independent phases give
	   no shell-to-shell alignment → SO(2) symmetric overall.  Optional small
	   radial jitter (fraction of Δr) further breaks the per-shell rotational
	   discreteness. */
	if(rank0_only && use_isoshell_ic){
		const char *seed_env   = getenv("SEDOV2D_ISOSHELL_SEED");
		const char *rjit_env   = getenv("SEDOV2D_ISOSHELL_RJITTER");
		unsigned int seed      = seed_env ? (unsigned int)atoi(seed_env) : 42u;
		postype rjit_frac      = rjit_env ? (postype)atof(rjit_env) : 0.0;
		srand(seed);
		postype dx_ic = Lx/(postype)nx;
		postype dr = dx_ic;
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
			/* per-shell random angular phase, uniform on [0, Δθ_k) */
			postype psi = ((postype)rand()/(postype)RAND_MAX) * dtheta;
			int m;
			for(m=0; m<Nk; m++){
				postype theta = psi + ((postype)m + 0.5) * dtheta;
				postype rjit = 0.0;
				if(rjit_frac > 0.0)
					rjit = dr * rjit_frac * (((postype)rand()/(postype)RAND_MAX) - 0.5);
				postype r_eff = rk + rjit;
				postype x = cx + r_eff*cos(theta);
				postype y = cy + r_eff*sin(theta);
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
		/* Anchor at exact center (cx,cy) */
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

		fprintf(stderr,"[Sedov2D-ISOSHELL] seed=%u rjitter=%.3f k_max=%d np=%zu (target nx*ny+1=%d, ratio=%.4f)\n",
			seed, (double)rjit_frac, k_max, np, nx*ny+1, (double)np/(double)(nx*ny+1)); fflush(stderr);
	}

	/* Sunflower (Vogel spiral) IC: r_i = c·√i, θ_i = i·φ_g where
	   φ_g = 2π/φ² is the golden angle.  No radial discreteness, no preferred
	   azimuthal axis (φ_g is the "most irrational" angle → broadband mode
	   amplitudes are minimized).  i=0 lands at origin → natural anchor.
	   Box-filling: oversample by π/4 ratio (fraction of circle inscribed in
	   the bounding box covered by the spiral disk). */
	if(rank0_only && use_sunflower_ic){
		postype dx_ic = Lx/(postype)nx;
		postype meanvol_ic = dx_ic*dx_ic;
		/* Scale c so that the area per particle is meanvol_ic.
		   For Vogel spiral, particle i sits in an annulus ~Δr_i ≈ c/(2√i),
		   with arc per particle ≈ 2π r_i / N_local.  Density is uniform when
		   c² = meanvol_ic/π, i.e. c = √(meanvol_ic/π). */
		postype c_sf = sqrt(meanvol_ic / M_PI);
		/* Largest radius needed: box corner farthest from center */
		postype r_corner = sqrt((cx-xmin)*(cx-xmin) + (cy-ymin)*(cy-ymin));
		postype r2;
		r2 = sqrt((xmax-cx)*(xmax-cx) + (cy-ymin)*(cy-ymin)); if(r2 > r_corner) r_corner = r2;
		r2 = sqrt((cx-xmin)*(cx-xmin) + (ymax-cy)*(ymax-cy)); if(r2 > r_corner) r_corner = r2;
		r2 = sqrt((xmax-cx)*(xmax-cx) + (ymax-cy)*(ymax-cy)); if(r2 > r_corner) r_corner = r2;
		size_t N_total = (size_t)((r_corner*r_corner)/(c_sf*c_sf)) + 4;
		postype phi_g = M_PI * (3.0 - sqrt(5.0));  /* golden angle = 2π(2-φ) */

		size_t i_sf;
		for(i_sf=1; i_sf<N_total; i_sf++){
			postype r_i  = c_sf * sqrt((postype)i_sf);
			postype th_i = (postype)i_sf * phi_g;
			postype x = cx + r_i*cos(th_i);
			postype y = cy + r_i*sin(th_i);
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
		/* Anchor at exact center (i=0 of the Vogel spiral) */
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

		fprintf(stderr,"[Sedov2D-SUNFLOWER] c=%g N_total=%zu np=%zu (target nx*ny+1=%d, ratio=%.4f)\n",
			(double)c_sf, N_total, np, nx*ny+1, (double)np/(double)(nx*ny+1)); fflush(stderr);
	}

	/* Glass IC: brick lattice + small uniform jitter as seed, then Lloyd.
	   Starting from a regular lattice gives Lloyd an O(1) initial V_max/V_min
	   (vs ~5-10 for pure random); convergence to V_max/V_min ~ 1 is then fast.
	   Jitter amplitude (fraction of Δx) controlled by SEDOV2D_GLASS_JITTER
	   (default 0.3 — large enough to break the lattice, small enough that
	   Lloyd doesn't have to swap neighbors). */
	if(rank0_only && use_glass_ic){
		const char *seed_env = getenv("SEDOV2D_GLASS_SEED");
		unsigned int seed    = seed_env ? (unsigned int)atoi(seed_env) : 42u;
		const char *jit_env  = getenv("SEDOV2D_GLASS_JITTER");
		postype jit_frac     = jit_env ? (postype)atof(jit_env) : (postype)0.3;
		srand(seed);

		postype dx_ic = Lx/(postype)nx;
		postype dy_ic = Ly/(postype)ny;
		int ix_g, iy_g;
		for(iy_g=0; iy_g<ny; iy_g++){
		  for(ix_g=0; ix_g<nx; ix_g++){
			postype rx = ((postype)rand()/(postype)RAND_MAX) - (postype)0.5;
			postype ry = ((postype)rand()/(postype)RAND_MAX) - (postype)0.5;
			postype x = xmin + (ix_g + (postype)0.5)*dx_ic + jit_frac*dx_ic*rx;
			postype y = ymin + (iy_g + (postype)0.5)*dy_ic + jit_frac*dy_ic*ry;
			if(x < xmin) x += Lx; else if(x >= xmax) x -= Lx;
			if(y < ymin) y += Ly; else if(y >= ymax) y -= Ly;

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

		fprintf(stderr,"[Sedov2D-GLASS] seed=%u jitter=%.3f np=%zu (target nx*ny+1=%d, ratio=%.4f) - Lloyd to follow\n",
			seed, (double)jit_frac, np, nx*ny+1, (double)np/(double)(nx*ny+1)); fflush(stderr);
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
	   With shell/sunflower ICs particles are already at uniform density and
	   isotropic; Lloyd would distort them into hex patches and break SO(2).
	   Glass IC starts from uniform random points and *requires* Lloyd to
	   converge to an amorphous packing — controlled by SEDOV2D_GLASS_NLLOYD. */
	{
		int n_glass = 0;
		if(use_brick_ic) n_glass = 100;
		if(use_glass_ic){
			const char *nlloyd_env = getenv("SEDOV2D_GLASS_NLLOYD");
			n_glass = nlloyd_env ? atoi(nlloyd_env) : 100;
		}
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
					fprintf(stderr,"[Sedov2D-LLOYD] it=%d Vmin=%g Vmax=%g ratio=%g\n",
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
		/* Anchor immobility is the position pin in walls_xy_postStage_blend.
		 * A 1000x mass makes rho=m/V a permanent density spike. LagMFM and
		 * the phase-1 HLL probe both keep the anchor at the ambient mass. */
		int av_mode_now = GAS_AVMODE(simpar);
		const char *phase1 = getenv("SEDOV_PHASE1");
		int phase1_on = (phase1 && phase1[0] == '1');
		postype M_anchor = (av_mode_now == 4 || phase1_on) ? rho_amb*meanvol
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
