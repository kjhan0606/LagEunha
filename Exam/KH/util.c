#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "voro.h"
//#include "ost.h"
#include "nnost.h"
#include "gnnost.h"
#include "exam.h"
#include "exam2d.h"
#include "kh.h"
#include "color.h"




static CellType *cells=NULL;

//postype cellsize = Ly/Ny;

//void colorizeit(float *, int, int, char *);

int kh_makemap(SimParameters *simpar, int icount){
	postype cellsize = KH_GridSize(simpar);
	postype Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	postype Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	int nximg = NX(simpar);
	int nyimg = NY(simpar);
	float *map = (float*)malloc(sizeof(float)*nximg*nyimg);
	float *img = (float*)malloc(sizeof(float)*nximg*nyimg);
	postype meanvol = Lx*Ly/nximg/nyimg;
//	float *ndist = (float*)malloc(sizeof(float)*nximg*nyimg);
	int i,j;
	int ii,jj;
	for(i=0;i<nximg*nyimg;i++) map[i] = KH_Rho1(simpar);
//	for(i=0;i<nximg*nyimg;i++) ndist[i] = 1.e20;
	postype pixsize = Lx/nximg;
	postype xmin,ymin,xmax,ymax;
	xmin = KH_XMIN(simpar);
	ymin = KH_YMIN(simpar);
	xmax = KH_XMAX(simpar);
	ymax = KH_YMAX(simpar);
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
		int index = ix+mx*iy;
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)bpi;
		cells[index].nmem ++;
		bpi->next = tmp;
	}

	for(j=0;j<nyimg;j++){
		postype yp = (j+0.5)*pixsize; // img pixel position
		if(yp < ymin || yp >= ymax) continue;
		for(i=0;i<nximg;i++){
			postype xp = (i+0.5)*pixsize; // img pixel position
			if(xp < xmin || xp >= xmax) continue;
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
								nearden = tt->mass/meanvol;
							}
						}
						tmp = tmp->next;
					}

				}
			}
//			map[i+nximg*j] = nearden;
			map[i+nximg*(nyimg-j-1)] = nearden;
//			ndist[i+nximg*j] = idist;

		}
	}
	/*
	int id;
	if(MYID(simpar)==0){
		float *idist = (float*)malloc(sizeof(float)*nximg*nyimg);
		MPI_Status status;

		for(id=1;id<NID(simpar);id++){
			MPI_Recv(idist, nximg*nyimg,MPI_FLOAT,id,id, MPI_COMM_WORLD, &status);
			MPI_Recv(img, nximg*nyimg,MPI_FLOAT,id,id, MPI_COMM_WORLD, &status);
			for(j=0;j<nximg*nyimg;j++){
				if(idist[j] < ndist[j]){
					map[j] = img[j];
					ndist[j] = idist[j];
				}
			}
		}
		free(idist);
	}
	else {
		MPI_Send(ndist, nximg*nyimg,MPI_FLOAT,0,MYID(simpar), MPI_COMM_WORLD);
		MPI_Send(map, nximg*nyimg,MPI_FLOAT,0,MYID(simpar), MPI_COMM_WORLD);
	}

	if(MYID(simpar)==0){
		float *img2 = (float*)malloc(sizeof(float)*nximg*nyimg);;
		for(j=0;j<nyimg;j++){
			for(i=0;i<nximg;i++){
				img2[i+nximg*(nyimg-j-1)] = map[i+nximg*j];
			}
		}
		char outfile[189]; 
		sprintf(outfile,"khmap.%.6d.ppm", icount); 
		colorizeit(7,img2, nximg,nyimg,"kh.sao", outfile, KH_Rho1(simpar), KH_Rho2(simpar));
		free(img2);
	}
	free(cells);
 	free(map); free(img);
	free(ndist);
	*/
	MPI_Reduce(map, img, nximg*nyimg, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);


    if(MYID(simpar)==0){
        char outfile[189];
        sprintf(outfile,"khmap.%.6d.ppm", icount);
        colorizeit(7, img, nximg,nyimg,"kh.sao", outfile, (double)KH_Rho1(simpar), (double)KH_Rho2(simpar));
    }
    free(cells);
    free(map); free(img);
	return 0;
}
int kh_maketscmap(SimParameters *simpar, int icount){
	int nx = NX(simpar);
	int ny = NY(simpar);
	int np = VORO_NP(simpar);
	int nximg, nyimg;
	nximg = nx;
	nyimg = ny;
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	float Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	float Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	float cellsize = Lx/nximg;
	int i;
	float *map = (float*)malloc(sizeof(float)*nximg*nyimg);
	float *tmap = (float*)malloc(sizeof(float)*nximg*nyimg);
	for(i=0;i<nximg*nyimg;i++) map[i] = 0;
	for(i=0;i<nximg*nyimg;i++) tmap[i] = 0;
	postype mscale = (nximg/Lx)*(nyimg/Ly);

	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		postype pmass = mscale*bpi->mass;
		postype xp = bpi->x/cellsize;
		postype yp = bpi->y/cellsize;
		int nearx = rint(xp);
		int neary = rint(yp);
		int ic = (nearx + nximg) %nximg;
		int jc = (neary + nyimg) %nyimg;
		postype xmin = xp -nearx;
		postype ymin = yp -neary;
		int icc = (nximg+nearx + (int)(copysign(1., xmin))) %nximg;
		int jcc = (nyimg+neary + (int)(copysign(1., ymin))) %nyimg;
		int iccc = (nximg+nearx - (int)(copysign(1., xmin))) %nximg;
		int jccc = (nyimg+neary - (int)(copysign(1., ymin))) %nyimg;
		float xd1 = fabs(xmin);
		float yd1 = fabs(ymin);
		float wx1 = (0.75-xd1*xd1)*pmass;
		float wy1 = (0.75-yd1*yd1);
		float wx3 = 0.5*pmass * (0.25+xd1*(xd1-1.));
		float wx2 = wx3 + pmass * xd1;
		float wy3 = 0.5 * (0.25+yd1*(yd1-1.));
		float wy2 = wy3 + yd1;

		map[ic+nximg*jc] += wx1*wy1;
		map[icc+nximg*jc] += wx2*wy1;
		map[iccc+nximg*jc] += wx3*wy1;
		map[ic+nximg*jcc] += wx1*wy2;
		map[icc+nximg*jcc] += wx2*wy2;
		map[iccc+nximg*jcc] += wx3*wy2;
		map[ic+nximg*jccc] += wx1*wy3;
		map[icc+nximg*jccc] += wx2*wy3;
		map[iccc+nximg*jccc] += wx3*wy3;

	}
	MPI_Reduce(map, tmap, nximg*nyimg, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	if(MYID(simpar)==0){
		char outfile[189];
		sprintf(outfile,"khmap.%.6d.ppm", icount);
		colorizeit(7,tmap, nximg,nyimg,"kh.sao", outfile, KH_Rho1(simpar), KH_Rho2(simpar));
	}
	free(map);
	free(tmap);

	return 0;
}

treevorork4particletype *kh_mkinitial(SimParameters *simpar, int *mp){
	int i,j,k;
	treevorork4particletype *res;
	postype rho1 = KH_Rho1(simpar), rho2 = KH_Rho2(simpar);
    postype deltay = KH_Deltay(simpar);  // tanh width 'a'
    postype dvy0 = KH_Vperturb(simpar);
    postype U1 = KH_Vel1(simpar);       // outer velocity
    postype U2 = KH_Vel2(simpar);       // inner velocity
    postype xmin,ymin,xmax,ymax;
    xmin = KH_XMIN(simpar);
    ymin = KH_YMIN(simpar);
    xmax = KH_XMAX(simpar);
    ymax = KH_YMAX(simpar);
    int nx = NX(simpar);
    int ny = NY(simpar);
    float Lx = KH_SIMBOX(simpar).x.max - KH_SIMBOX(simpar).x.min;
    float Ly = KH_SIMBOX(simpar).y.max - KH_SIMBOX(simpar).y.min;
	postype dmean = Lx/nx;
	GAS_dMean(simpar) = dmean;
	int av_mode = GAS_AVMODE(simpar);
#ifdef USE_CUDA
	int use_stress = 1;  /* GPU blend path needs stress-type for all av_modes */
#else
	int use_stress = (av_mode >= 1 || GAS_VISCOSITY(simpar) > 0);
#endif

	postype Gamma = GAS_GAMMA(simpar);

    if(use_stress)
        res = (treevorork4particletype*)my_malloc(sizeof(treevorostressrk4particletype)*nx*ny);
    else
        res = (treevorork4particletype*)my_malloc(sizeof(treevorork4particletype)*nx*ny);
    postype meanvol = Lx*Ly/nx/ny;
    size_t np = 0;
	// Lecoanet et al. (2016) Section 3.2 IC
	postype z1 = 0.5, z2 = 1.5;      // interface positions (fixed for [0,2] domain)
	postype sigma = 0.2;              // Gaussian envelope width for vy perturbation

	// Glass IC loading.  If a binary file at ./glass_NxxNy.bin (or pointed to
	// by env var LAGMFM_GLASS_FILE) exists with matching np=nx*ny, Lx, Ly,
	// use its positions; else fall back to regular grid + per-row McNally.
	double *glass_x = NULL, *glass_y = NULL;
	int use_glass = 0;
	{
		char path[512];
		const char *env = getenv("LAGMFM_GLASS_FILE");
		if(env && *env) snprintf(path, sizeof(path), "%s", env);
		else            snprintf(path, sizeof(path), "./glass_%dx%d.bin", nx, ny);
#if defined(IC_SOD2D)
		/* Sod with reflecting x-walls needs regular grid IC: glass irregularity
		 * puts some particles within << dx/2 of x=0/Lx, and their wall mirrors
		 * then create tiny Voronoi cells → density spikes. */
		FILE *fp = NULL;
#else
		FILE *fp = fopen(path, "rb");
#endif
		if(fp){
			int np_file;
			double Lx_file, Ly_file;
			if(fread(&np_file, sizeof(int), 1, fp)==1 &&
			   fread(&Lx_file, sizeof(double), 1, fp)==1 &&
			   fread(&Ly_file, sizeof(double), 1, fp)==1 &&
			   np_file == nx*ny &&
			   fabs(Lx_file-Lx) < 1e-9 && fabs(Ly_file-Ly) < 1e-9){
				glass_x = (double*)malloc(sizeof(double)*np_file);
				glass_y = (double*)malloc(sizeof(double)*np_file);
				if(fread(glass_x, sizeof(double), np_file, fp)==(size_t)np_file &&
				   fread(glass_y, sizeof(double), np_file, fp)==(size_t)np_file){
					use_glass = 1;
					DEBUGPRINT("KH: loaded glass IC from %s (np=%d)\n", path, np_file);
				} else {
					free(glass_x); free(glass_y); glass_x=glass_y=NULL;
				}
			} else {
				DEBUGPRINT("KH: glass file %s header mismatch, falling back to grid\n", path);
			}
			fclose(fp);
		} else {
			DEBUGPRINT("KH: glass file %s not found, using regular grid\n", path);
		}
	}

	size_t n_particles_total = (size_t)nx*(size_t)ny;
	for(size_t p=0; p<n_particles_total; p++){
		postype x, y;
		if(use_glass){
			x = (postype)glass_x[p];
			y = (postype)glass_y[p];
			i = (int)(p % (size_t)nx);
			j = (int)(p / (size_t)nx);
		} else {
			i = (int)(p % (size_t)nx);
			j = (int)(p / (size_t)nx);
			x = (postype)(i+0.5)*Lx/(postype)nx;
			y = (postype)(j+0.5)*Ly/(postype)ny;
		}
        postype rho,vx,press;
        postype vy_seed = 0;
        int keep = 1;   /* EUNHA_KEPLER_PROFILE=gizmo drops particles outside the ring */
#if defined(IC_SOD2D)
        /* 2D Sod shock tube: step at x=Lx/2.  For periodic BC, use wide box
         * (Lx ≥ 4) so waves don't wrap around to T_end=0.2.  Analytical
         * comparison valid in [Lx/2 - 1.2·T, Lx/2 + 1.8·T] roughly. */
        if(x < 0.5*Lx){ rho = 1.0;   press = 1.0; }
        else          { rho = 0.125; press = 0.1; }
        vx = 0;
        char iregion = (x < 0.5*Lx) ? 0 : 1;
#elif defined(IC_SEDOV2D)
        /* 2D Sedov (GIZMO/Wengen-style): top-hat energy over N_hot=64 central
         * particles, not a 1-cell spike.  P_hot = E(gamma-1)/(N_hot*h^2) so
         * the total deposited energy is E=1.  R_hot = sqrt(N_hot/pi)*h gives
         * a disk that contains exactly N_hot particles on a regular grid.
         *   N=64, gamma=1.4, h=L/nx  ->  P_hot ~= 102, ratio P_hot/P_amb=1e5.
         * The previous single-cell step (P_hot=500, ratio 5e5, N=12) made the
         * Laguerre cell construction degenerate and produced grid-anisotropy
         * in the blast. */
        rho   = 1.0;
        vx    = 0;
        {
            postype dx_c = x - 0.5*Lx, dy_c = y - 0.5*Ly;
            postype h0 = Lx/(postype)nx;
            postype R_hot = 4.51 * h0;    /* sqrt(64/pi) ~ 4.51 */
            /* P_hot = E*(gamma-1)/(N_hot*h^2) so total deposited E=1 for any gamma */
            postype P_hot = (Gamma - 1.0) / (64.0 * h0 * h0);
            if(dx_c*dx_c + dy_c*dy_c < R_hot*R_hot)
                press = P_hot;
            else
                press = 1.0e-3;
        }
        char iregion = 0;
#elif defined(IC_NOH2D)
        /* 2D Noh: radial inflow v = -r̂, ρ=1, P=1e-6. */
        rho = 1.0;
        press = 1.0e-6;
        postype dx_c = x - 0.5*Lx, dy_c = y - 0.5*Ly;
        postype r   = sqrt(dx_c*dx_c + dy_c*dy_c) + 1e-30;
        vx = -dx_c/r;
        vy_seed = -dy_c/r;
        char iregion = 0;
#elif defined(IC_DMR2D)
        /* Simplified DMR-like: Mach 10 planar oblique shock in periodic box.
         * Real DMR needs reflective bottom wall; this variant just tests
         * HLLC handling of Mach 10 shock state before wall interaction.
         *  Pre-shock (ahead of shock):  ρ=1.4, P=1, v=0
         *  Post-shock:  ρ=8, P=116.5, u=7.14, v=-4.13
         *  Initial shock front at x = 1/6 + y/√3 */
        {
            postype shock_x = Lx/6.0 + y / sqrt(3.0);
            if(x < shock_x){
                rho = 8.0; press = 116.5; vx = 7.145; vy_seed = -4.125;
            } else {
                rho = 1.4; press = 1.0;   vx = 0.0;   vy_seed = 0.0;
            }
        }
        char iregion = 0;
#elif defined(IC_RIEMANN2D)
        /* Schulz-Rinne config 3: all 4 shocks moving into UR quadrant.
         * Box [0,1]², partition at x=0.5, y=0.5.  (Lax & Liu 1998, Config 3.)
         *  UR (x>0.5, y>0.5): ρ=1.5,    u=0,     v=0,     P=1.5
         *  UL (x<0.5, y>0.5): ρ=0.5323, u=1.206, v=0,     P=0.3
         *  LL (x<0.5, y<0.5): ρ=0.138,  u=1.206, v=1.206, P=0.029
         *  LR (x>0.5, y<0.5): ρ=0.5323, u=0,     v=1.206, P=0.3
         * γ = 1.4, t_final ≈ 0.3 */
        {
            postype xc_ = 0.5*Lx, yc_ = 0.5*Ly;
            if(x > xc_ && y > yc_) {      /* UR */
                rho = 1.5;    vx = 0.0;    vy_seed = 0.0;   press = 1.5;
            } else if(x < xc_ && y > yc_){ /* UL */
                rho = 0.5323; vx = 1.206;  vy_seed = 0.0;   press = 0.3;
            } else if(x < xc_ && y < yc_){ /* LL */
                rho = 0.138;  vx = 1.206;  vy_seed = 1.206; press = 0.029;
            } else {                       /* LR */
                rho = 0.5323; vx = 0.0;    vy_seed = 1.206; press = 0.3;
            }
        }
        char iregion = 0;
#elif defined(IC_GRESHO2D)
        /* 2D Gresho vortex: stationary rotating vortex (Liska & Wendroff 2003).
         *  r < 0.2: u_φ = 5r      P = 5 + 25/2 r²
         *  0.2 < r < 0.4: u_φ = 2 - 5r  P = 9 + 25/2 r² - 20r + 4 ln(5r)
         *  r > 0.4: u_φ = 0       P = 3 + 4 ln 2
         *  ρ = 1, γ = 5/3 typically (low Mach). */
        rho = 1.0;
        {
            postype dx_c = x - 0.5*Lx, dy_c = y - 0.5*Ly;
            postype r = sqrt(dx_c*dx_c + dy_c*dy_c) + 1e-30;
            postype uphi;
            if(r < 0.2) {
                uphi = 5.0 * r;
                press = 5.0 + 12.5 * r * r;
            } else if(r < 0.4) {
                uphi = 2.0 - 5.0 * r;
                press = 9.0 + 12.5 * r * r - 20.0 * r + 4.0 * log(5.0 * r);
            } else {
                uphi = 0.0;
                press = 3.0 + 4.0 * log(2.0);
            }
            vx = -uphi * dy_c / r;      /* tangential: +φ̂ = (-y, x)/r */
            vy_seed = uphi * dx_c / r;
        }
        char iregion = 0;
#else
        {
        const char *eic = getenv("EUNHA_IC");
        if(eic && strcmp(eic, "gresho") == 0){
            rho = 1.0;
            postype dx_c = x - 0.5*Lx, dy_c = y - 0.5*Ly;
            postype r = sqrt(dx_c*dx_c + dy_c*dy_c) + 1e-30;
            postype uphi;
            if(r < 0.2) {
                uphi = 5.0 * r;
                press = 5.0 + 12.5 * r * r;
            } else if(r < 0.4) {
                uphi = 2.0 - 5.0 * r;
                press = 9.0 + 12.5 * r * r - 20.0 * r + 4.0 * log(5.0 * r);
            } else {
                uphi = 0.0;
                press = 3.0 + 4.0 * log(2.0);
            }
            vx = -uphi * dy_c / r;
            vy_seed = uphi * dx_c / r;
        } else if(eic && strcmp(eic, "noh") == 0){
            rho = 1.0;
            press = 1.0e-6;
            postype dx_c = x - 0.5*Lx, dy_c = y - 0.5*Ly;
            postype r = sqrt(dx_c*dx_c + dy_c*dy_c) + 1e-30;
            vx = -dx_c / r;
            vy_seed = -dy_c / r;
        } else if(eic && strcmp(eic, "tcftest") == 0){
            /* GFS_TCF unit tests (G2): rho = 1, P = P0 + bx x + by y,
             * v = (ux, uy) + random of amplitude va (deterministic in p). */
            const char *e_;
            double P0 = (e_ = getenv("TCF_T_P0")) && e_[0] ? atof(e_) : 1.0;
            double bx = (e_ = getenv("TCF_T_BX")) && e_[0] ? atof(e_) : 0.0;
            double by = (e_ = getenv("TCF_T_BY")) && e_[0] ? atof(e_) : 0.0;
            double ux = (e_ = getenv("TCF_T_UX")) && e_[0] ? atof(e_) : 0.0;
            double uy = (e_ = getenv("TCF_T_UY")) && e_[0] ? atof(e_) : 0.0;
            double va = (e_ = getenv("TCF_T_VA")) && e_[0] ? atof(e_) : 0.0;
            unsigned long long h_ = 0x9E3779B97F4A7C15ULL*((unsigned long long)p + 1);
            h_ ^= h_ >> 31; h_ *= 0xBF58476D1CE4E5B9ULL; h_ ^= h_ >> 29;
            double r1 = (double)(h_ & 0xFFFFFF)/16777216.0 - 0.5, r2 = (double)((h_ >> 24) & 0xFFFFFF)/16777216.0 - 0.5;
            rho = 1.0;
            press = P0 + bx*x + by*y;
            vx = ux + 2*va*r1; vy_seed = uy + 2*va*r2;
        } else if(eic && strcmp(eic, "sod") == 0){
            /* Sod shock tube along x in a periodic box: two Riemann problems
             * at x = Lx/4 and 3Lx/4. Inner (rho, P) = (1, 1), outer
             * (0.125, 0.1), v = 0; masses follow rho on the uniform lattice. */
            postype xr = x/Lx;
            if(xr >= 0.25 && xr < 0.75){ rho = 1.0; press = 1.0; }
            else { rho = 0.125; press = 0.1; }
            vx = 0; vy_seed = 0;
        } else if(eic && strcmp(eic, "vacslab") == 0){
            /* Gas slab into vacuum (particle version): lattice particles
             * kept for |x - Lx/2| < EUNHA_VAC_HALFW (default Lx/4), rho = 1,
             * P = 1, v = 0. Use with EUNHA_VORO_RHO. Each side is a 1D
             * rarefaction into vacuum; the front moves at 2c/(gamma-1). */
            /* EUNHA_VAC_P sets the pressure (default 1) and EUNHA_VAC_SHEAR
             * a uniform shear v_y = S (x - Lx/2) (default 0), periodic in y:
             * a cold sheared slab with two free surfaces, the Kepler inner
             * edge without orbit curvature. */
            const char *hws = getenv("EUNHA_VAC_HALFW");
            const char *vps = getenv("EUNHA_VAC_P");
            const char *vss = getenv("EUNHA_VAC_SHEAR");
            postype hw = (hws && hws[0]) ? atof(hws) : (postype)0.25*Lx;
            rho = 1.0; vx = 0;
            press = (vps && vps[0]) ? atof(vps) : (postype)1.0;
            vy_seed = ((vss && vss[0]) ? atof(vss) : (postype)0.0)*(x - (postype)0.5*Lx);
            if(fabs(x - 0.5*Lx) >= hw) keep = 0;
        } else if(eic && strcmp(eic, "kepler") == 0){
            const char *cxs = getenv("EUNHA_KEPLER_CX");
            const char *cys = getenv("EUNHA_KEPLER_CY");
            postype cx = (cxs && cxs[0]) ? atof(cxs) : (postype)0.5*Lx;
            postype cy = (cys && cys[0]) ? atof(cys) : (postype)0.5*Ly;
            postype dx_c = x - cx, dy_c = y - cy;
            postype R = sqrt(dx_c*dx_c + dy_c*dy_c);
            const char *es = getenv("EUNHA_KEPLER_EPS");
            const char *rins = getenv("EUNHA_KEPLER_RIN");
            const char *routs = getenv("EUNHA_KEPLER_ROUT");
            postype eps = (es && es[0]) ? atof(es) : (postype)0.0;
            /* Annulus Sigma falls as 1/R. Pressure is spatially constant and
             * tiny, as in the GIZMO cold Keplerian disk, so dP/dR = 0 and the
             * circular speed below is the radial equilibrium. */
            postype Rin = (rins && rins[0]) ? atof(rins) : (postype)2.0;
            postype Rout = (routs && routs[0]) ? atof(routs) : (postype)8.0;
            const char *profs = getenv("EUNHA_KEPLER_PROFILE");
            if(profs && strcmp(profs, "gizmo") == 0){
                /* Hopkins (2015, arXiv:1409.7395) Sec. 4.2.4, particle-code
                 * version: lattice particles are kept only for
                 * Rin <= R <= Rout (vacuum inside and outside), Sigma = 1.
                 * No density floor and no softened edges (those are for mesh
                 * codes, his footnote 24). P and v as below. */
                rho = 1;
                if(R < Rin || R > Rout) keep = 0;
            } else if(profs && strcmp(profs, "hopkinsp") == 0){
                /* Hopkins (2015) Sec. 4.2.4 particle-code IC: equal-mass
                 * lattice particles on -L < x, y < L about the centre
                 * (L = EUNHA_KEPLER_HALFBOX, default 2), thinned to match
                 * Eq. 34, Sigma = 0.01 + (R/Rin)^3 (R < Rin), 1, 0.01 +
                 * [1 + (R - Rout)/(0.05 Rout)]^-3 (R > Rout): a lattice point
                 * is kept when the R2 low-discrepancy value u(i, j) < Sigma
                 * (deterministic). Particle mass = meanvol (Sigma = 1). */
                const char *hbs = getenv("EUNHA_KEPLER_HALFBOX");
                postype hb = (hbs && hbs[0]) ? atof(hbs) : (postype)2.0;
                postype hw = (postype)0.05*Rout, shape;
                if(R < Rin) { postype q = R/Rin; shape = q*q*q; }
                else if(R <= Rout) shape = 1;
                else { postype q = 1 + (R-Rout)/hw; shape = 1/(q*q*q); }
                postype sig = (postype)0.01 + shape;
                double u = fmod(0.5 + 0.7548776662466927*(double)i + 0.5698402909980532*(double)j, 1.0);
                rho = 1;
                if(fabs(dx_c) >= hb || fabs(dy_c) >= hb) keep = 0;
                if(sig < 1 && u >= (double)sig) keep = 0;
            } else if(profs && strcmp(profs, "pakmor") == 0){
                /* Pakmor et al. (2016) App. B / Zier & Springel (2022): rho = 1 on Rin <= R <= Rout, else
                 * EUNHA_KEPLER_FLOOR (default 1e-5) AT REST; uniform pressure EUNHA_KEPLER_P
                 * (default (g-1)(5g/2)1e-5 = Pakmor's u). Freeze the background with EUNHA_KEPLER_FREEZE_RIN /
                 * EUNHA_KEPLER_EVOLVE_R. Gravity: EUNHA_KEPLER_SOFT=spline (exact point mass for R >= r_c). */
                const char *fls = getenv("EUNHA_KEPLER_FLOOR");
                postype hfloor = (fls && fls[0]) ? atof(fls) : (postype)1.0e-5;
                rho = (R >= Rin && R <= Rout) ? (postype)1.0 : hfloor;
            } else if(profs && strcmp(profs, "hopkins") == 0){
                /* Hopkins (2015, arXiv:1409.7395) Eq. 34, the mesh-code
                 * version of the cold Keplerian disk, scaled from his
                 * (0.5, 2) to (Rin, Rout):
                 *   rho = floor + (R/Rin)^3           R < Rin
                 *               + 1                   Rin <= R <= Rout
                 *               + [1+(R-Rout)/w]^-3   R > Rout
                 * floor 0.01 and w = 0.1*Rout/2 (0.4 for Rout = 8; the
                 * same as 0.1*Rin/0.5). EUNHA_KEPLER_FLOOR and
                 * EUNHA_KEPLER_WIDTH override them. P and v as below. */
                const char *fls = getenv("EUNHA_KEPLER_FLOOR");
                const char *ws = getenv("EUNHA_KEPLER_WIDTH");
                postype hfloor = (fls && fls[0]) ? atof(fls) : (postype)0.01;
                postype hw = (ws && ws[0]) ? atof(ws) : (postype)0.05*Rout;
                postype shape;
                if(R < Rin) { postype q = R/Rin; shape = q*q*q; }
                else if(R <= Rout) shape = 1;
                else { postype q = 1 + (R-Rout)/hw; shape = 1/(q*q*q); }
                rho = hfloor + shape;
                /* EUNHA_KEPLER_HALFBOX = L: keep only -L < x, y < L about the
                 * centre (the paper's -2 < x, y < 2 particle domain) */
                {
                    const char *hbs = getenv("EUNHA_KEPLER_HALFBOX");
                    if(hbs && hbs[0]){ postype hb = atof(hbs); if(fabs(dx_c) >= hb || fabs(dy_c) >= hb) keep = 0; }
                }
            } else {
            postype edge = (postype)0.25;
            postype inner = (postype)0.5*(1+tanh((R-Rin)/edge));
            postype outer = (postype)0.5*(1+tanh((Rout-R)/edge));
            postype floor = (postype)1.0e-3;
            rho = floor + (Rin / (R + (postype)1e-6)) * inner * outer;
            }
            postype Rsafe = (R > (postype)1e-8) ? R : (postype)1e-8;
            postype vamp = (eps > 0)
                ? R * pow(R*R + eps*eps, (postype)-0.75)
                : (postype)1.0 / sqrt(Rsafe);
            {
                /* EUNHA_KEPLER_SOFT=spline: v^2 = M(<R)/R with the cubic-spline enclosed mass (same as exam.c) */
                const char *ks = getenv("EUNHA_KEPLER_SOFT");
                if(ks && strcmp(ks, "spline") == 0){
                    const char *rcs = getenv("EUNHA_KEPLER_RC"); double rc = (rcs && rcs[0]) ? atof(rcs) : 0.25, u = (double)Rsafe/rc, M;
                    if(u >= 1) M = 1.0; else if(u <= 0.5) M = u*u*u*(32.0/3 - 192.0/5*u*u + 32*u*u*u); else M = -1.0/15 + u*u*u*(64.0/3 - 48*u + 192.0/5*u*u - 32.0/3*u*u*u);
                    vamp = (postype)sqrt(M/(double)Rsafe);
                }
                if(profs && strcmp(profs, "pakmor") == 0 && (R < Rin || R > Rout)){
                    /* background at rest (Pakmor) unless EUNHA_KEPLER_BG_ROTATE=1: co-rotating background (continuous v_c, no edge slip; global-dt variant) */
                    const char *bgr = getenv("EUNHA_KEPLER_BG_ROTATE");
                    if(!(bgr && atoi(bgr) > 0)) vamp = 0;
                }
            }
            {
                /* EUNHA_KEPLER_FORCE=gizmo: circular speed of that force,
                 * v^2 = g(r)/r with a_r = -g/r^2; at rest where g <= 0
                 * (the repulsive inner boundary, r < 0.216) */
                const char *kf = getenv("EUNHA_KEPLER_FORCE");
                if(kf && strcmp(kf, "gizmo") == 0){
                    double rr = (double)Rsafe, g;
                    if(rr <= 0.35) g = pow(rr/0.35, 2) - (0.35 - rr)/0.35;
                    else if(rr < 2.1) g = 1;
                    else g = 1 + (rr - 2.1)/0.1;
                    vamp = (postype)((g > 0) ? sqrt(g/rr) : 0.0);
                }
            }
            press = (postype)1.0e-6;
            if(profs && strcmp(profs, "pakmor") == 0){
                const char *ps = getenv("EUNHA_KEPLER_P");
                press = (ps && ps[0]) ? (postype)atof(ps) : (postype)((Gamma - 1)*(5.0*Gamma/2.0)*1.0e-5);
            }
            if(R < 1e-8){ vx = 0; vy_seed = 0; }
            else { vx = -dy_c/R*vamp; vy_seed = dx_c/R*vamp; }
        } else {
        postype profile = 0.5*(tanh((y-z1)/deltay) - tanh((y-z2)/deltay));
        rho = rho1 + (rho2-rho1)*profile;
        vx  = U1   + (U2-U1)*profile;
        press = KH_Pressure(simpar);
        vy_seed = dvy0*sin(2.0*M_PI*x)
                 *(exp(-(y-z1)*(y-z1)/(sigma*sigma))
                  +exp(-(y-z2)*(y-z2)/(sigma*sigma)));
        }
        }
        char iregion = (y > z1 && y < z2) ? 1 : 0;
#endif
        {
            postype vy = vy_seed;
            size_t indx = p;
            if(keep && x>=xmin && x < xmax && y>=ymin && y < ymax){
                res[np].u4if.indx  = indx;
                res[np].u4if.Flag[ENDIAN_OFFSET]  = iregion;
                res[np].x  = x;
                res[np].y  = y;
                res[np].vx  = vx;
                res[np].vy  = vy;
                /* EUNHA_IC_CELLVOL=<file>: per-particle initial cell volume (np_file doubles, glass order),
                 * e.g. the disc-clipped Voronoi area, so that m_i = rho V_i and U_i = P V_i/(g-1) match the
                 * target state cell by cell (Astra audit 16). Default: the mean volume. */
                postype vol_ic = meanvol;
                {
                    static double *cv_ = NULL; static long ncv_ = -1;
                    if(ncv_ < 0){
                        const char *cf_ = getenv("EUNHA_IC_CELLVOL"); ncv_ = 0;
                        if(cf_ && cf_[0]){
                            FILE *fc_ = fopen(cf_, "rb");
                            if(!fc_){ fprintf(stderr, "EUNHA_IC_CELLVOL: cannot open %s\n", cf_); exit(1); }
                            fseek(fc_, 0, SEEK_END); ncv_ = ftell(fc_)/sizeof(double); fseek(fc_, 0, SEEK_SET);
                            cv_ = (double*)malloc(ncv_*sizeof(double));
                            if(fread(cv_, sizeof(double), ncv_, fc_) != (size_t)ncv_){ fprintf(stderr, "EUNHA_IC_CELLVOL: short read\n"); exit(1); }
                            fclose(fc_);
                        }
                    }
                    if(ncv_ > 0){
                        if((long)p >= ncv_ || !(cv_[p] > 0)){ fprintf(stderr, "EUNHA_IC_CELLVOL: no volume for particle %ld\n", (long)p); exit(1); }
                        vol_ic = (postype)cv_[p];
                    }
                }
                res[np].mass  = rho*vol_ic;
                res[np].den  = rho;
                res[np].pressure  = press;
                res[np].ie  = res[np].pressure *vol_ic/(Gamma-1);
				/*
                res[np].ke  = 0.5*(res[np].vx*res[np].vx + res[np].vy*res[np].vy);
				*/
                res[np].csound  = sqrt(Gamma*res[np].pressure/res[np].den);
                res[np].w2  = (Lx/nx*GAS_Kappa(simpar))*(Lx/nx*GAS_Kappa(simpar));
                res[np].w2old  = res[np].w2;
                res[np].die  = 0;
				if(GAS_Kappa(simpar) <0) {
       		     res[np].w2 = - GAS_Kappa(simpar);
           		 res[np].w2old = - GAS_Kappa(simpar);
		        } 
				else {
			    	res[np].avgNeighboringPressure = press;
					res[np].w2 = getw2forHydroParticle(simpar, res+np,0);
					res[np].w2old = res[np].w2;
		            res[np].w2ceil = res[np].w2;
       		 	}
                np++;
            }
        }
    }
	if(use_glass){
		free(glass_x); free(glass_y);
	}
	if(use_stress){
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

	GAS_invw2Scale(simpar) = 1.L / KH_Pressure(simpar);
	/* w2Power is now read from params.dat ("GAS w2 power") */

	/* Specific-entropy normalization for w2_mode >= 1.
	 * s = ln(P/rho^gamma); KH has P common, two densities rho1, rho2.
	 * s_ref = (s1 + s2)/2, s_scale = |s1 - s2| = gamma*|ln(rho1/rho2)|.
	 * For mode 1 (entropy-A power-law), we also rescale invw2Scale to
	 * A_ref = P / rho_ref^gamma where rho_ref = sqrt(rho1*rho2) (geometric mean),
	 * so that the geometric mean of the two fluid w2 values ≈ base. */
	{
		postype P0 = KH_Pressure(simpar);
		postype rho1 = KH_Rho1(simpar), rho2 = KH_Rho2(simpar);
		postype gamma = GAS_GAMMA(simpar);
		postype s1 = log(P0 / pow(rho1, gamma));
		postype s2 = log(P0 / pow(rho2, gamma));
		GAS_W2SREF(simpar)   = 0.5*(s1 + s2);
		GAS_W2SSCALE(simpar) = fabs(s1 - s2);
		if(GAS_W2SSCALE(simpar) < 1e-12) GAS_W2SSCALE(simpar) = 1.0;
		if(GAS_W2MODE(simpar) == 1){
			postype rho_ref = sqrt(rho1 * rho2);
			GAS_invw2Scale(simpar) = pow(rho_ref, gamma) / P0;  /* 1/A_ref */
		}
	}


	int nbp = np;
	*mp = nbp;
	// voro-rk4 type tree-structure base particle.
	VORORK4_TBP(simpar) = res;
	VORORK4_BP(simpar) = (vorork4particletype*)res;
	// voro-rk4 type base particle.
	VORO_NP(simpar) = nbp;
	// declare the particles as not being boundary ghost particle
	{
		size_t p_size = (use_stress) ? sizeof(treevorostressrk4particletype)
		                             : sizeof(treevorork4particletype);
		char *bp_raw = (char*)res;
		for(i=0;i<np;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			bpi->u4if.Flag[ENDIAN_OFFSET] &= (~BoundaryGhostflag);
		}
	}
	DEBUGPRINT("P%d has np = %ld in box %g %g : %g %g\n", MYID(simpar), np, xmin,ymin,xmax,ymax);

    migrateTreeVorork4Particles(simpar);

	// measure the cell volume for each particle
	void kh_findVol(SimParameters *);
	kh_findVol(simpar);
	DEBUGPRINT("P%d passed initfindmass\n", MYID(simpar));

	res = VORORK4_TBP(simpar);
	nbp = VORO_NP(simpar);

	{
		size_t p_size = (use_stress) ? sizeof(treevorostressrk4particletype)
		                             : sizeof(treevorork4particletype);
		char *bp_raw = (char*)res;
		for(i=0;i<nbp;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			bpi->ie = bpi->pressure*bpi->volume/(Gamma-1);
			bpi->csound = sqrt(Gamma*bpi->pressure/bpi->den);
		}
	}
	VORO_NPAD(simpar) = 0;

	DEBUGPRINT("P%d has w2 values P= %g w2= %g w2ceil = %g den= %g meanV= %g vol= %g Cs= %g\n", 
			MYID(simpar), res[0].pressure, res[0].w2, res[0].w2ceil, res[0].den,
			meanvol, res[0].volume, res[0].csound);
    // update w2ceil and w2
    if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingTreeVorork4Particles);
	return res;
}


void kh_MkLinkedList(SimParameters *simpar){
	int np = VORO_NP(simpar);
	int i,j,k;
	postype xmin,ymin,xmax,ymax;
	int mx,my;
	postype pwidth = KH_GridSize(simpar);
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;

	BASICCELL_CELLWIDTH(simpar) = pwidth;
	BASICCELL_INVCELLWIDTH(simpar) = 1./pwidth;


	xmin = KH_XMIN(simpar)-pwidth;
	ymin = KH_YMIN(simpar)-pwidth;
	xmax = KH_XMAX(simpar)+pwidth;
	ymax = KH_YMAX(simpar)+pwidth;

	mx = ceil((xmax-xmin)/pwidth);
	my = ceil((ymax-ymin)/pwidth);

	BASICCELL_MX(simpar) = mx;
	BASICCELL_MY(simpar) = my;

	cells = VORO_BASICCELL(simpar)= (CellType*)realloc(VORO_BASICCELL(simpar),sizeof(CellType)*mx*my);


	postype cellsize = KH_GridSize(simpar);

	for(i=0;i<mx*my;i++) {
		cells[i].link = NULL;
		cells[i].nmem = 0;
	}

	DEBUGPRINT("P%d has np= %ld: %d %d: %g : %g %g\n", MYID(simpar), VORO_NP(simpar),
			mx,my, cellsize, xmin,ymin);

	char *bp_raw = (char*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		int ix,iy;
		ix = ((bpi->x-xmin)/cellsize);
		iy = ((bpi->y-ymin)/cellsize);
		size_t index = ix+mx*iy;
		if(index >=mx*my){
			printf("error detected: %ld : %d %d : %g %g < %g %g : %g %g >\n",
					index, mx,my, bpi->x, bpi->y, xmin,ymin, xmax,ymax);
			exit(0);
		}
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)bpi;
		cells[index].nmem ++;
		bpi->next = tmp;
	}
	DEBUGPRINT("P%d has %g %g : %g %g\n", MYID(simpar), KH_XMIN(simpar), KH_YMIN(simpar),
			KH_XMAX(simpar), KH_YMAX(simpar));
	paddingTreeVorork4Particles(simpar, pwidth);
	char *bpp_raw = (char*)VORORK4_TBPP(simpar);
	postype Xmin,Ymin,Xmax,Ymax;
	Xmin = Ymin = 1.e10;
	Xmax = Ymax = -1.e10;
	for(i=0;i<VORO_NPAD(simpar);i++){
		treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
		int ix,iy;
		ix = ((bppi->x-xmin)/cellsize);
		iy = ((bppi->y-ymin)/cellsize);
		Xmin = MIN(Xmin,bppi->x);
		Ymin = MIN(Ymin,bppi->y);
		Xmax = MAX(Xmax,bppi->x);
		Ymax = MAX(Ymax,bppi->y);
		size_t index = ix+mx*iy;
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)bppi;
		cells[index].nmem ++;
		bppi->next = tmp;
	}
#ifdef DEBUG
	DEBUGPRINT("P%d now after padding: %g %g : %g %g\n", MYID(simpar),Xmin,Ymin,Xmax,Ymax);
#endif
}

Voro2D_point *kh_Voro2D_FindNeighbor(SimParameters *simpar,
		int ix, int iy, int *nneigh){
	int i,j,k;
	int np,mp;
	Voro2D_point *neigh;
	int mx,my;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	cells = VORO_BASICCELL(simpar);
	/*
	float Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	float Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	*/

	np = 0;
	for(j=iy-1;j<=iy+1;j++){
		if(j <0 || j >=my) continue;
		for(i=ix-1;i<=ix+1;i++){
			if(i<0 || i>=mx) continue;
			np += cells[i+mx*j].nmem;
		}
	}
	neigh = (Voro2D_point*)malloc(sizeof(Voro2D_point)*np);
	np = 0;
	for(j=iy-1;j<=iy+1;j++){
		if(j<0 || j >=my)continue;
		for(i=ix-1;i<=ix+1;i++){
			if(i<0 || i>=mx) continue;
			struct linkedlisttype *tmp = cells[i+mx*j].link;
			while(tmp){
				treevorork4particletype *tt = (treevorork4particletype*)tmp;
				neigh[np].x = tt->x;
				neigh[np].y = tt->y;
				neigh[np].indx = PINDX(tt);
				neigh[np].csound = tt->csound;
				neigh[np].pressure = tt->pressure;
				neigh[np].w2 = tt->w2;
				neigh[np].bp = (void*)tt;
				np++;
				tmp = tmp->next;
			}

		}
	}
	*nneigh = np;
	return neigh;
}

treevorork4particletype *kh_Voro2D_FindCellBP(SimParameters *simpar,
		int ix, int iy, int *mp){
	int i,j,k;
	int np;
	treevorork4particletype *res;
	int mx,my;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	cells = VORO_BASICCELL(simpar);


	np = cells[ix+mx*iy].nmem;
	res = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*np);
	np = 0;
	size_t ipixel = ix+mx*iy;
	struct linkedlisttype *tmp = cells[ipixel].link;
	while(tmp){
		// exclude the boundary ghost particles.
		if(!IS_FLAG_ON(tmp,BoundaryGhostflag)){
			treevorork4particletype *tt = (treevorork4particletype*)tmp;
//			if(tt->u4if.Flag[ENDIAN_OFFSET]> 0)
			{
				res[np] = *tt;
				res[np].bp = tt;
				/*
				res[np].u4if.indx = tt-bp;
				res[np].x += iflag *Lx;
				res[np].y += jflag *Ly;
				*/
				np++;
			}
		}
		tmp = tmp->next;
	}
	*mp = np;
	return res;
}
Voro2D_point *searchCellRk4Neighbors2D_R(SimParameters *, int, int, int, int *);
int Voro2D_HasUnassignedCorner(Voro2D_Corner *, int);
int voro_expand_rmax(SimParameters *);
double voro_disc_radius(SimParameters *);
void voro2d_restrict_disc(Voro2D_Corner *, int, double);
void kh_findVol(SimParameters *simpar){
    treevorork4particletype *bp = VORORK4_TBP(simpar);
    int nbp = VORO_NP(simpar);
    /* The same cell as the hydro passes (Exam/exam.c): cap VORO_CELLBOX (the
     * old BOXSIZE/NX*5 used the unread cosmological Boxsize default 512),
     * widened neighbour search, disc restriction. The initial
     * ie = P V/(gamma-1) must use the volume the dynamics will see. */
    postype boxsize = VORO_CELLBOX(simpar);
    if(MYID(simpar) == 0)
        fprintf(stderr, "[KH_FINDVOL] cellbox=%g disc_R=%g\n", (double)boxsize, voro_disc_radius(simpar));
//  determine the minimum dpq for each particle by updating the w2ceil & w2
    det2d_dpqRK4(simpar,paddingTreeVorork4Particles);

//  DEBUGPRINT("P%d passed the initial minimum dpq %g %g %g\n", MYID(simpar), bp[0].w2ceil, bp[nbp/2].w2ceil, bp[nbp-1].w2ceil);


    postype xmin,ymin,zmin,xmax,ymax,zmax, cellsize;
    cellsize = BASICCELL_CELLWIDTH(simpar) = KH_GridSize(simpar);
    // xmin,ymin,xmax,ymax are the boundaries of the local domain.
    xmin = KH_XMIN(simpar)-cellsize;
    ymin = KH_YMIN(simpar)-cellsize;
    xmax = KH_XMAX(simpar)+cellsize;
    ymax = KH_YMAX(simpar)+cellsize;

//  DEBUGPRINT("P%d has cell info. %g : %g %g %g %g\n", MYID(simpar), cellsize, xmin,ymin,xmax,ymax);

    int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);
    // prepare the linked-list cells for mkLinkedList2D()  and tree findings below
    CellType *cells = VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*mx*my);
    // building the linked list with "paddingTreeVorork4Particles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
#if defined(IC_SOD2D)
    mkLinkedList2D_sod(simpar, cellsize,xmin,ymin,xmax,ymax, paddingTreeVorork4Particles);
#else
    mkLinkedList2D_oExam(simpar, cellsize,xmin,ymin,xmax,ymax, paddingTreeVorork4Particles);
#endif

    int i;
    /*
    for(i=0;i<nbp;i++){
        bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
    }
    */
    int iy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(iy=0;iy<my;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=0;ix<mx;ix++){
            int np;
            treevorork4particletype *p = findCellRk4BP2D(simpar,ix,iy,&np);
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
                ibp->dt = 1.e10; // initialization of dt
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
                {
                    Voro2D_point *ne_k = NULL, *nwe_k = NULL;
                    if(Voro2D_HasUnassignedCorner(vorocorner, ip)){
                        int R, nne;
                        for(R=2; R<=voro_expand_rmax(simpar); R++){
                            if(ne_k){ free(ne_k); free(nwe_k); }
                            ne_k = searchCellRk4Neighbors2D_R(simpar, ix, iy, R, &nne);
                            nwe_k = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nne);
                            ip = Voro2D_FindVC(&center, ne_k, nwe_k, nne, vorocorner, mp, boxsize);
                            if(!Voro2D_HasUnassignedCorner(vorocorner, ip)) break;
                        }
                    }
                    double vrR = voro_disc_radius(simpar);
                    if(vrR > 0) voro2d_restrict_disc(vorocorner, mp, vrR);
                    if(ne_k){ free(ne_k); free(nwe_k); }
                }
                ibp->volume = Area2DPolygon(vorocorner,mp);
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

}
