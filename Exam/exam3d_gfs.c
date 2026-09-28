/*
 * exam3d_gfs.c -- 3D GFS (geometric face scheme) hydro path, CPU first.
 *
 * Port of the 2D production path (exam.c: exam2d_vph_rk4_int_blend +
 * getAccVoro2DBlend_impl av_mode 5 + updateDenW2Pressure2DBlend) with the
 * dimensional changes only:
 *   edge length              -> face polygon area |S_f|
 *   edge midpoint            -> polygon area centroid c_f
 *   sqrt(V) (cell length)    -> cbrt(V)  (pair law, ghost CFL floor,
 *                                         acceleration CFL, entropy switch)
 *   2x2 gradients            -> 3x3 Green-Gauss gradients (no M(n,m) term;
 *                                the 2D four-point stencil has no 3D analogue)
 *   gfs_pair_work_limit nshare 6 -> GFS3D_PAIR_NSHARE (default 16, see below)
 * Geometry: Voro3D_FindVC (Voro/voro.c) per particle, w2 set explicitly on
 * the centre and every neighbour (0 = Voronoi). Exam/Sedov (the old 3D
 * prototype) is used only as the face-loop template; none of its physics is
 * reused (no ie >= 0 clamp, no pinned particle, dt from the full summed
 * accelerations, RK4 instead of its two-pass first-order update).
 *
 * Per face (i, j), av_mode 5:
 *   w = v_i + urad, urad = get3dUpqrad (0.5 (v_j - v_i) for equal w2)
 *       + Springel face rotation  -((v_j - v_i).(c_f - a)/d) e_r,
 *         a = midpoint (w = 0) or Laguerre anchor fact1 (x_j - x_i)
 *         (GFS_LAGUERRE_ROTATION=1)
 *   MUSCL {rho, P, v_n} at c_f with Barth-Jespersen limited gradients,
 *   HLLC in the face rest frame (gfs_riemann.h = exam.c solvers),
 *   extreme face (SEDOV_PHASE1, cell P ratio > 100): HLL star state along
 *   e_r from the cell-centred states, energy flux with v* e_r,
 *   optional Monaghan AV (GAS AlphaVis), capped pair pressure
 *   (gfs_pair_pressure_len with cbrt(V)) limited by gfs_pair_work_limit.
 *   dE/dt += -p* (u_a . S), dv/dt += -p* S / m.
 * Integrator: RK4 exactly as exam2d_vph_rk4_int_blend (dt from K1, stage
 * shifts x0+k1/2, x0+k2/2, x0+k3, final (k1+2k2+2k3+k4)/6), k_v includes
 * self-gravity and the external-force hook, k_ie = (dE/dt - m v.a_hydro) dt
 * under SEDOV_PHASE1. Stage floor P<=0 -> 1e-6 (or K rho^gamma with
 * GFS_DUAL_ENERGY) is booked (sfl_cum) in the ledger, final floor likewise
 * (floor_cum). GFS_DUAL_ENERGY, GFS_ENTROPY_SWITCH / GFS_HALF_LIMIT /
 * GFS_ES_COEF, GFS_FLOOR_LOG have the 2D meaning.
 *
 * Pair work-limit face count (nshare). In 2D nshare = 6 = the mean face
 * count of any planar Voronoi tessellation (Euler), so one face may spend at
 * most ie/6 of each cell. The 3D analogue is the mean face count of a 3D
 * Voronoi cell: 15.54 (Poisson), ~14.5 (relaxed glass), 14 (bcc), 12 (fcc),
 * 6 non-degenerate faces on a cubic lattice. 16 >= all of these, so the
 * per-face budget never exceeds ie/(mean face count), i.e. the sum over a
 * typical cell's faces stays below ie. Override with GFS3D_PAIR_NSHARE.
 *
 * Boundaries per axis (params "Hydro3D boundary", default from the IC
 * periodic flags: periodic, else reflect):
 *   periodic : periodic images through the link-cell stencil.
 *   reflect  : mirror images (v_n reversed). Cells are cut exactly by the
 *              wall plane; the wall face carries the reflective HLLC
 *              pressure and zero energy flux. Particles that cross a wall
 *              are reflected at the end of the step (not at RK stages).
 *   outflow  : mirror images with copied velocity (zero gradient);
 *              particles outside the box at the end of a step are removed
 *              and their energy is booked (out_cum).
 * Self-gravity ("Hydro3D gravity = 1"): Barnes-Hut octree, monopole,
 * opening angle LAG3D_THETA (default 0.5), Plummer softening from params
 * (or the IC header). LAG3D_GRAV_DIRECT=1 uses O(N^2) direct summation,
 * =0 forces the tree; unset: direct for N <= 20000, tree above.
 * Isolated only (error with periodic axes). Epot = 1/2 sum m phi.
 * External hook: LAG3D_PM_GM / LAG3D_PM_X,Y,Z / LAG3D_PM_EPS point mass,
 * LAG3D_ACC="ax,ay,az" uniform field; both enter k_v and Epot.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <ctype.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "../indxflag.h"   /* voro.h needs indxflag (voro.c defines VORO_MAIN) */
#include "voro.h"
#include "gfs_pair.h"
#include "gfs_riemann.h"
#include "exam3d_gfs.h"
#include "Tests3D/common/lag3d_ic.h"

__attribute__((used)) static const char lag3d_marker[] = "LAGEUNHA_3D_GFS_V1";

#define G3_MP 8192          /* Voro3D vertex buffer per thread */
#define BC_PERIODIC 0
#define BC_REFLECT  1
#define BC_OUTFLOW  2

typedef struct {
	int j;                  /* real particle index of the neighbour */
	int mf;                 /* bit a set: image mirrored across axis a */
	double rel[3];          /* x_j(image) - x_i */
	double S[3];            /* outward area vector */
	double c[3];            /* face centroid relative to x_i */
} g3face;

typedef struct { g3face *f; long n, cap; } g3fbuf;

typedef struct {
	/* parameters */
	double gamma, courant, tend, dumpdt, alphavis, betavis, etavis;
	double kappa, eps, G, theta, nshare;
	int av_mode, use_muscl, gravity, grav_direct, maxsteps;
	int bc[3];
	double lo[3], hi[3], L[3];
	char icfile[512];
	/* switches (environment, 2D names) */
	int phase1, floor_log, lag_rot, es_on, half_on;
	double de_eta, es_coef;
	/* external hook */
	double pm_gm, pm_x[3], pm_eps, acc[3];
	int pm_on, acc_on;
	/* particles */
	long n;
	double *x, *y, *z, *vx, *vy, *vz, *m, *ie, *vol, *den, *P, *cs, *K, *w2, *w2old;
	double *ax, *ay, *az, *gx, *gy, *gz, *phi, *phiext, *die, *dK, *ie0, *V0, *dtp;
	double *grad;           /* 15 per particle: drho[3] dP[3] dv[3][3] (dv[a][b] = d v_a / d x_b) */
	double *k;              /* 4 stages x 8 (x y z vx vy vz ie K) per particle */
	int64_t *id;
	/* faces */
	int nthreads;
	g3fbuf *tb;
	int *ftid, *fcnt;
	long *foff;
	int geom_valid, grav_valid;
	/* link cells */
	int nc[3];
	double cw[3], cwmin;
	int *head, *next;
	/* bookkeeping */
	double dtold, time;
	long step;
	double sfl_step, sfl_cum, fl_cum, de_cum, es_cum, hl_cum, out_cum;
	long long npair;
	int max_R;
	long n_expand, n_jitter;
	int sfl_nlog, stage_reset;
	double E0tot, E0hyd;
	int E0tot_set;
} g3sim;

static g3sim *G3S = NULL;   /* for the qsort-free helpers below */

/* ------------------------------------------------------------ params */
static void g3_trim(char *s){
	char *e;
	while(*s && isspace((unsigned char)*s)) memmove(s, s+1, strlen(s));
	e = s + strlen(s);
	while(e > s && isspace((unsigned char)e[-1])) *--e = 0;
}

/* value of "define <key> = value # comment", or NULL */
static int g3_param(const char *file, const char *key, char *out, int len){
	FILE *fp = fopen(file, "r");
	char line[1024];
	int found = 0;
	if(!fp) return 0;
	while(fgets(line, sizeof(line), fp)){
		char *p = line, *eq, *h;
		if(strncmp(line, "#End of Ascii Header", 20) == 0) break;
		while(*p && isspace((unsigned char)*p)) p++;
		if(strncmp(p, "define", 6) != 0) continue;
		p += 6;
		eq = strchr(p, '=');
		if(!eq) continue;
		*eq = 0;
		{
			char k[512];
			strncpy(k, p, sizeof(k)-1); k[sizeof(k)-1] = 0;
			g3_trim(k);
			if(strcmp(k, key) != 0) continue;
		}
		h = strchr(eq+1, '#');
		if(h) *h = 0;
		strncpy(out, eq+1, len-1); out[len-1] = 0;
		g3_trim(out);
		found = 1;
	}
	fclose(fp);
	return found;
}

static double g3_pd(const char *file, const char *key, double def){
	char v[512];
	if(g3_param(file, key, v, sizeof(v)) && v[0]) return atof(v);
	return def;
}

int exam3d_gfs_is_hydro3d(const char *paramfile){
	char v[512];
	if(!paramfile) return 0;
	if(!g3_param(paramfile, "Simulation Model", v, sizeof(v))) return 0;
	return strcmp(v, "Hydro3D") == 0;
}

static int g3_env_on(const char *name){
	const char *s = getenv(name);
	return (s && s[0] == '1') ? 1 : 0;
}

static double g3_env_d(const char *name, double def){
	const char *s = getenv(name);
	return (s && s[0]) ? atof(s) : def;
}

/* ------------------------------------------------------------ memory */
static void *g3_alloc(size_t n){
	void *p = calloc(n ? n : 1, 1);
	if(!p){ fprintf(stderr, "[3D] out of memory (%zu bytes)\n", n); exit(12); }
	return p;
}

static void g3_alloc_particles(g3sim *s, long n){
	double **f[] = {&s->x,&s->y,&s->z,&s->vx,&s->vy,&s->vz,&s->m,&s->ie,&s->vol,&s->den,
		&s->P,&s->cs,&s->K,&s->w2,&s->w2old,&s->ax,&s->ay,&s->az,&s->gx,&s->gy,&s->gz,
		&s->phi,&s->phiext,&s->die,&s->dK,&s->ie0,&s->V0,&s->dtp};
	size_t nf = sizeof(f)/sizeof(f[0]), q;
	for(q=0;q<nf;q++) *f[q] = (double*)g3_alloc(sizeof(double)*n);
	s->grad = (double*)g3_alloc(sizeof(double)*15*n);
	s->k = (double*)g3_alloc(sizeof(double)*32*n);
	s->id = (int64_t*)g3_alloc(sizeof(int64_t)*n);
	s->ftid = (int*)g3_alloc(sizeof(int)*n);
	s->fcnt = (int*)g3_alloc(sizeof(int)*n);
	s->foff = (long*)g3_alloc(sizeof(long)*n);
	s->next = (int*)g3_alloc(sizeof(int)*n);
}

/* ------------------------------------------------------------ link cells */
static void g3_build_cells(g3sim *s){
	long i;
	int a, nct;
	double hmean = cbrt(s->L[0]*s->L[1]*s->L[2] / (double)(s->n > 0 ? s->n : 1));
	for(a=0;a<3;a++){
		int n = (int)(s->L[a] / (2.0*hmean));
		if(n < 1) n = 1;
		if(n > 1024) n = 1024;
		s->nc[a] = n;
		s->cw[a] = s->L[a] / n;
	}
	s->cwmin = fmin(s->cw[0], fmin(s->cw[1], s->cw[2]));
	nct = s->nc[0]*s->nc[1]*s->nc[2];
	free(s->head);
	s->head = (int*)g3_alloc(sizeof(int)*nct);
	for(a=0;a<nct;a++) s->head[a] = -1;
	for(i=0;i<s->n;i++){
		double p[3] = {s->x[i], s->y[i], s->z[i]};
		int c[3];
		for(a=0;a<3;a++){
			c[a] = (int)floor((p[a] - s->lo[a]) / s->cw[a]);
			if(c[a] < 0) c[a] = 0;
			if(c[a] >= s->nc[a]) c[a] = s->nc[a]-1;
		}
		int ic = (c[2]*s->nc[1] + c[1])*s->nc[0] + c[0];
		s->next[i] = s->head[ic];
		s->head[ic] = (int)i;
	}
}

/* ------------------------------------------------------------ tessellation */
typedef struct {
	Voro3D_point *cand, *nw;
	int *cj, *cmf;
	double *crel;
	long cap;
	Voro3D_Vertex *vert;
	int *iwork;
	g3face *tmpf;
	long tmpcap;
} g3work;

static void g3_work_init(g3work *w){
	memset(w, 0, sizeof(*w));
	w->vert = (Voro3D_Vertex*)g3_alloc(sizeof(Voro3D_Vertex)*G3_MP);
	w->iwork = (int*)g3_alloc(sizeof(int)*G3_MP);
}

static void g3_work_free(g3work *w){
	free(w->cand); free(w->nw); free(w->cj); free(w->cmf); free(w->crel);
	free(w->vert); free(w->iwork); free(w->tmpf);
	memset(w, 0, sizeof(*w));
}

static void g3_work_grow(g3work *w, long need){
	if(need <= w->cap) return;
	long nc = need*2 + 256;
	w->cand = (Voro3D_point*)realloc(w->cand, sizeof(Voro3D_point)*nc);
	/* 8 guard entries in front: get3Dw2pCeil (end of Voro3D_FindVC) reads
	 * neighbors[related] also for the initial-cube walls (related -1..-6)
	 * before we get to reject such a cell. */
	free(w->nw);
	w->nw = (Voro3D_point*)calloc(nc + 8, sizeof(Voro3D_point));
	w->cj = (int*)realloc(w->cj, sizeof(int)*nc);
	w->cmf = (int*)realloc(w->cmf, sizeof(int)*nc);
	w->crel = (double*)realloc(w->crel, sizeof(double)*3*nc);
	if(!w->cand || !w->nw || !w->cj || !w->cmf || !w->crel){
		fprintf(stderr, "[3D] out of memory (candidates)\n"); exit(12);
	}
	w->cap = nc;
}

/* Area vector and area centroid of the face that starts on edge
 * start -> next. Same traversal as Voro3D_norm_polygon (marks considered). */
static int g3_polygon(Voro3D_Vertex *start, Voro3D_Vertex *next, double S[3], double c[3]){
	Voro3D_Vertex *now = start, *nnext;
	double M[3][3] = {{0}}, vs[3] = {0};
	int nv = 0, it = 0, i, a, b;
	S[0] = S[1] = S[2] = 0;
	do {
		i = 0;
		while(next->link[i] != now){
			i = (i+1)%3;
			if(++it > 100000) return -1;
		}
		i = (i-1+3)%3;
		nnext = next->link[i];
		double r1[3] = {next->x - start->x, next->y - start->y, next->z - start->z};
		double r2[3] = {nnext->x - next->x, nnext->y - next->y, nnext->z - next->z};
		double A[3] = {0.5*(r1[1]*r2[2]-r1[2]*r2[1]), 0.5*(r1[2]*r2[0]-r1[0]*r2[2]),
			0.5*(r1[0]*r2[1]-r1[1]*r2[0])};
		double ct[3] = {(start->x+next->x+nnext->x)/3.0, (start->y+next->y+nnext->y)/3.0,
			(start->z+next->z+nnext->z)/3.0};
		for(a=0;a<3;a++){
			S[a] += A[a];
			for(b=0;b<3;b++) M[a][b] += A[a]*ct[b];
		}
		vs[0] += next->x; vs[1] += next->y; vs[2] += next->z; nv++;
		next->considered[i] = Yes;
		now = next;
		next = nnext;
		if(nv > 10000) return -1;
	} while(now != start);
	double s2 = S[0]*S[0] + S[1]*S[1] + S[2]*S[2];
	if(s2 > 0){
		for(b=0;b<3;b++) c[b] = (S[0]*M[0][b] + S[1]*M[1][b] + S[2]*M[2][b]) / s2;
	} else {
		for(b=0;b<3;b++) c[b] = vs[b] / (nv > 0 ? nv : 1);
	}
	return nv;
}

/* Voronoi/Laguerre cell of particle i. Faces go to the thread buffer.
 * Returns the volume (<= 0 on failure). */
static double g3_tess_one(g3sim *s, long i, g3work *w, g3fbuf *tb, int tid){
	double xi[3] = {s->x[i], s->y[i], s->z[i]};
	int ci[3], a, R;
	for(a=0;a<3;a++){
		ci[a] = (int)floor((xi[a] - s->lo[a]) / s->cw[a]);
		if(ci[a] < 0) ci[a] = 0;
		if(ci[a] >= s->nc[a]) ci[a] = s->nc[a]-1;
	}
	int Rlimit = s->nc[0];
	if(s->nc[1] > Rlimit) Rlimit = s->nc[1];
	if(s->nc[2] > Rlimit) Rlimit = s->nc[2];
	Rlimit += 2;
	int why = 0, attempt = 0;
	for(R = 1; R <= Rlimit; R++){
		long nc = 0;
		int d[3];
		for(d[2]=-R; d[2]<=R; d[2]++) for(d[1]=-R; d[1]<=R; d[1]++) for(d[0]=-R; d[0]<=R; d[0]++){
			int cc[3], mir[3], skip = 0;
			double sh[3];
			for(a=0;a<3;a++){
				int c = ci[a] + d[a];
				sh[a] = 0; mir[a] = 0;
				if(s->bc[a] == BC_PERIODIC){
					int wv = (c >= 0) ? c / s->nc[a] : -((-c + s->nc[a] - 1) / s->nc[a]);
					cc[a] = c - wv*s->nc[a];
					sh[a] = wv * s->L[a];
				} else {
					if(c < 0){ cc[a] = -1 - c; mir[a] = 1; }
					else if(c >= s->nc[a]){ cc[a] = 2*s->nc[a] - 1 - c; mir[a] = 2; }
					else cc[a] = c;
					if(cc[a] < 0 || cc[a] >= s->nc[a]) skip = 1;
				}
			}
			if(skip) continue;
			int ic = (cc[2]*s->nc[1] + cc[1])*s->nc[0] + cc[0];
			int ident = (sh[0] == 0 && sh[1] == 0 && sh[2] == 0 && !mir[0] && !mir[1] && !mir[2]);
			int mf = (mir[0] ? 1 : 0) | (mir[1] ? 2 : 0) | (mir[2] ? 4 : 0);
			int j;
			for(j = s->head[ic]; j >= 0; j = s->next[j]){
				double p[3] = {s->x[j], s->y[j], s->z[j]};
				if(ident && j == i) continue;
				for(a=0;a<3;a++){
					p[a] += sh[a];
					if(mir[a] == 1) p[a] = 2*s->lo[a] - p[a];
					else if(mir[a] == 2) p[a] = 2*s->hi[a] - p[a];
				}
				g3_work_grow(w, nc+1);
				Voro3D_point *q = w->cand + nc;
				memset(q, 0, sizeof(*q));
				q->x = p[0] - xi[0]; q->y = p[1] - xi[1]; q->z = p[2] - xi[2];
				q->w2 = s->w2[j];
				q->indx = (int)nc;
				w->cj[nc] = j; w->cmf[nc] = mf;
				w->crel[3*nc] = q->x; w->crel[3*nc+1] = q->y; w->crel[3*nc+2] = q->z;
				nc++;
			}
		}
		if(nc < 4) continue;
		/* Retry after a topology failure: Voro3D_FindVC can build a broken
		 * vertex graph when many generators are exactly cospherical (cubic
		 * lattice + symmetric flow). Jitter the neighbour positions by
		 * 1e-9..1e-7 cell widths (deterministic in j) for the geometry only;
		 * rel[] used by the physics stays exact. */
		if(attempt > 0){
			double amp = 1e-9 * pow(10.0, attempt-1) * s->cwmin;
			long q;
			for(q=0;q<nc;q++){
				uint64_t hsh = (uint64_t)(w->cj[q]*8 + w->cmf[q]) * 0x9E3779B97F4A7C15ULL
					+ (uint64_t)attempt * 0xBF58476D1CE4E5B9ULL;
				double jit[3];
				for(a=0;a<3;a++){
					hsh ^= hsh >> 31; hsh *= 0x94D049BB133111EBULL; hsh ^= hsh >> 29;
					jit[a] = ((double)(hsh >> 11) / 9007199254740992.0) * 2.0 - 1.0;
				}
				w->cand[q].x += amp*jit[0]; w->cand[q].y += amp*jit[1]; w->cand[q].z += amp*jit[2];
			}
		}
		Voro3D_point ctr;
		memset(&ctr, 0, sizeof(ctr));
		ctr.w2 = s->w2[i];
		ctr.indx = -1;
		double boxsize = 2.0 * R * s->cwmin;
		/* mp only sizes Voro3D_InitializeVertex's "mark the rest Inactive"
		 * loop; FindVC reads vertices [0, ip) only, so a small value avoids
		 * touching the whole G3_MP buffer for every cell. */
		int ip = Voro3D_FindVC(&ctr, w->cand, w->nw + 8, (int)nc, w->vert, 16, boxsize, 0, w->iwork);
		if(ip <= 0 || ip >= G3_MP){ why = 1; if(attempt < 3){ attempt++; R--; } else attempt = 0; continue; }
		int bad = 0, iv, jj;
		double rmax2 = 0;
		for(iv=0; iv<ip; iv++){
			Voro3D_Vertex *v = w->vert + iv;
			if(v->status != Active) continue;
			double r2 = v->x*v->x + v->y*v->y + v->z*v->z;
			if(r2 > rmax2) rmax2 = r2;
			for(jj=0;jj<3;jj++) if(v->related[jj] < 0) bad = 1;
			v->considered[0] = v->considered[1] = v->considered[2] = No;
		}
		if(bad){ why = 2; attempt = 0; continue; }
		if(2.0*sqrt(rmax2) > R * s->cwmin){ why = 3; attempt = 0; continue; }
		/* faces */
		long nf = 0;
		double V = 0;
		for(iv=0; iv<ip && !bad; iv++){
			Voro3D_Vertex *v = w->vert + iv;
			if(v->status != Active) continue;
			for(jj=0;jj<3;jj++){
				if(v->considered[jj] != No) continue;
				int kk = v->related[(jj+2)%3];
				double S[3], c[3];
				if(g3_polygon(v, v->link[jj], S, c) < 0 || kk < 0 || kk >= (int)nc){ bad = 1; break; }
				int cid = w->nw[8 + kk].indx;
				if(cid < 0 || cid >= (int)nc){ bad = 1; break; }
				if(nf >= w->tmpcap){
					w->tmpcap = 2*w->tmpcap + 64;
					w->tmpf = (g3face*)realloc(w->tmpf, sizeof(g3face)*w->tmpcap);
					if(!w->tmpf){ fprintf(stderr, "[3D] out of memory (faces)\n"); exit(12); }
				}
				g3face *f = w->tmpf + nf;
				f->j = w->cj[cid]; f->mf = w->cmf[cid];
				for(a=0;a<3;a++){ f->rel[a] = w->crel[3*cid+a]; f->S[a] = S[a]; f->c[a] = c[a]; }
				V += (c[0]*S[0] + c[1]*S[1] + c[2]*S[2]) / 3.0;
				nf++;
			}
		}
		if(bad){ why = 4; if(attempt < 3){ attempt++; R--; } else attempt = 0; continue; }
		double sgn = 1.0;
		if(V < 0){ sgn = -1.0; V = -V; }
		if(!(V > 0)){ why = 5; if(attempt < 3){ attempt++; R--; } else attempt = 0; continue; }
		double amin = 1e-12 * pow(V, 2.0/3.0);
		long f0 = tb->n, q;
		for(q=0;q<nf;q++){
			g3face *f = w->tmpf + q;
			for(a=0;a<3;a++) f->S[a] *= sgn;
			double A = sqrt(f->S[0]*f->S[0] + f->S[1]*f->S[1] + f->S[2]*f->S[2]);
			if(!(A > amin)) continue;
			if(tb->n >= tb->cap){
				tb->cap = 2*tb->cap + 1024;
				tb->f = (g3face*)realloc(tb->f, sizeof(g3face)*tb->cap);
				if(!tb->f){ fprintf(stderr, "[3D] out of memory (face buffer)\n"); exit(12); }
			}
			tb->f[tb->n++] = *f;
		}
		s->ftid[i] = tid;
		s->foff[i] = f0;
		s->fcnt[i] = (int)(tb->n - f0);
		if(R > 1){
#ifdef _OPENMP
#pragma omp atomic
#endif
			s->n_expand++;
		}
		if(attempt > 0){
#ifdef _OPENMP
#pragma omp atomic
#endif
			s->n_jitter++;
		}
		if(R > s->max_R){
#ifdef _OPENMP
#pragma omp critical (g3maxr)
#endif
			{ if(R > s->max_R) s->max_R = R; }
		}
		return V;
	}
	/* why: 1 vertex buffer, 2 cell touches the initial cube, 3 security
	 * radius not covered, 4 face topology, 5 volume <= 0 */
	fprintf(stderr, "[3D] tessellation failed for particle %ld at (%g %g %g) up to R=%d (reason %d)\n",
			i, xi[0], xi[1], xi[2], Rlimit, why);
	return -1;
}

static inline g3face *g3_faces(g3sim *s, long i){
	return s->tb[s->ftid[i]].f + s->foff[i];
}

/* Tessellate all particles: faces + volumes. Returns 0 on success. */
static int g3_tessellate(g3sim *s){
	int t, fail = 0;
	g3_build_cells(s);
	for(t=0;t<s->nthreads;t++) s->tb[t].n = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+:fail)
#endif
	{
		int tid = 0;
		long i;
#ifdef _OPENMP
		tid = omp_get_thread_num();
#endif
		g3work w;
		g3_work_init(&w);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
		for(i=0;i<s->n;i++){
			double V = g3_tess_one(s, i, &w, &s->tb[tid], tid);
			if(!(V > 0)){ fail++; s->vol[i] = 0; s->fcnt[i] = 0; }
			else s->vol[i] = V;
		}
		g3_work_free(&w);
	}
	return fail;
}

/* ------------------------------------------------------------ images */
/* velocity sign r[a] and position sign sg[a] of the image described by mf */
static inline void g3_image_signs(const g3sim *s, int mf, double r[3], double sg[3]){
	int a;
	for(a=0;a<3;a++){
		int m = (mf >> a) & 1;
		sg[a] = m ? -1.0 : 1.0;
		r[a] = (m && s->bc[a] == BC_REFLECT) ? -1.0 : 1.0;
	}
}

static void g3_center(const g3sim *s, double c[3]);
/* ------------------------------------------------------------ primitives */
/* den, P, c from ie and V. The stage floor (P <= 0) is the 2D one
 * (updateDenW2Pressure2DBlend): K rho^gamma with GFS_DUAL_ENERGY and K > 0,
 * else P = 1e-6. The injected ie is summed in sfl_step (booked). */
static void g3_primitives(g3sim *s){
	long i;
	double inj = 0, gm1 = s->gamma - 1.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:inj)
#endif
	for(i=0;i<s->n;i++){
		double V = s->vol[i];
		if(!(V > 0)) continue;
		s->den[i] = s->m[i] / V;
		s->P[i] = s->ie[i] / V * gm1;
		if(s->P[i] <= 0){
			double ie_old = s->ie[i];
			if(s->de_eta > 0 && s->K[i] > 0 && s->den[i] > 0)
				s->P[i] = s->K[i] * pow(s->den[i], s->gamma);
			else
				s->P[i] = 1e-6;
			if(s->stage_reset){
				/* 2D behaviour: ie is reset, the difference persists
				 * through the RK4 combination (booked as sfl). */
				s->ie[i] = s->P[i] * V / gm1;
				inj += s->ie[i] - ie_old;
			}
			/* else: pressure-only floor, ie keeps the energy-equation
			 * value; the end-of-step floor (floor_cum) is the only refill. */
			if(s->floor_log){
				int nl;
#ifdef _OPENMP
#pragma omp atomic capture
#endif
				nl = s->sfl_nlog++;
				if(nl < 8){
					double c0[3];
					g3_center(s, c0);
					double dx = s->x[i]-c0[0], dy = s->y[i]-c0[1], dz = s->z[i]-c0[2];
					fprintf(stderr, "[RK4SF] step=%ld i=%ld x=%.5f y=%.5f z=%.5f R=%.4f den=%.4e ie_stage=%.4e vol=%.4e |v|=%.4e\n",
							s->step+1, i, s->x[i], s->y[i], s->z[i], sqrt(dx*dx+dy*dy+dz*dz),
							s->den[i], ie_old, V,
							sqrt(s->vx[i]*s->vx[i] + s->vy[i]*s->vy[i] + s->vz[i]*s->vz[i]));
				}
			}
		}
		s->cs[i] = sqrt(s->gamma * s->P[i] / s->den[i]);
	}
	s->sfl_step += inj;
}

/* Green-Gauss gradients and the av_mode 5 Barth-Jespersen limiter
 * (updateDenW2Pressure2DBlend, 3D). */
static void g3_gradients(g3sim *s){
	long i;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
	for(i=0;i<s->n;i++){
		double *g = s->grad + 15*i;
		int q, a, b, nf = s->fcnt[i];
		g3face *F = g3_faces(s, i);
		double V = s->vol[i];
		double vi[3] = {s->vx[i], s->vy[i], s->vz[i]};
		for(q=0;q<15;q++) g[q] = 0;
		if(!(V > 0) || nf == 0) continue;
		for(q=0;q<nf;q++){
			g3face *f = F + q;
			int j = f->j;
			double r[3], sg[3];
			g3_image_signs(s, f->mf, r, sg);
			double d2 = f->rel[0]*f->rel[0] + f->rel[1]*f->rel[1] + f->rel[2]*f->rel[2];
			double wf = (d2 > 0) ? 0.5 + 0.5*(s->w2[i] - s->w2[j])/d2 : 0.5;
			double vj[3] = {r[0]*s->vx[j], r[1]*s->vy[j], r[2]*s->vz[j]};
			double dr = wf*(s->den[j] - s->den[i]);
			double dp = wf*(s->P[j] - s->P[i]);
			for(b=0;b<3;b++){
				g[b]   += dr * f->S[b];
				g[3+b] += dp * f->S[b];
				for(a=0;a<3;a++) g[6+3*a+b] += wf*(vj[a] - vi[a]) * f->S[b];
			}
		}
		for(q=0;q<15;q++) g[q] /= V;
		if(s->av_mode == 5){
			double Qi[5] = {s->den[i], vi[0], vi[1], vi[2], s->P[i]};
			double Qmax[5], Qmin[5], alpha[5] = {1,1,1,1,1};
			for(q=0;q<5;q++) Qmax[q] = Qmin[q] = Qi[q];
			for(q=0;q<nf;q++){
				g3face *f = F + q;
				int j = f->j;
				if(f->mf) continue;         /* ghosts excluded, as in 2D */
				double Qj[5] = {s->den[j], s->vx[j], s->vy[j], s->vz[j], s->P[j]};
				for(a=0;a<5;a++){
					if(Qj[a] > Qmax[a]) Qmax[a] = Qj[a];
					if(Qj[a] < Qmin[a]) Qmin[a] = Qj[a];
				}
			}
			for(q=0;q<nf;q++){
				g3face *f = F + q;
				double dl[5];
				dl[0] = g[0]*f->c[0] + g[1]*f->c[1] + g[2]*f->c[2];
				for(a=0;a<3;a++) dl[1+a] = g[6+3*a]*f->c[0] + g[6+3*a+1]*f->c[1] + g[6+3*a+2]*f->c[2];
				dl[4] = g[3]*f->c[0] + g[4]*f->c[1] + g[5]*f->c[2];
				for(a=0;a<5;a++){
					double rr;
					if(dl[a] > 1e-30){
						rr = (Qmax[a] - Qi[a]) / dl[a];
						if(rr < 0) rr = 0;
						if(rr < alpha[a]) alpha[a] = rr;
					} else if(dl[a] < -1e-30){
						rr = (Qmin[a] - Qi[a]) / dl[a];
						if(rr < 0) rr = 0;
						if(rr < alpha[a]) alpha[a] = rr;
					}
				}
			}
			for(b=0;b<3;b++){
				g[b] *= alpha[0];
				g[3+b] *= alpha[4];
				for(a=0;a<3;a++) g[6+3*a+b] *= alpha[1+a];
			}
		}
	}
}

/* ------------------------------------------------------------ gravity */
static void g3_ext_accel(const g3sim *s, double x, double y, double z,
		double *ax, double *ay, double *az, double *phi){
	*ax = *ay = *az = 0; *phi = 0;
	if(s->pm_on){
		double dx = x - s->pm_x[0], dy = y - s->pm_x[1], dz = z - s->pm_x[2];
		double r2 = dx*dx + dy*dy + dz*dz + s->pm_eps*s->pm_eps;
		double r = sqrt(r2);
		if(r > 0){
			double f = -s->pm_gm / (r2*r);
			*ax += f*dx; *ay += f*dy; *az += f*dz;
			*phi += -s->pm_gm / r;
		}
	}
	if(s->acc_on){
		*ax += s->acc[0]; *ay += s->acc[1]; *az += s->acc[2];
		*phi += -(s->acc[0]*x + s->acc[1]*y + s->acc[2]*z);
	}
}

typedef struct { double c[3], h, com[3], m; int child[8]; int first, cnt; } g3node;
typedef struct { g3node *nd; int nn, cap; int *idx; } g3tree;

static int g3_tree_new(g3tree *t){
	if(t->nn >= t->cap){
		t->cap = 2*t->cap + 1024;
		t->nd = (g3node*)realloc(t->nd, sizeof(g3node)*t->cap);
		if(!t->nd){ fprintf(stderr, "[3D] out of memory (tree)\n"); exit(12); }
	}
	memset(t->nd + t->nn, 0, sizeof(g3node));
	return t->nn++;
}

static void g3_tree_build_node(g3sim *s, g3tree *t, int nodeid, int depth){
	g3node *nd = t->nd + nodeid;
	int k, a, first = nd->first, cnt = nd->cnt;
	double M = 0, cm[3] = {0,0,0};
	for(k=first;k<first+cnt;k++){
		int i = t->idx[k];
		M += s->m[i]; cm[0] += s->m[i]*s->x[i]; cm[1] += s->m[i]*s->y[i]; cm[2] += s->m[i]*s->z[i];
	}
	nd->m = M;
	for(a=0;a<3;a++) nd->com[a] = (M > 0) ? cm[a]/M : nd->c[a];
	for(k=0;k<8;k++) nd->child[k] = -1;
	if(cnt <= 8 || depth > 60) return;
	/* partition into octants */
	int cntc[8] = {0}, start[8], pos[8];
	double c0 = nd->c[0], c1 = nd->c[1], c2 = nd->c[2];
	for(k=first;k<first+cnt;k++){
		int i = t->idx[k];
		int o = (s->x[i] >= c0) | ((s->y[i] >= c1) << 1) | ((s->z[i] >= c2) << 2);
		cntc[o]++;
	}
	start[0] = first;
	for(k=1;k<8;k++) start[k] = start[k-1] + cntc[k-1];
	for(k=0;k<8;k++) pos[k] = start[k];
	int *tmp = (int*)g3_alloc(sizeof(int)*cnt);
	for(k=first;k<first+cnt;k++){
		int i = t->idx[k];
		int o = (s->x[i] >= c0) | ((s->y[i] >= c1) << 1) | ((s->z[i] >= c2) << 2);
		tmp[pos[o]++ - first] = i;
	}
	memcpy(t->idx + first, tmp, sizeof(int)*cnt);
	free(tmp);
	double hh = 0.5 * nd->h;
	for(k=0;k<8;k++){
		if(cntc[k] == 0) continue;
		int cid = g3_tree_new(t);
		nd = t->nd + nodeid;            /* realloc may move */
		g3node *ch = t->nd + cid;
		ch->h = hh;
		ch->c[0] = c0 + ((k & 1) ? hh : -hh);
		ch->c[1] = c1 + ((k & 2) ? hh : -hh);
		ch->c[2] = c2 + ((k & 4) ? hh : -hh);
		ch->first = start[k]; ch->cnt = cntc[k];
		nd->child[k] = cid;
		g3_tree_build_node(s, t, cid, depth+1);
	}
}

/* Self-gravity (Plummer softened) + external hook into gx, gy, gz;
 * phi = self potential, phiext = external potential. */
static void g3_gravity(g3sim *s){
	long i;
	double G = s->G, e2 = s->eps*s->eps, th2 = s->theta*s->theta;
	if(!s->gravity && !s->pm_on && !s->acc_on){
		memset(s->gx, 0, sizeof(double)*s->n); memset(s->gy, 0, sizeof(double)*s->n);
		memset(s->gz, 0, sizeof(double)*s->n); memset(s->phi, 0, sizeof(double)*s->n);
		memset(s->phiext, 0, sizeof(double)*s->n);
		return;
	}
	g3tree t = {0};
	if(s->gravity && !s->grav_direct){
		double lo[3] = {1e300,1e300,1e300}, hi[3] = {-1e300,-1e300,-1e300};
		int a;
		t.idx = (int*)g3_alloc(sizeof(int)*s->n);
		for(i=0;i<s->n;i++){
			double p[3] = {s->x[i], s->y[i], s->z[i]};
			t.idx[i] = (int)i;
			for(a=0;a<3;a++){ if(p[a] < lo[a]) lo[a] = p[a]; if(p[a] > hi[a]) hi[a] = p[a]; }
		}
		int root = g3_tree_new(&t);
		double h = 0;
		for(a=0;a<3;a++){ t.nd[root].c[a] = 0.5*(lo[a]+hi[a]); if(hi[a]-lo[a] > h) h = hi[a]-lo[a]; }
		t.nd[root].h = 0.5*h*1.000001 + 1e-300;
		t.nd[root].first = 0; t.nd[root].cnt = (int)s->n;
		g3_tree_build_node(s, &t, root, 0);
	}
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 128)
#endif
	for(i=0;i<s->n;i++){
		double xi = s->x[i], yi = s->y[i], zi = s->z[i];
		double ax = 0, ay = 0, az = 0, ph = 0;
		if(s->gravity){
			if(s->grav_direct){
				long j;
				for(j=0;j<s->n;j++){
					if(j == i) continue;
					double dx = s->x[j]-xi, dy = s->y[j]-yi, dz = s->z[j]-zi;
					double r2 = dx*dx + dy*dy + dz*dz + e2;
					double ir = 1.0/sqrt(r2), ir3 = ir*ir*ir*s->m[j];
					ax += dx*ir3; ay += dy*ir3; az += dz*ir3; ph -= s->m[j]*ir;
				}
			} else {
				int stack[512], sp = 0;
				stack[sp++] = 0;
				while(sp > 0){
					g3node *nd = t.nd + stack[--sp];
					double dx = nd->com[0]-xi, dy = nd->com[1]-yi, dz = nd->com[2]-zi;
					double r2 = dx*dx + dy*dy + dz*dz;
					double s2 = 4.0*nd->h*nd->h;
					int inside = fabs(xi - nd->c[0]) <= nd->h && fabs(yi - nd->c[1]) <= nd->h
						&& fabs(zi - nd->c[2]) <= nd->h;
					int leaf = (nd->child[0] < 0 && nd->child[1] < 0 && nd->child[2] < 0 && nd->child[3] < 0
						&& nd->child[4] < 0 && nd->child[5] < 0 && nd->child[6] < 0 && nd->child[7] < 0);
					if(!inside && s2 < th2*r2){
						double rr = r2 + e2, ir = 1.0/sqrt(rr), ir3 = ir*ir*ir*nd->m;
						ax += dx*ir3; ay += dy*ir3; az += dz*ir3; ph -= nd->m*ir;
					} else if(leaf){
						int k;
						for(k=nd->first;k<nd->first+nd->cnt;k++){
							int j = t.idx[k];
							if(j == i) continue;
							double ddx = s->x[j]-xi, ddy = s->y[j]-yi, ddz = s->z[j]-zi;
							double rr = ddx*ddx + ddy*ddy + ddz*ddz + e2;
							double ir = 1.0/sqrt(rr), ir3 = ir*ir*ir*s->m[j];
							ax += ddx*ir3; ay += ddy*ir3; az += ddz*ir3; ph -= s->m[j]*ir;
						}
					} else {
						int k;
						for(k=0;k<8;k++) if(nd->child[k] >= 0){
							if(sp >= 512){ fprintf(stderr, "[3D] tree stack overflow\n"); exit(13); }
							stack[sp++] = nd->child[k];
						}
					}
				}
			}
			ax *= G; ay *= G; az *= G; ph *= G;
		}
		double ex, ey, ez, pe;
		g3_ext_accel(s, xi, yi, zi, &ex, &ey, &ez, &pe);
		s->gx[i] = ax + ex; s->gy[i] = ay + ey; s->gz[i] = az + ez;
		s->phi[i] = ph; s->phiext[i] = pe;
	}
	free(t.nd); free(t.idx);
}

/* ------------------------------------------------------------ face forces */
/* get3dUpqrad: face (anchor) velocity relative to v_i. Same formula as
 * get2dUpqradRk4 / get3dUpqradRk4 (Voro/voro_eunha.c); 0.5 (v_j - v_i)
 * when w2_i = w2_j and the weights do not change. */
static inline void g3_upqrad(const g3sim *s, long i, long j, const double upq[3],
		const double er[3], double d, double out[3]){
	double wp2 = s->w2[i], wq2 = s->w2[j];
	double d2 = d*d;
	double fact1 = 0.5*(1.0 + (wp2 - wq2)/d2);
	double vpw = 0, vqw = 0, fact2;
	if(s->dtold > 0){
		double dwp = (sqrt(wp2) - sqrt(s->w2old[i]))/s->dtold;
		double dwq = (sqrt(wq2) - sqrt(s->w2old[j]))/s->dtold;
		vpw = (dwp > 0) ? fmin(s->cs[i], dwp) : fmax(-s->cs[i], dwp);
		vqw = (dwq > 0) ? fmin(s->cs[j], dwq) : fmax(-s->cs[j], dwq);
	}
	double erdu = er[0]*upq[0] + er[1]*upq[1] + er[2]*upq[2];
	fact2 = (sqrt(wp2)*vpw - sqrt(wq2)*vqw)/d - (wp2 - wq2)/d2 * erdu;
	out[0] = fact1*upq[0] + fact2*er[0];
	out[1] = fact1*upq[1] + fact2*er[1];
	out[2] = fact1*upq[2] + fact2*er[2];
}

/* Hydro accelerations, dE/dt (die) and the hydro time step
 * (getAccVoro2DBlend_impl, av_mode 5, 3D). g[xyz] must be current. */
static double g3_forces(g3sim *s){
	long i;
	double Dt = 1e30;
	long long npair = 0;
	const double C = s->courant, gam = s->gamma;
	const int phase1 = s->phase1, muscl = s->use_muscl;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 128) reduction(min:Dt) reduction(+:npair)
#endif
	for(i=0;i<s->n;i++){
		int q, a, b, nf = s->fcnt[i];
		g3face *F = g3_faces(s, i);
		const double vi[3] = {s->vx[i], s->vy[i], s->vz[i]};
		const double Pi = s->P[i], rhoi = s->den[i], csi = s->cs[i], Vi = s->vol[i];
		const double *gi = s->grad + 15*i;
		double Fx[3] = {0,0,0}, dte = 0, die = 0, dti = 1e30;
		for(q=0;q<nf;q++){
			g3face *f = F + q;
			long j = f->j;
			int ghost = (f->mf != 0);
			double r[3], sg[3];
			g3_image_signs(s, f->mf, r, sg);
			const double vj[3] = {r[0]*s->vx[j], r[1]*s->vy[j], r[2]*s->vz[j]};
			const double Pj = s->P[j], rhoj = s->den[j], csj = s->cs[j];
			double d = sqrt(f->rel[0]*f->rel[0] + f->rel[1]*f->rel[1] + f->rel[2]*f->rel[2]);
			double area = sqrt(f->S[0]*f->S[0] + f->S[1]*f->S[1] + f->S[2]*f->S[2]);
			if(!(d > 0) || !(area > 0)) continue;
			double er[3] = {f->rel[0]/d, f->rel[1]/d, f->rel[2]/d};
			double nh[3] = {f->S[0]/area, f->S[1]/area, f->S[2]/area};
			double upq[3] = {vj[0]-vi[0], vj[1]-vi[1], vj[2]-vi[2]};
			double urad[3];
			g3_upqrad(s, i, j, upq, er, d, urad);
			/* Springel face rotation (voronoi_face_rotation, 3D) */
			{
				double w2i = s->w2[i], w2j = s->w2[j], dx[3], delta;
				int use = 1;
				if(w2i == 0 && w2j == 0){
					for(a=0;a<3;a++) dx[a] = f->c[a] - 0.5*f->rel[a];
				} else if(s->lag_rot){
					double fact1 = 0.5*(1.0 + (w2i - w2j)/(d*d));
					for(a=0;a<3;a++) dx[a] = f->c[a] - fact1*f->rel[a];
				} else use = 0;
				if(use){
					delta = -(upq[0]*dx[0] + upq[1]*dx[1] + upq[2]*dx[2]) / d;
					for(a=0;a<3;a++) urad[a] += delta*er[a];
				}
			}
			double wn = (vi[0]+urad[0])*nh[0] + (vi[1]+urad[1])*nh[1] + (vi[2]+urad[2])*nh[2];
			double vnL = vi[0]*nh[0] + vi[1]*nh[1] + vi[2]*nh[2];
			double vnR = vj[0]*nh[0] + vj[1]*nh[1] + vj[2]*nh[2];
			double pL = Pi, pR = Pj, rhoL = rhoi, rhoR = rhoj;
			if(muscl){
				const double *gj = s->grad + 15*j;
				double xiF[3] = {f->c[0], f->c[1], f->c[2]};
				double xjF[3] = {f->c[0]-f->rel[0], f->c[1]-f->rel[1], f->c[2]-f->rel[2]};
				double drl = 0, drr = 0, dpl = 0, dpr = 0, dvl = 0, dvr = 0;
				for(b=0;b<3;b++){
					drl += gi[b]*xiF[b];
					drr += sg[b]*gj[b]*xjF[b];
					dpl += gi[3+b]*xiF[b];
					dpr += sg[b]*gj[3+b]*xjF[b];
					for(a=0;a<3;a++){
						dvl += nh[a]*gi[6+3*a+b]*xiF[b];
						dvr += nh[a]*r[a]*sg[b]*gj[6+3*a+b]*xjF[b];
					}
				}
				rhoL += drl; rhoR += drr; pL += dpl; pR += dpr; vnL += dvl; vnR += dvr;
			}
			if(pL < 1e-10) pL = 1e-10;
			if(pR < 1e-10) pR = 1e-10;
			if(rhoL < 1e-10) rhoL = 1e-10;
			if(rhoR < 1e-10) rhoR = 1e-10;
			double pst, vnst, pi_total;
			gfs_hllc_face_rest_frame(rhoL, pL, vnL, csi, rhoR, pR, vnR, csj, wn, gam, phase1, &pst, &vnst);
			pi_total = pst;
			int riem = 0;
			double rvn = 0;
			if(phase1 && d > 0){
				double pmin_c = fmin(Pi, Pj), pmax_c = fmax(Pi, Pj);
				if(pmax_c > 100.0*fmax(pmin_c, 1.0e-30)){
					double ps, vns;
					gfs_hll_star_state(rhoi, Pi, vi[0]*er[0]+vi[1]*er[1]+vi[2]*er[2], csi,
							rhoj, Pj, vj[0]*er[0]+vj[1]*er[1]+vj[2]*er[2], csj, gam, &ps, &vns);
					pi_total = ps; riem = 1; rvn = vns;
				}
			}
			double rvel = upq[0]*er[0] + upq[1]*er[1] + upq[2]*er[2];
			if(s->alphavis > 0 && rvel < 0){
				double meanden = 0.5*(rhoi + rhoj), meancs = 0.5*(csi + csj);
				pi_total += (-s->alphavis*meancs*rvel + s->betavis*rvel*rvel)*meanden;
			}
			if(!ghost){
				double vclose = -rvel;
				if(vclose < 0) vclose = 0;
				double pp = gfs_pair_pressure_len(d, cbrt(Vi), cbrt(s->vol[j]),
						rhoi, rhoj, csi, csj, vclose);
				pp = gfs_pair_work_limit(pp, area, s->dtold, s->m[i], s->m[j], rvel,
						s->ie[i], s->ie[j], s->nshare);
				if(pp > 0){ pi_total += pp; npair++; }
			}
			double ua[3];
			if(riem){ for(a=0;a<3;a++) ua[a] = rvn*er[a]; }
			else { for(a=0;a<3;a++) ua[a] = vi[a] + urad[a]; }
			die += -pi_total*(urad[0]*f->S[0] + urad[1]*f->S[1] + urad[2]*f->S[2]);
			dte += -pi_total*(ua[0]*f->S[0] + ua[1]*f->S[1] + ua[2]*f->S[2]);
			for(a=0;a<3;a++) Fx[a] += -pi_total*f->S[a];
			/* time step */
			double vsig = csi + csj - fmin(0.0, rvel);
			double heff = 0.25*cbrt(Vi);
			double dcfl = ghost ? fmax(d, heff) : d;
			double dt = 2.0*C*dcfl/vsig;
			double du = sqrt(upq[0]*upq[0] + upq[1]*upq[1] + upq[2]*upq[2]);
			if(du > 0){
				double dt3 = 0.1*dcfl/du;
				if(dt3 < dt) dt = dt3;
			}
			if(isnan(dt)) dt = 1e-10;
			if(dt < dti) dti = dt;
		}
		s->die[i] = phase1 ? dte : die;
		s->dK[i] = 0;             /* av_mode 5: no dissipative dK channel */
		s->ax[i] = Fx[0]/s->m[i]; s->ay[i] = Fx[1]/s->m[i]; s->az[i] = Fx[2]/s->m[i];
		/* acceleration CFL with the full summed acceleration (hydro + gravity) */
		{
			double atx = s->ax[i] + s->gx[i], aty = s->ay[i] + s->gy[i], atz = s->az[i] + s->gz[i];
			double amag = sqrt(atx*atx + aty*aty + atz*atz);
			if(amag > 0 && Vi > 0){
				double dta = 0.25*sqrt(cbrt(Vi)/amag);
				if(dta < dti) dti = dta;
			}
		}
		s->dtp[i] = dti;
		if(dti < Dt) Dt = dti;
	}
	s->npair += npair;
	return Dt;
}

/* One RK stage evaluation at the current positions. */
static int g3_eval(g3sim *s, double *Dt){
	if(!s->geom_valid){
		int nfail = g3_tessellate(s);
		if(nfail){ fprintf(stderr, "[3D] %d cells failed\n", nfail); return -1; }
	}
	s->geom_valid = 0;
	g3_primitives(s);
	g3_gradients(s);
	if(!s->grav_valid) g3_gravity(s);
	s->grav_valid = 0;
	double d = g3_forces(s);
	if(Dt) *Dt = d;
	return 0;
}

/* ------------------------------------------------------------ boundaries */
static void g3_wrap(g3sim *s){
	long i;
	int a;
	for(a=0;a<3;a++){
		if(s->bc[a] != BC_PERIODIC) continue;
		double *p = (a == 0) ? s->x : (a == 1) ? s->y : s->z;
		double lo = s->lo[a], L = s->L[a];
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(i=0;i<s->n;i++){
			if(p[i] < lo || p[i] >= lo + L){
				p[i] -= L*floor((p[i] - lo)/L);
				if(p[i] >= lo + L) p[i] -= L;
				if(p[i] < lo) p[i] = lo;
			}
		}
	}
}

/* End of a full step: periodic wrap, specular reflection at reflecting
 * walls (energy conserving), removal of particles outside outflow walls. */
static void g3_end_boundaries(g3sim *s){
	long i, k;
	int a;
	g3_wrap(s);
	for(a=0;a<3;a++){
		if(s->bc[a] != BC_REFLECT) continue;
		double *p = (a == 0) ? s->x : (a == 1) ? s->y : s->z;
		double *v = (a == 0) ? s->vx : (a == 1) ? s->vy : s->vz;
		for(i=0;i<s->n;i++){
			if(p[i] < s->lo[a]){ p[i] = 2*s->lo[a] - p[i]; v[i] = -v[i]; }
			else if(p[i] > s->hi[a]){ p[i] = 2*s->hi[a] - p[i]; v[i] = -v[i]; }
		}
	}
	if(s->bc[0] != BC_OUTFLOW && s->bc[1] != BC_OUTFLOW && s->bc[2] != BC_OUTFLOW) return;
	double eout = 0;
	for(i=0,k=0;i<s->n;i++){
		double p[3] = {s->x[i], s->y[i], s->z[i]};
		int out = 0;
		for(a=0;a<3;a++) if(s->bc[a] == BC_OUTFLOW && (p[a] < s->lo[a] || p[a] > s->hi[a])) out = 1;
		if(out){
			eout += s->ie[i] + 0.5*s->m[i]*(s->vx[i]*s->vx[i] + s->vy[i]*s->vy[i] + s->vz[i]*s->vz[i]);
			continue;
		}
		if(k != i){
			double **f[] = {&s->x,&s->y,&s->z,&s->vx,&s->vy,&s->vz,&s->m,&s->ie,&s->vol,&s->den,
				&s->P,&s->cs,&s->K,&s->w2,&s->w2old,&s->ax,&s->ay,&s->az,&s->gx,&s->gy,&s->gz,
				&s->phi,&s->phiext,&s->die,&s->dK,&s->ie0,&s->V0,&s->dtp};
			size_t q;
			for(q=0;q<sizeof(f)/sizeof(f[0]);q++) (*f[q])[k] = (*f[q])[i];
			memcpy(s->grad + 15*k, s->grad + 15*i, sizeof(double)*15);
			s->id[k] = s->id[i];
		}
		k++;
	}
	if(k != s->n){
		fprintf(stderr, "[3D] outflow removed %ld particles, E_out=%.6e\n", s->n - k, eout);
		s->out_cum += eout;
		s->n = k;
		s->grav_valid = 0;
	}
}

/* ------------------------------------------------------------ RK4 */
#define K3(s,i,st,f) ((s)->k[((size_t)(i)*4 + (st))*8 + (f)])

static void g3_derivs(g3sim *s, int st, double Dt){
	long i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<s->n;i++){
		K3(s,i,st,0) = s->vx[i]*Dt;
		K3(s,i,st,1) = s->vy[i]*Dt;
		K3(s,i,st,2) = s->vz[i]*Dt;
		/* self-gravity and the external hook enter v only (as kepler_accel) */
		K3(s,i,st,3) = (s->ax[i] + s->gx[i])*Dt;
		K3(s,i,st,4) = (s->ay[i] + s->gy[i])*Dt;
		K3(s,i,st,5) = (s->az[i] + s->gz[i])*Dt;
		if(s->phase1){
			/* phase1_ie_stage: die holds dE/dt; subtract the hydro kinetic power */
			double dke = s->m[i]*(s->vx[i]*s->ax[i] + s->vy[i]*s->ay[i] + s->vz[i]*s->az[i]);
			K3(s,i,st,6) = (s->die[i] - dke)*Dt;
		} else {
			K3(s,i,st,6) = s->die[i]*Dt;
		}
		K3(s,i,st,7) = s->dK[i]*Dt;
	}
}

/* x += a*k[sa] + b*k[sb] for all 8 variables */
static void g3_shift(g3sim *s, int sa, double ca, int sb, double cb){
	long i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<s->n;i++){
		double d[8];
		int f;
		for(f=0;f<8;f++) d[f] = ca*K3(s,i,sa,f) + (sb >= 0 ? cb*K3(s,i,sb,f) : 0.0);
		s->x[i] += d[0]; s->y[i] += d[1]; s->z[i] += d[2];
		s->vx[i] += d[3]; s->vy[i] += d[4]; s->vz[i] += d[5];
		s->ie[i] += d[6]; s->K[i] += d[7];
	}
}

static void g3_energies(const g3sim *s, double *ek, double *ei, double *ep){
	long i;
	double a = 0, b = 0, c = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:a,b,c)
#endif
	for(i=0;i<s->n;i++){
		a += 0.5*s->m[i]*(s->vx[i]*s->vx[i] + s->vy[i]*s->vy[i] + s->vz[i]*s->vz[i]);
		b += s->ie[i];
		c += 0.5*s->m[i]*s->phi[i] + s->m[i]*s->phiext[i];
	}
	*ek = a; *ei = b; *ep = c;
}

static void g3_center(const g3sim *s, double c[3]){
	int a;
	for(a=0;a<3;a++) c[a] = s->pm_on ? s->pm_x[a] : 0.5*(s->lo[a] + s->hi[a]);
}

/* One RK4 step (exam2d_vph_rk4_int_blend). Returns dt, < 0 on failure. */
static double g3_step(g3sim *s, double dtcap){
	long i;
	double Dt, dum;
	double gm1 = s->gamma - 1.0;
	s->sfl_step = 0;
	s->sfl_nlog = 0;
	s->npair = 0;
	for(i=0;i<s->n;i++) s->w2old[i] = s->w2[i];
	/* K1 */
	if(g3_eval(s, &Dt)) return -1;
	int es_step = s->es_on || s->half_on;
	if(es_step){
		for(i=0;i<s->n;i++){ s->ie0[i] = s->ie[i]; s->V0[i] = s->vol[i]; }
	}
	if(!(Dt > 0) || !isfinite(Dt)){ fprintf(stderr, "[3D] bad dt %g\n", Dt); return -1; }
	if(dtcap > 0 && Dt > dtcap) Dt = dtcap;
	g3_derivs(s, 0, Dt);
	/* K2 at x0 + k1/2 */
	g3_shift(s, 0, 0.5, -1, 0); g3_wrap(s);
	if(g3_eval(s, &dum)) return -1;
	g3_derivs(s, 1, Dt);
	/* K3 at x0 + k2/2 */
	g3_shift(s, 1, 0.5, 0, -0.5); g3_wrap(s);
	if(g3_eval(s, &dum)) return -1;
	g3_derivs(s, 2, Dt);
	/* K4 at x0 + k3 */
	g3_shift(s, 2, 1.0, 1, -0.5); g3_wrap(s);
	if(g3_eval(s, &dum)) return -1;
	g3_derivs(s, 3, Dt);
	/* undo K4 shift, final combination */
	g3_shift(s, 2, -1.0, -1, 0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<s->n;i++){
		double d[8];
		int f;
		for(f=0;f<8;f++) d[f] = (K3(s,i,0,f) + 2*K3(s,i,1,f) + 2*K3(s,i,2,f) + K3(s,i,3,f))/6.0;
		s->x[i] += d[0]; s->y[i] += d[1]; s->z[i] += d[2];
		s->vx[i] += d[3]; s->vy[i] += d[4]; s->vz[i] += d[5];
		s->ie[i] += d[6]; s->K[i] += d[7];
	}
	g3_end_boundaries(s);
	/* new volumes (exam2dUpdateVol); the faces are reused by the next K1 */
	{
		int nfail = g3_tessellate(s);
		if(nfail){ fprintf(stderr, "[3D] %d cells failed (final)\n", nfail); return -1; }
		s->geom_valid = 1;
	}
	/* GFS_ENTROPY_SWITCH / GFS_HALF_LIMIT on the final state */
	double es_de = 0, hl_de = 0;
	long es_n = 0, hl_n = 0;
	if(es_step){
		for(i=0;i<s->n;i++){
			double V = s->vol[i], V0 = s->V0[i], ie0 = s->ie0[i];
			double ieE = s->ie[i], ie1 = ieE;
			if(!(V > 0) || !(V0 > 0) || !(ie0 > 0)) continue;
			if(s->es_on){
				double g2 = s->gx[i]*s->gx[i] + s->gy[i]*s->gy[i] + s->gz[i]*s->gz[i];
				double eg = s->m[i]*sqrt(g2)*cbrt(V);
				double eth = (0.5*ie0 > ieE) ? 0.5*ie0 : ieE;
				if(s->es_coef*eg > eth){
					ie1 = ie0*pow(V0/V, gm1);
					es_de += ie1 - ieE; es_n++;
				}
			}
			if(s->half_on && ie1 < 0.5*ie0){
				hl_de += 0.5*ie0 - ie1; hl_n++;
				ie1 = 0.5*ie0;
			}
			s->ie[i] = ie1;
		}
	}
	/* GFS_DUAL_ENERGY and GFS_FLOOR_LOG */
	double de_inj = 0;
	long de_n = 0;
	if(s->de_eta > 0 || s->floor_log){
		int nlog = 0;
		double c0[3];
		g3_center(s, c0);
		for(i=0;i<s->n;i++){
			if(!(s->vol[i] > 0)) continue;
			double den = s->m[i]/s->vol[i];
			double ieK = (s->K[i] > 0) ? s->K[i]*pow(den, s->gamma)*s->vol[i]/gm1 : 0.0;
			int reset = (s->de_eta > 0 && ieK > 0 && s->ie[i] < s->de_eta*ieK);
			if(s->floor_log && (reset || s->ie[i] <= 0) && nlog < 8){
				double dx = s->x[i]-c0[0], dy = s->y[i]-c0[1], dz = s->z[i]-c0[2];
				fprintf(stderr, "[RK4F] r0 x=%.5f y=%.5f z=%.5f R=%.4f den=%.4e ie=%.4e ie_K=%.4e ie/ie_K=%.3e vol=%.4e reset=%d\n",
						s->x[i], s->y[i], s->z[i], sqrt(dx*dx+dy*dy+dz*dz), den, s->ie[i], ieK,
						(ieK > 0 ? s->ie[i]/ieK : 0.0), s->vol[i], reset);
				nlog++;
			}
			if(reset){ de_inj += ieK - s->ie[i]; de_n++; s->ie[i] = ieK; }
		}
	}
	/* final floor (booked) and primitives */
	double fl = 0;
	long nneg = 0;
	for(i=0;i<s->n;i++){
		if(s->ie[i] <= 0){ nneg++; fl += 1e-6*s->vol[i]/gm1 - s->ie[i]; }
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<s->n;i++){
		s->den[i] = s->m[i]/s->vol[i];
		s->P[i] = s->ie[i]/s->vol[i]*gm1;
		if(s->P[i] <= 0){ s->P[i] = 1e-6; s->ie[i] = s->P[i]*s->vol[i]/gm1; }
		s->cs[i] = sqrt(s->gamma*s->P[i]/s->den[i]);
	}
	/* gravity at the new positions: Epot for the log, reused by the next K1 */
	g3_gravity(s);
	s->grav_valid = 1;
	s->fl_cum += fl; s->de_cum += de_inj; s->es_cum += es_de; s->hl_cum += hl_de;
	s->sfl_cum += s->sfl_step;
	s->step++;
	s->time += Dt;
	s->dtold = Dt;
	{
		double ek, ei, ep;
		static double Et0 = 0, Eh0 = 0;
		static int first = 1;
		g3_energies(s, &ek, &ei, &ep);
		double eh = ek + ei, et = eh + ep;
		if(first){ Et0 = s->E0tot_set ? s->E0tot : et; Eh0 = s->E0tot_set ? s->E0hyd : eh; first = 0; }
		fprintf(stderr,
			"[RK4E] step=%ld dt=%.3e E_hyd=%.6e E_pot=%.6e E_tot=%.6e dEtot/|Etot0|=%.3e dEhyd/Ehyd0=%.3e n_ie_le0=%ld floor_inj=%.3e floor_cum=%.3e n_pair=%lld E_int=%.6e n_de=%ld de_inj=%.3e de_cum=%.3e sfl=%.3e sfl_cum=%.3e",
			s->step, Dt, eh, ep, et, (Et0 != 0 ? (et-Et0)/fabs(Et0) : 0), (Eh0 != 0 ? (eh-Eh0)/Eh0 : 0),
			nneg, fl, s->fl_cum, s->npair, ei, de_n, de_inj, s->de_cum, s->sfl_step, s->sfl_cum);
		if(es_step)
			fprintf(stderr, " n_es=%ld es_de=%.3e es_cum=%.3e n_half=%ld hl_de=%.3e hl_cum=%.3e",
					es_n, es_de, s->es_cum, hl_n, hl_de, s->hl_cum);
		if(s->out_cum != 0) fprintf(stderr, " out_cum=%.3e", s->out_cum);
		fprintf(stderr, "\n");
		printf("[E3D] step= %ld t= %.9e dt= %.6e Ekin= %.12e Eint= %.12e Epot= %.12e Etot= %.12e dE_rel= %.4e floor_cum= %.4e sfl_cum= %.4e de_cum= %.4e es_cum= %.4e hl_cum= %.4e out_cum= %.4e npair= %lld N= %ld Rmax= %d njit= %ld\n",
				s->step, s->time, Dt, ek, ei, ep, et, (Et0 != 0 ? (et-Et0)/fabs(Et0) : 0),
				s->fl_cum, s->sfl_cum, s->de_cum, s->es_cum, s->hl_cum, s->out_cum, s->npair, s->n, s->max_R, s->n_jitter);
		fflush(stdout);
	}
	return Dt;
}

/* ------------------------------------------------------------ output */
static int g3_write_snap(const g3sim *s, int idx, const lag3d_header *h0){
	char name[64];
	lag3d_data d;
	long i, n = s->n;
	memset(&d, 0, sizeof(d));
	d.h = *h0;
	d.h.version = 1;
	d.h.np = n;
	d.h.time = s->time;
	d.h.gamma = s->gamma;
	d.h.flags = LAG3D_FLAG_VOL | (s->gravity ? LAG3D_FLAG_POT : 0);
	d.x = s->x; d.y = s->y; d.z = s->z; d.vx = s->vx; d.vy = s->vy; d.vz = s->vz; d.mass = s->m;
	d.u = (double*)g3_alloc(sizeof(double)*n);
	d.rho = (double*)g3_alloc(sizeof(double)*n);
	d.vol = s->vol;
	d.pot = s->gravity ? s->phi : NULL;
	d.id = s->id;
	for(i=0;i<n;i++){ d.u[i] = s->ie[i]/s->m[i]; d.rho[i] = s->m[i]/s->vol[i]; }
	snprintf(name, sizeof(name), "snap_%06d.l3d", idx);
	int rc = lag3d_write(name, &d);
	printf("[3D] wrote %s t=%.6e N=%ld rc=%d\n", name, s->time, n, rc);
	fflush(stdout);
	free(d.u); free(d.rho);
	return rc;
}

static const char *g3_bcname(int b){
	return b == BC_PERIODIC ? "periodic" : b == BC_REFLECT ? "reflect" : "outflow";
}

static int g3_parse_bc(const char *v, int bc[3]){
	char buf[256], *tok, *save = NULL;
	int k = 0, a;
	strncpy(buf, v, sizeof(buf)-1); buf[sizeof(buf)-1] = 0;
	for(tok = strtok_r(buf, ", ", &save); tok && k < 3; tok = strtok_r(NULL, ", ", &save)){
		if(strcmp(tok, "periodic") == 0) bc[k++] = BC_PERIODIC;
		else if(strcmp(tok, "reflect") == 0 || strcmp(tok, "reflecting") == 0) bc[k++] = BC_REFLECT;
		else if(strcmp(tok, "outflow") == 0) bc[k++] = BC_OUTFLOW;
		else return -1;
	}
	if(k == 1){ bc[1] = bc[2] = bc[0]; k = 3; }
	for(a=k;a<3;a++) return -1;
	return 0;
}

/* ------------------------------------------------------------ driver */
int exam3d_gfs_run(const char *pf, int myid, int nranks){
	g3sim S, *s = &S;
	char v[512];
	int a;
	if(myid != 0) return 0;
	memset(s, 0, sizeof(*s));
	G3S = s;
	printf("LAGEUNHA_3D_GFS_V1: 3D GFS hydro (CPU/OpenMP), params %s\n", pf);
	if(nranks > 1)
		printf("[3D] note: the 3D path runs on rank 0 only; %d other ranks idle (use -np 1 and OMP_NUM_THREADS)\n", nranks-1);
	/* parameters */
	s->gamma = g3_pd(pf, "Gamma_hydro", 5.0/3.0);
	s->courant = g3_pd(pf, "Courant_hydro", 0.3);
	s->av_mode = (int)g3_pd(pf, "GAS av_mode", 5);
	s->use_muscl = (int)g3_pd(pf, "GAS use_muscl", 1);
	s->kappa = g3_pd(pf, "GAS kappa", 0);
	s->alphavis = g3_pd(pf, "GAS AlphaVis", 0);
	s->betavis = g3_pd(pf, "GAS BetaVis", 2.0*s->alphavis);
	s->tend = g3_pd(pf, "Hydro3D t_end", 0);
	s->dumpdt = g3_pd(pf, "Hydro3D dump dt", 0);
	s->gravity = (int)g3_pd(pf, "Hydro3D gravity", 0);
	s->eps = g3_pd(pf, "Hydro3D softening", 0);
	s->theta = g3_pd(pf, "Hydro3D theta", g3_env_d("LAG3D_THETA", 0.5));
	s->maxsteps = (int)g3_pd(pf, "Hydro3D max steps", g3_env_d("LAG3D_MAXSTEPS", 2000000000.0));
	if(getenv("HYDRO_TSTOP") && atof(getenv("HYDRO_TSTOP")) > 0) s->tend = atof(getenv("HYDRO_TSTOP"));
	if(getenv("EUNHA_DUMP_DT") && atof(getenv("EUNHA_DUMP_DT")) > 0) s->dumpdt = atof(getenv("EUNHA_DUMP_DT"));
	{
		int way = (int)g3_pd(pf, "Hydro time-stepping way", 1);
		int em = (int)g3_pd(pf, "GAS entropy_mode", 0);
		double fc = g3_pd(pf, "Voro centroid shift factor", 0);
		if(way != 1) printf("[3D] warning: Hydro time-stepping way = %d; the 3D path is RK4 only (way 1)\n", way);
		if(s->av_mode != 5){ fprintf(stderr, "[3D] GAS av_mode = %d not supported (3D path is av_mode 5)\n", s->av_mode); return 2; }
		if(em != 0){ fprintf(stderr, "[3D] GAS entropy_mode = %d not supported (ie only)\n", em); return 2; }
		if(s->kappa > 0){ fprintf(stderr, "[3D] GAS kappa > 0 (adaptive Laguerre weights) not supported in 3D\n"); return 2; }
		if(fc > 0) printf("[3D] warning: centroid shift %g ignored in 3D\n", fc);
		if(g3_pd(pf, "GAS gpu_enabled", 0) != 0) printf("[3D] warning: gpu_enabled ignored (no 3D GPU kernel)\n");
	}
	if(!g3_param(pf, "Hydro3D IC file", s->icfile, sizeof(s->icfile)) || !s->icfile[0]){
		fprintf(stderr, "[3D] missing 'Hydro3D IC file'\n"); return 2;
	}
	/* switches */
	s->phase1 = g3_env_on("SEDOV_PHASE1");
	s->floor_log = g3_env_on("GFS_FLOOR_LOG");
	s->lag_rot = g3_env_on("GFS_LAGUERRE_ROTATION");
	s->de_eta = g3_env_d("GFS_DUAL_ENERGY", 0); if(!(s->de_eta > 0)) s->de_eta = 0;
	{
		const char *e = getenv("GFS_ENTROPY_SWITCH");
		s->es_on = (e && e[0] && atoi(e) > 0) ? 1 : 0;
		s->es_coef = 0.01;
		e = getenv("GFS_ES_COEF"); if(e && e[0] && atof(e) > 0) s->es_coef = atof(e);
		e = getenv("GFS_HALF_LIMIT");
		s->half_on = (e && e[0]) ? (atoi(e) > 0 ? 1 : 0) : s->es_on;
	}
	s->nshare = g3_env_d("GFS3D_PAIR_NSHARE", 16.0);
	s->stage_reset = (int)g3_env_d("GFS3D_STAGE_FLOOR", 0);
	/* LAG3D_GRAV_DIRECT: 1 direct summation, 0 tree; unset = direct for
	 * N <= 20000 (exact pair symmetry, cheap), tree above. Set after the IC. */
	s->grav_direct = getenv("LAG3D_GRAV_DIRECT") && getenv("LAG3D_GRAV_DIRECT")[0]
		? atoi(getenv("LAG3D_GRAV_DIRECT")) : -1;
	s->pm_gm = g3_env_d("LAG3D_PM_GM", 0); s->pm_on = s->pm_gm != 0;
	s->pm_x[0] = g3_env_d("LAG3D_PM_X", 0); s->pm_x[1] = g3_env_d("LAG3D_PM_Y", 0); s->pm_x[2] = g3_env_d("LAG3D_PM_Z", 0);
	s->pm_eps = g3_env_d("LAG3D_PM_EPS", 0);
	if(getenv("LAG3D_ACC")){
		if(sscanf(getenv("LAG3D_ACC"), "%lf,%lf,%lf", s->acc, s->acc+1, s->acc+2) == 3)
			s->acc_on = (s->acc[0] != 0 || s->acc[1] != 0 || s->acc[2] != 0);
	}
	/* initial conditions */
	lag3d_data ic;
	int rc = lag3d_read(s->icfile, &ic);
	if(rc){ fprintf(stderr, "[3D] cannot read LAG3DV1 file %s (rc=%d)\n", s->icfile, rc); return 3; }
	if(fabs(ic.h.gamma - s->gamma) > 1e-5)
		printf("[3D] warning: IC gamma %.6f differs from Gamma_hydro %.6f; using Gamma_hydro\n", ic.h.gamma, s->gamma);
	for(a=0;a<3;a++){
		s->lo[a] = ic.h.box[2*a]; s->hi[a] = ic.h.box[2*a+1]; s->L[a] = s->hi[a] - s->lo[a];
		if(!(s->L[a] > 0)){ fprintf(stderr, "[3D] bad box in IC\n"); return 3; }
		s->bc[a] = ic.h.periodic[a] ? BC_PERIODIC : BC_REFLECT;
	}
	if(g3_param(pf, "Hydro3D boundary", v, sizeof(v)) && v[0]){
		if(g3_parse_bc(v, s->bc)){ fprintf(stderr, "[3D] bad 'Hydro3D boundary = %s'\n", v); return 2; }
	}
	if(getenv("LAG3D_BC") && getenv("LAG3D_BC")[0]){
		if(g3_parse_bc(getenv("LAG3D_BC"), s->bc)){ fprintf(stderr, "[3D] bad LAG3D_BC\n"); return 2; }
	}
	s->G = (ic.h.G > 0) ? ic.h.G : 1.0;
	if(!(s->eps > 0)) s->eps = ic.h.softening;
	if(s->gravity && (s->bc[0] == BC_PERIODIC || s->bc[1] == BC_PERIODIC || s->bc[2] == BC_PERIODIC)){
		fprintf(stderr, "[3D] self-gravity is isolated only; periodic axes are not supported with gravity\n");
		return 2;
	}
	s->n = ic.h.np;
	g3_alloc_particles(s, s->n);
	{
		long i;
		double w2c = (s->kappa < 0) ? -s->kappa : 0.0;
		for(i=0;i<s->n;i++){
			s->x[i] = ic.x[i]; s->y[i] = ic.y[i]; s->z[i] = ic.z[i];
			s->vx[i] = ic.vx[i]; s->vy[i] = ic.vy[i]; s->vz[i] = ic.vz[i];
			s->m[i] = ic.mass[i];
			s->ie[i] = ic.mass[i]*ic.u[i];     /* u is specific */
			s->w2[i] = s->w2old[i] = w2c;
			s->id[i] = ic.id[i];
		}
	}
	if(s->grav_direct < 0) s->grav_direct = (s->n <= 20000) ? 1 : 0;
	lag3d_header h0 = ic.h;
	s->time = ic.h.time;
	lag3d_free(&ic);
#ifdef _OPENMP
	s->nthreads = omp_get_max_threads();
#else
	s->nthreads = 1;
#endif
	s->tb = (g3fbuf*)g3_alloc(sizeof(g3fbuf)*s->nthreads);
	printf("[3D] N=%ld box=[%g,%g]x[%g,%g]x[%g,%g] bc=%s,%s,%s gamma=%.6f courant=%g t_end=%g dump_dt=%g\n",
			s->n, s->lo[0], s->hi[0], s->lo[1], s->hi[1], s->lo[2], s->hi[2],
			g3_bcname(s->bc[0]), g3_bcname(s->bc[1]), g3_bcname(s->bc[2]),
			s->gamma, s->courant, s->tend, s->dumpdt);
	printf("[3D] av_mode=%d muscl=%d w2=%g phase1=%d floor_log=%d dual_energy=%g entropy_switch=%d(coef %g) half_limit=%d lag_rot=%d alphavis=%g pair_nshare=%g stage_floor_reset=%d threads=%d\n",
			s->av_mode, s->use_muscl, s->w2[0], s->phase1, s->floor_log, s->de_eta, s->es_on, s->es_coef,
			s->half_on, s->lag_rot, s->alphavis, s->nshare, s->stage_reset, s->nthreads);
	printf("[3D] gravity=%d G=%g eps=%g theta=%g direct=%d point_mass=%d(GM %g at %g %g %g eps %g) uniform_acc=%d\n",
			s->gravity, s->G, s->eps, s->theta, s->grav_direct, s->pm_on, s->pm_gm,
			s->pm_x[0], s->pm_x[1], s->pm_x[2], s->pm_eps, s->acc_on);
	if(!(s->tend > s->time)){ fprintf(stderr, "[3D] t_end %g <= t0 %g\n", s->tend, s->time); return 2; }
	/* initial state */
	g3_wrap(s);
	if(g3_tessellate(s)){ fprintf(stderr, "[3D] initial tessellation failed\n"); return 4; }
	s->geom_valid = 1;
	g3_primitives(s);
	if(s->sfl_step != 0) printf("[3D] initial floor injected %.3e\n", s->sfl_step);
	s->sfl_step = 0;
	if(s->de_eta > 0){
		long i, nk = 0;
		for(i=0;i<s->n;i++) if(s->K[i] <= 0 && s->den[i] > 0 && s->P[i] > 0){
			s->K[i] = s->P[i]/pow(s->den[i], s->gamma); nk++;
		}
		fprintf(stderr, "[DUALE] eta=%g K set from P/rho^gamma on %ld particles\n", s->de_eta, nk);
	}
	g3_gravity(s);
	s->grav_valid = 1;
	{
		long i;
		double vt = 0, vbox = s->L[0]*s->L[1]*s->L[2], ek, ei, ep;
		for(i=0;i<s->n;i++) vt += s->vol[i];
		g3_energies(s, &ek, &ei, &ep);
		s->E0tot = ek + ei + ep; s->E0hyd = ek + ei; s->E0tot_set = 1;
		printf("[3D] initial: sum V=%.12e box V=%.12e (ratio-1=%.3e) cells expanded=%ld maxR=%d\n",
				vt, vbox, vt/vbox - 1.0, s->n_expand, s->max_R);
		printf("[E3D] step= 0 t= %.9e dt= 0 Ekin= %.12e Eint= %.12e Epot= %.12e Etot= %.12e dE_rel= 0 N= %ld\n",
				s->time, ek, ei, ep, ek+ei+ep, s->n);
	}
	int isnap = 0;
	g3_write_snap(s, isnap++, &h0);
	double tnext = (s->dumpdt > 0) ? s->time + s->dumpdt : s->tend;
	if(tnext > s->tend) tnext = s->tend;
	double wall0 = (double)clock();
	time_t tw0 = time(NULL);
	int status = 0;
	while(s->time < s->tend*(1 - 1e-12) && s->step < s->maxsteps){
		double cap = tnext - s->time;
		if(cap <= 0) cap = s->tend - s->time;
		double dt = g3_step(s, cap);
		if(dt < 0){ status = 5; fprintf(stderr, "[3D] step %ld failed at t=%g\n", s->step, s->time); break; }
		if(s->time >= tnext*(1 - 1e-12)){
			g3_write_snap(s, isnap++, &h0);
			tnext += s->dumpdt > 0 ? s->dumpdt : s->tend;
			if(tnext > s->tend) tnext = s->tend;
		}
		if(s->step % 20 == 0)
			printf("[3D] step %ld t=%.6e dt=%.3e wall=%lds cells_expanded(total)=%ld maxR=%d\n",
					s->step, s->time, dt, (long)(time(NULL) - tw0), s->n_expand, s->max_R);
	}
	(void)wall0;
	if(status == 0 && s->step >= s->maxsteps && s->time < s->tend*(1-1e-12)){
		printf("[3D] stopped at max steps %d, t=%g\n", s->maxsteps, s->time);
		g3_write_snap(s, isnap++, &h0);
	}
	printf("[3D] done: steps=%ld t=%.6e wall=%lds status=%d\n", s->step, s->time, (long)(time(NULL) - tw0), status);
	fflush(stdout);
	return status;
}
