#ifdef USE_CUDA
/* ================================================================
 *  exam_gpu_extract.c
 *
 *  CPU-side routines for GPU force loop:
 *    1. fillParticleSoA      — AoS → SoA conversion
 *    2. extractFaceCSR_2D    — Voronoi linked-list → CSR flat array
 *    3. writeBackResults      — GPU output → AoS particle struct
 *    4. getAccVoro2DBlend_GPU — top-level wrapper replacing getAccVoro2DBlend
 *
 *  Compiled with mpiicx (C, not CUDA). Calls GPU functions from exam_gpu.cu.
 * ================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "exam_gpu.h"
#include "../eunha.h"
#include "../Voro/voro.h"
#include "../Voro/voro_eunha.h"
/* exam.h/exam2d.h omitted: TStruct dependency.
   CellType = HydroTreeLinkedCell (from eunha.h). */

/* ================================================================
 *  ptr_to_soa_index: convert particle pointer → SoA index
 *
 *  Real particles:    [0, nbp)    from VORORK4_TBP
 *  Padding particles: [nbp, nbp+npad) from VORORK4_TBPP
 * ================================================================ */
static inline int ptr_to_soa_index(
    const void *bp,
    const char *real_base, const char *pad_base,
    int nbp, int npad, size_t p_size)
{
    const char *p = (const char *)bp;

    ptrdiff_t off = p - real_base;
    if (off >= 0 && off < (ptrdiff_t)(nbp * p_size))
        return (int)(off / (ptrdiff_t)p_size);

    if (pad_base != NULL) {
        off = p - pad_base;
        if (off >= 0 && off < (ptrdiff_t)(npad * p_size))
            return nbp + (int)(off / (ptrdiff_t)p_size);
    }
    return -1;
}

int gpu_ptr_to_soa_index(
    const void *bp,
    const char *real_base, const char *pad_base,
    int nbp, int npad, size_t p_size)
{
    return ptr_to_soa_index(bp, real_base, pad_base, nbp, npad, p_size);
}

/* ================================================================
 *  ThreadFaceBuf: per-thread growable face buffer for OMP extraction
 * ================================================================ */
typedef struct {
    int count, capacity;
    double *c1x, *c1y, *c2x, *c2y;
    int *neighbor, *kp, *km, *is_ghost, *owner;
} ThreadFaceBuf;

static void tfb_init(ThreadFaceBuf *b, int cap)
{
    b->count = 0;
    b->capacity = cap;
    b->c1x      = (double *)malloc(cap * sizeof(double));
    b->c1y      = (double *)malloc(cap * sizeof(double));
    b->c2x      = (double *)malloc(cap * sizeof(double));
    b->c2y      = (double *)malloc(cap * sizeof(double));
    b->neighbor  = (int *)malloc(cap * sizeof(int));
    b->kp        = (int *)malloc(cap * sizeof(int));
    b->km        = (int *)malloc(cap * sizeof(int));
    b->is_ghost  = (int *)malloc(cap * sizeof(int));
    b->owner     = (int *)malloc(cap * sizeof(int));
}

static void tfb_grow(ThreadFaceBuf *b)
{
    b->capacity *= 2;
    b->c1x      = (double *)realloc(b->c1x,      b->capacity * sizeof(double));
    b->c1y      = (double *)realloc(b->c1y,      b->capacity * sizeof(double));
    b->c2x      = (double *)realloc(b->c2x,      b->capacity * sizeof(double));
    b->c2y      = (double *)realloc(b->c2y,      b->capacity * sizeof(double));
    b->neighbor  = (int *)realloc(b->neighbor,    b->capacity * sizeof(int));
    b->kp        = (int *)realloc(b->kp,          b->capacity * sizeof(int));
    b->km        = (int *)realloc(b->km,          b->capacity * sizeof(int));
    b->is_ghost  = (int *)realloc(b->is_ghost,    b->capacity * sizeof(int));
    b->owner     = (int *)realloc(b->owner,       b->capacity * sizeof(int));
}

static void tfb_free(ThreadFaceBuf *b)
{
    free(b->c1x);  free(b->c1y);
    free(b->c2x);  free(b->c2y);
    free(b->neighbor);
    free(b->kp);   free(b->km);
    free(b->is_ghost);
    free(b->owner);
}

/* Overflow cache statics */
static ThreadFaceBuf *s_overflow_tbufs = NULL;
static int *s_overflow_face_count = NULL;
static int s_overflow_nthreads = 0;
static int s_overflow_total = 0;
static int s_overflow_nbp = 0;

/* ================================================================
 *  Fused tessellation+extraction context
 *
 *  updateDenW2Pressure2DBlend extracts CSR face data during its
 *  gradient face-ring traversal, storing it here.  getAccVoro2DBlend_GPU
 *  then scatters this data into FaceCSR without a second tessellation.
 * ================================================================ */
static struct {
    ThreadFaceBuf *tbufs;
    int *face_count;
    int nthreads;
    int nbp;
    int npad;
    size_t p_size;
    const char *real_base;
    const char *pad_base;
    int total_faces;
    int ready;
} g_fused = {0};

void gpu_fused_begin(int nbp, int npad, size_t p_size, int nthreads,
                     const char *real_base, const char *pad_base)
{
    g_fused.nbp = nbp;
    g_fused.npad = npad;
    g_fused.p_size = p_size;
    g_fused.nthreads = nthreads;
    g_fused.real_base = real_base;
    g_fused.pad_base = pad_base;
    g_fused.total_faces = 0;
    g_fused.ready = 0;

    int init_cap = (20 * nbp / nthreads) + 256;
    g_fused.tbufs = (ThreadFaceBuf *)malloc(nthreads * sizeof(ThreadFaceBuf));
    int t;
    for (t = 0; t < nthreads; t++)
        tfb_init(&g_fused.tbufs[t], init_cap);
    g_fused.face_count = (int *)calloc(nbp, sizeof(int));
}

void gpu_fused_add_face(int tid, int owner_idx,
    double c1x, double c1y, double c2x, double c2y,
    int neighbor_idx, int kp_idx, int km_idx, int is_ghost)
{
    ThreadFaceBuf *mybuf = &g_fused.tbufs[tid];
    if (mybuf->count >= mybuf->capacity)
        tfb_grow(mybuf);
    int f = mybuf->count++;
    mybuf->c1x[f]     = c1x;
    mybuf->c1y[f]     = c1y;
    mybuf->c2x[f]     = c2x;
    mybuf->c2y[f]     = c2y;
    mybuf->neighbor[f] = neighbor_idx;
    mybuf->kp[f]       = kp_idx;
    mybuf->km[f]       = km_idx;
    mybuf->is_ghost[f] = is_ghost;
    mybuf->owner[f]    = owner_idx;
    g_fused.face_count[owner_idx]++;
}

void gpu_fused_end(void)
{
    int t;
    g_fused.total_faces = 0;
    for (t = 0; t < g_fused.nthreads; t++)
        g_fused.total_faces += g_fused.tbufs[t].count;
    g_fused.ready = 1;
}

int gpu_fused_ready(void) { return g_fused.ready; }

static void gpu_fused_scatter(FaceCSR *faces, int max_capacity)
{
    int nbp = g_fused.nbp;
    int t;

    faces->n_particles = nbp;
    faces->n_faces_total = g_fused.total_faces;

    faces->face_offset[0] = 0;
    {
        int ii;
        for (ii = 0; ii < nbp; ii++)
            faces->face_offset[ii + 1] =
                faces->face_offset[ii] + g_fused.face_count[ii];
    }

    if (g_fused.total_faces > max_capacity) {
        s_overflow_tbufs = g_fused.tbufs;
        s_overflow_face_count = g_fused.face_count;
        s_overflow_nthreads = g_fused.nthreads;
        s_overflow_total = g_fused.total_faces;
        s_overflow_nbp = nbp;
        g_fused.tbufs = NULL;
        g_fused.face_count = NULL;
        g_fused.ready = 0;
        return;
    }

    {
        int *write_pos = (int *)malloc(nbp * sizeof(int));
        int ii;
        for (ii = 0; ii < nbp; ii++)
            write_pos[ii] = faces->face_offset[ii];

        for (t = 0; t < g_fused.nthreads; t++) {
            ThreadFaceBuf *b = &g_fused.tbufs[t];
            int f;
            for (f = 0; f < b->count; f++) {
                int owner = b->owner[f];
                int dst = write_pos[owner]++;
                faces->c1x[dst]          = b->c1x[f];
                faces->c1y[dst]          = b->c1y[f];
                faces->c2x[dst]          = b->c2x[f];
                faces->c2y[dst]          = b->c2y[f];
                faces->neighbor_idx[dst] = b->neighbor[f];
                faces->kp_idx[dst]       = b->kp[f];
                faces->km_idx[dst]       = b->km[f];
                faces->is_ghost[dst]     = b->is_ghost[f];
            }
        }
        free(write_pos);
    }

    free(g_fused.face_count);
    g_fused.face_count = NULL;
    for (t = 0; t < g_fused.nthreads; t++)
        tfb_free(&g_fused.tbufs[t]);
    free(g_fused.tbufs);
    g_fused.tbufs = NULL;
    g_fused.ready = 0;
}

/* ================================================================
 *  fillParticleSoA: copy particle data from AoS → SoA (pinned host)
 *
 *  Iterates contiguous real and padding arrays. No cell traversal needed.
 * ================================================================ */
static void fillParticleSoA(
    SimParameters *simpar,
    ParticleSoA *h_parts,
    int has_stress)
{
    size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
    char *real_base = (char *)VORORK4_TBP(simpar);
    char *pad_base  = (char *)VORORK4_TBPP(simpar);
    int nbp  = VORO_NP(simpar);
    int npad = VORO_NPAD(simpar);
    int n_total = nbp + npad;
    int i;

    h_parts->n_total = n_total;

#pragma omp parallel for schedule(static)
    for (i = 0; i < n_total; i++) {
        char *raw;
        if (i < nbp)
            raw = real_base + i * p_size;
        else
            raw = pad_base + (i - nbp) * p_size;

        treevorork4particletype *bp = (treevorork4particletype *)raw;

        h_parts->x[i]        = bp->x;
        h_parts->y[i]        = bp->y;
        h_parts->vx[i]       = bp->vx;
        h_parts->vy[i]       = bp->vy;
        h_parts->mass[i]     = bp->mass;
        h_parts->indx[i]     = (long long)PINDX(bp);
        h_parts->den[i]      = bp->den;
        h_parts->pressure[i] = bp->pressure;
        h_parts->csound[i]   = bp->csound;
        h_parts->volume[i]   = bp->volume;
        h_parts->ie[i]       = bp->ie;
        h_parts->w2[i]       = bp->w2;
        h_parts->w2old[i]    = bp->w2old;
        h_parts->w2ceil[i]   = bp->w2ceil;
        h_parts->avgNeighP[i]= bp->avgNeighboringPressure;

        if (has_stress) {
            treevorostressrk4particletype *sbp =
                (treevorostressrk4particletype *)raw;
            h_parts->gUxx[i]     = sbp->stress.gUxx;
            h_parts->gUxy[i]     = sbp->stress.gUxy;
            h_parts->gUyx[i]     = sbp->stress.gUyx;
            h_parts->gUyy[i]     = sbp->stress.gUyy;
            h_parts->dPdx[i]     = sbp->stress.dPdx;
            h_parts->dPdy[i]     = sbp->stress.dPdy;
            h_parts->dRhodx[i]   = sbp->stress.dRhodx;
            h_parts->dRhody[i]   = sbp->stress.dRhody;
            h_parts->tauxx[i]    = sbp->stress.tauxx;
            h_parts->tauxy[i]    = sbp->stress.tauxy;
            h_parts->tauyy[i]    = sbp->stress.tauyy;
            h_parts->alpha_cd[i] = sbp->stress.alpha_cd;
            h_parts->divv[i]     = sbp->stress.divv;
            h_parts->lap_vx[i]   = sbp->stress.lap_vx;
            h_parts->lap_vy[i]   = sbp->stress.lap_vy;
            h_parts->K[i]        = sbp->stress.K;
            /* LagMFM fields */
            h_parts->E_inv_xx[i] = sbp->stress.E_inv_xx;
            h_parts->E_inv_xy[i] = sbp->stress.E_inv_xy;
            h_parts->E_inv_yx[i] = sbp->stress.E_inv_yx;
            h_parts->E_inv_yy[i] = sbp->stress.E_inv_yy;
            h_parts->h_mfm[i]    = sbp->stress.h_mfm;
        }
    }
}

/* Overflow cache statics (declared above, before fused context) */

/* ================================================================
 *  extractFaceCSR_2D: Voronoi linked-list → CSR flat array
 *
 *  Mirrors the cell-iteration structure of getAccVoro2DBlend.
 *  Single-pass extraction: iterate cells, build Voronoi per particle,
 *  record faces into temporary buffer, then compact to CSR.
 *
 *  Returns total face count.
 * ================================================================ */
int extractFaceCSR_2D(
    SimParameters *simpar,
    FaceCSR *faces,
    ParticleSoA *parts,
    int *particle_map,
    int *particle_map_size,
    int max_faces_capacity,
    postype xmin, postype ymin, postype xmax, postype ymax,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
    treevorork4particletype *(*find2DCellBP)(SimParameters *, int, int, int *),
    void (*mkLinkedList2D)(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)))
{
    size_t p_size  = TVORORK4_DDINFO(simpar)[0].n_size;
    char *real_base = (char *)VORORK4_TBP(simpar);
    char *pad_base  = (char *)VORORK4_TBPP(simpar);
    int nbp  = VORO_NP(simpar);
    int npad = VORO_NPAD(simpar);
    postype boxsize = BOXSIZE(simpar) / NX(simpar) * 5;

    int mx = BASICCELL_MX(simpar);
    int my = BASICCELL_MY(simpar);

    faces->n_particles = nbp;

    /* --- Phase 1: parallel face extraction with per-thread buffers --- */
    int *face_count = (int *)calloc(nbp, sizeof(int));
    int total_faces = 0;
    int nthreads, t;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#else
    nthreads = 1;
#endif

    int init_cap = (20 * nbp / nthreads) + 256;
    ThreadFaceBuf *tbufs =
        (ThreadFaceBuf *)malloc(nthreads * sizeof(ThreadFaceBuf));
    for (t = 0; t < nthreads; t++)
        tfb_init(&tbufs[t], init_cap);

#pragma omp parallel
    {
        int tid;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#else
        tid = 0;
#endif
        ThreadFaceBuf *mybuf = &tbufs[tid];

        int mp = 1000;
        Voro2D_Corner *vorocorner =
            (Voro2D_Corner *)malloc(sizeof(Voro2D_Corner) * mp);
        int iy;

#pragma omp for schedule(dynamic, 4)
        for (iy = 1; iy < my - 1; iy++) {
            int ix;
            for (ix = 1; ix < mx - 1; ix++) {
                int np;
                treevorork4particletype *p =
                    find2DCellBP(simpar, ix, iy, &np);
                int nneigh;
                Voro2D_point *neighbors =
                    find2DNeighboringBP(simpar, ix, iy, &nneigh);
                Voro2D_point *neighwork =
                    (Voro2D_point *)malloc(sizeof(Voro2D_point) * nneigh);
                int i;

                for (i = 0; i < np; i++) {
                    Voro2D_point center;
                    center.x = p[i].x;
                    center.y = p[i].y;
                    center.indx = PINDX(p + i);
                    center.csound = p[i].csound;
                    center.w2 = p[i].w2;

                    treevorork4particletype *ibp_rk4 = p[i].bp;
                    int idx_i = ptr_to_soa_index(ibp_rk4, real_base,
                                                 pad_base, nbp, npad,
                                                 p_size);
                    if (idx_i < 0 || idx_i >= nbp) continue;

                    (void)Voro2D_FindVC(&center, neighbors, neighwork,
                                        nneigh, vorocorner, mp, boxsize);

                    Voro2D_Corner *tmp = vorocorner;
                    do {
                        if (tmp->upperrelated >= 0) {
                            void *jbp_ptr =
                                neighwork[tmp->upperrelated].bp;
                            int idx_j = ptr_to_soa_index(
                                jbp_ptr, real_base, pad_base,
                                nbp, npad, p_size);

                            int ghost =
                                (PINDX((treevorork4particletype *)jbp_ptr)
                                 == MAX_INDEX);

                            Voro2D_Corner *tmp2 = tmp->upperlink;

                            int idx_kp = -1, idx_km = -1;
                            if (tmp2->upperrelated >= 0) {
                                void *kp_ptr =
                                    neighwork[tmp2->upperrelated].bp;
                                idx_kp = ptr_to_soa_index(
                                    kp_ptr, real_base, pad_base,
                                    nbp, npad, p_size);
                            }
                            if (tmp->lowerrelated >= 0) {
                                void *km_ptr =
                                    neighwork[tmp->lowerrelated].bp;
                                idx_km = ptr_to_soa_index(
                                    km_ptr, real_base, pad_base,
                                    nbp, npad, p_size);
                            }

                            if (mybuf->count >= mybuf->capacity)
                                tfb_grow(mybuf);

                            {
                                int f = mybuf->count++;
                                mybuf->c1x[f]     = tmp->x;
                                mybuf->c1y[f]     = tmp->y;
                                mybuf->c2x[f]     = tmp2->x;
                                mybuf->c2y[f]     = tmp2->y;
                                mybuf->neighbor[f] = idx_j;
                                mybuf->kp[f]       = idx_kp;
                                mybuf->km[f]       = idx_km;
                                mybuf->is_ghost[f] = ghost;
                                mybuf->owner[f]    = idx_i;
                            }

                            face_count[idx_i]++;
                        }
                        tmp = tmp->upperlink;
                    } while (tmp != vorocorner);
                }
                free(neighwork);
                free(p);
                free(neighbors);
            }
        }

        free(vorocorner);
    } /* end omp parallel */

    /* Sum total faces from all thread buffers */
    for (t = 0; t < nthreads; t++)
        total_faces += tbufs[t].count;

    /* --- Phase 2: build CSR face_offset from face_count --- */
    faces->n_faces_total = total_faces;
    faces->face_offset[0] = 0;
    {
        int ii;
        for (ii = 0; ii < nbp; ii++)
            faces->face_offset[ii + 1] =
                faces->face_offset[ii] + face_count[ii];
    }

    /* --- Phase 3: scatter thread buffers into CSR order --- */
    if (total_faces > max_faces_capacity) {
        /* Cache thread buffers for scatter after caller reallocates */
        s_overflow_tbufs = tbufs;
        s_overflow_face_count = face_count;
        s_overflow_nthreads = nthreads;
        s_overflow_total = total_faces;
        s_overflow_nbp = nbp;
        if (particle_map_size) *particle_map_size = nbp + npad;
        return total_faces;
    }

    {
        int *write_pos = (int *)malloc(nbp * sizeof(int));
        int ii;
        for (ii = 0; ii < nbp; ii++)
            write_pos[ii] = faces->face_offset[ii];

        for (t = 0; t < nthreads; t++) {
            ThreadFaceBuf *b = &tbufs[t];
            int f;
            for (f = 0; f < b->count; f++) {
                int owner = b->owner[f];
                int dst = write_pos[owner]++;

                faces->c1x[dst]          = b->c1x[f];
                faces->c1y[dst]          = b->c1y[f];
                faces->c2x[dst]          = b->c2x[f];
                faces->c2y[dst]          = b->c2y[f];
                faces->neighbor_idx[dst] = b->neighbor[f];
                faces->kp_idx[dst]       = b->kp[f];
                faces->km_idx[dst]       = b->km[f];
                faces->is_ghost[dst]     = b->is_ghost[f];
            }
        }
        free(write_pos);
    }

    /* Cleanup */
    free(face_count);
    for (t = 0; t < nthreads; t++) tfb_free(&tbufs[t]);
    free(tbufs);

    if (particle_map_size) *particle_map_size = nbp + npad;

    return total_faces;
}

/* ================================================================
 *  extractFaceCSR_2D_scatter_cached: scatter from overflow cache
 * ================================================================ */
static void extractFaceCSR_2D_scatter_cached(FaceCSR *faces)
{
    int nbp = s_overflow_nbp;
    int t;

    faces->n_particles = nbp;
    faces->n_faces_total = s_overflow_total;
    faces->face_offset[0] = 0;
    {
        int ii;
        for (ii = 0; ii < nbp; ii++)
            faces->face_offset[ii + 1] =
                faces->face_offset[ii] + s_overflow_face_count[ii];
    }

    {
        int *write_pos = (int *)malloc(nbp * sizeof(int));
        int ii;
        for (ii = 0; ii < nbp; ii++)
            write_pos[ii] = faces->face_offset[ii];

        for (t = 0; t < s_overflow_nthreads; t++) {
            ThreadFaceBuf *b = &s_overflow_tbufs[t];
            int f;
            for (f = 0; f < b->count; f++) {
                int owner = b->owner[f];
                int dst = write_pos[owner]++;
                faces->c1x[dst]          = b->c1x[f];
                faces->c1y[dst]          = b->c1y[f];
                faces->c2x[dst]          = b->c2x[f];
                faces->c2y[dst]          = b->c2y[f];
                faces->neighbor_idx[dst] = b->neighbor[f];
                faces->kp_idx[dst]       = b->kp[f];
                faces->km_idx[dst]       = b->km[f];
                faces->is_ghost[dst]     = b->is_ghost[f];
            }
        }
        free(write_pos);
    }

    /* Free cached buffers */
    free(s_overflow_face_count);
    for (t = 0; t < s_overflow_nthreads; t++)
        tfb_free(&s_overflow_tbufs[t]);
    free(s_overflow_tbufs);
    s_overflow_tbufs = NULL;
    s_overflow_face_count = NULL;
}

/* ================================================================
 *  writeBackResults: copy GPU output back to particle AoS
 * ================================================================ */
static void writeBackResults_range(
    SimParameters *simpar,
    const ParticleSoA *h_parts,
    int has_stress,
    int i_start, int i_end)
{
    size_t p_size  = TVORORK4_DDINFO(simpar)[0].n_size;
    char *real_base = (char *)VORORK4_TBP(simpar);
    int i;

#pragma omp parallel for schedule(static)
    for (i = i_start; i < i_end; i++) {
        treevorork4particletype *bp =
            (treevorork4particletype *)(real_base + i * p_size);

        bp->ax  = h_parts->ax_out[i];
        bp->ay  = h_parts->ay_out[i];
        bp->die = h_parts->die_out[i];
        bp->dt  = h_parts->dt_out[i];

        if (has_stress) {
            treevorostressrk4particletype *sbp =
                (treevorostressrk4particletype *)bp;
            sbp->stress.vsig_max = h_parts->vsig_max_out[i];
            sbp->stress.v_smooth_x = h_parts->vsmoothx_out[i];
            sbp->stress.v_smooth_y = h_parts->vsmoothy_out[i];
            /* Entropy rate (used by RK4 iff entropy_mode==1; mirrors
             * exam.c:4546-4568 CPU path). dK_out is always written by the
             * GPU kernel — leaving it for the CPU integrator to consume. */
            sbp->stress.dK = h_parts->dK_out[i];
        }
    }
}

static void writeBackResults(
    SimParameters *simpar,
    const ParticleSoA *h_parts,
    int has_stress)
{
    writeBackResults_range(simpar, h_parts, has_stress, 0, VORO_NP(simpar));
}

/* ================================================================
 *  Global GPU context (persistent across RK4 stages)
 * ================================================================ */
static GPUContext g_gpu_ctx = {0};

GPUContext *gpu_get_context(void) { return &g_gpu_ctx; }

/* ================================================================
 *  getAccVoro2DBlend_GPU: top-level wrapper replacing getAccVoro2DBlend
 *
 *  Called from exam2d_vph_rk4_int_blend in place of getAccVoro2DBlend.
 *  All particles processed on GPU (single stream).
 *
 *  Sequence:
 *    1. Build cell linked list (or reuse from fused tessellation)
 *    2. Extract Voronoi faces to CSR (or reuse from GPU tessellation)
 *    3. Fill particle SoA + upload to GPU
 *    4. Launch force kernel → CFL reduction
 *    5. Download results + write back to particle struct
 *    6. MPI_Allreduce(min_dt) → return
 * ================================================================ */

static int gpu_call_count = 0;

double getAccVoro2DBlend_GPU(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    postype OrderOfAccuracy, postype Courant, postype Gamma,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
    treevorork4particletype *(*find2DCellBP)(SimParameters *, int, int, int *),
    void (*mkLinkedList2D)(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)))
{
    double t0_total = MPI_Wtime();
    int nbp  = VORO_NP(simpar);
    int npad;
    int n_total;
    int n_faces;
    int av_mode = GAS_AVMODE(simpar);
    int has_stress = (av_mode >= 1);
    int fused = gpu_fused_ready();

    postype Lx = SIMBOX(simpar).x.max;
    postype cellsize = BASICCELL_CELLWIDTH(simpar);

    if (fused) {
        /* Fused path: cells already freed by updateDenW2Pressure2DBlend,
           padding still alive from that call. */
        npad = VORO_NPAD(simpar);
    } else {
        /* Fallback: build cell linked list (identical to CPU path).
         * GPU-tess path kept TBPP alive in updateDenW2Pressure2DBlend; free it
         * here before paddingTreeVorork4Particles overwrites the pointer. */
        if (VORORK4_TBPP(simpar) != NULL) {
            my_free(VORORK4_TBPP(simpar));
            VORORK4_TBPP(simpar) = NULL;
        }
        int mx, my;
        BASICCELL_MX(simpar) = mx = (int)ceil((xmax - xmin) / cellsize);
        BASICCELL_MY(simpar) = my = (int)ceil((ymax - ymin) / cellsize);
        VORO_BASICCELL(simpar) =
            (HydroTreeLinkedCell *)my_malloc(sizeof(HydroTreeLinkedCell) * mx * my);
        mkLinkedList2D(simpar, cellsize, xmin, ymin, xmax, ymax,
                       paddingAllTreeParticles);
        npad = VORO_NPAD(simpar);
    }
    n_total = nbp + npad;

    /* --- Initialize dt to large value for all particles --- */
    {
        size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
        char *real_base = (char *)VORORK4_TBP(simpar);
        int i;
        for (i = 0; i < nbp; i++) {
            treevorork4particletype *bp =
                (treevorork4particletype *)(real_base + i * p_size);
            bp->dt = 1.e10f;
        }
    }

    /* --- Lazy GPU context initialization --- */
    if (!g_gpu_ctx.initialized) {
        int max_p = (int)(n_total * 1.2) + 1024;
        int max_f = max_p * 20;
        gpu_init(&g_gpu_ctx, max_p, max_f, MYID(simpar));
    }

    /* Ensure capacity */
    if (n_total > g_gpu_ctx.max_particles ||
        nbp * 20 > g_gpu_ctx.max_faces) {
        gpu_free(&g_gpu_ctx);
        int max_p = (int)(n_total * 1.5) + 1024;
        int max_f = max_p * 20;
        gpu_init(&g_gpu_ctx, max_p, max_f, MYID(simpar));
    }

    double t1_extract = MPI_Wtime();

    int tess_on_device = g_gpu_ctx.tess_faces_on_device;

    if (tess_on_device) {
        /* GPU tessellation already produced faces + SoA on device.
           No extraction, no SoA fill, no upload needed. */
        n_faces = g_gpu_ctx.d_faces.n_faces_total;
        g_gpu_ctx.tess_faces_on_device = 0;  /* consume the flag */
    } else if (fused) {
        /* --- Fused path: scatter pre-extracted faces into CSR --- */
        gpu_fused_scatter(&g_gpu_ctx.h_faces, g_gpu_ctx.max_faces);
        n_faces = g_gpu_ctx.h_faces.n_faces_total;

        if (n_faces > g_gpu_ctx.max_faces) {
            gpu_free(&g_gpu_ctx);
            int max_p = (int)(n_total * 1.5) + 1024;
            int max_f = (int)(n_faces * 1.5) + 1024;
            gpu_init(&g_gpu_ctx, max_p, max_f, MYID(simpar));
            extractFaceCSR_2D_scatter_cached(&g_gpu_ctx.h_faces);
        }
    } else {
        /* --- Fallback: full extraction with tessellation --- */
        n_faces = extractFaceCSR_2D(
            simpar,
            &g_gpu_ctx.h_faces,
            &g_gpu_ctx.h_parts,
            NULL, NULL,
            g_gpu_ctx.max_faces,
            xmin, ymin, xmax, ymax,
            paddingAllTreeParticles,
            find2DNeighboringBP,
            find2DCellBP,
            mkLinkedList2D);

        if (n_faces > g_gpu_ctx.max_faces) {
            gpu_free(&g_gpu_ctx);
            int max_p = (int)(n_total * 1.5) + 1024;
            int max_f = (int)(n_faces * 1.5) + 1024;
            gpu_init(&g_gpu_ctx, max_p, max_f, MYID(simpar));
            extractFaceCSR_2D_scatter_cached(&g_gpu_ctx.h_faces);
        }
    }

    double t2_soa = MPI_Wtime();
    double t3_upload;

    if (!tess_on_device) {
        /* --- Fill particle SoA --- */
        fillParticleSoA(simpar, &g_gpu_ctx.h_parts, has_stress);

        t3_upload = MPI_Wtime();
        /* --- Upload to GPU --- */
        gpu_upload_particles(&g_gpu_ctx, n_total, has_stress);
        gpu_upload_faces(&g_gpu_ctx, nbp, n_faces);
    } else {
        t3_upload = MPI_Wtime();
    }

    double t4_kernel = MPI_Wtime();
    /* --- Launch force kernel --- */
    GPUPhysicsParams params;
    params.av_mode     = av_mode;
    params.use_muscl   = GAS_USEMUSCL(simpar);
    params.OoA         = OrderOfAccuracy;
    params.Courant     = Courant;
    params.Gamma       = Gamma;
    params.alphavis    = GAS_AlphaVis(simpar);
    params.betavis     = GAS_BetaVis(simpar);
    params.etavis      = GAS_ETAVIS(simpar) * Lx / NX(simpar);
    params.epsvis      = GAS_EPSVIS(simpar);
    params.nu_phys     = GAS_VISCOSITY(simpar);
    params.prandtl     = GAS_PRANDTL(simpar);
    params.cd_amax     = GAS_CDAMAX(simpar);
    params.blend_theta = GAS_BLENDTHETA(simpar);
    params.dtold       = GAS_dtold(simpar);
    params.hyperv_alpha = GAS_HYPERVALPHA(simpar);
    params.hyperv_force_cap = GAS_HYPERVFORCECAP(simpar);

    /* --- Launch GPU kernel (all particles) --- */
    double Dtime_gpu = gpu_launch_force_kernel(&g_gpu_ctx, nbp, &params);

    double t5_download = MPI_Wtime();
    gpu_download_results(&g_gpu_ctx, nbp);

    double t6_writeback = MPI_Wtime();
    writeBackResults_range(simpar, &g_gpu_ctx.h_parts, has_stress, 0, nbp);

    double Dtime_local = Dtime_gpu;

    /* --- MPI reduction for global minimum dt --- */
    double Dtime;
    MPI_Allreduce(&Dtime_local, &Dtime, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM(simpar));

    double t7_done = MPI_Wtime();

    /* Free: fused path only has padding; fallback has cells + padding */
    if (!fused)
        my_free(VORO_BASICCELL(simpar));
    my_free(VORORK4_TBPP(simpar));
    VORORK4_TBPP(simpar) = NULL;

    gpu_call_count++;
    if (MYID(simpar) == 0 && gpu_call_count % 4 == 0) {
        printf("[GPU step %d] np=%d faces=%d %s | "
               "extract=%.1fms SoA=%.1fms upload=%.1fms "
               "kernel=%.1fms download=%.1fms writeback=%.1fms | "
               "total=%.1fms dt=%.3e\n",
               gpu_call_count / 4, nbp, n_faces,
               fused ? "(fused)" : "(full)",
               (t2_soa - t1_extract) * 1e3,
               (t3_upload - t2_soa) * 1e3,
               (t4_kernel - t3_upload) * 1e3,
               (t5_download - t4_kernel) * 1e3,
               (t6_writeback - t5_download) * 1e3,
               (t7_done - t6_writeback) * 1e3,
               (t7_done - t0_total) * 1e3,
               Dtime_gpu);
        fflush(stdout);
    }

    return Dtime;
}

/* ================================================================
 *  GPU_VALIDATE: CPU reference force on SAME CSR data as GPU kernel
 *
 *  Strategy: run getAccVoro2DBlend_GPU normally, then run an identical
 *  force loop on the CPU using the same host-side FaceCSR + ParticleSoA.
 *  This validates the GPU kernel arithmetic, not the Voronoi extraction.
 * ================================================================ */
#define N_VALIDATE_STEPS 4

/* --- CPU reference helpers (exact copies of GPU device functions) --- */

static inline void cpu_get2dUpqradRk4(
    double xi, double yi, double vxi, double vyi,
    float w2i, float w2oldi, float csi,
    double xj, double yj, double vxj, double vyj,
    float w2j, float w2oldj, float csj,
    double dtold,
    double *uradx, double *urady)
{
    double upqx = vxj - vxi, upqy = vyj - vyi;
    double erx = xj - xi, ery = yj - yi;
    double dpq2 = erx * erx + ery * ery;
    double dpq = sqrt(dpq2);
    double dpq_inv = 1.0 / dpq;
    erx *= dpq_inv;  ery *= dpq_inv;
    double er_dot_upq = erx * upqx + ery * upqy;

    double wp2 = (double)w2i, wq2 = (double)w2j;
    double wpold2 = (double)w2oldi, wqold2 = (double)w2oldj;
    double fact1 = 0.5 * (1.0 + (wp2 - wq2) / dpq2);

    double dwpdt = (sqrt(wp2) - sqrt(wpold2)) / dtold;
    double dwqdt = (sqrt(wq2) - sqrt(wqold2)) / dtold;
    double vpw = dwpdt > 0 ? fmin((double)csi, dwpdt) : fmax(-(double)csi, dwpdt);
    double vqw = dwqdt > 0 ? fmin((double)csj, dwqdt) : fmax(-(double)csj, dwqdt);

    double fact2 = (sqrt(wp2) * vpw - sqrt(wq2) * vqw) * dpq_inv;
    fact2 -= (wp2 - wq2) / dpq2 * er_dot_upq;
    *uradx = fact1 * upqx + fact2 * erx;
    *urady = fact1 * upqy + fact2 * ery;
}

static inline void cpu_hllc_face_2d(
    double rhoL, double pL, double vnL, double cL,
    double rhoR, double pR, double vnR, double cR,
    double Gamma, double *pstar, double *vnstar)
{
    double ZL = rhoL * cL, ZR = rhoR * cR;
    double GP1 = Gamma + 1.0;
    double p_pvrs = (ZR * pL + ZL * pR + ZL * ZR * (vnL - vnR)) / (ZL + ZR);
    if (p_pvrs < 0) p_pvrs = 0;
    double qL = 1.0, qR = 1.0;
    if (p_pvrs > pL) qL = sqrt(1.0 + GP1 / (2.0 * Gamma) * (p_pvrs / pL - 1.0));
    if (p_pvrs > pR) qR = sqrt(1.0 + GP1 / (2.0 * Gamma) * (p_pvrs / pR - 1.0));
    double WL = rhoL * cL * qL, WR = rhoR * cR * qR;
    double Ws = WL + WR;
    *pstar  = (WR * pL + WL * pR + WL * WR * (vnL - vnR)) / Ws;
    *vnstar = (WL * vnL + WR * vnR + pL - pR) / Ws;
}

static inline void cpu_hllc_face_2d_rest_frame(
    double rhoL, double pL, double vnL_lab, double cL,
    double rhoR, double pR, double vnR_lab, double cR,
    double wn, double Gamma, double *pstar, double *vnstar_lab)
{
    double pst, vnst;
    cpu_hllc_face_2d(rhoL, pL, vnL_lab - wn, cL,
                     rhoR, pR, vnR_lab - wn, cR, Gamma, &pst, &vnst);
    *pstar = pst;
    *vnstar_lab = vnst + wn;
}

/* --- CPU reference force loop on CSR data (mirrors GPU kernel exactly) --- */
static void cpu_reference_force_csr(
    const FaceCSR *faces, const ParticleSoA *parts,
    const GPUPhysicsParams *P, int n_particles,
    double *ref_ax, double *ref_ay,
    float *ref_die, float *ref_dt, double *ref_vsig)
{
    int i;
    for (i = 0; i < n_particles; i++) {
        double ibp_x = parts->x[i], ibp_y = parts->y[i];
        double ibp_vx = parts->vx[i], ibp_vy = parts->vy[i];
        double ibp_mass = parts->mass[i];
        double ibp_den = (double)parts->den[i];
        double ibp_pressure = (double)parts->pressure[i];
        double ibp_csound = (double)parts->csound[i];
        double ibp_volume = (double)parts->volume[i];
        float ibp_w2 = parts->w2[i], ibp_w2old = parts->w2old[i];

        double fx = 0, fy = 0, die = 0;
        float my_dt = 1.0e10f;
        double my_vsig_max = 0;
        int f_begin = faces->face_offset[i];
        int f_end = faces->face_offset[i + 1];
        int f;

        for (f = f_begin; f < f_end; f++) {
            int j = faces->neighbor_idx[f];
            if (j < 0) continue;
            int jbp_is_ghost = faces->is_ghost[f];

            double jbp_x = parts->x[j], jbp_y = parts->y[j];
            double jbp_vx = parts->vx[j], jbp_vy = parts->vy[j];
            double jbp_den = (double)parts->den[j];
            double jbp_pressure = (double)parts->pressure[j];
            double jbp_csound = (double)parts->csound[j];
            float jbp_w2 = parts->w2[j], jbp_w2old = parts->w2old[j];

            double line_x = faces->c2x[f] - faces->c1x[f];
            double line_y = faces->c2y[f] - faces->c1y[f];
            double dSx = line_y, dSy = -line_x;
            double facearea = sqrt(dSx * dSx + dSy * dSy);

            double drx = jbp_x - ibp_x, dry = jbp_y - ibp_y;
            double dramp = sqrt(drx * drx + dry * dry);
            double dramp_inv = 1.0 / (dramp + 1.0e-30);
            double erx = drx * dramp_inv, ery = dry * dramp_inv;

            double pi_total;
            double tau_dot_dS_x = 0, tau_dot_dS_y = 0;

            if (P->av_mode == 0 && P->nu_phys <= 0) {
                int kp = faces->kp_idx[f], km = faces->km_idx[f];
                double p_kp = (kp >= 0) ? (double)parts->pressure[kp] : ibp_pressure;
                double p_km = (km >= 0) ? (double)parts->pressure[km] : jbp_pressure;
                double w = P->OoA / 3.0;
                pi_total = (0.5 - w) * (ibp_pressure + jbp_pressure) + w * (p_kp + p_km);

                double uijx = jbp_vx - ibp_vx, uijy = jbp_vy - ibp_vy;
                double rvel = erx * uijx + ery * uijy;
                if (rvel < 0) {
                    double wcomp = sqrt((double)ibp_w2) + sqrt((double)jbp_w2);
                    double scaleFactor = (wcomp == 0 ? P->etavis : wcomp);
                    double drampScale = dramp / scaleFactor;
                    double mu = rvel / (drampScale + P->epsvis / drampScale);
                    double meanden = 0.5 * (ibp_den + jbp_den);
                    double meanCsound = 0.5 * (ibp_csound + jbp_csound);
                    pi_total += (-P->alphavis * meanCsound * mu + P->betavis * mu * mu) * meanden;
                }
            } else if (jbp_is_ghost) {
                int kp = faces->kp_idx[f], km = faces->km_idx[f];
                double p_kp = (kp >= 0) ? (double)parts->pressure[kp] : ibp_pressure;
                double p_km = (km >= 0) ? (double)parts->pressure[km] : jbp_pressure;
                double w = P->OoA / 3.0;
                pi_total = (0.5 - w) * (ibp_pressure + jbp_pressure) + w * (p_kp + p_km);

                double uijx = jbp_vx - ibp_vx, uijy = jbp_vy - ibp_vy;
                double rvel = erx * uijx + ery * uijy;
                if (rvel < 0) {
                    double wcomp = sqrt((double)ibp_w2) + sqrt((double)jbp_w2);
                    double scaleFactor = (wcomp == 0 ? P->etavis : wcomp);
                    double drampScale = dramp / scaleFactor;
                    double mu = rvel / (drampScale + P->epsvis / drampScale);
                    double meanden = 0.5 * (ibp_den + jbp_den);
                    double meanCsound = 0.5 * (ibp_csound + jbp_csound);
                    pi_total += (-P->alphavis * meanCsound * mu + P->betavis * mu * mu) * meanden;
                }
            } else if (P->av_mode == 5) {
                double ds_mag_inv = 1.0 / (facearea + 1.0e-30);
                double nx_hat = dSx * ds_mag_inv, ny_hat = dSy * ds_mag_inv;
                double uradx_t, urady_t;
                cpu_get2dUpqradRk4(ibp_x, ibp_y, ibp_vx, ibp_vy, ibp_w2, ibp_w2old, (float)ibp_csound,
                    jbp_x, jbp_y, jbp_vx, jbp_vy, jbp_w2, jbp_w2old, (float)jbp_csound,
                    P->dtold, &uradx_t, &urady_t);
                double wx = ibp_vx + uradx_t, wy = ibp_vy + urady_t;
                double wn = wx * nx_hat + wy * ny_hat;
                double vnL_lab = ibp_vx * nx_hat + ibp_vy * ny_hat;
                double vnR_lab = jbp_vx * nx_hat + jbp_vy * ny_hat;
                double pL = ibp_pressure, pR = jbp_pressure;
                double rhoL = ibp_den, rhoR = jbp_den;

                if (P->use_muscl) {
                    double xf_rel_i = 0.5 * (faces->c1x[f] + faces->c2x[f]);
                    double yf_rel_i = 0.5 * (faces->c1y[f] + faces->c2y[f]);
                    double dx_ij = jbp_x - ibp_x, dy_ij = jbp_y - ibp_y;
                    double dx_iF = xf_rel_i, dy_iF = yf_rel_i;
                    double dx_jF = xf_rel_i - dx_ij, dy_jF = yf_rel_i - dy_ij;
                    rhoL += parts->dRhodx[i]*dx_iF + parts->dRhody[i]*dy_iF;
                    rhoR += parts->dRhodx[j]*dx_jF + parts->dRhody[j]*dy_jF;
                    pL += parts->dPdx[i]*dx_iF + parts->dPdy[i]*dy_iF;
                    pR += parts->dPdx[j]*dx_jF + parts->dPdy[j]*dy_jF;
                    double dvxL = parts->gUxx[i]*dx_iF + parts->gUxy[i]*dy_iF;
                    double dvyL = parts->gUyx[i]*dx_iF + parts->gUyy[i]*dy_iF;
                    vnL_lab += nx_hat*dvxL + ny_hat*dvyL;
                    double dvxR = parts->gUxx[j]*dx_jF + parts->gUxy[j]*dy_jF;
                    double dvyR = parts->gUyx[j]*dx_jF + parts->gUyy[j]*dy_jF;
                    vnR_lab += nx_hat*dvxR + ny_hat*dvyR;
                }
                if (pL < 1.0e-10) pL = 1.0e-10;
                if (pR < 1.0e-10) pR = 1.0e-10;
                if (rhoL < 1.0e-10) rhoL = 1.0e-10;
                if (rhoR < 1.0e-10) rhoR = 1.0e-10;

                double pst, vnst_lab;
                cpu_hllc_face_2d_rest_frame(rhoL, pL, vnL_lab, ibp_csound,
                    rhoR, pR, vnR_lab, jbp_csound, wn, P->Gamma, &pst, &vnst_lab);
                pi_total = pst;

                if (P->alphavis > 0) {
                    double uijx = jbp_vx - ibp_vx, uijy = jbp_vy - ibp_vy;
                    double rvel = erx * uijx + ery * uijy;
                    if (rvel < 0) {
                        double wcomp = sqrt((double)ibp_w2) + sqrt((double)jbp_w2);
                        double scaleFactor = (wcomp == 0 ? P->etavis : wcomp);
                        double drampScale = dramp / scaleFactor;
                        double mu = rvel / (drampScale + P->epsvis / drampScale);
                        double meanden = 0.5 * (ibp_den + jbp_den);
                        double meanCsound = 0.5 * (ibp_csound + jbp_csound);
                        pi_total += (-P->alphavis*meanCsound*mu + P->betavis*mu*mu)*meanden;
                    }
                }
            } else if (P->av_mode == 3) {
                double ds_mag_inv = 1.0 / (facearea + 1.0e-30);
                double nx_hat = dSx * ds_mag_inv, ny_hat = dSy * ds_mag_inv;
                double vnL = ibp_vx*nx_hat + ibp_vy*ny_hat;
                double vnR = jbp_vx*nx_hat + jbp_vy*ny_hat;
                double pL = ibp_pressure, pR = jbp_pressure;

                if (P->use_muscl) {
                    double xf_rel_i = 0.5 * (faces->c1x[f] + faces->c2x[f]);
                    double yf_rel_i = 0.5 * (faces->c1y[f] + faces->c2y[f]);
                    double dx_ij = jbp_x - ibp_x, dy_ij = jbp_y - ibp_y;
                    double dx_iF = xf_rel_i, dy_iF = yf_rel_i;
                    double dx_jF = xf_rel_i - dx_ij, dy_jF = yf_rel_i - dy_ij;
                    double dp_i = parts->dPdx[i]*dx_iF + parts->dPdy[i]*dy_iF;
                    double dp_j = parts->dPdx[j]*dx_jF + parts->dPdy[j]*dy_jF;
                    pL += dp_i; pR += dp_j;
                    double dvnL0 = (parts->gUxx[i]*nx_hat+parts->gUxy[i]*ny_hat)*dx_iF
                                 + (parts->gUyx[i]*nx_hat+parts->gUyy[i]*ny_hat)*dy_iF;
                    double dvnR0 = (parts->gUxx[j]*nx_hat+parts->gUxy[j]*ny_hat)*dx_jF
                                 + (parts->gUyx[j]*nx_hat+parts->gUyy[j]*ny_hat)*dy_jF;
                    vnL += dvnL0; vnR += dvnR0;
                    double pmin = fmin(ibp_pressure, jbp_pressure);
                    double pmax = fmax(ibp_pressure, jbp_pressure);
                    pL = fmax(pmin, fmin(pmax, pL));
                    pR = fmax(pmin, fmin(pmax, pR));
                    double vnmin = fmin(vnL-dvnL0, vnR-dvnR0);
                    double vnmax = fmax(vnL-dvnL0, vnR-dvnR0);
                    vnL = fmax(vnmin, fmin(vnmax, vnL));
                    vnR = fmax(vnmin, fmin(vnmax, vnR));
                }
                if (pL < 1.0e-10) pL = 1.0e-10;
                if (pR < 1.0e-10) pR = 1.0e-10;

                double pst, vnst;
                cpu_hllc_face_2d(ibp_den, pL, vnL, ibp_csound,
                    jbp_den, pR, vnR, jbp_csound, P->Gamma, &pst, &vnst);
                pi_total = pst;

                if (P->cd_amax > 0) {
                    double dvn_hllc = vnR - vnL;
                    if (dvn_hllc < 0) {
                        double alpha_face = 0.5*(parts->alpha_cd[i]+parts->alpha_cd[j]);
                        double rho_mean = 0.5*(ibp_den+jbp_den);
                        double vsig_cd = ibp_csound+jbp_csound-fmin(0.0, dvn_hllc);
                        pi_total += 0.5*alpha_face*vsig_cd*rho_mean*(-dvn_hllc);
                    }
                }
                if (P->alphavis > 0) {
                    double uijx = jbp_vx-ibp_vx, uijy = jbp_vy-ibp_vy;
                    double rvel = erx*uijx + ery*uijy;
                    if (rvel < 0) {
                        double wcomp = sqrt((double)ibp_w2)+sqrt((double)jbp_w2);
                        double scaleFactor = (wcomp == 0 ? P->etavis : wcomp);
                        double drampScale = dramp/scaleFactor;
                        double mu = rvel/(drampScale+P->epsvis/drampScale);
                        double meanden = 0.5*(ibp_den+jbp_den);
                        double meanCsound = 0.5*(ibp_csound+jbp_csound);
                        pi_total += (-P->alphavis*meanCsound*mu+P->betavis*mu*mu)*meanden;
                    }
                }
            } else {
                /* av_mode==0 (with nu_phys>0), 1, or 2 */
                int kp = faces->kp_idx[f], km = faces->km_idx[f];
                double p_kp = (kp >= 0) ? (double)parts->pressure[kp] : ibp_pressure;
                double p_km = (km >= 0) ? (double)parts->pressure[km] : jbp_pressure;
                double w = P->OoA / 3.0;
                double p_mnm = (0.5-w)*(ibp_pressure+jbp_pressure) + w*(p_kp+p_km);

                if (P->av_mode == 1 || (P->av_mode == 0 && P->nu_phys > 0)) {
                    double txx_face, txy_face, tyy_face;
                    if (P->OoA > 0 && kp >= 0 && km >= 0) {
                        txx_face = (0.5-w)*(parts->tauxx[i]+parts->tauxx[j])
                                 + w*(parts->tauxx[kp]+parts->tauxx[km]);
                        txy_face = (0.5-w)*(parts->tauxy[i]+parts->tauxy[j])
                                 + w*(parts->tauxy[kp]+parts->tauxy[km]);
                        tyy_face = (0.5-w)*(parts->tauyy[i]+parts->tauyy[j])
                                 + w*(parts->tauyy[kp]+parts->tauyy[km]);
                    } else {
                        txx_face = 0.5*(parts->tauxx[i]+parts->tauxx[j]);
                        txy_face = 0.5*(parts->tauxy[i]+parts->tauxy[j]);
                        tyy_face = 0.5*(parts->tauyy[i]+parts->tauyy[j]);
                    }
                    tau_dot_dS_x = -(txx_face*dSx + txy_face*dSy);
                    tau_dot_dS_y = -(txy_face*dSx + tyy_face*dSy);
                }

                if (P->av_mode == 0 || P->av_mode == 1) {
                    pi_total = p_mnm;
                    if (P->alphavis > 0) {
                        double uijx = jbp_vx-ibp_vx, uijy = jbp_vy-ibp_vy;
                        double rvel = erx*uijx + ery*uijy;
                        if (rvel < 0) {
                            double wcomp = sqrt((double)ibp_w2)+sqrt((double)jbp_w2);
                            double scaleFactor = (wcomp == 0 ? P->etavis : wcomp);
                            double drampScale = dramp/scaleFactor;
                            double mu = rvel/(drampScale+P->epsvis/drampScale);
                            double meanden = 0.5*(ibp_den+jbp_den);
                            double meanCsound = 0.5*(ibp_csound+jbp_csound);
                            pi_total += (-P->alphavis*meanCsound*mu+P->betavis*mu*mu)*meanden;
                        }
                    }
                } else {
                    /* av_mode==2: two-tier blend */
                    double alpha_max_pq = fmax(parts->alpha_cd[i], parts->alpha_cd[j]);
                    double theta_pq = 0;
                    if (parts->divv[i] < -1.0e-10 || parts->divv[j] < -1.0e-10)
                        theta_pq = P->blend_theta;
                    double f_pq = fmin(1.0, alpha_max_pq/(P->cd_amax+1.0e-30)+theta_pq);

                    double p_hllc = p_mnm;
                    double Pi_cd10 = 0;

                    if (f_pq > 1.0e-6) {
                        double ds_mag_inv = 1.0/(facearea+1.0e-30);
                        double nx_hat = dSx*ds_mag_inv, ny_hat = dSy*ds_mag_inv;
                        double vnL = ibp_vx*nx_hat+ibp_vy*ny_hat;
                        double vnR = jbp_vx*nx_hat+jbp_vy*ny_hat;
                        double rho_mean_det = 0.5*(ibp_den+jbp_den);
                        double p_mean_det = 0.5*(ibp_pressure+jbp_pressure);
                        double drho_rel = fabs(ibp_den-jbp_den)/(rho_mean_det+1.0e-30);
                        double dp_rel = fabs(ibp_pressure-jbp_pressure)/(p_mean_det+1.0e-30);
                        int is_pure_contact = (drho_rel > 0.1) && (dp_rel < 0.02);

                        double pL = ibp_pressure, pR = jbp_pressure;
                        double rhoL = ibp_den, rhoR = jbp_den;
                        double cL = ibp_csound, cR = jbp_csound;

                        if (P->use_muscl) {
                            double xf_rel_i = 0.5*(faces->c1x[f]+faces->c2x[f]);
                            double yf_rel_i = 0.5*(faces->c1y[f]+faces->c2y[f]);
                            double dx_ij = jbp_x-ibp_x, dy_ij = jbp_y-ibp_y;
                            double dx_iF = xf_rel_i, dy_iF = yf_rel_i;
                            double dx_jF = xf_rel_i-dx_ij, dy_jF = yf_rel_i-dy_ij;
                            double dp_i = parts->dPdx[i]*dx_iF+parts->dPdy[i]*dy_iF;
                            double dp_j = parts->dPdx[j]*dx_jF+parts->dPdy[j]*dy_jF;
                            pL += dp_i; pR += dp_j;
                            double dvnL0 = 0, dvnR0 = 0;
                            if (!is_pure_contact) {
                                dvnL0 = (parts->gUxx[i]*nx_hat+parts->gUxy[i]*ny_hat)*dx_iF
                                       + (parts->gUyx[i]*nx_hat+parts->gUyy[i]*ny_hat)*dy_iF;
                                dvnR0 = (parts->gUxx[j]*nx_hat+parts->gUxy[j]*ny_hat)*dx_jF
                                       + (parts->gUyx[j]*nx_hat+parts->gUyy[j]*ny_hat)*dy_jF;
                                vnL += dvnL0; vnR += dvnR0;
                            }
                            double pmin = fmin(ibp_pressure, jbp_pressure);
                            double pmax_v = fmax(ibp_pressure, jbp_pressure);
                            pL = fmax(pmin, fmin(pmax_v, pL));
                            pR = fmax(pmin, fmin(pmax_v, pR));
                            double vnmin = fmin(vnL-dvnL0, vnR-dvnR0);
                            double vnmax = fmax(vnL-dvnL0, vnR-dvnR0);
                            vnL = fmax(vnmin, fmin(vnmax, vnL));
                            vnR = fmax(vnmin, fmin(vnmax, vnR));
                        }
                        pL = P->use_muscl ? pL : ibp_pressure;
                        pR = P->use_muscl ? pR : jbp_pressure;

                        double wx = 0.5*(ibp_vx+jbp_vx), wy = 0.5*(ibp_vy+jbp_vy);
                        double wn = wx*nx_hat+wy*ny_hat;
                        double pst, vnst_lab;
                        cpu_hllc_face_2d_rest_frame(rhoL, pL, vnL, cL,
                            rhoR, pR, vnR, cR, wn, P->Gamma, &pst, &vnst_lab);
                        p_hllc = pst;

                        if (f_pq <= 0.5) {
                            double dvn = vnR - vnL;
                            if (dvn < 0) {
                                double alpha_face = 0.5*(parts->alpha_cd[i]+parts->alpha_cd[j]);
                                double rho_mean = 0.5*(rhoL+rhoR);
                                double vsig_cd = cL+cR-fmin(0.0, dvn);
                                Pi_cd10 = 0.5*alpha_face*vsig_cd*rho_mean*(-dvn);
                            }
                        }
                    }
                    pi_total = (1.0-f_pq)*p_mnm + f_pq*(p_hllc+Pi_cd10);
                    tau_dot_dS_x *= (1.0-f_pq);
                    tau_dot_dS_y *= (1.0-f_pq);
                }
            }

            /* Face velocity */
            double uradx, urady;
            cpu_get2dUpqradRk4(ibp_x, ibp_y, ibp_vx, ibp_vy, ibp_w2, ibp_w2old, (float)ibp_csound,
                jbp_x, jbp_y, jbp_vx, jbp_vy, jbp_w2, jbp_w2old, (float)jbp_csound,
                P->dtold, &uradx, &urady);

            die += -pi_total * (uradx * dSx + urady * dSy);
            die += tau_dot_dS_x * uradx + tau_dot_dS_y * urady;

            if (P->nu_phys > 0 && P->prandtl > 0 && !jbp_is_ghost) {
                double chi = P->nu_phys / P->prandtl;
                double Ti = ibp_pressure / ibp_den;
                double Tj = jbp_pressure / jbp_den;
                double rho_face = 0.5 * (ibp_den + jbp_den);
                die += chi * rho_face * (Tj - Ti) / dramp * facearea;
            }

            fx += -pi_total * dSx + tau_dot_dS_x;
            fy += -pi_total * dSy + tau_dot_dS_y;

            double dvx = jbp_vx - ibp_vx, dvy = jbp_vy - ibp_vy;
            double VdotR = dvx * erx + dvy * ery;
            double vsig = jbp_csound + ibp_csound - fmin(0.0, VdotR);
            double heff = 0.25 * sqrt(ibp_volume);
            double dramp_cfl = fmax(dramp, heff);
            float dt_face = (float)(2.0 * P->Courant * dramp_cfl / vsig);
            if (dt_face < my_dt) my_dt = dt_face;

            if (P->av_mode >= 1) {
                double vsig_cd = jbp_is_ghost ? (ibp_csound+jbp_csound) : vsig;
                if (vsig_cd > my_vsig_max) my_vsig_max = vsig_cd;
            }
        } /* end face loop */

        if (P->av_mode >= 1 || P->nu_phys > 0) {
            double h_i = sqrt(ibp_volume);
            double nu_cd_i = (P->av_mode >= 1) ? parts->alpha_cd[i]*h_i*ibp_csound : 0;
            double chi = (P->nu_phys > 0 && P->prandtl > 0) ? P->nu_phys/P->prandtl : 0;
            double nu_eff = fmax(fmax(P->nu_phys, nu_cd_i), chi);
            if (nu_eff > 0) {
                float dt_visc = (float)(0.5 * h_i * h_i / nu_eff);
                if (dt_visc < my_dt) my_dt = dt_visc;
            }
        }

        ref_ax[i] = fx / ibp_mass;
        ref_ay[i] = fy / ibp_mass;
        ref_die[i] = (float)die;
        ref_dt[i] = my_dt;
        ref_vsig[i] = my_vsig_max;
    }
}

/* --- Validate function --- */
double getAccVoro2DBlend_GPU_validate(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    postype OrderOfAccuracy, postype Courant, postype Gamma,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
    treevorork4particletype *(*find2DCellBP)(SimParameters *, int, int, int *),
    void (*mkLinkedList2D)(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)))
{
    static int validate_count = 0;

    if (validate_count >= N_VALIDATE_STEPS) {
        return getAccVoro2DBlend_GPU(simpar, xmin, ymin, xmax, ymax,
            OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
            find2DNeighboringBP, find2DCellBP, mkLinkedList2D);
    }
    validate_count++;

    int nbp = VORO_NP(simpar);
    postype Lx = SIMBOX(simpar).x.max;

    /* --- Step 1: Run GPU path normally --- */
    double Dtime_gpu = getAccVoro2DBlend_GPU(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
        find2DNeighboringBP, find2DCellBP, mkLinkedList2D);

    /* After GPU wrapper returns, g_gpu_ctx.h_faces and h_parts
       still hold the CSR face data and SoA particle data (pinned memory).
       h_parts->ax_out etc. hold the GPU kernel output. */

    /* --- Step 2: Run CPU reference on the SAME CSR data --- */
    GPUPhysicsParams params;
    params.av_mode     = GAS_AVMODE(simpar);
    params.use_muscl   = GAS_USEMUSCL(simpar);
    params.OoA         = OrderOfAccuracy;
    params.Courant     = Courant;
    params.Gamma       = Gamma;
    params.alphavis    = GAS_AlphaVis(simpar);
    params.betavis     = GAS_BetaVis(simpar);
    params.etavis      = GAS_ETAVIS(simpar) * Lx / NX(simpar);
    params.epsvis      = GAS_EPSVIS(simpar);
    params.nu_phys     = GAS_VISCOSITY(simpar);
    params.prandtl     = GAS_PRANDTL(simpar);
    params.cd_amax     = GAS_CDAMAX(simpar);
    params.blend_theta = GAS_BLENDTHETA(simpar);
    params.dtold       = GAS_dtold(simpar);

    double *ref_ax  = (double *)malloc(nbp * sizeof(double));
    double *ref_ay  = (double *)malloc(nbp * sizeof(double));
    float  *ref_die = (float  *)malloc(nbp * sizeof(float));
    float  *ref_dt  = (float  *)malloc(nbp * sizeof(float));
    double *ref_vsig = (double *)malloc(nbp * sizeof(double));

    cpu_reference_force_csr(&g_gpu_ctx.h_faces, &g_gpu_ctx.h_parts,
                            &params, nbp,
                            ref_ax, ref_ay, ref_die, ref_dt, ref_vsig);

    /* --- Step 3: Compare GPU output vs CPU reference --- */
    double max_rel_ax = 0, max_rel_ay = 0, max_rel_die = 0, max_rel_dt = 0;
    double max_abs_ax = 0, max_abs_die = 0;
    int worst_rel_ax = 0, worst_abs_ax = 0;
    int n_mismatch_ax = 0, n_mismatch_die = 0;
    int i;

    for (i = 0; i < nbp; i++) {
        double gpu_ax = g_gpu_ctx.h_parts.ax_out[i];
        double gpu_ay = g_gpu_ctx.h_parts.ay_out[i];
        float  gpu_die = g_gpu_ctx.h_parts.die_out[i];
        float  gpu_dt = g_gpu_ctx.h_parts.dt_out[i];

        double scale, err;

        /* ax absolute error */
        double abs_ax = fabs(gpu_ax - ref_ax[i]);
        if (abs_ax > max_abs_ax) { max_abs_ax = abs_ax; worst_abs_ax = i; }

        /* ax relative error */
        scale = fmax(fabs(gpu_ax), fabs(ref_ax[i])) + 1.0e-30;
        err = fabs(gpu_ax - ref_ax[i]) / scale;
        if (err > max_rel_ax) { max_rel_ax = err; worst_rel_ax = i; }
        /* Count as mismatch only if BOTH relative and absolute exceed threshold */
        if (err > 1.0e-6 && abs_ax > 1.0e-11) n_mismatch_ax++;

        /* ay relative error */
        scale = fmax(fabs(gpu_ay), fabs(ref_ay[i])) + 1.0e-30;
        err = fabs(gpu_ay - ref_ay[i]) / scale;
        if (err > max_rel_ay) max_rel_ay = err;

        /* die relative error */
        scale = fmax(fabs(gpu_die), fabs(ref_die[i])) + 1.0e-30;
        err = fabs(gpu_die - ref_die[i]) / scale;
        if (err > max_rel_die) max_rel_die = err;
        if (err > 1.0e-6 && fabs(gpu_die - ref_die[i]) > 1.0e-8) n_mismatch_die++;

        /* dt relative error */
        scale = fmax(fabs(gpu_dt), fabs(ref_dt[i])) + 1.0e-30;
        err = fabs(gpu_dt - ref_dt[i]) / scale;
        if (err > max_rel_dt) max_rel_dt = err;
    }

    fprintf(stderr,
        "[GPU_VALIDATE %d] nbp=%d  (GPU kernel vs CPU reference on SAME CSR data)\n"
        "  ax:  max_rel=%g (i=%d) max_abs=%g  mismatches(>1e-6)=%d/%d\n"
        "  ay:  max_rel=%g\n"
        "  die: max_rel=%g  mismatches=%d\n"
        "  dt:  max_rel=%g\n",
        validate_count, nbp,
        max_rel_ax, worst_rel_ax, max_abs_ax, n_mismatch_ax, nbp,
        max_rel_ay,
        max_rel_die, n_mismatch_die,
        max_rel_dt);

    if (n_mismatch_ax > 0) {
        int wi = worst_rel_ax;
        fprintf(stderr,
            "  worst ax i=%d: gpu=%g ref=%g  diff=%g\n"
            "    gpu_ay=%g ref_ay=%g  gpu_die=%g ref_die=%g\n",
            wi, g_gpu_ctx.h_parts.ax_out[wi], ref_ax[wi],
            g_gpu_ctx.h_parts.ax_out[wi] - ref_ax[wi],
            g_gpu_ctx.h_parts.ay_out[wi], ref_ay[wi],
            g_gpu_ctx.h_parts.die_out[wi], ref_die[wi]);
    }

    /* PASS if no real mismatches (abs threshold filters out near-zero noise) */
    if (n_mismatch_ax == 0 && n_mismatch_die == 0)
        fprintf(stderr, "  >>> PASS: GPU kernel matches CPU reference\n");
    else
        fprintf(stderr, "  >>> FAIL: %d ax + %d die real mismatches\n",
                n_mismatch_ax, n_mismatch_die);

    free(ref_ax); free(ref_ay); free(ref_die); free(ref_dt); free(ref_vsig);

    /* GPU results already in particles from writeBackResults — keep them */
    return Dtime_gpu;
}

/* ================================================================
 *  gpu_build_cell_csr: traverse linked list → flat CellCSR
 *
 *  Called after mkLinkedList2D. Walks cell linked lists to build
 *  a CSR-format cell array for GPU tessellation.
 *  Includes BOTH real and padding particles (padding follow linked
 *  list via the `next` pointer, same as real particles).
 * ================================================================ */
void gpu_build_cell_csr(GPUContext *ctx, void *simpar_opaque,
                        int mx, int my, size_t p_size)
{
    SimParameters *simpar = (SimParameters *)simpar_opaque;
    HydroTreeLinkedCell *cells = VORO_BASICCELL(simpar);
    char *real_base = (char *)VORORK4_TBP(simpar);
    char *pad_base  = (char *)VORORK4_TBPP(simpar);
    int nbp  = VORO_NP(simpar);
    int npad = VORO_NPAD(simpar);
    int n_cells = mx * my;

    /* Ensure tessellation buffers are allocated */
    if (ctx->max_cells < n_cells || ctx->max_cell_entries < nbp + npad) {
        if (ctx->max_cells > 0)
            gpu_free_tess_buffers(ctx);
        int est_entries = (nbp + npad) * 2;  /* generous estimate */
        gpu_alloc_tess_buffers(ctx, ctx->max_particles, n_cells, est_entries);
    }

    /* First pass: count entries per cell */
    int running_count = 0;
    int cell;
    for (cell = 0; cell < n_cells; cell++) {
        ctx->h_cells.cell_offset[cell] = running_count;
        linkedlisttype *tmp = cells[cell].link;
        while (tmp) {
            running_count++;
            tmp = tmp->next;
        }
    }
    ctx->h_cells.cell_offset[n_cells] = running_count;

    /* Reallocate if needed */
    if (running_count > ctx->max_cell_entries) {
        fprintf(stderr, "[GPU] Reallocating cell entries: %d -> %d\n",
                ctx->max_cell_entries, running_count * 2);
        gpu_free_tess_buffers(ctx);
        gpu_alloc_tess_buffers(ctx, ctx->max_particles, n_cells, running_count * 2);
        /* Redo count pass */
        running_count = 0;
        for (cell = 0; cell < n_cells; cell++) {
            ctx->h_cells.cell_offset[cell] = running_count;
            linkedlisttype *tmp = cells[cell].link;
            while (tmp) {
                running_count++;
                tmp = tmp->next;
            }
        }
        ctx->h_cells.cell_offset[n_cells] = running_count;
    }

    /* Second pass: fill particle indices */
    running_count = 0;
    for (cell = 0; cell < n_cells; cell++) {
        linkedlisttype *tmp = cells[cell].link;
        while (tmp) {
            int idx = ptr_to_soa_index((void *)tmp, real_base, pad_base,
                                       nbp, npad, p_size);
            ctx->h_cells.particle_index[running_count] = idx;
            running_count++;
            tmp = tmp->next;
        }
    }

    ctx->h_cells.n_cells = n_cells;
    ctx->h_cells.n_entries = running_count;
}

/* ================================================================
 *  gpu_writeback_tess_results: pinned host SoA → AoS particle structs
 *
 *  After GPU tessellation + download, write results back to CPU
 *  particle arrays (real particles only, [0, nbp)).
 * ================================================================ */
void gpu_writeback_tess_results(void *simpar_opaque, GPUContext *ctx,
                                int nbp, int has_stress)
{
    SimParameters *simpar = (SimParameters *)simpar_opaque;
    size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
    char *real_base = (char *)VORORK4_TBP(simpar);
    ParticleSoA *h = &ctx->h_parts;
    int i;

#pragma omp parallel for schedule(static)
    for (i = 0; i < nbp; i++) {
        char *raw = real_base + i * p_size;
        treevorork4particletype *bp = (treevorork4particletype *)raw;

        bp->den                    = h->den[i];
        bp->pressure               = h->pressure[i];
        bp->csound                 = h->csound[i];
        bp->volume                 = h->volume[i];
        bp->w2ceil                 = h->w2ceil[i];
        bp->avgNeighboringPressure = h->avgNeighP[i];

        if (has_stress) {
            treevorostressrk4particletype *sbp =
                (treevorostressrk4particletype *)raw;
            sbp->stress.gUxx   = h->gUxx[i];
            sbp->stress.gUxy   = h->gUxy[i];
            sbp->stress.gUyx   = h->gUyx[i];
            sbp->stress.gUyy   = h->gUyy[i];
            sbp->stress.divv   = h->gUxx[i] + h->gUyy[i]; /* recalc from limited */
            sbp->stress.dPdx   = h->dPdx[i];
            sbp->stress.dPdy   = h->dPdy[i];
            sbp->stress.dRhodx = h->dRhodx[i];
            sbp->stress.dRhody = h->dRhody[i];
            sbp->stress.tauxx  = h->tauxx[i];
            sbp->stress.tauxy  = h->tauxy[i];
            sbp->stress.tauyy  = h->tauyy[i];
            sbp->stress.divv   = h->divv[i];
            sbp->stress.lap_vx = h->lap_vx[i];
            sbp->stress.lap_vy = h->lap_vy[i];
            /* K may have been clamped to K_floor by pressure_stress_kernel
             * when entropy_mode==1. When entropy_mode==0, GPU does not touch
             * K so this is a no-op (h->K[i] equals what we uploaded). */
            sbp->stress.K = h->K[i];
        }
    }
}

/* ================================================================
 *  det2d_dpqRK4_GPU: GPU replacement for det2d_dpqRK4.
 *  Uses CellCSR nearest-neighbor search instead of k-d tree.
 *  Computes w2ceil and clips w2 for all av_modes.
 * ================================================================ */
void det2d_dpqRK4_GPU(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    void mkLinkedList2D(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)))
{
    double t0 = MPI_Wtime();

    /* --- Build cell linked list --- */
    postype cellsize = BASICCELL_CELLWIDTH(simpar);
    int mx, my;
    BASICCELL_MX(simpar) = mx = (int)ceil((xmax - xmin) / cellsize);
    BASICCELL_MY(simpar) = my = (int)ceil((ymax - ymin) / cellsize);
    fprintf(stderr,"[DPQ] P%d enter mx=%d my=%d cellsize=%g\n", MYID(simpar), mx, my, cellsize); fflush(stderr);
    VORO_BASICCELL(simpar) =
        (HydroTreeLinkedCell *)my_malloc(sizeof(HydroTreeLinkedCell) * mx * my);
    fprintf(stderr,"[DPQ] P%d before mkLinkedList2D\n", MYID(simpar)); fflush(stderr);
    mkLinkedList2D(simpar, cellsize, xmin, ymin, xmax, ymax,
                   paddingAllTreeParticles);
    fprintf(stderr,"[DPQ] P%d after mkLinkedList2D npad=%ld\n", MYID(simpar), (long)VORO_NPAD(simpar)); fflush(stderr);

    int nbp  = VORO_NP(simpar);
    int npad = VORO_NPAD(simpar);
    int n_total = nbp + npad;
    size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
    fprintf(stderr,"[DPQ] P%d nbp=%d npad=%d n_total=%d p_size=%zu\n", MYID(simpar), nbp, npad, n_total, p_size); fflush(stderr);

    /* --- Lazy GPU context initialization --- */
    if (!g_gpu_ctx.initialized) {
        int max_p = (int)(n_total * 1.2) + 1024;
        int max_f = max_p * 20;
        fprintf(stderr,"[DPQ] P%d gpu_init max_p=%d\n", MYID(simpar), max_p); fflush(stderr);
        gpu_init(&g_gpu_ctx, max_p, max_f, MYID(simpar));
        fprintf(stderr,"[DPQ] P%d gpu_init done\n", MYID(simpar)); fflush(stderr);
    }
    if (n_total > g_gpu_ctx.max_particles) {
        gpu_free(&g_gpu_ctx);
        int max_p = (int)(n_total * 1.5) + 1024;
        int max_f = max_p * 20;
        gpu_init(&g_gpu_ctx, max_p, max_f, MYID(simpar));
    }

    /* --- Build CellCSR --- */
    fprintf(stderr,"[DPQ] P%d before gpu_build_cell_csr\n", MYID(simpar)); fflush(stderr);
    gpu_build_cell_csr(&g_gpu_ctx, simpar, mx, my, p_size);
    fprintf(stderr,"[DPQ] P%d after gpu_build_cell_csr\n", MYID(simpar)); fflush(stderr);

    /* --- Fill SoA and upload (only x, y, indx, w2, w2ceil needed) --- */
    int has_stress = (GAS_AVMODE(simpar) >= 1) ? 1 : 0;
    fillParticleSoA(simpar, &g_gpu_ctx.h_parts, has_stress);
    fprintf(stderr,"[DPQ] P%d after fillParticleSoA has_stress=%d\n", MYID(simpar), has_stress); fflush(stderr);
    gpu_upload_particles(&g_gpu_ctx, n_total, has_stress);
    gpu_upload_cells(&g_gpu_ctx, mx * my, g_gpu_ctx.h_cells.n_entries);
    fprintf(stderr,"[DPQ] P%d after upload\n", MYID(simpar)); fflush(stderr);

    /* --- Launch nearest-neighbor kernel --- */
    GPUPhysicsParams params;
    memset(&params, 0, sizeof(params));
    params.cellsize = cellsize;
    params.mx       = mx;
    params.my       = my;
    params.xmin     = xmin;
    params.ymin     = ymin;

    double t1 = MPI_Wtime();
    fprintf(stderr,"[DPQ] P%d before nearest_neighbor kernel\n", MYID(simpar)); fflush(stderr);
    gpu_launch_nearest_neighbor(&g_gpu_ctx, nbp, n_total, &params,
                                (float)GAS_Kappa(simpar));
    fprintf(stderr,"[DPQ] P%d after nearest_neighbor kernel\n", MYID(simpar)); fflush(stderr);

    /* --- Download w2ceil + w2 --- */
    gpu_download_w2ceil(&g_gpu_ctx, nbp);
    fprintf(stderr,"[DPQ] P%d after download\n", MYID(simpar)); fflush(stderr);
    double t2 = MPI_Wtime();

    /* --- Write back w2ceil + w2 to particle struct --- */
    char *bp_raw = (char *)VORORK4_TBP(simpar);
    fprintf(stderr,"[DPQ] P%d before writeback bp_raw=%p p_size=%zu nbp=%d w2ceil_ptr=%p w2_ptr=%p\n",
        MYID(simpar), (void*)bp_raw, p_size, nbp,
        (void*)g_gpu_ctx.h_parts.w2ceil, (void*)g_gpu_ctx.h_parts.w2); fflush(stderr);
    int i;
#pragma omp parallel for schedule(static)
    for (i = 0; i < nbp; i++) {
        treevorork4particletype *bpi =
            (treevorork4particletype *)(bp_raw + i * p_size);
        bpi->w2ceil = g_gpu_ctx.h_parts.w2ceil[i];
        bpi->w2     = g_gpu_ctx.h_parts.w2[i];
    }
    fprintf(stderr,"[DPQ] P%d after writeback\n", MYID(simpar)); fflush(stderr);

    /* --- Free cells + padding --- */
    fprintf(stderr,"[DPQ] P%d before free basiccell=%p tbpp=%p\n",
        MYID(simpar), (void*)VORO_BASICCELL(simpar), (void*)VORORK4_TBPP(simpar)); fflush(stderr);
    my_free(VORO_BASICCELL(simpar));
    my_free(VORORK4_TBPP(simpar));
    fprintf(stderr,"[DPQ] P%d after free\n", MYID(simpar)); fflush(stderr);

    double t3 = MPI_Wtime();
    if (MYID(simpar) == 0) {
        static int dpq_count = 0;
        dpq_count++;
        printf("[GPU dpq %d] np=%d kernel=%.1fms writeback=%.1fms total=%.1fms\n",
               dpq_count, nbp,
               (t2 - t1) * 1e3, (t3 - t2) * 1e3, (t3 - t0) * 1e3);
    }
}

/* ================================================================
 *  LagMFM (av_mode=4) GPU wrapper functions
 * ================================================================ */

/* Write back density kernel results to CPU particle struct */
static void writeBackLagmfmDensity(
    SimParameters *simpar,
    const ParticleSoA *h_parts)
{
    size_t p_size  = TVORORK4_DDINFO(simpar)[0].n_size;
    char *real_base = (char *)VORORK4_TBP(simpar);
    int nbp = VORO_NP(simpar);
    int i;

#pragma omp parallel for schedule(static)
    for (i = 0; i < nbp; i++) {
        treevorostressrk4particletype *sbp =
            (treevorostressrk4particletype *)(real_base + i * p_size);

        sbp->den      = h_parts->den[i];
        sbp->volume   = h_parts->volume[i];
        sbp->pressure = h_parts->pressure[i];
        sbp->csound   = h_parts->csound[i];

        sbp->stress.E_inv_xx = h_parts->E_inv_xx[i];
        sbp->stress.E_inv_xy = h_parts->E_inv_xy[i];
        sbp->stress.E_inv_yx = h_parts->E_inv_yx[i];
        sbp->stress.E_inv_yy = h_parts->E_inv_yy[i];
        sbp->stress.h_mfm    = h_parts->h_mfm[i];

        sbp->stress.gUxx = h_parts->gUxx[i];
        sbp->stress.gUxy = h_parts->gUxy[i];
        sbp->stress.gUyx = h_parts->gUyx[i];
        sbp->stress.gUyy = h_parts->gUyy[i];
        sbp->stress.dPdx = h_parts->dPdx[i];
        sbp->stress.dPdy = h_parts->dPdy[i];
        sbp->stress.divv = h_parts->divv[i];

        sbp->stress.tauxx = h_parts->tauxx[i];
        sbp->stress.tauxy = h_parts->tauxy[i];
        sbp->stress.tauyy = h_parts->tauyy[i];

        sbp->w2ceil = h_parts->w2ceil[i];
    }
}

/* ================================================================
 *  updateDenW2Pressure2D_LagMFM_GPU:
 *    GPU version of the LagMFM density pass.
 *    Called from exam2d_vph_rk4_int_lagmfm in place of
 *    updateDenW2Pressure2D_LagMFM.
 *
 *    Handles: cell list build, SoA fill, CellCSR build,
 *             GPU upload, density kernel, download, writeback.
 *    The cell list + padding are kept alive for the force call.
 * ================================================================ */
void updateDenW2Pressure2D_LagMFM_GPU(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    postype Gamma,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    void mkLinkedList2D(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)),
    postype Dtime)
{
    double t0 = MPI_Wtime();

    /* NOTE: w2 refresh (det2d_dpqRK4 + getw2forHydroParticle + applyW2Controls)
     * is done in exam.c before this call, since those functions are static/inline
     * in exam.c and not accessible from this file. */

    /* --- Build cell linked list --- */
    postype cellsize = BASICCELL_CELLWIDTH(simpar);
    int mx, my;
    BASICCELL_MX(simpar) = mx = (int)ceil((xmax - xmin) / cellsize);
    BASICCELL_MY(simpar) = my = (int)ceil((ymax - ymin) / cellsize);
    VORO_BASICCELL(simpar) =
        (HydroTreeLinkedCell *)my_malloc(sizeof(HydroTreeLinkedCell) * mx * my);
    mkLinkedList2D(simpar, cellsize, xmin, ymin, xmax, ymax,
                   paddingAllTreeParticles);

    int nbp  = VORO_NP(simpar);
    int npad = VORO_NPAD(simpar);
    int n_total = nbp + npad;
    size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;

    /* --- Lazy GPU context initialization --- */
    if (!g_gpu_ctx.initialized) {
        int max_p = (int)(n_total * 1.2) + 1024;
        int max_f = max_p * 20;
        gpu_init(&g_gpu_ctx, max_p, max_f, MYID(simpar));
    }

    if (n_total > g_gpu_ctx.max_particles) {
        gpu_free(&g_gpu_ctx);
        int max_p = (int)(n_total * 1.5) + 1024;
        int max_f = max_p * 20;
        gpu_init(&g_gpu_ctx, max_p, max_f, MYID(simpar));
    }

    /* --- Build CellCSR --- */
    gpu_build_cell_csr(&g_gpu_ctx, simpar, mx, my, p_size);

    /* --- Fill SoA and upload --- */
    fillParticleSoA(simpar, &g_gpu_ctx.h_parts, 1 /* has_stress */);
    gpu_upload_particles(&g_gpu_ctx, n_total, 1);
    gpu_upload_lagmfm_fields(&g_gpu_ctx, n_total);
    gpu_upload_cells(&g_gpu_ctx, mx * my,
                     g_gpu_ctx.h_cells.n_entries);

    /* --- Launch density kernel --- */
    GPUPhysicsParams params;
    memset(&params, 0, sizeof(params));
    params.av_mode       = GAS_AVMODE(simpar);
    params.Gamma         = Gamma;
    params.nu_phys       = GAS_VISCOSITY(simpar);
    params.lagmfm_eta      = 1.4;  /* Hopkins 2015 fiducial */
    params.lagmfm_h_iter   = 50;
    params.lagmfm_h_tol    = 0.001;
    params.lagmfm_pure_gizmo = 1;
    params.cellsize      = cellsize;
    params.mx            = mx;
    params.my            = my;
    params.xmin          = xmin;
    params.ymin          = ymin;

    double t1 = MPI_Wtime();
    gpu_launch_lagmfm_density_kernel(&g_gpu_ctx, nbp, n_total, &params);

    /* --- Download density results --- */
    gpu_download_lagmfm_density(&g_gpu_ctx, nbp);

    double t2 = MPI_Wtime();

    /* --- Write back to CPU particle struct --- */
    writeBackLagmfmDensity(simpar, &g_gpu_ctx.h_parts);

    /* --- Free cells but KEEP padding alive for force call --- */
    my_free(VORO_BASICCELL(simpar));

    double t3 = MPI_Wtime();
    if (MYID(simpar) == 0) {
        static int den_count = 0;
        den_count++;
        if (den_count % 4 == 0)
            printf("[GPU lagmfm den %d] np=%d kernel=%.1fms download+wb=%.1fms "
                   "total=%.1fms\n",
                   den_count/4, nbp,
                   (t2 - t1) * 1e3, (t3 - t2) * 1e3, (t3 - t0) * 1e3);
    }
}

/* ================================================================
 *  getAccVoro2D_LagMFM_GPU:
 *    GPU version of the LagMFM force computation.
 *    Reuses the SoA and CellCSR from the density call above.
 *    (Padding was kept alive; cells need to be rebuilt.)
 * ================================================================ */
double getAccVoro2D_LagMFM_GPU(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    postype OrderOfAccuracy, postype Courant, postype Gamma,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    void mkLinkedList2D(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)))
{
    double t0 = MPI_Wtime();

    int nbp  = VORO_NP(simpar);
    int npad = VORO_NPAD(simpar);
    int n_total = nbp + npad;
    size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;

    postype cellsize = BASICCELL_CELLWIDTH(simpar);
    int mx, my;
    BASICCELL_MX(simpar) = mx = (int)ceil((xmax - xmin) / cellsize);
    BASICCELL_MY(simpar) = my = (int)ceil((ymax - ymin) / cellsize);
    /* LagMFM density kept TBPP alive; free before rebuild to avoid leak. */
    if (VORORK4_TBPP(simpar) != NULL) {
        my_free(VORORK4_TBPP(simpar));
        VORORK4_TBPP(simpar) = NULL;
    }
    VORO_BASICCELL(simpar) =
        (HydroTreeLinkedCell *)my_malloc(sizeof(HydroTreeLinkedCell) * mx * my);
    mkLinkedList2D(simpar, cellsize, xmin, ymin, xmax, ymax,
                   paddingAllTreeParticles);
    npad = VORO_NPAD(simpar);
    n_total = nbp + npad;

    /* Rebuild CellCSR and SoA for force pass (padding may have changed) */
    gpu_build_cell_csr(&g_gpu_ctx, simpar, mx, my, p_size);
    fillParticleSoA(simpar, &g_gpu_ctx.h_parts, 1);
    gpu_upload_particles(&g_gpu_ctx, n_total, 1);
    gpu_upload_lagmfm_fields(&g_gpu_ctx, n_total);
    gpu_upload_cells(&g_gpu_ctx, mx * my,
                     g_gpu_ctx.h_cells.n_entries);

    /* Physics params */
    postype Lx = SIMBOX(simpar).x.max;
    GPUPhysicsParams params;
    memset(&params, 0, sizeof(params));
    params.av_mode       = GAS_AVMODE(simpar);
    params.Courant       = Courant;
    params.Gamma         = Gamma;
    params.alphavis      = GAS_AlphaVis(simpar);
    params.betavis       = GAS_BetaVis(simpar);
    params.etavis        = GAS_ETAVIS(simpar) * Lx / NX(simpar);
    params.epsvis        = GAS_EPSVIS(simpar);
    params.nu_phys       = GAS_VISCOSITY(simpar);
    params.lagmfm_eta      = 1.4;  /* Hopkins 2015 fiducial */
    params.lagmfm_pure_gizmo = 1;
    params.cellsize      = cellsize;
    params.mx            = mx;
    params.my            = my;
    params.xmin          = xmin;
    params.ymin          = ymin;

    double t1 = MPI_Wtime();
    double Dtime_local = gpu_launch_lagmfm_force_kernel(
        &g_gpu_ctx, nbp, n_total, &params);

    /* Download results */
    gpu_download_results(&g_gpu_ctx, nbp);

    double t2 = MPI_Wtime();

    /* Write back */
    writeBackResults(simpar, &g_gpu_ctx.h_parts, 1);

    /* MPI reduction */
    double Dtime;
    MPI_Allreduce(&Dtime_local, &Dtime, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM(simpar));

    /* Free */
    my_free(VORO_BASICCELL(simpar));
    my_free(VORORK4_TBPP(simpar));

    double t3 = MPI_Wtime();
    if (MYID(simpar) == 0) {
        static int force_count = 0;
        force_count++;
        if (force_count % 4 == 0)
            printf("[GPU lagmfm force %d] np=%d kernel=%.1fms "
                   "download+wb=%.1fms total=%.1fms Dtime=%g\n",
                   force_count/4, nbp,
                   (t2 - t1) * 1e3, (t3 - t2) * 1e3,
                   (t3 - t0) * 1e3, Dtime);
    }

    return Dtime;
}

#endif /* USE_CUDA */
