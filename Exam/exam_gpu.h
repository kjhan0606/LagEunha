#ifndef EXAM_GPU_H
#define EXAM_GPU_H

#ifdef USE_CUDA

#include <stddef.h>

/* ================================================================
 *  GPU Data Structures for Voronoi Force Loop + GPU Tessellation
 *
 *  Structures:
 *    FaceCSR      — face connectivity + geometry in CSR format
 *    ParticleSoA  — particle fields in Structure-of-Arrays layout
 *    CellCSR      — cell linked-list in CSR format (for GPU tessellation)
 *    PolygonOut   — per-particle polygon output (GPU tessellation temp)
 *
 *  GPU tessellation produces FaceCSR directly on device, eliminating
 *  CPU Voronoi + CPU→GPU face transfer.
 * ================================================================ */

/* Maximum polygon corners per particle (box starts with 4,
   each neighbor adds at most 1 → 48 handles ~44 neighbors) */
#define GPU_MAX_CORNERS 48

/* MAX_INDEX for ghost detection (must match eunha.h) */
#define MAX_INDEX_GPU 72057594037927935LL

/* ---- Cell linked-list in CSR format (for GPU tessellation) ---- */
typedef struct {
    int n_cells;              /* mx * my */
    int n_entries;            /* total particles in linked lists (real + padding) */
    int *cell_offset;         /* size: n_cells + 1 (CSR offsets) */
    int *particle_index;      /* size: n_entries (SoA index per entry) */
} CellCSR;

/* ---- Per-particle polygon output from GPU tessellation ---- */
typedef struct {
    double cx[GPU_MAX_CORNERS], cy[GPU_MAX_CORNERS];
    int    edge_neigh[GPU_MAX_CORNERS];  /* SoA index of face neighbor, <0 = box edge */
    int    n_corners;
    int    n_real_faces;                  /* count of edges with edge_neigh >= 0 */
} PolygonOut;

/* ---- CSR (Compressed Sparse Row) face topology ---- */
typedef struct {
    int n_particles;        /* number of real particles processed */
    int n_faces_total;      /* total face count across all particles */

    /* CSR offsets: particle i owns faces [face_offset[i], face_offset[i+1]) */
    int *face_offset;       /* size: n_particles + 1 */

    /* Per-face Voronoi corner coordinates (RELATIVE to owner particle) */
    double *c1x, *c1y;     /* corner 1 (tmp)  */
    double *c2x, *c2y;     /* corner 2 (tmp->upperlink) */

    /* Per-face connectivity */
    int *neighbor_idx;      /* index of neighbor j in ParticleSoA (-1 if invalid) */

    /* M(n,m) 4-point stencil neighbor indices */
    int *kp_idx;            /* (tmp->upperlink)->upperrelated neighbor index (-1 if invalid) */
    int *km_idx;            /* tmp->lowerrelated neighbor index (-1 if invalid) */

    /* Ghost flag: 1 = wall mirror (PINDX == MAX_INDEX) */
    int *is_ghost;
} FaceCSR;

/* ---- Structure-of-Arrays particle layout for GPU ---- */
typedef struct {
    int n_total;            /* real + padding particles */

    /* Position / velocity / mass (PosType = double with XYZDBL) */
    double *x, *y;
    double *vx, *vy;
    double *mass;

    /* Particle index for ghost detection (PINDX) */
    long long *indx;

    /* VoroQ fields (always float in CPU struct) */
    float *den, *pressure, *csound, *volume;
    float *ie, *w2, *w2old;
    float *w2ceil;                          /* Laguerre ceiling */
    float *avgNeighP;                       /* avg neighboring pressure */

    /* Stress fields (PosType = double) — only used when av_mode >= 1 */
    double *gUxx, *gUxy, *gUyx, *gUyy;   /* velocity gradient tensor */
    double *dPdx, *dPdy;                   /* pressure gradient (MUSCL) */
    double *dRhodx, *dRhody;               /* density gradient (Arepo MUSCL) */
    double *tauxx, *tauxy, *tauyy;         /* NS viscous stress */
    double *alpha_cd;                       /* CD10 viscosity coefficient */
    double *divv;                           /* divergence of velocity */

    /* LagMFM (av_mode=4) fields — only allocated when needed */
    double *E_inv_xx, *E_inv_xy, *E_inv_yx, *E_inv_yy;  /* inverse moment matrix */
    double *h_mfm;                          /* kernel bandwidth */

    /* Output arrays (written by GPU kernel, read back to CPU) */
    double *ax_out, *ay_out;
    float  *die_out;
    float  *dt_out;
    double *vsig_max_out;                   /* for CD10 */
} ParticleSoA;

/* ---- GPU memory manager (persistent across RK4 stages) ---- */
typedef struct {
    /* Device-side copies */
    FaceCSR    d_faces;     /* device pointers inside */
    ParticleSoA d_parts;    /* device pointers inside */

    /* Host-side staging (pinned memory for async transfer) */
    FaceCSR    h_faces;     /* pinned host pointers */
    ParticleSoA h_parts;    /* pinned host pointers */

    /* Allocation capacities */
    int max_particles;
    int max_faces;

    /* GPU tessellation buffers (persistent) */
    CellCSR  d_cells;       /* device cell CSR */
    CellCSR  h_cells;       /* pinned host cell CSR */
    int max_cells;           /* allocated n_cells capacity */
    int max_cell_entries;    /* allocated n_entries capacity */

    PolygonOut *d_polygons;  /* device: max_particles PolygonOut structs */
    int *d_face_count;       /* device: per-particle real face count */
    int *d_face_offset;      /* device: prefix-summed face offsets (n_particles+1) */

    /* CUB scratch for prefix sum (tessellation) */
    void  *d_cub_scan_temp;
    size_t cub_scan_temp_bytes;

    /* Flag: GPU tessellation produced faces on device, skip CPU extraction */
    int tess_faces_on_device;

    /* CUB scratch (persistent, avoids cudaMalloc per kernel call) */
    void *d_cub_temp;
    void *d_cub_min_out;     /* float* on device */
    void *h_cub_min_out;     /* float* on pinned host */
    size_t cub_temp_bytes;

    /* CUDA stream for async transfers */
    void *stream;            /* cudaStream_t, opaque to C code */

    int initialized;
} GPUContext;

/* ---- Physics parameters passed to GPU kernel ---- */
typedef struct {
    int    av_mode;
    int    use_muscl;
    double OoA;            /* OrderOfAccuracy */
    double Courant;
    double Gamma;
    double alphavis, betavis, etavis, epsvis;
    double nu_phys;
    double prandtl;
    double cd_amax;
    double blend_theta;
    double dtold;
    /* LagMFM (av_mode=4) */
    double lagmfm_eta;       /* 1.8: h = eta * sqrt(V) */
    int    lagmfm_h_iter;    /* 3: Newton iterations for adaptive h */
    double lagmfm_h_tol;     /* 0.005: convergence tolerance */
    double cellsize;         /* cell linked-list cell width */
    int    mx, my;           /* cell grid dimensions */
    double xmin, ymin;       /* domain lower-left corner */
} GPUPhysicsParams;

/* ================================================================
 *  Function declarations
 * ================================================================ */

#ifdef __cplusplus
extern "C" {
#endif

/* --- GPU memory management (exam_gpu.cu) --- */
void gpu_init(GPUContext *ctx, int max_particles, int max_faces, int mpi_rank);
void gpu_free(GPUContext *ctx);

/* --- Upload / download (exam_gpu.cu) --- */
void gpu_upload_particles(GPUContext *ctx, int n_total, int has_stress);
void gpu_upload_faces(GPUContext *ctx, int n_particles, int n_faces);
void gpu_download_results(GPUContext *ctx, int n_particles);

/* --- Force kernel launch (exam_gpu.cu) --- */
double gpu_launch_force_kernel(GPUContext *ctx, int n_particles,
                               const GPUPhysicsParams *params);

/* --- GPU tessellation (exam_gpu.cu) --- */
void gpu_alloc_tess_buffers(GPUContext *ctx, int max_particles,
                            int max_cells, int max_cell_entries);
void gpu_free_tess_buffers(GPUContext *ctx);
void gpu_upload_cells(GPUContext *ctx, int n_cells, int n_entries);
void gpu_tessellate_and_build_faces(GPUContext *ctx, int nbp,
    int mx, int my, double cellsize, double xmin, double ymin,
    double boxsize, double OoA, int av_mode, double Gamma,
    double nu_phys, double prandtl, double cd_amax);
void gpu_download_tess_results(GPUContext *ctx, int nbp, int has_stress);

/* --- CPU-side helpers for GPU tessellation (exam_gpu_extract.c) --- */
GPUContext *gpu_get_context(void);
void gpu_build_cell_csr(GPUContext *ctx, void *simpar_opaque,
                        int mx, int my, size_t p_size);
void gpu_writeback_tess_results(void *simpar_opaque, GPUContext *ctx,
                                int nbp, int has_stress);

/* --- Face extraction and GPU wrapper are declared in exam_gpu_extract.c ---
 *  getAccVoro2DBlend_GPU() is called from exam.c; its prototype is provided
 *  there via a local forward declaration to avoid pulling voro.h into this
 *  header (which is also included by the CUDA .cu file).
 */

/* --- Fused tessellation+extraction API (exam_gpu_extract.c) ---
 *  Allows updateDenW2Pressure2DBlend to extract CSR faces during
 *  the gradient face-ring traversal, eliminating the second tessellation
 *  that getAccVoro2DBlend_GPU would otherwise perform.
 */
void gpu_fused_begin(int nbp, int npad, size_t p_size, int nthreads,
                     const char *real_base, const char *pad_base);
void gpu_fused_add_face(int tid, int owner_idx,
    double c1x, double c1y, double c2x, double c2y,
    int neighbor_idx, int kp_idx, int km_idx, int is_ghost);
void gpu_fused_end(void);
int  gpu_fused_ready(void);

int gpu_ptr_to_soa_index(const void *bp,
    const char *real_base, const char *pad_base,
    int nbp, int npad, size_t p_size);

/* --- Download GPU face CSR to host for comparison --- */
void gpu_download_face_csr(GPUContext *ctx, int n_particles, int n_faces,
    int *h_face_offset, double *h_c1x, double *h_c1y, double *h_c2x, double *h_c2y,
    int *h_neighbor_idx, int *h_kp_idx, int *h_km_idx, int *h_is_ghost);

/* --- GPU nearest-neighbor (replaces det2d_dpqRK4) --- */
void gpu_launch_nearest_neighbor(GPUContext *ctx, int n_particles, int n_total,
                                  const GPUPhysicsParams *params, float kappa);
void gpu_download_w2ceil(GPUContext *ctx, int n_particles);

/* --- LagMFM (av_mode=4) GPU kernels (exam_gpu.cu) --- */
void gpu_launch_lagmfm_density_kernel(GPUContext *ctx, int n_particles,
                                       int n_total,
                                       const GPUPhysicsParams *params);
double gpu_launch_lagmfm_force_kernel(GPUContext *ctx, int n_particles,
                                       int n_total,
                                       const GPUPhysicsParams *params);
void gpu_upload_lagmfm_fields(GPUContext *ctx, int n_total);
void gpu_download_lagmfm_density(GPUContext *ctx, int n_particles);

#ifdef __cplusplus
}
#endif

#endif /* USE_CUDA */
#endif /* EXAM_GPU_H */
