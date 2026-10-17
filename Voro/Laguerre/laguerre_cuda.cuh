/**
 * @file laguerre_cuda.cuh
 * @brief CUDA header for parallel Laguerre diagram computation
 */

#ifndef LAGUERRE_CUDA_CUH
#define LAGUERRE_CUDA_CUH

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <math.h>

#define CUDA_BLOCK_SIZE 256
#define MAX_NEIGHBORS 64
#define MAX_VERTICES 128
#define MAX_CORNERS 64

#ifdef LAGUERRE_SINGLE
    typedef float real_t;
    #define EPSILON         (1.0e-6f)
    #define BOUNDARY_TOL    (1.0e-6f)
    #define CUDA_SQRT       sqrtf
    #define CUDA_FABS       fabsf
#else
    typedef double real_t;
    #define EPSILON         (1.0e-9)
    #define BOUNDARY_TOL    (1.0e-6)
    #define CUDA_SQRT       sqrt
    #define CUDA_FABS       fabs
#endif

#define CUDA_MIN(a, b)  ((a) < (b) ? (a) : (b))
#define CUDA_MAX(a, b)  ((a) > (b) ? (a) : (b))

#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
} while(0)

enum ClassifyResult { CLASSIFY_INSIDE = -1, CLASSIFY_BOUNDARY = -1, CLASSIFY_OUTSIDE = 1 };
enum VertexStatus { STATUS_INACTIVE = -1, STATUS_ACTIVE = 1 };

struct __align__(16) Point3D_GPU {
    real_t x, y, z, w2, dist2;
    int index, padding;
};

struct __align__(16) Point2D_GPU {
    real_t x, y, w2, dist2;
    int index, padding;
};

struct __align__(16) Vertex3D_GPU {
    real_t x, y, z;
    int status, neighbor_sites[3], adjacent[3];
};

struct __align__(16) Corner2D_GPU {
    real_t x, y;
    int status, lower_neighbor, upper_neighbor, lower_link, upper_link, padding;
};

struct CellResult3D {
    real_t volume, centroid_x, centroid_y, centroid_z;
    int num_vertices, num_faces;
};

struct CellResult2D {
    real_t area, centroid_x, centroid_y;
    int num_corners;
};

__device__ __forceinline__ real_t dot3d_gpu(real_t ax, real_t ay, real_t az, real_t bx, real_t by, real_t bz) {
    return ax * bx + ay * by + az * bz;
}

__device__ __forceinline__ real_t dot2d_gpu(real_t ax, real_t ay, real_t bx, real_t by) {
    return ax * bx + ay * by;
}

__device__ __forceinline__ real_t compute_weight_fraction_gpu(real_t w1_sq, real_t w2_sq, real_t dist_sq) {
    return 0.5 + 0.5 * (w1_sq - w2_sq) / dist_sq;
}

__device__ int classify_vertex_3d_gpu(const Point3D_GPU *neighbor, real_t vx, real_t vy, real_t vz, real_t wfrac, real_t *offset);
__device__ int classify_corner_2d_gpu(const Point2D_GPU *neighbor, real_t cx, real_t cy, real_t wfrac);
__device__ void find_edge_intersection_3d_gpu(real_t ax, real_t ay, real_t az, real_t bx, real_t by, real_t bz, const Point3D_GPU *neighbor, real_t wfrac, real_t *ix, real_t *iy, real_t *iz);
__device__ void find_edge_intersection_2d_gpu(real_t ax, real_t ay, real_t bx, real_t by, const Point2D_GPU *neighbor, real_t wfrac, real_t *ix, real_t *iy, real_t *t);
__device__ void intersect_three_planes_gpu(real_t n1x, real_t n1y, real_t n1z, real_t d1, real_t n2x, real_t n2y, real_t n2z, real_t d2, real_t n3x, real_t n3y, real_t n3z, real_t d3, real_t *ix, real_t *iy, real_t *iz);

__global__ void construct_cells_3d_kernel(const Point3D_GPU *sites, int num_sites, const Point3D_GPU *neighbors, const int *neighbor_counts, const int *neighbor_offsets, CellResult3D *results, real_t box_size);
__global__ void construct_cells_2d_kernel(const Point2D_GPU *sites, int num_sites, const Point2D_GPU *neighbors, const int *neighbor_counts, const int *neighbor_offsets, CellResult2D *results, real_t box_size);

#ifdef __cplusplus
extern "C" {
#endif

int laguerre_cuda_init(int device_id);
void laguerre_cuda_cleanup(void);
int laguerre_construct_cells_3d(const Point3D_GPU *h_sites, int num_sites, const Point3D_GPU *h_neighbors, const int *h_neighbor_counts, CellResult3D *h_results, real_t box_size);
int laguerre_construct_cells_2d(const Point2D_GPU *h_sites, int num_sites, const Point2D_GPU *h_neighbors, const int *h_neighbor_counts, CellResult2D *h_results, real_t box_size);
void laguerre_cuda_get_device_info(int *major, int *minor, size_t *mem, int *mp);

#ifdef __cplusplus
}
#endif

#endif
