/**
 * @file laguerre_cuda.cu
 * @brief CUDA implementation of parallel Laguerre diagram construction
 */

#include "laguerre_cuda.cuh"

/*============================================================================
 * DEVICE HELPER FUNCTIONS
 *============================================================================*/

__device__ void init_bounding_box_3d_local(Vertex3D_GPU *vertices, int *num_vertices, real_t box_size) {
    real_t half = 0.5 * box_size;
    for (int i = 0; i < 8; i++) {
        int j = i % 4, k = i / 4;
        vertices[i].x = box_size * (((j + 1) % 4) / 2) - half;
        vertices[i].y = box_size * (j / 2) - half;
        vertices[i].z = box_size * k - half;
        vertices[i].status = STATUS_ACTIVE;
        vertices[i].adjacent[0] = (j - 1 + 4) % 4 + k * 4;
        if (i < 4) {
            vertices[i].neighbor_sites[0] = -(j + 1);
            vertices[i].neighbor_sites[1] = -((j + 3) % 4 + 1);
            vertices[i].neighbor_sites[2] = -5;
            vertices[i].adjacent[1] = (j + 1) % 4 + k * 4;
            vertices[i].adjacent[2] = (i + 4) % 8;
        } else {
            vertices[i].neighbor_sites[0] = -(j + 1);
            vertices[i].neighbor_sites[1] = -6;
            vertices[i].neighbor_sites[2] = -((j + 3) % 4 + 1);
            vertices[i].adjacent[2] = (j + 1) % 4 + k * 4;
            vertices[i].adjacent[1] = (i + 4) % 8;
        }
    }
    *num_vertices = 8;
}

__device__ void init_bounding_box_2d_local(Corner2D_GPU *corners, int *num_corners, real_t box_size) {
    real_t half = 0.5 * box_size;
    for (int i = 0; i < 4; i++) {
        corners[i].x = box_size * (((i + 1) % 4) / 2) - half;
        corners[i].y = box_size * (i / 2) - half;
        corners[i].status = STATUS_ACTIVE;
        corners[i].lower_link = (i - 1 + 4) % 4;
        corners[i].upper_link = (i + 1) % 4;
        corners[i].upper_neighbor = -(i + 1);
        corners[i].lower_neighbor = -((i + 3) % 4 + 1);
    }
    *num_corners = 4;
}

__device__ int classify_vertex_3d_gpu(const Point3D_GPU *neighbor, real_t vx, real_t vy, real_t vz, real_t wfrac, real_t *offset) {
    real_t proj = dot3d_gpu(vx, vy, vz, neighbor->x, neighbor->y, neighbor->z);
    real_t n_sq = dot3d_gpu(neighbor->x, neighbor->y, neighbor->z, neighbor->x, neighbor->y, neighbor->z);
    real_t thresh = wfrac * n_sq;
    if (CUDA_FABS(thresh) < EPSILON) { *offset = -1; return CLASSIFY_INSIDE; }
    if (proj / thresh > 1 + BOUNDARY_TOL) { *offset = proj - thresh; return CLASSIFY_OUTSIDE; }
    *offset = -1; return CLASSIFY_INSIDE;
}

__device__ int classify_corner_2d_gpu(const Point2D_GPU *neighbor, real_t cx, real_t cy, real_t wfrac) {
    real_t proj = dot2d_gpu(cx, cy, neighbor->x, neighbor->y);
    real_t thresh = wfrac * dot2d_gpu(neighbor->x, neighbor->y, neighbor->x, neighbor->y);
    if (proj > thresh) return CLASSIFY_OUTSIDE;
    if (proj < thresh) return CLASSIFY_INSIDE;
    return CLASSIFY_BOUNDARY;
}

__device__ void find_edge_intersection_3d_gpu(real_t ax, real_t ay, real_t az, real_t bx, real_t by, real_t bz,
                                              const Point3D_GPU *neighbor, real_t wfrac, real_t *ix, real_t *iy, real_t *iz) {
    real_t dx = bx - ax, dy = by - ay, dz = bz - az;
    real_t e_dot_n = dot3d_gpu(dx, dy, dz, neighbor->x, neighbor->y, neighbor->z);
    real_t a_dot_n = dot3d_gpu(ax, ay, az, neighbor->x, neighbor->y, neighbor->z);
    real_t n_sq = dot3d_gpu(neighbor->x, neighbor->y, neighbor->z, neighbor->x, neighbor->y, neighbor->z);
    real_t t = (CUDA_FABS(e_dot_n) < EPSILON) ? wfrac : (wfrac * n_sq - a_dot_n) / e_dot_n;
    t = CUDA_MIN(CUDA_MAX(t, EPSILON), 1.0 - EPSILON);
    *ix = ax + t * dx; *iy = ay + t * dy; *iz = az + t * dz;
}

__device__ void find_edge_intersection_2d_gpu(real_t ax, real_t ay, real_t bx, real_t by,
                                              const Point2D_GPU *neighbor, real_t wfrac, real_t *ix, real_t *iy, real_t *t_out) {
    real_t dx = bx - ax, dy = by - ay;
    real_t e_dot_n = dot2d_gpu(dx, dy, neighbor->x, neighbor->y);
    real_t a_dot_n = dot2d_gpu(ax, ay, neighbor->x, neighbor->y);
    real_t n_sq = dot2d_gpu(neighbor->x, neighbor->y, neighbor->x, neighbor->y);
    real_t t = (CUDA_FABS(e_dot_n) < EPSILON) ? wfrac : (wfrac * n_sq - a_dot_n) / e_dot_n;
    t = CUDA_MIN(CUDA_MAX(t, EPSILON), 1.0 - EPSILON);
    *ix = ax + t * dx; *iy = ay + t * dy; *t_out = t;
}

__device__ void intersect_three_planes_gpu(real_t n1x, real_t n1y, real_t n1z, real_t d1,
                                           real_t n2x, real_t n2y, real_t n2z, real_t d2,
                                           real_t n3x, real_t n3y, real_t n3z, real_t d3,
                                           real_t *ix, real_t *iy, real_t *iz) {
    real_t det = n1x * (n2y * n3z - n3y * n2z) - n1y * (n2x * n3z - n3x * n2z) + n1z * (n2x * n3y - n3x * n2y);
    if (CUDA_FABS(det) < EPSILON) { *ix = *iy = *iz = 0; return; }
    real_t inv = 1.0 / det;
    *ix = (d1 * (n2y * n3z - n3y * n2z) - n1y * (d2 * n3z - d3 * n2z) + n1z * (d2 * n3y - d3 * n2y)) * inv;
    *iy = (n1x * (d2 * n3z - d3 * n2z) - d1 * (n2x * n3z - n3x * n2z) + n1z * (n2x * d3 - n3x * d2)) * inv;
    *iz = (n1x * (n2y * d3 - n3y * d2) - n1y * (n2x * d3 - n3x * d2) + d1 * (n2x * n3y - n3x * n2y)) * inv;
}

__device__ int compact_vertices_local(Vertex3D_GPU *vertices, int num_vertices, int *index_map) {
    int new_count = 0;
    for (int i = 0; i < num_vertices; i++)
        index_map[i] = (vertices[i].status == STATUS_ACTIVE) ? new_count++ : -1;
    int write_idx = 0;
    for (int i = 0; i < num_vertices; i++) {
        if (vertices[i].status == STATUS_ACTIVE) {
            if (write_idx != i) vertices[write_idx] = vertices[i];
            for (int j = 0; j < 3; j++)
                if (vertices[i].adjacent[j] >= 0 && vertices[i].adjacent[j] < num_vertices)
                    vertices[write_idx].adjacent[j] = index_map[vertices[i].adjacent[j]];
            write_idx++;
        }
    }
    return new_count;
}

__device__ real_t compute_volume_local(const Vertex3D_GPU *vertices, int num_vertices) {
    if (num_vertices < 4) return 0;
    real_t vol = 0;
    for (int v = 0; v < num_vertices; v++) {
        for (int e = 0; e < 3; e++) {
            int next = vertices[v].adjacent[e];
            if (next < 0 || next >= num_vertices || next <= v) continue;
            real_t sx = 0, sy = 0, sz = 0;
            int prev = v, curr = next, cnt = 0;
            do {
                sx += vertices[prev].y * vertices[curr].z - vertices[prev].z * vertices[curr].y;
                sy += vertices[prev].z * vertices[curr].x - vertices[prev].x * vertices[curr].z;
                sz += vertices[prev].x * vertices[curr].y - vertices[prev].y * vertices[curr].x;
                int temp = curr;
                for (int k = 0; k < 3; k++) {
                    if (vertices[curr].adjacent[k] == prev) { curr = vertices[curr].adjacent[(k + 2) % 3]; break; }
                }
                prev = temp;
            } while (curr != next && ++cnt < MAX_VERTICES);
            vol += CUDA_SQRT(sx*sx + sy*sy + sz*sz) / 6.0;
        }
    }
    return vol / 2.0;
}

__device__ real_t compute_area_local(const Corner2D_GPU *corners, int start, int num_active) {
    if (num_active < 3) return 0;
    real_t area = 0;
    int curr = start, cnt = 0;
    do {
        int next = corners[curr].upper_link;
        area += corners[curr].x * corners[next].y - corners[next].x * corners[curr].y;
        curr = next;
    } while (curr != start && ++cnt < MAX_CORNERS);
    return CUDA_FABS(area) * 0.5;
}

__device__ void compute_centroid_local(const Corner2D_GPU *corners, int start, real_t *cx, real_t *cy) {
    real_t area = 0, sx = 0, sy = 0;
    int curr = start, cnt = 0;
    do {
        int next = corners[curr].upper_link;
        real_t cross = corners[curr].x * corners[next].y - corners[next].x * corners[curr].y;
        area += cross;
        sx += (corners[curr].x + corners[next].x) * cross;
        sy += (corners[curr].y + corners[next].y) * cross;
        curr = next;
    } while (curr != start && ++cnt < MAX_CORNERS);
    area *= 0.5;
    *cx = (CUDA_FABS(area) > EPSILON) ? sx / (6 * area) : 0;
    *cy = (CUDA_FABS(area) > EPSILON) ? sy / (6 * area) : 0;
}

/*============================================================================
 * CUDA KERNELS
 *============================================================================*/

__global__ void construct_cells_3d_kernel(const Point3D_GPU *sites, int num_sites,
                                          const Point3D_GPU *neighbors, const int *neighbor_counts,
                                          const int *neighbor_offsets, CellResult3D *results, real_t box_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_sites) return;
    
    Vertex3D_GPU vertices[MAX_VERTICES];
    int index_map[MAX_VERTICES];
    int num_vertices;
    
    Point3D_GPU center = sites[idx];
    int num_neigh = neighbor_counts[idx];
    int neigh_start = neighbor_offsets[idx];
    
    init_bounding_box_3d_local(vertices, &num_vertices, box_size);
    
    for (int n = 0; n < num_neigh && n < MAX_NEIGHBORS; n++) {
        Point3D_GPU neighbor = neighbors[neigh_start + n];
        neighbor.x -= center.x; neighbor.y -= center.y; neighbor.z -= center.z;
        real_t dist_sq = neighbor.x * neighbor.x + neighbor.y * neighbor.y + neighbor.z * neighbor.z;
        if (dist_sq < EPSILON) continue;
        real_t wfrac = compute_weight_fraction_gpu(center.w2, neighbor.w2, dist_sq);
        
        for (int i = 0; i < num_vertices; i++) {
            if (vertices[i].status != STATUS_ACTIVE) continue;
            real_t offset;
            if (classify_vertex_3d_gpu(&neighbor, vertices[i].x, vertices[i].y, vertices[i].z, wfrac, &offset) == CLASSIFY_OUTSIDE)
                vertices[i].status = STATUS_INACTIVE;
        }
        num_vertices = compact_vertices_local(vertices, num_vertices, index_map);
        if (num_vertices < 4) break;
    }
    
    results[idx].num_vertices = num_vertices;
    results[idx].volume = compute_volume_local(vertices, num_vertices);
    real_t cx = 0, cy = 0, cz = 0;
    for (int i = 0; i < num_vertices; i++) { cx += vertices[i].x; cy += vertices[i].y; cz += vertices[i].z; }
    if (num_vertices > 0) {
        results[idx].centroid_x = cx / num_vertices + center.x;
        results[idx].centroid_y = cy / num_vertices + center.y;
        results[idx].centroid_z = cz / num_vertices + center.z;
    }
}

__global__ void construct_cells_2d_kernel(const Point2D_GPU *sites, int num_sites,
                                          const Point2D_GPU *neighbors, const int *neighbor_counts,
                                          const int *neighbor_offsets, CellResult2D *results, real_t box_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_sites) return;
    
    Corner2D_GPU corners[MAX_CORNERS];
    int num_corners, num_active = 4, first_active = 0;
    
    Point2D_GPU center = sites[idx];
    int num_neigh = neighbor_counts[idx];
    int neigh_start = neighbor_offsets[idx];
    
    init_bounding_box_2d_local(corners, &num_corners, box_size);
    
    for (int n = 0; n < num_neigh && n < MAX_NEIGHBORS && num_corners + 2 < MAX_CORNERS; n++) {
        Point2D_GPU neighbor = neighbors[neigh_start + n];
        neighbor.x -= center.x; neighbor.y -= center.y;
        real_t dist_sq = neighbor.x * neighbor.x + neighbor.y * neighbor.y;
        if (dist_sq < EPSILON) continue;
        real_t wfrac = compute_weight_fraction_gpu(center.w2, neighbor.w2, dist_sq);
        
        if (classify_corner_2d_gpu(&neighbor, corners[first_active].x, corners[first_active].y, wfrac) == CLASSIFY_INSIDE) {
            bool found = false;
            int curr = first_active;
            for (int i = 0; i < num_active; i++) {
                if (classify_corner_2d_gpu(&neighbor, corners[curr].x, corners[curr].y, wfrac) == CLASSIFY_OUTSIDE) {
                    first_active = curr; found = true; break;
                }
                curr = corners[curr].upper_link;
            }
            if (!found) continue;
        }
        
        int lower = corners[first_active].upper_link, cnt = 0;
        while (classify_corner_2d_gpu(&neighbor, corners[lower].x, corners[lower].y, wfrac) == CLASSIFY_OUTSIDE && ++cnt < num_active)
            lower = corners[lower].upper_link;
        int upper = corners[first_active].lower_link; cnt = 0;
        while (classify_corner_2d_gpu(&neighbor, corners[upper].x, corners[upper].y, wfrac) == CLASSIFY_OUTSIDE && ++cnt < num_active)
            upper = corners[upper].lower_link;
        
        int lower_prev = corners[lower].lower_link;
        real_t t1;
        find_edge_intersection_2d_gpu(corners[lower_prev].x, corners[lower_prev].y, corners[lower].x, corners[lower].y, &neighbor, wfrac, &corners[num_corners].x, &corners[num_corners].y, &t1);
        corners[num_corners].status = STATUS_ACTIVE;
        corners[num_corners].lower_link = lower_prev;
        
        int upper_next = corners[upper].upper_link;
        real_t t2;
        find_edge_intersection_2d_gpu(corners[upper].x, corners[upper].y, corners[upper_next].x, corners[upper_next].y, &neighbor, wfrac, &corners[num_corners + 1].x, &corners[num_corners + 1].y, &t2);
        corners[num_corners + 1].status = STATUS_ACTIVE;
        corners[num_corners + 1].upper_link = upper_next;
        
        corners[num_corners].upper_link = num_corners + 1;
        corners[num_corners + 1].lower_link = num_corners;
        corners[lower_prev].upper_link = num_corners;
        corners[upper_next].lower_link = num_corners + 1;
        
        int curr_mark = lower;
        while (curr_mark != upper_next) { corners[curr_mark].status = STATUS_INACTIVE; num_active--; curr_mark = corners[curr_mark].upper_link; }
        num_corners += 2;
        num_active += 2;
        first_active = num_corners - 2;
    }
    
    results[idx].num_corners = num_active;
    results[idx].area = compute_area_local(corners, first_active, num_active);
    compute_centroid_local(corners, first_active, &results[idx].centroid_x, &results[idx].centroid_y);
    results[idx].centroid_x += center.x;
    results[idx].centroid_y += center.y;
}

/*============================================================================
 * HOST API
 *============================================================================*/

static bool cuda_initialized = false;
static int cuda_device_id = 0;

int laguerre_cuda_init(int device_id) {
    if (cuda_initialized) return 0;
    int count;
    CUDA_CHECK(cudaGetDeviceCount(&count));
    if (count == 0) { fprintf(stderr, "No CUDA devices\n"); return -1; }
    if (device_id < 0) {
        size_t max_mem = 0;
        for (int i = 0; i < count; i++) {
            cudaDeviceProp prop;
            CUDA_CHECK(cudaGetDeviceProperties(&prop, i));
            if (prop.totalGlobalMem > max_mem) { max_mem = prop.totalGlobalMem; device_id = i; }
        }
    }
    CUDA_CHECK(cudaSetDevice(device_id));
    cuda_device_id = device_id;
    cuda_initialized = true;
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, device_id));
    printf("CUDA: %s (SM %d.%d, %.2f GB)\n", prop.name, prop.major, prop.minor, prop.totalGlobalMem / 1e9);
    return 0;
}

void laguerre_cuda_cleanup(void) {
    if (cuda_initialized) { cudaDeviceReset(); cuda_initialized = false; }
}

int laguerre_construct_cells_3d(const Point3D_GPU *h_sites, int num_sites, const Point3D_GPU *h_neighbors, const int *h_neighbor_counts, CellResult3D *h_results, real_t box_size) {
    if (!cuda_initialized && laguerre_cuda_init(-1) != 0) return -1;
    
    int total = 0;
    int *h_offsets = (int *)malloc(num_sites * sizeof(int));
    for (int i = 0; i < num_sites; i++) { h_offsets[i] = total; total += h_neighbor_counts[i]; }
    
    Point3D_GPU *d_sites, *d_neighbors;
    int *d_counts, *d_offsets;
    CellResult3D *d_results;
    
    CUDA_CHECK(cudaMalloc(&d_sites, num_sites * sizeof(Point3D_GPU)));
    CUDA_CHECK(cudaMalloc(&d_neighbors, total * sizeof(Point3D_GPU)));
    CUDA_CHECK(cudaMalloc(&d_counts, num_sites * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_offsets, num_sites * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_results, num_sites * sizeof(CellResult3D)));
    
    CUDA_CHECK(cudaMemcpy(d_sites, h_sites, num_sites * sizeof(Point3D_GPU), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_neighbors, h_neighbors, total * sizeof(Point3D_GPU), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_counts, h_neighbor_counts, num_sites * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_offsets, h_offsets, num_sites * sizeof(int), cudaMemcpyHostToDevice));
    
    int blocks = (num_sites + CUDA_BLOCK_SIZE - 1) / CUDA_BLOCK_SIZE;
    construct_cells_3d_kernel<<<blocks, CUDA_BLOCK_SIZE>>>(d_sites, num_sites, d_neighbors, d_counts, d_offsets, d_results, box_size);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    
    CUDA_CHECK(cudaMemcpy(h_results, d_results, num_sites * sizeof(CellResult3D), cudaMemcpyDeviceToHost));
    
    cudaFree(d_sites); cudaFree(d_neighbors); cudaFree(d_counts); cudaFree(d_offsets); cudaFree(d_results);
    free(h_offsets);
    return 0;
}

int laguerre_construct_cells_2d(const Point2D_GPU *h_sites, int num_sites, const Point2D_GPU *h_neighbors, const int *h_neighbor_counts, CellResult2D *h_results, real_t box_size) {
    if (!cuda_initialized && laguerre_cuda_init(-1) != 0) return -1;
    
    int total = 0;
    int *h_offsets = (int *)malloc(num_sites * sizeof(int));
    for (int i = 0; i < num_sites; i++) { h_offsets[i] = total; total += h_neighbor_counts[i]; }
    
    Point2D_GPU *d_sites, *d_neighbors;
    int *d_counts, *d_offsets;
    CellResult2D *d_results;
    
    CUDA_CHECK(cudaMalloc(&d_sites, num_sites * sizeof(Point2D_GPU)));
    CUDA_CHECK(cudaMalloc(&d_neighbors, total * sizeof(Point2D_GPU)));
    CUDA_CHECK(cudaMalloc(&d_counts, num_sites * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_offsets, num_sites * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_results, num_sites * sizeof(CellResult2D)));
    
    CUDA_CHECK(cudaMemcpy(d_sites, h_sites, num_sites * sizeof(Point2D_GPU), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_neighbors, h_neighbors, total * sizeof(Point2D_GPU), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_counts, h_neighbor_counts, num_sites * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_offsets, h_offsets, num_sites * sizeof(int), cudaMemcpyHostToDevice));
    
    int blocks = (num_sites + CUDA_BLOCK_SIZE - 1) / CUDA_BLOCK_SIZE;
    construct_cells_2d_kernel<<<blocks, CUDA_BLOCK_SIZE>>>(d_sites, num_sites, d_neighbors, d_counts, d_offsets, d_results, box_size);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    
    CUDA_CHECK(cudaMemcpy(h_results, d_results, num_sites * sizeof(CellResult2D), cudaMemcpyDeviceToHost));
    
    cudaFree(d_sites); cudaFree(d_neighbors); cudaFree(d_counts); cudaFree(d_offsets); cudaFree(d_results);
    free(h_offsets);
    return 0;
}

void laguerre_cuda_get_device_info(int *major, int *minor, size_t *mem, int *mp) {
    if (!cuda_initialized) { *major = *minor = *mp = 0; *mem = 0; return; }
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, cuda_device_id);
    *major = prop.major; *minor = prop.minor; *mem = prop.totalGlobalMem; *mp = prop.multiProcessorCount;
}
