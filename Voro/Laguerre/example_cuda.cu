/**
 * @file example_cuda.cu
 * @brief Example usage of the Laguerre diagram CUDA library
 * 
 * Compile: nvcc -O3 example_cuda.cu laguerre_cuda.cu -o example_cuda
 * Run: ./example_cuda
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime.h>
#include "laguerre_cuda.cuh"

#define NUM_SITES       10000
#define NUM_NEIGHBORS   50
#define BOX_SIZE        100.0

void generateRandomPoints3D(Point3D_GPU *pts, int n, real_t box) {
    real_t half = box / 2.0;
    for (int i = 0; i < n; i++) {
        pts[i].x = ((real_t)rand() / RAND_MAX) * box - half;
        pts[i].y = ((real_t)rand() / RAND_MAX) * box - half;
        pts[i].z = ((real_t)rand() / RAND_MAX) * box - half;
        pts[i].w2 = ((real_t)rand() / RAND_MAX) * 0.5;
        pts[i].dist2 = 0;
        pts[i].index = i;
    }
}

void generateRandomPoints2D(Point2D_GPU *pts, int n, real_t box) {
    real_t half = box / 2.0;
    for (int i = 0; i < n; i++) {
        pts[i].x = ((real_t)rand() / RAND_MAX) * box - half;
        pts[i].y = ((real_t)rand() / RAND_MAX) * box - half;
        pts[i].w2 = ((real_t)rand() / RAND_MAX) * 0.5;
        pts[i].dist2 = 0;
        pts[i].index = i;
    }
}

int compareByDist2(const void *a, const void *b) {
    const Point3D_GPU *pa = (const Point3D_GPU *)a;
    const Point3D_GPU *pb = (const Point3D_GPU *)b;
    return (pa->dist2 < pb->dist2) ? -1 : (pa->dist2 > pb->dist2) ? 1 : 0;
}

void prepareNeighborLists(const Point3D_GPU *sites, int ns, Point3D_GPU *neighbors, int *counts, int maxN) {
    for (int i = 0; i < ns; i++) {
        int cnt = 0;
        for (int j = 0; j < ns && cnt < maxN; j++) {
            if (i == j) continue;
            int off = i * maxN + cnt;
            neighbors[off] = sites[j];
            real_t dx = sites[j].x - sites[i].x;
            real_t dy = sites[j].y - sites[i].y;
            real_t dz = sites[j].z - sites[i].z;
            neighbors[off].dist2 = dx*dx + dy*dy + dz*dz;
            cnt++;
        }
        qsort(&neighbors[i * maxN], cnt, sizeof(Point3D_GPU), compareByDist2);
        counts[i] = cnt;
    }
}

int main(void) {
    printf("Laguerre Diagram Example (CUDA)\n");
    printf("===============================\n\n");
    
    if (laguerre_cuda_init(-1) != 0) { 
        fprintf(stderr, "CUDA initialization failed\n"); 
        return 1; 
    }
    
    int major, minor, mp; 
    size_t mem;
    laguerre_cuda_get_device_info(&major, &minor, &mem, &mp);
    printf("Device: SM %d.%d, %.2f GB, %d MPs\n\n", major, minor, mem/1e9, mp);
    
    srand((unsigned)time(NULL));
    
    /* 3D Test */
    printf("3D Test: %d sites, %d neighbors\n", NUM_SITES, NUM_NEIGHBORS);
    
    Point3D_GPU *sites = (Point3D_GPU *)malloc(NUM_SITES * sizeof(Point3D_GPU));
    Point3D_GPU *neighbors = (Point3D_GPU *)malloc(NUM_SITES * NUM_NEIGHBORS * sizeof(Point3D_GPU));
    int *counts = (int *)malloc(NUM_SITES * sizeof(int));
    CellResult3D *results = (CellResult3D *)malloc(NUM_SITES * sizeof(CellResult3D));
    
    generateRandomPoints3D(sites, NUM_SITES, BOX_SIZE);
    prepareNeighborLists(sites, NUM_SITES, neighbors, counts, NUM_NEIGHBORS);
    
    cudaEvent_t t1, t2;
    cudaEventCreate(&t1); 
    cudaEventCreate(&t2);
    cudaEventRecord(t1);
    
    laguerre_construct_cells_3d(sites, NUM_SITES, neighbors, counts, results, BOX_SIZE);
    
    cudaEventRecord(t2);
    cudaEventSynchronize(t2);
    float ms;
    cudaEventElapsedTime(&ms, t1, t2);
    
    real_t vol = 0;
    int valid = 0;
    for (int i = 0; i < NUM_SITES; i++) {
        if (results[i].num_vertices >= 4) { 
            vol += results[i].volume; 
            valid++; 
        }
    }
    
    printf("  Valid: %d, Volume: %.2f (expected: %.2f)\n", valid, vol, BOX_SIZE*BOX_SIZE*BOX_SIZE);
    printf("  Time: %.2f ms (%.0f cells/s)\n\n", ms, NUM_SITES/(ms/1000));
    
    /* 2D Test */
    printf("2D Test: %d sites\n", NUM_SITES);
    
    Point2D_GPU *sites2 = (Point2D_GPU *)malloc(NUM_SITES * sizeof(Point2D_GPU));
    Point2D_GPU *neighbors2 = (Point2D_GPU *)malloc(NUM_SITES * NUM_NEIGHBORS * sizeof(Point2D_GPU));
    CellResult2D *results2 = (CellResult2D *)malloc(NUM_SITES * sizeof(CellResult2D));
    
    generateRandomPoints2D(sites2, NUM_SITES, BOX_SIZE);
    for (int i = 0; i < NUM_SITES; i++) {
        int cnt = 0;
        for (int j = 0; j < NUM_SITES && cnt < NUM_NEIGHBORS; j++) {
            if (i != j) neighbors2[i * NUM_NEIGHBORS + cnt++] = sites2[j];
        }
        counts[i] = cnt;
    }
    
    cudaEventRecord(t1);
    laguerre_construct_cells_2d(sites2, NUM_SITES, neighbors2, counts, results2, BOX_SIZE);
    cudaEventRecord(t2);
    cudaEventSynchronize(t2);
    cudaEventElapsedTime(&ms, t1, t2);
    
    real_t area = 0;
    valid = 0;
    for (int i = 0; i < NUM_SITES; i++) {
        if (results2[i].num_corners >= 3) { 
            area += results2[i].area; 
            valid++; 
        }
    }
    
    printf("  Valid: %d, Area: %.2f (expected: %.2f)\n", valid, area, BOX_SIZE*BOX_SIZE);
    printf("  Time: %.2f ms (%.0f cells/s)\n", ms, NUM_SITES/(ms/1000));
    
    cudaEventDestroy(t1); 
    cudaEventDestroy(t2);
    free(sites); 
    free(neighbors); 
    free(counts); 
    free(results);
    free(sites2); 
    free(neighbors2); 
    free(results2);
    laguerre_cuda_cleanup();
    
    printf("\nDone!\n");
    return 0;
}
