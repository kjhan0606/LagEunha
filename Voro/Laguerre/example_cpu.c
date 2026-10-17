/**
 * @file example_cpu.c
 * @brief Example usage of the Laguerre diagram CPU library
 * 
 * Compile: gcc -O3 example_cpu.c laguerre.o -o example_cpu -lm
 * Run: ./example_cpu
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define LAGUERRE_MAIN
#include "laguerre.h"
#undef LAGUERRE_MAIN

#define NUM_SITES       100
#define NUM_NEIGHBORS   30
#define BOX_SIZE        10.0
#define MAX_VERTS       256

void generate_random_points(Point3D *pts, int n, real_t box) {
    real_t half = box / 2.0;
    for (int i = 0; i < n; i++) {
        pts[i].x = ((real_t)rand() / RAND_MAX) * box - half;
        pts[i].y = ((real_t)rand() / RAND_MAX) * box - half;
        pts[i].z = ((real_t)rand() / RAND_MAX) * box - half;
        pts[i].weightSquared = ((real_t)rand() / RAND_MAX) * 0.5;
        pts[i].index = i;
        pts[i].squaredDist = 0;
        pts[i].soundSpeed = 1.0;
    }
}

void compute_distances(Point3D *pts, int n, Point3D *ref) {
    for (int i = 0; i < n; i++) {
        real_t dx = pts[i].x - ref->x;
        real_t dy = pts[i].y - ref->y;
        real_t dz = pts[i].z - ref->z;
        pts[i].squaredDist = dx*dx + dy*dy + dz*dz;
    }
}

int main(void) {
    printf("Laguerre Diagram Example (CPU)\n");
    printf("==============================\n\n");
    
    srand((unsigned)time(NULL));
    
    Point3D *sites = malloc(NUM_SITES * sizeof(Point3D));
    Point3D *neighbors = malloc(NUM_SITES * sizeof(Point3D));
    Vertex3D *vertices = malloc(MAX_VERTS * sizeof(Vertex3D));
    int *work = malloc(MAX_VERTS * sizeof(int));
    
    printf("Generating %d random points...\n", NUM_SITES);
    generate_random_points(sites, NUM_SITES, BOX_SIZE);
    
    printf("Constructing cells...\n\n");
    
    real_t total_vol = 0;
    int degen = 0;
    clock_t start = clock();
    
    for (int i = 0; i < NUM_SITES; i++) {
        Point3D *center = &sites[i];
        
        int nc = 0;
        for (int j = 0; j < NUM_SITES && nc < NUM_NEIGHBORS; j++) {
            if (i == j) continue;
            neighbors[nc++] = sites[j];
        }
        
        compute_distances(neighbors, nc, center);
        qsort(neighbors, nc, sizeof(Point3D), compareByDistanceSquared3D);
        
        for (int j = 0; j < nc; j++) {
            neighbors[j].x -= center->x;
            neighbors[j].y -= center->y;
            neighbors[j].z -= center->z;
        }
        
        int nv = initBoundingBox3D(vertices, MAX_VERTS, BOX_SIZE);
        
        for (int n = 0; n < nc; n++) {
            Point3D *nb = &neighbors[n];
            real_t d2 = dot3D(nb, nb);
            if (d2 < EPSILON) continue;
            
            real_t wf = computeWeightFraction(center->weightSquared, nb->weightSquared, d2);
            
            for (int v = 0; v < nv; v++) {
                if (vertices[v].status != STATUS_ACTIVE) continue;
                real_t off;
                if (classifyVertexPosition3D(nb, &vertices[v], wf, &off) == POINT_OUTSIDE)
                    vertices[v].status = STATUS_INACTIVE;
            }
            
            real_t md2;
            nv = compactActiveVertices3D(vertices, nv, &md2, work);
            if (nv < 4) { degen++; break; }
        }
        
        if (nv >= 4) {
            real_t vol = computePolyhedronVolume3D(vertices, nv);
            total_vol += vol;
            if (i < 5) printf("  Site %d: %d vertices, vol=%.4f\n", i, nv, vol);
        }
    }
    
    double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
    
    printf("\nResults:\n");
    printf("  Total cells: %d, Degenerate: %d\n", NUM_SITES, degen);
    printf("  Total volume: %.4f (expected: %.4f)\n", total_vol, BOX_SIZE*BOX_SIZE*BOX_SIZE);
    printf("  Time: %.3f s (%.1f cells/s)\n", elapsed, NUM_SITES/elapsed);
    
    free(sites); free(neighbors); free(vertices); free(work);
    printf("\nDone!\n");
    return 0;
}
