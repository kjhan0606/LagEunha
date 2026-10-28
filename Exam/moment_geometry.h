#ifndef EXAM_MOMENT_GEOMETRY_H
#define EXAM_MOMENT_GEOMETRY_H
#include <stdint.h>
#include <stddef.h>
#define MOMENT_GEOMETRY_CORNERS 128
/* Separate double sidecar; no change to legacy ParticleSoA/PolygonOut bytes.
 * kind: 0 active moment cell, 1 physical frozen reservoir, 2 exterior registry.
 * MPI ownership is deliberately not encoded in kind. */
typedef struct {
 uint64_t epoch;
 int count;
 double area,centroid_relative[2],second_relative[4];
 double x[MOMENT_GEOMETRY_CORNERS],y[MOMENT_GEOMETRY_CORNERS];
 int64_t neighbor[MOMENT_GEOMETRY_CORNERS];
 int neighbor_kind[MOMENT_GEOMETRY_CORNERS];
 double nx[MOMENT_GEOMETRY_CORNERS],ny[MOMENT_GEOMETRY_CORNERS];
 double length[MOMENT_GEOMETRY_CORNERS],wn_midpoint[MOMENT_GEOMETRY_CORNERS];
} MomentGeometryPolygon;
#ifdef __cplusplus
extern "C" {
#endif
/* Host arrays in/out. backend=0 CPU/OpenMP,1 CUDA. Strict geometry only.
 * candidates are registry row indices in CSR, not stable global IDs.
 * Registry completeness is an independent required caller validation.
 * All polygons commit together after every cell succeeds; status is separate.
 * status: 1 invalid input,2 capacity,3 clipping,4 empty,5 box boundary,
 * 6 invalid geometry/ownership,7 CUDA failure,8 allocation failure.
 * Output vertices are relative to that cell's generator; physical quantities
 * are geometric integrals (density not included). wn_midpoint is geometric plane speed
 * evaluated at the edge midpoint; reconstruct the affine face speed at actual
 * quadrature points rather than treating this as constant. It is
 * not the physical wall boundary condition used by the research provider. */
int gpu_moment_geometry_export(int backend,int threads,int n,int registry_count,
 const int64_t *ids,const double *xy,const double *weights,const double *velocity,
 const double *weight_rates,const int *kind,const int64_t *selected,
 const int64_t *offsets,const int64_t *candidates,const double *bounds,
 uint64_t epoch,MomentGeometryPolygon *output,int *status);
size_t gpu_moment_geometry_polygon_size(void);
#ifdef __cplusplus
}
#endif
#endif
