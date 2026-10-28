#ifndef EXAM_MOMENT_GEOMETRY_V2_H
#define EXAM_MOMENT_GEOMETRY_V2_H
#include <stdint.h>
#include <stddef.h>
#define MG2_CAP 128
/* A vertex is an intersection of its two support planes, not the rounded x/y.
 * Plane IDs >=0 are stable registry IDs. -1/-2/-3/-4 are box bottom/right/top/left.
 * x/y, midpoint and second integrals are relative to the cell generator.
 * edge_dx/dy are rounded EXACT differences, never differences of rounded x/y.
 * geometry_wn is evaluated at the midpoint and is NOT the hydro boundary BC.
 * endpoint_rounding_difference includes ordinary last-bit subtraction error;
 * endpoint_collapse/reversal and max_relative_vector_error separately expose
 * actual failure of an endpoint-only representation. */
typedef struct {
 uint64_t epoch;
 int count,endpoint_rounding_difference,endpoint_collapse,endpoint_reversal;
 double max_relative_vector_error;
 uint64_t filtered_checks,exact_checks,added_candidates;
 double area,centroid[2],second[4];
 double x[MG2_CAP],y[MG2_CAP];
 int64_t support_a[MG2_CAP],support_b[MG2_CAP],neighbor[MG2_CAP];
 int neighbor_kind[MG2_CAP];
 double midpoint_x[MG2_CAP],midpoint_y[MG2_CAP];
 double edge_dx[MG2_CAP],edge_dy[MG2_CAP];
 double length[MG2_CAP],length_lo[MG2_CAP],length_hi[MG2_CAP];
 double nx[MG2_CAP],ny[MG2_CAP],geometry_wn[MG2_CAP];
} MomentGeometryV2;
#ifdef __cplusplus
extern "C" {
#endif
/* Explicit HOST native exact-resolution stage. Not a GPU geometry claim.
 * Input candidate CSR is in stable global ID order. Whole-registry completion
 * is mandatory in this audit API. statuses 1 input,2 cap,4 empty,5 box,6 geometry,
 * 8 allocation,9 shared topology,10 output representation/overflow,11 audit I/O.
 * Outputs transactional; status is diagnostic and may change on failure.
 * audit_path receives exact global rational endpoints and candidate witnesses. */
int moment_geometry_v2_resolve(int threads,int n,int registry_count,
 const int64_t *ids,const double *xy,const double *W,const double *velocity,
 const double *Wdot,const int *kind,const int64_t *selected,const int64_t *offsets,
 const int64_t *candidates,const double *bounds,uint64_t epoch,
 MomentGeometryV2 *output,int *status,const char *audit_path);
size_t moment_geometry_v2_size(void);
/* Predicate audit interfaces; they do not generate topology on the GPU.
 * Filter mode0 support rows[n,4], mode1 owner[n] and enclosures[n,4].
 * mode1 returns n*registry_count signs in vertex-major Cartesian order.
 * sign0 is unresolved. exact_queries sign2 denotes singular supports. */
int moment_geometry_v2_filter(int backend,int threads,int mode,int registry_count,
 const double *xy,const double *W,const double *bounds,int n,
 const int64_t *queries_or_owners,const double *enclosures,int8_t *signs);
int moment_geometry_v2_exact_queries(int registry_count,const double *xy,
 const double *W,const double *bounds,int n,const int64_t *queries,int *signs);
#ifdef __cplusplus
}
#endif
#endif
