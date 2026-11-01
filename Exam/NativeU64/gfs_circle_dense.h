#ifndef GFS_CIRCLE_DENSE_H
#define GFS_CIRCLE_DENSE_H
#include "gfs_wall_dense.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Candidate search on the supplied numerical cubic only. A CANDIDATE is
 * not a verified physical event: re-run the native RK trial and reproduce
 * the bracket/contact before applying any impulse. This is not fourth-order
 * event certification. No particle state is advanced or changed here.
 * On x86, the implementation scopes gradual underflow to the predicate and
 * restores the caller's complete MXCSR, including FTZ/DAZ and status flags.
 * Compile the implementation separately with strict floating-point flags.
 */
enum GfsCircleDenseStatus {
    GFS_CIRCLE_DENSE_UNRESOLVED=-1,
    GFS_CIRCLE_DENSE_NONE=0,
    GFS_CIRCLE_DENSE_CANDIDATE=1
};
enum GfsCircleDenseReason {
    GFS_CIRCLE_DENSE_EXCLUDED=0,
    GFS_CIRCLE_DENSE_START_CONTACT=1,
    GFS_CIRCLE_DENSE_SIGN_BRACKET=2,
    GFS_CIRCLE_DENSE_CONTACT_WINDOW=3,
    GFS_CIRCLE_DENSE_INVALID_INPUT=4,
    GFS_CIRCLE_DENSE_OUTSIDE_START=5,
    GFS_CIRCLE_DENSE_UNCERTAIN_START=6,
    GFS_CIRCLE_DENSE_NODE_LIMIT=7,
    GFS_CIRCLE_DENSE_DEPTH_LIMIT=8,
    GFS_CIRCLE_DENSE_ROUNDING_LIMIT=9,
    GFS_CIRCLE_DENSE_FP_ENVIRONMENT=10
};
typedef struct {
    double radius,step,radial_tolerance,minimum_normal_speed;
    int max_nodes,max_depth;
} GfsCircleDenseOptions;
typedef struct {
    int status,reason,nodes,max_depth,refinement_evaluations;
    int dense_sign_bracket_certified;
    double theta_lo,theta_hi;
} GfsCircleDenseResult;

GfsCircleDenseResult gfs_circle_dense_first_candidate(
    GfsWallCubic x,GfsWallCubic y,GfsCircleDenseOptions options);

#ifdef __cplusplus
}
#endif
#endif
