#ifndef GFS_MATERIAL_STAGE_TIME_H
#define GFS_MATERIAL_STAGE_TIME_H
/* Read-only time of the state currently being evaluated by native material
 * SSPRK2. phase: 1 = first RHS, 2 = predictor RHS, 3 = accepted-candidate
 * endpoint geometry. Return zero outside an evaluation; callers must not
 * silently substitute the step-start time for a missing stage context.
 * Host callbacks receive this before constructing boundary data for CUDA.
 * The existing RK step-start context remains unchanged for rollback/logs.
 */
#ifdef __cplusplus
extern "C" {
#endif
int gfs_material_native_evaluation_time(double *time, int *phase);
#ifdef __cplusplus
}
#endif
#endif
