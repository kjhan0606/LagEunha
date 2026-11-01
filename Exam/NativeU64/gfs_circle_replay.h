#ifndef GFS_CIRCLE_REPLAY_H
#define GFS_CIRCLE_REPLAY_H
#include "gfs_circle_dense.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct { double x,y,vx,vy; } GfsCircleReplayState;
typedef struct {
    GfsCircleReplayState end;
    double kx[4],ky[4]; /* Actual RK4 position increments, including h. */
} GfsCircleReplayTrial;
/* Every call must reconstruct the trial from the same immutable accepted
 * start. The callback must not publish state, ledgers, or an accepted time.
 * Return 0 on a valid native trial, nonzero on geometry/RHS/state failure.
 */
typedef int (*GfsCircleReplayCallback)(void *,double,GfsCircleReplayTrial *);
typedef struct {
    double cx,cy,radius,radial_tolerance,minimum_normal_speed;
    int max_trials;
    GfsCircleDenseOptions dense;
} GfsCircleReplayOptions;
enum GfsCircleReplayStatus { GFS_CIRCLE_REPLAY_UNRESOLVED=-1,
                            GFS_CIRCLE_REPLAY_CONTACT=1 };
enum GfsCircleReplayReason {
    GFS_CIRCLE_REPLAY_REFINED=0, GFS_CIRCLE_REPLAY_INITIAL_CONTACT=1,
    GFS_CIRCLE_REPLAY_INVALID=2, GFS_CIRCLE_REPLAY_OUTSIDE_START=3,
    GFS_CIRCLE_REPLAY_TRIAL_FAILED=4, GFS_CIRCLE_REPLAY_NO_OUTSIDE_BRACKET=5,
    GFS_CIRCLE_REPLAY_DENSE_UNRESOLVED=6, GFS_CIRCLE_REPLAY_EARLIER_CANDIDATE=7,
    GFS_CIRCLE_REPLAY_GRAZING=8, GFS_CIRCLE_REPLAY_TRIAL_LIMIT=9,
    GFS_CIRCLE_REPLAY_TIME_RESOLUTION=10, GFS_CIRCLE_REPLAY_INCONSISTENT_RK=11
};
typedef struct {
    int status,reason,trials;
    double h,lo,hi,glo,ghi,gap,normal_speed;
    GfsCircleReplayState contact;
} GfsCircleReplayResult;

/* Refine an actual numerical-RK inside/outside bracket. This does not prove
 * accuracy of the physical trajectory, apply an impulse, or commit a step.
 * A zero-time initial contact is explicit; the driver must ledger its impulse
 * before retrying a positive advance. Uncertainty means shrink/retry, never
 * accept an unchecked event or ignore an earlier dense candidate.
 */
GfsCircleReplayResult gfs_circle_replay_refine(GfsCircleReplayState start,
    double outside_h,GfsCircleReplayOptions options,
    GfsCircleReplayCallback callback,void *context);
#ifdef __cplusplus
}
#endif
#endif
