#ifndef GFS_CIRCLE_OBSERVER_H
#define GFS_CIRCLE_OBSERVER_H
#include <stddef.h>
#include <stdint.h>
#include "gfs_circle_replay.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint64_t id;
    GfsCircleReplayState start;
    GfsCircleReplayTrial trial;
} GfsCircleObservation;
/* Native caller must preserve the complete accepted ID/mass/payload set and
 * provide the same immutable accepted start through every retry. Its existing
 * identity/rollback validation remains authoritative. Initial scan reads only
 * id/start. Completed-trial scan reads the actual native increments/endpoint.
 */
typedef int (*GfsCircleObservationView)(void *,size_t,GfsCircleObservation *);
typedef struct {
    double cx,cy,radius,minimum_normal_speed;
    int max_trials,max_dense_nodes,max_dense_depth;
} GfsCircleObserverOptions;
enum GfsCircleObserverAction {
    GFS_CIRCLE_OBSERVER_STOP=-1, GFS_CIRCLE_OBSERVER_RETRY=0,
    GFS_CIRCLE_OBSERVER_QUIET=1, GFS_CIRCLE_OBSERVER_CONTACT=2,
    GFS_CIRCLE_OBSERVER_ZERO=3
};
enum GfsCircleObserverReason {
    GFS_CIRCLE_OBSERVER_INITIAL=0, GFS_CIRCLE_OBSERVER_NO_CANDIDATE=1,
    GFS_CIRCLE_OBSERVER_BRACKET=2, GFS_CIRCLE_OBSERVER_REFINED=3,
    GFS_CIRCLE_OBSERVER_OUT_AND_BACK=4, GFS_CIRCLE_OBSERVER_INWARD_START=5,
    GFS_CIRCLE_OBSERVER_UNCERTAIN=6, GFS_CIRCLE_OBSERVER_OVERLAP=7,
    GFS_CIRCLE_OBSERVER_GRAZING=8, GFS_CIRCLE_OBSERVER_INVALID=9,
    GFS_CIRCLE_OBSERVER_LIMIT=10, GFS_CIRCLE_OBSERVER_FAILED_STAGE=11,
    GFS_CIRCLE_OBSERVER_BASE_CHANGED=12, GFS_CIRCLE_OBSERVER_TIME_LIMIT=13
};
typedef struct {
    GfsCircleObserverOptions options;
    uint64_t epoch,selected_id;
    size_t count,zero_contacts;
    int action,reason,initial_checked,refining,trials,stalled,bisections,switches;
    double ceiling,lo,hi,glo,ghi;
    GfsCircleReplayState selected_start;
} GfsCircleObserver;

double gfs_circle_observer_tolerance(GfsCircleReplayState,GfsCircleObserverOptions);
void gfs_circle_observer_init(GfsCircleObserver *,GfsCircleObserverOptions,
    size_t count,uint64_t accepted_epoch,double initial_ceiling);
int gfs_circle_observer_initial(GfsCircleObserver *,GfsCircleObservationView,void *);
/* Every completed trial must be observed, even with all endpoints inside.
 * actual_h can be smaller than ceiling because of native CFL limits.
 * CONTACT/QUIET are requests for the caller, not state/ledger commits.
 */
int gfs_circle_observer_observe(GfsCircleObserver *,uint64_t accepted_epoch,
    double actual_h,GfsCircleObservationView,void *);
/* Failed geometry/RHS stages are never bracket observations. */
int gfs_circle_observer_reject_stage(GfsCircleObserver *,uint64_t accepted_epoch,
    double actual_h);

/* ZERO is separate from positive-time retry. The native adapter must validate
 * the whole initial batch, atomically ledger/apply its impulses once, publish
 * a restart-distinguishable post-impulse epoch, then reinitialize the observer
 * from the new baseline. This pure observer applies no impulse or mutation.
 */
#ifdef __cplusplus
}
#endif
#endif
