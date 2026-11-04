/* Event-resolved specular impacts for constant-velocity KDK drift segments.
 * The domain is the same regular fixed polygon used by the CUDA tessellator.
 * Corners are deliberately rejected: assigning a single normal there would
 * make the result depend on face-search order.
 */
#ifndef GFS_KDK_FIXEDPOLY_DRIFT_INC
#define GFS_KDK_FIXEDPOLY_DRIFT_INC
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "kepler_slip_wall.h"
#include "gfs_fixed_polygon_geometry.h"

typedef struct {
    unsigned long long impacts;
    double wall_jx, wall_jy, wall_torque, kinetic_roundoff;
    int failure_code, failure_side;
    double failure_remaining, failure_gap;
} GfsKdkFixedPolyWallStats;

static double gfs_kdk_wall_nx[GFS_FIXEDPOLY_SIDES];
static double gfs_kdk_wall_ny[GFS_FIXEDPOLY_SIDES];
static int gfs_kdk_wall_geometry_ready = 0;

static int gfs_kdk_fixedpoly_wall_events(void)
{
    const char *s = getenv("GFS_KEPLER_KDK_FIXEDPOLY_WALL_EVENTS");
    return s && strcmp(s, "1") == 0;
}

static void gfs_kdk_fixedpoly_wall_init(void)
{
    if (gfs_kdk_wall_geometry_ready) return;
    const double step = 2.0 * M_PI / (double)GFS_FIXEDPOLY_SIDES;
    for (int k = 0; k < GFS_FIXEDPOLY_SIDES; ++k) {
        const double angle = ((double)k + 0.5) * step;
        gfs_kdk_wall_nx[k] = cos(angle);
        gfs_kdk_wall_ny[k] = sin(angle);
    }
    gfs_kdk_wall_geometry_ready = 1;
}

static double gfs_kdk_fixedpoly_gap(double x, double y, int *side)
{
    const double step = 2.0 * M_PI / (double)GFS_FIXEDPOLY_SIDES;
    const double apothem = 2.75 * cos(M_PI / (double)GFS_FIXEDPOLY_SIDES);
    const double angle = atan2(y - 3.0, x - 3.0);
    int k = (int)floor(angle / step);
    k %= GFS_FIXEDPOLY_SIDES;
    if (k < 0) k += GFS_FIXEDPOLY_SIDES;
    if (side) *side = k;
    return gfs_kdk_wall_nx[k] * (x - 3.0) +
           gfs_kdk_wall_ny[k] * (y - 3.0) - apothem;
}

/* Return 1 on success, 0 on invalid/out-of-domain/corner trajectories. */
static int gfs_kdk_fixedpoly_wall_drift(
        double *x, double *y, double *vx, double *vy,
        double elapsed, double mass, GfsKdkFixedPolyWallStats *stats)
{
    if (stats) memset(stats, 0, sizeof(*stats));
    if (!x || !y || !vx || !vy || !stats)
        return 0;
    if (!isfinite(*x) || !isfinite(*y) || !isfinite(*vx) ||
        !isfinite(*vy) || !isfinite(elapsed) || elapsed < 0.0 ||
        !isfinite(mass) || !(mass > 0.0)) {
        stats->failure_code = 1;
        return 0;
    }

    gfs_kdk_fixedpoly_wall_init();
    const double apothem = 2.75 * cos(M_PI / (double)GFS_FIXEDPOLY_SIDES);
    double remaining = elapsed;
    double tol = 32.0 * DBL_EPSILON *
        fmax(2.75, fabs(*x - 3.0) + fabs(*y - 3.0));
#define GFS_KDK_WALL_FAIL(CODE,SIDE,GAP) do { \
        stats->failure_code=(CODE); stats->failure_side=(SIDE); \
        stats->failure_remaining=remaining; stats->failure_gap=(GAP); \
        return 0; \
    } while (0)
    int start_side = -1;
    const double start_gap = gfs_kdk_fixedpoly_gap(*x, *y, &start_side);
    if (start_gap > tol) GFS_KDK_WALL_FAIL(2,start_side,start_gap);

    for (int bounce = 0; remaining > 0.0; ++bounce) {
        if (bounce >= 16) GFS_KDK_WALL_FAIL(11,-1,remaining);
        const double time_tol = 64.0 * DBL_EPSILON *
            fmax(1.0, fabs(remaining));
        const double xe = *x + remaining * *vx;
        const double ye = *y + remaining * *vy;
        if (!isfinite(xe) || !isfinite(ye)) GFS_KDK_WALL_FAIL(10,-1,NAN);
        const double end_gap = gfs_kdk_fixedpoly_gap(xe, ye, NULL);

        /* Convexity makes an interior endpoint a complete no-crossing test.
         * Near contact, still solve the face event so an endpoint impact gets
         * its impulse in this drift rather than in the following step.
         */
        if (end_gap < -tol) {
            *x = xe;
            *y = ye;
            remaining = 0.0;
            break;
        }

        int first = -1;
        double first_t = DBL_MAX;
        double first_vn = 0.0;
        for (int k = 0; k < GFS_FIXEDPOLY_SIDES; ++k) {
            const double rx = *x - 3.0, ry = *y - 3.0;
            const double gap0 = gfs_kdk_wall_nx[k] * rx +
                gfs_kdk_wall_ny[k] * ry - apothem;
            if (gap0 > tol) GFS_KDK_WALL_FAIL(3,k,gap0);
            const double vn = gfs_kdk_wall_nx[k] * *vx +
                gfs_kdk_wall_ny[k] * *vy;
            if (!(vn > 0.0)) continue;
            double thit = -gap0 / vn;
            if (thit < -time_tol) continue;
            if (thit < 0.0) thit = 0.0;
            if (thit <= remaining + time_tol && thit < first_t) {
                first = k;
                first_t = thit;
                first_vn = vn;
            }
        }
        if (first < 0) {
            if (end_gap <= tol) {
                *x = xe;
                *y = ye;
                remaining = 0.0;
                break;
            }
            GFS_KDK_WALL_FAIL(4,-1,end_gap);
        }
        if (first_t > remaining) first_t = remaining;

        double hit_x = *x + first_t * *vx;
        double hit_y = *y + first_t * *vy;
        const double hit_gap = gfs_kdk_wall_nx[first] * (hit_x - 3.0) +
            gfs_kdk_wall_ny[first] * (hit_y - 3.0) - apothem;
        hit_x -= hit_gap * gfs_kdk_wall_nx[first];
        hit_y -= hit_gap * gfs_kdk_wall_ny[first];

        /* Reject a vertex contact, and any inconsistent face ordering. */
        const double hit_tol = 32.0 * DBL_EPSILON *
            fmax(2.75, fabs(hit_x - 3.0) + fabs(hit_y - 3.0));
        for (int k = 0; k < GFS_FIXEDPOLY_SIDES; ++k) {
            if (k == first) continue;
            const double other_gap = gfs_kdk_wall_nx[k] * (hit_x - 3.0) +
                gfs_kdk_wall_ny[k] * (hit_y - 3.0) - apothem;
            if (other_gap > hit_tol || fabs(other_gap) <= hit_tol)
                GFS_KDK_WALL_FAIL(5,k,other_gap);
        }

        GfsSlipImpulse impulse;
        if (!(first_vn > 0.0) ||
            gfs_slip_plane_impact(hit_x, hit_y, 3.0, 3.0,
                gfs_kdk_wall_nx[first], gfs_kdk_wall_ny[first],
                apothem, mass, *vx, *vy, &impulse) != GFS_SLIP_REFLECTED)
            GFS_KDK_WALL_FAIL(6,first,first_vn);
        *x = hit_x;
        *y = hit_y;
        *vx = impulse.vx_after;
        *vy = impulse.vy_after;
        stats->wall_jx += impulse.wall_jx;
        stats->wall_jy += impulse.wall_jy;
        stats->wall_torque += impulse.wall_torque_impulse;
        stats->kinetic_roundoff += impulse.fluid_kinetic_roundoff;
        ++stats->impacts;
        if (!isfinite(stats->wall_jx) || !isfinite(stats->wall_jy) ||
            !isfinite(stats->wall_torque) ||
            !isfinite(stats->kinetic_roundoff))
            GFS_KDK_WALL_FAIL(7,first,stats->kinetic_roundoff);

        remaining -= first_t;
        if (remaining < 0.0 && remaining >= -time_tol) remaining = 0.0;
        if (remaining < 0.0) GFS_KDK_WALL_FAIL(8,first,remaining);
    }

    tol = 32.0 * DBL_EPSILON *
        fmax(2.75, fabs(*x - 3.0) + fabs(*y - 3.0));
    const double final_gap = gfs_kdk_fixedpoly_gap(*x, *y, NULL);
    if (!isfinite(*x) || !isfinite(*y) || !isfinite(*vx) ||
        !isfinite(*vy) || !(final_gap <= tol))
        GFS_KDK_WALL_FAIL(9,-1,final_gap);
#undef GFS_KDK_WALL_FAIL
    return 1;
}
#endif
