#ifndef GFS_KDK_LTS_REGISTER_H
#define GFS_KDK_LTS_REGISTER_H

/*
 * Extensive conservative transfer storage for asynchronous KDK.
 *
 * A register stores increments to total momentum and total (internal+kinetic)
 * energy. Mass is deliberately absent: the initial W=0 method is fixed-mass.
 * One canonical face evaluation adds +transfer to one endpoint and the exact
 * negative to the other. Applying a register reconstructs internal energy
 * after the momentum change, so the kinetic-energy change is included exactly.
 *
 * This is a solver-independent primitive. The current GPU force kernel does
 * not yet emit canonical face transfers; this header does not activate LTS.
 */
#include <math.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double dpx;
    double dpy;
    double dEtot;
} GfsKdkLtsRegister;

typedef struct {
    double mass;
    double vx;
    double vy;
    double U;       /* extensive internal energy */
} GfsKdkLtsConservedCell;

typedef struct {
    double x;
    double y;
    double vx;
    double vy;
    double U;
} GfsKdkLtsPredictedState;

static inline void gfs_kdk_lts_register_clear(GfsKdkLtsRegister *r)
{
    if (!r) return;
    r->dpx = 0.0;
    r->dpy = 0.0;
    r->dEtot = 0.0;
}

/*
 * Add one time-integrated, oriented face transfer. dEtot_i is the change in
 * endpoint i's total energy (including kinetic); endpoint j receives its
 * exact negative. Returns 0 without mutation for invalid arguments/values.
 */
static inline int gfs_kdk_lts_register_pair(
    GfsKdkLtsRegister *reg, size_t ncell, size_t i, size_t j,
    double dpx_i, double dpy_i, double dEtot_i)
{
    if (!reg || i >= ncell || j >= ncell || i == j ||
        !isfinite(dpx_i) || !isfinite(dpy_i) || !isfinite(dEtot_i))
        return 0;
    const double ix = reg[i].dpx + dpx_i;
    const double iy = reg[i].dpy + dpy_i;
    const double ie = reg[i].dEtot + dEtot_i;
    const double jx = reg[j].dpx - dpx_i;
    const double jy = reg[j].dpy - dpy_i;
    const double je = reg[j].dEtot - dEtot_i;
    if (!isfinite(ix) || !isfinite(iy) || !isfinite(ie) ||
        !isfinite(jx) || !isfinite(jy) || !isfinite(je))
        return 0;
    reg[i].dpx = ix;
    reg[i].dpy = iy;
    reg[i].dEtot = ie;
    reg[j].dpx = jx;
    reg[j].dpy = jy;
    reg[j].dEtot = je;
    return 1;
}

/*
 * Accumulate one cell's integrated RHS from the current variational operator.
 * dEtot_rate is the cell contribution v·F + Udot, multiplied by dt here;
 * only its sum over all cells is constrained by the global work identity.
 * This deliberately does not impose pairwise momentum antisymmetry.
 */
static inline int gfs_kdk_lts_register_cell_rhs(
    GfsKdkLtsRegister *reg, size_t ncell, size_t i,
    double force_x, double force_y, double dEtot_rate, double dt)
{
    if (!reg || i >= ncell || !isfinite(force_x) || !isfinite(force_y) ||
        !isfinite(dEtot_rate) || !isfinite(dt) || dt < 0.0)
        return 0;
    const double px = reg[i].dpx + force_x*dt;
    const double py = reg[i].dpy + force_y*dt;
    const double e = reg[i].dEtot + dEtot_rate*dt;
    if (!isfinite(px) || !isfinite(py) || !isfinite(e)) return 0;
    reg[i].dpx = px;
    reg[i].dpy = py;
    reg[i].dEtot = e;
    return 1;
}

/*
 * Consume a cell's accumulated register without changing fixed mass. The
 * extensive total energy E=U+|p|^2/(2m) receives dEtot; the updated internal
 * energy is then E_new-|p_new|^2/(2m). On failure neither object is mutated.
 */
static inline int gfs_kdk_lts_consume_register(
    GfsKdkLtsConservedCell *c, GfsKdkLtsRegister *r)
{
    if (!c || !r || !(c->mass > 0.0) || !(c->U > 0.0) ||
        !isfinite(c->mass) || !isfinite(c->vx) || !isfinite(c->vy) ||
        !isfinite(c->U) || !isfinite(r->dpx) || !isfinite(r->dpy) ||
        !isfinite(r->dEtot))
        return 0;

    const double old_k = 0.5 * c->mass * (c->vx*c->vx + c->vy*c->vy);
    const double px = c->mass*c->vx + r->dpx;
    const double py = c->mass*c->vy + r->dpy;
    const double vx = px / c->mass;
    const double vy = py / c->mass;
    const double new_k = 0.5 * c->mass * (vx*vx + vy*vy);
    const double U = c->U + old_k + r->dEtot - new_k;
    if (!isfinite(old_k) || !isfinite(px) || !isfinite(py) ||
        !isfinite(vx) || !isfinite(vy) || !isfinite(new_k) ||
        !(U > 0.0) || !isfinite(U))
        return 0;

    c->vx = vx;
    c->vy = vy;
    c->U = U;
    gfs_kdk_lts_register_clear(r);
    return 1;
}

/* Second-order inactive-cell predictor for face states at local event times. */
static inline int gfs_kdk_lts_predict(
    const GfsKdkLtsPredictedState *base,
    double ax, double ay, double Udot, double elapsed,
    GfsKdkLtsPredictedState *pred)
{
    if (!base || !pred || !isfinite(ax) || !isfinite(ay) ||
        !isfinite(Udot) || !isfinite(elapsed) || elapsed < 0.0 ||
        !isfinite(base->x) || !isfinite(base->y) ||
        !isfinite(base->vx) || !isfinite(base->vy) || !(base->U > 0.0))
        return 0;
    const double h2 = 0.5 * elapsed * elapsed;
    pred->x = base->x + base->vx*elapsed + ax*h2;
    pred->y = base->y + base->vy*elapsed + ay*h2;
    pred->vx = base->vx + ax*elapsed;
    pred->vy = base->vy + ay*elapsed;
    pred->U = base->U + Udot*elapsed;
    return isfinite(pred->x) && isfinite(pred->y) &&
           isfinite(pred->vx) && isfinite(pred->vy) &&
           pred->U > 0.0 && isfinite(pred->U);
}

#ifdef __cplusplus
}
#endif
#endif
