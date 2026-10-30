#ifndef GFS_KEPLER_SLIP_WALL_H
#define GFS_KEPLER_SLIP_WALL_H

/* Accepted-impact impulse primitive for the native CUDA/MPI solver.
 * This header does NOT locate events, move positions, process RK predictors,
 * integrate pressure traction, or enable a boundary condition. The caller
 * must first commit an impact at the boundary used by its tessellator.
 * No thermal-energy or mass field is changed to repair kinetic roundoff.
 */
#include <float.h>
#include <math.h>

#ifdef __CUDACC__
#define GFS_SLIP_HD __host__ __device__
#else
#define GFS_SLIP_HD
#endif

enum GfsSlipStatus {
    GFS_SLIP_INVALID = -1,
    GFS_SLIP_NO_IMPULSE = 0,
    GFS_SLIP_REFLECTED = 1
};

typedef struct GfsSlipImpulse {
    double vx_after, vy_after;
    double wall_jx, wall_jy;
    double wall_torque_impulse;
    double fluid_kinetic_roundoff;
} GfsSlipImpulse;

GFS_SLIP_HD static inline int gfs_slip_impulse(
    double rx, double ry, double nx, double ny,
    double mass, double vx, double vy, GfsSlipImpulse *out)
{
    if (!out || !isfinite(rx) || !isfinite(ry) || !isfinite(nx) ||
        !isfinite(ny) || !isfinite(mass) || !(mass > 0.) ||
        !isfinite(vx) || !isfinite(vy) ||
        fabs(nx*nx + ny*ny - 1.) > 64.*DBL_EPSILON)
        return GFS_SLIP_INVALID;
    const double normal_speed = nx*vx + ny*vy;
    if (!isfinite(normal_speed)) return GFS_SLIP_INVALID;
    const double reflection = normal_speed > 0. ? 2.*normal_speed : 0.;
    const double ux = vx - reflection*nx;
    const double uy = vy - reflection*ny;
    const double jx = mass*reflection*nx;
    const double jy = mass*reflection*ny;
    /* Algebraically equals K_after-K_before without subtracting large K. */
    const double dk = .5*mass*((ux-vx)*(ux+vx) + (uy-vy)*(uy+vy));
    const double torque = rx*jy - ry*jx;
    if (!isfinite(ux) || !isfinite(uy) || !isfinite(jx) || !isfinite(jy) ||
        !isfinite(dk) || !isfinite(torque)) return GFS_SLIP_INVALID;
    out->vx_after = ux;
    out->vy_after = uy;
    out->wall_jx = jx;
    out->wall_jy = jy;
    out->wall_torque_impulse = torque;
    out->fluid_kinetic_roundoff = dk;
    return normal_speed > 0. ? GFS_SLIP_REFLECTED : GFS_SLIP_NO_IMPULSE;
}

/* Plane convention matches fixed_polygon_postStage_blend:
 * n dot (position-center) = apothem, with outward unit normal n.
 * The event controller must separately reject simultaneous corner impacts.
 */
GFS_SLIP_HD static inline int gfs_slip_plane_impact(
    double px, double py, double cx, double cy, double nx, double ny,
    double apothem, double mass, double vx, double vy, GfsSlipImpulse *out)
{
    const double rx = px-cx, ry = py-cy;
    const double scale = fmax(1., fmax(fabs(apothem), fabs(rx)+fabs(ry)));
    if (!isfinite(apothem) || !(apothem > 0.) ||
        !isfinite(rx) || !isfinite(ry) ||
        !isfinite(nx) || !isfinite(ny) ||
        fabs(nx*rx+ny*ry-apothem) > 64.*DBL_EPSILON*scale)
        return GFS_SLIP_INVALID;
    return gfs_slip_impulse(rx, ry, nx, ny, mass, vx, vy, out);
}

/* Only for a truly circular geometric boundary, never polygon clipping. */
GFS_SLIP_HD static inline int gfs_slip_circle_impact(
    double px, double py, double cx, double cy, double radius,
    double mass, double vx, double vy, GfsSlipImpulse *out)
{
    const double rx = px-cx, ry = py-cy;
    const double r = hypot(rx, ry);
    if (!isfinite(radius) || !(radius > 0.) || !isfinite(r) || !(r > 0.) ||
        fabs(r-radius) > 64.*DBL_EPSILON*fmax(1.,radius))
        return GFS_SLIP_INVALID;
    return gfs_slip_impulse(rx, ry, rx/r, ry/r, mass, vx, vy, out);
}

#undef GFS_SLIP_HD
#endif
