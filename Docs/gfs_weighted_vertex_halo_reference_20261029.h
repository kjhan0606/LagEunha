#ifndef GFS_WEIGHTED_VERTEX_HALO_REFERENCE_20261029_H
#define GFS_WEIGHTED_VERTEX_HALO_REFERENCE_20261029_H

#include "gfs_vertex_halo.h"

#ifdef __CUDACC__
#define GFS_WVH_INLINE static __host__ __device__ __forceinline__
#else
#define GFS_WVH_INLINE static inline
#endif

/* Experimental sufficient condition, not connected to the native solver.
 * Wmax must bound EVERY generator weight at this exact RK stage, including
 * omitted generators. W and Wmax have units of length squared. Complete
 * interior search, a valid convex polygon, and actual halo coverage are
 * separate preconditions. The existing native weighted search is unchanged.
 *
 * With delta=Wmax-Wp>=0, every omitted standard Laguerre plane is harmless
 * if the disk centered at p+(dx,dy), radius sqrt(dx^2+dy^2+delta), is inside
 * the covered rectangle. This also suffices for clamp(wfrac,.05,.95):
 * clamp(.5*(d^2+Wp-Wq),.05*d^2,.95*d^2) >= .5*(d^2-delta).
 * The strict directed bounds concern the represented inputs, not errors
 * in the polygon or floating point clipping that produced them. */
GFS_WVH_INLINE int gfs_weighted_vertex_halo_disk_inside_reference(
    double px, double py, double dx, double dy,
    double owner_weight, double global_weight_upper,
    double xmin, double ymin, double xmax, double ymax)
{
    if (!isfinite(owner_weight) || !isfinite(global_weight_upper)
        || global_weight_upper < owner_weight) return 0;
    /* Exact ordinary-Voronoi reduction, also for equal nonzero weights. */
    if (global_weight_upper == owner_weight)
        return gfs_vertex_halo_disk_inside(px,py,dx,dy,xmin,ymin,xmax,ymax);
    if (!isfinite(px) || !isfinite(py) || !isfinite(dx) || !isfinite(dy)
        || !isfinite(xmin) || !isfinite(ymin) || !isfinite(xmax) || !isfinite(ymax)
        || !(xmin<xmax) || !(ymin<ymax)) return 0;
    double delta_hi=gfs_vh_add_hi(global_weight_upper,-owner_weight);
    if (!(delta_hi>0) || !isfinite(delta_hi)) return 0;
    double vxlo=gfs_vh_add_lo(px,dx), vxhi=gfs_vh_add_hi(px,dx);
    double vylo=gfs_vh_add_lo(py,dy), vyhi=gfs_vh_add_hi(py,dy);
    double clearance=fmin(fmin(gfs_vh_sub_lo(vxlo,xmin),gfs_vh_sub_lo(xmax,vxhi)),
                          fmin(gfs_vh_sub_lo(vylo,ymin),gfs_vh_sub_lo(ymax,vyhi)));
    if (!(clearance>0) || !isfinite(clearance)) return 0;
    double radius2_hi=gfs_vh_add_hi(
        gfs_vh_add_hi(gfs_vh_mul_hi(dx,dx),gfs_vh_mul_hi(dy,dy)),delta_hi);
    double clearance2_lo=gfs_vh_mul_lo(clearance,clearance);
    return isfinite(radius2_hi) && clearance2_lo>radius2_hi;
}

#undef GFS_WVH_INLINE
#endif
