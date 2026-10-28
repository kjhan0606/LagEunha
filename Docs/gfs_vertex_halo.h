#ifndef GFS_VERTEX_HALO_H
#define GFS_VERTEX_HALO_H

#include <math.h>

#ifdef __CUDACC__
#define GFS_VH_INLINE static __host__ __device__ __forceinline__
#else
#define GFS_VH_INLINE static inline
#endif

/* Directed bounds for the represented binary inputs. These do not bound
 * errors in constructing the polygon itself; native geometry parity is a
 * separate gate. No negative geometric margin is accepted by a tolerance. */
GFS_VH_INLINE double gfs_vh_add_lo(double a, double b) {
#ifdef __CUDA_ARCH__
    return __dadd_rd(a,b);
#else
    return nextafter(a+b,-INFINITY);
#endif
}
GFS_VH_INLINE double gfs_vh_add_hi(double a, double b) {
#ifdef __CUDA_ARCH__
    return __dadd_ru(a,b);
#else
    return nextafter(a+b,INFINITY);
#endif
}
GFS_VH_INLINE double gfs_vh_sub_lo(double a, double b) {
#ifdef __CUDA_ARCH__
    return __dsub_rd(a,b);
#else
    return nextafter(a-b,-INFINITY);
#endif
}
GFS_VH_INLINE double gfs_vh_mul_lo(double a, double b) {
#ifdef __CUDA_ARCH__
    return __dmul_rd(a,b);
#else
    return nextafter(a*b,-INFINITY);
#endif
}
GFS_VH_INLINE double gfs_vh_mul_hi(double a, double b) {
#ifdef __CUDA_ARCH__
    return __dmul_ru(a,b);
#else
    return nextafter(a*b,INFINITY);
#endif
}

/* The disk centered at p+(dx,dy), radius hypot(dx,dy), must lie strictly
 * inside the rectangle actually covered by halo data. W=0 only. */
GFS_VH_INLINE int gfs_vertex_halo_disk_inside(
    double px, double py, double dx, double dy,
    double xmin, double ymin, double xmax, double ymax)
{
    if (!isfinite(px) || !isfinite(py) || !isfinite(dx) || !isfinite(dy)
        || !isfinite(xmin) || !isfinite(ymin) || !isfinite(xmax) || !isfinite(ymax)
        || !(xmin<xmax) || !(ymin<ymax)) return 0;
    double vxlo=gfs_vh_add_lo(px,dx), vxhi=gfs_vh_add_hi(px,dx);
    double vylo=gfs_vh_add_lo(py,dy), vyhi=gfs_vh_add_hi(py,dy);
    double clearance=fmin(fmin(gfs_vh_sub_lo(vxlo,xmin),gfs_vh_sub_lo(xmax,vxhi)),
                          fmin(gfs_vh_sub_lo(vylo,ymin),gfs_vh_sub_lo(ymax,vyhi)));
    if (!(clearance>0) || !isfinite(clearance)) return 0;
    double radius2_hi=gfs_vh_add_hi(gfs_vh_mul_hi(dx,dx),gfs_vh_mul_hi(dy,dy));
    double clearance2_lo=gfs_vh_mul_lo(clearance,clearance);
    return isfinite(radius2_hi) && clearance2_lo>radius2_hi;
}

#undef GFS_VH_INLINE
#endif
