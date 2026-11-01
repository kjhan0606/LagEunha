#ifndef GFS_FIXED_POLYGON_GEOMETRY_H
#define GFS_FIXED_POLYGON_GEOMETRY_H

/* One compile-time count for CUDA clipping and host wall/event geometry.
 * The default preserves the active 1024-side calculation. The 2048-side
 * build is an isolated nested-domain sensitivity experiment, not a new IC.
 */
#ifndef GFS_FIXEDPOLY_SIDES
#define GFS_FIXEDPOLY_SIDES 1024
#endif
#if GFS_FIXEDPOLY_SIDES != 1024 && GFS_FIXEDPOLY_SIDES != 2048
#error "This comparison supports only 1024 or 2048 polygon sides"
#endif

#endif
