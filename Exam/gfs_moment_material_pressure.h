#ifndef GFS_MOMENT_MATERIAL_PRESSURE_H
#define GFS_MOMENT_MATERIAL_PRESSURE_H
#include "gfs_moment.h"

/* Pressure part of an inviscid material-face flux in a shared normal frame.
 * This is valid for zero mass flux only when the physical normal interface
 * velocity and the face's actual geometric normal velocity agree to the
 * explicitly supplied absolute tolerance. It is not a Riemann solver and
 * must not be used to force a nonmaterial Voronoi face to be material. */
GFS_MOM_HD static inline int gfs_material_pressure_flux(
    double pressure,double nx,double ny,double material_vn,double geometric_vn,
    double velocity_tolerance,double flux[4]){
    if(!(pressure>=0) || !isfinite(pressure) || !isfinite(nx) || !isfinite(ny) ||
       !isfinite(material_vn) || !isfinite(geometric_vn) ||
       !(velocity_tolerance>=0) || !isfinite(velocity_tolerance) ||
       fabs(nx*nx+ny*ny-1)>1e-12 || fabs(material_vn-geometric_vn)>velocity_tolerance)return 0;
    double px=pressure*nx,py=pressure*ny,e=pressure*material_vn;
    if(!isfinite(px) || !isfinite(py) || !isfinite(e))return 0;
    flux[0]=0;flux[1]=px;flux[2]=py;flux[3]=e;return 1;
}
#endif
