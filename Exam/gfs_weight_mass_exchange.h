#ifndef GFS_WEIGHT_MASS_EXCHANGE_H
#define GFS_WEIGHT_MASS_EXCHANGE_H
#include "gfs_moment.h"

/* Contribution of changing power weights to a face's normal speed.
 * w2dot is the derivative of the signed power weight (length squared),
 * and distance is the positive generator separation. The normal points
 * i -> j. This matches the weight term of the power-face kinematics.
 * It is a proposed incremental exchange relative to a separately specified
 * fixed-mass baseline, not the full Euler ALE mass flux. */
GFS_MOM_HD static inline int gfs_weight_exchange_flux(
    double w2dot_i,double w2dot_j,double distance,double rho_face,
    double ux,double uy,double specific_total_energy,double flux[4]){
    if(!isfinite(w2dot_i) || !isfinite(w2dot_j) ||
       !(distance>0) || !isfinite(distance) || !(rho_face>0) || !isfinite(rho_face) ||
       !isfinite(ux) || !isfinite(uy) || !isfinite(specific_total_energy))return 0;
    double normal_weight_speed=.5*(w2dot_i-w2dot_j)/distance;
    double mflux=-rho_face*normal_weight_speed;
    flux[0]=mflux;flux[1]=mflux*ux;flux[2]=mflux*uy;flux[3]=mflux*specific_total_energy;
    return isfinite(flux[0]) && isfinite(flux[1]) && isfinite(flux[2]) && isfinite(flux[3]);
}

/* One canonical face point and common state. Apply once with opposite signs
 * in fixed reference frames. The caller supplies pressure work separately,
 * chooses a physically justified face state, controls positivity/timestep,
 * and owns race-free shared-face/MPI accumulation. */
GFS_MOM_HD static inline void gfs_weight_exchange_pair(
    GfsMomentState *ri,GfsMomentState *rj,double face_x,double face_y,
    double ox_i,double oy_i,double ox_j,double oy_j,double length,const double flux[4]){
    gfs_moment_face_rhs(ri,face_x-ox_i,face_y-oy_i,flux,length);
    gfs_moment_face_rhs(rj,face_x-ox_j,face_y-oy_j,flux,-length);
}
#endif
