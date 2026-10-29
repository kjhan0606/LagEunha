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

/* One face quadrature point with its own quadrature weight. This point
 * helper alone is not an exact integral of a segment's second moments. */
GFS_MOM_HD static inline void gfs_weight_exchange_point(
    GfsMomentState *ri,GfsMomentState *rj,double face_x,double face_y,
    double ox_i,double oy_i,double ox_j,double oy_j,double quadrature_weight,const double flux[4]){
    gfs_moment_face_rhs(ri,face_x-ox_i,face_y-oy_i,flux,quadrature_weight);
    gfs_moment_face_rhs(rj,face_x-ox_j,face_y-oy_j,flux,-quadrature_weight);
}

/* Two-point Gauss integration for a straight face with a constant shared
 * flux: exact up to quadratic geometric moments. The caller supplies
 * pressure work separately, chooses the face state, controls positivity
 * and timestep, and owns race-free shared-face/MPI accumulation. */
GFS_MOM_HD static inline int gfs_weight_exchange_segment(
    GfsMomentState *ri,GfsMomentState *rj,double ax,double ay,double bx,double by,
    double ox_i,double oy_i,double ox_j,double oy_j,const double flux[4]){
    if(!isfinite(ax) || !isfinite(ay) || !isfinite(bx) || !isfinite(by) ||
       !isfinite(ox_i) || !isfinite(oy_i) || !isfinite(ox_j) || !isfinite(oy_j))return 0;
    for(int k=0;k<4;k++)if(!isfinite(flux[k]))return 0;
    double dx=bx-ax,dy=by-ay,length=hypot(dx,dy);
    if(!(length>0) || !isfinite(length))return 0;
    double cx=ax+.5*dx,cy=ay+.5*dy;
    const double scale=.28867513459481288225; /* 1/(2 sqrt(3)) */
    gfs_weight_exchange_point(ri,rj,cx-scale*dx,cy-scale*dy,ox_i,oy_i,ox_j,oy_j,.5*length,flux);
    gfs_weight_exchange_point(ri,rj,cx+scale*dx,cy+scale*dy,ox_i,oy_i,ox_j,oy_j,.5*length,flux);
    return 1;
}

/* Two-point segment exchange with constant density/specific internal energy
 * and a single shared affine velocity at the midpoint. This represents only
 * weight-driven advection; pressure traction/work remain separate. No state
 * is modified if either quadrature state or resulting RHS is invalid. */
GFS_MOM_HD static inline int gfs_weight_exchange_segment_affine(
    GfsMomentState *ri,GfsMomentState *rj,double ax,double ay,double bx,double by,
    double ox_i,double oy_i,double ox_j,double oy_j,double wdot_i,double wdot_j,
    double distance,double rho,double ux,double uy,const double grad[4],double specific_internal_energy){
    if(ri==rj || !isfinite(ax) || !isfinite(ay) || !isfinite(bx) || !isfinite(by) ||
       !isfinite(ox_i) || !isfinite(oy_i) || !isfinite(ox_j) || !isfinite(oy_j) ||
       !isfinite(ux) || !isfinite(uy) || !(specific_internal_energy>0) || !isfinite(specific_internal_energy))return 0;
    for(int k=0;k<4;k++)if(!isfinite(grad[k]))return 0;
    double dx=bx-ax,dy=by-ay,length=hypot(dx,dy);
    if(!(length>0) || !isfinite(length))return 0;
    double cx=ax+.5*dx,cy=ay+.5*dy,fluxes[2][4];
    const double scale=.28867513459481288225;
    for(int k=0;k<2;k++){
        double sign=k?1.:-1.;double rx=sign*scale*dx,ry=sign*scale*dy;
        double vx=ux+grad[0]*rx+grad[1]*ry,vy=uy+grad[2]*rx+grad[3]*ry;
        double specific_energy=specific_internal_energy+.5*(vx*vx+vy*vy);
        if(!gfs_weight_exchange_flux(wdot_i,wdot_j,distance,rho,vx,vy,specific_energy,fluxes[k]))return 0;
    }
    GfsMomentState ai=*ri,aj=*rj;
    for(int k=0;k<2;k++){
        double sign=k?1.:-1.;
        gfs_weight_exchange_point(&ai,&aj,cx+sign*scale*dx,cy+sign*scale*dy,
            ox_i,oy_i,ox_j,oy_j,.5*length,fluxes[k]);
    }
    for(int q=0;q<GFM_NSTATE;q++)if(!isfinite(ai.a[q]) || !isfinite(aj.a[q]))return 0;
    *ri=ai;*rj=aj;return 1;
}
#endif
