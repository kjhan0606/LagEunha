#ifndef GFS_WEIGHT_MOTION_REMAP_20261027_H
#define GFS_WEIGHT_MOTION_REMAP_20261027_H

#include "gfs_spin_affine_face_exchange_20261027.h"

#if defined(__CUDACC__)
#define GFS_WEIGHT_HD __host__ __device__
#else
#define GFS_WEIGHT_HD
#endif

/* Caller-supplied face states. Energy is specific TOTAL hydro energy at the
 * face; nonrigid is its separately transported reservoir per unit mass.
 * The caller must supply consistent reconstruction, EOS and gravity work. */
struct gfs_weight_face_state {
    double density,vx,vy,specific_energy,specific_nonrigid;
};

struct gfs_weight_motion_input {
    double cx,cy,face_length,site_distance,weight_i,weight_j;
    double wdot_i,wdot_j,normal_separation_rate,dt;
    double inertia_i_new,inertia_j_new;
    gfs_weight_face_state state_i,state_j;
};

GFS_WEIGHT_HD static inline bool gfs_weight_face_state_valid(const gfs_weight_face_state &s) {
    return isfinite(s.density) && s.density>0. && isfinite(s.vx) && isfinite(s.vy) &&
           isfinite(s.specific_energy) && isfinite(s.specific_nonrigid) && s.specific_nonrigid>=0.;
}

/* For a power face, let d be generator separation, n point i->j, and
 * DeltaW=Wi-Wj. Relative to the Voronoi face at the same generator positions,
 * the exact normal face-speed correction is
 * (DeltaWdot)/(2*d)-DeltaW*(n dot (vj-vi))/(2*d^2). A positive speed expands
 * i and sweeps material from j into i. This is only the weight-relative
 * geometric mass flux, not the full ALE/Riemann, pressure or gravity flux. */
GFS_WEIGHT_HD static inline bool gfs_weight_motion_transfer(
    const gfs_weight_motion_input &in,gfs_spin_affine_transfer *out)
{
    if(!out || !isfinite(in.cx) || !isfinite(in.cy) || !isfinite(in.face_length) ||
       !(in.face_length>=0.) || !isfinite(in.site_distance) || !(in.site_distance>0.) ||
       !isfinite(in.weight_i) || !isfinite(in.weight_j) ||
       !isfinite(in.wdot_i) || !isfinite(in.wdot_j) ||
       !isfinite(in.normal_separation_rate) || !isfinite(in.dt) || !(in.dt>=0.) ||
       !isfinite(in.inertia_i_new) || !(in.inertia_i_new>0.) ||
       !isfinite(in.inertia_j_new) || !(in.inertia_j_new>0.) ||
       !gfs_weight_face_state_valid(in.state_i) || !gfs_weight_face_state_valid(in.state_j))return false;
    const double speed=(in.wdot_i-in.wdot_j)/(2.*in.site_distance)
        -(in.weight_i-in.weight_j)*in.normal_separation_rate/
          (2.*in.site_distance*in.site_distance);
    if(!isfinite(speed))return false;
    const gfs_weight_face_state &donor=(speed>0.)?in.state_j:in.state_i;
    const double dm=-donor.density*in.face_length*speed*in.dt;
    gfs_spin_affine_transfer f;
    f.mechanical.cx=in.cx;f.mechanical.cy=in.cy;f.mechanical.dm=dm;
    f.mechanical.dpx=dm*donor.vx;f.mechanical.dpy=dm*donor.vy;
    f.mechanical.denergy=dm*donor.specific_energy;
    f.mechanical.inertia_i_new=in.inertia_i_new;f.mechanical.inertia_j_new=in.inertia_j_new;
    f.nonrigid_delta_i=-dm*donor.specific_nonrigid;
    f.nonrigid_delta_j=dm*donor.specific_nonrigid;
    if(!isfinite(dm) || !isfinite(f.mechanical.dpx) || !isfinite(f.mechanical.dpy) ||
       !isfinite(f.mechanical.denergy) || !isfinite(f.nonrigid_delta_i) ||
       !isfinite(f.nonrigid_delta_j))return false;
    *out=f;
    return true;
}

#undef GFS_WEIGHT_HD
#endif
