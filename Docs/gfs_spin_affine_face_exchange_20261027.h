#ifndef GFS_SPIN_AFFINE_FACE_EXCHANGE_20261027_H
#define GFS_SPIN_AFFINE_FACE_EXCHANGE_20261027_H

#include "gfs_spin_face_exchange_20261027.h"

#if defined(__CUDACC__)
#define GFS_SPIN_AFFINE_HD __host__ __device__
#else
#define GFS_SPIN_AFFINE_HD
#endif

/* Explicit unresolved non-rigid affine kinetic reservoir. The face flux,
 * torque and both reservoir increments are supplied by the caller; this
 * accounting wrapper does not define their physical constitutive law. */
struct gfs_spin_affine_cell {
    gfs_spin_cell mechanical;
    double nonrigid_energy;
};

struct gfs_spin_affine_transfer {
    gfs_spin_face_transfer mechanical;
    double nonrigid_delta_i, nonrigid_delta_j;
};

GFS_SPIN_AFFINE_HD static inline double gfs_spin_affine_thermal(
    const gfs_spin_affine_cell &a)
{
    return gfs_spin_thermal_energy(a.mechanical)-a.nonrigid_energy;
}

/* Transactional update of pair mass, momentum, total energy, mass moments
 * and intrinsic angular momentum. The residual thermal energy excludes
 * both spin and non-rigid kinetic energy. Caller is responsible for the
 * actual physical fluxes, mass-moment volume term and shape evolution. */
GFS_SPIN_AFFINE_HD static inline int gfs_spin_affine_face_exchange(
    gfs_spin_affine_cell *i,gfs_spin_affine_cell *j,const gfs_spin_affine_transfer &f)
{
    if(!i || !j || i==j || !isfinite(i->nonrigid_energy) ||
       !isfinite(j->nonrigid_energy) || !isfinite(f.nonrigid_delta_i) ||
       !isfinite(f.nonrigid_delta_j))return -1;
    if(i->nonrigid_energy<0. || j->nonrigid_energy<0.)return -4;
    if(!gfs_spin_cell_finite(i->mechanical) || !gfs_spin_cell_finite(j->mechanical) ||
       !(i->mechanical.m>0.) || !(j->mechanical.m>0.) ||
       !(i->mechanical.inertia>0.) || !(j->mechanical.inertia>0.))return -2;
    if(!(gfs_spin_affine_thermal(*i)>0.) || !(gfs_spin_affine_thermal(*j)>0.))return -3;
    gfs_spin_affine_cell a=*i,b=*j;
    const int status=gfs_spin_face_exchange(&a.mechanical,&b.mechanical,f.mechanical);
    if(status)return status;
    a.nonrigid_energy+=f.nonrigid_delta_i;
    b.nonrigid_energy+=f.nonrigid_delta_j;
    if(!isfinite(a.nonrigid_energy) || !isfinite(b.nonrigid_energy))return -1;
    if(a.nonrigid_energy<0. || b.nonrigid_energy<0.)return -4;
    const double ui=gfs_spin_affine_thermal(a),uj=gfs_spin_affine_thermal(b);
    if(!isfinite(ui) || !isfinite(uj) || !(ui>0.) || !(uj>0.))return -3;
    *i=a;*j=b;
    return 0;
}

#undef GFS_SPIN_AFFINE_HD
#endif
