#ifndef GFS_AFFINE_SPIN_RECONSTRUCTION_20261027_H
#define GFS_AFFINE_SPIN_RECONSTRUCTION_20261027_H

#include <math.h>

#if defined(__CUDACC__)
#define GFS_AFFINE_HD __host__ __device__
#else
#define GFS_AFFINE_HD
#endif

/* Q is the mass central second-moment tensor about the mass centroid.
 * A is a reconstructed velocity gradient. S is an independently supplied
 * intrinsic angular momentum, not the generator's orbital momentum. */
struct gfs_affine_spin_input {
    double qxx,qxy,qyy;
    double axx,axy,ayx,ayy;
    double spin;
};

struct gfs_affine_spin_output {
    double axx,axy,ayx,ayy;
    double reference_spin,reference_omega,state_omega;
    double nonrigid_energy,spin_energy;
};

GFS_AFFINE_HD static inline double gfs_affine_energy(
    double axx,double axy,double ayx,double ayy,double qxx,double qxy,double qyy)
{
    return .5*((axx*axx+ayx*ayx)*qxx + 2.*(axx*axy+ayx*ayy)*qxy +
               (axy*axy+ayy*ayy)*qyy);
}

/* Remove the rigid projection already present in A before adding the
 * independently evolved spin. This avoids counting local vorticity twice.
 * The remaining Q-orthogonal affine component is retained, not turned into
 * heat. Returns false without changing *out for invalid/nonfinite inputs. */
GFS_AFFINE_HD static inline bool gfs_affine_spin_reconstruct(
    const gfs_affine_spin_input &in, gfs_affine_spin_output *out)
{
    if (!out || !isfinite(in.qxx) || !isfinite(in.qxy) || !isfinite(in.qyy) ||
        !isfinite(in.axx) || !isfinite(in.axy) || !isfinite(in.ayx) ||
        !isfinite(in.ayy) || !isfinite(in.spin)) return false;
    const double det=in.qxx*in.qyy-in.qxy*in.qxy;
    const double inertia=in.qxx+in.qyy;
    if (!(in.qxx>0.) || !(in.qyy>0.) || !(det>0.) ||
        !isfinite(det) || !(inertia>0.) || !isfinite(inertia)) return false;
    gfs_affine_spin_output a;
    a.reference_spin=in.ayx*in.qxx-in.axy*in.qyy+(in.ayy-in.axx)*in.qxy;
    a.reference_omega=a.reference_spin/inertia;
    a.state_omega=in.spin/inertia;
    const double bxy=in.axy+a.reference_omega;
    const double byx=in.ayx-a.reference_omega;
    a.axx=in.axx; a.axy=bxy-a.state_omega;
    a.ayx=byx+a.state_omega; a.ayy=in.ayy;
    a.nonrigid_energy=gfs_affine_energy(in.axx,bxy,byx,in.ayy,in.qxx,in.qxy,in.qyy);
    a.spin_energy=in.spin*in.spin/(2.*inertia);
    if (!isfinite(a.reference_spin) || !isfinite(a.reference_omega) ||
        !isfinite(a.state_omega) || !isfinite(a.axy) || !isfinite(a.ayx) ||
        !isfinite(a.nonrigid_energy) || !isfinite(a.spin_energy) ||
        a.nonrigid_energy<0. || a.spin_energy<0.) return false;
    *out=a;
    return true;
}

/* Optional excess-spin parameterization: in.spin is delta_S relative to
 * the rigid spin already contained in A, not the total intrinsic spin.
 * delta_S=0 leaves the input gradient bitwise unchanged. Its time evolution
 * requires d(delta_S)/dt = physical_torque - d(S_A)/dt; this function is
 * only a local algebraic reconstruction and does not compute that RHS. */
GFS_AFFINE_HD static inline bool gfs_affine_spin_reconstruct_excess(
    const gfs_affine_spin_input &in,gfs_affine_spin_output *out)
{
    if(!out || !isfinite(in.spin))return false;
    gfs_affine_spin_input probe=in;
    probe.spin=0.;
    gfs_affine_spin_output reference;
    if(!gfs_affine_spin_reconstruct(probe,&reference))return false;
    const double total_spin=reference.reference_spin+in.spin;
    if(!isfinite(total_spin))return false;
    probe.spin=total_spin;
    gfs_affine_spin_output result;
    if(!gfs_affine_spin_reconstruct(probe,&result))return false;
    if(in.spin==0.) {
        result.axx=in.axx;result.axy=in.axy;result.ayx=in.ayx;result.ayy=in.ayy;
    }
    *out=result;
    return true;
}

#undef GFS_AFFINE_HD
#endif
