#ifndef GFS_FACE_WORK_H
#define GFS_FACE_WORK_H
/* Shared C/CUDA traction work. The caller must supply the actual boundary
 * velocity. This helper does not enforce material-boundary kinematics.
 * S is the outward oriented face-area vector; tauS is the stress traction
 * integrated over the face. Units: energy/time for work, volume/time for sweep.
 */
#ifdef __CUDACC__
#define GFS_FACE_HD __host__ __device__
#else
#define GFS_FACE_HD
#endif
static GFS_FACE_HD inline double gfs_face_sweep2d(
    double sx, double sy, double wx, double wy)
{
    return wx*sx + wy*sy;
}
static GFS_FACE_HD inline double gfs_face_work2d(
    double pressure, double sx, double sy, double taux, double tauy,
    double wx, double wy)
{
    return -pressure*gfs_face_sweep2d(sx, sy, wx, wy)
           + taux*wx + tauy*wy;
}
/* Difference in power from using a work velocity distinct from geometry.
 * Zero net sweep alone does not make the sum of this defect vanish when
 * face pressures differ. Diagnostic only: never an energy correction.
 */
static GFS_FACE_HD inline double gfs_face_work_defect2d(
    double pressure, double sx, double sy, double taux, double tauy,
    double workx, double worky, double meshx, double meshy)
{
    return gfs_face_work2d(pressure, sx, sy, taux, tauy,
                         workx-meshx, worky-meshy);
}
#undef GFS_FACE_HD
#endif
