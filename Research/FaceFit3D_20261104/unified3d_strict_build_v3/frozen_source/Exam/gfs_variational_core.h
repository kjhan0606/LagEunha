#ifndef GFS_VARIATIONAL_CORE_H
#define GFS_VARIATIONAL_CORE_H

/*
 * Dimension-independent INTERIOR algebra for unweighted Voronoi cells.
 * Not a complete shock solver or a boundary/time-integration prescription.
 * Coordinates and face centroids below are relative to generator i.
 * Conservative claims require reciprocal faces and closed cell surfaces.
 */
#include <math.h>
#ifdef __CUDACC__
#define GFS_CORE_INLINE static __host__ __device__ __forceinline__
#else
#define GFS_CORE_INLINE static inline
#endif

typedef struct {
    double bi[3], bj[3], normal[3], area, separation;
    int dim;
} GfsVariationalFace;

GFS_CORE_INLINE int gfs_variational_face(int dim, const double r[3],
        const double centroid[3], double area, GfsVariationalFace *f)
{
    if (!f || dim < 1 || dim > 3 || !(area > 0) || !isfinite(area))
        return 0;
    double d2=0;
    for (int a=0;a<dim;a++) {
        if (!isfinite(r[a]) || !isfinite(centroid[a])) return 0;
        d2+=r[a]*r[a];
    }
    if (!(d2>0) || !isfinite(d2)) return 0;
    const double invd=1.0/sqrt(d2);
    f->dim=dim; f->area=area; f->separation=1.0/invd;
    for (int a=0;a<3;a++) {
        f->normal[a]=a<dim?r[a]*invd:0;
        const double off=a<dim?centroid[a]-.5*r[a]:0;
        f->bi[a]=.5*area*f->normal[a]+area*off*invd;
        f->bj[a]=.5*area*f->normal[a]-area*off*invd;
    }
    return 1;
}

/* Relative face sweep. Summing a closed cell gives its actual volume rate. */
GFS_CORE_INLINE double gfs_variational_sweep(const GfsVariationalFace *f,
        const double vi[3], const double vj[3])
{
    double result=0;
    for (int a=0;a<f->dim;a++) result+=f->bj[a]*(vj[a]-vi[a]);
    return result;
}

GFS_CORE_INLINE void gfs_variational_pressure(const GfsVariationalFace *f,
        double pi, double pj, double force[3])
{
    for (int a=0;a<3;a++) force[a]=(pi-pj)*f->bi[a];
}

GFS_CORE_INLINE int gfs_variational_adiabatic(double U0, double V0,
        double V1, double gamma, double *U1)
{
    if (!U1 || !(U0>0) || !(V0>0) || !(V1>0) || !(gamma>1) ||
        !isfinite(U0) || !isfinite(V0) || !isfinite(V1) || !isfinite(gamma))
        return 0;
    const double value=U0*pow(V0/V1,gamma-1);
    if (!(value>0) || !isfinite(value)) return 0;
    *U1=value;
    return 1;
}

/*
 * Candidate reconstructed Monaghan-type pair dissipation, not promoted.
 * delta_v_face is reconstructed v_j-v_i at the same face.
 * Actual approach is also required, ensuring nonnegative pair heating.
 * One conventional alpha=1, beta=2 linear/quadratic closure. Reconstruction
 * removes an exactly affine velocity difference. The unreconstructed
 * approach sets the nonlinear signal speed. No restoring force is added.
 */
GFS_CORE_INLINE double gfs_variational_pair_pressure(
        const GfsVariationalFace *f, const double vi[3], const double vj[3],
        const double delta_v_face[3], double rhoi, double ci,
        double rhoj, double cj, double compression_fraction,
        double length_i,double length_j)
{
    double actual=0,reconstructed=0;
    for (int a=0;a<f->dim;a++) {
        actual+=(vj[a]-vi[a])*f->normal[a];
        reconstructed+=delta_v_face[a]*f->normal[a];
    }
    if (!(actual<0)) return 0;
    /* Reconstruction is a bounded viscosity coefficient, never a finite
     * force at zero actual approach. This also bounds q/|dv.n|. */
    const double closing=fmin(-actual,fmax(fmax(0.0,-reconstructed),
                              fmin(1.0,fmax(0.0,compression_fraction))*(-actual)));
    if (!(closing>0)) return 0;
    if (!(rhoi>0) || !(rhoj>0) || !(ci>0) || !(cj>0) ||
        !isfinite(rhoi) || !isfinite(rhoj) || !isfinite(ci) || !isfinite(cj) ||
        !(length_i>0) || !(length_j>0) || !isfinite(length_i) || !isfinite(length_j))
        return 0;
    const double rho=.5*rhoi+.5*rhoj,c=.5*ci+.5*cj;
    const double scale=(.5*length_i+.5*length_j)/f->separation;
    const double mu=scale*(-actual),mu_effective=scale*closing;
    return rho*(c+2.0*mu)*mu_effective;
}

/* One compression/rotation switch, evaluated from the velocity fit.
 * Geometric Voronoi volume changes are deliberately not called div(v).
 * The reconstruction residual still detects a checkerboard with zero fit.
 */
GFS_CORE_INLINE double gfs_variational_compression_fraction(
        double divi,double curli,double divj,double curlj)
{
    double compression=fmax(0.0,-divi)+fmax(0.0,-divj);
    double rotation=fabs(curli)+fabs(curlj);
    return compression+rotation>0?compression/(compression+rotation):0;
}

#undef GFS_CORE_INLINE
#endif
