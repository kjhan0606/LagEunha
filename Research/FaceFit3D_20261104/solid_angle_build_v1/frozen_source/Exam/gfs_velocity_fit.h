#ifndef GFS_VELOCITY_FIT_H
#define GFS_VELOCITY_FIT_H
#include <float.h>
#include <math.h>
#ifdef __CUDACC__
#define GFS_FIT_INLINE static __host__ __device__ __forceinline__
#else
#define GFS_FIT_INLINE static inline
#endif

/* A single small least-squares velocity fit for dimensions 1, 2, and 3.
 * Samples are generator displacements and velocity DIFFERENCES.
 * The fit is affine exact, rotation covariant, and Galilean invariant.
 * It is used for numerical viscosity, not primitive-state extrapolation.
 * No fit regularization or silent zero-gradient fallback is permitted.
 * Native hydrodynamics uses face solid-angle weights A/distance^(dim-1).
 * Thus a vanishing face cannot contribute a finite angular sample. The
 * 1D unit-area face reduces exactly to the original fit.
 */
typedef struct { double M[3][3],R[3][3];int dim,count; } GfsVelocityFit;
GFS_FIT_INLINE int gfs_velocity_fit_init(GfsVelocityFit *f,int dim)
{
    if(!f||dim<1||dim>3)return 0;
    f->dim=dim;f->count=0;
    for(int a=0;a<3;a++)for(int b=0;b<3;b++)f->M[a][b]=f->R[a][b]=0;
    return 1;
}
GFS_FIT_INLINE int gfs_velocity_fit_add_weight(GfsVelocityFit *f,
        const double r[3],const double dv[3],double weight)
{
    if(!(weight>=0)||!isfinite(weight))return 0;
    if(weight==0)return 1;
    double d2=0;
    for(int a=0;a<f->dim;a++){
        if(!isfinite(r[a])||!isfinite(dv[a]))return 0;
        d2+=r[a]*r[a];
    }
    if(!(d2>0)||!isfinite(d2))return 0;
    for(int a=0;a<f->dim;a++)for(int b=0;b<f->dim;b++){
        f->M[a][b]+=weight*(r[a]*r[b]/d2);
        f->R[a][b]+=weight*(dv[a]*r[b]/d2);
    }
    f->count++;return 1;
}
GFS_FIT_INLINE int gfs_velocity_fit_add(GfsVelocityFit *f,
        const double r[3],const double dv[3])
{
    return gfs_velocity_fit_add_weight(f,r,dv,1);
}
GFS_FIT_INLINE int gfs_velocity_fit_add_face(GfsVelocityFit *f,
        const double r[3],const double dv[3],double area)
{
    if(!(area>=0)||!isfinite(area))return 0;
    double d2=0;for(int a=0;a<f->dim;a++)d2+=r[a]*r[a];
    if(!(d2>0)||!isfinite(d2))return 0;
    double weight=area;
    if(f->dim==2)weight/=sqrt(d2);
    else if(f->dim==3)weight/=d2;
    return gfs_velocity_fit_add_weight(f,r,dv,weight);
}
GFS_FIT_INLINE int gfs_velocity_fit_solve(const GfsVelocityFit *f,
        double G[3][3],double *condition_bound)
{
    const int d=f->dim;
    double L[3][3]={{0}},inverse[3][3]={{0}},trace=0;
    if(f->count<d)return 0;
    for(int a=0;a<d;a++)trace+=f->M[a][a];
    if(!(trace>0)||!isfinite(trace))return 0;
    for(int a=0;a<d;a++)for(int b=0;b<=a;b++){
        double q=f->M[a][b];
        for(int k=0;k<b;k++)q-=L[a][k]*L[b][k];
        if(a==b){
            if(!(q>64*DBL_EPSILON*trace)||!isfinite(q))return 0;
            L[a][b]=sqrt(q);
        }else L[a][b]=q/L[b][b];
    }
    for(int column=0;column<d;column++){
        double y[3]={0},x[3]={0};
        for(int a=0;a<d;a++){
            double q=a==column?1:0;
            for(int k=0;k<a;k++)q-=L[a][k]*y[k];
            y[a]=q/L[a][a];
        }
        for(int a=d-1;a>=0;a--){
            double q=y[a];
            for(int k=a+1;k<d;k++)q-=L[k][a]*x[k];
            x[a]=q/L[a][a];inverse[a][column]=x[a];
        }
    }
    double normM2=0,normI2=0;
    for(int a=0;a<d;a++)for(int b=0;b<d;b++){
        normM2+=f->M[a][b]*f->M[a][b];
        normI2+=inverse[a][b]*inverse[a][b];
    }
    double condition=sqrt(normM2*normI2);
    if(!isfinite(condition)||condition>1/sqrt(DBL_EPSILON))return 0;
    if(condition_bound)*condition_bound=condition;
    for(int a=0;a<3;a++)for(int b=0;b<3;b++){
        G[a][b]=0;
        for(int k=0;k<d&&a<d&&b<d;k++)G[a][b]+=f->R[a][k]*inverse[k][b];
        if(!isfinite(G[a][b]))return 0;
    }
    return 1;
}
#undef GFS_FIT_INLINE
#endif
