#ifndef GFS_MATERIAL_MAP_H
#define GFS_MATERIAL_MAP_H
#include <math.h>
#ifdef __CUDACC__
#define GFS_MAT_HD __host__ __device__
#else
#define GFS_MAT_HD
#endif
#define GFS_MATERIAL_N 5
/* Compatible material-map algebra only. Reference basis contains constants
 * and coordinates; x(X),v(X) use the SAME basis, fixed mass rho0 dX, and mass
 * matrix. e is physical specific internal energy, not total-minus-affine K.
 * No Euler remap, timestep, viscosity, or production solver is implied.
 */
typedef struct {
    int n;
    double x[2][GFS_MATERIAL_N],v[2][GFS_MATERIAL_N];
} GfsMaterialMap;
typedef struct {
    double position[2],velocity[2],F[4],J,grad_v[4];
} GfsMaterialPoint;
GFS_MAT_HD static inline int gfs_material_point(const GfsMaterialMap *map,
    const double *N,const double *Nx,const double *Ny,GfsMaterialPoint *point){
    if(map->n<1 || map->n>GFS_MATERIAL_N)return 0;
    for(int k=0;k<2;k++){point->position[k]=0;point->velocity[k]=0;}
    for(int k=0;k<4;k++)point->F[k]=point->grad_v[k]=0;
    double vX[4]={0};
    for(int a=0;a<map->n;a++)for(int k=0;k<2;k++){
        point->position[k]+=map->x[k][a]*N[a];point->velocity[k]+=map->v[k][a]*N[a];
        point->F[2*k]+=map->x[k][a]*Nx[a];point->F[2*k+1]+=map->x[k][a]*Ny[a];
        vX[2*k]+=map->v[k][a]*Nx[a];vX[2*k+1]+=map->v[k][a]*Ny[a];
    }
    point->J=point->F[0]*point->F[3]-point->F[1]*point->F[2];
    if(!(point->J>0) || !isfinite(point->J))return 0;
    double inv[4]={point->F[3]/point->J,-point->F[1]/point->J,-point->F[2]/point->J,point->F[0]/point->J};
    for(int k=0;k<2;k++){
        point->grad_v[2*k]=vX[2*k]*inv[0]+vX[2*k+1]*inv[2];
        point->grad_v[2*k+1]=vX[2*k]*inv[1]+vX[2*k+1]*inv[3];
    }return 1;
}
/* sigmaxx,sigmaxy,sigmayy is a symmetric Cauchy stress, including -P I.
 * Same-point traction, virtual work and physical heat are paired exactly.
 * Numerical shock stress must be derived separately; none is invented here.
 */
GFS_MAT_HD static inline int gfs_material_stress_point(const GfsMaterialMap *map,
    const double *N,const double *Nx,const double *Ny,double rho0,double reference_weight,
    const double sigma[3],double gx,double gy,double force[2][GFS_MATERIAL_N],
    double *internal_rate,double *gravity_power){
    GfsMaterialPoint p;if(!(rho0>0) || !(reference_weight>=0) ||
        !gfs_material_point(map,N,Nx,Ny,&p))return 0;
    double inv[4]={p.F[3]/p.J,-p.F[1]/p.J,-p.F[2]/p.J,p.F[0]/p.J};
    double dV=reference_weight*p.J,dm=reference_weight*rho0;
    for(int a=0;a<map->n;a++){
        double dx=Nx[a]*inv[0]+Ny[a]*inv[2],dy=Nx[a]*inv[1]+Ny[a]*inv[3];
        force[0][a]+= -dV*(sigma[0]*dx+sigma[1]*dy)+dm*N[a]*gx;
        force[1][a]+= -dV*(sigma[1]*dx+sigma[2]*dy)+dm*N[a]*gy;
    }
    *internal_rate+=dV*(sigma[0]*p.grad_v[0]+sigma[1]*(p.grad_v[1]+p.grad_v[2])+sigma[2]*p.grad_v[3]);
    *gravity_power+=dm*(p.velocity[0]*gx+p.velocity[1]*gy);return 1;
}
GFS_MAT_HD static inline void gfs_material_mass_point(int n,const double *N,
    double dm,double M[GFS_MATERIAL_N][GFS_MATERIAL_N]){
    for(int a=0;a<n;a++)for(int b=0;b<n;b++)M[a][b]+=dm*N[a]*N[b];
}
GFS_MAT_HD static inline double gfs_material_kinetic(const GfsMaterialMap *map,
    const double M[GFS_MATERIAL_N][GFS_MATERIAL_N]){
    double K=0;for(int k=0;k<2;k++)for(int a=0;a<map->n;a++)for(int b=0;b<map->n;b++)
        K+=.5*map->v[k][a]*M[a][b]*map->v[k][b];return K;
}
typedef struct {
    double mass,center[2],momentum[2],kinetic,internal_kinetic;
    double intrinsic_angular,polar_inertia,rigid_rotation_kinetic,resolved_remainder;
} GfsMaterialMoments;
/* Cell-local basis has N0=1. These are moments of the ACTUAL represented
 * velocity/mass field, not an independently transported spin/kinetic tensor.
 * The rigid rotational projection and remaining kinetic energy are diagnostics
 * only. Neither is added to or removed from physical internal energy e.
 */
GFS_MAT_HD static inline int gfs_material_moments(const GfsMaterialMap *map,
    const double M[GFS_MATERIAL_N][GFS_MATERIAL_N],GfsMaterialMoments *out){
    out->mass=M[0][0];if(!(out->mass>0))return 0;
    double q[2]={0},Q[2]={0},Dxy=0,Dyx=0;out->momentum[0]=out->momentum[1]=0;
    for(int k=0;k<2;k++)for(int a=0;a<map->n;a++){
        q[k]+=M[0][a]*map->x[k][a];out->momentum[k]+=M[0][a]*map->v[k][a];
        for(int b=0;b<map->n;b++)Q[k]+=map->x[k][a]*M[a][b]*map->x[k][b];
    }
    for(int a=0;a<map->n;a++)for(int b=0;b<map->n;b++){
        Dxy+=map->x[0][a]*M[a][b]*map->v[1][b];Dyx+=map->x[1][a]*M[a][b]*map->v[0][b];
    }
    out->center[0]=q[0]/out->mass;out->center[1]=q[1]/out->mass;
    out->kinetic=gfs_material_kinetic(map,M);
    out->internal_kinetic=out->kinetic-.5*(out->momentum[0]*out->momentum[0]+out->momentum[1]*out->momentum[1])/out->mass;
    out->intrinsic_angular=Dxy-Dyx-(q[0]*out->momentum[1]-q[1]*out->momentum[0])/out->mass;
    out->polar_inertia=Q[0]+Q[1]-(q[0]*q[0]+q[1]*q[1])/out->mass;
    if(!(out->polar_inertia>0))return 0;
    out->rigid_rotation_kinetic=.5*out->intrinsic_angular*out->intrinsic_angular/out->polar_inertia;
    out->resolved_remainder=out->internal_kinetic-out->rigid_rotation_kinetic;
    return isfinite(out->resolved_remainder);
}
#undef GFS_MAT_HD
#endif
