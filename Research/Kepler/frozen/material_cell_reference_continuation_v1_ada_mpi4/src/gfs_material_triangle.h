#ifndef GFS_MATERIAL_TRIANGLE_H
#define GFS_MATERIAL_TRIANGLE_H
#include <math.h>
#ifdef __CUDACC__
#define GFS_TRI_HD __host__ __device__
#else
#define GFS_TRI_HD
#endif
/* Shared H1 P1 reference mesh with two continuous INPUT-velocity residual
 * functions per connected component. Residual coefficients are ordinary
 * evolving kinematic unknowns, not an orbit constraint or residual source.
 * The caller supplies global DOF ownership and input function/derivatives.
 */
typedef struct {double x[3][2],grad[3][2],v0[3][2],area;} GfsMaterialTriangle;
GFS_TRI_HD static inline int gfs_material_triangle_init(GfsMaterialTriangle *t){
    double ax=t->x[1][0]-t->x[0][0],ay=t->x[1][1]-t->x[0][1];
    double bx=t->x[2][0]-t->x[0][0],by=t->x[2][1]-t->x[0][1],det=ax*by-ay*bx;
    if(!(det>0) || !isfinite(det))return 0;t->area=.5*det;
    t->grad[1][0]=by/det;t->grad[1][1]=-bx/det;t->grad[2][0]=-ay/det;t->grad[2][1]=ax/det;
    t->grad[0][0]=-t->grad[1][0]-t->grad[2][0];t->grad[0][1]=-t->grad[1][1]-t->grad[2][1];return 1;
}
GFS_TRI_HD static inline void gfs_material_triangle_shape(const GfsMaterialTriangle *t,
    double x,double y,const double v0[2],const double B0[4],const double scale[2],
    int enriched,double N[5],double Nx[5],double Ny[5]){
    double rx=x-t->x[0][0],ry=y-t->x[0][1];
    N[1]=t->grad[1][0]*rx+t->grad[1][1]*ry;N[2]=t->grad[2][0]*rx+t->grad[2][1]*ry;N[0]=1-N[1]-N[2];
    for(int a=0;a<3;a++){Nx[a]=t->grad[a][0];Ny[a]=t->grad[a][1];}
    if(!enriched)return;
    for(int k=0;k<2;k++){
        double r=v0[k]-t->v0[0][k]-N[1]*(t->v0[1][k]-t->v0[0][k])-N[2]*(t->v0[2][k]-t->v0[0][k]);
        N[3+k]=r/scale[k];Nx[3+k]=B0[2*k];Ny[3+k]=B0[2*k+1];
        for(int a=0;a<3;a++){Nx[3+k]-=t->grad[a][0]*t->v0[a][k];Ny[3+k]-=t->grad[a][1]*t->v0[a][k];}
        Nx[3+k]/=scale[k];Ny[3+k]/=scale[k];
    }
}
#undef GFS_TRI_HD
#endif
