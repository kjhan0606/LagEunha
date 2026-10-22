#ifndef GFS_MATERIAL_P3_H
#define GFS_MATERIAL_P3_H
#include "gfs_material_triangle.h"
#ifdef __CUDACC__
#define GFS_P3_HD __host__ __device__
#else
#define GFS_P3_HD
#endif
#define GFS_MATERIAL_P3_POLY 10
#define GFS_MATERIAL_P3_N 12
/* Complete cubic Lagrange kinematics on a fixed material triangle, plus two
 * initial-velocity interpolation residuals. No source/gravity enrichment.
 * Global ownership/trace constraints and evolution are caller responsibilities.
 */
typedef struct {GfsMaterialTriangle geometry;double v0[10][2],remainder[10][2];} GfsMaterialP3Triangle;
GFS_P3_HD static inline void gfs_material_p3_index(int a,int index[3]){
    const int order[10][3]={{3,0,0},{0,3,0},{0,0,3},{2,1,0},{1,2,0},{0,2,1},{0,1,2},{1,0,2},{2,0,1},{1,1,1}};
    for(int k=0;k<3;k++)index[k]=order[a][k];
}
GFS_P3_HD static inline void gfs_material_p3_node(const GfsMaterialP3Triangle *t,int a,double x[2]){
    int index[3];gfs_material_p3_index(a,index);for(int k=0;k<2;k++){x[k]=0;for(int j=0;j<3;j++)x[k]+=index[j]*t->geometry.x[j][k]/3.;}
}
GFS_P3_HD static inline void gfs_material_p3_factor(int n,double l,double *p,double *d){
    if(n==0){*p=1;*d=0;}else if(n==1){*p=3*l;*d=3;}
    else if(n==2){*p=1.5*l*(3*l-1);*d=9*l-1.5;}
    else{*p=.5*l*(3*l-1)*(3*l-2);*d=(13.5*l-9)*l+1;}
}
GFS_P3_HD static inline void gfs_material_p3_polynomial(const GfsMaterialP3Triangle *t,
    double x,double y,double N[10],double Nx[10],double Ny[10]){
    const GfsMaterialTriangle *g=&t->geometry;double rx=x-g->x[0][0],ry=y-g->x[0][1],l[3];
    l[1]=g->grad[1][0]*rx+g->grad[1][1]*ry;l[2]=g->grad[2][0]*rx+g->grad[2][1]*ry;l[0]=1-l[1]-l[2];
    for(int a=0;a<10;a++){
        int index[3];double p[3],d[3];gfs_material_p3_index(a,index);
        for(int k=0;k<3;k++)gfs_material_p3_factor(index[k],l[k],&p[k],&d[k]);
        N[a]=p[0]*p[1]*p[2];Nx[a]=Ny[a]=0;
        for(int k=0;k<3;k++){double v=d[k]*p[(k+1)%3]*p[(k+2)%3];Nx[a]+=v*g->grad[k][0];Ny[a]+=v*g->grad[k][1];}
    }
}
GFS_P3_HD static inline void gfs_material_p3_shape(const GfsMaterialP3Triangle *t,
    double x,double y,const double v0[2],const double B0[4],const double scale[2],
    double N[12],double Nx[12],double Ny[12]){
    gfs_material_p3_polynomial(t,x,y,N,Nx,Ny);
    for(int k=0;k<2;k++){
        double r=v0[k]-t->v0[0][k],dx=B0[2*k],dy=B0[2*k+1];
        for(int a=1;a<10;a++){double dv=t->v0[a][k]-t->v0[0][k];r-=N[a]*dv;dx-=Nx[a]*dv;dy-=Ny[a]*dv;}
        N[10+k]=r/scale[k];Nx[10+k]=dx/scale[k];Ny[10+k]=dy/scale[k];
    }
}
/* Same residual from any accurate input remainder q=v0-T3(v0):
 * q-I3(q)=v0-I3(v0). T3 is an INITIAL-input polynomial; no force correction.
 */
GFS_P3_HD static inline void gfs_material_p3_shape_remainder(const GfsMaterialP3Triangle *t,
    double x,double y,const double q[2],const double Dq[4],const double scale[2],
    double N[12],double Nx[12],double Ny[12]){
    gfs_material_p3_polynomial(t,x,y,N,Nx,Ny);
    for(int k=0;k<2;k++){
        double r=q[k]-t->remainder[0][k],dx=Dq[2*k],dy=Dq[2*k+1];
        for(int a=1;a<10;a++){double dq=t->remainder[a][k]-t->remainder[0][k];r-=N[a]*dq;dx-=Nx[a]*dq;dy-=Ny[a]*dq;}
        N[10+k]=r/scale[k];Nx[10+k]=dx/scale[k];Ny[10+k]=dy/scale[k];
    }
}
#undef GFS_P3_HD
#endif
