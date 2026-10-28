#ifndef LAGFORCE_AFFINE_RELAXATION_RATE_H
#define LAGFORCE_AFFINE_RELAXATION_RATE_H
#include <math.h>
#ifdef __CUDACC__
#define ARR_INLINE static __host__ __device__ __forceinline__
#else
#define ARR_INLINE static inline
#endif
/* Experimental 17-moment affine contact relaxation. No advection/gravity.
 * gamma includes the positive face quadrature weight. Unit tangent required.
 * Does not certify spatial realizability, choose viscosity or repair a state.
 * Outputs remain untouched on failure; input/output buffers must not alias.
 */
ARR_INLINE int arr_velocity(const double *q,const double *f,double *v) {
    for(int k=0;k<17;k++)if(!isfinite(q[k]))return 1;
    if(!(q[0]>0.0) || !(q[15]-0.5*(q[12]+q[14])>0.0))return 1;
    const double cx=q[1]/q[0],cy=q[2]/q[0],vx=q[6]/q[0],vy=q[7]/q[0];
    const double xx=q[3]-cx*q[1],xy=q[4]-cx*q[2],yy=q[5]-cy*q[2];
    const double det=xx*yy-xy*xy;
    if(!(xx>0.0 && yy>0.0 && det>0.0) || !isfinite(det))return 2;
    const double rx=f[0]-cx,ry=f[1]-cy;
    const double ax=(yy*rx-xy*ry)/det,ay=(xx*ry-xy*rx)/det;
    v[0]=vx+(q[8]-q[1]*vx)*ax+(q[10]-q[2]*vx)*ay;
    v[1]=vy+(q[9]-q[1]*vy)*ax+(q[11]-q[2]*vy)*ay;
    return isfinite(v[0]) && isfinite(v[1]) ? 0 : 1;
}
ARR_INLINE int arr_pair(const double *a,const double *b,const double *f,
                        const double *t,double gamma,double *da,double *db) {
    if(!isfinite(gamma) || gamma<0.0)return 1;
    for(int k=0;k<2;k++)if(!isfinite(f[k]) || !isfinite(t[k]))return 1;
    if(fabs(t[0]*t[0]+t[1]*t[1]-1.0)>1e-12)return 1;
    double va[2],vb[2];int code=arr_velocity(a,f,va);if(code)return code;
    code=arr_velocity(b,f,vb);if(code)return code;
    const double s=(va[0]-vb[0])*t[0]+(va[1]-vb[1])*t[1];
    const double F[2]={-gamma*s*t[0],-gamma*s*t[1]},heat=gamma*s*s;
    double out[2][17]={{0.0}};
    for(int side=0;side<2;side++) {
        const double sign=side ? -1.0 : 1.0;
        const double fx=sign*F[0],fy=sign*F[1];const double *v=side ? vb : va;
        double *d=out[side];d[6]=fx;d[7]=fy;
        d[8]=f[0]*fx;d[9]=f[0]*fy;d[10]=f[1]*fx;d[11]=f[1]*fy;
        d[12]=2.0*v[0]*fx;d[13]=v[0]*fy+v[1]*fx;d[14]=2.0*v[1]*fy;
        d[15]=v[0]*fx+v[1]*fy+0.5*heat;
        for(int k=0;k<17;k++)if(!isfinite(d[k]))return 1;
    }
    for(int k=0;k<17;k++){da[k]=out[0][k];db[k]=out[1][k];}
    return 0;
}
#undef ARR_INLINE
#endif
