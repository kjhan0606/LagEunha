#ifndef LAGFORCE_MOMENT_AFFINE_KICK_H
#define LAGFORCE_MOMENT_AFFINE_KICK_H
#include "moment_affine_recovery.h"
#ifdef __CUDACC__
#define MK_INLINE static __host__ __device__ __forceinline__
#else
#define MK_INLINE static inline
#endif
/* Experimental frozen-position affine acceleration substep, not full hydro.
 * Input q[17], F=dp/dt[2], J=dD/dt[4 row-major], dt, delta_Ehydro.
 * Find a(x)=F/m+A(x-c) matching F,J. Apply v'=v+dt*a exactly to
 * p,D,Pi, including the dt^2 acceleration variance. q,Q,m,Phi unchanged.
 * R is unchanged algebraically: adding an affine velocity cannot change
 * the residual orthogonal to affine functions. Caller supplies shared
 * midpoint pressure work, and must check thermal positivity after update.
 * No pressure/gravity closure or production integration is implied.
 */
MK_INLINE int mk_kick(const double *in,double *out) {
    double o[21];if(ma_recover(in,o))return 1;
    double m=in[0],fx=in[17],fy=in[18],dt=in[23];
    double sx=o[11],sxy=o[12],sy=o[13],det=sx*sy-sxy*sxy;
    double j00=in[19]-o[0]*fx,j01=in[20]-o[0]*fy;
    double j10=in[21]-o[1]*fx,j11=in[22]-o[1]*fy;
    double A[2][2]={{(j00*sy-j10*sxy)/det,(j10*sx-j00*sxy)/det},
                    {(j01*sy-j11*sxy)/det,(j11*sx-j01*sxy)/det}};
    double S[2][2]={{sx,sxy},{sxy,sy}},C[2][2]={{o[14],o[15]},{o[16],o[17]}};
    double acc[2]={fx/m,fy/m},p[2]={in[6],in[7]},Pi[2][2]={{in[12],in[13]},{in[13],in[14]}};
    for(int i=0;i<17;i++)out[i]=in[i];
    for(int i=0;i<2;i++)for(int j=i;j<2;j++) {
        double cross=p[i]*acc[j]+acc[i]*p[j],square=m*acc[i]*acc[j];
        for(int k=0;k<2;k++) {
            cross+=A[j][k]*C[k][i]+A[i][k]*C[k][j];
            for(int l=0;l<2;l++)square+=A[i][k]*S[k][l]*A[j][l];
        }
        double value=Pi[i][j]+dt*cross+dt*dt*square;
        out[(i==0 && j==0)?12:((i==1 && j==1)?14:13)]=value;
    }
    out[6]+=dt*fx;out[7]+=dt*fy;
    for(int k=0;k<4;k++)out[8+k]+=dt*in[19+k];
    out[15]+=in[24];
    return 0;
}
#undef MK_INLINE
#endif
