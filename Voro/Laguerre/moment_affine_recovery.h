#ifndef LAGFORCE_MOMENT_AFFINE_RECOVERY_H
#define LAGFORCE_MOMENT_AFFINE_RECOVERY_H
#ifdef __CUDACC__
#define MA_INLINE static __host__ __device__ __forceinline__
#else
#define MA_INLINE static inline
#endif
/* Diagnostic implementation of the 2026-10-09 affine moment projection.
 * Input: 17 raw moments from weight_remap_moments.h.
 * Output: c[2],vbar[2],B[4 row-major],R[xx,xy,yy],S[xx,xy,yy],
 * C[4 row-major],U,K_affine,K_residual (21 doubles).
 * S=Q-qq'/m, C=D-qp'/m, B=C' S^-1, R=Pi-pp'/m-B C.
 * No clipping of R, no conversion of residual motion to heat, no force.
 * Positive S and U are required; R realizability is checked by the caller.
 */
MA_INLINE int ma_recover(const double *q,double *o) {
    double m=q[0]; if(!(m>0.0)) return 1;
    double cx=q[1]/m,cy=q[2]/m,vx=q[6]/m,vy=q[7]/m;
    double sx=q[3]-q[1]*cx,sxy=q[4]-q[1]*cy,sy=q[5]-q[2]*cy;
    double det=sx*sy-sxy*sxy,U=q[15]-0.5*(q[12]+q[14]);
    if(!(sx>0.0 && sy>0.0 && det>0.0 && U>0.0)) return 1;
    double c00=q[8]-q[1]*vx,c01=q[9]-q[1]*vy;
    double c10=q[10]-q[2]*vx,c11=q[11]-q[2]*vy;
    double b00=(c00*sy-c10*sxy)/det,b01=(c10*sx-c00*sxy)/det;
    double b10=(c01*sy-c11*sxy)/det,b11=(c11*sx-c01*sxy)/det;
    double ax=b00*c00+b01*c10,ay=b10*c01+b11*c11;
    double rx=q[12]-q[6]*vx-ax,ry=q[14]-q[7]*vy-ay;
    double rxy=q[13]-q[6]*vy-0.5*(b00*c01+b01*c11+b10*c00+b11*c10);
    o[0]=cx;o[1]=cy;o[2]=vx;o[3]=vy;
    o[4]=b00;o[5]=b01;o[6]=b10;o[7]=b11;
    o[8]=rx;o[9]=rxy;o[10]=ry;
    o[11]=sx;o[12]=sxy;o[13]=sy;
    o[14]=c00;o[15]=c01;o[16]=c10;o[17]=c11;
    o[18]=U;o[19]=0.5*(ax+ay);o[20]=0.5*(rx+ry);
    return 0;
}
#undef MA_INLINE
#endif
