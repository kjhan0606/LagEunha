#ifndef LAGFORCE_MOMENT_HYDRO_CELL_SOURCE_H
#define LAGFORCE_MOMENT_HYDRO_CELL_SOURCE_H
#include "moment_state_layout.h"
#ifdef __CUDACC__
#define LHC_INLINE static __host__ __device__ __forceinline__
#else
#define LHC_INLINE static inline
#endif
/* Input raw18, pressure force F[2], centered face pressure moment Hface[4],
 * PV=P*area. Source contains kinematic volume and compatible pressure Pi;
 * pressure momentum/work ALREADY belong to the face flux. No double count.
 * diagnostic[8]: A row-major, a0[2], scaled backward residual, cond_inf(S).
 * Return 0 success,1 invalid input,2 inadmissible spatial covariance,
 * 3 nonfinite intermediate/output. No clipping. Both outputs transactional.
 * Output arrays must not overlap each other; inputs may alias output arrays.
 */
LHC_INLINE int lmh_cell_source(const double *r,const double *F,const double *hf,
                              double pv,double *out,double *diagnostic){
    if(!r || !F || !hf || !out || !diagnostic || !isfinite(pv))return 1;
    for(int k=0;k<18;k++)if(!isfinite(r[k]))return 1;
    for(int k=0;k<2;k++)if(!isfinite(F[k]))return 1;
    for(int k=0;k<4;k++)if(!isfinite(hf[k]))return 1;
    if(r[0]<=0. || r[4]!=r[5] || r[14]!=r[15])return 1;
    const double m=r[0],c[2]={r[1]/m,r[2]/m};
    const double S[4]={r[3]-(r[1]*r[1])/m,r[4]-(r[1]*r[2])/m,
                       r[5]-(r[2]*r[1])/m,r[6]-(r[2]*r[2])/m};
    // Equal raw Q offdiagonals can acquire unequal last bits in q*c.
    // Form one common qx*qy/m term, without averaging the stored raw tensor.
    const double sxy=r[4]-(r[1]*r[2])/m;
    double scale=fmax(fmax(fabs(S[0]),fabs(S[3])),fabs(sxy));
    if(!isfinite(scale) || !isfinite(sxy) || !isfinite(c[0]) || !isfinite(c[1]))return 3;
    if(scale<=0.)return 2;
    const double s00=S[0]/scale,s01=sxy/scale,s11=S[3]/scale;
    if(s00<=0.)return 2;
    const double l00=sqrt(s00),l10=s01/l00,rem=s11-l10*l10;
    if(!isfinite(l10) || !isfinite(rem))return 3;
    if(rem<=0.)return 2;
    const double l11=sqrt(rem),H[4]={hf[0]+pv,hf[1],hf[2],hf[3]+pv};
    double A[4];
    for(int j=0;j<2;j++){
        const double y0=(H[j]/scale)/l00;
        const double y1=(H[2+j]/scale-l10*y0)/l11;
        const double b1=y1/l11,b0=(y0-l10*b1)/l00;
        A[2*j]=b0;A[2*j+1]=b1;
    }
    const double a0[2]={F[0]/m-A[0]*c[0]-A[1]*c[1],
                         F[1]/m-A[2]*c[0]-A[3]*c[1]};
    double src[18]={0.},diag[8];
    src[1]=r[7];src[2]=r[8];src[3]=2*r[9];
    src[4]=src[5]=r[10]+r[11];src[6]=2*r[12];
    src[9]=r[13]+pv;src[10]=r[14];src[11]=r[15];src[12]=r[16]+pv;
    for(int i=0;i<2;i++)for(int j=i;j<2;j++){
        const double value=a0[i]*r[7+j]+r[7+i]*a0[j]
            +A[2*i]*r[9+j]+A[2*i+1]*r[11+j]
            +r[9+i]*A[2*j]+r[11+i]*A[2*j+1];
        src[13+2*i+j]=value;src[13+2*j+i]=value;
    }
    double residual=0.,anorm=0.,hnorm=0.;
    for(int i=0;i<2;i++){
        double rr=0.,hh=0.;
        for(int j=0;j<2;j++){
            const double v=(i==0?s00:s01)*A[2*j]+(i==0?s01:s11)*A[2*j+1]-H[2*i+j]/scale;
            rr+=fabs(v);hh+=fabs(H[2*i+j]/scale);
        }
        residual=fmax(residual,rr);hnorm=fmax(hnorm,hh);
        anorm=fmax(anorm,fabs(A[i])+fabs(A[2+i])); // norm of A^T
    }
    const double snorm=fmax(fabs(s00)+fabs(s01),fabs(s01)+fabs(s11));
    const double denom=snorm*anorm+hnorm;
    const double det=(l00*l11)*(l00*l11);
    for(int k=0;k<4;k++)diag[k]=A[k];
    diag[4]=a0[0];diag[5]=a0[1];diag[6]=denom>0.?residual/denom:residual;
    diag[7]=snorm*fmax(fabs(s11)+fabs(s01),fabs(s00)+fabs(s01))/det;
    for(int k=0;k<18;k++)if(!isfinite(src[k]))return 3;
    for(int k=0;k<8;k++)if(!isfinite(diag[k]))return 3;
    for(int k=0;k<18;k++)out[k]=src[k];
    for(int k=0;k<8;k++)diagnostic[k]=diag[k];
    return 0;
}
#undef LHC_INLINE
#endif
