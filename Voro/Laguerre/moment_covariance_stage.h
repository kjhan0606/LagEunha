#ifndef LAGFORCE_MOMENT_COVARIANCE_STAGE_H
#define LAGFORCE_MOMENT_COVARIANCE_STAGE_H
#include "moment_state_layout.h"
#ifdef __CUDACC__
#define MCS_INLINE static __host__ __device__ __forceinline__
#else
#define MCS_INLINE static inline
#endif
/* One receiver incidence. extra20 = xiL[8],xiR[8],affineL[2],affineR[2].
 * side0 owns left i; side1 owns right j. Q stores four columns as [4][2]. */
MCS_INLINE int lmc_transport_node(const double *a,const double *extra,int side,double mass,
                                double *Kout,double *Qout){
    if(!a||!extra||!Kout||!Qout||(side!=0&&side!=1)||!isfinite(mass)||mass<=0.)return 1;
    for(int k=0;k<33;k++)if(!isfinite(a[k]))return 1;
    for(int k=0;k<20;k++)if(!isfinite(extra[k]))return 1;
    if(a[4]<0.)return 1;
    double K[4]={0.},Q[8];
    for(int s=0;s<4;s++){
        if(a[21+s]<0.||a[25+s]>0.)return 1;
        const double loss=side==0?a[21+s]:-a[25+s];
        const double incoming=side==0?-a[25+s]:a[21+s];
        const double x=extra[8*side+2*s],y=extra[8*side+2*s+1];
        const double weight=a[4]/mass*loss;
        K[0]+=weight*x*x;K[1]+=weight*x*y;K[3]+=weight*y*y;
        const double amplitude=sqrt(a[4]*incoming);
        for(int d=0;d<2;d++)Q[2*s+d]=amplitude*(a[(side==0?13:5)+2*s+d]-extra[16+2*side+d]);
    }
    K[2]=K[1];
    for(int k=0;k<4;k++)if(!isfinite(K[k]))return 3;
    for(int k=0;k<8;k++)if(!isfinite(Q[k]))return 3;
    for(int k=0;k<4;k++)Kout[k]=K[k];
    for(int k=0;k<8;k++)Qout[k]=Q[k];return 0;
}
/* Gravity point input: x[2],positive dm,vmean[2],acceleration[2]. */
MCS_INLINE int lmc_gravity_node(const double *a,double *out){
    if(!a||!out)return 1;
    for(int k=0;k<7;k++)if(!isfinite(a[k]))return 1;
    if(a[2]<=0.)return 1;
    double t[18]={0.};
    for(int i=0;i<2;i++){
        t[7+i]=a[2]*a[5+i];
        for(int j=0;j<2;j++)t[9+2*i+j]=a[2]*a[i]*a[5+j];
    }
    t[13]=2*a[2]*a[3]*a[5];t[16]=2*a[2]*a[4]*a[6];
    t[14]=t[15]=a[2]*(a[3]*a[6]+a[5]*a[4]);
    t[17]=a[2]*(a[3]*a[5]+a[4]*a[6]);
    for(int k=0;k<18;k++)if(!isfinite(t[k]))return 3;
    for(int k=0;k<18;k++)out[k]=t[k];return 0;
}
/* Symmetric SPD 2x2 S solve S X=C, output B=X^T. */
MCS_INLINE int lmc_affine_B(const double *r,double *B){
    if(!r||!B)return 1;
    for(int k=0;k<18;k++)if(!isfinite(r[k]))return 1;
    if(r[0]<=0.||r[4]!=r[5])return 1;
    const double m=r[0],s00=r[3]-r[1]*r[1]/m,s01=r[4]-r[1]*r[2]/m,s11=r[6]-r[2]*r[2]/m;
    const double scale=fmax(fmax(fabs(s00),fabs(s01)),fabs(s11));
    if(!isfinite(scale))return 3;if(scale<=0.||s00<=0.)return 2;
    const double l0=sqrt(s00/scale),l10=(s01/scale)/l0,d=s11/scale-l10*l10;
    if(!isfinite(d))return 3;if(d<=0.)return 2;
    const double l1=sqrt(d);double tmp[4];
    for(int j=0;j<2;j++){
        const double c0=r[9+j]-r[1]*r[7+j]/m,c1=r[11+j]-r[2]*r[7+j]/m;
        const double y0=(c0/scale)/l0,y1=(c1/scale-l10*y0)/l1;
        tmp[2*j+1]=y1/l1;tmp[2*j]=(y0-l10*tmp[2*j+1])/l0;
    }
    for(int k=0;k<4;k++)if(!isfinite(tmp[k]))return 3;
    for(int k=0;k<4;k++)B[k]=tmp[k];return 0;
}
/* Covariance acceleration source and unique symmetric Sylvester generator.
 * Outputs B[4],G[4],Ag[4], diagnostics[3]: G antisymmetry before averaging,
 * normwise Sylvester backward residual, cond_inf of scaled 3x3 operator.
 * Strict SPD R=L L^T, no clipping or regularization, transactional outputs.
 * R is formed before scaling; extreme L overflow/underflow is rejected.
 * Inputs and output arrays must not overlap one another.
 */
MCS_INLINE int lmc_gravity_covariance(const double *r,const double *g,const double *L,
                                     double *Bout,double *Gout,double *Agout,double *diagout){
    if(!r||!g||!L||!Bout||!Gout||!Agout||!diagout)return 1;
    for(int k=0;k<18;k++)if(!isfinite(g[k]))return 1;
    for(int k=0;k<7;k++)if(g[k]!=0.)return 1;
    if(g[14]!=g[15])return 1;
    for(int k=0;k<4;k++)if(!isfinite(L[k]))return 1;
    double B[4];int status=lmc_affine_B(r,B);if(status)return status;
    double dc[4],dt[4],G[4];const double m=r[0];
    for(int i=0;i<2;i++)for(int j=0;j<2;j++){
        dc[2*i+j]=g[9+2*i+j]-r[1+i]*g[7+j]/m;
        dt[2*i+j]=g[13+2*i+j]-(g[7+i]*r[7+j]+r[7+i]*g[7+j])/m;
    }
    for(int i=0;i<2;i++)for(int j=0;j<2;j++){
        G[2*i+j]=dt[2*i+j]-B[2*i]*dc[j]-B[2*i+1]*dc[2+j]
            -dc[i]*B[2*j]-dc[2+i]*B[2*j+1];
    }
    const double antisym=G[1]-G[2];G[1]=G[2]=.5*(G[1]+G[2]);
    const double R[4]={L[0]*L[0]+L[1]*L[1],L[0]*L[2]+L[1]*L[3],
                       L[0]*L[2]+L[1]*L[3],L[2]*L[2]+L[3]*L[3]};
    const double scale=fmax(R[0],R[3]);if(!isfinite(scale))return 3;if(scale<=0.)return 2;
    const double aa=R[0]/scale,bb=R[1]/scale,cc=R[3]/scale,rt=sqrt(2.);
    if(aa<=0.)return 2;
    const double r10=bb/sqrt(aa),rrem=cc-r10*r10;
    if(!isfinite(rrem))return 3;if(rrem<=0.)return 2;
    // SPD in orthonormal symmetric-matrix coordinates (xx,sqrt(2)xy,yy).
    const double M[9]={2*aa,rt*bb,0.,rt*bb,aa+cc,rt*bb,0.,rt*bb,2*cc};
    double chol[9]={0.};
    for(int i=0;i<3;i++)for(int j=0;j<=i;j++){
        double v=M[3*i+j];for(int k=0;k<j;k++)v-=chol[3*i+k]*chol[3*j+k];
        if(!isfinite(v))return 3;
        if(i==j){if(v<=0.)return 2;chol[3*i+j]=sqrt(v);}else chol[3*i+j]=v/chol[3*j+j];
    }
    const double rhs[3]={G[0]/scale,rt*G[1]/scale,G[3]/scale};double z[3];
    for(int i=0;i<3;i++){double v=rhs[i];for(int j=0;j<i;j++)v-=chol[3*i+j]*z[j];z[i]=v/chol[3*i+i];}
    for(int i=2;i>=0;i--){double v=z[i];for(int j=i+1;j<3;j++)v-=chol[3*j+i]*z[j];z[i]=v/chol[3*i+i];}
    double Ag[4]={z[0],z[1]/rt,z[1]/rt,z[2]},inv[9]={0.};
    for(int column=0;column<3;column++){
        double v[3];for(int i=0;i<3;i++){double b=i==column?1.:0.;for(int j=0;j<i;j++)b-=chol[3*i+j]*v[j];v[i]=b/chol[3*i+i];}
        for(int i=2;i>=0;i--){double b=v[i];for(int j=i+1;j<3;j++)b-=chol[3*j+i]*v[j];v[i]=b/chol[3*i+i];}
        for(int i=0;i<3;i++)inv[3*i+column]=v[i];
    }
    double mn=0.,in=0.,rn=0.,zn=0.,bn=0.;
    for(int i=0;i<3;i++){
        double mm=0.,ii=0.,res=-rhs[i];for(int j=0;j<3;j++){mm+=fabs(M[3*i+j]);ii+=fabs(inv[3*i+j]);res+=M[3*i+j]*z[j];}
        mn=fmax(mn,mm);in=fmax(in,ii);rn=fmax(rn,fabs(res));zn=fmax(zn,fabs(z[i]));bn=fmax(bn,fabs(rhs[i]));
    }
    const double denom=mn*zn+bn;double diagnostics[3]={antisym,denom>0.?rn/denom:rn,mn*in};
    for(int k=0;k<4;k++)if(!isfinite(B[k])||!isfinite(G[k])||!isfinite(Ag[k]))return 3;
    for(int k=0;k<3;k++)if(!isfinite(diagnostics[k]))return 3;
    for(int k=0;k<4;k++){Bout[k]=B[k];Gout[k]=G[k];Agout[k]=Ag[k];}
    for(int k=0;k<3;k++)diagout[k]=diagnostics[k];return 0;
}
#undef MCS_INLINE
#endif
