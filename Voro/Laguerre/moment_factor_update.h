#ifndef LAGFORCE_MOMENT_FACTOR_UPDATE_H
#define LAGFORCE_MOMENT_FACTOR_UPDATE_H
#include "moment_covariance_stage.h"
#ifdef __CUDACC__
#define MF_INLINE static __host__ __device__ __forceinline__
#else
#define MF_INLINE static inline
#endif
MF_INLINE void lmf_multiply(const double *a,const double *b,double *out){
    double t[4];for(int i=0;i<2;i++)for(int j=0;j<2;j++)t[2*i+j]=a[2*i]*b[j]+a[2*i+1]*b[2+j];
    for(int k=0;k<4;k++)out[k]=t[k];
}
MF_INLINE double lmf_norm(const double *a){return fmax(fabs(a[0])+fabs(a[1]),fabs(a[2])+fabs(a[3]));}
/* exp(X) for general real 2x2 X. Explicit bounded scaling/Taylor18/squaring.
 * diagnostics = scaling count, scaled norm, last Taylor-term norm. */
MF_INLINE int lmf_exponential(const double *X,double *out,double *diagnostics){
    for(int k=0;k<4;k++)if(!isfinite(X[k]))return 1;
    double norm=lmf_norm(X);if(!isfinite(norm))return 3;
    int s=0;while(norm>.5 && s<60){norm*=.5;s++;}if(norm>.5)return 4;
    double scaled[4],term[4]={1.,0.,0.,1.},value[4]={1.,0.,0.,1.};
    for(int k=0;k<4;k++)scaled[k]=ldexp(X[k],-s);
    for(int n=1;n<=18;n++){
        double next[4];lmf_multiply(term,scaled,next);
        for(int k=0;k<4;k++){term[k]=next[k]/n;value[k]+=term[k];}
    }
    const double last=lmf_norm(term);
    for(int n=0;n<s;n++)lmf_multiply(value,value,value);
    for(int k=0;k<4;k++)if(!isfinite(value[k]))return 3;
    for(int k=0;k<4;k++)out[k]=value[k];
    diagnostics[0]=s;diagnostics[1]=norm;diagnostics[2]=last;return 0;
}
/* Symmetric T R+R T=J, R=L L^T strictly SPD. J is constructed symmetric.
 * A 3x3 SPD system in the (xx,sqrt(2)xy,yy) basis; no regularization.
 * Outputs T[4], diagnostics[2]=backward residual,operator cond_inf. */
MF_INLINE int lmf_loss(const double *L,const double *K,double *out,double *diagnostics){
    if(K[1]!=K[2])return 1;
    double R[4],J[4];
    for(int i=0;i<2;i++)for(int j=i;j<2;j++){
        R[2*i+j]=R[2*j+i]=L[2*i]*L[2*j]+L[2*i+1]*L[2*j+1];
        const double v=(L[2*i]*K[0]+L[2*i+1]*K[2])*L[2*j]
            +(L[2*i]*K[1]+L[2*i+1]*K[3])*L[2*j+1];
        J[2*i+j]=J[2*j+i]=v;
    }
    const double scale=fmax(R[0],R[3]),rt=sqrt(2.);
    if(!isfinite(scale))return 3;if(scale<=0.||R[0]<=0.)return 2;
    const double a=R[0]/scale,b=R[1]/scale,c=R[3]/scale;
    if(c-(b/sqrt(a))*(b/sqrt(a))<=0.)return 2;
    const double M[9]={2*a,rt*b,0.,rt*b,a+c,rt*b,0.,rt*b,2*c};
    double chol[9]={0.};
    for(int i=0;i<3;i++)for(int j=0;j<=i;j++){
        double t=M[3*i+j];for(int k=0;k<j;k++)t-=chol[3*i+k]*chol[3*j+k];
        if(!isfinite(t))return 3;
        if(i==j){if(t<=0.)return 2;chol[3*i+j]=sqrt(t);}else chol[3*i+j]=t/chol[3*j+j];
    }
    double rhs[3]={J[0]/scale,rt*J[1]/scale,J[3]/scale},z[3];
    for(int i=0;i<3;i++){double t=rhs[i];for(int j=0;j<i;j++)t-=chol[3*i+j]*z[j];z[i]=t/chol[3*i+i];}
    for(int i=2;i>=0;i--){double t=z[i];for(int j=i+1;j<3;j++)t-=chol[3*j+i]*z[j];z[i]=t/chol[3*i+i];}
    double inverse[9];
    for(int col=0;col<3;col++){
        double v[3];for(int i=0;i<3;i++){double t=i==col?1.:0.;for(int j=0;j<i;j++)t-=chol[3*i+j]*v[j];v[i]=t/chol[3*i+i];}
        for(int i=2;i>=0;i--){double t=v[i];for(int j=i+1;j<3;j++)t-=chol[3*j+i]*v[j];v[i]=t/chol[3*i+i];}
        for(int i=0;i<3;i++)inverse[3*i+col]=v[i];
    }
    double mn=0.,in=0.,res=0.,zn=0.,bn=0.;
    for(int i=0;i<3;i++){
        double mm=0.,ii=0.,rr=-rhs[i];for(int j=0;j<3;j++){mm+=fabs(M[3*i+j]);ii+=fabs(inverse[3*i+j]);rr+=M[3*i+j]*z[j];}
        mn=fmax(mn,mm);in=fmax(in,ii);res=fmax(res,fabs(rr));zn=fmax(zn,fabs(z[i]));bn=fmax(bn,fabs(rhs[i]));
    }
    const double denom=mn*zn+bn;double t[4]={z[0],z[1]/rt,z[1]/rt,z[2]};
    for(int k=0;k<4;k++)if(!isfinite(t[k]))return 3;
    for(int k=0;k<4;k++)out[k]=t[k];
    diagnostics[0]=denom>0.?res/denom:res;diagnostics[1]=mn*in;return 0;
}
/* Append one row of U^T to its thin QR using actual orthogonal Givens
 * rotations. R is upper triangular with nonnegative diagonal. No UU^T root. */
MF_INLINE void lmf_qr_row(double x,double y,double *r){
    const double norm=hypot(r[0],x);
    if(norm>0.){
        const double c=r[0]/norm,s=x/norm;
        const double next=c*r[1]+s*y;y=-s*r[1]+c*y;r[0]=norm;r[1]=next;
    }
    r[2]=hypot(r[2],y);
}
/* U column: first two are Efull*Lbase, then sqrt(h)*Ehalf*Q columns. */
MF_INLINE void lmf_u_column(int column,const double *base,const double *ef,
                           const double *eh,const double *q,double root_h,double *u){
    if(column<2){
        for(int i=0;i<2;i++)u[i]=ef[2*i]*base[column]+ef[2*i+1]*base[2+column];
    }else{
        const double x=q[2*(column-2)],y=q[2*(column-2)+1];
        for(int i=0;i<2;i++)u[i]=root_h*(eh[2*i]*x+eh[2*i+1]*y);
    }
}
/* Factor map + Pi rebuild. Lbase and Llinear are DISTINCT stage inputs.
 * q stores [columns][2] in original incoming order. No Psi or ledgers here.
 * frozen/h0 exact-copy bypass (before SPD checks); other failures are
 * transactional. Output arrays must be separate; inputs may alias outputs.
 * Return1 invalid,2 unsupported rank/SPD,3 nonfinite,4 exp scale limit,
 * 5 nonpositive final thermal. diagnostic40 layout documented in harness.
 */
MF_INLINE int lmf_update(const double *base,const double *linear,const double *A,
                        const double *K,const double *q,int columns,double h,
                        const double *provisional,int frozen,double *outL,
                        double *outRaw,double *outDiagnostic){
    if(!base||!linear||!A||!K||!provisional||!outL||!outRaw||!outDiagnostic||columns<0||(!q&&columns)||!isfinite(h)||h<0.)return 1;
    for(int k=0;k<4;k++)if(!isfinite(base[k]))return 1;
    for(int k=0;k<18;k++)if(!isfinite(provisional[k]))return 1;
    double diag[40]={0.},raw[18],newL[4];for(int k=0;k<18;k++)raw[k]=provisional[k];
    if(frozen||h==0.){
        for(int k=0;k<4;k++)newL[k]=base[k];
    }else{
        for(int k=0;k<4;k++)if(!isfinite(linear[k])||!isfinite(A[k])||!isfinite(K[k]))return 1;
        for(int k=0;k<2*columns;k++)if(!isfinite(q[k]))return 1;
        double loss[4],ld[2],ac[4],arg[4],eh[4],ef[4],ed[3];
        int status=lmf_loss(linear,K,loss,ld);if(status)return status;
        for(int k=0;k<4;k++){ac[k]=A[k]-loss[k];arg[k]=.5*h*ac[k];}
        status=lmf_exponential(arg,eh,ed);if(status)return status;
        lmf_multiply(eh,eh,ef);
        double qr[3]={0.},gram[4]={0.};const double rh=sqrt(h);
        for(int k=0;k<columns+2;k++){
            double u[2];lmf_u_column(k,base,ef,eh,q,rh,u);
            for(int i=0;i<2;i++)for(int j=0;j<2;j++)gram[2*i+j]+=u[i]*u[j];
            lmf_qr_row(u[0],u[1],qr);
        }
        if(qr[0]<=0.||qr[2]<=0.)return 2;
        newL[0]=qr[0];newL[1]=0.;newL[2]=qr[1];newL[3]=qr[2];
        double ortho[4]={0.},reconstruction=0.;
        for(int k=0;k<columns+2;k++){
            double u[2];lmf_u_column(k,base,ef,eh,q,rh,u);
            const double v[2]={u[0]/qr[0],(u[1]-(u[0]/qr[0])*qr[1])/qr[2]};
            for(int i=0;i<2;i++)for(int j=0;j<2;j++)ortho[2*i+j]+=v[i]*v[j];
            reconstruction=fmax(reconstruction,fabs(v[0]*qr[0]-u[0]));
            reconstruction=fmax(reconstruction,fabs(v[0]*qr[1]+v[1]*qr[2]-u[1]));
        }
        double R[4];for(int i=0;i<2;i++)for(int j=0;j<2;j++)R[2*i+j]=newL[2*i]*newL[2*j]+newL[2*i+1]*newL[2*j+1];
        double B[4];status=lmc_affine_B(raw,B);if(status)return status;
        const double m=raw[0],S[4]={raw[3]-raw[1]*raw[1]/m,raw[4]-raw[1]*raw[2]/m,raw[5]-raw[2]*raw[1]/m,raw[6]-raw[2]*raw[2]/m};
        diag[28]=raw[17]-.5*(raw[13]+raw[16]);
        for(int i=0;i<2;i++)for(int j=i;j<2;j++){
            const double resolved=(B[2*i]*S[0]+B[2*i+1]*S[2])*B[2*j]
                +(B[2*i]*S[1]+B[2*i+1]*S[3])*B[2*j+1];
            const double pi=raw[7+i]*raw[7+j]/m+resolved+R[2*i+j];
            raw[13+2*i+j]=raw[13+2*j+i]=pi;
        }
        diag[29]=raw[17]-.5*(raw[13]+raw[16]);diag[30]=diag[29]-diag[28];
        if(!isfinite(diag[29]))return 3;if(diag[29]<=0.)return 5;
        for(int k=0;k<4;k++){
            diag[k]=loss[k];diag[4+k]=ac[k];diag[8+k]=eh[k];diag[12+k]=ef[k];diag[16+k]=R[k];diag[20+k]=B[k];diag[24+k]=raw[13+k]-provisional[13+k];
        }
        diag[31]=ld[0];diag[32]=ld[1];for(int k=0;k<3;k++)diag[33+k]=ed[k];
        ortho[0]-=1.;ortho[3]-=1.;diag[36]=lmf_norm(ortho);diag[37]=reconstruction;
        double delta[4];for(int k=0;k<4;k++)delta[k]=R[k]-gram[k];diag[38]=lmf_norm(delta);
        diag[39]=lmf_norm(gram);
    }
    for(int k=0;k<4;k++)if(!isfinite(newL[k]))return 3;
    for(int k=0;k<18;k++)if(!isfinite(raw[k]))return 3;
    for(int k=0;k<40;k++)if(!isfinite(diag[k]))return 3;
    for(int k=0;k<4;k++)outL[k]=newL[k];for(int k=0;k<18;k++)outRaw[k]=raw[k];
    for(int k=0;k<40;k++)outDiagnostic[k]=diag[k];return 0;
}
#undef MF_INLINE
#endif
