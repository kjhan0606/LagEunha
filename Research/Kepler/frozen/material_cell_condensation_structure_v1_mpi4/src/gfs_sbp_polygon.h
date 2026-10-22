#ifndef GFS_SBP_POLYGON_H
#define GFS_SBP_POLYGON_H
#include <math.h>
#include <float.h>
/* Host construction of a P2 positive diagonal-norm generalized SBP matrix.
 * The caller supplies positive cubature M, boundary diagonal E, and P2 values
 * V and exact derivatives Vx,Vy at the SAME nodes. Boundary nodes are a subset
 * of these nodes (selection trace). Volume degree 3 / face degree 4 exactness
 * is required; compatibility is measured, never assumed or repaired using a
 * fluid state. Q matrices are reusable by native CPU or CUDA volume kernels.
 * No Euler state projection, evolution, viscosity, or residual forcing here.
 */
typedef struct {
    double compatibility,skew_error,rowsum_error,polynomial_error;
    double minimum_mass,gram_minimum_pivot;
} GfsSbpMatrixAudit;

static inline int gfs_sbp_polygon_p2(int n,const double *mass,
    const double *ex,const double *ey,const double *V,const double *Vx,
    const double *Vy,double *Qx,double *Qy,double *scratch,GfsSbpMatrixAudit *out){
    if(n<6 || !mass || !ex || !ey || !V || !Vx || !Vy || !Qx || !Qy || !scratch || !out)return 0;
    out->compatibility=out->skew_error=out->rowsum_error=out->polynomial_error=0;
    out->minimum_mass=DBL_MAX;out->gram_minimum_pivot=DBL_MAX;
    double L[6][6]={{0}},inverse[6][6]={{0}};
    for(int i=0;i<n;i++){
        if(!(mass[i]>0) || !isfinite(mass[i]))return 0;
        out->minimum_mass=fmin(out->minimum_mass,mass[i]);
        for(int a=0;a<6;a++)for(int b=0;b<=a;b++)L[a][b]+=V[6*i+a]*V[6*i+b];
    }
    for(int a=0;a<6;a++)for(int b=0;b<=a;b++){
        double z=L[a][b];for(int c=0;c<b;c++)z-=L[a][c]*L[b][c];
        if(a==b){if(!(z>0) || !isfinite(z))return 0;
            out->gram_minimum_pivot=fmin(out->gram_minimum_pivot,z);L[a][b]=sqrt(z);
        }else L[a][b]=z/L[b][b];
    }
    for(int col=0;col<6;col++){
        double z[6];for(int a=0;a<6;a++){
            z[a]=a==col?1:0;for(int b=0;b<a;b++)z[a]-=L[a][b]*z[b];z[a]/=L[a][a];
        }
        for(int a=5;a>=0;a--){for(int b=a+1;b<6;b++)z[a]-=L[b][a]*z[b];z[a]/=L[a][a];inverse[a][col]=z[a];}
    }
    /* P=(V^T V)^-1 V^T is an operator-construction right inverse only.
     * Euler states keep all n degrees of freedom; none are compressed to P2.
     * scratch has 12*n doubles: row-major P^T then A. */
    double *P=scratch,*A=scratch+6*n;
    for(int i=0;i<n;i++)for(int a=0;a<6;a++){
        double z=0;for(int b=0;b<6;b++)z+=inverse[a][b]*V[6*i+b];P[6*i+a]=z;
    }
    for(int axis=0;axis<2;axis++){
        const double *E=axis?ey:ex,*D=axis?Vy:Vx;double *Q=axis?Qy:Qx;
        for(int i=0;i<n;i++)for(int a=0;a<6;a++)A[6*i+a]=mass[i]*D[6*i+a]-.5*E[i]*V[6*i+a];
        double C[6][6]={{0}};
        for(int a=0;a<6;a++)for(int b=0;b<6;b++)for(int i=0;i<n;i++)C[a][b]+=A[6*i+a]*V[6*i+b];
        double scale=0;for(int i=0;i<n;i++)scale+=fabs(E[i]);if(!(scale>0))return 0;
        for(int a=0;a<6;a++)for(int b=0;b<6;b++)
            out->compatibility=fmax(out->compatibility,fabs(C[a][b]+C[b][a])/scale);
        /* Skew symmetrization removes floating point cubature roundoff only.
         * The unmodified compatibility residual above is a required gate. */
        double K[6][6];for(int a=0;a<6;a++)for(int b=0;b<6;b++)K[a][b]=.5*(C[a][b]-C[b][a]);
        for(int i=0;i<n;i++){
            Q[i*n+i]=.5*E[i];
            for(int j=i+1;j<n;j++){
                double s=0;for(int a=0;a<6;a++){
                    s+=A[6*i+a]*P[6*j+a]-A[6*j+a]*P[6*i+a];
                    for(int b=0;b<6;b++)s+=P[6*i+a]*K[a][b]*P[6*j+b];
                }Q[i*n+j]=s;Q[j*n+i]=-s;
            }
        }
        for(int i=0;i<n;i++){
            double row=0;for(int j=0;j<n;j++){
                row+=Q[i*n+j];out->skew_error=fmax(out->skew_error,fabs(Q[i*n+j]+Q[j*n+i]-(i==j?E[i]:0))/scale);
            }out->rowsum_error=fmax(out->rowsum_error,fabs(row)/scale);
            for(int a=0;a<6;a++){
                double z=-mass[i]*D[6*i+a];for(int j=0;j<n;j++)z+=Q[i*n+j]*V[6*j+a];
                out->polynomial_error=fmax(out->polynomial_error,fabs(z)/scale);
            }
        }
    }return 1;
}
#endif
