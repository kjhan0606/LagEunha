#ifndef LAGFORCE_REVERSIBLE_MEMORY_H
#define LAGFORCE_REVERSIBLE_MEMORY_H
#include <math.h>
#ifdef __CUDACC__
#define RM_INLINE static __host__ __device__ __forceinline__
#else
#define RM_INLINE static inline
#endif

/* Experimental owner-centered residual, not a complete hydro solver.
 * dx/dy: periodic minimum-image offsets of ALL candidates inside a certified
 * search radius. Caller must guarantee no missing point inside coverage.
 * owner_slot identifies the unique zero offset. work needs 7*n doubles.
 * Outputs are usable only on return 0. meta = radius, max moment residual.
 * No state reset, thermal conversion, force scatter or timestep occurs here.
 */
RM_INLINE int rm_quadratic_rows_impl(int n, const double *dx, const double *dy,
    int owner_slot, double coverage, double target,
    double *row, double *work, double *meta, double *gradient_x, double *gradient_y)
{
    if((gradient_x==0)!=(gradient_y==0))return 1;
    if(n<7 || owner_slot<0 || owner_slot>=n || !(coverage>0)
       || !(target>6 && target<n) || !isfinite(coverage)) return 1;
    if(dx[owner_slot]!=0 || dy[owner_slot]!=0) return 1;
    for(int i=0;i<n;i++) if(!isfinite(dx[i]) || !isfinite(dy[i])) return 1;
    double high=coverage,low=0,sum=0;
    for(int i=0;i<n;i++) {
        double t=fmax(1-(dx[i]*dx[i]+dy[i]*dy[i])/(high*high),0.0);
        sum+=t*t*t*t;
    }
    if(!(sum>target)) return 2;
    for(int iter=0;iter<55;iter++) {
        double mid=(low+high)*.5;sum=0;
        for(int i=0;i<n;i++) {
            double t=fmax(1-(dx[i]*dx[i]+dy[i]*dy[i])/(mid*mid),0.0);
            sum+=t*t*t*t;
        }
        if(sum<target)low=mid;else high=mid;
    }
    double radius=(low+high)*.5;
    double *Q=work,*sqrtw=work+6*n,R[36]={0},y[6]={0};
    for(int i=0;i<n;i++) {
        double xx=dx[i]/radius,yy=dy[i]/radius;
        double t=fmax(1-xx*xx-yy*yy,0.0),p[6]={1,xx,yy,xx*xx,xx*yy,yy*yy};
        sqrtw[i]=t*t;
        for(int j=0;j<6;j++)Q[j*n+i]=sqrtw[i]*p[j];
    }
    /* Twice-reorthogonalized modified Gram-Schmidt, six columns. */
    for(int j=0;j<6;j++) {
        for(int pass=0;pass<2;pass++)for(int k=0;k<j;k++) {
            double dot=0;
            for(int i=0;i<n;i++)dot+=Q[k*n+i]*Q[j*n+i];
            R[k*6+j]+=dot;
            for(int i=0;i<n;i++)Q[j*n+i]-=dot*Q[k*n+i];
        }
        double norm=0;
        for(int i=0;i<n;i++)norm+=Q[j*n+i]*Q[j*n+i];
        norm=sqrt(norm);
        if(!(norm>1e-12*sqrt(target)) || !isfinite(norm))return 3;
        R[j*6+j]=norm;
        for(int i=0;i<n;i++)Q[j*n+i]/=norm;
    }
    /* R^T y = p(0), prediction=sqrt(W)*Q*y. */
    for(int j=0;j<6;j++) {
        double b=(j==0)?1:0;
        for(int k=0;k<j;k++)b-=R[k*6+j]*y[k];
        y[j]=b/R[j*6+j];
    }
    double moments[6]={0},maximum=0;
    for(int i=0;i<n;i++) {
        double a=0;
        for(int j=0;j<6;j++)a+=Q[j*n+i]*y[j];
        row[i]=(i==owner_slot?1.:0.)-sqrtw[i]*a;
        if(!isfinite(row[i]))return 4;
        double xx=dx[i]/radius,yy=dy[i]/radius;
        double p[6]={1,xx,yy,xx*xx,xx*yy,yy*yy};
        for(int j=0;j<6;j++)moments[j]+=p[j]*row[i];
    }
    for(int j=0;j<6;j++)maximum=fmax(maximum,fabs(moments[j]));
    meta[0]=radius;meta[1]=maximum;
    if(!(maximum<1e-10))return 4;
    if(gradient_x){
        /* Derivative prediction rows in physical coordinates. Reuse the
         * identical factorization: R^T y=e_axis/radius, D=sqrt(W)Qy.
         * The residual row above is not recomputed or modified. */
        for(int axis=0;axis<2;axis++){
            double gy[6]={0},gm[6]={0};
            double *g=axis==0?gradient_x:gradient_y;
            for(int j=0;j<6;j++){
                double b=j==axis+1?1./radius:0.;
                for(int k=0;k<j;k++)b-=R[k*6+j]*gy[k];
                gy[j]=b/R[j*6+j];
            }
            for(int i=0;i<n;i++){
                double a=0;
                for(int j=0;j<6;j++)a+=Q[j*n+i]*gy[j];
                g[i]=sqrtw[i]*a;
                if(!isfinite(g[i]))return 4;
                double xx=dx[i]/radius,yy=dy[i]/radius;
                double p[6]={1,xx,yy,xx*xx,xx*yy,yy*yy};
                for(int j=0;j<6;j++)gm[j]+=p[j]*(radius*g[i]);
            }
            for(int j=0;j<6;j++)
                if(!(fabs(gm[j]-(j==axis+1?1.:0.))<1e-10))return 4;
        }
    }
    return 0;
}
/* Original API: derivative work is omitted. */
RM_INLINE int rm_quadratic_row(int n,const double *dx,const double *dy,
    int owner_slot,double coverage,double target,double *row,double *work,double *meta)
{
    return rm_quadratic_rows_impl(n,dx,dy,owner_slot,coverage,target,row,work,meta,0,0);
}
/* Optional gradient rows need n doubles each, separate from row and work.
 * sum_j gradient_x[j]*v[j] gives dv/dx at owner (and similarly y).
 * Spin in 2D is .5*(gradient_x dot vy - gradient_y dot vx).
 * This helper does not change the production memory transport equation. */
RM_INLINE int rm_quadratic_row_gradient(int n,const double *dx,const double *dy,
    int owner_slot,double coverage,double target,double *row,double *work,double *meta,
    double *gradient_x,double *gradient_y)
{
    if(!gradient_x || !gradient_y)return 1;
    return rm_quadratic_rows_impl(n,dx,dy,owner_slot,coverage,target,row,work,meta,gradient_x,gradient_y);
}
#undef RM_INLINE
#endif
