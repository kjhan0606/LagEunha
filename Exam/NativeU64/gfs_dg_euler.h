#ifndef GFS_DG_EULER_H
#define GFS_DG_EULER_H
#include "gfs_moment_initial.h"
/* Euler spatial representation, not an independent kinetic reservoir.
 * All four conserved components use the SAME mean-free polynomial basis.
 * Projection/primitive/weak-form helpers are native C and CUDA compatible.
 * No coefficient, density, pressure, or total-energy clipping occurs here.
 */
#define GFS_DG_MAX_ORDER 6
#define GFS_DG_MAX_BASIS 28
#define GFS_DG_MAX_GAUSS 12
typedef struct {
    GfsMomentIC ic;
    int order,ng[3];
    double node[3][GFS_DG_MAX_GAUSS],weight[3][GFS_DG_MAX_GAUSS];
    double density_slope[2]; /* exact affine controls only, not Kepler */
} GfsDgConfig;
typedef struct {
    double area,gx,gy,l00,l10,l11;
    double mean[GFS_DG_MAX_BASIS];
} GfsDgBasis;
typedef struct { double c[4][GFS_DG_MAX_BASIS]; } GfsDgPolynomial;
enum {
    DG_MIN_RHO,DG_MIN_P,DG_MAX_P,DG_MAX_PERR,DG_L1_PERR,DG_MEAN_P,
    DG_MAX_P_QUAD_DELTA,DG_WORST_X,DG_WORST_Y,
    DG_QUAD_MASS_REL,DG_QUAD_PX_REL,DG_QUAD_PY_REL,DG_QUAD_E_OVER_THERMAL,
    DG_CONSERV_RHO,DG_CONSERV_PX,DG_CONSERV_PY,DG_CONSERV_E,
    DG_NMETRICS
};
typedef struct {
    GfsDgBasis basis;
    GfsDgPolynomial polynomial;
    double integral[4],coarse_integral[4],metrics[DG_NMETRICS];
    int status,nonpositive_rho,nonpositive_p,worst_kind;
} GfsDgProjectionAudit;

GFS_MOM_HD static inline int gfs_dg_nb(int p){return (p+1)*(p+2)/2;}

/* Positive tensor Gauss quadrature on a centroid-fan triangle, Duffy map.
 * q indexes edge*ng^2 + iu*ng+iv. ng>=p+1 integrates a Pp Gram matrix
 * exactly in exact arithmetic. The three independent rules are for
 * projection, its integration control, and spatial-error assessment.
 */
GFS_MOM_HD static inline void gfs_dg_quad(GfsMomentPolygon poly,
    const GfsDgBasis *b,const GfsDgConfig *cfg,int rule,int q,
    double *x,double *y,double *w){
    int n=cfg->ng[rule],edge=q/(n*n),uv=q%(n*n),iu=uv/n,iv=uv%n;
    int j=(edge+1)%poly.n;
    double u=cfg->node[rule][iu],v=cfg->node[rule][iv];
    double ax=b->gx,ay=b->gy,bx=poly.x[edge],by=poly.y[edge],cx=poly.x[j],cy=poly.y[j];
    /* The determinant of two identical vectors must be exactly zero. A
     * contracted multiply-subtract can otherwise return one product's
     * rounding residual, giving an inactive edge a signed fake measure. */
    if(bx==cx && by==cy){*x=bx;*y=by;*w=0;return;}
    double area2=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax);
    *x=(1-u)*(1-v)*ax+u*bx+(1-u)*v*cx;
    *y=(1-u)*(1-v)*ay+u*by+(1-u)*v*cy;
    *w=area2*(1-u)*cfg->weight[rule][iu]*cfg->weight[rule][iv];
}

GFS_MOM_HD static inline void gfs_dg_monomial(const GfsDgBasis *b,int p,
    double x,double y,double *phi,double *dx,double *dy){
    double u=(x-b->gx)/b->l00,v=(y-b->gy-b->l10*u)/b->l11;
    double up[GFS_DG_MAX_ORDER+1],vp[GFS_DG_MAX_ORDER+1];up[0]=vp[0]=1;
    for(int k=1;k<=p;k++){up[k]=up[k-1]*u;vp[k]=vp[k-1]*v;}
    int k=0;for(int degree=0;degree<=p;degree++)for(int j=0;j<=degree;j++,k++){
        int i=degree-j;phi[k]=up[i]*vp[j];
        if(dx || dy){
            double pu=i?i*up[i-1]*vp[j]:0,pv=j?j*up[i]*vp[j-1]:0;
            if(dx)dx[k]=pu/b->l00-pv*b->l10/(b->l00*b->l11);
            if(dy)dy[k]=pv/b->l11;
        }
    }
}

GFS_MOM_HD static inline void gfs_dg_evaluate(const GfsDgBasis *b,
    const GfsDgPolynomial *poly,int order,double x,double y,double u[4]){
    double phi[GFS_DG_MAX_BASIS];gfs_dg_monomial(b,order,x,y,phi,0,0);
    int n=gfs_dg_nb(order);
    for(int c=0;c<4;c++){
        u[c]=poly->c[c][0];for(int k=1;k<n;k++)u[c]=fma(poly->c[c][k],phi[k]-b->mean[k],u[c]);
    }
}

GFS_MOM_HD static inline double gfs_dg_pressure(const double u[4],double gamma){
    if(!(u[0]>0))return -DBL_MAX;
    double eth=fma(-.5*u[1],u[1]/u[0],u[3]);eth=fma(-.5*u[2],u[2]/u[0],eth);
    return (gamma-1)*eth;
}

GFS_MOM_HD static inline void gfs_dg_ic_state(const GfsDgConfig *cfg,
    double x,double y,double u[4]){
    double rho,vx,vy;gfs_moment_ic_value(&cfg->ic,x,y,&rho,&vx,&vy);
    if(cfg->ic.kind==1)rho*=1+cfg->density_slope[0]*(x-cfg->ic.cx)+cfg->density_slope[1]*(y-cfg->ic.cy);
    u[0]=rho;u[1]=rho*vx;u[2]=rho*vy;
    u[3]=cfg->ic.pressure/(cfg->ic.gamma-1)+.5*rho*(vx*vx+vy*vy);
}

/* Single factorization for all Euler components. Mean-free high modes
 * decouple the constant: c0=integral/area is preserved structurally.
 * No regularization or projection onto an admissible cone is hidden here. */
GFS_MOM_HD static inline int gfs_dg_factor(double *H,int n){
    for(int i=1;i<n;i++)for(int j=1;j<=i;j++){
        double x=H[i*GFS_DG_MAX_BASIS+j];
        for(int k=1;k<j;k++)x=fma(-H[i*GFS_DG_MAX_BASIS+k],H[j*GFS_DG_MAX_BASIS+k],x);
        if(i==j){if(!(x>0) || !isfinite(x))return 0;H[i*GFS_DG_MAX_BASIS+j]=sqrt(x);}
        else H[i*GFS_DG_MAX_BASIS+j]=x/H[j*GFS_DG_MAX_BASIS+j];
    }
    return 1;
}

GFS_MOM_HD static inline void gfs_dg_solve(const double *L,int n,double *rhs){
    for(int i=1;i<n;i++){
        for(int j=1;j<i;j++)rhs[i]=fma(-L[i*GFS_DG_MAX_BASIS+j],rhs[j],rhs[i]);
        rhs[i]/=L[i*GFS_DG_MAX_BASIS+i];
    }
    for(int i=n-1;i>=1;i--){
        for(int j=i+1;j<n;j++)rhs[i]=fma(-L[j*GFS_DG_MAX_BASIS+i],rhs[j],rhs[i]);
        rhs[i]/=L[i*GFS_DG_MAX_BASIS+i];
    }
}

GFS_MOM_HD static inline int gfs_dg_project(GfsMomentPolygon polygon,
    const GfsDgConfig *cfg,double generator_x,double generator_y,
    GfsDgProjectionAudit *out,GfsDgPolynomial *coarse){
    out->status=-1;
    if(cfg->order<1 || cfg->order>GFS_DG_MAX_ORDER)return 0;
    GfsMomentDensity geom={0};if(!gfs_moment_geometry(polygon,&geom))return 0;
    GfsDgBasis *b=&out->basis;
    b->area=geom.area;b->gx=geom.gx;b->gy=geom.gy;
    b->l00=geom.l00;b->l10=geom.l10;b->l11=geom.l11;
    int n=gfs_dg_nb(cfg->order);
    double gram[GFS_DG_MAX_BASIS*GFS_DG_MAX_BASIS]={0};
    double rhs[2][4][GFS_DG_MAX_BASIS]={0};
    for(int k=0;k<GFS_DG_MAX_BASIS;k++)b->mean[k]=0;
    for(int rule=0;rule<2;rule++){
        int nq=polygon.n*cfg->ng[rule]*cfg->ng[rule];
        for(int q=0;q<nq;q++){
            double x,y,w,phi[GFS_DG_MAX_BASIS],u[4];gfs_dg_quad(polygon,b,cfg,rule,q,&x,&y,&w);
            if(w==0)continue;if(!(w>0)){out->status=-3;return 0;}
            gfs_dg_monomial(b,cfg->order,x,y,phi,0,0);
            gfs_dg_ic_state(cfg,generator_x+x,generator_y+y,u);
            if(!(u[0]>0) || !(gfs_dg_pressure(u,cfg->ic.gamma)>0)){out->status=-4;return 0;}
            for(int i=0;i<n;i++){
                if(rule==0){
                    b->mean[i]+=w*phi[i];
                    for(int j=1;j<=i;j++)gram[i*GFS_DG_MAX_BASIS+j]+=w*phi[i]*phi[j];
                }
                for(int c=0;c<4;c++)rhs[rule][c][i]+=w*u[c]*phi[i];
            }
        }
    }
    for(int k=0;k<n;k++)b->mean[k]/=b->area;b->mean[0]=1;
    for(int i=1;i<n;i++)for(int j=1;j<=i;j++)
        gram[i*GFS_DG_MAX_BASIS+j]-=b->area*b->mean[i]*b->mean[j];
    out->status=-2;if(!gfs_dg_factor(gram,n))return 0;
    for(int rule=0;rule<2;rule++)for(int c=0;c<4;c++){
        double mean=rhs[rule][c][0];
        for(int k=1;k<n;k++)rhs[rule][c][k]-=b->mean[k]*mean;
        gfs_dg_solve(gram,n,rhs[rule][c]);
        GfsDgPolynomial *p=rule?&out->polynomial:coarse;
        p->c[c][0]=mean/b->area;
        for(int k=1;k<n;k++)p->c[c][k]=rhs[rule][c][k];
        if(rule)out->integral[c]=mean;else out->coarse_integral[c]=mean;
    }
    out->status=1;return 1;
}

/* Standard Euler point flux. Kinetic energy (including resolved rotation
 * and strain) stays inside the single total-energy state and Riemann flux. */
GFS_MOM_HD static inline int gfs_dg_euler_flux(const double u[4],double gamma,
    double fx[4],double fy[4]){
    double p=gfs_dg_pressure(u,gamma);if(!(u[0]>0) || !(p>0))return 0;
    double vx=u[1]/u[0],vy=u[2]/u[0];
    fx[0]=u[1];fx[1]=u[1]*vx+p;fx[2]=u[2]*vx;fx[3]=(u[3]+p)*vx;
    fy[0]=u[2];fy[1]=u[1]*vy;fy[2]=u[2]*vy+p;fy[3]=(u[3]+p)*vy;
    return 1;
}

/* Fixed-in-time test functions for a stage: integrate phi*(Euler source)
 * plus grad(phi).Euler flux. Face term uses one actual ALE quadrature
 * flux shared by the two cells. A future moving basis needs its explicit
 * time derivative; these helpers do not silently omit it for evolution. */
GFS_MOM_HD static inline int gfs_dg_volume_point(const GfsDgBasis *b,int order,
    const double u[4],double gamma,double x,double y,double gx,double gy,
    double weight,double rhs[4][GFS_DG_MAX_BASIS]){
    double phi[GFS_DG_MAX_BASIS],dx[GFS_DG_MAX_BASIS],dy[GFS_DG_MAX_BASIS],fx[4],fy[4];
    if(!gfs_dg_euler_flux(u,gamma,fx,fy))return 0;
    gfs_dg_monomial(b,order,x,y,phi,dx,dy);
    double source[4]={0,u[0]*gx,u[0]*gy,u[1]*gx+u[2]*gy};
    int n=gfs_dg_nb(order);
    for(int k=0;k<n;k++)for(int c=0;c<4;c++){
        double test=k?phi[k]-b->mean[k]:1;
        rhs[c][k]+=weight*(fx[c]*dx[k]+fy[c]*dy[k]+source[c]*test);
    }
    return 1;
}

GFS_MOM_HD static inline void gfs_dg_face_point(const GfsDgBasis *b,int order,
    double x,double y,const double ale_flux[4],double signed_weight,
    double rhs[4][GFS_DG_MAX_BASIS]){
    double phi[GFS_DG_MAX_BASIS];gfs_dg_monomial(b,order,x,y,phi,0,0);
    int n=gfs_dg_nb(order);
    for(int k=0;k<n;k++)for(int c=0;c<4;c++)
        rhs[c][k]-=signed_weight*ale_flux[c]*(k?phi[k]-b->mean[k]:1);
}
#endif
