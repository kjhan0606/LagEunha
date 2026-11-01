#ifndef GFS_ENRICHED_EULER_H
#define GFS_ENRICHED_EULER_H
#include "gfs_dg_euler.h"
/* Common scalar Euler space P2 + the four DISCLOSED reference fields.
 * Coefficients are ordinary unknowns, not a stationary-state constraint.
 * This header supplies representation and full fixed-basis Euler weak form;
 * it does NOT supply a moving-mesh time integrator or weight/mass solve.
 * There is no kinetic reservoir, residual subtraction or energy correction.
 */
#define GFS_ENRICH_N 10
#define GFS_ENRICH_NA 64
#define GFS_ENRICH_NR 32
#define GFS_ENRICH_SEGMENTS 32
typedef struct {
    GfsDgConfig field;
    int na,nr;
    double an[GFS_ENRICH_NA],aw[GFS_ENRICH_NA];
    double rn[GFS_ENRICH_NR],rw[GFS_ENRICH_NR];
    double boost[2];
} GfsEnrichConfig;
typedef struct {
    GfsDgBasis geometry;
    double generator[2],scale[4];
    /* f_k = sum_j L_kj q_j. A zero diagonal explicitly excludes q_k.
     * Recursive residual evaluation avoids subtracting nearly equal Gram
     * matrices of a smooth reference field and its polynomial projection. */
    double L[GFS_ENRICH_N][GFS_ENRICH_N],mean[GFS_ENRICH_N];
    double minimum_norm,omitted_norm,orthogonality_error,area_error;
    int rank,status;
} GfsEnrichBasis;
typedef struct { double c[4][GFS_ENRICH_N]; } GfsEnrichState;
typedef struct {
    double anchor[2],cut[GFS_ENRICH_SEGMENTS+1];
    int ns,anchor_edge;
} GfsEnrichQuadrature;

/* Algebraically identical to .5*(1+tanh(z)), evaluated without subtracting
 * a number nearly equal to one. This changes no IC function or parameter.
 * Derivative is with respect to r for z=(r-r0)/.25. */
GFS_MOM_HD static inline double gfs_enrich_taper(double z,double *dr){
    double e=exp(-2*fabs(z)),a=1/(1+e);
    if(dr)*dr=8*e*a*a;
    return z>=0?a:e*a;
}

GFS_MOM_HD static inline double gfs_enrich_density_contrast(const GfsEnrichConfig *cfg,
    double x,double y){
    const GfsMomentIC *ic=&cfg->field.ic;
    double dx=x-ic->cx,dy=y-ic->cy;
    if(ic->kind==1)return ic->rho*(cfg->field.density_slope[0]*dx+cfg->field.density_slope[1]*dy);
    double r=hypot(dx,dy);
    return ic->rin/(r+1.e-6)*gfs_enrich_taper((r-ic->rin)/.25,0)*gfs_enrich_taper((ic->rout-r)/.25,0);
}

/* Conservative reference value and its actual analytic derivatives.
 * The radial density derivative is undefined exactly at r=0; requesting it
 * there fails. Values at that point remain valid. No center smoothing. */
GFS_MOM_HD static inline int gfs_enrich_reference(const GfsEnrichConfig *cfg,
    double x,double y,double U[4],double *Ux,double *Uy){
    const GfsMomentIC *ic=&cfg->field.ic;
    double dx=x-ic->cx,dy=y-ic->cy,rho,vx,vy,rx=0,ry=0;
    double vxx,vxy,vyx,vyy;
    if(ic->kind==1){
        double a=cfg->field.density_slope[0],b=cfg->field.density_slope[1];
        rho=ic->rho*(1+a*dx+b*dy);rx=ic->rho*a;ry=ic->rho*b;
        vx=ic->vx+ic->B[0]*dx+ic->B[1]*dy;
        vy=ic->vy+ic->B[2]*dx+ic->B[3]*dy;
        vxx=ic->B[0];vxy=ic->B[1];vyx=ic->B[2];vyy=ic->B[3];
    }else{
        double r=hypot(dx,dy),di,do_;
        double inner=gfs_enrich_taper((r-ic->rin)/.25,&di),outer=gfs_enrich_taper((ic->rout-r)/.25,&do_),den=r+1.e-6;
        rho=.001+ic->rin/den*inner*outer;
        if(Ux || Uy){
            if(r==0)return 0;
            double dr=ic->rin*((di*outer-inner*do_)/den-inner*outer/(den*den));
            rx=dr*dx/r;ry=dr*dy/r;
        }
        double z=r*r+ic->eps*ic->eps,omega=pow(z,-.75);
        double ox=-1.5*dx*omega/z,oy=-1.5*dy*omega/z;
        vx=-dy*omega;vy=dx*omega;
        vxx=-dy*ox;vxy=-omega-dy*oy;vyx=omega+dx*ox;vyy=dx*oy;
    }
    vx+=cfg->boost[0];vy+=cfg->boost[1];
    U[0]=rho;U[1]=rho*vx;U[2]=rho*vy;
    double v2=vx*vx+vy*vy;U[3]=ic->pressure/(ic->gamma-1)+.5*rho*v2;
    if(Ux){Ux[0]=rx;Ux[1]=rx*vx+rho*vxx;Ux[2]=rx*vy+rho*vyx;
        Ux[3]=.5*rx*v2+rho*(vx*vxx+vy*vyx);}
    if(Uy){Uy[0]=ry;Uy[1]=ry*vx+rho*vxy;Uy[2]=ry*vy+rho*vyy;
        Uy[3]=.5*ry*v2+rho*(vx*vxy+vy*vyy);}
    return rho>0 && isfinite(U[3]);
}

/* Fan about the closest point of the convex polygon to the physical
 * reference origin (the origin itself if inside). This resolves the actual
 * r+1e-6 field rather than moving generators or changing the density. */
GFS_MOM_HD static inline int gfs_enrich_quadrature_init(GfsMomentPolygon p,
    const GfsEnrichConfig *cfg,const GfsEnrichBasis *b,GfsEnrichQuadrature *q){
    double ox=cfg->field.ic.cx-b->generator[0],oy=cfg->field.ic.cy-b->generator[1];
    int inside=1,best_edge=-1;double best=DBL_MAX,ax=0,ay=0;
    for(int i=0;i<p.n;i++){
        int j=(i+1)%p.n;double ex=p.x[j]-p.x[i],ey=p.y[j]-p.y[i],e2=ex*ex+ey*ey;
        if(e2==0)continue;
        if(ex*(oy-p.y[i])-ey*(ox-p.x[i])<0)inside=0;
        double t=((ox-p.x[i])*ex+(oy-p.y[i])*ey)/e2;t=fmax(0.,fmin(1.,t));
        /* A projected endpoint is the stored endpoint itself. Rebuilding
         * it as a+(b-a) can shift it by one ulp and place the fan anchor
         * just outside the adjacent support plane. No geometric tolerance
         * or finite-length edge cutoff is used. */
        double x=t==0?p.x[i]:(t==1?p.x[j]:p.x[i]+t*ex);
        double y=t==0?p.y[i]:(t==1?p.y[j]:p.y[i]+t*ey),d=hypot(x-ox,y-oy);
        if(d<best){best=d;ax=x;ay=y;best_edge=i;}
    }
    if(inside){ax=ox;ay=oy;best=0;}
    if(!(best<DBL_MAX))return 0;
    q->anchor[0]=ax;q->anchor[1]=ay;
    q->anchor_edge=inside?-1:best_edge;
    double h=0;for(int i=0;i<p.n;i++)h=fmax(h,hypot(p.x[i]-ax,p.y[i]-ay));
    if(!(h>0))return 0;
    double first=fmin(1.,(best+1.e-6)/h);q->cut[0]=0;q->cut[1]=first;q->ns=1;
    while(q->cut[q->ns]<1 && q->ns<GFS_ENRICH_SEGMENTS){
        q->cut[q->ns+1]=fmin(1.,4*q->cut[q->ns]);q->ns++;
    }
    return q->cut[q->ns]==1;
}

GFS_MOM_HD static inline int gfs_enrich_nquad(GfsMomentPolygon p,
    const GfsEnrichConfig *cfg,const GfsEnrichQuadrature *rule){
    return p.n*rule->ns*cfg->na*cfg->nr;
}

GFS_MOM_HD static inline int gfs_enrich_quad(GfsMomentPolygon p,
    const GfsEnrichConfig *cfg,const GfsEnrichQuadrature *r,int q,
    double *x,double *y,double *w){
    int a=q%cfg->na;q/=cfg->na;int radial=q%cfg->nr;q/=cfg->nr;
    int s=q%r->ns,edge=q/r->ns,j=(edge+1)%p.n;
    double ax=r->anchor[0],ay=r->anchor[1],bx=p.x[edge]-ax,by=p.y[edge]-ay;
    double cx=p.x[j]-ax,cy=p.y[j]-ay;
    if(p.x[edge]==p.x[j] && p.y[edge]==p.y[j]){*x=p.x[edge];*y=p.y[edge];*w=0;return 1;}
    double jac=bx*cy-by*cx;
    /* Anchor can lie on an actual polygon edge. Exact collinearity has
     * zero measure. Negative resolved measure is an explicit failure. */
    if(edge==r->anchor_edge || (bx==0 && by==0) || (cx==0 && cy==0))jac=0;
    if(jac<0)return 0;
    double t=r->cut[s]+(r->cut[s+1]-r->cut[s])*cfg->rn[radial],a0=cfg->an[a];
    *x=ax+t*((1-a0)*bx+a0*cx);*y=ay+t*((1-a0)*by+a0*cy);
    *w=jac*t*(r->cut[s+1]-r->cut[s])*cfg->rw[radial]*cfg->aw[a];
    return isfinite(*w);
}

GFS_MOM_HD static inline int gfs_enrich_raw(const GfsEnrichBasis *b,
    const GfsEnrichConfig *cfg,double x,double y,double f[GFS_ENRICH_N],
    double *fx,double *fy){
    double U[4],Ux[4],Uy[4];
    gfs_dg_monomial(&b->geometry,2,x,y,f,fx,fy);
    if(!gfs_enrich_reference(cfg,b->generator[0]+x,b->generator[1]+y,U,
        fx?Ux:0,fy?Uy:0))return 0;
    for(int c=0;c<4;c++){
        f[6+c]=U[c]/b->scale[c];
        if(fx)fx[6+c]=Ux[c]/b->scale[c];if(fy)fy[6+c]=Uy[c]/b->scale[c];
    }
    /* P2 already contains constants. Removing the known constant density
     * floor BEFORE adding it, and the constant thermal E BEFORE adding it,
     * leaves the SAME function space while avoiding a tiny halo density
     * feature being computed by subtracting two nearly equal full densities.
     * Derivatives are unchanged; initialized means still use full rho/E. */
    f[6]=gfs_enrich_density_contrast(cfg,b->generator[0]+x,b->generator[1]+y)/b->scale[0];
    f[9]=.5*(U[1]*U[1]+U[2]*U[2])/(U[0]*b->scale[3]);
    return 1;
}

GFS_MOM_HD static inline void gfs_enrich_transform(const GfsEnrichBasis *b,
    int n,const double *f,const double *fx,const double *fy,double *q,double *qx,double *qy){
    for(int i=0;i<n;i++){
        q[i]=f[i];if(qx)qx[i]=fx[i];if(qy)qy[i]=fy[i];
        for(int j=0;j<i;j++){
            q[i]=fma(-b->L[i][j],q[j],q[i]);
            if(qx)qx[i]=fma(-b->L[i][j],qx[j],qx[i]);
            if(qy)qy[i]=fma(-b->L[i][j],qy[j],qy[i]);
        }
        if(b->L[i][i]>0){q[i]/=b->L[i][i];if(qx)qx[i]/=b->L[i][i];if(qy)qy[i]/=b->L[i][i];}
        else{q[i]=0;if(qx)qx[i]=0;if(qy)qy[i]=0;}
    }
}

GFS_MOM_HD static inline int gfs_enrich_q(const GfsEnrichBasis *b,
    const GfsEnrichConfig *cfg,int n,double x,double y,double q[GFS_ENRICH_N],
    double *qx,double *qy){
    double f[GFS_ENRICH_N],fx[GFS_ENRICH_N],fy[GFS_ENRICH_N];
    if(!gfs_enrich_raw(b,cfg,x,y,f,qx?fx:0,qy?fy:0))return 0;
    gfs_enrich_transform(b,n,f,fx,fy,q,qx,qy);return 1;
}

/* Twice reorthogonalized MGS using direct residual values. No Gram
 * subtraction, diagonal ridge, or silently forced positive pivot. Excluded
 * functions remain recorded with zero diagonal; downstream pressure and
 * conserved-field error gates must judge whether exclusion is acceptable. */
GFS_MOM_HD static inline int gfs_enrich_build(GfsMomentPolygon p,
    const GfsEnrichConfig *cfg,double gx,double gy,GfsEnrichBasis *b,
    GfsEnrichState *state,double integral[4]){
    b->status=-1;b->rank=0;b->minimum_norm=DBL_MAX;b->omitted_norm=0;
    b->orthogonality_error=0;b->area_error=0;
    for(int i=0;i<GFS_ENRICH_N;i++){
        b->mean[i]=0;for(int j=0;j<GFS_ENRICH_N;j++)b->L[i][j]=0;
    }
    GfsMomentDensity geom={0};if(!gfs_moment_geometry(p,&geom))return 0;
    b->geometry.area=geom.area;b->geometry.gx=geom.gx;b->geometry.gy=geom.gy;
    b->geometry.l00=geom.l00;b->geometry.l10=geom.l10;b->geometry.l11=geom.l11;
    b->generator[0]=gx;b->generator[1]=gy;
    double U[4];if(!gfs_enrich_reference(cfg,gx+geom.gx,gy+geom.gy,U,0,0))return 0;
    double velocity=sqrt((U[1]*U[1]+U[2]*U[2])/(U[0]*U[0])+
        cfg->field.ic.gamma*cfg->field.ic.pressure/U[0]);
    b->scale[0]=U[0];b->scale[1]=b->scale[2]=U[0]*velocity;
    b->scale[3]=U[3];
    GfsEnrichQuadrature rule;if(!gfs_enrich_quadrature_init(p,cfg,b,&rule))return 0;
    int nq=gfs_enrich_nquad(p,cfg,&rule);b->L[0][0]=1;b->rank=1;
    for(int k=1;k<GFS_ENRICH_N;k++){
        double H[GFS_ENRICH_N][GFS_ENRICH_N]={0};
        for(int pass=0;pass<2;pass++){
            double a[GFS_ENRICH_N]={0};
            for(int n=0;n<nq;n++){
                double x,y,w,q[GFS_ENRICH_N],f[GFS_ENRICH_N];
                if(!gfs_enrich_quad(p,cfg,&rule,n,&x,&y,&w)){b->status=-2;return 0;}
                if(w==0)continue;
                if(!gfs_enrich_raw(b,cfg,x,y,f,0,0))return 0;
                gfs_enrich_transform(b,k,f,0,0,q,0,0);
                double residual=f[k];for(int j=0;j<k;j++)residual=fma(-b->L[k][j],q[j],residual);
                for(int j=0;j<k;j++){
                    a[j]+=w*residual*q[j];
                    if(pass==0)for(int l=0;l<=j;l++)H[j][l]+=w*q[j]*q[l]/geom.area;
                }
            }
            /* Solve with the ACTUAL previous-basis Gram, rather than
             * assuming preceding near-polynomial residuals are exactly
             * orthogonal in floating point. No diagonal regularization. */
            if(pass==0)for(int i=0;i<k;i++)if(b->L[i][i]>0)for(int j=0;j<=i;j++)if(b->L[j][j]>0){
                double v=H[i][j];for(int l=0;l<j;l++)if(b->L[l][l]>0)v=fma(-H[i][l],H[j][l],v);
                if(i==j){if(!(v>0)){b->status=-5;return 0;}H[i][i]=sqrt(v);}
                else H[i][j]=v/H[j][j];
            }
            for(int i=0;i<k;i++)if(b->L[i][i]>0){
                a[i]/=geom.area;for(int j=0;j<i;j++)if(b->L[j][j]>0)a[i]=fma(-H[i][j],a[j],a[i]);
                a[i]/=H[i][i];
            }
            for(int i=k-1;i>=0;i--)if(b->L[i][i]>0){
                for(int j=i+1;j<k;j++)if(b->L[j][j]>0)a[i]=fma(-H[j][i],a[j],a[i]);a[i]/=H[i][i];
            }
            for(int j=0;j<k;j++)if(b->L[j][j]>0)b->L[k][j]+=a[j];
        }
        double norm=0;
        for(int n=0;n<nq;n++){
            double x,y,w,q[GFS_ENRICH_N],f[GFS_ENRICH_N];
            if(!gfs_enrich_quad(p,cfg,&rule,n,&x,&y,&w)){b->status=-2;return 0;}
            if(w==0)continue;
            if(!gfs_enrich_raw(b,cfg,x,y,f,0,0))return 0;
            gfs_enrich_transform(b,k,f,0,0,q,0,0);
            double residual=f[k];for(int j=0;j<k;j++)residual=fma(-b->L[k][j],q[j],residual);
            norm+=w*residual*residual;
        }
        norm=sqrt(norm/geom.area);if(!isfinite(norm)){b->status=-3;return 0;}
        if(norm>64*DBL_EPSILON){b->L[k][k]=norm;b->rank++;b->minimum_norm=fmin(b->minimum_norm,norm);}
        else{b->L[k][k]=0;b->omitted_norm=fmax(b->omitted_norm,norm);}
    }
    double gram[GFS_ENRICH_N][GFS_ENRICH_N]={0},area=0;
    for(int c=0;c<4;c++)integral[c]=0;
    for(int n=0;n<nq;n++){
        double x,y,w,q[GFS_ENRICH_N];
        if(!gfs_enrich_quad(p,cfg,&rule,n,&x,&y,&w)){b->status=-2;return 0;}
        if(w==0)continue;
        if(!gfs_enrich_q(b,cfg,GFS_ENRICH_N,x,y,q,0,0) ||
           !gfs_enrich_reference(cfg,gx+x,gy+y,U,0,0))return 0;
        area+=w;
        for(int i=0;i<GFS_ENRICH_N;i++){
            b->mean[i]+=w*q[i];
            for(int j=0;j<=i;j++)gram[i][j]+=w*q[i]*q[j];
        }
        for(int c=0;c<4;c++)integral[c]+=w*U[c];
    }
    b->area_error=fabs(area/geom.area-1);
    for(int i=0;i<GFS_ENRICH_N;i++){
        b->mean[i]/=geom.area;
        for(int j=0;j<=i;j++){
            double target=(i==j && b->L[i][i]>0)?1.:0;
            b->orthogonality_error=fmax(b->orthogonality_error,fabs(gram[i][j]/geom.area-target));
        }
    }
    for(int c=0;c<4;c++){
        state->c[c][0]=integral[c]/geom.area;
        for(int j=1;j<GFS_ENRICH_N;j++)state->c[c][j]=j<=6+c?b->scale[c]*b->L[6+c][j]:0;
    }
    b->status=1;return 1;
}

GFS_MOM_HD static inline int gfs_enrich_evaluate(const GfsEnrichBasis *b,
    const GfsEnrichConfig *cfg,const GfsEnrichState *s,double x,double y,double U[4]){
    double q[GFS_ENRICH_N];if(!gfs_enrich_q(b,cfg,GFS_ENRICH_N,x,y,q,0,0))return 0;
    for(int c=0;c<4;c++){
        U[c]=s->c[c][0];
        for(int j=1;j<GFS_ENRICH_N;j++)U[c]=fma(s->c[c][j],q[j]-b->mean[j],U[c]);
    }
    return 1;
}

/* Physical-coordinate test functions frozen in time. This is the full Euler
 * volume term, not equilibrium-residual cancellation. The caller supplies
 * actual acceleration and one shared numerical ALE face flux. */
GFS_MOM_HD static inline int gfs_enrich_volume_point(const GfsEnrichBasis *b,
    const GfsEnrichConfig *cfg,const GfsEnrichState *s,double x,double y,
    double gx,double gy,double w,double rhs[4][GFS_ENRICH_N]){
    double q[GFS_ENRICH_N],qx[GFS_ENRICH_N],qy[GFS_ENRICH_N],U[4],fx[4],fy[4];
    if(!gfs_enrich_q(b,cfg,GFS_ENRICH_N,x,y,q,qx,qy) ||
       !gfs_enrich_evaluate(b,cfg,s,x,y,U) ||
       !gfs_dg_euler_flux(U,cfg->field.ic.gamma,fx,fy))return 0;
    double source[4]={0,U[0]*gx,U[0]*gy,U[1]*gx+U[2]*gy};
    for(int k=0;k<GFS_ENRICH_N;k++)for(int c=0;c<4;c++){
        double test=k?q[k]-b->mean[k]:1;
        rhs[c][k]+=w*(fx[c]*qx[k]+fy[c]*qy[k]+source[c]*test);
    }
    return 1;
}

GFS_MOM_HD static inline int gfs_enrich_face_point(const GfsEnrichBasis *b,
    const GfsEnrichConfig *cfg,double x,double y,const double flux[4],double w,
    double rhs[4][GFS_ENRICH_N]){
    double q[GFS_ENRICH_N];if(!gfs_enrich_q(b,cfg,GFS_ENRICH_N,x,y,q,0,0))return 0;
    for(int k=0;k<GFS_ENRICH_N;k++)for(int c=0;c<4;c++)
        rhs[c][k]-=w*flux[c]*(k?q[k]-b->mean[k]:1);
    return 1;
}

enum {
    ER_MIN_RHO,ER_MIN_P,ER_MAX_PERR,ER_L1_PERR,ER_FIELD_ERROR,ER_BOOST_PERR,
    ER_CONSERV_RHO,ER_CONSERV_PX,ER_CONSERV_PY,ER_CONSERV_E,
    ER_QUAD_RHO,ER_QUAD_PX,ER_QUAD_PY,ER_QUAD_E_OVER_THERMAL,ER_LITERAL_RHO_REL,ER_NMETRICS
};
typedef struct {
    GfsEnrichBasis basis;
    GfsEnrichState state;
    double integral[4],metrics[ER_NMETRICS];
    int status,nonpositive_rho,nonpositive_p,reserved;
} GfsEnrichAudit;
#endif
