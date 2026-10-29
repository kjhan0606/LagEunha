#ifndef GFS_MOMENT_DENSITY_H
#define GFS_MOMENT_DENSITY_H
#include "gfs_moment.h"
/* Positive exponential-quadratic density on the actual convex polygon.
 * Realizability is a nonlinear solve for all five normalized moments, NOT
 * a reset to a uniform cell. Failure rejects the state; no density clipping,
 * Hessian regularization, moment projection or added thermal energy.
 */
typedef struct { const double *x,*y; int n,refine; } GfsMomentPolygon;
typedef struct {
    double area,gx,gy,l00,l10,l11,lambda[5],logZ,residual;
    int iterations,status;
} GfsMomentDensity;

GFS_MOM_HD static inline int gfs_moment_geometry(GfsMomentPolygon p,
    GfsMomentDensity *d){
    if(p.n<3 || (p.refine!=0 && p.refine!=1)) return 0;
    double a=0,x=0,y=0,xx=0,xy=0,yy=0;
    for(int i=0;i<p.n;i++){
        int j=(i+1)%p.n; double xi=p.x[i],yi=p.y[i],xj=p.x[j],yj=p.y[j];
        double z=xi*yj-xj*yi;
        a+=z; x+=(xi+xj)*z; y+=(yi+yj)*z;
        xx+=(xi*xi+xi*xj+xj*xj)*z;
        xy+=(2*xi*yi+xi*yj+xj*yi+2*xj*yj)*z;
        yy+=(yi*yi+yi*yj+yj*yj)*z;
    }
    if(!(a>0) || !isfinite(a)) return 0;
    d->area=.5*a; d->gx=x/(3*a); d->gy=y/(3*a);
    double vxx=xx/(6*a)-d->gx*d->gx;
    double vxy=xy/(12*a)-d->gx*d->gy;
    double vyy=yy/(6*a)-d->gy*d->gy;
    if(!(vxx>0)) return 0;
    d->l00=sqrt(vxx); d->l10=vxy/d->l00;
    double pivot=fma(-d->l10,d->l10,vyy);
    if(!(pivot>0) || !isfinite(pivot)) return 0;
    d->l11=sqrt(pivot);
    /* Native raw polygons may repeat an exactly identical adjacent vertex.
     * This is an inactive edge of zero measure, not a negative-area triangle.
     * All nonzero edges require strictly positive centroid-fan areas. */
    for(int i=0;i<p.n;i++){
        int j=(i+1)%p.n;
        if(p.x[i]==p.x[j] && p.y[i]==p.y[j])continue;
        if(!((p.x[i]-d->gx)*(p.y[j]-d->gy)-
             (p.y[i]-d->gy)*(p.x[j]-d->gx)>0)) return 0;
    }
    return 1;
}

GFS_MOM_HD static inline int gfs_moment_nquad(GfsMomentPolygon p){
    return 7*p.n*(p.refine?4:1);
}

/* Positive 7-point degree-five Dunavant rule. Optional one-level triangle
 * subdivision is a quadrature check/choice, never a change to the mesh. */
GFS_MOM_HD static inline void gfs_moment_quad(GfsMomentPolygon p,
    const GfsMomentDensity *d,int k,double *x,double *y,double *weight){
    int node=k%7,t=k/7,sub=0;
    if(p.refine){sub=t%4;t/=4;}
    int j=(t+1)%p.n;
    double ax=d->gx,ay=d->gy,bx=p.x[t],by=p.y[t],cx=p.x[j],cy=p.y[j];
    double area=.5*((bx-ax)*(cy-ay)-(by-ay)*(cx-ax));
    if(p.refine){
        double abx=.5*(ax+bx),aby=.5*(ay+by),bcx=.5*(bx+cx),bcy=.5*(by+cy);
        double cax=.5*(cx+ax),cay=.5*(cy+ay);
        if(sub==0){bx=abx;by=aby;cx=cax;cy=cay;}
        else if(sub==1){ax=abx;ay=aby;cx=bcx;cy=bcy;}
        else if(sub==2){ax=cax;ay=cay;bx=bcx;by=bcy;}
        else {ax=abx;ay=aby;bx=bcx;by=bcy;cx=cax;cy=cay;}
        area*=.25;
    }
    double a,b,c,w;
    if(node==0){a=b=c=1./3.;w=.225;}
    else {
        double special=node<=3?.05971587178976982046:.79742698535308732240;
        double other=node<=3?.47014206410511508977:.10128650732345633880;
        int v=(node-1)%3;
        a=v==0?special:other;b=v==1?special:other;c=v==2?special:other;
        w=node<=3?.13239415278850618074:.12593918054482715260;
    }
    *x=a*ax+b*bx+c*cx;*y=a*ay+b*by+c*cy;*weight=area*w;
}

GFS_MOM_HD static inline void gfs_moment_features(const GfsMomentDensity *d,
    double x,double y,double f[5]){
    double u=(x-d->gx)/d->l00,v=(y-d->gy-d->l10*u)/d->l11;
    f[0]=u;f[1]=v;f[2]=u*u;f[3]=u*v;f[4]=v*v;
}

GFS_MOM_HD static inline double gfs_moment_dot5(const double a[5],const double b[5]){
    double z=0;for(int k=0;k<5;k++)z+=a[k]*b[k];return z;
}

GFS_MOM_HD static inline int gfs_moment_partition(GfsMomentPolygon p,
    const GfsMomentDensity *d,const double lambda[5],double *logZ,
    double mean[5],double H[25]){
    int n=gfs_moment_nquad(p);double peak=-DBL_MAX;
    for(int k=0;k<n;k++){
        double x,y,w,f[5];gfs_moment_quad(p,d,k,&x,&y,&w);
        if(w==0)continue;
        gfs_moment_features(d,x,y,f);double e=gfs_moment_dot5(lambda,f);
        if(!isfinite(e))return 0;if(e>peak)peak=e;
    }
    double z=0;for(int a=0;a<5;a++)mean[a]=0;
    if(H)for(int a=0;a<25;a++)H[a]=0;
    for(int k=0;k<n;k++){
        double x,y,w,f[5];gfs_moment_quad(p,d,k,&x,&y,&w);gfs_moment_features(d,x,y,f);
        if(w==0)continue;
        w*=exp(gfs_moment_dot5(lambda,f)-peak);z+=w;
        for(int a=0;a<5;a++){
            mean[a]+=w*f[a];
            if(H)for(int b=0;b<=a;b++)H[5*a+b]+=w*f[a]*f[b];
        }
    }
    if(!(z>0) || !isfinite(z))return 0;
    *logZ=peak+log(z);
    for(int a=0;a<5;a++)mean[a]/=z;
    if(H)for(int a=0;a<5;a++)for(int b=0;b<=a;b++){
        double v=H[5*a+b]/z-mean[a]*mean[b];H[5*a+b]=H[5*b+a]=v;
    }
    return isfinite(*logZ);
}

GFS_MOM_HD static inline int gfs_moment_solve5(double H[25],double b[5]){
    for(int i=0;i<5;i++)for(int j=0;j<=i;j++){
        double v=H[5*i+j];for(int k=0;k<j;k++)v-=H[5*i+k]*H[5*j+k];
        if(i==j){if(!(v>0) || !isfinite(v))return 0;H[5*i+j]=sqrt(v);}
        else H[5*i+j]=v/H[5*j+j];
    }
    for(int i=0;i<5;i++){
        for(int j=0;j<i;j++)b[i]-=H[5*i+j]*b[j];b[i]/=H[5*i+i];
    }
    for(int i=4;i>=0;i--){
        for(int j=i+1;j<5;j++)b[i]-=H[5*j+i]*b[j];b[i]/=H[5*i+i];
        if(!isfinite(b[i]))return 0;
    }
    return 1;
}

/* State q,Q are relative to fixed O; polygon vertices relative to generator,
 * whose offset from the appropriate image of O is (rx,ry). */
GFS_MOM_HD static inline int gfs_moment_density_fit_seeded(GfsMomentPolygon p,
    const GfsMomentState *s,double rx,double ry,const double seed[5],GfsMomentDensity *d){
    d->status=0;d->iterations=0;d->residual=DBL_MAX;
    if(!gfs_moment_geometry(p,d)){d->status=-1;return 0;}
    GfsMomentPrimitive r;if(!gfs_moment_primitive(s,&r)){d->status=-2;return 0;}
    double m=s->a[GFM_M],a=1/d->l00,b=-d->l10/(d->l00*d->l11),c=1/d->l11;
    double dx=r.c[0]-rx-d->gx,dy=r.c[1]-ry-d->gy;
    double mx=a*dx,my=b*dx+c*dy;
    double target[5]={mx,my,a*a*r.S[0]/m+mx*mx,
        a*(b*r.S[0]+c*r.S[1])/m+mx*my,
        (b*b*r.S[0]+2*b*c*r.S[1]+c*c*r.S[2])/m+my*my};
    for(int k=0;k<5;k++)d->lambda[k]=seed?seed[k]:0;
    for(int it=0;it<32;it++){
        double mean[5],H[25],grad[5];
        if(!gfs_moment_partition(p,d,d->lambda,&d->logZ,mean,H)){d->status=-3;return 0;}
        double err=0;for(int k=0;k<5;k++){
            grad[k]=mean[k]-target[k];err=fmax(err,fabs(grad[k]));
        }
        d->iterations=it;d->residual=err;
        if(err<=5.e-12){d->status=1;return 1;}
        double obj=d->logZ-gfs_moment_dot5(d->lambda,target),step[5];
        for(int k=0;k<5;k++)step[k]=grad[k];
        if(!gfs_moment_solve5(H,step)){d->status=-4;return 0;}
        double descent=gfs_moment_dot5(grad,step);if(!(descent>0)){d->status=-5;return 0;}
        int accepted=0;double alpha=1;
        for(int ls=0;ls<24;ls++,alpha*=.5){
            double trial[5],tm[5],tz;
            for(int k=0;k<5;k++)trial[k]=d->lambda[k]-alpha*step[k];
            if(!gfs_moment_partition(p,d,trial,&tz,tm,0))continue;
            double te=0;for(int k=0;k<5;k++)te=fmax(te,fabs(tm[k]-target[k]));
            double tobj=tz-gfs_moment_dot5(trial,target);
            if(tobj<=obj-1.e-4*alpha*descent || te<err*(1-.1*alpha)){
                for(int k=0;k<5;k++)d->lambda[k]=trial[k];accepted=1;break;
            }
        }
        if(!accepted){d->status=-6;return 0;}
    }
    d->status=-7;return 0;
}

/* Original zero-initialized solve remains the default. A warm seed is only
 * an initial guess, expressed in this polygon's current normalized basis.
 * The same residual and positivity criteria apply; a failed or nonfinite
 * guess retries the original cold solve without projecting the state. */
GFS_MOM_HD static inline int gfs_moment_density_fit(GfsMomentPolygon p,
    const GfsMomentState *s,double rx,double ry,GfsMomentDensity *d){
    return gfs_moment_density_fit_seeded(p,s,rx,ry,0,d);
}
GFS_MOM_HD static inline int gfs_moment_density_fit_warm(GfsMomentPolygon p,
    const GfsMomentState *s,double rx,double ry,const double seed[5],GfsMomentDensity *d){
    int finite=seed!=0;
    if(seed)for(int k=0;k<5;k++)if(!isfinite(seed[k]))finite=0;
    if(finite && gfs_moment_density_fit_seeded(p,s,rx,ry,seed,d))return 1;
    return gfs_moment_density_fit(p,s,rx,ry,d);
}

GFS_MOM_HD static inline double gfs_moment_density_value(const GfsMomentDensity *d,
    double mass,double x,double y){
    double f[5];gfs_moment_features(d,x,y,f);
    return mass*exp(gfs_moment_dot5(d->lambda,f)-d->logZ);
}
#endif
