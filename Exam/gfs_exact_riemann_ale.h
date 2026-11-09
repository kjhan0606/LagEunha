#ifndef GFS_EXACT_RIEMANN_ALE_H
#define GFS_EXACT_RIEMANN_ALE_H

#ifdef __CUDACC__
#define GFS_EXACT_HD static __host__ __device__ inline
#else
#define GFS_EXACT_HD static inline
#endif
/* Reference Godunov ALE flux for ideal-gas Euler. The same normal 1D exact
 * Riemann solution is embedded in dimensions 1--3. This iterative host/device
 * reference is not yet a production/performance-qualified solver. */
#include "gfs_fv_ale.h"
#include <float.h>

typedef struct {double p,u;int iterations,vacuum;} GfsExactStar;
typedef struct {
    GfsFvPrimitive left,right;
    GfsExactStar star;
    double normal[3],gamma,ul,ur,cl,cr;
    double rho_star_left,rho_star_right;
    double left_shock,left_head,left_tail;
    double right_shock,right_tail,right_head;
    double vacuum_left_edge,vacuum_right_edge;
    int dim,left_is_shock,right_is_shock,valid;
} GfsExactFan;

/* Inverse of the same exact ideal-gas wave curve: pressure at a piston
 * moving with outward normal speed w. s=u_n-w is positive in compression.
 * Used by shared material nodes; no empirical viscosity coefficient.
 */
GFS_EXACT_HD int gfs_exact_piston(double rho,double p,double gamma,double s,
        double *star,double *dpds){
    if(!(rho>0)||!(p>0)||!(gamma>1)||!isfinite(s))return 0;
    double c=sqrt(gamma*p/rho);
    if(s<=0){
        double y=1+(gamma-1)*s/(2*c);
        if(y<=0){*star=0;*dpds=0;return 1;}
        *star=p*pow(y,2*gamma/(gamma-1));
        *dpds=rho*c*pow(y,(gamma+1)/(gamma-1));
    }else{
        double a=.5*(gamma+1),b=a*s,r=hypot(b,2*c);
        *star=p+.5*rho*s*(b+r);
        *dpds=.5*rho*(2*b+r+s*a*b/r);
    }
    return isfinite(*star)&&isfinite(*dpds)&&*star>=0&&*dpds>=0;
}

GFS_EXACT_HD double gfs_exact_wave(double p,GfsFvPrimitive q,double c,double g){
    double gm1=g-1,gp1=g+1;
    if(p>q.p){double A=2/(gp1*q.rho),B=gm1/gp1*q.p;return (p-q.p)*sqrt(A/(p+B));}
    /* Same weak-wave arithmetic as gfs_exact_wave_impedance: pow(x,a)-1
     * loses the wave entirely for adjacent representable pressures.
     * This changes no wave curve, EOS, floor or convergence tolerance. */
    if(p==0)return -2*c/gm1;
    double r=(p-q.p)/q.p;
    double ell=r>-.5?log1p(r):log(p)-log(q.p);
    return 2*c/gm1*expm1((gm1/(2*g))*ell);
}
GFS_EXACT_HD double gfs_exact_wave_derivative(double p,GfsFvPrimitive q,double c,double g){
    double gm1=g-1,gp1=g+1;
    if(p>q.p){double A=2/(gp1*q.rho),B=gm1/gp1*q.p,s=sqrt(A/(p+B));return s*(1-.5*(p-q.p)/(p+B));}
    return 1/(q.rho*c)*pow(p/q.p,-gp1/(2*g));
}

GFS_EXACT_HD int gfs_exact_star(GfsFvPrimitive l,GfsFvPrimitive r,
        double ul,double ur,double g,GfsExactStar *s){
    double cl=sqrt(g*l.p/l.rho),cr=sqrt(g*r.p/r.rho),gm1=g-1;
    /* Form the velocity jump before adding wave speeds. Adding each wave
     * to an absolute bulk velocity first loses it in a boosted frame. */
    double jump=ur-ul;
    double vacuum_limit=2*(cl+cr)/gm1;
    if(jump>=vacuum_limit){s->p=0;s->u=.5*(ul+ur)+ (cl-cr)/gm1;s->iterations=0;s->vacuum=1;return 1;}
    /* Both waves rarefy when the root lies below both input pressures.
     * Solve that branch analytically: avoids hundreds of bisections near
     * vacuum without introducing a pressure floor or changing the EOS.
     * Normalize by min(pL,pR) so the expression is pressure-unit invariant. */
    double pm=fmin(l.p,r.p);
    double fm=gfs_exact_wave(pm,l,cl,g)+gfs_exact_wave(pm,r,cr,g)+jump;
    if(fm>=0){
        double a=gm1/(2*g);
        double remaining=.5*gm1*(vacuum_limit-jump);
        double denom=cl*pow(pm/l.p,a)+cr*pow(pm/r.p,a);
        /* In the weak-wave limit remaining and denom nearly coincide.
         * Use F(pm) to evaluate their difference without subtracting sound
         * speeds. The deep-rarefaction form retains the resolved distance
         * from vacuum. Both are the same analytic two-rarefaction root. */
        double decrement=.5*gm1*fm/denom;
        double p;
        if(decrement<.5){
            double exponent=log1p(-decrement)/a;
            /* For gamma close to one a modest decrement can still imply
             * a tiny, representable pressure. Adding expm1 to one would
             * cancel it to zero; retain that root multiplicatively. */
            p=exponent>-.5?pm+pm*expm1(exponent):pm*exp(exponent);
        }else p=pm*pow(remaining/denom,1/a);
        if(!(p>0)||!isfinite(p))return 0;
        s->p=p;s->iterations=0;s->vacuum=0;
        s->u=.5*(ul+ur)+.5*(gfs_exact_wave(p,r,cr,g)-gfs_exact_wave(p,l,cl,g));
        return isfinite(s->u);
    }
    double lo=0,hi=fmax(l.p,r.p),flo=jump-2*(cl+cr)/gm1;
    double fhi=gfs_exact_wave(hi,l,cl,g)+gfs_exact_wave(hi,r,cr,g)+jump;
    for(int k=0;fhi<=0&&k<1024;k++){
        if(hi>DBL_MAX/4)return 0;
        hi*=2;fhi=gfs_exact_wave(hi,l,cl,g)+gfs_exact_wave(hi,r,cr,g)+jump;
    }
    if(!(fhi>0)||!isfinite(fhi)||!(hi>0))return 0;
    double pp=.5*(l.p+r.p)-.125*jump*(l.rho+r.rho)*(cl+cr);
    /* hi is a valid positive-pressure endpoint. For a weak wave pp may
     * round exactly to hi; discarding it for hi/2 loses a near-exact guess
     * and lets the sound-speed residual tolerance dominate the signal. */
    double p=(pp>lo&&pp<=hi)?pp:.5*(lo+hi);
    for(int it=1;it<=64;it++){
        double f=gfs_exact_wave(p,l,cl,g)+gfs_exact_wave(p,r,cr,g)+jump;
        double fp=gfs_exact_wave_derivative(p,l,cl,g)+gfs_exact_wave_derivative(p,r,cr,g);
        if(!isfinite(f)||!(fp>0)||!isfinite(fp))return 0;
        if(f>0)hi=p;else{lo=p;flo=f;}
        /* The stopping scale must be invariant to a common bulk velocity. */
        if(fabs(f)<=2e-14*(2*(cl+cr)+fabs(jump))||
           /* Pressure units must not set an absolute convergence floor.
            * Near vacuum, a 1-unit floor can accept a bracket much larger
            * than p_star while the wave residual is still appreciable. */
           hi-lo<=4e-14*fmax(DBL_MIN,hi)){
            s->p=p;
            s->u=.5*(ul+ur)+.5*(gfs_exact_wave(p,r,cr,g)-gfs_exact_wave(p,l,cl,g));
            s->iterations=it;s->vacuum=0;
            return s->p>0&&isfinite(s->p)&&isfinite(s->u);
        }
        double pn=p-f/fp;
        if(!(pn>lo&&pn<hi)||!isfinite(pn))pn=.5*(lo+hi);
        p=pn;
    }
    (void)flo;
    return 0;
}

GFS_EXACT_HD double gfs_exact_star_density(GfsFvPrimitive q,double p,double g){
    double gm1=g-1,gp1=g+1;
    if(p>q.p){double ratio=p/q.p,mu=gm1/gp1;return q.rho*(ratio+mu)/(mu*ratio+1);}
    return q.rho*pow(p/q.p,1/g);
}

GFS_EXACT_HD int gfs_exact_riemann_prepare(int dim,GfsFvPrimitive l,GfsFvPrimitive r,
        const double normal[3],double gamma,GfsExactFan *fan){
    if(!fan||!gfs_fv_input(dim,l)||!gfs_fv_input(dim,r)||!(gamma>1)||!isfinite(gamma))return 0;
    GfsExactFan f={0};f.dim=dim;f.left=l;f.right=r;f.gamma=gamma;
    double norm=0;
    for(int a=0;a<dim;a++){
        if(!isfinite(normal[a]))return 0;
        f.normal[a]=normal[a];norm+=normal[a]*normal[a];
        f.ul+=l.u[a]*normal[a];f.ur+=r.u[a]*normal[a];
    }
    if(fabs(norm-1)>128*DBL_EPSILON)return 0;
    f.cl=sqrt(gamma*l.p/l.rho);f.cr=sqrt(gamma*r.p/r.rho);
    if(!(f.cl>0)||!(f.cr>0)||!isfinite(f.cl)||!isfinite(f.cr))return 0;
    if(!gfs_exact_star(l,r,f.ul,f.ur,gamma,&f.star))return 0;
    if(f.star.vacuum){
        f.vacuum_left_edge=f.ul+2*f.cl/(gamma-1);
        f.vacuum_right_edge=f.ur-2*f.cr/(gamma-1);
        if(!isfinite(f.vacuum_left_edge)||!isfinite(f.vacuum_right_edge))return 0;
    }else{
        f.rho_star_left=gfs_exact_star_density(l,f.star.p,gamma);
        f.rho_star_right=gfs_exact_star_density(r,f.star.p,gamma);
        f.left_is_shock=f.star.p>l.p;f.right_is_shock=f.star.p>r.p;
        if(f.left_is_shock)
            f.left_shock=f.ul-f.cl*sqrt((gamma+1)/(2*gamma)*(f.star.p/l.p)+(gamma-1)/(2*gamma));
        else{
            f.left_head=f.ul-f.cl;
            f.left_tail=f.star.u-f.cl*pow(f.star.p/l.p,(gamma-1)/(2*gamma));
        }
        if(f.right_is_shock)
            f.right_shock=f.ur+f.cr*sqrt((gamma+1)/(2*gamma)*(f.star.p/r.p)+(gamma-1)/(2*gamma));
        else{
            f.right_head=f.ur+f.cr;
            f.right_tail=f.star.u+f.cr*pow(f.star.p/r.p,(gamma-1)/(2*gamma));
        }
        if(!(f.rho_star_left>0)||!(f.rho_star_right>0)||
           !isfinite(f.rho_star_left)||!isfinite(f.rho_star_right))return 0;
    }
    f.valid=1;*fan=f;return 1;
}

GFS_EXACT_HD int gfs_exact_riemann_sample(const GfsExactFan *fan,double wn,
        double flux[5],GfsExactStar *star){
    if(!fan||!fan->valid||!flux||!isfinite(wn))return 0;
    int d=fan->dim;double gamma=fan->gamma,xi=wn;
    GfsFvPrimitive q;
    if(fan->star.vacuum){
        if(xi>=fan->vacuum_left_edge&&xi<=fan->vacuum_right_edge){
            for(int c=0;c<5;c++)flux[c]=0;
            if(star)*star=fan->star;
            return 1;
        }
        if(xi<fan->vacuum_left_edge){
            if(xi<=fan->ul-fan->cl)q=fan->left;
            else{
                double u=2/(gamma+1)*(fan->cl+.5*(gamma-1)*fan->ul+xi);
                double c=2/(gamma+1)*(fan->cl+.5*(gamma-1)*(fan->ul-xi));
                if(!(c>0))return 0;
                q=fan->left;q.rho*=pow(c/fan->cl,2/(gamma-1));
                q.p*=pow(c/fan->cl,2*gamma/(gamma-1));
                for(int a=0;a<d;a++)q.u[a]+=(u-fan->ul)*fan->normal[a];
            }
        }else{
            if(xi>=fan->ur+fan->cr)q=fan->right;
            else{
                double u=2/(gamma+1)*(-fan->cr+.5*(gamma-1)*fan->ur+xi);
                double c=2/(gamma+1)*(fan->cr-.5*(gamma-1)*(fan->ur-xi));
                if(!(c>0))return 0;
                q=fan->right;q.rho*=pow(c/fan->cr,2/(gamma-1));
                q.p*=pow(c/fan->cr,2*gamma/(gamma-1));
                for(int a=0;a<d;a++)q.u[a]+=(u-fan->ur)*fan->normal[a];
            }
        }
    }else if(xi<=fan->star.u){
        if(fan->left_is_shock){
            q=fan->left;
            if(xi>=fan->left_shock){
                q.rho=fan->rho_star_left;q.p=fan->star.p;
                for(int a=0;a<d;a++)q.u[a]+=(fan->star.u-fan->ul)*fan->normal[a];
            }
        }else if(xi<=fan->left_head)q=fan->left;
        else if(xi>=fan->left_tail){
            q=fan->left;q.rho=fan->rho_star_left;q.p=fan->star.p;
            for(int a=0;a<d;a++)q.u[a]+=(fan->star.u-fan->ul)*fan->normal[a];
        }else{
            double u=2/(gamma+1)*(fan->cl+.5*(gamma-1)*fan->ul+xi);
            double c=2/(gamma+1)*(fan->cl+.5*(gamma-1)*(fan->ul-xi));
            if(!(c>0))return 0;
            q=fan->left;q.rho*=pow(c/fan->cl,2/(gamma-1));
            q.p*=pow(c/fan->cl,2*gamma/(gamma-1));
            for(int a=0;a<d;a++)q.u[a]+=(u-fan->ul)*fan->normal[a];
        }
    }else{
        if(fan->right_is_shock){
            q=fan->right;
            if(xi<=fan->right_shock){
                q.rho=fan->rho_star_right;q.p=fan->star.p;
                for(int a=0;a<d;a++)q.u[a]+=(fan->star.u-fan->ur)*fan->normal[a];
            }
        }else if(xi>=fan->right_head)q=fan->right;
        else if(xi<=fan->right_tail){
            q=fan->right;q.rho=fan->rho_star_right;q.p=fan->star.p;
            for(int a=0;a<d;a++)q.u[a]+=(fan->star.u-fan->ur)*fan->normal[a];
        }else{
            double u=2/(gamma+1)*(-fan->cr+.5*(gamma-1)*fan->ur+xi);
            double c=2/(gamma+1)*(fan->cr-.5*(gamma-1)*(fan->ur-xi));
            if(!(c>0))return 0;
            q=fan->right;q.rho*=pow(c/fan->cr,2/(gamma-1));
            q.p*=pow(c/fan->cr,2*gamma/(gamma-1));
            for(int a=0;a<d;a++)q.u[a]+=(u-fan->ur)*fan->normal[a];
        }
    }
    if(!(q.rho>0)||!(q.p>0)||!isfinite(q.rho)||!isfinite(q.p))return 0;
    double U[5],F[5];gfs_fv_state_flux(d,q,fan->normal,gamma,U,F);
    for(int c=0;c<5;c++){flux[c]=F[c]-wn*U[c];if(!isfinite(flux[c]))return 0;}
    if(star)*star=fan->star;
    return 1;
}

GFS_EXACT_HD int gfs_exact_riemann_ale_flux(int dim,GfsFvPrimitive l,GfsFvPrimitive r,
        const double normal[3],double gamma,double wn,double flux[5],GfsExactStar *star){
    GfsExactFan fan;
    return gfs_exact_riemann_prepare(dim,l,r,normal,gamma,&fan)&&
           gfs_exact_riemann_sample(&fan,wn,flux,star);
}
#undef GFS_EXACT_HD
#endif
