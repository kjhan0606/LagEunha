#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "gfs_moment_initial.h"
#include "gfs_ale_flux.h"

static void require(bool ok,const char *message){
    if(!ok){fprintf(stderr,"FAIL: %s\n",message);exit(1);}
}
static double maxabs(const GfsMomentState &s){
    double z=0;for(int k=0;k<GFM_NSTATE;k++)z=fmax(z,fabs(s.a[k]));return z;
}
int main(){
    double xx[6]={-.8,.45,1.,.7,-.4,-.9},yy[6]={-.4,-.7,-.1,.7,.8,.2};
    GfsMomentPolygon poly={xx,yy,6,1};GfsMomentDensity geom={};
    require(gfs_moment_geometry(poly,&geom),"polygon geometry");
    double true_lambda[5]={-.15,.1,.07,-.03,-.04};
    GfsMomentState s={};double physical_U=0;
    const double u=1.e-5,gamma=5./3.,B[4]={.03,-.7,.9,-.04};
    for(int k=0;k<gfs_moment_nquad(poly);k++){
        double x,y,w,f[5];gfs_moment_quad(poly,&geom,k,&x,&y,&w);
        gfs_moment_features(&geom,x,y,f);
        double rho=exp(gfs_moment_dot5(true_lambda,f));
        double vx=.2+B[0]*x+B[1]*y,vy=-.1+B[2]*x+B[3]*y;
        gfs_moment_accumulate(&s,w*rho,x,y,vx,vy,w*rho*u);physical_U+=w*rho*u;
    }
    GfsMomentPrimitive r={};require(gfs_moment_primitive(&s,&r),"admissible affine state");
    double berr=0;for(int k=0;k<4;k++)berr=fmax(berr,fabs(r.B[k]-B[k]));
    require(berr<2.e-14,"affine B exactness");
    require(fabs(r.U/physical_U-1)<2.e-10,"affine thermal closure");
    GfsMomentDensity fit={};require(gfs_moment_density_fit(poly,&s,0,0,&fit),"nonuniform density realization");
    double density_error=0;
    for(int k=0;k<gfs_moment_nquad(poly);k++){
        double x,y,w,f[5];gfs_moment_quad(poly,&geom,k,&x,&y,&w);
        gfs_moment_features(&geom,x,y,f);double exact=exp(gfs_moment_dot5(true_lambda,f));
        double actual=gfs_moment_density_value(&fit,s.a[GFM_M],x,y);
        density_error=fmax(density_error,fabs(actual/exact-1));
    }
    require(density_error<2.e-10,"positive density pointwise recovery");
    GfsMomentState shifted=gfs_moment_shift(s,.2,-.1),boosted=gfs_moment_boost(s,3,-2);
    GfsMomentPrimitive sr={},br={};
    require(gfs_moment_primitive(&shifted,&sr) && gfs_moment_primitive(&boosted,&br),"origin/boost admissibility");
    require(fabs(sr.U/r.U-1)<2.e-10 && fabs(br.U/r.U-1)<2.e-9,"origin/boost thermal invariance");
    require(fabs(gfs_moment_Lz(&s,1.2,-.3)-gfs_moment_Lz(&shifted,1.4,-.4))<2.e-14,"origin angular invariance");
    GfsMomentState source={},reference={};gfs_moment_euler_source(&s,&r,gamma,&source);
    reference.a[GFM_QX]=s.a[GFM_PX];reference.a[GFM_QY]=s.a[GFM_PY];
    reference.a[GFM_QXX]=2*s.a[GFM_DXX];
    reference.a[GFM_QXY]=s.a[GFM_DXY]+s.a[GFM_DYX];reference.a[GFM_QYY]=2*s.a[GFM_DYY];
    for(int k=0;k<gfs_moment_nquad(poly);k++){
        double x,y,w,f[5];gfs_moment_quad(poly,&geom,k,&x,&y,&w);gfs_moment_features(&geom,x,y,f);
        double rho=exp(gfs_moment_dot5(true_lambda,f));
        double vx=.2+B[0]*x+B[1]*y,vy=-.1+B[2]*x+B[3]*y;
        reference.a[GFM_DXX]+=w*rho*(vx*vx+(gamma-1)*u);
        reference.a[GFM_DXY]+=w*rho*vx*vy;reference.a[GFM_DYX]+=w*rho*vy*vx;
        reference.a[GFM_DYY]+=w*rho*(vy*vy+(gamma-1)*u);
    }
    double source_error=0;for(int k=0;k<GFM_NSTATE;k++)source_error=fmax(source_error,fabs(source.a[k]-reference.a[k]));
    require(source_error<2.e-14,"Euler volume source against direct quadrature");
    GfsMomentState fi={},fj={};double flux[4]={.2,.3,-.7,.4};
    gfs_moment_face_rhs(&fi,.7-1.2,.8+.4,flux,.31);
    gfs_moment_face_rhs(&fj,.7+.1,.8-.3,flux,-.31);
    double pair_angular=fabs(gfs_moment_Lz(&fi,1.2,-.4)+gfs_moment_Lz(&fj,-.1,.3));
    require(pair_angular<2.e-16,"canonical pair orbital+intrinsic angular balance");
    GfsMomentState pair=gfs_moment_shift(fi,-1.2,.4),other=gfs_moment_shift(fj,.1,-.3);
    for(int k=0;k<GFM_NSTATE;k++)pair.a[k]+=other.a[k];
    require(maxabs(pair)<2.e-16,"all origin-transformed face moments conservative");
    GfsMomentState impossible=s;impossible.a[GFM_QXX]=100*s.a[GFM_M];
    GfsMomentDensity fail={};require(!gfs_moment_density_fit(poly,&impossible,0,0,&fail),"support-unrealizable SPD moments rejected");
    /* A generator outside its convex power cell still has a valid centroid
     * fan; no generator-based negative quadrature weights are used. */
    double tx[3]={2,3,2},ty[3]={0,0,1};GfsMomentPolygon outside={tx,ty,3,1};
    GfsMomentDensity outside_geom={};require(gfs_moment_geometry(outside,&outside_geom),"outside-generator polygon");
    double area_sum=0;for(int k=0;k<gfs_moment_nquad(outside);k++){
        double x,y,w;gfs_moment_quad(outside,&outside_geom,k,&x,&y,&w);
        require(w>0,"positive power-cell quadrature");area_sum+=w;
    }
    require(fabs(area_sum-.5)<2.e-15,"quadrature area");
    double repeated_x[4]={2,3,3,2},repeated_y[4]={0,0,0,1};
    GfsMomentPolygon repeated={repeated_x,repeated_y,4,1};
    GfsMomentDensity repeated_geom={};require(gfs_moment_geometry(repeated,&repeated_geom),"inactive repeated native vertex");
    double repeated_area=0;for(int k=0;k<gfs_moment_nquad(repeated);k++){
        double x,y,w;gfs_moment_quad(repeated,&repeated_geom,k,&x,&y,&w);
        require(w>=0,"nonnegative quadrature with an inactive edge");repeated_area+=w;
    }
    require(fabs(repeated_area-.5)<2.e-15,"inactive edge preserves exact measure");
    /* Actual HLLC/HLL flux, not an assumed SPD matrix. At frozen primitives,
     * mass flux must decrease monotonically with outward mesh normal speed.
     * Derivative checks do not yet constitute the MPI nonlinear weight solve. */
    double worst_slope=-DBL_MAX,min_negative_derivative=DBL_MAX;int modes=0,samples=0;
    for(int variant=0;variant<6;variant++){
        double pressure=variant%2?1.e-6:1.;
        GfsAlePrimitive l={1.2,2*pressure,.3,.8},rr={.7,pressure,-.2,-.1};
        if(variant>=2){rr.vx=.31;rr.vy=.79;}
        if(variant>=4){rr=l;}
        double prev=DBL_MAX;
        for(int k=0;k<=200;k++){
            double wn=-4+.04*k;
            GfsAleFlux z=gfs_ale_flux(l,rr,.6,.8,.6*wn,.8*wn,gamma);
            require(z.valid,"ALE monotonicity valid flux");modes|=1<<z.mode;samples++;
            if(k)require(z.flux[0]<prev,"actual massflux strictly mesh-speed monotone");
            prev=z.flux[0];double dw=1.e-6;
            GfsAleFlux lo=gfs_ale_flux(l,rr,.6,.8,.6*(wn-dw),.8*(wn-dw),gamma);
            GfsAleFlux hi=gfs_ale_flux(l,rr,.6,.8,.6*(wn+dw),.8*(wn+dw),gamma);
            require(lo.valid && hi.valid,"finite-difference flux validity");
            double slope=(hi.flux[0]-lo.flux[0])/(2*dw);
            worst_slope=fmax(worst_slope,slope);min_negative_derivative=fmin(min_negative_derivative,-slope);
            require(slope<0,"finite-difference massflux derivative negative");
        }
    }
    printf("{\"status\":\"pass\",\"payload_bytes\":%zu,\"state_components\":%d,\"B_error\":%.17g,\"density_pointwise_relative_error\":%.17g,\"moment_residual\":%.17g,\"density_iterations\":%d,\"source_error\":%.17g,\"pair_angular_residual\":%.17g,\"flux_samples\":%d,\"mode_bitmask\":%d,\"largest_massflux_derivative\":%.17g,\"smallest_minus_derivative\":%.17g}\n",
        sizeof(GfsMomentPayload),GFM_NSTATE,berr,density_error,fit.residual,fit.iterations,source_error,pair_angular,samples,modes,worst_slope,min_negative_derivative);
    return 0;
}
