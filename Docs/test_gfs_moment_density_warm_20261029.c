#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Exam/gfs_moment_density.h"

static void check(int condition,const char *message){if(!condition){fprintf(stderr,"%s\n",message);exit(1);}}
static GfsMomentState make_state(GfsMomentPolygon poly,const double lambda[5]){
    GfsMomentDensity geom;
    check(gfs_moment_geometry(poly,&geom),"geometry");
    double logz,mean[5];
    check(gfs_moment_partition(poly,&geom,lambda,&logz,mean,0),"partition");
    GfsMomentState state;gfs_moment_zero(&state);
    for(int k=0;k<gfs_moment_nquad(poly);k++){
        double x,y,w,f[5];gfs_moment_quad(poly,&geom,k,&x,&y,&w);
        gfs_moment_features(&geom,x,y,f);
        double dm=w*exp(gfs_moment_dot5(lambda,f)-logz);
        gfs_moment_accumulate(&state,dm,x,y,.3+.2*x-.1*y,-.4+.3*x+.1*y,.4*dm);
    }
    return state;
}
int main(void){
    const double x[]={-1.0,.7,1.1,.4,-.9},y[]={-.4,-.8,.2,.9,.7};
    GfsMomentPolygon poly={x,y,5,1};
    const double first[]={.35,-.23,.09,.04,-.08};
    const double second[]={.351,-.229,.091,.039,-.081};
    GfsMomentState a=make_state(poly,first),b=make_state(poly,second);
    GfsMomentDensity old,cold,warm,bad;
    check(gfs_moment_density_fit(poly,&a,0,0,&old),"first cold fit");
    check(gfs_moment_density_fit(poly,&b,0,0,&cold),"second cold fit");
    check(gfs_moment_density_fit_warm(poly,&b,0,0,old.lambda,&warm),"warm fit");
    check(warm.iterations<cold.iterations,"warm start did not reduce Newton iterations");
    check(warm.residual<=5e-12 && cold.residual<=5e-12,"residual tolerance");
    for(int k=0;k<5;k++)check(fabs(warm.lambda[k]-cold.lambda[k])<1e-9,"cold/warm mismatch");
    double invalid[]={NAN,0,0,0,0};
    check(gfs_moment_density_fit_warm(poly,&b,0,0,invalid,&bad),"invalid seed fallback");
    for(int k=0;k<5;k++)check(bad.lambda[k]==cold.lambda[k],"fallback did not use cold solution");
    double extreme[]={1e6,-1e6,1e6,1e6,-1e6};
    check(!gfs_moment_density_fit_seeded(poly,&b,0,0,extreme,&bad),"extreme seed should exercise retry");
    check(gfs_moment_density_fit_warm(poly,&b,0,0,extreme,&bad),"finite bad seed fallback");
    for(int k=0;k<5;k++)check(bad.lambda[k]==cold.lambda[k],"finite fallback differs from cold");
    GfsMomentDensity in_place=old;
    check(gfs_moment_density_fit_warm(poly,&b,0,0,in_place.lambda,&in_place),"in-place seed");
    for(int k=0;k<5;k++)check(fabs(in_place.lambda[k]-cold.lambda[k])<1e-9,"in-place mismatch");
    double sx[5],sy[5];
    for(int k=0;k<5;k++){sx[k]=1.02*x[k]+.03*y[k]+.2;sy[k]=.01*x[k]+.98*y[k]-.1;}
    GfsMomentPolygon moved={sx,sy,5,1};
    GfsMomentState moved_state=make_state(moved,second);
    GfsMomentDensity moved_cold,moved_warm;
    check(gfs_moment_density_fit(moved,&moved_state,0,0,&moved_cold),"moved polygon cold fit");
    check(gfs_moment_density_fit_warm(moved,&moved_state,0,0,old.lambda,&moved_warm),"moved polygon warm fit");
    for(int k=0;k<5;k++)check(fabs(moved_warm.lambda[k]-moved_cold.lambda[k])<1e-9,"moved polygon mismatch");
    GfsMomentState impossible=b;impossible.a[GFM_QXX]=-1;
    check(!gfs_moment_density_fit_warm(poly,&impossible,0,0,old.lambda,&bad),"invalid state accepted");
    printf("moment density fit: cold=%d warm=%d iterations; residual %.3g %.3g; finite/nonfinite fallback, alias and moved polygon pass\n",
           cold.iterations,warm.iterations,cold.residual,warm.residual);
    return 0;
}
