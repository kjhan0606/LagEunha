#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Exam/gfs_moment_material_pressure.h"
static void ok(int a,const char *s){if(!a){fprintf(stderr,"%s\n",s);exit(1);}}
static void near(double a,double b,const char *s){if(!isfinite(a)||fabs(a-b)>2e-12*(1+fabs(b))){fprintf(stderr,"%s %.17g %.17g\n",s,a,b);exit(2);}}
static void check_square(double omega){
    const double gamma=5./3.,pressure=.6;
    GfsMomentState s,rhs;gfs_moment_zero(&s);gfs_moment_zero(&rhs);
    s.a[GFM_M]=1;s.a[GFM_QXX]=s.a[GFM_QYY]=1./12.;
    s.a[GFM_DXY]=omega/12.;s.a[GFM_DYX]=-omega/12.;
    s.a[GFM_E]=pressure/(gamma-1)+omega*omega/12.;
    GfsMomentPrimitive primitive;ok(gfs_moment_primitive(&s,&primitive),"primitive");
    near(primitive.U,pressure/(gamma-1),"thermal");
    gfs_moment_euler_source(&s,&primitive,gamma,&rhs);
    double abscissa=1./(2*sqrt(3.));
    for(int edge=0;edge<4;edge++)for(int q=0;q<2;q++){
        double t=(q?1:-1)*abscissa,x,y,nx,ny;
        if(edge==0){x=.5;y=t;nx=1;ny=0;}
        else if(edge==1){x=-.5;y=t;nx=-1;ny=0;}
        else if(edge==2){x=t;y=.5;nx=0;ny=1;}
        else{x=t;y=-.5;nx=0;ny=-1;}
        double un=(-omega*y)*nx+(omega*x)*ny,flux[4];
        ok(gfs_material_pressure_flux(pressure,nx,ny,un,un,1e-14,flux),"material face");
        gfs_moment_face_rhs(&rhs,x,y,flux,.5);
    }
    near(rhs.a[GFM_M],0,"mass");near(rhs.a[GFM_PX],0,"px");near(rhs.a[GFM_PY],0,"py");
    near(rhs.a[GFM_E],0,"pressure energy");near(gfs_moment_Lz(&rhs,0,0),0,"angular momentum");
    near(rhs.a[GFM_QXX],0,"Qxx");near(rhs.a[GFM_QXY],0,"Qxy");near(rhs.a[GFM_QYY],0,"Qyy");
    near(rhs.a[GFM_DXX],omega*omega/12.,"Dxx");near(rhs.a[GFM_DYY],omega*omega/12.,"Dyy");
    near(rhs.a[GFM_DXY],0,"Dxy");near(rhs.a[GFM_DYX],0,"Dyx");
}
int main(void){
    check_square(0);check_square(2);
    double flux[4]={7,8,9,10},before[4];memcpy(before,flux,sizeof(flux));
    ok(!gfs_material_pressure_flux(1,1,0,.2,.1,1e-6,flux),"nonmaterial face accepted");
    ok(!memcmp(before,flux,sizeof(flux)),"rejection modified flux");
    ok(!gfs_material_pressure_flux(1,2,0,.2,.2,1e-6,flux),"invalid normal accepted");
    ok(gfs_material_pressure_flux(2,1,0,.3,.3,0,flux),"valid face rejected");
    GfsMomentState a,b;gfs_moment_zero(&a);gfs_moment_zero(&b);
    gfs_moment_face_rhs(&a,.5,.2,flux,.7);gfs_moment_face_rhs(&b,-.5,.2,flux,-.7);
    near(a.a[GFM_PX]+b.a[GFM_PX],0,"pair px");near(a.a[GFM_E]+b.a[GFM_E],0,"pair energy");
    near(gfs_moment_Lz(&a,0,0)+gfs_moment_Lz(&b,1,0),0,"pair angular momentum");
    puts("uniform pressure, rigid rotation, material-face and pair balance checks pass");
}
