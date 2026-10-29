/* Independent point-sum checks of affine moment volume and face balances. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Exam/gfs_moment.h"

static void eq(double a,double b,const char *what){
    if(!isfinite(a)||!isfinite(b)||fabs(a-b)>2e-12*(1+fabs(a)+fabs(b))){
        fprintf(stderr,"%s: %.17g != %.17g\n",what,a,b);exit(1);
    }
}
int main(void){
    /* Unequal positive quadrature masses and non-isotropic coordinates avoid
     * relying on a square or special symmetry. Velocities form affine shear
     * plus rotation and translation. Each point carries thermal energy. */
    const double x[]={-.7,.9,.4,-.2,.1},y[]={-.3,-.4,.8,.5,-.1};
    const double dm[]={.4,.7,.6,.8,.5};
    double vx[5],vy[5],p=0,q=0,xx=0,xy=0,yy=0,uu=0,uv=0,vv=0;
    GfsMomentState s,rhs;gfs_moment_zero(&s);gfs_moment_zero(&rhs);
    for(int i=0;i<5;i++){
        vx[i]=1.1+.31*x[i]-.74*y[i];
        vy[i]=-.2+.27*x[i]-.16*y[i];
        gfs_moment_accumulate(&s,dm[i],x[i],y[i],vx[i],vy[i],.2*dm[i]);
        p+=dm[i]*vx[i];q+=dm[i]*vy[i];
        xx+=2*dm[i]*x[i]*vx[i];
        xy+=dm[i]*(x[i]*vy[i]+y[i]*vx[i]);
        yy+=2*dm[i]*y[i]*vy[i];
        uu+=dm[i]*vx[i]*vx[i];uv+=dm[i]*vx[i]*vy[i];vv+=dm[i]*vy[i]*vy[i];
    }
    GfsMomentPrimitive primitive;
    if(!gfs_moment_primitive(&s,&primitive))return 2;
    eq(primitive.B[0],.31,"grad xx");eq(primitive.B[1],-.74,"grad xy");
    eq(primitive.B[2],.27,"grad yx");eq(primitive.B[3],-.16,"grad yy");
    eq(primitive.U,.6,"thermal energy");
    const double gamma=5./3.;gfs_moment_euler_source(&s,&primitive,gamma,&rhs);
    eq(rhs.a[GFM_QX],p,"first x");eq(rhs.a[GFM_QY],q,"first y");
    eq(rhs.a[GFM_QXX],xx,"second xx");eq(rhs.a[GFM_QXY],xy,"second xy");eq(rhs.a[GFM_QYY],yy,"second yy");
    eq(rhs.a[GFM_DXX],uu+(gamma-1)*primitive.U,"momentum moment xx");
    eq(rhs.a[GFM_DXY],uv,"momentum moment xy");eq(rhs.a[GFM_DYX],uv,"momentum moment yx");
    eq(rhs.a[GFM_DYY],vv+(gamma-1)*primitive.U,"momentum moment yy");
    eq(gfs_moment_Lz(&rhs,0,0),0,"isotropic pressure internal torque");

    /* General nonzero mass flux at one common physical face point. Opposite
     * signs, distinct cell reference origins, exactly one shared Euler flux. */
    const double flux[]={.4,-.7,.9,1.3},X=.8,Y=-.6,ox1=-.2,oy1=.3,ox2=1.7,oy2=-.4;
    GfsMomentState f1,f2;gfs_moment_zero(&f1);gfs_moment_zero(&f2);
    gfs_moment_face_rhs(&f1,X-ox1,Y-oy1,flux,2.3);
    gfs_moment_face_rhs(&f2,X-ox2,Y-oy2,flux,-2.3);
    eq(f1.a[GFM_M]+f2.a[GFM_M],0,"face mass");
    eq(f1.a[GFM_PX]+f2.a[GFM_PX],0,"face px");
    eq(f1.a[GFM_PY]+f2.a[GFM_PY],0,"face py");
    eq(f1.a[GFM_E]+f2.a[GFM_E],0,"face energy");
    eq(gfs_moment_Lz(&f1,ox1,oy1)+gfs_moment_Lz(&f2,ox2,oy2),0,"face angular momentum");
    eq(f1.a[GFM_QX]+ox1*f1.a[GFM_M]+f2.a[GFM_QX]+ox2*f2.a[GFM_M],0,"face first moment x");
    eq(f1.a[GFM_QY]+oy1*f1.a[GFM_M]+f2.a[GFM_QY]+oy2*f2.a[GFM_M],0,"face first moment y");

    /* Central gravity, including a nonzero origin, has zero total torque. */
    GfsMomentState grav;gfs_moment_zero(&grav);
    gfs_moment_gravity_source(&grav,.7,.8-ox1,-.6-oy1,.1,.2,-.8,.6);
    eq(gfs_moment_Lz(&grav,ox1,oy1),0,"central gravity torque");
    eq(grav.a[GFM_E],.7*(.1*(-.8)+.2*.6),"gravity work");
    puts("affine volume, shared ALE face, and central-gravity moment balances pass");
    return 0;
}
