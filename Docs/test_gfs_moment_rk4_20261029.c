#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Exam/gfs_moment_rk4.h"
#include "Exam/gfs_weight_mass_exchange.h"
static void ok(int c,const char *s){if(!c){fprintf(stderr,"%s\n",s);exit(1);}}
static void near(double x,double y){ok(isfinite(x)&&fabs(x-y)<1e-12*(1+fabs(y)),"numeric mismatch");}
static GfsMomentState initial(double mass){
    GfsMomentState s;gfs_moment_zero(&s);s.a[GFM_M]=mass;
    s.a[GFM_QXX]=mass/12.;s.a[GFM_QYY]=mass/12.;s.a[GFM_E]=10*mass;return s;
}
static double oscillator(int n){
    GfsMomentState state=initial(1);state.a[GFM_PX]=1;
    for(int step=0;step<n;step++){
        GfsMomentRK4 r={0};ok(gfs_moment_rk4_begin(&r,&state,1./n,1),"begin oscillator");
        for(int k=0;k<4;k++){
            GfsMomentState rhs;gfs_moment_zero(&rhs);
            rhs.a[GFM_PX]=-r.stage.a[GFM_PY];rhs.a[GFM_PY]=r.stage.a[GFM_PX];
            ok(gfs_moment_rk4_add_rhs(&r,&rhs),"oscillator stage");
        }
        ok(gfs_moment_rk4_candidate(&r,&state),"oscillator candidate");
    }
    return hypot(state.a[GFM_PX]-cos(1.),state.a[GFM_PY]-sin(1.));
}
int main(void){
    double e5=oscillator(5),e10=oscillator(10),ratio=e5/e10;
    ok(ratio>15 && ratio<17,"fourth-order convergence");
    GfsMomentState base=initial(1),rhs,candidate;GfsMomentRK4 r={0},before;
    ok(gfs_moment_rk4_begin(&r,&base,.2,1),"begin fixed mass");
    ok(!gfs_moment_rk4_candidate(&r,&candidate),"premature candidate");
    gfs_moment_zero(&rhs);rhs.a[GFM_M]=1;before=r;
    ok(!gfs_moment_rk4_add_rhs(&r,&rhs) && !memcmp(&r,&before,sizeof(r)),"reject mass flux and retain scratch");
    gfs_moment_zero(&rhs);rhs.a[GFM_E]=-1000;before=r;
    ok(!gfs_moment_rk4_add_rhs(&r,&rhs) && !memcmp(&r,&before,sizeof(r)),"reject negative thermal state");
    gfs_moment_zero(&rhs);rhs.a[GFM_PX]=NAN;before=r;
    ok(!gfs_moment_rk4_add_rhs(&r,&rhs) && !memcmp(&r,&before,sizeof(r)),"reject nonfinite derivative");
    near(base.a[GFM_M],1);near(base.a[GFM_E],10);

    /* Two cells, one shared conservative weight-flux contribution at every
     * stage; test both fixed coordinate origins and tentative candidate. */
    GfsMomentState a=initial(5),b=initial(7),fa,fb,outa,outb;
    double flux[4];ok(gfs_weight_exchange_flux(.01,-.02,1,1,.1,-.2,1.5,flux),"flux");
    gfs_moment_zero(&fa);gfs_moment_zero(&fb);
    ok(gfs_weight_exchange_segment(&fa,&fb,.5,-.2,.5,.2,0,0,1,0,flux),"face segment");
    GfsMomentRK4 ra={0},rb={0};
    ok(gfs_moment_rk4_begin(&ra,&a,.01,0) && gfs_moment_rk4_begin(&rb,&b,.01,0),"begin pair");
    for(int k=0;k<4;k++)ok(gfs_moment_rk4_add_rhs(&ra,&fa) && gfs_moment_rk4_add_rhs(&rb,&fb),"pair stage");
    ok(gfs_moment_rk4_candidate(&ra,&outa) && gfs_moment_rk4_candidate(&rb,&outb),"pair candidate");
    near(outa.a[GFM_M]+outb.a[GFM_M],a.a[GFM_M]+b.a[GFM_M]);
    near(outa.a[GFM_PX]+outb.a[GFM_PX],0);near(outa.a[GFM_PY]+outb.a[GFM_PY],0);
    near(outa.a[GFM_E]+outb.a[GFM_E],a.a[GFM_E]+b.a[GFM_E]);
    near(gfs_moment_Lz(&outa,0,0)+gfs_moment_Lz(&outb,1,0),0);
    printf("moment RK4: error ratio %.6f, guarded stages and conservative pair pass\n",ratio);
    return 0;
}
