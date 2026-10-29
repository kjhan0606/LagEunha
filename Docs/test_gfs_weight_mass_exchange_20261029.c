#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Exam/gfs_weight_mass_exchange.h"
static void eq(double a,double b){if(!isfinite(a)||fabs(a-b)>1e-12*(1+fabs(b)))exit(1);}
int main(void){
    double f[4],reverse[4],zero[4];
    if(!gfs_weight_exchange_flux(.12,-.04,2.,1.5,.3,-.7,4.,f) ||
       !gfs_weight_exchange_flux(-.04,.12,2.,1.5,.3,-.7,4.,reverse) ||
       !gfs_weight_exchange_flux(0,0,2.,1.5,.3,-.7,4.,zero))return 2;
    eq(f[0],-1.5*.16/4.);for(int k=0;k<4;k++){eq(f[k]+reverse[k],0);eq(zero[k],0);}
    GfsMomentState a,b;gfs_moment_zero(&a);gfs_moment_zero(&b);
    double ox1=-.2,oy1=.3,ox2=1.7,oy2=-.4;
    gfs_weight_exchange_pair(&a,&b,.8,-.6,ox1,oy1,ox2,oy2,2.3,f);
    eq(a.a[GFM_M]+b.a[GFM_M],0);eq(a.a[GFM_PX]+b.a[GFM_PX],0);
    eq(a.a[GFM_PY]+b.a[GFM_PY],0);eq(a.a[GFM_E]+b.a[GFM_E],0);
    eq(gfs_moment_Lz(&a,ox1,oy1)+gfs_moment_Lz(&b,ox2,oy2),0);
    eq(a.a[GFM_QX]+ox1*a.a[GFM_M]+b.a[GFM_QX]+ox2*b.a[GFM_M],0);
    eq(a.a[GFM_QY]+oy1*a.a[GFM_M]+b.a[GFM_QY]+oy2*b.a[GFM_M],0);
    eq(a.a[GFM_QXX]+2*ox1*a.a[GFM_QX]+ox1*ox1*a.a[GFM_M]+
       b.a[GFM_QXX]+2*ox2*b.a[GFM_QX]+ox2*ox2*b.a[GFM_M],0);
    eq(a.a[GFM_QYY]+2*oy1*a.a[GFM_QY]+oy1*oy1*a.a[GFM_M]+
       b.a[GFM_QYY]+2*oy2*b.a[GFM_QY]+oy2*oy2*b.a[GFM_M],0);
    eq(a.a[GFM_QXY]+ox1*a.a[GFM_QY]+oy1*a.a[GFM_QX]+ox1*oy1*a.a[GFM_M]+
       b.a[GFM_QXY]+ox2*b.a[GFM_QY]+oy2*b.a[GFM_QX]+ox2*oy2*b.a[GFM_M],0);
    double vx=.41,vy=-.2,energy=4.+vx*.3+vy*(-.7)+.5*(vx*vx+vy*vy);
    double boosted[4];if(!gfs_weight_exchange_flux(.12,-.04,2.,1.5,.3+vx,-.7+vy,energy,boosted))return 3;
    eq(boosted[0],f[0]);eq(boosted[1],f[1]+vx*f[0]);eq(boosted[2],f[2]+vy*f[0]);
    eq(boosted[3],f[3]+vx*f[1]+vy*f[2]+.5*(vx*vx+vy*vy)*f[0]);
    if(gfs_weight_exchange_flux(0,0,0,1,0,0,1,f) ||
       gfs_weight_exchange_flux(NAN,0,1,1,0,0,1,f) ||
       gfs_weight_exchange_flux(0,0,1,-1,0,0,1,f))return 4;
    puts("weight-linked pair exchange: Wdot=0, reciprocity, conservation, boost and invalid-input checks pass");
    return 0;
}
