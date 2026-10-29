#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Exam/gfs_weight_mass_exchange.h"
static void ok(int b,const char *s){if(!b){fprintf(stderr,"%s\n",s);exit(1);}}
static void near(double a,double b,const char *s){if(!isfinite(a)||fabs(a-b)>2e-12*(1+fabs(b))){fprintf(stderr,"%s %.17g %.17g\n",s,a,b);exit(2);}}
int main(void){
    double ax=-.3,ay=.2,bx=1.9,by=-1.4,dx=bx-ax,dy=by-ay,L=hypot(dx,dy);
    double ox=-.2,oy=.3,ox2=1.7,oy2=-.4,cx=.5*(ax+bx)-ox,cy=.5*(ay+by)-oy;
    double wi=.12,wj=-.04,d=2.,rho=1.5,ux=.3,uy=-.7,e=4.,G[]={.2,-.4,.6,.1};
    double fm=-rho*(wi-wj)/(2*d),dvx=G[0]*dx+G[1]*dy,dvy=G[2]*dx+G[3]*dy;
    GfsMomentState a,b;gfs_moment_zero(&a);gfs_moment_zero(&b);
    ok(gfs_weight_exchange_segment_affine(&a,&b,ax,ay,bx,by,ox,oy,ox2,oy2,wi,wj,d,rho,ux,uy,G,e),"affine exchange");
    near(a.a[GFM_M],-fm*L,"mass");near(a.a[GFM_PX],-fm*L*ux,"px");near(a.a[GFM_PY],-fm*L*uy,"py");
    near(a.a[GFM_QXX],-fm*L*(cx*cx+dx*dx/12),"qxx");
    near(a.a[GFM_QXY],-fm*L*(cx*cy+dx*dy/12),"qxy");
    near(a.a[GFM_QYY],-fm*L*(cy*cy+dy*dy/12),"qyy");
    near(a.a[GFM_DXX],-fm*L*(cx*ux+dx*dvx/12),"dxx");
    near(a.a[GFM_DXY],-fm*L*(cx*uy+dx*dvy/12),"dxy");
    near(a.a[GFM_DYX],-fm*L*(cy*ux+dy*dvx/12),"dyx");
    near(a.a[GFM_DYY],-fm*L*(cy*uy+dy*dvy/12),"dyy");
    near(a.a[GFM_E],-fm*L*(e+.5*(ux*ux+uy*uy)+(dvx*dvx+dvy*dvy)/24),"energy");
    near(a.a[GFM_M]+b.a[GFM_M],0,"mass sum");near(a.a[GFM_PX]+b.a[GFM_PX],0,"px sum");
    near(a.a[GFM_PY]+b.a[GFM_PY],0,"py sum");near(a.a[GFM_E]+b.a[GFM_E],0,"energy sum");
    near(gfs_moment_Lz(&a,ox,oy)+gfs_moment_Lz(&b,ox2,oy2),0,"angular sum");
    double vx=.41,vy=-.2;GfsMomentState shifted_a,shifted_b;gfs_moment_zero(&shifted_a);gfs_moment_zero(&shifted_b);
    ok(gfs_weight_exchange_segment_affine(&shifted_a,&shifted_b,ax,ay,bx,by,ox,oy,ox2,oy2,wi,wj,d,rho,ux+vx,uy+vy,G,e),"boosted exchange");
    GfsMomentState expected=gfs_moment_boost(a,vx,vy);
    for(int k=0;k<GFM_NSTATE;k++)near(shifted_a.a[k],expected.a[k],"boost");
    GfsMomentState before_a=a,before_b=b;G[0]=NAN;
    ok(!gfs_weight_exchange_segment_affine(&a,&b,ax,ay,bx,by,ox,oy,ox2,oy2,wi,wj,d,rho,ux,uy,G,e),"invalid gradient accepted");
    ok(!memcmp(&a,&before_a,sizeof(a))&&!memcmp(&b,&before_b,sizeof(b)),"invalid input changed state");
    G[0]=.2;GfsMomentState za,zb;gfs_moment_zero(&za);gfs_moment_zero(&zb);
    ok(gfs_weight_exchange_segment_affine(&za,&zb,ax,ay,bx,by,ox,oy,ox2,oy2,0,0,d,rho,ux,uy,G,e),"W0 reduction");
    for(int k=0;k<GFM_NSTATE;k++){near(za.a[k],0,"zero rate i");near(zb.a[k],0,"zero rate j");}
    ok(!gfs_weight_exchange_segment_affine(&za,&zb,ax,ay,bx,by,ox,oy,ox2,oy2,wi,wj,d,rho,ux,uy,G,-1),"negative internal energy accepted");
    puts("affine segment mass, moment, energy, boost and transactional checks pass");
}
