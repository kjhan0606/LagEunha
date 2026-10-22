/* Native instantaneous Euler modal/angular contract; no evolution. */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "gfs_dg_euler.h"
#include "gfs_ale_flux.h"
static double angular(const GfsDgBasis &b,double ox,double oy,
    const double rhs[4][GFS_DG_MAX_BASIS]){
    double x=ox+b.gx+b.l00*b.mean[1],y=oy+b.gy+b.l10*b.mean[1]+b.l11*b.mean[2];
    return x*rhs[2][0]+b.l00*rhs[2][1]-y*rhs[1][0]-b.l10*rhs[1][1]-b.l11*rhs[1][2];
}
int main(){
    GfsDgBasis a={},b={};a.area=b.area=1;
    a.gx=.02;a.gy=-.03;a.l00=.1;a.l10=.03;a.l11=.08;a.mean[1]=.1;a.mean[2]=-.2;
    b.gx=-.03;b.gy=.05;b.l00=.08;b.l10=-.02;b.l11=.12;b.mean[1]=-.05;b.mean[2]=.15;
    const double ax=1.2,ay=-.4,bx=-.1,by=.3,gamma=5./3.;
    double ra[4][GFS_DG_MAX_BASIS]={},rb[4][GFS_DG_MAX_BASIS]={};
    double nodes[5]={-.906179845938664,-.538469310105683,0,.538469310105683,.906179845938664};
    double weights[5]={.236926885056189,.478628670499366,.568888888888889,.478628670499366,.236926885056189};
    for(int k=0;k<5;k++){
        double t=nodes[k],x=.7+.2*t,y=.8+.3*t;
        GfsAlePrimitive l={1.2+.1*t,1.e-6*(1+.1*t),.3-.02*t,.1+.1*t};
        GfsAlePrimitive r={.7+.02*t,2.e-6*(1-.05*t),-.1+.07*t,.2-.1*t};
        GfsAleFlux f=gfs_ale_flux(l,r,.6,.8,.03*t,-.02*t,gamma);if(!f.valid)return 1;
        gfs_dg_face_point(&a,2,x-ax,y-ay,f.flux,weights[k],ra);
        gfs_dg_face_point(&b,2,x-bx,y-by,f.flux,-weights[k],rb);
    }
    double pair=fabs(angular(a,ax,ay,ra)+angular(b,bx,by,rb));
    double conserved=0;for(int c=0;c<4;c++)conserved=fmax(conserved,fabs(ra[c][0]+rb[c][0]));
    double orbital_only=(ax*ra[2][0]-ay*ra[1][0])+(bx*rb[2][0]-by*rb[1][0]);
    double volume[4][GFS_DG_MAX_BASIS]={};double u[4]={1.2,.3,-.2,1.5+.5*(.3*.3+.2*.2)/1.2};
    if(!gfs_dg_volume_point(&a,2,u,gamma,.04,-.01,0,0,.2,volume))return 1;
    double volume_L=fabs(angular(a,ax,ay,volume));
    if(pair>1.e-14 || conserved>1.e-14 || volume_L>1.e-14 || fabs(orbital_only)<1.e-5)return 1;
    printf("{\"status\":\"pass\",\"canonical_pair_total_angular_residual\":%.17g,\"mean_Euler_pair_residual\":%.17g,\"orbital_only_pair_torque\":%.17g,\"symmetric_Euler_volume_angular_residual\":%.17g,\"scope\":\"instantaneous modal face/source helpers, not time-integrator angular conservation\"}\n",pair,conserved,orbital_only,volume_L);
}
