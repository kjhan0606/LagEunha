/* Frozen periodic velocity fields, actual project Voronoi geometry.
 * Not a PDE time integrator. No gravity or Kepler-specific construction.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../indxflag.h"
#include "../../Voro/voro.h"
#include "../gfs_variational_core.h"
#include "../gfs_velocity_fit.h"
typedef struct {int sx,sy;double r[3],c[3],A;} LocalFace;
static LocalFace faces[32];static int nf;
static double template(double h){
    Voro2D_point in[8]={0},guard[16]={0},ctr={0},*sorted=guard+8;
    Voro2D_Corner corners[128];int count=0;ctr.indx=-1;
    for(int y=-1;y<=1;y++)for(int x=-1;x<=1;x++)if(x||y){
        in[count].x=x*h;in[count].y=y*h;in[count].indx=count;count++;
    }
    int nv=Voro2D_FindVC(&ctr,in,sorted,count,corners,16,4);
    if(nv<3||nv>32)exit(2);
    double V=0;nf=0;Voro2D_Corner *a=corners;
    for(int k=0;k<nv;k++){
        Voro2D_Corner *b=a->upperlink;int j=a->upperrelated;if(j<0)exit(2);
        LocalFace *f=faces+nf++;f->r[0]=sorted[j].x;f->r[1]=sorted[j].y;f->r[2]=0;
        f->c[0]=.5*(a->x+b->x);f->c[1]=.5*(a->y+b->y);f->c[2]=0;
        f->A=hypot(b->x-a->x,b->y-a->y);f->sx=(int)llround(f->r[0]/h);f->sy=(int)llround(f->r[1]/h);
        V+=.5*f->A*(f->c[0]*f->r[0]+f->c[1]*f->r[1])/hypot(f->r[0],f->r[1]);a=b;
    }
    if(fabs(V/(h*h)-1)>1e-12)exit(2);return V;
}
static int run(int n,int checkerboard,int exact){
    int N=n*n;double h=1./n,V=template(h),gamma=5./3.,cs=sqrt(gamma*1e-6);
    double *v=calloc(2*N,sizeof(double)),*G=calloc(4*N,sizeof(double));if(!v||!G)return 2;
    for(int y=0;y<n;y++)for(int x=0;x<n;x++){
        int i=y*n+x;double a=2*M_PI*(x+.5)*h,b=4*M_PI*(y+.5)*h;
        if(checkerboard){v[2*i]=x%2?1:-1;v[2*i+1]=0;}
        else{v[2*i]=sin(a)*cos(b);v[2*i+1]=-.5*cos(a)*sin(b);}
    }
    double max_condition=0,maxdiv=0;
    for(int y=0;y<n;y++)for(int x=0;x<n;x++){
        int i=y*n+x;GfsVelocityFit fit;gfs_velocity_fit_init(&fit,2);double g[3][3],condition;
        for(int k=0;k<nf;k++){
            LocalFace *f=faces+k;int j=((y+f->sy+n)%n)*n+(x+f->sx+n)%n;
            double dv[3]={v[2*j]-v[2*i],v[2*j+1]-v[2*i+1],0};
            if(!gfs_velocity_fit_add_face(&fit,f->r,dv,f->A))return 2;
        }
        if(!gfs_velocity_fit_solve(&fit,g,&condition))return 2;
        max_condition=fmax(max_condition,condition);
        G[4*i]=g[0][0];G[4*i+1]=g[0][1];G[4*i+2]=g[1][0];G[4*i+3]=g[1][1];
        if(exact&&!checkerboard){
            double a=2*M_PI*(x+.5)*h,b=4*M_PI*(y+.5)*h;
            G[4*i]=2*M_PI*cos(a)*cos(b);G[4*i+1]=-4*M_PI*sin(a)*sin(b);
            G[4*i+2]=M_PI*sin(a)*sin(b);G[4*i+3]=-G[4*i];
        }
        maxdiv=fmax(maxdiv,fabs(G[4*i]+G[4*i+3]));
    }
    long double H=0,K=0;
    for(int y=0;y<n;y++)for(int x=0;x<n;x++){
        int i=y*n+x;K+=.5L*V*(v[2*i]*v[2*i]+v[2*i+1]*v[2*i+1]);
        for(int k=0;k<nf;k++){
            LocalFace *f=faces+k;int j=((y+f->sy+n)%n)*n+(x+f->sx+n)%n;
            double vi[3]={v[2*i],v[2*i+1],0},vj[3]={v[2*j],v[2*j+1],0};
            double jump[3]={vj[0]-vi[0]-.5*((G[4*i]+G[4*j])*f->r[0]+(G[4*i+1]+G[4*j+1])*f->r[1]),
                            vj[1]-vi[1]-.5*((G[4*i+2]+G[4*j+2])*f->r[0]+(G[4*i+3]+G[4*j+3])*f->r[1]),0};
            GfsVariationalFace face;if(!gfs_variational_face(2,f->r,f->c,f->A,&face))return 2;
            double sf=gfs_variational_compression_fraction(G[4*i]+G[4*i+3],G[4*i+2]-G[4*i+1],G[4*j]+G[4*j+3],G[4*j+2]-G[4*j+1]);
            double q=gfs_variational_pair_pressure(&face,vi,vj,jump,1,cs,1,cs,sf,h,h);
            double closing=(vi[0]-vj[0])*face.normal[0]+(vi[1]-vj[1])*face.normal[1];
            H+=.5*q*f->A*fmax(0,closing);
        }
    }
    printf("{\"scope\":\"frozen periodic field, actual project Voro geometry, shared production closure\",\"dim\":2,\"Nside\":%d,\"field\":\"%s\",\"gradient\":\"%s\",\"faces_per_cell\":%d,\"max_condition_bound\":%.17g,\"max_abs_fit_divergence\":%.17g,\"kinetic_energy\":%.17Lg,\"heat_rate\":%.17Lg,\"heat_rate_over_kinetic\":%.17Lg,\"crossing_time_normalized_heat\":%.17Lg}\n",n,checkerboard?"grid_checkerboard":"fixed_wavelength_divergence_free_streamfunction",exact?"analytic_oracle":"common_LSQ",nf,max_condition,maxdiv,K,H,H/K,H*h/(2*K));
    free(v);free(G);return 0;
}
int main(void){int failed=0;for(int n=16;n<=64;n*=2){failed|=run(n,0,0);failed|=run(n,0,1);failed|=run(n,1,0);}return failed;}
