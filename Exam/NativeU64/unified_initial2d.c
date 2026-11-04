/* Smooth test ICs for actual rectangular Voronoi cells. Mass and thermal
 * energy are cell integrals; velocities are generator samples. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <mpi.h>
#include "eunha.h"
#include "gfs_weight_checkpoint.h"
static double log_cosh(double x){x=fabs(x);return x+log1p(exp(-2*x))-log(2.);}
static double kh_profile(double y){return .5*(tanh((y-.5)/.05)-tanh((y-1.5)/.05));}
static double kh_average(double y,double hy){
    double a=y-.5*hy,b=y+.5*hy;
    return .025*(log_cosh((b-.5)/.05)-log_cosh((a-.5)/.05)
                   -log_cosh((b-1.5)/.05)+log_cosh((a-1.5)/.05))/hy;
}
static double gresho_pressure(double x,double y){
    double r=hypot(x-.5,y-.5);
    if(r<=.2)return 5+12.5*r*r;
    if(r<=.4)return 9+12.5*r*r-20*r+4*log(5*r);
    return 3+4*log(2.);
}
static double gresho_average(double x,double y,double hx,double hy,int fine){
    const double nodes8[]={-.9602898564975363,-.7966664774136267,-.5255324099163290,-.1834346424956498,.1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363};
    const double weights8[]={.1012285362903763,.2223810344533745,.3137066458778873,.3626837833783620,.3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763};
    const double nodes4[]={-.8611363115940526,-.3399810435848563,.3399810435848563,.8611363115940526};
    const double weights4[]={.3478548451374539,.6521451548625461,.6521451548625461,.3478548451374539};
    const double *q=fine?nodes8:nodes4,*w=fine?weights8:weights4;int n=fine?8:4;double value=0;
    for(int j=0;j<n;j++)for(int i=0;i<n;i++)value+=.25*w[i]*w[j]*gresho_pressure(x+.5*hx*q[i],y+.5*hy*q[j]);
    return value;
}
int main(int argc,char **argv){
    if(argc!=8){fprintf(stderr,"usage: initial2d Nx Ny Lx Ly gresho|kh gamma output\n");return 2;}
    int nx=atoi(argv[1]),ny=atoi(argv[2]);double lx=atof(argv[3]),ly=atof(argv[4]),gamma=(double)(float)atof(argv[6]);
    int gresho=!strcmp(argv[5],"gresho"),kh=!strcmp(argv[5],"kh");
    if(nx<16||ny<16||nx>2048||ny>4096||!(gamma>1)||(!gresho&&!kh)||lx!=1||(gresho?ly!=1:ly!=2))return 2;
    _Static_assert(sizeof(treevorostressrk4particletype)==696,"native IC ABI");
    uint64_t count=(uint64_t)nx*ny;double hx=lx/nx,hy=ly/ny,V=hx*hy;
    treevorostressrk4particletype *p=calloc(count,sizeof(*p));if(!p)return 2;
    long double M=0,E=0;double max_quad_difference=0,max_P_error=0;
    for(int j=0;j<ny;j++)for(int i=0;i<nx;i++){
        uint64_t k=(uint64_t)j*nx+i;treevorostressrk4particletype *a=p+k;
        double x=(i+.5)*hx,y=(j+.5)*hy,rho=1,P,vx,vy;
        if(gresho){
            double dx=x-.5,dy=y-.5,r=hypot(dx,dy),speed=r<=.2?5*r:r<=.4?2-5*r:0;
            vx=r>0?-speed*dy/r:0;vy=r>0?speed*dx/r:0;P=gresho_average(x,y,hx,hy,1);
            max_quad_difference=fmax(max_quad_difference,fabs(P-gresho_average(x,y,hx,hy,0)));
        }else{
            rho=1+kh_average(y,hy);P=10;vx=-1+2*kh_profile(y);
            vy=.01*sin(2*M_PI*x)*(exp(-pow((y-.5)/.2,2))+exp(-pow((y-1.5)/.2,2)));
        }
        a->u4if.indx=k;a->x=x;a->y=y;a->mass=rho*V;a->vx=vx;a->vy=vy;
        a->weight=(GfsWeightState){0,0,rho,GFS_WEIGHT_STATE_TAG};
        a->weight_thermal=(GfsWeightThermal){P*V/(gamma-1),0,V};
        a->volume=(float)V;a->ie=(float)a->weight_thermal.U;a->den=(float)rho;a->pressure=(float)P;
        a->csound=(float)sqrt(gamma*P/rho);a->avgNeighboringPressure=a->pressure;a->stress.K=P/pow(rho,gamma);
        M+=a->mass;E+=a->weight_thermal.U+.5L*a->mass*(vx*vx+vy*vy);
        max_P_error=fmax(max_P_error,fabs((gamma-1)*a->weight_thermal.U/V-P));
    }
    unsigned char footer[52]={0};double t=0,dt=0;float dtold=0,alpha=0,beta=0;
    memcpy(footer,&t,8);memcpy(footer+8,&dt,8);memcpy(footer+16,&dtold,4);
    memcpy(footer+20,&nx,4);memcpy(footer+24,&ny,4);memcpy(footer+28,&lx,8);memcpy(footer+36,&ly,8);
    memcpy(footer+44,&alpha,4);memcpy(footer+48,&beta,4);
    char error[GFS_WEIGHT_ERROR_BYTES]={0};FILE *out=fopen(argv[7],"w+bx");
    uint64_t features=GFS_WEIGHT_FEATURE_INERTIA|GFS_WEIGHT_FEATURE_FROZEN|GFS_WEIGHT_FEATURE_THERMAL64;
    if(!out||gfs_weight_snapshot_write(out,p,count,features,footer,error)||fclose(out)){fprintf(stderr,"IC write failed %s\n",error);return 2;}
    printf("{\"problem\":\"%s\",\"Nx\":%d,\"Ny\":%d,\"cells\":%llu,\"gamma\":%.17g,\"mass\":%.17Lg,\"energy\":%.17Lg,\"max_pressure_projection_roundoff\":%.17g,\"max_pressure_quadrature_8_vs_4\":%.17g,\"generator_velocity_samples\":true}\n",argv[5],nx,ny,(unsigned long long)count,gamma,M,E,max_P_error,max_quad_difference);
    free(p);return 0;
}
