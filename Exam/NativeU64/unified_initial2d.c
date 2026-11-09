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
/* IC-only control: exact polygon overlap of a rotated physical square source
 * with an initial rectangular Voronoi cell. Coordinates are relative to the
 * source center. This does not rotate the mesh or modify hydro operators. */
static double sedov_source_overlap(double x,double y,double hx,double hy,
                                   double width,double angle){
    double p[12][2]={{x-hx/2,y-hy/2},{x+hx/2,y-hy/2},
                     {x+hx/2,y+hy/2},{x-hx/2,y+hy/2}},q[12][2];
    double c=cos(angle),s=sin(angle);
    const double normal[4][2]={{c,s},{-c,-s},{-s,c},{s,-c}};
    int n=4;
    for(int k=0;k<4&&n;k++){
        int m=0;
        for(int i=0;i<n;i++){
            int j=(i+1)%n;
            double a=normal[k][0]*p[i][0]+normal[k][1]*p[i][1]-width/2;
            double b=normal[k][0]*p[j][0]+normal[k][1]*p[j][1]-width/2;
            if(a<=0){if(m>=12)return NAN;q[m][0]=p[i][0];q[m++][1]=p[i][1];}
            if((a<=0)!=(b<=0)){
                double t=a/(a-b);if(m>=12)return NAN;
                q[m][0]=p[i][0]+t*(p[j][0]-p[i][0]);
                q[m++][1]=p[i][1]+t*(p[j][1]-p[i][1]);
            }
        }
        n=m;memcpy(p,q,(size_t)n*sizeof(p[0]));
    }
    if(n<3)return 0;
    double area=0;
    for(int i=1;i+1<n;i++)
        area+=(p[i][0]-p[0][0])*(p[i+1][1]-p[0][1])
             -(p[i][1]-p[0][1])*(p[i+1][0]-p[0][0]);
    return area/2;
}


/* Rectangle/disk intersection via an analytic quadrant primitive. The disk
 * is a source-geometry control only. Long-double primitive differences keep
 * boundary-cell deposition accurate without angular sampling or smoothing. */
static long double sedov_disk_primitive(long double x,long double y,long double r){
    long double sign=(x<0?-1.L:1.L)*(y<0?-1.L:1.L);
    x=fminl(fabsl(x),r);y=fminl(fabsl(y),r);
    if(x*x+y*y<=r*r)return sign*x*y;
    long double a=sqrtl(fmaxl(0.L,r*r-y*y));
    long double fx=.5L*(x*sqrtl(fmaxl(0.L,r*r-x*x))+r*r*asinl(x/r));
    long double fa=.5L*(a*y+r*r*asinl(fminl(1.L,a/r)));
    return sign*(a*y+fx-fa);
}
static double sedov_disk_overlap(double x,double y,double hx,double hy,double radius){
    long double r=radius,x0=(long double)x-hx/2.L,x1=(long double)x+hx/2.L;
    long double y0=(long double)y-hy/2.L,y1=(long double)y+hy/2.L;
    long double nearx=fmaxl(0.L,fabsl(x)-hx/2.L),neary=fmaxl(0.L,fabsl(y)-hy/2.L);
    if(nearx*nearx+neary*neary>=r*r)return 0;
    long double farx=fabsl(x)+hx/2.L,fary=fabsl(y)+hy/2.L;
    if(farx*farx+fary*fary<=r*r)return hx*hy;
    long double area=sedov_disk_primitive(x1,y1,r)-sedov_disk_primitive(x0,y1,r)
                    -sedov_disk_primitive(x1,y0,r)+sedov_disk_primitive(x0,y0,r);
    return (double)area;
}

int main(int argc,char **argv){
    if(argc<8||argc>10){fprintf(stderr,"usage: initial2d Nx Ny Lx Ly gresho|kh|sod_periodic|rt_modal|sedov2d|shear_periodic|noh2d gamma output [frozen|dynamic [seed_amplitude]]\n");return 2;}
    int dynamic=argc>=9&&!strcmp(argv[8],"dynamic");
    if(argc>=9&&!dynamic&&strcmp(argv[8],"frozen"))return 2;
    double kh_seed=!strcmp(argv[5],"rt_modal")?.001:.01;
    if(argc==10){
        char *end=NULL;kh_seed=strtod(argv[9],&end);
        if(end==argv[9]||*end||!isfinite(kh_seed)||kh_seed<0||(strcmp(argv[5],"kh")&&strcmp(argv[5],"rt_modal")))return 2;
    }
    int nx=atoi(argv[1]),ny=atoi(argv[2]);double lx=atof(argv[3]),ly=atof(argv[4]),gamma=(double)(float)atof(argv[6]);
    int gresho=!strcmp(argv[5],"gresho"),kh=!strcmp(argv[5],"kh"),sod=!strcmp(argv[5],"sod_periodic"),rt=!strcmp(argv[5],"rt_modal"),sedov=!strcmp(argv[5],"sedov2d"),shear=!strcmp(argv[5],"shear_periodic"),noh=!strcmp(argv[5],"noh2d");
    if(nx<16||ny<16||nx>2048||ny>4096||!(gamma>1)||(!gresho&&!kh&&!sod&&!rt&&!sedov&&!shear&&!noh)||!(lx>0)||!isfinite(lx)||(!noh&&lx!=1)||!(ly>0)||((gresho||shear)&&ly!=1)||(kh&&ly!=2)||((rt||sedov)&&(ly!=1||nx!=ny||ny%2))||(noh&&(ly!=lx||nx!=ny||nx%2)))return 2;
    /* Optional physical square-source width: IC definition only, not a hydro
     * parameter. Default retains the original four-cell deposition arithmetic. */
    double source_width=2.0/nx;int fixed_source=0;
    if(sedov){const char *env=getenv("GFS_SEDOV_SOURCE_WIDTH");if(env&&*env){
        char *end=NULL;source_width=strtod(env,&end);
        if(end==env||*end||!isfinite(source_width)||!(source_width>0)||source_width>.5)return 2;
        fixed_source=1;
    }}
    /* Optional source orientation in degrees; zero keeps the original IC
     * arithmetic exactly, including default four-cell deposition. */
    double source_angle=0;
    if(sedov){const char *env=getenv("GFS_SEDOV_SOURCE_ANGLE_DEG");if(env&&*env){
        char *end=NULL;double degrees=strtod(env,&end);
        if(end==env||*end||!isfinite(degrees))return 2;
        source_angle=fmod(degrees,90.0)*(M_PI/180.0);
    }}
    /* Equal second radial moment to the square of source_width: <r^2>
     * = width^2/6 = radius^2/2. No new hydro parameter is introduced. */
    int disk_source=0;double source_radius=source_width/sqrt(3.0);
    if(sedov){const char *env=getenv("GFS_SEDOV_SOURCE_SHAPE");if(env&&*env){
        if(!strcmp(env,"disk"))disk_source=1;
        else if(strcmp(env,"square"))return 2;
    }}
    if(disk_source&&(!fixed_source||source_angle!=0||!(M_PI*source_radius*source_radius>0)))return 2;
    _Static_assert(sizeof(treevorostressrk4particletype)==696,"native IC ABI");
    uint64_t count=(uint64_t)nx*ny;double hx=lx/nx,hy=ly/ny,V=hx*hy;
    treevorostressrk4particletype *p=calloc(count,sizeof(*p));if(!p)return 2;
    double *rt_mode=NULL;
    if(rt){
        const char *path=getenv("GFS_RT_MODE_FILE");FILE *mode=path?fopen(path,"r"):NULL;int rows;
        if(!mode||fscanf(mode,"%d",&rows)!=1||rows!=ny){free(p);if(mode)fclose(mode);return 2;}
        rt_mode=malloc((size_t)2*ny*sizeof(double));if(!rt_mode){fclose(mode);free(p);return 2;}
        for(int j=0;j<ny;j++){double y;
            if(fscanf(mode,"%lf%lf%lf",&y,&rt_mode[2*j],&rt_mode[2*j+1])!=3||
               fabs(y-(j+.5)*hy)>1e-12||!isfinite(rt_mode[2*j])||!isfinite(rt_mode[2*j+1])){
                fclose(mode);free(rt_mode);free(p);return 2;
            }
        }fclose(mode);
    }
    long double M=0,E=0;double max_quad_difference=0,max_P_error=0;
    for(int j=0;j<ny;j++)for(int i=0;i<nx;i++){
        uint64_t k=(uint64_t)j*nx+i;treevorostressrk4particletype *a=p+k;
        double x=(i+.5)*hx,y=(j+.5)*hy,rho=1,P,vx,vy;
        if(noh){
            /* Existing NativeU64/util.c Noh convention: rho=1, P=1e-6,
             * radial speed one. Even square grids avoid the central singularity.
             * This writer supplies an IC, not analytic inflow boundary data. */
            double dx=x-.5*lx,dy=y-.5*ly,r=hypot(dx,dy);
            if(!(r>0)){free(p);return 2;}
            rho=1;P=1e-6;vx=-dx/r;vy=-dy/r;
        }else if(gresho){
            double dx=x-.5,dy=y-.5,r=hypot(dx,dy),speed=r<=.2?5*r:r<=.4?2-5*r:0;
            vx=r>0?-speed*dy/r:0;vy=r>0?speed*dx/r:0;P=gresho_average(x,y,hx,hy,1);
            max_quad_difference=fmax(max_quad_difference,fabs(P-gresho_average(x,y,hx,hy,0)));
        }else if(rt){
            /* Flat aligned interface, g_y=-.5, hydrostatic pressure. Linear
             * pressure within each rectangle integrates to its center value.
             * Use the same float-representable gamma as the native driver. */
            rho=y<.5?1:2;P=y<.5?2-.5*y:1.75-(y-.5);
            vx=kh_seed*rt_mode[2*j+1]*sin(2*M_PI*x);
            vy=kh_seed*rt_mode[2*j]*cos(2*M_PI*x);
        }else if(shear){
            /* Exact steady periodic shear: rho=1,P=1e-6,vx=.1 sin(2pi y),vy=0.
             * Generator velocity samples follow the existing IC convention.
             * No added physical viscosity, forcing, or subcell energy model. */
            P=1e-6;vx=.1*sin(2*M_PI*y);vy=0;
        }else if(sedov){
            /* Unit explosion energy in the four central cells, above a
             * uniform P=1e-6 background. Fixed mass, no velocity/weight kick.
             * This finite deposition converges to a point source as h->0;
             * compare only before the blast reaches the periodic boundary. */
            vx=vy=0;P=1e-6;
            if(disk_source){
                double area=sedov_disk_overlap(x-.5,y-.5,hx,hy,source_radius);
                if(!isfinite(area)||area<0||area>V*(1+1e-12)){free(p);return 2;}
                P+=(gamma-1)*area/(M_PI*source_radius*source_radius*V);
            }else if(source_angle!=0){
                double area=sedov_source_overlap(x-.5,y-.5,hx,hy,source_width,source_angle);
                if(!isfinite(area)||area<0||area>V*(1+1e-12)){free(p);return 2;}
                P+=(gamma-1)*area/(source_width*source_width*V);
            }else if(fixed_source){
                double lo=.5-.5*source_width,hi=.5+.5*source_width;
                double ox=fmax(0,fmin(x+.5*hx,hi)-fmax(x-.5*hx,lo));
                double oy=fmax(0,fmin(y+.5*hy,hi)-fmax(y-.5*hy,lo));
                P+=(gamma-1)*ox*oy/(source_width*source_width*V);
            }else if((i==nx/2-1||i==nx/2)&&(j==ny/2-1||j==ny/2))P+=(gamma-1)*.25/V;
        }else if(sod){
            /* Aligned discontinuities: exact initial rectangular cell integrals.
             * Periodicity adds a second Riemann problem at x=0; compare
             * the central wave before the two solutions interact. */
            if(nx%2){free(p);return 2;}
            rho=x<.5?1:.125;P=x<.5?1:.1;vx=vy=0;
        }else{
            rho=1+kh_average(y,hy);P=10;vx=-1+2*kh_profile(y);
            vy=kh_seed*sin(2*M_PI*x)*(exp(-pow((y-.5)/.2,2))+exp(-pow((y-1.5)/.2,2)));
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
    uint64_t features=GFS_WEIGHT_FEATURE_INERTIA|GFS_WEIGHT_FEATURE_THERMAL64|(dynamic?0:GFS_WEIGHT_FEATURE_FROZEN);
    if(!out||gfs_weight_snapshot_write(out,p,count,features,footer,error)||fclose(out)){fprintf(stderr,"IC write failed %s\n",error);return 2;}
    printf("{\"problem\":\"%s\",\"Nx\":%d,\"Ny\":%d,\"cells\":%llu,\"gamma\":%.17g,\"mass\":%.17Lg,\"energy\":%.17Lg,\"max_pressure_projection_roundoff\":%.17g,\"max_pressure_quadrature_8_vs_4\":%.17g,\"generator_velocity_samples\":true}\n",argv[5],nx,ny,(unsigned long long)count,gamma,M,E,max_P_error,max_quad_difference);
    free(rt_mode);free(p);return 0;
}
