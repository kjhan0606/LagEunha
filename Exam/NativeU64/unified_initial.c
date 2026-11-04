/* Versioned native IC writer. Uniform generators have exact rectangular
 * periodic Voronoi cells. Discontinuous conserved states are integrated
 * over those cells, including cells cut by an oblique interface.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <mpi.h>
#include "eunha.h"
#include "gfs_weight_checkpoint.h"
typedef struct {double x,y;} Point;
static int clip(Point *p,int n,double nx,double ny,double b){
    Point q[32];int m=0;
    for(int k=0;k<n;k++){
        Point a=p[k],c=p[(k+1)%n];double fa=nx*a.x+ny*a.y-b,fc=nx*c.x+ny*c.y-b;
        if(fa<=0)q[m++]=a;
        if((fa<=0)!=(fc<=0)){double t=fa/(fa-fc);q[m++]=(Point){a.x+t*(c.x-a.x),a.y+t*(c.y-a.y)};}
    }
    memcpy(p,q,m*sizeof(Point));return m;
}
static double stripe_fraction(double x,double y,double h,int diagonal,double L){
    double phase=x/L+(diagonal?y/L:0),area=0,ny=diagonal?1./L:0;
    int a=(int)floor(phase-2*h/L)-1,b=(int)floor(phase+2*h/L)+1;
    for(int k=a;k<=b;k++){
        Point p[32]={{-.5*h,-.5*h},{.5*h,-.5*h},{.5*h,.5*h},{-.5*h,.5*h}};
        int n=clip(p,4,1./L,ny,k+.75-phase);
        n=clip(p,n,-1./L,-ny,phase-k-.25);
        double A=0;for(int i=0;i<n;i++)A+=p[i].x*p[(i+1)%n].y-p[i].y*p[(i+1)%n].x;
        area+=.5*fabs(A);
    }
    return fmin(1,fmax(0,area/(h*h)));
}
int main(int argc,char **argv){
    if(argc!=7){fprintf(stderr,"usage: unified_initial N L contact|moving_contact|sod|shear axis|diagonal gamma output\n");return 2;}
    int n=atoi(argv[1]),diagonal=!strcmp(argv[4],"diagonal");double L=atof(argv[2]);
    const char *problem=argv[3];double gamma=(double)(float)atof(argv[5]);
    if(n<8||n>2048||!(L>0)||!(gamma>1)||
       (strcmp(argv[4],"axis")&&!diagonal))return 2;
    int contact=!strcmp(problem,"contact")||!strcmp(problem,"moving_contact"),sod=!strcmp(problem,"sod"),shear=!strcmp(problem,"shear");
    if(!contact&&!sod&&!shear)return 2;
    if(shear&&diagonal)return 2;
    uint64_t count=(uint64_t)n*n;double h=L/n,V=h*h;
    treevorostressrk4particletype *p=calloc(count,sizeof(*p));if(!p)return 2;
    long double M=0,E=0;double maxPressureError=0;
    for(int j=0;j<n;j++)for(int i=0;i<n;i++){
        uint64_t k=(uint64_t)j*n+i;treevorostressrk4particletype *a=p+k;
        double x=(i+.5)*h,y=(j+.5)*h,f=contact||sod?stripe_fraction(x,y,h,diagonal,L):1;
        double rho=contact||sod?.125+.875*f:1,pressure=sod?.1+.9*f:1;
        double vx=!strcmp(problem,"moving_contact")?1:0,vy=!strcmp(problem,"moving_contact")?.5:0;
        /* The variational method evolves generator velocities. Sampling
         * the analytic shear at the generator avoids reclassifying an
         * unresolved velocity variance as initial thermal energy. */
        if(shear)vx=sin(2*M_PI*y/L);
        a->u4if.indx=k;a->x=x;a->y=y;a->mass=rho*V;a->vx=vx;a->vy=vy;
        a->weight=(GfsWeightState){0,0,a->mass/(h*h),GFS_WEIGHT_STATE_TAG};
        a->weight_thermal=(GfsWeightThermal){pressure*V/(gamma-1),0,V};
        a->volume=(float)V;a->ie=(float)a->weight_thermal.U;
        a->den=(float)rho;a->pressure=(float)((gamma-1)*a->weight_thermal.U/V);
        a->csound=(float)sqrt(gamma*a->pressure/rho);a->avgNeighboringPressure=a->pressure;
        a->stress.K=a->pressure/pow(rho,gamma);
        M+=a->mass;E+=a->weight_thermal.U+.5L*a->mass*(vx*vx+vy*vy);
        maxPressureError=fmax(maxPressureError,fabs((gamma-1)*a->weight_thermal.U/V-pressure));
    }
    unsigned char footer[52]={0};double t=0,dt=0;float dtold=0,alpha=0,beta=0;
    memcpy(footer,&t,8);memcpy(footer+8,&dt,8);memcpy(footer+16,&dtold,4);
    memcpy(footer+20,&n,4);memcpy(footer+24,&n,4);memcpy(footer+28,&L,8);memcpy(footer+36,&L,8);
    memcpy(footer+44,&alpha,4);memcpy(footer+48,&beta,4);
    char error[GFS_WEIGHT_ERROR_BYTES];FILE *out=fopen(argv[6],"w+bx");
    uint64_t features=GFS_WEIGHT_FEATURE_INERTIA|GFS_WEIGHT_FEATURE_FROZEN|GFS_WEIGHT_FEATURE_THERMAL64;
    if(!out||gfs_weight_snapshot_write(out,p,count,features,footer,error)){fprintf(stderr,"IC write failed %s\n",error);return 2;}
    if(fclose(out))return 2;
    printf("{\"problem\":\"%s\",\"Nside\":%d,\"cells\":%llu,\"L\":%.17g,\"gamma\":%.17g,\"diagonal\":%d,\"mass\":%.17Lg,\"total_energy\":%.17Lg,\"max_pressure_projection_shift\":%.17g,\"mass_internal_energy_cell_integrals\":true,\"generator_velocity_samples\":%s}\n",problem,n,(unsigned long long)count,L,gamma,diagonal,M,E,maxPressureError,shear?"true":"false");
    free(p);return 0;
}
