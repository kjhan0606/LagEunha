/* Uniform cubic Voronoi ICs. Discontinuous mass and internal energy are
 * actual cell integrals. Smooth velocities are generator samples. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../Tests3D/common/lag3d_ic.h"
static double sum_cdf(double t,int d){
    if(t<=0)return 0;if(t>=d)return 1;
    double f=0;int choose=1;
    for(int k=0;k<=d;k++){
        if(k)choose=choose*(d-k+1)/k;
        double u=fmax(0,t-k);f+=(k%2?-1:1)*choose*pow(u,d);
    }
    if(d==3)f/=6;return f;
}
static double fraction(double phase,int n,int d){
    double f=0;
    for(int k=(int)floor(phase)-2;k<=(int)floor(phase)+2;k++)
        f+=sum_cdf((k+.75-phase)*n+.5*d,d)-sum_cdf((k+.25-phase)*n+.5*d,d);
    return fmin(1,fmax(0,f));
}
int main(int argc,char **argv){
    if(argc!=7){fprintf(stderr,"usage: initial3d N L contact|moving_contact|sod|shear axis|diagonal gamma output\n");return 2;}
    int n=atoi(argv[1]),diagonal=!strcmp(argv[4],"diagonal");double L=atof(argv[2]),gamma=atof(argv[5]);
    const char *problem=argv[3];int contact=!strcmp(problem,"contact")||!strcmp(problem,"moving_contact"),sod=!strcmp(problem,"sod"),shear=!strcmp(problem,"shear");
    if(n<4||n>256||!(L>0)||!(gamma>1)||(!contact&&!sod&&!shear)||(strcmp(argv[4],"axis")&&!diagonal))return 2;
    FILE *reserved=fopen(argv[6],"wx");if(!reserved)return 2;fclose(reserved);
    lag3d_data d={0};long count=(long)n*n*n;d.h=(lag3d_header){.version=1,.flags=LAG3D_FLAG_VOL,.np=count,.time=0,.gamma=gamma,.box={0,L,0,L,0,L},.periodic={1,1,1}};
    double **arrays[]={&d.x,&d.y,&d.z,&d.vx,&d.vy,&d.vz,&d.mass,&d.u,&d.rho,&d.vol};
    for(int a=0;a<10;a++){*arrays[a]=calloc(count,sizeof(double));if(!*arrays[a])return 2;}
    d.id=calloc(count,sizeof(*d.id));if(!d.id)return 2;
    double h=L/n,V=h*h*h;long double M=0,U=0;
    for(int k=0;k<n;k++)for(int j=0;j<n;j++)for(int i=0;i<n;i++){
        long at=((long)k*n+j)*n+i;double x=(i+.5)*h,y=(j+.5)*h,z=(k+.5)*h;
        double phase=(x+(diagonal?y+z:0))/L,f=contact||sod?fraction(phase,n,diagonal?3:1):1;
        double rho=contact||sod?.125+.875*f:1,P=sod?.1+.9*f:1;
        d.x[at]=x;d.y[at]=y;d.z[at]=z;d.id[at]=at;d.mass[at]=rho*V;d.rho[at]=rho;d.vol[at]=V;d.u[at]=P/((gamma-1)*rho);
        if(!strcmp(problem,"moving_contact")){d.vx[at]=1;d.vy[at]=.5;d.vz[at]=.25;}
        if(shear)d.vx[at]=sin(2*M_PI*(y+(diagonal?z:0))/L);
        M+=d.mass[at];U+=d.mass[at]*d.u[at];
    }
    int rc=lag3d_write(argv[6],&d);
    printf("{\"problem\":\"%s\",\"Nside\":%d,\"cells\":%ld,\"diagonal\":%d,\"mass\":%.17Lg,\"internal_energy\":%.17Lg,\"rc\":%d}\n",problem,n,count,diagonal,M,U,rc);
    lag3d_free(&d);return rc?2:0;
}
