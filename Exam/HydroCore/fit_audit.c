#include <stdio.h>
#include <math.h>
#include "../gfs_velocity_fit.h"
int main(void){
    int pass=1;double worst=0;
    const double A[3][3]={{.2,1.5,-.1},{-.7,.4,.6},{.8,-.3,-.6}};
    for(int d=1;d<=3;d++){
        GfsVelocityFit f;gfs_velocity_fit_init(&f,d);
        for(int j=1;j<=12;j++){
            double r[3]={0},dv[3]={0};
            for(int a=0;a<d;a++)r[a]=sin(.7*j*(a+1))+.2*cos(.37*j*(a+2));
            for(int a=0;a<d;a++)for(int b=0;b<d;b++)dv[a]+=A[a][b]*r[b];
            if(!gfs_velocity_fit_add(&f,r,dv))return 2;
        }
        double G[3][3],cond;
        if(!gfs_velocity_fit_solve(&f,G,&cond))return 2;
        for(int a=0;a<d;a++)for(int b=0;b<d;b++)worst=fmax(worst,fabs(G[a][b]-A[a][b]));
    }
    GfsVelocityFit f;gfs_velocity_fit_init(&f,3);
    double r[3]={1,0,0},dv[3]={1,0,0},G[3][3];
    for(int i=0;i<4;i++)gfs_velocity_fit_add(&f,r,dv);
    int rejects_rank_deficient=!gfs_velocity_fit_solve(&f,G,NULL);
    pass=worst<1e-13&&rejects_rank_deficient;
    printf("{\"affine_dimensions\":[1,2,3],\"max_gradient_error\":%.17g,\"rejects_rank_deficient\":%s,\"pass\":%s}\n",
        worst,rejects_rank_deficient?"true":"false",pass?"true":"false");
    return pass?0:1;
}
