#include <stdio.h>
#include <math.h>
#include "../gfs_velocity_fit.h"
int main(void){
    int pass=1;double worst=0;
    const double A[3][3]={{.2,1.5,-.1},{-.7,.4,.6},{.8,-.3,-.6}};
    for(int d=1;d<=3;d++)for(int weighted=0;weighted<=1;weighted++){
        GfsVelocityFit f;gfs_velocity_fit_init(&f,d);
        for(int j=1;j<=12;j++){
            double r[3]={0},dv[3]={0};
            for(int a=0;a<d;a++)r[a]=sin(.7*j*(a+1))+.2*cos(.37*j*(a+2));
            for(int a=0;a<d;a++)for(int b=0;b<d;b++)dv[a]+=A[a][b]*r[b];
            double area=.2+j*.17;
            if(!(weighted?gfs_velocity_fit_add_face(&f,r,dv,area):gfs_velocity_fit_add(&f,r,dv)))return 2;
        }
        double G[3][3],cond;
        if(!gfs_velocity_fit_solve(&f,G,&cond))return 2;
        for(int a=0;a<d;a++)for(int b=0;b<d;b++)worst=fmax(worst,fabs(G[a][b]-A[a][b]));
    }
    GfsVelocityFit f;gfs_velocity_fit_init(&f,3);
    double r[3]={1,0,0},dv[3]={1,0,0},G[3][3];
    for(int i=0;i<4;i++)gfs_velocity_fit_add(&f,r,dv);
    int rejects_rank_deficient=!gfs_velocity_fit_solve(&f,G,NULL);
    double point_limit_error=0,equal_sample_error=0;
    for(int d=2;d<=3;d++){
        GfsVelocityFit geom,uniform;gfs_velocity_fit_init(&geom,d);gfs_velocity_fit_init(&uniform,d);
        /* Nonlinear planar velocity. Cardinal faces preserve transverse
         * symmetry. A diagonal touching only a point must disappear. */
        for(int a=0;a<d;a++)for(int sign=-1;sign<=1;sign+=2){
            double r[3]={0},dv[3]={0};r[a]=sign;dv[0]=r[0]*r[0];
            if(!gfs_velocity_fit_add_face(&geom,r,dv,1)||!gfs_velocity_fit_add(&uniform,r,dv))return 2;
        }
        double r[3]={1,1,0},dv[3]={1,0,0};
        if(!gfs_velocity_fit_add_face(&geom,r,dv,1e-16)||!gfs_velocity_fit_add(&uniform,r,dv))return 2;
        if(!gfs_velocity_fit_solve(&geom,G,NULL))return 2;point_limit_error=fmax(point_limit_error,fabs(G[0][1]));
        if(!gfs_velocity_fit_solve(&uniform,G,NULL))return 2;equal_sample_error=fmax(equal_sample_error,fabs(G[0][1]));
    }
    pass=worst<1e-13&&rejects_rank_deficient&&point_limit_error<1e-14&&equal_sample_error>.1;
    printf("{\"affine_dimensions\":[1,2,3],\"max_gradient_error\":%.17g,\"rejects_rank_deficient\":%s,\"vanishing_face_transverse_error\":%.17g,\"equal_sample_control_error\":%.17g,\"pass\":%s}\n",
        worst,rejects_rank_deficient?"true":"false",point_limit_error,equal_sample_error,pass?"true":"false");
    return pass?0:1;
}
