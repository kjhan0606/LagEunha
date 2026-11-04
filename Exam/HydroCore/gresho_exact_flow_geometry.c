/* GEOMETRY-ONLY control. Move birth generators by the analytic stationary
 * Gresho flow. Mass, entropy inputs, and weights are not numerically evolved.
 * The native solver is subsequently used only to measure actual Voronoi V.
 * This file is not an accepted hydrodynamic checkpoint or restart authority. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "eunha.h"
#include "gfs_weight_checkpoint.h"
int main(int argc,char **argv){
    if(argc!=4)return 2;char *end;double t=strtod(argv[2],&end);
    if(*end||!(t>=0)||!isfinite(t))return 2;
    FILE *in=fopen(argv[1],"rb");GfsWeightSnapshotHeader header;
    treevorostressrk4particletype *p=NULL;unsigned char footer[52];char error[GFS_WEIGHT_ERROR_BYTES]={0};
    if(!in||gfs_weight_snapshot_read(in,&header,&p,footer,error)||fclose(in))return 2;
    double initial_time;memcpy(&initial_time,footer,sizeof(initial_time));if(initial_time!=0)return 2;
    double max_radial_error=0;
    for(uint64_t i=0;i<header.count;i++){
        if(p[i].weight.W!=0||p[i].weight.pi!=0)return 2;
        double x=p[i].x-.5,y=p[i].y-.5,r=hypot(x,y),speed=r<=.2?5*r:r<=.4?2-5*r:0;
        double omega=r>0?speed/r:0,angle=omega*t,c=cos(angle),s=sin(angle);
        double dx=c*x-s*y,dy=s*x+c*y;
        p[i].x=.5+dx;p[i].y=.5+dy;p[i].vx=-omega*dy;p[i].vy=omega*dx;
        max_radial_error=fmax(max_radial_error,fabs(hypot(dx,dy)-r));
    }
    memcpy(footer,&t,sizeof(t));FILE *out=fopen(argv[3],"w+bx");
    if(!out||gfs_weight_snapshot_write(out,p,header.count,header.features,footer,error)||fclose(out))return 2;
    printf("{\"scope\":\"analytic incompressible flow geometry only, no PDE evolution; thermal input left at birth\",\"cells\":%llu,\"mapping_time\":%.17g,\"max_radial_mapping_error\":%.17g}\n",(unsigned long long)header.count,t,max_radial_error);
    free(p);return 0;
}
