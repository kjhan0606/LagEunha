/* Diagnostic state variant, not a trajectory: replace only generator velocity
 * by the exact stationary Gresho field. Preserve geometry, mass, U and time.
 * The native checkpoint codec validates and rewrites its checksum. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "eunha.h"
#include "gfs_weight_checkpoint.h"
int main(int argc,char **argv){
    if(argc!=3)return 2;
    FILE *in=fopen(argv[1],"rb");GfsWeightSnapshotHeader h;
    treevorostressrk4particletype *p=NULL;unsigned char footer[52];char error[GFS_WEIGHT_ERROR_BYTES]={0};
    if(!in||gfs_weight_snapshot_read(in,&h,&p,footer,error)||fclose(in)){
        fprintf(stderr,"read failed: %s\n",error);return 2;
    }
    double max_change=0;
    for(uint64_t i=0;i<h.count;i++){
        if(p[i].weight.W!=0||p[i].weight.pi!=0)return 2;
        double dx=p[i].x-.5,dy=p[i].y-.5,r=hypot(dx,dy),v=r<=.2?5*r:r<=.4?2-5*r:0;
        double vx=r>0?-v*dy/r:0,vy=r>0?v*dx/r:0;
        max_change=fmax(max_change,hypot(p[i].vx-vx,p[i].vy-vy));p[i].vx=vx;p[i].vy=vy;
    }
    FILE *out=fopen(argv[2],"w+bx");
    if(!out||gfs_weight_snapshot_write(out,p,h.count,h.features,footer,error)||fclose(out)){
        fprintf(stderr,"write failed: %s\n",error);return 2;
    }
    printf("{\"scope\":\"exact generator velocity diagnostic variant; no evolution\",\"cells\":%llu,\"max_velocity_replacement\":%.17g}\n",(unsigned long long)h.count,max_change);
    free(p);return 0;
}
