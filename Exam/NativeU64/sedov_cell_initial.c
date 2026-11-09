/* Attach certified actual-cell geometry/source integrals to a pristine native
 * Sedov IC. Hydro and the existing regular-grid IC writer stay unchanged. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "eunha.h"
#include "gfs_weight_checkpoint.h"
int main(int argc,char **argv){
 if(argc!=5)return 2;
 char *end=NULL;double gamma=(double)(float)strtod(argv[4],&end);
 if(end==argv[4]||*end||!isfinite(gamma)||gamma<=1)return 2;
 FILE *f=fopen(argv[1],"rb");unsigned char foot[52];char err[256]={0};
 GfsWeightSnapshotHeader h;treevorostressrk4particletype *p=NULL;
 if(!f||gfs_weight_snapshot_read(f,&h,&p,foot,err))return 2;fclose(f);
 double t,dt,lx,ly;int nx,ny;
 memcpy(&t,foot,8);memcpy(&dt,foot+8,8);memcpy(&nx,foot+20,4);memcpy(&ny,foot+24,4);
 memcpy(&lx,foot+28,8);memcpy(&ly,foot+36,8);
 if(t!=0||dt!=0||!h.count||lx!=1||ly!=1||nx<=0||ny<=0||(uint64_t)nx*ny!=h.count){free(p);return 2;}
 f=fopen(argv[2],"r");unsigned long long count=0;
 if(!f||fscanf(f,"%llu",&count)!=1||count!=h.count){if(f)fclose(f);free(p);return 2;}
 long double volume=0,excess=0;
 for(uint64_t k=0;k<h.count;k++){
  double x,y,V,E;
  if(fscanf(f,"%lf%lf%lf%lf",&x,&y,&V,&E)!=4||!isfinite(x)||!isfinite(y)||
     !isfinite(V)||!isfinite(E)||x<0||x>=1||y<0||y>=1||!(V>0)||E<0||
     p[k].u4if.indx!=k||p[k].vx!=0||p[k].vy!=0||p[k].weight.W!=0||p[k].weight.pi!=0){
   fclose(f);free(p);return 2;
  }
  volume+=V;excess+=E;double U=1e-6*V/(gamma-1)+E,P=(gamma-1)*U/V;
  p[k].x=x;p[k].y=y;p[k].mass=V;
  p[k].weight_thermal=(GfsWeightThermal){U,0,V};
  p[k].volume=(float)V;p[k].ie=(float)U;p[k].den=1;p[k].pressure=(float)P;
  p[k].csound=(float)sqrt(gamma*P);p[k].avgNeighboringPressure=p[k].pressure;
  p[k].stress.K=P;p[k].stress.dK=0;
 }
 char extra;int trailing=fscanf(f," %c",&extra);fclose(f);
 if(trailing!=EOF||fabsl(volume-1)>1e-11L||fabsl(excess-1)>1e-11L){free(p);return 2;}
 f=fopen(argv[3],"w+bx");if(!f){free(p);return 2;}
 int rc=gfs_weight_snapshot_write(f,p,h.count,h.features,foot,err);
 if(fclose(f))rc=2;free(p);
 if(!rc)printf("{\"cells\":%llu,\"volume\":%.17Lg,\"excess_energy\":%.17Lg}\n",count,volume,excess);
 return rc;
}
