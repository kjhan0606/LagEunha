/* Diagnostic only: same strict shared clip, no state or geometry commit. */
#include "power_polygon_clip.h"
#include <stdio.h>
#include <stdint.h>
int main(int argc,char **argv){
 if(argc!=2)return 1;FILE *f=fopen(argv[1],"rb");if(!f)return 2;
 int64_t nr,owner;if(fread(&nr,8,1,f)!=1||fread(&owner,8,1,f)!=1)return 3;
 double center[3];if(fread(center,8,3,f)!=3)return 3;
 double x[128]={-3-center[0],3-center[0],3-center[0],-3-center[0]},y[128]={-3-center[1],-3-center[1],3-center[1],3-center[1]},tx[128],ty[128];
 int64_t tag[128]={-1,-2,-3,-4},tt[128];int nc=4;
 for(int k=0;k<nr;k++){
  int64_t j;double a[3];if(fread(&j,8,1,f)!=1||fread(a,8,3,f)!=3)return 3;
  double dx=a[0]-center[0],dy=a[1]-center[1],threshold=.5*(dx*dx+dy*dy+center[2]-a[2]);
  printf("cut=%lld owner=%lld delta=(%.17g,%.17g) threshold=%.17g\n",(long long)j,(long long)owner,dx,dy,threshold);
  for(int e=0;e<nc;e++){
   int en=(e+1)%nc;double da=x[e]*dx+y[e]*dy,db=x[en]*dx+y[en]*dy;
   if((da<=threshold)!=(db<=threshold))printf(" crossing e=%d tag=%lld A=(%.17g,%.17g) B=(%.17g,%.17g) da=%.17g db=%.17g denominator=%.17g t=%.17g\n",e,(long long)tag[e],x[e],y[e],x[en],y[en],da,db,db-da,(threshold-da)/(db-da));
  }
  int nn=0,err=power_polygon_clip(nc,x,y,tag,dx,dy,threshold,j,128,tx,ty,tt,1,0.,&nn);if(err)return err+10;
  nc=nn;for(int e=0;e<nc;e++){x[e]=tx[e];y[e]=ty[e];tag[e]=tt[e];printf(" vertex=%d xy=(%.17g,%.17g) outgoingID=%lld\n",e,x[e],y[e],(long long)tag[e]);}
 }
 fclose(f);return 0;
}
