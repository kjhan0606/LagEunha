#ifndef GFS_MATERIAL_POSITIVE_FAN2D_H
#define GFS_MATERIAL_POSITIVE_FAN2D_H
#include "gfs_material_trace_csr2d.h"
#include <stdio.h>
/* Research-only periodic 2D continuous piecewise-affine reversible trace.
 * Positive barycentric coefficients on owner/adjacent-angle-neighbor triangles.
 * Failure is explicit; there is no extrapolation or negative-weight fallback. */
typedef struct {int k;double angle;} FanAngle;
static int fan_cmp(const void*a,const void*b){double x=((const FanAngle*)a)->angle,y=((const FanAngle*)b)->angle;return (x>y)-(x<y);}
static int gfs_material_positive_fan2d(int n,int nd,const GfsTraceDirected2D*edge,GfsTraceCsr2D*t){
 GfsTraceDirected2D*fit=gfs_trace_centroid_coordinates2d(n,nd,edge,NULL);if(!fit)return 0;
 int*off=calloc(n+1,sizeof(int));FanAngle*order=malloc(nd*sizeof(*order));double *fan_weights=calloc(t->entries,sizeof(double));
 if(!off||!order||!fan_weights){free(off);free(order);free(fan_weights);free(fit);return 0;}
 for(int k=0;k<nd;k++)off[edge[k].i+1]++;for(int k=0;k<n;k++)off[k+1]+=off[k];
 for(int o=0;o<n;o++){for(int k=off[o];k<off[o+1];k++)order[k]=(FanAngle){k,atan2(fit[k].dy,fit[k].dx)};qsort(order+off[o],off[o+1]-off[o],sizeof(*order),fan_cmp);}
 int row=0;
 for(int k=0;k<nd;k++)if(edge[k].i<edge[k].j){
  const GfsTraceFace2D*f=t->face+row++;
  for(int side=0;side<2;side++){
   int o=side?f->j:f->i,start=side?f->sr:f->sl,found=0;double best_slack=-HUGE_VAL;
   double x=fit[k].cx-(side?fit[k].dx:0),y=fit[k].cy-(side?fit[k].dy:0);
   for(int z=off[o];z<off[o+1]&&!found;z++){
    int a=order[z].k,b=order[z+1==off[o+1]?off[o]:z+1].k;
    double ax=fit[a].dx,ay=fit[a].dy,bx=fit[b].dx,by=fit[b].dy,det=ax*by-ay*bx;
    if(det<=1e-12*hypot(ax,ay)*hypot(bx,by))continue;
    double u=(x*by-y*bx)/det,v=(ax*y-ay*x)/det,w=1-u-v;
    best_slack=fmax(best_slack,fmin(w,fmin(u,v)));
    if(u>=0&&v>=0&&w>=0){fan_weights[start]=w;fan_weights[start+1+a-off[o]]=u;fan_weights[start+1+b-off[o]]=v;found=1;}
   }
   if(!found){fprintf(stderr,"[MATERIAL_POSITIVE_FAN_FAIL] face=%d side=%d owner=%d target=%.17g,%.17g best_slack=%.17g neighbors=%d\n",row-1,side,o,x,y,best_slack,off[o+1]-off[o]);free(off);free(order);free(fit);free(fan_weights);return 0;}
  }
 }
 for(int j=0;j<t->entries;j++)t->entry[j].reversible_weight=fan_weights[j];
 free(off);free(order);free(fit);free(fan_weights);return 1;
}

#endif
