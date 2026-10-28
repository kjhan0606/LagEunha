#ifndef EXAM_POWER_POLYGON_CLIP_H
#define EXAM_POWER_POLYGON_CLIP_H
#include <math.h>
#ifdef __CUDACC__
#define PPC_HD __host__ __device__
#else
#define PPC_HD
#endif
/* Production Sutherland-Hodgman edge convention: tag[e] owns e -> e+1.
 * strict=0 preserves legacy endpoint clamp and capacity truncation.
 * strict=1 has neither; caller discards the entire polygon on failure. */
template<class Tag> PPC_HD int power_polygon_clip(int nc,const double *x,const double *y,
 const Tag *tag,double nx,double ny,double threshold,Tag neighbor,int capacity,
 double *ox,double *oy,Tag *ot,int strict,double fallback_fraction,int *count){
 int n=0;
 for(int e=0;e<nc;e++){
  int en=(e+1)%nc;double ax=x[e],ay=y[e],bx=x[en],by=y[en];
  double da=ax*nx+ay*ny,db=bx*nx+by*ny;int ai=da<=threshold,bi=db<=threshold;
  if(ai){
   if(n>=capacity){if(strict)return 1;}
   else{ox[n]=ax;oy[n]=ay;ot[n]=tag[e];n++;}
  }
  if(ai!=bi){
   double den=db-da,t;
   if(strict){if(!isfinite(den)||den==0.)return 2;t=(threshold-da)/den;if(!isfinite(t)||t<0.||t>1.)return 2;}
   else{t=fabs(den)>1e-30?(threshold-da)/den:fallback_fraction;t=fmax(1e-9,fmin(1.-1e-9,t));}
   if(n>=capacity){if(strict)return 1;}
   else{ox[n]=ax+t*(bx-ax);oy[n]=ay+t*(by-ay);ot[n]=ai?neighbor:tag[e];n++;}
  }
 }
 if(strict&&n){
  /* Exact duplicate points only. Keep the outgoing edge of the last copy.
   * No distance threshold, coordinate movement, or finite face deletion. */
  int k=0;
  for(int e=0;e<n;e++){
   if(k&&ox[e]==ox[k-1]&&oy[e]==oy[k-1])ot[k-1]=ot[e];
   else{ox[k]=ox[e];oy[k]=oy[e];ot[k]=ot[e];k++;}
  }
  n=k;if(n>1&&ox[n-1]==ox[0]&&oy[n-1]==oy[0])n--;
 }
 *count=n;return 0;
}
#endif
