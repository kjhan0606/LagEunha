#ifndef GFS_VOLUME_ANCHOR_GRAPH_H
#define GFS_VOLUME_ANCHOR_GRAPH_H
#include <math.h>
#include <stdlib.h>
#include <limits.h>
/* J x = sum c_ij (x_i-x_j) + a_i*x_i, c_ij,a_i >= 0.
 * Every row must have a directed positive-coefficient path to an anchor.
 * Reverse BFS checks this sufficient nonsingularity condition without
 * altering coefficients, dropping small edges or modifying the RHS.
 * Self faces (including stationary walls) contribute exactly zero.
 * This is an anchored-domain certificate, not a periodic gauge solver.
 */
static inline int gfs_volume_anchor_graph(int n,const int*off,const int*owner,
 const double*c,int*reached){
 if(reached)*reached=0;
 if(n<1||!off||!owner||!c||off[0]!=0||off[n]<0||off[n]>INT_MAX-n)return -1;
 int nf=off[n],answer=-1,tail=0;
 int *reverse_off=calloc((size_t)n+1,sizeof(int));
 int *cursor=malloc((size_t)n*sizeof(int)),*queue=malloc((size_t)n*sizeof(int));
 int *reverse=malloc((size_t)(nf?nf:1)*sizeof(int));
 unsigned char *seen=calloc((size_t)n,1);
 if(!reverse_off||!cursor||!queue||!reverse||!seen)goto done;
 for(int i=0;i<n;i++){
  if(off[i+1]<off[i]||off[i+1]>nf)goto done;
  for(int k=off[i];k<off[i+1];k++){
   int j=owner[k];double a=c[k];
   if(j < -1||j>=n||!(a>=0)||!isfinite(a))goto done;
   if(a==0||j==i)continue;
   if(j==-1){if(!seen[i]){seen[i]=1;queue[tail++]=i;}}
   else reverse_off[j+1]++;
  }
 }
 for(int i=0;i<n;i++){reverse_off[i+1]+=reverse_off[i];cursor[i]=reverse_off[i];}
 for(int i=0;i<n;i++)for(int k=off[i];k<off[i+1];k++){
  int j=owner[k];if(j>=0&&j!=i&&c[k]>0)reverse[cursor[j]++]=i;
 }
 for(int head=0;head<tail;head++){
  int j=queue[head];
  for(int k=reverse_off[j];k<reverse_off[j+1];k++){
   int i=reverse[k];if(!seen[i]){seen[i]=1;queue[tail++]=i;}
  }
 }
 if(reached)*reached=tail;
 answer=tail==n;
done:
 free(seen);free(reverse);free(queue);free(cursor);free(reverse_off);return answer;
}
#endif
