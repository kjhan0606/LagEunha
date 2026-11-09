#ifndef GFS_MATERIAL_EXTERIOR_ALIAS2D_H
#define GFS_MATERIAL_EXTERIOR_ALIAS2D_H
#include "gfs_material_exterior2d.h"
/* Resolve only certified periodic exterior samples to their real owners.
 * Geometric displacement/centroid remains that of the actual image face.
 * Caller certifies ID, translated position, equal weight and scalar state.
 * The downstream common CSR builder certifies reciprocal interior faces.
 * Other exterior sites retain their physical boundary-data slots.
 * Input is the unchanged geometry/pair list from exterior2d_build (or a
 * previous alias call). Remapping indices does not change A/(2 distance).
 * Preserve the validated coefficients and canonical summation order. */
static inline int gfs_material_exterior_alias2d(GfsMaterialExterior2D*g,const int*alias){
 if(!g||!alias||g->n<1||g->nt<g->n||g->directed_count<1||g->pair_count<1||
    !g->directed||!g->pair||!g->pair_directed)return 0;
 for(int j=g->n;j<g->nt;j++)if(alias[j]<-1||alias[j]>=g->n)return 0;
 int pairs=0;
 for(int k=0;k<g->directed_count;k++){
  GfsTraceDirected2D*e=g->directed+k;
#ifdef GFS_MATERIAL_FIXED_SOLID
  if(e->i==e->j){
   if(e->i<0||e->i>=g->n||pairs>=g->pair_count||g->pair_directed[pairs]!=k||g->pair[pairs].i!=e->i||g->pair[pairs].j!=e->j||g->pair[pairs].coefficient!=0)return 0;
   pairs++;continue;
  }
#endif
  if(e->i<0||e->i>=g->n||e->j<0||e->j>=g->nt||e->i==e->j)return 0;
  if(e->j>=g->n&&alias[e->j]>=0&&alias[e->j]==e->i)return 0;
  if(e->i<e->j){
   if(pairs>=g->pair_count||g->pair_directed[pairs]!=k||
      g->pair[pairs].i!=e->i||g->pair[pairs].j!=e->j||
      !(g->pair[pairs].coefficient>0)||!isfinite(g->pair[pairs].coefficient))return 0;
   pairs++;
  }
 }
 if(pairs!=g->pair_count)return 0;
 for(int k=0;k<g->directed_count;k++){
  GfsTraceDirected2D*e=g->directed+k;
  if(e->j>=g->n&&alias[e->j]>=0)e->j=alias[e->j];
 }
 int count=0;
 for(int k=0;k<pairs;k++){
  int z=g->pair_directed[k];const GfsTraceDirected2D*e=g->directed+z;
  if(e->i<e->j
#ifdef GFS_MATERIAL_FIXED_SOLID
     || e->i==e->j
#endif
  ){
   double coefficient=g->pair[k].coefficient;
   g->pair[count]=(GfsMaterialExteriorEdge){e->i,e->j,coefficient};
   g->pair_directed[count++]=z;
  }
 }
 g->pair_count=count;return count>0;
}
#endif
