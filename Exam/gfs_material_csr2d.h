#ifndef GFS_MATERIAL_CSR2D_H
#define GFS_MATERIAL_CSR2D_H
/* Host-side initialization from the production FaceCSR layout. This is a
 * one-time topology import, not per-step retessellation or an evolution mode.
 * Bounded single-rank meshes only; periodic/MPI image ownership is not inferred.
 */
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <string.h>
#include "gfs_material_mesh2d.h"
typedef struct GfsMaterialCsr2D {
 int cells,nodes,entries;
 int *offset,*vertex;
 double (*x)[2];
 double max_merge_distance;
} GfsMaterialCsr2D;
static inline void gfs_material_csr_free(GfsMaterialCsr2D *m){
 free(m->offset);free(m->vertex);free(m->x);memset(m,0,sizeof(*m));
}
typedef struct {int a,b,owner,sign,boundary;} GfsMaterialEdgeKey;
static inline int gfs_material_edge_compare(const void *a,const void *b){
 const GfsMaterialEdgeKey *x=(const GfsMaterialEdgeKey*)a,*y=(const GfsMaterialEdgeKey*)b;
 if(x->a!=y->a)return x->a<y->a?-1:1;
 return (x->b>y->b)-(x->b<y->b);
}
static inline uint64_t gfs_material_bin_hash(int64_t x,int64_t y){
 uint64_t a=(uint64_t)x*UINT64_C(11400714819323198485),b=(uint64_t)y*UINT64_C(14029467366897019727);
 return a^(b+(a<<6)+(a>>2));
}
/* Coordinates c1/c2 are owner-relative, exactly as in FaceCSR. boundary[f]
 * explicitly marks an outer face. Imported loops are oriented CCW. Failure
 * leaves output empty; no incomplete or ambiguous topology is accepted.
 */
static inline int gfs_material_csr_import(int nc,const int *offset,
 const double *qx,const double *qy,const double *c1x,const double *c1y,
 const double *c2x,const double *c2y,const int *boundary,GfsMaterialCsr2D *out){
 GfsMaterialCsr2D m={0};int *head=NULL,*chain=NULL,*end1=NULL,*end2=NULL,*used=NULL;
 int64_t *bx=NULL,*by=NULL;GfsMaterialEdgeKey *keys=NULL;int ok=0;
 if(!out||nc<1||!offset||offset[0]!=0)return 0;
 memset(out,0,sizeof(*out));int nf=offset[nc];if(nf<3*nc||nf>100000000)return 0;
 for(int i=0;i<nc;i++)if(offset[i+1]-offset[i]<3)return 0;
 double xmin=DBL_MAX,ymin=DBL_MAX,xmax=-DBL_MAX,ymax=-DBL_MAX;
 for(int i=0;i<nc;i++)for(int f=offset[i];f<offset[i+1];f++)for(int k=0;k<2;k++){
  double x=qx[i]+(k?c2x[f]:c1x[f]),y=qy[i]+(k?c2y[f]:c1y[f]);
  if(!isfinite(x)||!isfinite(y))return 0;
  xmin=fmin(xmin,x);xmax=fmax(xmax,x);ymin=fmin(ymin,y);ymax=fmax(ymax,y);
 }
 double scale=fmax(xmax-xmin,ymax-ymin),tol=128*DBL_EPSILON*scale;
 if(!(tol>0)||!isfinite(tol))return 0;
 size_t cap=1;while(cap<(size_t)4*nf)cap*=2;
 head=(int*)malloc(cap*sizeof(int));chain=(int*)malloc(2*nf*sizeof(int));
 bx=(int64_t*)malloc(2*nf*sizeof(int64_t));by=(int64_t*)malloc(2*nf*sizeof(int64_t));
 end1=(int*)malloc(nf*sizeof(int));end2=(int*)malloc(nf*sizeof(int));used=(int*)calloc(nf,sizeof(int));
 keys=(GfsMaterialEdgeKey*)malloc(nf*sizeof(*keys));
 m.x=(double(*)[2])malloc(2*nf*sizeof(*m.x));m.offset=(int*)malloc((nc+1)*sizeof(int));m.vertex=(int*)malloc(nf*sizeof(int));
 if(!head||!chain||!bx||!by||!end1||!end2||!used||!keys||!m.x||!m.offset||!m.vertex)goto done;
 for(size_t k=0;k<cap;k++)head[k]=-1;
 m.cells=nc;m.entries=nf;memcpy(m.offset,offset,(nc+1)*sizeof(int));
 for(int i=0;i<nc;i++)for(int f=offset[i];f<offset[i+1];f++)for(int k=0;k<2;k++){
  double x=qx[i]+(k?c2x[f]:c1x[f]),y=qy[i]+(k?c2y[f]:c1y[f]);
  int64_t ix=(int64_t)floor((x-xmin)/tol),iy=(int64_t)floor((y-ymin)/tol);int found=-1;
  for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++){
   size_t h=gfs_material_bin_hash(ix+dx,iy+dy)&(cap-1);
   for(int v=head[h];v>=0;v=chain[v])if(bx[v]==ix+dx&&by[v]==iy+dy){
    double distance=hypot(m.x[v][0]-x,m.x[v][1]-y);
    if(distance<=tol){if(found>=0&&found!=v)goto done;found=v;m.max_merge_distance=fmax(m.max_merge_distance,distance);}
   }
  }
  if(found<0){
   found=m.nodes++;m.x[found][0]=x;m.x[found][1]=y;bx[found]=ix;by[found]=iy;
   size_t h=gfs_material_bin_hash(ix,iy)&(cap-1);chain[found]=head[h];head[h]=found;
  }
  if(k)end2[f]=found;else end1[f]=found;
 }
 for(int i=0;i<nc;i++){
  int first=offset[i],last=offset[i+1],f=first;
  for(int k=first;k<last;k++){
   if(used[f]||end1[f]==end2[f])goto done;
   used[f]=1;m.vertex[k]=end1[f];int next=-1;
   for(int j=first;j<last;j++)if(end1[j]==end2[f]){if(next>=0)goto done;next=j;}
   if(next<0)goto done;
   f=next;
  }
  if(f!=first)goto done;
  double area=gfs_mm_area(last-first,m.vertex+first,m.x);
  if(area<0){for(int a=first,b=last-1;a<b;a++,b--){int t=m.vertex[a];m.vertex[a]=m.vertex[b];m.vertex[b]=t;}}
  else if(!(area>0))goto done;
  for(int k=first;k<last;k++){
   int a=m.vertex[k],b=m.vertex[k+1==last?first:k+1];
   for(int j=first;j<last;j++)if(m.vertex[j]!=a&&m.vertex[j]!=b){
    int q=m.vertex[j];if(gfs_mm_cross(m.x[b][0]-m.x[a][0],m.x[b][1]-m.x[a][1],m.x[q][0]-m.x[a][0],m.x[q][1]-m.x[a][1])<=0)goto done;
   }
   int bf=-1;for(int j=first;j<last;j++)if((end1[j]==a&&end2[j]==b)||(end1[j]==b&&end2[j]==a)){bf=j;break;}
   if(bf<0)goto done;
   keys[k].a=a<b?a:b;keys[k].b=a<b?b:a;keys[k].owner=i;keys[k].sign=a<b?1:-1;keys[k].boundary=boundary[bf]!=0;
  }
 }
 qsort(keys,nf,sizeof(*keys),gfs_material_edge_compare);
 for(int a=0;a<nf;){
  int b=a+1;while(b<nf&&keys[b].a==keys[a].a&&keys[b].b==keys[a].b)b++;
  if(b-a==1){if(!keys[a].boundary)goto done;}
  else if(b-a==2){if(keys[a].boundary||keys[a+1].boundary||keys[a].owner==keys[a+1].owner||keys[a].sign==keys[a+1].sign)goto done;}
  else goto done;
  a=b;
 }
 *out=m;memset(&m,0,sizeof(m));ok=1;
 done:
 free(head);free(chain);free(bx);free(by);free(end1);free(end2);free(used);free(keys);gfs_material_csr_free(&m);return ok;
}
#endif
