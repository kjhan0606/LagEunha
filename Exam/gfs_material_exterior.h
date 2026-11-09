#ifndef GFS_MATERIAL_EXTERIOR_H
#define GFS_MATERIAL_EXTERIOR_H
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include "gfs_material_power_face_motion.h"
#ifdef __CUDACC__
#define GFS_EXTERIOR_HD static __host__ __device__ inline
#else
#define GFS_EXTERIOR_HD static inline
#endif
/* Power-face volume derivative, valid in dimensions 1..3. d and c are
 * neighbor position and face centroid relative to owner i. c is the AREA
 * centroid in 3D. W has units length^2. No mass transfer is introduced.
 * Vdot_i = base + coefficient*(Wdot_i-Wdot_j).
 * This is geometry, not a prescription for a fluid velocity or pressure. */
GFS_EXTERIOR_HD int gfs_material_power_face_rate(int dim,const double*d,
 const double*c,const double*vi,const double*vj,double area,
 double*base,double*coefficient){
 if(!base||!coefficient)return 0;
 GfsMaterialPowerFaceMotion face;
 if(!gfs_material_power_face_motion(dim,d,c,area,&face)||
    !gfs_material_power_face_motion_rate(dim,&face,vi,vj,0,0,base))return 0;
 *coefficient=face.weight;return 1;
}
/* Undirected graph: i is owned; j<n is another owned cell, j>=n a
 * prescribed exterior site. The latter's Wdot contributes to RHS, not to
 * the unknown vector. An anchored graph must NOT have its mean removed. */
typedef struct {int i,j;double coefficient;} GfsMaterialExteriorEdge;
GFS_EXTERIOR_HD double gfs_material_exterior_edge_action(
 double coefficient,double owner_value,double neighbor_value){
 return coefficient*(owner_value-neighbor_value);
}
#undef GFS_EXTERIOR_HD

static inline void gfs_material_exterior_matvec(int n,int nf,
 const GfsMaterialExteriorEdge*edge,const double*x,double*y){
 memset(y,0,(size_t)n*sizeof(double));
 for(int k=0;k<nf;k++){
  int i=edge[k].i,j=edge[k].j;
  double a=gfs_material_exterior_edge_action(edge[k].coefficient,x[i],j<n?x[j]:0.);
  y[i]+=a;if(j<n)y[j]-=a;
 }
}
static inline int gfs_material_exterior_root(int*parent,int i){
 while(parent[i]!=i){parent[i]=parent[parent[i]];i=parent[i];}return i;
}
/* Same graph certificate for CPU and CUDA. Each connected component must
 * have a prescribed exterior anchor; no diagonal regularization is added. */
static inline int gfs_material_exterior_diagonal(int n,int nt,int nf,
 const GfsMaterialExteriorEdge*edge,double*diag){
 if(n<1||nt<=n||nf<1||!edge||!diag)return 0;
 int ok=0,*parent=(int*)malloc((size_t)n*sizeof(int));
 unsigned char*anchor=(unsigned char*)calloc((size_t)n,1);
 if(!parent||!anchor)goto done;
 for(int i=0;i<n;i++){parent[i]=i;diag[i]=0;}
 for(int k=0;k<nf;k++){
  int i=edge[k].i,j=edge[k].j;double c=edge[k].coefficient;
  if(i<0||i>=n||j<0||j>=nt||i==j||!(c>0)||!isfinite(c))goto done;
  diag[i]+=c;
  if(j<n){diag[j]+=c;int a=gfs_material_exterior_root(parent,i),b=gfs_material_exterior_root(parent,j);parent[a]=b;}
 }
 for(int k=0;k<nf;k++)if(edge[k].j>=n)anchor[gfs_material_exterior_root(parent,edge[k].i)]=1;
 for(int i=0;i<n;i++)if(!anchor[gfs_material_exterior_root(parent,i)]||!(diag[i]>0)||!isfinite(diag[i]))goto done;
 ok=1;
 done:free(parent);free(anchor);return ok;
}
/* Jacobi PCG with the same 1e-10 true-residual criterion as the existing
 * closed-domain solver. Every connected component must have an exterior
 * anchor. Invalid/singular graphs fail; no regularization or mean repair.
 * rhs already contains known-exterior contributions +c*Wdot_external.
 * Shared by native open-boundary evolution and the CPU reference audits. */
static inline int gfs_material_exterior_solve(int n,int nt,int nf,
 const GfsMaterialExteriorEdge*edge,const double*rhs,double*q,
 int*iterations,double*relative){
 if(n<1||nt<=n||nf<1||!edge||!rhs||!q||!iterations||!relative)return 0;
 int ok=0;
 double*buf=(double*)calloc((size_t)n*5,sizeof(double));
 if(!buf)goto done;
 {
 double*r=buf,*z=r+n,*direction=z+n,*Ap=direction+n,*diag=Ap+n;
 *iterations=0;*relative=INFINITY;
 for(int i=0;i<n;i++){q[i]=0;if(!isfinite(rhs[i]))goto done;}
 if(!gfs_material_exterior_diagonal(n,nt,nf,edge,diag))goto done;
 double norm=0,rz=0;
 for(int i=0;i<n;i++){
  r[i]=rhs[i];z[i]=r[i]/diag[i];direction[i]=z[i];norm+=r[i]*r[i];rz+=r[i]*z[i];
 }
 if(norm==0){*relative=0;ok=1;goto done;}
 if(!isfinite(norm)||!(rz>0)||!isfinite(rz))goto done;
 int limit=n>256?n:256;
 for(int it=0;it<limit;it++){
  gfs_material_exterior_matvec(n,nf,edge,direction,Ap);
  double curvature=0;
  for(int i=0;i<n;i++)curvature+=direction[i]*Ap[i];
  if(!(curvature>0)||!isfinite(curvature))goto done;
  double alpha=rz/curvature,res=0;
  for(int i=0;i<n;i++){q[i]+=alpha*direction[i];r[i]-=alpha*Ap[i];res+=r[i]*r[i];}
  *iterations=it+1;
  if(!isfinite(res))goto done;
  if(sqrt(res/norm)<=1e-10){
   gfs_material_exterior_matvec(n,nf,edge,q,Ap);res=0;
   for(int i=0;i<n;i++){r[i]=rhs[i]-Ap[i];res+=r[i]*r[i];}
   *relative=sqrt(res/norm);
   if(isfinite(*relative)&&*relative<=1e-10){ok=1;goto done;}
   rz=0;for(int i=0;i<n;i++){z[i]=r[i]/diag[i];direction[i]=z[i];rz+=r[i]*z[i];}
   continue;
  }
  double next=0;for(int i=0;i<n;i++){z[i]=r[i]/diag[i];next+=r[i]*z[i];}
  if(!(next>0)||!isfinite(next))goto done;
  double beta=next/rz;rz=next;
  for(int i=0;i<n;i++)direction[i]=z[i]+beta*direction[i];
 }
 gfs_material_exterior_matvec(n,nf,edge,q,Ap);
 double res=0;for(int i=0;i<n;i++){double e=rhs[i]-Ap[i];res+=e*e;}
 *relative=sqrt(res/norm);
 }
 done:free(buf);return ok;
}
/* Backend selection changes only where the identical SPD system is solved.
 * Keep the CPU entry point available as an independent arithmetic reference. */
#ifdef GFS_MATERIAL_CUDA_EXTERIOR
#ifdef __cplusplus
extern "C" {
#endif
int gfs_material_cuda_exterior_pcg(int,int,int,const GfsMaterialExteriorEdge*,
 const double*,double*,int*,double*);
#ifdef __cplusplus
}
#endif
#endif
static inline int gfs_material_exterior_solve_selected(int n,int nt,int nf,
 const GfsMaterialExteriorEdge*edge,const double*rhs,double*q,int*iterations,double*relative){
#ifdef GFS_MATERIAL_CUDA_EXTERIOR
 return gfs_material_cuda_exterior_pcg(n,nt,nf,edge,rhs,q,iterations,relative);
#else
 return gfs_material_exterior_solve(n,nt,nf,edge,rhs,q,iterations,relative);
#endif
}
#endif
