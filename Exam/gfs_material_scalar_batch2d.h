#ifndef GFS_MATERIAL_SCALAR_BATCH2D_H
#define GFS_MATERIAL_SCALAR_BATCH2D_H
#include <stdint.h>
#include "gfs_material_scalar_project2d.h"
#ifdef GFS_MATERIAL_CUDA_PROJECT
#include "gfs_material_cuda_project.h"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

typedef int (*GfsMaterialProjectScalar2D)(int,int,const int*,
 const GfsTraceDirected2D*,const GfsTraceDirected2D*,const GfsFvPrimitive*,
 const unsigned char*,double*,double*,double*);

/* Same per-cell scalar projector for periodic and prescribed-exterior data.
 * Only n_owned gradients are written; neighbor values may use n_sample
 * physical samples. Riemann boundary slots are not gradient samples.
 * The callback must write only its cell's gradient and supplied workspace.
 * A failed batch must be discarded before any face trace is published. */
static inline int gfs_material_project_scalar_cpu_batch2d(int n_owned,
 int n_sample,int degree,const int*off,const GfsTraceDirected2D*edge,
 const GfsTraceDirected2D*fit,const GfsFvPrimitive*state,
 const unsigned char*wall,double*grad,GfsMaterialProjectScalar2D project){
 if(n_owned<1||n_sample<n_owned||degree<1||!off||!edge||!fit||!state||
    !grad||!project||off[0]!=0)return 0;
 for(int i=0;i<n_owned;i++){
  if(off[i+1]<off[i]||off[i+1]-off[i]>degree)return 0;
  for(int k=off[i];k<off[i+1];k++)
   if(edge[k].i!=i||edge[k].j<0||edge[k].j>=n_sample)return 0;
 }
 int threads=1,allok=1;
#ifdef _OPENMP
 threads=omp_get_max_threads();
#endif
 if((size_t)threads>SIZE_MAX/((size_t)degree*6*sizeof(double)))return 0;
 double*A=(double*)malloc((size_t)threads*degree*6*sizeof(double));
 double*b=(double*)malloc((size_t)threads*degree*3*sizeof(double));
 if(!A||!b){free(A);free(b);return 0;}
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(&:allok)
#endif
 for(int i=0;i<n_owned;i++){
  int tid=0;
#ifdef _OPENMP
  tid=omp_get_thread_num();
#endif
  double*localA=A+(size_t)tid*degree*6,*localb=b+(size_t)tid*degree*3;
  for(int scalar=0;scalar<2;scalar++)
   allok &= project(i,scalar,off,edge,fit,state,wall,localA,localb,grad);
 }
 free(A);free(b);return allok;
}

/* A CUDA batch uses the same fixed scalar function. Preserve arbitrary CPU
 * callbacks, and never hide a failed CUDA solve by retrying it on the CPU. */
static inline int gfs_material_project_scalar_batch2d(int n_owned,
 int n_sample,int degree,const int*off,const GfsTraceDirected2D*edge,
 const GfsTraceDirected2D*fit,const GfsFvPrimitive*state,
 const unsigned char*wall,double*grad,GfsMaterialProjectScalar2D project){
#ifdef GFS_MATERIAL_CUDA_PROJECT
 if(project==gfs_project_one_scalar){
  if(n_owned<1||n_sample<n_owned||degree<1||!off||!grad||off[0]!=0)return 0;
  for(int i=0;i<n_owned;i++)
   if(off[i+1]<off[i]||off[i+1]-off[i]>degree)return 0;
  if(!gfs_cuda_project_gradients_samples(n_owned,n_sample,off[n_owned],
       off,edge,fit,state,wall,grad))return 0;
  if(getenv("GFS_CUDA_PROJECT_AUDIT")){
   double*reference=calloc((size_t)n_owned*4,sizeof(double));
   if(!reference)return 0;
   int pass=gfs_material_project_scalar_cpu_batch2d(n_owned,n_sample,degree,
     off,edge,fit,state,wall,reference,project);double worst=0;
   if(pass)for(size_t k=0;k<(size_t)4*n_owned;k++){
    double error=fabs(grad[k]-reference[k])/(1+fabs(reference[k]));
    worst=fmax(worst,error);if(!isfinite(grad[k])||error>1e-10)pass=0;
   }
   fprintf(stderr,"[CUDA_PROJECT_AUDIT] cells=%d samples=%d directed_faces=%d max_scaled_error=%.17g pass=%d\n",
     n_owned,n_sample,off[n_owned],worst,pass);
   free(reference);if(!pass)return 0;
  }
  return 1;
 }
#endif
 return gfs_material_project_scalar_cpu_batch2d(n_owned,n_sample,degree,
   off,edge,fit,state,wall,grad,project);
}
#endif
