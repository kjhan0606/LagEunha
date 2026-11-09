#ifndef GFS_MATERIAL_PROJECT_SHARED2D_H
#define GFS_MATERIAL_PROJECT_SHARED2D_H
#if !defined(GFS_MATERIAL_DIRECT_FV) || !defined(GFS_MATERIAL_CENTROID_TRACE) || !defined(GFS_MATERIAL_FACE_NORMAL_LIMIT)
#error "Projected gradients require direct FV centroid traces"
#endif
#include "gfs_material_gradient_project.h"
#ifdef GFS_MATERIAL_PRESSURE_ACCELERATION
#include "gfs_material_pressure_offset.h"
#endif
#include "gfs_material_scalar_batch2d.h"
/* One limited gradient per scalar per cell. Shared prepared thermodynamic
 * traces feed both C and CUDA, preserving a single force/work Riemann state. */
static inline int gfs_material_project_shared2d(GfsTraceCsr2D*t,int nf,
 const GfsTraceDirected2D*edge,const GfsFvPrimitive*c,const unsigned char*walls,const GfsTraceDirected2D*fit){
 if(!t||!edge||!c||nf<1)return 0;
 int n=t->cells,ok=0,*off=calloc((size_t)n+1,sizeof(int)),maxdegree=0;
 double *grad=calloc((size_t)n*4,sizeof(double));

 if(!off||!grad||!fit)goto done;
 for(int k=0;k<nf;k++){if(edge[k].i<0||edge[k].i>=n||edge[k].j<0||edge[k].j>=n||(k&&edge[k].i<edge[k-1].i))goto done;off[edge[k].i+1]++;}
 for(int i=0;i<n;i++){if(off[i+1]>maxdegree)maxdegree=off[i+1];off[i+1]+=off[i];}
 if(!gfs_material_project_scalar_batch2d(n,n,maxdegree,off,edge,fit,c,
    walls,grad,gfs_project_one_scalar))goto done;

 {
 int row=0;
 for(int k=0;k<nf;k++)if(edge[k].i<edge[k].j||(walls&&walls[k])){
  if(row>=t->faces)goto done;GfsTraceFace2D*f=t->face+row++;
  if(f->i!=edge[k].i||f->j!=edge[k].j)goto done;
  for(int side=0;side<2;side++){
   int i=side?f->j:f->i;double x=fit[k].cx-(side?fit[k].dx:0),y=fit[k].cy-(side?fit[k].dy:0);
   if(f->i==f->j){x=fit[k].cx;y=fit[k].cy;}
   for(int scalar=0;scalar<2;scalar++){
    double value=(scalar?c[i].p:c[i].rho)+grad[4*i+2*scalar]*x+grad[4*i+2*scalar+1]*y;
    if(!(value>0)||!isfinite(value))goto done;f->thermo[2*side+scalar]=value;
   }
  }
 }
 if(row!=t->faces)goto done;
 }
 ok=1;
 done:free(off);free(grad);return ok;
}
#endif
