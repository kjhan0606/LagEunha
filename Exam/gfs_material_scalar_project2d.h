#ifndef GFS_MATERIAL_SCALAR_PROJECT2D_H
#define GFS_MATERIAL_SCALAR_PROJECT2D_H
/* Existing qualified pressure-offset / neighbor-bounded scalar projection.
 * Moved unchanged from native_chart_original_csr/shared_fit.h so periodic
 * and prescribed-exterior candidates use one implementation. */
#include "gfs_material_gradient_project.h"
#include <stdio.h>
#ifdef GFS_MATERIAL_PRESSURE_ACCELERATION
#include "gfs_material_pressure_offset.h"
#endif
static inline
#ifdef __CUDACC__
__host__ __device__
#endif
int gfs_project_one_scalar(int i,int scalar,const int *off,
 const GfsTraceDirected2D *edge,const GfsTraceDirected2D *fit,
 const GfsFvPrimitive *c,const unsigned char *walls,double *A,double *b,double *grad){

  double qi=scalar?c[i].p:c[i].rho,lo=qi,hi=qi,M[3][3]={{0}},v[3]={0},g0[3]={0};
  if(!(qi>0)||!isfinite(qi))return 0;
  for(int k=off[i];k<off[i+1];k++){
   int j=edge[k].j;double qj=scalar?c[j].p:c[j].rho;lo=fmin(lo,qj);hi=fmax(hi,qj);
   if(walls&&walls[k])continue;
   double x[2]={fit[k].dx,fit[k].dy},d2=x[0]*x[0]+x[1]*x[1];if(!(d2>0))return 0;
   double w=gfs_trace_metric2d(fit+k,d2)/d2;
#ifdef GFS_MATERIAL_PRESSURE_ACCELERATION
   if(scalar){double ri[2]={fit[k].cx,fit[k].cy},dij[2]={fit[k].dx,fit[k].dy};
    if(!gfs_material_pressure_offset(2,c[i].rho,c[j].rho,ri,dij,x))return 0;}
#endif
   for(int a=0;a<2;a++){v[a]+=w*x[a]*(qj-qi);for(int z=0;z<2;z++)M[a][z]+=w*x[a]*x[z];}
  }
  /* Scale the normal matrix before the generic small solve. */
  double ms=fmax(M[0][0],M[1][1]);if(!(ms>0))return 0;
  for(int a=0;a<2;a++){v[a]/=ms;for(int z=0;z<2;z++)M[a][z]/=ms;}
  if(!mg_solve(2,M,v,g0))return 0;
  int nc=0,it=0;
  for(int k=off[i];k<off[i+1];k++){
   double x=fit[k].cx,y=fit[k].cy,rr=hypot(x,y);if(!(rr>0))return 0;
   /* Always retain positivity at the physical evaluation point. */
   A[2*nc]=-x/rr;A[2*nc+1]=-y/rr;b[nc++]=qi/rr;
   if(!(walls&&walls[k])){
    double dx=fit[k].dx,dy=fit[k].dy;
#ifdef GFS_MATERIAL_PRESSURE_ACCELERATION
    if(scalar){double ri[2]={x,y},dij[2]={dx,dy},effective[2];if(!gfs_material_pressure_offset(2,c[i].rho,c[edge[k].j].rho,ri,dij,effective))return 0;dx=effective[0];dy=effective[1];}
#endif
    double dr=hypot(dx,dy);if(!(dr>0))return 0;
    A[2*nc]=dx/dr;A[2*nc+1]=dy/dr;b[nc++]=(hi-qi)/dr;
    A[2*nc]=-dx/dr;A[2*nc+1]=-dy/dr;b[nc++]=(qi-lo)/dr;
   }
  }
  if(!mg_project(2,nc,A,b,g0,grad+4*i+2*scalar,&it)){
#ifndef __CUDA_ARCH__
   fprintf(stderr,"[GRADIENT_PROJECT_FAIL] cell=%d scalar=%d constraints=%d iterations=%d\n",i,scalar,nc,it);
   FILE *fp=fopen("gradient_project_failure.csv","wx");if(fp){fprintf(fp,"gx,gy\n%.17g,%.17g\nax,ay,b\n",g0[0],g0[1]);for(int j=0;j<nc;j++)fprintf(fp,"%.17g,%.17g,%.17g\n",A[2*j],A[2*j+1],b[j]);fclose(fp);}
#endif
   return 0;
  }
 
 /* Keep a single gradient strictly inside the positive face domain.
    Margin tied to the existing projection solver's 8*512 eps tolerance. */
 double theta=1,*g=grad+4*i+2*scalar;
 for(int k=off[i];k<off[i+1];k++){
  double dx=g[0]*fit[k].cx,dy=g[1]*fit[k].cy,delta=dx+dy;
  double margin=4096*DBL_EPSILON*(qi+fabs(dx)+fabs(dy));
  if(qi+delta<=margin&&delta<0)theta=fmin(theta,nextafter((qi-margin)/(-delta),0.));
 }
 if(!(theta>0)||!isfinite(theta))return 0;
 if(theta<1){g[0]*=theta;g[1]*=theta;}

 return 1;
}
#endif
