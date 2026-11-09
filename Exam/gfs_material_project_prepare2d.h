#ifndef GFS_MATERIAL_PROJECT_PREPARE2D_H
#define GFS_MATERIAL_PROJECT_PREPARE2D_H
#if !defined(GFS_MATERIAL_DIRECT_FV) || !defined(GFS_MATERIAL_CENTROID_TRACE) || !defined(GFS_MATERIAL_FACE_NORMAL_LIMIT)
#error "Projected gradients require direct FV centroid traces"
#endif
#include "gfs_material_gradient_project.h"
#ifdef GFS_MATERIAL_PRESSURE_ACCELERATION
#include "gfs_material_pressure_offset.h"
#endif
/* One limited gradient per scalar per cell. Shared prepared thermodynamic
 * traces feed both C and CUDA, preserving a single force/work Riemann state. */
static inline int gfs_material_project_prepare2d(GfsTraceCsr2D*t,int nf,
 const GfsTraceDirected2D*edge,const GfsFvPrimitive*c,const unsigned char*walls){
 if(!t||!edge||!c||nf<1)return 0;
 int n=t->cells,ok=0,*off=calloc((size_t)n+1,sizeof(int)),maxdegree=0;
 double *grad=calloc((size_t)n*4,sizeof(double)),*A=NULL,*b=NULL;
 GfsTraceDirected2D*fit=gfs_trace_centroid_coordinates2d(n,nf,edge,walls);
 if(!off||!grad||!fit)goto done;
 for(int k=0;k<nf;k++){if(edge[k].i<0||edge[k].i>=n||edge[k].j<0||edge[k].j>=n||(k&&edge[k].i<edge[k-1].i))goto done;off[edge[k].i+1]++;}
 for(int i=0;i<n;i++){if(off[i+1]>maxdegree)maxdegree=off[i+1];off[i+1]+=off[i];}
 A=malloc((size_t)maxdegree*4*sizeof(double));b=malloc((size_t)maxdegree*2*sizeof(double));if(!A||!b)goto done;
 for(int i=0;i<n;i++)for(int scalar=0;scalar<2;scalar++){
  double qi=scalar?c[i].p:c[i].rho,lo=qi,hi=qi,M[3][3]={{0}},v[3]={0},g0[3]={0};
  if(!(qi>0)||!isfinite(qi))goto done;
  for(int k=off[i];k<off[i+1];k++){
   int j=edge[k].j;double qj=scalar?c[j].p:c[j].rho;lo=fmin(lo,qj);hi=fmax(hi,qj);
   if(walls&&walls[k])continue;
   double x[2]={fit[k].dx,fit[k].dy},d2=x[0]*x[0]+x[1]*x[1];if(!(d2>0))goto done;
   double w=gfs_trace_metric2d(fit+k,d2)/d2;
#ifdef GFS_MATERIAL_PRESSURE_ACCELERATION
   if(scalar){double ri[2]={fit[k].cx,fit[k].cy},dij[2]={fit[k].dx,fit[k].dy};
    if(!gfs_material_pressure_offset(2,c[i].rho,c[j].rho,ri,dij,x))goto done;}
#endif
   for(int a=0;a<2;a++){v[a]+=w*x[a]*(qj-qi);for(int z=0;z<2;z++)M[a][z]+=w*x[a]*x[z];}
  }
  /* Scale the normal matrix before the generic small solve. */
  double ms=fmax(M[0][0],M[1][1]);if(!(ms>0))goto done;
  for(int a=0;a<2;a++){v[a]/=ms;for(int z=0;z<2;z++)M[a][z]/=ms;}
  if(!mg_solve(2,M,v,g0))goto done;
  int nc=0,it=0;
  for(int k=off[i];k<off[i+1];k++){
   double x=fit[k].cx,y=fit[k].cy,rr=hypot(x,y);if(!(rr>0))goto done;
   if(!(walls&&walls[k])){A[2*nc]=x/rr;A[2*nc+1]=y/rr;b[nc++]=(hi-qi)/rr;}
   A[2*nc]=-x/rr;A[2*nc+1]=-y/rr;b[nc++]=(qi-((walls&&walls[k])?0:lo))/rr;
  }
  if(!mg_project(2,nc,A,b,g0,grad+4*i+2*scalar,&it)){
   fprintf(stderr,"[GRADIENT_PROJECT_FAIL] cell=%d scalar=%d constraints=%d iterations=%d\n",i,scalar,nc,it);
   FILE *fp=fopen("gradient_project_failure.csv","wx");if(fp){fprintf(fp,"gx,gy\n%.17g,%.17g\nax,ay,b\n",g0[0],g0[1]);for(int j=0;j<nc;j++)fprintf(fp,"%.17g,%.17g,%.17g\n",A[2*j],A[2*j+1],b[j]);fclose(fp);}goto done;
  }
 }
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
 done:free(off);free(grad);free(fit);free(A);free(b);return ok;
}
#endif
