#ifndef MATERIAL_GRADIENT_PROJECT_H
#define MATERIAL_GRADIENT_PROJECT_H
#include <math.h>
#include <float.h>
#ifdef __CUDACC__
#define GFS_GRADIENT_HD __host__ __device__
#else
#define GFS_GRADIENT_HD
#endif
/* Candidate: Euclidean closest gradient under normalized halfspaces A g <= b.
 * 1 <= dim <= 3, b >= 0 so the zero slope is feasible. No physical parameter.
 * Caller supplies geometry and pressure/density bounds in consistent units.
 * Active set has at most dim independent constraints. Failure is explicit.
 */
static GFS_GRADIENT_HD int mg_solve(int n,double m[3][3],double rhs[3],double x[3]){
 double a[3][4]={{0}};
 for(int i=0;i<n;i++){for(int j=0;j<n;j++)a[i][j]=m[i][j];a[i][n]=rhs[i];}
 for(int k=0;k<n;k++){
  int p=k;for(int i=k+1;i<n;i++)if(fabs(a[i][k])>fabs(a[p][k]))p=i;
  if(fabs(a[p][k])<128*DBL_EPSILON)return 0;
  for(int j=k;j<=n;j++){double q=a[k][j];a[k][j]=a[p][j];a[p][j]=q;}
  double v=a[k][k];for(int j=k;j<=n;j++)a[k][j]/=v;
  for(int i=0;i<n;i++)if(i!=k){v=a[i][k];for(int j=k;j<=n;j++)a[i][j]-=v*a[k][j];}
 }
 for(int i=0;i<n;i++)x[i]=a[i][n];
 return 1;
}
static GFS_GRADIENT_HD int mg_project(int dim,int n,const double *A,const double*b,const double*g0,double*g,int*iterations){
 if(dim<1||dim>3||n<0||!A||!b||!g0||!g||!iterations)return 0;
 double scale=0;for(int j=0;j<dim;j++){if(!isfinite(g0[j]))return 0;scale=fmax(scale,fabs(g0[j]));g[j]=0;}
 for(int i=0;i<n;i++){
  double nn=0;for(int j=0;j<dim;j++){if(!isfinite(A[dim*i+j]))return 0;nn+=A[dim*i+j]*A[dim*i+j];}
  if(!isfinite(b[i])||b[i]<0||fabs(nn-1)>512*DBL_EPSILON)return 0;
 }
 if(scale==0){*iterations=0;return 1;}
 double tol=512*DBL_EPSILON*scale;int active[3],na=0;
 for(int it=0;it<32*(n+1);it++){
  *iterations=it+1;double h[3]={0},d[3]={0},Q[3][3]={{0}},L[3][3]={{0}},t[3]={0},lambda[3]={0};
  for(int j=0;j<dim;j++)h[j]=g0[j]-g[j];
  /* QR of active normals avoids squaring the condition number in A A^T.
   * Nearly opposite faces can be independent even when their Gram matrix
   * loses that distinction to floating-point rounding. */
  for(int i=0;i<na;i++){
   for(int j=0;j<dim;j++)Q[i][j]=A[dim*active[i]+j];
   for(int pass=0;pass<2;pass++)for(int k=0;k<i;k++){
    double v=0;for(int j=0;j<dim;j++)v+=Q[i][j]*Q[k][j];L[i][k]+=v;
    for(int j=0;j<dim;j++)Q[i][j]-=v*Q[k][j];
   }
   double nn=0;for(int j=0;j<dim;j++)nn+=Q[i][j]*Q[i][j];nn=sqrt(nn);
   if(nn<128*DBL_EPSILON)return 0;L[i][i]=nn;
   for(int j=0;j<dim;j++){Q[i][j]/=nn;t[i]+=Q[i][j]*h[j];}
  }
  for(int i=na-1;i>=0;i--){lambda[i]=t[i];for(int k=i+1;k<na;k++)lambda[i]-=L[k][i]*lambda[k];lambda[i]/=L[i][i];}
  double norm=0;for(int j=0;j<dim;j++){d[j]=h[j];for(int i=0;i<na;i++)d[j]-=t[i]*Q[i][j];norm=fmax(norm,fabs(d[j]));}
  /* Full-rank active constraints leave no feasible tangent direction. */
  if(na==dim)norm=0;
  if(norm<=tol){
   int remove=-1;double worst=-tol;for(int i=0;i<na;i++)if(lambda[i]<worst){worst=lambda[i];remove=i;}
   if(remove>=0){active[remove]=active[--na];continue;}
   for(int i=0;i<n;i++){double v=-b[i];for(int j=0;j<dim;j++)v+=A[dim*i+j]*g[j];if(v>8*tol)return 0;}
   return 1;
  }
  double step=1;int blocker=-1;
  for(int i=0;i<n;i++){
   int exists=0;for(int k=0;k<na;k++)exists|=active[k]==i;if(exists)continue;
   double v=0,ad=0;for(int j=0;j<dim;j++){v+=A[dim*i+j]*g[j];ad+=A[dim*i+j]*d[j];}
   if(ad>tol){double s=fmax(0.0,(b[i]-v)/ad);if(s<step){step=s;blocker=i;}}
  }
  for(int j=0;j<dim;j++)g[j]+=step*d[j];
  if(blocker>=0){if(na==dim)return 0;active[na++]=blocker;}
 }
 return 0;
}
#undef GFS_GRADIENT_HD
#endif
