#ifndef GFS_MATERIAL_SMALL_GENERAL_H
#define GFS_MATERIAL_SMALL_GENERAL_H
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
/* Pivoted nonsymmetric solve; original-system residual certification.
 * No diagonal regularization, silent singular fallback or SPD assumption. */
static inline int gfs_material_small_general(int n,const double*A,const double*b,double*x){
 if(n<1||!A||!b||!x)return 0;
 double*M=malloc((size_t)n*(n+1)*sizeof(double));if(!M)return 0;
 int ok=0;double scale=0;
 for(int i=0;i<n;i++){for(int j=0;j<n;j++){M[i*(n+1)+j]=A[i*n+j];scale=fmax(scale,fabs(A[i*n+j]));}M[i*(n+1)+n]=b[i];}
 if(!(scale>0)||!isfinite(scale))goto done;
 for(int k=0;k<n;k++){
  int p=k;for(int i=k+1;i<n;i++)if(fabs(M[i*(n+1)+k])>fabs(M[p*(n+1)+k]))p=i;
  if(!(fabs(M[p*(n+1)+k])>64*DBL_EPSILON*scale))goto done;
  if(p!=k)for(int j=k;j<=n;j++){double v=M[k*(n+1)+j];M[k*(n+1)+j]=M[p*(n+1)+j];M[p*(n+1)+j]=v;}
  for(int i=k+1;i<n;i++){double a=M[i*(n+1)+k]/M[k*(n+1)+k];for(int j=k+1;j<=n;j++)M[i*(n+1)+j]-=a*M[k*(n+1)+j];}
 }
 for(int i=n-1;i>=0;i--){double a=M[i*(n+1)+n];for(int j=i+1;j<n;j++)a-=M[i*(n+1)+j]*x[j];x[i]=a/M[i*(n+1)+i];if(!isfinite(x[i]))goto done;}
 {double err=0,size=0;for(int i=0;i<n;i++){double r=-b[i];size+=fabs(b[i]);for(int j=0;j<n;j++){r+=A[i*n+j]*x[j];size+=fabs(A[i*n+j]*x[j]);}err+=fabs(r);}ok=isfinite(err)&&err<=1e-10*size;}
 done:free(M);return ok;
}
#endif
