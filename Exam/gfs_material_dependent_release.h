#ifndef GFS_MATERIAL_DEPENDENT_RELEASE_H
#define GFS_MATERIAL_DEPENDENT_RELEASE_H
#include <stdlib.h>
#include <math.h>
#include <float.h>
/* Pivoted Schur test: reject indefinite matrices; never regularize a null
 * direction. The tolerance measures roundoff against the original scale. */
static inline int gfs_material_release_psd(int m,const double*S,const double*bp){
 double *a=malloc(((size_t)m*m+2*(size_t)m)*sizeof(double));if(!a)return 0;
 double *rhs=a+m*m,*bs=rhs+m;for(int i=0;i<m;i++){rhs[i]=bp[i];bs[i]=fabs(bp[i]);}
 double scale=0;int ok=0;for(int i=0;i<m*m;i++){a[i]=S[i];if(!isfinite(a[i]))goto done;scale=fmax(scale,fabs(a[i]));}
 double tol=64*m*DBL_EPSILON*scale;
 for(int k=0;k<m;k++){
  int p=k;for(int i=k+1;i<m;i++)if(a[i*m+i]>a[p*m+p])p=i;
  if(a[p*m+p]<=tol){ok=1;for(int i=k;i<m;i++)for(int j=k;j<m;j++)if(fabs(a[i*m+j])>tol)ok=0;for(int i=k;i<m;i++)if(!isfinite(rhs[i])||fabs(rhs[i])>1e-10*bs[i])ok=0;goto done;}
  if(p!=k){double v=rhs[k];rhs[k]=rhs[p];rhs[p]=v;v=bs[k];bs[k]=bs[p];bs[p]=v;for(int j=0;j<m;j++){double z=a[k*m+j];a[k*m+j]=a[p*m+j];a[p*m+j]=z;}for(int i=0;i<m;i++){double z=a[i*m+k];a[i*m+k]=a[i*m+p];a[i*m+p]=z;}}
  for(int i=k+1;i<m;i++){double v=a[i*m+k]*rhs[k]/a[k*m+k];rhs[i]-=v;bs[i]+=fabs(v);}
  for(int i=k+1;i<m;i++)for(int j=k+1;j<m;j++)a[i*m+j]-=a[i*m+k]*a[k*m+j]/a[k*m+k];
 }
 ok=1;
 done:free(a);return ok;
}
/* Solve the small convex dual only to identify a releasable active inequality.
 * Coordinate minimization discovers the positive support; the original SPD
 * solver then certifies ALL original inequalities and pressure equalities on
 * the retained support. The Schur test also checks full pressure-drive range
 * compatibility. No hydrodynamic coefficient or regularization is added.
 * The caller removes ONE zero multiplier and repeats its existing physical
 * feasibility checks. A finite work budget rejects rather than accepting an
 * unconverged support. Algebra is independent of spatial dimension. */
static inline int gfs_material_dependent_release(int m,const double*S,
 const double*b,const double*bp,int (*solve)(int,const double*,const double*,double*)){
 if(m<2||!S||!b||!bp||!solve||!gfs_material_release_psd(m,S,bp))return -1;
 double *mem=calloc((size_t)m*m+6*(size_t)m,sizeof(double));int *ids=malloc(m*sizeof(int));
 if(!mem||!ids){free(mem);free(ids);return -1;}
 double *A=mem,*x=A+m*m,*r=x+m,*rhs=r+m,*xp=rhs+m,*rp=xp+m,*z=rp+m;int answer=-1;
 for(int i=0;i<m;i++){if(!(S[i*m+i]>0)||!isfinite(b[i])||!isfinite(bp[i]))goto done;r[i]=b[i];}
 for(int sweep=0;sweep<100*m+256;sweep++){
  for(int i=0;i<m;i++){
   double next=fmax(0,x[i]+r[i]/S[i*m+i]),delta=next-x[i];if(!isfinite(next))goto done;x[i]=next;
   for(int j=0;j<m;j++)r[j]-=S[j*m+i]*delta;
  }
  int converged=1;
  for(int i=0;i<m;i++){
   double v=b[i],scale=fabs(b[i]);for(int j=0;j<m;j++){v-=S[i*m+j]*x[j];scale+=fabs(S[i*m+j]*x[j]);}r[i]=v;
   if((x[i]>0?fabs(v):fmax(v,0))>1e-10*scale)converged=0;
  }
  if(!converged)continue;
  int n=0;for(int i=0;i<m;i++)if(x[i]>0)ids[n++]=i;
  if(n==0||n==m)goto done;
  for(int i=0;i<n;i++){rhs[i]=b[ids[i]];rp[i]=bp[ids[i]];for(int j=0;j<n;j++)A[i*n+j]=S[ids[i]*m+ids[j]];}
  if(!solve(n,A,rhs,z)||!solve(n,A,rp,xp))goto done;
  for(int j=0;j<n;j++)if(!isfinite(z[j])||z[j]<0)goto done;
  for(int i=0;i<m;i++){
   double v=b[i],pv=bp[i],scale=fabs(b[i]),ps=fabs(bp[i]);
   for(int j=0;j<n;j++){double a=S[i*m+ids[j]]*z[j],pa=S[i*m+ids[j]]*xp[j];v-=a;pv-=pa;scale+=fabs(a);ps+=fabs(pa);}
   if(!isfinite(v)||!isfinite(pv)||(x[i]>0?fabs(v):fmax(v,0))>1e-10*scale||(x[i]>0&&fabs(pv)>1e-10*ps))goto done;
   r[i]=v;
  }
  for(int i=0;i<m;i++)if(x[i]==0&&(answer<0||r[i]<r[answer]))answer=i;
  break;
 }
 done:free(mem);free(ids);return answer;
}
#endif
