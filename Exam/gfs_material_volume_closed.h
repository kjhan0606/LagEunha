#ifndef GFS_MATERIAL_VOLUME_CLOSED_H
#define GFS_MATERIAL_VOLUME_CLOSED_H
#include "gfs_material_volume_inverse.h"
static int gfs_volume_component_root(int*p,int i){while(p[i]!=i){p[i]=p[p[i]];i=p[i];}return i;}
/* Fix one constant-weight gauge per connected component. Only the auxiliary
 * solve replaces a row. Certify ALL original equations afterwards, including
 * the gauge rows. Never subtract a RHS mean or change a target volume. */
static int gfs_material_volume_closed_correction_scaled(int n,const long long*ids,
 const int*off,const int*owner,const double*c,const double*rhs,double*x,
 GfsVolumeLinearStats*stats,const double*target_volume,double volume_tolerance){
 if(!stats)return 0;memset(stats,0,sizeof(*stats));
 if(n<1||!ids||!off||!owner||!c||!rhs||!x||off[0]!=0||off[n]<0||off[n]>INT_MAX-n)return 0;
 if(target_volume){
  if(!(volume_tolerance>0)||!isfinite(volume_tolerance))return 0;
  for(int i=0;i<n;i++)if(!(target_volume[i]>0)||!isfinite(target_volume[i]))return 0;
 }

 int nf=off[n],ok=0,components=0;
 int *parent=malloc(n*sizeof(int)),*pin=malloc(n*sizeof(int)),*goff=malloc(((size_t)n+1)*sizeof(int));
 int *gowner=malloc(((size_t)nf+n)*sizeof(int));double *gc=malloc(((size_t)nf+n)*sizeof(double)),*b=malloc(n*sizeof(double));
 if(!parent||!pin||!goff||!gowner||!gc||!b)goto done;
 for(int i=0;i<n;i++){parent[i]=i;pin[i]=-1;if(off[i+1]<off[i]||off[i+1]>nf||!isfinite(rhs[i]))goto done;}
 for(int i=0;i<n;i++)for(int k=off[i];k<off[i+1];k++){
  int j=owner[k];if(j<0||j>=n||!(c[k]>=0)||!isfinite(c[k]))goto done;
  if(c[k]>0&&j!=i){int a=gfs_volume_component_root(parent,i),z=gfs_volume_component_root(parent,j);if(a!=z)parent[z]=a;}
 }
 for(int i=0;i<n;i++){int r=gfs_volume_component_root(parent,i);if(pin[r]<0||ids[i]<ids[pin[r]])pin[r]=i;}
 int count=0;
 for(int i=0;i<n;i++){
  goff[i]=count;b[i]=rhs[i];
  if(pin[gfs_volume_component_root(parent,i)]==i){gowner[count]=-1;gc[count++]=1;b[i]=0;components++;}
  else for(int k=off[i];k<off[i+1];k++){gowner[count]=owner[k];gc[count++]=c[k];}
 }
 goff[n]=count;
 int reached=0;
 if(gfs_volume_anchor_graph(n,goff,gowner,gc,&reached)!=1)goto done;
 if(!gfs_lu_prepare_basis(n,ids,goff,gowner,gc))goto done;
 if(!gfs_material_volume_linear_correction(n,goff,gowner,gc,b,x,stats))goto done;
 stats->relative_residual=gfs_lu_true_residual(n,off,owner,c,rhs,x);
 /* The pinned auxiliary solve retains its 1e-10 algebra certificate.
  * Closed geometry targets have finite precision and need not sum exactly.
  * Certify EVERY original equation in the units of the requested geometric
  * accuracy, rather than dividing by a correction that tends to zero.
  * No RHS projection, mean removal, or target-volume change is performed.
  * The caller must recompute geometry and enforce the SAME physical tolerance. */
 double physical_residual=0;
 if(target_volume){
  ok=1;
  for(int i=0;i<n;i++){
   long double ax=0;
   for(int k=off[i];k<off[i+1];k++)ax+=(long double)c[k]*((long double)x[i]-x[owner[k]]);
   long double r=fabsl(ax-(long double)rhs[i])/(long double)target_volume[i];
   if(!isfinite(r)||r>volume_tolerance)ok=0;
   if(r>physical_residual)physical_residual=(double)r;
  }
 }else ok=isfinite(stats->relative_residual)&&stats->relative_residual<=1e-10;
 fprintf(stderr,"[MATERIAL_VOLUME_CLOSED_LINEAR] components=%d original_residual=%.17g physical_residual=%.17g physical_tolerance=%.17g accepted=%d\n",components,stats->relative_residual,physical_residual,target_volume?volume_tolerance:0,ok);

done:
 free(b);free(gc);free(gowner);free(goff);free(pin);free(parent);return ok;
}
static int gfs_material_volume_closed_correction(int n,const long long*ids,
 const int*off,const int*owner,const double*c,const double*rhs,double*x,
 GfsVolumeLinearStats*stats){
 return gfs_material_volume_closed_correction_scaled(n,ids,off,owner,c,rhs,x,stats,NULL,0);
}
#endif
