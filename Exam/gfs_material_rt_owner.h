#ifndef GFS_MATERIAL_RT_OWNER_H
#define GFS_MATERIAL_RT_OWNER_H
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <stdlib.h>
/* Sidecar order is exactly the padding-array order used by fillParticleSoA.
 * A reflected image retains its native ghost sentinel in the particle ID. */
typedef struct {long long owner_id;int reflected_y;} GfsMaterialRtOwner;
static inline int gfs_material_rt_resolve_owners(int n,int nt,
 const long long *id,const double *x,const double *y,
 int count,const GfsMaterialRtOwner *pad,double lx,double ymin,double ymax,
 int *alias,unsigned char *reflection){
 if(n<1||nt<n||count!=nt-n||!id||!x||!y||!alias||!reflection||
    (count&&!pad)||!(lx>0)||!isfinite(lx)||!isfinite(ymin)||!isfinite(ymax)||!(ymax>ymin))return 0;
 size_t cap=1;while(cap<2*(size_t)n)cap*=2;
 int *slots=malloc(cap*sizeof(*slots));if(!slots)return 0;
 for(size_t k=0;k<cap;k++)slots[k]=-1;
 int ok=0;
 for(int i=0;i<n;i++){
  if(id[i]<0||!isfinite(x[i])||!isfinite(y[i]))goto done;
  size_t k=((uint64_t)id[i]*UINT64_C(11400714819323198485))&(cap-1);
  while(slots[k]>=0){if(id[slots[k]]==id[i])goto done;k=(k+1)&(cap-1);}
  slots[k]=i;alias[i]=i;reflection[i]=0;
 }
 for(int i=n;i<nt;i++){
  const GfsMaterialRtOwner *p=pad+i-n;
  if(p->owner_id<0||(p->reflected_y!=0&&p->reflected_y!=1)||!isfinite(x[i])||!isfinite(y[i]))goto done;
  size_t k=((uint64_t)p->owner_id*UINT64_C(11400714819323198485))&(cap-1);
  while(slots[k]>=0&&id[slots[k]]!=p->owner_id)k=(k+1)&(cap-1);
  if(slots[k]<0)goto done;
  int j=slots[k];double dx=x[i]-x[j];
  double tol=128*DBL_EPSILON*(lx+ymax-ymin+fabs(x[i])+fabs(y[i])+fabs(x[j])+fabs(y[j]));
  if(fabs(dx-lx*nearbyint(dx/lx))>tol)goto done;
  if(p->reflected_y){
   if(fmin(fabs(y[i]+y[j]-2*ymin),fabs(y[i]+y[j]-2*ymax))>tol)goto done;
  }else if(id[i]!=id[j]||fabs(y[i]-y[j])>tol)goto done;
  alias[i]=j;reflection[i]=(unsigned char)p->reflected_y;
 }
 ok=1;
 done:free(slots);return ok;
}
#endif
