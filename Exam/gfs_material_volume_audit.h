#ifndef GFS_MATERIAL_VOLUME_AUDIT_H
#define GFS_MATERIAL_VOLUME_AUDIT_H
#include <math.h>
#include <stdlib.h>

/* Dimension-independent, read-only finite-step geometric conservation audit.
 * Match immutable material IDs; neither particle order nor cell volume is
 * advanced here. Rebegin after a rejected trial to discard both stage samples.
 */
typedef struct {long long id;double volume,rate[2];int second;} GfsMaterialVolumeRow;
typedef struct {int n,phase;GfsMaterialVolumeRow*row;} GfsMaterialVolumeAudit;
static inline int gfs_material_volume_order(const void*a,const void*b){
 long long x=((const GfsMaterialVolumeRow*)a)->id,y=((const GfsMaterialVolumeRow*)b)->id;
 return (x>y)-(x<y);
}
static inline void gfs_material_volume_audit_free(GfsMaterialVolumeAudit*a){
 if(a){free(a->row);a->row=NULL;a->n=a->phase=0;}
}
static inline int gfs_material_volume_audit_start(GfsMaterialVolumeAudit*a,int n){
 if(!a||n<1)return 0;
 gfs_material_volume_audit_free(a);
 a->row=calloc((size_t)n,sizeof(*a->row));if(!a->row)return 0;
 a->n=n;return 1;
}
/* Preserve a certified first evaluation of the identical step-start state.
 * Discard every second-stage sample, including a partially failed stage.
 * The caller must separately certify identical state, time and geometry. */
static inline int gfs_material_volume_audit_rewind_first(GfsMaterialVolumeAudit*a,int n){
 if(!a||!a->row||a->n!=n||n<1||a->phase<1||a->phase>2)return 0;
 for(int i=0;i<n;i++){a->row[i].rate[1]=0;a->row[i].second=0;}
 a->phase=1;return 1;
}
static inline int gfs_material_volume_audit_stage(GfsMaterialVolumeAudit*a,int phase,
 int n,const long long*ids,const double*volume,const double*rate){
 if(!a||!a->row||n!=a->n||!ids||!volume||!rate||phase!=a->phase+1||phase>2)return 0;
 for(int i=0;i<n;i++)if(!(volume[i]>0)||!isfinite(volume[i])||!isfinite(rate[i]))return 0;
 if(phase==1){
  for(int i=0;i<n;i++)a->row[i]=(GfsMaterialVolumeRow){ids[i],volume[i],{rate[i],0},0};
  qsort(a->row,n,sizeof(*a->row),gfs_material_volume_order);
  for(int i=1;i<n;i++)if(a->row[i-1].id==a->row[i].id)return 0;
 }else{
  for(int i=0;i<n;i++){
   GfsMaterialVolumeRow key={0};key.id=ids[i];
   GfsMaterialVolumeRow*r=bsearch(&key,a->row,n,sizeof(*a->row),gfs_material_volume_order);
   if(!r||r->second)return 0;
   r->rate[1]=rate[i];r->second=1;
  }
 }
 a->phase=phase;return 1;
}
static inline int gfs_material_volume_audit_cell(const GfsMaterialVolumeAudit*a,
 long long id,double endpoint,double dt,double*defect,double*relative){
 if(!a||a->phase!=2||!a->row||!(endpoint>0)||!isfinite(endpoint)||
    !(dt>0)||!isfinite(dt)||!defect||!relative)return 0;
 GfsMaterialVolumeRow key={0};key.id=id;
 const GfsMaterialVolumeRow*r=bsearch(&key,a->row,a->n,sizeof(*a->row),gfs_material_volume_order);
 if(!r||!r->second)return 0;
 *defect=endpoint-r->volume-.5*dt*(r->rate[0]+r->rate[1]);
 *relative=fabs(*defect)/r->volume;
 return isfinite(*defect)&&isfinite(*relative);
}
#endif
