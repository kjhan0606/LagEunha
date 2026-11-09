#include "gfs_material_stage_time.h"
#include "../gfs_material_volume_audit.h"
static GfsMaterialVolumeAudit gfs_material_native_volume_audit;
int gfs_material_native_volume_audit_begin(int n){
 return gfs_material_volume_audit_start(&gfs_material_native_volume_audit,n);
}
int gfs_material_native_volume_audit_rewind_first(int n){
 return gfs_material_volume_audit_rewind_first(&gfs_material_native_volume_audit,n);
}
void gfs_material_native_volume_audit_end(void){
 gfs_material_volume_audit_free(&gfs_material_native_volume_audit);
}
int gfs_material_native_volume_audit_cell(long long id,double V,double dt,double*defect,double*relative){
 return gfs_material_volume_audit_cell(&gfs_material_native_volume_audit,id,V,dt,defect,relative);
}
static int gfs_material_native_volume_sample(int n,const long long*ids,
 const double*volume,const double*rate){
 if(!gfs_material_native_volume_audit.row)return 1;
 double time;int phase;
 if(!gfs_material_native_evaluation_time(&time,&phase))return 0;
 /* RT evaluates the endpoint for its gravitational-potential diagnostic.
  * It is not a third SSP stage and must not overwrite either stored rate. */
 if(phase==3)return gfs_material_native_volume_audit.phase==2&&gfs_material_native_volume_audit.n==n;
 if(!gfs_material_volume_audit_stage(&gfs_material_native_volume_audit,phase,n,ids,volume,rate))return 0;
 long double V=0,Q=0,absolute=0;
 for(int i=0;i<n;i++){V+=volume[i];Q+=rate[i];absolute+=fabsl(rate[i]);}
 fprintf(stderr,"[MATERIAL_VOLUME_STAGE_SUM] t=%.17g phase=%d cells=%d volume=%.21Lg rate=%.21Lg absolute_rate=%.21Lg\n",time,phase,n,V,Q,absolute);
 return 1;
}
