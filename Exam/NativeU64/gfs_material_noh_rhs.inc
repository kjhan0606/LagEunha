/* Native routing only: reuses the already qualified exterior face operator. */
static int gfs_material_boundary_power_valid;
static double gfs_material_boundary_power_value;
int gfs_material_native_boundary_power(double *power){
 if(!power||!gfs_material_boundary_power_valid)return 0;
 *power=gfs_material_boundary_power_value;return 1;
}
#ifdef GFS_MATERIAL_NATIVE_EXTERIOR
#include "gfs_material_exterior_native2d.h"
#include "../gfs_material_scalar_project2d.h"
static int gfs_material_noh_prescribe(double x,double y,void*data,GfsFvPrimitive*out){
 double point[3]={x,y,0};return gfs_noh_exterior_state(data,point,out);
}
#endif
static int gfs_material_noh_native_rhs(SimParameters*s,GPUContext*ctx,int n,int nt){
 gfs_material_boundary_power_valid=0;
#ifdef GFS_MATERIAL_NATIVE_EXTERIOR
 GfsNohExterior b;GfsMaterialExteriorNativeStats2D stats;
 if(!gfs_material_noh_configuration(s,&b)||nt!=4*n)return 0;
 double *offset=calloc((size_t)2*nt,sizeof(double)),*known=calloc(nt,sizeof(double));
 double*vdot=gfs_material_native_volume_audit.row?malloc((size_t)n*sizeof(double)):NULL;
 if(!offset||!known||(gfs_material_native_volume_audit.row&&!vdot)){free(offset);free(known);free(vdot);return 0;}
 int ok=gfs_material_exterior_native_rhs2d(ctx,n,nt,b.gamma,gfs_material_noh_prescribe,&b,
  offset,known,gfs_project_one_scalar,
#ifdef GFS_MATERIAL_CUDA_RHS
  1,
#else
  0,
#endif
  vdot,&stats);
 if(ok&&vdot)ok=gfs_material_native_volume_sample(n,ctx->h_parts.indx,ctx->h_parts.weight_volume64,vdot);
 free(vdot);free(offset);free(known);if(!ok)return 0;
 gfs_material_last_face_crossing=stats.face_crossing_time;gfs_material_last_damping_bound=0;
 gfs_material_last_geometry_bound=stats.geometry_timestep;
 gfs_material_last_potential_valid=0;memcpy(gfs_material_last_ledger,stats.thermal_ledger,sizeof(stats.thermal_ledger));
 gfs_material_boundary_power_value=stats.boundary_ledger[3];gfs_material_boundary_power_valid=1;
#ifdef GFS_MATERIAL_SWEPT_CFL
 if(stats.geometry_sweep_coefficient>0){
  double time;int phase;if(!gfs_material_native_evaluation_time(&time,&phase))return 0;
  int i=stats.geometry_limiting_cell;
  fprintf(stderr,"[MATERIAL_GEOMETRIC_CFL] t=%.17g phase=%d coefficient=%.17g dt=%.17g id=%lld\n",time,phase,
   stats.geometry_sweep_coefficient,stats.geometry_timestep,i>=0?(long long)ctx->h_parts.indx[i]:-1LL);
 }
#endif
 if(getenv("GFS_MATERIAL_STAGE_TIME_AUDIT"))fprintf(stderr,
  "[NOH_NATIVE_RHS] t=%.17g owned=%d exterior=%d faces=%d boundary_power=%.17g energy_rate_defect=%.17g weight_residual=%.17g\n",
  b.time,n,nt-n,stats.faces,stats.boundary_ledger[3],stats.energy_rate_defect,stats.weight_residual);
 return 1;
#else
 (void)s;(void)ctx;(void)n;(void)nt;
 fprintf(stderr,"[NOH_NATIVE_RHS] requires explicit GFS_MATERIAL_NATIVE_EXTERIOR build\n");return 0;
#endif
}
