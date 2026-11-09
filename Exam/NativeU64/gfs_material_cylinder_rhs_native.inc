/* Native Cylinder routing: same certified mixed-boundary C/CUDA operator. */
static int gfs_material_wall_power_valid;
static double gfs_material_wall_power_value;
int gfs_material_native_wall_power(double *power){
 if(!power||!gfs_material_wall_power_valid)return 0;
 *power=gfs_material_wall_power_value;return 1;
}
#ifdef GFS_MATERIAL_NATIVE_EXTERIOR
#include "gfs_material_cylinder_rhs2d.h"
#endif
static int gfs_material_cylinder_native_rhs(SimParameters*s,GPUContext*ctx,int n,int nt){
#ifdef GFS_MATERIAL_NATIVE_EXTERIOR
 const GfsCylinderMaterialSites*b=gfs_material_cylinder_sites(n,nt,ctx->h_parts.indx);
 GfsMaterialCylinderNativeStats2D stats;
 if(!b||NID(s)!=1)return 0;
 double *vdot=malloc((size_t)n*sizeof(double));if(!vdot)return 0;
 int ok=gfs_material_cylinder_native_rhs2d(ctx,n,nt,GAS_GAMMA(s),b->center,b->radius,
  b->alias,b->image,b->lowx,b->lowy,gfs_project_one_scalar,
#ifdef GFS_MATERIAL_CUDA_RHS
  gfs_material_cuda_rhs_moving_walls2d,
#else
  gfs_material_rhs_moving_walls2d,
#endif
  vdot,&stats);
 if(ok&&gfs_material_native_volume_audit.row){
  if(!gfs_material_native_volume_sample(n,ctx->h_parts.indx,ctx->h_parts.weight_volume64,vdot)){
   fprintf(stderr,"[MATERIAL_VOLUME_AUDIT] invalid stage/ID sampling\n");
   free(vdot);MPI_Abort(MPI_COMM(s),99);return 0;
  }
 }
#ifdef GFS_MATERIAL_FIXED_SOLID
 if(ok&&getenv("GFS_MATERIAL_FIXED_EXPORT")){
  FILE*ff=fopen("fixed_faces.bin","wb"),*sf=fopen("fixed_state.bin","wb");
  if(!ff||!sf)MPI_Abort(MPI_COMM(s),99);
  int header[2]={n,ctx->h_faces.face_offset[n]};
  fwrite("FIXEDB01",1,8,ff);fwrite(header,sizeof(header),1,ff);
  for(int i=0;i<n;i++){
   ParticleSoA*p=&ctx->h_parts;FaceCSR*f=&ctx->h_faces;
   double state[9]={p->x[i],p->y[i],p->vx[i],p->vy[i],p->mass[i],p->weight_U64[i],p->weight_volume64[i],p->weight_W[i],(double)p->indx[i]};
   fwrite(state,sizeof(state),1,sf);
   for(int k=f->face_offset[i];k<f->face_offset[i+1];k++){
    int id[2]={i,f->neighbor_idx[k]};double edge[4]={f->c1x[k],f->c1y[k],f->c2x[k],f->c2y[k]};
    fwrite(id,sizeof(id),1,ff);fwrite(edge,sizeof(edge),1,ff);
   }
  }
  if(fclose(ff)||fclose(sf))MPI_Abort(MPI_COMM(s),99);
  fprintf(stderr,"[FIXED_SOLID_EXPORT_COMPLETE] cells=%d faces=%d wall_power=%.17g evolved_steps=0\n",n,header[1],stats.exchange.wall_power);
  fflush(NULL);MPI_Finalize();exit(0);
 }
#endif
 free(vdot);
 if(!ok)return 0;
 if(stats.exchange.volume_error>1e-8||fabs(stats.exchange.energy_defect)>1e-8)return 0;
 gfs_material_last_face_crossing=stats.laboratory_crossing_time;
 gfs_material_last_geometry_bound=stats.geometry_timestep;
 gfs_material_last_damping_bound=0;gfs_material_last_potential_valid=0;
 memcpy(gfs_material_last_ledger,stats.thermal_ledger,sizeof(stats.thermal_ledger));
 gfs_material_boundary_power_value=stats.exchange.wall_power+stats.exchange.exterior_ledger[3];
 gfs_material_boundary_power_valid=1;
 gfs_material_wall_power_value=stats.exchange.wall_power;gfs_material_wall_power_valid=1;
 double t;int phase;if(!gfs_material_native_evaluation_time(&t,&phase))return 0;
 fprintf(stderr,"[CYLINDER_NATIVE_RHS] t=%.17g phase=%d owned=%d exterior=%d faces=%d volume_error=%.17g energy_defect=%.17g wall_power=%.17g exterior_power=%.17g laboratory_crossing=%.17g iterations=%d\n",
  t,phase,n,nt-n,stats.faces,stats.exchange.volume_error,stats.exchange.energy_defect,
  stats.exchange.wall_power,stats.exchange.exterior_ledger[3],stats.laboratory_crossing_time,stats.exchange.weight_iterations);
#ifdef GFS_MATERIAL_SWEPT_CFL
 if(stats.geometry_sweep_coefficient>0){
  int i=stats.geometry_limiting_cell;
  fprintf(stderr,"[MATERIAL_GEOMETRIC_CFL] t=%.17g phase=%d coefficient=%.17g dt=%.17g id=%lld\n",t,phase,
   stats.geometry_sweep_coefficient,stats.geometry_timestep,i>=0?(long long)ctx->h_parts.indx[i]:-1LL);
 }
#endif
 return 1;
#else
 (void)s;(void)ctx;(void)n;(void)nt;
 fprintf(stderr,"[CYLINDER_NATIVE_RHS] requires GFS_MATERIAL_NATIVE_EXTERIOR build\n");return 0;
#endif
}
