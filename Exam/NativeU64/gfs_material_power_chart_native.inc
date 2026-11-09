extern void gpu_power_chart_domain(double,double);

static double chart_ymin,chart_ymax;
/* Existing RT sidecar resolves mirror sentinel IDs before any GPU geometry.
 * Do not treat a mirror sentinel as an ordinary periodic particle ID. */
const int *gfs_power_chart_owners(GPUContext *ctx,int n,int nt,double lx,double ly){
 if(gfs_material_cylinder_mode()){
  const GfsCylinderMaterialSites*b=gfs_material_cylinder_sites(n,nt,ctx->h_parts.indx);
  return b?b->alias:NULL;
 }
 if(ly>0)return gpu_unified_owner_aliases(ctx,n,nt,lx,ly);
 static int *alias=NULL;static unsigned char *reflection=NULL;static int capacity=0;
 if(nt>capacity){free(alias);free(reflection);alias=malloc(nt*sizeof(int));reflection=malloc(nt);capacity=nt;}
 if(!alias||!reflection)return NULL;
 int count=0;const GfsMaterialRtOwner *pad=gfs_material_rt_owner_map(&count);ParticleSoA *a=&ctx->h_parts;
 if(!gfs_material_rt_resolve_owners(n,nt,a->indx,a->x,a->y,count,pad,lx,chart_ymin,chart_ymax,alias,reflection))return NULL;
 for(int i=0;i<nt;i++){
  if(reflection[i])alias[i]=-1;
  else if(a->weight_W[i]!=a->weight_W[alias[i]])return NULL;
 }
 return alias;
}

static void gfs_material_power_chart_config(SimParameters*simpar){
    if(NID(simpar)!=1){fprintf(stderr,"[POWER_CHART] multi-rank not yet qualified\n");MPI_Abort(MPI_COMM(simpar),99);}
    const char *chart_rt=getenv("GFS_MATERIAL_RT_NATIVE");
    int chart_wall=chart_rt&&!strcmp(chart_rt,"1");
    chart_ymin=Ymin_HydroExam(simpar);chart_ymax=Ymax_HydroExam(simpar);
    if(gfs_material_noh_native_mode())gpu_power_chart_domain(0,0);
    else gpu_power_chart_domain(chart_wall?Xmax_HydroExam(simpar)-Xmin_HydroExam(simpar):SIMBOX(simpar).x.max-SIMBOX(simpar).x.min,
        chart_wall?0:SIMBOX(simpar).y.max-SIMBOX(simpar).y.min);

}
