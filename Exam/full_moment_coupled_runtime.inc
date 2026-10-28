/* Private adapter: native MPI owns the nonlinear time transaction. Providers
 * remain geometry/fit/audit services. Included after Runtime helpers so the
 * legacy explicit control path retains its existing implementation. */
#include <limits.h>
typedef struct {
 Runtime *r;int step;double started,gamma;double *base_area;
 FullMomentState midpoint;
 char folder[512],base_metadata[512],provider_capture[512],mapped_path[512];
 FullMomentPotentialNorm *anderson_metric;
 size_t metric_cells;
 int metric_epoch,metric_initialized;
} CoupledRuntime;
static int coupled_read(const char *root,const char *name,void *out,size_t bytes){
 char file[1024];snprintf(file,sizeof file,"%s/%s.bin",root,name);FILE *f=fopen(file,"rb");
 if(!f)return 1;int status=fread(out,1,bytes,f)!=bytes;status|=fgetc(f)!=EOF;status|=fclose(f)!=0;return status;
}
static int coupled_metadata(Runtime *r,char *out){
 int status=0;if(!r->rank){
  status=path(r->request.fits_path,r->fits);FullMomentProviderResponse response;
  if(!status)status=rpc(r,FMP_COUPLED_METADATA,&response);
  if(!status)status=path(out,response.output_path);
 }
 if(collective(r,status))return 1;MPI_Bcast(out,512,MPI_CHAR,0,r->comm);return 0;
}
static uint64_t coupled_hash_bytes(uint64_t hash,const void *data,size_t size){
 const unsigned char *p=data;
 for(size_t i=0;i<size;i++){hash^=p[i];hash*=1099511628211ULL;}
 return hash;
}
static uint64_t coupled_hash_long_double(uint64_t hash,long double value){
 char text[128];int n=snprintf(text,sizeof text,"%.21La",value);
 if(n<0||(size_t)n>=sizeof text)return 0;
 return coupled_hash_bytes(hash,text,(size_t)n+1);
}
static uint64_t coupled_metric_fingerprint(const FullMomentPotentialNorm *metric,size_t n){
 uint64_t hash=14695981039346656037ULL;
 for(size_t i=0;i<n;i++){
  hash=coupled_hash_bytes(hash,&metric[i].dependent,sizeof metric[i].dependent);
  for(int j=0;j<2;j++)hash=coupled_hash_long_double(hash,metric[i].origin[j]);
  hash=coupled_hash_long_double(hash,metric[i].scale);
  for(int j=0;j<6;j++)hash=coupled_hash_long_double(hash,metric[i].reference[j]);
  if(!metric[i].dependent){
   for(int j=0;j<6;j++)hash=coupled_hash_long_double(hash,metric[i].projection[j]);
   hash=coupled_hash_long_double(hash,metric[i].sigma);
   hash=coupled_hash_long_double(hash,metric[i].rank_error);
  }
 }
 return hash;
}
static void coupled_metric_summary(const FullMomentPotentialNorm *metric,size_t n,
 FullMomentCoupledResidual *residual){
 long double sigma_min=INFINITY,sigma_max=0.,margin_min=INFINITY;int dependent=0;
 for(size_t i=0;i<n;i++){
  if(metric[i].dependent){dependent++;continue;}
  if(metric[i].sigma<sigma_min)sigma_min=metric[i].sigma;
  if(metric[i].sigma>sigma_max)sigma_max=metric[i].sigma;
  long double margin=metric[i].sigma-metric[i].rank_error;
  if(margin<margin_min)margin_min=margin;
 }
 residual->dependent_cells=dependent;
 residual->actual_sigma_min=dependent==(int)n?0.:sigma_min;
 residual->actual_sigma_max=sigma_max;
 residual->actual_rank_margin_min=dependent==(int)n?0.:margin_min;
}
static int coupled_type2_fit(void *context,size_t count,int maximum_depth,
 const long double *current,const long double *const history[3],double gamma[3],
 int *used_depth,int *qr_rank,long double *predicted_norm2,long double *condition){
 CoupledRuntime *c=context;Runtime *r=c->r;int fit=0;
 if(!r->rank)fit=full_moment_coupled_type2_fit(count,maximum_depth,current,history,
   gamma,used_depth,qr_rank,predicted_norm2,condition);
 int status=0;
 status|=MPI_Bcast(&fit,1,MPI_INT,0,r->comm)!=MPI_SUCCESS;
 status|=MPI_Bcast(gamma,3,MPI_DOUBLE,0,r->comm)!=MPI_SUCCESS;
 status|=MPI_Bcast(used_depth,1,MPI_INT,0,r->comm)!=MPI_SUCCESS;
 status|=MPI_Bcast(qr_rank,1,MPI_INT,0,r->comm)!=MPI_SUCCESS;
 status|=MPI_Bcast(predicted_norm2,1,MPI_LONG_DOUBLE,0,r->comm)!=MPI_SUCCESS;
 status|=MPI_Bcast(condition,1,MPI_LONG_DOUBLE,0,r->comm)!=MPI_SUCCESS;
 if(status)return -1;
 return fit;
}
static int coupled_collective_status(void *context,int local_status){
 return collective(((CoupledRuntime*)context)->r,local_status);
}
static int coupled_midpoint_write(CoupledRuntime *c,const FullMomentState *base,
 const FullMomentState *trial,uint64_t epoch,char *out){
 Runtime *r=c->r;int status=0;
 snprintf(out,512,"%s/states/e%llu",r->memory,(unsigned long long)epoch);
 if(!r->rank){
  size_t n=base->header.cells;double *raw=malloc((20*n+4)*sizeof(double));
  double *factors=malloc(4*n*sizeof(double));long double *Psi=calloc(n,sizeof(long double));
  status=!raw||!factors||!Psi;int cell=-1;
  if(!status)status=full_moment_coupled_midpoint((int)n,base->state,trial->state,base->mass,base->Psi,trial->Psi,raw,Psi,factors,&cell);
  if(!status){
   FullMomentState mid={0};mid.header=base->header;mid.header.epoch=epoch;mid.header.time=(c->step-.5)*r->h;
   mid.ids=base->ids;mid.mass=base->mass;mid.state=raw;mid.factors=factors;mid.Psi=Psi;
   mid.weight_rate_seed=trial->weight_rate_seed;
   status=full_moment_state_write(out,&mid);
  }
  if(status)fprintf(stderr,"[COUPLED_MIDPOINT_REJECTED] epoch=%llu cell=%d code=%d\n",(unsigned long long)epoch,cell,status);
  free(raw);free(factors);free(Psi);
 }
 return collective(r,status);
}
static int coupled_residual_record(CoupledRuntime *c,int iteration,int code,
 const FullMomentCoupledResidual *residual,int failed){
 Runtime *r=c->r;int status=0;if(!r->rank){
  char file[1024];snprintf(file,sizeof file,"%s/iteration_%02d.json",c->folder,iteration);
  FILE *f=fopen(file,"wx");if(!f)status=1;else{
   fprintf(f,"{\"iteration\":%d,\"physical_provider_count\":%d,\"provider_epoch\":%llu,\"map_epoch\":%llu,\"residual_status\":%d,\"failed_cell\":%d,\"converged\":%s,\"maximum\":[",iteration,r->providers,(unsigned long long)r->epoch,(unsigned long long)r->allocated_epoch,code,failed,!code&&full_moment_coupled_converged(residual)?"true":"false");
   for(int i=0;i<FM_RES_COUNT;i++)fprintf(f,"%s\"%.21Lg\"",i?",":"",residual->maximum[i]);
   fprintf(f,"],\"cell\":[");for(int i=0;i<FM_RES_COUNT;i++)fprintf(f,"%s%lld",i?",":"",residual->cell[i]);
   fprintf(f,"],\"dependent_cells\":%d,\"fixed_metric_norm2\":\"%.21Lg\",\"actual_psi_remainder_numerator_max\":\"%.21Lg\",\"actual_sigma_min\":\"%.21Lg\",\"actual_sigma_max\":\"%.21Lg\",\"actual_rank_margin_min\":\"%.21Lg\",\"metric_epoch\":%d,\"metric_reset_cells\":%d,\"actual_metric_fingerprint\":\"%016llx\",\"anderson_metric_fingerprint\":\"%016llx\",\"accepted\":false}\n",residual->dependent_cells,residual->vector_norm2,residual->psi_remainder_numerator_max,residual->actual_sigma_min,residual->actual_sigma_max,residual->actual_rank_margin_min,residual->metric_epoch,residual->metric_reset_cells,(unsigned long long)residual->actual_metric_fingerprint,(unsigned long long)residual->anderson_metric_fingerprint);
   status=fflush(f)!=0;if(!status)status=fsync(fileno(f))!=0;status|=fclose(f)!=0;
  }
 }
 return collective(r,status);
}
static int coupled_evaluate(void *context,const FullMomentState *base,
 const FullMomentState *trial,int iteration,FullMomentState *mapped,FullMomentCoupledResidual *residual,
 long double *solver_residual,size_t solver_count){
 CoupledRuntime *c=context;Runtime *r=c->r;int status=0;memset(residual,0,sizeof *residual);
 if(!solver_residual||solver_count!=(size_t)r->cells*21)return 1;
 if(solver_count>(size_t)INT_MAX)return 1;
 memset(solver_residual,0,solver_count*sizeof(*solver_residual));
 if(release(r))return 1;
 full_moment_state_close(&c->midpoint);
 uint64_t epoch=++r->allocated_epoch;char mid_path[512],mid_geometry[512]={0};
 if(coupled_midpoint_write(c,base,trial,epoch,mid_path))return 1;
 if(geometry(r,epoch,mid_path,mid_geometry)||configure(r,epoch,mid_path,mid_geometry,0))return 1;
 char label[96];snprintf(label,sizeof label,"coupled_midpoint_%04d_%02d",c->step,iteration);
 if(admissibility(r,c->step,0,label))return 1;
 char metadata[512]={0};if(coupled_metadata(r,metadata))return 1;
 status=full_moment_state_load(mid_path,&c->midpoint);if(collective(r,status))return 1;
 r->request.step=c->step;r->request.stage=iteration;
 if(rhs(r))return 1;
 uint64_t mapped_epoch=++r->allocated_epoch;
 snprintf(c->mapped_path,sizeof c->mapped_path,"%s/states/e%llu",r->memory,(unsigned long long)mapped_epoch);
 if(full_moment_native_factor(r->comm,&r->service,&r->request,base,&c->midpoint,r->native_stage,r->root,r->rates,c->mapped_path,mapped_epoch,r->h,c->step*r->h,r->threads,1))return 1;
 if(save_factor_audit(r,mapped_epoch))return 1;
 status=full_moment_state_load(c->mapped_path,mapped);if(collective(r,status))return 1;
 int failed=-1;
 if(!r->rank){
  FullMomentPotentialNorm *potential=malloc((size_t)r->cells*sizeof(*potential));status=!potential;
  if(!status)status=coupled_read(metadata,"potential_norms",potential,(size_t)r->cells*sizeof(*potential));
  int reset_cells=0,new_metric=!c->metric_initialized,metric_epoch=c->metric_epoch;
  uint64_t actual_fingerprint=0,fixed_fingerprint=0;
  if(!status&&(!c->anderson_metric||c->metric_cells!=(size_t)r->cells))status=1;
  if(!status&&new_metric){memcpy(c->anderson_metric,potential,(size_t)r->cells*sizeof(*potential));c->metric_epoch=0;c->metric_initialized=1;}
  if(!status&&!new_metric){
   for(int i=0;i<r->cells;i++)if(c->anderson_metric[i].dependent!=potential[i].dependent)reset_cells++;
   if(reset_cells){memcpy(c->anderson_metric,potential,(size_t)r->cells*sizeof(*potential));c->metric_epoch++;}
  }
  if(!status){
   metric_epoch=c->metric_epoch;
   actual_fingerprint=coupled_metric_fingerprint(potential,(size_t)r->cells);
   fixed_fingerprint=coupled_metric_fingerprint(c->anderson_metric,(size_t)r->cells);
   char source[1024],destination[1024];
   snprintf(source,sizeof source,"%s/potential_norms.bin",metadata);
   snprintf(destination,sizeof destination,"%s/actual_metric_%02d.bin",c->folder,iteration);
   status=copy_file(source,destination);
   if(!status&&(new_metric||reset_cells)){
    snprintf(destination,sizeof destination,"%s/anderson_metric_epoch_%02d.bin",c->folder,c->metric_epoch);
    status=copy_file(source,destination);
   }
  }
  if(!status)status=full_moment_coupled_residual(r->cells,c->gamma,base->state,base->mass,c->base_area,trial->state,mapped->state,trial->Psi,mapped->Psi,potential,residual,&failed);
  if(!status){
   coupled_metric_summary(potential,(size_t)r->cells,residual);
   residual->metric_epoch=metric_epoch;residual->metric_reset_cells=reset_cells;
   residual->actual_metric_fingerprint=actual_fingerprint;
   residual->anderson_metric_fingerprint=fixed_fingerprint;
  }
  if(!status)status=full_moment_coupled_residual_vector(r->cells,c->gamma,base->state,base->mass,c->base_area,trial->state,mapped->state,trial->Psi,mapped->Psi,c->anderson_metric,solver_residual,solver_count,&failed);
  if(!status){
   long double norm2=0.;for(size_t k=0;k<solver_count;k++)norm2+=solver_residual[k]*solver_residual[k];
   if(!isfinite(norm2))status=1;else residual->vector_norm2=norm2;
  }
  free(potential);
 }
 MPI_Bcast(residual,sizeof *residual,MPI_BYTE,0,r->comm);
 MPI_Bcast(solver_residual,(int)solver_count,MPI_LONG_DOUBLE,0,r->comm);
 if(coupled_residual_record(c,iteration,status,residual,failed))return 1;
 return collective(r,status);
}
static int coupled_finalize(void *context,const FullMomentState *base,
 const FullMomentState *trial,const FullMomentState *mapped,
 const FullMomentCoupledResidual *residual,int iteration){
 (void)trial;(void)residual;(void)iteration;
 CoupledRuntime *c=context;Runtime *r=c->r;int status=0;
 if(!r->rank){
  r->request.outer=0;
  char base_path[512];snprintf(base_path,sizeof base_path,"%s/states/e%llu",r->memory,(unsigned long long)base->header.epoch);
  status=path(r->request.fits_path,base_path)|path(r->request.bounds_path,c->base_metadata)|path(r->request.aux_path,r->initial)|path(r->request.output_dir,c->provider_capture);
  FullMomentProviderResponse response;if(!status)status=rpc(r,FMP_MIDPOINT_OBSERVE,&response);
 }
 if(collective(r,status))return 1;
 /* Preserve the small actual owned factor records before provider RELEASE.
  * Incoming columns have already been audited; no giant face dump is kept. */
 const char *parts[]={"header","deck","factor","incoming_offsets"};
 for(int k=0;k<4;k++){
  char source[1024],destination[1024];snprintf(source,sizeof source,"%s/factor_%llu/%s.bin",r->root,(unsigned long long)mapped->header.epoch,parts[k]);
  snprintf(destination,sizeof destination,"%s/factor_%s_rank%d.bin",c->folder,parts[k],r->rank);
  status|=copy_file(source,destination);
 }
 if(collective(r,status)||release(r))return 1;
 /* Geometry-only proof of the midpoint defined by the VERIFIED map output.
  * It performs no provider/RHS and is never configured as a physical stage. */
 uint64_t proof_epoch=++r->allocated_epoch;char proof_path[512],proof_geometry[512]={0};
 if(coupled_midpoint_write(c,base,mapped,proof_epoch,proof_path)||geometry(r,proof_epoch,proof_path,proof_geometry))return 1;
 if(!r->rank){
  r->request.epoch=proof_epoch;r->request.outer=1;
  char folder[512];snprintf(folder,sizeof folder,"%s/agreement",c->folder);
  status=path(r->request.state_path,proof_path)|path(r->request.geometry_path,proof_geometry)|path(r->request.bounds_path,c->provider_capture)|path(r->request.fits_path,r->initial)|path(r->request.aux_path,r->physics)|path(r->request.output_dir,folder);
  FullMomentProviderResponse response;if(!status)status=rpc(r,FMP_MIDPOINT_OBSERVE,&response);
 }
 if(collective(r,status))return 1;
 char mapped_geometry[512]={0};
 if(geometry(r,mapped->header.epoch,c->mapped_path,mapped_geometry)||configure(r,mapped->header.epoch,c->mapped_path,mapped_geometry,1))return 1;
 char label[96];snprintf(label,sizeof label,"candidate_audit_%04d",c->step);
 if(admissibility(r,c->step,1,label))return 1;
 if(!r->rank){
  char folder[512];snprintf(folder,sizeof folder,"%s/accepted_journal",c->folder);
  status=path(r->request.bounds_path,c->provider_capture)|path(r->request.output_dir,folder);
  memcpy(r->request.reserved,&r->h,sizeof r->h);FullMomentProviderResponse response;
  if(!status)status=rpc(r,FMP_COUPLED_JOURNAL,&response);
 }
 return collective(r,status);
}
static int coupled_count(void *context){return ((CoupledRuntime*)context)->r->providers;}
static double coupled_elapsed(void *context){return MPI_Wtime()-((CoupledRuntime*)context)->started;}
static int coupled_step(Runtime *r,const FullMomentState *base,int step,double started,FullMomentState *out){
 CoupledRuntime c={0};c.r=r;c.step=step;c.started=started;int status=0;
 snprintf(c.folder,sizeof c.folder,"%s/coupled_%04d",r->output,step);
 snprintf(c.base_metadata,sizeof c.base_metadata,"%s/base_metadata",c.folder);
 snprintf(c.provider_capture,sizeof c.provider_capture,"%s/provider",c.folder);
 if(!r->rank)status=mkdir(c.folder,0700)!=0||mkdir(c.base_metadata,0700)!=0;
 if(collective(r,status))return 1;
 char metadata[512]={0};if(coupled_metadata(r,metadata))return 1;
 if(!r->rank){
  c.metric_cells=(size_t)r->cells;
  c.base_area=malloc(c.metric_cells*sizeof(double));
  c.anderson_metric=malloc(c.metric_cells*sizeof(*c.anderson_metric));
  status=!c.base_area||!c.anderson_metric;
  if(!status)status=coupled_read(metadata,"areas",c.base_area,(size_t)r->cells*sizeof(double))||coupled_read(metadata,"gamma",&c.gamma,sizeof c.gamma);
  const char *parts[]={"areas","gamma","ids"};
  for(int k=0;k<3&&!status;k++){
   char source[1024],destination[1024];snprintf(source,sizeof source,"%s/%s.bin",metadata,parts[k]);snprintf(destination,sizeof destination,"%s/%s.bin",c.base_metadata,parts[k]);status=copy_file(source,destination);
  }
 }
 if(collective(r,status)){free(c.base_area);free(c.anderson_metric);return 1;}
 FullMomentCoupledCallbacks callbacks={0};callbacks.context=&c;callbacks.evaluate=coupled_evaluate;callbacks.finalize=coupled_finalize;
 callbacks.type2_fit=coupled_type2_fit;
 callbacks.collective_status=coupled_collective_status;
 callbacks.provider_count=coupled_count;callbacks.elapsed_seconds=coupled_elapsed;
 callbacks.maximum_step_providers=10;callbacks.maximum_total_providers=r->h==.005?50:100;
 /* Fixed disjoint 12h branch allowances enforce the combined 24h ceiling. */
 callbacks.wall_limit_seconds=12*3600.;FullMomentCoupledControlReport result={0};
 int code=full_moment_coupled_control(base,&callbacks,out,&result);
 if(!r->rank){
  char file[1024];snprintf(file,sizeof file,"%s/control.json",c.folder);FILE *f=fopen(file,"wx");
  if(!f)status=1;else{
   fprintf(f,"{\"status_code\":%d,\"iterations\":%d,\"converged\":%s,\"final_checks_passed\":%s,\"provider_calls\":%d,\"branch_provider_total\":%d,\"compute_elapsed_seconds\":%.17g,\"anderson_depth\":[",code,result.iterations,result.converged?"true":"false",result.final_checks_passed?"true":"false",result.providers_after-result.providers_before,r->providers,coupled_elapsed(&c));
   for(int i=0;i<result.iterations;i++)fprintf(f,"%s%d",i?",":"",result.anderson_depth[i]);
   fprintf(f,"],\"anderson_qr_rank\":[");for(int i=0;i<result.iterations;i++)fprintf(f,"%s%d",i?",":"",result.anderson_qr_rank[i]);
   fprintf(f,"],\"anderson_gamma\":[");for(int i=0;i<result.iterations;i++)fprintf(f,"%s[%.17g,%.17g,%.17g]",i?",":"",result.anderson_gamma[i][0],result.anderson_gamma[i][1],result.anderson_gamma[i][2]);
   fprintf(f,"],\"anderson_lambda\":[");for(int i=0;i<result.iterations;i++)fprintf(f,"%s%.17g",i?",":"",result.anderson_lambda[i]);
   fprintf(f,"],\"anderson_condition\":[");for(int i=0;i<result.iterations;i++){if(i)fputc(',',f);if(isfinite(result.anderson_condition[i]))fprintf(f,"%.17g",result.anderson_condition[i]);else fprintf(f,"null");}
   fprintf(f,"],\"anderson_predicted_norm2\":[");for(int i=0;i<result.iterations;i++)fprintf(f,"%s\"%.21Lg\"",i?",":"",result.anderson_predicted_norm2[i]);
   fprintf(f,"],\"prediction_miss_streak\":[");for(int i=0;i<result.iterations;i++)fprintf(f,"%s%d",i?",":"",result.prediction_miss_streak[i]);
   fprintf(f,"],\"anderson_fallback\":[");for(int i=0;i<result.iterations;i++)fprintf(f,"%s%d",i?",":"",result.anderson_fallback[i]);
   fprintf(f,"],\"solver_history\":[");
   for(int i=0;i<result.iterations;i++){
    const FullMomentCoupledResidual *q=&result.iteration_residual[i];
    fprintf(f,"%s{\"iteration\":%d,\"fixed_metric_norm2\":\"%.21Lg\",\"actual_maximum\":[",i?",":"",i,q->vector_norm2);
    for(int k=0;k<FM_RES_COUNT;k++)fprintf(f,"%s\"%.21Lg\"",k?",":"",q->maximum[k]);
    fprintf(f,"],\"raw_psi_remainder_numerator_max\":\"%.21Lg\",\"dependent_cells\":%d,\"sigma_min\":\"%.21Lg\",\"sigma_max\":\"%.21Lg\",\"rank_margin_min\":\"%.21Lg\",\"metric_epoch\":%d,\"metric_reset_cells\":%d,\"actual_metric_fingerprint\":\"%016llx\",\"anderson_metric_fingerprint\":\"%016llx\",\"qr_depth\":%d,\"qr_rank\":%d,\"qr_condition\":",q->psi_remainder_numerator_max,q->dependent_cells,q->actual_sigma_min,q->actual_sigma_max,q->actual_rank_margin_min,q->metric_epoch,q->metric_reset_cells,(unsigned long long)q->actual_metric_fingerprint,(unsigned long long)q->anderson_metric_fingerprint,result.anderson_depth[i],result.anderson_qr_rank[i]);
    if(isfinite(result.anderson_condition[i]))fprintf(f,"%.17g",result.anderson_condition[i]);else fprintf(f,"null");
    fprintf(f,",\"gamma\":[%.17g,%.17g,%.17g],\"lambda\":%.17g,\"predicted_norm2\":\"%.21Lg\",\"prediction_miss_streak\":%d,\"fallback_reason\":%d}",result.anderson_gamma[i][0],result.anderson_gamma[i][1],result.anderson_gamma[i][2],result.anderson_lambda[i],result.anderson_predicted_norm2[i],result.prediction_miss_streak[i],result.anderson_fallback[i]);
   }
   fprintf(f,"],\"accepted\":false,\"outer_native_commit_required\":true}\n");
   status=fflush(f)!=0;if(!status)status=fsync(fileno(f))!=0;status|=fclose(f)!=0;
  }
 }
 full_moment_state_close(&c.midpoint);free(c.base_area);free(c.anderson_metric);
 return collective(r,status||code);
}
