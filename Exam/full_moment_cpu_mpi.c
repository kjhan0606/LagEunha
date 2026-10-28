/* Native main-thread MPI control of persistent, spawned rank-local providers.
 * This gate verifies owned CPU preparation, not an accepted evolution step.
 */
#define _POSIX_C_SOURCE 200809L
#include "full_moment_cpu_mpi.h"
#include "full_moment_provider.h"
#include "full_moment_native_stage.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>

static int all_status(MPI_Comm c,int status){int out;MPI_Allreduce(&status,&out,1,MPI_INT,MPI_MAX,c);return out;}
static int set_path(char *dest,const char *source){size_t n=strlen(source);if(n>=512)return 1;memcpy(dest,source,n+1);return 0;}
static char *read_text(const char *path,int *length){
 FILE *f=fopen(path,"rb");if(!f)return NULL;
 if(fseek(f,0,SEEK_END)){fclose(f);return NULL;}long n=ftell(f);rewind(f);
 if(n<2||n>INT_MAX){fclose(f);return NULL;}char *p=malloc(n+1);
 if(!p||fread(p,1,n,f)!=(size_t)n){free(p);fclose(f);return NULL;}
 fclose(f);p[n]=0;*length=(int)n;return p;
}
int full_moment_join_rows(MPI_Comm c,int rank,int size,const char *local,const char *global){
 int count=0;char *data=read_text(local,&count);int status=!data;
 if(!status)status=data[0]!='['||data[count-1]!=']';
 if(all_status(c,status)){free(data);return 1;}
 int *counts=malloc(size*sizeof(int)),*offsets=malloc(size*sizeof(int));
 if(all_status(c,!counts||!offsets)){free(data);free(counts);free(offsets);return 1;}
 count-=2;MPI_Allgather(&count,1,MPI_INT,counts,1,MPI_INT,c);
 int64_t total=0;for(int r=0;r<size;r++){offsets[r]=(int)total;total+=counts[r];}
 status=total>INT_MAX;char *out=rank?NULL:malloc((size_t)total+1);
 if(!rank&&!out)status=1;
 if(all_status(c,status)){free(data);free(counts);free(offsets);free(out);return 1;}
 MPI_Gatherv(data+1,count,MPI_CHAR,out,counts,offsets,MPI_CHAR,0,c);
 if(!rank){FILE *f=fopen(global,"wx");status=!f;if(f){
  status|=fputc('[',f)==EOF;
  for(int r=0;r<size;r++){
   if(r)status|=fputc(',',f)==EOF;
   status|=fwrite(out+offsets[r],1,counts[r],f)!=(size_t)counts[r];
  }
  status|=fputc(']',f)==EOF;status|=fclose(f)!=0;
 }}
 free(data);free(counts);free(offsets);free(out);return all_status(c,status);
}
int full_moment_cpu_mpi_validation(const char *path,MPI_Comm c){
 int rank,size;MPI_Comm_rank(c,&rank);MPI_Comm_size(c,&size);
 char python[512]="",script[512]="",state[512]="",geometry[512]="",output[512]="",physics_file[512]="",physics_sha[65]="";
 char line[4096];int workers=7,ncells=0,q48=0,status=0,actual_stage=0,threads=8;
 FILE *f=fopen(path,"r");if(!f)status=1;
 if(f){while(fgets(line,sizeof line,f)){
  sscanf(line,"moment_python = %511s",python);sscanf(line,"moment_provider_script = %511s",script);
  sscanf(line,"moment_state = %511s",state);sscanf(line,"moment_geometry = %511s",geometry);
  sscanf(line,"moment_output = %511s",output);sscanf(line,"moment_physics_file = %511s",physics_file);
  sscanf(line,"moment_physics_sha = %64s",physics_sha);sscanf(line,"moment_workers = %d",&workers);
  sscanf(line,"moment_cells = %d",&ncells);sscanf(line,"moment_q48 = %d",&q48);
  sscanf(line,"moment_actual_stage_validation = %d",&actual_stage);sscanf(line,"moment_threads = %d",&threads);
 }fclose(f);}
 status|=!python[0]||!script[0]||!state[0]||!geometry[0]||!output[0]||!physics_file[0]||strlen(physics_sha)!=64||workers<1||ncells<size;
 /* Initial bulk-buffer transport is deliberately limited to one host. */
 char host[256]={0};gethostname(host,sizeof host);char *hosts=calloc(size,256);
 if(all_status(c,status||!hosts)){free(hosts);return 90;}
 MPI_Allgather(host,256,MPI_CHAR,hosts,256,MPI_CHAR,c);
 for(int r=1;r<size;r++)status|=strcmp(hosts,hosts+256*r)!=0;free(hosts);
 if(all_status(c,status||threads<1))return 90;
 uint64_t token=rank?0:(((uint64_t)time(NULL)<<24)^(uint64_t)getpid());MPI_Bcast(&token,1,MPI_UINT64_T,0,c);
 char root[512],root0[512],log[1024],report[1024],local[1024],bounds[512],fits[512];
 snprintf(root,sizeof root,"/dev/shm/eunha_fmp_%llu_r%d",(unsigned long long)token,rank);
 snprintf(root0,sizeof root0,"/dev/shm/eunha_fmp_%llu_r0",(unsigned long long)token);
 snprintf(bounds,sizeof bounds,"%s/global_bounds.json",root0);snprintf(fits,sizeof fits,"%s/global_fits.json",root0);
 snprintf(log,sizeof log,"%s/provider_rank%d.log",output,rank);
 FullMomentProvider service={0};FullMomentProviderRequest req={0};FullMomentProviderResponse res={0};
 status=full_moment_provider_start(&service,python,script,workers,log);
 if(all_status(c,status)){if(service.pid>0)full_moment_provider_close(&service);return 91;}
 req.epoch=1;req.run_token=token;req.rank=rank;req.nranks=size;req.ncells=ncells;req.outer=q48;
 memcpy(req.physics_sha,physics_sha,64);
 status=set_path(req.state_path,state)|set_path(req.geometry_path,geometry)|set_path(req.output_dir,root)|set_path(req.aux_path,physics_file);
 double start=MPI_Wtime();req.command=FMP_CONFIGURE;
 if(!status)status=full_moment_provider_call(&service,&req,&res);
 if(all_status(c,status)){full_moment_provider_close(&service);return 92;}
 double phase[3]={0};uint64_t gn=0;
 for(int command=FMP_BOUNDS;command<=FMP_VELOCITY;command++){
  double begin=MPI_Wtime();req.command=command;
  if(command>=FMP_DENSITY)status|=set_path(req.bounds_path,bounds);
  if(command==FMP_VELOCITY)status|=set_path(req.fits_path,fits);
  if(!status)status=full_moment_provider_call(&service,&req,&res);
  if(all_status(c,status)){full_moment_provider_close(&service);return 93;}
  snprintf(local,sizeof local,"%s/rows_%d.json",root,command);
  if(command<FMP_VELOCITY)status=full_moment_join_rows(c,rank,size,local,command==FMP_BOUNDS?bounds:fits);
  if(all_status(c,status)){full_moment_provider_close(&service);return 94;}
  double elapsed=MPI_Wtime()-begin;MPI_Allreduce(&elapsed,phase+command-FMP_BOUNDS,1,MPI_DOUBLE,MPI_MAX,c);
  if(command==FMP_VELOCITY)gn=res.gravity_nodes;
  printf("[FULL_MOMENT_CPU_OWNED] rank=%d ranks=%d command=%d owned_start=%llu owned_cells=%llu gravity_nodes=%llu max_seconds=%.9g\n",rank,size,command,(unsigned long long)res.owned_start,(unsigned long long)res.owned_count,(unsigned long long)res.gravity_nodes,phase[command-FMP_BOUNDS]);fflush(stdout);
 }
 double elapsed=MPI_Wtime()-start,max_elapsed=0;MPI_Allreduce(&elapsed,&max_elapsed,1,MPI_DOUBLE,MPI_MAX,c);
 snprintf(report,sizeof report,"%s/provider_rank%d.json",output,rank);f=fopen(report,"wx");
 status=!f;if(f){fprintf(f,"{\"status\":\"passed_owned_cpu_MPI_preparation\",\"rank\":%d,\"ranks\":%d,\"cells\":%d,\"workers\":%d,\"host\":\"%s\",\"buffer_root\":\"%s\",\"owned_start\":%llu,\"owned_count\":%llu,\"gravity_nodes\":%llu,\"bounds_seconds\":%.17g,\"density_seconds\":%.17g,\"velocity_seconds\":%.17g,\"max_wall_seconds\":%.17g,\"evolution_steps\":0}\n",rank,size,ncells,workers,host,root,(unsigned long long)res.owned_start,(unsigned long long)res.owned_count,(unsigned long long)gn,phase[0],phase[1],phase[2],max_elapsed);status|=fclose(f)!=0;}
 if(all_status(c,status)){full_moment_provider_close(&service);return 95;}
 if(actual_stage){
  char velocity[512],native_root[512]={0};snprintf(velocity,sizeof velocity,"%s/global_velocity.json",root0);
  status=full_moment_join_rows(c,rank,size,local,velocity);
  if(all_status(c,status)){full_moment_provider_close(&service);return 96;}
  if(!rank){
   req.command=FMP_FACE_STAGE;status=set_path(req.aux_path,velocity);
   if(!status)status=full_moment_provider_call(&service,&req,&res);
   if(status==100+FMP_NATIVE_STAGE_REQUIRED){memcpy(native_root,res.output_path,512);status=0;}else status=1;
  }
  if(all_status(c,status)){full_moment_provider_close(&service);return 97;}
  MPI_Bcast(native_root,512,MPI_CHAR,0,c);
  status=full_moment_native_stage(c,native_root,root,req.epoch,req.run_token,threads);
  if(all_status(c,status)){full_moment_provider_close(&service);return 98;}
  if(!rank){
   char native_outputs[1024];snprintf(native_outputs,sizeof native_outputs,"%s/outputs",native_root);
   req.command=FMP_RESUME_STAGE;status=set_path(req.output_dir,native_outputs);
   if(!status)status=full_moment_provider_call(&service,&req,&res);
   if(!status){
    char *data=NULL;int length=0;data=read_text(res.output_path,&length);
    char final[1024];snprintf(final,sizeof final,"%s/actual_stage_report.json",output);
    FILE *report_file=data?fopen(final,"wx"):NULL;status=!report_file;
    if(report_file){status=fwrite(data,1,length,report_file)!=(size_t)length;status|=fclose(report_file)!=0;}
    free(data);
   }
  }
  if(all_status(c,status)){full_moment_provider_close(&service);return 99;}
 }
 req.command=FMP_STOP;status|=full_moment_provider_call(&service,&req,&res);
 status|=full_moment_provider_close(&service);
 return all_status(c,status)?95:0;
}
