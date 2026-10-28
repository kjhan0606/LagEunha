/* Authoritative midpoint control for the original MPI executable.
 * Providers supply geometry/reconstruction/root/audits only. Native ranks
 * own stage kernels, factor maps, x/W/raw/L/Psi updates and collective commit.
 */
#define _POSIX_C_SOURCE 200809L
#include "full_moment_evolution.h"
#include "full_moment_provider.h"
#include "full_moment_provider_mpi.h"
#include "full_moment_cpu_mpi.h"
#include "full_moment_native_stage.h"
#include "full_moment_native_factor.h"
#include "full_moment_state_io.h"
#include "full_moment_evolution_budget.h"
#include "full_moment_coupled_control.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/statvfs.h>

typedef struct {
 MPI_Comm comm;int rank,ranks,cells,workers,threads,steps,validation,providers;
 double h,coupled_max_time;uint64_t token,epoch,allocated_epoch;
 char python[512],script[512],input[512],initial[512],initial_geometry[512],output[512],physics[512],physics_sha[65];
 char restart_manifest[512];
 char fork_manifest[512],integrator[128];int coupled;
 char memory[512],state_path[512],root[512],root0[512],geometry[512],bounds[512],fits[512],native_stage[512],rates[512];
 FullMomentProvider service;FullMomentProviderRequest request;
} Runtime;
static int collective(Runtime *r,int status){int out;MPI_Allreduce(&status,&out,1,MPI_INT,MPI_MAX,r->comm);return out;}
static int path(char *out,const char *value){size_t n=strlen(value);if(n>=512)return 1;memcpy(out,value,n+1);return 0;}
static int rpc(Runtime *r,int command,FullMomentProviderResponse *response){r->request.command=command;return full_moment_provider_call(&r->service,&r->request,response);}
static int copy_file(const char *source,const char *destination){
 FILE *in=fopen(source,"rb");if(!in)return 1;FILE *out=fopen(destination,"wbx");if(!out){fclose(in);return 1;}
 char buffer[16384];size_t count;int status=0;
 while((count=fread(buffer,1,sizeof buffer,in))>0)if(fwrite(buffer,1,count,out)!=count){status=1;break;}
 status|=ferror(in);status|=fclose(in)!=0;status|=fflush(out)!=0;if(!status)status=fsync(fileno(out))!=0;status|=fclose(out)!=0;return status;
}
static void preserve_failure(Runtime *r){
 if(r->rank)return;char folder[1024];snprintf(folder,sizeof folder,"%s/failure",r->output);
 if(mkdir(folder,0700))return;
 int input_status=1;struct stat st;
 if(!stat(r->request.state_path,&st)){
  char target[1200];
  if(S_ISDIR(st.st_mode)){
   FullMomentState failed={0};input_status=full_moment_state_load(r->request.state_path,&failed);
   if(!input_status){snprintf(target,sizeof target,"%s/input_state",folder);input_status=full_moment_state_publish(target,&failed);full_moment_state_close(&failed);}
  }else if(S_ISREG(st.st_mode)){
   snprintf(target,sizeof target,"%s/input_state.npz",folder);input_status=copy_file(r->request.state_path,target);
  }
 }
 int saved_provider_state=0,saved_node_witness=0;
 const char *names[]={"failure_state.npz","bounded_node_witness.npz"};
 for(int k=0;k<2;k++){
  char source[1200],target[1200];snprintf(source,sizeof source,"%s/actual_stage_audit/%s",r->root,names[k]);snprintf(target,sizeof target,"%s/%s",folder,names[k]);
  if(!stat(source,&st)&&S_ISREG(st.st_mode)&&st.st_size<8*1024*1024){int copied=copy_file(source,target)==0;if(k)saved_node_witness=copied;else saved_provider_state=copied;}
 }
 char file[1200];snprintf(file,sizeof file,"%s/manifest.json",folder);FILE *f=fopen(file,"wx");
 if(f){fprintf(f,"{\"status\":\"failed_native_transaction_preserved\",\"input_copy_status\":%d,\"provider_failure_state_saved\":%s,\"bounded_node_witness_saved\":%s,\"epoch\":%llu,\"request_epoch\":%llu,\"step\":%u,\"stage\":%u,\"last_request\":%u,\"provider_attempts\":%d,\"h\":%.17g,\"run_token\":%llu,\"source_state\":\"%s\",\"source_geometry\":\"%s\",\"node_local_root\":\"%s\",\"physics_sha\":\"%s\",\"previous_accepted_preserved\":true,\"retry\":false,\"timestep_reduced\":false}\n",input_status,saved_provider_state?"true":"false",saved_node_witness?"true":"false",(unsigned long long)r->epoch,(unsigned long long)r->request.epoch,r->request.step,r->request.stage,r->request.command,r->providers,r->h,(unsigned long long)r->token,r->request.state_path,r->geometry,r->memory,r->physics_sha);
  fflush(f);fsync(fileno(f));fclose(f);
 }
 fprintf(stderr,"[FULL_MOMENT_FAILURE_PRESERVED] %s input_status=%d\n",folder,input_status);
}
static int save_factor_audit(Runtime *r,uint64_t epoch){
 char source[1024],destination[1024];snprintf(source,sizeof source,"%s/factor_%llu/audit.json",r->root,(unsigned long long)epoch);
 snprintf(destination,sizeof destination,"%s/factor_%llu_rank%d.json",r->output,(unsigned long long)epoch,r->rank);
 return collective(r,copy_file(source,destination));
}
static int configure(Runtime *r,uint64_t epoch,const char *state,const char *geometry,int independent_q48){
 r->epoch=epoch;snprintf(r->root,sizeof r->root,"/dev/shm/eunha_fmp_%llu_e%llu_r%d",(unsigned long long)r->token,(unsigned long long)epoch,r->rank);
 snprintf(r->root0,sizeof r->root0,"/dev/shm/eunha_fmp_%llu_e%llu_r0",(unsigned long long)r->token,(unsigned long long)epoch);
 snprintf(r->bounds,sizeof r->bounds,"%s/global_bounds.json",r->root0);snprintf(r->fits,sizeof r->fits,"%s/global_fits.json",r->root0);
 int status=path(r->state_path,state)|path(r->geometry,geometry);
 FullMomentProviderRequest *q=&r->request;q->epoch=epoch;q->run_token=r->token;q->rank=r->rank;q->nranks=r->ranks;q->ncells=r->cells;q->outer=independent_q48;
 memcpy(q->physics_sha,r->physics_sha,64);
 status|=path(q->state_path,state)|path(q->geometry_path,geometry)|path(q->output_dir,r->root)|path(q->aux_path,r->physics);
 FullMomentProviderResponse response;if(!status)status=rpc(r,FMP_CONFIGURE,&response);
 if(collective(r,status))return 1;
 for(int command=FMP_BOUNDS;command<=FMP_DENSITY;command++){
  if(command==FMP_DENSITY)status|=path(q->bounds_path,r->bounds);
  if(!status)status=rpc(r,command,&response);
  if(collective(r,status))return 1;
  char local[1024];snprintf(local,sizeof local,"%s/rows_%d.json",r->root,command);
  if(full_moment_join_rows(r->comm,r->rank,r->ranks,local,command==FMP_BOUNDS?r->bounds:r->fits))return 1;
 }
 return 0;
}
static int geometry(Runtime *r,uint64_t epoch,const char *state,char *out){
 int status=0;if(!r->rank){
  char destination[512];snprintf(destination,sizeof destination,"%s/geometry/e%llu",r->memory,(unsigned long long)epoch);
  r->request.epoch=epoch;status=path(r->request.state_path,state)|path(r->request.output_dir,destination);
  FullMomentProviderResponse response;if(!status)status=rpc(r,FMP_GEOMETRY,&response);
  if(!status)status=path(out,response.output_path);
 }
 if(collective(r,status))return 1;MPI_Bcast(out,512,MPI_CHAR,0,r->comm);return 0;
}
static int release(Runtime *r){
 r->request.epoch=r->epoch;
 return full_moment_collective_release(r->comm,&r->service,&r->request,r->geometry);
}
static int rhs(Runtime *r){
 if(r->providers>=400)return 1;
 double start=MPI_Wtime();FullMomentProviderResponse response;int status=path(r->request.fits_path,r->fits)|path(r->request.bounds_path,r->bounds);
 if(!status)status=rpc(r,FMP_VELOCITY,&response);
 if(collective(r,status))return 1;
 char local[1024],global[512];snprintf(local,sizeof local,"%s/rows_4.json",r->root);snprintf(global,sizeof global,"%s/global_velocity.json",r->root0);
 if(full_moment_join_rows(r->comm,r->rank,r->ranks,local,global))return 1;
 r->providers++;
 if(!r->rank){
  status=path(r->request.aux_path,global);
  if(!status)status=rpc(r,FMP_FACE_STAGE,&response);
  if(status==100+FMP_NATIVE_STAGE_REQUIRED){status=path(r->native_stage,response.output_path);}else status=1;
 }
 if(collective(r,status))return 1;MPI_Bcast(r->native_stage,512,MPI_CHAR,0,r->comm);
 if(full_moment_native_stage(r->comm,r->native_stage,r->root,r->epoch,r->token,r->threads))return 1;
 if(!r->rank){
  char outputs[512];snprintf(outputs,sizeof outputs,"%s/outputs",r->native_stage);status=path(r->request.output_dir,outputs);
  if(!status)status=rpc(r,FMP_RESUME_STAGE,&response);
  if(!status){char destination[1024];snprintf(destination,sizeof destination,"%s/rhs_%llu.json",r->output,(unsigned long long)r->epoch);status=copy_file(response.output_path,destination);}
 }
 if(collective(r,status))return 1;
 snprintf(r->rates,sizeof r->rates,"%s/rates",r->root0);
 printf("[FULL_MOMENT_RHS] rank=%d epoch=%llu provider=%d seconds=%.9g accepted=false\n",r->rank,(unsigned long long)r->epoch,r->providers,MPI_Wtime()-start);fflush(stdout);return 0;
}
static int admissibility(Runtime *r,int step,int checkpoint,const char *label){
 int status=0;if(!r->rank){
  char output[512];snprintf(output,sizeof output,"%s/%s",r->output,label);
  r->request.outer=checkpoint;r->request.step=step;
  r->request.reserved[8]=r->coupled!=0;
  status=path(r->request.output_dir,output)|path(r->request.aux_path,r->initial)|path(r->request.fits_path,r->fits);
  FullMomentProviderResponse response;if(!status)status=rpc(r,FMP_ADMISSIBILITY,&response);
 }
 return collective(r,status);
}
static int checkpoint(Runtime *r,int step){
 if(r->validation)return 1;
 if(step==r->steps)return 1;
 double t=step*r->h;double points[]={.01,.025,.05,.1,.25,.5,1.,2.,4.,6.,8.,10.};
 for(size_t i=0;i<sizeof points/sizeof points[0];i++)if(fabs(t-points[i])<32*DBL_EPSILON*fmax(1.,t))return 1;
 return 0;
}
static int report(Runtime *r,int step,const FullMomentState *accepted,double seconds,int final){
 int status=0;if(!r->rank){
  char name[1024];snprintf(name,sizeof name,"%s/%s_%04d.json",r->output,final?"completed":"accepted",step);
  char temporary[1100];snprintf(temporary,sizeof temporary,"%s.pending",name);
  if(access(name,F_OK)==0)status=1;FILE *f=status?NULL:fopen(temporary,"wx");status|=!f;
  if(f){fprintf(f,"{\"status\":\"%s\",\"native_time_control\":true,\"ranks\":%d,\"cells\":%d,\"step\":%d,\"time\":%.17g,\"h\":%.17g,\"physical_providers\":%d,\"step_wall_seconds\":%.17g,\"model\":\"FullMomentKepler\",\"integrator\":\"%s\",\"fixed_mass\":true,\"Re_infinity\":true,\"accepted_commit\":true,\"long_time_stability_verified\":false}\n",final?"completed_native_evolution_batch":r->coupled?"accepted_native_coupled_midpoint_step":"accepted_native_midpoint_step",r->ranks,r->cells,step,accepted->header.time,r->h,r->providers,seconds,r->coupled?"self-consistent midpoint + Sylvester/exponential/positive QR":"explicit midpoint + Sylvester/exponential/positive QR");
   status|=fflush(f)!=0;if(!status)status=fsync(fileno(f))!=0;status|=fclose(f)!=0;
   if(!status)status=rename(temporary,name)!=0;
   if(!status){int directory=open(r->output,O_RDONLY|O_DIRECTORY);if(directory<0)status=1;else{status=fsync(directory)!=0;close(directory);}}
  }
 }
 return collective(r,status);
}
#include "full_moment_coupled_runtime.inc"
int full_moment_evolution_run(const char *config,MPI_Comm comm){
 Runtime r={0};r.comm=comm;MPI_Comm_rank(comm,&r.rank);MPI_Comm_size(comm,&r.ranks);r.cells=19945;r.workers=7;r.threads=8;
 FullMomentState initial={0},accepted={0},midpoint={0},candidate={0};int status=0,code=1;double started=MPI_Wtime();
 FILE *f=fopen(config,"r");char line[4096];status=!f;
 if(f){while(fgets(line,sizeof line,f)){
  sscanf(line,"moment_python = %511s",r.python);sscanf(line,"moment_provider_script = %511s",r.script);
  sscanf(line,"moment_state = %511s",r.input);sscanf(line,"moment_initial = %511s",r.initial);
  sscanf(line,"moment_restart_manifest = %511s",r.restart_manifest);
  sscanf(line,"moment_fork_manifest = %511s",r.fork_manifest);
  sscanf(line,"moment_integrator = %127s",r.integrator);
  sscanf(line,"moment_geometry = %511s",r.initial_geometry);sscanf(line,"moment_output = %511s",r.output);
  sscanf(line,"moment_physics_file = %511s",r.physics);sscanf(line,"moment_physics_sha = %64s",r.physics_sha);
  sscanf(line,"moment_workers = %d",&r.workers);sscanf(line,"moment_threads = %d",&r.threads);
  sscanf(line,"moment_steps = %d",&r.steps);sscanf(line,"moment_validation = %d",&r.validation);sscanf(line,"moment_h = %lf",&r.h);
  sscanf(line,"moment_coupled_max_time = %lf",&r.coupled_max_time);
 }fclose(f);}
 if(!r.initial[0])path(r.initial,r.input);
 r.coupled=!strcmp(r.integrator,"self_consistent_midpoint");
 status|=r.integrator[0]&&!r.coupled&&strcmp(r.integrator,"explicit_midpoint");
 if(r.coupled){
  status|=(!r.fork_manifest[0]&&!r.restart_manifest[0])||(r.fork_manifest[0]&&r.restart_manifest[0]);
  status|=!isfinite(r.coupled_max_time)||r.coupled_max_time<.03-DBL_EPSILON||r.coupled_max_time>10.+DBL_EPSILON;
  status|=(r.fork_manifest[0]&&fabs(r.coupled_max_time-.03)>DBL_EPSILON);
  status|=r.steps*r.h>r.coupled_max_time+DBL_EPSILON;
 }else status|=r.fork_manifest[0]!=0||r.coupled_max_time!=0;
 status|=!r.python[0]||!r.script[0]||!r.input[0]||(!r.initial_geometry[0]&&!r.restart_manifest[0]&&!r.fork_manifest[0])||!r.output[0]||!r.physics[0]||strlen(r.physics_sha)!=64||r.ranks<1||r.ranks>2||r.workers<1||r.threads<1||r.steps<1||r.steps>1000000||!isfinite(r.h)||r.h<=0||r.h>.005;
 if(r.coupled)status|=(r.h!=.005&&r.h!=.0025);
 struct stat input_stat;if(stat(r.input,&input_stat))status=1;
 else if(r.restart_manifest[0])status|=!S_ISDIR(input_stat.st_mode)||!strcmp(r.input,r.initial);
 else if(r.fork_manifest[0])status|=!S_ISREG(input_stat.st_mode)||!strcmp(r.input,r.initial);
 else status|=!S_ISREG(input_stat.st_mode)||strcmp(r.input,r.initial);
 char host[256]={0};gethostname(host,sizeof host);char *hosts=calloc(r.ranks,256);
 if(collective(&r,status||!hosts)){free(hosts);return 111;}
 MPI_Allgather(host,256,MPI_CHAR,hosts,256,MPI_CHAR,comm);
 for(int k=1;k<r.ranks;k++)status|=strcmp(hosts,hosts+256*k)!=0;free(hosts);
 if(collective(&r,status))return 112;
 if(!r.rank)r.token=((uint64_t)time(NULL)<<24)^(uint64_t)getpid();MPI_Bcast(&r.token,1,MPI_UINT64_T,0,comm);
 snprintf(r.memory,sizeof r.memory,"/dev/shm/eunha_run_%llu",(unsigned long long)r.token);
 if(!r.rank){
  status=mkdir(r.memory,0700)!=0;char states[1024],geo[1024];snprintf(states,sizeof states,"%s/states",r.memory);snprintf(geo,sizeof geo,"%s/geometry",r.memory);
  status|=mkdir(states,0700)!=0||mkdir(geo,0700)!=0;
  struct statvfs space;if(statvfs("/dev/shm",&space))status=1;else status|=(uint64_t)space.f_bavail*space.f_frsize<UINT64_C(32)*1024*1024*1024;
 }
 if(collective(&r,status))return 113;
 char log[1024];snprintf(log,sizeof log,"%s/provider_rank%d.log",r.output,r.rank);
 if(r.coupled){setenv("LAGFORCE_NATIVE_FACE_STAGE","1",1);setenv("LAGFORCE_COUPLED_JOURNAL","1",1);}
 status=full_moment_provider_start(&r.service,r.python,r.script,r.workers,log);
 if(collective(&r,status))goto close;
 printf("[FULL_MOMENT_NATIVE_RUNTIME] rank=%d ranks=%d host=%s pid=%ld h=%.17g steps=%d binary_sha=%s\n",r.rank,r.ranks,host,(long)getpid(),r.h,r.steps,getenv("EUNHA_MOMENT_BINARY_SHA")?getenv("EUNHA_MOMENT_BINARY_SHA"):"unrecorded");fflush(stdout);
 /* Restart payload never passes CONFIGURE with the previous run's token.
  * Verify its model/baseline/audits first, then native-copy into a new epoch.
  * Geometry is regenerated from the imported moving sites and weights. */
 char base_path[512];snprintf(base_path,sizeof base_path,"%s/states/e1",r.memory);
 if(r.restart_manifest[0]||r.fork_manifest[0]){
  r.request.epoch=1;r.request.run_token=r.token;r.request.rank=r.rank;r.request.nranks=r.ranks;r.request.ncells=r.cells;
  memcpy(r.request.physics_sha,r.physics_sha,64);memcpy(r.request.reserved,&r.h,sizeof r.h);
  if(!r.rank){
   r.request.step=r.steps;r.request.outer=r.coupled;
   memcpy(r.request.reserved,&r.h,sizeof r.h);
   memcpy(r.request.reserved+8,&r.coupled_max_time,sizeof r.coupled_max_time);
   int is_fork=r.fork_manifest[0]!=0;
   char verified[512];snprintf(verified,sizeof verified,"%s/%s",r.memory,is_fork?"fork_verified":"restart_verified");
   status=path(r.request.state_path,r.input)|path(r.request.bounds_path,is_fork?r.fork_manifest:r.restart_manifest)|path(r.request.fits_path,r.initial)|path(r.request.aux_path,r.physics)|path(r.request.output_dir,verified);
   FullMomentProviderResponse response;if(!status)status=rpc(&r,is_fork?FMP_VERIFY_FORK:FMP_VERIFY_RESTART,&response);
   FullMomentStateHeader expected={0};
   if(!status){
    char file[1024];snprintf(file,sizeof file,"%s/expected_header.bin",response.output_path);FILE *header=fopen(file,"rb");
    if(!header)status=1;else{status=fread(&expected,1,sizeof expected,header)!=sizeof expected;status|=fgetc(header)!=EOF;fclose(header);}
   }
   char import_source[1024];snprintf(import_source,sizeof import_source,"%s/serialized_source",response.output_path);
   if(!status)status=full_moment_state_import(r.coupled?import_source:r.input,base_path,&expected,r.token,&initial);
   if(!status){
    char source[1024],destination[1024];snprintf(source,sizeof source,"%s/receipt.json",response.output_path);snprintf(destination,sizeof destination,"%s/%s_import.json",r.output,r.coupled?"fork":"restart");status=copy_file(source,destination);
   }
  }
  if(collective(&r,status))goto close;
  char imported_geometry[512]={0};if(geometry(&r,1,base_path,imported_geometry))goto close;
  if(configure(&r,1,base_path,imported_geometry,r.fork_manifest[0]?0:1))goto close;
 }else{
 if(configure(&r,1,r.input,r.initial_geometry,1))goto close;
 char serialized[512]={0};if(!r.rank){FullMomentProviderResponse response;status=rpc(&r,FMP_EXPORT_STATE,&response);if(!status)status=path(serialized,response.output_path);}
 if(collective(&r,status))goto close;MPI_Bcast(serialized,512,MPI_CHAR,0,comm);
 status=full_moment_state_load(serialized,&initial);if(collective(&r,status))goto close;
 /* Keep a native-owned initial/base path independent of provider buffers. */
 if(!r.rank)status=full_moment_state_write(base_path,&initial);
 if(collective(&r,status))goto close;
 }
 status=full_moment_state_load(base_path,&accepted);if(collective(&r,status))goto close;
 int start_step=0;
 status=full_moment_evolution_budget(accepted.header.time,r.h,r.steps,&start_step);
 if(collective(&r,status))goto close;
 if(r.restart_manifest[0]){
  if(admissibility(&r,start_step,1,"restart_audit"))goto close;
  printf("[FULL_MOMENT_IMPORTED] rank=%d time=%.17g start_step=%d; payload unchanged, new epoch/token, no accepted time step\n",r.rank,accepted.header.time,start_step);fflush(stdout);
 }
 r.allocated_epoch=r.epoch;
 for(int step=start_step+1;step<=r.steps;step++){
  double step_start=MPI_Wtime();r.request.step=step;r.request.stage=0;
  if(r.coupled){
   if(coupled_step(&r,&accepted,step,started,&candidate))goto close;
  }else{
  if(rhs(&r))goto close;
  uint64_t half_epoch=r.epoch+1;char half_path[512];snprintf(half_path,sizeof half_path,"%s/states/e%llu",r.memory,(unsigned long long)half_epoch);
  if(full_moment_native_factor(comm,&r.service,&r.request,&accepted,&accepted,r.native_stage,r.root,r.rates,half_path,half_epoch,r.h/2,(step-.5)*r.h,r.threads,r.validation))goto close;
  if(save_factor_audit(&r,half_epoch))goto close;
  if(release(&r))goto close;
  char half_geometry[512]={0};if(geometry(&r,half_epoch,half_path,half_geometry))goto close;
  if(configure(&r,half_epoch,half_path,half_geometry,0))goto close;
  char half_label[64];snprintf(half_label,sizeof half_label,"half_audit_%04d",step);
  if(admissibility(&r,step,0,half_label))goto close;
  status=full_moment_state_load(half_path,&midpoint);if(collective(&r,status))goto close;
  r.request.stage=1;if(rhs(&r))goto close;
  uint64_t full_epoch=r.epoch+1;char full_path[512];snprintf(full_path,sizeof full_path,"%s/states/e%llu",r.memory,(unsigned long long)full_epoch);
  int check=checkpoint(&r,step);
  if(full_moment_native_factor(comm,&r.service,&r.request,&accepted,&midpoint,r.native_stage,r.root,r.rates,full_path,full_epoch,r.h,step*r.h,r.threads,r.validation||check))goto close;
  if(save_factor_audit(&r,full_epoch))goto close;
  full_moment_state_close(&midpoint);
  if(release(&r))goto close;
  char full_geometry[512]={0};if(geometry(&r,full_epoch,full_path,full_geometry))goto close;
  if(configure(&r,full_epoch,full_path,full_geometry,check))goto close;
  char label[64];snprintf(label,sizeof label,"candidate_audit_%04d",step);
  if(admissibility(&r,step,check,label))goto close;
  status=full_moment_state_load(full_path,&candidate);if(collective(&r,status))goto close;
  }
  /* Every rank has passed factor/reconstruction/geometry/conservation gates.
   * Native commit occurs only here; the previous accepted snapshot survives. */
  if(!r.rank){char saved[512];snprintf(saved,sizeof saved,"%s/accepted_%04d",r.output,step);status=full_moment_state_publish(saved,&candidate);}
  if(collective(&r,status))goto close;
  full_moment_state_close(&accepted);accepted=candidate;memset(&candidate,0,sizeof candidate);
  if(report(&r,step,&accepted,MPI_Wtime()-step_start,0))goto close;
  printf("[FULL_MOMENT_ACCEPTED] rank=%d step=%d time=%.17g h=%.17g providers=%d seconds=%.9g\n",r.rank,step,accepted.header.time,r.h,r.providers,MPI_Wtime()-step_start);fflush(stdout);
 }
 if(report(&r,r.steps,&accepted,MPI_Wtime()-started,1))goto close;
 if(release(&r))goto close;
 code=0;
close:
 if(code&&!r.rank)fprintf(stderr,"[FULL_MOMENT_FAILED] epoch=%llu providers=%d; previous accepted files retained, no retry/dt change/legacy fallback\n",(unsigned long long)r.epoch,r.providers);
 full_moment_state_close(&candidate);full_moment_state_close(&midpoint);full_moment_state_close(&accepted);full_moment_state_close(&initial);
 if(r.service.pid>0){if(!code){FullMomentProviderResponse response;status=rpc(&r,FMP_STOP,&response);if(status)code=1;}if(full_moment_provider_close(&r.service))code=1;}
 if(code)preserve_failure(&r);
 return collective(&r,code)?114:0;
}
