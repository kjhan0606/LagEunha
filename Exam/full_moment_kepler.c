/* Opt-in native executable integration for the full moment model.
 * The first integration gate exercises real MPI-owned CUDA stages on an
 * immutable capture. It is explicitly not an accepted evolution step.
 * The normal legacy RK4 path is selected when this model is absent.
 */
#include "full_moment_kepler.h"
#include "moment_owned_backend.h"
#include "full_moment_cpu_mpi.h"
#include "full_moment_evolution.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>

extern int lagmoment_stage(int,int,int,int,int64_t,int64_t,const double*,const double*,const int64_t*,const int64_t*,const int64_t*,const int64_t*,const double*,const double*,const double*,const int64_t*,double*,double*,double*,double*);
int full_moment_kepler_is_config(const char *path){
 FILE *f=fopen(path,"r");if(!f)return 0;char line[2048];int found=0;
 while(fgets(line,sizeof line,f))if(strstr(line,"define Simulation Model")&&strstr(line,"FullMomentKepler"))found=1;
 fclose(f);return found;
}
static int collective(MPI_Comm c,int status){int result;MPI_Allreduce(&status,&result,1,MPI_INT,MPI_MAX,c);return result;}
static void *read_array(const char *root,const char *name,size_t width,int64_t *count){
 char path[4096];snprintf(path,sizeof path,"%s/%s",root,name);FILE *f=fopen(path,"rb");if(!f)return NULL;
 if(fseek(f,0,SEEK_END)||ftell(f)<0){fclose(f);return NULL;}
 long size=ftell(f);rewind(f);if(size%(long)width){fclose(f);return NULL;}
 void *p=malloc((size_t)size);if(!p||fread(p,1,(size_t)size,f)!=(size_t)size){free(p);fclose(f);return NULL;}
 fclose(f);*count=size/(int64_t)width;return p;
}
static int share(MPI_Comm c,void **p,int64_t elements,MPI_Datatype type,int rank){
 int bytes=0;MPI_Type_size(type,&bytes);if(rank)*p=malloc((size_t)elements*bytes);
 int status=collective(c,!*p);if(status)return status;
 for(int64_t lo=0;lo<elements;){int chunk=(int)((elements-lo)>(1<<25)?(1<<25):elements-lo);
  MPI_Bcast((char*)*p+(size_t)lo*bytes,chunk,type,0,c);lo+=chunk;}
 return 0;
}
static int write_array(const char *root,const char *name,const void *a,size_t bytes){
 char path[4096];snprintf(path,sizeof path,"%s/%s",root,name);FILE *f=fopen(path,"wb");if(!f)return 1;
 int error=fwrite(a,1,bytes,f)!=bytes;error|=fclose(f)!=0;return error;
}
int full_moment_kepler_run(const char *path,MPI_Comm c){
 int rank,size;MPI_Comm_rank(c,&rank);MPI_Comm_size(c,&size);
 char capture[2048]="",prep[2048]="",stage[2048]="",output[2048]="",line[4096];int validation=0,threads=1;
 FILE *f=fopen(path,"r");int status=!f;int cpu_validation=0,native_evolution=0;
 if(f){while(fgets(line,sizeof line,f)){
  sscanf(line,"moment_capture_root = %2047s",capture);sscanf(line,"moment_prepare_root = %2047s",prep);
  sscanf(line,"moment_stage_root = %2047s",stage);sscanf(line,"moment_output = %2047s",output);
  sscanf(line,"moment_bridge_validation = %d",&validation);sscanf(line,"moment_threads = %d",&threads);
  sscanf(line,"moment_cpu_validation = %d",&cpu_validation);
  sscanf(line,"moment_native_evolution = %d",&native_evolution);
 }fclose(f);}
 if(native_evolution==1)return full_moment_evolution_run(path,c);
 if(cpu_validation==1)return full_moment_cpu_mpi_validation(path,c);
 status|=validation!=1||!capture[0]||!prep[0]||!stage[0]||!output[0]||threads<1;
 if(collective(c,status)){if(!rank)fprintf(stderr,"[FULL_MOMENT] Native bridge validation configuration required; evolution runtime not yet enabled. No legacy fallback.\n");return 77;}
 char host[256],gpu[256];unsigned char uuid[16];int visible=0;
 MPI_Comm local_communicator;int local_rank=0;
 MPI_Comm_split_type(c,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&local_communicator);
 MPI_Comm_rank(local_communicator,&local_rank);MPI_Comm_free(&local_communicator);
 gethostname(host,sizeof host);int ordinal=local_rank;
 /* Scheduler may expose one bound device per rank. Try ordinal zero first
  * in that case; the global UUID check below rejects accidental sharing. */
 status=lagmoment_select_visible_device(0,&visible,gpu,uuid);
 if(!status&&visible>1)status=lagmoment_select_visible_device(ordinal%visible,&visible,gpu,uuid);
 if(collective(c,status))return 78;
 unsigned char *uuids=malloc(16*size);MPI_Allgather(uuid,16,MPI_BYTE,uuids,16,MPI_BYTE,c);
 for(int i=0;i<size;i++)for(int j=0;j<i;j++)if(!memcmp(uuids+16*i,uuids+16*j,16))status=1;
 free(uuids);if(collective(c,status)){if(!rank)fprintf(stderr,"[FULL_MOMENT] GPU UUID collision across ranks\n");return 79;}
 char uuid_hex[33];for(int i=0;i<16;i++)snprintf(uuid_hex+2*i,3,"%02x",uuid[i]);
 printf("[FULL_MOMENT_ID] rank=%d ranks=%d host=%s pid=%ld gpu=%s UUID=%s visible=%d binary_sha=%s\n",rank,size,host,(long)getpid(),gpu,uuid_hex,visible,getenv("EUNHA_MOMENT_BINARY_SHA")?getenv("EUNHA_MOMENT_BINARY_SHA"):"unrecorded");fflush(stdout);
 int64_t N=0,NN=0,NI=0,NG=0,t=0;
 double *a=NULL,*raw=NULL,*extra=NULL,*L=NULL,*gn=NULL,*nf=NULL;
 int64_t *owners=NULL,*off=NULL,*inc=NULL,*fr=NULL,*goff=NULL;
 if(!rank){
  a=read_array(capture,"inputs.bin",33*sizeof(double),&NN);
  raw=read_array(prep,"cells.bin",19*sizeof(double),&N);
  owners=read_array(prep,"owners.bin",2*sizeof(int64_t),&t);status|=!owners||t!=NN;
  off=read_array(prep,"offsets.bin",sizeof(int64_t),&t);status|=!off||t!=N+1;
  inc=read_array(prep,"incidence.bin",sizeof(int64_t),&NI);
  fr=read_array(prep,"frozen.bin",sizeof(int64_t),&t);status|=!fr||t!=N;
  extra=read_array(stage,"node_extra.bin",20*sizeof(double),&t);status|=!extra||t!=NN;
  L=read_array(stage,"factors.bin",4*sizeof(double),&t);status|=!L||t!=N;
  gn=read_array(stage,"gravity_nodes.bin",7*sizeof(double),&NG);
  goff=read_array(stage,"gravity_offsets.bin",sizeof(int64_t),&t);status|=!goff||t!=N+1;
  status|=!a||!raw||!inc||!gn||N<size||N>INT_MAX/106||NN>INT_MAX/33||NI>INT_MAX/16;
 }
 if(collective(c,status))return 80;
 MPI_Bcast(&N,1,MPI_INT64_T,0,c);MPI_Bcast(&NN,1,MPI_INT64_T,0,c);MPI_Bcast(&NI,1,MPI_INT64_T,0,c);MPI_Bcast(&NG,1,MPI_INT64_T,0,c);
 double begin=MPI_Wtime();nf=rank?NULL:malloc(19*NN*sizeof(double));if(!rank&&!nf)status=1;
 if(!rank&&!status)status=lagmoment_canonical_nodes(1,threads,(int)NN,a,nf);
 if(collective(c,status))return 81;
 if(share(c,(void**)&a,33*NN,MPI_DOUBLE,rank)||share(c,(void**)&nf,19*NN,MPI_DOUBLE,rank)
  ||share(c,(void**)&raw,19*N,MPI_DOUBLE,rank)||share(c,(void**)&extra,20*NN,MPI_DOUBLE,rank)
  ||share(c,(void**)&L,4*N,MPI_DOUBLE,rank)||share(c,(void**)&gn,7*NG,MPI_DOUBLE,rank)
  ||share(c,(void**)&owners,2*NN,MPI_INT64_T,rank)||share(c,(void**)&off,N+1,MPI_INT64_T,rank)
  ||share(c,(void**)&inc,NI,MPI_INT64_T,rank)||share(c,(void**)&fr,N,MPI_INT64_T,rank)
  ||share(c,(void**)&goff,N+1,MPI_INT64_T,rank))return 82;
 int64_t lo=N*rank/size,hi=N*(rank+1)/size,n=hi-lo,ni=off[hi]-off[lo],ng=goff[hi]-goff[lo];
 int64_t *owned=malloc(n*sizeof(int64_t)),*node_ids=malloc(NN*sizeof(int64_t)),*of=malloc((n+1)*sizeof(int64_t)),*gf=malloc((n+1)*sizeof(int64_t));
 double *out=malloc(n*94*sizeof(double)),*q=malloc(ni*16*sizeof(double)),*hydro=malloc(n*106*sizeof(double));
 status=!owned||!node_ids||!of||!gf||!out||!q||!hydro;
 if(collective(c,status))return 83;
 for(int64_t i=0;i<n;i++)owned[i]=lo+i;for(int64_t i=0;i<NN;i++)node_ids[i]=i;
 for(int64_t i=0;i<=n;i++){of[i]=off[lo+i]-off[lo];gf[i]=goff[lo+i]-goff[lo];}
 double kernel_begin=MPI_Wtime();
 status=lagmoment_owned_stage(1,threads,1,1,(int)N,(int)n,(int)NN,ni,ng,owned,node_ids,a,nf,raw+19*lo,owners,of,inc+off[lo],fr+lo,extra,L+4*lo,gn+7*goff[lo],gf,out,q,hydro);
 double kernel_seconds=MPI_Wtime()-kernel_begin;
 if(collective(c,status))return 84;
 printf("[FULL_MOMENT_OWNED] rank=%d owned_cells=%lld incidences=%lld gravity_nodes=%lld actual_GPU_stage_calls=1 stage_seconds=%.9g\n",rank,(long long)n,(long long)ni,(long long)ng,kernel_seconds);fflush(stdout);
 int *counts=malloc(size*sizeof(int)),*displacements=malloc(size*sizeof(int));
 double *allout=rank?NULL:malloc(N*94*sizeof(double)),*allq=rank?NULL:malloc(NI*16*sizeof(double)),*allhydro=rank?NULL:malloc(N*106*sizeof(double));
 status=!counts||!displacements||(!rank&&(!allout||!allq||!allhydro));if(collective(c,status))return 85;
 for(int r=0;r<size;r++){int64_t l=N*r/size,h=N*(r+1)/size;counts[r]=(int)((h-l)*94);displacements[r]=(int)(l*94);}
 MPI_Gatherv(out,(int)(n*94),MPI_DOUBLE,allout,counts,displacements,MPI_DOUBLE,0,c);
 for(int r=0;r<size;r++){int64_t l=N*r/size,h=N*(r+1)/size;counts[r]=(int)((off[h]-off[l])*16);displacements[r]=(int)(off[l]*16);}
 MPI_Gatherv(q,(int)(ni*16),MPI_DOUBLE,allq,counts,displacements,MPI_DOUBLE,0,c);
 for(int r=0;r<size;r++){int64_t l=N*r/size,h=N*(r+1)/size;counts[r]=(int)((h-l)*106);displacements[r]=(int)(l*106);}
 MPI_Gatherv(hydro,(int)(n*106),MPI_DOUBLE,allhydro,counts,displacements,MPI_DOUBLE,0,c);
 double elapsed=MPI_Wtime()-begin,max_elapsed=0.;MPI_Reduce(&elapsed,&max_elapsed,1,MPI_DOUBLE,MPI_MAX,0,c);
 if(!rank){
  double *oracle=malloc(N*94*sizeof(double)),*oq=malloc(NI*16*sizeof(double)),*oh=malloc(N*106*sizeof(double)),ledger[12];
  status=!oracle||!oq||!oh;
  if(!status)status=lagmoment_stage(1,threads,(int)N,(int)NN,NI,NG,a,raw,owners,off,inc,fr,extra,L,gn,goff,oracle,oq,oh,ledger);
  if(!status)status=memcmp(allout,oracle,N*94*sizeof(double))||memcmp(allq,oq,NI*16*sizeof(double))||memcmp(allhydro,oh,N*106*sizeof(double));
  if(!status)status=write_array(output,"stage.bin",allout,N*94*sizeof(double))||write_array(output,"Q.bin",allq,NI*16*sizeof(double))||write_array(output,"hydro.bin",allhydro,N*106*sizeof(double));
  char file[4096];snprintf(file,sizeof file,"%s/native_result.json",output);FILE *report=fopen(file,"w");
  if(report){fprintf(report,"{\"status\":\"%s\",\"ranks\":%d,\"cells\":%lld,\"canonical_nodes\":%lld,\"production_node_passes\":1,\"extra_oracle_node_passes\":1,\"owned_stage_calls\":%d,\"MPI_stage_wall_seconds\":%.17g,\"evolution_steps\":0,\"MPI_evolution_verified\":false}\n",status?"failed_native_owned_MPI_component":"passed_native_owned_MPI_component",size,(long long)N,(long long)NN,size,max_elapsed);fclose(report);}else status=1;
  free(oracle);free(oq);free(oh);
 }
 status=collective(c,status);
 free(a);free(nf);free(raw);free(extra);free(L);free(gn);free(owners);free(off);free(inc);free(fr);free(goff);
 free(owned);free(node_ids);free(of);free(gf);free(out);free(q);free(hydro);free(counts);free(displacements);free(allout);free(allq);free(allhydro);
 return status?86:0;
}
