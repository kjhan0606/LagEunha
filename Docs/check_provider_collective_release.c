#define _POSIX_C_SOURCE 200809L
#include "full_moment_provider_mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
static int collective(MPI_Comm c,int value){int out;MPI_Allreduce(&value,&out,1,MPI_INT,MPI_MAX,c);return out;}
int main(int argc,char **argv){
 int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);MPI_Comm c=MPI_COMM_WORLD;
 int rank,size;MPI_Comm_rank(c,&rank);MPI_Comm_size(c,&size);int status=argc!=8||size!=2||provided<MPI_THREAD_FUNNELED;
 if(collective(c,status)){MPI_Finalize();return 2;}
 uint64_t token=rank?0:(((uint64_t)time(NULL)<<24)^(uint64_t)getpid());MPI_Bcast(&token,1,MPI_UINT64_T,0,c);
 char memory[512],geo_parent[512],epoch_root[512],geometry[512],owned[512],log[1024];
 snprintf(memory,sizeof memory,"/dev/shm/eunha_run_%llu",(unsigned long long)token);
 snprintf(geo_parent,sizeof geo_parent,"%s/geometry",memory);snprintf(epoch_root,sizeof epoch_root,"%s/e1",geo_parent);snprintf(geometry,sizeof geometry,"%s/geometry",epoch_root);
 snprintf(owned,sizeof owned,"/dev/shm/eunha_fmp_release_test_%llu_r%d",(unsigned long long)token,rank);snprintf(log,sizeof log,"%s/rank%d.log",argv[7],rank);
 if(!rank){
  status=mkdir(memory,0700)||mkdir(geo_parent,0700)||mkdir(epoch_root,0700)||mkdir(geometry,0700);
  for(int k=0;k<2&&!status;k++){char src[1024],dest[1024];const char *name=k?"full_exact_supports.json":"mesh.npz";
   snprintf(src,sizeof src,"%s/%s",argv[4],name);snprintf(dest,sizeof dest,"%s/%s",geometry,name);status=symlink(src,dest)!=0;}
 }
 if(collective(c,status)){MPI_Finalize();return 3;}
 FullMomentProvider provider={0};FullMomentProviderRequest request={0};FullMomentProviderResponse response;
 status=full_moment_provider_start(&provider,argv[1],argv[2],1,log);
 if(collective(c,status)){if(provider.pid>0)full_moment_provider_close(&provider);MPI_Finalize();return 4;}
 request.command=FMP_CONFIGURE;request.epoch=1;request.run_token=token;request.rank=rank;request.nranks=size;request.ncells=19945;
 memcpy(request.physics_sha,argv[6],64);strcpy(request.state_path,argv[3]);strcpy(request.geometry_path,geometry);strcpy(request.output_dir,owned);strcpy(request.aux_path,argv[5]);
 status=full_moment_provider_call(&provider,&request,&response);
 if(!collective(c,status)){
  /* Rank0 must keep shared geometry while rank1 is still before its hash. */
  if(rank==1){struct timespec delay={1,500000000};nanosleep(&delay,NULL);}
  status=full_moment_collective_release(c,&provider,&request,geometry);
  status|=access(geometry,F_OK)==0;
 }
 status=collective(c,status);
 if(!status){request.command=FMP_STOP;status=full_moment_provider_call(&provider,&request,&response);}
 status|=full_moment_provider_close(&provider);status=collective(c,status);
 if(!status&&!rank){
  char file[1024];snprintf(file,sizeof file,"%s/report.json",argv[7]);FILE *f=fopen(file,"wx");status=!f;
  if(f){fprintf(f,"{\"status\":\"passed_collective_two_phase_release\",\"ranks\":2,\"delayed_rank\":1,\"delay_seconds\":1.5,\"late_input_hash_passed\":true,\"shared_geometry_deleted_after_all_ACK\":true,\"actual_providers\":0,\"accepted_steps\":0}\n");status|=fclose(f)!=0;}
  status|=rmdir(geo_parent)!=0||rmdir(memory)!=0;
 }
 status=collective(c,status);MPI_Finalize();return status?5:0;
}
