/* Exercise the real typed provider verifier + native rebind, no RHS/MPI. */
#define _POSIX_C_SOURCE 200809L
#include "full_moment_provider.h"
#include "full_moment_state_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <float.h>
#include <sys/stat.h>
int main(int argc,char **argv){
 if(argc!=8)return 2;
 FullMomentState original={0},imported={0};FullMomentProvider service={0};
 if(full_moment_state_load(argv[3],&original))return 3;
 uint64_t token=((uint64_t)time(NULL)<<24)^(uint64_t)getpid();
 char memory[512],output[512],log[1024];snprintf(memory,sizeof memory,"/dev/shm/eunha_run_%llu",(unsigned long long)token);
 snprintf(output,sizeof output,"%s/restart_verified",memory);snprintf(log,sizeof log,"%s/provider.log",argv[7]);
 if(mkdir(argv[7],0700)||mkdir(memory,0700))return 4;
 if(full_moment_provider_start(&service,argv[1],argv[2],1,log))return 5;
 FullMomentProviderRequest q={0};FullMomentProviderResponse response;
 q.epoch=1;q.run_token=token;q.rank=0;q.nranks=1;q.ncells=original.header.cells;
 memcpy(q.physics_sha,"5128e9d2475f60389ba6b59aaeab34a64d5f75d96bbaf91fa9d13ee9f8d55efc",64);
 snprintf(q.state_path,512,"%s",argv[3]);snprintf(q.bounds_path,512,"%s",argv[4]);
 snprintf(q.fits_path,512,"%s",argv[5]);snprintf(q.aux_path,512,"%s",argv[6]);snprintf(q.output_dir,512,"%s",output);
 double h=.001;memcpy(q.reserved,&h,sizeof h);q.command=FMP_VERIFY_RESTART;
 int status=full_moment_provider_call(&service,&q,&response);
 char header_file[1024];snprintf(header_file,sizeof header_file,"%s/expected_header.bin",response.output_path);
 FullMomentStateHeader expected={0};FILE *f=status?NULL:fopen(header_file,"rb");
 if(!f)status=1;else{status|=fread(&expected,1,sizeof expected,f)!=sizeof expected;status|=fgetc(f)!=EOF;fclose(f);}
 if(!status)status=memcmp(&expected,&original.header,sizeof expected)!=0;
 char imported_path[1024];snprintf(imported_path,sizeof imported_path,"%s/imported",argv[7]);
 if(!status)status=full_moment_state_import(argv[3],imported_path,&expected,token,&imported);
 if(!status){
  for(int k=0;k<6;k++)status|=original.bytes[k]!=imported.bytes[k]||memcmp(original.mapping[k],imported.mapping[k],original.bytes[k]);
  status|=imported.header.epoch!=1||imported.header.run_token!=token||imported.header.time!=original.header.time;
 }
 q.command=FMP_STOP;status|=full_moment_provider_call(&service,&q,&response)!=0;
 status|=full_moment_provider_close(&service)!=0;
 if(!status){
  char report[1024];snprintf(report,sizeof report,"%s/report.json",argv[7]);FILE *out=fopen(report,"wx");
  if(!out)status=1;else{fprintf(out,"{\"status\":\"passed_typed_restart_verification_and_native_import\",\"payload_bytes_exact\":true,\"time\":%.17g,\"new_epoch\":1,\"physical_providers\":0,\"accepted_time_steps\":0,\"runtime_physical_restart_verified\":false}\n",imported.header.time);fclose(out);}
 }
 full_moment_state_close(&original);full_moment_state_close(&imported);
 /* Remove only our two read-only verification artifacts after all checks. */
 if(!status){char file[1024];snprintf(file,sizeof file,"%s/expected_header.bin",output);unlink(file);snprintf(file,sizeof file,"%s/receipt.json",output);unlink(file);rmdir(output);rmdir(memory);}
 return status?6:0;
}
