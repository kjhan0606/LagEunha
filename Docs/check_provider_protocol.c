#define _POSIX_C_SOURCE 200809L
#include "full_moment_provider.h"
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>

int main(int argc,char **argv){
 if(argc!=9)return 2;
 int result[3]={0};uint64_t token=((uint64_t)time(NULL)<<24)^(uint64_t)getpid();
 for(int k=0;k<3;k++){
  FullMomentProvider service={0};FullMomentProviderRequest req={0};FullMomentProviderResponse res={0};
  char log[1024],root[512];snprintf(log,sizeof log,"%s/case%d.log",argv[8],k);
  snprintf(root,sizeof root,"/dev/shm/eunha_fmp_protocol_%llu_case%d",(unsigned long long)token,k);
  if(full_moment_provider_start(&service,argv[1],argv[2],1,log))return 3;
  req.command=FMP_CONFIGURE;req.epoch=1;req.run_token=token+k;req.rank=0;req.nranks=1;req.ncells=19945;
  snprintf(req.state_path,sizeof req.state_path,"%s",argv[3]);snprintf(req.geometry_path,sizeof req.geometry_path,"%s",argv[4]);
  snprintf(req.aux_path,sizeof req.aux_path,"%s",argv[5]);memcpy(req.physics_sha,argv[6],64);
  snprintf(req.output_dir,sizeof req.output_dir,"%s",root);
  if(k==2)memset(req.physics_sha,'0',64);
  int status=full_moment_provider_call(&service,&req,&res);
  if(k==2){if(status!=101||full_moment_provider_close(&service)!=3)return 4;result[k]=1;continue;}
  if(status)return 5;
  if(k==0){
   uint64_t sequence=service.sequence;req.command=FMP_STOP;req.run_token++;
   memset(&res,0x6d,sizeof res);FullMomentProviderResponse before=res;
   if(full_moment_provider_call(&service,&req,&res)!=2||memcmp(&before,&res,sizeof res)||service.sequence!=sequence)return 6;
   req.run_token--;if(full_moment_provider_call(&service,&req,&res)||full_moment_provider_close(&service))return 7;
   char context[1024];snprintf(context,sizeof context,"%s/context.json",root);if(unlink(context)||rmdir(root))return 8;
   result[k]=1;
  }else{
   req.command=FMP_BOUNDS;req.epoch=0;
   if(full_moment_provider_call(&service,&req,&res)!=101||full_moment_provider_close(&service)!=3)return 9;
   result[k]=1;
  }
 }
 char path[1024];snprintf(path,sizeof path,"%s/report.json",argv[8]);FILE *f=fopen(path,"wx");if(!f)return 10;
 fprintf(f,"{\"status\":\"passed_native_provider_protocol\",\"request_bytes\":%zu,\"response_bytes\":%zu,\"wrong_token_no_IO_or_sequence_change\":%s,\"stale_epoch_rejected\":%s,\"wrong_physics_rejected\":%s,\"provider_calls\":0,\"evolution_steps\":0}\n",sizeof(FullMomentProviderRequest),sizeof(FullMomentProviderResponse),result[0]?"true":"false",result[1]?"true":"false",result[2]?"true":"false");
 return fclose(f)?11:0;
}
