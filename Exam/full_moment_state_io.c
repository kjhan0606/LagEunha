#define _POSIX_C_SOURCE 200809L
#include "full_moment_state_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
static const char *names[6]={"ids","mass","state","factors","Psi","weight_rate_seed"};
static int valid_header(const FullMomentStateHeader *h){
 uint16_t endian=1;if(*(unsigned char*)&endian!=1)return 0;
 if(h->magic!=FM_STATE_MAGIC||h->version!=1||h->cells<1||h->cells>INT_MAX/20
  ||h->state_elements!=20*h->cells+4||h->ld_bytes!=sizeof(long double)||h->ld_mantissa!=LDBL_MANT_DIG||!isfinite(h->time)||h->time<0)return 0;
 for(int i=0;i<7;i++)if(h->reserved[i])return 0;return 1;
}
void full_moment_state_close(FullMomentState *s){
 if(!s)return;for(int k=0;k<6;k++)if(s->mapping[k])munmap(s->mapping[k],s->bytes[k]);memset(s,0,sizeof *s);
}
int full_moment_state_load(const char *root,FullMomentState *out){
 if(!root||!out)return 1;FullMomentState s={0};char path[4096];snprintf(path,sizeof path,"%s/header.bin",root);
 FILE *f=fopen(path,"rb");if(!f)return 2;int bad=fread(&s.header,1,sizeof s.header,f)!=sizeof s.header;
 bad|=fgetc(f)!=EOF;fclose(f);if(bad||!valid_header(&s.header))return 3;
 uint64_t n=s.header.cells;size_t bytes[6]={n*8,n*8,(20*n+4)*8,n*4*8,n*sizeof(long double),n*8};
 for(int k=0;k<6;k++){
  snprintf(path,sizeof path,"%s/%s.bin",root,names[k]);int fd=open(path,O_RDONLY);struct stat st;
  if(fd<0){full_moment_state_close(&s);return 4;}
  if(fstat(fd,&st)||st.st_size!=(off_t)bytes[k]){close(fd);full_moment_state_close(&s);return 5;}
  void *p=mmap(NULL,bytes[k],PROT_READ,MAP_SHARED,fd,0);close(fd);
  if(p==MAP_FAILED){full_moment_state_close(&s);return 6;}s.mapping[k]=p;s.bytes[k]=bytes[k];
 }
 s.ids=s.mapping[0];s.mass=s.mapping[1];s.state=s.mapping[2];s.factors=s.mapping[3];s.Psi=s.mapping[4];s.weight_rate_seed=s.mapping[5];
 for(uint64_t i=0;i<n;i++){
  if((i&&s.ids[i]<=s.ids[i-1])||!isfinite(s.mass[i])||s.mass[i]<=0||!isfinite(s.Psi[i])||!isfinite(s.weight_rate_seed[i])){full_moment_state_close(&s);return 7;}
 }
 for(uint64_t i=0;i<20*n+4;i++)if(!isfinite(s.state[i])){full_moment_state_close(&s);return 7;}
 for(uint64_t i=0;i<4*n;i++)if(!isfinite(s.factors[i])){full_moment_state_close(&s);return 7;}
 *out=s;return 0;
}
static int write_file(const char *root,const char *name,const void *data,size_t bytes){
 char path[4096];snprintf(path,sizeof path,"%s/%s.bin",root,name);FILE *f=fopen(path,"wbx");if(!f)return 1;
 int bad=fwrite(data,1,bytes,f)!=bytes;bad|=fflush(f)!=0;
 if(!bad)bad=fsync(fileno(f))!=0;bad|=fclose(f)!=0;return bad;
}
int full_moment_state_write(const char *root,const FullMomentState *s){
 if(!root||!s||!valid_header(&s->header)||!s->ids||!s->mass||!s->state||!s->factors||!s->Psi||!s->weight_rate_seed)return 1;
 uint64_t n=s->header.cells;
 if(mkdir(root,0700))return 2;
 const void *arrays[6]={s->ids,s->mass,s->state,s->factors,s->Psi,s->weight_rate_seed};size_t bytes[6]={n*8,n*8,(20*n+4)*8,n*4*8,n*sizeof(long double),n*8};
 /* Header written last: an incomplete directory is never a readable state. */
 for(int k=0;k<6;k++)if(write_file(root,names[k],arrays[k],bytes[k]))return 3;
 if(write_file(root,"header",&s->header,sizeof s->header))return 3;
 int directory=open(root,O_RDONLY|O_DIRECTORY);if(directory<0)return 4;
 int status=fsync(directory)!=0;close(directory);if(status)return 4;
 return 0;
}
int full_moment_state_publish(const char *root,const FullMomentState *s){
 if(!root||strlen(root)>3500)return 1;struct stat st;if(lstat(root,&st)==0)return 2;
 char temporary[4096];snprintf(temporary,sizeof temporary,"%s.pending",root);
 int status=full_moment_state_write(temporary,s);if(status)return status;
 /* One authoritative rank0 writer; a complete directory becomes visible in
  * one rename. Existing accepted states are never opened for writing. */
 if(rename(temporary,root))return 5;
 char parent[4096];snprintf(parent,sizeof parent,"%s",root);char *slash=strrchr(parent,'/');
 if(!slash)return 6;*slash=0;int directory=open(parent,O_RDONLY|O_DIRECTORY);if(directory<0)return 6;
 status=fsync(directory)!=0;close(directory);return status?6:0;
}
int full_moment_state_import(const char *source,const char *destination,
 const FullMomentStateHeader *expected,uint64_t new_token,FullMomentState *out){
 if(!source||!destination||!expected||!out||!valid_header(expected)
  ||!expected->epoch||!expected->run_token||!new_token||new_token==expected->run_token)return 1;
 FullMomentState loaded={0};int status=full_moment_state_load(source,&loaded);
 if(status)return status;
 if(memcmp(&loaded.header,expected,sizeof *expected)){
  full_moment_state_close(&loaded);return 8;
 }
 /* Read-only mappings retain every payload bit, including long-double
  * padding bytes. Header rebind is explicit and confined to this new copy. */
 loaded.header.epoch=1;loaded.header.run_token=new_token;
 status=full_moment_state_publish(destination,&loaded);
 full_moment_state_close(&loaded);
 if(status)return status;
 return full_moment_state_load(destination,out);
}
