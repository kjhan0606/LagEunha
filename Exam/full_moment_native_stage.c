#define _POSIX_C_SOURCE 200809L
#include "full_moment_native_stage.h"
#include "moment_owned_backend.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>

typedef struct {void *p;size_t bytes;} Buffer;
static int collective(MPI_Comm c,int status){int out;MPI_Allreduce(&status,&out,1,MPI_INT,MPI_MAX,c);return out;}
static Buffer map_file(const char *root,const char *name,size_t elements,size_t width,int create){
 Buffer out={0};if(!elements||!width||elements>SIZE_MAX/width)return out;
 char path[4096];snprintf(path,sizeof path,"%s/%s.bin",root,name);
 int flags=create==1?O_CREAT|O_EXCL|O_RDWR:create==2?O_RDWR:O_RDONLY;
 int fd=open(path,flags,0600);if(fd<0)return out;
 size_t bytes=elements*width;struct stat st;
 if((create==1&&ftruncate(fd,(off_t)bytes))||fstat(fd,&st)||st.st_size!=(off_t)bytes){close(fd);return out;}
 void *p=mmap(NULL,bytes,create?PROT_READ|PROT_WRITE:PROT_READ,MAP_SHARED,fd,0);close(fd);
 if(p==MAP_FAILED)return out;out.p=p;out.bytes=bytes;return out;
}
static int device(MPI_Comm c,int rank,int size){
 int local_rank,count=0;MPI_Comm local;MPI_Comm_split_type(c,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&local);
 MPI_Comm_rank(local,&local_rank);MPI_Comm_free(&local);
 char name[256];unsigned char uuid[16];int status=lagmoment_select_visible_device(0,&count,name,uuid);
 if(!status&&count>1)status=lagmoment_select_visible_device(local_rank%count,&count,name,uuid);
 if(collective(c,status))return 1;
 unsigned char *uuids=malloc(16*size);if(collective(c,!uuids)){free(uuids);return 1;}
 MPI_Allgather(uuid,16,MPI_BYTE,uuids,16,MPI_BYTE,c);
 for(int i=0;i<size;i++)for(int j=0;j<i;j++)status|=memcmp(uuids+16*i,uuids+16*j,16)==0;
 free(uuids);if(collective(c,status))return 1;
 char hex[33];for(int k=0;k<16;k++)snprintf(hex+2*k,3,"%02x",uuid[k]);
 printf("[FULL_MOMENT_ACTUAL_GPU] rank=%d ranks=%d GPU=%s UUID=%s\n",rank,size,name,hex);fflush(stdout);return 0;
}
int full_moment_native_select_device(MPI_Comm c){int rank,size;MPI_Comm_rank(c,&rank);MPI_Comm_size(c,&size);return device(c,rank,size);}
int full_moment_native_stage(MPI_Comm c,const char *global_root,const char *owned_root,uint64_t epoch,uint64_t token,int threads){
 int rank,size;MPI_Comm_rank(c,&rank);MPI_Comm_size(c,&size);int status=0,result=1;
 Buffer buffers[20]={{0}};int nb=0;Buffer header=map_file(global_root,"header",6,8,0);buffers[nb++]=header;
 uint64_t *h=header.p;int64_t N=0,NN=0,NI=0,NG=0;
 if(!h)status=1;else{
  status=h[0]!=epoch||h[1]!=token||h[2]<(uint64_t)size||h[2]>INT_MAX/106||h[3]<1||h[3]>INT_MAX||h[4]<1||h[4]>INT_MAX/16||h[5]<h[2]||h[5]>INT64_MAX/7;
  if(!status){N=h[2];NN=h[3];NI=h[4];NG=h[5];}
 }
 if(collective(c,status))goto cleanup;
 if(device(c,rank,size))goto cleanup;
 double start=MPI_Wtime();
#define MAP(type,name,count,root,file,create) Buffer name##_buffer=map_file(root,file,(size_t)(count),sizeof(type),create); buffers[nb++]=name##_buffer;type *name=name##_buffer.p;status|=!name
 MAP(double,a,33*NN,global_root,"node_inputs",0);
 MAP(double,raw,19*N,global_root,"cells",0);
 MAP(double,extra,20*NN,global_root,"extra",0);
 MAP(double,L,4*N,global_root,"factors",0);
 MAP(int64_t,owners,2*NN,global_root,"owners",0);
 MAP(int64_t,off,N+1,global_root,"offsets",0);
 MAP(int64_t,inc,NI,global_root,"incidence",0);
 MAP(int64_t,fr,N,global_root,"frozen",0);
 MAP(int64_t,goff,N+1,global_root,"gravity_offsets",0);
 if(collective(c,status))goto cleanup;
 int64_t lo=N*rank/size,hi=N*(rank+1)/size,n=hi-lo;
 status=off[0]!=0||off[N]!=NI||goff[0]!=0||goff[N]!=NG;
 for(int64_t i=0;i<N;i++)status|=off[i+1]<off[i]||goff[i+1]<=goff[i];
 if(collective(c,status))goto cleanup;
 int64_t ni=off[hi]-off[lo],ng=goff[hi]-goff[lo];
 MAP(double,gn,7*ng,owned_root,"gravity",0);
 if(collective(c,status))goto cleanup;
 Buffer flux={0};if(!rank)flux=map_file(global_root,"canonical_flux",19*NN,8,1);
 status=!rank&&!flux.p;if(collective(c,status)){if(flux.p)munmap(flux.p,flux.bytes);goto cleanup;}
 if(rank)flux=map_file(global_root,"canonical_flux",19*NN,8,0);
 buffers[nb++]=flux;double *nf=flux.p;
 if(collective(c,!nf))goto cleanup;
 if(!rank)status=lagmoment_canonical_nodes(1,threads,(int)NN,a,nf);
 if(collective(c,status))goto cleanup;
 /* Both ranks map the same node-local immutable flux file after this fence. */
 MPI_Barrier(c);
 int64_t *owned=malloc(n*sizeof(int64_t)),*nodes=malloc(NN*sizeof(int64_t));
 int64_t *of=malloc((n+1)*sizeof(int64_t)),*gf=malloc((n+1)*sizeof(int64_t));
 double *stage=malloc(n*94*sizeof(double)),*Q=malloc(ni*16*sizeof(double)),*hydro=malloc(n*106*sizeof(double));
 status=!owned||!nodes||!of||!gf||!stage||!Q||!hydro;
 if(collective(c,status))goto free_local;
 for(int64_t i=0;i<n;i++)owned[i]=lo+i;for(int64_t k=0;k<NN;k++)nodes[k]=k;
 for(int64_t i=0;i<=n;i++){of[i]=off[lo+i]-off[lo];gf[i]=goff[lo+i]-goff[lo];}
 double kernel_start=MPI_Wtime();
 status=lagmoment_owned_stage(1,threads,epoch,epoch,(int)N,(int)n,(int)NN,ni,ng,owned,nodes,a,nf,raw+19*lo,owners,of,inc+off[lo],fr+lo,extra,L+4*lo,gn,gf,stage,Q,hydro);
 if(status)fprintf(stderr,"[FULL_MOMENT_NATIVE_STAGE_ERROR] rank=%d epoch=%llu status=%d\n",rank,(unsigned long long)epoch,status);
 double kernel_seconds=MPI_Wtime()-kernel_start;
 if(collective(c,status))goto free_local;
 char output[4096];snprintf(output,sizeof output,"%s/outputs",global_root);
 Buffer allstage={0},allQ={0},allhydro={0},ledgers={0};
 if(!rank){
  allstage=map_file(output,"stage",N*94,8,1);
  allQ=map_file(output,"Q",NI*16,8,1);
  allhydro=map_file(output,"hydro",N*106,8,1);
  ledgers=map_file(output,"ledgers",12,8,1);
  buffers[nb++]=allstage;buffers[nb++]=allQ;buffers[nb++]=allhydro;buffers[nb++]=ledgers;
  status=!allstage.p||!allQ.p||!allhydro.p||!ledgers.p;
 }
 if(collective(c,status))goto free_local;
 /* Every rank is already constrained to one host. Create the shared outputs
  * once, then let each owner write its disjoint canonical slice directly. */
 if(rank){
  allstage=map_file(output,"stage",N*94,8,2);
  allQ=map_file(output,"Q",NI*16,8,2);
  allhydro=map_file(output,"hydro",N*106,8,2);
  buffers[nb++]=allstage;buffers[nb++]=allQ;buffers[nb++]=allhydro;
 }
 if(collective(c,!allstage.p||!allQ.p||!allhydro.p))goto free_local;
 memcpy((double *)allstage.p+94*lo,stage,(size_t)n*94*sizeof(double));
 memcpy((double *)allQ.p+16*off[lo],Q,(size_t)ni*16*sizeof(double));
 memcpy((double *)allhydro.p+106*lo,hydro,(size_t)n*106*sizeof(double));
 MPI_Barrier(c);
 if(!rank)status=lagmoment_ordered_ledgers(1,threads,(int)N,(int)NN,nf,owners,fr,allstage.p,allhydro.p,ledgers.p);
 if(collective(c,status))goto free_local;
 double wall=MPI_Wtime()-start,max_wall=0;MPI_Allreduce(&wall,&max_wall,1,MPI_DOUBLE,MPI_MAX,c);
 printf("[FULL_MOMENT_ACTUAL_OWNED] rank=%d owned_cells=%lld gravity_nodes=%lld stage_seconds=%.9g MPI_GPU_seconds=%.9g\n",rank,(long long)n,(long long)ng,kernel_seconds,max_wall);fflush(stdout);
 result=0;
free_local:
 free(owned);free(nodes);free(of);free(gf);free(stage);free(Q);free(hydro);
cleanup:
 for(int i=0;i<nb;i++)if(buffers[i].p)munmap(buffers[i].p,buffers[i].bytes);
 return result;
#undef MAP
}
