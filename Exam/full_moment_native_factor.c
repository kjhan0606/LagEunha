#define _POSIX_C_SOURCE 200809L
#include "full_moment_native_factor.h"
#include "full_moment_update.h"
#include "full_moment_native_stage.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
extern int lagmoment_factor(int,int,int,int64_t,const double*,const double*,const double*,const double*,const double*,const int64_t*,double,const double*,const int64_t*,double*);
typedef struct {void *p;size_t bytes;} Map;
static int all_status(MPI_Comm c,int status){int out;MPI_Allreduce(&status,&out,1,MPI_INT,MPI_MAX,c);return out;}
static Map map(const char *root,const char *name,size_t elements,size_t width){
 Map out={0};if(elements>SIZE_MAX/width||!elements)return out;
 char path[4096];snprintf(path,sizeof path,"%s/%s.bin",root,name);int fd=open(path,O_RDONLY);if(fd<0)return out;
 struct stat st;size_t bytes=elements*width;
 if(fstat(fd,&st)||st.st_size!=(off_t)bytes){close(fd);return out;}
 void *p=mmap(NULL,bytes,PROT_READ,MAP_SHARED,fd,0);close(fd);if(p==MAP_FAILED)return out;
 out.p=p;out.bytes=bytes;return out;
}
static int write_array(const char *root,const char *name,const void *data,size_t bytes){
 char path[4096];snprintf(path,sizeof path,"%s/%s.bin",root,name);FILE *f=fopen(path,"wbx");if(!f)return 1;
 int status=fwrite(data,1,bytes,f)!=bytes;status|=fclose(f)!=0;return status;
}
int full_moment_native_factor(MPI_Comm c,FullMomentProvider *service,FullMomentProviderRequest *req,
 const FullMomentState *base,const FullMomentState *linear,const char *stage_root,
 const char *owned_root,const char *rates_root,const char *candidate_root,
 uint64_t candidate_epoch,double h,double candidate_time,int threads,int shadow_audit){
 int rank,size;MPI_Comm_rank(c,&rank);MPI_Comm_size(c,&size);int result=1,status=0;
 int64_t N=base->header.cells,lo=N*rank/size,hi=N*(rank+1)/size,n=hi-lo;
 status=N<size||N>INT_MAX/20||linear->header.cells!=(uint64_t)N||base->header.run_token!=req->run_token||linear->header.run_token!=req->run_token
  ||linear->header.epoch!=req->epoch||candidate_epoch<=req->epoch||!isfinite(h)||h<=0||h>.005||!isfinite(candidate_time)||candidate_time<=base->header.time||threads<1;
 if(!status)status=memcmp(base->ids,linear->ids,N*8)||memcmp(base->mass,linear->mass,N*8);
 if(all_status(c,status))return 1;
 if(full_moment_native_select_device(c))return 1;
 Map maps[8]={{0}};int nm=0;
#define MAP(type,var,elements,root,name) Map var##_map=map(root,name,(size_t)(elements),sizeof(type));maps[nm++]=var##_map;type *var=var##_map.p;status|=!var
 MAP(uint64_t,header,6,stage_root,"header");
 if(all_status(c,status))goto cleanup;
 status=header[0]!=req->epoch||header[1]!=req->run_token||header[2]!=(uint64_t)N||header[4]<1||header[4]>INT_MAX/16;
 if(all_status(c,status))goto cleanup;
 int64_t NI=header[4];char outputs[4096];snprintf(outputs,sizeof outputs,"%s/outputs",stage_root);
 MAP(int64_t,csr,N+1,stage_root,"offsets");MAP(double,stage,94*N,outputs,"stage");
 MAP(double,allQ,16*NI,outputs,"Q");MAP(double,gen,2*N,rates_root,"generator");
 MAP(double,z,N,rates_root,"weight_rates");MAP(double,psi_rate,N,rates_root,"Psi_rates");
 MAP(double,ledger,12,rates_root,"ledgers");
 if(all_status(c,status))goto cleanup;
 status=csr[0]!=0||csr[N]!=NI;
 for(int64_t i=0;i<N;i++)status|=csr[i+1]<csr[i];
 if(all_status(c,status))goto cleanup;
 int64_t ni=csr[hi]-csr[lo],columns=4*ni;
 double *A=malloc(n*4*8),*K=malloc(n*4*8),*raw=malloc(n*18*8),*incoming=malloc(ni*8*8),*out=malloc(n*63*8),*deck=malloc(n*36*8);
 int64_t *offset=malloc((n+1)*8),*frozen=calloc(n,8);
 double *candidate=malloc(n*20*8),*factors=malloc(n*4*8);long double *Psi=calloc(n,sizeof(long double));
 int *counts=malloc(size*sizeof(int)),*displacements=malloc(size*sizeof(int));
 double *global_state=rank?NULL:malloc((20*N+4)*8),*global_factors=rank?NULL:malloc(N*4*8);
 long double *global_Psi=rank?NULL:calloc(N,sizeof(long double));
 status=!A||!K||!raw||!incoming||!out||!deck||!offset||!frozen||!candidate||!factors||!Psi||!counts||!displacements||(!rank&&(!global_state||!global_factors||!global_Psi));
 if(all_status(c,status))goto free_local;
 for(int64_t i=0;i<n;i++){
  int64_t g=lo+i;raw[18*i]=base->mass[g];
  for(int k=1;k<18;k++)raw[18*i+k]=base->state[20*g+k+2]+h*stage[94*g+k];
  for(int k=0;k<4;k++){A[4*i+k]=stage[94*g+56+k];K[4*i+k]=stage[94*g+22+k];}
  double *d=deck+36*i;
  memcpy(d,base->factors+4*g,4*8);memcpy(d+4,linear->factors+4*g,4*8);
  memcpy(d+8,A+4*i,4*8);memcpy(d+12,K+4*i,4*8);memcpy(d+16,raw+18*i,18*8);d[34]=h;d[35]=0.;
 }
 for(int64_t i=0;i<=n;i++)offset[i]=4*(csr[lo+i]-csr[lo]);
 for(int64_t t=0;t<ni;t++)memcpy(incoming+8*t,allQ+16*(csr[lo]+t)+8,8*8);
 double start=MPI_Wtime();
 status=lagmoment_factor(1,threads,(int)n,columns,base->factors+4*lo,linear->factors+4*lo,A,K,incoming,offset,h,raw,frozen,out);
 if(status)fprintf(stderr,"[FULL_MOMENT_NATIVE_FACTOR_ERROR] rank=%d epoch=%llu status=%d h=%.17g\n",rank,(unsigned long long)req->epoch,status,h);
 if(all_status(c,status))goto free_local;
 char audit[1024];snprintf(audit,sizeof audit,"%s/factor_%llu",owned_root,(unsigned long long)candidate_epoch);
 status=mkdir(audit,0700)!=0;uint64_t factor_header[6]={UINT64_C(0x31544341464d464c),req->epoch,req->run_token,(uint64_t)n,(uint64_t)columns,(uint64_t)rank};
 if(!status)status=write_array(audit,"header",factor_header,sizeof factor_header)||write_array(audit,"deck",deck,n*36*8)
  ||write_array(audit,"factor",out,n*63*8)||write_array(audit,"incoming_Q",incoming,columns*2*8)||write_array(audit,"incoming_offsets",offset,(n+1)*8);
 if(all_status(c,status))goto free_local;
 FullMomentProviderResponse response;req->command=FMP_FACTOR_AUDIT;uint32_t previous_outer=req->outer;req->outer=shadow_audit!=0;
 if(strlen(audit)>=sizeof req->output_dir)status=1;else strcpy(req->output_dir,audit);
 if(!status)status=full_moment_provider_call(service,req,&response);
 req->outer=previous_outer;
 if(all_status(c,status))goto free_local;
 status=full_moment_build_candidate((int)n,h,base->state+20*lo,base->mass+lo,stage+94*lo,gen+2*lo,z+lo,base->Psi+lo,psi_rate+lo,out,candidate,Psi,factors);
 if(all_status(c,status))goto free_local;
 for(int r=0;r<size;r++){int64_t l=N*r/size,e=N*(r+1)/size;counts[r]=(int)((e-l)*20);displacements[r]=(int)(l*20);}
 MPI_Gatherv(candidate,(int)(n*20),MPI_DOUBLE,global_state,counts,displacements,MPI_DOUBLE,0,c);
 for(int r=0;r<size;r++){int64_t l=N*r/size,e=N*(r+1)/size;counts[r]=(int)((e-l)*4);displacements[r]=(int)(l*4);}
 MPI_Gatherv(factors,(int)(n*4),MPI_DOUBLE,global_factors,counts,displacements,MPI_DOUBLE,0,c);
 for(int r=0;r<size;r++){int64_t l=N*r/size,e=N*(r+1)/size;counts[r]=(int)(e-l);displacements[r]=(int)l;}
 MPI_Gatherv(Psi,(int)n,MPI_LONG_DOUBLE,global_Psi,counts,displacements,MPI_LONG_DOUBLE,0,c);
 if(!rank){
  for(int k=0;k<4;k++){double rate=0.;for(int row=0;row<3;row++)rate+=ledger[4*row+k];global_state[20*N+k]=base->state[20*N+k]+h*rate;}
  FullMomentState tentative={0};tentative.header=base->header;tentative.header.epoch=candidate_epoch;tentative.header.time=candidate_time;
  tentative.ids=base->ids;tentative.mass=base->mass;tentative.state=global_state;tentative.factors=global_factors;tentative.Psi=global_Psi;tentative.weight_rate_seed=z;
  status=full_moment_state_write(candidate_root,&tentative);
 }
 if(all_status(c,status))goto free_local;
 printf("[FULL_MOMENT_TENTATIVE] rank=%d cells=%lld candidate_time=%.17g h=%.17g seconds=%.9g accepted=false\n",rank,(long long)n,candidate_time,h,MPI_Wtime()-start);fflush(stdout);
 result=0;
free_local:
 free(A);free(K);free(raw);free(incoming);free(out);free(deck);free(offset);free(frozen);free(candidate);free(factors);free(Psi);free(counts);free(displacements);free(global_state);free(global_factors);free(global_Psi);
cleanup:
 for(int i=0;i<nm;i++)if(maps[i].p)munmap(maps[i].p,maps[i].bytes);
 return result;
#undef MAP
}
