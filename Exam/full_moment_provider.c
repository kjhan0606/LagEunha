#define _POSIX_C_SOURCE 200809L
#include "full_moment_provider.h"
#include <errno.h>
#include <float.h>
#include <spawn.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <signal.h>
#include <unistd.h>
#include <sys/wait.h>
extern char **environ;
static int transfer(int fd,void *buffer,size_t length,int writing){
 char *p=buffer;while(length){ssize_t n=writing?write(fd,p,length):read(fd,p,length);
  if(n<0&&errno==EINTR)continue;if(n<=0)return 1;p+=n;length-=(size_t)n;}
 return 0;
}
int full_moment_provider_start(FullMomentProvider *s,const char *python,const char *script,int workers,const char *log){
 if(!s||!python||!script||workers<1||!log)return 1;
 /* A failed child must become a collective error, not kill one MPI rank. */
 struct sigaction action={0};action.sa_handler=SIG_IGN;sigemptyset(&action.sa_mask);
 if(sigaction(SIGPIPE,&action,NULL))return 2;
 int in[2],out[2];if(pipe(in))return 2;if(pipe(out)){close(in[0]);close(in[1]);return 2;}
 for(int i=0;i<2;i++){fcntl(in[i],F_SETFD,FD_CLOEXEC);fcntl(out[i],F_SETFD,FD_CLOEXEC);}
 posix_spawn_file_actions_t a;posix_spawn_file_actions_init(&a);
 posix_spawn_file_actions_adddup2(&a,in[0],STDIN_FILENO);posix_spawn_file_actions_adddup2(&a,out[1],STDOUT_FILENO);
 posix_spawn_file_actions_addopen(&a,STDERR_FILENO,log,O_CREAT|O_EXCL|O_WRONLY,0600);
 posix_spawn_file_actions_addclose(&a,in[1]);posix_spawn_file_actions_addclose(&a,out[0]);
 char count[32];snprintf(count,sizeof count,"%d",workers);char *args[]={(char*)python,(char*)script,"--workers",count,NULL};
 int result=posix_spawnp(&s->pid,python,&a,NULL,args,environ);posix_spawn_file_actions_destroy(&a);
 close(in[0]);close(out[1]);
 if(result){close(in[1]);close(out[0]);return 3;}
 s->request_fd=in[1];s->response_fd=out[0];s->sequence=0;s->run_token=0;return 0;
}
int full_moment_provider_call(FullMomentProvider *s,FullMomentProviderRequest *r,FullMomentProviderResponse *out){
 if(!s||!r||!out||s->pid<=0)return 1;
 if(!s->run_token)s->run_token=r->run_token;
 if(!r->run_token||r->run_token!=s->run_token)return 2;
 r->magic=FMP_MAGIC;r->version=FMP_ABI;r->request_id=++s->sequence;
 r->ld_bytes=sizeof(long double);r->ld_mantissa=LDBL_MANT_DIG;
 FullMomentProviderResponse response;
 if(transfer(s->request_fd,r,sizeof *r,1)||transfer(s->response_fd,&response,sizeof response,0))return 3;
 if(response.magic!=FMP_MAGIC||response.version!=FMP_ABI||response.epoch!=r->epoch||response.run_token!=r->run_token||response.request_id!=r->request_id
  ||response.ld_bytes!=sizeof(long double)||response.ld_mantissa!=LDBL_MANT_DIG)return 4;
 *out=response;return response.status?100+(int)response.status:0;
}
int full_moment_provider_close(FullMomentProvider *s){
 if(!s||s->pid<=0)return 1;
 close(s->request_fd);close(s->response_fd);int status;
 while(waitpid(s->pid,&status,0)<0){if(errno!=EINTR)return 2;}
 s->pid=0;return WIFEXITED(status)&&WEXITSTATUS(status)==0?0:3;
}
