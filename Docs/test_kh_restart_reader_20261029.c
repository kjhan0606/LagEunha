/* Exercise the actual kh_readdata with singleton MPI and isolated fixtures. */
#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "eunha.h"
#include "voro.h"
#include "Exam/gfs_moment.h"
#include "Exam/gfs_thermal64.h"

void kh_readdata(SimParameters *,postype *,postype *,int);
void migrateTreeVorork4Particles(SimParameters *s){(void)s;}
void *my_malloc(size_t n){void *p=malloc(n);if(!p) abort();return p;}
void my_free(void *p){free(p);}

int main(int argc,char **argv){
    if(argc!=2) return 2;
    int old=!strcmp(argv[1],"old");
    int legacy=old || !strcmp(argv[1],"legacy"), mismatch=!strcmp(argv[1],"mismatch");
    int thermal=!strcmp(argv[1],"thermal"),moment=!strcmp(argv[1],"moment");
    int trunc_header=!strcmp(argv[1],"trunc_header"),trunc_record=!strcmp(argv[1],"trunc_record");
    int trunc_footer=!strcmp(argv[1],"trunc_footer"),tiny=!strcmp(argv[1],"tiny");
    int huge=!strcmp(argv[1],"huge"),badflag=!strcmp(argv[1],"badflag");
    if(!legacy && !mismatch && !thermal && !moment && !trunc_header && !trunc_record &&
       !trunc_footer && !tiny && !huge && !badflag && strcmp(argv[1],"plain")) return 2;
    char dir[]="/tmp/kh_guard_reader_XXXXXX";
    if(!mkdtemp(dir) || chdir(dir)) return 3;
    FILE *f=fopen("khout.000001.dat","wb");if(!f)return 3;
    int magic=-1,np=1,has_stress=!legacy;
    int stride=legacy?(int)sizeof(treevorork4particletype):(int)sizeof(treevorostressrk4particletype);
    if(mismatch)stride++;
    if(tiny){stride=1;has_stress=0;}
    if(huge)np=2000000000;
    if(badflag)has_stress=2;
    treevorostressrk4particletype record;
    memset(&record,0,sizeof(record));
    record.mass=2.5;record.ie=0.125;
    if(thermal){record.stress.E_inv_yy=GFS_T64_TAG;record.stress.E_inv_xx=0.123456789012345;}
    if(moment){GfsMomentPayload m={0};m.tag=GFS_MOMENT_TAG;memcpy(&record.rk4,&m,sizeof(m));}
    if(!old)fwrite(&magic,sizeof(magic),1,f);
    fwrite(&np,sizeof(np),1,f);
    if(!old && !trunc_header){fwrite(&stride,sizeof(stride),1,f);fwrite(&has_stress,sizeof(has_stress),1,f);}
    if(!trunc_header){
        size_t bytes=legacy?sizeof(treevorork4particletype):sizeof(record);
        if(trunc_record)bytes/=2;
        fwrite(&record,bytes,1,f);
        if(mismatch)fputc(0,f);
    }
    postype expected_t=1.25,expected_dt=0.01;float expected_old=0.02f;
    if(!trunc_header && !trunc_record){
        fwrite(&expected_t,sizeof(expected_t),1,f);
        if(!trunc_footer){fwrite(&expected_dt,sizeof(expected_dt),1,f);fwrite(&expected_old,sizeof(expected_old),1,f);}
    }
    fclose(f);
    MPI_Init(&argc,&argv);
    SimParameters sim;
    memset(&sim,0,sizeof(sim));
    MPI_COMM(&sim)=MPI_COMM_WORLD;MYID(&sim)=0;NID(&sim)=1;
    postype t=0,dt=0;
    kh_readdata(&sim,&t,&dt,1);
    treevorork4particletype *loaded=VORORK4_TBP(&sim);
    if(VORO_NP(&sim)!=1 || t!=expected_t || dt!=expected_dt || GAS_dtold(&sim)!=expected_old ||
       loaded->mass!=record.mass || loaded->ie!=record.ie) return 4;
#ifdef USE_CUDA
    if(thermal){
        treevorostressrk4particletype *full=(void*)loaded;
        if(!GFS_T64_TAGGED(*full) || GFS_T64_U(*full)!=GFS_T64_U(record)) return 5;
    }
    if(moment){
        treevorostressrk4particletype *full=(void*)loaded;
        GfsMomentPayload m;memcpy(&m,&full->rk4,sizeof(m));
        if(m.tag!=GFS_MOMENT_TAG) return 6;
    }
#endif
    my_free(loaded);MPI_Finalize();
    printf("reader accepted %s fixture\n",argv[1]);
    return 0;
}
