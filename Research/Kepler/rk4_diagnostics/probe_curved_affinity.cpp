#include <mpi.h>
#include <omp.h>
#include <sched.h>
#include <cstdio>
int main(int argc,char **argv){int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int rank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    #pragma omp parallel
    {cpu_set_t mask;CPU_ZERO(&mask);sched_getaffinity(0,sizeof(mask),&mask);
        #pragma omp critical
        printf("rank=%d thread=%d cpu=%d allowed=%d proc_bind=%d places=%d\n",rank,omp_get_thread_num(),sched_getcpu(),CPU_COUNT(&mask),(int)omp_get_proc_bind(),omp_get_num_places());}
    MPI_Finalize();return 0;}
