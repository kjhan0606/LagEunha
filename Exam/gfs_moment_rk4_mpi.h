#ifndef GFS_MOMENT_RK4_MPI_H
#define GFS_MOMENT_RK4_MPI_H
#include <mpi.h>
#include <stddef.h>
#include <float.h>
#include "gfs_moment_rk4.h"

/* All ranks must enter in the same order, including ranks with no cells.
 * Returns 1 only on unanimous acceptance, 0 on physical rejection, and -1 on
 * an MPI error when the communicator is configured to return errors.
 * This synchronizes only these moment arrays. An integrated solver must
 * separately stage/roll back geometry, face registers, heat and all other
 * coupled state, and vote consistently at every stage before collectives. */
static inline int gfs_moment_rk4_collective_commit(
    MPI_Comm comm,GfsMomentRK4 *scratch,GfsMomentState *accepted,size_t n,int local_stage_ok){
    int local=local_stage_ok!=0;
    double local_dt_min=DBL_MAX,local_dt_max=0;
    for(size_t i=0;local && i<n;i++){
        if(scratch[i].phase!=4){local=0;break;}
        if(!(scratch[i].dt>0) || !isfinite(scratch[i].dt)){local=0;break;}
        if(scratch[i].dt<local_dt_min)local_dt_min=scratch[i].dt;
        if(scratch[i].dt>local_dt_max)local_dt_max=scratch[i].dt;
        for(int q=0;q<GFM_NSTATE;q++)if(scratch[i].base.a[q]!=accepted[i].a[q])local=0;
        if(!local)break;
        GfsMomentPrimitive primitive;
        if(!gfs_moment_primitive(&scratch[i].stage,&primitive))local=0;
        if(scratch[i].fixed_mass && scratch[i].stage.a[GFM_M]!=scratch[i].base.a[GFM_M])local=0;
    }
    double global_dt_min,global_dt_max;
    if(MPI_Allreduce(&local_dt_min,&global_dt_min,1,MPI_DOUBLE,MPI_MIN,comm)!=MPI_SUCCESS ||
       MPI_Allreduce(&local_dt_max,&global_dt_max,1,MPI_DOUBLE,MPI_MAX,comm)!=MPI_SUCCESS)return -1;
    if(global_dt_min!=global_dt_max || global_dt_max==0)local=0;
    int global=0;
    if(MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,comm)!=MPI_SUCCESS)return -1;
    if(!global)return 0;
    for(size_t i=0;i<n;i++)accepted[i]=scratch[i].stage;
    return 1;
}
#endif
