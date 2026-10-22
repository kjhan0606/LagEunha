// Explicit preparation branch in the original native integrator lifetime.
// It never falls through to legacy particle RK scratch.
#include "material_full_reference.h"
static void gfs_material_reference_gate(SimParameters *sp,
    void (*padding)(SimParameters*,postype),
    Voro2D_point *(*neighbors)(SimParameters*,int,int,int*),
    treevorork4particletype *(*cell)(SimParameters*,int,int,int*),
    void (*linked)(SimParameters*,postype,postype,postype,postype,postype,void(*)(SimParameters*,postype)),
    void (*post)(SimParameters*)){
#ifdef USE_CUDA
    const char *raw=getenv("GFS_GPU_RAW_FACES"),*boundary=getenv("GFS_MATERIAL_BOUNDARY");
    const char *ic=getenv("EUNHA_IC");
    if(!GAS_GPU_ENABLED(sp)||gfs_rk4_context_time!=0||!raw||atoi(raw)<=0||
       !boundary||strcmp(boundary,"material_constant_Pext")||!ic||strcmp(ic,"kepler")||
       GAS_Kappa(sp)!=0||GAS_FCENTROID(sp)!=0||GAS_XSPHEPS(sp)!=0||
       GAS_VISCOSITY(sp)!=0||GAS_AlphaVis(sp)!=0||GAS_BetaVis(sp)!=0||
       GAS_CDAMAX(sp)!=0||GAS_HYPERVALPHA(sp)!=0||GAS_ENTROPY_MODE(sp)!=0||
       TVORORK4_DDINFO(sp)[0].n_size!=sizeof(treevorostressrk4particletype)){
        fprintf(stderr,"[MATERIAL_REFERENCE] incompatible fresh IC/boundary/Re-infinity contract\n");MPI_Abort(MPI_COMM(sp),99);
    }
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(sp);
    for(ptrdiff_t k=0;k<VORO_NP(sp);k++)if(p[k].w2!=0||gfs_moment_tagged(p+k))MPI_Abort(MPI_COMM(sp),99);
    const char *retained=getenv("GFS_MATERIAL_REFERENCE_INPUT");int result;
    if(retained&&retained[0]){
        result=gfs_material_initialize_full(retained,MPI_Comm_c2f(MPI_COMM(sp)),
            GAS_GAMMA(sp),1e-6,gfs_moment_env("EUNHA_KEPLER_EPS",0),SIMBOX(sp).x.max);
    }else{
        exam2dUpdateVol(sp,padding,neighbors,cell,linked);
        result=gfs_material_prepare_reference(gpu_get_context(),VORO_NP(sp),
            MPI_Comm_c2f(MPI_COMM(sp)),GAS_GAMMA(sp),1e-6,
            gfs_moment_env("EUNHA_KEPLER_EPS",0),SIMBOX(sp).x.max);
    }
    MPI_Barrier(MPI_COMM(sp));
    if(MYID(sp)==0){const char *evolve=getenv("GFS_MATERIAL_EVOLVE");if(evolve&&atoi(evolve)>0)fprintf(stderr,"[MATERIAL_FULL_EVOLUTION] result=%d; consult actual accepted-step/restart reports; t6 stability unproved\n",result);else fprintf(stderr,"[MATERIAL_INITIAL_ONLY] result=%d; reference/full initial operator only; no accepted step or t6 stability claim\n",result);}
    MPI_Finalize();exit(result);
#else
    fprintf(stderr,"[MATERIAL_REFERENCE] native CUDA required\n");MPI_Abort(MPI_COMM(sp),99);
#endif
}
