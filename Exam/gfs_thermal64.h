#ifndef GFS_THERMAL64_H
#define GFS_THERMAL64_H
/* Opt-in thermal state in otherwise unused smooth-av5 double payload slots.
 * Native particle ABI and full-stride checkpoint format remain unchanged.
 * These slots conflict with MFM, entropy repair and DIRECT_U: dispatch guards
 * must refuse those combinations and refuse tagged state in legacy mode.
 */
#include <stdlib.h>
#include <math.h>
#define GFS_T64_TAG 2026102603.0
#define GFS_T64_U(p) ((p).stress.E_inv_xx)
#define GFS_T64_V(p) ((p).stress.E_inv_xy)
#define GFS_T64_RHS(p) ((p).stress.E_inv_yx)
#define GFS_T64_TAGGED(p) ((p).stress.E_inv_yy==GFS_T64_TAG)
#define GFS_T64_ANCHOR(p) ((p).stress.h_mfm)
/* Optional process-constant environment cache. When enabled, callers must
 * keep GFS_NATIVE_CONST_CONFIG, GFS_THERMAL64, EUNHA_IC and the Kepler center/
 * softening settings fixed for this process. Default-off preserves the legacy
 * dynamic environment reads. */
static inline int gfs_native_const_config_on(void){
    static int mode=-1;
    int value=__atomic_load_n(&mode,__ATOMIC_RELAXED);
    if(value<0){
        const char *s=getenv("GFS_NATIVE_CONST_CONFIG");
        value=s && atoi(s)>0;
        __atomic_store_n(&mode,value,__ATOMIC_RELAXED);
    }
    return value;
}
static inline int gfs_thermal64_requested(void){
    if(!gfs_native_const_config_on()){
        const char *s=getenv("GFS_THERMAL64");
        return s && atoi(s)>0;
    }
    static int cached=-1;
    int value=__atomic_load_n(&cached,__ATOMIC_RELAXED);
    if(value<0){
        const char *s=getenv("GFS_THERMAL64");
        value=s && atoi(s)>0;
        __atomic_store_n(&cached,value,__ATOMIC_RELAXED);
    }
    return value;
}
#ifndef __CUDACC__
static inline double gfs_t64_U(const treevorostressrk4particletype *p){
    return gfs_thermal64_requested() && GFS_T64_TAGGED(*p)?GFS_T64_U(*p):(double)p->ie;
}
static inline double gfs_t64_V(const treevorostressrk4particletype *p){
    return gfs_thermal64_requested() && GFS_T64_TAGGED(*p)?GFS_T64_V(*p):(double)p->volume;
}
static inline double gfs_t64_Edot(const treevorostressrk4particletype *p){
    return gfs_thermal64_requested() && GFS_T64_TAGGED(*p)?
        GFS_T64_RHS(*p)+p->mass*(p->vx*p->ax+p->vy*p->ay):(double)p->die;
}
#endif
#endif
