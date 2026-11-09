#ifndef GFS_UNIFIED_RUNTIME_H
#define GFS_UNIFIED_RUNTIME_H
#include <stdlib.h>
#include <string.h>
static inline int gfs_unified_periodic_mode(void){
    const char *method=getenv("GFS_UNIFIED_HYDRO");
    const char *boundary=getenv("GFS_UNIFIED_BOUNDARY");
    return method&&!strcmp(method,"1")&&boundary&&!strcmp(boundary,"periodic");
}
/* Allow periodic images for the native variational K1 inspection without
 * switching on the separate unified force model. Evolution is still gated. */
static inline int gfs_native_periodic_rhs_mode(void){
    const char *p=getenv("GFS_WEIGHT_PERIODIC_RHS_ONLY");
    const char *r=getenv("GFS_WEIGHT_NATIVE_RHS_ONLY");
    const char *k=getenv("GFS_WALL_RHS_PROBE_ONLY");
    return p&&!strcmp(p,"1")&&r&&!strcmp(r,"1")&&k&&!strcmp(k,"1");
}
static inline int gfs_native_periodic_weight_kdk_mode(void){
    const char *p=getenv("GFS_NATIVE_WEIGHT_KDK_VALIDATE");
    return p&&!strcmp(p,"1");
}
static inline int gfs_native_periodic_weight_mode(void){
    return gfs_native_periodic_rhs_mode() || gfs_native_periodic_weight_kdk_mode();
}
#endif
