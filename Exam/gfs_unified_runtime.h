#ifndef GFS_UNIFIED_RUNTIME_H
#define GFS_UNIFIED_RUNTIME_H
#include <stdlib.h>
#include <string.h>
static inline int gfs_unified_periodic_mode(void){
    const char *method=getenv("GFS_UNIFIED_HYDRO");
    const char *boundary=getenv("GFS_UNIFIED_BOUNDARY");
    return method&&!strcmp(method,"1")&&boundary&&!strcmp(boundary,"periodic");
}
#endif
