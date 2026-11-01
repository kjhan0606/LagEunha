#ifndef GFS_WALL_VALIDATION_SCAN_H
#define GFS_WALL_VALIDATION_SCAN_H

#include <stddef.h>

#if !defined(GFS_WALL_SCAN_PARTICLE_TYPE) || !defined(GFS_WALL_SCAN_ID) || \
    !defined(GFS_WALL_SCAN_MASS) || !defined(GFS_WALL_SCAN_STATE_VALID)
#error "Define the wall-scan particle type, ID, mass, and state accessors before including this header."
#endif

static int gfs_wall_parallel_validation_scan(ptrdiff_t n,size_t extent,
        unsigned int generation,GfsWallBase *base,
        GFS_WALL_SCAN_PARTICLE_TYPE *particles,int rhs){
    int invalid=0;
    #pragma omp parallel for schedule(static) reduction(|:invalid)
    for(ptrdiff_t i=0;i<n;i++){
        GFS_WALL_SCAN_PARTICLE_TYPE *p=particles+i;
        size_t id=(size_t)GFS_WALL_SCAN_ID(p);
        if(id>=extent || !base[id].present || GFS_WALL_SCAN_MASS(p)!=base[id].m){
            invalid=1;continue;
        }
        unsigned int previous;
        #pragma omp atomic capture
        {previous=base[id].seen;base[id].seen=generation;}
        if(previous==generation || !GFS_WALL_SCAN_STATE_VALID(p,rhs))invalid=1;
    }
    return !invalid;
}

#endif
