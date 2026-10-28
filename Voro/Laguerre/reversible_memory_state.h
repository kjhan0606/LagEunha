#ifndef LAGFORCE_REVERSIBLE_MEMORY_STATE_H
#define LAGFORCE_REVERSIBLE_MEMORY_STATE_H
#include <stdint.h>
#include <math.h>
#ifdef __CUDACC__
#define RMS_INLINE static __host__ __device__ __forceinline__
#else
#define RMS_INLINE static inline
#endif

/* Dedicated auxiliary state: never alias thermal energy, entropy or weights.
 * Persistent identity is particle id; the caller migrates this record by id.
 * Checkpoints carry id,z only when phase==0. Trial fields are reconstructed.
 * This header does not alter the production particle/snapshot layout.
 */
typedef struct {
    uint64_t id;
    double z[2], start[2], k[4][2], dt;
    int phase; /* 0 accepted, 1..4 awaiting that RHS, 5 trial endpoint */
} rm_owner_state;

RMS_INLINE int rms_init(rm_owner_state *s,uint64_t id,double zx,double zy) {
    if(!isfinite(zx) || !isfinite(zy))return 1;
    s->id=id;s->z[0]=zx;s->z[1]=zy;s->dt=0;s->phase=0;
    for(int a=0;a<2;a++) {
        s->start[a]=s->z[a];
        for(int j=0;j<4;j++)s->k[j][a]=0;
    }
    return 0;
}
RMS_INLINE int rms_begin(rm_owner_state *s,double dt) {
    if(s->phase!=0 || !(dt>0) || !isfinite(dt))return 1;
    for(int a=0;a<2;a++) {
        if(!isfinite(s->z[a]))return 1;
        s->start[a]=s->z[a];
        for(int j=0;j<4;j++)s->k[j][a]=0;
    }
    s->dt=dt;s->phase=1;return 0;
}
/* Supply RHS evaluated at the current stage's x,v,z and acoustic coefficient.
 * On stages 1..3 expose the next z stage; on stage 4 expose trial endpoint.
 * A false/nonfinite input leaves state untouched.
 */
RMS_INLINE int rms_stage(rm_owner_state *s,int stage,double rx,double ry) {
    if(stage<1 || stage>4 || s->phase!=stage || !isfinite(rx) || !isfinite(ry))return 1;
    double rate[2]={rx,ry},next[2],increment[2];
    for(int a=0;a<2;a++) {
        increment[a]=s->dt*rate[a];
        if(stage<4)next[a]=s->start[a]+(stage==3?1.:.5)*increment[a];
        else next[a]=s->start[a]+(s->k[0][a]+2*s->k[1][a]+2*s->k[2][a]+increment[a])/6.;
        if(!isfinite(increment[a]) || !isfinite(next[a]))return 1;
    }
    for(int a=0;a<2;a++){s->k[stage-1][a]=increment[a];s->z[a]=next[a];}
    s->phase=stage+1;return 0;
}
RMS_INLINE int rms_accept(rm_owner_state *s) {
    if(s->phase!=5)return 1;
    /* Clear scratch deterministically: restart and uninterrupted records agree. */
    return rms_init(s,s->id,s->z[0],s->z[1]);
}
RMS_INLINE int rms_reject(rm_owner_state *s) {
    if(s->phase<1 || s->phase>5)return 1;
    return rms_init(s,s->id,s->start[0],s->start[1]);
}
RMS_INLINE double rms_energy(const rm_owner_state *s) {
    return .5*(s->z[0]*s->z[0]+s->z[1]*s->z[1]);
}
/* Checkpoint representation is defined fieldwise, not fwrite(struct).
 * Caller must verify matching hydro snapshot time, particle ID set and hash,
 * version/endianness, file length, and atomic pair publication before loading.
 */
RMS_INLINE int rms_pack(const rm_owner_state *s,uint64_t *id,double values[2]) {
    if(s->phase!=0 || !isfinite(s->z[0]) || !isfinite(s->z[1]))return 1;
    *id=s->id;values[0]=s->z[0];values[1]=s->z[1];return 0;
}
#undef RMS_INLINE
#endif
