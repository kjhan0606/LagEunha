#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <sys/types.h>
#include <mpi.h>
#include "eunha.h"
#include "gfs_weight_checkpoint.h"

#ifndef GFS_WEIGHT_DYNAMICS
#error The weight codec requires the explicit native weight particle layout
#endif
_Static_assert(sizeof(GfsWeightSnapshotHeader)==64,"Versioned header layout");
_Static_assert(sizeof(GfsWeightState)==32,"Persistent weight state size");
_Static_assert(sizeof(GfsWeightRK4)==64,"Particle-owned weight RK scratch");
_Static_assert(sizeof(PosType)==8,"Double coordinate build required");
_Static_assert(offsetof(treevorostressrk4particletype,weight)==576,"Legacy prefix must survive intact");
_Static_assert(offsetof(treevorostressrk4particletype,weight_rk4)==608,"Explicit RK extension");
_Static_assert(sizeof(treevorostressrk4particletype)==696,"Native weighted particle size");
_Static_assert(offsetof(treevorostressrk4particletype,weight_thermal)==672,"Explicit thermal state offset");
_Static_assert(sizeof(GfsWeightThermal)==24,"Explicit U/Udot/V state");

static int fail(char *error,const char *message){
    if(error)snprintf(error,GFS_WEIGHT_ERROR_BYTES,"%s",message);
    return 1;
}
static uint64_t fnv(uint64_t h,const void *input,size_t n){
    const unsigned char *p=input;
    for(size_t k=0;k<n;k++){h^=p[k];h*=UINT64_C(1099511628211);}
    return h;
}
static int header_valid(const GfsWeightSnapshotHeader *h,char *error){
    if(h->magic!=GFS_WEIGHT_SNAPSHOT_MAGIC || h->version!=2 || h->header_bytes!=64 ||
       h->endian_tag!=UINT32_C(0x01020304))return fail(error,"unsupported weight checkpoint header/version/endianness");
    if(h->record_bytes!=696 || h->prefix_bytes!=576 || h->layout_tag!=GFS_WEIGHT_LAYOUT_TAG)
        return fail(error,"unsupported native weighted particle layout");
    if((h->features!=(GFS_WEIGHT_FEATURE_INERTIA|GFS_WEIGHT_FEATURE_THERMAL64)) &&
       (h->features!=(GFS_WEIGHT_FEATURE_INERTIA|GFS_WEIGHT_FEATURE_THERMAL64|GFS_WEIGHT_FEATURE_FROZEN)))
        return fail(error,"unsupported weight checkpoint feature set");
    if(!h->count || h->count>INT_MAX || h->count>SIZE_MAX/696 ||
       h->payload_bytes!=h->count*UINT64_C(696))return fail(error,"invalid weight checkpoint count/payload length");
    return 0;
}
static int footer_valid(const unsigned char *footer){
    double t,dt,lx,ly;float dtold,alpha,beta;int32_t nx,ny;
    memcpy(&t,footer,8);memcpy(&dt,footer+8,8);memcpy(&dtold,footer+16,4);
    memcpy(&nx,footer+20,4);memcpy(&ny,footer+24,4);
    memcpy(&lx,footer+28,8);memcpy(&ly,footer+36,8);
    memcpy(&alpha,footer+44,4);memcpy(&beta,footer+48,4);
    return isfinite(t)&&t>=0&&isfinite(dt)&&dt>=0&&isfinite(dtold)&&dtold>=0&&
        nx>0&&ny>0&&isfinite(lx)&&lx>0&&isfinite(ly)&&ly>0&&
        isfinite(alpha)&&alpha>=0&&isfinite(beta)&&beta>=0;
}
int gfs_weight_particle_valid(const treevorostressrk4particletype *p,uint64_t features){
    if(!(features&GFS_WEIGHT_FEATURE_THERMAL64) ||
       !isfinite(p->weight_thermal.U) || !(p->weight_thermal.U>0) ||
       !isfinite(p->weight_thermal.Udot) || !isfinite(p->weight_thermal.V) || !(p->weight_thermal.V>0) ||
       p->ie!=(float)p->weight_thermal.U || p->volume!=(float)p->weight_thermal.V)return 0;
    if(!gfs_weight_state_valid(&p->weight) || !isfinite(p->mass) || !(p->mass>0) ||
       !isfinite(p->volume) || !(p->volume>0) || !isfinite(p->ie) || !(p->ie>0) ||
       !isfinite(p->x) || !isfinite(p->y) || !isfinite(p->vx) || !isfinite(p->vy))return 0;
    /* A new dynamics checkpoint must not silently select legacy float W2. */
    if(p->w2!=0 || p->w2old!=0)return 0;
    if((features&GFS_WEIGHT_FEATURE_FROZEN) && (p->weight.W!=0 || p->weight.pi!=0))return 0;
    double stage[8];memcpy(stage,&p->weight_rk4,sizeof(stage));
    for(int k=0;k<8;k++)if(!isfinite(stage[k]))return 0;
    return 1;
}
int gfs_weight_snapshot_begin(FILE *f,uint64_t n,uint64_t features,char *error){
    GfsWeightSnapshotHeader h;memset(&h,0,sizeof(h));
    h.magic=GFS_WEIGHT_SNAPSHOT_MAGIC;h.version=2;h.count=n;
    h.record_bytes=696;h.prefix_bytes=576;h.features=features;
    h.layout_tag=GFS_WEIGHT_LAYOUT_TAG;h.payload_bytes=n*UINT64_C(696);
    h.header_bytes=64;h.endian_tag=UINT32_C(0x01020304);
    if(header_valid(&h,error))return 1;
    if(fseeko(f,0,SEEK_SET) || fwrite(&h,sizeof(h),1,f)!=1)return fail(error,"cannot write weight checkpoint header");
    return 0;
}
int gfs_weight_snapshot_finish(FILE *f,char *error){
    GfsWeightSnapshotHeader h;
    if(fflush(f) || fseeko(f,0,SEEK_SET) || fread(&h,sizeof(h),1,f)!=1)return fail(error,"cannot finalize weight checkpoint header");
    if(header_valid(&h,error))return 1;
    if(fseeko(f,0,SEEK_END))return fail(error,"cannot size weight checkpoint");
    off_t size=ftello(f);
    if(size<0 || (uint64_t)size!=64+h.payload_bytes+52)return fail(error,"incomplete weight checkpoint body/footer");
    if(fseeko(f,64,SEEK_SET))return fail(error,"cannot seek weight checkpoint body");
    unsigned char buffer[65536];uint64_t left=h.payload_bytes+52;
    uint64_t hash=UINT64_C(14695981039346656037);
    while(left){
        size_t take=left>sizeof(buffer)?sizeof(buffer):(size_t)left;
        if(fread(buffer,1,take,f)!=take)return fail(error,"short weight checkpoint checksum read");
        hash=fnv(hash,buffer,take);left-=take;
    }
    h.body_fnv64=hash;
    if(fseeko(f,0,SEEK_SET) || fwrite(&h,sizeof(h),1,f)!=1 || fflush(f))
        return fail(error,"cannot commit weight checkpoint checksum");
    return 0;
}
int gfs_weight_snapshot_write(FILE *f,const treevorostressrk4particletype *p,uint64_t n,
                              uint64_t features,const unsigned char *footer,char *error){
    if(!footer_valid(footer))return fail(error,"invalid weight checkpoint time/domain footer");
    for(uint64_t k=0;k<n;k++)if(!gfs_weight_particle_valid(p+k,features))return fail(error,"invalid native weight particle");
    if(gfs_weight_snapshot_begin(f,n,features,error))return 1;
    if(fwrite(p,sizeof(*p),(size_t)n,f)!=(size_t)n || fwrite(footer,1,52,f)!=52)
        return fail(error,"cannot write weight checkpoint payload/footer");
    return gfs_weight_snapshot_finish(f,error);
}
int gfs_weight_snapshot_read(FILE *f,GfsWeightSnapshotHeader *header,
                             treevorostressrk4particletype **out,unsigned char *footer,char *error){
    *out=NULL;
    GfsWeightSnapshotHeader h;
    if(fseeko(f,0,SEEK_SET) || fread(&h,sizeof(h),1,f)!=1)return fail(error,"short weight checkpoint header");
    if(header_valid(&h,error))return 1;
    if(fseeko(f,0,SEEK_END))return fail(error,"cannot size weight checkpoint");
    off_t size=ftello(f);
    if(size<0 || (uint64_t)size!=64+h.payload_bytes+52)return fail(error,"weight checkpoint size mismatch/truncation");
    if(fseeko(f,64,SEEK_SET))return fail(error,"cannot seek weight checkpoint payload");
    treevorostressrk4particletype *p=malloc((size_t)h.payload_bytes);
    if(!p)return fail(error,"weight checkpoint allocation failed");
    if(fread(p,sizeof(*p),(size_t)h.count,f)!=(size_t)h.count || fread(footer,1,52,f)!=52){
        free(p);return fail(error,"short weight checkpoint payload/footer");
    }
    uint64_t hash=fnv(UINT64_C(14695981039346656037),p,(size_t)h.payload_bytes);
    hash=fnv(hash,footer,52);
    if(hash!=h.body_fnv64){free(p);return fail(error,"weight checkpoint body checksum mismatch");}
    if(!footer_valid(footer)){free(p);return fail(error,"invalid weight checkpoint time/domain footer");}
    for(uint64_t k=0;k<h.count;k++)if(!gfs_weight_particle_valid(p+k,h.features)){
        free(p);return fail(error,"invalid persistent W/pi/inertia/tag/stage state");
    }
    *header=h;*out=p;return 0;
}
