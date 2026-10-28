#ifndef LAGFORCE_REVERSIBLE_MEMORY_CHECKPOINT_H
#define LAGFORCE_REVERSIBLE_MEMORY_CHECKPOINT_H
/* Host-only sidecar codec. The caller supplies a fingerprint of the associated
 * hydro snapshot and its exact time. This codec does not compute that hash,
 * publish the hydro snapshot or migrate distributed particle state.
 * Records are sorted by unique persistent ID. Existing files are not replaced.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "reversible_memory_state.h"
static inline uint64_t rms_hash_bytes(uint64_t hash,const void *data,size_t n) {
    const unsigned char *p=(const unsigned char*)data;
    for(size_t i=0;i<n;i++){hash^=p[i];hash*=UINT64_C(1099511628211);}
    return hash;
}
static inline uint64_t rms_record_hash(uint64_t hash,uint64_t id,const double z[2]) {
    hash=rms_hash_bytes(hash,&id,sizeof(id));
    return rms_hash_bytes(hash,z,2*sizeof(double));
}
static inline int rms_checkpoint_write(const char *path,const rm_owner_state *s,
        size_t n,uint64_t hydro_hash,double time) {
    if(!n || !isfinite(time) || sizeof(double)!=8)return 1;
    uint64_t hash=UINT64_C(14695981039346656037);
    for(size_t i=0;i<n;i++) {
        uint64_t id;double z[2];
        if(rms_pack(s+i,&id,z) || (i && id<=s[i-1].id))return 1;
        hash=rms_record_hash(hash,id,z);
    }
    char temporary[4096];
    int length=snprintf(temporary,sizeof(temporary),"%s.tmp.%ld",path,(long)getpid());
    if(length<0 || (size_t)length>=sizeof(temporary))return 1;
    FILE *f=fopen(temporary,"wx");if(!f)return 1;
    uint64_t header[6]={UINT64_C(0x524d454d53544154),1,
        UINT64_C(0x0102030405060708),(uint64_t)n,hydro_hash,hash};
    int bad=fwrite(header,sizeof(header),1,f)!=1 || fwrite(&time,sizeof(time),1,f)!=1;
    for(size_t i=0;i<n && !bad;i++) {
        bad=fwrite(&s[i].id,sizeof(s[i].id),1,f)!=1 || fwrite(s[i].z,sizeof(double),2,f)!=2;
    }
    if(fclose(f)!=0)bad=1;
    if(bad)return 1;
    /* Atomic no-replace publication on the same filesystem. Interrupted temp
     * files remain distinguishable from accepted sidecars. */
    if(link(temporary,path)!=0)return 1;
    (void)unlink(temporary);return 0;
}
static inline int rms_checkpoint_read(const char *path,rm_owner_state *s,size_t n,
        uint64_t hydro_hash,double time) {
    if(!n || !isfinite(time) || sizeof(double)!=8 || n>SIZE_MAX/sizeof(*s))return 1;
    FILE *f=fopen(path,"rb");if(!f)return 1;
    uint64_t header[6];double saved_time;
    int bad=fread(header,sizeof(header),1,f)!=1 || fread(&saved_time,sizeof(saved_time),1,f)!=1;
    if(bad || header[0]!=UINT64_C(0x524d454d53544154) || header[1]!=1
       || header[2]!=UINT64_C(0x0102030405060708) || header[3]!=(uint64_t)n
       || header[4]!=hydro_hash || saved_time!=time){fclose(f);return 1;}
    rm_owner_state *copy=(rm_owner_state*)malloc(n*sizeof(*copy));
    if(!copy){fclose(f);return 1;}
    uint64_t hash=UINT64_C(14695981039346656037);
    for(size_t i=0;i<n && !bad;i++) {
        uint64_t id;double z[2];
        bad=fread(&id,sizeof(id),1,f)!=1 || fread(z,sizeof(double),2,f)!=2;
        if(bad)break;
        if(id!=s[i].id || (i && id<=copy[i-1].id) || rms_init(copy+i,id,z[0],z[1])){bad=1;break;}
        hash=rms_record_hash(hash,id,z);
    }
    if(!bad && (hash!=header[5] || fgetc(f)!=EOF || ferror(f)))bad=1;
    if(fclose(f)!=0)bad=1;
    if(!bad)memcpy(s,copy,n*sizeof(*s));
    free(copy);return bad?1:0;
}
#endif
