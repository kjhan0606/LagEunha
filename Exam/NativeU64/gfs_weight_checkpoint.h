#ifndef GFS_WEIGHT_CHECKPOINT_H
#define GFS_WEIGHT_CHECKPOINT_H
#include <stdint.h>
#include <stdio.h>
#include <stddef.h>

/* Deliberately distinct from the legacy -1 format. The current original
 * unpatched reader does NOT have an unknown-negative-magic guard. */
#define GFS_WEIGHT_SNAPSHOT_MAGIC INT32_C(-2147483607)
#define GFS_WEIGHT_SNAPSHOT_VERSION UINT32_C(2)
#define GFS_WEIGHT_FEATURE_INERTIA UINT64_C(1)
#define GFS_WEIGHT_FEATURE_FROZEN UINT64_C(2)
#define GFS_WEIGHT_FEATURE_THERMAL64 UINT64_C(4)
#define GFS_WEIGHT_LAYOUT_TAG UINT64_C(0x4746535736393632)
#define GFS_WEIGHT_PREFIX_BYTES 576u
#define GFS_WEIGHT_RECORD_BYTES 696u
#define GFS_WEIGHT_FOOTER_BYTES 52u
#define GFS_WEIGHT_ERROR_BYTES 256u

typedef struct {
    int32_t magic;
    uint32_t version;
    uint64_t count;
    uint32_t record_bytes,prefix_bytes;
    uint64_t features,layout_tag,payload_bytes;
    uint32_t header_bytes,endian_tag;
    uint64_t body_fnv64;
} GfsWeightSnapshotHeader;

/* Include eunha.h before this header when GFS_WEIGHT_DYNAMICS is enabled. */
int gfs_weight_particle_valid(const treevorostressrk4particletype*,uint64_t);
int gfs_weight_snapshot_begin(FILE*,uint64_t,uint64_t,char*);
int gfs_weight_snapshot_finish(FILE*,char*);
int gfs_weight_snapshot_write(FILE*,const treevorostressrk4particletype*,uint64_t,
                              uint64_t,const unsigned char*,char*);
int gfs_weight_snapshot_read(FILE*,GfsWeightSnapshotHeader*,
                             treevorostressrk4particletype**,unsigned char*,char*);
#endif
