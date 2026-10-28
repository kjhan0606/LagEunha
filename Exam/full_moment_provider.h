#ifndef EXAM_FULL_MOMENT_PROVIDER_H
#define EXAM_FULL_MOMENT_PROVIDER_H
#include <stdint.h>
#include <sys/types.h>

#define FMP_MAGIC UINT64_C(0x32564352534d464c)
#define FMP_ABI 2
enum FullMomentProviderCommand {FMP_CONFIGURE=1,FMP_BOUNDS=2,FMP_DENSITY=3,FMP_VELOCITY=4,FMP_RELEASE=5,FMP_STOP=6,FMP_FACE_STAGE=7,FMP_RESUME_STAGE=8,FMP_EXPORT_STATE=9,FMP_FACTOR_AUDIT=10,FMP_GEOMETRY=11,FMP_ADMISSIBILITY=12,FMP_DISPOSE_GEOMETRY=13,FMP_VERIFY_RESTART=14,FMP_COUPLED_METADATA=15,FMP_VERIFY_FORK=16,FMP_MIDPOINT_OBSERVE=17,FMP_COUPLED_JOURNAL=18};
/* Optional coupled-method metadata export: actual fit scales/geometric areas,
 * no provider physics evaluation or state update. ABI packet sizes unchanged. */
/* VERIFY_RESTART is read-only evidence verification before first CONFIGURE.
 * reserved[0:8] holds the requested native IEEE binary64 timestep (memcpy).
 * bounds_path=manifest, fits_path=original IC, aux_path=physics manifest.
 * Returned directory contains authenticated expected_header.bin + receipt. */
/* FMP_FACE_STAGE returns this validated event before awaiting RESUME_STAGE.
 * It requests the native owned GPU stage, not permission or a fallback. */
#define FMP_NATIVE_STAGE_REQUIRED 2
#pragma pack(push,1)
typedef struct {
 uint64_t magic,epoch,run_token,request_id;
 uint32_t version,command,step,stage,outer,rank,nranks,ncells,ld_bytes,ld_mantissa;
 char physics_sha[64];
 char state_path[512],geometry_path[512],bounds_path[512],fits_path[512],output_dir[512],aux_path[512];
 unsigned char reserved[888];
} FullMomentProviderRequest;
typedef struct {
 uint64_t magic,epoch,run_token,request_id;
 uint32_t version,status,ncells,ld_bytes,ld_mantissa,reserved;
 uint64_t payload_bytes,owned_start,owned_count,gravity_nodes,token;
 char output_path[512],message[256];unsigned char padding[160];
} FullMomentProviderResponse;
#pragma pack(pop)
_Static_assert(sizeof(FullMomentProviderRequest)==4096,"Provider request wire size");
_Static_assert(sizeof(FullMomentProviderResponse)==1024,"Provider response wire size");
typedef struct {pid_t pid;int request_fd,response_fd;uint64_t sequence,run_token;} FullMomentProvider;
int full_moment_provider_start(FullMomentProvider*,const char *python,const char *script,int workers,const char *log);
int full_moment_provider_call(FullMomentProvider*,FullMomentProviderRequest*,FullMomentProviderResponse*);
int full_moment_provider_close(FullMomentProvider*);
#endif
