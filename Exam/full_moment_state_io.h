#ifndef EXAM_FULL_MOMENT_STATE_IO_H
#define EXAM_FULL_MOMENT_STATE_IO_H
#include <stdint.h>
#include <stddef.h>
#define FM_STATE_MAGIC UINT64_C(0x31544154534d464c)
typedef struct {
 uint64_t magic,version,cells,epoch,run_token,state_elements,ld_bytes,ld_mantissa;
 double time;uint64_t reserved[7];
} FullMomentStateHeader;
_Static_assert(sizeof(FullMomentStateHeader)==128,"Native state header");
typedef struct {
 FullMomentStateHeader header;
 int64_t *ids;double *mass,*state,*factors,*weight_rate_seed;long double *Psi;
 void *mapping[6];size_t bytes[6];
} FullMomentState;
int full_moment_state_load(const char*,FullMomentState*);
void full_moment_state_close(FullMomentState*);
/* Writes a new tentative state directory only; never overwrites or commits. */
int full_moment_state_write(const char*,const FullMomentState*);
/* Durable, rank0-only atomic publication into a previously absent path. */
int full_moment_state_publish(const char*,const FullMomentState*);
/* Import a previously authenticated checkpoint into a new run. The caller
 * must first verify model/physics provenance and all source file digests.
 * expected is the authenticated ORIGINAL header, not a replacement header.
 * Only epoch/run_token change; every numerical payload byte is preserved.
 * out is unchanged on failure. This does not accept a physical time step. */
int full_moment_state_import(const char *source,const char *destination,
 const FullMomentStateHeader *expected,uint64_t new_token,FullMomentState *out);
#endif
