#ifndef LAGFORCE_MOMENT_OWNED_BACKEND_H
#define LAGFORCE_MOMENT_OWNED_BACKEND_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* FullMomentKepler bridge, independent of legacy particle storage.
 * registry_index denotes the stable-ID-sorted global registry ROW, not the
 * persistent particle ID value. Node owner rows keep their global meaning.
 * Only authoritative owned cells are advanced. Rank boundaries are not walls.
 * All arrays are contiguous binary64/int64. Counts are element counts.
 * Nonzero status leaves every caller output unchanged.
 */
int lagmoment_canonical_nodes(int backend,int threads,int nodes,
    const double *node_inputs_33,double *node_flux_19);
int lagmoment_select_visible_device(int ordinal,int *visible_count,
    char *name_256,unsigned char *uuid_16);
/* Canonical global order, using the already computed node flux. Does not
 * compute any node primitive or cell stage again. Transactional output. */
int lagmoment_ordered_ledgers(int backend,int threads,int cells,int nodes,
    const double *canonical_node_flux_19,const int64_t *global_node_owners_2,
    const int64_t *physical_frozen,const double *global_stage_94,
    const double *global_hydro_106,double *ledgers_12);
int lagmoment_owned_stage(int backend,int threads,uint64_t epoch,uint64_t expected_epoch,
    int registry_count,int owned_cells,int nodes,int64_t incidences,int64_t gravity_nodes,
    const int64_t *owned_registry_index,const int64_t *node_global_index,
    const double *node_inputs_33,const double *canonical_node_flux_19,
    const double *owned_cells_19,const int64_t *global_node_owners_2,
    const int64_t *owned_offsets,const int64_t *signed_incidence,
    const int64_t *physical_frozen,const double *node_extra_20,
    const double *owned_factors_4,const double *owned_gravity_nodes_7,
    const int64_t *owned_gravity_offsets,double *owned_stage_94,
    double *ordered_incoming_16,double *owned_hydro_106);
#ifdef __cplusplus
}
#endif
#endif
