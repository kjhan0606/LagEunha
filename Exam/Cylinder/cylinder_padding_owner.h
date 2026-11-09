#ifndef CYLINDER_PADDING_OWNER_H
#define CYLINDER_PADDING_OWNER_H
/* In the padding array's actual order. source_index is local+MPI-source
 * order at generation time; source_id must be checked before dereferencing.
 * No assumption that a remote source is locally owned is encoded here. */
enum { CYL_IMAGE_MPI=0, CYL_IMAGE_RADIAL=1, CYL_IMAGE_PERIODIC=2,
       CYL_IMAGE_INFLOW=3, CYL_IMAGE_OUTFLOW=4 };
typedef struct {
    long long source_id;
    int source_index, kind, y_shift;
} CylinderPaddingOwner;
#endif
