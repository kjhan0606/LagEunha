typedef HydroTreeLinkedCell CellType;

treevorork4particletype *cyl_mkinitial(SimParameters *, int *);
void paddingCylinderParticles(SimParameters *, postype);
/* Optional ownership sidecar. source_index indexes the local particles followed
 * by the filtered MPI padding prefix, never the newly generated wall ghosts.
 * source_id preserves the source ID before the legacy ghost sentinel is set.
 * y_shift counts additional periodic box translations applied to that source.
 * Remote sources remain remote; no local-owner assumption is made here. */
#include "cylinder_padding_owner.h"
/* Returns a newly malloc-allocated sidecar in padding-array order. Caller frees
 * it. Valid only for the padding produced by this call, before any reordering.
 * This is metadata only, not an activation of material curved-wall dynamics. */
void paddingCylinderParticlesMapped(SimParameters *, postype,
                                   CylinderPaddingOwner **, int *);
int cyl_makemap(SimParameters *, int);
void cyl_outdata(SimParameters *, int, postype, postype);
void cyl_readdata(SimParameters *, postype *, postype *, int);
int cyl_evolBP(treevorork4particletype *, postype, postype);
int cyl_skipBP(treevorork4particletype *);
void cyl_postStage(SimParameters *);
void cyl_prepareWallContext(SimParameters *);
void cyl_outflowCleanup(SimParameters *);

#define XOFP(simpar,bp) (bp->x)
#define YOFP(simpar,bp) (bp->y)
