typedef HydroTreeLinkedCell CellType;

treevorork4particletype *sedov2d_mkinitial(SimParameters *, int *);
int sedov2d_makemap(SimParameters *, int);

#define XOFP(simpar,bp) (bp->x)
#define YOFP(simpar,bp) (bp->y)
