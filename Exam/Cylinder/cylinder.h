typedef HydroTreeLinkedCell CellType;

treevorork4particletype *cyl_mkinitial(SimParameters *, int *);
void paddingCylinderParticles(SimParameters *, postype);
int cyl_makemap(SimParameters *, int);
void cyl_outdata(SimParameters *, int, postype, postype);
void cyl_readdata(SimParameters *, postype *, postype *, int);
int cyl_evolBP(treevorork4particletype *, postype, postype);

#define XOFP(simpar,bp) (bp->x)
#define YOFP(simpar,bp) (bp->y)
