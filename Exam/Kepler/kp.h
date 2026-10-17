/*
#define Nx 64
#define Ny 64
#define Nximg 256
#define Nyimg 256
#define Nxp 256
#define Nyp 256
*/


/*
#define Lx (1.L)
#define Ly (1.L)
#define boxsize (Lx/Nx*5.L)
*/


#define Gamma 1.4L
#define Zeta (0.125L/100.L)

#define Courant 0.2L
/*
#define alphavis 0.2
#define betavis 0.0
*/


typedef HydroTreeLinkedCell CellType;

void kp_MkLinkedList(SimParameters*);
Voro2D_point *kp_Voro2D_FindNeighbor(SimParameters *, int , int , int *, treevorork4particletype *);
treevorork4particletype *kp_Voro2D_FindCellBP(SimParameters *, int , int , int *, treevorork4particletype *);
Voro2D_point *Voro2D_FindNeighbor(SimParameters *, int , int , int *, treevorork4particletype *);
treevorork4particletype *Voro2D_FindCellBP(SimParameters *, int , int , int *, treevorork4particletype *);
double kp_vph2D(SimParameters *);
double kp_vph2D_rk4(SimParameters *);
treevorork4particletype *kp_mkinitial(SimParameters *, int *);
treevorork4particletype *old_kp_mkinitial(SimParameters *, int *);
void KP_TreeAllParticlePadding(SimParameters *, float);
void KP_TreeAllParticleMigrate(SimParameters *simpar);
int kp_makemap(SimParameters *, int);

#define XOFP(simpar,bp) (bp->x)
#define YOFP(simpar,bp) (bp->y)
