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

void rt_MkLinkedList(SimParameters*);
Voro2D_point *rt_Voro2D_FindNeighbor(SimParameters *, int , int , int *, treevoroparticletype *);
treevoroparticletype *rt_Voro2D_FindCellBP(SimParameters *, int , int , int *, treevoroparticletype *);
Voro2D_point *Voro2D_FindNeighbor(SimParameters *, int , int , int *, treevoroparticletype *);
treevoroparticletype *Voro2D_FindCellBP(SimParameters *, int , int , int *, treevoroparticletype *);
double rt_vph2D(SimParameters *);
treevoroparticletype *rt_mkinitial(SimParameters *, int *);
void RT_TreeAllParticlePadding(SimParameters *, float);
void RT_TreeAllParticleMigrate(SimParameters *simpar);
int rt_makemap(SimParameters *, int);

#define XOFP(simpar,bp) (bp->x)
#define YOFP(simpar,bp) (bp->y)
