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


#define Gamma 1.666666L
#define Zeta (0.125L/100.L)

#define Courant 0.2L
/*
#define alphavis 0.2
#define betavis 0.0
*/


/*
typedef struct CellType{
	int np;
	struct linkedlisttype *bp;
} CellType;
*/
typedef HydroTreeLinkedCell CellType;

void MkLinkedList(SimParameters*);
Voro2D_point *Voro2D_FindNeighbor(SimParameters *, int , int , int *, treevoroparticletype *);
treevoroparticletype *Voro2D_FindCellBP(SimParameters *, int , int , int *, treevoroparticletype *);
double vph2D(SimParameters *);
treevoroparticletype *mkinitial(SimParameters *, int *);
void KH_TreeAllParticlePadding(SimParameters *, float);
void KH_TreeAllParticleMigrate(SimParameters *simpar);
int makemap(SimParameters *, int);
int maketscmap(SimParameters *, int);

#define XOFP(simpar,bp) (bp->x)
#define YOFP(simpar,bp) (bp->y)
