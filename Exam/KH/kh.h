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

void kh_MkLinkedList(SimParameters*);
Voro2D_point *kh_Voro2D_FindNeighbor(SimParameters *, int , int , int *);
treevorork4particletype *kh_Voro2D_FindCellBP(SimParameters *,int,int,int *);
Voro2D_point *Voro2D_FindNeighbor(SimParameters *, int , int , int *, 
		treevorork4particletype *);
treevorork4particletype *Voro2D_FindCellBP(SimParameters *, int , int , 
		int *, treevorork4particletype *);
double kh_vph2D_rk4(SimParameters *);
double kh_vph2D(SimParameters *);
treevorork4particletype *kh_mkinitial(SimParameters *, int *);
void KH_TreeAllParticlePadding(SimParameters *, float);
void KH_TreeAllParticleMigrate(SimParameters *simpar);
int kh_makemap(SimParameters *, int);
void KH_AllParticlePadding(SimParameters *, postype );
double kh_w2Measure(SimParameters *, double , double ,double );

#define XOFP(simpar,bp) (bp->x)
#define YOFP(simpar,bp) (bp->y)
