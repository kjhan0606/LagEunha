#define Nx 64
#define Ny 64
#define Nximg 256
#define Nyimg 256
#define Nxp 256
#define Nyp 256


#define Lx (1.L)
#define Ly (1.L)
#define boxsize (Lx/Nx*5.L)


#define Gamma 1.666666L
#define Zeta (0.125L/100.L)

#define Courant 0.2L
#define alphavis 0.2
#define betavis 0.0

typedef struct CellType{
	int np;
	Voro2D_GasParticle *bp;
} CellType;

void MkLinkedList(Voro2D_GasParticle *, int );
Voro2D_point *Voro2D_FindNeighbor(int , int , int *, Voro2D_GasParticle *);
Voro2D_GasParticle *Voro2D_FindCellBP(int , int , int *, Voro2D_GasParticle *);
double vph2D(Voro2D_GasParticle *, int );
Voro2D_GasParticle *mkinitial(int *);

