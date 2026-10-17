#define Nx 16
#define Ny 16
#define Nximg 64
#define Nyimg 64
#define Nxp 64
#define Nyp 64
#define Nzp 64


#define Lx (1.L)
#define Ly (1.L)
#define boxsize (Lx/Nx*5.L)


#define Gamma 1.666666L
#define Zeta (0.125L/100.L)

#define Courant 0.1L
#define alphavis 0.2
#define betavis 0.4

typedef struct CellType{
	int np;
	Voro3D_GasParticle *bp;
} CellType;

void MkLinkedList(Voro3D_GasParticle *, int );

