#define Nx 32
#define Ny 32
#define Nz 32
#define Nximg 128
#define Nyimg 128
#define Nzimg 128
#define Nxp 128
#define Nyp 128
#define Nzp 128


#define Lx (1.L)
#define Ly (1.L)
#define Lz (1.L)
#define boxsize (Lx/Nx*5.L)

#ifndef _OPENMP
#define omp_get_thread_num() (0)
#define omp_get_num_threads() (1)
#endif

#define eBlast 1.
#define bckPressure 1.e-6

#define Gamma 1.666666L
#define Zeta (0.125L/100.L)

#define Courant 0.1L
/*
#define alphavis 2
#define betavis 4
*/
#define alphavis 2
#define betavis 4

typedef struct CellType{
	int np;
	Voro3D_GasParticle *bp;
} CellType;

void MkLinkedList(Voro3D_GasParticle *, int );
void Voro3D_FindNeighbor(Voro3D_point *,int , int , int , int *, Voro3D_GasParticle *);
void Voro3D_FindCellBP(Voro3D_GasParticle *,int , int , int , int *, Voro3D_GasParticle *);
double vph3D(Voro3D_GasParticle **, int *, ptype );
Voro3D_GasParticle *mkinitial(int *);
int get_npmax();

