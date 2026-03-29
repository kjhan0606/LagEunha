typedef long indxtype ;

typedef struct HaloInfo{
	int nsub,ndm,nstar,nsink,ngas,npall;
	dptype totm,mdm,mgas,msink,mstar;
	dptype x,y,z,vx,vy,vz;
}HaloInfo;

typedef struct Pos{
	dptype x,y,z;
}Pos;

/*
typedef struct SubInfo{
	int npdm,npgas,npsink,npstar,npall;
	Pos starpeak, dmpeak;
	dptype totm,mdm,mgas,msink,mstar;
	dptype x,y,z,vx,vy,vz;
}SubInfo;
*/
typedef struct SubInfo{
	int npdm,npgas,npsink,npstar,npall;
	dptype totm,mdm,mgas,msink,mstar;
	dptype x,y,z,vx,vy,vz;
}SubInfo;



#define NOT_HALO_MEMBER -9999
#define MPI_INDX MPI_LONG_LONG
#define nblur 3
#define pixmin D_pixmin
#define xhalo D_xhalo
#define yhalo D_yhalo
#define nendhalo D_nendhalo
#define numhalo D_numhalo
#define amax D_amax
#define a D_a
#define rng D_rng
#define size D_size
#define hubble D_hubble
#define ng D_ng
#define nspace D_nspace
#define omep D_omep
#define omeplam D_omeplam
#define epsilon D_epsilon
#define nx D_nx
#define ny D_ny
#define nz D_nz
#define onesolarmass D_onesolarmass
#define pntmass D_pntmass
#define r1kineticfact D_r1kineticfact
#define r2kineticfact D_r2kineticfact
#define com2real D_com2real
#define real2com D_real2com
#define potentfact D_potentfact
#define mpeak D_mpeak
#define ptl2halonum D_ptl2halonum
#define halo D_halo
#define yhalo D_yhalo
#define nendhalo D_nendhalo
#define numhalo D_numhalo
#define ulinked_halo D_ulinked_halo
#define dlinked_halo D_dlinked_halo
#define numulinked_halo D_numulinked_halo
#define myrank D_myrank
#define nrank D_nrank
#define m_tidal D_m_tidal
#define r_tidal D_r_tidal
#define den D_den
#define den2ptl D_den2ptl
#define magnitude D_magnitude
#define acoeff D_acoeff
/*
find_halo.c:extern float pixmin;
find_halo.c:extern int *xhalo[];
find_halo.c:extern int *yhalo[];
find_halo.c:extern int *nendhalo[];
find_halo.c:extern int numhalo[];
hdf_disk_IO.c:extern float amax,a,rng,size,hubble;
hdf_disk_IO.c:extern int ng,nspace;
hdf_disk_IO.c:extern float omep,omeplam,epsilon;
hdf_disk_IO.c:extern int nx,ny,nz;
hdf_disk_IO.c:extern double onesolarmass,pntmass,r1kineticfact,r2kineticfact;
hdf_disk_IO.c:extern double com2real,real2com,potentfact;
hdf_disk_IO.c:extern lint mpeak;
hdf_disk_IO.c:extern lint *ptl2halonum;
hfind.c:extern int physical_parameters(void);
link_slice_halo.c:extern int *halo[];
link_slice_halo.c:extern int *yhalo[];
link_slice_halo.c:extern int *nendhalo[];
link_slice_halo.c:extern int numhalo[];
link_slice_halo.c:extern int *ulinked_halo[];
link_slice_halo.c:extern int *dlinked_halo[];
link_slice_halo.c:extern int numulinked_halo[];
link_slice_halo.c:extern int numdlinked_halo[];
subhalo.den.c:extern int nx,ny,nz;
subhalo.den.c:extern int myrank,nrank;
subhalo.den.c:extern double onesolarmass;
subhalo.den.c:extern double com2real,real2com,potentfact;
subhalo.den.c:extern double pntmass;
subhalo.den.c:extern double r1kineticfact,r2kineticfact;
subhalo.den.c:extern float amax,a,rng,size,hubble;
subhalo.den.c:extern int ng,nspace;
subhalo.den.c:extern float omep,omeplam,acoeff[8];
subhalo.den.c:extern float epsilon;
subhalo.den.c:extern float m_tidal[NUM_MASS],r_tidal[NUM_MASS];
subhalo.den.c://extern lint *ptl2halonum;
subhalo.den.c:extern float *den;
subhalo.den.c:extern lint *den2ptl;
subhalo.den.c:extern int magnitude;
test.c:extern int physical_parameters(void);
*/
