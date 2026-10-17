#define OneGB 1073741824L
typedef float DenType;
#define MPI_DenType MPI_FLOAT
typedef struct GridInfo{ 
	int ix,iy,iz; 
	int jx,jy,jz; 
	int nx,ny,nz; 
	ptrdiff_t npix; 
	DenType *den;
} GridInfo;

