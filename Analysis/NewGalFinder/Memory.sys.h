#define MAX_MALLOC 5000000
#define MPI_INT8	MPI_LONG_LONG
typedef  long long INT8;
#define Malloc(a,b) malloc(a)

#define Free(a) free(a)

#define Realloc(a,b) realloc(a,b)
#define Make_Total_Memory(a) (a)
#define Calloc(a,b) calloc(a,b)

#define MEGABYTE 1048576L
#ifndef NMEG
#define NMEG 100L
#endif
#define MYINFINITY -1L

