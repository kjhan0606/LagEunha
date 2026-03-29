#define MAX_MALLOC 5000000
#define PPTR(A) ((void **)(&A))
/* for particle num > 2000000000 */
#define MPI_INT8	MPI_LONG_LONG
typedef  long long INT8;
typedef struct memorystrcut {
	INT8 Size;
	void *Starting;
	void **PtrToVariable;
} memorystruct;
INT8 Make_Total_Memory(INT8);
INT8 CheckAvailableMemory();
void *Calloc(INT8,INT8,void **);
void *Malloc(INT8, void **);
void *Realloc(void *, INT8,void **);
INT8 freespace();
void freelast(void *);
INT8 ptrsize(void *);
void *resizelast(void *,INT8);
void NumMemStack();
void FreeRightNumMemStack();
void LastSwitchPointer(void **);
void MemSwitchPointer(void **,void **);
INT8 CurMemStack();
void InitialOldMemStack(INT8);
void StackPosition(void *a);
void Free(void *a);
#define MEGABYTE 1048576L
#ifndef NMEG
#define NMEG 100L
#endif
#define MYINFINITY -1L

#ifndef MEMMAIN
#define Realloc(A,B) Realloc(A,B,PPTR(A))
#define MemSwitchPointer(A,B) MemSwitchPointer(PPTR(A),PPTR(B))
#endif

#ifdef MEMMAIN
#define EXTERN static
#else
#define EXTERN extern
#endif
EXTERN memorystruct Memory[MAX_MALLOC];
EXTERN void *FREE;
EXTERN INT8 tsize;
EXTERN INT8 Current_Stack,Current_Stack_org;
EXTERN void *START_TOTAL_MEMORY;
EXTERN void *LastEnd;
#undef EXTERN
