#define MAX_MALLOC 1000
#define PPTR(A) ((void **)(&A))
typedef struct memorystrcut {
	size_t Size;
	void *Starting;
	void **PtrToVariable;
} memorystruct;
size_t Make_Total_Memory();
size_t CheckAvailableMemory();
void *CALLOC(size_t,size_t,void **);
void *MALLOC(size_t, void **);
void *REALLOC(void *, size_t,void **);
size_t freespace();
void freelast(void *);
size_t ptrsize(void *);
void *resizelast(void *,size_t);
void NumMemStack();
void FreeRightNumMemStack();
void LastSwitchPointer(void **);
void MemSwitchPointer(void **,void **);
size_t CurMemStack();
void InitialOldMemStack(size_t);
void StackPosition(void *a);
void Free(void *a);

#define MEGABYTE 1048576L
#ifndef NMEG
#define NMEG 100L
#endif
#define MYINFINITY -1L



#define Malloc(a,b) MALLOC(b, ((void **)(&(a))))
#define Realloc(a,b) REALLOC(a, b,((void **)(&(a))))
#define Calloc(a,b,c) REALLOC(b, c,((void **)(&(a))))
