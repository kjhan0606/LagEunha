#define MAX_MALLOC 10000
#define PPTR(A) ((void **)(&A))
int Make_Total_Memory();
void *Calloc(size_t,size_t,void **);
void *Malloc(size_t, void **);
void *Realloc(void *, size_t,void **);
size_t freespace();
void freelast(void *);
size_t ptrsize(void *);
void *resizelast(void *,size_t);
void NumMemStack();
void FreeRightNumMemStack();
void LastSwitchPointer(void **);
void MemSwitchPointer(void **,void **);
int CurMemStack();
void InitialOldMemStack(int);
#define MEGABYTE 1048576
#ifndef NMEG
#define NMEG 100
#endif
#define INFINITY -1

#ifndef MEMMAIN
#define Realloc(A,B) Realloc(A,B,PPTR(A))
#define MemSwitchPointer(A,B) MemSwitchPointer(PPTR(A),PPTR(B))
#endif
