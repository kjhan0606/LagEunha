#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#define MEMMAIN
#include "Memory.h"
typedef struct memorystrcut {
	size_t Size;
	void *Starting;
	void **PtrToVariable;
} memorystruct;
static memorystruct Memory[MAX_MALLOC];
static void *FREE;
static size_t tsize;
static int Current_Stack,Current_Stack_org;
static void *START_TOTAL_MEMORY;
static void *LastEnd;
size_t CheckAvailableMemory(){
	int i;
	size_t size;
	size = 0;
	for(i=0;i<Current_Stack;i++) size += Memory[i].Size;
	size = tsize - size;
	return size;
}
int Make_Total_Memory(){
	size_t i,size;
	Current_Stack = 0;
	FREE = NULL;
	for(i=0;i<MAX_MALLOC;i++){
		Memory[i].Starting = NULL;
		Memory[i].Size = 0;
		Memory[i].PtrToVariable = NULL;
	}
	size = NMEG * MEGABYTE;
	FREE = (void *) malloc(size);
	LastEnd = (void*)((char *)FREE + size);
	tsize = size;
	START_TOTAL_MEMORY = FREE;
	if(FREE == NULL){
		size /= MEGABYTE;
		fprintf(stderr,"Error when initializing %d Mbytes memory space\n",size);
		return 0;
	}
	else {
		return 1;
	}

}
size_t ptrsize(void *p){
	size_t i;
	for(i=0;i<Current_Stack;i++){
		if(Memory[i].Starting == p) return Memory[i].Size;
	}
}
size_t freespace(){
	if(Current_Stack == 0){
		return tsize;
	}
	else{
		return tsize - ((size_t)((char *)Memory[Current_Stack-1].Starting
				-(char *)START_TOTAL_MEMORY) + Memory[Current_Stack-1].Size);
	}
}
void *Malloc(size_t size, void **src){
	void *value;
	if(size == 0) return NULL;
	if(size == INFINITY){
		if(Current_Stack == 0){
			size = tsize;
		}
		else {
			size = tsize - ( (size_t)((char *)Memory[Current_Stack-1].Starting-
					(char *)START_TOTAL_MEMORY) + Memory[Current_Stack-1].Size);
		}
	}
	Memory[Current_Stack].Size = size;
	Memory[Current_Stack].Starting = FREE;
	Memory[Current_Stack].PtrToVariable = src;
	value = FREE;
	FREE = (void *)((char *)FREE + size);
	if(FREE > (void *)((char *)START_TOTAL_MEMORY + tsize)) {
		FREE = (void *)((char *)FREE-size); /* restore the original mem */
		fprintf(stderr,"Error allocating memory %d bytes in Malloc()\n",size);
		fprintf(stderr,"need more %d bytes \n",(int)(size - 
					(tsize-((char *)FREE-(char *)START_TOTAL_MEMORY))));
		value = NULL;
		value = (void *) malloc(size);
		if(value == NULL) {
			fprintf(stderr,"Cannot allocate memory any more\n");
			exit(0);
		}
		else {
			return value;
		}
	}
	Current_Stack++;
	return value;
}
void *Calloc(size_t num, size_t size,void **src){
	size_t i;
	char *value;
	size = size *num;
	value = (char *)Malloc(size,src);
	for(i=0;i<size;i++){
		value[i] = 0;
	}
	return (void *) value;
}
void EraseMemory(size_t srcn){
	size_t size=0,esize;
	size_t i,j,k;
	void *from;
	FREE = (void *)((char *)FREE - Memory[srcn].Size);
	*(Memory[srcn].PtrToVariable) = NULL;
	esize = Memory[srcn].Size;
	if(srcn == Current_Stack-1){
		Current_Stack--;
		return;
	}
	else {
		from=Memory[srcn+1].Starting;
		for(i=srcn+1;i<Current_Stack;i++){
			Memory[i-1].Starting = (void *)((char *)Memory[i].Starting-esize);
			Memory[i-1].Size = Memory[i].Size;
			Memory[i-1].PtrToVariable = Memory[i].PtrToVariable;
			*(Memory[i-1].PtrToVariable) = Memory[i-1].Starting;
			size += Memory[i].Size;
		}
		Current_Stack--;
		memmove(Memory[srcn].Starting,from,size);
	}
}
void Free(void *a){
	size_t i;
	if(a == NULL){
/*
		fprintf(stderr,"Attempt to free NULL pointer\n");
*/
	}
	else {
		for(i=0;i<Current_Stack;i++){
			if((char *)a == Memory[i].Starting){
				EraseMemory(i);
				return;
			}
		}
		if(i==Current_Stack) {
			free(a);
			a = NULL;
		}
	}
}
void *Realloc(void *a, size_t size,void **src){
	size_t i,j,diffsize,nmove;
	void *from,*to;
	void *ptr,*ptr2;
	void *value;
	char *tmpptr;
	/* freeing the memory if size == 0 */
	if(size == 0) {
		Free(a);
		return NULL;
	}
	if(a == NULL){
		printf("Strange in Realloc ");
		printf("%d  size of %d\n",src,size);fflush(stdout);
		a = (void *)Malloc(size,src);
		return a;
	}
	for(i=0;i<Current_Stack;i++){
		if((char *)a == Memory[i].Starting){
			break;
		}
	}
	/* If this memory was allocated using malloc */
	if(i == Current_Stack) {
		return (void *) realloc(a,size);
	}
	diffsize = size-Memory[i].Size;
	if(diffsize == 0) return a;
	/* if there is no more free space */
	/*
	if(tsize <  (size_t)((char *)SMalloc[i]-
		(char *)START_TOTAL_MEMORY + size)) {
	*/
	if(tsize < (size_t)((char *)FREE-(char *)START_TOTAL_MEMORY + diffsize)) {
		fprintf(stderr,"Error allocating memory %d bytes in Realloc()\n",size);
		fprintf(stderr,"need more %d bytes \n",(int)(size - 
					(tsize-((char *)FREE-(char *)START_TOTAL_MEMORY))));
		fflush(stderr);
		value = NULL;
		value = (void *)malloc(size);
		if(value == NULL) {
			fprintf(stderr,"Cannot reallocate any more memory\n");
			fprintf(stderr,"Finish this program !!!!!!!!!!!!!!!! \n");
			exit(0);
		}
		ptr2 = Memory[i].Starting;
		memcpy(value,ptr2,Memory[i].Size);
		/*
		Current_Stack--;
		*/
		Free(a);
		return value;
	}
	if(i != Current_Stack-1){
		if(size > Memory[i].Size) {
			from = (void *)((char *)FREE-1);
			to = Memory[i+1].Starting;
			j = 0;
			tmpptr = (char *)Memory[Current_Stack-1].Starting+
				Memory[Current_Stack-1].Size +diffsize;
			for(ptr=from;ptr>=to;ptr=(void *)((char *)ptr-1)){
				*(tmpptr+j-1) = *((char *)ptr);
				j--;
			}
		}
		else if( size < Memory[i].Size) {
			nmove = 0;
			for(j=i+1;j<Current_Stack;j++){
				nmove += Memory[j].Size;
			}
			memmove((void *)((char *)Memory[i].Starting+size),
					Memory[i+1].Starting,nmove);
		}
		for(j=i+1;j<Current_Stack;j++){
			*(Memory[j].PtrToVariable) = 
				(void **) ((char *)*(Memory[j].PtrToVariable) + diffsize);
			/*
			SMalloc[j] += diffsize;
			*/
			Memory[j].Starting = (void *)((char *) Memory[j].Starting+diffsize);
		}
	}
	FREE = (void *)((char *)FREE + diffsize);
	Memory[i].Size = size;
	return a;
}
void *resizelast(void *p, size_t size){
	p = (void *)Realloc(p,size,Memory[Current_Stack-1].PtrToVariable);
	return;
}
void *ReallocLast(void *p,size_t size){
	void *result;
	size_t diffsize;
	diffsize = size - Memory[Current_Stack-1].Size;
	FREE = (void *)((char *)FREE + diffsize);
	if((char *)LastEnd <=(char *) FREE) 
		return NULL;
	else {
		Memory[Current_Stack-1].Size = size; 
		return Memory[Current_Stack-1].PtrToVariable;
	}
}
void freelast(void *p){
	if(Current_Stack == 0){
		fprintf(stderr,"No current allocated space\n");
	}
	else {
		/*
		Free(*P2Array[Current_Stack-1]);
		*/
		Free(p);
	}
	return;
}
void dumpptr(){
	return;
}
void NumMemStack(){
	Current_Stack_org = Current_Stack;
	return ;
}
void FreeRightNumMemStack(){
	int i;
	for(i=Current_Stack_org;i<Current_Stack;i++){
		*(Memory[i].PtrToVariable) = NULL;
	}
	FREE = (void *)((char *) Memory[Current_Stack_org-1].Starting
			+Memory[Current_Stack_org-1].Size);
	Current_Stack = Current_Stack_org;
	return;
}
int CurMemStack(){
	return Current_Stack;
}
void InitialOldMemStack(int oldstacknum){
	int i;
	for(i=oldstacknum;i<Current_Stack;i++){
		*(Memory[i].PtrToVariable) = NULL;
	}
	FREE = (void *)((char *) Memory[oldstacknum-1].Starting+
			Memory[oldstacknum-1].Size);
	Current_Stack = oldstacknum;
	return;
}
void LastSwitchPointer(void **den){
	Memory[Current_Stack-1].PtrToVariable = den;
}
void MemSwitchPointer(void **den, void **den2){
	int i;
	for(i=0;i<Current_Stack;i++){
		if(Memory[i].PtrToVariable == den){
			Memory[i].PtrToVariable = den2;
			return;
		}
	}
}
