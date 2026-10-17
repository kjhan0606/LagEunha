#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#define MEMMAIN
#include "Memory.h"
static void *SMalloc[MAX_MALLOC];
static void *FREE;
static size_t tsize,SizeMalloc[MAX_MALLOC];
static int Current_Stack,Current_Stack_org;
static void **P2Array[MAX_MALLOC];
static void *START_TOTAL_MEMORY;
int Make_Total_Memory(){
	size_t i,size;
	Current_Stack = 0;
	FREE = NULL;
	for(i=0;i<MAX_MALLOC;i++){
		SMalloc[i] = NULL;
		SizeMalloc[i] = 0;
		P2Array[i] = NULL;
	}
	size = NMEG * MEGABYTE;
	FREE = (void *) malloc(size);
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
		if(SMalloc[i] == p) return SizeMalloc[i];
	}
}
size_t freespace(){
	if(Current_Stack == 0){
		return tsize;
	}
	else{
		return tsize - ((size_t)((char *)SMalloc[Current_Stack-1]
				-(char *)START_TOTAL_MEMORY) + SizeMalloc[Current_Stack-1]);
	}
}
void *Malloc(size_t size, void **src){
	void *value;
	char *cvalue;
	if(size == 0) return NULL;
	if(size == INFINITY){
		if(Current_Stack == 0){
			size = tsize;
		}
		else {
			size = tsize - ( (size_t)((char *)SMalloc[Current_Stack-1]-
					(char *)START_TOTAL_MEMORY) + SizeMalloc[Current_Stack-1]);
		}
	}
	SizeMalloc[Current_Stack] = size;
	SMalloc[Current_Stack] = FREE;
	P2Array[Current_Stack] = src;
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
	FREE = (void *)((char *)FREE - SizeMalloc[srcn]);
	*P2Array[srcn] = NULL;
	esize = SizeMalloc[srcn];
	if(srcn == Current_Stack-1){
		Current_Stack--;
		return;
	}
	else {
		from=SMalloc[srcn+1];
		for(i=srcn+1;i<Current_Stack;i++){
			SMalloc[i-1] = (void *)((char *)SMalloc[i]-esize);
			SizeMalloc[i-1] = SizeMalloc[i];
			P2Array[i-1] = P2Array[i];
			*P2Array[i-1] = SMalloc[i-1];
			size += SizeMalloc[i];
		}
		Current_Stack--;
		memmove(SMalloc[srcn],from,size);
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
			if((char *)a == SMalloc[i]){
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
		if((char *)a == SMalloc[i]){
			break;
		}
	}
	/* If this memory was allocated using malloc */
	if(i == Current_Stack) {
		return (void *) realloc(a,size);
	}
	diffsize = size-SizeMalloc[i];
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
		ptr2 = SMalloc[i];
		memcpy(value,ptr2,SizeMalloc[i]);
		/*
		Current_Stack--;
		*/
		Free(a);
		return value;
	}
	if(i != Current_Stack-1){
		if(size > SizeMalloc[i]) {
			from = (void *)((char *)FREE-1);
			to = SMalloc[i+1];
			j = 0;
			tmpptr = (char *)SMalloc[Current_Stack-1]+
				SizeMalloc[Current_Stack-1] +diffsize;
			for(ptr=from;ptr>=to;ptr=(void *)((char *)ptr-1)){
				*(tmpptr+j-1) = *((char *)ptr);
				j--;
			}
		}
		else if( size < SizeMalloc[i]) {
			nmove = 0;
			for(j=i+1;j<Current_Stack;j++){
				nmove += SizeMalloc[j];
			}
			memmove((void *)((char *)SMalloc[i]+size),SMalloc[i+1],nmove);
		}
		for(j=i+1;j<Current_Stack;j++){
			*P2Array[j] = (void **) ((char *)*P2Array[j] + diffsize);
			/*
			SMalloc[j] += diffsize;
			*/
			SMalloc[j] = (void *)((char *) SMalloc[j]+diffsize);
		}
	}
	FREE = (void *)((char *)FREE + diffsize);
	SizeMalloc[i] = size;
	return a;
}
void *resizelast(void *p, size_t size){
	p = (void *)Realloc(p,size,P2Array[Current_Stack-1]);
	return;
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
		*P2Array[i] = NULL;
	}
	FREE = (void *)((char *) SMalloc[Current_Stack_org-1]
			+SizeMalloc[Current_Stack_org-1]);
	Current_Stack = Current_Stack_org;
	return;
}
int CurMemStack(){
	return Current_Stack;
}
void InitialOldMemStack(int oldstacknum){
	int i;
	if(Current_Stack == 0) return;
	for(i=oldstacknum;i<Current_Stack;i++){
		*P2Array[i] = NULL;
	}
	FREE = (void *)((char *) SMalloc[oldstacknum-1]+SizeMalloc[oldstacknum-1]);
	Current_Stack = oldstacknum;
	return;
}
void LastSwitchPointer(void **den){
	P2Array[Current_Stack-1] = den;
}
void MemSwitchPointer(void **den, void **den2){
	int i;
	for(i=0;i<Current_Stack;i++){
		if(P2Array[i] == den){
			P2Array[i] = den2;
			return;
		}
	}
}
