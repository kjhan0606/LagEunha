#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#define MEMMAIN
#include "Memory.h"
#undef MEMMAIN
INT8 CheckAvailableMemory(){
	INT8 i;
	INT8 size;
	size = 0;
	for(i=0;i<Current_Stack;i++) size += Memory[i].Size;
	size = tsize - size;
	return size;
}
INT8 SSSSMake_Total_Memory(){
	INT8 i,size;
	Current_Stack = 0;
	FREE = NULL;
	size = NMEG * MEGABYTE;
	for(i=0;i<5;i++){
			FREE = (void*)malloc(i*300L*MEGABYTE);
			if(FREE) printf("Suceed in %lld\n",i*600L*MEGABYTE);
			free(FREE);
	}
	return 1;
}
INT8 Make_Total_Memory(INT8 size){
	INT8 i;
	Current_Stack = 0;
	FREE = NULL;
	for(i=0;i<MAX_MALLOC;i++){
		Memory[i].Starting = NULL;
		Memory[i].Size = 0;
		Memory[i].PtrToVariable = NULL;
	}
//	size = NMEG * MEGABYTE;
	/*
	printf("Trying to allocating %ld memory\n",size);
	*/
	while((FREE = (void *) malloc(size))==NULL){
		fprintf(stderr,"Error: Not enough memory space\n");
		system("hostname");
		return 0L;
		exit(999);
		size -= 10*MEGABYTE;
		fprintf(stderr,"Warning: Being able to initialize only %lld Mbytes memory space\n",size);
		if(size<0) {
			fprintf(stderr,"Error: No available memory space\n");
			exit(999);
		}
	}
	/*
	{
			INT8 tempsize;
			tempsize = size/5L;
			for(i=1;i<=5;i++){
				if(FREE) free(FREE);
				size = tempsize*i;
				FREE = (void*)malloc(size);
				printf("total %ld size is allocated\n",size);
			}

	}
	*/
	/*
	printf("Allocated %ld memory %p\n",size,FREE);
	*/
	LastEnd = (void*)((char *)FREE + size);
	tsize = size;
	START_TOTAL_MEMORY = FREE;
	if(FREE == NULL){
		size /= MEGABYTE;
		fprintf(stderr,"Error when initializing %lld Mbytes memory space\n",size);
		return 0;
	}
	else {
		size /= MEGABYTE;
		fprintf(stderr,"Initializing %lld Mbytes memory space\n",size);
		return 1;
	}
}
INT8 ptrsize(void *p){
	INT8 i;
	for(i=0;i<Current_Stack;i++){
		if(Memory[i].Starting == p) return Memory[i].Size;
	}
	return 0L;
}
INT8 freespace(){
	INT8 free;
	if(Current_Stack == 0){
		return tsize;
	}
	else{
		free= tsize - ((INT8)((char *)Memory[Current_Stack-1].Starting
				-(char *)START_TOTAL_MEMORY) + Memory[Current_Stack-1].Size);
		/*
		printf("Total size =%ld and current_free_stack=%ld\n",tsize,free);
		*/
		return free;
	}
}
void *Malloc(INT8 size, void **src){
	void *value;
	char *a;
	if(size == 0) {
		printf("Warning : size %lld is being allocated\n",size);
		return NULL;
	}
	if(size == MYINFINITY){
		if(Current_Stack == 0){
			size = tsize;
		}
		else {
			size = tsize - ( (INT8)((char *)Memory[Current_Stack-1].Starting-
					(char *)START_TOTAL_MEMORY) + Memory[Current_Stack-1].Size);
		}
		printf("change allocated memory size from infinity to %lld\n",size);
	}
	Memory[Current_Stack].Size = size;
	Memory[Current_Stack].Starting = FREE;
	Memory[Current_Stack].PtrToVariable = src;
	value = FREE;
	FREE = (void *)((char *)FREE + size);
/*
	{
		a = (char *)malloc(size/20);
		for(INT8 i=0;i<size/20;i++) a[i] = (i%111)^0xa5;
	}
*/
	if(FREE > (void *)((char *)START_TOTAL_MEMORY + tsize)) {
		INT8 tmptsize;
		FREE = (void *)((char *)FREE-size); /* restore the original mem */
		fprintf(stderr,"Error allocating memory %lld bytes in Malloc()\n",size);
		fprintf(stderr,"need more %d bytes \n",(int)(size - 
					(tsize-((char *)FREE-(char *)START_TOTAL_MEMORY))));
                tmptsize = NMEG;
		fprintf(stderr,"total memory is %lld Mbytes\n",tmptsize);
		value = NULL;
		value = (void *) malloc(size);

/*
		for(INT8 i=0;i<size/20;i++) *((char*)value+i) = a[i];
*/

		if(value == NULL) {
			fprintf(stderr,"Cannot allocate memory any more\n");
			exit(0);
		}
		else {
			return value;
		}
	}
/*
	{
		free(a);
	} 
*/
	Current_Stack++;
	return value;
}
void *Calloc(INT8 num, INT8 size,void **src){
	INT8 i;
	char *value;
	size = size *num;
	value = (char *)Malloc(size,src);
	for(i=0;i<size;i++){
		value[i] = 0;
	}
	return (void *) value;
}
void EraseMemory(INT8 srcn){
	INT8 size=0,esize;
	INT8 i,j,k;
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
	INT8 i;
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
void *Realloc(void *a, INT8 size,void **src){
	INT8 i,j,diffsize,nmove;
	void *from,*to;
	void *ptr,*ptr2;
	void *value;
	char *tmpptr;
	/* freeing the memory if size == 0 */
	if(size <0) {
		fprintf(stderr,"Error,,, minus size reallocation \n");
		fflush(stderr);
		exit(9999);
	}
	if(size == 0) {
		Free(a);
		return NULL;
	}
	if(a == NULL){
		printf("Strange in Realloc ");
		printf("%p  size of %lld\n",src,size);fflush(stdout);
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
	if(tsize <  (INT8)((char *)SMalloc[i]-
		(char *)START_TOTAL_MEMORY + size)) {
	*/
	if(tsize < (INT8)((char *)FREE-(char *)START_TOTAL_MEMORY + diffsize)) {
		fprintf(stderr,"Error allocating memory %lld bytes in Realloc()\n",size);
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
void *resizelast(void *p, INT8 size){
	p = (void *)Realloc(p,size,Memory[Current_Stack-1].PtrToVariable);
	return p;
}
void *ReallocLast(void *p,INT8 size){
	void *result;
	INT8 diffsize;
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
	INT8 i;
	for(i=Current_Stack_org;i<Current_Stack;i++){
		*(Memory[i].PtrToVariable) = NULL;
	}
	FREE = (void *)((char *) Memory[Current_Stack_org-1].Starting
			+Memory[Current_Stack_org-1].Size);
	Current_Stack = Current_Stack_org;
	return;
}
INT8 CurMemStack(){
	return Current_Stack;
}
void InitialOldMemStack(INT8 oldstacknum){
	INT8 i;
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
	INT8 i;
	for(i=0;i<Current_Stack;i++){
		if(Memory[i].PtrToVariable == den){
			Memory[i].PtrToVariable = den2;
			return;
		}
	}
}
void StackPosition(void *a){
	INT8 i;
	for(i=0;i<Current_Stack;i++){
		if((char *)a == Memory[i].Starting){
			printf("Now at stack # %lld in %lld\n",i,Current_Stack);fflush(stdout);
		}
	}
}
