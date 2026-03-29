#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<mpi.h>
#include<math.h>

#define MIN(a,b) ( (a)<(b)? (a):(b))
#define MAX(a,b) ( (a)>(b)? (a):(b))

#define mpi_nchunk 500000000L

void BIG_MPI_Send(const void *a, size_t ncount, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
	size_t dsize;
	if(datatype == MPI_CHAR){
		dsize = sizeof(char);
	}
	else if(datatype == MPI_BYTE){
		dsize = sizeof(char);
	}
	else if(datatype == MPI_INT){
		dsize = sizeof(int);
	}
	else if(datatype == MPI_FLOAT){
		dsize = sizeof(float);
	}
	else if(datatype == MPI_LONG){
		dsize = sizeof(long);
	}
	else if(datatype == MPI_DOUBLE){
		dsize = sizeof(double);
	}
	else if(datatype == MPI_LONG_LONG){
		dsize = sizeof(long long);
	}
	else if(datatype == MPI_LONG_DOUBLE){
		dsize = sizeof(long double);
	}
	else {
		fprintf(stderr,"Unidentified types of data in mpi_send\n");
	}
	size_t worktodo = ncount;
	size_t nchunk;
	void *b = (void*)a;
	while(worktodo >0){
		nchunk = MIN(mpi_nchunk, worktodo);
		MPI_Send(b, nchunk, datatype, dest, tag, comm);
		worktodo = worktodo - nchunk;
		b = (void*)((char*)b + nchunk*dsize);
	}
}

void BIG_MPI_Recv(void *a, size_t ncount, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
		MPI_Status *status){
	size_t dsize;
	if(datatype == MPI_CHAR){
		dsize = sizeof(char);
	}
	else if(datatype == MPI_BYTE){
		dsize = sizeof(char);
	}
	else if(datatype == MPI_INT){
		dsize = sizeof(int);
	}
	else if(datatype == MPI_FLOAT){
		dsize = sizeof(float);
	}
	else if(datatype == MPI_LONG){
		dsize = sizeof(long);
	}
	else if(datatype == MPI_DOUBLE){
		dsize = sizeof(double);
	}
	else if(datatype == MPI_LONG_LONG){
		dsize = sizeof(long long);
	}
	else if(datatype == MPI_LONG_DOUBLE){
		dsize = sizeof(long double);
	}
	else {
		fprintf(stderr,"Unidentified types of data in mpi_send\n");
	}
	size_t worktodo = ncount;
	size_t nchunk;
	void *b = (void*)a;
	while(worktodo >0){
		nchunk = MIN(mpi_nchunk, worktodo);
		MPI_Recv(b, nchunk, datatype, dest, tag, comm, status);
		worktodo = worktodo - nchunk;
		b = (void*)((char*)b + nchunk*dsize);
	}
}
