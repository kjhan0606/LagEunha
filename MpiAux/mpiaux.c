#include <mpi.h>
#include <limits.h>
#include "eunha.h"
#include "mpiaux.h"

void Mpi_Basic_Set(SimParameters *simpar, MPI_Comm com){
	MPICOM(simpar) = MPI_COMM(simpar) = com;
	int myid,nid;
	MPI_Comm_rank(com,&myid);
	MPI_Comm_size(com,&nid);
	MYID(simpar) = myid;
	NID(simpar) = nid;
	MPI_COMM(simpar) = com;
	WGROUPSIZE(simpar) = (nid + NWGROUP(simpar) - 1) /NWGROUP(simpar);
}

void my_MPI_Sendrecv(void *sbase, long smem, MPI_Datatype datatype1, int dest, int stag,
	void *rbase, long rmem, MPI_Datatype datatype2, int src, int rtag, MPI_Comm comm, MPI_Status *status){
	long nmax, nchunk = INT_MAX/4, max;
	long i;
	MPI_Reduce(&smem, &nmax, 1,MPI_LONG, MPI_MAX, 0, comm);
	max = nmax;
	MPI_Reduce(&rmem, &nmax, 1,MPI_LONG, MPI_MAX, 0, comm);
	if(nmax > max) max = nmax;
	MPI_Bcast(&max, 1, MPI_LONG,0,comm);
	for(i=0;i<max;i+=nchunk){
		int ismem, irmem;
		long lsmem, lrmem;
		lsmem = smem -i;
		lrmem = rmem -i;
		if(lsmem > nchunk) ismem = nchunk;
		else if(lsmem <=0 ) ismem = 0;
		else ismem = lsmem;
		if(lrmem > nchunk) irmem = nchunk;
		else if(lrmem <=0 ) irmem = 0;
		else irmem = lrmem;
		MPI_Sendrecv((char*)sbase+i, ismem, datatype1,dest,stag, 
			(char*)rbase+i,irmem, datatype2, src, rtag, comm, status);
	}
}
void my_MPI_Send(void *sbase, long smem, MPI_Datatype datatype, int dest, int stag, MPI_Comm comm){
	long nchunk = INT_MAX/4;
	long i;
	for(i=0;i<smem;i+=nchunk){
		int ismem;
		long lsmem;
		lsmem = smem - i;
		if(lsmem > nchunk) ismem = nchunk;
		else if (lsmem <=0) ismem = 0;
		else ismem = lsmem;
		MPI_Send((char*) sbase+i, ismem, datatype, dest, stag, comm);
	}
}
void my_MPI_Recv(void *rbase, long smem, MPI_Datatype datatype, int src, int rtag, MPI_Comm comm, MPI_Status *status){
	long nchunk = INT_MAX/4;
	long i;
	for(i=0;i<smem;i+=nchunk){
		int ismem;
		long lsmem;
		lsmem = smem - i;
		if(lsmem > nchunk) ismem = nchunk;
		else if (lsmem <=0) ismem = 0;
		else ismem = lsmem;
		MPI_Recv((char*) rbase+i, ismem, datatype, src, rtag, comm, status);
	}
}

#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))
void StartParallelRW(int myid, int WGroupSize, MPI_Comm com){ 
	int iget; 
	int src,tgt,itag=129; 
	MPI_Status status; 
	iget = 1; 
	src = myid-1; 
	tgt = myid+1; 
	if(WGroupSize==0) WGroupSize = 1;
	if(RANKINGROUP(myid,WGroupSize) != 0 ) 
		MPI_Recv(&iget,1,MPI_INT,src,itag,com,&status);
}
void CloseParallelRW(int myid, int nid, int WGroupSize, MPI_Comm com){ 
	int isend,iget; 
	int src,tgt,itag=129; 
	isend = iget = 1; 
	src = myid-1; 
	tgt = myid+1; 
	if(WGroupSize==0) WGroupSize = 1;
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) 
		MPI_Send(&isend,1,MPI_INT,tgt,itag,com);
}
#undef GROUPID
#undef RANKINGROUP
