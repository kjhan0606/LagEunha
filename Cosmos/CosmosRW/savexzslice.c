#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include "eunha.h"

void savexzslice(SimParameters *simpar, DenType *den){
	MPI_Status status;
	float *tden;
	FILE *fp;
	char filename[190];
	ptrdiff_t mx,offset, i,j,k,ii, slocal_nz;
	ptrdiff_t nx,ny,nz;
	ptrdiff_t local_nz = LOCAL_NZ(simpar);

	nx = NX(simpar);
	ny = NY(simpar);
	nz = NZ(simpar);
	mx = 2*(nx/2+1);
	if(local_nz) tden = (DenType*)malloc(sizeof(DenType)*nx*local_nz);
	for(i=0;i<nx*local_nz;i++) tden[i] = 0;
	for(k=0;k<local_nz;k++) for(j=0;j<3;j++) for(i=0;i<nx;i++){
		tden[i+nx*k] += den[i+mx*(j+ny*k)];
	}
	if(MYID(simpar)==0){
		DenType *allden = (DenType*)malloc(sizeof(DenType)*nx*nz);
		for(i=0;i<nx*local_nz;i++) allden[i] = tden[i];
		offset = local_nz*nx;
		for(i=1;i<NID(simpar);i++){
			MPI_Recv(&slocal_nz,1,MPI_LONG,i,0,MPI_COMM(simpar), &status);
			MPI_Recv(allden+offset, slocal_nz*nx, MPI_DEN_T,i,1,MPI_COMM(simpar),&status);
			offset += slocal_nz *nx;
		}
		sprintf(filename,"xzslice.%.5d", STEPCOUNT(simpar));
		fp = fopen(filename,"w");
		int nnx,nnz;
		nnx = nx;
		nnz = nz;
		fwrite(&nnx,sizeof(int), 1,fp);
		fwrite(&nnz,sizeof(int), 1,fp);
		fwrite(allden,sizeof(DenType), nx*nz,fp);
		fclose(fp);
		free(allden);
	}
	else {
		MPI_Send(&local_nz,1, MPI_LONG, 0,0,MPI_COMM(simpar));
		MPI_Send(tden,local_nz*nx, MPI_DEN_T, 0,1,MPI_COMM(simpar));
	}
	if(local_nz) free(tden);
}
