#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
//#include "mpirks.h"
#include "parallelIO.h"
#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))

void pwrite(void *a, size_t lsize, size_t nmem, long offset, char *filename, SimParameters *simpar){
	int myid, nid;
	myid = MYID(simpar);
	nid = NID(simpar);
	int WGroupSize = WGROUPSIZE(simpar);
	int src = myid-1;
	int tgt = myid+1;
	int iget, isend,itag=1;
	MPI_Status status;

	if(RANKINGROUP(myid, WGroupSize) !=0) 
		MPI_Recv(&iget, 1, MPI_INT, src, itag, MPI_COMM(simpar), &status);

	{
		FILE *wp = fopen(filename,"r+");
		fseek(wp, offset, SEEK_SET);
		fwrite(a, lsize, nmem,  wp);
		fclose(wp);
	}
	if(GROUPID(myid, WGroupSize) == GROUPID(tgt, WGroupSize) && tgt <nid) 
		MPI_Send(&isend, 1, MPI_INT, tgt, itag, MPI_COMM(simpar));
}
void pread(void *a, size_t lsize, size_t nmem, long offset, char *filename, SimParameters *simpar){
	int myid, nid;
	myid = MYID(simpar);
	nid = NID(simpar);
	int WGroupSize = WGROUPSIZE(simpar);
	int src = myid-1;
	int tgt = myid+1;
	int iget,  isend,itag=1;
	MPI_Status status;

	if(RANKINGROUP(myid, WGroupSize) !=0) 
		MPI_Recv(&iget, 1, MPI_INT, src, itag, MPI_COMM(simpar), &status);

	{
		FILE *fp = fopen(filename,"r");
		fseek(fp, offset, SEEK_SET);
		fread(a, lsize, nmem,  fp);
		fclose(fp);
	}
	if(GROUPID(myid, WGroupSize) == GROUPID(tgt, WGroupSize) && tgt <nid) 
		MPI_Send(&isend, 1, MPI_INT, tgt, itag, MPI_COMM(simpar));
}

void pGraficRead(GridInfo *a,  SimParameters *simpar, int mtype, int ivec){
	int myid, nid;
	myid = MYID(simpar);
	nid = NID(simpar);
	int WGroupSize = WGROUPSIZE(simpar);
	int src = myid-1;
	int tgt = myid+1;
	int iget,  isend,itag=1;
	MPI_Status status;
	DenType *den = (DenType*)(a+1);

	if(RANKINGROUP(myid, WGroupSize) !=0) 
		MPI_Recv(&iget, 1, MPI_INT, src, itag, MPI_COMM(simpar), &status);

	{
		char filename[190];
		char it,iv;
		if(mtype==0) it = 't';
		else if(mtype==1) it = 'c';
		else if(mtype==2) it = 'b';
		if(ivec%3==0) iv = 'x';
		else if(ivec%3==1) iv = 'y';
		else if(ivec%3==2) iv = 'z';
		if(ivec <3) sprintf(filename,"%s./ic_pos%c%c",GRAFIC_DIRECTORY(simpar), it,iv);
		else sprintf(filename,"%s./ic_vel%c%c",GRAFIC_DIRECTORY(simpar), it,iv);
		printf("P%d is trying to read %s\n", myid, filename);
		FILE *fp = fopen(filename,"r");
		long offset;
		int chip;
		size_t i,nx,ny,nz;
		nx = NX(simpar);
		ny = NY(simpar);
		nz = NZ(simpar);
		fread(&chip, sizeof(int),1,fp);
		{
			int mx,my,mz;
			float dx,lx,ly,lz,astart, omegam, omegav, h0;
			fread(&mx,sizeof(int), 1, fp);
			fread(&my,sizeof(int), 1, fp);
			fread(&mz,sizeof(int), 1, fp);
			fread(&dx,sizeof(float), 1, fp);
			fread(&lx,sizeof(float), 1, fp);
			fread(&ly,sizeof(float), 1, fp);
			fread(&lz,sizeof(float), 1, fp);
			fread(&astart,sizeof(float), 1, fp);
			fread(&omegam,sizeof(float), 1, fp);
			fread(&omegav,sizeof(float), 1, fp);
			fread(&h0,sizeof(float), 1, fp);
			DEBUGPRINT("P%d mx/my/mz ... %d %d %d %g %g %g %g %g %g %g %g : with chip= %d\n",
					myid, mx,my,mz,dx,lx,ly,lz,astart, omegam,omegav,h0, chip);
			if(mx!=nx || my!=ny|| mz != nz || omegam != OMEP(simpar) || omegav != OMEPLAM(simpar)){
				fprintf(stderr,"ERROR: mismatch between grafic file and the input simulation parameters\n");
				exit(999);
			}
		}
		fread(&chip, sizeof(int), 1, fp);

		offset = 2*sizeof(int) + chip;
		offset += (2*sizeof(int) + sizeof(DenType)*nx*ny)*LOCAL_Z_START(simpar);
		fseek(fp, offset, SEEK_SET);
		for(i=0;i<LOCAL_NZ(simpar);i++){
			fread(&chip, sizeof(int), 1,fp);
			if(chip != sizeof(DenType)*nx*ny)DEBUGPRINT("P%d has error in reading data %d %ld\n", myid, chip, nx*ny*sizeof(DenType));
			size_t iread = fread(den, sizeof(DenType), nx*ny,fp);
			if(iread != nx*ny)DEBUGPRINT("P%d has error in reading data %ld %ld\n", myid, nx*ny, iread);
			fread(&chip, sizeof(int), 1,fp);
			DEBUGPRINT("P%d has value %g %g\n", MYID(simpar), den[0], den[nx*ny-1]);
			den += nx*ny;
		}
		fclose(fp);
	}
	if(GROUPID(myid, WGroupSize) == GROUPID(tgt, WGroupSize) && tgt <nid) 
		MPI_Send(&isend, 1, MPI_INT, tgt, itag, MPI_COMM(simpar));
}
