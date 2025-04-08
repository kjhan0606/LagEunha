#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
//#include "mpirms.h"
#include "denrw.h"

#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))

#define MAX(a,b) ( (a)>(b)?(a):(b))
#define MIN(a,b) ( (a)<(b)?(a):(b))


void voldenout( SimParameters *simpar, GridInfo *dengrid, GridInfo outgrid, char *outfile){
	ptrdiff_t i,j,k;
	int myid = MYID(simpar);
	int nid = NID(simpar);
	FILE *wp;
	int src = myid-1;
	int tgt = myid+1;
	int iget, WGroupSize, isend,itag=1;
	MPI_Status status;
	int headercount = 0;




	int  nnx,nny,nnz;
	nnx = outgrid.jx-outgrid.ix+1;
	nny = outgrid.jy-outgrid.iy+1;
	nnz = outgrid.jz-outgrid.iz+1;

#ifdef DEBUG
	if(1){
		DenType *inden = (DenType*)(dengrid+1);
		MPI_Barrier(MPI_COMM(simpar));
		DEBUGPRINT("*P%d has %d %d %d :: %g \n", MYID(simpar), nnx,nny,nnz,inden[nnx*nny*nnz-1]);
	}
#endif

	if(myid == 0){
		long npixels = 1;
		wp = fopen(outfile,"w");
		if(nnx>1){
			fwrite(&nnx,sizeof(int), 1, wp);
			headercount ++;
			npixels *= nnx;
		} 
		if(nny>1){
			fwrite(&nny,sizeof(int), 1, wp);
			headercount ++;
			npixels *= nny;
		}
		if(nnz>1) {
			fwrite(&nnz,sizeof(int), 1, wp);
			headercount ++;
			npixels *= nnz;
		}
		fseek(wp, sizeof(DenType)*(npixels-1L), SEEK_CUR);
		DenType null = 100;
		fwrite(&null, sizeof(DenType), 1, wp);
		fclose(wp);
	}
	MPI_Bcast(&headercount, 1, MPI_INT, 0, MPI_COMM(simpar));
	MPI_Barrier(MPI_COMM(simpar));

	ptrdiff_t ix,iy,iz,jx,jy,jz;
	ix = MAX(outgrid.ix,dengrid->ix);
	iy = MAX(outgrid.iy,dengrid->iy);
	iz = MAX(outgrid.iz,dengrid->iz);
	jx = MIN(outgrid.jx,dengrid->jx);
	jy = MIN(outgrid.jy,dengrid->jy);
	jz = MIN(outgrid.jz,dengrid->jz);
	DenType *outvolume;
	ptrdiff_t mx = (jx-ix + 1);
	ptrdiff_t my = (jy-iy + 1);
	ptrdiff_t mz = (jz-iz + 1);
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));
	DEBUGPRINT("-P%d has %ld %ld %ld ::: %ld %ld %ld : %ld %ld %ld :: %d %d %d\n", MYID(simpar), ix,iy,iz,jx,jy,jz, mx,my,mz,nnx,nny,nnz);
#endif
	if(mx>0 && my >0 && mz >0){
		int nx = dengrid->jx-dengrid->ix + 1;
		int ny = dengrid->jy-dengrid->iy + 1;
		int nz = dengrid->jz-dengrid->iz + 1;
		outvolume = (DenType*)malloc(sizeof(DenType)*mx*my*mz);
		DenType *inden = (DenType*)(dengrid+1);
		for(k=iz;k<=jz;k++) for(j=iy;j<=jy;j++) for(i=ix;i<=jx;i++){
			outvolume[i-ix + mx*(j - iy + my*(k-iz))] = inden[i-dengrid->ix + nx*(j-dengrid->iy+ny*(k-dengrid->iz))];
		}
	}
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));
	DEBUGPRINT("P%d has %ld %ld %ld ::: %ld %ld %ld : %ld %ld %ld :: %d %d %d\n", MYID(simpar), ix,iy,iz,jx,jy,jz, mx,my,mz,nnx,nny,nnz);
#endif

	WGroupSize = WGROUPSIZE(simpar);

	if(RANKINGROUP(myid, WGroupSize) !=0) 
		MPI_Recv(&iget, 1, MPI_INT, src, itag, MPI_COMM(simpar), &status);

	if(mx>0 && my >0 && mz>0){
		wp = fopen(outfile,"r+");
		for(k=iz;k<=jz;k++) for(j=iy;j<=jy;j++){
			long offset = sizeof(int)*headercount + sizeof(DenType)*(ix - outgrid.ix + nnx*(j-outgrid.iy+nny*(k-outgrid.iz)));
			fseek(wp, offset, SEEK_SET);
			fwrite(outvolume + mx*(j-iy + my*(k-iz)), sizeof(DenType), mx, wp);
			fflush(wp);
		}
		fclose(wp);
		free(outvolume);

	}
	if(GROUPID(myid, WGroupSize) == GROUPID(tgt, WGroupSize) && tgt <nid) 
		MPI_Send(&isend, 1, MPI_INT, tgt, itag, MPI_COMM(simpar));
}



void volumeout(SimParameters *simpar, GridInfo *ingrid, enum xyzslicetype xyz){ 
	GridInfo outgrid; 
	char outfile[190];
	sprintf(outfile,"%s.%.5d",outdenname[xyz],STEPCOUNT(simpar));
	if(xyz ==xzslice){ 
		outgrid.ix = outgrid.iy = outgrid.iz = 0; 
		outgrid.jx = NX(simpar)-1; outgrid.jy = 0; outgrid.jz = NZ(simpar)-1; 
	}   
	if(xyz ==xyslice){ 
		outgrid.ix = outgrid.iy = outgrid.iz = 0; 
		outgrid.jx = NX(simpar)-1; outgrid.jy = NY(simpar)-1; outgrid.jz = 0;
	}   
	if(xyz ==yzslice){ 
		outgrid.ix = outgrid.iy = outgrid.iz = 0; 
		outgrid.jx = 0; outgrid.jy = NY(simpar)-1; outgrid.jz = NZ(simpar)-1;
	}   
	if(xyz ==volume){ 
		outgrid.ix = outgrid.iy = outgrid.iz = 0; 
		outgrid.jx = NX(simpar)-1; outgrid.jy = NY(simpar)-1; outgrid.jz = NZ(simpar)-1;
	}   
	voldenout(simpar, ingrid, outgrid, outfile);
}

