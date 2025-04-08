/* alpha version of parallel domain deocomposition of the Recursive MultiSection (RMS).
   It now allows only the one-dimensional data.
   Further implementation is needed to allow for arbitrary dimensional data. */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<limits.h>
#include<math.h>
#include<mpi.h>
#include "mpirms.h"
/*
#include "../pmheader.h"
*/
#include "../MpiAux/mpiaux.h"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

#define FindOverLap(a,b,ix,iy,iz,jx,jy,jz) do{\
	ix = MAX(a->ix,b->ix);\
	iy = MAX(a->iy,b->iy);\
	iz = MAX(a->iz,b->iz);\
	jx = MIN(a->jx,b->jx);\
	jy = MIN(a->jy,b->jy);\
	jz = MIN(a->jz,b->jz);\
} while(0)


void BuildGridHierarchyFromCurrentGridTopology(DoDeInfo *ddinfo, int nddinfo, 
		GridInfo *outgrid){
	int i;
	GridInfo localgrid;
	ddinfo[nddinfo-1].sublocalgrid = (localgrid = *outgrid);
#ifdef DEBUG
	printf("%d : P%d has %d %d :: %d %d :: %d %d\n",nddinfo-1,ddinfo[0].myid,
				localgrid.ix,localgrid.jx, localgrid.iy,localgrid.jy, localgrid.iz,localgrid.jz);
#endif
	for(i=nddinfo-2;i>=0;i--){
		GridInfo tlocalgrid[ddinfo[i+1].nid];
		MPI_Comm com = ddinfo[i+1].com;
		MPI_Gather(&localgrid,sizeof(GridInfo),MPI_BYTE,
				tlocalgrid,sizeof(GridInfo), MPI_BYTE, 0, com);
		if(ddinfo[i+1].myid==0) {
			int j;
			for(j=0;j<ddinfo[i+1].nid;j++){
				localgrid.ix = MIN(localgrid.ix,tlocalgrid[j].ix);
				localgrid.iy = MIN(localgrid.iy,tlocalgrid[j].iy);
				localgrid.iz = MIN(localgrid.iz,tlocalgrid[j].iz);
				localgrid.jx = MAX(localgrid.jx,tlocalgrid[j].jx);
				localgrid.jy = MAX(localgrid.jy,tlocalgrid[j].jy);
				localgrid.jz = MAX(localgrid.jz,tlocalgrid[j].jz);
			}
		}
		MPI_Bcast(&localgrid,sizeof(GridInfo), MPI_BYTE, 0, com);
		ddinfo[i].sublocalgrid = localgrid;
#ifdef DEBUG
		printf("%d : P%d has %d %d :: %d %d :: %d %d\n",i,ddinfo[0].myid,
				localgrid.ix,localgrid.jx, localgrid.iy,localgrid.jy, localgrid.iz,localgrid.jz);
#endif
	}
}
int OverLapGrid(const void *a, const void *b){
	GridInfo *grid = (GridInfo*)a;
	GridInfo *tgrid = (GridInfo*)b;
	int ix,iy,iz;
	int jx,jy,jz;
	/*
	ix = MAX(grid->ix, tgrid->ix);
	iy = MAX(grid->iy, tgrid->iy);
	iz = MAX(grid->iz, tgrid->iz);
	jx = MIN(grid->jx, tgrid->jx);
	jy = MIN(grid->jy, tgrid->jy);
	jz = MIN(grid->jz, tgrid->jz);
	*/
	FindOverLap(grid,tgrid,ix,iy,iz,jx,jy,jz);
	if(jx >=ix && jy >= iy && jz >= iz)  return 1;
	else return 0;
}


void *ExtractGrid2Send(GridInfo *grid, GridInfo *targetlocalgrid, GridInfo *write){
	int ix,iy,iz;
	int jx,jy,jz;
	/*
	ix = MAX(grid->ix, targetlocalgrid->ix);
	iy = MAX(grid->iy, targetlocalgrid->iy);
	iz = MAX(grid->iz, targetlocalgrid->iz);
	jx = MIN(grid->jx, targetlocalgrid->jx);
	jy = MIN(grid->jy, targetlocalgrid->jy);
	jz = MIN(grid->jz, targetlocalgrid->jz);
	*/
	FindOverLap(grid,targetlocalgrid,ix,iy,iz,jx,jy,jz);
	if(jx >=ix && jy >= iy && jz >= iz) {
		long mx = (jx-ix+1);
		long my = (jy-iy+1);
		long mz = (jz-iz+1);

		write->npix = (long)mx*(long)my*(long)mz;
		write->ix = ix;
		write->iy = iy;
		write->iz = iz;
		write->jx = jx;
		write->jy = jy;
		write->jz = jz;
		long wx = grid->jx-grid->ix +1;
		long wy = grid->jy-grid->iy +1;
		long wz = grid->jz-grid->iz +1;
		long i,j,k;
		/*
		DenType *den = write->den;
		*/
		DenType *den = (DenType*)(write+1);
		DenType *gden = (DenType*)(grid + 1);
		long ni=0;
		for(k=iz;k<=jz;k++){
			long kk = k - grid->iz;
			for(j=iy;j<=jy;j++){
				long jj = j - grid->iy;
				for(i=ix;i<=jx;i++){
					long ii = i - grid->ix;
					/*
					den[ni++] = grid->den[ii+wx*(jj+wy*kk)];
					*/
					den[ni++] = gden[ii+wx*(jj+wy*kk)];
				}
			}
		}
		return (den +ni);
	}
	else {
		targetlocalgrid->npix = 0;
		return write;
	}
}
long  CountGrid2Send(GridInfo *grid, GridInfo *targetlocalgrid){
	int ix,iy,iz;
	int jx,jy,jz;
	/*
	ix = MAX(grid->ix, targetlocalgrid->ix);
	iy = MAX(grid->iy, targetlocalgrid->iy);
	iz = MAX(grid->iz, targetlocalgrid->iz);
	jx = MIN(grid->jx, targetlocalgrid->jx);
	jy = MIN(grid->jy, targetlocalgrid->jy);
	jz = MIN(grid->jz, targetlocalgrid->jz);
	*/
	FindOverLap(grid,targetlocalgrid,ix,iy,iz,jx,jy,jz);
	if(jx >=ix && jy >= iy && jz >= iz) {
		long mx = (jx-ix+1);
		long my = (jy-iy+1);
		long mz = (jz-iz+1);
		return mx*my*mz;
	}
	else {
		return 0;
	}
}
void *ExtractGrids2Send(void *nowmydata, GridInfo targetlocalgrid, long *grid_data_size){
	int i;
	long npixels=0;
	void *sendbuff;
	int ngrids = ((RGridInfo*)nowmydata)->ngrids;
	char *myrun = (char *)nowmydata + sizeof(int);
	int igrid = 0;
	for(i=0;i<ngrids;i++){
		GridInfo *nowgrid = (GridInfo *)myrun;
		long mpixels = CountGrid2Send(nowgrid, &targetlocalgrid);
		if(mpixels>0) igrid++;
		npixels += mpixels;
		myrun += sizeof(GridInfo) + sizeof(DenType)*nowgrid->npix; 
	}
	long sendsize = sizeof(int) + npixels*sizeof(DenType) + igrid*sizeof(GridInfo);

	sendbuff = (void*)malloc(sendsize);
	((RGridInfo*)sendbuff)->ngrids = igrid;
	myrun = (char *)nowmydata + sizeof(int);
	char *sendrun = (char*)sendbuff + sizeof(int);
	for(i=0;i<ngrids;i++){
		GridInfo *nowgrid = (GridInfo *)myrun;
		/*
		((GridInfo*)sendrun)->den = (DenType*)(sendrun + sizeof(GridInfo));
		*/
		sendrun = (char *)ExtractGrid2Send(nowgrid, &targetlocalgrid, (GridInfo*)sendrun);
		if(targetlocalgrid.npix >0) myrun += sizeof(GridInfo) + sizeof(DenType)*nowgrid->npix;
	}
	*grid_data_size = (char *)sendrun - (char*)sendbuff;
	return sendbuff;
}
void Dump2Outgrid(void *nowmydata, GridInfo *outgrid){
	int ngrids = ((RGridInfo*)nowmydata)->ngrids;
	char *run = (char*)nowmydata + sizeof(int);

	long mx = outgrid->jx - outgrid->ix + 1;
	long my = outgrid->jy - outgrid->iy + 1;
	long mz = outgrid->jz - outgrid->iz + 1;

	long ijk;
	for(ijk=0;ijk<ngrids;ijk++){
		GridInfo *nowgrid = (GridInfo*) run;
		int ix,iy,iz;
		int jx,jy,jz;
		/*
		ix = MAX(nowgrid->ix, outgrid->ix);
		iy = MAX(nowgrid->iy, outgrid->iy);
		iz = MAX(nowgrid->iz, outgrid->iz);
		jx = MIN(nowgrid->jx, outgrid->jx);
		jy = MIN(nowgrid->jy, outgrid->jy);
		jz = MIN(nowgrid->jz, outgrid->jz);
		*/
		FindOverLap(nowgrid,outgrid,ix,iy,iz,jx,jy,jz);
		if(jx >=ix && jy >= iy && jz >=iz) {
			long nx = nowgrid->jx - nowgrid->ix + 1;
			long ny = nowgrid->jy - nowgrid->iy + 1;
			long nz = nowgrid->jz - nowgrid->iz + 1;
			DenType *oden = (DenType*)(outgrid + 1);
			DenType *nden = (DenType*)(nowgrid + 1);

			int i,j,k;
			for(k=iz;k<=jz;k++){
				long ko = k-outgrid->iz;
				long ki = k-nowgrid->iz;
				for(j=iy;j<=jy;j++){
					long jo = j-outgrid->iy;
					long ji = j-nowgrid->iy;
					for(i=ix;i<=jx;i++){
						long io = i-outgrid->ix;
						long ii = i-nowgrid->ix;
						/*
						outgrid->den[io+mx*(jo+my*ko)] += nowgrid->den[ii+nx*(ji+ny*ki)];
						*/
						oden[io+mx*(jo+my*ko)] += nden[ii+nx*(ji+ny*ki)];
					}
				}
			}
		}
		run += sizeof(GridInfo) + sizeof(DenType)*nowgrid->npix;
	}
}


void *UpdateNowMyData(void *nowmydata, void *buff2get){
	RGridInfo *a, *b;
	int i;
	a = (RGridInfo *)nowmydata;
	b = (RGridInfo *)buff2get;

	GridInfo *aa = (GridInfo*)((char*)nowmydata + sizeof(int));
	GridInfo *bb = (GridInfo*)((char*)buff2get + sizeof(int));

	char *cc = (char *)aa;
	for(i=0;i<a->ngrids;i++){
		cc += sizeof(GridInfo) + sizeof(DenType)*( ((GridInfo*)cc)->npix);
	}

	long orgsize = sizeof(int) + (cc - (char*)aa);
	a->ngrids += b->ngrids;


	/* */
	cc = (char*)bb;
	for(i=0;i<b->ngrids;i++){
		cc +=  sizeof(GridInfo) + sizeof(DenType)*( ((GridInfo*)cc)->npix);
	}
	long addsize = (cc - (char*)bb);
	nowmydata = (void *)realloc(nowmydata, addsize + orgsize);
	memmove( (char*)nowmydata+orgsize, bb, addsize);



	return nowmydata;
}

void *WrongUpdateNowMyData(void *nowmydata, void *buff2get, GridInfo mylocalgrid){
	int igrids=0;
	int i;

	int ngrids = ( (RGridInfo*)nowmydata)->ngrids;
	char *running = (char*)nowmydata + sizeof(int);
	char *overlapped = running;
	for(i=0;i<ngrids;i++){
		long npix = ((GridInfo*)running)->npix;
		long size = sizeof(GridInfo)+sizeof(DenType)* npix;
		if(OverLapGrid(running, &mylocalgrid))
		{
			long j;
			for(j=0;j<size;j++) *(overlapped++) = *(running++);
			igrids ++;
		}
		else running += size;
	}


	long nowsize = overlapped - (char*)nowmydata;

	ngrids = ( (RGridInfo*)buff2get)->ngrids;
	running = (char*)buff2get + sizeof(int);
	for(i=0;i<ngrids;i++){
		long npix = ((GridInfo*)running)->npix;
		running += sizeof(GridInfo)+sizeof(DenType)* npix;
	}
	long addsize = running - (char*)buff2get;
	nowmydata = (void *)realloc(nowmydata, nowsize + addsize);

	memmove((char*)nowmydata+nowsize, buff2get, addsize);






	free(buff2get);
	igrids += ngrids;
	( (RGridInfo*)nowmydata)->ngrids = igrids;

	return nowmydata;
}


void *dogmigrate(void *nowmydata, DoDeInfo *ddinfo, int nbuff){
	MPI_Status status;
	MPI_Comm com;
	int src,dest;
	long nrecv,nsend;
	long i,j,k;
	int myid,nid;
	long n_size = ddinfo->n_size;
	int igroup,subgroupsize,subgroupid;
	GridInfo mylocalgrid = ddinfo->sublocalgrid;
	GridInfo yourlocalgrid;

	myid = ddinfo->myid;
	nid = ddinfo->nid;
	com = ddinfo->com;

	int nsubgroup = ddinfo->nsubgroup;
	subgroupsize = ddinfo->subgroupsize;
	subgroupid = ddinfo->subgroupid;

	for(igroup=1;igroup<nsubgroup;igroup++){
		dest = (myid + subgroupsize*igroup + nid)%nid;
		src =  (myid - subgroupsize*igroup + nid)%nid;
		int targetsubgid = (subgroupid + igroup+nsubgroup)%nsubgroup;

		MPI_Sendrecv(&mylocalgrid,sizeof(GridInfo),MPI_BYTE,src,0,
				&yourlocalgrid,sizeof(GridInfo),MPI_BYTE,dest,0,
				com,&status);

		long sizeinbyte2send, sizeinbyte2recv;

		void *buff2send = ExtractGrids2Send(nowmydata, yourlocalgrid, &sizeinbyte2send);
		MPI_Sendrecv(&sizeinbyte2send,1,MPI_LONG,dest,0,
				&sizeinbyte2recv,1,MPI_LONG,src,0,
				com,&status);

		void *buff2get = (void *)malloc(sizeinbyte2recv);
		MPI_Sendrecv(buff2send,sizeinbyte2send,MPI_BYTE,dest,0,
				buff2get,sizeinbyte2recv,MPI_BYTE,src,0,
				com,&status);
		if(sizeinbyte2recv>0) 
			nowmydata = UpdateNowMyData(nowmydata,buff2get);


		if(buff2send) free(buff2send);
	}

	if(subgroupsize > 1) nowmydata = dogmigrate(nowmydata, ddinfo+1, nbuff);
	return nowmydata;

} 

RGridInfo *Grid2RGrid(GridInfo *ingrid){
	long size = sizeof(int) + sizeof(GridInfo) + sizeof(DenType)*ingrid->npix;
	long size2 = sizeof(GridInfo) + sizeof(DenType)*ingrid->npix;
	long i;
	char *nowmydata = (char*)realloc(ingrid,size);
	char *p = nowmydata + size;
	char *q = nowmydata + size2;
	for(i=0;i<size2;i++){
		*(--p) = *(--q);
	}
	RGridInfo *a = (RGridInfo*)nowmydata;
	a->ngrids = 1;
	return (RGridInfo*)nowmydata;
}


void getingridnpix(GridInfo *ingrid){
	ingrid->npix = (long)(ingrid->jx-ingrid->ix+1) * 
		(long)(ingrid->jy-ingrid->iy+1) * 
		(long)(ingrid->jz-ingrid->iz+1);
}
void *getoutgridnpix(GridInfo *outgrid, DenType initialval){
	outgrid->npix = (long)(outgrid->jx-outgrid->ix+1) * 
		(long)(outgrid->jy-outgrid->iy+1) * 
		(long)(outgrid->jz-outgrid->iz+1);
	outgrid = realloc(outgrid, sizeof(GridInfo) + sizeof(DenType)*outgrid->npix);
	long i;
	DenType *oden = (DenType *) (outgrid + 1);
	for(i=0;i<outgrid->npix;i++) {
		oden[i] = initialval;
	}
	return outgrid;
}


void gmigrate(GridInfo *ingrid, GridInfo *outgrid,DoDeInfo *ddinfo, int nddinfo, int nbuff){

	BuildGridHierarchyFromCurrentGridTopology(ddinfo, nddinfo, outgrid);

	RGridInfo *nowmydata = Grid2RGrid(ingrid);

	nowmydata = dogmigrate(nowmydata, ddinfo, nbuff);


	Dump2Outgrid(nowmydata, outgrid);

	free(nowmydata);
}
