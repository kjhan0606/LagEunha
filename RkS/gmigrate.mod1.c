/* 
   Not yet tested. PLEASE DO NOT USE THIS CODE */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<limits.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "mpirms.h"
#include "mpiaux.h"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))



void BuildGridHierarchyFromCurrentGridTopology(DoDeInfo *ddinfo, int nddinfo, 
		GridInfo *outgrid){
	int i;
	GridInfo localgrid;
	ddinfo[nddinfo-1].sublocalgrid = (localgrid = *outgrid);
	for(i=nddinfo-1;i>0;i--){
		GridInfo tlocalgrid[ddinfo[i].nid];
		MPI_Comm com = ddinfo[i].com;
		MPI_Gather(&localgrid,sizeof(GridInfo),MPI_BYTE,
				tlocalgrid,sizeof(GridInfo), MPI_BYTE, 0, com);
		if(ddinfo[i].myid==0) {
			int j;
			for(j=0;j<ddinfo[i].nid;j++){
				localgrid.ix = MIN(localgrid.ix,tlocalgrid[j].ix);
				localgrid.iy = MIN(localgrid.iy,tlocalgrid[j].iy);
				localgrid.iz = MIN(localgrid.iz,tlocalgrid[j].iz);
				localgrid.jx = MAX(localgrid.jx,tlocalgrid[j].jx);
				localgrid.jy = MAX(localgrid.jy,tlocalgrid[j].jy);
				localgrid.jz = MAX(localgrid.jz,tlocalgrid[j].jz);
			}
		}
		MPI_Bcast(&localgrid,sizeof(GridInfo), MPI_BYTE, 0, com);
		ddinfo[i-1].sublocalgrid = localgrid;
	}
}

#define FindOverLap(a,b,ix,iy,iz,jx,jy,jz,ip,lx,ly,lz) do{\
	lx = a->nx *(ip%3-1);\
	ly = a->ny *((ip%9/3)-1);\
	lz = a->nz *(ip/9-1);\
	ptrdiff_t bx = a->ix + lx;\
	ptrdiff_t by = a->iy + ly;\
	ptrdiff_t bz = a->iz + lz;\
	ptrdiff_t ux = a->jx + lx;\
	ptrdiff_t uy = a->jy + ly;\
	ptrdiff_t uz = a->jz + lz;\
	ix = MAX(bx,b->ix);\
	iy = MAX(by,b->iy);\
	iz = MAX(bz,b->iz);\
	jx = MIN(ux,b->jx);\
	jy = MIN(uy,b->jy);\
	jz = MIN(uz,b->jz);\
} while(0)

int OverLapGrid(const void *a, const void *b, int ip){
	GridInfo *grid = (GridInfo*)a;
	GridInfo *tgrid = (GridInfo*)b;
	int ix,iy,iz;
	int jx,jy,jz;
	int lx,ly,lz;
	FindOverLap(grid,tgrid,ix,iy,iz,jx,jy,jz,ip,lx,ly,lz);
	if(jx >=ix && jy >= iy && jz >= iz)
		return 1;
	else 
		return 0;
}


void *ExtractGrid2Send(GridInfo *grid, GridInfo *targetlocalgrid, GridInfo *write, int ip){
	int ix,iy,iz;
	int jx,jy,jz;
	int lx,ly,lz;
	FindOverLap(grid,targetlocalgrid,ix,iy,iz,jx,jy,jz,ip,lx,ly,lz);
	if(jx >=ix && jy >= iy && jz >= iz) 
	{
		ptrdiff_t mx = (jx-ix+1);
		ptrdiff_t my = (jy-iy+1);
		ptrdiff_t mz = (jz-iz+1);

		write->npix = (ptrdiff_t)mx*(ptrdiff_t)my*(ptrdiff_t)mz;
		write->nx = grid->nx;
		write->ny = grid->ny;
		write->nz = grid->nz;
		write->ix = ix;
		write->iy = iy;
		write->iz = iz;
		write->jx = jx;
		write->jy = jy;
		write->jz = jz;
		ptrdiff_t wx = grid->jx-grid->ix +1;
		ptrdiff_t wy = grid->jy-grid->iy +1;
		ptrdiff_t wz = grid->jz-grid->iz +1;
		ptrdiff_t i,j,k;
		DenType *den = (DenType*)(write+1);
		DenType *gden = (DenType*)(grid + 1);
		ptrdiff_t ni=0;
		for(k=iz;k<=jz;k++){
			ptrdiff_t kk = k -lz - grid->iz;
			for(j=iy;j<=jy;j++){
				ptrdiff_t jj = j -ly - grid->iy;
				for(i=ix;i<=jx;i++){
					ptrdiff_t ii = i -lx - grid->ix;
					den[ni++] = gden[ii+wx*(jj+wy*kk)];
				}
			}
		}
		return (den +ni);
	}
	else {
		/* to notify that there is nothing to send */
		targetlocalgrid->npix = 0;
		return write;
	}
}
ptrdiff_t  NPixInAGrid2Send(GridInfo *grid, GridInfo *targetlocalgrid, int ip){
	int ix,iy,iz;
	int jx,jy,jz;
	int lx,ly,lz;

	FindOverLap(grid,targetlocalgrid,ix,iy,iz,jx,jy,jz,ip,lx,ly,lz);
	if(jx >=ix && jy >= iy && jz >= iz) 
	{
		ptrdiff_t mx = (jx-ix+1);
		ptrdiff_t my = (jy-iy+1);
		ptrdiff_t mz = (jz-iz+1);
		return (mx*my*mz);
	}
	else return 0L;
}
void *ExtractGrids2Send(void *nowmydata, GridInfo targetlocalgrid, ptrdiff_t *grid_data_size, int ngrids_in_this_level){
	int i;
	ptrdiff_t npixels=0;
	void *sendbuff;
	char *myrun = (char *)nowmydata + sizeof(int);
	int igrid = 0;
	int ip=13;
	/* ngrids_in_this_level : To send grids segments excluding the grids transferred in this level of communication */
	for(i=0;i<ngrids_in_this_level;i++)
	{
		GridInfo *nowgrid = (GridInfo *)myrun;
		if(i==0) {
			for(ip=0;ip<27;ip++){
				ptrdiff_t mpixels = NPixInAGrid2Send(nowgrid, &targetlocalgrid,ip);
				if(mpixels>0) igrid++;
				npixels += mpixels;
			}
		}
		else {
			ip = 13;
			ptrdiff_t mpixels = NPixInAGrid2Send(nowgrid, &targetlocalgrid,ip);
			if(mpixels>0) igrid++;
			npixels += mpixels;
		}

		myrun += sizeof(GridInfo) + sizeof(DenType)*nowgrid->npix; 
	}
	ptrdiff_t sendsize = sizeof(int) + npixels*sizeof(DenType) + igrid*sizeof(GridInfo);

	sendbuff = (void*)malloc(sendsize);
	((RGridInfo*)sendbuff)->ngrids = igrid;
	myrun = (char *)nowmydata + sizeof(int);
	char *sendrun = (char*)sendbuff + sizeof(int);
	for(i=0;i<ngrids_in_this_level;i++)
	{
		GridInfo *nowgrid = (GridInfo *)myrun;
		if(i==0){
			for(ip=0;ip<27;ip++){
				sendrun = (char *)ExtractGrid2Send(nowgrid, &targetlocalgrid, (GridInfo*)sendrun,ip);
			}
		}
		else {
			ip = 13;
			sendrun = (char *)ExtractGrid2Send(nowgrid, &targetlocalgrid, (GridInfo*)sendrun,ip);
		}
		myrun += sizeof(GridInfo) + sizeof(DenType)*nowgrid->npix;
	}
	*grid_data_size = (char *)sendrun - (char*)sendbuff;
	return sendbuff;
}
void Dump2Outgrid(void *nowmydata, GridInfo *outgrid){
	int ngrids = ((RGridInfo*)nowmydata)->ngrids;
	char *run = (char*)nowmydata + sizeof(int);

	ptrdiff_t mx = outgrid->jx - outgrid->ix + 1;
	ptrdiff_t my = outgrid->jy - outgrid->iy + 1;
	ptrdiff_t mz = outgrid->jz - outgrid->iz + 1;
	outgrid->npix = mx*my*mz;

	ptrdiff_t ijk;
	for(ijk=0;ijk<ngrids;ijk++){
		GridInfo *nowgrid = (GridInfo*) run;
		int ix,iy,iz;
		int jx,jy,jz;
		int lx,ly,lz;
		int ip,is;
		int step;
		if(ijk==0) { /* First grid contains the original ingrid data and therefore it should be periodically shifted 
					  to check the overlap with the requested grids*/
			is = 0;
			step = 1;
		}
		else { /* Rest grids are the trasferred grids and already reflect the peridodic conditions */ 
			is = 13;
			step = 10000;
		}
		for(ip=is;ip<27;ip+=step){
			FindOverLap(nowgrid,outgrid,ix,iy,iz,jx,jy,jz,ip,lx,ly,lz);
			if(jx >=ix && jy >= iy && jz >=iz) 
			{
				ptrdiff_t nx = nowgrid->jx - nowgrid->ix + 1;
				ptrdiff_t ny = nowgrid->jy - nowgrid->iy + 1;
				ptrdiff_t nz = nowgrid->jz - nowgrid->iz + 1;
				DenType *oden = (DenType*)(outgrid + 1);
				DenType *nden = (DenType*)(nowgrid + 1);
	
				int i,j,k;
				for(k=iz;k<=jz;k++){
					ptrdiff_t ko = k-outgrid->iz;
					ptrdiff_t ki = k-lz-nowgrid->iz;
					for(j=iy;j<=jy;j++){
						ptrdiff_t jo = j-outgrid->iy;
						ptrdiff_t ji = j-ly-nowgrid->iy;
						for(i=ix;i<=jx;i++){
							ptrdiff_t io = i-outgrid->ix;
							ptrdiff_t ii = i-lx-nowgrid->ix;
							oden[io+mx*(jo+my*ko)] += nden[ii+nx*(ji+ny*ki)];
						}
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

	ptrdiff_t orgsize = sizeof(int) + (cc - (char*)aa);
	a->ngrids += b->ngrids;


	/* */
	cc = (char*)bb;
	for(i=0;i<b->ngrids;i++){
		cc +=  sizeof(GridInfo) + sizeof(DenType)*( ((GridInfo*)cc)->npix);
	}
	ptrdiff_t addsize = (cc - (char*)bb);
	nowmydata = (void *)realloc(nowmydata, addsize + orgsize);
	memmove( (char*)nowmydata+orgsize, bb, addsize);



	return nowmydata;
}




RGridInfo *Grid2RGrid(GridInfo *ingrid){
	ptrdiff_t size = sizeof(int) + sizeof(GridInfo) + sizeof(DenType)*ingrid->npix;
	ptrdiff_t size2 = sizeof(GridInfo) + sizeof(DenType)*ingrid->npix;
	ptrdiff_t i;
	char *nowmydata = (char*)realloc(ingrid,size);

	{
		char *temp = (char*)malloc(size2);
		memmove(temp,       nowmydata, size2);
		memmove(nowmydata+sizeof(int) ,  temp, size2);
		free(temp);
	}

	RGridInfo *a = (RGridInfo*)nowmydata;
	a->ngrids = 1;
	return (RGridInfo*)nowmydata;
}


void initializeGrid(GridInfo *ingrid, DenType initialval){
	ingrid->npix = (ptrdiff_t)(ingrid->jx-ingrid->ix+1) * 
		(ptrdiff_t)(ingrid->jy-ingrid->iy+1) * 
		(ptrdiff_t)(ingrid->jz-ingrid->iz+1);
	ptrdiff_t i;
	DenType *den = (DenType*)(ingrid+1);
	for(i=0;i<ingrid->npix;i++){
		den[i] = initialval;
	}
}
void *getoutgridnpix(GridInfo *outgrid, DenType initialval){
	outgrid->npix = (ptrdiff_t)(outgrid->jx-outgrid->ix+1) * (ptrdiff_t)(outgrid->jy-outgrid->iy+1) * (ptrdiff_t)(outgrid->jz-outgrid->iz+1);
	outgrid = realloc(outgrid, sizeof(GridInfo) + sizeof(DenType)*outgrid->npix);
	DenType *oden = (DenType *) (outgrid + 1);
	ptrdiff_t i;
	for(i=0;i<outgrid->npix;i++) {
		oden[i] = initialval;
	}
	return outgrid;
}

void FillVolume(GridInfo *ingrid, GridInfo *outgrid){
	ptrdiff_t i,j,k;
	ptrdiff_t mx,my;
	ptrdiff_t nx,ny;
	ptrdiff_t ix,iy,iz,jx,jy,jz;
	ptrdiff_t lx,ly,lz;
	int ip;
	nx = outgrid->jx - outgrid->ix;
	ny = outgrid->jy - outgrid->iy;
	mx = ingrid->jx - ingrid->ix;
	my = ingrid->jy - ingrid->iy;
	ip= 13;
	{
		FindOverLap(ingrid,outgrid,ix,iy,iz,jx,jy,jz,ip,lx,ly,lz);
		if(jx >=ix && jy >= iy && jz >= iz){
			for(k=iz;k<=jz;k++) {
				ptrdiff_t shiftz1 = nx*ny*(k-outgrid->iz);
				ptrdiff_t shiftz2 = mx*my*(k-ingrid->iz-lz);
				for(j=iy;j<=jy;j++) { 
					ptrdiff_t shifty1 = nx*(j-outgrid->iy);
					ptrdiff_t shifty2 = mx*(j-ingrid->iy-ly);
					outgrid[ix-outgrid->ix + nx*(j-outgrid->iy+ny*(k-outgrid->iz))] += 
						ingrid[ix-ingrid->ix-lx+mx*(j-ingrid->iy-ly + my*(k-ingrid->iz-lz))];
				}
			}
		}
	}
}
void FillPeriodicVolume(GridInfo *ingrid, GridInfo *outgrid){
	ptrdiff_t i,j,k;
	ptrdiff_t mx,my;
	ptrdiff_t nx,ny;
	ptrdiff_t ix,iy,iz,jx,jy,jz;
	ptrdiff_t lx,ly,lz;
	int ip;
	nx = outgrid->jx - outgrid->ix;
	ny = outgrid->jy - outgrid->iy;
	mx = ingrid->jx - ingrid->ix;
	my = ingrid->jy - ingrid->iy;
	for(ip=0;ip<27;ip++){
		FindOverLap(ingrid,outgrid,ix,iy,iz,jx,jy,jz,ip,lx,ly,lz);
		if(jx >=ix && jy >= iy && jz >= iz){
			for(k=iz;k<=jz;k++) {
				ptrdiff_t shiftz1 = nx*ny*(k-outgrid->iz);
				ptrdiff_t shiftz2 = mx*my*(k-ingrid->iz-lz);
				for(j=iy;j<=jy;j++) { 
					ptrdiff_t shifty1 = nx*(j-outgrid->iy);
					ptrdiff_t shifty2 = mx*(j-ingrid->iy-ly);
					outgrid[ix-outgrid->ix + nx*(j-outgrid->iy+ny*(k-outgrid->iz))] += 
						ingrid[ix-ingrid->ix-lx+mx*(j-ingrid->iy-ly + my*(k-ingrid->iz-lz))];
				}
			}
		}
	}
}
void *dogmigrate(void *nowmydata, DoDeInfo *ddinfo){
	MPI_Status status;
	MPI_Comm com;
	int src,dest;
	ptrdiff_t nrecv,nsend;
	ptrdiff_t i,j,k;
	int myid,nid;
	ptrdiff_t n_size = ddinfo->n_size;
	int igroup,subgroupsize,subgroupid;
	GridInfo mylocalgrid = ddinfo->sublocalgrid;
	GridInfo yourlocalgrid;

	myid = ddinfo->myid;
	nid = ddinfo->nid;
	com = ddinfo->com;

	int nsubgroup = ddinfo->nsubgroup;
	subgroupsize = ddinfo->subgroupsize;
	subgroupid = ddinfo->subgroupid;


	int ngrids_in_this_level = ((RGridInfo*)nowmydata)->ngrids;

	for(igroup=1;igroup<nsubgroup;igroup++){
		dest = (myid + subgroupsize*igroup + nid)%nid;
		src =  (myid - subgroupsize*igroup + nid)%nid;
		int targetsubgid = (subgroupid + igroup+nsubgroup)%nsubgroup;

		MPI_Sendrecv(&mylocalgrid,sizeof(GridInfo),MPI_BYTE,src,0,
				&yourlocalgrid,sizeof(GridInfo),MPI_BYTE,dest,0,
				com,&status);

		ptrdiff_t sizeinbyte2send, sizeinbyte2recv;

		void *buff2send = ExtractGrids2Send(nowmydata, yourlocalgrid, &sizeinbyte2send, ngrids_in_this_level);
		MPI_Sendrecv(&sizeinbyte2send,1,MPI_INT64_T,dest,0,
				&sizeinbyte2recv,1,MPI_INT64_T,src,0,
				com,&status);

		void *buff2get = (void *)malloc(sizeinbyte2recv);
		MPI_Sendrecv(buff2send,sizeinbyte2send,MPI_BYTE,dest,0,
				buff2get,sizeinbyte2recv,MPI_BYTE,src,0,
				com,&status);
		if(sizeinbyte2recv>0) 
			nowmydata = UpdateNowMyData(nowmydata,buff2get);


		if(buff2send != NULL) free(buff2send);
		if(buff2get != NULL) free(buff2get);
	}

	if(subgroupsize > 1) nowmydata = dogmigrate(nowmydata, ddinfo+1);
	return nowmydata;
} 


void gmigrate(GridInfo *ingrid, GridInfo *outgrid,DoDeInfo *ddinfo, int nddinfo){
	BuildGridHierarchyFromCurrentGridTopology(ddinfo, nddinfo, outgrid);

	initializeGrid(outgrid, 0.L);

	FillPeriodicVolume(ingrid, outgrid);

	RGridInfo *nowmydata = Grid2RGrid(ingrid);

	nowmydata = dogmigrate(nowmydata, ddinfo);

	Dump2Outgrid(nowmydata, outgrid);

	free(nowmydata);
}
