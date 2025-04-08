/* alpha version of parallel domain deocomposition with the Recursive MultiSection (RMS).
   It now allows only the three-dimensional data.
   Further implementation is needed to allow for arbitrary dimensional data. */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<limits.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "mpirks.h"
#include "mpiaux.h"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

#define findOverlap(a,b,ix,iy,iz,jx,jy,jz,ip,lx,ly,lz) do{\
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


void BuildGridHierarchyFromCurrentGridTopology(DoDeInfo *ddinfo, int nddinfo, 
		GridInfo *outgrid){
	int i;
	GridInfo localgrid;
	ddinfo[nddinfo-1].sublocalgrid = (localgrid = *outgrid);
#ifdef DEBUG2
	DEBUGPRINT("%d : P%d has %d %d :: %d %d :: %d %d\n",nddinfo-1,ddinfo[0].myid,
				localgrid.ix,localgrid.jx, localgrid.iy,localgrid.jy, localgrid.iz,localgrid.jz);
#endif
	for(i=nddinfo-2;i>=0;i--){
		GridInfo tlocalgrid[ddinfo[i+1].nid];
		MPI_Comm com = ddinfo[i+1].com;
		MPI_Gather(&localgrid,sizeof(GridInfo),MPI_BYTE,
				tlocalgrid,sizeof(GridInfo), MPI_BYTE, 0, com);
		if(ddinfo[i+1].myid==0) {
			int j;
			for(j=1;j<ddinfo[i+1].nid;j++){
				localgrid.ix = MIN(localgrid.ix,tlocalgrid[j].ix);
				localgrid.iy = MIN(localgrid.iy,tlocalgrid[j].iy);
				localgrid.iz = MIN(localgrid.iz,tlocalgrid[j].iz);
				localgrid.jx = MAX(localgrid.jx,tlocalgrid[j].jx);
				localgrid.jy = MAX(localgrid.jy,tlocalgrid[j].jy);
				localgrid.jz = MAX(localgrid.jz,tlocalgrid[j].jz);

				localgrid.npix += tlocalgrid[j].npix;
			}
		}
		MPI_Bcast(&localgrid,sizeof(GridInfo), MPI_BYTE, 0, com);
		ddinfo[i].sublocalgrid = localgrid;
#ifdef DEBUG2
		DEBUGPRINT("%d : P%d has %d %d :: %d %d :: %d %d\n",i,ddinfo[0].myid,
				localgrid.ix,localgrid.jx, localgrid.iy,localgrid.jy, localgrid.iz,localgrid.jz);
#endif
	}
}
int isGridOverlapped(const void *a, const void *b, int ip){
	GridInfo *grid = (GridInfo*)a;
	GridInfo *tgrid = (GridInfo*)b;
	ptrdiff_t ix,iy,iz;
	ptrdiff_t jx,jy,jz;
	ptrdiff_t lx,ly,lz;

	findOverlap(grid,tgrid,ix,iy,iz,jx,jy,jz,ip,lx,ly,lz);
	if(jx >=ix && jy >= iy && jz >= iz)
		return 1;
	else 
		return 0;
}


void *ExtractGrid2Send(GridInfo *grid, GridInfo *targetlocalgrid, GridInfo *write, int ip){
	ptrdiff_t ix,iy,iz;
	ptrdiff_t jx,jy,jz;
	ptrdiff_t lx,ly,lz;
	findOverlap(grid,targetlocalgrid,ix,iy,iz,jx,jy,jz,ip,lx,ly,lz);
	if(jx >=ix && jy >= iy && jz >= iz) {
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
				ptrdiff_t yzoffset = wx*(jj+wy*kk) -lx -grid->ix;
				for(i=ix;i<=jx;i++){
					/*
					ptrdiff_t ii = i -lx - grid->ix;
					den[ni++] = gden[ii+wx*(jj+wy*kk)];
					*/
					den[ni++] = gden[i+yzoffset];
				}
			}
		}
		return (den +ni);
	}
	else {
		/* to notify that there is nothing to send */
		/*
		targetlocalgrid->npix = 0;
		*/
		return write;
	}
}
ptrdiff_t  CountPixels2Send(GridInfo *grid, GridInfo *targetlocalgrid, int ip){
	ptrdiff_t ix,iy,iz;
	ptrdiff_t jx,jy,jz;
	ptrdiff_t lx,ly,lz;

	findOverlap(grid,targetlocalgrid,ix,iy,iz,jx,jy,jz,ip,lx,ly,lz);
	if(jx >=ix && jy >= iy && jz >= iz) {
		ptrdiff_t mx = (jx-ix+1);
		ptrdiff_t my = (jy-iy+1);
		ptrdiff_t mz = (jz-iz+1);
		/*
#ifdef DEBUG
		DEBUGPRINT("countpixels2send %ld %ld %ld : %ld %ld %ld : %ld %ld %ld\n", 
				ix,iy,iz,jx,jy,jz,mx,my,mz );
#endif
*/
		return (mx*my*mz);
	}
	else return 0L;
}
void *extractGrids2Send(void *nowmydata, GridInfo targetlocalgrid, ptrdiff_t *grid_data_size, int ngrids_in_this_level){
	int i;
	ptrdiff_t npixels=0;
	void *sendbuff;
	char *myrun = (char *)nowmydata + sizeof(int);
	int igrid = 0;
	int ip=13; /* self */

	/* ngrids_in_this_level : To send the grid part except for the grids 
	 * already got in this OST level of communication */
	for(i=0;i<ngrids_in_this_level;i++) {
		GridInfo *nowgrid = (GridInfo *)myrun;
		if(i==0) {
			for(ip=0;ip<27;ip++){
				ptrdiff_t mpixels = CountPixels2Send(nowgrid, &targetlocalgrid,ip);
				if(mpixels>0) igrid++;
				npixels += mpixels;
			}
		}
		else {
			ip = 13;
			ptrdiff_t mpixels = CountPixels2Send(nowgrid, &targetlocalgrid,ip);
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
	for(i=0;i<ngrids_in_this_level;i++) {
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
#ifdef DEBUG2
	DEBUGPRINT("the number of grids to send= %d and npixels= %ld with tsize= %ld\n", ((RGridInfo*)sendbuff)->ngrids, npixels, *grid_data_size);
#endif

	return sendbuff;
}
void Dump2Outgrid(void *nowmydata, GridInfo *outgrid){
	int ngrids = ((RGridInfo*)nowmydata)->ngrids;
	char *run = (char*)nowmydata + sizeof(int);

	ptrdiff_t mx = outgrid->jx - outgrid->ix + 1;
	ptrdiff_t my = outgrid->jy - outgrid->iy + 1;
	ptrdiff_t mz = outgrid->jz - outgrid->iz + 1;

	if(mx<0 || my <0 || mz <0) {
		outgrid->npix = 0;
		return;
	}
	outgrid->npix = mx*my*mz;

	ptrdiff_t ijk;
	for(ijk=0;ijk<ngrids;ijk++){
		GridInfo *nowgrid = (GridInfo*) run;
		int ix,iy,iz;
		int jx,jy,jz;
		int lx,ly,lz;
		int ip,is;
		int step;
		if(ijk==0) { /* First grid contains the original ingrid data and 
						therefore it should be periodically shifted 
					  to check the overlap with the requested grids*/
			is = 0;
			step = 1;
		}
		else { /* The rests are the trasferred grids and already reflect the peridodic conditions */ 
			is = 13;
			step = 10000;
		}
		for(ip=is;ip<27;ip+=step){
			findOverlap(nowgrid,outgrid,ix,iy,iz,jx,jy,jz,ip,lx,ly,lz);
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
						ptrdiff_t offset1, offset2;
						offset1 = mx*(jo+my*ko) - outgrid->ix;
						offset2 = nx*(ji+ny*ki) -lx-nowgrid->ix;
						for(i=ix;i<=jx;i++){
							oden[i+offset1] += nden[i+offset2];
						}
					}
				}
			}
		}
		run += sizeof(GridInfo) + sizeof(DenType)*nowgrid->npix;
	}
}


void *UpdateNowMyData(DoDeInfo *ddinfo, void *nowmydata, void *buff2get){
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

	ptrdiff_t orgsizeinbyte = sizeof(int) + (cc - (char*)aa);
	a->ngrids += b->ngrids;

	/* */
	cc = (char*)bb;
	for(i=0;i<b->ngrids;i++){
#ifdef DEBUG2
		DEBUGPRINT("P%d/%d has G%d/%d npix= %ld to send\n", 
				ddinfo->myid, ddinfo->nid, i, b->ngrids,((GridInfo*)cc)->npix );
#endif
		cc +=  sizeof(GridInfo) + sizeof(DenType)*( ((GridInfo*)cc)->npix);
	}
	ptrdiff_t addsizeinbyte = (cc - (char*)bb);
	nowmydata = (void *)realloc(nowmydata, addsizeinbyte + orgsizeinbyte);
	memcpy( (char*)nowmydata+orgsizeinbyte, bb, addsizeinbyte);

	return nowmydata;
}
ptrdiff_t PickUpOverLapGrid(GridInfo *grid, GridInfo *tgrid){
	ptrdiff_t ix,iy,iz;
	ptrdiff_t jx,jy,jz;
	ptrdiff_t mx,my,mz,nx,ny;
	ptrdiff_t i,j,k;
	ptrdiff_t nnx,nny,nnz;
	nnx = grid->nx;
	nny = grid->ny;
	nnz = grid->nz;
	ix = MAX(grid->ix,tgrid->ix);
	iy = MAX(grid->iy,tgrid->iy);
	iz = MAX(grid->iz,tgrid->iz);
	jx = MIN(grid->jx,tgrid->jx);
	jy = MIN(grid->jy,tgrid->jy);
	jz = MIN(grid->jz,tgrid->jz);

	nx = grid->jx-grid->ix + 1;
	ny = grid->jy-grid->iy + 1;

	mx = (jx-ix+1);
	my = (jy-iy+1);
	mz = (jz-iz+1);
	ptrdiff_t npixels = mx*my*mz;
	DenType *den, *gden;
	den = (DenType*)(grid+1);
	gden = den;
	for(k=iz;k<=jz;k++){
		ptrdiff_t k1 = k - grid->iz;
		ptrdiff_t k2 = k - iz;
		for(j=iy;j<=jy;j++){
			ptrdiff_t j1 = j - grid->iy;
			ptrdiff_t j2 = j - iy;
			/*
			memmove(den+mx*(j2+my*k2), gden+ix-grid->ix+nx*(j1+ny*k1), mx * sizeof(DenType));
			*/
			ptrdiff_t i1, i2;
			i1 = mx*(j2+my*k2);
			i2 = ix-grid->ix + nx*(j1+ny*k1);
			for(i=0;i<mx;i++){
				den[i+i1] = gden[i+i2];
			}
		}
	}
	grid->npix = npixels;
	grid->ix = ix;
	grid->jx = jx;
	grid->iy = iy;
	grid->jy = jy;
	grid->iz = iz;
	grid->jz = jz;
	grid->nx = nnx;
	grid->ny = nny;
	grid->nz = nnz;
	grid->nxny = nnx*nny;
	return (npixels*sizeof(DenType) + sizeof(GridInfo));
}
void CopyGrid(GridInfo *left, GridInfo *right){
	*left = *right;
	ptrdiff_t i,npix = right->npix;
	DenType *a,*b;
	a = (DenType*)(left+1);
	b = (DenType*)(right+1);
	for(i=0;i<npix;i++){
		*(a++) = *(b++);
	}
}

void *DeleteOutofSubLocalGrid(void *nowmydata, DoDeInfo *ddinfo){
	GridInfo *mylocalgrid = &(ddinfo->sublocalgrid);
	GridInfo *nowgrid;
	int ngrids = ((RGridInfo*)nowmydata)->ngrids;
	int igrids,sgrids=0;
	char *right = (char *)nowmydata + sizeof(int);
	char *left = right;
	for(igrids=0;igrids<ngrids;igrids++){
		nowgrid = (GridInfo*) right;
		ptrdiff_t size = sizeof(GridInfo)+ sizeof(DenType)*((GridInfo*)right)->npix;
		if(igrids==0) {
			left += size;
			sgrids++;
		}
		else {
			if(isGridOverlapped( (const void*)right, (const void*) mylocalgrid, 13)){
				ptrdiff_t msize = PickUpOverLapGrid((GridInfo*)right, mylocalgrid);
				CopyGrid((GridInfo*)left, (GridInfo*)right);
				/*
				memmove(left, right, msize);
				*/
				left += msize;
				sgrids++;
			}
		}
		right += size;
	}
	((RGridInfo*)nowmydata)->ngrids = sgrids;
	nowmydata = (void *)realloc(nowmydata, left-(char*)nowmydata);
	return nowmydata;
}
/* grid migration */

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


#ifdef DEBUG2
	DEBUGPRINT("p%d/%d subgroupid/nsubgroup/subgroupsize= %d %d %d\n", myid, nid,
			subgroupid, nsubgroup, subgroupsize);
#endif
	int ngrids_in_this_level = ((RGridInfo*)nowmydata)->ngrids;

	for(igroup=1;igroup<nsubgroup;igroup++){
		dest = (myid + subgroupsize*igroup + nid)%nid;
		src =  (myid - subgroupsize*igroup + nid)%nid;
//		int targetsubgid = (subgroupid + igroup+nsubgroup)%nsubgroup;

		MPI_Sendrecv(&mylocalgrid,sizeof(GridInfo),MPI_BYTE,src,0,
				&yourlocalgrid,sizeof(GridInfo),MPI_BYTE,dest,0,
				com,&status);

		ptrdiff_t sizeinbyte2send, sizeinbyte2recv;

		void *buff2send = NULL;
		buff2send = extractGrids2Send(nowmydata, yourlocalgrid, &sizeinbyte2send, 
				ngrids_in_this_level);
		MPI_Sendrecv(&sizeinbyte2send,1,MPI_INT64_T,dest,1,
				&sizeinbyte2recv,1,MPI_INT64_T,src,1,
				com,&status);
#ifdef DEBUG2
		{
			RGridInfo *a = (RGridInfo*)buff2send;
			DEBUGPRINT("P%d/%d has ngrids to send %d with sizes = %ld\n", ddinfo->myid, ddinfo->nid,
					a->ngrids, sizeinbyte2send);
		}
#endif

		void *buff2get = NULL;
		buff2get = (void *)malloc(sizeinbyte2recv);

		my_MPI_Sendrecv(buff2send,sizeinbyte2send,MPI_BYTE,dest,10,
				buff2get,sizeinbyte2recv,MPI_BYTE,src,10,
				com,&status);
#ifdef DEBUG2
		{
			RGridInfo *a = (RGridInfo *)buff2send;
			RGridInfo *b = (RGridInfo *)buff2get;
			DEBUGPRINT("P%d/%d has ngrids of buff2send and buff2get %d %d with sizes= %ld %ld\n", 
					ddinfo->myid, ddinfo->nid, 
					a->ngrids, b->ngrids,
					sizeinbyte2send, sizeinbyte2recv);
		}
#endif
		MPI_Barrier(MPI_COMM_WORLD);
		if(sizeinbyte2recv>4)  // This is for the case of ngrids >=1
			nowmydata = UpdateNowMyData(ddinfo, nowmydata,buff2get);

#ifdef DEBUG2
		MPI_Barrier(com);
		ptrdiff_t size = 0;
		GridInfo *run = (GridInfo*)((int*)nowmydata + 1);
		for(i=0;i<((RGridInfo*)nowmydata)->ngrids;i++){
			size += run->npix;
			run = (GridInfo*) ( (char*) run + sizeof(GridInfo) + sizeof(DenType)*run->npix);
		}
		DEBUGPRINT("p%d has sizeof in bytes: %ld / %ld to recv/send from/to P%d P%d  with nblocks of grids %d %d on the current total size of %ld::: igroup/nsubgroup/subgroupsize = %d / %d / %d\n",
				myid, sizeinbyte2recv , sizeinbyte2send, src, dest,
				((RGridInfo*)buff2get)->ngrids,
				((RGridInfo*)buff2send)->ngrids,
				size,
				igroup, nsubgroup,subgroupsize
				);
#endif

		if(buff2send != NULL) free(buff2send);
		if(buff2get != NULL) free(buff2get);
	}
	nowmydata = DeleteOutofSubLocalGrid(nowmydata, ddinfo);

	if(subgroupsize > 1) nowmydata = dogmigrate(nowmydata, ddinfo+1);
	return nowmydata;
} 

RGridInfo *Grid2RGrid(GridInfo *ingrid){
	ptrdiff_t size = sizeof(int) + sizeof(GridInfo) + sizeof(DenType)*ingrid->npix;
	ptrdiff_t size2 = sizeof(GridInfo) + sizeof(DenType)*ingrid->npix;
	ptrdiff_t i;
	char *nowmydata = (char*)realloc(ingrid,size);

	if(0){
		char *p = nowmydata + size;
		char *q = nowmydata + size2;
		for(i=0;i<size2;i++){
			*(--p) = *(--q);
		}
	}
	else {
		char *temp = (char*)malloc(size2);
		memcpy(temp,       nowmydata, size2);
		memcpy(nowmydata+sizeof(int) ,  temp, size2);
		free(temp);
	}

	RGridInfo *a = (RGridInfo*)nowmydata;
	a->ngrids = 1;
	return (RGridInfo*)nowmydata;
}


void getingridnpix(GridInfo *ingrid){
	ingrid->npix = (ptrdiff_t)(ingrid->jx-ingrid->ix+1) * 
		(ptrdiff_t)(ingrid->jy-ingrid->iy+1) * 
		(ptrdiff_t)(ingrid->jz-ingrid->iz+1);
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

void gmigrate(GridInfo *ingrid, GridInfo *outgrid,DoDeInfo *ddinfo, int nddinfo){

#ifdef DEBUG2
	MPI_Barrier(ddinfo->com);
	DEBUGPRINT("P%d is here 1\n",ddinfo->myid);
#endif

	BuildGridHierarchyFromCurrentGridTopology(ddinfo, nddinfo, outgrid);
#ifdef DEBUG2
	MPI_Barrier(ddinfo->com);
	DEBUGPRINT("P%d is here 2\n",ddinfo->myid);
#endif

	RGridInfo *nowmydata = Grid2RGrid(ingrid);
#ifdef DEBUG2
	MPI_Barrier(ddinfo->com);
	DEBUGPRINT("P%d is here 3\n",ddinfo->myid);
#endif

	nowmydata = dogmigrate(nowmydata, ddinfo);

#ifdef DEBUG2
	MPI_Barrier(ddinfo->com);
	DEBUGPRINT("P%d is here 4\n",ddinfo->myid);
#endif
	Dump2Outgrid(nowmydata, outgrid);
#ifdef DEBUG2
	MPI_Barrier(ddinfo->com);
	DEBUGPRINT("P%d is here 5\n",ddinfo->myid);
#endif

	free(nowmydata);
}
/*
RGridInfo *dogmigrate0(RGridInfo *nowmydata, DoDeInfo *ddinfo, int istep){
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


	int ngrids_in_this_level = nowmydata->ngrids;

	for(igroup=1;igroup<nsubgroup;igroup++){
		dest = (myid + subgroupsize*igroup + nid)%nid;
		src =  (myid - subgroupsize*igroup + nid)%nid;
		int targetsubgid = (subgroupid + igroup+nsubgroup)%nsubgroup;

		MPI_Sendrecv(&mylocalgrid,sizeof(GridInfo),MPI_BYTE,src,0,
				&yourlocalgrid,sizeof(GridInfo),MPI_BYTE,dest,0,
				com,&status);

		ptrdiff_t sizeinbyte2send, sizeinbyte2recv;

		void *buff2send = extractGrids2Send(nowmydata, yourlocalgrid, &sizeinbyte2send, 
				ngrids_in_this_level);
		MPI_Sendrecv(&sizeinbyte2send,1,MPI_INT64_T,dest,0,
				&sizeinbyte2recv,1,MPI_INT64_T,src,0,
				com,&status);

		void *buff2get = (void *)malloc(sizeinbyte2recv);
		MPI_Sendrecv(buff2send,sizeinbyte2send,MPI_BYTE,dest,10,
				buff2get,sizeinbyte2recv,MPI_BYTE,src,10,
				com,&status);
		if(sizeinbyte2recv>0) nowmydata = UpdateNowMyData(ddinfo,nowmydata,buff2get);

		ptrdiff_t size = 0;
		GridInfo *run = (GridInfo*)((int*)nowmydata + 1);
		for(i=0;i<nowmydata->ngrids;i++){
			size += run->npix;
			run = (GridInfo*) ( (char*) run + sizeof(GridInfo) + sizeof(DenType)*run->npix);
		}
#ifdef DEBUG
		printf("%d :: P%d has npixels %ld / %ld to recv/send from/to P%d P%d  with ngrids %d %d with current total size of %ld\n",
				istep,
				myid, sizeinbyte2recv , sizeinbyte2send, src, dest,
				((RGridInfo*)buff2get)->ngrids,
				((RGridInfo*)buff2send)->ngrids,
				size
				);
#endif

		free(buff2send);
		free(buff2get);
	}
	nowmydata = DeleteOutofSubLocalGrid(nowmydata, ddinfo);
	if(istep ==0) {
		MPI_Finalize();
		exit(123);
	}

	if(subgroupsize > 1) nowmydata = dogmigrate0(nowmydata, ddinfo+1,istep+1);
	return nowmydata;
} 

void gmigrate0(GridInfo *ingrid, GridInfo *outgrid,DoDeInfo *ddinfo, int nddinfo){

#ifdef DEBUG
	DEBUGPRINT("P%d is here 1\n",ddinfo->myid);
#endif
	BuildGridHierarchyFromCurrentGridTopology(ddinfo, nddinfo, outgrid);
#ifdef DEBUG
	DEBUGPRINT("P%d is here 2\n",ddinfo->myid);
#endif

	RGridInfo *nowmydata = Grid2RGrid(ingrid);
#ifdef DEBUG
	DEBUGPRINT("P%d is here 3\n",ddinfo->myid);
#endif
	nowmydata = dogmigrate0(nowmydata, ddinfo, 0);

#ifdef DEBUG
	DEBUGPRINT("P%d is here 4\n",ddinfo->myid);
#endif
	Dump2Outgrid(nowmydata, outgrid);
#ifdef DEBUG
	DEBUGPRINT("P%d is here 5\n",ddinfo->myid);
#endif

	free(nowmydata);
}
*/
