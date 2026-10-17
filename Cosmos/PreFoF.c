#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<omp.h>
#include<mpi.h>
#include "eunha.h"
#include "params.h"
#include "ost.h"
#include "BasicCell.h"
#include "PreFoF.h"


#define MIN(a,b) ( (a)<(b) ? (a):(b) )
#define MAX(a,b) ( (a)>(b) ? (a):(b) )


typedef struct HaloBound{
	size_t nmem;
	PosType xmin,ymin,zmin,xmax,ymax,zmax;
	int dumpflag;
	FoFTPtlStruct *sibling;
} HaloBound;




ptrdiff_t Dump2FoFPtl(SimParameters *simpar, size_t ix,
		size_t iy, size_t iz, size_t mx, size_t my, size_t mz, FoFTPtlStruct **point,
		PosType xmin, PosType ymin, PosType zmin){
	size_t npoint = 0;
	TreeLinkedCell *BasicCell = BASICCELL(simpar);

#ifdef XYZDBL
	double cellxpos, cellypos, cellzpos;
#else
	PosType cellxpos, cellypos, cellzpos;
#endif
	cellxpos = CellWidth*ix + xmin;
	cellypos = CellWidth*iy + ymin;
	cellzpos = CellWidth*iz + zmin;

	size_t i,j,k;

	npoint = BasicCell[ix+mx*(iy+my*iz)].nmem;
	*point = (FoFTPtlStruct*)malloc(sizeof(FoFTPtlStruct)*npoint);
	FoFTPtlStruct *g = *point;
	linkedlisttype *tmp = BasicCell[ix+mx*(iy+my*iz)].link;
	npoint = 0;
	while(tmp){
		g[npoint].x = XofP(simpar, tmp) - cellxpos;
		g[npoint].y = YofP(simpar, tmp) - cellypos;
		g[npoint].z = ZofP(simpar, tmp) - cellzpos;
		g[npoint].bp = tmp;
		g[npoint].haloindx = -1;
		tmp = tmp->next;
		npoint ++;
	}
	return npoint;
}

#define DumpP2Disk_TRAD(simpar, TYPE, type, nfiles, nsave, tnsave, headeroffset) do{\
	for(i=0;i<TYPE##_NP(simpar);i++) UNSET_P_FLAG(simpar, TYPE, i, FoFflag);\
	ptrdiff_t offset = sizeof(struct linkedlisttype*);\
	ptrdiff_t isave=0;\
	MPI_Status status;\
	int itag = 1;\
	int isend, iget, WGroupSize = WGROUPSIZE(simpar);\
	int src = MYID(simpar) - 1;\
	int tgt = MYID(simpar) + 1;\
	type##particletype *tmp;\
	int NSAVE = (NSave+NID(simpar)-1)/NID(simpar);\
	type##particletype pbuff[NSAVE];\
	if(RANKINGROUP(MYID(simpar),WGroupSize) != 0 ) MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM(simpar),&status);\
	char outfile[190];\
	sprintf(outfile,"PreFoF."#TYPE".%.5d.%.5d",nowstep, NID(simpar));\
	FILE *wp = fopen(outfile,"w");\
	for(i=0;i<TYPE##_NP(simpar);i++){\
		if(IS_P_FLAG(simpar, TYPE, i, FoFflag)){\
			type##particletype *tmp = (type##particletype *) ((char*)(TYPE##_TBP(simpar)+i) + offset);\
			pbuff[isave++] = *tmp;\
			if(isave>= NSAVE){\
				fwrite(pbuff, sizeof(type##particletype), isave, wp);\
				isave = 0;\
			}\
		}\
	}\
	if(isave >0){\
		fwrite(pbuff, sizeof(type##particletype), isave, wp);\
	}\
	fclose(wp);\
	if(GROUPID(MYID(simpar),WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < NID(simpar))\
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);\
} while(0)

#define DumpP2Disk(simpar, TYPE, type, nfiles, nsave, tnsave, headeroffset) do{\
	ptrdiff_t _i, _j;\
	int NSAVE = (NSave+NID(simpar)-1)/NID(simpar);\
	for(_i=0;_i<nfiles;_i++) nsave[_i] = 0;\
	ptrdiff_t offset = sizeof(struct linkedlisttype*);\
	for(_i=0;_i<TYPE##_NP(simpar);_i++) if(IS_P_FLAG(simpar, TYPE, _i, FoFflag)) {\
		type##particletype *tmp = (type##particletype *) ((char*)(TYPE##_TBP(simpar)+_i) + offset);\
		ptrdiff_t ifile = ZofP(simpar, tmp)/zspacing;\
		nsave[ifile]++;\
	}\
	MPI_Gather(nsave, nfiles, MPI_LONG, tnsave, nfiles, MPI_LONG, 0, MPI_COMM(simpar));\
	MPI_Bcast(tnsave, nfiles*nid, MPI_LONG, 0, MPI_COMM(simpar));\
	if(MYID(simpar) == 0){\
		for(_i=0;_i<nfiles;_i++){\
			char outfile[190];\
			sprintf(outfile,"PreFoF."#TYPE".%.5d.%.5d",nowstep, (int)_i);\
			FILE *wp = fopen(outfile,"w");\
			float zmin = _i*zspacing; \
			float zmax = (_i+1)*zspacing;\
			SimParameters tsim = *simpar;\
			MYID(&tsim) = _i;\
			NID(&tsim) = nfiles;\
			int ncnt=0;\
			FILE_HEADER_ZMINZMAX(fprintf,wp, ,(&tsim),zmin,zmax);\
			headeroffset[_i] = ftell(wp);\
			fclose(wp);\
		}\
	}\
	MPI_Bcast(headeroffset, nfiles, MPI_LONG, 0, MPI_COMM(simpar));\
	MPI_Status status;\
	for(_i=0;_i<nfiles;_i++){\
		nsave[_i] = 0;\
		for(_j=0;_j<MYID(simpar);_j++){\
			nsave[_i] += tnsave[_i +nfiles*_j];\
		}\
		headeroffset[_i] += nsave[_i] * sizeof(type##particletype);\
	}\
	if(MYID(simpar) == NID(simpar)-1){\
		int null = 0;\
		for(_i=0;_i<nfiles;_i++){\
			char outfile[190];\
			sprintf(outfile,"PreFoF."#TYPE".%.5d.%.5d",nowstep, (int)_i);\
			FILE *wp = fopen(outfile,"r+");\
			fseek(wp, (headeroffset[_i]-sizeof(int)), SEEK_SET);\
			fwrite(&null, sizeof(int), 1, wp);\
			fclose(wp);\
		}\
	}\
	MPI_Barrier(MPI_COMM(simpar));\
	int itag = 1;\
	int isend, iget, WGroupSize = WGROUPSIZE(simpar);\
	int src = MYID(simpar) - 1;\
	int tgt = MYID(simpar) + 1;\
	type##particletype *tmp;\
	for(_i=0;_i<nfiles;_i++) nsave[_i] = 0;\
	type##particletype pbuff[nfiles][NSAVE];\
	if(RANKINGROUP(MYID(simpar),WGroupSize) != 0 ) MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM(simpar),&status);\
	for(_i=0;_i<TYPE##_NP(simpar);_i++){\
		if(IS_P_FLAG(simpar, TYPE, _i, FoFflag)){\
			type##particletype *tmp = (type##particletype *) ((char*)(TYPE##_TBP(simpar)+_i) + offset);\
			ptrdiff_t ifile = ZofP(simpar, tmp)/zspacing;\
			pbuff[ifile][nsave[ifile]++] = *tmp;\
			if(nsave[ifile]>= NSAVE){\
				char outfile[190];\
				sprintf(outfile,"PreFoF."#TYPE".%.5d.%.5d",nowstep, (int)ifile);\
				FILE *wp = fopen(outfile,"r+");\
				fseek(wp, headeroffset[ifile], SEEK_SET);\
				fwrite(pbuff[ifile], sizeof(type##particletype), nsave[ifile], wp);\
				fclose(wp);\
				headeroffset[ifile] += nsave[ifile]*sizeof(type##particletype);\
				nsave[ifile] = 0;\
			}\
		}\
	}\
	for(_i=0;_i<nfiles;_i++){\
		if(nsave[_i] >0){\
			char outfile[190];\
			sprintf(outfile,"PreFoF."#TYPE".%.5d.%.5d",nowstep, (int)_i);\
			FILE *wp = fopen(outfile,"r+");\
			fseek(wp, headeroffset[_i], SEEK_SET);\
			fwrite(pbuff[_i], sizeof(type##particletype), nsave[_i], wp);\
			fclose(wp);\
			headeroffset[_i] += nsave[_i]*sizeof(type##particletype);\
			nsave[_i] = 0;\
		}\
	}\
	if(GROUPID(MYID(simpar),WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < NID(simpar))\
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);\
} while(0)

void DumpFoFP2Disk(SimParameters *simpar){
	ptrdiff_t i,j,k;
	int nfiles;
	nfiles = MIN(NZ(simpar)/4, NID(simpar));
	nfiles = ceil(log(nfiles)/log(2.L));
	nfiles = (1<<nfiles);
	double zspacing = NZ(simpar)/nfiles;
	int nid = NID(simpar);
	long nsave[nfiles], tnsave[nfiles*nid], headeroffset[nfiles];
	int nowstep = STEPCOUNT(simpar);
	
	DumpP2Disk(simpar, DM,dm,nfiles, nsave, tnsave, headeroffset);
#ifndef GOTPM
	DumpP2Disk(simpar, SPH, sph,nfiles, nsave, tnsave, headeroffset);

	DumpP2Disk(simpar, STAR, star,nfiles, nsave, tnsave, headeroffset);

	DumpP2Disk(simpar , AGN, agn,nfiles, nsave, tnsave, headeroffset);
#endif
}

void CountHaloCandidates(SimParameters *simpar, size_t nhalo, HaloBound *halo,
		FoFTPtlStruct *point, size_t npoint){
	float fof_link = FOF_LINK(simpar);
	PosType xf,yf,zf;
	xf = CellWidth;
	yf = CellWidth;
	zf = CellWidth;
	size_t i, mh;
	FoFTPtlStruct *tmp;
	size_t localnp;
	for(i=0;i<nhalo;i++){
		halo[i].xmin = halo[i].ymin = halo[i].zmin = 1.E20;
		halo[i].xmax = halo[i].ymax = halo[i].zmax = -1.E20;
		halo[i].dumpflag = 0;
		halo[i].sibling = NULL;
		halo[i].nmem = 0;
	}
	for(i=0;i<npoint; i++){
		mh = point[i].haloindx;
		halo[mh].nmem ++;
		tmp = halo[mh].sibling;
		halo[mh].sibling = &(point[i]);
		point[i].sibling = tmp;
		halo[mh].xmin = MIN(halo[mh].xmin, point[i].x);
		halo[mh].ymin = MIN(halo[mh].ymin, point[i].y);
		halo[mh].zmin = MIN(halo[mh].zmin, point[i].z);
		halo[mh].xmax = MAX(halo[mh].xmax, point[i].x);
		halo[mh].ymax = MAX(halo[mh].ymax, point[i].y);
		halo[mh].zmax = MAX(halo[mh].zmax, point[i].z);
	}
	localnp = 0;
	for(i=0;i<nhalo;i++){
		if(halo[i].xmin >fof_link && halo[i].ymin > fof_link && halo[i].zmin >fof_link &&
				halo[i].xmax < xf-fof_link && halo[i].ymax < yf-fof_link && halo[i].zmax <zf-fof_link){
			if(halo[i].nmem >= FOF_MINNUM(simpar)){
				halo[i].dumpflag = 1;
				localnp += halo[i].nmem;
			}
		}
		else {
			halo[i].dumpflag = 1;
			localnp += halo[i].nmem;
		}
	}
	for(i=0;i<nhalo;i++){
		if(halo[i].dumpflag){
			tmp=halo[i].sibling;
			while(tmp){
				linkedlisttype *p;
				p = (linkedlisttype*)(tmp->bp);
				SET_FLAG(p, FoFflag);
				tmp = tmp->sibling;
			}
		}
	}
	/*
	return localnp;
	*/
}

void PreFoF(SimParameters *simpar){
	ptrdiff_t nx = NX(simpar);
	ptrdiff_t ny = NY(simpar);
	ptrdiff_t nz = NZ(simpar);

	PosType xmin,ymin,zmin,xmax,ymax,zmax;


	FOF_MINNUM(simpar) = 20;


	/*
	xmin = COS_XMIN(simpar);
	ymin = COS_YMIN(simpar);
	zmin = COS_ZMIN(simpar);
	xmax = COS_XMAX(simpar);
	ymax = COS_YMAX(simpar);
	zmax = COS_ZMAX(simpar);
	*/
	xmin = SIM_LXMIN(simpar,dm);
	ymin = SIM_LYMIN(simpar,dm);
	zmin = SIM_LZMIN(simpar,dm);
	xmax = SIM_LXMAX(simpar,dm);
	ymax = SIM_LYMAX(simpar,dm);
	zmax = SIM_LZMAX(simpar,dm);


#ifdef DEBUG
	DEBUGPRINT("P%d has %g %g : %g %g: %g %g\n",MYID(simpar), xmin,xmax, ymin, ymax, zmin,zmax);
	DEBUGPRINT("P%d has %ld %ld %ld %ld\n",MYID(simpar), DM_NP(simpar), SPH_NP(simpar),STAR_NP(simpar), AGN_NP(simpar));
#endif

	RESET_FOFFLAG(simpar);



	float fof_link = FOF_LINK(simpar);
	ptrdiff_t mx,my,mz;

	mx = ceil((xmax-xmin)/CellWidth);
	my = ceil((ymax-ymin)/CellWidth);
	mz = ceil((zmax-zmin)/CellWidth);
	if(MYID(simpar) == -1)
	{
		int kkk = 1;
		while(kkk) {
			kkk = 1;
		}
	}


	BASICCELL(simpar) = (TreeLinkedCell*)malloc(sizeof(TreeLinkedCell)*mx*my*mz);

	/*
	BuildLinkedList(simpar,mx,my,mz, xmin,ymin,zmin);
	*/
	{ 
		TreeLinkedCell *BasicCell = BASICCELL(simpar);
		ptrdiff_t _i;
		for(_i=0;_i<mx*my*mz;_i++) {
			BasicCell[_i].link = NULL;
			BasicCell[_i].nmem = 0;
		}
		treedmparticletype *bp = DM_TBP(simpar);
		for(_i=0;_i<DM_NP(simpar);_i++){
			ptrdiff_t ix,iy,iz;
			ix = (XofP(simpar,bp+_i)-xmin)/CellWidth;
			iy = (YofP(simpar,bp+_i)-ymin)/CellWidth;
			iz = (ZofP(simpar,bp+_i)-zmin)/CellWidth;
			ptrdiff_t ipos = ix+mx*(iy+my*iz);
			linkedlisttype *tmp = BasicCell[ipos].link;
			BasicCell[ipos].link = (linkedlisttype*)(bp+_i);
			BasicCell[ipos].nmem ++;
			(bp+_i)->next = tmp;
		}
	}



#ifdef DEBUG
	DEBUGPRINT("P%d is now at PreFoF   %ld %ld %ld\n",MYID(simpar), mx,my,mz);
#endif

	size_t ii;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 10)
#endif
	for(ii=0;ii<mx*my*mz;ii++){
		size_t iz = ii/(mx*my); 
		size_t iy = (ii % (mx*my)) / mx;
		size_t ix = (ii % mx);
		TreeLinkedCell *BasicCell = BASICCELL(simpar);
		FoFTPtlStruct *point;
		ptrdiff_t npoint = Dump2FoFPtl(simpar, ix,iy,iz, mx,my,mz, &point,
				xmin,ymin,zmin);
		if(npoint > 0){
			FoFPosition *linked = (FoFPosition*)malloc(sizeof(FoFPosition)*npoint);
			size_t i,nnode = npoint;
			FoFTStruct *TREECELL = (FoFTStruct*)malloc(sizeof(FoFTStruct)*nnode);
			FoF_Make_Tree(TREECELL, nnode, point, npoint, RECURSIVE);
			HaloBound *halo = (HaloBound*)malloc(sizeof(HaloBound)*npoint);
			size_t nhalo = 0;
			for(i=0;i<npoint;i++){
				if(point[i].included == NO){
					size_t num;
					FoFPosition p;
					p.x = point[i].x;
					p.y = point[i].y;
					p.z = point[i].z;
					num = pnew_fof_link(&p, TREECELL, fof_link, linked, nhalo);
					nhalo ++;
				}
			}
			free(TREECELL);
			CountHaloCandidates(simpar, nhalo, halo, point, npoint);
			free(halo);
			free(linked);
			free(point);
		}
	}
#ifdef DEBUG
	DEBUGPRINT("P%d passed with \n",MYID(simpar));
#endif
	DumpFoFP2Disk(simpar);
}

