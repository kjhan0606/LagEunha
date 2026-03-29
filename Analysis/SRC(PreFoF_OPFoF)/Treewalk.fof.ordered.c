/* This makes Tree structure and walks along tree structures.
 *
 * */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "eunha.h"
//#include "pmheader.h"
#include "mpi.h"
#include "fof.h"
#define IMOD(A,B) ((A) - ((A)/(B))*(B))
#define MIN(A,B) ((A)<(B) ? (A):(B))
#define MAX(A,B) ((A)>(B) ? (A):(B))
POSTYPE Lx2, Ly2,Lz2;
POSTYPE Lx, Ly,Lz;
/* open node in periodic boundary conditions */
enum boolean pfof_open(particle p,FoFTStruct *tree, POSTYPE fof_link){
	POSTYPE tmpx,tmpy,tmpz; 
	POSTYPE dist2,dist,r,diff; 
	POSTYPE ratio; 
	tmpx = fabs(p.x-tree->mono[0]); 
	tmpy = fabs(p.y-tree->mono[1]); 
	tmpz = fabs(p.z-tree->mono[2]); 
	if(tmpx > Lx2) tmpx = Lx-tmpx;
	if(tmpy > Ly2) tmpy = Ly-tmpy;
	if(tmpz > Lz2) tmpz = Lz-tmpz;
	r = tree->dist;
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz; 
	dist = sqrt(dist2);
	diff = dist - r;
	if(diff <= fof_link) return YES;
	else return NO;
}
enum boolean fof_open(particle p,FoFTStruct *tree, POSTYPE fof_link){
	POSTYPE tmpx,tmpy,tmpz; 
	POSTYPE dist2,dist,r,diff; 
	POSTYPE ratio; 
	tmpx = p.x-tree->mono[0]; 
	tmpy = p.y-tree->mono[1]; 
	tmpz = p.z-tree->mono[2]; 
	r = tree->dist;
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz; 
	dist = sqrt(dist2);
	diff = dist - r;
	if(diff <= fof_link) return YES;
	else return NO;
}
static int ntmp;   
static POSTYPE tmpx,tmpy,tmpz,dist2;   
static POSTYPE xx,yy,zz,xy,xz,yz,tmpxx; 
static POSTYPE qxx1,qxx2,qxx3,qxy,tmpdist;   
static POSTYPE fplmf1,fplmf2,fplmf3,fplmf4,tmptmp;   
static POSTYPE fplmf;  
/*
static TStruct *pointer;   
static TPtlStruct *ppointer;   
*/
static POSTYPE ptlmass;
static POSTYPE idist2,sqrtdist2,isqrtdist2; 

/* 
 * *p is the position at which you want to calculate force using tree
 * theta2 is the opening angle for tree walk
 * *tree is the tree structure.
 */
/* cross? is needed to check whether the nearby particle searching crosses 
 * the boundary face. */
enum boolean crossx,crossy,crossz;
#ifdef OLD
void WriteIsolatedHalo(size_t nhalo, HaloBound *halobound,FoFTPtlStruct *ptl,
		particle *linked, char *halofile , char *memparticlefile){
	int i;
	size_t j;
	size_t np;
	HaloQ haloq;
	FILE *hfp,*pfp;
	FoFTPtlStruct *ptr;
	int myid,nid;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	fflush(stdout);
	int mid;
#ifdef USE_MASTER
	mid = nid -1;
#else
	mid = nid;
#endif
	for(i=0;i<mid;i++){
		if(i==myid){
			printf("P%d is preparing to write isolated %ld halo data ",myid,nhalo);fflush(stdout);
			hfp=fopen(halofile,"a");
			pfp=fopen(memparticlefile,"a");
			for(j=0;j<nhalo;j++){
				if(halobound[j].boundflag ==0){
					ptr = halobound[j].sibling;
					np = 0;
					while(ptr){
						/*
						linked[np].x = xofP(ptr);
						linked[np].y = yofP(ptr);
						linked[np].z = zofP(ptr);
						*/
						linked[np].x = ptr->r[0];
						linked[np].y = ptr->r[1];
						linked[np].z = ptr->r[2];
						linked[np].vx = ptr->rv[0];
						linked[np].vy = ptr->rv[1];
						linked[np].vz = ptr->rv[2];
						linked[np].indx = ptr->indx;
						ptr = ptr->sibling;
						np++;
					}
					if(np >= MinNumMem){
						haloq=haloproperty(linked,np);
						fwrite(&haloq,sizeof(HaloQ),1,hfp);
						fwrite(linked,sizeof(particle),np,pfp);
					}
				}
			}
			printf("P%d have written isolated halo data\n",myid);fflush(stdout);
			fclose(hfp);fclose(pfp);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}
#elif USE_MASTER
void WriteIsolatedHalo(size_t nhalo, HaloBound *halobound,FoFTPtlStruct *ptl,
		particle *linked, char *halofile , char *memparticlefile){
	int i;
	size_t j;
	size_t np;
	HaloQ haloq;
	FILE *hfp,*pfp;
	FoFTPtlStruct *ptr;
	int rktag=99,sktag=987;
	int myid,nid;
	int masterid;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	masterid = nid -1;
	fflush(stdout);
	if(myid==masterid){
		printf("P%d is preparing to write isolated %ld halo data ",myid,nhalo);fflush(stdout);
		hfp=fopen(halofile,"a");
		pfp=fopen(memparticlefile,"a");
		int ncount = 1;
		while(ncount != nid){
			int ok = 1;
			int src;
			MPI_Status rstatus;
			MPI_Probe(MPI_ANY_SOURCE, sktag,MPI_COMM_WORLD,&rstatus);
			src = rstatus.MPI_SOURCE;
			MPI_Recv(&ok,1,MPI_INT,src,sktag,MPI_COMM_WORLD,&rstatus);
			MPI_Recv(&haloq,sizeof(HaloQ),MPI_BYTE,src,rktag,MPI_COMM_WORLD,&rstatus);
			while(haloq.np !=0){
				MPI_Recv(linked,haloq.np*sizeof(particle),MPI_BYTE,src,rktag,MPI_COMM_WORLD,&rstatus);
				fwrite(&haloq,sizeof(HaloQ),1,hfp);
				fwrite(linked,sizeof(particle),haloq.np,pfp);
				MPI_Recv(&haloq,sizeof(HaloQ),MPI_BYTE,src,rktag,MPI_COMM_WORLD,&rstatus);
			}
			ncount ++;
		}
		fclose(hfp);
		fclose(pfp);
		printf("P%d have written isolated halo data\n",myid);fflush(stdout);
	}
	else {
		int ok = 1;
		MPI_Send(&ok,1,MPI_INT,0,sktag,MPI_COMM_WORLD);
		for(j=0;j<nhalo;j++){
			if(halobound[j].boundflag ==0){
				ptr = halobound[j].sibling;
				np = 0;
				while(ptr){
					linked[np].x = xofP(ptr);
					linked[np].y = yofP(ptr);
					linked[np].z = zofP(ptr);
					linked[np].vx = ptr->rv[0];
					linked[np].vy = ptr->rv[1];
					linked[np].vz = ptr->rv[2];
					linked[np].indx = ptr->indx;
					ptr = ptr->sibling;
					np++;
				}
				if(np >= MinNumMem){
					haloq=haloproperty(linked,np);
					MPI_Send(&haloq,sizeof(HaloQ),MPI_BYTE,0,rktag,MPI_COMM_WORLD);
					MPI_Send(linked,haloq.np*sizeof(particle),MPI_BYTE,0,rktag,MPI_COMM_WORLD);
				}
			}
		}
		haloq.np = 0;
		MPI_Send(&haloq,sizeof(HaloQ),MPI_BYTE,0,rktag,MPI_COMM_WORLD);
		printf("P%d have sent all isolated halo data to P0\n",myid);fflush(stdout);
	}
}
#else 
void WriteIsolatedHalo(size_t nhalo, HaloBound *halobound,FoFTPtlStruct *ptl,
		particle *linked, char *halofile , char *memparticlefile){
	int i;
	size_t j;
	size_t np;
	HaloQ haloq;
	FILE *hfp,*pfp;
	FoFTPtlStruct *ptr;
	int rktag=99,sktag=987;
	int myid,nid;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	fflush(stdout);
	if(myid==0){
		printf("P%d is preparing to write isolated %ld halo data ",myid,nhalo);fflush(stdout);
		hfp=fopen(halofile,"a");
		pfp=fopen(memparticlefile,"a");
		{
			for(j=0;j<nhalo;j++){
				if(halobound[j].boundflag ==0){
					ptr = halobound[j].sibling;
					np = 0;
					while(ptr){
						linked[np].x = xofP(ptr);
						linked[np].y = yofP(ptr);
						linked[np].z = zofP(ptr);
						linked[np].vx = ptr->rv[0];
						linked[np].vy = ptr->rv[1];
						linked[np].vz = ptr->rv[2];
						linked[np].indx = ptr->indx;
						ptr = ptr->sibling;
						np++;
					}
					if(np >= MinNumMem){
						haloq=haloproperty(linked,np);
						fwrite(&haloq,sizeof(HaloQ),1,hfp);
						fwrite(linked,sizeof(particle),np,pfp);
					}
				}
			}
		}

		int ncount = 1;
		while(ncount != nid){
			int ok = 1;
			int src;
			MPI_Status rstatus;
			MPI_Probe(MPI_ANY_SOURCE, sktag,MPI_COMM_WORLD,&rstatus);
			src = rstatus.MPI_SOURCE;
			MPI_Recv(&ok,1,MPI_INT,src,sktag,MPI_COMM_WORLD,&rstatus);
			MPI_Recv(&haloq,sizeof(HaloQ),MPI_BYTE,src,rktag,MPI_COMM_WORLD,&rstatus);
			while(haloq.np !=0){
				MPI_Recv(linked,haloq.np*sizeof(particle),MPI_BYTE,src,rktag,MPI_COMM_WORLD,&rstatus);
				fwrite(&haloq,sizeof(HaloQ),1,hfp);
				fwrite(linked,sizeof(particle),haloq.np,pfp);
				MPI_Recv(&haloq,sizeof(HaloQ),MPI_BYTE,src,rktag,MPI_COMM_WORLD,&rstatus);
			}
			ncount ++;
		}
		fclose(hfp);
		fclose(pfp);
		printf("P%d have written isolated halo data\n",myid);fflush(stdout);
	}
	else {
		int ok = 1;
		MPI_Send(&ok,1,MPI_INT,0,sktag,MPI_COMM_WORLD);
		for(j=0;j<nhalo;j++){
			if(halobound[j].boundflag ==0){
				ptr = halobound[j].sibling;
				np = 0;
				while(ptr){
					linked[np].x = xofP(ptr);
					linked[np].y = yofP(ptr);
					linked[np].z = zofP(ptr);
					linked[np].vx = ptr->rv[0];
					linked[np].vy = ptr->rv[1];
					linked[np].vz = ptr->rv[2];
					linked[np].indx = ptr->indx;
					ptr = ptr->sibling;
					np++;
				}
				if(np >= MinNumMem){
					haloq=haloproperty(linked,np);
					MPI_Send(&haloq,sizeof(HaloQ),MPI_BYTE,0,rktag,MPI_COMM_WORLD);
					MPI_Send(linked,haloq.np*sizeof(particle),MPI_BYTE,0,rktag,MPI_COMM_WORLD);
				}
			}
		}
		haloq.np = 0;
		MPI_Send(&haloq,sizeof(HaloQ),MPI_BYTE,0,rktag,MPI_COMM_WORLD);
		printf("P%d have sent all isolated halo data to P0\n",myid);fflush(stdout);
	}
}

#endif
void WriteFinalHalo(size_t nhalo, HaloBound *halobound,FoFTPtlStruct *ptl,
		particle *linked, char *halofile , char *memparticlefile){
	int i;
	size_t j;
	size_t np;
	HaloQ haloq;
	FILE *hfp,*pfp;
	int myid,nid,mid;
	FoFTPtlStruct *ptr;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
#ifdef USE_MASTER
	mid = nid -1;
#else
	mid = nid;
#endif
	for(i=0;i<mid;i++){
		if(i==myid){
			hfp=fopen(halofile,"a");
			pfp=fopen(memparticlefile,"a");
			for(j=0;j<nhalo;j++){
				ptr = halobound[j].sibling;
				np = 0;
				while(ptr){
					/*
					linked[np].x = xofP(ptr);
					linked[np].y = yofP(ptr);
					linked[np].z = zofP(ptr);
					*/
					linked[np].x = ptr->r[0];
					linked[np].y = ptr->r[1];
					linked[np].z = ptr->r[2];
					linked[np].vx = ptr->rv[0];
					linked[np].vy = ptr->rv[1];
					linked[np].vz = ptr->rv[2];
					linked[np].indx = ptr->indx;
					ptr = ptr->sibling;
					np++;
				}
				if(np >= MinNumMem){
					haloq=haloproperty(linked,np);
					fwrite(&haloq,sizeof(HaloQ),1,hfp);
					fwrite(linked,sizeof(particle),np,pfp);
				}
			}
			printf("written\n");fflush(stdout);
			fclose(hfp);fclose(pfp);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

void ReadBottomFaceContact(FoFTPtlStruct *p,size_t npread, particle *linked,int src,
		int nstep, int nz){
	FoFTPtlStruct *tmp;
	size_t i;
	FILE *fp;
	char infile[190];
	size_t nread;
	int myid,nid;
	int j;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);

	int mid;
#ifdef USE_MASTER
	mid = nid-1;
#else
	mid = nid;
#endif

	sprintf(infile,"BottomFaceContactHalo%.5d%.5d.dat",src,nstep);
	if((fp=fopen(infile,"r")) == NULL)
#ifdef PREVIOUS
	{
		fprintf(stderr,"error opening %s\n",infile);
		exit(0);
	}
#else
	{ }
	else {
#endif
#ifdef OLD
	for(j=0;j<mid;j++){
		if(myid==j)
#endif
		{
			POSTYPE zoffset;
			if(myid==mid-1) zoffset = nz;
			else zoffset = 0;
			tmp = p;
			while((nread=fread(linked,sizeof(particle),MaxLinkedParticles,fp)) > 0 ){
				for(i=0;i<nread;i++){
					tmp->r[0] = linked[i].x;
					tmp->r[1] = linked[i].y;
					tmp->r[2] = linked[i].z + zoffset;
					tmp->rv[0] = linked[i].vx;
					tmp->rv[1] = linked[i].vy;
					tmp->rv[2] = linked[i].vz;
					tmp->indx = linked[i].indx;
					tmp++;
				}
			}
			printf("Now read %ld particles from %ld particles\n",tmp-p,npread);
		}
#ifdef OLD
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif
	fclose(fp);
#ifndef PREVIOUS
	}
#endif
}

size_t WriteBottomFaceContact(size_t nhalo, HaloBound *halobound,
		FoFTPtlStruct *ptl, particle *linked,int nowfile,int nstep){
	int i;
	size_t j;
	size_t np;
	HaloQ haloq;
	FILE *hfp,*pfp;
	char bottomboxhalo[80];
	size_t npwrite;
	FoFTPtlStruct *ptr;
	int myid,nid;
	
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);

#ifdef OLD
	int mid;
#ifdef USE_MASTER
	mid = nid -1;
#else
	mid = nid;
#endif
	for(i=0;i<mid;i++){
		if(i==myid)
#endif
		{
			sprintf(bottomboxhalo,"BottomFaceContactHalo%.5d%.5d.dat",myid,nstep);
			if(nowfile == 0) pfp=fopen(bottomboxhalo,"w");
			else  pfp=fopen(bottomboxhalo,"a");
			npwrite = 0;
			for(j=0;j<nhalo;j++){
				if(halobound[j].boundflag ==1){
					ptr = halobound[j].sibling;
					np = 0;
					while(ptr){
						linked[np].x = ptr->r[0];
						linked[np].y = ptr->r[1];
						linked[np].z = ptr->r[2];
						linked[np].vx = ptr->rv[0];
						linked[np].vy = ptr->rv[1];
						linked[np].vz = ptr->rv[2];
						linked[np].indx = ptr->indx;
						ptr = ptr->sibling;
						np++;
					}
					fwrite(linked,sizeof(particle),np,pfp);
					npwrite += np;
				}
			}
			fclose(pfp);
			printf("writing total %ld bottom FoF particles\n",npwrite);
		}
#ifdef OLD
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif
	return (npwrite);
}

void WriteAllHalo(size_t nhalo, HaloBound *halobound,FoFTPtlStruct *ptl,
		size_t np, particle *linked, char *halofile , char *memparticlefile){
	int i;
	size_t j;
	size_t mh,mp;
	HaloQ haloq;
	FILE *hfp,*pfp;
	FoFTPtlStruct *ptr;
	for(j=0;j<nhalo;j++){
		halobound[j].sibling = NULL;
	}
	for(j=0;j<np;j++){
		mh = ptl[j].haloindx;
		ptr = halobound[mh].sibling;
		halobound[mh].sibling = &(ptl[j]);
		ptl[j].sibling = ptr;
	}
			hfp=fopen(halofile,"a");
			pfp=fopen(memparticlefile,"a");
			for(j=0;j<nhalo;j++){
				ptr = halobound[j].sibling;
				mp = 0;
				while(ptr){
					/*
					linked[mp].x = xofP(ptr);
					linked[mp].y = yofP(ptr);
					linked[mp].z = zofP(ptr);
					*/
					linked[mp].x = ptr->r[0];
					linked[mp].y = ptr->r[1];
					linked[mp].z = ptr->r[2];
					linked[mp].vx = ptr->rv[0];
					linked[mp].vy = ptr->rv[1];
					linked[mp].vz = ptr->rv[2];
					linked[mp].indx = ptr->indx;
					ptr = ptr->sibling;
					mp++;
				}
				if(mp >= MinNumMem){
					haloq=haloproperty(linked,mp);
					fwrite(&haloq,sizeof(HaloQ),1,hfp);
					fwrite(linked,sizeof(particle),mp,pfp);
				}
			}
			fclose(hfp);fclose(pfp);
}
size_t StackUpContactParticleLeftWard(size_t nhalo,HaloBound *halobound,
		FoFTPtlStruct *ptl,size_t np){
	size_t nowp,i;
	FoFTPtlStruct *ptr;
	size_t mh;
	nowp = 0;
	for(i=0;i<np;i++){
		if(halobound[ptl[i].haloindx].boundflag ==2 ||
				halobound[ptl[i].haloindx].boundflag ==3){
			ptl[nowp] = ptl[i];
			nowp ++;
		}
	}
	return (nowp);
}

void CheckHaloBound(size_t nhalo,HaloBound *halo,FoFTPtlStruct *ptl,size_t np,
		POSTYPE fof_link,POSTYPE zdown,POSTYPE zup, POSTYPE zminlocal){
	size_t i,mh;
	FoFTPtlStruct *tmp;
	size_t maxnmem;

	for(i=0;i<nhalo;i++){
		halo[i].zmin = 2.E23;
		halo[i].zmax = -2.E23;
		halo[i].boundflag = 0;
		halo[i].sibling = NULL;
		halo[i].nmem = 0;
	}
	for(i=0;i<np;i++){
		mh = ptl[i].haloindx;
		halo[mh].nmem++;
		tmp = halo[mh].sibling;
		halo[mh].sibling = &(ptl[i]);
		ptl[i].sibling = tmp;
		halo[mh].zmin = MIN(halo[mh].zmin,ptl[i].r[2]);
		halo[mh].zmax = MAX(halo[mh].zmax,ptl[i].r[2]);
	}
	maxnmem = 0;
	for(i=0;i<nhalo;i++){
		if(halo[i].zmin <= zminlocal+fof_link) halo[i].boundflag ++;
		if(halo[i].zmax >= zup-fof_link)  halo[i].boundflag += 2;
		maxnmem = MAX(maxnmem,halo[i].nmem);
	}
	printf("Maximum FoF member of particles %ld\n",maxnmem);
	{
		size_t num1,num2,num3,num0;
		num1=num2=num3=num0=0;
		for(i=0;i<nhalo;i++){
			if(halo[i].boundflag == 0) num0 ++;
			else if(halo[i].boundflag ==1) num1 ++;
			else if(halo[i].boundflag ==2) num2 ++;
			else if(halo[i].boundflag ==3) num3 ++;
		}
		printf("num0=%ld num1=%ld num2=%ld num3=%ld\n",num0,num1,num2,num3);
	}
}
/* lx,ly, lz are the limits of the box */
size_t pnew_fof_link(particle *p,POSTYPE fof_link,FoFTStruct *tree,
		FoFTPtlStruct *ptl,particle *linked,size_t nhalo,
		POSTYPE lx,POSTYPE ly, POSTYPE lz){
	size_t ncount, now;
	void *ptr,*optr,*nptr;
	POSTYPE tmpx,tmpy,tmpz;
	POSTYPE fof_link2,dist2;
	POSTYPE Lpx,Lpy,Lpz;
	particle point;
	Lx = lx; Ly = ly; Lz = lz;
	Lx2 = Lx*0.5;
	Ly2 = Ly*0.5;
	Lz2 = Lz*0.5;
	fof_link2 = fof_link*fof_link;
	ncount = now = 0;
	point.x = p->x;
	point.y = p->y;
	point.z = p->z;
	do {
		optr = (void *) tree;
		ptr = (void*) tree;
		while(ptr != NULL){
			switch(((TYPE*)ptr)->type){
				case TYPE_TREE:
					if(((FoFTStruct *)ptr)->sibling ==
							((FoFTStruct *)ptr)->daughter){
						EraseFromTree(optr,ptr,((FoFTStruct *)ptr)->sibling);
						ptr = ((FoFTStruct *)ptr)->sibling;
					}
					else
					switch(pfof_open(point,ptr,fof_link)){
						case YES:
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->daughter);
							break;
						default :
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->sibling);
					}
					break;
				default :
					if(((FoFTPtlStruct*)ptr)->included == YES){
						nptr = ((FoFTPtlStruct *)ptr)->sibling;
						EraseFromTree(optr,ptr,nptr);
						ptr = nptr;
					}
					else {
						tmpx = fabs(point.x - ((FoFTPtlStruct*)ptr)->r[0]);
						if(tmpx>Lx2) tmpx = Lx-tmpx;

						tmpy = fabs(point.y - ((FoFTPtlStruct*)ptr)->r[1]);
						if(tmpy>Ly2) tmpy = Ly-tmpy;

						tmpz = fabs(point.z - ((FoFTPtlStruct*)ptr)->r[2]);
						if(tmpz>Lz2) tmpz = Lz-tmpz;

						dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
						dist2 = sqrt(dist2);
						if(dist2 <= fof_link){
							linked[ncount].x = ((FoFTPtlStruct*)ptr)->r[0];
							linked[ncount].y = ((FoFTPtlStruct*)ptr)->r[1];
							linked[ncount].z = ((FoFTPtlStruct*)ptr)->r[2];
							linked[ncount].vx = ((FoFTPtlStruct*)ptr)->rv[0];
							linked[ncount].vy = ((FoFTPtlStruct*)ptr)->rv[1];
							linked[ncount].vz = ((FoFTPtlStruct*)ptr)->rv[2];
							linked[ncount].indx = ((FoFTPtlStruct*)ptr)->indx;
							((FoFTPtlStruct*)ptr)->haloindx = nhalo;
							((FoFTPtlStruct*)ptr)->included = YES;
							ncount ++;
							nptr = ((FoFTPtlStruct *)ptr)->sibling;
							EraseFromTree(optr,ptr,nptr);
						}
						else optr = ptr;
						ptr = (void*)(((FoFTPtlStruct*)ptr)->sibling);
					}
			}
		}
		point = linked[now];
		/* for periodic boundary conditions in x,y, and z*/
		now ++;
	} while( now <= ncount);
	return (ncount);
}
int new_fof_link(particle *p,POSTYPE fof_link,FoFTStruct *tree,
		FoFTPtlStruct *ptl,particle *linked){
	int ncount, now;
	void *ptr,*optr,*nptr;
	POSTYPE tmpx,tmpy,tmpz;
	POSTYPE fof_link2,dist2;
	particle point;
	fof_link2 = fof_link*fof_link;
	ncount = now = 0;
	point.x = p->x;
	point.y = p->y;
	point.z = p->z;
	do {
		optr = (void *) tree;
		ptr = (void*) tree;
		while(ptr != NULL){
			switch(((TYPE*)ptr)->type){
				case TYPE_TREE:
					if(((FoFTStruct *)ptr)->sibling ==
							((FoFTStruct *)ptr)->daughter){
						EraseFromTree(optr,ptr,((FoFTStruct *)ptr)->sibling);
						ptr = ((FoFTStruct *)ptr)->sibling;
					}
					else
					switch(fof_open(point,ptr,fof_link)){
						case YES:
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->daughter);
							break;
						default :
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->sibling);
					}
					break;
				default :
					if(((FoFTPtlStruct*)ptr)->included == YES){
						nptr = ((FoFTPtlStruct *)ptr)->sibling;
						EraseFromTree(optr,ptr,nptr);
						ptr = nptr;
					}
					else {
						tmpx = point.x - ((FoFTPtlStruct*)ptr)->r[0];
						tmpy = point.y - ((FoFTPtlStruct*)ptr)->r[1];
						tmpz = point.z - ((FoFTPtlStruct*)ptr)->r[2];
						dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
						dist2 = sqrt(dist2);
						if(dist2 <= fof_link){
							linked[ncount].x = ((FoFTPtlStruct*)ptr)->r[0];
							linked[ncount].y = ((FoFTPtlStruct*)ptr)->r[1];
							linked[ncount].z = ((FoFTPtlStruct*)ptr)->r[2];
							linked[ncount].vx = ((FoFTPtlStruct*)ptr)->rv[0];
							linked[ncount].vy = ((FoFTPtlStruct*)ptr)->rv[1];
							linked[ncount].vz = ((FoFTPtlStruct*)ptr)->rv[2];
							linked[ncount].indx = ((FoFTPtlStruct*)ptr)->indx;
							((FoFTPtlStruct*)ptr)->included = YES;
							ncount ++;
							nptr = ((FoFTPtlStruct *)ptr)->sibling;
							EraseFromTree(optr,ptr,nptr);
						}
						else optr = ptr;
						ptr = (void*)(((FoFTPtlStruct*)ptr)->sibling);
					}
			}
		}
		point = linked[now];
		now ++;
	} while( now <= ncount);
	return (ncount);
}
void FoF_Make_Tree(FoFTStruct *TREE_START,FoFTPtlStruct *ptl,size_t np,Box box){
	FoFBeginEndTree beginend;
	FoFTStruct *NewTree;
	FoFTPtlStruct *ptr;
	ptr = ptl;
	while(ptr != NULL){
		ptr->included = NO;
		ptr = ptr->sibling;
	}
	TREE_START->sibling = NULL;
	beginend = FoF_divide_node(TREE_START,TREE_START+1,ptl,box,TREE_START);
}
FoFBeginEndTree FoF_divide_node(FoFTStruct *TREE_START,FoFTStruct *NewTree, 
		FoFTPtlStruct *ptl, Box box,FoFTStruct *ThisTree){ 
	FoFBeginEndTree beginend;
	FoFTStruct *p2tree,tmpnode[8];
	FoFTStruct *NowCal;
	void *from_sibling;
	FoFTPtlStruct *p2ptl,*tmpptr,*tmpptr2;
	Box tmpbox[8];
	int i,j,k,mnode,mx,my,mz;
	POSTYPE x0,y0,z0,inv_halfw,halfw;
	POSTYPE tmpx,tmpy,tmpz,tmpdist2,distmax;
	POSTYPE ptlmass;
	int count;
	ThisTree->type = TYPE_TREE;
	/*
	ThisTree->sibling = NULL;
	*/
	ThisTree->daughter = NULL;
	ThisTree->L = box.width;
	ThisTree->r0[0] = box.x;
	ThisTree->r0[1] = box.y;
	ThisTree->r0[2] = box.z;
	ThisTree->Nparticle = 0;
	ThisTree->mono[0]= ThisTree->mono[1]= ThisTree->mono[2]= 0.;

	p2ptl = ptl;
	while(p2ptl != NULL){
		ThisTree->Nparticle ++;
		ThisTree->mono[0] += p2ptl->r[0];
		ThisTree->mono[1] += p2ptl->r[1];
		ThisTree->mono[2] += p2ptl->r[2];
		p2ptl = p2ptl->sibling;
	}
	ThisTree->mono[0] /= ThisTree->Nparticle;
	ThisTree->mono[1] /= ThisTree->Nparticle;
	ThisTree->mono[2] /= ThisTree->Nparticle;
	distmax = -1.E20;
	p2ptl = ptl;
	while(p2ptl != NULL){
		tmpx = p2ptl->r[0] - ThisTree->mono[0];
		tmpy = p2ptl->r[1] - ThisTree->mono[1];
		tmpz = p2ptl->r[2] - ThisTree->mono[2];
		tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		distmax = MAX(distmax,tmpdist2);
		p2ptl = p2ptl->sibling;
	}
	ThisTree->dist2 = distmax;
	ThisTree->dist = sqrt(distmax);
	x0 = box.x;
	y0 = box.y;
	z0 = box.z;
	halfw = box.width*0.5;
	inv_halfw = 1./halfw;
	/* initialize temporary tree array */
	for(i=0;i<8;i++) {
		tmpnode[i].sibling = tmpnode[i].daughter = NULL;
		tmpnode[i].Nparticle = 0;
		tmpbox[i].width = halfw;
		/*
		tmpbox[i].x = x0+IMOD(i,2)*halfw;
		tmpbox[i].y = y0+(IMOD(i,4)/2)*halfw;
		*/
		tmpbox[i].x = x0+(i%2)*halfw;
		tmpbox[i].y = y0+((i%4)/2)*halfw;
		tmpbox[i].z = z0+(i/4)*halfw;
	}
	p2ptl = ptl;
	while(p2ptl != NULL){
		mx = (int)((p2ptl->r[0] - x0)*inv_halfw);
		my = (int)((p2ptl->r[1] - y0)*inv_halfw);
		mz = (int)((p2ptl->r[2] - z0)*inv_halfw);
		mx = MIN(mx,1);
		my = MIN(my,1);
		mz = MIN(mz,1);
		mx = MAX(mx,0);
		my = MAX(my,0);
		mz = MAX(mz,0);
		mnode = mx + 2*my + 4*mz;
		mnode = MIN(mnode,7);
		mnode = MAX(mnode,0);
		tmpnode[mnode].Nparticle ++; 
		tmpptr = tmpnode[mnode].daughter;
		tmpptr2 = p2ptl->sibling;
		tmpnode[mnode].daughter = p2ptl;
		p2ptl->sibling = tmpptr;
		p2ptl = tmpptr2;
	}
	/* Making link from Mother Node */
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle > 0) break;
	}
	if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
		ThisTree->daughter = (void *)(NewTree);
	}
	else {
		ThisTree->daughter = (void *) tmpnode[i].daughter ;
	}
	count = 0;
	NowCal = NewTree;
	from_sibling = NULL;
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
			NewTree->daughter = tmpnode[i].daughter;
			if(from_sibling != NULL) 
				((GENERAL_TPtl_POINTER*)from_sibling)->sibling = NewTree;
			from_sibling = NewTree;
			NewTree++;
			count ++;
		}
		else if(tmpnode[i].Nparticle > 0 ){
			tmpptr = tmpnode[i].daughter;
			if(from_sibling != NULL) 
				((GENERAL_TPtl_POINTER*)from_sibling)->sibling = tmpptr;
			while(tmpptr != NULL){
				from_sibling = tmpptr;
				tmpptr = tmpptr->sibling;
			}
		}
	}
	((GENERAL_TPtl_POINTER*)from_sibling)->sibling = ThisTree->sibling;
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
			beginend = FoF_divide_node(TREE_START,NewTree,
					tmpnode[i].daughter, tmpbox[i],NowCal);
			NewTree = beginend.start;
			NowCal ++;
		}
	}
	beginend.start = NewTree;
	return beginend;
}
