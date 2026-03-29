#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<mpi.h>
#include "fof.h"
#define max(A,B) ( (A) > (B) ? (A) : (B) )
#define min(A,B) ( (A) < (B) ? (A) : (B) )
#define TICKS 100.
#define GETCOMM(mp,mmp) {\
	MPI_Reduce(&mp,&mmp,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);\
		MPI_Bcast(&mmp,1,MPI_INT,0,MPI_COMM_WORLD);\
}
void migrate(size_t *mp,FoFTPtlStruct **Bp, int nz){
  	FoFTPtlStruct *p,tmp,*ptr,*bp,*sp;
	FoFTPtlStruct *sndp,*lp,*rp;
	MPI_Status status;
	int nid,myid;
	int src,dest;
	int nrecv,nsend,commsize;
	size_t i,j,k,np,sndnp;
	float zmin,zmax,zp,zup,zdown;

	bp = *Bp;
	np = *mp;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	zmin = ((float)nz/(float)nid)*myid;
	zmax = ((float)nz/(float)nid)*(myid+1);
	printf("P%d has zmin=%g zmax=%g\n",myid,zmin,zmax);

	commsize = 100;
	sp = bp;

	while(commsize != 0){
		rp = bp + np - 1;
		for(lp=sp;lp<=rp;){
			zp = lp->r[2];
			zup = fabs(zmax - zp);
			zup = min(zup,nz-zup);
			zdown = fabs(zmin-zp);
			zdown = min(zdown,nz-zdown);
			if((zp>=zmax||zp<zmin)&& zup<=zdown){
				tmp = *lp;
				*lp = *rp;
				*rp = tmp;
				rp --;
			}
			else lp ++;
		}
		sndnp = np-(rp-bp+1);
		sndp = (FoFTPtlStruct *)malloc(sizeof(FoFTPtlStruct)*sndnp);
		for(p=rp+1,ptr=sndp;p<bp+np;p++,ptr++) *ptr = *p;
		nsend = sndnp;
		dest = (myid+1+nid)%nid;
		src = (myid-1+nid)%nid;
		MPI_Sendrecv(&nsend,1,MPI_INT,dest,0,
				&nrecv,1,MPI_INT,src,0,MPI_COMM_WORLD,&status);
		bp = (FoFTPtlStruct *)realloc(bp,
				sizeof(FoFTPtlStruct)*(np-nsend+nrecv));
		*Bp = bp;
		MPI_Sendrecv(sndp,nsend*sizeof(FoFTPtlStruct),MPI_BYTE,dest,0,
				&(rp[1]),nrecv*sizeof(FoFTPtlStruct),MPI_BYTE,src,0,
				MPI_COMM_WORLD,&status);
		free(sndp);
		sp = bp+np-nsend;
		np = np-nsend+nrecv;
		GETCOMM(nsend,commsize);
		printf("+P%d has sent %d particle to P%d in total=%d \n",
				myid,sndnp,dest,commsize);
		fflush(stdout);
	}

	commsize =100;
	sp = bp;

	while(commsize != 0){
		rp = bp + np - 1;
		for(lp=sp;lp<=rp;){
			zp = lp->r[2];
			if(zp < zmin || zp >=zmax){
				tmp = *lp;
				*lp = *rp;
				*rp = tmp;
				rp --;
			}
			else lp ++;
		}
		sndnp = np-(rp-bp+1);
		sndp = (FoFTPtlStruct *)malloc(sizeof(FoFTPtlStruct)*sndnp);
		for(p=rp+1,ptr=sndp;p<bp+np;p++,ptr++) *ptr = *p;
		nsend = sndnp;
		dest = (myid-1+nid)%nid;
		src = (myid+1+nid)%nid;
		MPI_Sendrecv(&nsend,1,MPI_INT,dest,0,
				&nrecv,1,MPI_INT,src,0,MPI_COMM_WORLD,&status);
		bp = (FoFTPtlStruct *)realloc(bp,
				sizeof(FoFTPtlStruct)*(np-nsend+nrecv));
		*Bp = bp;
		MPI_Sendrecv(sndp,nsend*sizeof(FoFTPtlStruct),MPI_BYTE,dest,0,
				&(rp[1]),nrecv*sizeof(FoFTPtlStruct),MPI_BYTE,src,0,
				MPI_COMM_WORLD,&status);
		free(sndp);
		sp = bp+np-nsend;
		np = np-nsend+nrecv;
		GETCOMM(nsend,commsize);
		printf("-P%d has sent %d particles to P%d in total=%d \n",
				myid,sndnp,dest,commsize);
		fflush(stdout);
	}
	*mp = np;
}
