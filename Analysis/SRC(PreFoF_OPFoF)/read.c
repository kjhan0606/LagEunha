#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include <mpi.h>
#include "eunha.h"
#include "fof.h"
#define DEFINE_SIM_PARA
//#include "pmheader.h"
#undef DEFINE_SIM_PARA
#include "params.h"
#define MP 1000000


#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))



#define MIN(a,b) (a)<(b)? (a):(b)
#define MAX(a,b) (a)>(b)? (a):(b)
READTYPE rp[MP];
POSTYPE zstart,zwidth;
size_t np;
size_t  Fread(void *a,size_t b,size_t c, FILE *fp){
        char *A;
        char t1,t2,t3,t4;
        size_t i,nmem;
        nmem = fread(a,b,c,fp);
        A = (char *)a;
        for(i=0;i<b*nmem;i+=4){
                t1 = A[i];
                t2 = A[i+1];
                t3 = A[i+2];
                t4 = A[i+3];
                A[i] =t4;
                A[i+1] =t3;
                A[i+2] =t2;
                A[i+3] =t1;
        }
        return nmem;
}

#define Fread(a,b,c,d) fread(a,b,c,d)


#ifdef OLD
size_t readparticle(FoFTPtlStruct **Bp,size_t np,int nstep,int nfile,int nz,
		char *infile, int mp, SimParameters *simpar){
	FoFTPtlStruct *bp,*p;
	size_t j,k;
	int i;
	int next,before;
	FILE *fp;
	char outfile[190];
	FILE *wp;
	int nid,myid;
	size_t npadd;
	MPI_Status status;
	int itag;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);

	next = (myid+1+nid)%nid;
	before = (myid-1+nid)%nid;
	itag = 0;
	if(myid != 0) MPI_Recv(&i,1,MPI_INT,before, itag,MPI_COMM_WORLD,&status);
	{
			double zmin,zmax;
			SimParameters *rsim = simpar;
			fp=fopen(infile,"r");
			read_slab_head(fp, simpar);
			printf("P%d read the header of %s : %g %g : %g %g : %d\n", 
					myid, infile, ZMIN(simpar),ZMAX(simpar),HUBBLE(simpar),OMEP(simpar),XYZSHIFT(simpar));
			zstart = ZMIN(rsim);
			zwidth = ZMAX(rsim) - ZMIN(rsim);
			mp = DM_NP(rsim);
				
			bp = *Bp;
			bp = (FoFTPtlStruct *)realloc(bp,sizeof(FoFTPtlStruct)*(np+mp));
			*Bp = bp;
			p = bp + np;

			npadd = 0;
			printf("Now opening %s with %d particles from np= %d",infile,mp,np);
			fflush(stdout);
			zmin = 1e20;
			zmax = -1e20;
			while((mp=Fread(rp,sizeof(READTYPE),MP,fp))>0){
				for(j=0;j<mp;j++){
					p->r[0]=XofP(rsim,rp+j);
					p->r[1]=YofP(rsim,rp+j);
					p->r[2]=ZofP(rsim,rp+j);
					p->rv[0]=rp[j].vx;
					p->rv[1]=rp[j].vy;
					p->rv[2]=rp[j].vz;
					p->indx=PINDX(rp+j);
					zmin = MIN(zmin,p->r[2]);
					zmax = MAX(zmax,p->r[2]);
					npadd++;
					p++;
				}
			}
			printf("Well read %ld particles zmin: zmax   %lg : %lg\n",npadd,zmin,zmax);
			fflush(stdout);
			fclose(fp);
	}
	if(myid != nid-1) MPI_Send(&i,1,MPI_INT,next, itag,MPI_COMM_WORLD);
	return npadd;
}
#else
size_t readparticle(FoFTPtlStruct **Bp,size_t np,int nstep,int nfile,int nz,
		char *infile, int mp, SimParameters *simpar){
	FoFTPtlStruct *bp,*p;
	size_t j,k;
	int i;
	FILE *fp;
	char outfile[190];
	FILE *wp;
	int nid,myid;
	size_t npadd;
	MPI_Status status;
	int itag=1;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);

	int isend,iget;
	isend = iget = 1;
	int WGroupSize = WGSIZE;
	int src = myid -1;
	int tgt = myid + 1;
	if(RANKINGROUP(myid,WGroupSize) !=0)
		MPI_Recv(&i,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);
	npadd = 0;
	fp=fopen(infile,"r");
	if(fp){
			double zmin,zmax;
			SimParameters *rsim = simpar;
			printf("P%d tries to open %s\n", myid, infile);
			read_slab_head(fp, simpar);
			zstart = ZMIN(rsim);
			zwidth = ZMAX(rsim) - ZMIN(rsim);
			mp = simpar->np;
				
			bp = *Bp;
			bp = (FoFTPtlStruct *)realloc(bp,sizeof(FoFTPtlStruct)*(np+mp));
			*Bp = bp;
			p = bp + np;

			printf("P%d Now opening %s with %d particles from np= %d",myid,infile,mp,(int)np);
			fflush(stdout);
			zmin = 1e20;
			zmax = -1e20;
			while((mp=Fread(rp,sizeof(READTYPE),MP,fp))>0){
				for(j=0;j<mp;j++){
					p->r[0]=XofP(rsim,rp+j);
					p->r[1]=YofP(rsim,rp+j);
					p->r[2]=ZofP(rsim,rp+j);
					p->rv[0]=rp[j].vx;
					p->rv[1]=rp[j].vy;
					p->rv[2]=rp[j].vz;
					p->indx=rp[j].indx;
					zmin = MIN(zmin,p->r[2]);
					zmax = MAX(zmax,p->r[2]);
					npadd++;
					p++;
				}
			}
			printf("Well read %ld particles zmin: zmax   %lg : %lg\n",npadd,zmin,zmax);
			fflush(stdout);
			fclose(fp);
	}
	else {
		printf("P%d is now allocated with a file = %s . However, maybe, it is out of the file boundary\n",myid, infile);
	}
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid)
		MPI_Send(&i,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);
	return npadd;
}
#endif
