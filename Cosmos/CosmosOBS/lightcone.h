/* 여기서  theta, phi 는 모두 0 ` 90 도 이내 이어야 한다. */
/* 만약 바꿀려면 프로그램중에서 min,max 의 값과 acos(xp/sqrt(xp*xp+yp*yp)) 
 * 을 바꿔야한다. */
#define pi 3.1415926535L
#define X0 0.
#define Y0 0.
#define Z0 0.
#define NPERIODIC 8 /*Number of saving bits is 8*sizeof(int)*8bits */
typedef struct InBox {
    int nx,ny,nz;
} InBox;
typedef struct slcparticletype{
    treepmparticletype *bp;
#ifdef INDEX
    indextype indx;
#endif
#ifdef XYZDBL
    double x,y,z;
#else
    float x,y,z;
#endif
    float vx,vy,vz;
    float ovx,ovy,ovz,nvx,nvy,nvz;
    float dist;
    int type;
} slcparticletype;
typedef struct lightconetype{
#ifdef XYZDBL
    double x,y,z;
#else
    float x,y,z;
#endif
    float vx,vy,vz;
#ifdef INDEX
    indextype indx;
#endif
} lightconetype;

#ifdef OLD_OBS_SWAP
#define OBSTEMWRITE(bp,np,myid,nid,obsid) {\
	FILE *wp; int i,j,k;\
	char tempfile[190];\
	for(i=0;i<nid;i++){\
		sprintf(tempfile,"tempobs%.5d-id%.5d.dat",obsid,myid);\
		if(myid==i){\
			wp = fopen(tempfile,"w");\
			fwrite(&np,sizeof(int),1,wp);\
			fwrite(bp,sizeof(slcparticletype),np,wp);\
			fclose(wp);\
		}\
		MPI_Barrier(MPI_COMM_WORLD);\
	}\
	Free(bp);\
}
#define OBSTEMREAD(bp,np,myid,nid,obsid) {\
	FILE *fp; int i,j,k;\
	char tempfile[190];\
	for(i=0;i<nid;i++){\
		sprintf(tempfile,"tempobs%.5d-id%.5d.dat",obsid,myid);\
		if(myid==i){\
			np=0;\
			fp = fopen(tempfile,"r");\
			fread(&np,sizeof(int),1,fp);\
			if(np>0){\
				bp = (slcparticletype *)Malloc(sizeof(slcparticletype)*np,PPTR(bp));\
				fread(bp,sizeof(slcparticletype),np,fp);\
			}\
			fclose(fp);\
		}\
		MPI_Barrier(MPI_COMM_WORLD);\
	}\
}
#else
#define OBSTEMWRITE(bp,np,myid,nid,obsid) {\
	FILE *wp; int i,j,k,mp;\
	char tempfile[190];\
	MPI_Status status;\
	if(myid==0){\
		sprintf(tempfile,"tempobs%.5d.dat",obsid);\
		wp = fopen(tempfile,"w");\
		fwrite(&np,sizeof(int),1,wp);\
		fwrite(bp,sizeof(slcparticletype),np,wp);\
		Free(bp);\
		for(i=1;i<nid;i++){\
			MPI_Recv(&np,1,MPI_INT,i,i,MPI_COMM_WORLD,&status);\
			fwrite(&np,sizeof(int),1,wp);\
			if(np>0) {\
				bp = (slcparticletype *)Malloc(sizeof(slcparticletype)*np,PPTR(bp));\
				MPI_Recv(bp,np*sizeof(slcparticletype),MPI_BYTE,i,i,MPI_COMM_WORLD,&status);\
				fwrite(bp,sizeof(slcparticletype),np,wp);\
			}\
			Free(bp);\
		}\
		fclose(wp);\
	}\
	else {\
		MPI_Send(&np,1,MPI_INT,0,myid,MPI_COMM_WORLD);\
		if(np>0) MPI_Send(bp,np*sizeof(slcparticletype),MPI_BYTE,0,myid,MPI_COMM_WORLD);\
		Free(bp);\
	}\
}
#define OBSTEMREAD(bp,np,myid,nid,obsid) {\
	FILE *fp; int i,j,k,mp;\
	char tempfile[190];\
	MPI_Status status;\
	if(myid==nid-1){\
		sprintf(tempfile,"tempobs%.5d.dat",obsid);\
		fp = fopen(tempfile,"r");\
		for(i=0;i<nid;i++){\
			fread(&np,sizeof(int),1,fp);\
			if(np>0) {\
				bp = (slcparticletype*)Malloc(sizeof(slcparticletype)*np,PPTR(bp));\
				fread(bp,sizeof(slcparticletype),np,fp);\
			}\
			if(i<nid-1){\
				MPI_Send(&np,1,MPI_INT,i,i,MPI_COMM_WORLD);\
				if(np>0) {\
					MPI_Send(bp,np*sizeof(slcparticletype),MPI_BYTE,i,i,MPI_COMM_WORLD);\
					Free(bp);\
				}\
			}\
		}\
		fclose(fp);\
	}\
	else {\
		MPI_Recv(&np,1,MPI_INT,nid-1,myid,MPI_COMM_WORLD,&status);\
		if(np>0){\
			bp = (slcparticletype*)Malloc(sizeof(slcparticletype)*np,PPTR(bp));\
			MPI_Recv(bp,np*sizeof(slcparticletype),MPI_BYTE,nid-1,myid,MPI_COMM_WORLD,&status);\
		}\
	}\
}
#endif


#define writedownlightconedata(A,B,C) {\
    FILE *wp;\
    float x0,y0,z0;\
    int nmlcp;\
    lightconetype *ttmp;\
    x0 = SX0;\
    y0 = SY0;\
    z0 = SZ0;\
    {\
        if(myid==0) {\
            wp=fopen(A,"w");\
            fwrite(&nx,sizeof(int),1,wp);\
            fwrite(&ny,sizeof(int),1,wp);\
            fwrite(&nz,sizeof(int),1,wp);\
            fwrite(&omep,sizeof(float),1,wp);\
            fwrite(&omeplam,sizeof(float),1,wp);\
            fwrite(&omepb,sizeof(float),1,wp);\
            fwrite(&hubble,sizeof(float),1,wp);\
            fwrite(&boxsize,sizeof(float),1,wp);\
            fwrite(&amax,sizeof(float),1,wp);\
            fwrite(&a,sizeof(float),1,wp);\
            fwrite(&astep,sizeof(float),1,wp);\
            fwrite(&x0,sizeof(float),1,wp);\
            fwrite(&y0,sizeof(float),1,wp);\
            fwrite(&z0,sizeof(float),1,wp);\
            fwrite(&rmin2,sizeof(double),1,wp);\
            fwrite(&rmax2,sizeof(double),1,wp);\
            fwrite(&C,sizeof(int),1,wp);\
            fwrite(B,sizeof(lightconetype),C,wp);\
            for(i=1;i<nid;i++){ \
                MPI_Recv(&nmlcp,1,MPI_INT,i,i,MPI_COMM_WORLD,&status);\
                if(nmlcp>0)ttmp=(lightconetype*)Malloc(sizeof(lightconetype)*nmlcp,PPTR(ttmp));\
                MPI_Recv(ttmp,sizeof(lightconetype)*nmlcp,MPI_BYTE,i,i,\
                        MPI_COMM_WORLD,&status);\
                fwrite(&nmlcp,sizeof(int),1,wp);\
                fwrite(ttmp,sizeof(lightconetype),nmlcp,wp);\
                if(nmlcp>0)Free(ttmp);\
            }\
            fclose(wp);\
        }\
        else  {\
            MPI_Send(&C,1,MPI_INT,0,myid,MPI_COMM_WORLD);\
            MPI_Send(B,sizeof(lightconetype)*C,MPI_BYTE,0,myid,MPI_COMM_WORLD);\
        }\
    }\
}
