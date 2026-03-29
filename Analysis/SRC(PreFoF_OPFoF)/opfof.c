#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include <mpi.h>
#include "eunha.h"
//#include "pmheader.h"
#include "fof.h"
#include "Time.h"
#define MIN(a,b) a > b? b: a
#define MAX(a,b) a > b? a: b
#define pow2(a) ((a)*(a))
#define pow3(a) ((a)*(a)*(a))
particle p;
size_t readparticle(FoFTPtlStruct **,size_t ,int ,int ,int ,char *,int,SimParameters *);
POSTYPE fof_link = 0.2;
POSTYPE lx,ly,lz;

	float size,hubble,npower,omep,omepb,omeplam,bias,smooth;
	int nx,ny,nz,nspace;
	float ntree,ntree1,theta;
	float zinit,amax,astep,anow,a;
	int iseed;
	float vscale,rscale;
	FILE *hfp,*pfp;
	char halofile[190],memparticlefile[190];
void xhpsort(unsigned long n, particle ra[])
{
    unsigned long i,ir,j,l;
    particle rra;

    if (n < 2) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l];
        } else {
            rra=ra[ir];
            ra[ir]=ra[1];
            if (--ir == 1) {
                ra[1]=rra;
                break;
            }
        }
        i=l;
        j=l+l;
        while (j <= ir) {
            if (j < ir && ra[j].x < ra[j+1].x) j++;
            if (rra.x < ra[j].x) {
                ra[i]=ra[j];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i]=rra;
    }
}
void yhpsort(unsigned long n, particle ra[])
{
    unsigned long i,ir,j,l;
    particle rra;

    if (n < 2) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l];
        } else {
            rra=ra[ir];
            ra[ir]=ra[1];
            if (--ir == 1) {
                ra[1]=rra;
                break;
            }
        }
        i=l;
        j=l+l;
        while (j <= ir) {
            if (j < ir && ra[j].y < ra[j+1].y) j++;
            if (rra.y < ra[j].y) {
                ra[i]=ra[j];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i]=rra;
    }
}
void zhpsort(unsigned long n, particle ra[])
{
    unsigned long i,ir,j,l;
    particle rra;

    if (n < 2) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l];
        } else {
            rra=ra[ir];
            ra[ir]=ra[1];
            if (--ir == 1) {
                ra[1]=rra;
                break;
            }
        }
        i=l;
        j=l+l;
        while (j <= ir) {
            if (j < ir && ra[j].z < ra[j+1].z) j++;
            if (rra.z < ra[j].z) {
                ra[i]=ra[j];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i]=rra;
    }
}

/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */

HaloQ haloproperty(particle *member,size_t nmem){
	HaloQ halo;
	particle *pp;
	size_t i,j,k;
	double cx,cy,cz;
	double vx,vy,vz;
	POSTYPE xmin,xmax,ymin,ymax,zmin,zmax;
	POSTYPE x1,x2,y1,y2,z1,z2;
	halo.np=nmem;
	cx=cy=cz=vx=vy=vz=0;
	pp=member;
	xmax = ymax = zmax = -1.E25;
	xmin = ymin = zmin = +1.E25;
	for(i=0;i<nmem;i++){
		xmax = MAX(xmax,(pp->x));
		xmin = MIN(xmin,(pp->x));
		ymax = MAX(ymax,(pp->y));
		ymin = MIN(ymin,(pp->y));
		zmax = MAX(zmax,(pp->z));
		zmin = MIN(zmin,(pp->z));
		pp++;
	}
	if(xmin <= fof_link && xmax >=lx-fof_link){
		xhpsort(nmem,member-1);
		pp=member+1;
		for(i=1;i<nmem;i++){
			if(((pp->x)-(pp-1)->x) > 1.5*fof_link) {
				break;
			}
			pp++;
		}
		for(j=i;j<nmem;j++){
			member[j].x -= lx;
		}
	}
    if(ymin <= fof_link && ymax >=ly-fof_link){
                yhpsort(nmem,member-1);
                pp=member+1;
                for(i=1;i<nmem;i++){
                        if(((pp->y)-(pp-1)->y) > 1.5*fof_link) {
                                break;
                        }
						pp++;
                }
                for(j=i;j<nmem;j++){
                        member[j].y -= ly;
                }
     }
     if(zmin <= fof_link && zmax >=lz-fof_link){
                zhpsort(nmem,member-1);
                pp=member+1;
                for(i=1;i<nmem;i++){
                        if(((pp->z)-(pp-1)->z) > 1.5*fof_link) {
                                break;
                        }
						pp++;
                }
                for(j=i;j<nmem;j++){
                        member[j].z -= lz;
                }
     }


	pp=member;
	for(i=0;i<nmem;i++){
		cx+=(pp->x);
		cy+=(pp->y);
		cz+=(pp->z);
		vx+=pp->vx;
		vy+=pp->vy;
		vz+=pp->vz;
		pp++;
	}
	/* position in h^-1 Mpc */
	halo.x = cx/(POSTYPE)nmem *rscale;
	halo.y = cy/(POSTYPE)nmem *rscale;
	halo.z = cz/(POSTYPE)nmem *rscale;
	/* velocity in km/sec */
	halo.vx = vx/(float)nmem*vscale;
	halo.vy = vy/(float)nmem*vscale;
	halo.vz = vz/(float)nmem*vscale;
	return halo;
}
int myid,nid,mid;
SimParameters simpar;
int main(int argc,char *argv[]){
	double std,mean;
	int ntmp;
	size_t num,ii;
	POSTYPE tmpx,tmpy,tmpz,dist2;
	POSTYPE fplmf,ptlmass;
	size_t i,j,k;
	int si,snp;
	int nowfile;
	int N,M;
	size_t N3;
	size_t np,addhere;
	int nstep;
	FoFTPtlStruct *ptl;
	FoFTPtlStruct *ptr;
	particle *linked;
	int nfof;
	FoFBeginEndTree beginend;
	Box box;
	FoFTStruct *TREE;
	long long ntreemax = 1000000L;
	float wtime;
	int nfile;
	FILE *fp;
	size_t nhalo;
	POSTYPE zmin,zmax,zminlocal;
	size_t npadd,npwrite,npread;
	MPI_Status status;
	int initfile,finalfile;
	char infile[190];

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);



	if(argc != 3) {
		if(myid==0){
			fprintf(stderr,"Error in # of arguments\n");
			fprintf(stderr,"%%fof fileleadingheadername nstep\n");
		}
		MPI_Finalize();
		exit(199);
	}
	else {
		FILE *ffp;
		double r2kineticfact,HSUB;
		nstep = atoi(argv[2]);
		sprintf(infile,"%s.%.5d.%.5d",argv[1],nstep,myid);
		if(myid==0){
			printf("Now trying to open %s\n", infile);
			ffp = fopen(infile,"r");
			void read_slab_head(FILE *, SimParameters *);
			read_slab_head(ffp, &simpar);
			fclose(ffp);
		}
		MPI_Bcast(&simpar,sizeof(SimParameters),MPI_BYTE,0,MPI_COMM_WORLD);
		SimParameters *rsim = &simpar;
		nfile = NID(rsim);
		size = BOXSIZE(rsim);
		hubble = HUBBLE(rsim); 
		npower = NPOW(rsim); 
		omep = OMEP(rsim);
		omepb = OMEPB(rsim); 
		omeplam = OMEPLAM(rsim);
		bias = BIAS8(rsim);
		smooth = 0;
		nx = NX(rsim);
		ny = NY(rsim);
		nz = NZ(rsim);
		nspace = NSPACE(rsim);
		theta = THETA(rsim); 
		zinit = AMAX(rsim) - 1.; 
		astep = ASTEP(rsim); 
		anow = ANOW(rsim); 
        ny=nz=nx;
		a = anow;
		amax = AMAX(rsim);
		nspace =1;
		N3 = (size_t)(nx/nspace)*(size_t)(ny/nspace)*(size_t)(nz/nspace)/nfile;
		N3 = N3*0.5;
		rscale = size/nx;
		HSUB = sqrt(omep*pow3(amax/a)+omeplam+(1.-omep-omeplam)*pow2(amax/a));
		vscale = size/amax*a*a*100.*HSUB;
           
		fof_link = fof_link * nspace;
		if(myid==0){
			printf("P%d size = %g hubble = %g\n",myid,size,hubble);
			printf("P%d npower = %g omep = %g omeplam = %g bias = %g smooth=%g\n",
					myid,npower,omep,omeplam,bias,smooth);
			printf("P%d zinit = %g astep = %g anow = %g \n",myid,zinit,astep,anow);
			printf("P%d nx = %d ny= %d nz= %d nspace=%d\n",myid,nx,ny,nz,nspace);
			printf("P%d rscale = %g vscale= %g\n",myid,rscale,vscale);
		}
	}

	if((ptl = (FoFTPtlStruct *) malloc(sizeof(FoFTPtlStruct)*100)) == NULL){
		fprintf(stderr,"Error allocating ptl %ld\n",N3);
		exit(99);
	}
	M = 1;
	linked = (particle *) malloc(sizeof(particle)*MaxLinkedParticles);
	if((TREE = (FoFTStruct *) malloc(sizeof(FoFTStruct)*ntreemax)) == NULL){
                fprintf(stderr,"Error allocating TREE\n");
                exit(99);
        }
	wtime = WALLCLOCK();
	(void)WALLCLOCK();
	box.x = box.y = box.z = 0.;
	box.width = nx;
	sprintf(halofile,"FoF_halo_cat.%.5d",nstep);
	sprintf(memparticlefile,"FoF_member_particle.%.5d",nstep);
	if(myid==0)
	{
		hfp=fopen(halofile,"w");
		pfp=fopen(memparticlefile,"w");
		fwrite(&size,sizeof(float),1,hfp);
		fwrite(&hubble,sizeof(float),1,hfp);
		fwrite(&npower,sizeof(float),1,hfp);
		fwrite(&omep,sizeof(float),1,hfp);
		fwrite(&omepb,sizeof(float),1,hfp);
		fwrite(&omeplam,sizeof(float),1,hfp);
		fwrite(&bias,sizeof(float),1,hfp);
		fwrite(&nx,sizeof(int),1,hfp);
		fwrite(&nspace,sizeof(int),1,hfp);
		fwrite(&amax,sizeof(float),1,hfp);
		fwrite(&astep,sizeof(float),1,hfp);
		fwrite(&anow,sizeof(float),1,hfp);
		fclose(hfp);
	}
	lx = nx;
	ly = ny;
	lz = nz;
	printf("P%d has Lx=%g Ly=%g Lz=%g\n",myid,lx,ly,lz);
	np = 0;
	npwrite = 0;
#ifdef USE_MASTER
	mid = nid -1;
#error Not yet implmented
#else
	mid = nid;
#endif
#ifdef PREVIOUS_VERSION
	initfile = (nfile/mid)*myid;
	finalfile = (nfile/mid)*(myid+1);
#else
	/* realintifile & realfinalfile are the real file range
	 * but (initfile, finalfile) could be a fake-file range to be set to synchronize the 
	 * parallel process */
	int realinitfile,realfinalfile;
	int LocalNFiles = (nfile-1 +mid)/mid;
	int n_pes = (nfile + LocalNFiles - 1)/LocalNFiles;
	{
		if(myid >= n_pes) {
			/* These are for a fake  to redundant ranks to whom no files are allocated */
			initfile = nfile;
			finalfile = nfile + LocalNFiles;
			realinitfile = realfinalfile = -1;
		}
		else {
			realinitfile = myid * LocalNFiles;
			if(myid == n_pes-1) realfinalfile = nfile;
			else realfinalfile = realinitfile + LocalNFiles;
			initfile = realinitfile;
			finalfile = initfile + LocalNFiles;
		}
	}
#endif
	for(nowfile=initfile;nowfile<finalfile;nowfile++){
		sprintf(infile,"%s.%.5d.%.5d",argv[1],nstep,nowfile);
		printf("P%d is opening %s\n", myid, infile);
		npadd = readparticle(&ptl,np,nstep,nowfile,nz,infile,N3, &simpar);
		np += npadd;
		if(nowfile == finalfile-1){
			int src,dest;
			src = (myid+1+mid)%mid;
			dest = (myid-1+mid)%mid;
#ifndef PREVIOUS_VERSION
			if(myid == n_pes-1) { src = 0;}
			if(myid == 0) { dest = n_pes-1;}
			if(myid < n_pes)
			{
#endif
				MPI_Sendrecv(&npwrite,sizeof(size_t),MPI_BYTE,dest,0,
						&npread,sizeof(size_t),MPI_BYTE,src,0,MPI_COMM_WORLD,&status);
				{
					ptl = (FoFTPtlStruct*)realloc(ptl,
							sizeof(FoFTPtlStruct)*(np+npread));
					ReadBottomFaceContact(ptl+np,npread,linked,src,nstep,nz);
					np += npread;
				}
			}
		}
		HaloBound *halobound;
		nhalo = 0;
#ifndef PREVIOUS_VERSION
		if((nowfile < nfile || nowfile == finalfile-1) && np >0)
#endif
		{
			zmin = 2.E23;
			zmax = -2.E23;
			for(i=0;i<np;i++){
				zmin = MIN(zmin,ptl[i].r[2]);
				zmax = MAX(zmax,ptl[i].r[2]);
			}
	
			if(nowfile==initfile) zminlocal = zmin;

			printf("Now we have zmin=%g zmax=%g\n",zmin,zmax);
			ptl = (FoFTPtlStruct *)realloc(ptl,sizeof(FoFTPtlStruct)*np);
			if(np /NODE_HAVE_PARTICLE *1.9 > ntreemax){
				ntreemax = np/NODE_HAVE_PARTICLE *1.5;
				TREE = (FoFTStruct *)realloc(TREE,sizeof(FoFTStruct)*ntreemax);
			}
			for(i=0;i<np;i++){
				ptl[i].type = TYPE_PTL;
				ptl[i].sibling = &ptl[i+1];
				ptl[i].haloindx = -1;
			}
			ptl[np-1].sibling = NULL;
			FoF_Make_Tree(TREE,ptl,np,box);
			/* Local FoF and tag with a new halo number nhalo */
			for(i=0;i<np;i++){
				if(ptl[i].included == NO){
					p.x = ptl[i].r[0];
					p.y = ptl[i].r[1];
					p.z = ptl[i].r[2];
					num=pnew_fof_link(&p,fof_link,TREE,ptl,linked,nhalo,lx,ly,lz);
					nhalo ++;
				}
			}
			printf("P%d has %ld halos from %ld particles\n",myid,nhalo,np);
			halobound = (HaloBound *)malloc(sizeof(HaloBound)*nhalo);
			CheckHaloBound(nhalo,halobound,ptl,np,fof_link,
				zmin,zmax, zminlocal);
			printf("P%d has checked halo boundary\n",myid);fflush(stdout);
		}
		{
			if(nowfile != finalfile-1) {
				int flag;
				WriteIsolatedHalo(nhalo,halobound,ptl,linked,halofile,
						memparticlefile);
				printf("P%d has saved isolated halos\n",myid);fflush(stdout);
				if(nowfile == initfile) flag = 0;
				else flag = 1;
				npwrite += WriteBottomFaceContact(nhalo,halobound,ptl,linked,
						flag,nstep);
#ifndef PREVIOUS_VERSION
				if(nowfile < nfile)
#endif
				{
					np = StackUpContactParticleLeftWard(nhalo,halobound,ptl,np);
					printf("P%d has saved Bottom Faced Halo %ld\n",myid,npwrite);fflush(stdout);
					free(halobound);
				}
			}
			if(nowfile == finalfile-1){
				WriteFinalHalo(nhalo,halobound,ptl,linked,halofile,
						memparticlefile);
				MPI_Finalize();
				return 0;
			}
		}
	}
}
