#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
#include "voro.h"
#include "kh.h"
#include "DD2d.h"

#define Nneigh 50
Voro2D_point neigh[Nneigh];


void outdata(SimParameters *simpar, int nstep, postype t, postype dt){
	treevoroparticletype *bp = VORO_TBP(simpar);
	int np = VORO_NP(simpar);
	int myid, nid;
	myid = MYID(simpar);
	nid = NID(simpar);
	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);

	char outfile[190];
	sprintf(outfile,"khold.%.6d.dat",nstep);
	int i;
	int snp=0;
	FILE *wp;
	if(myid ==0){
		wp = fopen(outfile,"w");
		int tnp = VORO_TNP(simpar);
		fwrite(&tnp, sizeof(int),1,wp);
		fwrite(bp, sizeof(treevoroparticletype), np, wp);
		long offset = sizeof(int)+sizeof(treevoroparticletype)*VORO_TNP(simpar);
		fseek(wp, offset, SEEK_SET);
		fwrite(&t, sizeof(postype), 1, wp);
		fwrite(&dt, sizeof(postype), 1, wp);
		fwrite(&GAS_dtold(simpar), sizeof(float), 1, wp);
		fwrite(&nx, sizeof(int), 1, wp);
		fwrite(&ny, sizeof(int), 1, wp);
		fwrite(&(SIMBOX(simpar).x.max), sizeof(postype), 1, wp);
		fwrite(&(SIMBOX(simpar).y.max), sizeof(postype), 1, wp);
		fwrite(&GAS_AlphaVis(simpar), sizeof(float), 1, wp);
		fwrite(&GAS_BetaVis(simpar), sizeof(float), 1, wp);
		fclose(wp);
		int rnp;
		rnp = np;
		MPI_Send(&rnp, 1, MPI_INT, myid+1, 0, MPI_COMM_WORLD);
	}
	else{
		int rnp;
		MPI_Status status;
		MPI_Recv(&snp, 1, MPI_INT, myid-1, 0, MPI_COMM_WORLD,&status);
		long offset = sizeof(int)+sizeof(treevoroparticletype)*snp;
		wp = fopen(outfile,"r+");
		fseek(wp, offset, SEEK_SET);
		fwrite(bp, sizeof(treevoroparticletype), np, wp);
		fclose(wp);
		rnp = snp + np;
		if(myid != nid-1){
			MPI_Send(&rnp, 1, MPI_INT, myid+1, 0, MPI_COMM_WORLD);
		}
	}
}


void readdata(SimParameters *simpar, postype *t, postype *dt,int nstep){
	char infile[190];
	sprintf(infile,"khout.%.6d.dat", nstep);
	int myid = MYID(simpar);
	float dtold;
	if(myid ==0){
		FILE *fp = fopen(infile,"r");
		int np;
		fread(&np, sizeof(int), 1, fp);
		VORO_TNP(simpar) = VORO_NP(simpar) = np;
		VORO_TBP(simpar) = (treevoroparticletype*)malloc(sizeof(treevoroparticletype)*np);
		if(0){
			fread(VORO_TBP(simpar),sizeof(treevoroparticletype), np,fp);
		}
		else {
			treevorork4particletype *aa = malloc(sizeof(treevorork4particletype)*np);
			fread(aa,sizeof(treevorork4particletype), np,fp);
			treevoroparticletype *bb = VORO_TBP(simpar);
			int i;
			for(i=0;i<np;i++){
				bb[i].u4if = aa[i].u4if;
				bb[i].mass = aa[i].mass;
				bb[i].x = aa[i].x;
				bb[i].y = aa[i].y;
				bb[i].vx = aa[i].vx;
				bb[i].vy = aa[i].vy;
				bb[i].ie = aa[i].ie;
				bb[i].csound = aa[i].csound;
				bb[i].volume = aa[i].volume;
				bb[i].w2 = aa[i].w2;
				bb[i].w2old = aa[i].w2old;
				bb[i].w2ceil = aa[i].w2ceil;
				bb[i].den = aa[i].den;
				bb[i].pressure = aa[i].pressure;
			}
			free(aa);
		}
		fread(t,sizeof(postype), 1,fp);
		fread(dt,sizeof(postype), 1,fp);
		fread(&dtold,sizeof(float), 1,fp);
		DEBUGPRINT("P%d is detecting total np= %d t= %g dt= %g\n", MYID(simpar), np, *t, *dt);
		fclose(fp);
	}
	else {
		VORO_NP(simpar) = 0;
		VORO_TBP(simpar) = (treevoroparticletype*)malloc(sizeof(treevoroparticletype)*100);;
	}
	MPI_Bcast(&VORO_TNP(simpar), 1, MPI_INT64_T,0, MPI_COMM(simpar));
	MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(dt, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(&dtold, 1, MPI_FLOAT,0, MPI_COMM(simpar));
	KH_TreeAllParticleMigrate(simpar);
	treevoroparticletype *a = VORO_TBP(simpar);
	GAS_dtold(simpar) = dtold;
	/*
	DEBUGPRINT("P%d has read data p0: mass/x/y/vx/vy/ie/csound/volume/w2/den/pressure= %g %g %g %g %g"
			" %g %g %g %g %g %g\n", MYID(simpar),a->mass,a->x,a->y,a->vx,a->vy,a->ie,a->csound,
			a->volume,a->w2, a->den, a->pressure);
			*/
}


int RunKHOLD(SimParameters *simpar, int icont){
	postype t,dt;
	int np;
	treevoroparticletype *bp;
	int icount = 0;;
	int iflag,jflag;


	{
		GridInfo *grid = &(GRIDINFO(simpar));
		PosNX(grid) = NX(simpar);
		PosNY(grid) = NY(simpar);
		PosNZ(grid) = NZ(simpar);
		PosNXNY(grid) = NX(simpar)* NY(simpar);
	}
	void StartKHRkSDD(SimParameters *);
	StartKHRkSDD(simpar); /* This makes all the ddinfo including pivot */


	printf("P%d has volume: rmin= %g %g rmax= %g %g\n", MYID(simpar),
			SIM_LXMIN(simpar,voro), SIM_LYMIN(simpar,voro),
			SIM_LXMAX(simpar,voro), SIM_LYMAX(simpar,voro));
	KH_XMIN(simpar) = SIM_LXMIN(simpar,voro);
	KH_YMIN(simpar) = SIM_LYMIN(simpar,voro);
	KH_XMAX(simpar) = SIM_LXMAX(simpar,voro);
	KH_YMAX(simpar) = SIM_LYMAX(simpar,voro);

	VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*100);



	if(icont ==0) {
		t = 0;
		dt = 1.e-7;
		bp = mkinitial(simpar, &np);
		printf("made initial conditions\n");
		GAS_dtold(simpar) = dt;
	}
	else{
		postype dt;
		int nstep = icont;
		readdata(simpar, &t,&dt, nstep);
		icount = nstep;
	}

	do {
//		MkLinkedList(bp,np);
		postype dt = vph2D(simpar);
		t += dt;
		if(MYID(simpar)==0){
			printf("Time is %g with icount = %d\n",t, icount);
			fflush(stdout);
		}
		GAS_dtold(simpar) = dt;


		iflag = t * 10.;
		jflag = (t-dt) * 10.;

		icount ++;

		if(iflag != jflag || icount%100 ==0 || icount == 1 || icount == 2)
		{
			outdata(simpar,  icount, t, dt);
			makemap2D(simpar, icount);
		}




	} while(t<10.);
//	} while(icount < 30);
	return 1;
}
