#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
#include "voro.h"
#include "rt.h"
#include "DD2d.h"

#define Nneigh 50
Voro2D_point neigh[Nneigh];


void rt_outdata(SimParameters *simpar, int nstep, postype t, postype dt){
	treevoroparticletype *bp = VORO_TBP(simpar);
	int np = VORO_NP(simpar);
	int myid, nid;
	myid = MYID(simpar);
	nid = NID(simpar);
	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);


	char outfile[190];
	sprintf(outfile,"rtout.%.6d.dat",nstep);
	int i;
	int snp=0;
	FILE *wp;
	int tnp = VORO_TNP(simpar);
	for(i=0;i<nid;i++){
		if(myid ==i){
			if(myid==0){
				wp = fopen(outfile,"w"); 
				fwrite(&tnp, sizeof(int),1,wp);
			}
			else wp = fopen(outfile,"a");
			fwrite(bp, sizeof(treevoroparticletype), np, wp);
			fclose(wp);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if(myid==0){
		wp = fopen(outfile,"a");
		fwrite(&t, sizeof(postype), 1, wp);
		fwrite(&nx, sizeof(int), 1, wp);
		fwrite(&ny, sizeof(int), 1, wp);
		fwrite(&(SIMBOX(simpar).x.max), sizeof(postype), 1, wp);
		fwrite(&(SIMBOX(simpar).y.max), sizeof(postype), 1, wp);
		fwrite(&GAS_AlphaVis(simpar), sizeof(float), 1, wp);
		fwrite(&GAS_BetaVis(simpar), sizeof(float), 1, wp);
		fclose(wp);
	}
	/*
	MPI_Barrier(MPI_COMM_WORLD);
	DEBUGPRINT("P%d is now exiting the dump\n", myid);
	*/
}


void rt_readdata(SimParameters *simpar, postype *t, int nstep){
	char infile[190];
	sprintf(infile,"rtout.%.6d.dat", nstep);
	int myid = MYID(simpar);
	if(myid ==0){
		FILE *fp = fopen(infile,"r");
		int np;
		fread(&np, sizeof(int), 1, fp);
		VORO_TNP(simpar) = VORO_NP(simpar) = np;
		VORO_TBP(simpar) = (treevoroparticletype*)malloc(sizeof(treevoroparticletype)*np);
		fread(VORO_TBP(simpar),sizeof(treevoroparticletype), np,fp);
		fread(t,sizeof(postype), 1,fp);
		fclose(fp);
	}
	else {
		VORO_NP(simpar) = 0;
		VORO_TBP(simpar) = (treevoroparticletype*)malloc(sizeof(treevoroparticletype)*100);;

	}
	RT_TreeAllParticleMigrate(simpar);
	MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Barrier(MPI_COMM(simpar));
}


int RunRT(SimParameters *simpar, int icont){
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
	void StartRTRkSDD(SimParameters *);
	StartRTRkSDD(simpar); /* This makes all the ddinfo including pivot */


	printf("P%d has volume: rmin= %g %g rmax= %g %g\n", MYID(simpar),
			SIM_LXMIN(simpar,voro), SIM_LYMIN(simpar,voro),
			SIM_LXMAX(simpar,voro), SIM_LYMAX(simpar,voro));
	RT_XMIN(simpar) = SIM_LXMIN(simpar,voro);
	RT_YMIN(simpar) = SIM_LYMIN(simpar,voro);
	RT_XMAX(simpar) = SIM_LXMAX(simpar,voro);
	RT_YMAX(simpar) = SIM_LYMAX(simpar,voro);

	VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*100);



	if(icont ==0) {
		t = 0;
		bp = rt_mkinitial(simpar, &np);
		printf("made initial conditions\n");
	}
	else{
		postype dt;
		int nstep = icont;
		rt_readdata(simpar, &t, nstep);
		icount = nstep+1;
	}

	do {
//		MkLinkedList(bp,np);
		postype dt = rt_vph2D(simpar);
		t += dt;
		if(MYID(simpar)==0){
			printf("Time is %g with icount = %d\n",t, icount);
			fflush(stdout);
		}
		DEBUGPRINT("P%d is now at t = %g with dt= %g\n", MYID(simpar), t,dt);

		/*
		if(icount == 152 && MYID(simpar) == 2){
			FILE *wp = fopen("checkic.dat","w");
			treevoroparticletype *bp = VORO_TBP(simpar);
			int i;
			for(i=0;i<VORO_NP(simpar);i++){
				fprintf(wp,"%g %g %g %g\n", bp[i].x, bp[i].y, bp[i].vx, bp[i].vy);
			}
			fclose(wp);
		}
		*/

		iflag = t * 10.;
		jflag = (t-dt) * 10.;


		if(iflag != jflag)
		{
			rt_outdata(simpar,  icount, t, dt);
			rt_makemap(simpar, icount);
		}


		icount ++;

	} while(t<10.);
	return 1;
}
