#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include "eunha.h"
#include "kjhrw.h"
#include "params.h"
#include "timerutil.h"
//#include "mpirms.h"
#include "domaindecomp.h"


/*
float cputime0[10];
float cputime1[10];
*/


#define MAX(a,b) ((a) > (b) ? (a): (b))

#define MIN(a,b) ((a) < (b) ? (a): (b))




#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))


#include<sys/times.h>

/*
static SimParameters simpar_backup;

void BackUpbeforeReadHead(SimParameters *simpar){
	simpar_backup = *simpar;
}
void RecoverAfterReadHead(SimParameters *simpar){
	*PTR2MPIINFO(simpar) = *PTR2MPIINFO(simpar_backup);
	DM_DDFUNC(simpar) = DM_DDFUNC(simpar_backup);
	SPH_DDFUNC(simpar) = SPH_DDFUNC(simpar_backup);
	STAR_DDFUNC(simpar) = STAR_DDFUNC(simpar_backup);
	AGN_DDFUNC(simpar) = AGN_DDFUNC(simpar_backup);
	TDM_DDFUNC(simpar) = TDM_DDFUNC(simpar_backup);
	TSPH_DDFUNC(simpar) = TSPH_DDFUNC(simpar_backup);
	TSTAR_DDFUNC(simpar) = TSTAR_DDFUNC(simpar_backup);
	TAGN_DDFUNC(simpar) = TAGN_DDFUNC(simpar_backup);
}
*/

void jread(SimParameters *simpar){
	FILE *rvfile;
	int myid, nid;
	MPI_Status status;
	char infile[190];
	float cputime0[10];
	float cputime1[10];

	myid = MYID(simpar);
	nid = NID(simpar);


	int src, tgt, WGroupSize, isend, iget, itag = 1;
	WGroupSize = WGROUPSIZE(simpar);
	src = myid -1;
	tgt = myid +1;

#ifdef GOTPM
	sprintf(infile,"%s.%.5d", RV_FILE(simpar),MYID(simpar));
#else
	sprintf(infile,"%s.%.5d.%.5d", RV_FILE(simpar), IndT_ISUBSTEP(simpar), MYID(simpar));
#endif
	SimParameters testpar;
	if(myid ==0){
		rvfile = fopen(infile,"r");
		read_head(rvfile, &testpar);
		fclose(rvfile);
	}

	if(RANKINGROUP(myid, WGroupSize) !=0) 
		MPI_Recv(&iget, 1, MPI_INT, src, itag, MPI_COMM(simpar), &status);
	fprintf(stderr,"P%d starts reading data: %s : sum= %ld : %ld, %ld, %ld, %ld   ",myid,
			infile, NPSUM(simpar), DM_NP(simpar), SPH_NP(simpar), STAR_NP(simpar), AGN_NP(simpar) );
	rvfile = fopen(infile,"r");

	/*
	Basic_MPI mpi_backup = *PTR2MPIINFO(simpar);
	*simpar = read_head(rvfile);
	*PTR2MPIINFO(simpar) = mpi_backup;
	*/
	read_head(rvfile, simpar);

	fflush(stderr);
	TIMER_START(5);
	/*
	fread(simpar, sizeof(SimParameters), 1, rvfile);
	*/
	if(PTYPE(simpar) == TREETYPE){
		DM_TBP(simpar) = (treedmparticletype*)malloc(sizeof(treedmparticletype)*DM_NP(simpar));
		SPH_TBP(simpar) = (treesphparticletype*)malloc(sizeof(treesphparticletype)*SPH_NP(simpar));
		STAR_TBP(simpar) = (treestarparticletype*)malloc(sizeof(treestarparticletype)*STAR_NP(simpar));
		AGN_TBP(simpar) = (treeagnparticletype*)malloc(sizeof(treeagnparticletype)*AGN_NP(simpar));
		fread(DM_TBP(simpar), sizeof(treedmparticletype), DM_NP(simpar), rvfile);
		fread(SPH_TBP(simpar), sizeof(treesphparticletype), SPH_NP(simpar), rvfile);
		fread(STAR_TBP(simpar), sizeof(treestarparticletype), STAR_NP(simpar), rvfile);
		fread(AGN_TBP(simpar), sizeof(treeagnparticletype), AGN_NP(simpar), rvfile);
	}
	else if(PTYPE(simpar) == PMTYPE){
		DM_BP(simpar) = (dmparticletype*)malloc(sizeof(dmparticletype)*DM_NP(simpar));
		SPH_BP(simpar) = (sphparticletype*)malloc(sizeof(sphparticletype)*SPH_NP(simpar));
		STAR_BP(simpar) = (starparticletype*)malloc(sizeof(starparticletype)*STAR_NP(simpar));
		AGN_BP(simpar) = (agnparticletype*)malloc(sizeof(agnparticletype)*AGN_NP(simpar));
		fread(DM_BP(simpar), sizeof(dmparticletype), DM_NP(simpar), rvfile);
		fread(SPH_BP(simpar), sizeof(sphparticletype), SPH_NP(simpar), rvfile);
		fread(STAR_BP(simpar), sizeof(starparticletype), STAR_NP(simpar), rvfile);
		fread(AGN_BP(simpar), sizeof(agnparticletype), AGN_NP(simpar), rvfile);
	}
	int iflag = fread(simpar, sizeof(SimParameters), 1, rvfile);
	fclose(rvfile);

	if(GROUPID(myid, WGroupSize) == GROUPID(tgt, WGroupSize) && tgt <nid)
		MPI_Send(&isend, 1, MPI_INT, tgt, itag, MPI_COMM(simpar));




	/* To build RMS communcation topology */
	buildrmscom(&DM_DDFUNC(simpar), DM_DDINFO(simpar), MPI_COMM(simpar));
	/* To do the domain decomposiiton */
	if(iflag ==0) {
		DomainDecomp(simpar, 1);
	}


#ifdef DEBUG
	float xmin,ymin,zmin, xmax,ymax, zmax;
	xmin = ymin = zmin = 1.e20;
	xmax = ymax = zmax = -1.e20; 
	if(PTYPE(simpar) == TREETYPE){
		ptrdiff_t j;
		for(j=0;j<DM_NP(simpar);j++){
			xmin = MIN(xmin, XofP(simpar,DM_TBP(simpar)+j));
			ymin = MIN(ymin, YofP(simpar,DM_TBP(simpar)+j));
			zmin = MIN(zmin, ZofP(simpar,DM_TBP(simpar)+j));
			xmax = MAX(xmax, XofP(simpar,DM_TBP(simpar)+j));
			ymax = MAX(ymax, YofP(simpar,DM_TBP(simpar)+j));
			zmax = MAX(zmax, ZofP(simpar,DM_TBP(simpar)+j));
		}
	}
	else if(PTYPE(simpar) == PMTYPE){
		ptrdiff_t j;
		for(j=0;j<DM_NP(simpar);j++){
			xmin = MIN(xmin, XofP(simpar,DM_BP(simpar)+j));
			ymin = MIN(ymin, YofP(simpar,DM_BP(simpar)+j));
			zmin = MIN(zmin, ZofP(simpar,DM_BP(simpar)+j));
			xmax = MAX(xmax, XofP(simpar,DM_BP(simpar)+j));
			ymax = MAX(ymax, YofP(simpar,DM_BP(simpar)+j));
			zmax = MAX(zmax, ZofP(simpar,DM_BP(simpar)+j));
		}
	}
	TIMER_STOP(5);
	fprintf(stderr," -- P%d finishes reading data %ld/ %ld/ %ld/ %ld in %g sec. with xyz= %g %g : %g %g : %g %g\n",
			myid, DM_NP(simpar), SPH_NP(simpar), STAR_NP(simpar), AGN_NP(simpar), ELAPSED_TIME(5),
			xmin, xmax, ymin,ymax, zmin, zmax);
#endif

}

void jwrite(SimParameters *simpar){
	FILE *rvfile;
	int myid, nid;
	MPI_Status status;
	char outfile[190];
	float cputime0[10];
	float cputime1[10];


	sprintf(RV_FILE(simpar),"%s.%.5d", RV_FILEPREFIX(simpar),STEPCOUNT(simpar));

	myid = MYID(simpar);
	nid = NID(simpar);


	int src, tgt, WGroupSize, isend, iget, itag = 1;
	WGroupSize = WGROUPSIZE(simpar);
	src = myid -1;
	tgt = myid +1;


	if(RANKINGROUP(myid, WGroupSize) !=0) 
		MPI_Recv(&iget, 1, MPI_INT, src, itag, MPI_COMM(simpar), &status);

#ifdef GOTPM
	sprintf(outfile,"%s.%.5d", RV_FILE(simpar),MYID(simpar));
#else
	sprintf(outfile,"%s.%.5d.%.5d", RV_FILE(simpar),IndT_ISUBSTEP(simpar), MYID(simpar));
#endif

	rvfile = fopen(outfile,"w");
	write_head(rvfile, simpar);
	fprintf(stderr,"P%d starts writing data: %s : sum= %ld : %ld, %ld, %ld, %ld  \n ",myid,
			outfile, NPSUM(simpar), DM_NP(simpar), SPH_NP(simpar), STAR_NP(simpar), AGN_NP(simpar) );
	fflush(stderr);
	TIMER_START(5);
	/*
	fwrite(simpar, sizeof(SimParameters), 1, rvfile);
	*/
	if(PTYPE(simpar) == TREETYPE){
		fwrite(DM_TBP(simpar), sizeof(treedmparticletype), DM_NP(simpar), rvfile);
		fwrite(SPH_TBP(simpar), sizeof(treesphparticletype), SPH_NP(simpar), rvfile);
		fwrite(STAR_TBP(simpar), sizeof(treestarparticletype), STAR_NP(simpar), rvfile);
		fwrite(AGN_TBP(simpar), sizeof(treeagnparticletype), AGN_NP(simpar), rvfile);
	}
	else  if(PTYPE(simpar) == PMTYPE){
		fwrite(DM_BP(simpar), sizeof(dmparticletype), DM_NP(simpar), rvfile);
		fwrite(SPH_BP(simpar), sizeof(sphparticletype), SPH_NP(simpar), rvfile);
		fwrite(STAR_BP(simpar), sizeof(starparticletype), STAR_NP(simpar), rvfile);
		fwrite(AGN_BP(simpar), sizeof(agnparticletype), AGN_NP(simpar), rvfile);
	}
	fwrite(simpar, sizeof(SimParameters), 1, rvfile);
	fclose(rvfile);
	TIMER_STOP(5);

	if(GROUPID(myid, WGroupSize) == GROUPID(tgt, WGroupSize) && tgt <nid)
		MPI_Send(&isend, 1, MPI_INT, tgt, itag, MPI_COMM(simpar));

	if(myid==0){
		FILE *paramfile;
		char paramname[190];
#ifdef GOTPM
		sprintf(paramname,"params.%.5d", STEPCOUNT(simpar));
#else
		sprintf(paramname,"params.%.5d_%.5d", STEPCOUNT(simpar), IndT_ISUBSTEP(simpar));
#endif
		paramfile = fopen(paramname,"w");
		{
			void write_sim_parameter_file(FILE*, SimParameters *);
			write_sim_parameter_file(paramfile, simpar);
		}
		fclose(paramfile);
	}

	fprintf(stderr," -- P%d finishes writing data %ld/ %ld/ %ld/ %ld in %g sec.\n",
			myid, DM_NP(simpar), SPH_NP(simpar), STAR_NP(simpar), AGN_NP(simpar), ELAPSED_TIME(5));
#ifdef DEBUG
	MPI_Barrier(MPI_COMM(simpar));
	DEBUGPRINT("P%d Now Ending the dump\n",MYID(simpar));
#endif

}
