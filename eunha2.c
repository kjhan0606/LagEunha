#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<mpi.h>
#include<unistd.h>

#include "eunha.h"
#include "fft.h"
#include "mpiaux.h"
#include "cosmology.h"
#include "params.h"

#ifdef PGCC
#define FUNC(MAIN) (MAIN_)
#elif INTEL
#define FUNC(MAIN) (main)
#elif CRAY
#define FUNC(MAIN) (main)
#else
#error wrong compiler setting
#endif

void checkarg(int argc){
	if(argc !=2 && argc !=3){
		DEBUGPRINT0("Error in the argument\n");
		DEBUGPRINT0("usage: eunha.exe paramsfile\n");
		MPI_Finalize();
		exit(100);
	}
	if(sizeof(ptrdiff_t) != sizeof(long)){
		DEBUGPRINT("Error in different size of ptrdiff_t & long:  %ld : %ld\n",sizeof(ptrdiff_t), sizeof(long));
		MPI_Finalize();
		exit(999);
	}
}
void HAMB(SimParameters *simpar){
	char hostname[190];
	gethostname(hostname,190);
	if(MYID(simpar)==0) {
		printf("###############################################\n");
		printf("EUNHA2: compiled at %s %s \n",__TIME__, __DATE__);
		printf("###############################################\n");
	}
	DEBUGPRINT("P%d is on %s with pid=%d\n",MYID(simpar),hostname,getpid());
	if(MYID(simpar) == -1){
		long kkk = 1;
		while(kkk) {
			kkk = 1;
		}
	}
}

int FUNC(MAIN)(int argc, char **argv)
{
	int icont;
	SimParameters simpar;
	MPI_Comm com = MPI_COMM_WORLD;
#ifdef LAM_MPI
	printf("%d\n",MPI_Init(NULL,NULL));
#else
	mpi_fftw_initialize(argc, argv);
#endif
	checkarg(argc);

	FILE *fp = fopen(argv[1],"r");



	/* Make Default SimParameter for Cosmological Simulation */
	mk_default_param(&simpar, "WMAP5");

	/* Initial Setting of MPI Communicator */
	Mpi_Basic_Set(&simpar, com);

	HAMB(&simpar);


	/* Read Simulation Input Parameters */
	ReadSimulationParameters(fp, &icont, &simpar);

	if(SIMMODEL( (&simpar) ) == Cosmos) {
		RunCosmos(&simpar, icont); // icont from ReadSimulationParameters is an input here.
	}
	else if(SIMMODEL( (&simpar) ) == KH) {
		if(argc ==3) icont = atoi(argv[2]);
		int RunKH(SimParameters *, int);
		RunKH(&simpar, icont);
	}
	/*
	else if(SIMMODEL( (&simpar) ) == RT) {
		if(argc ==3) icont = atoi(argv[2]);
		int RunRT(SimParameters *, int);
		RunRT(&simpar, icont);
	}
	else if(SIMMODEL( (&simpar) ) == Kepler) {
		if(argc ==3) icont = atoi(argv[2]);
		int RunKepler(SimParameters *, int);
		RunKepler(&simpar, icont);
	}
	*/
	else if(SIMMODEL( (&simpar) ) == MkGlass2D) {
		if(argc ==3) icont = atoi(argv[2]);
		int Make2DGlass(SimParameters *, int);
		Make2DGlass(&simpar, icont);
	}

	MPI_Finalize();

	
	return 0;
}
