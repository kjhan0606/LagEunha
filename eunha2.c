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
		DEBUGPRINT("Error in different size of ptrdiff_t & long:  %zu : %zu\n",sizeof(ptrdiff_t), sizeof(long));
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

void DEBUGGING(SimParameters *simpar){
	char hostname[190];
	gethostname(hostname,190);
	if(MYID(simpar)==0) {
		printf("###############################################\n");
		printf("EUNHA2: compiled at %s %s \n",__TIME__, __DATE__);
		printf("###############################################\n");
	}
	DEBUGPRINT("P%d is on %s with pid=%d\n",MYID(simpar),hostname,getpid());
	{
		long kkk = 1;
		while(kkk) {
			kkk = 1;
		}
	}
}
// Global counter for malloc calls
static int malloc_calls = 0;

    // Wrapper for malloc
void* my_malloc(size_t size) {
    malloc_calls++; // Increment the counter
	return malloc(size); // Call the original malloc
}
void* my_realloc(void *a, size_t size) {
	if(size==0) {
		my_free(a);
		return NULL;
	}
	else {
		return realloc(a, size); // Call the original realloc
	}
}

// Function to get the count
int get_malloc_call_count() {
	return malloc_calls;
}
void my_free(void *a) {
	--malloc_calls;
	free(a);
}


int FUNC(MAIN)(int argc, char **argv)
{
	int icont;
	SimParameters simpar;
	MPI_Comm com = MPI_COMM_WORLD;

	for(int i = 0; i < argc; i++) if(argv[i] == NULL) argv[i] = "";

#ifdef LAM_MPI
	printf("%d\n",MPI_Init(NULL,NULL));
#else
	mpi_fftw_initialize(argc, argv);
#endif
	checkarg(argc);

	FILE *fp = fopen(argv[1],"r");



	/* Make Default SimParameter for Cosmological Simulation */
	mk_default_param(&simpar, "WMAP5");
	fprintf(stderr,"[DBG] after mk_default_param simmodel=%d NWGROUP=%d\n",
		SIMMODEL((&simpar)), NWGROUP((&simpar))); fflush(stderr);

	/* Initial Setting of MPI Communicator */
	Mpi_Basic_Set(&simpar, com);
	fprintf(stderr,"[DBG] after Mpi_Basic_Set myid=%d\n", MYID((&simpar))); fflush(stderr);

	HAMB(&simpar);
	fprintf(stderr,"[DBG] after HAMB\n"); fflush(stderr);


	/* Read Simulation Input Parameters */
	ReadSimulationParameters(fp, &icont, &simpar);
	fprintf(stderr,"[DBG] after ReadSimulationParameters simmodel=%d\n", SIMMODEL((&simpar))); fflush(stderr);

	if(SIMMODEL( (&simpar) ) == Cosmos) {
		RunCosmos(&simpar, icont); // icont from ReadSimulationParameters is an input here.
	}
	/*
	else if(SIMMODEL( (&simpar) ) == KH_LF) {
		if(argc ==3) icont = atoi(argv[2]);
//		int RunKHOLD(SimParameters *, int);
//		RunKHOLD(&simpar, icont);
		int RunKH(SimParameters *, int, int);
		RunKH(&simpar, icont, 0);
	}
	*/
	else if(SIMMODEL( (&simpar) ) == KH) {
		if(argc ==3) icont = atoi(argv[2]);
		int RunKH(SimParameters *, int);
		RunKH(&simpar, icont);
	}
	else if(SIMMODEL( (&simpar) ) == RT) {
		if(argc ==3) icont = atoi(argv[2]);
		int RunRT(SimParameters *, int);
		RunRT(&simpar, icont);
	}
	else if(SIMMODEL( (&simpar) ) == RT_LF) {
		if(argc ==3) icont = atoi(argv[2]);
		int RunRT_LF(SimParameters *, int);
		RunRT_LF(&simpar, icont);
	}
	/*
	else if(SIMMODEL( (&simpar) ) == Kepler) {
		if(argc ==3) icont = atoi(argv[2]);
		int RunKepler(SimParameters *, int);
		RunKepler(&simpar, icont);
	}
	*/
	/*
	else if(SIMMODEL( (&simpar) ) == MkGlass2D) {
		if(argc ==3) icont = atoi(argv[2]);
		int Make2DGlass(SimParameters *, int);
		Make2DGlass(&simpar, icont);
	}
	*/
	else if(SIMMODEL( (&simpar) ) == MkGlass2D) {
		if(argc ==3) icont = atoi(argv[2]);
		int Make2DGlass(SimParameters *, int, int);
		Make2DGlass(&simpar, icont,1);
	}
	else if(SIMMODEL( (&simpar) ) == Cylinder) {
		if(argc ==3) icont = atoi(argv[2]);
		int RunCylinder(SimParameters *, int);
		RunCylinder(&simpar, icont);
	}
	else if(SIMMODEL( (&simpar) ) == Sedov2D) {
		if(argc ==3) icont = atoi(argv[2]);
		int RunSedov2D(SimParameters *, int);
		RunSedov2D(&simpar, icont);
	}

	MPI_Finalize();

	
	return 0;
}
