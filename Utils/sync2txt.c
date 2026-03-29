#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<mpi.h>
#include<math.h>

#include "eunha.h"
#include "fft.h"
#include "mpiaux.h"
#include "cosmology.h"
#include "params.h"



int main(int argc, char **argv){
	int nstep;
	SimParameters simpar;
	FILE *fp;
	int nid,myid;
	MPI_Comm com;


	mpi_fftw_initialize(argc, argv);

	/* Make Default SimParameter for Cosmological Simulation */
    mk_default_param(&simpar, "WMAP5");


	com = MPI_COMM_WORLD;
	Mpi_Basic_Set(&simpar, com);

	char infile[190];

	sprintf(infile,"%s.%.5d",argv[1],MYID(&simpar));

	fp = fopen(infile,"r");
	int icont = 0;

	printf("P%d has %p :: %s\n", MYID(&simpar), fp, infile);

	ReadSimulationParameters(fp, &icont, &simpar);

	InitializeReadStart(&simpar);

	void jread(SimParameters*);
	jread(&simpar);
	int i,j;
	float ratio,vratio;
	float Hsub;

	Hsub = sqrt(OMEP(&simpar)*pow(AMAX(&simpar)/ANOW(&simpar),3) + OMEPLAM(&simpar)
			+(1.-OMEP(&simpar)-OMEPLAM(&simpar))*pow(AMAX(&simpar)/ANOW(&simpar),2));
	ratio = BOXSIZE(&simpar)/NX(&simpar);
	vratio = BOXSIZE(&simpar)*ANOW(&simpar)*ANOW(&simpar)/AMAX(&simpar)*100.*Hsub;

	for(i=0;i<NID(&simpar);i++){
		FILE *wp;
		if(i == MYID(&simpar)){
			if(i==0) wp= fopen(argv[2],"w");
			else wp = fopen(argv[2],"a");
			dmparticletype *bp = DM_BP(&simpar);
			for(j=0;j<DM_NP(&simpar);j++){
				fprintf(wp,"%g %g %g %g %g %g\n", 
						ratio*XofP(&simpar,bp+j),ratio*YofP(&simpar,bp+j),ratio*ZofP(&simpar,bp+j),
						(bp+j)->vx*vratio, (bp+j)->vy*vratio,(bp+j)->vz*vratio);
			}
			fclose(wp);
		}
		MPI_Barrier(MPI_COMM(&simpar));
	}

	MPI_Finalize();
}
