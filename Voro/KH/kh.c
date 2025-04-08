#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "eunha.sub.h"
#include "kjhrw.h"
#include "parallelIO.h"


int KH_2D(SimParameters *simpar, int icont){

	if(MYID(simpar)==0){
		printf("Kelvin Helmholtz Simulation in 2D is Starting\n");
	}

	SetEnvGridInfo(simpar);
	StartRkSDD(simpar);
}
