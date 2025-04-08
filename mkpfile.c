#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<math.h>

#include "eunha.h"
#include "params.h"


int main (int argc, char **argv){
	int i,j,k;
	FILE *wp;
	SimParameters sim;


	if(argc !=3){
		DEBUGPRINT("Please input as mkpfile.exe [outputfilename] [cosmology:(B)WMAP3/(B)WMAP5/NoExpand/Blast/KH/BowShock/RT/Kepler/MkGlass2D]\n now narg= %d\n", argc);
		exit(99);
	}
	wp = fopen(argv[1],"w");
	{
		void mk_default_param(SimParameters *,char *);
		mk_default_param(&sim,argv[2]);
	}
	{
		void write_default_sim_parameter_file (FILE *, SimParameters *);
		write_default_sim_parameter_file(wp,&sim);
	}
	fclose(wp);


}
