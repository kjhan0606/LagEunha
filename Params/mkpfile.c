#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<math.h>

#include "eunha.h"
#include "params.h"


int main(int argc, char **argv){
	int i,j,k;
	FILE *wp;
	SimParameters sim;


	if(argc !=3){
		fprintf(stderr,"Please input as mkpfile.exe [outputfilename] [cosmology:WMAP3/WMAP5/NoExpand/Blast/KH/BowShock]\n");
		exit(99);
	}
	wp = fopen(argv[1],"w");
	GAS_SPHFLAG((&sim)) = 'Y';
	{
		void mk_default_param(SimParameters *,char *);
		if(strcmp(argv[2],"WMAP3")==0 || strcmp(argv[2],"WMAP5") ==0){
			printf("With Hydro or not? 1/0\n");
			char yesno=fgetc(stdin);
			if(yesno=='0') GAS_SPHFLAG((&sim)) = 'N';
		}
		mk_default_param(&sim,argv[2]);
	}
	{
		void write_default_sim_parameter_file (FILE *, SimParameters);
		write_default_sim_parameter_file(wp,sim);
	}
	fclose(wp);


}
