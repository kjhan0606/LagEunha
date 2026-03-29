#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>

#include "eunha.h"
#include "params.h"

int read_viewers(char *paramfile,SimParameters *simpar){
	int isize,ncnt=0;
	int nviewer, mview = 0;
	char line[MAX_LINE_LENGTH];
	FILE *fp;
	fp  = fopen(paramfile,"r");
   	while(strcmp(P_Closing,fgets(line, MAX_LINE_LENGTH,fp)) != 0){ /* Loop until it reaches P_Closing */ 
		int ncnt = 0;
		if(line[0] != '#'){
			if(strstr(line,"define") != NULL ) {
				VIEWER2(sscanf,line, &, simpar,ncnt);
				if(ncnt == 1) {
					nviewer = ANIM_NVIEWER(simpar);
					break;
				}
			}
		}
	}
	fclose(fp);
	if(mview != nviewer){
		fprintf(stderr,"Strange in the number of viewers\n");
	}
	return nviewer;
}


void getsimparm_(char *paramfile,float *boxsize,float *hubble,float *npower, 
		float *omep, float *omepb, float*omeplam, float *wlam,float *bias, float *rng, 
		float *amax, float *da, float *anow, int *iflag, char *outfile){ 
	SimParameters rsimpar; 
	SimParameters *simpar = &rsimpar; 
	char cpar[190];
	FILE *fp; 
	sprintf(cpar,"%s",paramfile);
	fp = fopen(cpar,"r");
   	rsimpar=read_sim_parameter_file(rsimpar, fp);
   	*boxsize = BOXSIZE(simpar);
   	*hubble = HUBBLE(simpar);
   	*npower = NPOW(simpar);
   	*omep = OMEP(simpar);
   	*omepb = OMEPB(simpar);
   	*omeplam = OMEPLAM(simpar);
   	*bias = BIAS8(simpar);
   	*rng = NX(simpar);
   	*amax = AMAX(simpar);
   	*da = ASTEP(simpar);
   	*anow = ANOW(simpar);
   	*iflag = POWREADFLAG(simpar);
	*wlam = WLAM0(simpar);
   	sprintf(outfile,"%s", POWFILENAME(simpar));
   	fclose(fp);
}

