#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include "eunha.h"
#include "params.h"
#include "voro.h"
#include "kh.h"


int main(int argc, char **argv){


	int i,j,k;

	int nstep;
	char fheader[190], fname[190];
	SimParameters simpar;
	int icont;
	postype time;

	nstep = atoi(argv[2]);


	FILE *fp = fopen(argv[1],"r");
	/* Make Default SimParameter for Cosmological Simulation */ 
	mk_default_param(&simpar, "WMAP5");
	/* Read Simulation Input Parameters */ 
	read_head(fp,  &simpar);
	fclose(fp);

	if(SIMMODEL(&simpar)==KH){
		sprintf(fheader,"khout");
	}
	else if(SIMMODEL(&simpar)==RT){
		sprintf(fheader,"rtout");
	}

	sprintf(fname,"%s.%.6d.dat", fheader,nstep);

	fp = fopen(fname,"r");

	int np;
	fread(&np, sizeof(int), 1, fp); 
	treevorork4particletype *bp = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*(np));

	fread(bp, sizeof(treevorork4particletype), np, fp);
	fread(&time, sizeof(postype), 1, fp);
	fclose(fp);

	double ke, pe, te, ie;
	ke = pe = ie = te = 0;

	for(i=0;i<np;i++){
		ke += 0.5*bp[i].mass*(bp[i].x*bp[i].x + bp[i].y*bp[i].y);
		ie += bp[i].ie;
	}
	if(SIMMODEL(&simpar) == RT){
		double acc = RT_ACC(&simpar);
		for(i=0;i<np;i++){
			pe += bp[i].mass*acc*bp[i].y;
		}
	}
	te = ke + pe + ie;

	printf("t= %g ke= %g pe= %g ie= %g te= %g\n", time, ke, pe, ie,te);
	free(bp);

}
