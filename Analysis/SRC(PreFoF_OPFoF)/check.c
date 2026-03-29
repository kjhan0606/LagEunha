#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "pmheader.h"
#include "params.h"


#define MIN(a,b) (a)<(b)? (a):(b)
#define MAX(a,b) (a)>(b)? (a):(b)


SimParameters simpar;

pmparticletype *bp;

int main(){
	int istep;
	long long mp;


	for(istep=0;istep<1000;istep++){ 
		char infile[190]; 
		sprintf(infile,"PreFoF.01328%.5d",istep); 
		FILE *fp = fopen(infile,"r"); 
		read_head(fp, &simpar); 
		bp = (pmparticletype*)malloc(sizeof(pmparticletype)*simpar.np); 
		mp = fread(bp,sizeof(pmparticletype),simpar.np,fp); 
		fclose(fp); 
		double zmin,zmax; 
		zmin = 1.e12; 
		zmax = -1e12; 
		long long i; 
		for(i=0;i<simpar.np;i++){ 
			zmin = MIN(zmin,ZofP(bp+i)); 
			zmax = MAX(zmax,ZofP(bp+i));
			/*
			if(zmin < 10){
				printf("%ld\n",i);
			}
			*/
		}
		printf("%s : %g %g\n",infile,zmin,zmax);
		free(bp);
	}

}
