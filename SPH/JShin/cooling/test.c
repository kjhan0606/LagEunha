#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>



int MAIN_(int argc, char **argv){
	void readcoolingdata_();
	readcoolingdata_();

	void getcoolheatredshift_(double*);
	double redshift = 128;
	getcoolheatredshift_(&redshift);

	void getcoolheatfinal_(double*,double*,double*,int*,double*,double*);
	double gastemp,gasdens,gasmetal,coolrate,heatrate;
	int Uvshield = 0;

	gastemp = 259.0000000000000;
	gasmetal = 	1.9999999552965164E-002;
	gasdens = 9.9999996826552220E-021;
	getcoolheatfinal_(&gastemp,&gasdens,&gasmetal,&Uvshield, &coolrate,&heatrate);

	printf("%g %g\n",coolrate,heatrate);

}
