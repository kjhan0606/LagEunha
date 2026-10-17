#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "cosmology.h"
#include "indT.h"
#include "flags.h"
int flagpsmeasure(SimParameters *simpar){
	float amax; float anow; float astep;int step;
	amax = AMAX(simpar);
	anow = ANOW(simpar);
#ifdef GOTPM
	astep = ASTEP(simpar);
#else
	astep = IndT_DAMIN(simpar);
#endif
	step = STEPCOUNT(simpar);
	FILE *fp;
	int saveflag;
	double red,redn;
	float redi;
	double redb,red1,red2;
	red = amax/(double)(anow)-1L;
	redn = amax/(double)(anow+astep)-1L;
	redb = amax/(double)(anow-astep)-1L;
	red1 = 0.5L*(red+redn); 
	red2 = 0.5L*(red+redb);
	if(MYID(simpar)==0){
		if((fp=fopen("WriteSync&WholeDen.flag","r"))){
			while(fscanf(fp,"%g",&redi)!= EOF){
				if( redi >=red1 && redi<red2){
					fclose(fp);
					saveflag = 1;
					goto out;
				}
			}
			fclose(fp);
		}
	}
	saveflag = 0;
out: 
	MPI_Bcast(&saveflag,1,MPI_INT,0,MPI_COMM(simpar));
	if(step%20==1 || step<10 || fabs(amax-anow) < 1.e-4) saveflag = 1;
	if(saveflag ==0) CONT_FLAGPS(simpar) = 'N';
	else CONT_FLAGPS(simpar) = 'Y';
	return saveflag;
}
int flagPreFoF(SimParameters *simpar){
	float amax; float anow; float astep;
	amax = AMAX(simpar);
	anow = ANOW(simpar);
#ifdef GOTPM
	astep = ASTEP(simpar);
#else
	astep = IndT_DAMIN(simpar);
#endif
	FILE *fp;
	int istep;
	int saveflag;
	double red,redn;
	float redi;
	double redb,red1,red2;
	red = amax/(double)(anow)-1L;
	redn = amax/(double)(anow+astep)-1L;
	redb = amax/(double)(anow-astep)-1L;
	red1 = 0.5L*(red+redn);
	red2 = 0.5L*(red+redb);
	if(MYID(simpar)==0){
		if((fp=fopen("PreFoF.flag","r"))){
			while(fscanf(fp,"%g",&redi)!= EOF){
				if( redi >=red1 && redi<red2){
					fclose(fp);
					saveflag= 1;
					DEBUGPRINT("Turning on the PreFoF with %g %g ::: %g<= %g<= %g\n",anow, astep, red1, redi, red2);
					goto out;
				}
			}
			fclose(fp);
		}
	}
	saveflag = 0;
out: 
	MPI_Bcast(&saveflag,1,MPI_INT,0,MPI_COMM(simpar));
	if(saveflag ==0) CONT_FLAGPREFOF(simpar) = 'N';
	else CONT_FLAGPREFOF(simpar) = 'Y';
	return saveflag;
}
int flagsyncpdata(SimParameters *simpar){
	float amax; float anow; float astep;
	amax = AMAX(simpar);
	anow = ANOW(simpar);
#ifdef GOTPM
	astep = ASTEP(simpar);
#else
	astep = IndT_DAMIN(simpar);
#endif
	FILE *fp;
	int istep;
	int saveflag;
	double red,redn;
	float redi;
	double redb,red1,red2;
	red = amax/(double)(anow)-1L;
	redn = amax/(double)(anow+astep)-1L;
	redb = amax/(double)(anow-astep)-1L;
	red1 = 0.5L*(red+redn);
	red2 = 0.5L*(red+redb);

	if(MYID(simpar)==0){
		if((fp=fopen("WriteSyncPData.flag","r"))){
			while(fscanf(fp,"%g",&redi)!= EOF){
				if( redi >=red1 && redi<red2){
					fclose(fp);
					DEBUGPRINT("Turning on the SyncWrite with %g %g ::: %g<= %g<= %g\n",anow, astep, red1, redi, red2);
					saveflag= 1;
					goto out;
				}
			}
			fclose(fp);
		}
	}
	saveflag = 0;
out: 
	MPI_Bcast(&saveflag,1,MPI_INT,0,MPI_COMM(simpar));
	if(saveflag ==0) CONT_FLAGSYNCP(simpar) = 'N';
	else CONT_FLAGSYNCP(simpar) = 'Y';
	return saveflag;
}
int flagwholeden(SimParameters *simpar){
	float amax; float anow; float astep;
	amax = AMAX(simpar);
	anow = ANOW(simpar);
#ifdef GOTPM
	astep = ASTEP(simpar);
#else
	astep = IndT_DAMIN(simpar);
#endif
	FILE *fp;
	int istep;
	int saveflag;
	double red,redn;
	float redi;
	double redb,red1,red2;
    red = amax/(double)(anow)-1L;
	redn = amax/(double)(anow+astep)-1L;
	redb = amax/(double)(anow-astep)-1L;
	red1 = 0.5L*(red+redn);
	red2 = 0.5L*(red+redb);

	if(MYID(simpar)==0){
		if((fp=fopen("WriteWholeDen.flag","r"))){
			while(fscanf(fp,"%g",&redi)!= EOF){
				if( redi >=red1 && redi<red2){
					fclose(fp);
					saveflag= 1;
					goto out;
				}
			}
			fclose(fp);
		}
	}
	saveflag = 0;
out: 
	MPI_Bcast(&saveflag,1,MPI_INT,0,MPI_COMM(simpar));
	if(saveflag ==0) CONT_FLAGWHOLEDEN(simpar) = 'N';
	else CONT_FLAGWHOLEDEN(simpar) = 'Y';
	return saveflag;
}

int flagsuddenstopInd(SimParameters *simpar){
	int nowstep; int nowsubtime;
	nowstep = STEPCOUNT(simpar);
	nowsubtime = IndT_NSUBSTEP(simpar);
	int stopstep,stopsubtime;
	int flag;
	FILE *fp;
	char flagcontinue[100]; 
	flag= 0;
	if(MYID(simpar)==0){
		if((fp=fopen("SuddenstopInd.flag","r"))){
			while(fscanf(fp,"%d",&stopstep)!=EOF){
				fscanf(fp,"%d",&stopsubtime);
				if(stopstep == nowstep && stopsubtime == nowsubtime) {
					if(fscanf(fp,"%s",flagcontinue)!=EOF){
						if(strcmp(flagcontinue,"+")==0) flag = 1;
						else if(strcmp(flagcontinue,"-")==0) flag = 2;
						else flag = 3;
					}
					else flag= 0;
					goto out;
				}
				else {
					flag= 4;
					fscanf(fp,"%s",flagcontinue);
				}
			}
			fclose(fp);
		}
	}
out:
	MPI_Bcast(&flag,1,MPI_INT,0,MPI_COMM(simpar));
	if(flag == 0) CONT_FLAGCONTINUE(simpar) = 'P';
	else if(flag == 1) CONT_FLAGCONTINUE(simpar) = 'D';
	else if(flag == 2) CONT_FLAGCONTINUE(simpar) = 'S';
	else if(flag == 3) CONT_FLAGCONTINUE(simpar) = 'E';
	else CONT_FLAGCONTINUE(simpar) = 'P';
	return flag;
}


int flagsuddenstopGOTPM(SimParameters *simpar){
	int nowstep; int nowsubtime;
	nowstep = STEPCOUNT(simpar);
	nowsubtime = IndT_NSUBSTEP(simpar);
	int stopstep,stopsubtime;
	int flag;
	FILE *fp;
	char flagcontinue[100]; 
	flag= 0;
	if(MYID(simpar)==0){
		if((fp=fopen("Suddenstop.flag","r"))){
			while(fscanf(fp,"%d",&stopstep)!=EOF){
				if(stopstep == nowstep){
					if(fscanf(fp,"%s",flagcontinue)!=EOF){
						if(strcmp(flagcontinue,"+")==0) flag = 1;
						else if(strcmp(flagcontinue,"-")==0) flag = 2;
						else flag = 3;
					}
					else flag= 0;
					goto out;
				}
				else {
					flag= 4;
					fscanf(fp,"%s",flagcontinue);
				}
			}
			fclose(fp);
		}
	}
out:
	MPI_Bcast(&flag,1,MPI_INT,0,MPI_COMM(simpar));
	if(flag == 0) CONT_FLAGCONTINUE(simpar) = 'P';
	else if(flag == 1) CONT_FLAGCONTINUE(simpar) = 'D';
	else if(flag == 2) CONT_FLAGCONTINUE(simpar) = 'S';
	else if(flag == 3) CONT_FLAGCONTINUE(simpar) = 'E';
	else CONT_FLAGCONTINUE(simpar) = 'P';
	return flag;
}

