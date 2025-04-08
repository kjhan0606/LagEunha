#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include<sys/times.h>
#include "eunha.h"
#include "nnost.h"
#include "hydroBasicCell.h"
#include "indT.h"
#include "timerutil.h"
#include "sph.h"
#include "CosmosEvolFactor.h"
#include "recfast.h"

//jhshin1
#include "mu.h"
//jhshin2

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

float cputime1[1], cputime0[1];






//jhshin1
//jhshin2
void FindSphDensity(SimParameters *);
void CompleteSph(SimParameters *);

float simparqsimp(float (*func)(SimParameters *, float), float , float , SimParameters *);

double GetTime(SimParameters *, float ,float );
float HofEz(SimParameters *, float );

static int istellarmodel=1;

void DetermineSphFactor(SimParameters *simpar){
	float zp1,ellb,gamma,alpha;
	float Hz;
	double red;
	zp1 = AMAX(simpar)/ANOW(simpar);
	red = zp1 - 1;
	ellb = BOXSIZE(simpar)/NX(simpar);
	gamma = GAS_GAMMA(simpar);
	alpha = GAS_ALPHA(simpar);
	GAS_VISFORCEFACTOR(simpar) = 4*M_PI*alpha*gamma/
		(3*OMEP(simpar)*H0*H0*ellb*ellb*zp1) *kBovermH / (NX(simpar)*NX(simpar));
	GAS_SPHFORCEFACTOR(simpar) = 8*M_PI/
		(3*OMEP(simpar)*H0*H0*ellb*ellb*zp1) *kBovermH / (NX(simpar)*NX(simpar));
	GAS_TEMPEVOLFACTOR(simpar) = -0.25*alpha*gamma*(gamma-1);
	Hz = H0*HofEz(simpar, red);
	GAS_G1(simpar)= Hz*ellb/zp1;
	GAS_G2(simpar) =GAS_G1(simpar)*ANOW(simpar)*NX(simpar);
	GAS_G12RATIO(simpar) = ANOW(simpar)*NX(simpar);
	GAS_G1EXPANSION(simpar)= 1;
	{
		float redshift, Hzzp1(SimParameters *, float);
		double H0incgs;
		redshift = red;
		H0incgs = H0/Mpccgs*HUBBLE(simpar);
		LOOKBACKTIME(simpar) = simparqsimp(Hzzp1,0,redshift,simpar)/(H0incgs*Myr); /* real time */ 
	}
	{
		double pntmass;
		pntmass = 3.L/8.L/M_PI/P_G*pow(100.L*100000.L*HUBBLE(simpar),2)*Mpccgs; 
		pntmass = pntmass*pow(BOXSIZE(simpar)/NX(simpar)*NSPACE(simpar)/HUBBLE(simpar),3)*
			OMEP(simpar);
		SPH_MASSINSOLARMASS(simpar) = pntmass/onesolarmass;
	}

	if(0){
		double nHnow = 3*H0/Mpccgs*H0/Mpccgs/(8.*M_PI*P_G)*OMEPB(simpar)*(1-GAS_YP(simpar))/mH;
		double nHenow = nHnow/not4*GAS_YP(simpar)/(1-GAS_YP(simpar));
		double nHpHenow = nHnow + nHenow;
		GAS_SIMDEN2HNDEN(simpar) = OMEP(simpar)/OMEPB(simpar) *zp1*zp1*zp1 * nHnow;
		GAS_SIMDEN2HPHENDEN(simpar) = OMEP(simpar)/OMEPB(simpar) *zp1*zp1*zp1 * nHpHenow;
	}
	{
		double H0incgs;
		H0incgs  = H0/Mpccgs*HUBBLE(simpar);
		GAS_RHOC0(simpar) =  3.L/(8.L*M_PI*P_G)* H0incgs*H0incgs;
		GAS_MEANRHO(simpar) =  GAS_RHOC0(simpar)* OMEP(simpar) * zp1*zp1*zp1;
		GAS_RHOS2RHOR(simpar)= OMEP(simpar) * GAS_RHOC0(simpar) *zp1*zp1*zp1;
		/*
		GAS_ENTROPYFACT(simpar) = 1./pow(OMEP(simpar)*GAS_MEANRHO(simpar)c0*zp1*zp1*zp1,1.-simpar.sph.gamma);
		*/
		/* This should be inverse of the mean density in cgs */
		GAS_ENTROPYFACT(simpar) = 1./(OMEP(simpar)*GAS_RHOC0(simpar)*zp1*zp1*zp1);
		GAS_ENTROPYEVOLFACT(simpar) = 0.25*alpha*(gamma-1)*gamma;
		GAS_SFRHO577(simpar)= GAS_MEANRHO(simpar)/GAS_RHOS2RHOR(simpar)*GAS_SFVIRIALDEN(simpar);
	}
	{
		if(GAS_COOLFLAG(simpar) == 'Y') GAS_COOLFACT(simpar) = 1;
		else GAS_COOLFACT(simpar) = 0;
		if(GAS_BGHEATFLAG(simpar) == 'Y') GAS_HEATFACT(simpar)= 1;
		else GAS_HEATFACT(simpar)= 0;
	}
	GAS_MINENTROPYFACTOR(simpar) = kBovermH * GAS_MINTEMP(simpar)/1.23;

	if(istellarmodel){
		void readstellarmodel_(int *);
		if(MYID(simpar)==0){
			printf("Now start reading stellar model\n");
		}
		readstellarmodel_(&(MYID(simpar)));
		istellarmodel = 0;
	}
	if(MYID(simpar) ==0) {
		printf("Expanding medium Sph Factors %g %g\n",
				GAS_VISFORCEFACTOR(simpar), GAS_SPHFORCEFACTOR(simpar));
		fflush(stdout);
	}
}

void StaticDetermineSphFactor(SimParameters *simpar){
	float zp1,ellb,gamma,alpha;
	float Hz;
	double red=0;
	if(MYID(simpar) ==0) printf("Static Medium Sph Factors\n");
	zp1 = 1;
	ellb = BOXSIZE(simpar)/NX(simpar);
	gamma = GAS_GAMMA(simpar);
	alpha = GAS_ALPHA(simpar);

	GAS_VISFORCEFACTOR(simpar) = alpha/2.*FFTELLB(simpar) *kBovermH *GAS_GAMMA(simpar)/ (NX(simpar)*NX(simpar));
	GAS_SPHFORCEFACTOR(simpar) = FFTELLB(simpar) *kBovermH / (NX(simpar)*NX(simpar));
	GAS_ENTROPYEVOLFACT(simpar) = 0.25*GAS_ALPHA(simpar)*(GAS_GAMMA(simpar)-1)*GAS_GAMMA(simpar) *FFTIME(simpar)/ (ellb*Mpccgs);
	GAS_TEMPEVOLFACTOR(simpar) = -0.25*alpha*gamma*(gamma-1);

	GAS_G1(simpar) = 1;
	GAS_G2(simpar) = ellb*Mpccgs/FFTIME(simpar) * NX(simpar);
	GAS_G12RATIO(simpar)= GAS_G2(simpar);
	GAS_G1EXPANSION(simpar)= 0;
	{
		float redshift, Hzzp1(SimParameters *, float);
		/* lookbacktime in Myr */
		LOOKBACKTIME(simpar) = (AMAX(simpar) - ANOW(simpar))*FFTIME(simpar)/Myr;
		
	}
	{
		double pntmass;
		pntmass = TOTMASS(simpar);
        pntmass = pntmass/(NX(simpar)*NX(simpar)*NX(simpar)/NSPACE(simpar)/NSPACE(simpar)/NSPACE(simpar));
		SPH_MASSINSOLARMASS(simpar) = pntmass;
	}


	void InitializeMu(SimParameters *, double); InitializeMu(simpar, red);

	if(0){
		double nHnow = 3*H0/Mpccgs*H0/Mpccgs/(8.*M_PI*P_G)*OMEPB(simpar)*(1-GAS_YP(simpar))/mH;
		double nHenow = nHnow/not4*GAS_YP(simpar)/(1-GAS_YP(simpar));
		double nHpHenow = nHnow + nHenow;
		GAS_SIMDEN2HNDEN(simpar) = OMEP(simpar)/OMEPB(simpar) *zp1*zp1*zp1 * nHnow;
		GAS_SIMDEN2HPHENDEN(simpar) = OMEP(simpar)/OMEPB(simpar) *zp1*zp1*zp1 * nHpHenow;
	}
	{
		double meanrho,boxsizecgs;
		boxsizecgs = BOXSIZE(simpar)*Mpccgs;
		meanrho = TOTMASS(simpar)*onesolarmass/(pow(boxsizecgs,3.L));
		GAS_MEANRHO(simpar) =  meanrho;
		GAS_RHOS2RHOR(simpar) = meanrho;
		/* This should be inverse of the mean density in cgs */
		GAS_ENTROPYFACT(simpar) = 1./meanrho;
		GAS_SFRHO577(simpar) = GAS_MEANRHO(simpar)/GAS_RHOS2RHOR(simpar)*GAS_SFVIRIALDEN(simpar);
	}
	{
		if(GAS_COOLFLAG(simpar) == 'Y') GAS_COOLFACT(simpar) = 1;
		else GAS_COOLFACT(simpar) = 0;
		if(GAS_BGHEATFLAG(simpar) == 'Y') GAS_HEATFACT(simpar) = 1;
		else GAS_HEATFACT(simpar) = 0;
	}
	GAS_MINENTROPYFACTOR(simpar) = kBovermH * GAS_MINTEMP(simpar)/1.23;

	if(istellarmodel){
		void readstellarmodel_(int *);
		readstellarmodel_(&(MYID(simpar)));
		istellarmodel = 0;
	}
}


void HydroBuildLinkedList(SimParameters *simpar){
	PosType CellWidth = BASICCELL_CELLWIDTH(simpar);
	size_t mx,my,mz;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	mz = BASICCELL_MZ(simpar);
	PosType xmin,ymin,zmin;
	xmin = SIM_LXMIN(simpar,sph);
	ymin = SIM_LYMIN(simpar,sph);
	zmin = SIM_LZMIN(simpar,sph);

	HydroTreeLinkedCell *SPH_BasicCell = SPH_BASICCELL(simpar);
	HydroTreeLinkedCell *STAR_BasicCell = STAR_BASICCELL(simpar);
	HydroTreeLinkedCell *AGN_BasicCell = AGN_BASICCELL(simpar);

	ptrdiff_t i;
	for(i=0;i<mx*my*mz;i++){
		SPH_BasicCell[i].link = NULL;
		SPH_BasicCell[i].nmem = 0;
		SPH_BasicCell[i].calflag = 0;

		STAR_BasicCell[i].nmem = 0;
		STAR_BasicCell[i].link = NULL;
		STAR_BasicCell[i].calflag = 0;

		AGN_BasicCell[i].nmem = 0;
		AGN_BasicCell[i].link = NULL;
		AGN_BasicCell[i].calflag = 0;
	}
	LinkParticles(simpar, SPH, sph, SPH_BasicCell, mx,my,xmin,ymin,zmin,CellWidth);
	LinkParticles(simpar, STAR, star, STAR_BasicCell, mx,my,xmin,ymin,zmin,CellWidth);
	LinkParticles(simpar, AGN, agn, AGN_BasicCell, mx,my,xmin,ymin,zmin,CellWidth);
}

void DestroyHydroLinkedCell(SimParameters *simpar){
	free(SPH_BASICCELL(simpar));
	free(STAR_BASICCELL(simpar));
	free(AGN_BASICCELL(simpar));
}
/* 
 *  J. Shin's algorithm. Please check it.
 *
 *
 *
 *
 */

//jhshin1
void Init_dAs(SimParameters *simpar){
	long i;
#ifdef _OPENMP
#pragma omp parallel  private(i)
#endif
	for(i = 0;i<SPH_NP(simpar);i++){
		(SPH_TBP(simpar)+i)->dAs = 0;
		(SPH_TBP(simpar)+i)->delta_e= 0;
	}
}

float GetTempGivenMu(SimParameters *, treesphparticletype*);

void Calc_mutempnden(SimParameters *simpar){
 	long i;
	if(GAS_CONSTMU(simpar) !='Y'){
#ifdef _OPENMP
#pragma omp parallel  private(i)
#endif
		for(i=0;i<SPH_NP(simpar);i++){
			treesphparticletype *p = SPH_TBP(simpar)+i;
			if(IsNowOneOfTwoSteps(IndT_NOWTSUBDIV(simpar),GetTsubPower(p),GetSphTsubPower(p),IndT_MAXTSUBPOWER(simpar)))
			{
				p->temp = GetT(simpar, p);
				p->mu   = getmu(simpar, p->temp,p->rho);
				p->nden = GAS_RHOS2RHOR(simpar)*(p->rho)*(1.0-GAS_YP(simpar)-(p->metallicity))/mH; 
				p->temp = GetTempGivenMu(simpar, p);
			}
		}
	}
	else {
#ifdef _OPENMP
#pragma omp parallel  private(i)
#endif
		for(i=0;i<SPH_NP(simpar);i++){
			treesphparticletype *p = SPH_TBP(simpar)+i;
			if(IsNowOneOfTwoSteps(IndT_NOWTSUBDIV(simpar),GetTsubPower(p),GetSphTsubPower(p),IndT_MAXTSUBPOWER(simpar)))
			{
				p->nden = GAS_RHOS2RHOR(simpar)*(p->rho)*(1.0-GAS_YP(simpar)-(p->metallicity))/mH; 
			}
		}
	}
}

static int readcalcevoltable=1;
static int past_stepcount=0;
//jhshin2

static int readcooltable=1;
void CoolingAndHeating(SimParameters *simpar){
	double H0incgs,Hincgs,zp1,dTcgs,gamma,gammam1,mgamma;
	double redshift=0;
	gamma = GAS_GAMMA(simpar);
	gammam1 = gamma-1;
	mgamma = -gamma;
	if(BGEXPAND(simpar) == 'Y'){
		H0incgs  = H0/Mpccgs*HUBBLE(simpar);
		zp1 = AMAX(simpar)/ANOW(simpar);
		Hincgs = H0incgs * sqrt(OMEP(simpar)*zp1*zp1*zp1+(1-OMEP(simpar)-OMEPLAM(simpar))*zp1*zp1+OMEPLAM(simpar));
		dTcgs = 1/Hincgs;
	}
	else {
		dTcgs = FFTIME(simpar);
	}

	if(MYID(simpar)==0){
		printf("Now entering into cooling and heating routine %g %g : %g \n",
				GAS_ENTROPYFACT(simpar),Hincgs,
				GAS_ENTROPYFACT(simpar)*(gammam1)/Hincgs);
	}
	if(readcooltable){
		void readcoolingdata_(int *);
		readcoolingdata_(&(MYID(simpar)));
		readcooltable = 0;
        if(MYID(simpar)==0) printf("Finish to read all cooling&heating data. \n");
	}
//jhshin1
	if(BGEXPAND(simpar) == 'Y') {
		redshift = AMAX(simpar)/ANOW(simpar)-1.0;
		if( readcalcevoltable==1 || STEPCOUNT(simpar)>past_stepcount ){
			double zp1 = 1+redshift;
		    void getcoolheatredshift_(double *);
            getcoolheatredshift_(&redshift);
			void calcevoldata_(double *, int *, int *, int *, float *);
			float a1,a2;
			a1 = ANOW(simpar);
			a2 = a1+ASTEP(simpar);
			double globaltimestep=GetTime(simpar, a1,a2)/(365.L*24.L*3600.L); //yr
			int include_shielding;					//1:yes, 0:no
            {   
				double redshift = AMAX(simpar)/ANOW(simpar)-1;
               	void InitializeMu(SimParameters *, double); InitializeMu(simpar, redshift);
				//read Mu data for fortran
				void readmu_(int *, double *);
                readmu_(&MYID(simpar),aMU);
            }

		 	include_shielding=0;
		 	calcevoldata_(&globaltimestep,&include_shielding,&MYID(simpar),&NID(simpar),&GAS_YP(simpar));
			if(redshift <= 8.9){
		  		include_shielding=1;
		  		calcevoldata_(&globaltimestep,&include_shielding,&MYID(simpar),&NID(simpar),&GAS_YP(simpar));
			}
			if(MYID(simpar)==0){ 
				printf("Reading enetropy evolution data. %g %g\n", redshift, globaltimestep);
			}
			past_stepcount=STEPCOUNT(simpar);
			readcalcevoltable=0;
		}
	}
	else{
		 if(readcalcevoltable){
			redshift = 0; 
			void getcoolheatredshift_(double*);
			void calcevoldata_(double *, int *, int *, int *,float *); 
			getcoolheatredshift_(&redshift);
            if(MYID(simpar)==0) printf("Finish to read c&h(redshift). \n");
			double globaltimestep=ASTEP(simpar)*FFTIME(simpar)/(365.0L*24.0L*3600.0L);    //yr
            int include_shielding;                                                 // 1: yes, 0:no
			{
                //read Mu data for fortran
				void readmu_(int *, double *);
                readmu_(&MYID(simpar),aMU);
            }
            include_shielding=0;
			calcevoldata_(&globaltimestep,&include_shielding,&MYID(simpar),&NID(simpar),&GAS_YP(simpar));
            if(MYID(simpar)==0) printf("Finish to calc evoldata. \n");
			readcalcevoltable=0;
			if(MYID(simpar)==0){ 
				printf("Reading enetropy evolution data. \n");
			}
		}
	}
	MPI_Barrier(MPI_COMM(simpar));
	
	{	
		void getdeltae_(double*,double*,double*,int*,int*,double*);
		void getcoolheatfinal_(double*,double*,double*,int*,double*,double*);
		double coolrate,heatrate;
		float UVShieldDen = GAS_UVSHIELDDEN(simpar) * mH /(1-GAS_YP(simpar));
		long long i;
#ifdef _OPENMP
#pragma omp parallel 
#endif
		for(i=0;i<SPH_NP(simpar);i++){
			int Uvshield=0;
			int dUvshield=0;
			treesphparticletype *p = SPH_TBP(simpar)+i;
			double gastemp,gasnden,gasmetal,gasener;
			int sub_nstep;
       		int ii = GetTsubPower(p);
       		int jj = GetSphTsubPower(p);
       		int kk = GetSphNextTsubPower(p);
			if(ii>jj) jj=ii;
			/*
			if(isnowdmstep(IndT_NOWTSUBDIV(simpar),jj,maxTsubpower)){
			*/
			if(IsNowOneOfTwoSteps(IndT_NOWTSUBDIV(simpar),ii,jj,IndT_MAXTSUBPOWER(simpar))){
				gastemp  = p->temp;
				gasnden  = p->nden;
				gasmetal = MAX(1.E-4,p->metallicity)/0.02;
				double delta_e_real,delta_e_sim,delta_e_sim_1,delta_e_sim_2,mid_entropy,mid_temp;
				if(BGEXPAND(simpar)=='Y'){
					if((redshift<8.9&&gasnden>GAS_UVSHIELDDEN(simpar))) Uvshield=1;
					if((redshift<8.9&&gasnden<(GAS_UVSHIELDDEN(simpar)+GAS_DUVSHIELDDEN(simpar)))) dUvshield=1;
				}
				if(CONT_HALFSTEP(simpar)==HALFPULL){
	           		sub_nstep= jj+1;
					if(dUvshield==1){
						int iflag=0,jflag=1;
						double delta_e_reali,delta_e_realj;
						float fraci,fracj;
						fraci = (gasnden-GAS_UVSHIELDDEN(simpar))/GAS_DUVSHIELDDEN(simpar);
						fracj = 1-fraci;
						getdeltae_(&gastemp,&gasnden,&gasmetal,&sub_nstep,&iflag,&delta_e_reali);
        	        	getdeltae_(&gastemp,&gasnden,&gasmetal,&sub_nstep,&jflag,&delta_e_realj);
						delta_e_real = fraci*delta_e_reali + fracj*delta_e_realj;
					}
					else {
        	        	getdeltae_(&gastemp,&gasnden,&gasmetal,&sub_nstep,&Uvshield,&delta_e_real);
					}
                    delta_e_sim=GAS_ENTROPYFACT(simpar)*gammam1*delta_e_real*pow(p->rho,mgamma);
				}else if(CONT_HALFSTEP(simpar)==HALFPUSH){ 
					if(kk>jj) jj=kk; 
					sub_nstep= jj+1; 
					if(dUvshield==1){
						int iflag=0,jflag=1;
						double delta_e_reali,delta_e_realj;
						float fraci,fracj;
						fraci = (gasnden-GAS_UVSHIELDDEN(simpar))/GAS_DUVSHIELDDEN(simpar);
						fracj = 1-fraci;
						getdeltae_(&gastemp,&gasnden,&gasmetal,&sub_nstep,&iflag,&delta_e_reali);
        	        	getdeltae_(&gastemp,&gasnden,&gasmetal,&sub_nstep,&jflag,&delta_e_realj);
						delta_e_real = fraci*delta_e_reali + fracj*delta_e_realj;
					}
					else { 
						getdeltae_(&gastemp,&gasnden,&gasmetal,&sub_nstep,&Uvshield,&delta_e_real); 
					}
					delta_e_sim=GAS_ENTROPYFACT(simpar)*gammam1*delta_e_real*pow(p->rho,mgamma); 
				}else{ 
					sub_nstep= jj+1; 
					if(dUvshield==1){
						int iflag=0,jflag=1;
						double delta_e_reali,delta_e_realj;
						float fraci,fracj;
						fraci = (gasnden-GAS_UVSHIELDDEN(simpar))/GAS_DUVSHIELDDEN(simpar);
						fracj = 1-fraci;
						getdeltae_(&gastemp,&gasnden,&gasmetal,&sub_nstep,&iflag,&delta_e_reali);
        	        	getdeltae_(&gastemp,&gasnden,&gasmetal,&sub_nstep,&jflag,&delta_e_realj);
						delta_e_real = fraci*delta_e_reali + fracj*delta_e_realj;
					}
					else { 
						getdeltae_(&gastemp,&gasnden,&gasmetal,&sub_nstep,&Uvshield,&delta_e_real); 
					}

					delta_e_sim_1=GAS_ENTROPYFACT(simpar)*gammam1*delta_e_real*pow(p->rho,mgamma);
                    gasener = gammam1*kBovermH/(p->mu)*(p->rho)*(p->temp);
                    mid_temp = gastemp*(1.0+delta_e_real/gasener);
					if(kk>=jj) {
						sub_nstep = kk+1;
					}
					else if(IsNowSphStep(IndT_NOWTSUBDIV(simpar),jj-1,IndT_MAXTSUBPOWER(simpar))){
						sub_nstep = jj-1+1;
					}
					else {
						sub_nstep = jj+1;
					}
					if(dUvshield==1){
						int iflag=0,jflag=1;
						double delta_e_reali,delta_e_realj;
						float fraci,fracj;
						fraci = (gasnden-GAS_UVSHIELDDEN(simpar))/GAS_DUVSHIELDDEN(simpar);
						fracj = 1-fraci;
						getdeltae_(&mid_temp,&gasnden,&gasmetal,&sub_nstep,&iflag,&delta_e_reali);
        	        	getdeltae_(&mid_temp,&gasnden,&gasmetal,&sub_nstep,&jflag,&delta_e_realj);
						delta_e_real = fraci*delta_e_reali + fracj*delta_e_realj;
					}
					else { 
						getdeltae_(&mid_temp,&gasnden,&gasmetal,&sub_nstep,&Uvshield,&delta_e_real); 
					}

					delta_e_sim_2=GAS_ENTROPYFACT(simpar)*gammam1*delta_e_real*pow(p->rho,mgamma);
					delta_e_sim=delta_e_sim_1+delta_e_sim_2;
				}
				p->delta_e=delta_e_sim;
				if(isnan(delta_e_sim)){
					printf("P%d index=%ld: %g %g %g %g \n",MYID(simpar),i,gastemp,gasnden,gasmetal,delta_e_real); 
				} 
			}
		}
	}
//jhshin2
}

void DeleteHydroPadding(SimParameters *simpar){
	free(SPH_TBPP(simpar)); SPH_NPAD(simpar) = 0;
	free(STAR_TBPP(simpar)); STAR_NPAD(simpar) = 0;
	free(AGN_TBPP(simpar)); AGN_NPAD(simpar) = 0;
}

void MakeStar(SimParameters *simpar);

int nsubsphpow,nsubfixedpow,nsubnbodypow;


void sphmain(SimParameters *simpar){
	int Numnear = SPH_NUMNEAR(simpar);
	float astep = ASTEP(simpar)/pow(2.,IndT_MAXTSUBPOWER(simpar));


	if(BGEXPAND(simpar)== 'Y') {
		DetermineSphFactor(simpar);
		GetEvolFactorArr(simpar, ANOW(simpar),astep); 
	}
	else {
		StaticDetermineSphFactor(simpar);
        StaticGetEvolFactorArr(simpar,ANOW(simpar),astep); 
	}

	Init_dAs(simpar);

	float width = 4;

	ppadding(SPH_TBP(simpar), SPH_NP(simpar), (void**)(&SPH_TBPP(simpar)), &SPH_NPAD(simpar), TSPH_DDINFO(simpar), 
			NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar), 3);
	ppadding(STAR_TBP(simpar), STAR_NP(simpar), (void**)(&STAR_TBPP(simpar)), &STAR_NPAD(simpar), TSTAR_DDINFO(simpar), 
			NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar), 3);
	ppadding(AGN_TBP(simpar), AGN_NP(simpar), (void**)(&AGN_TBPP(simpar)), &AGN_NPAD(simpar), TAGN_DDINFO(simpar), 
			NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar), 3);


	IndT_IFLAGFIXEDLIST(simpar) = 0;
	if(0){
		/* To get more information on whether starting the fixedneighbor time step or not */
		void GetFixedSubStepInfo(SimParameters *);
		GetFixedSubStepInfo(simpar);
		IndT_NSUBFIXED(simpar) = MAX(IndT_NSUBFIXED(simpar),IndT_NSUBNBODY(simpar));
		if(IndT_NSUBFIXED(simpar) > IndT_NSUBGAS(simpar)) IndT_IFLAGFIXEDLIST(simpar) = 0;
		else if (isnowstep(IndT_NOWTSUBDIV(simpar), IndT_NSUBFIXED(simpar),IndT_MAXTSUBPOWER(simpar))) IndT_IFLAGFIXEDLIST(simpar) = 0;
		else IndT_IFLAGFIXEDLIST(simpar) = 1;
		if(MYID(simpar)==0) printf("Now SubStep %d : %d %d %d\n",IndT_IFLAGFIXEDLIST(simpar),
				IndT_NSUBNBODY(simpar),IndT_NSUBGAS(simpar),IndT_NSUBFIXED(simpar));
	}

	if(IndT_IFLAGFIXEDLIST(simpar) == 1){
#ifdef NOTYET
		void GetSphParticleAddress(SimParameters *); GetSphParticleAddress(simpar);
		void GetInfoOfNeighborWorkList(SimParameters *); GetInfoOfNeighborWorkList(simpar);
		MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) printf("Now herereereer \n");
		void SendNeighborList(SimParameters *); SendNeighborList(simpar);
		MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) printf("Now herereereer \n");
		void DoSphDenFixedNeighbor(SimParameters *); DoSphDenFixedNeighbor(simpar);
		MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) printf("Now herereereer \n");
		MPI_Finalize();exit(0);
#endif
	}
	else {
		void FreeFixedNeighborList(SimParameters *); FreeFixedNeighborList(simpar);
		void InitFixedNeighborList(SimParameters *); InitFixedNeighborList(simpar);
		size_t SPHBuildLinkedList(SimParameters *);    SPHBuildLinkedList(simpar); /*Making the linked list for the SPH calculation */
		FindSphDensity(simpar);
	}



//jhshin1 
	if(SPH_INITFLAG(simpar) == 1){ /* Initializing Entropy as T=1e4K */ 
		double redshift = AMAX(simpar)/ANOW(simpar)-1; 
		void InitializeMu(SimParameters *, double); InitializeMu(simpar, redshift);			
		if(BGEXPAND(simpar)=='Y'){ 
			SPH_INITMASS(simpar) = SPH_MASSINSOLARMASS(simpar);
			treesphparticletype *sp = SPH_TBP(simpar);
			float temp;
			long ii;
			for(ii=0;ii<SPH_NP(simpar);ii++){
				temp       = GetT(simpar, sp);
				sp->mu     = getmu(simpar, temp,sp->rho);
				sp++;
			}
		}else if (SIMMODEL(simpar)==Static){
			if(STEPCOUNT(simpar)==1){
				long long ii;
				if(MYID(simpar)==0) printf("jhshin's Entropy initializing. \n");
				treesphparticletype *sp = SPH_TBP(simpar);
				for(ii=0;ii<SPH_NP(simpar);ii++){
					sp->mu = 1;
					sp->Entropy= kBovermH/(sp->mu)*1e4*pow(sp->rho,1-GAS_GAMMA(simpar));
					sp->temp = GetT(simpar, sp); 
					sp++;
				}
			}
			else {
				long long ii;
				treesphparticletype *sp = SPH_TBP(simpar);
				for(ii=0;ii<SPH_NP(simpar);ii++){
					sp->mu = 1;
					sp++;
				}
			}
		}
		SPH_INITFLAG(simpar) = 0;
	}
	Calc_mutempnden(simpar);
//jhshin2
	if(GAS_BGHEATFLAG(simpar) == 'Y' || GAS_COOLFLAG(simpar)=='Y') CoolingAndHeating(simpar);




	if(IndT_IFLAGFIXEDLIST(simpar) == 1){
		if(MYID(simpar)==0) printf("Before 127 \n");
		void DoSphComFixedNeighbor(SimParameters *); DoSphComFixedNeighbor(simpar);
		if(MYID(simpar)==0) printf("Before 127 \n");
		void UpdateTimeLimiter4Ghosts(SimParameters *); UpdateTimeLimiter4Ghosts(simpar);
		if(MYID(simpar)==0) printf("Before 127 \n");
	}
	else {
		CompleteSph(simpar);
		DestroyHydroLinkedCell(simpar); /*Destroying the linked list for the SPH calculation */
	}

	if(GAS_SFFLAG(simpar) == 'Y') MakeStar(simpar);




	DeleteHydroPadding(simpar);




	/*
	int determine_fixed_timestep(SimParameters *);
	int ii = determine_fixed_timestep(simpar);
	*/


}
void CHECK1(SimParameters *simpar){
	int ix,iy,iz;
	long i;
	float xmin,ymin,zmin,xmax,ymax,zmax;
	long invCellWidth = BASICCELL_INVCELLWIDTH(simpar);
	xmin = ymin = zmin = 1.e20;
	xmax = ymax = zmax =-1.e20;
	for(i=0;i<SPH_NP(simpar);i++){
		ix = ((SPH_TBP(simpar)+i)->x - SIM_LXMIN(simpar,sph))*invCellWidth;
		iy = ((SPH_TBP(simpar)+i)->y - SIM_LYMIN(simpar,sph))*invCellWidth;
		iz = ((SPH_TBP(simpar)+i)->z - SIM_LZMIN(simpar,sph))*invCellWidth;
			xmin = MIN(xmin,(SPH_TBP(simpar)+i)->x);
			ymin = MIN(ymin,(SPH_TBP(simpar)+i)->y);
			zmin = MIN(zmin,(SPH_TBP(simpar)+i)->z);
			xmax = MAX(xmax,(SPH_TBP(simpar)+i)->x);
			ymax = MAX(ymax,(SPH_TBP(simpar)+i)->y);
			zmax = MAX(zmax,(SPH_TBP(simpar)+i)->z);
	}
	printf("P%d has %g %g %g : %g %g %g    ::: %ld\n",MYID(simpar),xmin,ymin,zmin,xmax,ymax,zmax,SPH_NP(simpar));
	MPI_Finalize();exit(99);
}
void SetTemp2Sph(SimParameters *simpar){
	long long i;
	float red,Temp,OmegaB,OmegaC,OmegaL,HO,Tnow,Yp;
	float ai = 1.;
	float meanrhoB;
	red = AMAX(simpar)/ai - 1;
	OmegaB = OMEPB(simpar);
	OmegaC = OMEP(simpar);
	OmegaL = OMEPLAM(simpar);
	HO = 100.*HUBBLE(simpar);
	Tnow = GAS_TRAD0(simpar);
	Yp = GAS_YP(simpar);
	myrecfast_(&red,&Temp,&OmegaB,&OmegaC,&OmegaL,&HO,&Tnow,&Yp);
	if(MYID(simpar)==0) {
		printf("######################################\n");
		printf("From RecFast: The mean gas temperature is %g at z= %g\n",Temp,red);
		printf("######################################\n");
	}
	meanrhoB = OMEPB(simpar)/OMEP(simpar);
#ifdef _OPENMP
#pragma omp parallel 
#endif
	for(i=0;i<SPH_NP(simpar);i++){
		float dTemp = (GAS_GAMMA(simpar)-1)*Temp*((SPH_TBP(simpar)+i)->rho/meanrhoB -1 );
		float TEMP = Temp+dTemp;
		float MU = getmu(simpar, TEMP,(SPH_TBP(simpar)+i)->rho);
		(SPH_TBP(simpar)+i)->Entropy = kBovermH * TEMP/MU * pow((SPH_TBP(simpar)+i)->rho,1-GAS_GAMMA(simpar));
	}
	if(MYID(simpar)==0){
		printf("################ Setting Initial Condition for the SPH particles. ##########\n\n");
		printf("P%d has Entropy_s=%g for first particle\n\n\n",
				MYID(simpar),SPH_TBP(simpar)->Entropy);
	}
}
void FindInitialGasTemp(SimParameters *simpar){
	if(MYID(simpar)==0) TIMER_START(0);
	if(BGEXPAND(simpar) == 'Y') DetermineSphFactor(simpar);
	HydroBuildLinkedList(simpar); /*Making the linked list for the SPH calculation */
	FindSphDensity(simpar);
	DestroyHydroLinkedCell(simpar); /*Destroying the linked list for the SPH calculation */
	SetTemp2Sph(simpar);
	float Temp = GetT(simpar, SPH_TBP(simpar));
	float Tmu = SPH_TBP(simpar)->Entropy *pow(SPH_TBP(simpar)->rho,GAS_GAMMA(simpar)-1)*mHoverkB;
	if(MYID(simpar)==0) {
		printf("P%d has sph den/temp/hsml/Tmu %g %g %g %g ith numnear= %d\n",MYID(simpar), 
			SPH_TBP(simpar)->rho,Temp,
			SPH_TBP(simpar)->hsml,Tmu,SPH_NUMNEAR(simpar));
		TIMER_STOP(0);
		printf("Time needed for the Initial Gas Temperature estimation %g second\n\n",ELAPSED_TIME(0));
	}
}
