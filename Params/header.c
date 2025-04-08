#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "params.h"
//#include "mpirms.h"

void Mpi_Basic_Set(SimParameters *, MPI_Comm );


void CheckConsistencyParam(SimParameters *simpar){
	/*
	if(SIMMODEL(simpar) != Cosmos){
		printf("Error in setting the total mass in the non-expanding simulation box\n");
		printf("You don't need to set the value in the cosmological simulation \n");
		printf("But you have to set the value in the non-expanding simulation \n");
		printf("It should be positive one. The current value is %g\n",TOTMASS(simpar));
		exit(999);
	}
	*/
}



void write_head(FILE *wp, SimParameters *wsimpar){
	int isize,ncnt = 0;
	FILE_HEADER(fprintf,wp, ,wsimpar);
}
void DetermineFreeFallTime(SimParameters *simpar){
		double psize = (COSBOXSIZE(simpar)/NX(simpar)) * Mpccgs;
		double pmass = TOTMASS(simpar)*onesolarmass/ pow(NX(simpar), 3.L);
		double meanrhor = pmass/psize/psize/psize;
		FFTIME(simpar) = 1./sqrt(P_G*meanrhor);
		FFTELLB(simpar) = (FFTIME(simpar)*FFTIME(simpar))/(psize*psize);
}
void read_slab_head(FILE *fp, SimParameters *simpar){
	SimParameters rsimpar;
	int isize,ncnt=0;
	char line[MAX_LINE_LENGTH];

	char defaultcosmos[190];

	sprintf(defaultcosmos,"WMAP5");


	AGN_NP(simpar) = SPH_NP(simpar) = VORO_NP(simpar) = STAR_NP(simpar) = DM_NP(simpar) = NPSUM(simpar) = 0;
	NSPACE(simpar) = 1;

	SIMMODEL(simpar) = Cosmos;

	/* Loop until it reaches P_Closing */
	while(strcmp(P_Closing,fgets(line, MAX_LINE_LENGTH,fp)) != 0){
		int ncnt=0;
		if(line[0] != '#'){
			if(strstr(line,"define") != NULL) FILE_SLAB_HEADER(sscanf,line,&,simpar);
			if(strstr(line,"=") != NULL && strstr(line,"define")!=NULL && ncnt ==0) {
				fprintf(stderr,"Warning: the following parameter line is unknown:\n %s\n\n\n",line);
			}
		}
	}
	NXNY(simpar) = NX(simpar)*NY(simpar);
	/* Now change the stepnum to this step */
	STEPNUM(simpar) = STEPCOUNT(simpar);
	ZINIT(simpar) = AMAX(simpar) + 1.;
	if(SIMMODEL(simpar) != Cosmos) DetermineFreeFallTime(simpar);
	CheckConsistencyParam(simpar);
}

void read_head(FILE *fp, SimParameters *simpar){
	SimParameters rsimpar;
	int isize,ncnt=0;
	char line[MAX_LINE_LENGTH];

	char defaultcosmos[190];

	sprintf(defaultcosmos,"WMAP5");
	/*
	mk_default_param(simpar,defaultcosmos);
	*/


	AGN_NP(simpar) = SPH_NP(simpar) = VORO_NP(simpar) = STAR_NP(simpar) = DM_NP(simpar) = NPSUM(simpar) = 0;

	SIMMODEL(simpar) = Cosmos;

	/* Loop until it reaches P_Closing */
	while(strcmp(P_Closing,fgets(line, MAX_LINE_LENGTH,fp)) != 0){
		int ncnt=0;
		/* check whether it is a comment */
		if(line[0] != '#'){
			if(strstr(line,"define") != NULL) FILE_HEADER(sscanf,line,&,simpar);
			if(strstr(line,"=") != NULL && strstr(line,"define")!=NULL && ncnt ==0) {
				fprintf(stderr,"Warning: the following parameter line is unknown:\n %s\n\n\n",line);
			}
		}
	}
	/* Now change the stepnum to this step */
	STEPNUM(simpar) = STEPCOUNT(simpar);
	ZINIT(simpar) = AMAX(simpar) + 1.;
	if(SIMMODEL(simpar) != Cosmos) DetermineFreeFallTime(simpar);
	CheckConsistencyParam(simpar);
}

void write_default_sim_parameter_file(FILE *wp, SimParameters *wsimpar){
	int isize,ncnt = 0;
	DEFAULT_PARAMS(fprintf,wp,, wsimpar);
}

void write_sim_parameter_file(FILE *wp, SimParameters *wsimpar){
	int isize,ncnt = 0;
	PARAMS(fprintf,wp,,wsimpar);
}
void MkDefaultCosParam(SimParameters *defsim, char *cosmology){
	double zi;
	
	PORDER(defsim) = 2;
	/* Simulation Parameters for WMAP 5-year */
	if(strcmp(cosmology,"WMAP3")==0){
		OMEP(defsim) = 0.238;
		OMEPLAM(defsim) = 0.762;
		OMEPB(defsim) = 0.042;
		HUBBLE(defsim) = 0.732;
		NPOW(defsim) = 0.958;
		BIAS8(defsim) = 1.314;
	}
	else {
		OMEP(defsim) = 0.26;
		OMEPLAM(defsim) = 0.74;
		OMEPB(defsim) = 0.044;
		HUBBLE(defsim) = 0.72;
		NPOW(defsim) = 0.96;
		BIAS8(defsim) = 1.26;
	}
	sprintf(POWFILENAME(defsim),"camb.z=47.dat");
	WLAM0(defsim) = -1;
	WLAM1(defsim) = 0;
	FNL(defsim) = 0.;
	GNL(defsim) = 0.;
	COSBOXSIZE(defsim) = 1024.;
	AMAX(defsim) = 48.;
	ANOW(defsim) = 1.;
	NX(defsim) = NY(defsim) = NZ(defsim) = 1024;
	NXNY(defsim) = NX(defsim) * NY(defsim);
	ZINIT(defsim) = AMAX(defsim)/ANOW(defsim)-1;
	THETA(defsim) = 0.3;
	RSPHERE(defsim) = 4;
	ASTEP(defsim) = 0.025;
	NSTEP(defsim) = 1881;
	STEPCOUNT(defsim) = 1.;
	FLAG_IC(defsim) = 0;
	sprintf(GRAFIC_DIRECTORY(defsim)," ");
	ISEED(defsim) = -56;
	POWREADFLAG(defsim) = 1;
	NDDINFO(defsim) = 0;
#ifdef XYZDBL
	XYZSHIFT(defsim) = 1;
#else
	XYZSHIFT(defsim) = 0;
#endif
	sprintf(RV_FILE(defsim),"INITIAL");
	sprintf(RV_FILEPREFIX(defsim),"INITIAL");
	zi = ZINIT(defsim);
	OMEI(defsim) = OMEP(defsim) * pow(1+zi,3)/(OMEP(defsim)*pow(1+zi,3) +
			OMEPLAM(defsim)+(1-OMEP(defsim)-OMEPLAM(defsim))*pow(1+zi,2));
	sprintf(ANIM_VIEWFILE(defsim),"\0");
	GAS_VISCOSITY(defsim) = 0;
}






void ReadSimulationParameters(FILE *simfile, int *icont, SimParameters *simpar){
	char filename[100];
	FILE *rvfile;
	int abortflag;
	if(MYID(simpar)==0){
		MPI_Comm com = MPI_COMM(simpar);
		int myid, nid;
		myid = MYID(simpar);
		nid = NID(simpar);

		read_head(simfile, simpar);

		if(SIMMODEL(simpar)==Cosmos){
			printf("*** INPUT COSMOLOGICAL PARAMETER LIST ***\n\n");
			printf("Box size = %g Mpc/h, Hubble parameter =  %g\n",COSBOXSIZE(simpar),HUBBLE(simpar));
			printf("Power index = %g Omega_m0= %g Omega_b0= %g Omega_L0= %g Bias_8Mpc= %g \n",
					NPOW(simpar), OMEP(simpar), OMEPB(simpar), OMEPLAM(simpar), BIAS8(simpar));
			printf("\n **** Time Parameters **** \n");
			printf("Initial Redshift = %g Delta_a =  %g Current a = %g\n\n",
					AMAX(simpar)-1,ASTEP(simpar), ANOW(simpar));
			printf("RV_FILE= %s\n",RV_FILE(simpar));
		}
		else if(SIMMODEL(simpar)== Static){
			printf("*** INPUT STATIC PARAMETER LIST ***\n\n");
			printf("Box size = %g Mpc/h, Total mass in Msun/h =  %g\n",STATBOXSIZE(simpar),TOTMASS(simpar));
			printf("Free Fall Time = %g \n", FFTIME(simpar)/365.E6/24./3600.);
			printf("Final Time = %g Delta_T =  %g Current Time = %g\n\n",
					AMAX(simpar)-1,ASTEP(simpar), ANOW(simpar)-1);
		}
		else if(SIMMODEL(simpar)== KH){
			printf("Kelvin Helmholtz Instability Test is set\n");
			printf("box size nx/ny/nz= %ld %ld %ld\n", NX(simpar), NY(simpar), NZ(simpar));
			KH_XMIN(simpar) = KH_YMIN(simpar) = 0;
		}
		else if(SIMMODEL(simpar) == RT){
			printf("Rayleigh Tayloer Instability Test is set\n");
			printf("box size nx/ny/nz= %ld %ld %ld\n", NX(simpar), NY(simpar), NZ(simpar));
			RT_XMIN(simpar) = RT_YMIN(simpar) = 0;
		}
		else if(SIMMODEL(simpar) == Kepler){
			printf("Kepler L-mom. Test is set\n");
			printf("box size nx/ny/nz= %ld %ld %ld\n", NX(simpar), NY(simpar), NZ(simpar));
			KP_XMIN(simpar) = KP_YMIN(simpar) = 0;
		}
		else if(SIMMODEL(simpar) == MkGlass2D){
			printf("Making 2D Glass IC is set\n");
			printf("box size nx/ny/nz= %ld %ld %ld\n", NX(simpar), NY(simpar), NZ(simpar));
			GL2D_XMIN(simpar) = GL2D_YMIN(simpar) = 0;
		}
		MPI_COMM(simpar) = com;
	}
	MPI_Bcast(simpar, sizeof(SimParameters), MPI_BYTE, 0, COM(simpar));

	Mpi_Basic_Set(simpar, MPI_COMM(simpar));

	abortflag = 0;
	if(MYID(simpar) ==0){
		if(strcmp(RV_FILE(simpar),"INITIAL")==0) {
			*icont = 0;
		}
		else {
			sprintf(filename,"%s.%.5d",RV_FILE(simpar), 0);
			rvfile = fopen(filename,"r");
			if(rvfile == NULL){
				fprintf(stderr,"Can't open initial rvfile %s -exiting \n",RV_FILE(simpar));
				fprintf(stderr,"Can't open initial rvfile %s -exiting \n",filename);
				abortflag = 1;
			}
			else {
				fprintf(stderr,"Successfully opened rvfile %s -starting simulation\n",RV_FILE(simpar));
				fclose(rvfile);
				*icont = 1;
			}
		}
	}
	MPI_Bcast(&abortflag, 1, MPI_INT, 0, MPI_COMM(simpar));
	if(abortflag){
		fprintf(stderr,"Aborting can't find initial rvfile\n");
		MPI_Finalize();
		exit(0);
	}
	MPI_Bcast(icont, 1, MPI_INT, 0, MPI_COMM(simpar));
	if(*icont ==1) ANOW(simpar) = 1.;
	NXNY(simpar) = NX(simpar)*NY(simpar);
}



SimParameters readcheck_head(FILE *fp){
	SimParameters rsimpar, *simpar;
	simpar = &rsimpar;
	int isize,ncnt=0;
	char line[MAX_LINE_LENGTH];

	char defaultcosmos[190];

	sprintf(defaultcosmos,"WMAP5");

	mk_default_param(simpar,defaultcosmos);
	SPH_NP(simpar) = VORO_NP(simpar) = STAR_NP(simpar) = DM_NP(simpar) = AGN_NP(simpar) = NPSUM(simpar) = 0;
	SIMMODEL(simpar) = Cosmos;
	/* Loop until it reaches P_Closing */
	while(strcmp(P_Closing,fgets(line, MAX_LINE_LENGTH,fp)) != 0){
		int ncnt=0;
		/* check whether it is a comment */
		if(line[0] != '#'){
			if(strstr(line,"define") != NULL) FILE_HEADER(sscanf,line,&,simpar);
			if(strstr(line,"=") != NULL && strstr(line,"define")!=NULL && ncnt ==0) {
				fprintf(stderr,"Warning: the following parameter line is unknown:\n %s\n\n\n",line);
			}
			else {
				printf(line);
			}
		}
	}
	/* Now change the stepnum to this step */
	STEPNUM(simpar) = STEPCOUNT(simpar);
	ZINIT(simpar) = AMAX(simpar) + 1;
	if(SIMMODEL(simpar) != Cosmos) DetermineFreeFallTime(simpar);
	CheckConsistencyParam(simpar);
	return rsimpar;
}

SimParameters  read_sim_parameter_file(SimParameters presimpar,FILE *fp){
	SimParameters rsimpar, *simpar;
	int isize,ncnt=0;
	char line[MAX_LINE_LENGTH];
	rsimpar = presimpar;
	simpar = &rsimpar;
	while(strcmp(P_Closing,fgets(line, MAX_LINE_LENGTH,fp)) != 0){ /* Loop until it reaches P_Closing */
		int ncnt=0;
		if(line[0] != '#'){ /* check whether it is a comment */
			if(strstr(line,"define") != NULL) FILE_HEADER(sscanf,line,&,simpar);
			if(strstr(line,"=") != NULL && strstr(line,"define")!=NULL && ncnt ==0) {
				DEBUGPRINT("Warning: the following parameter line is unknown:\n %s\n\n\n",line);
			}
		}
	}
	STEPNUM(simpar) = STEPCOUNT(simpar);
	if(SIMMODEL(simpar) != Cosmos) DetermineFreeFallTime(simpar);
	CheckConsistencyParam(simpar);
	return rsimpar;
}
void mk_default_kh_param(SimParameters *defsim){
		BGEXPAND(defsim) = 'N';
		GAS_TYPE(defsim) = 'V';
		GAS_VISCOSITY(defsim) = 0;
		SIMMODEL(defsim) = KH;
		GAS_CONSTMU(defsim) = 'Y';
		HUBBLE(defsim) = 1;
		STATBOXSIZE(defsim) = NX(defsim);
		TOTMASS(defsim) = 2.7755E11*0.26*pow(STATBOXSIZE(defsim),3.L);

		AMAX(defsim) = 94;
		ANOW(defsim) = 1.;

		ASTEP(defsim) = 0.005;
		NSTEP(defsim) = 94000;
		NX(defsim) = NY(defsim) = 128;
		NZ(defsim) = 1;
		NXNY(defsim) = NX(defsim) * NY(defsim);


		KH_Rho1(defsim) = 1;
		KH_Rho2(defsim) = 2;
		KH_Vel1(defsim) = 0.5;
		KH_Vel2(defsim) = -0.5;
		KH_Deltay(defsim) = 0.025;
		KH_Pressure(defsim) = 5./2.;
		KH_Vperturb(defsim) = 0.01;
		KH_XMIN(defsim) = 0;
		KH_XMAX(defsim) = 1;
		KH_YMIN(defsim) = 0;
		KH_YMAX(defsim) = 1;
		KH_GridSize(defsim) = (KH_XMAX(defsim)-KH_XMIN(defsim))/NX(defsim)*4;
		KH_OA(defsim) = 0.33333333;

		GRAVITY(defsim) = 'N';
		GAS_SFFLAG(defsim) = 'N';
		GAS_COOLFLAG(defsim) = 'N';
		GAS_SNFBFLAG(defsim) = 'N';
		GAS_BGHEATFLAG(defsim) = 'N';
	}
void mk_default_rt_param(SimParameters *defsim){
		BGEXPAND(defsim) = 'N';
		GAS_TYPE(defsim) = 'V';
		GAS_VISCOSITY(defsim) = 0;
		SIMMODEL(defsim) = RT;
		TOTMASS(defsim) = 2.7755E11*0.26*pow(STATBOXSIZE(defsim),3.L);
		RT_XMAX(defsim) = 0.5;
		RT_YMAX(defsim) = 1;
		RT_XMIN(defsim) = 0;
		RT_YMIN(defsim) = 0;
		NX(defsim) = 128;
		NY(defsim) = 256;
		NZ(defsim) = 1;
		RT_GridSize(defsim) = (RT_XMAX(defsim)-RT_XMIN(defsim))/NX(defsim)*4;
		DetermineFreeFallTime(defsim);
		RT_DEN1(defsim) = 1;
		RT_DEN2(defsim) =  2;
		RT_ACC(defsim) = -0.5;
		GRAVITY(defsim) = 'N';
		RT_Deltay(defsim) = 0.025;
		RT_Vperturb(defsim) = 0.025;
		RT_Phalf(defsim) = 10./7.;
		RT_OA(defsim) = 0.;
}
void mk_default_kepler_param(SimParameters *defsim){
	BGEXPAND(defsim) = 'N'; 
	GAS_TYPE(defsim) = 'V'; 
	GAS_VISCOSITY(defsim) = 0;
	SIMMODEL(defsim) = Kepler; 
	TOTMASS(defsim) = 2.7755E11*0.26*pow(STATBOXSIZE(defsim),3.L); 
	KP_XMAX(defsim) = 6.;
	KP_YMAX(defsim) = 6.; 
	KP_XMIN(defsim) = 0; 
	KP_YMIN(defsim) = 0; 
	NX(defsim) = 256; 
	NY(defsim) = 256; 
	NZ(defsim) = 1; 
	KP_GridSize(defsim) = (KP_XMAX(defsim)-KP_XMIN(defsim))/NX(defsim)*4; 
	DetermineFreeFallTime(defsim); 
	GRAVITY(defsim) = 'N'; 
	KP_OA(defsim) = 0.;
	KP_EPS(defsim) = 0.1;

	GAS_AlphaVis(defsim) = .2;
	GAS_BetaVis(defsim) = 0.1;
	GAS_ETAVIS(defsim) = 1;
	GAS_EPSVIS(defsim) = 0.01;
	GAS_VISCOSITY(defsim) = 0;
}
void mk_default_2dglass_param(SimParameters *defsim){
    BGEXPAND(defsim) = 'N';
    GAS_TYPE(defsim) = 'V';
    GAS_VISCOSITY(defsim) = 0;
    SIMMODEL(defsim) = MkGlass2D;
    TOTMASS(defsim) = 2.7755E11*0.26*pow(STATBOXSIZE(defsim),3.L);
    GL2D_XMAX(defsim) = 1.;
    GL2D_YMAX(defsim) = 1.;
    GL2D_XMIN(defsim) = 0;
    GL2D_YMIN(defsim) = 0;
    NX(defsim) = 256;
    NY(defsim) = 256;
    NZ(defsim) = 1;
    GL2D_GridSize(defsim) = (GL2D_XMAX(defsim)-GL2D_XMIN(defsim))/NX(defsim)*4;
    DetermineFreeFallTime(defsim);
    GRAVITY(defsim) = 'N';
    GL2D_OA(defsim) = 0.;
    GL2D_EPS(defsim) = 0.1;
    GL2D_Kappa(defsim) = 0.55;

    GAS_AlphaVis(defsim) = .2;
    GAS_BetaVis(defsim) = 0.1;
    GAS_ETAVIS(defsim) = 1;
    GAS_EPSVIS(defsim) = 0.01;
    GAS_VISCOSITY(defsim) = 0;
}


void mk_default_bowshock_param(SimParameters *defsim){
		BGEXPAND(defsim) = 'N';
		GAS_CONSTMU(defsim) = 'Y';
		HUBBLE(defsim) = 1;
		GAS_TYPE(defsim) = 'V';
		GAS_VISCOSITY(defsim) = 0;
		SIMMODEL(defsim) = BowShock;
		STATBOXSIZE(defsim) = 2;
		TOTMASS(defsim) = 2.7755E11*0.26*pow(STATBOXSIZE(defsim),3.L);
		AMAX(defsim) = 94;
		ANOW(defsim) = 1.;
		NX(defsim) = NY(defsim) = NZ(defsim) = 128;
		NXNY(defsim) = NX(defsim) * NY(defsim);
		DetermineFreeFallTime(defsim);
		BS_MU(defsim) = 1;
		BS_INITTEMP(defsim) = 1000;
		BS_HIGHDEN(defsim) = 10;
		BS_LOWDEN(defsim) = 10;
		BS_VELSHOCK(defsim) = 1;
		GRAVITY(defsim) = 'N';
		GAS_SFFLAG(defsim) = 'N';
		GAS_COOLFLAG(defsim) = 'N';
		GAS_SNFBFLAG(defsim) = 'N';
		GAS_BGHEATFLAG(defsim) = 'N';
		ASTEP(defsim) = 1./BS_VELSHOCK(defsim)/NX(defsim);
		NSTEP(defsim) = (AMAX(defsim)-ANOW(defsim))/ASTEP(defsim);
}

void mk_default_blast_param(SimParameters *defsim){
		SIMMODEL(defsim) = Blast; 
		GAS_TYPE(defsim) = 'V';
		GAS_VISCOSITY(defsim) = 0;
		BGEXPAND(defsim) = 'N';
		GAS_CONSTMU(defsim) = 'Y';
		HUBBLE(defsim) = 1;
		STATBOXSIZE(defsim) = 0.1;
		TOTMASS(defsim) = 2.7755E11*0.26*pow(STATBOXSIZE(defsim),3.L);
		AMAX(defsim) = 94.;
		ANOW(defsim) = 1.;
		ASTEP(defsim) = 0.0001;
		NSTEP(defsim) = 940000;
		NX(defsim) = NY(defsim) = NZ(defsim) = 256;
		NXNY(defsim) = NX(defsim) * NY(defsim);
		BLAST_MU(defsim) = 1;
		BLAST_INITTEMP(defsim) = 10;

		BLAST_EXPLODINGTEMP(defsim) = 1.e6;
		GRAVITY(defsim) = 'N';
		GAS_TYPE(defsim) = 'V';
		GAS_SFFLAG(defsim) = 'N';
		GAS_COOLFLAG(defsim) = 'N';
		GAS_SNFBFLAG(defsim) = 'N';
		GAS_BGHEATFLAG(defsim) = 'N';
}


void mk_default_param(SimParameters *defsim, char *cosmology){
	double zi;
	/* Simulation Parameters for WMAP 5-year */
	SIMMODEL(defsim) = Cosmos;
	GAS_TYPE(defsim) = 'N';
	EXTERNALFORCE(defsim) = 'N';
	PORDER(defsim) = 2;
	GRAVITY(defsim) = 'Y'; 
#ifdef PMonly
	GRV_EFOLD(defsim) = 0.5;
#else
	GRV_EFOLD(defsim) = 0.9;
#endif
    

	COSBOXSIZE(defsim) = 512;
	AMAX(defsim) = 48;
	ANOW(defsim) = 1.;
#ifdef GOTPM
	ASTEP(defsim) = 0.25;
	NSTEP(defsim) = 940;
#else
	ASTEP(defsim) = 0.5;
	NSTEP(defsim) = 94;
#endif
	NX(defsim) = NY(defsim) = NZ(defsim) = 512;
	NXNY(defsim) = NX(defsim)*NY(defsim);

	NWGROUP(defsim) = 4;
	ANIM_FLAG(defsim) = 'N';


	IndT_FLAG(defsim) = 'N';
	GAS_TYPE(defsim) = 'N';

	FOF_LINK(defsim) = 0.2;

	sprintf(GRAFIC_DIRECTORY(defsim),"./");
	/*
	SIM_LXMIN(defsim,dm) = SIM_LYMIN(defsim,dm) = SIM_LZMIN(defsim,dm) = 0;
	SIM_LXMAX(defsim,dm) = SIM_LYMAX(defsim,dm) = SIM_LZMAX(defsim,dm) = COSBOXSIZE(defsim);
	*/

	GAS_AlphaVis(defsim) = .2;
	GAS_BetaVis(defsim) = 0;
	GAS_ETAVIS(defsim) = 0;
	GAS_EPSVIS(defsim) = 0;
	GAS_VISCOSITY(defsim) = 0;
	GAS_COURANT(defsim) = 0.2;

	if(strcmp(cosmology,"WMAP3")==0 || strcmp(cosmology,"BWMAP3")==0){
		SIMMODEL(defsim) = Cosmos;
		BGEXPAND(defsim) = 'Y';
		OMEP(defsim) = 0.238;
		OMEPLAM(defsim) = 0.762;
		OMEPB(defsim) = 0.042;
		HUBBLE(defsim) = 0.732;
		NPOW(defsim) = 0.958;
		BIAS8(defsim) = 1.314;
		sprintf(POWFILENAME(defsim),"camb.z=47.dat");
		WLAM0(defsim) = -1;
		WLAM1(defsim) = 0;
#ifdef GOTPM
		IndT_FLAG(defsim) = 'N';
#else
		IndT_FLAG(defsim) = 'Y';
#endif
		if(strcmp(cosmology,"BWMAP3") !=0) GAS_TYPE(defsim) = 'N';
		else GAS_TYPE(defsim) = 'V';
	}
	else if(strcmp(cosmology, "WMAP5") ==0 || strcmp(cosmology, "BWMAP5") ==0 ) {
		SIMMODEL(defsim) = Cosmos;
		BGEXPAND(defsim) = 'Y';
		OMEP(defsim) = 0.26;
		OMEPLAM(defsim) = 0.74;
		OMEPB(defsim) = 0.044;
		HUBBLE(defsim) = 0.72;
		NPOW(defsim) = 0.96;
		BIAS8(defsim)= 1.26;
		sprintf(POWFILENAME(defsim),"camb.z=47.dat");
		WLAM0(defsim) = -1;
		WLAM1(defsim) = 0;
#ifdef GOTPM
		IndT_FLAG(defsim) = 'N';
#else
		IndT_FLAG(defsim) = 'Y';
#endif
		if(strcmp(cosmology, "BWMAP5") ==0) GAS_TYPE(defsim) = 'Y';
		else GAS_TYPE(defsim) = 'V';
	}
	else if(strcmp(cosmology,"NoExpand")==0 || strcmp(cosmology,"BNoExpand")==0){
		SIMMODEL(defsim) = Static;
		GAS_TYPE(defsim) = 'V';
		BGEXPAND(defsim) = 'N';
		OMEP(defsim) = 0.26;
		HUBBLE(defsim) = 0.732;
		GAS_CONSTMU(defsim) = 'Y';
		if(strcmp(cosmology,"BNoExpand")==0) GAS_TYPE(defsim) = 'Y';
		else GAS_TYPE(defsim) = 'N';
	}
	else if(strcmp(cosmology,"Blast")==0){
		mk_default_blast_param(defsim);
	}
	else if(strcmp(cosmology,"KH")==0){
		mk_default_kh_param(defsim);
	}
	else if(strcmp(cosmology,"BowShock")==0){
		mk_default_bowshock_param(defsim);
	}
	else if(strcmp(cosmology,"RT")==0){
		mk_default_rt_param(defsim);
	}
	else if(strcmp(cosmology,"Kepler")==0){
		mk_default_kepler_param(defsim);
	}
	else if(strcmp(cosmology,"MkGlass2D")==0){
		mk_default_2dglass_param(defsim);
	}
	else {
		exit(999);
	}

	if(GAS_TYPE(defsim) != 'N'){
		GAS_TYPE(defsim) = 'V';
		GAS_SFFLAG(defsim) = 'Y';
		GAS_SNFBFLAG(defsim) = 'Y';
		GAS_COOLFLAG(defsim) = 'Y';
		GAS_BGHEATFLAG(defsim) = 'Y';
		GAS_FNFLAG(defsim) = 'Y';
	
		GAS_DURANT(defsim) = 0.1;
		GAS_MINTEMP(defsim) = 10;
		GAS_CONSTMU(defsim) = 'N';
		SPH_TIMESTEPLIMITER(defsim) = 2;
		SPH_CONSTNEIGHPOW(defsim) = 32;
		SPH_INTERACTIONSPHERE(defsim) = 4;
	
		GAS_YP(defsim) = 0.25;
		GAS_MSTAR(defsim) = 1./3.;
		GAS_CSTAR(defsim) = 1./3.;
		SPH_NUMNEAR(defsim) = 30;
		GAS_TRAD0(defsim) = 2.725;
		GAS_GAMMA(defsim) = 1+ 2./3.;
		GAS_ALPHA(defsim) = 1.;
		GAS_COURANT(defsim) = 0.2;
		GAS_SFVIRDEN(defsim) = 57.7;
		GAS_SFGASDEN(defsim) = 0.1;
		GAS_SFTEMP(defsim) = 1.E4;
		GAS_UVSHIELDDEN(defsim) = 1.4E-2;
		GAS_DUVSHIELDDEN(defsim) = 1.4E-2;
		GAS_INITMETAL(defsim) = 1.E-4;
		SPH_INITFLAG(defsim) = 1;
	}
	FNL(defsim) = 0;
	GNL(defsim) = 0.;
	NSPACE(defsim) = 1;
	NXNY(defsim) = NX(defsim)*NY(defsim);
	ZINIT(defsim) = AMAX(defsim)/ANOW(defsim) - 1;
	THETA(defsim) = 0.3;
	RSPHERE(defsim) = 4;
	GRV_EPSILON(defsim) = 0.1;

	STEPCOUNT(defsim) = STEPNUM(defsim) = 1;
	FLAG_IC(defsim) = 0;
	sprintf(GRAFIC_DIRECTORY(defsim)," ");
	ISEED(defsim) = -56;
	NSKIP(defsim) = 0;
	PTYPE(defsim) = PMTYPE;
	POWREADFLAG(defsim) = 1;


#ifdef XYZDBL
	XYZSHIFT(defsim) = 1;
#else
	XYZSHIFT(defsim) = 0;
#endif
	SPH_NP(defsim) = 0;
	VORO_NP(defsim) = 0;
	STAR_NP(defsim) = 0;
	AGN_NP(defsim) = 0;
	DM_NP(defsim) = 0;
	NPSUM(defsim) = 0;

	if(IndT_FLAG(defsim) == 'Y'){
		IndT_NSUBSTEP(defsim) = 0;
		IndT_NSUBGAS(defsim) = 0;
		IndT_ISUBSTEP(defsim) = 0;
		IndT_NOWTSUBDIV(defsim) = IndT_MAXTSUBPOWER(defsim)  = 0;
	}


	sprintf(RV_FILE(defsim),"INITIAL");
	sprintf(RV_FILEPREFIX(defsim),"INITIAL");
	sprintf(GLACIALHEADER(defsim),"");
	float zp1 = ZINIT(defsim)+1;
	/*
	OMEI(defsim) = OMEP(defsim) *pow(zp1,3)/(OMEP(defsim)*pow(zp1,3) +
			OMEPLAM(defsim)+(1-OMEP(defsim)-OMEPLAM(defsim))*pow(zp1,2));
			*/
	float HofEz(SimParameters *, float);
	OMEI(defsim) = OMEP(defsim) *zp1*zp1*zp1/HofEz(defsim, ZINIT(defsim))/HofEz(defsim, ZINIT(defsim));
	ANIM_FLAG(defsim) = 0;
	PMPreFoFFLAG(defsim)= 'Y';
	printf("%s is detected\n",cosmology);
}

/*
#include "mpi.h"
void determine_mpi_misc_param(SimParameters *simpar){
	int nid,myid;
	float zi;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	simpar->nid = nid;
	simpar->myid = myid;
#ifdef XYZDBL
	simpar->xyzshiftflag = 1;
#else
	simpar->xyzshiftflag = 0;
#endif
	zi = simpar->amax-1;
	simpar->omei = simpar->omep*pow(1+zi,3)/(simpar->omep*pow(1+zi,3) +
			simpar->omeplam+(1-simpar->omep-simpar->omeplam)*pow(1+zi,2));
	simpar->mx = simpar->nx/simpar->nspace;
	simpar->my = simpar->ny/simpar->nspace;
	simpar->mz = simpar->nz/simpar->nspace;
	simpar->mxmy = simpar->mx * simpar->my;
	simpar->lnx = simpar->nx;
	simpar->lny = simpar->ny;
	simpar->lnz = simpar->nz;
	simpar->sphere_radius = 4;
	simpar->particle_radius = 4;
	simpar->rth = 8/(simpar->boxsize/simpar->nx);
	simpar->zinit = simpar->amax -1;
}
*/
