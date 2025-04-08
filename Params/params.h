void write_head(FILE *, SimParameters *);
void ReadSimulationParameters(FILE *, int *, SimParameters *);
void mk_default_param(SimParameters *, char *);
void CheckConsistencyParam(SimParameters *);
void write_head(FILE *, SimParameters *);
void DetermineFreeFallTime(SimParameters *);
void write_default_sim_parameter_file(FILE *, SimParameters *);
void write_sim_parameter_file(FILE *, SimParameters *);
void MkDefaultCosParam(SimParameters *, char *);
SimParameters readcheck_head(FILE *);
SimParameters  read_sim_parameter_file(SimParameters ,FILE *);
void  read_head(FILE*, SimParameters *);
void  read_slab_head(FILE*, SimParameters *);


#define Space ""

#define MAX_LINE_LENGTH 200
#define S_FLOAT   "%f"Space
#define S_DOUBLE   "%lg"Space
#define S_INT     "%d"Space
#define S_LONG     "%ld"Space
#define S_ULONG     "%lu"Space
#define S_LONG_LONG     "%lld"Space
#define S_STRING  "%s"Space
#define S_CHAR  "%c"Space
#define SET     "define "

/*
#define Cosmos 0
#define Galaxy 1
#define KH     3
#define RT     6
#define Kepler     7
*/

#ifdef GOTPM
#define P_Starting  "#Start of the Ascii Header of the GOTPM Simulation\n"
#else
#define P_Starting  "#Start of the Ascii Header of the EUNHA Simulation\n"
#endif

#define COMM_MODEL "# <<SIM_MODEL value>> 0(Cosmos), 1(Galaxy) 3(KH:2D), 4(Blast:3D),\n#    5(Bow Shock:3D), 6(RT:2D), 7(Keplyer:2D), 30(MkGlass2D)\n"
#define P_MODEL     SET"SIM_MODEL           = "S_INT" # SIMULATION MODEL \n"
#define P_GASTYPE   SET"flag for Hydro      = "S_CHAR" # N: no hydro, S: SPH hydro, V: Voro hydro\n"
#define P_INITIC    SET"INITIAL CONDITION   = "S_INT" # How to generate initial condtion, 1: Zeldovich, 2: 2nd order LPT\n"
#define P_EXTFORCE  SET"External Force      = "S_CHAR" # External Force N:No, Y:Yes\n"
#define P_BGCOSMOS  SET"BG_Cosmology        = "S_CHAR" # Background Expansion Y/N\n"
#define P_Nid       SET"Nid                 = "S_INT" # Number of Processors\n"
#define P_Myid      SET"Myid                = "S_INT" # Processor id\n"
#define P_Omep      SET"OmegaMatter0        = "S_FLOAT" # Current Omega Matter\n"
#define P_Omepb     SET"OmegaBaryon0        = "S_FLOAT" # Current Omega Baryon\n"
#define P_Omeplamb  SET"OmegaLambda0        = "S_FLOAT" # Current Omega Lambda\n"
#define P_Wlam0     SET"w0(eos of DE)       = "S_FLOAT" # w0 of Lambda\n"
#define P_Wlam1     SET"wa(eos of DE)       = "S_FLOAT" # w1 of Lambda\n"
#define P_Omei      SET"OmegaMatterI        = "S_FLOAT" # Initial Omega Matter\n"
#define P_fNL       SET"fNL                 = "S_FLOAT"\n"
#define P_gNL       SET"gNL                 = "S_FLOAT"\n"
#define P_Hubble    SET"Hubble              = "S_FLOAT" # Hubble expansion parameter\n"
#define P_nPS       SET"nPS                 = "S_FLOAT" # Spectral index\n"
#define P_Boxsize   SET"Boxsize(Mpc/h)      = "S_FLOAT" # Simulation Box size in Mpc/h\n"
#define P_Bias      SET"Bias                = "S_FLOAT" # Bias factor, inverse of sigma_8\n"
#define P_Amax      SET"Amax                = "S_FLOAT" # Maximum expansion factor/ final time\n"
#define P_Anow      SET"Anow                = "S_DOUBLE" # Current expansion factor/ current time\n"
#define P_Nx        SET"Nx                  = "S_LONG" # Number of grids along x axis\n"
#define P_Ny        SET"Ny                  = "S_LONG" # Number of grids along y axis\n"
#define P_Nz        SET"Nz                  = "S_LONG" # Number of grids along z axis\n"
#define P_Nspace    SET"Nspace              = "S_INT" # On every grids (1), on every other grids in 3d (2)\n"
#define P_Mx        SET"Mx                  = "S_LONG"\n"
#define P_My        SET"My                  = "S_LONG"\n"
#define P_Mz        SET"Mz                  = "S_LONG"\n"
#define P_MxMy      SET"MxMy                = "S_LONG"\n"
#define P_Lx        SET"Lx                  = "S_DOUBLE"\n"
#define P_Ly        SET"Ly                  = "S_DOUBLE"\n"
#define P_Lz        SET"Lz                  = "S_DOUBLE"\n"
#define P_Theta     SET"Theta               = "S_FLOAT" # Current cell opening angle in Tree correction\n"
#define P_EPSILON   SET"EPSILON             = "S_FLOAT" # Force Smoothing Length\n"
#define P_Rsphere   SET"Rsphere             = "S_FLOAT"\n"
#define P_Pradius   SET"Pradius             = "S_FLOAT"\n"
#define P_ken       SET"Ken                 = "S_FLOAT"\n"
#define P_ktot      SET"Ktot                = "S_FLOAT"\n"
#define P_const0    SET"Const0              = "S_FLOAT"\n"
#define P_poten0    SET"Poten0              = "S_FLOAT"\n"
#define P_fact1     SET"Fact1               = "S_FLOAT"\n"
#define P_fact2     SET"Fact2               = "S_FLOAT"\n"
#define P_pfact     SET"Pfact               = "S_FLOAT"\n"
#define P_stepnum   SET"Stepnum             = "S_INT"\n"
#define P_nskip     SET"Nskip               = "S_INT"\n"
#define P_local_nz  SET"Local_nz            = "S_LONG"\n"
#define P_local_zs  SET"Local_z_start       = "S_LONG"\n"
#define P_Zinit     SET"Z_init              = "S_FLOAT"\n"
#define P_Rsmooth   SET"Rsmooth             = "S_FLOAT"\n"
#define P_Rthooth   SET"Rth                 = "S_FLOAT"\n"
#define P_Astep     SET"Astep               = "S_DOUBLE"\n"
#define P_Nstep     SET"Nstep               = "S_INT" # Number of time steps\n"
#define P_nsub      SET"Nsub                = "S_INT"\n"
#define P_isub      SET"Isub                = "S_INT"\n"
#define P_Istep     SET"Stepcount           = "S_INT" # Current time step number\n"
#define P_Iseed     SET"Iseed               = "S_LONG"\n"
#define P_flagIC    SET"flagIC              = "S_INT" #0:self generate, 1:grafic1 read\n"
#define P_grafic    SET"grafic1 dir. name   = "S_STRING"\n"
#define P_rvfile    SET"rvfilename          = "S_STRING"\n"
#define P_rvprefix  SET"rvprefix            = "S_STRING"\n"
#define P_npsum     SET"Np                  = "S_LONG"\n"
#define P_npdm      SET"DNp                 = "S_LONG"\n"
#define P_npbm      SET"BNp                 = "S_LONG"\n"
#define P_npstar    SET"SNp                 = "S_LONG"\n"
#define P_npagn     SET"ANp                 = "S_LONG"\n"
#define P_PTYPE     SET"Particle Type       = "S_INT"\n"
#define P_glacial   SET"GlacialFileName     = "S_STRING"\n"

#ifdef XYZDBL
#define P_Xmin      SET"Xmin                = "S_DOUBLE" # Bottom x value of local domain slab\n"
#define P_Xmax      SET"Xmax                = "S_DOUBLE" # Top x value of local domain slab\n"
#define P_Ymin      SET"Ymin                = "S_DOUBLE" # Bottom y value of local domain slab\n"
#define P_Ymax      SET"Ymax                = "S_DOUBLE" # Top y value of local domain slab\n"
#define P_Zmin      SET"Zmin                = "S_DOUBLE" # Bottom z value of local domain slab\n"
#define P_Zmax      SET"Zmax                = "S_DOUBLE" # Top z value of local domain slab\n"
#else
#define P_Xmin      SET"Xmin                = "S_FLOAT" # Bottom x value of local domain slab\n"
#define P_Xmax      SET"Xmax                = "S_FLOAT" # Top x value of local domain slab\n"
#define P_Ymin      SET"Ymin                = "S_FLOAT" # Bottom y value of local domain slab\n"
#define P_Ymax      SET"Ymax                = "S_FLOAT" # Top y value of local domain slab\n"
#define P_Zmin      SET"Zmin                = "S_FLOAT" # Bottom z value of local domain slab\n"
#define P_Zmax      SET"Zmax                = "S_FLOAT" # Top z value of local domain slab\n"
#endif

#define P_NumNear   SET"Num_nearby_SPH      = "S_INT" # Number of neighbor particles (SPH)\n"
#define P_Yp        SET"Yp_SPH              = "S_FLOAT" # Mass ratio of H to He (SPH)\n"
#define P_Trad0     SET"Trad0_SPH           = "S_FLOAT" # Current background radiation temperature (SPH)\n"
#define P_Gamma     SET"Gamma_hydro         = "S_FLOAT" # Polytropic index (hydro)\n"
#define P_Alpha     SET"Alpha_SPH           = "S_FLOAT"\n"
#define P_Courant   SET"Courant_hydro       = "S_FLOAT" # Courant number in hydro\n"
#define P_Powerflag SET"Powerflag           = "S_INT" #0: self generate, 1: read oCAMB, 2: read nCAMB\n"
#define P_Powfile   SET"Powerfile           = "S_STRING" # 0: No input 1: CAMB 2: Ascii\n"
#define P_Apowfile  SET"AsciiPowerfile      = "S_STRING"\n"
#define P_PMPreFoF  SET"PM_PreFoF_flag      = "S_CHAR"\n"
#define P_XYZSHIFT  SET"XYZ_Shift_Flag      = "S_INT"\n"
#define P_PTSIZE    SET"Particle_TYPE_Size  = "S_ULONG"\n"
#define P_Viewer    SET"ViewerTrackFile     = "S_STRING"\n"
#define P_DDINFOSIZE  SET"DDINFO SIZE         = "S_INT"\n"
/* individual time step */
#define P_INDTFLAG  SET"Ind Time Flag       = "S_CHAR"\n"
#define P_INDABEFO  SET"Ind A Before        = "S_DOUBLE"\n"
#define P_INDASTEP  SET"Ind astep           = "S_DOUBLE"\n"
#define P_INDDA     SET"Ind da              = "S_DOUBLE"\n"
#define P_INDNTSUB  SET"Ind nowTsubdiv      = "S_INT"\n"
#define P_INDMTSUB  SET"Ind maxTsubpower    = "S_INT"\n"
#define P_INDPreFoF SET"Ind flag 4 PreFoF   = "S_CHAR"\n"
#define P_INDsyncp  SET"Ind flag 4 syncp    = "S_CHAR"\n"
#define P_INDnsubstepcount  SET"Ind nsubstepcount   = "S_INT"\n"
#define P_INDtnumcount  SET"Ind tnumcount       = "S_INT"\n"
#define P_INDanext  SET"Ind anext           = "S_DOUBLE"\n"

/* Substep control */
#define P_SUBiflagfixed  SET"SubStep fixed Flag  = "S_CHAR"\n"
#define P_SUBsph    SET"SubStep sph power   = "S_INT"\n"
#define P_SUBnbody  SET"SubStep Nbody power = "S_INT"\n"
#define P_SUBfixed  SET"SubStep Fixed power = "S_INT"\n"

/* Flags for the sph */
#define P_SPH       SET"Flag for SPH        = "S_CHAR" # flag for basic SPH (Y/N)\n"
#define P_CONSTMU   SET"Flag CONST MU       = "S_CHAR" # flag for Constant Mu in sph  (Y/N)\n"
#define P_STARFORM  SET"FlagStarformation   = "S_CHAR" # flag for star formation in SPH (Y/N)\n"
#define P_SNFB      SET"FlagSNFeedBack      = "S_CHAR" # flag for supernovae feedback in SPH (Y/N)\n"
#define P_COOL      SET"FlagCooling         = "S_CHAR" # flag for cooling in SPH (Y/N)\n"
#define P_BGHEAT    SET"FlagBackgroundHeat  = "S_CHAR" # flag for background heating in SPH (Y/N)\n"

/* SPH parameters */
#define P_SFVIRIAL  SET"SF_Virial_Den       = "S_FLOAT" # Cosmic virial density criteria for SF\n"
#define P_SFGASDEN  SET"SF_Gas_Den          = "S_FLOAT" # Gas density criterion for SF\n"
#define P_SFTEMP    SET"SF_Temperature      = "S_FLOAT" # Temperature criterion for SF\n"
#define P_MSTAR     SET"SF_Mstar            = "S_FLOAT" # Maximum Number of stars per gas particle (SF)\n"
#define P_CSTAR     SET"SF_Cstar            = "S_FLOAT" # Star formation efficiency (SF)\n"
#define P_MINTEMP   SET"MIN_TEMP_SPH        = "S_FLOAT" # Minimum temperature of gas particle\n"
#define P_HEATDEN   SET"UV Shielding GasDen = "S_FLOAT" # UV shielding gas density\n"
#define P_DHEATDEN  SET"UV Shielding dGasDen= "S_FLOAT" # Transition Width of UV shielding gas density\n"
#define P_IMETAL    SET"Init Metallicity    = "S_FLOAT" # Initial metallicity\n"


#define P_TSLIMITER SET"Time Step Limiter   = "S_INT" # Step number in the timestep limiter\n"

#define P_TOTALMASS SET"TMass_Box_inMsun    = "S_FLOAT" # Total mass in the simulation box (No expansion only)\n"
#define P_FFTIME    SET"Free Fall Time      = "S_DOUBLE" # Free Fall time in sec. (No expansion only)\n"
#define P_CONNEIGH  SET"Pow. of Const Neigh = "S_INT" # subPower of Constant Neighbor Approximation\n"

#define P_FIXNEIGH  SET"Fixed Neighbor List = "S_CHAR" # Fixed Neighbor List (Y/N)\n"
#define P_DURANT    SET"Durant Factor       = "S_FLOAT" # Durant factor (float)\n"

#define P_GRAVITY   SET"Self Gravity        = "S_CHAR" # Whether to include the self gravity\n"
/* Blast parameters */
#define P_BlastInitTemp  SET"Blast Init Temp     = "S_FLOAT" # Initial Temperature of the Blast Test\n"
#define P_BlastMU   SET"Blast Constat MU    = "S_FLOAT" # Constant value of MU of the Blast Test\n"
#define P_BlastTemp SET"Blast Exploding Temp= "S_FLOAT" # Initial Temperature of the Explosion of the Blast Test\n"
/* KH parameters */
#define P_KH_Rho1      SET"KH Den in R1     = "S_FLOAT" # Initial Density in Region 1\n"
#define P_KH_Rho2      SET"KH Den in R2     = "S_FLOAT" # Initial Density in Region 2\n"
#define P_KH_Vel1      SET"KH Vel in R1     = "S_FLOAT" # Initial Velosity in Region 1\n"
#define P_KH_Vel2      SET"KH Vel in R2     = "S_FLOAT" # Initial Velosity in Region 2\n"
#define P_KH_Deltay    SET"KH Transition size  = "S_FLOAT" # Initial Transition region of the KH Test\n"
#define P_KH_Pressure  SET"KH Global Pressure  = "S_FLOAT" # Initial global pressure of the KH Test\n"
#define P_KH_Vperturb  SET"KH Vel Perturb Amp  = "S_FLOAT" # Initial Perturbation of Vel in Y-direction\n"
#define P_KH_Xmin      SET"KH Xmin             = "S_DOUBLE" # Minimum X of the simulation space\n"
#define P_KH_Ymin      SET"KH Ymin             = "S_DOUBLE" # Minimum Y of the simulation space\n"
#define P_KH_Xmax      SET"KH Xmax             = "S_DOUBLE" # Maximum X of the simulation space\n"
#define P_KH_Ymax      SET"KH Ymax             = "S_DOUBLE" # Maximum Y of the simulation space\n"
#define P_KH_GridSize  SET"KH Grid Size        = "S_DOUBLE" # Grid size for Linked List\n"
#define P_KH_OrderAcc  SET"KH Order of Accuracy= "S_FLOAT" # Order of Accuracy in Voro mode. 0<=w<0.5\n"
/* BlowShock parameters */
#define P_BSInitTemp  SET"BS    Init Temp     = "S_FLOAT" # Initial Temperature of the BowShock Test\n"
#define P_BSMU  	SET"BS    Constat MU    = "S_FLOAT" # Constant value of MU of the BowShcok Test\n"
#define P_BSHden 	SET"BS High Density     = "S_FLOAT" # High Density in the BowShcok Test\n"
#define P_BSLden 	SET"BS Low Density      = "S_FLOAT" # Low Density in the BowShcok Test\n"
#define P_BSVel  	SET"BS Supersonic Vel.  = "S_FLOAT" # Supersonic Velocity in the BowShcok Test\n"

/* Rayleigh Taylor Instibility parameters */
#define P_RT_XMAX  SET"RT maximum x (box)  = "S_DOUBLE" # Maximum X in the RT Test\n"
#define P_RT_YMAX   SET"RT maximum y (box)  = "S_DOUBLE" # Maximum Y in the RT Test\n"
#define P_RT_DEN1 SET"RT Lower density     = "S_FLOAT" # Low Density in the RT Test\n"
#define P_RT_DEN2 SET"RT Higher density         = "S_FLOAT" # High Density in the RT Test\n"
#define P_RT_extG  SET"RT External Force   = "S_FLOAT" # External Gravity in ydirection in the RT Test\n"
#define P_RT_OrderAcc  SET"RT Order of Accuracy= "S_FLOAT" # Order of Accuracy in Voro mode. 0<=w<0.5\n"
#define P_RT_GridSize  SET"RT Grid Size        = "S_DOUBLE" # Grid size for Linked List\n"
#define P_RT_Deltay    SET"RT Transition size  = "S_FLOAT" # Initial Transition region of the RT Test\n"
#define P_RT_Vperturb  SET"RT Vel Perturb Amp  = "S_FLOAT" # Initial Perturbation ampl. of Vel. in Y-direction\n"
#define P_RT_Phalf  SET"RT Pressure at Ly/2  = "S_FLOAT" # Hydrostatic pressure at Ly/2\n"

#define P_KP_EPS    SET"Keplyer EPS of g    = "S_FLOAT" # g-epsilon\n"
#define P_KP_XMAX   SET"Kepler xmax (box)   = "S_DOUBLE" # Maximum X in the Kepler Test\n"
#define P_KP_YMAX   SET"Kepler ymax (box)   = "S_DOUBLE" # Maximum Y in the Kepler Test\n"
#define P_KP_OrderAcc  SET"Kepler Order of Accuracy = "S_FLOAT" # Order of Accuracy in Voro mode. 0<=w<0.5\n"
#define P_KP_GridSize  SET"Kepler Grid Size = "S_DOUBLE" # Grid size for Linked List\n"

#define P_GL2D_XMAX   SET"MkGL2D xmax (box)   = "S_DOUBLE" # Maximum X in makding 2D glass\n"
#define P_GL2D_YMAX   SET"MkGL2D ymax (box)   = "S_DOUBLE" # Maximum Y in making 2D glass\n"
#define P_GL2D_OrderAcc  SET"MkGL2D Order of Accuracy = "S_FLOAT" # Order of Accuracy in Voro mode. 0<=w<0.5\n"
#define P_GL2D_GridSize  SET"MkGL2D Grid Size    = "S_DOUBLE" # Grid size for Linked List\n"
#define P_GL2D_Kappa     SET"MkGL2D Kappa        = "S_FLOAT" # coefficient to determine w2 in the 2d glass\n"


#define P_VORO_AlphaVis  SET"Alpha parameter in Voro  = "S_FLOAT" # Alpha factor of Voronoi AV\n"
#define P_VORO_BetaVis  SET"Beta parameter in Voro  = "S_FLOAT" # Beta factor of Voronoi AV\n"
#define P_VORO_EtaVis  SET"Eta parameter in Voro  = "S_FLOAT" # Eta factor of Voronoi AV\n"
#define P_VORO_EpsVis  SET"Eps parameter in Voro  = "S_FLOAT" # Eps factor of Voronoi AV\n"
#define P_VORO_Viscosity  SET"Re parameter in Voro   = "S_FLOAT" # Reynolds number of Voronoi\n"

#define P_WGRPS  SET"FILE I/O Group Size = "S_INT" # FILE I/O Group Size\n"

#define P_Closing   "#End of Ascii Header\n"
#define P_NULL      "##################################################\n"
#define P_VOROVISINFO      "# INFORMATION ON THE VORONOI VISCOSITY (AV:artificial viscosity)\n"

#define MPI_DATA(frw,wp,sp,simpar)do{\
 	ncnt += frw(wp,P_Myid,sp MYID(simpar));\
	ncnt += frw(wp,P_Nid,sp NID(simpar));\
}while(0)
#define ENERGY(frw,wp,sp,simpar){\
}
#define CORE(frw,wp,sp,simpar) do{\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,COMM_MODEL);\
	ncnt += frw(wp,P_MODEL,sp SIMMODEL(simpar));\
	ncnt += frw(wp,P_GASTYPE,sp GAS_TYPE(simpar));\
	/*\
	ncnt += frw(wp,P_EXTFORCE,sp EXTERNALFORCE(simpar));\
	*/\
	ncnt += frw(wp,P_GRAVITY,sp GRAVITY(simpar));\
	if(SIMMODEL(simpar) == Cosmos){\
		ncnt += frw(wp,P_flagIC,sp FLAG_IC(simpar));\
		ncnt += frw(wp,P_grafic,   GRAFIC_DIRECTORY(simpar));\
		ncnt += frw(wp,P_INITIC,sp PORDER(simpar));\
		if(POWREADFLAG(simpar)!=0){\
			ncnt += frw(wp,P_Powfile,   POWFILENAME(simpar));\
			if(POWREADFLAG(simpar)==2) ncnt += frw(wp,P_Apowfile,   INPAPKFILENAME(simpar));\
		}\
		ncnt += frw(wp,"### Hubble parameter is in 100km/sec/Mpc.\n");\
		ncnt += frw(wp,P_Hubble,sp HUBBLE(simpar));\
		ncnt += frw(wp,P_Omep,sp OMEP(simpar));\
		ncnt += frw(wp,P_Omepb,sp OMEPB(simpar));\
		ncnt += frw(wp,P_Omeplamb,sp OMEPLAM(simpar));\
		ncnt += frw(wp,P_Wlam0,sp WLAM0(simpar));\
		ncnt += frw(wp,P_Wlam1,sp WLAM1(simpar));\
		ncnt += frw(wp,P_fNL,sp FNL(simpar));\
		ncnt += frw(wp,P_gNL,sp GNL(simpar));\
		ncnt += frw(wp,P_NULL);\
		ncnt += frw(wp,"### nPS is the power spectral index.\n");\
		ncnt += frw(wp,"### Bias factor is inverse of sigma_8.\n");\
		ncnt += frw(wp,P_nPS,sp NPOW(simpar));\
		ncnt += frw(wp,P_Bias,sp BIAS8(simpar));\
		ncnt += frw(wp,P_Iseed,sp ISEED(simpar));\
		ncnt += frw(wp,P_Powerflag,sp POWREADFLAG(simpar));\
		ncnt += frw(wp,P_Boxsize,sp COSBOXSIZE(simpar));\
		ncnt += frw(wp,P_NULL);\
		ncnt += frw(wp,P_Amax,sp AMAX(simpar));\
		ncnt += frw(wp,P_Anow,sp ANOW(simpar));\
		ncnt += frw(wp,P_Astep,sp ASTEP(simpar));\
		ncnt += frw(wp,P_Theta,sp THETA(simpar));\
		ncnt += frw(wp,P_EPSILON,sp GRV_EPSILON(simpar));\
	}\
	else if(SIMMODEL(simpar) == Static){\
		ncnt += frw(wp,P_Boxsize,sp STATBOXSIZE(simpar));\
		ncnt += frw(wp,P_TOTALMASS,sp TOTMASS(simpar));\
		ncnt += frw(wp,P_FFTIME,sp FFTIME(simpar));\
	}\
	else if(SIMMODEL(simpar)  == Blast){\
		ncnt += frw(wp,P_Courant,sp GAS_COURANT(simpar));\
		ncnt += frw(wp,P_Gamma,sp GAS_GAMMA(simpar));\
		ncnt += frw(wp,P_BlastInitTemp,sp BLAST_INITTEMP(simpar));\
		ncnt += frw(wp,P_BlastMU,sp GAS_MU(simpar));\
		ncnt += frw(wp,P_BlastTemp,sp BLAST_EXPLODINGTEMP(simpar));\
	}\
	else if(SIMMODEL(simpar) == KH){\
		ncnt += frw(wp,P_Courant,sp GAS_COURANT(simpar));\
		ncnt += frw(wp,P_Gamma,sp GAS_GAMMA(simpar));\
		ncnt += frw(wp,P_KH_Rho1,sp KH_Rho1(simpar));\
		ncnt += frw(wp,P_KH_Rho2,sp KH_Rho2(simpar));\
		ncnt += frw(wp,P_KH_Vel1,sp KH_Vel1(simpar));\
		ncnt += frw(wp,P_KH_Vel2,sp KH_Vel2(simpar));\
		ncnt += frw(wp,P_KH_Deltay,sp KH_Deltay(simpar));\
		ncnt += frw(wp,P_KH_Pressure,sp KH_Pressure(simpar));\
		ncnt += frw(wp,P_KH_Vperturb,sp KH_Vperturb(simpar));\
		ncnt += frw(wp,P_KH_Xmax,sp KH_XMAX(simpar));\
		ncnt += frw(wp,P_KH_Ymax,sp KH_YMAX(simpar));\
		ncnt += frw(wp,P_KH_Xmin,sp KH_XMIN(simpar));\
		ncnt += frw(wp,P_KH_Ymin,sp KH_YMIN(simpar));\
		ncnt += frw(wp,P_KH_GridSize,sp KH_GridSize(simpar));\
		ncnt += frw(wp,P_KH_OrderAcc,sp KH_OA(simpar));\
	}\
	else if(SIMMODEL(simpar) == BowShock){\
		ncnt += frw(wp,P_Courant,sp GAS_COURANT(simpar));\
		ncnt += frw(wp,P_Gamma,sp GAS_GAMMA(simpar));\
		ncnt += frw(wp,P_BSInitTemp,sp BS_INITTEMP(simpar));\
		ncnt += frw(wp,P_BSMU,sp GAS_MU(simpar));\
		ncnt += frw(wp,P_BSHden,sp BS_HIGHDEN(simpar));\
		ncnt += frw(wp,P_BSLden,sp BS_LOWDEN(simpar));\
		ncnt += frw(wp,P_BSVel,sp BS_VELSHOCK(simpar));\
	}\
	else if(SIMMODEL(simpar) == RT){\
		ncnt += frw(wp,P_Courant,sp GAS_COURANT(simpar));\
		ncnt += frw(wp,P_Gamma,sp GAS_GAMMA(simpar));\
		ncnt += frw(wp,P_RT_XMAX,sp RT_XMAX(simpar));\
		ncnt += frw(wp,P_RT_YMAX,sp RT_YMAX(simpar));\
		ncnt += frw(wp,P_RT_DEN1,sp RT_DEN1(simpar));\
		ncnt += frw(wp,P_RT_DEN2,sp RT_DEN2(simpar));\
		ncnt += frw(wp,P_RT_extG,sp RT_ACC(simpar));\
		ncnt += frw(wp,P_RT_GridSize,sp RT_GridSize(simpar));\
		ncnt += frw(wp,P_RT_OrderAcc,sp RT_OA(simpar));\
		ncnt += frw(wp,P_RT_Deltay,sp RT_Deltay(simpar));\
		ncnt += frw(wp,P_RT_Vperturb,sp RT_Vperturb(simpar));\
		ncnt += frw(wp,P_RT_Phalf,sp RT_Phalf(simpar));\
	}\
	else if(SIMMODEL(simpar) == Kepler){\
		ncnt += frw(wp,P_Courant,sp GAS_COURANT(simpar));\
		ncnt += frw(wp,P_Gamma,sp GAS_GAMMA(simpar));\
		ncnt += frw(wp,P_KP_XMAX,sp KP_XMAX(simpar));\
		ncnt += frw(wp,P_KP_YMAX,sp KP_YMAX(simpar));\
		ncnt += frw(wp,P_KP_EPS,sp KP_EPS(simpar));\
		ncnt += frw(wp,P_KP_GridSize,sp KP_GridSize(simpar));\
		ncnt += frw(wp,P_KP_OrderAcc,sp KP_OA(simpar));\
	}\
	else if(SIMMODEL(simpar) == MkGlass2D){\
		ncnt += frw(wp,P_Courant,sp GAS_COURANT(simpar));\
		ncnt += frw(wp,P_Gamma,sp GAS_GAMMA(simpar));\
		ncnt += frw(wp,P_GL2D_XMAX,sp GL2D_XMAX(simpar));\
		ncnt += frw(wp,P_GL2D_YMAX,sp GL2D_YMAX(simpar));\
		ncnt += frw(wp,P_GL2D_GridSize,sp GL2D_GridSize(simpar));\
		ncnt += frw(wp,P_GL2D_OrderAcc,sp GL2D_OA(simpar));\
		ncnt += frw(wp,P_GL2D_Kappa,sp GL2D_Kappa(simpar));\
	}\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,"### Nx, Ny, and Nz should be same & nspace should be int.\n");\
	ncnt += frw(wp,P_Nx,sp NX(simpar));\
	ncnt += frw(wp,P_Ny,sp NY(simpar));\
	ncnt += frw(wp,P_Nz,sp NZ(simpar));\
	ncnt += frw(wp,P_NULL);\
	if(ANIM_FLAG(simpar) == 'Y') {\
		ncnt += frw(wp,P_Viewer,   ANIM_VIEWFILE(simpar));\
	}\
	if(SIMMODEL(simpar) == Cosmos){\
		ncnt += frw(wp,P_NULL);\
		ncnt += frw(wp,P_Nstep,sp NSTEP(simpar));\
		ncnt += frw(wp,P_Istep,sp STEPCOUNT(simpar));\
		ncnt += frw(wp,P_stepnum,sp STEPNUM(simpar));\
		ncnt += frw(wp,P_NULL);\
		ncnt += frw(wp,P_rvfile,   RV_FILE(simpar));\
		ncnt += frw(wp,P_rvprefix,   RV_FILEPREFIX(simpar));\
		ncnt += frw(wp,P_PMPreFoF, sp  PMPreFoFFLAG(simpar));\
		ncnt += frw(wp,P_NULL);\
		ncnt += frw(wp,P_glacial,   GLACIALHEADER(simpar));\
		ncnt += frw(wp,P_npsum,sp NPSUM(simpar));\
		ncnt += frw(wp,P_npdm,sp DM_NP(simpar));\
		ncnt += frw(wp,P_npbm,sp SPH_NP(simpar));\
		ncnt += frw(wp,P_npstar,sp STAR_NP(simpar));\
		ncnt += frw(wp,P_npagn,sp AGN_NP(simpar));\
		ncnt += frw(wp,P_PTYPE,sp PTYPE(simpar));\
		ncnt += frw(wp,P_INDTFLAG, sp IndT_FLAG(simpar));\
		{ \
			int ddinfosize=sizeof(DoDeInfo);\
			ncnt += frw(wp,P_DDINFOSIZE,sp ddinfosize);\
			if(ddinfosize != sizeof(DoDeInfo)){\
				DEBUGPRINT("P%d has error in the size of DoDeInfo %d : %ld\n",MYID(simpar), ddinfosize, sizeof(DoDeInfo));\
				exit(999);\
			}\
		}\
	}\
	if(IndT_FLAG(simpar)=='Y')\
	{\
		ncnt += frw(wp,P_nsub,sp IndT_NSUBSTEP(simpar));\
		ncnt += frw(wp,P_isub,sp IndT_ISUBSTEP(simpar));\
		ncnt += frw(wp,"################ INDIVIDUAL TIMESTEP  ############\n");\
		ncnt += frw(wp,P_INDABEFO,sp IndT_ABEFORE(simpar));\
		ncnt += frw(wp,P_INDASTEP,sp IndT_ASTEP(simpar));\
		ncnt += frw(wp,P_INDDA,sp IndT_DA(simpar));\
		ncnt += frw(wp,P_INDNTSUB,sp IndT_NOWTSUBDIV(simpar));\
		ncnt += frw(wp,P_INDMTSUB,sp IndT_MAXTSUBPOWER(simpar));\
		ncnt += frw(wp,P_INDPreFoF,sp IndT_IFLAGPREFOF(simpar));\
		ncnt += frw(wp,P_INDsyncp,sp IndT_IFLAGSYNCPDATA(simpar));\
		ncnt += frw(wp,P_INDnsubstepcount,sp IndT_NSUBSTEPCOUNT(simpar));\
		ncnt += frw(wp,P_INDtnumcount,sp IndT_TNUMCOUNT(simpar));\
		ncnt += frw(wp,P_INDanext,sp IndT_ANEXT(simpar));\
		ncnt += frw(wp,"################ SubStep Power  ############\n");\
		ncnt += frw(wp,P_SUBiflagfixed,sp IndT_IFLAGFIXEDLIST(simpar));\
		ncnt += frw(wp,P_SUBsph,sp IndT_NSUBGAS(simpar));\
		ncnt += frw(wp,P_SUBnbody,sp IndT_NSUBNBODY(simpar));\
		ncnt += frw(wp,P_SUBfixed,sp IndT_NSUBFIXED(simpar));\
	}\
	if(SIMMODEL(simpar) == Cosmos){\
		ncnt += frw(wp,"################    GAS FLAGS   ##################\n");\
		ncnt += frw(wp,P_GASTYPE,sp GAS_TYPE(simpar));\
		ncnt += frw(wp,P_WGRPS,sp WGROUPSIZE(simpar));\
	}\
	if(GAS_TYPE(simpar)=='S'){\
		ncnt += frw(wp,P_CONSTMU,sp GAS_CONSTMU(simpar));\
		ncnt += frw(wp,P_STARFORM,sp GAS_SFFLAG(simpar));\
		ncnt += frw(wp,P_COOL,sp GAS_COOLFLAG(simpar));\
		ncnt += frw(wp,P_SNFB,sp GAS_SNFBFLAG(simpar));\
		ncnt += frw(wp,P_BGHEAT,sp GAS_BGHEATFLAG(simpar));\
		ncnt += frw(wp,P_FIXNEIGH ,sp GAS_FNFLAG(simpar));\
		ncnt += frw(wp,P_DURANT ,sp GAS_DURANT(simpar));\
		ncnt += frw(wp,"################ SPH parameters ##################\n");\
		ncnt += frw(wp,P_Yp,sp GAS_YP(simpar));\
		ncnt += frw(wp,P_NumNear,sp SPH_NUMNEAR(simpar));\
		ncnt += frw(wp,P_Trad0,sp GAS_TRAD0(simpar));\
		ncnt += frw(wp,P_Gamma,sp GAS_GAMMA(simpar));\
		ncnt += frw(wp,P_Alpha,sp GAS_ALPHA(simpar));\
		ncnt += frw(wp,P_Courant,sp GAS_COURANT(simpar));\
		ncnt += frw(wp,P_HEATDEN,sp GAS_UVSHIELDDEN(simpar));\
		ncnt += frw(wp,P_DHEATDEN,sp GAS_DUVSHIELDDEN(simpar));\
		ncnt += frw(wp,P_IMETAL,sp GAS_INITMETAL(simpar));\
		ncnt += frw(wp,P_MINTEMP,sp GAS_MINTEMP(simpar));\
		ncnt += frw(wp,P_CONNEIGH,sp SPH_CONSTNEIGHPOW(simpar));\
		ncnt += frw(wp,P_NULL);\
		if(GAS_SFFLAG(simpar) == 'Y') {\
			ncnt += frw(wp,P_SFVIRIAL,sp GAS_SFVIRDEN(simpar));\
			ncnt += frw(wp,P_SFGASDEN,sp GAS_SFGASDEN(simpar));\
			ncnt += frw(wp,P_SFTEMP,sp GAS_SFTEMP(simpar));\
			ncnt += frw(wp,P_MSTAR ,sp GAS_MSTAR(simpar));\
			ncnt += frw(wp,P_CSTAR ,sp GAS_CSTAR(simpar));\
		}\
		ncnt += frw(wp,P_TSLIMITER ,sp SPH_TIMESTEPLIMITER(simpar));\
	}\
	else if(GAS_TYPE(simpar)=='V'){\
		ncnt += frw(wp,P_Courant,sp GAS_COURANT(simpar));\
		ncnt += frw(wp,P_VOROVISINFO);\
		ncnt += frw(wp,P_VORO_AlphaVis,sp GAS_AlphaVis(simpar));\
		ncnt += frw(wp,P_VORO_BetaVis,sp GAS_BetaVis(simpar));\
		ncnt += frw(wp,P_VORO_EtaVis,sp GAS_ETAVIS(simpar));\
		ncnt += frw(wp,P_VORO_EpsVis,sp GAS_EPSVIS(simpar));\
		ncnt += frw(wp,P_VORO_Viscosity,sp GAS_VISCOSITY(simpar));\
	}\
	ncnt += frw(wp,P_NULL);\
}while(0)
#define AUX1(frw,wp,sp,simpar) do{\
	ncnt += frw(wp,"### Miscellaneous parameters\n");\
	ncnt += frw(wp,P_Omei,sp OMEI(simpar));\
}while(0)
#define DEFAULT_PARAMS(frw,wp,sp,simpar) do{\
	ncnt += frw(wp,P_Starting);\
	CORE(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_Closing);\
}while(0)

#define VIEWER(frw,wp,sp,simpar,ncnt)do{\
	ncnt += frw(wp,P_Viewer,    ANIM_VIEWFILE(simpar));\
}while(0)
#define VIEWER2(frw,wp,sp,simpar,ncnt){\
	ncnt += frw(wp,P_Viewer,    ANIM_VIEWFILE((simpar)));\
}


#define PARAMS(frw,wp,sp,simpar) do{\
	ncnt += frw(wp,P_Starting);\
	CORE(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_XYZSHIFT,sp XYZSHIFT(simpar));\
	PTYPESIZE(simpar) = sizeof(dmparticletype);\
	ncnt += frw(wp,P_PTSIZE,sp PTYPESIZE(simpar));\
	ncnt += frw(wp,P_Closing);\
}while(0)

#define FILE_HEADER(frw,wp,sp,simpar) do{\
	ncnt += frw(wp,P_Starting);\
	ncnt += frw(wp,P_NULL);\
	MPI_DATA(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_NULL);\
	CORE(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_XYZSHIFT,sp XYZSHIFT(simpar));\
	PTYPESIZE(simpar) = sizeof(dmparticletype);\
	ncnt += frw(wp,P_PTYPE,sp PTYPE(simpar));\
	AUX1(frw,wp,sp,simpar);\
	ENERGY(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_fact1,sp Evol_FACT1(simpar));\
	ncnt += frw(wp,P_fact2,sp Evol_FACT2(simpar));\
	ncnt += frw(wp,P_pfact,sp Evol_PFACT(simpar));\
	ncnt += frw(wp,P_local_nz,sp LOCAL_NZ(simpar));\
	ncnt += frw(wp,P_local_zs,sp LOCAL_Z_START(simpar));\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,P_Closing);\
} while(0)

#define FILE_SLAB_HEADER(frw,wp,sp,simpar) do{\
	ncnt += frw(wp,P_Starting);\
	ncnt += frw(wp,P_NULL);\
	MPI_DATA(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_NULL);\
	CORE(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_NULL);\
	AUX1(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_Zmin,sp ZMIN(simpar));\
	ncnt += frw(wp,P_Zmax,sp ZMAX(simpar));\
	ncnt += frw(wp,P_XYZSHIFT,sp XYZSHIFT(simpar));\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,P_Closing);\
}while(0)


#define FILE_HEADER_ZMINZMAX(frw,wp,sp,simpar,zmin,zmax) do{\
	ncnt += frw(wp,P_Starting);\
    ncnt += frw(wp,P_NULL);\
    MPI_DATA(frw,wp,sp,simpar);\
    CORE(frw,wp,sp,simpar);\
	AUX1(frw,wp,sp,simpar);\
	ZMIN(simpar) = zmin;\
	ZMAX(simpar) = zmax;\
	ncnt += frw(wp,P_Zmin,sp ZMIN(simpar));\
	ncnt += frw(wp,P_Zmax,sp ZMAX(simpar));\
}while(0)
