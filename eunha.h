#include "astrounits.h"

#ifdef XYZDBL

typedef double PosType;
#define MPI_POSTYPE MPI_DOUBLE

#else

typedef float PosType;
#define MPI_POSTYPE MPI_FLOAT

#endif

#ifndef DENTYPE
typedef float DenType;
#define DENTYPE
#endif

#define MPI_PTRDIFF_T MPI_LONG

#define MPI_DEN_T MPI_FLOAT
#define MAX_ZOOM_LEVEL 5

#define WHERESTR  "[file %s, line %d]: "
#define WHEREARG  __FILE__, __LINE__
#define DEBUGPRINT2(...)       fprintf(stderr, __VA_ARGS__);
#define DEBUGPRINT(_fmt, ...)  do{if(DEBUG) DEBUGPRINT2(WHERESTR _fmt, WHEREARG, __VA_ARGS__)}while(0)
#define DEBUGPRINT0(_fmt) do{if(DEBUG) DEBUGPRINT2(WHERESTR _fmt, WHEREARG)} while(0)


#ifndef EUNHA_H
#define EUNHA_H

#ifndef FFTW3_H
#include "fftw3.h"
#include "fftw3-mpi.h"
#define FFTW3_H
#endif 

#ifndef MPIRMS_H
#include "mpirks.h"
#define MPIRMS_H
#endif 


enum pmtreetype {PMTYPE=0, TREETYPE=1};
#define ntype 5
enum mtype {totm=0,cdm=1,sph=2, star=3, voro=4};
enum boolean {YES=01, NO=00};
enum InitCondition {ZELD=1, TWOLPT=2};
enum dimension {X=1, Y=2, Z=3, VX = 1, VY = 2, VZ = 3};



#include <sys/types.h>
enum SimulationModels {Cosmos=0, Static=1,ZoomedCosmos=2,KH=3,Blast=4,BowShock=5, RT=6,Kepler=7,MkGlass2D=30};
/*
#define Cosmos 0
#define Static 1
#define ZoomedCosmos 2
#define KH 3
#define Blast 4
#define BowShock 5
#define RT 6
*/

/*
enum PMSTATUS {PUSH=0, PULL=1};
enum HalfStep { halfpush=0, halfpull=1, kick=2};
*/
#define PUSH 0
#define PULL 1
#define FULL 2
#define HALFPUSH 0
#define HALFPULL 1
#define KICK 2


typedef struct FFTW_INFO{ 
	fftwf_plan p, ip;
	fftwf_plan np, inp;
	ptrdiff_t local_grid_size, local_z_start, local_nz;
	PosType zstart, zfinal;
	ptrdiff_t local_grid_size_after_transpose,local_y_start_after_transpose, local_ny_after_transpose;
	/*
	MPI_Comm com;
	*/
}FFTW_INFO;

#define LOCAL_GRID_SIZE(simpar) ((simpar)->fftwgrid.fftw_info.local_grid_size)
#define LOCAL_GRID_SIZE_AFTER_TRANSPOSE(simpar) ((simpar)->fftwgrid.fftw_info.local_grid_size_after_transpose)
#define LOCAL_Y_START_AFTER_TRANSPOSE(simpar) ((simpar)->fftwgrid.fftw_info.local_y_start_after_transpose)
#define LOCAL_NY_AFTER_TRANSPOSE(simpar) ((simpar)->fftwgrid.fftw_info.local_ny_after_transpose)
#define FFTWINFO(simpar) ((simpar)->fftwgrid.fftw_info)
#define ZSTART(simpar) ((simpar)->fftwgrid.fftw_info.zstart)
#define ZFINAL(simpar) ((simpar)->fftwgrid.fftw_info.zfinal)
#define LOCAL_NZ(simpar) ((simpar)->fftwgrid.fftw_info.local_nz)
#define LOCAL_Z_START(simpar) ((simpar)->fftwgrid.fftw_info.local_z_start)
#define FFTW_F_PLAN(simpar) ((simpar)->fftwgrid.fftw_info.p)
#define FFTW_B_PLAN(simpar) ((simpar)->fftwgrid.fftw_info.ip)
#define FFTW_NF_PLAN(simpar) ((simpar)->fftwgrid.fftw_info.np)
#define FFTW_NB_PLAN(simpar) ((simpar)->fftwgrid.fftw_info.inp)
/*
#define FFTW_COMM(simpar) ((simpar)->fftwgrid.fftw_info.com)
*/
#define FFTW_COMM(simpar) MPI_COMM(simpar)

typedef struct FFTWGridInfo{
	/*
	ptrdiff_t nx,ny,nz,nspace;
	*/
	FFTW_INFO fftw_info;
	GridInfo gridinfo;
	union{
		DenType *rden;
		fftwf_complex *cden;
	} den;
} FFTWGridInfo;

#define FFTW_RDEN(simpar) ( (simpar)->fftwgrid.den.rden)
#define FFTW_CDEN(simpar) ( (simpar)->fftwgrid.den.cden)


typedef struct Basic_MPI{
	int myid,nid;
	MPI_Comm com;
}Basic_MPI;

typedef struct HydroExamBox{
	PosType xmin,ymin,zmin,xmax,ymax,zmax;
} HydroExamBox;

typedef struct SimBox{
	float boxsize;
	SimBoxRange simbox;
	HydroExamBox hydroExamBox;
	int nddinfo;
	DoDeInfo *dm_ddinfo, *sph_ddinfo, *voro_ddinfo, *star_ddinfo, *agn_ddinfo;
	DoDeInfo *tdm_ddinfo, *tsph_ddinfo, *tvoro_ddinfo, *tstar_ddinfo, *tagn_ddinfo;
	DoDeFunc dm_ddfunc, sph_ddfunc, voro_ddfunc, star_ddfunc, agn_ddfunc;
	DoDeFunc tdm_ddfunc, tsph_ddfunc, tvoro_ddfunc, tstar_ddfunc, tagn_ddfunc;

	DoDeInfo *vorork4_ddinfo;
	DoDeInfo *tvorork4_ddinfo;
	DoDeFunc tvorork4_ddfunc;
	DoDeFunc vorork4_ddfunc;
} Simbox;

#define Xmin_HydroExam(simpar) ((simpar)->simmodel.simbox.hydroExamBox.xmin)
#define Ymin_HydroExam(simpar) ((simpar)->simmodel.simbox.hydroExamBox.ymin)
#define Zmin_HydroExam(simpar) ((simpar)->simmodel.simbox.hydroExamBox.zmin)
#define Xmax_HydroExam(simpar) ((simpar)->simmodel.simbox.hydroExamBox.xmax)
#define Ymax_HydroExam(simpar) ((simpar)->simmodel.simbox.hydroExamBox.ymax)
#define Zmax_HydroExam(simpar) ((simpar)->simmodel.simbox.hydroExamBox.zmax)

/*
#define COS_XMIN(simpar) ( (simpar)->simmodel.simbox.simbox.x.min)
#define COS_XMAX(simpar) ( (simpar)->simmodel.simbox.simbox.x.max)
#define COS_YMIN(simpar) ( (simpar)->simmodel.simbox.simbox.y.min)
#define COS_YMAX(simpar) ( (simpar)->simmodel.simbox.simbox.y.max)
#define COS_ZMIN(simpar) ( (simpar)->simmodel.simbox.simbox.z.min)
#define COS_ZMAX(simpar) ( (simpar)->simmodel.simbox.simbox.z.max)
*/

/*
#define SIM_XMIN(simpar) ( (simpar)->simmodel.simbox.simbox.x.min)
#define SIM_XMAX(simpar) ( (simpar)->simmodel.simbox.simbox.x.max)
#define SIM_YMIN(simpar) ( (simpar)->simmodel.simbox.simbox.y.min)
#define SIM_YMAX(simpar) ( (simpar)->simmodel.simbox.simbox.y.max)
#define SIM_ZMIN(simpar) ( (simpar)->simmodel.simbox.simbox.z.min)
#define SIM_ZMAX(simpar) ( (simpar)->simmodel.simbox.simbox.z.max)
*/

#define SIM_LXMIN(simpar,type) \
	( (simpar)->simmodel.simbox.type##_ddinfo[NDDINFO(simpar)-1].lgroup.r.rmin[0])
#define SIM_LXMAX(simpar,type) \
	( (simpar)->simmodel.simbox.type##_ddinfo[NDDINFO(simpar)-1].lgroup.r.rmax[0])
#define SIM_LYMIN(simpar,type) \
	( (simpar)->simmodel.simbox.type##_ddinfo[NDDINFO(simpar)-1].lgroup.r.rmin[1])
#define SIM_LYMAX(simpar,type) \
	( (simpar)->simmodel.simbox.type##_ddinfo[NDDINFO(simpar)-1].lgroup.r.rmax[1])
#define SIM_LZMIN(simpar,type) \
	( (simpar)->simmodel.simbox.type##_ddinfo[NDDINFO(simpar)-1].lgroup.r.rmin[2])
#define SIM_LZMAX(simpar,type) \
	( (simpar)->simmodel.simbox.type##_ddinfo[NDDINFO(simpar)-1].lgroup.r.rmax[2])


typedef struct IndTInfo {
	char indTflag;
	char iflagfixedlist; 
	int nowTsubdiv, maxTsubpower; 
	int nsubstep,tnumcount; 
	int nsubgas,nsubfixed,nsubnbody;
	int nsubstepcount,isubstep;
	double anext;
	double abefore,astep, da,damin; 
}IndTInfo;


typedef struct EvolInfo{
	float pfact, fact1, fact2;
} EvolInfo;

typedef struct TimeInfo {
	int nstep,nskip,stepcount,stepnum;
	float ai,amax;
	float zinit, redshift;
	double anow, astep;
	double lookbacktime;
	int pmstatus;
	/*
	int iflagPreFoF,iflagsyncpdata; 
	*/
	IndTInfo indt;
	EvolInfo evolinfo;
} TimeInfo;
#define LOOKBACKTIME(simpar) ( (simpar)->timeinfo.lookbacktime)

typedef struct IC{
	int porder;
	int flagIC;
	float growthfactor;
	float pamp[ntype], damp1[ntype], damp2[ntype], vamp1[ntype], vamp2[ntype];
	long iseed;
	char graficDir[64];
} IC;

typedef struct DE{
	float wlam0, wlama, Pwcorection;
//	char constw;
} DE;
typedef struct Cosmology{
	float omei,omep, omepb,omeplam,omepk,cosconx; 
	float fNL, gNL;
	float hubble, npow;
	float bias8;
	char powfilename[64], inpapkfilename[64];
	char GlacialHeader[64];
	int powreadflag;
	DE de;
	IC ic;
} Cosmology;


typedef struct StaticBG{ 
	char externalforce;
	float TotMass;
	double ffTime, ffTimeellb;
} StaticBG;



typedef struct KHINS{
	float rho1,rho2, vel1, vel2, deltay,pressure, vperturb;
	float xmin,ymin,xmax,ymax;
	float OrderofAccuracy,Kappa;
	int Nxp, Nyp;
}KHINS;

#define VoroAccuracyOrder(simpar) ( (simpar)->physics.gasinfo.AccuracyOrder)
#define KH_OA(simpar) VoroAccuracyOrder(simpar)

typedef struct RTI{
	float rho1,rho2, gravity, deltay, vperturb, Phalf;
	float xmin,ymin,xmax,ymax;
	float OrderofAccuracy,kappa;
	int Nxp, Nyp;
}RTI;
#define RT_OA(simpar) VoroAccuracyOrder(simpar)

typedef struct KeplyerI{
    float eps;
    float xmin,ymin,xmax,ymax;
    float OrderofAccuracy,kappa;
    int Nxp, Nyp;
}KPI;
#define KP_OA(simpar) VoroAccuracyOrder(simpar)
#define KP_EPS(simpar) ((simpar)->simmodel.kp.eps)

typedef struct Glass2D{
    float eps;
    float xmin,ymin,xmax,ymax;
    float OrderofAccuracy, Kappa;
    int Nxp, Nyp;
}Glass2D;
#define GL2D_OA(simpar) VoroAccuracyOrder(simpar)
#define GL2D_EPS(simpar) ((simpar)->simmodel.gl2d.eps)








/* gastype: sph=S, voro=V, nogas= N */
typedef struct gas_flag{
	char gastype, flagSF, flagCOOL, flagSNFB, flagBGHEAT, flagFixedneigh;
	char constmu;
}gas_flag;

typedef struct GasInfo{
	float mu;
	float visforcefactor; 
	float alphavis, betavis; // artificial viscosity factor of Voronoi
	float etavis, epsvis; // artificial viscosity factor of Voronoi 
	float viscosity; // viscosity in the Voronoi */
	float Tempevolfactor,sphforcefactor; 
	float Yp,alpha,gamma,rhos2rhor; 
	float initMetal, coolfact, heatfact;
	float entropyfact,entropyevolfact;
	float g12ratio,g1expansion,Courant,minTemp,minEntropyfactor,Durant; 
	float Trad0,SFvirialDen,SFrho577,SFgasden,SFtemp,SFrealgasdenfactor; 
	float HmassFrac,UVShieldDen,dUVShieldDen; 
	float Mstar;/* fractional mass of sph particle which goes into star */ 
	float xoffset,yoffset,zoffset; // offset of gas particle w.r.t. dm particle
	float Cstar;
	float simden2Hnden,simden2HpHenden;
	float AccuracyOrder;
	double g1,g2, rhoc0,meanrho; 
	struct gas_flag flag;
}GasInfo;
#define GAS_MEANRHO(simpar) ((simpar)->physics.gasinfo.meanrho)
#define GAS_RHOS2RHOR(simpar) ((simpar)->physics.gasinfo.rhos2rhor)
#define GAS_RHOC0(simpar) ((simpar)->physics.gasinfo.rhoc0)
#define GAS_RHOC(simpar) ((simpar)->physics.gasinfo.rhoc)
#define GAS_RHO(simpar) ((simpar)->physics.gasinfo.rho)
#define GAS_VISFORCEFACTOR(simpar) ((simpar)->physics.gasinfo.visforcefactor)
#define GAS_AlphaVis(simpar) ((simpar)->physics.gasinfo.alphavis)
#define GAS_BetaVis(simpar) ((simpar)->physics.gasinfo.betavis)
#define GAS_ETAVIS(simpar) ((simpar)->physics.gasinfo.etavis)
#define GAS_EPSVIS(simpar) ((simpar)->physics.gasinfo.epsvis)
#define GAS_VISCOSITY(simpar) ((simpar)->physics.gasinfo.viscosity)
#define GAS_SPHFORCEFACTOR(simpar) ((simpar)->physics.gasinfo.sphforcefactor)
#define GAS_TEMPEVOLFACTOR(simpar) ((simpar)->physics.gasinfo.Tempevolfactor)
#define GAS_G1(simpar) ((simpar)->physics.gasinfo.g1)
#define GAS_G2(simpar) ((simpar)->physics.gasinfo.g2)
#define GAS_G12RATIO(simpar) ((simpar)->physics.gasinfo.g12ratio)
#define GAS_G1EXPANSION(simpar) ((simpar)->physics.gasinfo.g1expansion)
#define GAS_ENTROPYFACT(simpar) ((simpar)->physics.gasinfo.entropyfact)
#define GAS_ENTROPYEVOLFACT(simpar) ((simpar)->physics.gasinfo.entropyevolfact)
#define GAS_SFRHO577(simpar) ((simpar)->physics.gasinfo.SFrho577)
#define GAS_SFVIRIALDEN(simpar) ((simpar)->physics.gasinfo.SFvirialDen)
#define GAS_COOLFACT(simpar) ((simpar)->physics.gasinfo.coolfact)
#define GAS_HEATFACT(simpar) ((simpar)->physics.gasinfo.heatfact)
#define GAS_MINENTROPYFACTOR(simpar) ((simpar)->physics.gasinfo.minEntropyfactor)
#define GAS_SIMDEN2HNDEN(simpar) ((simpar)->physics.gasinfo.simden2Hnden)
#define GAS_SIMDEN2HPHENDEN(simpar) ((simpar)->physics.gasinfo.simden2HpHenden)
typedef struct Gravity{
	char gravflag;
	float epsilon, efold;
	float theta, rsphere;
}Gravity;


typedef struct BlastInfo{
	float InitTemp, BlastExplodingTemp;
	int mu;
}BlastInfo;

typedef struct BowShockInfo{
	float InitTemp, highden, lowden, velshock;
	float mu;
}BowShockInfo;
typedef struct RTInfo{
	float InitTemp, highmass, lowmass, force;
	float mu;
}RTInfo;


typedef union indxflag{
	size_t indx;
	char Flag[8];
} indxflag;

typedef struct linkedlisttype{
	struct linkedlisttype *next;
	indxflag u4if;
#ifndef GOTPM
	char indt[4];
#endif
	float x,y,z;
	float mass;
} linkedlisttype;


#define LINKEDTYPE_MASS(simpar, tmp) ( (tmp)->mass )

typedef struct dmparticletype{
	indxflag u4if;
#ifndef GOTPM
	char indt[4];
#endif
	float x,y,z;
	float vx,vy,vz;
#ifndef GOTPM
	float ax,ay,az;
#endif
}dmparticletype;
typedef struct treedmparticletype{
	struct linkedlisttype *next;
	indxflag u4if;
#ifndef GOTPM
	char indt[4];
#endif
	float x,y,z;
	float vx,vy,vz;
#ifndef GOTPM
	float ax,ay,az;
#endif
}treedmparticletype;

/* NOTE: sph particle could be used for other hydro particle type */

#define SphQ float Entropy,rho,fi,F1,hsml,metallicity,dAs,divVel,mu,delta_e,nden,temp; int nstar
#define BPSphQ SphQ;int pid;float ometallicity,odAs,omass; treesphparticletype *bp
#define StarQ float sftime,metallicity, swlasttime; int SNmark


typedef struct sphparticletype{
	indxflag u4if;
	char indt[4];
	float x,y,z;
	float mass;
	float vx,vy,vz;
	float ax,ay,az;
	SphQ;
	struct treesphparticletype *bp;
}sphparticletype;

typedef struct treesphparticletype{
	struct linkedlisttype *next;
	indxflag u4if;
	char indt[4];
	float x,y,z;
	float mass;
	float vx,vy,vz;
	float ax,ay,az;
	SphQ;
	struct treesphparticletype *bp;
}treesphparticletype;

#define VoroQ float \
	den,\
	pressure,\
	te,\
	ie,\
	ke, \
	die,\
	dke,\
	dte, \
	csound,\
	dt, \
	volume,\
	w2, w2ceil // Laguerre circular radius

typedef struct voroparticletype{
    indxflag u4if;
    char indt[4];
    float x,y,z;
    float mass;
    float vx,vy,vz;
    float ax,ay,az;
    VoroQ;
    struct treevoroparticletype *bp;
}voroparticletype;

typedef struct treevoroparticletype{
    struct linkedlisttype *next;
    indxflag u4if;
    char indt[4];
    float x,y,z;
    float mass;
    float vx,vy,vz;
    float ax,ay,az;
    VoroQ;
    struct treevoroparticletype *bp;
}treevoroparticletype;

// RK4 coefficients
typedef struct RK4{
	float k1x,k1y,k1vx,k1vy;
	float k2x,k2y,k2vx,k2vy;
	float k3x,k3y,k3vx,k3vy;
	float k4x,k4y,k4vx,k4vy;
	float w2backup;
} RK4;


typedef struct vorork4particletype{
    indxflag u4if;
    char indt[4];
    float x,y,z;
    float mass;
    float vx,vy,vz;
    float ax,ay,az;
    VoroQ;
    RK4 rk4;
    struct treevorork4particletype *bp;
}vorork4particletype;
typedef struct treevorork4particletype{
    struct linkedlisttype *next;
    indxflag u4if;
    char indt[4];
    float x,y,z;
    float mass;
    float vx,vy,vz;
    float ax,ay,az;
    VoroQ;
    RK4 rk4;
    struct treevorork4particletype *bp;
}treevorork4particletype;



typedef struct starparticletype{
	indxflag u4if;
	char indt[4];
	float x,y,z;
	float mass;
	float vx,vy,vz;
	float ax,ay,az;
	StarQ;
}starparticletype;
typedef struct treestarparticletype{
	struct linkedlisttype *next;
	indxflag u4if;
	char indt[4];
	float x,y,z;
	float mass;
	float vx,vy,vz;
	float ax,ay,az;
	StarQ;
}treestarparticletype;

typedef struct agnparticletype{
	indxflag u4if;
	char indt[4];
	float x,y,z;
	float mass;
	float vx,vy,vz;
	float ax,ay,az;
}agnparticletype;
typedef struct treeagnparticletype{
	struct linkedlisttype *next;
	indxflag u4if;
	char indt[4];
	float x,y,z;
	float mass;
	float vx,vy,vz;
	float ax,ay,az;
}treeagnparticletype;

typedef struct DM{
	float mass;
	ptrdiff_t np,tnp;
	ptrdiff_t npad;
	union{ 
		dmparticletype *bp;
		treedmparticletype *tbp;
	} p;
	union{
		dmparticletype *bp;
		treedmparticletype *tbp;
	}padding;
} DM;
typedef struct SPH{
	float initmass;
	float massinsolarmass, interactionsphere;
	int NumNear,timesteplimiter, constneighborpower,init_flag;
	ptrdiff_t np,tnp;
	ptrdiff_t npad;
	union{ 
		sphparticletype *bp;
		treesphparticletype *tbp;
	} p;
	union{ 
		sphparticletype *bp;
		treesphparticletype *tbp;
	} padding;
} SPH;

typedef struct VORO{
    float initmass;
    float massinsolarmass, interactionsphere;
    int NumNear,timesteplimiter, constneighborpower,init_flag;
    ptrdiff_t np,tnp;
    ptrdiff_t npad;
    union{
        voroparticletype *bp;
        vorork4particletype *bprk4;
        treevoroparticletype *tbp;
        treevorork4particletype *tbprk4;
    } p;
    union{
        voroparticletype *bp;
        vorork4particletype *bprk4;
        treevoroparticletype *tbp;
        treevorork4particletype *tbprk4;
    } padding;
} VORO;



typedef struct STAR{
	ptrdiff_t np,tnp;
	ptrdiff_t npad;
	union{ 
		starparticletype *bp;
		treestarparticletype *tbp;
	}p;
	union{ 
		starparticletype *bp;
		treestarparticletype *tbp;
	} padding;
} STAR;
typedef struct AGN{
	ptrdiff_t np,tnp;
	ptrdiff_t npad;
	union{ 
		agnparticletype *bp;
		treeagnparticletype *tbp;
	}p;
	union{ 
		agnparticletype *bp;
		treeagnparticletype *tbp;
	}padding;
} AGN;


typedef struct Particles{
	int xyzshiftflag;
	enum pmtreetype ptype;
	size_t npsum, ptypesize;
	DM dm;
	SPH sph;
	VORO voro;
	STAR star;
	AGN agn;
} Particles;

typedef struct SimModels{
	char BGExpand;
	char flagHydro;
	int SimModel;
	Cosmology cosmos;
	StaticBG staticworld;
	KHINS kh;
	RTI rt;
	BlastInfo blast;
	BowShockInfo bowshock;
	RTInfo rtinfo;
	KPI kp;
	Glass2D gl2d;
	Simbox simbox;
}SimModels;

typedef struct FILE_IO{
	char rvfilename[100], rvprefix[100];
	int WGroupSize, NWGroup;
} FILE_IO;
typedef struct HALF{
	char first, second;
} PMHALF;
typedef struct Step_ControlFlag{
	char PMPreFoFflag;
	char savexzslice;
	char flagpsmeasure, flagPreFoF, flagsyncpdata,flagwholeden;
	char flagcontinue;
	int halfstep;
	PMHALF pmhalf;
} Step_ControlFlag;



#define CONT_FLAGPREFOF(simpar) ((simpar)->control.flagPreFoF)
#define CONT_FLAGPS(simpar) ((simpar)->control.flagpsmeasure)
#define CONT_FLAGSYNCP(simpar) ((simpar)->control.flagsyncpdata)
#define CONT_FLAGWHOLEDEN(simpar) ((simpar)->control.flagwholeden)
/* E: Emergency stop
 * D: Dump & Continue
 * P: Pass
 * S: Stop after Dump
 * */
#define CONT_FLAGCONTINUE(simpar) ((simpar)->control.flagcontinue)
#define CONT_HALFSTEP(simpar) ((simpar)->control.halfstep)

#define PM_FIRST_HALF(simpar) ((simpar)->control.pmhalf.first)
#define PM_SECOND_HALF(simpar) ((simpar)->control.pmhalf.second)

typedef struct Animation{
	char flag;
	int nviewer;
	char viewfilename[100];
}Animation;
typedef struct ObserverMode{
	int flag;
	char observer[100];
}ObserverMode;


typedef struct Physics{
	GasInfo gasinfo;
	Gravity gravinfo;
}Physics;

typedef struct Virialization{
	PosType fof_link;
	size_t minNum;
} Virialization;
#define FOF_LINK(simpar) ( (simpar)->virial.fof_link )
#define FOF_MINNUM(simpar) ( (simpar)->virial.minNum )

typedef struct HydroTreeLinkedCell{ 
	linkedlisttype *link; 
	size_t nmem; 
	int calflag;
} HydroTreeLinkedCell;
typedef struct TreeLinkedCell{ 
	linkedlisttype *link; 
	size_t nmem; 
} TreeLinkedCell;


typedef struct LinkedCell{
	PosType CellWidth, hydroGridSize;
	PosType invCellWidth;
	size_t mx,my,mz;
	TreeLinkedCell *BasicCell;
	HydroTreeLinkedCell *SPH_BasicCell;
	HydroTreeLinkedCell *VORO_BasicCell;
	HydroTreeLinkedCell *VORORK4_BasicCell;
	HydroTreeLinkedCell *STAR_BasicCell;
	HydroTreeLinkedCell *AGN_BasicCell;
} LinkedCell;


typedef struct Analysis{
	PosType zmin,zmax;
} Analysis;

typedef struct SimParameters{
	int izoom,nzoom;
	struct Basic_MPI mpi_info;
	struct FFTWGridInfo fftwgrid;
	struct TimeInfo timeinfo;
	struct SimModels simmodel;
	struct LinkedCell linkedcell;

	struct Physics physics;
	struct Particles bp;

	struct FILE_IO fileio;
	struct Step_ControlFlag control;
	struct Animation anim;
	struct ObserverMode obs;
	struct Virialization virial;
	struct Analysis anal;
} SimParameters;

#endif

// the domain boundary in z direction for a given simpar
#define ZMIN(simpar) ( (simpar)->anal.zmin) 
#define ZMAX(simpar) ( (simpar)->anal.zmax) 



// extract the information on the basic cell for the tree and sph force calculations
#define BASICCELL_CELLWIDTH(simpar) ((simpar)->linkedcell.CellWidth)
#define BASICCELL_INVCELLWIDTH(simpar) ( (simpar)->linkedcell.invCellWidth )
#define HydroGridSize(simpar) ((simpar)->linkedcell.hydroGridSize)

#define BASICCELL_MX(simpar) ( (simpar)->linkedcell.mx )
#define BASICCELL_MY(simpar) ( (simpar)->linkedcell.my )
#define BASICCELL_MZ(simpar) ( (simpar)->linkedcell.mz )

#define BASICCELL(simpar) ( (simpar)->linkedcell.BasicCell)
#define SPH_BASICCELL(simpar) ( (simpar)->linkedcell.SPH_BasicCell)
#define VORO_BASICCELL(simpar) ( (simpar)->linkedcell.VORO_BasicCell)
#define VORORK4_BASICCELL(simpar) ( (simpar)->linkedcell.VORORK4_BasicCell)
#define STAR_BASICCELL(simpar) ( (simpar)->linkedcell.STAR_BasicCell)
#define AGN_BASICCELL(simpar) ( (simpar)->linkedcell.AGN_BasicCell)

// returning the pointer to the mpi_info
#define PTR2MPIINFO(simpar) ( &((simpar)->mpi_info))

// a flag to do the gravity calculation
#define GRAVITY(simpar) ((simpar)->physics.gravinfo.gravflag)

#define GAS_TYPE(simpar) ((simpar)->physics.gasinfo.flag.gastype)
// mean molecular weight of gas
#define GAS_MU(simpar) ((simpar)->physics.gasinfo.mu)


// a flag to use the constant mean molecular weight
#define GAS_CONSTMU(simpar) ((simpar)->physics.gasinfo.flag.constmu)

// a flag to use the Star formation algorithm 
#define GAS_SFFLAG(simpar) ((simpar)->physics.gasinfo.flag.flagSF)

// a flag to cool the gas
#define GAS_COOLFLAG(simpar) ((simpar)->physics.gasinfo.flag.flagCOOL)

// a flag to use the SN feedbacks
#define GAS_SNFBFLAG(simpar) ((simpar)->physics.gasinfo.flag.flagSNFB)
// a flag to use the fixed number of neighbors in the SPH
#define GAS_FNFLAG(simpar) ((simpar)->physics.gasinfo.flag.flagFixedneigh)
// a flag to use the background heating
#define GAS_BGHEATFLAG(simpar) ((simpar)->physics.gasinfo.flag.flagBGHEAT)



#define GAS_DURANT(simpar) ((simpar)->physics.gasinfo.Durant)
// input Yp value
#define GAS_YP(simpar) ((simpar)->physics.gasinfo.Yp)
#define GAS_TRAD0(simpar) ((simpar)->physics.gasinfo.Trad0)
// input Gamma value of gas
#define GAS_GAMMA(simpar) ((simpar)->physics.gasinfo.gamma)
#define GAS_ALPHA(simpar) ((simpar)->physics.gasinfo.alpha)
// input Courant number
#define GAS_COURANT(simpar) ((simpar)->physics.gasinfo.Courant)
// UV shielding density threshold
#define GAS_UVSHIELDDEN(simpar) ((simpar)->physics.gasinfo.UVShieldDen)
#define GAS_DUVSHIELDDEN(simpar) ((simpar)->physics.gasinfo.dUVShieldDen)

#define GAS_INITMETAL(simpar) ((simpar)->physics.gasinfo.initMetal)
#define GAS_MINTEMP(simpar) ((simpar)->physics.gasinfo.minTemp)
#define GAS_SFVIRDEN(simpar) ((simpar)->physics.gasinfo.SFvirialDen)
#define GAS_SFGASDEN(simpar) ((simpar)->physics.gasinfo.SFgasden)
#define GAS_SFTEMP(simpar) ((simpar)->physics.gasinfo.SFtemp)
#define GAS_MSTAR(simpar) ((simpar)->physics.gasinfo.Mstar)
#define GAS_CSTAR(simpar) ((simpar)->physics.gasinfo.Cstar)


#define SPH_NUMNEAR(simpar) ((simpar)->bp.sph.NumNear)
#define SPH_INITFLAG(simpar) ((simpar)->bp.sph.init_flag)
#define SPH_MASSINSOLARMASS(simpar) ((simpar)->bp.sph.massinsolarmass)
#define SPH_INIT_MASS(simpar) ((simpar)->bp.sph.initmass)
#define SPH_INTERACTIONSPHERE(simpar) ((simpar)->bp.sph.interactionsphere)
#define SPH_INITMASS(simpar) ((simpar)->bp.sph.initmass)
#define SPH_CONSTNEIGHPOW(simpar) ((simpar)->bp.sph.constneighborpower)
#define SPH_TIMESTEPLIMITER(simpar) ((simpar)->bp.sph.timesteplimiter)

#define VORO_NUMNEAR(simpar) ((simpar)->bp.voro.NumNear)
#define VORO_INITFLAG(simpar) ((simpar)->bp.voro.init_flag)
#define VORO_MASSINSOLARMASS(simpar) ((simpar)->bp.voro.massinsolarmass)
#define VORO_INIT_MASS(simpar) ((simpar)->bp.voro.initmass)
#define VORO_INTERACTIONSPHERE(simpar) ((simpar)->bp.voro.interactionsphere)
#define VORO_INITMASS(simpar) ((simpar)->bp.voro.initmass)
#define VORO_CONSTNEIGHPOW(simpar) ((simpar)->bp.voro.constneighborpower)
#define VORO_TIMESTEPLIMITER(simpar) ((simpar)->bp.voro.timesteplimiter)



#define DM_NPP(simpar) ((simpar)->bp.dm.npad)
#define STAR_NPP(simpar) ((simpar)->bp.star.npad)
#define SPH_NPP(simpar) ((simpar)->bp.sph.npad)
#define VORO_NPP(simpar) ((simpar)->bp.voro.npad)
#define AGN_NPP(simpar) ((simpar)->bp.agn.npad)

#define DM_BPP(simpar) ((simpar)->bp.dm.padding.bp)
#define SPH_BPP(simpar) ((simpar)->bp.sph.padding.bp)
#define VORO_BPP(simpar) ((simpar)->bp.voro.padding.bp)
#define VORORK4_BPP(simpar) ((simpar)->bp.voro.padding.bprk4)
#define STAR_BPP(simpar) ((simpar)->bp.star.padding.bp)
#define AGN_BPP(simpar) ((simpar)->bp.agn.padding.bp)

#define DM_TBPP(simpar) ((simpar)->bp.dm.padding.tbp)
#define SPH_TBPP(simpar) ((simpar)->bp.sph.padding.tbp)
#define VORO_TBPP(simpar) ((simpar)->bp.voro.padding.tbp)
#define VORORK4_TBPP(simpar) ((simpar)->bp.voro.padding.tbprk4)
#define STAR_TBPP(simpar) ((simpar)->bp.star.padding.tbp)
#define AGN_TBPP(simpar) ((simpar)->bp.agn.padding.tbp)

#define DM_TNP(simpar) ((simpar)->bp.dm.tnp)
#define SPH_TNP(simpar) ((simpar)->bp.sph.tnp)
#define VORO_TNP(simpar) ((simpar)->bp.voro.tnp)
#define STAR_TNP(simpar) ((simpar)->bp.star.tnp)
#define AGN_TNP(simpar) ((simpar)->bp.agn.tnp)



#define BLAST_INITTEMP(simpar) ((simpar)->simmodel.blast.InitTemp)
#define BLAST_EXPLODINGTEMP(simpar) ((simpar)->simmodel.blast.BlastExplodingTemp)
#define BLAST_MU(simpar) ((simpar)->simmodel.blast.mu)

#define KH_Rho1(simpar) ((simpar)->simmodel.kh.rho1)
#define KH_Rho2(simpar) ((simpar)->simmodel.kh.rho2)
#define KH_Vel1(simpar) ((simpar)->simmodel.kh.vel1)
#define KH_Vel2(simpar) ((simpar)->simmodel.kh.vel2)
#define KH_Deltay(simpar) ((simpar)->simmodel.kh.deltay)
#define KH_Vperturb(simpar) ((simpar)->simmodel.kh.vperturb)
#define KH_Pressure(simpar) ((simpar)->simmodel.kh.pressure)
/* local domain range */
#define KH_XMAX(simpar) Xmax_HydroExam(simpar)
#define KH_YMAX(simpar) Ymax_HydroExam(simpar)
#define KH_XMIN(simpar) Xmin_HydroExam(simpar)
#define KH_YMIN(simpar) Ymin_HydroExam(simpar)
#define KH_GridSize(simpar) HydroGridSize(simpar)
#define KH_Kappa(simpar) ((simpar)->simmodel.kh.Kappa)
/* RT instability */
#define RT_XMAX(simpar) Xmax_HydroExam(simpar)
#define RT_YMAX(simpar) Ymax_HydroExam(simpar)
#define RT_XMIN(simpar) Xmin_HydroExam(simpar)
#define RT_YMIN(simpar) Ymin_HydroExam(simpar)
#define RT_DEN1(simpar) ((simpar)->simmodel.rt.rho1)
#define RT_DEN2(simpar) ((simpar)->simmodel.rt.rho2)
#define RT_ACC(simpar) ((simpar)->simmodel.rt.gravity)
#define RT_GridSize(simpar) HydroGridSize(simpar)
#define RT_Deltay(simpar) ((simpar)->simmodel.rt.deltay)
#define RT_Vperturb(simpar) ((simpar)->simmodel.rt.vperturb)
#define RT_Phalf(simpar) ((simpar)->simmodel.rt.Phalf)
#define RT_Kappa(simpar) ((simpar)->simmodel.rt.Kappa)
/* Kepler test */
#define KP_XMAX(simpar) Xmax_HydroExam(simpar)
#define KP_YMAX(simpar) Ymax_HydroExam(simpar)
#define KP_XMIN(simpar) Xmin_HydroExam(simpar)
#define KP_YMIN(simpar) Ymin_HydroExam(simpar)
#define KP_ACC(simpar) ((simpar)->simmodel.kp.gravity)
#define KP_GridSize(simpar) HydroGridSize(simpar)
#define KP_Phalf(simpar) ((simpar)->simmodel.kp.Phalf)
#define KP_Kappa(simpar) ((simpar)->simmodel.kp.Kappa)
/* 2D Glass Mading */
#define GL2D_XMAX(simpar) Xmax_HydroExam(simpar)
#define GL2D_YMAX(simpar) Ymax_HydroExam(simpar)
#define GL2D_XMIN(simpar) Xmin_HydroExam(simpar)
#define GL2D_YMIN(simpar) Ymin_HydroExam(simpar)
#define GL2D_ACC(simpar) ((simpar)->simmodel.gl2d.gravity)
#define GL2D_GridSize(simpar) HydroGridSize(simpar)
#define GL2D_Phalf(simpar) ((simpar)->simmodel.gl2d.Phalf)
#define GL2D_Kappa(simpar) ((simpar)->simmodel.gl2d.Kappa)



#define BS_INITTEMP(simpar) ((simpar)->simmodel.bowshock.InitTemp)
#define BS_HIGHDEN(simpar) ((simpar)->simmodel.bowshock.highden)
#define BS_LOWDEN(simpar) ((simpar)->simmodel.bowshock.lowden)
#define BS_VELSHOCK(simpar) ((simpar)->simmodel.bowshock.velshock)
#define BS_MU(simpar) ((simpar)->simmodel.bowshock.mu)
#define RT_MU(simpar) ((simpar)->simmodel.rt.mu)
#define RT_INITTEMP(simpar) ((simpar)->simmodel.rt.InitTemp)
#define RT_HIGHMASS(simpar) ((simpar)->simmodel.rt.highmass)
#define RT_LOWMASS(simpar) ((simpar)->simmodel.rt.lowmass)
#define RT_FORCE(simpar) ((simpar)->simmodel.rt.force)

#define RV_FILE(simpar) ((simpar)->fileio.rvfilename)
#define RV_FILEPREFIX(simpar) ((simpar)->fileio.rvprefix)
#define WGROUPSIZE(simpar) ((simpar)->fileio.WGroupSize)
#define NWGROUP(simpar) ((simpar)->fileio.NWGroup)

#define ANIM_VIEWFILE(simpar) ((simpar)->anim.viewfilename)
#define ANIM_FLAG(simpar) ((simpar)->anim.flag)
#define ANIM_NVIEWER(simpar) ((simpar)->anim.nviewer)

#define OBS_FLAG(simpar) ((simpar)->obs.flag)

#define PMPreFoFFLAG(simpar) ((simpar)->control.PMPreFoFflag)
#define FLAG4SAVEXZSLICE(simpar) ((simpar)->control.savexzslice)
#define GLACIALHEADER(simpar) ((simpar)->simmodel.cosmos.GlacialHeader)


#define NPSUM(simpar) ((simpar)->bp.npsum)

#define DM_NP(simpar) ((simpar)->bp.dm.np)
#define STAR_NP(simpar) ((simpar)->bp.star.np)
#define SPH_NP(simpar) ((simpar)->bp.sph.np)
#define VORO_NP(simpar) ((simpar)->bp.voro.np)
#define AGN_NP(simpar) ((simpar)->bp.agn.np)

#define DM_NPAD(simpar) ((simpar)->bp.dm.npad)
#define SPH_NPAD(simpar) ((simpar)->bp.sph.npad)
#define VORO_NPAD(simpar) ((simpar)->bp.voro.npad)
#define STAR_NPAD(simpar) ((simpar)->bp.star.npad)
#define AGN_NPAD(simpar) ((simpar)->bp.agn.npad)

#define DM_BP(simpar) ((simpar)->bp.dm.p.bp)
#define DM_TBP(simpar) ((simpar)->bp.dm.p.tbp)
#define SPH_BP(simpar) ((simpar)->bp.sph.p.bp)
#define SPH_TBP(simpar) ((simpar)->bp.sph.p.tbp)
#define VORO_BP(simpar) ((simpar)->bp.voro.p.bp)
#define VORORK4_BP(simpar) ((simpar)->bp.voro.p.bprk4)
#define VORO_TBP(simpar) ((simpar)->bp.voro.p.tbp)
#define VORORK4_TBP(simpar) ((simpar)->bp.voro.p.tbprk4)
#define STAR_BP(simpar) ((simpar)->bp.star.p.bp)
#define STAR_TBP(simpar) ((simpar)->bp.star.p.tbp)
#define AGN_BP(simpar) ((simpar)->bp.agn.p.bp)
#define AGN_TBP(simpar) ((simpar)->bp.agn.p.tbp)

/*
#define IFLAGSYNCPDATA(simpar) ((simpar)->timeinfo.iflagsyncpdata)
#define IFLAGPREFOF(simpar) ((simpar)->timeinfo.iflagPreFoF)
#define IFLAGSYNCPDATA(simpar) ((simpar)->control.flagsyncpdata)
#define IFLAGPREFOF(simpar) ((simpar)->control.flagPreFoF)
*/

#define XYZSHIFT(simpar) ((simpar)->bp.xyzshiftflag)

#define PTYPE(simpar) ((simpar)->bp.ptype)
#define PTYPESIZE(simpar) ((simpar)->bp.ptypesize)

#define Evol_PFACT(simpar) ((simpar)->timeinfo.evolinfo.pfact)
#define Evol_FACT1(simpar) ((simpar)->timeinfo.evolinfo.fact1)
#define Evol_FACT2(simpar) ((simpar)->timeinfo.evolinfo.fact2)

#define IndT_FLAG(simpar) ((simpar)->timeinfo.indt.indTflag)
#define IndT_ABEFORE(simpar) ((simpar)->timeinfo.indt.abefore)
#define IndT_DAMIN(simpar) ((simpar)->timeinfo.indt.damin)
#define IndT_ASTEP(simpar) ((simpar)->timeinfo.indt.astep)
#define IndT_DA(simpar) ((simpar)->timeinfo.indt.da)
#define IndT_NOWTSUBDIV(simpar) ((simpar)->timeinfo.indt.nowTsubdiv)
#define IndT_MAXTSUBPOWER(simpar) ((simpar)->timeinfo.indt.maxTsubpower)

/*
#define IndT_IFLAGPREFOF(simpar) ((simpar)->timeinfo.iflagPreFoF)
#define IndT_IFLAGSYNCPDATA(simpar) ((simpar)->timeinfo.iflagsyncpdata)
*/
#define IndT_IFLAGPREFOF(simpar) ((simpar)->control.flagPreFoF)
#define IndT_IFLAGSYNCPDATA(simpar) ((simpar)->control.flagsyncpdata)

#define IndT_NSUBSTEP(simpar) ((simpar)->timeinfo.indt.nsubstep)
#define IndT_NSUBSTEPCOUNT(simpar) ((simpar)->timeinfo.indt.nsubstepcount)
#define IndT_TNUMCOUNT(simpar) ((simpar)->timeinfo.indt.tnumcount)
#define IndT_ANEXT(simpar) ((simpar)->timeinfo.indt.anext)
#define IndT_NSUBGAS(simpar) ((simpar)->timeinfo.indt.nsubgas)
#define IndT_NSUBNBODY(simpar) ((simpar)->timeinfo.indt.nsubnbody)
#define IndT_NSUBFIXED(simpar) ((simpar)->timeinfo.indt.nsubfixed)
#define IndT_IFLAGFIXEDLIST(simpar) ((simpar)->timeinfo.indt.iflagfixedlist)



#define GRV_EPSILON(simpar) ((simpar)->physics.gravinfo.epsilon)
#define GRV_EFOLD(simpar) ((simpar)->physics.gravinfo.efold)


#define SIMMODEL(simpar) ((simpar)->simmodel.SimModel)
#define flagHydro(simpar) ((simpar)->simmodel.flagHydro)
#define BGEXPAND(simpar) ((simpar)->simmodel.BGExpand)

#define MYID(simpar) 	((simpar)->mpi_info.myid)
#define NID(simpar) 	((simpar)->mpi_info.nid)
#define MPICOM(simpar) 	((simpar)->mpi_info.com)
#define COM(simpar) 	((simpar)->mpi_info.com)
#define MPI_COMM(simpar) ((simpar)->mpi_info.com)

#define NZOOM(simpar) 	((simpar)->nzoom)
#define IZOOM(simpar) 	((simpar)->izoom)

#define COSBOXSIZE(simpar) ((simpar)->simmodel.simbox.boxsize)
#define STATBOXSIZE(simpar) ((simpar)->simmodel.simbox.boxsize)
#define BOXSIZE(simpar) ((simpar)->simmodel.simbox.boxsize)

#define HUBBLE(simpar) 	((simpar)->simmodel.cosmos.hubble)
#define COSCONX(simpar) 	((simpar)->simmodel.cosmos.cosconx)
#define OMEP(simpar) 	((simpar)->simmodel.cosmos.omep)
#define OMEPB(simpar) 	((simpar)->simmodel.cosmos.omepb)
#define OMEPLAM(simpar) 	((simpar)->simmodel.cosmos.omeplam)

#define WLAM0(simpar) 	((simpar)->simmodel.cosmos.de.wlam0)
#define WLAM1(simpar) 	((simpar)->simmodel.cosmos.de.wlama)
#define ConstW(simpar) 	((simpar)->simmodel.cosmos.de.constw)
#define PWCOR(simpar) 	((simpar)->simmodel.cosmos.de.Pwcorrection)
#define OMEI(simpar) 	((simpar)->simmodel.cosmos.omei)
#define FNL(simpar) 	((simpar)->simmodel.cosmos.fNL)
#define GNL(simpar) 	((simpar)->simmodel.cosmos.gNL)
#define NPOW(simpar) 	((simpar)->simmodel.cosmos.npow)
#define BIAS8(simpar) 	((simpar)->simmodel.cosmos.bias8)

#define AI(simpar) 	((simpar)->timeinfo.ai)
#define REDSHIFT(simpar) 	((simpar)->timeinfo.redshift)
#define AMAX(simpar) 	((simpar)->timeinfo.amax)
#define ANOW(simpar) 	((simpar)->timeinfo.anow)
#define ASTEP(simpar) 	((simpar)->timeinfo.astep)
#define ZINIT(simpar) 	((simpar)->timeinfo.zinit)
#define THETA(simpar) 	((simpar)->physics.gravinfo.theta)
#define RSPHERE(simpar) 	((simpar)->physics.gravinfo.rsphere)
#define NSTEP(simpar) 	((simpar)->timeinfo.nstep)
#define PMSTATUS(simpar) 	((simpar)->timeinfo.pmstatus)
#define NSKIP(simpar) 	((simpar)->timeinfo.nskip)
#define STEPCOUNT(simpar) 	((simpar)->timeinfo.stepcount)
#define STEPNUM(simpar) 	((simpar)->timeinfo.stepnum)
#define IndT_ISUBSTEP(simpar) 	((simpar)->timeinfo.indt.isubstep)
#define AGN_BP(simpar) ((simpar)->bp.agn.p.bp)

#define TOTMASS(simpar) 	((simpar)->simmodel.staticworld.TotMass)
#define FFTIME(simpar) 	((simpar)->simmodel.staticworld.ffTime)
#define FFTELLB(simpar) 	((simpar)->simmodel.staticworld.ffTimeellb)

#define NX(simpar) 	((simpar)->fftwgrid.gridinfo.nx)
#define NY(simpar) 	((simpar)->fftwgrid.gridinfo.ny)
#define NZ(simpar) 	((simpar)->fftwgrid.gridinfo.nz)
#define NXNY(simpar) 	((simpar)->fftwgrid.gridinfo.nxny)
#define NSPACE(simpar) 	((simpar)->fftwgrid.gridinfo.nspace)

#define ISEED(simpar) 	((simpar)->simmodel.cosmos.ic.iseed)
#define FLAG_IC(simpar) 	((simpar)->simmodel.cosmos.ic.flagIC)
#define GRAFIC_DIRECTORY(simpar) 	((simpar)->simmodel.cosmos.ic.graficDir)
#define PORDER(simpar) 	((simpar)->simmodel.cosmos.ic.porder)

#define COS_PAMP(simpar,itype) ((simpar)->simmodel.cosmos.ic.pamp[itype])
#define COS_DAMP1(simpar,itype) ((simpar)->simmodel.cosmos.ic.damp1[itype])
#define COS_DAMP2(simpar,itype) ((simpar)->simmodel.cosmos.ic.damp2[itype])
#define COS_VAMP1(simpar,itype) ((simpar)->simmodel.cosmos.ic.vamp1[itype])
#define COS_VAMP2(simpar,itype) ((simpar)->simmodel.cosmos.ic.vamp2[itype])

#define COS_GROWTH(simpar) ((simpar)->simmodel.cosmos.ic.growthfactor)

#define POWREADFLAG(simpar) 	((simpar)->simmodel.cosmos.powreadflag)
#define POWFILENAME(simpar) 	((simpar)->simmodel.cosmos.powfilename)

#define INPAPKFILENAME(simpar) 	((simpar)->simmodel.cosmos.inpapkfilename)


#define EXTERNALFORCE(simpar) ((simpar)->simmodel.staticworld.externalforce)



#define GRIDINFO(simpar) ((simpar)->fftwgrid.gridinfo)

#define NPIVOT(simpar,j) ( ((simpar)->simmodel.simbox.dm_ddinfo)[j].npivot )
#define PIVOT(simpar,j) ( ((simpar)->simmodel.simbox.dm_ddinfo)[j].pivot )
#define NDDINFO(simpar) ((simpar)->simmodel.simbox.nddinfo)

#define DM_DDINFO(simpar) ((simpar)->simmodel.simbox.dm_ddinfo)
#define SPH_DDINFO(simpar) ((simpar)->simmodel.simbox.sph_ddinfo)
#define VORO_DDINFO(simpar) ((simpar)->simmodel.simbox.voro_ddinfo)
#define VORORK4_DDINFO(simpar) ((simpar)->simmodel.simbox.vorork4_ddinfo)
#define STAR_DDINFO(simpar) ((simpar)->simmodel.simbox.star_ddinfo)
#define AGN_DDINFO(simpar) ((simpar)->simmodel.simbox.agn_ddinfo)
#define TDM_DDINFO(simpar) ((simpar)->simmodel.simbox.tdm_ddinfo)
#define TSPH_DDINFO(simpar) ((simpar)->simmodel.simbox.tsph_ddinfo)
#define TVORO_DDINFO(simpar) ((simpar)->simmodel.simbox.tvoro_ddinfo)
#define TVORORK4_DDINFO(simpar) ((simpar)->simmodel.simbox.tvorork4_ddinfo)
#define TSTAR_DDINFO(simpar) ((simpar)->simmodel.simbox.tstar_ddinfo)
#define TAGN_DDINFO(simpar) ((simpar)->simmodel.simbox.tagn_ddinfo)

#define DM_DDFUNC(simpar) ((simpar)->simmodel.simbox.dm_ddfunc)
#define SPH_DDFUNC(simpar) ((simpar)->simmodel.simbox.sph_ddfunc)
#define VORO_DDFUNC(simpar) ((simpar)->simmodel.simbox.voro_ddfunc)
#define VORORK4_DDFUNC(simpar) ((simpar)->simmodel.simbox.vorork4_ddfunc)
#define STAR_DDFUNC(simpar) ((simpar)->simmodel.simbox.star_ddfunc)
#define AGN_DDFUNC(simpar) ((simpar)->simmodel.simbox.agn_ddfunc)

#define TDM_DDFUNC(simpar) ((simpar)->simmodel.simbox.tdm_ddfunc)
#define TSPH_DDFUNC(simpar) ((simpar)->simmodel.simbox.tsph_ddfunc)
#define TVORO_DDFUNC(simpar) ((simpar)->simmodel.simbox.tvoro_ddfunc)
#define TVORORK4_DDFUNC(simpar) ((simpar)->simmodel.simbox.tvorork4_ddfunc)
#define TSTAR_DDFUNC(simpar) ((simpar)->simmodel.simbox.tstar_ddfunc)
#define TAGN_DDFUNC(simpar) ((simpar)->simmodel.simbox.tagn_ddfunc)

#define COS_SIMBOX(simpar) ((simpar)->simmodel.simbox.simbox)
#define KH_SIMBOX(simpar) ((simpar)->simmodel.simbox.simbox)
#define RT_SIMBOX(simpar) ((simpar)->simmodel.simbox.simbox)
#define KP_SIMBOX(simpar) ((simpar)->simmodel.simbox.simbox)
#define GL2D_SIMBOX(simpar) ((simpar)->simmodel.simbox.simbox)
#define SIMBOX(simpar) ((simpar)->simmodel.simbox.simbox)
#define STAT_SIMBOX(simpar) ((simpar)->simmodel.simbox.simbox)



#define DM_MASS(simpar,j) ((simpar)->bp.dm.mass)
#define SPH_MASS(simpar,j) ((simpar)->bp.sph.p.bp[j].mass)
#define VORO_MASS(simpar,j) ((simpar)->bp.voro.p.bp[j].mass)
#define VORORK4_MASS(simpar,j) ((simpar)->bp.voro.p.bp[j].mass)
#define STAR_MASS(simpar,j) ((simpar)->bp.star.p.bp[j].mass)
#define AGN_MASS(simpar,j) ((simpar)->bp.agn.p.bp[j].mass)

#define TDM_MASS(simpar,j) ((simpar)->bp.dm.mass)
#define TSPH_MASS(simpar,j) ((simpar)->bp.sph.p.tbp[j].mass)
#define TVORO_MASS(simpar,j) ((simpar)->bp.voro.p.tbp[j].mass)
#define TVORORK4_MASS(simpar,j) ((simpar)->bp.voro.p.tbprk4[j].mass)
#define TSTAR_MASS(simpar,j) ((simpar)->bp.star.p.tbp[j].mass)
#define TAGN_MASS(simpar,j) ((simpar)->bp.agn.p.tbp[j].mass)

/***************************************************************************************************/
/* How to handle index and flag union */
/* index range  -4611686018427387903 < range(indx) <= 4611686018427387904 */
#define FLAG_RESET (0x00)
#define DMflag (1<<0)
#define SPHflag (1<<1)
#define STARflag (1<<2)
#define AGNflag (1<<3)
#define FoFflag (1<<4)
#define BoundaryGhostflag (1<<5)
#define VOROflag (1<<6)
#define NULL2  (1<<7)

/*
#define ENDIAN_OFFSET 0
*/
#define ENDIAN_OFFSET 7

#define CLEAR_P_FLAG(simpar, TYPE, i,flag)    (TYPE##_TBP(simpar)[i].u4if.Flag[ENDIAN_OFFSET] =   FLAG_RESET )
#define SET_P_FLAG(simpar, TYPE, i,flag)    (TYPE##_TBP(simpar)[i].u4if.Flag[ENDIAN_OFFSET] |=   flag )
#define UNSET_P_FLAG(simpar, TYPE, i,flag)  (TYPE##_TBP(simpar)[i].u4if.Flag[ENDIAN_OFFSET]  &= (~flag) )
#define IS_P_FLAG(simpar, TYPE, i,flag)     (TYPE##_TBP(simpar)[i].u4if.Flag[ENDIAN_OFFSET]  &    flag )
#define TOGGLE_P_FLAG(simpar, TYPE, i,flag) (TYPE##_TBP(simpar)[i].u4if.Flag[ENDIAN_OFFSET]  ^=   flag )

#define CLEAR_FLAG(p)    ((p)->u4if.Flag[ENDIAN_OFFSET] =   FLAG_RESET )
#define SET_FLAG(p, flag)  ((p)->u4if.Flag[ENDIAN_OFFSET] |= flag)
#define UNSET_FLAG(p, flag)  ((p)->u4if.Flag[ENDIAN_OFFSET] &= (~flag))

#define IS_FLAG_ON(p, flag)  ((p)->u4if.Flag[ENDIAN_OFFSET] &  flag)
#define IS_FLAG(a, b)  IS_FLAG_ON(a,b)


#define TOGGLE_FLAG(p, flag)  ((p)->u4if.Flag[ENDIAN_OFFSET] ^= flag)


#define INDX(simpar, T, TYPE, i) ((TYPE##_##T##BP(simpar)[i].u4if.indx << 8) >> 8)
#define PINDX(p) (((p)->u4if.indx <<8) >>8)

#define CHANGEINDX(p, a) do{\
	char _b = ((p)->u4if.Flag)[ENDIAN_OFFSET];\
	(p)->u4if.indx = (a);\
	((p)->u4if.Flag)[ENDIAN_OFFSET] = _b;\
}while(0)



/***************************************************************************************************/

#define XofBp(p) ((p)->x)
#define YofBp(p) ((p)->y)
#define ZofBp(p) ((p)->z)

#ifdef XYZDBL
#		define XofP(simpar,p) fmod((double)((p)->x)+(double)((PINDX(p)%NX(simpar)))+NX(simpar),NX(simpar))
#		define YofP(simpar,p) fmod((double)((p)->y)+(double)(((PINDX(p)%NXNY(simpar))\
				/NX(simpar)))+NY(simpar),NY(simpar))
#		define ZofP(simpar,p) fmod((double)((p)->z)+(double)((PINDX(p)/NXNY(simpar)))+NZ(simpar),NZ(simpar))

#		define Xoftype(simpar,type,p) fmod((double)(((type*)p)->x)+(double)((PINDX((type*)p)%NX(simpar)))+NX(simpar),NY(simpar))
#		define Yoftype(simpar,type,p) fmod((double)(((type*)p)->y)+(double)(((PINDX((type*)p)%NXNY(simpar))/NX(simpar)))+NY(simpar),NY(simpar))
#		define Zoftype(simpar,type,p) fmod((double)(((type*)p)->z)+(double)((PINDX((type*)p)/NXNY(simpar)))+NZ(simpar),NZ(simpar))
#else
#	define XofP(simpar, p) ((p)->x)
#	define YofP(simpar, p) ((p)->y)
#	define ZofP(simpar, p) ((p)->z)
#	define Xoftype(simpar, type,p) (((type*)p)->x)
#	define Yoftype(simpar, type,p) (((type*)p)->y)
#	define Zoftype(simpar, type,p) (((type*)p)->z)
#endif

#define MofP(p,lvsimpar) (lvsimpar.mass)




#ifndef MPI_INT64_T
#define MPI_INT64_T MPI_LONG
#endif



void HAMB(SimParameters *);
int nullfct00(void);
int nullfct01(void);
