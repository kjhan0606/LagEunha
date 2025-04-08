#define OneGB 1073741824L
#define quarterGB 268435456L




#define SWAP(a,b,stmp) do{\
	memcpy(stmp,a,n_size);\
	memcpy(a,b,n_size);\
	memcpy(b,stmp,n_size);\
}while(0)

#define PosNX(gridinfo) ( (gridinfo)->nx )
#define PosNY(gridinfo) ( (gridinfo)->ny )
#define PosNZ(gridinfo) ( (gridinfo)->nz )
#define PosNXNY(gridinfo) ( (gridinfo)->nxny )

#ifdef XYZDBL
#define xPos(gridinfo,p) fmod((double)((p)->x)+(double)((PINDX(p)%PosNX(gridinfo)))+PosNX(gridinfo),PosNX(gridinfo))
#define yPos(gridinfo,p) fmod((double)((p)->y)+(double)(((PINDX(p)%PosNXNY(gridinfo))\
				                                                /PosNX(gridinfo)))+PosNY(gridinfo),PosNY(gridinfo))
#define zPos(gridinfo,p) fmod((double)((p)->z)+(double)((PINDX(p)/PosNXNY(gridinfo)))+PosNZ(gridinfo),PosNZ(gridinfo))
#else
#define xPos(gridinfo,p) ( (p)->x)
#define yPos(gridinfo,p) ( (p)->y)
#define zPos(gridinfo,p) ( (p)->z)
#endif



#define DomainError 0.01
#define MAX_ITER 50
#define MPI_DenType MPI_FLOAT
#ifndef DENTYPE
typedef float DenType;
#define DENTYPE
#endif

#ifndef MPIRKS_H
#define MPIRKS_H
typedef struct GridInfo{ 
	int ix,iy,iz; 
	int jx,jy,jz; 
	ptrdiff_t nx,ny,nz,nxny,nspace; 
	ptrdiff_t npix; 
	/*
	DenType *den;
	*/
} GridInfo;

typedef struct RGridInfo{
	int ngrids;
} RGridInfo;


typedef struct PrimeNumber{
	int prime, factor;
} PrimeNumber;


enum DDdirection {newway=0, oldway = 1};
enum RMSComBuild {NewBuildCom=0, OldCom= 1};
enum Boundary {Upper=1, Lower=0};

enum PadPtl {PadPtl=1, SimPtl=0};

/*
#define min(a,b) ((a)<(b)? (a):(b))
#define max(a,b) ((a)>(b)? (a):(b))

*/



typedef struct Range{
	PosType min,max;
}Range;
typedef struct SimBoxRange{
	Range x,y,z,w;
} SimBoxRange;
typedef struct BoxMinMax{
	PosType xmin,ymin,zmin,wmin;
	PosType xmax,ymax,zmax,wmax;
}BoxMinMax;

typedef struct DoDeFunc{
	int (*xcompare)(GridInfo *,const void *, const void *);
	int (*ycompare)(GridInfo *,const void *, const void *);
	int (*zcompare)(GridInfo *,const void *, const void *);
	int (*wcompare)(GridInfo *,const void *, const void *);
	int (*xpinner)(GridInfo *,const void *, const void *, const void *, PosType, int,Range *);
	int (*ypinner)(GridInfo *,const void *, const void *, const void *, PosType, int,Range *);
	int (*zpinner)(GridInfo *,const void *, const void *, const void *, PosType, int,Range *);
	int (*wdist)(GridInfo *,const void *, const void *, const void *, PosType, int, Range *);
	int (*muladd)(GridInfo *,const void *, const void *, float,char, int);
	int (*divdir)(GridInfo *, const void *, ptrdiff_t, MPI_Comm);
	int (*insidebox)(GridInfo *,const void *, BoxMinMax *, PosType *, SimBoxRange *, int, int, enum PadPtl);
	int (*edgeptl)(GridInfo *,const void *, SimBoxRange *, PosType *);
	void (*shift2pos)(GridInfo *,const void *);
	char xyzchip;
}DoDeFunc;


int getprimenumber(int, PrimeNumber *);


typedef struct Offset{
	size_t x,y,z,w;
}Offset;


typedef struct DoDeInfo{
	MPI_Comm com;
	ptrdiff_t n_size;
	int myid,nid;
	int nsubgroup,npivot;
	int idirection;
	int subgroupsize,subgroupid;
	int (*xcompare)(GridInfo *,const void *, const void *);
	int (*ycompare)(GridInfo *,const void *, const void *);
	int (*zcompare)(GridInfo *,const void *, const void *);
	int (*wcompare)(GridInfo *,const void *, const void *);
	int (*xpinner)(GridInfo *,const void *, const void *, const void *, PosType, int,Range *);
	int (*ypinner)(GridInfo *,const void *, const void *, const void *, PosType, int,Range *);
	int (*zpinner)(GridInfo *,const void *, const void *, const void *, PosType, int,Range *);
	int (*wdist)(GridInfo *,const void *, const void *, const void *, PosType, int, Range *);
	int (*muladd)(GridInfo *,const void *, const void *, float,char, int);
	int (*insidebox)(GridInfo *,const void *, BoxMinMax *,PosType *, SimBoxRange *, int, int, enum PadPtl);
	int (*edgeptl)(GridInfo *,const void *, SimBoxRange *, PosType *);
	void (*shift2pos)(GridInfo *,const void *);
	GridInfo sublocalgrid;
	char xyzchip;
	union {
		Offset xyz;
		size_t r[4];
	}memoffset;
	char *pivot;
	union {
		BoxMinMax xyz;
		struct {
			PosType rmin[4];
			PosType rmax[4];
		}r;
	}lgroup;
} DoDeInfo;



typedef struct Padding {
	ptrdiff_t npad;
	ptrdiff_t psize;
	void *ppad;
}PadType;


#endif
/* end of MPIRKS_H */

void old_mpirks(void **, ptrdiff_t *, ptrdiff_t , DoDeFunc *, DoDeInfo *, MPI_Comm, GridInfo *, enum DDdirection,
		enum RMSComBuild);

void mpirks(void **, ptrdiff_t *, ptrdiff_t , DoDeFunc *, DoDeInfo *, MPI_Comm, GridInfo *, enum DDdirection,
		enum RMSComBuild);
void pmigrate(void **, ptrdiff_t *, DoDeInfo *, GridInfo *);
void ExtractLocalDomainVolume(DoDeInfo *, int , SimBoxRange);

void *ExtractPtls2Send(void *, ptrdiff_t , void *, ptrdiff_t , DoDeInfo *, BoxMinMax , PosType , SimBoxRange , 
		ptrdiff_t *, ptrdiff_t , GridInfo *);

void ppadding(void *, ptrdiff_t, void **, ptrdiff_t *, DoDeInfo *, int, SimBoxRange, PosType, GridInfo *, int );
void buildrmscom(DoDeFunc *, DoDeInfo *, MPI_Comm );
void gmigrate(GridInfo *, GridInfo *, DoDeInfo *, int);
void gmigrate0(GridInfo *, GridInfo *, DoDeInfo *, int);
void getingridnpix(GridInfo *);
void *getoutgridnpix(GridInfo *, DenType);
int getNextPrimeNumber(int );


void getgridnpix(GridInfo *);

int getdirection(SimBoxRange );
int getdirection2d(SimBoxRange );




#define ISUNSIGNED(a) (((typeof(a))-1) >= 0)


#define MakeDoDeFunc3D(type,ddfunc,xcompare3d,ycompare3d,zcompare3d,xpinner3d,ypinner3d,zpinner3d,muladd3d,divdir3d,InSideBox,Edgeptl,Dshift2pos,ddinfo) do{\
	ddfunc.xcompare = type##xcompare3d;\
	ddfunc.ycompare = type##ycompare3d;\
	ddfunc.zcompare = type##zcompare3d;\
	ddfunc.xpinner = type##xpinner3d;\
	ddfunc.ypinner = type##ypinner3d;\
	ddfunc.zpinner = type##zpinner3d;\
	ddfunc.muladd = type##muladd3d;\
	ddfunc.divdir = type##divdir3d;\
	ddfunc.insidebox = type##InSideBox;\
	ddfunc.edgeptl = type##Edgeptl;\
	ddfunc.shift2pos = type##Dshift2pos;\
	ddfunc.xyzchip = 'z';\
} while(0)
#define MakeDoDeFunc2D(type,ddfunc,xcompare2d,ycompare2d,xpinner2d,ypinner2d,dmuladd,\
		divdir2d,InSideBox,Edgeptl,Dshift2pos,ddinfo) do{\
	ddfunc.xcompare = type##xcompare2d;\
	ddfunc.ycompare = type##ycompare2d;\
	ddfunc.xpinner = type##xpinner2d;\
	ddfunc.ypinner = type##ypinner2d;\
	ddfunc.muladd = type##dmuladd;\
	ddfunc.divdir = type##divdir2d;\
	ddfunc.insidebox = type##InSideBox;\
	ddfunc.edgeptl = type##Edgeptl;\
	ddfunc.shift2pos = type##Dshift2pos;\
	ddfunc.xyzchip = 'y';\
} while(0)
#define MakeDoDeFunc4D(type,ddfunc,xcompare4d,ycompare4d,zcompare4d,wcompare4d,xpinner4d,ypinner4d,zpinner4d,wdist,muladd4d,divdir4d,InSideBox,\
		Edgeptl,Dshift2pos,ddinfo) do{\
	ddfunc.xcompare = type##xcompare4d;\
	ddfunc.ycompare = type##ycompare4d;\
	ddfunc.zcompare = type##zcompare4d;\
	ddfunc.wcompare = type##wcompare4d;\
	ddfunc.xpinner = type##xpinner4d;\
	ddfunc.ypinner = type##ypinner4d;\
	ddfunc.zpinner = type##zpinner4d;\
	ddfunc.muladd = type##muladd4d;\
	ddfunc.divdir = type##divdir4d;\
	ddfunc.insidebox = type##InSideBox;\
	ddfunc.edgeptl = type##Edgeptl;\
	ddfunc.shift2pos = type##Dshift2pos;\
	ddfunc.xyzchip = 'w';\
} while(0)

#define GenerateShift2Pos(Shift2Pos,ptype) void ptype##shift2pos(GridInfo *gridinfo,const void *a){\
	ptype *aa = (ptype*)a;\
	aa->x = xPos(gridinfo, aa);\
	aa->y = yPos(gridinfo, aa);\
	aa->z = zPos(gridinfo, aa);\
	aa->u4if.indx = 0;\
}

#define GenerateShift2Pos2D(Shift2Pos,ptype) void ptype##shift2pos(GridInfo *gridinfo,const void *a){\
	ptype *aa = (ptype*)a;\
	aa->x = xPos(gridinfo, aa);\
	aa->y = yPos(gridinfo, aa);\
	aa->u4if.indx = 0;\
}

#define GenerateMemberCompare(compare,ptype,member,memtype) int ptype##compare(GridInfo *gridinfo,const void *a, const void *b){\
	memtype _a,_b;\
	ptype *aa = (ptype*)a;\
	ptype *bb = (ptype*)b;\
	_a = member##Pos(gridinfo, aa);\
	_b = member##Pos(gridinfo, bb);\
	if(_a>_b) return 1;\
	else if(_a == _b) return 0;\
	else return -1;\
}

#define GenerateMemberDistance(dist,ptype,member,memtype) int ptype##dist(GridInfo *gridinfo, const void *a, const void *b, \
		PosType width, Range *L){\
	memtype _a,_b;\
	ptype *aa = (ptype*)a;\
	ptype *bb = (ptype*)b;\
	_a = member##Pos(gridinfo, aa);\
	_b = member##Pos(gridinfo, bb);\
	memtype d = fabs(_a-_b);\
	if(d <width) return 1;\
	else if(d ==width) return  0;\
	else return -1;\
} 


#define GenerateInsideRange(dist,ptype,member,memtype) int ptype##dist(GridInfo *gridinfo, \
		const void *a, const void *l, const void *u, PosType width, int iperiod, Range *L){\
	memtype _a,_u,_l;\
	ptype *aa = (ptype*)a;\
	_a = member##Pos(gridinfo, aa);\
	ptype *ll = (ptype*)l;\
	ptype *uu = (ptype*)u;\
	_l = member##Pos(gridinfo,ll) - width + iperiod*(L->max);\
	_u = member##Pos(gridinfo,uu) + width + iperiod*(L->max);\
	if(_a >= _l && _a <_u) return 1;\
	else return 0;\
} 


#define MakeDoDeInfo2D(number,prime, ptype,mx,my,Info, ninfo) do {\
    int nprime = getprimenumber(number,prime);\
    int _i,_j,_k;\
    _j=0;\
    for(_i=0;_i<nprime;_i++) _j += prime[_i].factor;\
    Info = (DoDeInfo*)malloc(sizeof(DoDeInfo)*_j);\
    _k=0;\
    for(_i=0;_i<nprime;_i++) for(_j=0;_j<prime[_i].factor;_j++){\
        Info[_k].pivot = (char *)malloc(sizeof(ptype)*(prime[_i].prime-1));\
        Info[_k].memoffset.xyz.x = offsetof(ptype,mx);\
        Info[_k].memoffset.xyz.y = offsetof(ptype,my);\
        _k++;\
    }\
    ninfo = _k;\
} while(0)


#define MakeDoDeInfo3D(number,prime, ptype,mx,my,mz,Info, ninfo) do {\
	int nprime = getprimenumber(number,prime);\
	int _i,_j,_k;\
	_j=0;\
	for(_i=0;_i<nprime;_i++) _j += prime[_i].factor;\
	Info = (DoDeInfo*)malloc(sizeof(DoDeInfo)*_j);\
	_k=0;\
	for(_i=0;_i<nprime;_i++) for(_j=0;_j<prime[_i].factor;_j++){\
		Info[_k].pivot = (char *)malloc(sizeof(ptype)*(prime[_i].prime-1));\
		Info[_k].memoffset.xyz.x = offsetof(ptype,mx);\
		Info[_k].memoffset.xyz.y = offsetof(ptype,my);\
		Info[_k].memoffset.xyz.z = offsetof(ptype,mz);\
		_k++;\
	}\
	ninfo = _k;\
} while(0)

#define MakeDoDeInfo4D(number,prime, ptype,mx,my,mz,mw,Info, ninfo) do {\
	int nprime = getprimenumber(number,prime);\
	int _i,_j,_k;\
	_j=0;\
	for(_i=0;_i<nprime;_i++) _j += prime[_i].factor;\
	Info = (DoDeInfo*)malloc(sizeof(DoDeInfo)*_j);\
	_k=0;\
	for(_i=0;_i<nprime;_i++) for(_j=0;_j<prime[_i].factor;_j++){\
		Info[_k].pivot = (char *)malloc(sizeof(ptype)*(prime[_i].prime-1));\
		Info[_k].memoffset.xyz.x = offsetof(ptype,mx);\
		Info[_k].memoffset.xyz.y = offsetof(ptype,my);\
		Info[_k].memoffset.xyz.z = offsetof(ptype,mz);\
		Info[_k].memoffset.xyz.w = offsetof(ptype,mw);\
		_k++;\
	}\
	ninfo = _k;\
} while(0)






