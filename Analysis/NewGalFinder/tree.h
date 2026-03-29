#define sqr(A) ((A)*(A))
#define NODE 0
#define PARTICLE 1
#define NODE_HAVE_PARTICLE 12
#define EPSILON (0.0)
#define MinNumMem 20
#define MaxLinkedParticles 1000000

#define minstarmass 1.e6
#define INDXTYPE long long
#ifdef XYZDBL
#	define POSTYPE double
#else
#	define POSTYPE float
#endif

// imported from nnost.h
#define MINCELLWIDTH 1.E-6
//#define MAX_NUM_NEAR 128
#define MAX_NUM_NEAR 32
#define MIN_CELL_PARTICLE_NUM 10
#define RECURSIVE 1
#define PTHREAD 0
#define SERIALIZED -1
#define yes 1
#define no 0
#define NEW 2
/*
#define ENCLOSE 2
*/

enum where {OUT=0, IN=1, CROSS=2};
enum boolean {YES,NO, ENCLOSE};


#define divideThisNode(thisNode,nparticles) ((thisNode->nodesize > 0.5*MINCELLWIDTH ? yes:no) && (nparticles>=MIN_CELL_PARTICLE_NUM ? yes:no))
typedef struct FoFPosition {
    POSTYPE x,y,z;
} FoFPosition;







enum {TYPE_TREE = 0,TYPE_PTL = 1, TYPE_STAR = 1, TYPE_DM=2, TYPE_GAS=3, 
	TYPE_SINK=4, TYPE_AGN=4, TYPE_ALL=5};


typedef struct Box{
	dptype x,y,z,width;
} Box;

typedef struct HaloQ{
    size_t np,npstar,npgas,npdm,npsink;
    POSTYPE x,y,z;
    float mass, mstar,mgas,mdm,msink;
    float vx,vy,vz;
}HaloQ;

/*      */
#define GENERAL_TYPE unsigned int type
#define GENERALTPtlPOINTER GENERAL_TYPE;void *sibling
/*      */
typedef struct TYPE {
	GENERAL_TYPE;
} TYPE;
typedef struct GENERAL_TPtl_POINTER {
	GENERALTPtlPOINTER;
} GENERAL_TPtl_POINTER;
typedef struct TPtlStruct{
	GENERALTPtlPOINTER;
	dptype x,y,z;
	dptype mass;
	int indx;
} TPtlStruct;

typedef struct TStruct{
	GENERALTPtlPOINTER;
	void *daughter;
	dptype L;
	dptype dist,dist2;
	dptype x0,y0,z0;
	dptype mass;
	int Nparticle;
	dptype monox,monoy,monoz;
	dptype trQ;
	dptype quad[6];
	dptype dist_over_thetasq;
} TStruct;

typedef struct BeginEndTree{
	TStruct *start;
	void *EndTree;
	void *EndPtl;
} BeginEndTree;
/*      */
typedef struct DMParticle{
	dptype x,y,z;
	dptype vx,vy,vz;
} DMParticle;


typedef struct FoFTStruct{
	unsigned int type;
	void *sibling;
	POSTYPE L;
	POSTYPE dist2;
	POSTYPE dist;
	POSTYPE x0,y0,z0;
	POSTYPE monox,monoy,monoz;
	int Nparticle;
	void *daughter;
	float maxlink02, minlink02;
} FoFTStruct;
typedef struct FoFBeginEndTree{
	FoFTStruct *start;
	void *EndTree;
	void *EndPtl;
} FoFBeginEndTree;



/* for Fof */
typedef struct FoFTPtlStruct{
	unsigned int type;
	void *sibling;
	enum boolean included;
	dptype x,y,z,vx,vy,vz,mass;
	float link02;
	union {
		DmType dm;
		StarType star;
		SinkType sink;
		GasType gas;
	}p;
	size_t haloindx;
	int indx;
} FoFTPtlStruct;

/*
typedef struct FoFTPtlStruct particle;
*/
typedef struct particle{
	dptype x,y,z;
	dptype link02;
	int indx;
}particle;
/*
typedef particle FoFTPtlStruct;
*/


typedef struct SimpleBasicParticleType{
	unsigned int type;
	dptype mass,x,y,z;
	dptype vx,vy,vz;
	dptype link02;
	idtype indx;
	struct SimpleBasicParticleType *bp;
} SimpleBasicParticleType;

typedef struct LinkedListGrid{
	SimpleBasicParticleType *bp;
	int np;
} LinkedListGrid;

typedef struct ompFoFParticleType{
	dptype x,y,z, link02;
	idtype indx,imother;
} ompFoFParticleType;

typedef struct pforce {
	dptype x,y,z;
} pforce;

#define EraseFromTree(optr,ptr,nptr) do {\
	switch(((TYPE*)optr)->type) {\
		case TYPE_TREE:\
			if(((FoFTStruct*)optr)->daughter == ptr) \
				((FoFTStruct*)optr)->daughter = nptr;\
			else ((FoFTStruct*)optr)->sibling = nptr;\
			break;\
		default :\
			((FoFTPtlStruct*)optr)->sibling = nptr;\
	}\
}while(0)


void treeforce(particle*,float, TStruct *,TPtlStruct *,pforce *);
void Make_Tree(TStruct *,size_t, TPtlStruct *,size_t ,float, int );
float treeplumpotential(particle*,float, TStruct *,TPtlStruct *);
TStruct *divide_node(TStruct *,TStruct *, float , int ); 
BeginEndTree divide_node_Near(TStruct *,TStruct *, TPtlStruct *, Box ,TStruct *); 

typedef struct HaloBound{
	size_t nmem;
	POSTYPE zmin,zmax;
	int boundflag;
	FoFTPtlStruct *sibling;
} HaloBound;

typedef struct Coretype{
    int peak;
    int nummem, numstar;
    float starmass;
    float coredensity;
    float cx,cy,cz;
    float cvx,cvy,cvz;
    float Rtidal,density;
    unsigned char flag;
}Coretype;



FoFTStruct *FoF_divide_node(FoFTStruct *,FoFTStruct *, int);
void FoF_Make_Tree(FoFTStruct *, size_t, FoFTPtlStruct *,size_t ,int );
int new_fof_link(particle*,POSTYPE, FoFTStruct *, FoFTPtlStruct *,particle *);
int destroy_new_fof_link(particle*,POSTYPE, FoFTStruct *, FoFTPtlStruct *,particle *);
void destroy_omp_fof_link(particle *,POSTYPE ,FoFTStruct *, FoFTPtlStruct *);
size_t pnew_fof_link(particle*,POSTYPE, FoFTStruct *, FoFTPtlStruct *,particle *,
		size_t nhalo, POSTYPE,POSTYPE,POSTYPE);
TStruct *calThreadFreeNodeStart(size_t , size_t , size_t ,size_t , TStruct *, TStruct *);
FoFTStruct *calFoFThreadFreeNodeStart(size_t , size_t , size_t ,size_t , FoFTStruct *, FoFTStruct *);

void old_Make_Tree_Near(TStruct *, TPtlStruct *, int , Box);
