#define sqr(A) ((A)*(A))
#define NODE 0
#define PARTICLE 1
#define MIN_CELL_PARTICLE_NUM 10
#define MINCELLWIDTH 1.E-5
#define EPSILON (0.1)
enum where {OUT=0, IN=1, CROSS=2};
#define RECURSIVE 1
#define PTHREAD 0
#define SERIALIZED -1
#define YES 1
#define NO 0
/*      */
#define GENERAL_TYPE unsigned int type: 1
#define GENERALTPtlPOINTER GENERAL_TYPE;void *sibling
/*      */
enum {TYPE_TREE = 0,TYPE_PTL = 1};
#ifdef XYZDBL
#define PosType  double
#else
#define PosType  float
#endif

#define NodePosType float



typedef struct TYPE {
	GENERAL_TYPE;
} TYPE;
typedef struct GENERAL_TPtl_POINTER {
	GENERALTPtlPOINTER;
} GENERAL_TPtl_POINTER;

typedef struct TPtlStruct{
	GENERALTPtlPOINTER;
	PosType x,y,z;
	float mass;
	size_t indx;
} TPtlStruct;

typedef struct TStruct{
	GENERALTPtlPOINTER;
	void *daughter;
	PosType nodesize;
	float mass;
	int Nparticle;
	float dist_over_thetasq;
	PosType monox,monoy,monoz;
	float quad[6],trQ;
} TStruct;



/* for Fof */
typedef struct FoFTPtlStruct{
	GENERALTPtlPOINTER;
	PosType x,y,z;
	float mass;
	size_t haloindx;
	linkedlisttype *bp;
	int included;
} FoFTPtlStruct;

typedef struct FoFTStruct{
	GENERALTPtlPOINTER;
	void *daughter;
	PosType nodesize;
	/*
	float mass;
	PosType L;
	PosType dist;
	PosType x0,y0,z0;
	*/
	PosType monox,monoy,monoz;
	int Nparticle;
} FoFTStruct;

/*      */
typedef struct DMParticle{
	PosType x,y,z;
	float vx,vy,vz;
} DMParticle;
typedef struct Position {
	PosType x,y,z;
	float mass;
	float *axyz;
} Position;

typedef struct FoFPosition {
	PosType x,y,z;
	/*
	float mass;
	float vx,vy,vz;
	*/
} FoFPosition;

typedef struct pforce {
	PosType x,y,z;
} pforce;

#define EraseFromTree(optr,ptr,nptr) {\
	switch(((TYPE*)optr)->type) {\
		case TYPE_TREE:\
			if(((FoFTStruct*)optr)->daughter == ptr) \
				((FoFTStruct*)optr)->daughter == nptr;\
			else ((FoFTStruct*)optr)->sibling = nptr;\
			break;\
		default :\
			((FoFTPtlStruct*)ptr)->sibling = nptr;\
	}\
}\
	
void treeforce(SimParameters *, Position *, TStruct *, size_t ,float , float*);
void treeforce2(SimParameters *, Position *, TStruct *, size_t ,float , void*);
void treeforce3(SimParameters *, Position *, TStruct *, size_t ,float , void*);
void Make_Tree(TStruct *,size_t, TPtlStruct *,size_t , float, int);
void Make_Tree2(TStruct *,size_t, TPtlStruct *,size_t , float, int);
void Make_Tree3(TStruct *,size_t, TPtlStruct *,size_t , float, int);
float treeplumpotential(Position*,float, TStruct *,TPtlStruct *);
TStruct *divide_node(TStruct *,TStruct *, float, int );
TStruct *divide_node2(TStruct *,TStruct *, float, int );
TStruct *divide_node3(TStruct *,TStruct *, float, int );
size_t pnew_fof_link(FoFPosition *, FoFTStruct *, PosType , FoFPosition *, size_t );
void FoF_Make_Tree(FoFTStruct *,size_t, FoFTPtlStruct *,size_t , int);
int Find_Near(Position *,int ,TStruct *, TPtlStruct *,float *,float);
ptrdiff_t insidecount(Position , TStruct *, float );
ptrdiff_t insidecount2(Position , TStruct *, float );
ptrdiff_t insidecount3(Position , TStruct *, float );
