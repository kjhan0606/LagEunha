#define sqr(A) ((A)*(A))
#define NODE 0
#define PARTICLE 1
#define MIN_CELL_PARTICLE_NUM 10
#define MINCELLWIDTH 1.E-5
#define EPSILON (0.1)
#define RECURSIVE 1
#define PTHREAD 0
#define ForLoop -1
/*      */
#define GENERAL_TYPE unsigned int type: 1
#define GENERALTPtlPOINTER GENERAL_TYPE;void *sibling
/*      */
enum {TYPE_TREE = 0,TYPE_PTL = 1};

#define NodePosType float



typedef struct TYPE {
	GENERAL_TYPE;
} TYPE;
typedef struct GENERAL_TPtl_POINTER {
	GENERALTPtlPOINTER;
} GENERAL_TPtl_POINTER;
/* for Fof */
typedef struct FoFTPtlStruct{
	GENERALTPtlPOINTER;
	PosType x,y,z;
	int indx;
	float mass;
	enum boolean included;
	size_t haloindx;
} FoFTPtlStruct;

typedef struct FoFTStruct{
	GENERALTPtlPOINTER;
	void *daughter;
	PosType monox,monoy,monoz;
	int Nparticle;
	float mass;
	PosType nodesize;
} FoFTStruct;

typedef struct Position {
	PosType x,y,z;
} Position;

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
	
/*
void treeforce(Position*,float, TStruct *,TPtlStruct *,pforce *);
void Make_Tree(TStruct *,int, TPtlStruct *,int , float, int);
float treeplumpotential(Position*,float, TStruct *,TPtlStruct *);
TStruct *divide_node(TStruct *,TStruct *, float, int );
*/
FoFTStruct *FoF_divide_node(FoFTStruct *,FoFTStruct *,  int); 
void FoF_Make_Tree(FoFTStruct *, int, FoFTPtlStruct *, int,  int );
size_t pnew_fof_link(Position*,FoFTStruct *, PosType, Position *, size_t);
