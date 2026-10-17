/* It prohibits to include "ost.h" and "nnost.h" at the same time */
#define MAX_NUM_NEAR 128
#define MINWIDTH 1.e-6
#define sqr(A) ((A)*(A))
#define NODE 0
#define PARTICLE 1
#define NODE_HAVE_PARTICLE 5
#define MIN_CELL_PARTICLE_NUM 10
#define MINCELLWIDTH 1.E-6

#define RECURSIVE 1
#define PTHREAD 0
#define SERIALIZED -1
/*
#define YES 1
#define NO 0
*/

typedef float PosType;

/*      */
#define GENERAL_TYPE unsigned int type: 1
#define GENERALTPtlPOINTER GENERAL_TYPE;void *sibling
/*      */

enum {TYPE_TREE = 0,TYPE_PTL = 1};
typedef struct TYPE {
	GENERAL_TYPE;
} TYPE;
typedef struct GENERAL_TPtl_POINTER {
	GENERALTPtlPOINTER;
} GENERAL_TPtl_POINTER;



typedef struct TPtlStruct{
	GENERALTPtlPOINTER;
	PosType mass,x,y,z;
	float vx,vy,vz;
	size_t indx; // This is to identify self referencing.
} TPtlStruct;
typedef struct TStruct{
	GENERALTPtlPOINTER;
	void *daughter;
	float dist,mass,nodesize;
	int Nparticle;
	PosType monox,monoy,monoz;
} TStruct;





/*      */
typedef struct DMParticle{
	float x,y,z;
	float vx,vy,vz;
} DMParticle;

typedef struct particle {
	dptype x,y,z;
} particle;
typedef struct pforce {
	float x,y,z;
} pforce;

void Make_NN_Tree(TStruct *,size_t, TPtlStruct *,size_t ,int );
TStruct *divide_nn_node(TStruct *,TStruct *, int); 
int Find_Near(particle *,int ,TStruct *, TPtlStruct *,PosType *, long long *, float *, float *);
int Find_Near2(particle *, int , TStruct *, TPtlStruct *, size_t , PosType *, long long *, float *, float*);
