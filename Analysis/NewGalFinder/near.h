#ifdef NNNNNEW
#define sqr(A) ((A)*(A))
#define NODE 0
#define PARTICLE 1
#define NODE_HAVE_PARTICLE 5
#define EPSILON (0.0L)


#define INDXTYPE long long
#ifdef XYZDBL
#	define POSTYPE double
#else
#	define POSTYPE float
#endif

/*
enum boolean {YES=01, NO=02};
*/
typedef struct Box{
	double x,y,z,width;
} Box;
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
	float r[3];
	float mass;
	float vx,vy,vz;
	int indx;
} TPtlStruct;
typedef struct TStruct{
	GENERALTPtlPOINTER;
	void *daughter;
	float L;
	float dist2;
	float r0[3];
	float mass;
	int Nparticle;
	float mono[3];
	float mrr;
	/*
	float quad[6];
	*/
} TStruct;
typedef struct BeginEndTree{
	TStruct *start;
	void *EndTree;
	void *EndPtl;
} BeginEndTree;


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
}
	
#endif


void Make_Tree_Near(TStruct *,TPtlStruct *,int ,Box );
BeginEndTree divide_node_Near(TStruct *,TStruct *, TPtlStruct *, Box ,TStruct *); 
int Find_Near(particle *,int ,TStruct *, TPtlStruct *,float *, int*,float*, float *);
int Limited_Find_Near(particle *,int ,TStruct *, TPtlStruct *,float *, int*,float*, float *);
Box findsphbox(TPtlStruct *, int );
void findsphdensity(SimpleBasicParticleType *,int ,long long *, int , float *);
void starfindsphdensity(SimpleBasicParticleType *,int ,long long *, int , float *);
