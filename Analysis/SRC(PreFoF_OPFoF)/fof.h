#define MinNumMem 20
#define MaxLinkedParticles 100000000 
#define sqr(A) ((A)*(A))
#define NODE 0
#define PARTICLE 1
#define NODE_HAVE_PARTICLE 20
#define EPSILON (0.1)
#define INDXTYPE long long
#ifdef XYZDBL
#	define POSTYPE double
#else
#	define POSTYPE float
#endif
/*
#define YES '0'
#define NO '1'
*/
typedef struct Box{
	POSTYPE x,y,z;
	float width;
} Box;
typedef struct HaloQ{
	size_t np;
	POSTYPE x,y,z;
	float vx,vy,vz;
}HaloQ;
/*      */
/*
#define GENERAL_TYPE unsigned int type: 1
#define GENERALTPtlPOINTER GENERAL_TYPE;void *sibling
#define GENERAL_TYPE unsigned int type
*/
/*      */
enum {TYPE_TREE = 0,TYPE_PTL = 1};
typedef struct TYPE {
	unsigned int type;
} TYPE;
typedef struct GENERAL_TPtl_POINTER {
	unsigned int type;
	void *sibling;
} GENERAL_TPtl_POINTER;
typedef struct TPtlStruct{
	unsigned int type;
	void *sibling;
	POSTYPE r[3];
	float mass;
	INDXTYPE indx;
} TPtlStruct;
/* for Fof */
typedef struct FoFTPtlStruct{
	unsigned int type;
	void *sibling;
	enum boolean included;
	POSTYPE r[3];
	float rv[3];
	INDXTYPE indx;
	/*
	TPtlStruct *pointer;
	*/
	size_t haloindx;
} FoFTPtlStruct;
typedef struct HaloBound{
	size_t nmem;
	POSTYPE zmin,zmax;
	int boundflag;
	FoFTPtlStruct *sibling;
} HaloBound;
typedef struct FoFTStruct{
	unsigned int type;
	void *sibling;
	POSTYPE L;
	POSTYPE dist2;
	POSTYPE dist;
	POSTYPE r0[3];
	POSTYPE mono[3];
	int Nparticle;
	void *daughter;
} FoFTStruct;
typedef struct FoFBeginEndTree{
	FoFTStruct *start;
	void *EndTree;
	void *EndPtl;
} FoFBeginEndTree;

typedef struct particle {
	POSTYPE x,y,z;
	float vx,vy,vz;
	INDXTYPE indx;
} particle;

/*
typedef struct readtype{
	float x,y,z;
	float vx,vy,vz;
	INDXTYPE indx;
}READTYPE;
*/
typedef dmparticletype READTYPE;

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
	
FoFBeginEndTree FoF_divide_node(FoFTStruct *,FoFTStruct *, FoFTPtlStruct *, 
		Box ,FoFTStruct *); 
void FoF_Make_Tree(FoFTStruct *,FoFTPtlStruct *,size_t ,Box );
int new_fof_link(particle*,POSTYPE, FoFTStruct *, FoFTPtlStruct *,particle *);
size_t pnew_fof_link(particle*,POSTYPE, FoFTStruct *, FoFTPtlStruct *,particle *,
		size_t nhalo, POSTYPE,POSTYPE,POSTYPE);
void CheckHaloBound(size_t ,HaloBound *,FoFTPtlStruct *,size_t , POSTYPE ,
		POSTYPE , POSTYPE, POSTYPE);
HaloQ haloproperty(particle *,size_t );
void WriteIsolatedHalo(size_t , HaloBound *,FoFTPtlStruct *, particle *, char *,
		char *);
void WriteFinalHalo(size_t , HaloBound *,FoFTPtlStruct *, particle *, char *,
		char *);
void WriteAllHalo(size_t , HaloBound *,FoFTPtlStruct *, size_t, 
		particle *, char *, char *);

size_t StackUpContactParticleLeftWard(size_t ,HaloBound *, FoFTPtlStruct *,
		size_t );
void ReadBottomFaceContact(FoFTPtlStruct *,size_t,particle *,int ,int,int);
size_t WriteBottomFaceContact(size_t , HaloBound *, FoFTPtlStruct *, particle *,
		int,int);


