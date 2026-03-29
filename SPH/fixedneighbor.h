#define No_More_Neigh -1
typedef struct NeighborList{
	treesphparticletype *pi,*pj;
	long long nj, njnow;
	int pid;
}NeighborList;
NeighborList *neigh;

#define insert_local_neighbor(simpar,p,i,pj) do { \
	NeighborList *bp = neighborlist+i;\
	bp->pid = MYID(simpar);\
	bp->nj = pj-SPH_TBP(simpar);\
	bp->njnow = pj-SPH_TBP(simpar);\
	bp->pj = pj;\
	bp->pi = p;\
} while(0)
#define insert_other_neighbor(simpar,p,i,pj) do { \
	NeighborList *bp = neighborlist +i;\
	bp->pid = pj->pid;\
	bp->pj = pj->bp;\
	bp->pi = p;\
} while(0)

int determine_fixed_timestep(SimParameters *);
void InitFixedNeighborList(SimParameters *);
void FreeFixedNeighborList(SimParameters *);
void InsertFixedNeighborList(SimParameters *, treesphparticletype *,int , TPtlStruct **);
void GetSphParticleAddress(SimParameters *);
void GetInfoOfNeighborWorkList(SimParameters *);
void SendNeighborList(SimParameters *);
void DoSphDenFixedNeighbor(SimParameters *);
void DoSphComFixedNeighbor(SimParameters *);





