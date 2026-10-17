#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<float.h>
#include<unistd.h>
#ifdef _OPENMP
#include<omp.h>
#endif

#include "eunha.h"
#include "nnost.h"


int nullfct0();
int nullfct1();


#ifndef _OPENMP
#define Omp_get_thread_num() nullfct0()
#define Omp_get_num_threads() nullfct1()
#else
#define Omp_get_thread_num() omp_get_thread_num()
#define Omp_get_num_threads() omp_get_num_threads()

#endif


#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

/*
#define YES 1
#define NO 0
*/


#define SubCellDicision(a,b) ((a)>(b)? 1:0)
#define DivideNode(ThisNode,nparticles) ((ThisNode->nodesize > 0.5*MINCELLWIDTH ? YES:NO) && (nparticles>=MIN_CELL_PARTICLE_NUM ? YES:NO))



TStruct *divide_nn_node(TStruct *ThisNode, TStruct *SpareNode, int recursiveflag)
{
	TStruct *p2tree, tmpnode[8];
	void *LeftSibling;
	TPtlStruct *p2ptl, *tmpptr, *tmpptr2, *nodeparticles;
	size_t i, j, k, mnode, mx,my,mz;
	PosType tmpx,tmpy,tmpz,tmpdist2, distmax;
	float ptlmass;

	nodeparticles = ThisNode->daughter;
	ThisNode->type = TYPE_TREE;
	ThisNode->mass = 0;
	ThisNode->monox = ThisNode->monoy = ThisNode->monoz = 0;

	for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
		ptlmass = p2ptl->mass;
		ThisNode->mass += ptlmass;
		ThisNode->monox += ptlmass * p2ptl->x;
		ThisNode->monoy += ptlmass * p2ptl->y;
		ThisNode->monoz += ptlmass * p2ptl->z;
	}
	ThisNode->monox /= ThisNode->mass;
	ThisNode->monoy /= ThisNode->mass;
	ThisNode->monoz /= ThisNode->mass;
	distmax = -1.E20;
	for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
		tmpx = p2ptl->x - ThisNode->monox;
		tmpy = p2ptl->y - ThisNode->monoy;
		tmpz = p2ptl->z - ThisNode->monoz;
		tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		distmax = MAX(distmax, tmpdist2);
	}
	ThisNode->nodesize = sqrt(distmax);

	for(i=0;i<8;i++) {
		tmpnode[i].sibling = tmpnode[i].daughter = NULL;
		tmpnode[i].Nparticle = 0;
	}
	for(p2ptl=nodeparticles;p2ptl;) {
		mx = SubCellDicision(p2ptl->x,ThisNode->monox);
		my = SubCellDicision(p2ptl->y,ThisNode->monoy);
		mz = SubCellDicision(p2ptl->z,ThisNode->monoz);
		mnode = mx + 2*my + 4*mz;
		tmpnode[mnode].Nparticle ++;
		tmpptr = tmpnode[mnode].daughter;
		tmpptr2 = p2ptl->sibling;
		tmpnode[mnode].daughter = p2ptl;
		p2ptl->sibling = tmpptr;
		p2ptl = tmpptr2;
	}
	for(i=0;i<8;i++) if(tmpnode[i].Nparticle >0) break;

	TStruct *FirstDaughter,*NowDaughter;
	FirstDaughter = NowDaughter =  SpareNode;
	/* Link to the first daughter or first particle */
	if( DivideNode(ThisNode, tmpnode[i].Nparticle) )
		ThisNode->daughter = (void*)FirstDaughter;
	else 
		ThisNode->daughter = (void*) tmpnode[i].daughter;

	LeftSibling = NULL;

	for(i=0;i<8;i++) {
		if( DivideNode(ThisNode, tmpnode[i].Nparticle) ){
			NowDaughter->daughter = tmpnode[i].daughter;
			NowDaughter->Nparticle = tmpnode[i].Nparticle;
			if(LeftSibling) ((GENERAL_TPtl_POINTER*) LeftSibling)->sibling = NowDaughter;
			LeftSibling = NowDaughter;
			NowDaughter++;
		}
		else if(tmpnode[i].Nparticle >0) {
			tmpptr = tmpnode[i].daughter;
			if(LeftSibling) ((GENERAL_TPtl_POINTER*)LeftSibling)->sibling = tmpptr;
			for(;tmpptr;tmpptr=tmpptr->sibling) LeftSibling = tmpptr;
		}
	}
	((GENERAL_TPtl_POINTER *)LeftSibling)->sibling = ThisNode->sibling;
	SpareNode = NowDaughter;


	if(recursiveflag == RECURSIVE){
		TStruct *NextJobNode = FirstDaughter;
		for(i=0;i<8;i++){
			if(DivideNode(ThisNode, tmpnode[i].Nparticle) ){
				SpareNode = divide_nn_node(NextJobNode,SpareNode, recursiveflag);
				NextJobNode ++;
			}
		}
	}
	return SpareNode;
}

void Make_NN_Tree(TStruct *TREE_START, size_t nnode, TPtlStruct *ptl, size_t np, int recursiveflag){
	size_t i;
	TStruct *SpareNode = TREE_START+1;
	TREE_START->sibling = NULL;
	for(i=0;i<np;i++) {
		ptl[i].sibling = ptl + i +1;
		ptl[i].type = TYPE_PTL;
	}
	ptl[np-1].sibling = NULL;
	TREE_START->daughter = &(ptl[0]);
	TREE_START->Nparticle = np;
	if(recursiveflag == RECURSIVE) SpareNode = divide_nn_node(TREE_START, SpareNode, recursiveflag);
	else if(recursiveflag == SERIALIZED){
		TStruct *work;
		for(work=TREE_START;SpareNode-work >0; work++){
			SpareNode = divide_nn_node(work, SpareNode, SERIALIZED);
		}
	}
	else if(recursiveflag == PTHREAD) 
	{
		TStruct *work = TREE_START;
		do{
			SpareNode = divide_nn_node(work, SpareNode, PTHREAD);
			work ++;
		}
		while( work < SpareNode && (SpareNode-work) <= 50);
		size_t twork = (SpareNode-work);

#ifdef _OPENMP
#pragma omp parallel 
#endif
		{
			int pid = Omp_get_thread_num();
			int npid = Omp_get_num_threads();
			size_t mys, myf;
			size_t worksize = (twork + npid - 1)/npid;
			size_t j;
			mys = worksize *pid;
			myf = MIN(worksize *(pid+1), twork);
			TStruct *threadSpareNode = SpareNode + pid * ( (nnode - (SpareNode-TREE_START))/ npid);
			for(j=mys;j<myf;j++)
			{
				threadSpareNode = divide_nn_node(work+j,threadSpareNode,RECURSIVE);
			}
		}
	}
}



int near_open(particle *point, TStruct *tree, int npneigh, PosType maxdist , int Num_neighbor){
	PosType tmpx,tmpy,tmpz,dist2, r2, dist, r, sortdist;
	if(npneigh >= Num_neighbor) {
		tmpx = point->x - tree->monox;
		tmpy = point->y - tree->monoy;
		tmpz = point->z - tree->monoz;
		r = tree->nodesize;
		dist = sqrt(tmpx*tmpx+ tmpy*tmpy + tmpz*tmpz);
		if(dist-r > maxdist) return NO;
		else return YES;
	}
	else return YES;
}

/**********************************************************************************************************************/
/**********************************************************************************************************************/
/*******                NEAREST NEIGHBOR SEARCHING ALGORITHM         **************************************************/
/**********************************************************************************************************************/
/**********************************************************************************************************************/

typedef struct nearestneighbor{
	PosType dist2;
	TPtlStruct *bp;
} Neighbor;

typedef struct gpu_nearestneighbor{
	PosType dist2;
	int ibp;
} GPU_Neighbor;



#define INSERT(dist2, maxdist, maxdist2, npneigh, neighbor) do{\
	for(i=0;i<npneigh;i++) if(neighbor[i].dist2>dist2) break;\
	for(j=npneigh-1;j>=i;j--) neighbor[j+1] = neighbor[j];\
	neighbor[i].dist2 = dist2;\
	neighbor[i].bp = (TPtlStruct*)ptr;\
	npneigh ++;\
	npneigh=MIN(Num_neighbor,npneigh);\
	maxdist2 = (neighbor[npneigh-1].dist2);\
	maxdist = sqrt(maxdist2);\
}while(0)


int Find_Near(particle *point, int Num_neighbor, TStruct *tree, TPtlStruct *ptl, PosType *maxr, TPtlStruct **bpneighbor){
	int i,j,k;
	PosType dist2, tmpx,tmpy,tmpz;
	PosType maxdist, maxdist2;
	void *ptr;
	int npneigh=0;
	Neighbor neighbor[MAX_NUM_NEAR]={0};
	maxdist = maxdist2 = 1.E23;
	ptr = (void*)tree;
	while(ptr != NULL){
		switch( ((TYPE*)ptr)->type) {
			case TYPE_TREE:
				switch( near_open(point, ptr, npneigh, maxdist, Num_neighbor)){
					case YES:
						ptr = (void *)(((TStruct*)ptr)->daughter);
						break;
					default:
						ptr = (void *)(((TStruct*)ptr)->sibling);
				}
				break;
			default:
				tmpx = point->x - ((TPtlStruct *)ptr)->x;
				tmpy = point->y - ((TPtlStruct *)ptr)->y;
				tmpz = point->z - ((TPtlStruct *)ptr)->z;
				dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
				if(npneigh < Num_neighbor || dist2<maxdist2){
					INSERT(dist2, maxdist, maxdist2, npneigh, neighbor);
				}
				ptr = (void *) (((TPtlStruct*)ptr)->sibling);
		}
	}
	*maxr = maxdist;
	for(i=0;i<npneigh;i++) bpneighbor[i] = neighbor[i].bp;
	return npneigh;
}

#define DIRECT_INSERT(dist2,maxdist,maxdist2,npneigh,neighbor,grv) do{\
	for(i=0;i<npneigh;i++) if(neighbor[i].dist2>dist2) break;\
	for(j=npneigh-1;j>=i;j--) neighbor[j+1] = neighbor[j];\
	neighbor[i].dist2 = dist2;\
	neighbor[i].bp = grv;\
	npneigh ++;\
	npneigh=MIN(Num_neighbor,npneigh);\
	maxdist2 = (neighbor[npneigh-1].dist2);\
	maxdist = sqrt(maxdist2);\
}while(0)



int Direct_Find_Near(particle *point, int Num_neighbor, TPtlStruct *grv, int ngrv, PosType *maxr, TPtlStruct **bpneighbor){
	int i,j,ii;
	int npneigh = 0;
	PosType maxdist, maxdist2;
	Neighbor neighbor[MAX_NUM_NEAR]={0};
	maxdist = maxdist2 = 1.E23;
	for(ii=0;ii<ngrv;ii++){
		PosType tmpx,tmpy,tmpz, dist2;
		tmpx = point->x - grv[ii].x;
		tmpy = point->y - grv[ii].y;
		tmpz = point->z - grv[ii].z;
		dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		if(npneigh < Num_neighbor || dist2 < maxdist2){
			DIRECT_INSERT(dist2,maxdist,maxdist2,npneigh,neighbor,grv+ii);
		}
	}
	*maxr = maxdist;
	for(i=0;i<npneigh;i++) bpneighbor[i] = neighbor[i].bp;
	return npneigh;
}

int Find_Near2(particle *point, int Num_neighbor, TStruct *tree, TPtlStruct *ptl, size_t nptl, 
		PosType *maxr, TPtlStruct **bpneighbor){
	size_t i,j,k;
	size_t iptl = 0;
	PosType dist2, tmpx,tmpy,tmpz;
	PosType maxdist, maxdist2;
	void *ptr;
	int npneigh=0;
	Neighbor neighbor[nptl];
	maxdist = maxdist2 = 1.E23;
	ptr = (void*)tree;
	while(ptr != NULL){
		switch( ((TYPE*)ptr)->type) {
			case TYPE_TREE:
				switch( near_open(point, ptr, npneigh, maxdist, Num_neighbor)){
					case YES:
						ptr = (void *)(((TStruct*)ptr)->daughter);
						break;
					default:
						ptr = (void *)(((TStruct*)ptr)->sibling);
				}
				break;
			default:
				tmpx = point->x - ((TPtlStruct *)ptr)->x;
				tmpy = point->y - ((TPtlStruct *)ptr)->y;
				tmpz = point->z - ((TPtlStruct *)ptr)->z;
				dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
				if(iptl < Num_neighbor || dist2<maxdist2){
					neighbor[iptl].dist2 = dist2;
					neighbor[iptl].bp = ptl;
					iptl ++;
					if(iptl>Num_neighbor) maxdist2 = MIN(maxdist2, dist2);
					else maxdist2 = MAX(maxdist2, dist2);
				}
				ptr = (void *) (((TPtlStruct*)ptr)->sibling);
		}
	}
	npneigh = 0;
	for(i=0;i<iptl;i++){
		if(neighbor[i].dist2 <= maxdist2) bpneighbor[npneigh++] = neighbor[i].bp;
	}
	return npneigh;
}
