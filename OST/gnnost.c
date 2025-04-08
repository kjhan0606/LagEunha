/* This is a generalized version of the near neighbor searching with the OST. */
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
#include "gnnost.h"


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



#define SubCellDicision(a,b) ((a)>(b)? YES:NO)
#define furtherDivision(thisNode,nparticles) \
	((thisNode->nodesize > 0.5*MINCELLWIDTH ? YES:NO) \
	 && (nparticles>=MIN_CELL_PARTICLE_NUM ? YES:NO))


TStruct *divide_gnn_node(
		TStruct *thisNode, 
		TStruct *freeNode, 
		void (*rule2Divide)(TStruct *, TStruct *),
		void (*findCentroid)(TStruct *),
		void (*findCellSize)(TStruct *),
		int recursiveflag
		) {
	TStruct *p2tree, tmpNode[8];
	void *leftSibling;
	TPtlStruct *p2ptl, *tmpptr, *tmpptr2, *nodeparticles;
	size_t i, j, k, mnode, mx,my,mz;
	PosType tmpx,tmpy,tmpz,tmpdist2, distmax;
	float particleMass;

	thisNode->type = TYPE_TREE;
	thisNode->mass = 0;


	findCentroid(thisNode);

	/*
	nodeparticles = thisNode->daughter;
	thisNode->monox = thisNode->monoy = thisNode->monoz = 0;
	for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
		particleMass = p2ptl->mass;
		thisNode->mass += particleMass;
		thisNode->monox += particleMass * p2ptl->x;
		thisNode->monoy += particleMass * p2ptl->y;
		thisNode->monoz += particleMass * p2ptl->z;
	}
	thisNode->monox /= thisNode->mass;
	thisNode->monoy /= thisNode->mass;
	thisNode->monoz /= thisNode->mass;
	*/


	findCellSize(thisNode);
	/*
	nodeparticles = thisNode->daughter;
	distmax = -1.E20;
	for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
		tmpx = p2ptl->x - thisNode->monox;
		tmpy = p2ptl->y - thisNode->monoy;
		tmpz = p2ptl->z - thisNode->monoz;
		tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		distmax = MAX(distmax, tmpdist2);
	}
	thisNode->nodesize = sqrt(distmax);
	*/

	for(i=0;i<8;i++) {
		tmpNode[i].sibling = tmpNode[i].daughter = NULL;
		tmpNode[i].Nparticle = 0;
	}
	rule2Divide(thisNode, tmpNode);
	/*
	nodeparticles = thisNode->daughter;
	for(p2ptl=nodeparticles;p2ptl;) {
		mx = SubCellDicision(p2ptl->x,thisNode->monox);
		my = SubCellDicision(p2ptl->y,thisNode->monoy);
		mz = SubCellDicision(p2ptl->z,thisNode->monoz);
		mnode = mx + 2*my + 4*mz;
		tmpNode[mnode].Nparticle ++;
		tmpptr = tmpNode[mnode].daughter;
		tmpptr2 = p2ptl->sibling;
		tmpNode[mnode].daughter = p2ptl;
		p2ptl->sibling = tmpptr;
		p2ptl = tmpptr2;
	}
	*/



	for(i=0;i<8;i++) if(tmpNode[i].Nparticle >0) break;

	TStruct *firstDaughter,*nowDaughter;
	firstDaughter = nowDaughter =  freeNode;
	/* Link to the first daughter or first particle */
	if( furtherDivision(thisNode, tmpNode[i].Nparticle) )
		thisNode->daughter = (void*)firstDaughter;
	else 
		thisNode->daughter = (void*) tmpNode[i].daughter;

	leftSibling = NULL;

	for(i=0;i<8;i++) {
		if( furtherDivision(thisNode, tmpNode[i].Nparticle) ){
			nowDaughter->daughter = tmpNode[i].daughter;
			nowDaughter->Nparticle = tmpNode[i].Nparticle;
			if(leftSibling) ((GENERAL_TPtl_POINTER*) leftSibling)->sibling = nowDaughter;
			leftSibling = nowDaughter;
			nowDaughter++;
		}
		else if(tmpNode[i].Nparticle >0) {
			tmpptr = tmpNode[i].daughter;
			if(leftSibling) ((GENERAL_TPtl_POINTER*)leftSibling)->sibling = tmpptr;
			for(;tmpptr;tmpptr=tmpptr->sibling) leftSibling = tmpptr;
		}
	}
	((GENERAL_TPtl_POINTER *)leftSibling)->sibling = thisNode->sibling;
	freeNode = nowDaughter;


	if(recursiveflag == RECURSIVE){
		TStruct *nextJobNode = firstDaughter;
		for(i=0;i<8;i++){
			if(furtherDivision(thisNode, tmpNode[i].Nparticle) ){
				freeNode = divide_gnn_node(nextJobNode,freeNode, rule2Divide, findCentroid, 
						findCellSize, recursiveflag);
				nextJobNode ++;
			}
		}
	}
	return freeNode;
}

void Make_GNN_Tree(
		TStruct *treeStarting, 
		size_t nnode, 
		void (*rule2Divide)(TStruct *, TStruct *), /* It defines the division rule. */
		void (*findCentroid)(TStruct *), /* It calculates the center of mass . */
		void (*findCellSize)(TStruct *), /* It calculates the cell size . */
		int recursiveflag
		){
	size_t i;
	TStruct *freeNode = treeStarting+1;
	treeStarting->sibling = NULL;
	/*
	for(i=0;i<np;i++) {
		ptl[i].sibling = ptl + i +1;
		ptl[i].type = TYPE_PTL;
	}
	ptl[np-1].sibling = NULL;
	treeStarting->daughter = ptl;
	treeStarting->Nparticle = np;
	*/
	/* This is the recursive sequential code. */
	if(recursiveflag == RECURSIVE)  
		freeNode = divide_gnn_node(treeStarting, freeNode, 
				rule2Divide, findCentroid, findCellSize,
				recursiveflag);
	/* This is a code which is sequential but not recursive. */
	else if(recursiveflag == SERIALIZED){ 
		TStruct *work;
		for(work=treeStarting;freeNode-work >0; work++){
			freeNode = divide_gnn_node(work, freeNode, 
					rule2Divide, findCentroid, findCellSize,
					SERIALIZED);
		}
	}
	/* This is an OMP code triggered when 50 nodes are in the waiting list. */
	else if(recursiveflag == PTHREAD)  
	{
		TStruct *work = treeStarting;
		do{
			freeNode = divide_gnn_node(work, freeNode, 
					rule2Divide,findCentroid, findCellSize,
					PTHREAD);
			work ++;
		}
		while( work<freeNode && (freeNode-work) <= 50);
		size_t twork = (freeNode-work);

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
			TStruct *threadfreeNode = freeNode + 
				pid * ( (nnode - (freeNode-treeStarting))/ npid);
			for(j=mys;j<myf;j++)
			{
				threadfreeNode = divide_gnn_node(work+j,threadfreeNode,
						rule2Divide, findCentroid, findCellSize,
						RECURSIVE);
			}
		}
	}
}



/*
int near_open(particle *point, TStruct *tree, int npneigh, PosType maxdist , int Num_neighbor){
	PosType tmpx,tmpy,tmpz,dist2, r2, dist, r, sortdist;
	if(npneigh >= Num_neighbor) {
		tmpx = point->x = tree->monox;
		tmpy = point->y = tree->monoy;
		tmpz = point->z = tree->monoz;
		r = tree->nodesize;
		dist = sqrt(tmpx*tmpx+ tmpy*tmpy + tmpz*tmpz);
		if(dist-r > maxdist) return NO;
		else return YES;
	}
	else return YES;
}
*/

/**********************************************************************************************************************/
/**********************************************************************************************************************/
/*******                NEAREST NEIGHBOR SEARCHING ALGORITHM         **************************************************/
/**********************************************************************************************************************/
/**********************************************************************************************************************/




int insGnear(void *ptr, PosType dist2, PosType *maxdist2,  
		int npneigh, Neighbor *neighbor, int Num_neighbor) 
{
	int i,j;
	for(i=0;i<npneigh;i++) if(neighbor[i].dist2>dist2) break;
	for(j=npneigh-1;j>=i;j--) neighbor[j+1] = neighbor[j];
	neighbor[i].dist2 = dist2;
	neighbor[i].bp = ptr;
	npneigh ++;
	npneigh=MIN(Num_neighbor,npneigh);
	*maxdist2 = (neighbor[npneigh-1].dist2);
	return npneigh;
}


int find_GNear(
		void *point, 
		int Num_neighbor, 
		TStruct *tree, 
		/*
		TPtlStruct *ptl, 
		*/
		PosType *maxr, 
		void **bpneighbor,
		int (*nearOpen)(void *, TStruct *, int, PosType , int),
		PosType (*getDistance)(void *, void *),
		int (*insGnear)( void *, PosType, PosType *, int, Neighbor *, int)
		){
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
				switch( nearOpen(point, ptr, npneigh, maxdist, Num_neighbor)){
					case YES:
						ptr = (void *)(((TStruct*)ptr)->daughter);
						break;
					default:
						ptr = (void *)(((TStruct*)ptr)->sibling);
				}
				break;
			default:
				dist2 = getDistance(point, ptr);
				if(npneigh < Num_neighbor || dist2<maxdist2){
					npneigh = insGnear(ptr, dist2,  &maxdist2, npneigh, neighbor, Num_neighbor);
					maxdist = sqrt(maxdist2);
				}
				ptr = (void *) (((TPtlStruct*)ptr)->sibling);
		}
	}
	*maxr = maxdist;
	for(i=0;i<npneigh;i++) bpneighbor[i] = (void*)(neighbor[i].bp);
	return npneigh;
}

PosType find_GNearest(
        void *point,
        TStruct *tree,
        int (*nearOpen)(void *, TStruct *, PosType ),
        PosType (*getDistance)(void *, void *)
        ){
    int i,j,k;
    PosType dist2, tmpx,tmpy,tmpz;
    PosType mindist, mindist2;
    void *ptr;
    int npneigh=0;
    Neighbor neighbor[MAX_NUM_NEAR]={0};
    mindist = mindist2 = 1.E23;
    ptr = (void*)tree;
    while(ptr != NULL){
        switch( ((TYPE*)ptr)->type) {
            case TYPE_TREE:
                switch( nearOpen(point, ptr, mindist )){
                    case YES:
                        ptr = (void *)(((TStruct*)ptr)->daughter);
                        break;
                    default:
                        ptr = (void *)(((TStruct*)ptr)->sibling);
                }
                break;
            default:
				if( ((TPtlStruct*)point)->indx != ((TPtlStruct*)ptr)->indx) {
	                dist2 = getDistance(point, ptr);
					if(dist2 < mindist2){
						mindist2 = dist2;
						mindist = sqrt(mindist2);
					}
				}
                ptr = (void *) (((TPtlStruct*)ptr)->sibling);
        }
    }
    return mindist2;
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



