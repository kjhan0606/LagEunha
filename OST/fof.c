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
#include "ost.h"


int nullfct0() { return 0;}
int nullfct1() { return 1;}
#ifndef _OPENMP
#define omp_get_thread_num() nullfct0()
#define omp_get_num_threads() nullfct1()
#endif



#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define YES 1
#define NO 0


#define SubCellDicision(a,b) ((a)>(b)? 1:0)
#define DivideNode(ThisNode,nparticles) ((ThisNode->nodesize > 0.5*MINCELLWIDTH ? YES:NO) & (nparticles>=MIN_CELL_PARTICLE_NUM ? YES:NO))

#define EraseFromTree(optr,ptr,nptr) {\
	switch(((TYPE*)optr)->type) {\
		case TYPE_TREE:\
			if(((FoFTStruct*)optr)->daughter == ptr) ((FoFTStruct*)optr)->daughter == nptr;\
			else ((FoFTStruct*)optr)->sibling = nptr;\
			break;\
		default :\
			((FoFTPtlStruct*)ptr)->sibling = nptr;\
	}\
}\




FoFTStruct *FoF_divide_node(FoFTStruct *ThisNode, FoFTStruct *SpareNode, int recursiveflag)
{
	FoFTStruct *p2tree, tmpnode[8];
	void *LeftSibling;
	FoFTPtlStruct *p2ptl, *tmpptr, *tmpptr2, *nodeparticles;
	size_t i, j, k, mnode, mx,my,mz;
	PosType tmpx,tmpy,tmpz,tmpdist2, distmax;
	float ptlmass;

	nodeparticles = ThisNode->daughter;
	ThisNode->type = TYPE_TREE;
	ThisNode->Nparticle = 0;
	ThisNode->monox = ThisNode->monoy = ThisNode->monoz = 0;

	for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
		ThisNode->Nparticle += 1;
		ThisNode->monox += p2ptl->x;
		ThisNode->monoy += p2ptl->y;
		ThisNode->monoz += p2ptl->z;
	}
	ThisNode->monox /= ThisNode->Nparticle;
	ThisNode->monoy /= ThisNode->Nparticle;
	ThisNode->monoz /= ThisNode->Nparticle;
	distmax = -1.E20;
	float q0,q1,q2,q3,q4,q5;
	q0=q1=q2=q3=q4=q5=0;
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

	FoFTStruct *FirstDaughter,*NowDaughter;
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
		FoFTStruct *NextJobNode = FirstDaughter;
		for(i=0;i<8;i++){
			if(DivideNode(ThisNode, tmpnode[i].Nparticle) ){
				SpareNode = FoF_divide_node(NextJobNode,SpareNode, recursiveflag);
				NextJobNode ++;
			}
		}
	}
	return SpareNode;
}

void FoF_Make_Tree(FoFTStruct *TREE_START, size_t nnode, FoFTPtlStruct *ptl, size_t np, int recursiveflag){
	size_t i;
	FoFTStruct *SpareNode = TREE_START+1;
	TREE_START->sibling = NULL;
	for(i=0;i<np;i++) {
		ptl[i].sibling = ptl + i +1;
		ptl[i].type = TYPE_PTL;
		ptl[i].included = NO;
	}
	ptl[np-1].sibling = NULL;
	TREE_START->daughter = &(ptl[0]);
	TREE_START->Nparticle = np;
	if(recursiveflag == RECURSIVE) SpareNode = FoF_divide_node(TREE_START, SpareNode, recursiveflag);
	else if(recursiveflag == SERIALIZED){
		FoFTStruct *work;
		for(work=TREE_START;SpareNode-work >0; work++){
			SpareNode = FoF_divide_node(work, SpareNode, SERIALIZED);
		}
	}
	else if(recursiveflag == PTHREAD) 
	{
		FoFTStruct *work = TREE_START;
		do{
			SpareNode = FoF_divide_node(work, SpareNode, PTHREAD);
			work ++;
		}
		while( (SpareNode-work) <= 50);
		size_t twork = (SpareNode-work);

#ifdef _OPENMP
#pragma omp parallel 
#endif
		{
			int pid = omp_get_thread_num();
			int npid = omp_get_num_threads();
			size_t mys, myf;
			size_t worksize = (twork + npid - 1)/npid;
			size_t j;
			mys = worksize *pid;
			myf = MIN(worksize *(pid+1), twork);
			FoFTStruct *threadSpareNode = SpareNode + pid * ( (nnode - (SpareNode-TREE_START))/ npid);
			for(j=mys;j<myf;j++)
			{
				threadSpareNode = FoF_divide_node(work+j,threadSpareNode,RECURSIVE);
			}
		}
	}
}

inline enum where FoF_InsideOpen(FoFPosition p, FoFTStruct *tree, float maxdist) {
	PosType tmpx = p.x-tree->monox;
	PosType tmpy = p.y-tree->monoy;
	PosType tmpz = p.z-tree->monoz;
	PosType dist = sqrt(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);
	if(dist+tree->nodesize <= maxdist) return IN;
	else if(dist-tree->nodesize > maxdist) return OUT;
	else return CROSS;

}



int pfof_open(FoFPosition p, FoFTStruct *tree, PosType fof_link){
	PosType tmpx, tmpy, tmpz;
	PosType dist2, dist, r, diff;
	PosType ratio;
	tmpx = fabs(p.x-tree->monox);
	tmpy = fabs(p.y-tree->monoy);
	tmpz = fabs(p.z-tree->monoz);

	r = tree->nodesize;
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
	dist = sqrt(dist2);
	diff = dist - r;
	if(diff <=fof_link) return YES;
	else return NO;
}


size_t pnew_fof_link(FoFPosition *p, FoFTStruct *tree, PosType fof_link, FoFPosition *linked, size_t nhalo){
	size_t now, ncount = 0;
	void *ptr, *optr, *nptr;
	PosType tmpx,tmpy,tmpz;
	PosType fof_link2, dist2;
	PosType Lpx, Lpy, Lpz;
	FoFPosition point;
	fof_link2 = fof_link * fof_link;
	ncount = now = 0;
	point = *p;
	/*
	point.x = p->x;
	point.y = p->y;
	point.z = p->z;
	*/

	ptr = (void*)tree;
	do{
		optr = (void *) tree;
		ptr = (void *) tree;
		while(ptr){
			switch(((TYPE*)ptr)->type) {
				case TYPE_TREE:
					if(((FoFTStruct *)ptr)->sibling == ((FoFTStruct *)ptr)->daughter){
						EraseFromTree(optr,ptr,((FoFTStruct *)ptr)->sibling);
						ptr = ((FoFTStruct *)ptr)->sibling;
					}
					else
						switch(pfof_open(point, ptr, fof_link)){
							case YES:
								optr = ptr;
								ptr = (void *)(((FoFTStruct *)ptr)->daughter);
								break;
							default:
								optr = ptr;
								ptr = (void *)(((FoFTStruct *)ptr)->sibling);
						}
					break;
				default:
					if(((FoFTPtlStruct*)ptr)->included == YES){
						nptr = ((FoFTPtlStruct*)ptr)->sibling;
						EraseFromTree(optr,ptr,nptr);
						ptr = nptr;
					}
					else {
						tmpx = fabs(point.x - ((FoFTPtlStruct*)ptr)->x);
						tmpy = fabs(point.y - ((FoFTPtlStruct*)ptr)->y);
						tmpz = fabs(point.z - ((FoFTPtlStruct*)ptr)->z);

						dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
						dist2 = sqrt(dist2);
						if(dist2<=fof_link){
							linked[ncount].x = ((FoFTPtlStruct*)ptr)->x;
							linked[ncount].y = ((FoFTPtlStruct*)ptr)->y;
							linked[ncount].z = ((FoFTPtlStruct*)ptr)->z;

							((FoFTPtlStruct*)ptr)->haloindx = nhalo;
							((FoFTPtlStruct*)ptr)->included = YES;
							ncount ++;

							nptr = ((FoFTPtlStruct *)ptr)->sibling;
							EraseFromTree(optr,ptr,nptr);
						}
						else optr = ptr;
						ptr = (void*)(((FoFTPtlStruct*)ptr)->sibling);
					}
			}
		}
		point = linked[now];
		now ++;
	} while(now <=ncount);
	return ncount;

}

