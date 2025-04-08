#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<string.h>
#include<float.h>
#include<unistd.h>
#ifdef _OPENMP
#include<omp.h>
#endif

#include "eunha.h"
#include "ost.h"
#include "Model4TreeForce.h"


int nullfct00() { return 0;}
int nullfct10() { return 1;}
#ifndef _OPENMP
#define omp_get_thread_num() nullfct00()
#define omp_get_num_threads() nullfct10()
#endif



#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define YES 1
#define NO 0


#define SubCellDicision(a,b) ((a)>(b)? 1:0)
#define ShallIDivideNode(ThisNode,nparticles) ((ThisNode->nodesize > 0.5*MINCELLWIDTH ? YES:NO) && (nparticles>=MIN_CELL_PARTICLE_NUM ? YES:NO))



TStruct *divide_node(TStruct *ThisNode, TStruct *FreeNode, float thetasq,int recursiveflag)
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
	float q0,q1,q2,q3,q4,q5;
	q0=q1=q2=q3=q4=q5=0;
	for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
		ptlmass = p2ptl->mass;
		tmpx = p2ptl->x - ThisNode->monox;
		tmpy = p2ptl->y - ThisNode->monoy;
		tmpz = p2ptl->z - ThisNode->monoz;
		tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		distmax = MAX(distmax, tmpdist2);
		q0 += ptlmass * tmpx*tmpx;
		q1 += ptlmass * tmpy*tmpy;
		q2 += ptlmass * tmpz*tmpz;
		q3 += ptlmass * tmpx*tmpy;
		q4 += ptlmass * tmpx*tmpz;
		q5 += ptlmass * tmpy*tmpx;
	}
	ThisNode->nodesize = sqrt(distmax);
	ThisNode->quad[0] = q0; ThisNode->quad[1] = q1; ThisNode->quad[2] = q2; ThisNode->quad[3] = q3; ThisNode->quad[4] = q4; ThisNode->quad[5] = q5;
	ThisNode->trQ = q0+q1+q2;
	ThisNode->dist_over_thetasq = distmax/thetasq;

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
	FirstDaughter = NowDaughter =  FreeNode;
	/* Link to the first daughter or first particle */
	if( ShallIDivideNode(ThisNode, tmpnode[i].Nparticle) )
		ThisNode->daughter = (void*)FirstDaughter;
	else 
		ThisNode->daughter = (void*) tmpnode[i].daughter;

	LeftSibling = NULL;

	for(i=0;i<8;i++) {
		if( ShallIDivideNode(ThisNode, tmpnode[i].Nparticle) ){
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
	FreeNode = NowDaughter;


	if(recursiveflag == RECURSIVE){
		TStruct *NextJobNode = FirstDaughter;
		for(i=0;i<8;i++){
			if(ShallIDivideNode(ThisNode, tmpnode[i].Nparticle) ){
				FreeNode = divide_node(NextJobNode,FreeNode, thetasq,recursiveflag);
				NextJobNode ++;
			}
		}
	}
	return FreeNode;
}

void Make_Tree(TStruct *TREE_START, size_t nnode, TPtlStruct *ptl, size_t np, float thetasq,int recursiveflag){
	size_t i;
	TStruct *FreeNode = TREE_START+1;
	TREE_START->sibling = NULL;
	for(i=0;i<np;i++) {
		ptl[i].sibling = ptl + (i +1);
		ptl[i].type = TYPE_PTL;
	}
	ptl[np-1].sibling = NULL;
	TREE_START->daughter = &(ptl[0]);
	TREE_START->Nparticle = np;
	if(recursiveflag == RECURSIVE) FreeNode = divide_node(TREE_START, FreeNode, thetasq,recursiveflag);
	else if(recursiveflag == SERIALIZED){
		TStruct *work;
		for(work=TREE_START;FreeNode-work >0; work++){
			FreeNode = divide_node(work, FreeNode, thetasq,SERIALIZED);
		}
	}
	else if(recursiveflag == PTHREAD) 
	{
		TStruct *work = TREE_START;
		do{
			FreeNode = divide_node(work, FreeNode, thetasq,PTHREAD);
			work ++;
		}
		while( work < FreeNode && (FreeNode-work) <= 50);
		size_t twork = (FreeNode-work);

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
			TStruct *threadFreeNode = FreeNode + pid * ( (nnode - (FreeNode-TREE_START))/ npid);
			for(j=mys;j<myf;j++)
			{
				threadFreeNode = divide_node(work+j,threadFreeNode,thetasq,RECURSIVE);
			}
		}
	}
}




inline enum where InsideOpen(Position p, TStruct *tree, float maxdist) {
	PosType tmpx = p.x-tree->monox;
	PosType tmpy = p.y-tree->monoy;
	PosType tmpz = p.z-tree->monoz;
	PosType dist = sqrt(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);
	if(dist+tree->nodesize <= maxdist) return IN;
	else if(dist-tree->nodesize > maxdist) return OUT;
	else return CROSS;

}


ptrdiff_t insidecount(Position point, TStruct *tree, float maxdist){
	ptrdiff_t count = 0;
	PosType tmpx,tmpy,tmpz;
	PosType maxdist2 = maxdist*maxdist;

	void *ptr = (void*)tree;
	while(ptr) {
		switch(((TYPE*)ptr)->type) {
			case TYPE_TREE:
				switch(InsideOpen(point, ptr, maxdist)){
					case OUT:
						ptr = (void*)(((TStruct*)ptr)->sibling);
						break;
					case IN:
						count += ((TStruct*) ptr)->Nparticle;
						ptr = (void*)(((TStruct*)ptr)->sibling);
						break;

					default:
						ptr = (void*)(((TStruct *)ptr)->daughter);
				}
				break;
			default:
				tmpx = point.x - ((TPtlStruct*)ptr)->x;
				tmpy = point.y - ((TPtlStruct*)ptr)->y;
				tmpz = point.z - ((TPtlStruct*)ptr)->z;
				PosType dist2 = (tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);
				if(dist2 < maxdist2) count ++;
				ptr= (void*)(((TPtlStruct*)ptr)->sibling);
		}
	}
	return count;
}

#define AddNode(p,ptr) do {\
	TStruct *pointer= (TStruct*)ptr;\
	float tmpx,tmpy,tmpz,dist2;\
	tmpx = (p->x)-(pointer->monox);\
	tmpy = (p->y)-(pointer->monoy);\
	tmpz = (p->z)-(pointer->monoz);\
	dist2 = (tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);\
	if(dist2 <=rspheresq){\
		xnode[*nnode] = tmpx;\
		ynode[*nnode] = tmpy;\
		znode[*nnode] = tmpz;\
		xxnode[*nnode] = pointer->quad[0];\
		yynode[*nnode] = pointer->quad[1];\
		zznode[*nnode] = pointer->quad[2];\
		xynode[*nnode] = pointer->quad[3];\
		xznode[*nnode] = pointer->quad[4];\
		yznode[*nnode] = pointer->quad[5];\
		trnode[*nnode] = pointer->trQ;\
		mnode[*nnode] = pointer->mass;\
		(*nnode) ++;\
	}\
} while(0)

#define AddPar(ptr,npar,mpar,xpar,ypar,zpar) do{\
	TPtlStruct *pointer= (TPtlStruct*)ptr;\
	float tmpx,tmpy,tmpz,dist2;\
	tmpx = p->x - pointer->x;\
	tmpy = p->y - pointer->y;\
	tmpz = p->z - pointer->z;\
	dist2 = (tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);\
	if(dist2 <= rspheresq){\
		xpar[*npar] = tmpx;\
		ypar[*npar] = tmpy;\
		zpar[*npar] = tmpz;\
		mpar[*npar] = pointer->mass;\
		(*npar) ++;\
	}\
} while(0)


inline int nodeopen(Position *p, TStruct *node){
	NodePosType tmpx,tmpy,tmpz;
	NodePosType dist2;
	NodePosType ratio;
	tmpx = p->x - node->monox;
	tmpy = p->y - node->monoy;
	tmpz = p->z - node->monoz;
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
	if(dist2 <node->dist_over_thetasq) return YES;
	else return NO;
}

void inline treeforcesearch(p, tree, rspheresq, npar, xpar,ypar,zpar, mpar,
		nnode,xnode,ynode,znode,xxnode,yynode,zznode,xynode,xznode,yznode,trnode,mnode)
Position *p;
TStruct *tree;
size_t *npar, *nnode;
float rspheresq;
float *xpar,*ypar,*zpar,*mpar;
float *xnode,*ynode,*znode;
float *xxnode,*yynode, *zznode, *xynode,*xznode,*yznode,*trnode,*mnode;
{
	void *ptr, *optr;
	size_t ntmp;

	ptr = (void *) tree;
	*nnode = *npar = 0;
	while(ptr){
		switch (((TYPE*)ptr)->type) {
			case TYPE_TREE:
				switch (nodeopen(p,ptr)){
					case YES:
						ptr = (void *)(((TStruct *)ptr)->daughter);
						break;
					default:
						AddNode(p,ptr);
						ptr = (void *)(((TStruct *)ptr)->sibling);
				}
				break;
			default:
				AddPar(ptr,npar,mpar,xpar,ypar,zpar);
				ptr = (void *)(((TPtlStruct *)ptr)->sibling);
		}
	}
	return;
}
inline void ParticleGravForce(SimParameters *simpar, Position *pos, size_t np, float *xpar, 
		float *ypar, float *zpar, float *mpar){
	size_t i;
	float ax,ay,az;
	ax = ay = az = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: ax,ay,az)
#endif
	for(i=0;i<np;i++){
		float dist = sqrt(xpar[i]*xpar[i] + ypar[i]*ypar[i] + zpar[i]*zpar[i]);
		float fplmf = GetPointGrav(dist) * mpar[i];
		ax += xpar[i]*fplmf;
		ay += ypar[i]*fplmf;
		az += zpar[i]*fplmf;
	}
#ifndef GOTPM
	(pos->axyz[0]) += ax;
	(pos->axyz[1]) += ay;
	(pos->axyz[2]) += az;
#else
	float vfact2 = Evol_FACT2(simpar);
	(pos->axyz[0]) += vfact2*ax;
	(pos->axyz[1]) += vfact2*ay;
	(pos->axyz[2]) += vfact2*az;
		
#endif
}
inline void CellGravForce(SimParameters *simpar, Position *pos,
		size_t nnode, float *xnode, float *ynode, float *znode, 
		float *xxnode, float *yynode, float *zznode, float *xynode,
		float *xznode, float *yznode, float *trnode, float *mnode){
	size_t i;
	float ax,ay,az;
	ax = ay = az = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: ax,ay,az)
#endif
	for(i=0;i<nnode;i++){
		float dist = sqrt(xnode[i]*xnode[i] + ynode[i]*ynode[i] + znode[i]*znode[i]);
		float qxx1 = xxnode[i]*xnode[i] + xynode[i]*ynode[i] + xznode[i]*znode[i];
		float qxx2 = xynode[i]*xnode[i] + yynode[i]*ynode[i] + yznode[i]*znode[i];
		float qxx3 = xznode[i]*xnode[i] + yznode[i]*ynode[i] + zznode[i]*znode[i];
		float qxy = xxnode[i]*xnode[i]*xnode[i] + yynode[i]*ynode[i]*ynode[i] + 
			zznode[i]*znode[i]*znode[i] + 
			2*( xynode[i]*xnode[i]*ynode[i] + xznode[i]*xnode[i]*znode[i] + 
					yznode[i]*ynode[i]*znode[i] );
		float fplmf1, fplmf2, fplmf3;
		GetCellGrav(dist, &fplmf1, &fplmf2, &fplmf3);
		float tmptmp = mnode[i]*fplmf1 + trnode[i]*fplmf2 + qxy*fplmf3;
		float twofplmf2 = 2*fplmf2;
		ax += xnode[i]*tmptmp + qxx1 * twofplmf2;
		ay += ynode[i]*tmptmp + qxx2 * twofplmf2;
		az += znode[i]*tmptmp + qxx3 * twofplmf2;
	}

#ifndef GOTPM
	(pos->axyz[0]) += ax;
	(pos->axyz[1]) += ay;
	(pos->axyz[2]) += az;
#else
	float vfact2 = Evol_FACT2(simpar);
	(pos->axyz[0]) += vfact2*ax;
	(pos->axyz[1]) += vfact2*ay;
	(pos->axyz[2]) += vfact2*az;
#endif
}

void treeforce(SimParameters *simpar,
		Position *pos, TStruct *TREECELL, size_t linesize,float rspheresq, float *work){
	float *xpar,*ypar,*zpar,*mpar;
	float *xnode,*ynode,*znode;
	float *xxnode,*yynode, *zznode, *xynode,*xznode,*yznode,*trnode,*mnode;
	xpar = work;
	ypar = work+1*linesize;
	zpar = work+2*linesize;
	mpar = work+3*linesize;
	xnode = work+4*linesize;
	ynode = work+5*linesize;
	znode = work+6*linesize;
	xxnode = work+7*linesize;
	yynode = work+8*linesize;
	zznode = work+9*linesize;
	xynode = work+10*linesize;
	xznode = work+11*linesize;
	yznode = work+12*linesize;
	trnode = work+13*linesize;
	mnode = work+14*linesize;
	size_t np, nnode;
	treeforcesearch(pos, TREECELL, rspheresq, &np, xpar,ypar,zpar,mpar, &nnode, xnode,ynode,znode,
			xxnode, yynode, zznode, xynode, xznode, yznode, trnode, mnode);
	ParticleGravForce(simpar, pos,np,xpar,ypar,zpar,mpar);
	CellGravForce(simpar, pos, nnode, xnode,ynode,znode, xxnode, yynode, zznode, xynode, xznode,
			yznode, trnode, mnode);
}
