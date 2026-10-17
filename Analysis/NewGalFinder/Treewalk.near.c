/* This makes Tree structure and walks along tree structures.
 *
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "ramses.h"
#include "tree.h"
//#include "near.h"
#define MIN(A,B) ((A)<(B) ? (A):(B))
#define MAX(A,B) ((A)>(B) ? (A):(B))
void old_Make_Tree_Near(TStruct *TREE_START,TPtlStruct *ptl,int np,Box box){
	BeginEndTree beginend;
	long i;
	TStruct *NewTree;
	TREE_START->sibling = NULL;
	for(i=0;i<np;i++) ptl[i].sibling = &ptl[i+1];
	ptl[np-1].sibling = NULL;
	beginend = divide_node_Near(TREE_START,TREE_START+1,ptl,box,TREE_START);
}
BeginEndTree divide_node_Near(TStruct *TREE_START,TStruct *NewTree, 
		TPtlStruct *ptl, Box box,TStruct *ThisTree){ 
	BeginEndTree beginend;
	TStruct *p2tree, tmpnode[8];
	TStruct *NowCal;
	void *from_sibling;
	TPtlStruct *p2ptl,*tmpptr,*tmpptr2;
	Box tmpbox[8];
	int i,j,k,mnode,mx,my,mz;
	dptype x0,y0,z0,inv_halfw,halfw;
	dptype tmpx,tmpy,tmpz,tmpdist2,distmax;
	dptype ptlmass;
	int count;

	/*   */
	ThisTree->type = TYPE_TREE;
	/*
	ThisTree->sibling = NULL;
	*/
	ThisTree->daughter = NULL;
	ThisTree->L = box.width;
	ThisTree->x0 = box.x;
	ThisTree->y0 = box.y;
	ThisTree->z0 = box.z;
	ThisTree->mass = 0;
	ThisTree->Nparticle = 0;
	ThisTree->monox= ThisTree->monoy= ThisTree->monoz= 0.;
	/*
	ThisTree->quad[0] = ThisTree->quad[1] = ThisTree->quad[2] = 0.;
	ThisTree->quad[3] = ThisTree->quad[4] = ThisTree->quad[5] = 0.;
	ThisTree->mrr = 0.;
	*/

	p2ptl = ptl;
	while(p2ptl != NULL){
		ptlmass = p2ptl->mass;
		ThisTree->mass += ptlmass;
		ThisTree->Nparticle ++;
		ThisTree->monox += ptlmass * p2ptl->x;
		ThisTree->monoy += ptlmass * p2ptl->y;
		ThisTree->monoz += ptlmass * p2ptl->z;
		p2ptl = p2ptl->sibling;
	}
	ThisTree->monox /= ThisTree->mass;
	ThisTree->monoy /= ThisTree->mass;
	ThisTree->monoz /= ThisTree->mass;
	distmax = -1.E20;
	p2ptl = ptl;
	while(p2ptl != NULL){
		ptlmass = p2ptl->mass;
		tmpx = p2ptl->x - ThisTree->monox;
		tmpy = p2ptl->y - ThisTree->monoy;
		tmpz = p2ptl->z - ThisTree->monoz;
		tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		distmax = MAX(distmax,tmpdist2);
		/*
		ThisTree->mrr += ptlmass * tmpdist2;
		ThisTree->quad[0] += ptlmass * tmpx * tmpx;
		ThisTree->quad[1] += ptlmass * tmpy * tmpy;
		ThisTree->quad[2] += ptlmass * tmpz * tmpz;
		ThisTree->quad[3] += ptlmass * tmpx * tmpy;
		ThisTree->quad[4] += ptlmass * tmpx * tmpz;
		ThisTree->quad[5] += ptlmass * tmpy * tmpz;
		*/
		p2ptl = p2ptl->sibling;
	}
	ThisTree->dist2 = distmax;
	ThisTree->dist = sqrt(distmax);
	p2ptl = ptl;
	x0 = box.x;
	y0 = box.y;
	z0 = box.z;
	halfw = box.width*0.5;
	inv_halfw = 1./halfw;
	/* initialize temporary tree array */
	for(i=0;i<8;i++) {
		tmpnode[i].sibling = tmpnode[i].daughter = NULL;
		tmpnode[i].Nparticle = 0;
		tmpbox[i].width = halfw;
		tmpbox[i].x = x0+(i%2)*halfw;
		tmpbox[i].y = y0+((i%4)/2)*halfw;
		tmpbox[i].z = z0+(i/4)*halfw;
	}
	dptype x2,y2,z2;
	x2 = box.x + halfw;
	y2 = box.y + halfw;
	z2 = box.z + halfw;

	while(p2ptl != NULL){
		/*
		mx = (int)((p2ptl->x - x0)*inv_halfw);
		my = (int)((p2ptl->y - y0)*inv_halfw);
		mz = (int)((p2ptl->z - z0)*inv_halfw);
		mx = MIN(1,MAX(0,mx));
		my = MIN(1,MAX(0,my));
		mz = MIN(1,MAX(0,mz));
		*/
		if(p2ptl->x > x2)mx = 1;
		else mx = 0;
		if(p2ptl->y > y2)my = 1;
		else my = 0;
		if(p2ptl->z > z2)mz = 1;
		else mz = 0;
		mnode = mx + 2*my + 4*mz;
		tmpnode[mnode].Nparticle ++; 
		tmpptr = tmpnode[mnode].daughter;
		tmpptr2 = p2ptl->sibling;
		tmpnode[mnode].daughter = p2ptl;
		p2ptl->sibling = tmpptr;
		p2ptl = tmpptr2;
	}
	/* Making link from Mother Node */
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle > 0) break;
	}
	if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
		ThisTree->daughter = (void *)(NewTree);
	}
	else {
		ThisTree->daughter = (void *) tmpnode[i].daughter ;
	}
	count = 0;
	NowCal = NewTree;
	from_sibling = NULL;
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
			NewTree->daughter = tmpnode[i].daughter;
			if(from_sibling != NULL) 
				((GENERAL_TPtl_POINTER*)from_sibling)->sibling = NewTree;
			from_sibling = NewTree;
			NewTree++;
			count ++;
		}
		else if(tmpnode[i].Nparticle > 0 ){
			tmpptr = tmpnode[i].daughter;
			if(from_sibling != NULL) 
				((GENERAL_TPtl_POINTER*)from_sibling)->sibling = tmpptr;
			while(tmpptr != NULL){
				from_sibling = tmpptr;
				tmpptr = tmpptr->sibling;
			}
		}
	}
	((GENERAL_TPtl_POINTER*)from_sibling)->sibling = ThisTree->sibling;
	for(i=0;i<8;i++){
		if(tmpnode[i].Nparticle >= NODE_HAVE_PARTICLE){
			beginend = divide_node_Near(TREE_START,NewTree,tmpnode[i].daughter, 
					tmpbox[i],NowCal);
			NewTree = beginend.start;
			NowCal ++;
		}
	}
	beginend.start = NewTree;
	return beginend;
}
#define OFFSET 0.001
Box findsphbox(TPtlStruct *ptl, int nend){
	float xmin,ymin,zmin,xmax,ymax,zmax;
	float width;
	Box box;
	int i;
	box.x = box.y = box.z = 0.;
	width = xmax=ymax=zmax = -1.E25;
	xmin=ymin=zmin = 1.E25;
	for(i=0;i<nend;i++){
		xmin = MIN(xmin,ptl[i].x);
		ymin = MIN(ymin,ptl[i].y);
		zmin = MIN(zmin,ptl[i].z);
		xmax = MAX(xmax,ptl[i].x);
		ymax = MAX(ymax,ptl[i].y);
		zmax = MAX(zmax,ptl[i].z);
	}
	width = MAX(width,xmax-xmin);
	width = MAX(width,ymax-ymin);
	width = MAX(width,zmax-zmin);
	box.x = xmin;
	box.y = ymin;
	box.z = zmin;
	box.width = width;
	box.x = box.x -box.width*OFFSET;
	box.y = box.y -box.width*OFFSET;
	box.z = box.z -box.width*OFFSET;
	box.width = box.width + box.width*OFFSET*2.;
	return box;
}
#undef OFFSET
