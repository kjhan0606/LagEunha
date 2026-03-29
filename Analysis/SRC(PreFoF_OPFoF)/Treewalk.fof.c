/* This makes Tree structure and walks along tree structures.
 *
 * */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "fof.h"
#define IMOD(A,B) ((A) - ((A)/(B))*(B))
#define MIN(A,B) ((A)<(B) ? (A):(B))
#define MAX(A,B) ((A)>(B) ? (A):(B))
float Lx2, Ly2,Lz2;
float Lx, Ly,Lz;
/* open node in periodic boundary conditions */
enum boolean pfof_open(particle p,FoFTStruct *tree, float fof_link){
	float tmpx,tmpy,tmpz; 
	float dist2,dist,r,diff; 
	float ratio; 
	tmpx = fabs(p.x-tree->mono[0]); 
	tmpy = fabs(p.y-tree->mono[1]); 
	tmpz = fabs(p.z-tree->mono[2]); 
	if(tmpx > Lx2) tmpx = Lx-tmpx;
	if(tmpy > Ly2) tmpy = Ly-tmpy;
	if(tmpz > Lz2) tmpz = Lz-tmpz;
	r = tree->dist;
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz; 
	dist = sqrt(dist2);
	diff = dist - r;
	if(diff <= fof_link) return YES;
	else return NO;
}
enum boolean fof_open(particle p,FoFTStruct *tree, float fof_link){
	float tmpx,tmpy,tmpz; 
	float dist2,dist,r,diff; 
	float ratio; 
	tmpx = p.x-tree->mono[0]; 
	tmpy = p.y-tree->mono[1]; 
	tmpz = p.z-tree->mono[2]; 
	r = tree->dist;
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz; 
	dist = sqrt(dist2);
	diff = dist - r;
	if(diff <= fof_link) return YES;
	else return NO;
}
static int ntmp;   
static float tmpx,tmpy,tmpz,dist2;   
static float xx,yy,zz,xy,xz,yz,tmpxx; 
static float qxx1,qxx2,qxx3,qxy,tmpdist;   
static float fplmf1,fplmf2,fplmf3,fplmf4,tmptmp;   
static float fplmf;  
/*
static TStruct *pointer;   
static TPtlStruct *ppointer;   
*/
static float ptlmass;
static float idist2,sqrtdist2,isqrtdist2; 

/* 
 * *p is the position at which you want to calculate force using tree
 * theta2 is the opening angle for tree walk
 * *tree is the tree structure.
 */
/* cross? is needed to check whether the nearby particle searching crosses 
 * the boundary face. */
enum boolean crossx,crossy,crossz;

/* lx,ly, lz are the limits of the box */
int pnew_fof_link(particle *p,float fof_link,FoFTStruct *tree,
		FoFTPtlStruct *ptl,particle *linked,float lx,float ly, float lz){
	int ncount, now;
	void *ptr,*optr,*nptr;
	float tmpx,tmpy,tmpz;
	float fof_link2,dist2;
	float Lpx,Lpy,Lpz;
	particle point;
	Lx = lx; Ly = ly; Lz = lz;
	Lx2 = Lx*0.5;
	Ly2 = Ly*0.5;
	Lz2 = Lz*0.5;
	fof_link2 = fof_link*fof_link;
	ncount = now = 0;
	point.x = p->x;
	point.y = p->y;
	point.z = p->z;
	do {
		optr = (void *) tree;
		ptr = (void*) tree;
		while(ptr != NULL){
			switch(((TYPE*)ptr)->type){
				case TYPE_TREE:
					if(((FoFTStruct *)ptr)->sibling ==
							((FoFTStruct *)ptr)->daughter){
						EraseFromTree(optr,ptr,((FoFTStruct *)ptr)->sibling);
						ptr = ((FoFTStruct *)ptr)->sibling;
					}
					else
					switch(pfof_open(point,ptr,fof_link)){
						case YES:
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->daughter);
							break;
						default :
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->sibling);
					}
					break;
				default :
					if(((FoFTPtlStruct*)ptr)->included == YES){
						nptr = ((FoFTPtlStruct *)ptr)->sibling;
						EraseFromTree(optr,ptr,nptr);
						ptr = nptr;
					}
					else {
						tmpx = fabs(point.x - ((FoFTPtlStruct*)ptr)->r[0]);
						if(tmpx>Lx2) {
							tmpx = Lx-tmpx;
							crossx = YES;
						}
						else crossy = NO;

						tmpy = fabs(point.y - ((FoFTPtlStruct*)ptr)->r[1]);
						if(tmpy>Ly2) {
							tmpy = Ly-tmpy;
							crossy = YES;
						}
						else crossy = NO;

						tmpz = fabs(point.z - ((FoFTPtlStruct*)ptr)->r[2]);
						if(tmpz>Lz2) {
							tmpz = Lz-tmpz;
							crossz = YES;
						}
						else crossz = NO;

						dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
						if(dist2 <= fof_link2){
							if(crossx == YES) 
								linked[ncount].x = Lpx+((FoFTPtlStruct*)ptr)->r[0];
							else 
								linked[ncount].x = ((FoFTPtlStruct*)ptr)->r[0];
							if(crossy == YES) 
								linked[ncount].y = Lpy+((FoFTPtlStruct*)ptr)->r[1];
							else 
								linked[ncount].y = ((FoFTPtlStruct*)ptr)->r[1];
							if(crossz == YES) 
								linked[ncount].z = Lpz+((FoFTPtlStruct*)ptr)->r[2];
							else 
								linked[ncount].z = ((FoFTPtlStruct*)ptr)->r[2];

							linked[ncount].vx = ((FoFTPtlStruct*)ptr)->rv[0];
							linked[ncount].vy = ((FoFTPtlStruct*)ptr)->rv[1];
							linked[ncount].vz = ((FoFTPtlStruct*)ptr)->rv[2];
							linked[ncount].indx = ((FoFTPtlStruct*)ptr)->indx;
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
		/* for periodic boundary conditions in x,y, and z*/
		if(point.x <= fof_link) Lpx = -Lx;
		else if(point.x >= Lx-fof_link) Lpx = Lx;
		if(point.y <= fof_link) Lpy = -Ly;
		else if(point.y >= Ly-fof_link) Lpy = Ly;
		if(point.z <= fof_link) Lpz = -Lz;
		else if(point.z >= Lz-fof_link) Lpz = Lz;
		now ++;
	} while( now <= ncount);
	return (ncount);
}
int new_fof_link(particle *p,float fof_link,FoFTStruct *tree,
		FoFTPtlStruct *ptl,particle *linked){
	int ncount, now;
	void *ptr,*optr,*nptr;
	float tmpx,tmpy,tmpz;
	float fof_link2,dist2;
	particle point;
	fof_link2 = fof_link*fof_link;
	ncount = now = 0;
	point.x = p->x;
	point.y = p->y;
	point.z = p->z;
	do {
		optr = (void *) tree;
		ptr = (void*) tree;
		while(ptr != NULL){
			switch(((TYPE*)ptr)->type){
				case TYPE_TREE:
					if(((FoFTStruct *)ptr)->sibling ==
							((FoFTStruct *)ptr)->daughter){
						EraseFromTree(optr,ptr,((FoFTStruct *)ptr)->sibling);
						ptr = ((FoFTStruct *)ptr)->sibling;
					}
					else
					switch(fof_open(point,ptr,fof_link)){
						case YES:
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->daughter);
							break;
						default :
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->sibling);
					}
					break;
				default :
					if(((FoFTPtlStruct*)ptr)->included == YES){
						nptr = ((FoFTPtlStruct *)ptr)->sibling;
						EraseFromTree(optr,ptr,nptr);
						ptr = nptr;
					}
					else {
						tmpx = point.x - ((FoFTPtlStruct*)ptr)->r[0];
						tmpy = point.y - ((FoFTPtlStruct*)ptr)->r[1];
						tmpz = point.z - ((FoFTPtlStruct*)ptr)->r[2];
						dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
						if(dist2 <= fof_link2){
							linked[ncount].x = ((FoFTPtlStruct*)ptr)->r[0];
							linked[ncount].y = ((FoFTPtlStruct*)ptr)->r[1];
							linked[ncount].z = ((FoFTPtlStruct*)ptr)->r[2];
							linked[ncount].vx = ((FoFTPtlStruct*)ptr)->rv[0];
							linked[ncount].vy = ((FoFTPtlStruct*)ptr)->rv[1];
							linked[ncount].vz = ((FoFTPtlStruct*)ptr)->rv[2];
							linked[ncount].indx = ((FoFTPtlStruct*)ptr)->indx;
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
	} while( now <= ncount);
	return (ncount);
}
void FoF_Make_Tree(FoFTStruct *TREE_START,FoFTPtlStruct *ptl,int np,Box box){
	FoFBeginEndTree beginend;
	FoFTStruct *NewTree;
	FoFTPtlStruct *ptr;
	ptr = ptl;
	while(ptr != NULL){
		ptr->included = NO;
		ptr = ptr->sibling;
	}
	TREE_START->sibling = NULL;
	beginend = FoF_divide_node(TREE_START,TREE_START+1,ptl,box,TREE_START);
}
FoFBeginEndTree FoF_divide_node(FoFTStruct *TREE_START,FoFTStruct *NewTree, 
		FoFTPtlStruct *ptl, Box box,FoFTStruct *ThisTree){ 
	FoFBeginEndTree beginend;
	FoFTStruct *p2tree,tmpnode[8];
	FoFTStruct *NowCal;
	void *from_sibling;
	FoFTPtlStruct *p2ptl,*tmpptr,*tmpptr2;
	Box tmpbox[8];
	int i,j,k,mnode,mx,my,mz;
	float x0,y0,z0,inv_halfw,halfw;
	float tmpx,tmpy,tmpz,tmpdist2,distmax;
	float ptlmass;
	int count;
	ThisTree->type = TYPE_TREE;
	/*
	ThisTree->sibling = NULL;
	*/
	ThisTree->daughter = NULL;
	ThisTree->L = box.width;
	ThisTree->r0[0] = box.x;
	ThisTree->r0[1] = box.y;
	ThisTree->r0[2] = box.z;
	ThisTree->Nparticle = 0;
	ThisTree->mono[0]= ThisTree->mono[1]= ThisTree->mono[2]= 0.;

	p2ptl = ptl;
	while(p2ptl != NULL){
		ThisTree->Nparticle ++;
		ThisTree->mono[0] += p2ptl->r[0];
		ThisTree->mono[1] += p2ptl->r[1];
		ThisTree->mono[2] += p2ptl->r[2];
		p2ptl = p2ptl->sibling;
	}
	ThisTree->mono[0] /= ThisTree->Nparticle;
	ThisTree->mono[1] /= ThisTree->Nparticle;
	ThisTree->mono[2] /= ThisTree->Nparticle;
	distmax = -1.E20;
	p2ptl = ptl;
	while(p2ptl != NULL){
		tmpx = p2ptl->r[0] - ThisTree->mono[0];
		tmpy = p2ptl->r[1] - ThisTree->mono[1];
		tmpz = p2ptl->r[2] - ThisTree->mono[2];
		tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		distmax = MAX(distmax,tmpdist2);
		p2ptl = p2ptl->sibling;
	}
	ThisTree->dist2 = distmax;
	ThisTree->dist = sqrt(distmax);
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
		tmpbox[i].x = x0+IMOD(i,2)*halfw;
		tmpbox[i].y = y0+(IMOD(i,4)/2)*halfw;
		tmpbox[i].z = z0+(i/4)*halfw;
	}
	p2ptl = ptl;
	while(p2ptl != NULL){
		mx = (int)((p2ptl->r[0] - x0)*inv_halfw);
		my = (int)((p2ptl->r[1] - y0)*inv_halfw);
		mz = (int)((p2ptl->r[2] - z0)*inv_halfw);
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
			beginend = FoF_divide_node(TREE_START,NewTree,
					tmpnode[i].daughter, tmpbox[i],NowCal);
			NewTree = beginend.start;
			NowCal ++;
		}
	}
	beginend.start = NewTree;
	return beginend;
}
