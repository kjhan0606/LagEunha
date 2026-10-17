#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include "ramses.h"
#include "tree.h"
#ifdef _OPENMP
#include<omp.h>
#endif

/*
#include "fof.h"
*/

#include "mathsub.h"
#include "force_spline.h"
extern float diff[NSPLINE][3],slope[NSPLINE][3],ran2nran;


#define nodesize dist
int nullfct0() { return 0;}
int nullfct1() { return 1;}
#ifndef _OPENMP
#define Omp_get_thread_num() nullfct0()
#define Omp_get_num_threads() nullfct1()
#else
#define Omp_get_thread_num() omp_get_thread_num()
#define Omp_get_num_threads() omp_get_num_threads()
#endif

#define SubCellDicision(a,b) ((a)>(b)? 1:0)
#define DivideNode(ThisNode,nparticles) ((ThisNode->nodesize > 0.5*MINCELLWIDTH ? yes:no) & (nparticles>=MIN_CELL_PARTICLE_NUM ? yes:no))

FoFTStruct *calFoFThreadFreeNodeStart(size_t navail, size_t twork, size_t mys,size_t myf, 
		FoFTStruct *work, FoFTStruct *nextFreeNode){
	FoFTStruct *threadFreeNode, *tmp;
	size_t i, tnp;
	tnp = 0;
	tmp = work;
	for(i=0;i<twork;i++) {
		tnp += tmp->Nparticle;
		tmp ++;
	}
	float ratio = (float)navail/(float)tnp;
	size_t rnp = 0;
	tmp = work;
	for(i=0;i<mys;i++){
		rnp += tmp->Nparticle;
		tmp++;
	}
	threadFreeNode = nextFreeNode + (size_t)(rnp*ratio);
	return threadFreeNode;
}

enum boolean fof_open(particle p,FoFTStruct *tree, float fof_link){
    dptype tmpx,tmpy,tmpz;
    dptype dist2,dist,r,diffr;
    dptype ratio, link02;
    tmpx = p.x-tree->monox;
    tmpy = p.y-tree->monoy;
    tmpz = p.z-tree->monoz;
    r = tree->dist;
    dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
    dist = sqrtf(dist2);
    diffr = dist - r;
	link02 = 0.5*(p.link02 + tree->maxlink02);
    if(diffr <= link02) return YES;
    else return NO;
}

enum boolean treeopen(particle *p,TStruct *tree, float theta2){
	dptype tmpx,tmpy,tmpz; 
	dptype dist2; 
	dptype ratio; 
	tmpx = p->x-tree->monox;
	tmpy = p->y-tree->monoy;
	tmpz = p->z-tree->monoz;
	dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz; 
	ratio = tree->dist2/theta2; 
	if(dist2 < ratio) return YES;
	else return NO;
}
#define CELLFORCE(p,ptr,force,ran2nran) \
{\
	int ntmp;   \
	dptype tmpx,tmpy,tmpz,dist2,dist;   \
	dptype xx,yy,zz,xy,xz,yz,tmpxx; \
	dptype qxx1,qxx2,qxx3,qxy,tmpdist;   \
	dptype fplmf1,fplmf2,fplmf3,fplmf4,tmptmp;   \
	TStruct *pointer;   \
	pointer = (TStruct *) ptr;   \
	tmpx = pointer->monox - p->x;   \
	tmpy = pointer->monoy - p->y; \
	tmpz = pointer->monoz - p->z;   \
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;   \
	xx = pointer->quad[0];   \
	yy = pointer->quad[1];   \
	zz = pointer->quad[2];   \
	xy = pointer->quad[3];   \
	xz = pointer->quad[4];   \
	yz = pointer->quad[5];   \
	tmpxx = xx + yy + zz;   \
	qxx1 = xx*tmpx+xy*tmpy+xz*tmpz;   \
	qxx2 = xy*tmpx+yy*tmpy+yz*tmpz;   \
	qxx3 = xz*tmpx+yz*tmpy+zz*tmpz;   \
	qxy = xx*tmpx*tmpx+yy*tmpy*tmpy+zz*tmpz*tmpz+   \
		2.*(xy*tmpx*tmpy+xz*tmpx*tmpz+yz*tmpy*tmpz); \
	dist = sqrtf(dist2);\
	ntmp = dist/ran2nran;   \
	tmpdist = (dist-ntmp*ran2nran);   \
	fplmf1 = diff[ntmp][0] + slope[ntmp][0]*tmpdist;   \
	fplmf2 = diff[ntmp][1] + slope[ntmp][1]*tmpdist;   \
	fplmf3 = diff[ntmp][2] + slope[ntmp][2]*tmpdist;   \
	tmptmp = pointer->mass * fplmf1 + tmpxx*fplmf2 +   \
		qxy*fplmf3;   \
	force->x += tmpx*tmptmp+2.*qxx1*fplmf2; \
	force->y += tmpy*tmptmp+2.*qxx2*fplmf2; \
	force->z += tmpz*tmptmp+2.*qxx3*fplmf2; \
}

#define PARTICLEFORCE(p,ptr,force,ran2nran) \
{ \
	int ntmp;  \
	dptype tmpx,tmpy,tmpz,dist2,dist;  \
	dptype fplmf;  \
	dptype ptlmass; \
	TPtlStruct *pointer; \
	pointer = (TPtlStruct *) ptr; \
	ptlmass = pointer->mass; \
	tmpx = pointer->x - p->x;  \
	tmpy = pointer->y - p->y;  \
	tmpz = pointer->z - p->z;  \
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;  \
	dist = sqrtf(dist2);\
	ntmp = dist/ran2nran;  \
	fplmf = diff[ntmp][0] + slope[ntmp][0] * (dist-ntmp*ran2nran);  \
	fplmf *= ptlmass; \
	force->x +=  tmpx * fplmf;  \
	force->y +=  tmpy * fplmf;  \
	force->z +=  tmpz * fplmf;  \
}

#define CELLPLUMPOTENT(p,ptr,potent,epsilon2) \
{\
	int ntmp;   \
	dptype tmpx,tmpy,tmpz,dist2,dist;   \
	dptype xx,yy,zz,xy,xz,yz,tmpxx; \
	dptype qxx1,qxx2,qxx3,qxy,tmpdist;   \
	dptype fplmf1,fplmf2,fplmf3,fplmf4,tmptmp;   \
	dptype idist2,sqrtdist2,isqrtdist2; \
	TStruct *pointer;   \
	pointer = (TStruct *) ptr;   \
	tmpx = p->x - pointer->monox;   \
	tmpy = p->y - pointer->monoy;   \
	tmpz = p->z - pointer->monoz;   \
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;   \
	xx = pointer->quad[0];   \
	yy = pointer->quad[1];   \
	zz = pointer->quad[2];   \
	xy = pointer->quad[3];   \
	xz = pointer->quad[4];   \
	yz = pointer->quad[5];   \
	tmpxx = xx + yy + zz;   \
	qxy = xx*tmpx*tmpx+yy*tmpy*tmpy+zz*tmpz*tmpz+   \
		2.*(xy*tmpx*tmpy+xz*tmpx*tmpz+yz*tmpy*tmpz); \
	dist = sqrtf(dist2);\
	ntmp = dist/ran2nran;   \
	ntmp = MIN(ntmp,NSPLINE-1);\
	tmpdist = (dist-ntmp*ran2nran);   \
	fplmf1 = diff[ntmp][0] + slope[ntmp][0]*tmpdist;   \
	fplmf2 = diff[ntmp][1] + slope[ntmp][1]*tmpdist;   \
	fplmf3 = diff[ntmp][2] + slope[ntmp][2]*tmpdist;   \
	potent += pointer->mass * fplmf1 + tmpxx*fplmf2 + \
		qxy*fplmf3; \
}

#define PARTICLEPLUMPOTENT(p,ptr,potent,epsilon2) \
{\
	float tmpx,tmpy,tmpz;   \
	TPtlStruct *pointer; \
	float fplmf,ptlmass,dist2,dist; \
	int ntmp; \
	pointer = (TPtlStruct *) ptr; \
	ptlmass = pointer->mass; \
	tmpx = p->x -  pointer->x; \
	tmpy = p->y -  pointer->y; \
	tmpz = p->z -  pointer->z; \
	dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;  \
	if(dist2 >0) fplmf = -ptlmass/sqrt(dist2+epsilon2);\
	else fplmf = 0;\
	potent += fplmf; \
}


/* 
 * *p is the position at which you want to calculate force using tree
 * theta2 is the opening angle for tree walk
 * *tree is the tree structure.
 */

void treeforce(particle *p,float theta2,TStruct *tree, TPtlStruct *ptl,
		pforce *force){
	void *ptr,*optr;
	ptr = (void *)tree;
	force->x = force->y = force->z = 0.;
	while(ptr != NULL){
		switch(((TYPE*)ptr)->type ) {
			case TYPE_TREE:
				switch(treeopen(p,ptr,theta2)){
					case YES: 
						ptr = (void *)(((TStruct*)ptr)->daughter);
						break;
					default :
						CELLFORCE(p,ptr,force,ran2nran);
						ptr = (void *)(((TStruct*)ptr)->sibling);
				}
				break;
			default: 
				PARTICLEFORCE(p,ptr,force,ran2nran);
				ptr = (void*)(((TPtlStruct*)ptr)->sibling);
		}
	}
	return;
}
float treeplumpotential(particle *p,float theta2,TStruct *tree,
		TPtlStruct *ptl){
	void *ptr;
	float potent=0.;
	float epsilon2;
	epsilon2 = EPSILON * EPSILON;
	ptr = (void *) tree;
	while(ptr != NULL){
		switch(((TYPE*)ptr)->type ) {
			case TYPE_TREE:
				switch(treeopen(p,ptr,theta2)){
					case YES: 
						ptr = (void *)(((TStruct*)ptr)->daughter);
						break;
					default :
						CELLPLUMPOTENT(p,ptr,potent,epsilon2);
						ptr = (void *)(((TStruct*)ptr)->sibling);
				}
				break;
			default: 
				PARTICLEPLUMPOTENT(p,ptr,potent,epsilon2);
				ptr = (void*)(((TPtlStruct*)ptr)->sibling);
		}
	}
	return potent;
}
typedef struct pthreadStruct{
    POSTYPE fof_link;
    FoFTStruct *tree;
    FoFTPtlStruct *ptl;
	particle point;
    int nptl, id, nthreads, haloid;
} pthreadStruct;

void pthreadFoFGroup(void *args){
    pthreadStruct *p = (pthreadStruct *)args;
    POSTYPE fof_link = p->fof_link;
    FoFTStruct *tree = p->tree;
    FoFTPtlStruct *ptl = p->ptl;
    int nptl = p->nptl;
    int pid = p->id;
    int nthreads = p->nthreads;
    int haloid = p->haloid;
    particle point = p->point;
    {
        void *ptr = (void *) tree;
        while(ptr){
			FoFTPtlStruct *tmp;
			int ioff;
            switch(((TYPE*)ptr)->type){
                case TYPE_TREE:
                    switch(fof_open(point,ptr,fof_link)){
                        case NO:
                            ptr = (void *)(((FoFTStruct*)ptr)->sibling);
                            break;
                        default:
                            ptr = (void *)(((FoFTStruct*)ptr)->daughter);
                    }
                    break;
                default:
                    tmp = (FoFTPtlStruct*)ptr;
                    ioff = tmp - ptl;
                    if(ioff%nthreads == pid && tmp->included == NO) {
                        POSTYPE tmpx = point.x - tmp->x;
                        POSTYPE tmpy = point.y - tmp->y;
                        POSTYPE tmpz = point.z - tmp->z;
                        POSTYPE dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
                        dist2 = sqrt(dist2);
                        if(dist2 <= 0.5*(point.link02+tmp->link02)){
                            tmp->included =NEW;
                            tmp->haloindx =haloid;
                        }
                    }
                    ptr = (void*)(tmp->sibling);
            }

        }
    }
}

int new_fof_link(particle *p,POSTYPE fof_link,FoFTStruct *tree,
        FoFTPtlStruct *ptl,particle *linked){
    int ncount, now;
    void *ptr,*optr,*nptr;
    POSTYPE fof_link2;
    particle point;
    ncount = now = 0;
	point = *p;
    do {
        ptr = (void*) tree;
        while(ptr != NULL){
            switch(((TYPE*)ptr)->type){
                case TYPE_TREE:
                    switch(fof_open(point,ptr,fof_link)){
                        case YES:
                            ptr = (void *)(((FoFTStruct*)ptr)->daughter);
                            break;
                        default :
                            ptr = (void *)(((FoFTStruct*)ptr)->sibling);
                    }
                    break;
                default :
                    if(((FoFTPtlStruct*)ptr)->included == NO) {
                        POSTYPE tmpx = point.x - ((FoFTPtlStruct*)ptr)->x;
                        POSTYPE tmpy = point.y - ((FoFTPtlStruct*)ptr)->y;
                        POSTYPE tmpz = point.z - ((FoFTPtlStruct*)ptr)->z;
                        POSTYPE dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
						POSTYPE dist = sqrtf(dist2);
                        if(dist <= 0.5*(point.link02+((FoFTPtlStruct*)ptr)->link02))
						{
                            linked[ncount].x = ((FoFTPtlStruct*)ptr)->x;
                            linked[ncount].y = ((FoFTPtlStruct*)ptr)->y;
                            linked[ncount].z = ((FoFTPtlStruct*)ptr)->z;
                            linked[ncount].link02 = ((FoFTPtlStruct*)ptr)->link02;
                            ((FoFTPtlStruct*)ptr)->included = YES;
                            ncount ++;
                        }
                    }
                    ptr = (void*)(((FoFTPtlStruct*)ptr)->sibling);
            }
        }
        point = linked[now];
        now ++;
    } while( now <= ncount);
    return ncount;
}
int ScoopUpParticles(void *optr, particle *linked, int oncount){
    int ncount=oncount;
    void *ptr = optr;
    void *terminal = (void*)(((FoFTPtlStruct*)ptr)->sibling);
    while(ptr != terminal){
		FoFTPtlStruct *tmp;
        switch(((TYPE*)ptr)->type){
            case TYPE_TREE:
                ptr = ((FoFTStruct*)ptr)->daughter;
                break;
            default :
                tmp = (FoFTPtlStruct*)ptr;
                if(((FoFTPtlStruct*)ptr)->included == NO){
                    linked[ncount].x = tmp->x;
                    linked[ncount].y = tmp->y;
                    linked[ncount].z = tmp->z;
                    linked[ncount].link02 = tmp->link02;
                    tmp->included = YES;
                    ncount ++;
                }
                ptr = (void*)(tmp->sibling);
        }
    }
    ((FoFTStruct *)ptr)->sibling = ((FoFTStruct *)ptr)->daughter;
    return ncount;
}
void ScoopUpParticles2New(void *optr){
    void *ptr = optr;
    void *terminal = (void*)(((FoFTPtlStruct*)ptr)->sibling);
    while(ptr != terminal){
        FoFTPtlStruct *tmp;
        switch(((TYPE*)ptr)->type){
            case TYPE_TREE:
                ptr = ((FoFTStruct*)ptr)->daughter;
                break;
            default :
                tmp = (FoFTPtlStruct*)ptr;
                if(((FoFTPtlStruct*)ptr)->included == NO){
                    tmp->included = NEW;
                }
                ptr = (void*)(tmp->sibling);
        }
    }
    ((FoFTStruct *)ptr)->sibling = ((FoFTStruct *)ptr)->daughter;
    return ;
}

void destroy_omp_fof_link(particle *p,POSTYPE fof_link,FoFTStruct *tree, FoFTPtlStruct *ptl){
    void *ptr,*optr,*nptr;
	POSTYPE fof_link2;
    particle point;
    point = *p;
    {
        optr = ptr = (void*) tree;
        while(ptr != NULL){
            switch(((TYPE*)ptr)->type){
                case TYPE_TREE:
					if(((FoFTStruct *)ptr)->sibling ==
							((FoFTStruct *)ptr)->daughter){
						EraseFromTree(optr,ptr,((FoFTStruct *)ptr)->sibling);
						ptr =  ((FoFTStruct *)ptr)->sibling;
					}
					else
					switch(fof_open(point,ptr,fof_link)){ 
						case ENCLOSE:
							ScoopUpParticles2New(ptr);
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->sibling); 
							break;
						case NO : 
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->sibling); 
							break;
						default :
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->daughter); 
					}
                    break;
                default :
                    if(((FoFTPtlStruct*)ptr)->included == YES ||
							((FoFTPtlStruct*)ptr)->included == NEW ) {
						nptr = ((FoFTPtlStruct *)ptr)->sibling;
						EraseFromTree(optr,ptr,nptr);
						ptr = nptr;
					}
					else
					{
						FoFTPtlStruct *tmp = (FoFTPtlStruct*)ptr;
                        POSTYPE tmpx = point.x - tmp->x;
                        POSTYPE tmpy = point.y - tmp->y;
                        POSTYPE tmpz = point.z - tmp->z;
                        POSTYPE dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
                        POSTYPE dist = sqrtf(dist2);
                        if(dist <= 0.5*(point.link02+tmp->link02))
                        {
                            tmp->included = NEW;

							nptr = tmp->sibling;
							EraseFromTree(optr,ptr,nptr);
                        }
						else optr = ptr;
                    	ptr = (void*)(tmp->sibling);
                    }
            }
        }
    }
}

int destroy_new_fof_link(particle *p,POSTYPE fof_link,FoFTStruct *tree,
        FoFTPtlStruct *ptl,particle *linked){
    int ncount, now;
    void *ptr,*optr,*nptr;
	POSTYPE fof_link2;
    particle point;
    ncount = now = 0;
    point = *p;
    do {
        optr = ptr = (void*) tree;
//        optr = (void*) tree;
        while(ptr != NULL){
            switch(((TYPE*)ptr)->type){
                case TYPE_TREE:
					if(((FoFTStruct *)ptr)->sibling ==
							((FoFTStruct *)ptr)->daughter){
						EraseFromTree(optr,ptr,((FoFTStruct *)ptr)->sibling);
						ptr =  ((FoFTStruct *)ptr)->sibling;
					}
					else
					switch(fof_open(point,ptr,fof_link)){ 
						case ENCLOSE:
							ncount = ScoopUpParticles(ptr,linked,ncount);
							optr = ptr;
							ptr = (void *)(((FoFTStruct*)ptr)->sibling); 
							break;
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
                    if(((FoFTPtlStruct*)ptr)->included == YES) {
						nptr = ((FoFTPtlStruct *)ptr)->sibling;
						EraseFromTree(optr,ptr,nptr);
						ptr = nptr;
					}
					else
					{
                        POSTYPE tmpx = point.x - ((FoFTPtlStruct*)ptr)->x;
                        POSTYPE tmpy = point.y - ((FoFTPtlStruct*)ptr)->y;
                        POSTYPE tmpz = point.z - ((FoFTPtlStruct*)ptr)->z;
                        POSTYPE dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
                        POSTYPE dist = sqrtf(dist2);
                        if(dist <= 0.5*(point.link02+((FoFTPtlStruct*)ptr)->link02))
                        {
                            linked[ncount].x = ((FoFTPtlStruct*)ptr)->x;
                            linked[ncount].y = ((FoFTPtlStruct*)ptr)->y;
                            linked[ncount].z = ((FoFTPtlStruct*)ptr)->z;
                            linked[ncount].link02 = ((FoFTPtlStruct*)ptr)->link02;
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
    return ncount;
}


TStruct *divide_node(TStruct *motherNode, TStruct *nextFreeNode, float thetasq,int recursiveflag)
{
    TStruct *p2tree, tmpnode[8];
    void *lastOffspring;
    TPtlStruct *p2ptl, *tmpptr, *tmpptr2, *nodeparticles;
    size_t i, j, k, mnode, mx,my,mz;
    float tmpx,tmpy,tmpz,tmpdist2, distmax;
    float ptlmass;

    nodeparticles = motherNode->daughter;
    motherNode->type = TYPE_TREE;
    motherNode->mass = 0;
    motherNode->monox = motherNode->monoy = motherNode->monoz = 0;

    for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
        ptlmass = p2ptl->mass;
        motherNode->mass += ptlmass;
        motherNode->monox += ptlmass * p2ptl->x;
        motherNode->monoy += ptlmass * p2ptl->y;
        motherNode->monoz += ptlmass * p2ptl->z;
    }
    motherNode->monox /= motherNode->mass;
    motherNode->monoy /= motherNode->mass;
    motherNode->monoz /= motherNode->mass;
    distmax = -1.E20;
    float q0,q1,q2,q3,q4,q5;
    q0=q1=q2=q3=q4=q5=0;
    for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
        ptlmass = p2ptl->mass;
        tmpx = p2ptl->x - motherNode->monox;
        tmpy = p2ptl->y - motherNode->monoy;
        tmpz = p2ptl->z - motherNode->monoz;
        tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
        distmax = MAX(distmax, tmpdist2);
        q0 += ptlmass * tmpx*tmpx;
        q1 += ptlmass * tmpy*tmpy;
        q2 += ptlmass * tmpz*tmpz;
        q3 += ptlmass * tmpx*tmpy;
        q4 += ptlmass * tmpx*tmpz;
        q5 += ptlmass * tmpy*tmpx;
    }
    motherNode->dist2 = distmax; 
    motherNode->nodesize = sqrt(distmax); 
	motherNode->quad[0] = q0; 
	motherNode->quad[1] = q1; 
	motherNode->quad[2] = q2; 
	motherNode->quad[3] = q3; 
	motherNode->quad[4] = q4; 
	motherNode->quad[5] = q5;
    motherNode->trQ = q0+q1+q2;
    motherNode->dist_over_thetasq = distmax/thetasq;

    for(i=0;i<8;i++) {
        tmpnode[i].sibling = tmpnode[i].daughter = NULL;
        tmpnode[i].Nparticle = 0;
    }
    for(p2ptl=nodeparticles;p2ptl;) {
        mx = SubCellDicision(p2ptl->x,motherNode->monox);
        my = SubCellDicision(p2ptl->y,motherNode->monoy);
        mz = SubCellDicision(p2ptl->z,motherNode->monoz);
        mnode = mx + 2*my + 4*mz;
        tmpnode[mnode].Nparticle ++;
        tmpptr = tmpnode[mnode].daughter;
        tmpptr2 = p2ptl->sibling;
        tmpnode[mnode].daughter = p2ptl;
        p2ptl->sibling = tmpptr;
        p2ptl = tmpptr2;
    }
    for(i=0;i<8;i++) if(tmpnode[i].Nparticle >0) break;

    TStruct *firstDaughter,*newDaughter;
    firstDaughter = newDaughter =  nextFreeNode;
    /* Link to the first daughter or first particle */
    if( divideThisNode(motherNode, tmpnode[i].Nparticle) )
        motherNode->daughter = (void*)firstDaughter;
    else
        motherNode->daughter = (void*) tmpnode[i].daughter;

    lastOffspring = NULL;

    for(i=0;i<8;i++) {
        if( divideThisNode(motherNode, tmpnode[i].Nparticle) ){
            newDaughter->daughter = tmpnode[i].daughter;
            newDaughter->Nparticle = tmpnode[i].Nparticle;
            if(lastOffspring) ((GENERAL_TPtl_POINTER*) lastOffspring)->sibling = newDaughter;
            lastOffspring = newDaughter;
            newDaughter++;
        }
        else if(tmpnode[i].Nparticle >0) {
            tmpptr = tmpnode[i].daughter;
            if(lastOffspring) ((GENERAL_TPtl_POINTER*)lastOffspring)->sibling = tmpptr;
            for(;tmpptr;tmpptr=tmpptr->sibling) lastOffspring = tmpptr;
        }
    }
	// This is to link the last offspring to mother's sibling.
    ((GENERAL_TPtl_POINTER *)lastOffspring)->sibling = motherNode->sibling;
    nextFreeNode = newDaughter;

    if(recursiveflag == RECURSIVE){
        TStruct *NextJobNode = firstDaughter;
        for(i=0;i<8;i++){
            if(divideThisNode(motherNode, tmpnode[i].Nparticle) ){
                nextFreeNode = divide_node(NextJobNode,nextFreeNode, thetasq,recursiveflag);
                NextJobNode ++;
            }
        }
    }
    return nextFreeNode;
}

void Make_Tree(
		TStruct *TREE_START, 
		size_t nnode, 
		TPtlStruct *ptl, 
		size_t np, 
		float thetasq,
		int recursiveflag
		){
    size_t i;
    TStruct *nextFreeNode = TREE_START+1;
    TREE_START->sibling = NULL;
    for(i=0;i<np;i++) {
        ptl[i].sibling = ptl + (i +1);
        ptl[i].type = TYPE_PTL;
    }
    ptl[np-1].sibling = NULL;
    TREE_START->daughter = &(ptl[0]);
    TREE_START->Nparticle = np;
    if(recursiveflag == RECURSIVE) nextFreeNode = divide_node(TREE_START, nextFreeNode, 
			thetasq,recursiveflag);
    else if(recursiveflag == SERIALIZED){
        TStruct *work;
        for(work=TREE_START;nextFreeNode-work >0; work++){
            nextFreeNode = divide_node(work, nextFreeNode, thetasq,SERIALIZED);
        }
    }
    else if(recursiveflag == PTHREAD)
    {
        TStruct *work = TREE_START;
        do{
            nextFreeNode = divide_node(work, nextFreeNode, thetasq,PTHREAD);
            work ++;
        }
        while( work < nextFreeNode && (nextFreeNode-work) < 64);
        size_t twork = (nextFreeNode-work);
		if(twork <=0) return;

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
            TStruct *threadnextFreeNode = nextFreeNode + pid * ( (nnode - (nextFreeNode-TREE_START))/ npid);
			size_t navail =  nnode - (nextFreeNode-TREE_START);
			threadnextFreeNode = calThreadFreeNodeStart(navail, twork, mys,myf, work, nextFreeNode);
            for(j=mys;j<myf;j++)
            {
                threadnextFreeNode = divide_node(work+j,threadnextFreeNode,thetasq,RECURSIVE);
            }
        }
    }
}


FoFTStruct *FoF_divide_node(
		FoFTStruct *motherNode, 
		FoFTStruct *nextFreeNode, 
		int recursiveflag
		){
    FoFTStruct *p2tree, tmpnode[8];
    void *lastOffspring;
    FoFTPtlStruct *p2ptl, *tmpptr, *tmpptr2, *nodeparticles;
    size_t i, j, k, mnode, mx,my,mz;
    float tmpx,tmpy,tmpz,tmpdist2, distmax;
    float ptlmass;

    nodeparticles = motherNode->daughter;
    motherNode->type = TYPE_TREE;
    motherNode->Nparticle = 0;
    motherNode->monox = motherNode->monoy = motherNode->monoz = 0;
	motherNode->maxlink02 = 0;

    for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
        motherNode->Nparticle += 1;
        motherNode->monox += p2ptl->x;
        motherNode->monoy += p2ptl->y;
        motherNode->monoz += p2ptl->z;
    }
    motherNode->monox /= motherNode->Nparticle;
    motherNode->monoy /= motherNode->Nparticle;
    motherNode->monoz /= motherNode->Nparticle;
    distmax = -1.E20;
    float q0,q1,q2,q3,q4,q5;
    q0=q1=q2=q3=q4=q5=0;
    for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
        tmpx = p2ptl->x - motherNode->monox;
        tmpy = p2ptl->y - motherNode->monoy;
        tmpz = p2ptl->z - motherNode->monoz;
        tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
        distmax = MAX(distmax, tmpdist2);
		motherNode->maxlink02 = MAX(motherNode->maxlink02,p2ptl->link02);
    }
    motherNode->nodesize = sqrt(distmax);

    for(i=0;i<8;i++) {
        tmpnode[i].sibling = tmpnode[i].daughter = NULL;
        tmpnode[i].Nparticle = 0;
    }
    for(p2ptl=nodeparticles;p2ptl;) {
        mx = SubCellDicision(p2ptl->x,motherNode->monox);
        my = SubCellDicision(p2ptl->y,motherNode->monoy);
        mz = SubCellDicision(p2ptl->z,motherNode->monoz);
        mnode = mx + 2*my + 4*mz;
        tmpnode[mnode].Nparticle ++;
        tmpptr = tmpnode[mnode].daughter;
        tmpptr2 = p2ptl->sibling;
        tmpnode[mnode].daughter = p2ptl;
        p2ptl->sibling = tmpptr;
        p2ptl = tmpptr2;
    }
    for(i=0;i<8;i++) if(tmpnode[i].Nparticle >0) break;

    FoFTStruct *firstDaughter,*newDaughter;
    firstDaughter = newDaughter =  nextFreeNode;
    /* Link to the first daughter or first particle */
    if( divideThisNode(motherNode, tmpnode[i].Nparticle) )
        motherNode->daughter = (void*)firstDaughter;
    else
        motherNode->daughter = (void*) tmpnode[i].daughter;

    lastOffspring = NULL;

    for(i=0;i<8;i++) {
        if( divideThisNode(motherNode, tmpnode[i].Nparticle) ){
            newDaughter->daughter = tmpnode[i].daughter;
            newDaughter->Nparticle = tmpnode[i].Nparticle;
            if(lastOffspring) ((GENERAL_TPtl_POINTER*) lastOffspring)->sibling = newDaughter;
            lastOffspring = newDaughter;
            newDaughter++;
        }
        else if(tmpnode[i].Nparticle >0) {
            tmpptr = tmpnode[i].daughter;
            if(lastOffspring) ((GENERAL_TPtl_POINTER*)lastOffspring)->sibling = tmpptr;
            for(;tmpptr;tmpptr=tmpptr->sibling) lastOffspring = tmpptr;
        }
    }
    ((GENERAL_TPtl_POINTER *)lastOffspring)->sibling = motherNode->sibling;
    nextFreeNode = newDaughter;


    if(recursiveflag == RECURSIVE){
        FoFTStruct *NextJobNode = firstDaughter;
        for(i=0;i<8;i++){
            if(divideThisNode(motherNode, tmpnode[i].Nparticle) ){
                nextFreeNode = FoF_divide_node(NextJobNode,nextFreeNode, recursiveflag);
                NextJobNode ++;
            }
        }
    }
    return nextFreeNode;
}

void FoF_Make_Tree(FoFTStruct *TREE_START, size_t nnode, FoFTPtlStruct *ptl, size_t np, int recursiveflag){
    size_t i;
    FoFTStruct *nextFreeNode = TREE_START+1;
    TREE_START->sibling = NULL;
	/*
    for(i=0;i<np;i++) {
        ptl[i].sibling = ptl + i +1;
        ptl[i].type = TYPE_PTL;
        ptl[i].included = NO;
    }
    ptl[np-1].sibling = NULL;
	*/
	FoFTPtlStruct *tmp = ptl;
	while(tmp){
		tmp->type = TYPE_PTL;
		tmp->included = NO;
		tmp = tmp->sibling;
	}

    TREE_START->daughter = &(ptl[0]);
    TREE_START->Nparticle = np;
    if(recursiveflag == RECURSIVE) nextFreeNode = FoF_divide_node(TREE_START, nextFreeNode, recursiveflag);
    else if(recursiveflag == SERIALIZED){
        FoFTStruct *work;
        for(work=TREE_START;nextFreeNode-work >0; work++){
            nextFreeNode = FoF_divide_node(work, nextFreeNode, SERIALIZED);
        }
    }
    else if(recursiveflag == PTHREAD)
    {
        FoFTStruct *work = TREE_START;
        do{
            nextFreeNode = FoF_divide_node(work, nextFreeNode, PTHREAD);
            work ++;
        }
        while( work < nextFreeNode && (nextFreeNode-work) < 64);
        size_t twork = (nextFreeNode-work);
		if(twork <=0) return;

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
            FoFTStruct *threadnextFreeNode = nextFreeNode + 
				pid * ( (nnode - (nextFreeNode-TREE_START)+npid-1)/ npid);
			size_t navail =  nnode - (nextFreeNode-TREE_START);
			threadnextFreeNode = calFoFThreadFreeNodeStart(navail, twork, 
					mys,myf, work, nextFreeNode);
            for(j=mys;j<myf;j++)
            {
                threadnextFreeNode = FoF_divide_node(work+j,threadnextFreeNode,RECURSIVE);
            }
        }
    }
}

/*
inline enum where FoF_InsideOpen(FoFPosition p, FoFTStruct *tree, float maxdist) {
    float tmpx = p.x-tree->monox;
    float tmpy = p.y-tree->monoy;
    float tmpz = p.z-tree->monoz;
    float dist = sqrt(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);
    if(dist+tree->nodesize <= maxdist) return IN;
    else if(dist-tree->nodesize > maxdist) return OUT;
    else return CROSS;

}
int pfof_open(FoFPosition p, FoFTStruct *tree, float fof_link){
    float tmpx, tmpy, tmpz;
    float dist2, dist, r, diff;
    float ratio;
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
size_t pnew_fof_link(FoFPosition *p, FoFTStruct *tree, float fof_link, FoFPosition *linked, size_t nhalo){
    size_t now, ncount = 0;
    void *ptr, *optr, *nptr;
    float tmpx,tmpy,tmpz;
    float fof_link2, dist2;
    float Lpx, Lpy, Lpz;
    FoFPosition point;
    fof_link2 = fof_link * fof_link;
    ncount = now = 0;
    point = *p;

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


*/
