/******************************************************************************/
/******************************************************************************/
/*******                NEAREST NEIGHBOR SEARCHING ALGORITHM         **********/
/******************************************************************************/
/******************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<float.h>
#include<unistd.h>
#ifdef _OPENMP
#include<omp.h>
#endif

#include "Memory.h"

#include "ramses.h"
#include "tree.h"
#include "params.h"


#define nodesize dist


int nullfct0();
int nullfct1();

float W4(float ,float );

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
#define OFFSET 0.001
Box findbox(TPtlStruct *ptl, int nend){
    float xmin,ymin,zmin,xmax,ymax,zmax;
    float width;
    Box box;
    int i;
    box.x = box.y = box.z = 0.;
    width=xmax=ymax=zmax = -1.E25;
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



#define SubCellDecision(a,b) ((a)>(b)? 1:0)
#define DivideNode(ThisNode,nparticles) ((ThisNode->nodesize > 0.5*MINCELLWIDTH ? yes:no) && (nparticles>=MIN_CELL_PARTICLE_NUM ? yes:no))

TStruct *calThreadFreeNodeStart(size_t navail, size_t twork, size_t mys,size_t myf, 
		TStruct *work, TStruct *nextFreeNode){
	TStruct *threadFreeNode, *tmp;
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


TStruct *new_divide_node_Near(TStruct *motherNode, TStruct *nextFreeNode, int recursiveflag)
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
	for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling) {
		tmpx = p2ptl->x - motherNode->monox;
		tmpy = p2ptl->y - motherNode->monoy;
		tmpz = p2ptl->z - motherNode->monoz;
		tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
		distmax = MAX(distmax, tmpdist2);
	}
	motherNode->nodesize = sqrt(distmax);

	for(i=0;i<8;i++) {
		tmpnode[i].sibling = tmpnode[i].daughter = NULL;
		tmpnode[i].Nparticle = 0;
	}
	for(p2ptl=nodeparticles;p2ptl;) {
		mx = SubCellDecision(p2ptl->x,motherNode->monox);
		my = SubCellDecision(p2ptl->y,motherNode->monoy);
		mz = SubCellDecision(p2ptl->z,motherNode->monoz);
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
	// Link mother node to the first daughter or particle (first offsprig) 
	if( divideThisNode(motherNode, tmpnode[i].Nparticle) )
		motherNode->daughter = (void*)firstDaughter;
	else 
		motherNode->daughter = (void*) tmpnode[i].daughter;

	// make links between the current offspring
	// lastOffsping: the terminal point of the current offsprings
	lastOffspring = NULL;
	for(i=0;i<8;i++) {
		if( divideThisNode(motherNode, tmpnode[i].Nparticle)){
			newDaughter->daughter = tmpnode[i].daughter;
			newDaughter->Nparticle = tmpnode[i].Nparticle;
			if(lastOffspring) ((GENERAL_TPtl_POINTER*) lastOffspring)->sibling = newDaughter;
			lastOffspring = (void *)newDaughter; // update lastOffspring
			newDaughter++;
		}
		else if(tmpnode[i].Nparticle >0) {
			tmpptr = tmpnode[i].daughter;
			if(lastOffspring) ((GENERAL_TPtl_POINTER*)lastOffspring)->sibling = tmpptr;
			for(;tmpptr;tmpptr=tmpptr->sibling) lastOffspring = (void*) tmpptr;
		}
	}
	// close by setting the sibling of the last offspring to the mother's sibling.
	((GENERAL_TPtl_POINTER *)lastOffspring)->sibling = motherNode->sibling;
	nextFreeNode = newDaughter;


	if(recursiveflag == RECURSIVE){
		TStruct *NextJobNode = firstDaughter;
		for(i=0;i<8;i++){
			if(divideThisNode(motherNode, tmpnode[i].Nparticle)){
				nextFreeNode = new_divide_node_Near(NextJobNode,nextFreeNode, recursiveflag);
				NextJobNode ++;
			}
		}
	}
	return nextFreeNode;
}

void new_Make_Tree_Near(
		TStruct *TREE_START, 
		size_t nnode, 
		TPtlStruct *ptl, 
		size_t np, 
		int recursiveflag){
	size_t i;
	TStruct *nextFreeNode = TREE_START+1;
	TREE_START->sibling = NULL;
	for(i=0;i<np;i++) {
		ptl[i].sibling = ptl + i +1;
		ptl[i].type = TYPE_PTL;
	}
	ptl[np-1].sibling = NULL;
	TREE_START->daughter = &(ptl[0]);
	TREE_START->Nparticle = np;
	if(recursiveflag == RECURSIVE) nextFreeNode = new_divide_node_Near(TREE_START, nextFreeNode, recursiveflag);
	else if(recursiveflag == SERIALIZED){
		TStruct *work;
		for(work=TREE_START;nextFreeNode-work >0; work++){
			nextFreeNode = new_divide_node_Near(work, nextFreeNode, SERIALIZED);
		}
	}
	else if(recursiveflag == PTHREAD) 
	{
		TStruct *work = TREE_START;
		do{
			nextFreeNode = new_divide_node_Near(work, nextFreeNode, PTHREAD);
			work ++;
		}
		while( work < nextFreeNode &&  (nextFreeNode-work) < 64);
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
			TStruct *threadnextFreeNode = nextFreeNode + 
				pid * ( (nnode - (nextFreeNode-TREE_START)+npid-1)/ npid);
			size_t navail =  nnode - (nextFreeNode-TREE_START);
			threadnextFreeNode = calThreadFreeNodeStart(navail, twork, 
					mys,myf, work, nextFreeNode);
			for(j=mys;j<myf;j++) {
				threadnextFreeNode = new_divide_node_Near(work+j,threadnextFreeNode,RECURSIVE);
			}
		}
	}
}



enum boolean near_open(particle *point, TStruct *tree, int npneigh, float maxdist , int Num_neighbor){
	dptype tmpx,tmpy,tmpz,dist2, r2, dist, r;
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


typedef struct nearestneighbor{
	dptype dist2;
	float mass;
	long long indx;
} Neighbor;



#define INSERT(dist2, maxdist, maxdist2, npneigh, neighbor) do{\
	for(i=0;i<npneigh;i++) if(neighbor[i].dist2>dist2) break;\
	for(j=MIN(Num_neighbor,npneigh-1)-1;j>=i;j--) {\
		neighbor[j+1] = neighbor[j];\
	}\
	neighbor[i].dist2 = dist2;\
	neighbor[i].mass = ( (TPtlStruct*)ptr)->mass;\
	neighbor[i].indx = ( (TPtlStruct*)ptr)-ptl;\
	npneigh ++;\
	npneigh=MIN(Num_neighbor,npneigh);\
	maxdist2 = (neighbor[npneigh-1].dist2);\
	maxdist = sqrtf(maxdist2);\
}while(0)

int sortneighdist2(const void *a, const void *b){
	Neighbor *aa = (Neighbor*)a;
	Neighbor *bb = (Neighbor*)b;
	if(aa->dist2 < bb->dist2) return -1;
	else if(aa->dist2 > bb->dist2) return 1;
	else return 0;
}

#define INSERT_PreCopied(ptr, dist2, maxdist, maxdist2, npneigh, neighbor) do{\
	int _iflag = 0;\
	int _i,_j;\
	long long _ioffset = ( (TPtlStruct*)ptr)-ptl;\
	for(_i=0;_i<npneigh;_i++) {\
		if(neighbor[_i].indx == _ioffset) {\
			_iflag = 1;\
			break;\
		}\
	}\
	if(_iflag == 0) \
	{\
		for(_i=0;_i<npneigh;_i++) if(neighbor[_i].dist2>dist2) break;\
		for(_j=MIN(npneigh,Num_neighbor-1)-1;_j>=_i;_j--) {\
			neighbor[_j+1] = neighbor[_j];\
		}\
		neighbor[_i].dist2 = dist2;\
		neighbor[_i].indx = _ioffset;\
		npneigh ++;\
		npneigh=MIN(Num_neighbor,npneigh);\
		maxdist2 = (neighbor[npneigh-1].dist2);\
		maxdist = sqrt(maxdist2);\
	}\
}while(0)



int Find_Near_PreCopied(
		particle *point, 
		int Num_neighbor, 
		TStruct *tree, 
		TPtlStruct *ptl, 
		Neighbor *neighbor,
		int npneigh
		){
	int i,j,k;
	dptype dist2, tmpx,tmpy,tmpz;
	void *ptr;
	ptr = (void*)tree;
	dptype maxdist2, maxdist;
	if(npneigh==0) {
		maxdist2 = maxdist = 1.e20;
	}
	else{
		maxdist2 = neighbor[npneigh-1].dist2;
		maxdist = sqrt(maxdist2);
	}
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
					INSERT_PreCopied(ptr, dist2, maxdist, maxdist2, npneigh, neighbor);
				}
				ptr = (void *) (((TPtlStruct*)ptr)->sibling);
		}
	}
	return npneigh;
}

int Find_Near(
		particle *point, 
		int Num_neighbor, 
		TStruct *tree, 
		TPtlStruct *ptl, 
		float *maxr, 
		long long *nearindx,
		dptype *DIST2,
		float *MASS
		){
	int i,j,k;
	dptype dist2, tmpx,tmpy,tmpz;
	dptype maxdist, maxdist2;
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
	for(i=0;i<npneigh;i++) {
		nearindx[i] = neighbor[i].indx;
		DIST2[i] = neighbor[i].dist2;
		MASS[i] = neighbor[i].mass;
	}
	return npneigh;
}



int nearConstOpen(particle *point, TStruct *tree, dptype constR){
	dptype tmpx = point->x - tree->monox;
	dptype tmpy = point->y - tree->monoy; 
	dptype tmpz = point->z - tree->monoz;
	dptype r = tree->nodesize;
    dptype dist = sqrtf(tmpx*tmpx+ tmpy*tmpy + tmpz*tmpz); 
	if(dist-r > constR) return NO; 
	else return YES;
}
dptype getDenConstR (
        particle *point,
        TStruct *tree,
        TPtlStruct *ptl,
		dptype constR,
		int *count
        ){
	int icount = 0;
	dptype tmpx,tmpy,tmpz,dist2,mass;
    dptype den = 0;
    void *ptr = (void*)tree;
    while(ptr != NULL){
        switch( ((TYPE*)ptr)->type) {
            case TYPE_TREE:
                switch( nearConstOpen(point, ptr, constR)){
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
                mass = ((TPtlStruct *)ptr)->mass;
                if(dist2<constR){
                    den += W4(dist2,constR/2.) * mass;
					icount ++;
                }
                ptr = (void *) (((TPtlStruct*)ptr)->sibling);
        }
    }
	*count = icount;
    return den;
}

int Find_Near2(
		particle *point, 
		int Num_neighbor, 
		TStruct *tree, 
		TPtlStruct *ptl, 
		size_t nptl, 
		float *maxr, 
		long long *nearindx,
		dptype *DIST2,
		float *MASS
		){
	size_t i,j,k;
	size_t iptl = 0;
	dptype dist2, tmpx,tmpy,tmpz;
	dptype maxdist, maxdist2;
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
					neighbor[iptl].mass = ptl->mass;
					neighbor[iptl].indx = ptl->indx;
					iptl ++;
					if(iptl>Num_neighbor) maxdist2 = MIN(maxdist2, dist2);
					else maxdist2 = MAX(maxdist2, dist2);
				}
				ptr = (void *) (((TPtlStruct*)ptr)->sibling);
		}
	}
	npneigh = 0;
	for(i=0;i<iptl;i++){
		if(neighbor[i].dist2 <= maxdist2) {
			nearindx[npneigh] = neighbor[i].indx;
			DIST2[npneigh] = neighbor[i].dist2;
			MASS[npneigh] = neighbor[i].mass;
			npneigh ++;
		}
	}
	return npneigh;
}

double f(double);
#define AMPW41 (-0.2387324)
#define AMPW42 (0.07957747)
float W4(float dist,float h){
    float twomr,r;
    r = dist/h;
    if(r<1.){
        /*
        return 1/PI*(1.-1.5*dist*dist+0.75*dist*dist*dist);
        */
        return AMPW41*((2.-r)*r*r-1.333333)/(h*h*h);
    }
    else if(r < 2.){
        twomr = (2.-r)/h;
        /*
        return 1/PI*0.25*twomr*twomr*twomr;
        */
        return AMPW42*twomr*twomr*twomr;
    }
    else {
        return 0.;
    }
}
#undef AMPW41
#undef AMPW42


void findsphdensity(SimpleBasicParticleType *bp,int np,long long *nearindex, int Numnear,
        float *densph){
    float *h;
    double std,mean;
    int ntmp;
    float tmpx,tmpy,tmpz;
    dptype *dist2;
	float *mass;
    float fplmf,ptlmass;
    long long i,j,k;
    int N,M;
    TPtlStruct *ptl;
    TStruct *TREE;
    float theta = 1.;
    float wtime;
    long iseed=-9;
    float x0,y0,z0,pscale;
    float Rmin,Rmax,size;
    int index;
    particle *p;
	p = (particle *) Malloc(sizeof(particle)*np,PPTR(p));
    for(i=0;i<np;i++){
        p[i].x = bp[i].x;
        p[i].y = bp[i].y;
        p[i].z = bp[i].z;
    }
    ptl = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
    for(i=0;i<np;i++){
        ptl[i].type = TYPE_PTL;
        ptl[i].x = bp[i].x;
        ptl[i].y = bp[i].y;
        ptl[i].z = bp[i].z;
        ptl[i].sibling = ptl+(i+1);
        ptl[i].mass = bp[i].mass;
    }
    ptl[np-1].sibling = NULL;
	size_t nnode = MAX(np/2,65*10000);
    TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE));
	int recursiveflag;
	if(np > 65*10000) recursiveflag = PTHREAD;
	else recursiveflag = RECURSIVE;

	recursiveflag = RECURSIVE;
	new_Make_Tree_Near(TREE,nnode, ptl,np,recursiveflag);
    float mindist;
    mindist = 2.e25;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k) schedule(guided)
#endif
    for(i=0;i<np;i++){
        long long tmpindx[Numnear];
        dptype dist2[Numnear];
        float mass[Numnear];
        int res;
        float neardist;
        k = i*Numnear;
        res = Find_Near(p+i,Numnear,TREE,ptl,&neardist,tmpindx,dist2, mass);
        densph[i] = 0;
        for(j=0;j<res;j++){
            nearindex[k+j] = tmpindx[j];
            dist2[j] = sqrtf(dist2[j]);
            densph[i] += W4(dist2[j],neardist/2.);
        }
    }
    Free(TREE);Free(ptl);Free(p);
    /*
    printf("ins=", tmpindx);
    fprintf(stdout,"Terminated with success\n");fflush(stdout);
    */
    DEBUGPRINT("exiting findsphdensity\n");
}

void starfindsphdensity(SimpleBasicParticleType *bp,int np,long long *nearindex, int Numnear,
        float *densph){
    float *h;
    double std,mean;
    int ntmp;
    float tmpx,tmpy,tmpz;
    dptype *dist2;
	float *mass;
    float fplmf,ptlmass;
    long long i,j,k;
    int N,M;
    TPtlStruct *ptl;
    TStruct *TREE;
    float theta = 1.;
    float wtime;
    long iseed=-9;
    float x0,y0,z0,pscale;
    float Rmin,Rmax,size;
    int index;
    int mp;
    particle *p;


    p = (particle *) Malloc(sizeof(particle)*np,PPTR(p));
    for(i=0;i<np;i++){
        p[i].x = bp[i].x;
        p[i].y = bp[i].y;
        p[i].z = bp[i].z;
    }

    ptl = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
    mp = 0;
    for(i=0;i<np;i++){
        if(bp[i].type == TYPE_STAR){
            ptl[mp].type = TYPE_PTL;
            ptl[mp].x = bp[i].x;
            ptl[mp].y = bp[i].y;
            ptl[mp].z = bp[i].z;
            ptl[mp].sibling = ptl+(mp+1);
            ptl[mp].mass = bp[i].mass;
            mp++;
        }
    }
    ptl[mp-1].sibling = NULL;
	size_t nnode = MAX(mp/2,65*10000);
    TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE));

	int recursiveflag;
	if(np > 65*10000) recursiveflag = PTHREAD;
	else recursiveflag = RECURSIVE;

	recursiveflag = RECURSIVE;
    new_Make_Tree_Near(TREE,nnode, ptl,mp,recursiveflag);
    DEBUGPRINT("after Star_Make_Tree_near for nstar = %d\n",mp);
    float mindist;
    mindist = 2.e25;
//    int NumNearDen = NUMNEARDEN;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k) schedule(dynamic)
#endif
    for(i=0;i<np;i++){
        int res;
        long long tmpindx[NUMNEARDEN];
        dptype tmpd2[NUMNEARDEN];
        float tmpmass[NUMNEARDEN];
        float neardist;
        res = Find_Near(p+i,NUMNEARDEN,TREE,ptl,&neardist,tmpindx,tmpd2, tmpmass);
        densph[i] = 0;
        for(j=0;j<res;j++){
            tmpd2[j] = sqrtf(tmpd2[j]);
            densph[i] += W4(tmpd2[j],neardist/2.) * tmpmass[j];
        }
    }
    DEBUGPRINT("End of allocating the stellar density for np = %d\n",np);

	// now find the k-nearest neighbors.
    for(i=0;i<np;i++){
        ptl[i].type = TYPE_PTL;
        ptl[i].x = bp[i].x;
        ptl[i].y = bp[i].y;
        ptl[i].z = bp[i].z;
        ptl[i].sibling = ptl+i+1;
        ptl[i].mass = bp[i].mass;
    }
    ptl[np-1].sibling = NULL;
	recursiveflag = RECURSIVE;
	nnode = MAX(np/2, 65*10000);
	TREE = Realloc(TREE, sizeof(TStruct)*nnode);
    new_Make_Tree_Near(TREE,nnode, ptl,np,recursiveflag);
    DEBUGPRINT("before Make_Tree_near2 with np= %d\n", np);
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k) schedule(dynamic)
#endif
    for(i=0;i<np;i++){
        int res;
        long long tmpindx[NUMNEIGHBOR];
        dptype dist2[NUMNEIGHBOR];
        float mass[NUMNEIGHBOR];
        float neardist;
        k = i*Numnear;
        res = Find_Near(p+i,Numnear,TREE,ptl,&neardist,tmpindx,dist2, mass);
        if(res != Numnear){
            printf("error occurred %d %d : %d\n", res, Numnear, np);
            exit(99);
        }
        for(j=0;j<res;j++){
            nearindex[k+j] = tmpindx[j];
        }
    }
//    Free(mass);Free(dist2); Free(h);
    Free(TREE);Free(ptl);Free(p);
    DEBUGPRINT("end of starfinddensity()\n");
}

void findStellarCore(
		SimpleBasicParticleType *bp,int np,
		int Numnear,
        float *densph,
		Coretype **Core,
		int *NumCore,  // return value for number of core
		int maxnumcore,
		long long **Nearindex){
	Coretype *core = *Core;
	int numcore;
//    float *h;
    double std,mean;
    int ntmp;
    float tmpx,tmpy,tmpz;
    float fplmf,ptlmass;
    long long i,j,k;
    int N,M;
    TPtlStruct *ptl;
    TStruct *TREE;
    float theta = 1.;
    float wtime;
    long iseed=-9;
    float x0,y0,z0,pscale;
    float Rmin,Rmax,size;
    int index;
    int mp;
	int nstar;
    particle *p;


    p = (particle *) Malloc(sizeof(particle)*np,PPTR(p));
    for(i=0;i<np;i++){
        p[i].x = bp[i].x;
        p[i].y = bp[i].y;
        p[i].z = bp[i].z;
    }

    ptl = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
    nstar = 0;
    for(i=0;i<np;i++){
        if(bp[i].type == TYPE_STAR){
            ptl[nstar].type = TYPE_PTL;
            ptl[nstar].x = bp[i].x;
            ptl[nstar].y = bp[i].y;
            ptl[nstar].z = bp[i].z;
            ptl[nstar].sibling = ptl+(nstar+1);
            ptl[nstar].mass = bp[i].mass;
            ptl[nstar].indx = i;
            nstar++;
        }
    }
    ptl[nstar-1].sibling = NULL;
	size_t nnode = MAX(nstar/2,65*10000);
    DEBUGPRINT("nstar = %d, nnode= %zu\n",nstar, nnode);
    TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE));


	int recursiveflag;
	if(nnode > 65*10000) recursiveflag = PTHREAD;
	else recursiveflag = RECURSIVE;
	recursiveflag = RECURSIVE;
    new_Make_Tree_Near(TREE,nnode, ptl,nstar,recursiveflag);
    DEBUGPRINT("after Star_Make_Tree_near for nstar = %d\n",nstar);
//  all to star
#ifdef _OPENMP
#pragma omp parallel for private(i,j) schedule(guided)
#endif
    for(i=0;i<np;i++){
        int res, icount;
        long long tmpindx[NUMNEARDEN];
        dptype tmpd2[NUMNEARDEN];
        float tmpmass[NUMNEARDEN];
        float neardist;
#if defined mADV   || defined ADV
		dptype constR = MIN_CONST_R_SMOOTHING;
		densph[i] = getDenConstR(p+i, TREE,ptl, constR, &icount);
		if(icount < NUMNEARDEN){
	        res = Find_Near(p+i,NUMNEARDEN,TREE,ptl,&neardist,tmpindx,tmpd2, tmpmass);
       		densph[i] = 0;
	        for(j=0;j<res;j++){ 
				tmpd2[j] = sqrtf(tmpd2[j]); 
				densph[i] += W4(tmpd2[j],neardist/2.) * tmpmass[j];
			}
		}
#else
        res = Find_Near(p+i,NUMNEARDEN,TREE,ptl,&neardist,tmpindx,tmpd2, tmpmass);
        densph[i] = 0;
        for(j=0;j<res;j++){
            tmpd2[j] = sqrtf(tmpd2[j]);
            densph[i] += W4(tmpd2[j],neardist/2.) * tmpmass[j];
        }
#endif
    }
    DEBUGPRINT("End of getDenConstR()\n");




	// now find the k-nearest neighbors.
	mp = 0;
    for(i=0;i<np;i++){
		if(bp[i].type == TYPE_STAR){
			p[mp].x = bp[i].x;
			p[mp].y = bp[i].y;
			p[mp].z = bp[i].z;
			p[mp].indx = i;
			mp++;
		}
    }
	int numSNeigh = MIN(nstar, NUMSTELLARNEIGHBORS);
	long long *nearindex = *Nearindex 
		= (long long*)Malloc(sizeof(long long)*numSNeigh*nstar, PPTR(*Nearindex));
//  star to star
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k) schedule(guided)
#endif
    for(i=0;i<nstar;i++){
        int res;
        long long tmpindx[NUMSTELLARNEIGHBORS];
        float neardist;
		dptype dist2[NUMSTELLARNEIGHBORS];
		float mass[NUMSTELLARNEIGHBORS];
        k = i*numSNeigh;
        res = Find_Near(p+i,numSNeigh,TREE,ptl,&neardist,tmpindx,
				dist2, mass);
        if(res != numSNeigh){
            DEBUGPRINT("error occurred %d %d : %d\n", res, numSNeigh, np);
            exit(99);
        }
        for(j=0;j<res;j++){
            nearindex[k+j] = tmpindx[j];
        }
    }

//	Free(h);
    Free(TREE);
    DEBUGPRINT("end of stellarneighbor()\n");
	int istar = 0;
	numcore  = 0;
	for(i=0;i<np;i++){
		if(bp[i].type == TYPE_STAR && densph[i]>PEAKTHRESHOLD){
			k = istar*numSNeigh;
			int iflag = 0;
			for(j=0;j<numSNeigh;j++){
				if(densph[i] < densph[ptl[nearindex[k+j]].indx]) {
					iflag = 1;
					break;
				}
			}
			if(iflag ==0){
				long long ibp = p[istar].indx;
				core[numcore].peak = ibp;
				core[numcore].cx = bp[ibp].x;
				core[numcore].cy = bp[ibp].y;
				core[numcore].cz = bp[ibp].z;
				core[numcore].density = densph[i];
				numcore ++;
				if(numcore > maxnumcore-10){
					maxnumcore += MAXNUMCORE;
					*Core = Realloc(*Core, sizeof(Coretype)*maxnumcore);
					core = *Core;
				}
			}
		}
		if(bp[i].type == TYPE_STAR) istar ++;
	}
	Free(ptl); Free(p);
	DEBUGPRINT("The number of cores : %d and before MergingPeak\n", numcore);
	for(i=0;i<numcore;i++){
		DEBUGPRINT("p%lld: %g %g %g\n", i,core[i].cx, core[i].cy,core[i].cz);
	}




	int MergingPeak(SimpleBasicParticleType *, int, Coretype *, int, int);
	if(numcore >10) numcore = MergingPeak(bp,np, core, numcore,0);
	core = *Core = Realloc(*Core, sizeof(Coretype)*numcore);

	DEBUGPRINT("The number of cores : %d and after MergingPeak\n", numcore);

	nearindex = *Nearindex = Realloc(*Nearindex, sizeof(long long)*np*Numnear);
    p = (particle *) Malloc(sizeof(particle)*np,PPTR(p));
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<np;i++){
        p[i].x = bp[i].x;
        p[i].y = bp[i].y;
        p[i].z = bp[i].z;
    }
    ptl = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<np;i++){ 
		ptl[i].type = TYPE_PTL; 
		ptl[i].x = bp[i].x; 
		ptl[i].y = bp[i].y; 
		ptl[i].z = bp[i].z; 
		ptl[i].sibling = ptl+(i+1); 
		ptl[i].mass = bp[i].mass; 
    }
    ptl[np-1].sibling = NULL;
	nnode = MAX(np/2,65*10000);
    TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE));

	if(nnode > 65*10000) recursiveflag = PTHREAD;
	else recursiveflag = RECURSIVE;
//  all to all
	recursiveflag = RECURSIVE;
    new_Make_Tree_Near(TREE,nnode, ptl,np,recursiveflag);
    DEBUGPRINT("after All_Make_Tree_near for np = %d\n",np);
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k) schedule(guided)
#endif
    for(i=0;i<np;i++){
//		DEBUGPRINT("p%d  is now being processed: %d\n", i, recursiveflag);
        int res;
        long long tmpindx[NUMNEIGHBOR];
        float mass[NUMNEIGHBOR];
        dptype dist2[NUMNEIGHBOR];
        float neardist;
        k = i*Numnear;
        res = Find_Near(p+i,Numnear,TREE,ptl,&neardist,tmpindx,
				dist2, mass);
        if(res != Numnear){
            DEBUGPRINT("error occurred %d %d : %d\n", res, Numnear, np);
            exit(99);
        }
        for(j=0;j<res;j++){
            nearindex[k+j] = tmpindx[j];
        }
    }

	Free(TREE);
	Free(ptl);
	Free(p);


	*NumCore = numcore;

	return;
}

void lagFindStellarCore(
		SimpleBasicParticleType *bp,
		int np,
		int Numnear,
        float *densph,
		Coretype **Core,
		int *NumCore,  // return value for number of core
		int maxnumcore,
		long long *nearindex, // index array for neighbor network
		int ptype
		){
	Coretype *core = *Core;
	int numcore;
//    float *h;
    double std,mean;
    int ntmp;
    float tmpx,tmpy,tmpz;
    float fplmf,ptlmass;
    long long i,j,k;
    int N,M;
    TPtlStruct *ptl;
    TStruct *TREE;
    float theta = 1.;
    float wtime;
    long iseed=-9;
    float x0,y0,z0,pscale;
    float Rmin,Rmax,size;
    int index;
    int mp;
    particle *p;





	size_t nnode;
	int recursiveflag; 


	int nx,ny,nz;
	long long mx,ncells;
	double xmin,ymin,zmin,xmax,ymax,zmax;
	double cellsize = TSC_CELL_SIZE;

	{
		float RG = Gaussian_Smoothing_Length;
		int nbuff = NCELLBUFF;
		float *denGrid;

		xmin = ymin = zmin = 1e20;
		xmax = ymax = zmax = -1e20;
		for(i=0;i<np;i++){
			xmin = MIN(xmin, bp[i].x);
			ymin = MIN(ymin, bp[i].y);
			zmin = MIN(zmin, bp[i].z);
			xmax = MAX(xmax, bp[i].x);
			ymax = MAX(ymax, bp[i].y);
			zmax = MAX(zmax, bp[i].z);
		}
		xmax = xmax + (xmax-xmin)*0.0001; // making a buffer to avoid the numerical exception
		ymax = ymax + (ymax-ymin)*0.0001;
		zmax = zmax + (zmax-zmin)*0.0001;
		nx = (xmax-xmin)/cellsize + nbuff; // buffer for the periodic boundary condition
		ny = (ymax-ymin)/cellsize + nbuff;
		nz = (zmax-zmin)/cellsize + nbuff;
		xmin -= cellsize*nbuff/2.; ymin -= cellsize*nbuff/2.; zmin -= cellsize*nbuff/2.;

		mx = 2*(nx/2+1);
		ncells = mx*ny*nz;
		denGrid = (float*)Malloc(sizeof(float)*ncells,PPTR(denGrid));
		{	
			/*
			LOGPRINT("np= %d nx/y/z= %d %d %d xmin/ymin/zmin = %g %g %g / %g %g %g cellsize= %g\n",
					np, nx,ny,nz,xmin,ymin,zmin,xmax,ymax,zmax,cellsize);
					*/
			void assign_density_TSC(SimpleBasicParticleType *, int, float *, int, int, int,
				double, double, double, double, int);
			assign_density_TSC(bp, np, denGrid, nx,ny,nz,xmin,ymin,zmin,cellsize,
					ptype);
		}
		if(0){

			FILE *wp = fopen("denmap.out","w");
			fwrite(&mx, sizeof(int), 1, wp);
			fwrite(&ny, sizeof(int), 1, wp);
			fwrite(&nz, sizeof(int), 1, wp);
			fwrite(denGrid, sizeof(float), mx*ny*nz,wp);
			fclose(wp);
		}
		{
			void gaussian_Smoothing(float *,int ,int ,int , double , float );
			gaussian_Smoothing(denGrid,nx,ny,nz, cellsize, RG);
		}
		if(0){
			FILE *wp = fopen("gS.denmap.out","w");
			fwrite(&mx, sizeof(int), 1, wp);
			fwrite(&ny, sizeof(int), 1, wp);
			fwrite(&nz, sizeof(int), 1, wp);
			fwrite(denGrid, sizeof(float), mx*ny*nz,wp);
			fclose(wp);
		}
		{
			void findDen(SimpleBasicParticleType *, int, float *, double, double, double,
					double, float *, int, int, int);
			findDen(bp,np, densph, xmin,ymin,zmin, cellsize, denGrid, nx,ny,nz);
		}
		if(0){
			float *checkden = (float*)Malloc(sizeof(float)*ncells,PPTR(checkden));
			for(i=0;i<ncells;i++) checkden[i] = 0;
			for(i=0;i<np;i++){
				int ix = (bp[i].x-xmin)/cellsize;
				int iy = (bp[i].y-ymin)/cellsize;
				if(densph[i] > checkden[ix + mx*(iy)]) checkden[ix + mx*(iy)] = densph[i];
			}
			FILE *wp = fopen("checkden.out","w");
            fwrite(&mx, sizeof(int), 1, wp);
            fwrite(&ny, sizeof(int), 1, wp);
            fwrite(&nz, sizeof(int), 1, wp);
            fwrite(checkden, sizeof(float), mx*ny*nz,wp);
            fclose(wp);
			exit(1);
		}


		LinkedListGrid *linkedListGrid= (LinkedListGrid*)Malloc(sizeof(LinkedListGrid)
				*ncells, PPTR(linkedListGrid));
		for(i=0;i<ncells;i++) {
			linkedListGrid[i].bp = NULL;
			linkedListGrid[i].np = 0;
		}
		int halfnbuff = NCELLBUFF/2;
		if(ptype == TYPE_STAR){
			int nstar = 0;
			for(i=0;i<np;i++){
				if(bp[i].type == TYPE_STAR){
					long long ir = rint((bp[i].x-xmin)/cellsize);
					long long jr = rint((bp[i].y-ymin)/cellsize);
					long long kr = rint((bp[i].z-zmin)/cellsize);
					long long ioff = ir+mx*(jr+ny*kr);
					SimpleBasicParticleType *tmp = linkedListGrid[ioff].bp;
					linkedListGrid[ioff].bp = bp+i;
					linkedListGrid[ioff].np ++;
					bp[i].bp = tmp;
					nstar ++;
				}
			}
			DEBUGPRINT("After building LinkedList with nstar= %d\n", nstar);
		}
		else if (ptype == TYPE_ALL) {
			int nparticles = 0;
			for(i=0;i<np;i++){
//				if(bp[i].type == TYPE_STAR)
				{
					long long ir = rint((bp[i].x-xmin)/cellsize);
					long long jr = rint((bp[i].y-ymin)/cellsize);
					long long kr = rint((bp[i].z-zmin)/cellsize);
					long long ioff = ir+mx*(jr+ny*kr);
					SimpleBasicParticleType *tmp = linkedListGrid[ioff].bp;
					linkedListGrid[ioff].bp = bp+i;
					linkedListGrid[ioff].np ++;
					bp[i].bp = tmp;
					nparticles ++;
				}
			}
			DEBUGPRINT("After building LinkedList with all_type = %d\n", nparticles);
		}
		numcore = 0;
		int nthreads;
#ifdef _OPENMP
#pragma omp parallel shared(nthreads)
#endif
		{
#ifdef _OPENMP
#pragma omp master
#endif
			{
				nthreads = Omp_get_num_threads();
			}
		}
		int num_cores[nthreads];
		for(i=0;i<nthreads;i++) num_cores[i] = 0;
		DEBUGPRINT("%d threads are ready\n", nthreads);
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
		for(k=nbuff/2;k<nz-nbuff/2;k++){
			int thread_id = Omp_get_thread_num();
			for(j=nbuff/2;j<ny-nbuff/2;j++){
				for(i=nbuff/2;i<nx-nbuff/2;i++){
					long long ioff = i + mx*(long long)(j+ny*k);
					if(denGrid[ioff] > PEAKTHRESHOLD){
						int peakflag=1;
						int i1,j1,k1;
						for(k1=-1;k1<2;k1++) for(j1=-1;j1<2;j1++) for(i1=-1;i1<2;i1++){
							long long joff = ioff+i1 + mx*(long long)(j1 + ny*k1);
							if(joff != ioff && denGrid[ioff]<=denGrid[joff] ){
								peakflag = 0;
								break;
							}
						}
						if(peakflag == 1 && linkedListGrid[ioff].np >0 ){
							SimpleBasicParticleType *tmp = linkedListGrid[ioff].bp;
							float maxden = -1.e20;
							int ibp=-1;
							while(tmp){
								ibp = tmp-bp;
								if(densph[ibp] > maxden){
									maxden = densph[ibp];
								}
								tmp = tmp->bp;
							}
							if(ibp>0) num_cores[thread_id] ++; 
						}
					}
				}
			}
		}
		numcore = 0;
		for(i=0;i<nthreads;i++) {
//			DEBUGPRINT("T%d has numcores= %d\n", i, num_cores[i]);
			numcore += num_cores[i];
		}
		DEBUGPRINT("The number of cores: %d in initial finding\n", numcore);

		int ioff_cores[nthreads];
		ioff_cores[0] =0;
		for(i=1;i<nthreads;i++) {
			ioff_cores[i] = ioff_cores[i-1] + num_cores[i-1];
//			DEBUGPRINT("T%d has core offset= %d\n", i, ioff_cores[i]);
		}
		numcore =0;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k) firstprivate(numcore)
#endif
        for(k=nbuff/2;k<nz-nbuff/2;k++){
            int thread_id = Omp_get_thread_num();
			int koff = ioff_cores[thread_id];
            for(j=nbuff/2;j<ny-nbuff/2;j++){
                for(i=nbuff/2;i<nx-nbuff/2;i++){
                    long long ioff = i + mx*(long long)(j+ny*k);
                    if(denGrid[ioff] > PEAKTHRESHOLD){
                        int peakflag=1;
                        int i1,j1,k1;
                        for(k1=-1;k1<2;k1++) for(j1=-1;j1<2;j1++) for(i1=-1;i1<2;i1++){
                            long long joff = ioff+i1 + mx*(long long)(j1 + ny*k1);
                            if(joff != ioff && denGrid[ioff]<=denGrid[joff]) {
                                peakflag = 0;
                                break;
                            }
                        }
                        if(peakflag == 1 && linkedListGrid[ioff].np >0 ){
                            SimpleBasicParticleType *tmp = linkedListGrid[ioff].bp;
                            float maxden = -1.e20;
                            long long ibp=-1;
							long long jbp;
                            while(tmp){
                                ibp = tmp-bp;
                                if(densph[ibp] > maxden){
                                    maxden = densph[ibp];
									jbp = ibp;
                                }
                                tmp = tmp->bp;
                            }
							if(ibp > 0){
	                            core[koff+numcore].peak = jbp; 
								core[koff+numcore].cx = bp[jbp].x; 
								core[koff+numcore].cy = bp[jbp].y; 
								core[koff+numcore].cz = bp[jbp].z; 
								core[koff+numcore].density = maxden;
								/*
								LOGPRINT("p%d has c%d with koff= %d xyz= %g %g %g rho= %g\n",
										thread_id, koff+numcore, koff, bp[jbp].x,
										bp[jbp].y, bp[jbp].z,maxden);
										*/
								numcore ++;
							}
                        }
                    }
                }
            }
        }
		numcore = 0;
		for(i=0;i<nthreads;i++) numcore += num_cores[i];
		if(ptype == TYPE_STAR){
			DEBUGPRINT("The number of stellar cores: %d before MergingPeak\n", numcore);
		}
		else {
			DEBUGPRINT("The number of All-type  cores: %d before MergingPeak\n", numcore);
		}
		Free(linkedListGrid);
		Free(denGrid);

	}


//	DEBUGPRINT("C70 has peak id= %d  den= %g before merging\n", core[70].peak, core[70].density);

	int MergingPeak(SimpleBasicParticleType *, int, Coretype *, int, int);
	if(numcore >10) numcore = MergingPeak(bp,np, core, numcore,0);
	core = *Core = Realloc(*Core, sizeof(Coretype)*numcore);

//	DEBUGPRINT("C70 has peak id= %d  den= %g after merging\n", core[70].peak, core[70].density);



	DEBUGPRINT("The number of cores : %d and after MergingPeak\n", numcore);

    p = (particle *) Malloc(sizeof(particle)*np,PPTR(p));
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<np;i++){
        p[i].x = bp[i].x;
        p[i].y = bp[i].y;
        p[i].z = bp[i].z;
    }
    ptl = (TPtlStruct *)Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<np;i++){ 
		ptl[i].type = TYPE_PTL; 
		ptl[i].x = bp[i].x; 
		ptl[i].y = bp[i].y; 
		ptl[i].z = bp[i].z; 
		ptl[i].sibling = ptl+(i+1); 
		ptl[i].mass = 1;
    }
    ptl[np-1].sibling = NULL;
	nnode = MAX(np/2,65*10000);
    TREE = (TStruct *) Malloc(sizeof(TStruct)*nnode,PPTR(TREE));

//  all to all
	if(1){
		if(nnode > 65*10000) recursiveflag = PTHREAD;
		else recursiveflag = RECURSIVE;
		recursiveflag = SERIALIZED;
	    new_Make_Tree_Near(TREE,nnode, ptl,np,recursiveflag);
	}
	else { 
		Box box;
		Box findbox(TPtlStruct *, int);
		box = findbox(ptl, np);
		DEBUGPRINT("box size is %g %g %g width= %g\n", box.x,box.y,box.z,box.width);
		old_Make_Tree_Near(TREE, ptl,np,box);
	}
    DEBUGPRINT("after All_Make_Tree_near for np = %d\n",np);
	
	if(np < 100000)
	{ // an old and slow version
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k) schedule(guided)
#endif
	    for(i=0;i<np;i++){ 
            long long tmpindx[NUMNEIGHBOR];
            float mass[NUMNEIGHBOR];
            dptype dist2[NUMNEIGHBOR];
            float neardist;
			int res = Find_Near(p+i,Numnear,TREE,ptl,&neardist,tmpindx,
                    dist2, mass);
            if(res != Numnear){
                DEBUGPRINT("error occurred %d %d : %d\n", res, Numnear, np);
                exit(99);
            }
            k = i*Numnear;
            for(j=0;j<res;j++){
                nearindex[k+j] = tmpindx[j];
            }
		}
	}
	else { // a hacked version precopying the neighbor from the linked list
		LinkedListGrid *linkedListGrid = (LinkedListGrid*)Malloc(sizeof(LinkedListGrid)
				*ncells, PPTR(linkedListGrid));
		for(i=0;i<ncells;i++) {
			linkedListGrid[i].bp = NULL;
			linkedListGrid[i].np = 0;
		}
		for(i=0;i<np;i++){
			long long ir = rint((bp[i].x-xmin)/cellsize);
			long long jr = rint((bp[i].y-ymin)/cellsize);
			long long kr = rint((bp[i].z-zmin)/cellsize);
			long long ioff = ir+mx*(jr+ny*kr);
			SimpleBasicParticleType *tmp = linkedListGrid[ioff].bp;
			linkedListGrid[ioff].bp = bp+i;
			bp[i].bp = tmp;
			linkedListGrid[ioff].np ++;
		}
		int nthreads=1;
#ifdef _OPENMP
#pragma omp parallel
		{
			int it = Omp_get_thread_num();
			if(it ==0) nthreads = Omp_get_num_threads();
		}
#endif
		Neighbor *tneighbor = (Neighbor*)Malloc(sizeof(Neighbor)*nthreads*(NUMNEIGHBOR),PPTR(tneighbor));
		int npneigh=0;
		long long ioff;
		int nchunk = ncells/nthreads/mx;
		nchunk = MAX(1, nchunk);
		nchunk = MIN(1, nchunk);

		DEBUGPRINT("nthreads= %d nchunk = %d\n",nthreads, nchunk);

#ifdef _OPENMP
#pragma omp parallel for private(ioff,i,j,k) firstprivate(npneigh) schedule(dynamic,nchunk) num_threads(nthreads)
#endif
		for(ioff=0;ioff<ncells;ioff++){
			int threadID = Omp_get_thread_num();
			Neighbor *neighbor = tneighbor + threadID*(NUMNEIGHBOR);
			if(linkedListGrid[ioff].np>0 && linkedListGrid[ioff].np < DEEPSIZE){
				SimpleBasicParticleType *tmp = linkedListGrid[ioff].bp;
				/*
				// arrange templete from the linked list
				int inext = 0; // last neighbor index not-overlapped with target particles of current linkedlist cell.
				while(tmp) {
					long long indx = tmp-bp;
					int ialready = -1;
					for(j=0;j<npneigh;j++){
						if(neighbor[j].indx == indx){
							ialready = j;
							break;
						}
					}
					if(ialready ==-1) {
						neighbor[inext].indx = indx;
					}
					else {
						neighbor[inext].indx = neighbor[ialready].indx;
						neighbor[ialready].indx = neighbor[npneigh-1].indx;
						npneigh --;
					}
					inext ++;
					if(inext == Numnear) break;
					tmp = tmp->bp;
				}
				npneigh = MAX(inext,npneigh);
				*/
				// do the k-nearest finding
				tmp = linkedListGrid[ioff].bp;
				while(tmp){
					long long ibp = tmp-bp;
					particle pp;
					pp.x = tmp->x; 
					pp.y = tmp->y; 
					pp.z = tmp->z;
					for(i=0;i<npneigh;i++){ // update neighbor.dist2 for tmp
						long long indx = neighbor[i].indx;
						dptype tmpx = pp.x-bp[indx].x;
						dptype tmpy = pp.y-bp[indx].y;
						dptype tmpz = pp.z-bp[indx].z;
						dptype dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
						neighbor[i].dist2 = dist2;
					}
					if(npneigh) qsort(neighbor, npneigh, sizeof(Neighbor), sortneighdist2);


        			npneigh = Find_Near_PreCopied(&pp,Numnear,TREE,ptl, neighbor, npneigh);
					if(npneigh != Numnear){ 
						DEBUGPRINT("error occurred %d %d : %d && ioff= %lld\n", 
								npneigh, Numnear, np, ioff); 
						exit(99); 
					} 
					k = ibp*Numnear;
					for(j=0;j<npneigh;j++){ 
						nearindex[k+j] = neighbor[j].indx; 
					}
					tmp = tmp->bp;
				}
			}
		}

		DEBUGPRINT("shallow linked list part is done\n");


		for(ioff=0;ioff<ncells;ioff++){
			if(linkedListGrid[ioff].np>=DEEPSIZE) {// for a deeper linked list. 
											// This value should be larger than Num_neighbor.
				{
					SimpleBasicParticleType *tmp = linkedListGrid[ioff].bp;
					// arrange templete from the linked list
					int inext = 0;
					while(tmp && inext < Numnear) {
						long long indx = tmp-bp;
						tneighbor[inext].indx = indx;
						inext ++;
						tmp = tmp->bp;
					}
					npneigh = inext;
					for(i=1;i<nthreads;i++){
						for(j=0;j<npneigh;j++){
							tneighbor[i*NUMNEIGHBOR+j] = tneighbor[j];
						}
					}
				}
				// do the k-nearest finding
#ifdef _OPENMP
#pragma omp parallel firstprivate(npneigh)  private(i,j,k)
#endif
				{
					int threadID = Omp_get_thread_num();
					Neighbor *neighbor = tneighbor + threadID*(NUMNEIGHBOR);
					int icount=0;
					SimpleBasicParticleType *tmp = linkedListGrid[ioff].bp;
					while(tmp){
						if(icount%nthreads == threadID) { // allocate the job using id
							long long ibp = tmp-bp;
							particle pp;
							pp.x = tmp->x;
							pp.y = tmp->y;
							pp.z = tmp->z;
							for(i=0;i<npneigh;i++){ // update neighbor.dist2 for tmp
								long long indx = neighbor[i].indx;
								dptype tmpx = pp.x-bp[indx].x;
								dptype tmpy = pp.y-bp[indx].y;
								dptype tmpz = pp.z-bp[indx].z;
								dptype dist2 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
								neighbor[i].dist2 = dist2;
							}
							if(npneigh) qsort(neighbor, npneigh, sizeof(Neighbor), sortneighdist2);
		
		
   		   	  				npneigh = Find_Near_PreCopied(&pp,Numnear,TREE,ptl, neighbor, npneigh);
							if(npneigh != Numnear){ 
								DEBUGPRINT("error occurred %d %d : %d && ioff= %lld\n", 
										npneigh, Numnear, np, ioff); 
								exit(99); 
							} 
							k = ibp*Numnear;
							for(j=0;j<npneigh;j++){ 
								nearindex[k+j] = neighbor[j].indx; 
							}
						}
						icount ++;
						tmp = tmp->bp;
					}
				}
			}
		}




		Free(tneighbor);
		Free(linkedListGrid);
    }

	Free(TREE);
	Free(ptl);
	Free(p);

	DEBUGPRINT("making nearbor list is completed\n");


	*NumCore = numcore;

	return;
}

