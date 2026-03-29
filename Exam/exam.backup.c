#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>

#ifdef _OPENMP
#include<omp.h>
#endif

#include "eunha.h"

#include "voro.h"
#include "voro_eunha.h"
#include "nnost.h"
#include "gnnost.h"
#include "exam.h"
#include "exam2d.h"

double kh_w2Measure2D(SimParameters *simpar, double dmean, double pressure,double rho){
#define p0_w2 2.5
#define rho0_w2 1.
    double res = GAS_Kappa(simpar)*dmean*pow(pressure/p0_w2*rho0_w2/rho,0.5);
    double res2 = res*res;
    return res2;
#undef p0_w2
#undef rho0_w2
}

double kp_w2Measure2D(SimParameters *simpar, double dmean, double pressure,double rho){
#define p0_w2 2.5
#define rho0_w2 1.
    double res = GAS_Kappa(simpar)*dmean*pow(pressure/p0_w2*rho0_w2/rho,0.5);
    double res2 = res*res;
    return res2;
#undef p0_w2
#undef rho0_w2
}

double rt_w2Measure2D(SimParameters *simpar, double dmean, double pressure,double rho){
#define p0_w2 2.5
#define rho0_w2 1.
    double res = GAS_Kappa(simpar)*dmean*pow(pressure/p0_w2*rho0_w2/rho,0.5);
    double res2 = res*res;
    return res2;
#undef p0_w2
#undef rho0_w2
}

double gl2d_w2Measure(SimParameters *simpar, double dmean, double pressure,double rho){
#define p0_w2 1.
#define rho0_w2 1.
    double res = GAS_Kappa(simpar)*dmean*pow(pressure/p0_w2*rho0_w2/rho,0.5);
    double res2 = res*res;
    return res2;
#undef p0_w2
#undef rho0_w2
}




void ex2d_findCentroid(TStruct *thisNode){
    TPtlStruct *nodeparticles = (TPtlStruct*)thisNode->daughter;
    thisNode->monox = thisNode->monoy = 0;
    thisNode->mass = 0;
    TPtlStruct *p2ptl;
    for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling){
		postype mass = p2ptl->mass;
        thisNode->mass += mass;
        thisNode->monox += mass* p2ptl->x;
        thisNode->monoy += mass* p2ptl->y;
    }
    thisNode->monox = thisNode->monox/thisNode->mass;
    thisNode->monoy = thisNode->monoy/thisNode->mass;
}

void ex3d_findCentroid(TStruct *thisNode){
    TPtlStruct *nodeparticles = (TPtlStruct*)thisNode->daughter;
    thisNode->monox = thisNode->monoy = thisNode->monoz = 0;
    thisNode->mass = 0;
    TPtlStruct *p2ptl;
    for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling){
		postype mass = p2ptl->mass;
        thisNode->mass += mass;
        thisNode->monox += mass * p2ptl->x;
        thisNode->monoy += mass * p2ptl->y;
        thisNode->monoz += mass * p2ptl->z;
    }
    thisNode->monox = thisNode->monox/thisNode->mass;
    thisNode->monoy = thisNode->monoy/thisNode->mass;
    thisNode->monoz = thisNode->monoz/thisNode->mass;
}


void ex2d_findCellSize(TStruct *thisNode){
    TPtlStruct *nodeparticles = (TPtlStruct *) thisNode->daughter;
    PosType distmax = 1.e20;
    TPtlStruct *p2ptl;
    for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling){
        PosType tmpx = p2ptl->x - thisNode->monox;
        PosType tmpy = p2ptl->y - thisNode->monoy;
        PosType tmpdist2 = tmpx*tmpx + tmpy*tmpy;
        distmax = MAX(distmax, tmpdist2);
    }
    thisNode->nodesize = sqrt(distmax);
}

void ex3d_findCellSize(TStruct *thisNode){
    TPtlStruct *nodeparticles = (TPtlStruct *) thisNode->daughter;
    PosType distmax = 1.e20;
    TPtlStruct *p2ptl;
    for(p2ptl=nodeparticles;p2ptl;p2ptl=p2ptl->sibling){
        PosType tmpx = p2ptl->x - thisNode->monox;
        PosType tmpy = p2ptl->y - thisNode->monoy;
        PosType tmpz = p2ptl->z - thisNode->monoz;
        PosType tmpdist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
        distmax = MAX(distmax, tmpdist2);
    }
    thisNode->nodesize = sqrt(distmax);
}



void ex2d_idivision(TStruct *thisNode, TStruct *tmpNode){
    TPtlStruct *nodeparticles = thisNode->daughter;
    TPtlStruct *p2ptl;
    for(p2ptl=nodeparticles;p2ptl;) {
        int mx = (p2ptl->x>thisNode->monox ? 1:0);
        int my = (p2ptl->y>thisNode->monoy ? 1:0);
        int mnode = mx + 2*my;
        tmpNode[mnode].Nparticle ++;
        TPtlStruct *tmpptr = (TPtlStruct *)tmpNode[mnode].daughter;
        TPtlStruct *tmpptr2 = p2ptl->sibling;
        tmpNode[mnode].daughter = (void *)p2ptl;
        p2ptl->sibling = tmpptr;
        p2ptl = tmpptr2;
    }
}

void ex3d_idivision(TStruct *thisNode, TStruct *tmpNode){
    TPtlStruct *nodeparticles = thisNode->daughter;
    TPtlStruct *p2ptl;
    for(p2ptl=nodeparticles;p2ptl;) {
        int mx = (p2ptl->x>thisNode->monox ? 1:0);
        int my = (p2ptl->y>thisNode->monoy ? 1:0);
        int mz = (p2ptl->z>thisNode->monoz ? 1:0);
        int mnode = mx + 2*(my + 2*mz);
        tmpNode[mnode].Nparticle ++;
        TPtlStruct *tmpptr = (TPtlStruct *)tmpNode[mnode].daughter;
        TPtlStruct *tmpptr2 = p2ptl->sibling;
        tmpNode[mnode].daughter = (void *)p2ptl;
        p2ptl->sibling = tmpptr;
        p2ptl = tmpptr2;
    }
}

int nearest2dOpen(void *p, TStruct *tree, postype mindist){
	TPtlStruct *point = (TPtlStruct*)p;
	postype tmpx = point->x - tree->monox;
	postype tmpy = point->y - tree->monoy;
	postype dist2 = tmpx*tmpx + tmpy * tmpy;
	postype dist = sqrt(dist2);
	if( (dist- tree->nodesize) > mindist) return NO;
	else return YES;
}

int nearest3dOpen(void *p, TStruct *tree, postype mindist){
    TPtlStruct *point = (TPtlStruct*)p;
    postype tmpx = point->x - tree->monox;
    postype tmpy = point->y - tree->monoy;
    postype tmpz = point->z - tree->monoz;
    postype dist2 = tmpx*tmpx + tmpy * tmpy + tmpz *tmpz;
    postype dist = sqrt(dist2);
    if( dist- tree->nodesize > mindist) return NO;
    else return YES;
}


postype ex2d_dist(void *p, void *pp){
	TPtlStruct *point = (TPtlStruct*)p;
	TPtlStruct *neigh = (TPtlStruct *)pp;
	postype tmpx = point->x - neigh->x;
	postype tmpy = point->y - neigh->y;
	postype d2 = tmpx*tmpx + tmpy*tmpy;
	return d2;
}

postype ex3d_dist(void *p, void *pp){
    TPtlStruct *point = (TPtlStruct*)p;
    TPtlStruct *neigh = (TPtlStruct *)pp;
    postype tmpx = point->x - neigh->x;
    postype tmpy = point->y - neigh->y;
    postype tmpz = point->z - neigh->z;
    postype d2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
    return d2;
}



void det2d_dpq(
		SimParameters *simpar, 
		void (*paddingAllTreeParticles)(SimParameters *, postype) 
		){
	treevoroparticletype *bp = VORO_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;

	{
		postype cellsize = HydroGridSize(simpar);
    	paddingAllTreeParticles(simpar, cellsize);
	}
	int mp = VORO_NPAD(simpar);

	TPtlStruct *ptl = (TPtlStruct*)my_malloc(sizeof(TPtlStruct)*(np+mp));
	int nnode = np+mp;
	TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode);

	postype xmin,ymin,xmax,ymax;
	xmin = ymin = 1.e20;
	xmax = ymax = -1.e20;
	// this is for the local-domain particles.
	for(i=0;i<np;i++){
		ptl[i].x = bp[i].x;
		ptl[i].y = bp[i].y;
		ptl[i].mass = bp[i].mass; /* in this case mass */
		ptl[i].sibling = ptl + i+1;
		ptl[i].type = TYPE_PTL;
		ptl[i].indx = i;
		xmin = MIN(xmin, bp[i].x);
		ymin = MIN(ymin, bp[i].y);
		xmax = MAX(xmax, bp[i].x);
		ymax = MAX(ymax, bp[i].y);
	}

	// this is for the boundary ghost particles.
	bp = VORO_TBPP(simpar);
	for(i=0;i<mp;i++){
		ptl[np+i].x = bp[i].x;
		ptl[np+i].y = bp[i].y;
		ptl[np+i].mass = bp[i].mass; /* in this case mass */
		ptl[np+i].sibling = ptl + np+i+1;
		ptl[np+i].type = TYPE_PTL;
		ptl[np+i].indx = np+i;
	}
	// terminate the particle sibling pointer
	ptl[np+mp-1].sibling = NULL;
	// tree node has a pointer (daughter) linked to the first particle.
	tree->daughter = ptl;
	// total number of particles linked to the tree node is (np+mp).
	tree->Nparticle = np+mp;

	Make_GNN_Tree(tree, nnode, ex2d_idivision, ex2d_findCentroid,
			ex2d_findCellSize, RECURSIVE);

	// back to the local-domain particles
	bp = VORO_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<np;i++){ 
		bp[i].w2ceil = find_GNearest(ptl+i,  tree, nearest2dOpen, ex2d_dist);
		bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
	}

	my_free(ptl);
	my_free(tree);
	my_free(VORO_TBPP(simpar));
	if(0){
		postype mindpq,maxdpq;
		mindpq = 1.e20;
		maxdpq = -1.e20;
		for(i=0;i<np;i++){
			mindpq = MIN(mindpq, bp[i].w2ceil);
			maxdpq = MAX(maxdpq, bp[i].w2ceil);
		}
		DEBUGPRINT("P%d passed the finding of dpq: min/max= %g %g\n",MYID(simpar), mindpq,maxdpq);
	}
}

void det3d_dpq(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype)
		){
    treevoroparticletype *bp = VORO_TBP(simpar);
    int np = VORO_NP(simpar);

	{
		postype cellsize = HydroGridSize(simpar);
    	paddingAllTreeParticles(simpar, cellsize);

	}
	int mp = VORO_NPAD(simpar);

    TPtlStruct *ptl = (TPtlStruct*)my_malloc(sizeof(TPtlStruct)*(np+mp));
    int nnode = np+mp;
    TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode);
    int i;

    for(i=0;i<np;i++){
        ptl[i].x = bp[i].x;
        ptl[i].y = bp[i].y;
        ptl[i].z = bp[i].z;
		ptl[i].mass = bp[i].mass;
        ptl[i].sibling = ptl + i+1;
        ptl[i].type = TYPE_PTL;
        ptl[i].indx = i;
    }

	bp = VORO_TBPP(simpar);
	for(i=0;i<mp;i++){
		ptl[np+i].x = bp[i].x;
		ptl[np+i].y = bp[i].y;
		ptl[np+i].z = bp[i].z;
		ptl[np+i].mass = bp[i].mass; /* in this case mass */
		ptl[np+i].sibling = ptl + np+i+1;
		ptl[np+i].type = TYPE_PTL;
		ptl[np+i].indx = np+i;
	}
	ptl[np+mp-1].sibling = NULL;
    tree->daughter = ptl;
    tree->Nparticle = np+mp;

    Make_GNN_Tree(tree, nnode, ex3d_idivision, ex3d_findCentroid,
            ex3d_findCellSize, RECURSIVE);

	// back to the local-domain particles
	bp = VORO_TBP(simpar);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(i=0;i<np;i++){
        bp[i].w2ceil = find_GNearest(ptl+i,  tree, nearest3dOpen, ex3d_dist);
    }
    my_free(ptl);
    my_free(tree);
	my_free(VORO_TBPP(simpar));
}
void det2d_dpqRK4(
		SimParameters *simpar, 
		void (*paddingAllTreeParticles)(SimParameters *, postype) ){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;

	int mp;
	{
		postype cellsize = HydroGridSize(simpar);
    	paddingAllTreeParticles(simpar, cellsize);
		mp = VORO_NPAD(simpar);
	}

	TPtlStruct *ptl = (TPtlStruct*)my_malloc(sizeof(TPtlStruct)*(np+mp));
	int nnode = np+mp;
	TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode);


	postype xmin,ymin,xmax,ymax;
	xmin = ymin = 1.e20;
	xmax = ymax = -1.e20;
	// this is for the local-domain particles.
	for(i=0;i<np;i++){
		ptl[i].x = bp[i].x;
		ptl[i].y = bp[i].y;
		ptl[i].mass = bp[i].mass; /* in this case mass */
		ptl[i].sibling = ptl + i+1;
		ptl[i].type = TYPE_PTL;
		ptl[i].indx = i;
		xmin = MIN(xmin, bp[i].x);
		ymin = MIN(ymin, bp[i].y);
		xmax = MAX(xmax, bp[i].x);
		ymax = MAX(ymax, bp[i].y);
	}

	// this is for the boundary ghost particles.
	bp = VORORK4_TBPP(simpar);
	for(i=0;i<mp;i++){
		ptl[np+i].x = bp[i].x;
		ptl[np+i].y = bp[i].y;
		ptl[np+i].mass = bp[i].mass; /* in this case mass */
		ptl[np+i].sibling = ptl + np+i+1;
		ptl[np+i].type = TYPE_PTL;
		ptl[np+i].indx = np+i;
	}
	// terminate the particle sibling pointer
	ptl[np+mp-1].sibling = NULL;
	// tree node has a pointer (daughter) linked to the first particle.
	tree->daughter = ptl;
	// total number of particles linked to the tree node is (np+mp).
	tree->Nparticle = np+mp;

	Make_GNN_Tree(tree, nnode, ex2d_idivision, ex2d_findCentroid,
			ex2d_findCellSize, RECURSIVE);

	// back to the local-domain particles
	bp = VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<np;i++){ 
		bp[i].w2ceil = find_GNearest(ptl+i,  tree, nearest2dOpen, ex2d_dist);
		bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
	}

	my_free(ptl);
	my_free(tree);
	my_free(VORORK4_TBPP(simpar));
	if(0){
		postype mindpq,maxdpq;
		mindpq = 1.e20;
		maxdpq = -1.e20;
		for(i=0;i<np;i++){
			mindpq = MIN(mindpq, bp[i].w2ceil);
			maxdpq = MAX(maxdpq, bp[i].w2ceil);
		}
		DEBUGPRINT("P%d passed the finding of dpq: min/max= %g %g\n",MYID(simpar), mindpq,maxdpq);
	}
}

void det3d_dpqRK4(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype)
		){
    treevorork4particletype *bp = VORORK4_TBP(simpar);
    int np = VORO_NP(simpar);

	{
		postype cellsize = HydroGridSize(simpar);
    	paddingAllTreeParticles(simpar, cellsize);

	}
	int mp = VORO_NPAD(simpar);

    TPtlStruct *ptl = (TPtlStruct*)my_malloc(sizeof(TPtlStruct)*(np+mp));
    int nnode = np+mp;
    TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode);
    int i;

    for(i=0;i<np;i++){
        ptl[i].x = bp[i].x;
        ptl[i].y = bp[i].y;
        ptl[i].z = bp[i].z;
		ptl[i].mass = bp[i].mass;
        ptl[i].sibling = ptl + i+1;
        ptl[i].type = TYPE_PTL;
        ptl[i].indx = i;
    }

	bp = VORORK4_TBPP(simpar);
	for(i=0;i<mp;i++){
		ptl[np+i].x = bp[i].x;
		ptl[np+i].y = bp[i].y;
		ptl[np+i].z = bp[i].z;
		ptl[np+i].mass = bp[i].mass; /* in this case mass */
		ptl[np+i].sibling = ptl + np+i+1;
		ptl[np+i].type = TYPE_PTL;
		ptl[np+i].indx = np+i;
	}
	ptl[np+mp-1].sibling = NULL;
    tree->daughter = ptl;
    tree->Nparticle = np+mp;

    Make_GNN_Tree(tree, nnode, ex3d_idivision, ex3d_findCentroid,
            ex3d_findCellSize, RECURSIVE);

	// back to the local-domain particles
	bp = VORORK4_TBP(simpar);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(i=0;i<np;i++){
        bp[i].w2ceil = find_GNearest(ptl+i,  tree, nearest3dOpen, ex3d_dist);
    }
    my_free(ptl);
    my_free(tree);
	my_free(VORORK4_TBPP(simpar));
}

Voro2D_point *searchCellRk4Neighbors2D(
		SimParameters *simpar,
        int ix, 
		int iy, 
		int *nneigh
		){
    int i,j,k;
    int np,mp;
    Voro2D_point *neigh;
    int mx,my;
    mx = BASICCELL_MX(simpar);
    my = BASICCELL_MY(simpar);
    CellType *cells = VORO_BASICCELL(simpar);
    float Lx = KH_SIMBOX(simpar).x.max - KH_SIMBOX(simpar).x.min;
    float Ly = KH_SIMBOX(simpar).y.max - KH_SIMBOX(simpar).y.min;

    np = 0;
    for(j=iy-1;j<=iy+1;j++){
		if(j<0 || j >= my) continue;
        for(i=ix-1;i<=ix+1;i++){
			if(i<0 || i >= mx) continue;
            np += cells[i+mx*j].nmem;
        }
    }
	// twice the size for the boundary conditions
    neigh = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*np);
    np = 0;
    for(j=iy-1;j<=iy+1;j++){
		if( j<0 || j >= my) continue;
        for(i=ix-1;i<=ix+1;i++){
			if(i<0 || i >= mx) continue;
            struct linkedlisttype *tmp = cells[i+mx*j].link;
            while(tmp){
                treevorork4particletype *tt = (treevorork4particletype*)tmp;
                neigh[np].x = tt->x;
                neigh[np].y = tt->y;
                neigh[np].indx = PINDX(tt);
                neigh[np].csound = tt->csound;
                neigh[np].pressure = tt->pressure;
                neigh[np].w2 = tt->w2;
                neigh[np].bp = (void*)tt;
                neigh[np].u4if = tt->u4if;
                np++;
                tmp = tmp->next;
            }

        }
    }
    *nneigh = np;
    return neigh;
}

treevorork4particletype *findCellRk4BP2D(SimParameters *simpar, int ix, int iy, int *mp){
    int i,j,k;
    int np;
    treevorork4particletype *res;
    int mx,my;
    mx = BASICCELL_MX(simpar);
    my = BASICCELL_MY(simpar);
    CellType *cells = VORO_BASICCELL(simpar);

    np = cells[ix+mx*iy].nmem;
	// twice the size for the boundary conditions
    res = (treevorork4particletype*)my_malloc(sizeof(treevorork4particletype)*np);
    np = 0;
    int ipixel = ix+mx*iy;
    struct linkedlisttype *tmp = cells[ipixel].link;
    while(tmp){
        if(!IS_FLAG_ON(tmp,BoundaryGhostflag)){
            treevorork4particletype *tt = (treevorork4particletype*)tmp;
            {
                res[np] = *tt;
                res[np].bp = tt;
                np++;
            }
        }
        tmp = tmp->next;
    }
    *mp = np;
    return res;
}

void mkLinkedList2D_rt(
		SimParameters *simpar, 
		postype cellsize, 
		postype xmin, 
		postype ymin,
		postype xmax, 
		postype ymax, 
		void (*paddingTreeAllParticle)(SimParameters *, postype)
		){ 
	treevorork4particletype *bp = VORORK4_TBP(simpar); 
	int np = VORO_NP(simpar);
    int i,j,k;
    int mx,my;
	postype Gamma = GAS_GAMMA(simpar);
	postype Ly = SIMBOX(simpar).y.max;

    mx = BASICCELL_MX(simpar);
    my = BASICCELL_MY(simpar);

	CellType *cells = VORO_BASICCELL(simpar); 

	for(i=0;i<mx*my;i++) {
        cells[i].link = NULL;
        cells[i].nmem = 0;
    }
	for(i=0;i<VORO_NP(simpar);i++){
        int ix,iy;
        ix = ((bp[i].x-xmin)/cellsize);
        iy = ((bp[i].y-ymin)/cellsize);
        int index = ix+mx*iy;
//		if(PINDX(bp+i) ==0)DEBUGPRINT("P%d has p%ld has pressure= %g vol= %g\n", MYID(simpar), PINDX(bp+i), bp[i].pressure, bp[i].volume);
        if(index >=mx*my || index <0){
            DEBUGPRINT("error detected: %d : %d %d : %g %g %g %g\n",
                    index, mx,my, bp[i].x, bp[i].y, xmin,ymin);
        }
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)(bp+i);
        cells[index].nmem ++;
        bp[i].next = tmp;
    }
	// now for the particle padding from other ranks
    paddingTreeAllParticle(simpar, cellsize);

	{
	    treevorork4particletype *bpp = VORORK4_TBPP(simpar); 
		for(i=0;i<VORO_NPAD(simpar);i++){ 
			int ix,iy; 
			ix = ((bpp[i].x-xmin)/cellsize); 
			iy = ((bpp[i].y-ymin)/cellsize); 
			int index = ix+mx*iy; 
			struct linkedlisttype *tmp = cells[index].link; 
			cells[index].link = (struct linkedlisttype*)(bpp+i); 
			cells[index].nmem ++; 
			bpp[i].next = tmp; 
		}
	}



	{
		postype rho1 = RT_DEN1(simpar);
		postype rho2 = RT_DEN2(simpar);
		postype deltay = RT_Deltay(simpar);
		postype ycen = 0.5*Ly;
		postype getDen(postype , postype , postype , postype ,postype );
		int ny = NY(simpar);
		treevorork4particletype *bpp = VORORK4_TBPP(simpar);
		postype dy = Ly/(postype)ny;
//		DEBUGPRINT("P%d has %g %g deltay= %g pcenter = %g\n", MYID(simpar), rho1, rho2, deltay, RT_Phalf(simpar));
		// top pressure
		postype pup = RT_Phalf(simpar);
		for(j=ny/2+1;j<ny+1;j++){
			postype y = (postype)(j+0.5)*Ly/(postype)ny; 
			postype rho = getDen(rho1, rho2, deltay, ycen, y); 
			postype dr = dy;
			pup += GAS_ACCY(simpar) * rho * dr;
		}
		for(i=0;i<VORO_NPAD(simpar);i++){
			postype getDen(postype , postype , postype , postype ,postype );
			if(bpp[i].y >= Ly) {
				bpp[i].den = getDen(rho1, rho2, deltay, ycen, bpp[i].y);
				postype rho = getDen(rho1, rho2, deltay, ycen, 0.5*(bpp[i].y+Ly));
				postype dr = (bpp[i].y- Ly);
				bpp[i].pressure = pup + GAS_ACCY(simpar)* rho * dr;
				bpp[i].ie = bpp[i].pressure*bpp[i].volume/(Gamma-1);
        		bpp[i].csound = sqrt(Gamma*bpp[i].pressure/bpp[i].den);
		        if(GAS_Kappa(simpar) <0) { 
					bpp[i].w2 = -GAS_Kappa(simpar);
		        }
	       		else if (GAS_Kappa(simpar) >0){ 
					bpp[i].w2 = getw2forHydroParticle(simpar,(bpp+i)); 
        		}
				bpp[i].w2old = bpp[i].w2; // static boundary conditions.
				bpp[i].vx = bpp[i].vy = 0;
			}
		}
		// bottom pressure
		postype pbottom = RT_Phalf(simpar);
		for(j=ny/2-1;j>=0;j--){
			postype y = (postype)(j+0.5)*Ly/(postype)ny;
			postype rho = getDen(rho1,rho2, deltay,ycen,y);
			postype dr = -dy;
			pbottom += GAS_ACCY(simpar)*rho*dr;
		}
		for(i=0;i<VORO_NPAD(simpar);i++){
			postype getDen(postype , postype , postype , postype ,postype );
			if(bpp[i].y < 0) {
				bpp[i].den = getDen(rho1, rho2, deltay, ycen, bpp[i].y);
				postype rho = getDen(rho1, rho2, deltay, ycen, bpp[i].y*0.5);
				postype dr = bpp[i].y;
				bpp[i].pressure =  pbottom +  GAS_ACCY(simpar)* rho * dr;
				bpp[i].ie = bpp[i].pressure*bpp[i].volume/(Gamma-1);
        		bpp[i].csound = sqrt(Gamma*bpp[i].pressure/bpp[i].den);
		        if(GAS_Kappa(simpar) <0) { 
					bpp[i].w2 = -GAS_Kappa(simpar);
		        }
	       		else if (GAS_Kappa(simpar) >0){ 
					bpp[i].w2 = getw2forHydroParticle(simpar,(bpp+i)); 
        		}
				bpp[i].w2old = bpp[i].w2;// static boundary conditions
				bpp[i].vx = bpp[i].vy = 0;
			}
		}
	}
}

void old_mkLinkedList2D_rt(
		SimParameters *simpar, 
		postype cellsize, 
		postype xmin, 
		postype ymin,
		postype xmax, 
		postype ymax, 
		void (*paddingTreeAllParticle)(SimParameters *, postype)
		){ 
	treevorork4particletype *bp = VORORK4_TBP(simpar); 
	int np = VORO_NP(simpar);
    int i,j,k;
    int mx,my;
	postype Gamma = GAS_GAMMA(simpar);
	postype Ly = SIMBOX(simpar).y.max;
	postype rho1 = RT_DEN1(simpar), rho2 = RT_DEN2(simpar);
	postype deltay = RT_Deltay(simpar);
	postype ycen = 0.5*Ly;

    mx = BASICCELL_MX(simpar);
    my = BASICCELL_MY(simpar);

	CellType *cells = VORO_BASICCELL(simpar); 

	for(i=0;i<mx*my;i++) {
        cells[i].link = NULL;
        cells[i].nmem = 0;
    }
	for(i=0;i<VORO_NP(simpar);i++){
        int ix,iy;
        ix = ((bp[i].x-xmin)/cellsize);
        iy = ((bp[i].y-ymin)/cellsize);
        int index = ix+mx*iy;
        if(index >=mx*my || index <0){
            DEBUGPRINT("error detected: %d : %d %d : %g %g %g %g\n",
                    index, mx,my, bp[i].x, bp[i].y, xmin,ymin);
        }
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)(bp+i);
        cells[index].nmem ++;
        bp[i].next = tmp;
    }
	// now for the particle padding from other ranks
    paddingTreeAllParticle(simpar, cellsize);


	int ny = NY(simpar);

	if(Ymax_HydroExam(simpar) >= (Ly - EPS)){
		treevorork4particletype *bpp = VORORK4_TBPP(simpar);
		int mpad = 0;
		for(i=0;i<VORO_NPAD(simpar);i++){
			if(bpp[i].y < Ly) {
				bpp[mpad] = bpp[i]; mpad ++;
			}
		}
		int nmirror = 0;
//	    pre-estimate the space for the mirror particles
		for(i=0;i<VORO_NP(simpar);i++){
			if((Ly-bp[i].y) < cellsize) nmirror++;
		}
		for(i=0;i<mpad;i++){
			if((Ly-bpp[i].y) < cellsize) nmirror++;
		}
		if(VORO_NPAD(simpar) != (mpad+nmirror)) {
			DEBUGPRINT("P%d is resizing the boundary ghost particles: mpad= %d nmirror= %d : org= %ld\n", 
					MYID(simpar), mpad,nmirror, VORO_NPAD(simpar));
			bpp = (VORORK4_TBPP(simpar) = (treevorork4particletype *)
					realloc(VORORK4_TBPP(simpar), sizeof(treevorork4particletype)*(mpad+nmirror)));
		}
		// starting from the empty slot
		treevorork4particletype *tmp = VORORK4_TBPP(simpar) + mpad; 
		for(i=0;i<VORO_NP(simpar);i++){
			postype getDen(postype , postype , postype , postype ,postype );
			if((Ly-bp[i].y) < cellsize) {
				*tmp = bp[i];
				postype dy = 2*(Ly-bp[i].y);
				tmp->u4if.indx = MAX_INDEX;
				tmp->y += dy;
				tmp->vy = -bp[i].vy;
				tmp->den = getDen(rho1, rho2, deltay, ycen, tmp->y);
				tmp->pressure += GAS_ACCY(simpar)* getDen(rho1,rho2,deltay,ycen,Ly) * dy;
				tmp->ie = tmp->pressure*tmp->volume/(Gamma-1);
				CLEAR_FLAG(tmp); 
				SET_FLAG(tmp,BoundaryGhostflag);
				SET_FLAG(tmp,VOROflag);
				if(PINDX(tmp)==130822) DEBUGPRINT("p%lu is mirrored xy= %g %g\n", PINDX(tmp), tmp->x,tmp->y);
				tmp++;
			}
		}
		for(i=0;i<mpad;i++){
			postype getDen(postype , postype , postype , postype ,postype );
			if((Ly-bpp[i].y) < cellsize) {
				*tmp = bpp[i];
				tmp->u4if.indx = MAX_INDEX;
				postype dy = 2*(Ly-bpp[i].y);
				tmp->y += dy;
				tmp->vy = -bpp[i].vy;
				tmp->den = getDen(rho1, rho2, deltay, ycen, tmp->y);
				tmp->pressure += GAS_ACCY(simpar)* getDen(rho1,rho2,deltay,ycen,Ly) * dy;
				tmp->ie = tmp->pressure*tmp->volume/(Gamma-1);
				CLEAR_FLAG(tmp); 
				SET_FLAG(tmp,BoundaryGhostflag);
				SET_FLAG(tmp,VOROflag);
				tmp++;
			}
		}
		VORO_NPAD(simpar) = mpad + nmirror;
	}
	if(Ymin_HydroExam(simpar)< EPS)
	{
		treevorork4particletype *bpp = VORORK4_TBPP(simpar);
		int mpad = 0;
		for(i=0;i<VORO_NPAD(simpar);i++){
			if(bpp[i].y >= 0.) {
				bpp[mpad] = bpp[i]; mpad ++;
			}
		}
		int nmirror = 0;
		for(i=0;i<VORO_NP(simpar);i++){
			if(bp[i].y <= cellsize) nmirror++;
		}
		for(i=0;i<mpad;i++){
			if(bpp[i].y <= cellsize) nmirror++;
		}
		if(VORO_NPAD(simpar) != (mpad+nmirror)) {
			DEBUGPRINT("P%d is resizing the boundary ghost particles: mpad= %d nmirror= %d : org= %ld\n", 
					MYID(simpar), mpad,nmirror, VORO_NPAD(simpar));
			bpp = (VORORK4_TBPP(simpar) = (treevorork4particletype *) realloc(VORORK4_TBPP(simpar), 
						sizeof(treevorork4particletype)*(mpad+nmirror)));
		}
		// starting from the empty slot
		treevorork4particletype *tmp = VORORK4_TBPP(simpar) + mpad; 
		for(i=0;i<VORO_NP(simpar);i++){
			postype getDen(postype , postype , postype , postype ,postype );
			if(bp[i].y <= cellsize) {
				*tmp = bp[i];
				tmp->u4if.indx = MAX_INDEX;
				postype dy = -bp[i].y;
				tmp->y = dy;
				tmp->vy = -bp[i].vy;
				tmp->den = getDen(rho1, rho2, deltay, ycen, tmp->y);
				tmp->pressure -= GAS_ACCY(simpar)* getDen(rho1,rho2,deltay,ycen,0.) * 2*(-dy);
				tmp->ie = tmp->pressure*tmp->volume/(Gamma-1);
				CLEAR_FLAG(tmp); 
				SET_FLAG(tmp,BoundaryGhostflag);
				SET_FLAG(tmp,VOROflag);
				tmp++;
			}
		}
		for(i=0;i<mpad;i++){
			postype getDen(postype , postype , postype , postype ,postype );
			if(bpp[i].y <= cellsize) {
				*tmp = bpp[i];
				tmp->u4if.indx = MAX_INDEX;
				postype dy = -bpp[i].y;
				tmp->y = dy;
				tmp->vy = -bpp[i].vy;
				tmp->den = getDen(rho1, rho2, deltay, ycen, tmp->y);
				tmp->pressure -= GAS_ACCY(simpar)* getDen(rho1,rho2,deltay,ycen,0.) * 2*(-dy);
				tmp->ie = tmp->pressure*tmp->volume/(Gamma-1);
				CLEAR_FLAG(tmp); 
				SET_FLAG(tmp,BoundaryGhostflag);
				SET_FLAG(tmp,VOROflag);
				tmp++;
			}
		}
		VORO_NPAD(simpar) = mpad + nmirror;
	}
	{
	    treevorork4particletype *bpp = VORORK4_TBPP(simpar); 
		for(i=0;i<VORO_NPAD(simpar);i++){ 
			int ix,iy; 
			ix = ((bpp[i].x-xmin)/cellsize); 
			iy = ((bpp[i].y-ymin)/cellsize); 
			int index = ix+mx*iy; 
			struct linkedlisttype *tmp = cells[index].link; 
			cells[index].link = (struct linkedlisttype*)(bpp+i); 
			cells[index].nmem ++; 
			bpp[i].next = tmp; 
		}
	}
}
void mkLinkedList2D_oExam(
		SimParameters *simpar, 
		postype cellsize, 
		postype xmin, 
		postype ymin,
		postype xmax, 
		postype ymax, 
		void (*paddingTreeAllParticle)(SimParameters *, postype)
		){ 
	treevorork4particletype *bp = VORORK4_TBP(simpar); 
	int np = VORO_NP(simpar);
    int i,j,k;
    int mx,my;
	postype Gamma = GAS_GAMMA(simpar);
	postype Ly = SIMBOX(simpar).y.max;
	postype rho1 = RT_DEN1(simpar), rho2 = RT_DEN2(simpar);
	postype deltay = RT_Deltay(simpar);
	postype ycen = 0.5*Ly;

    mx = BASICCELL_MX(simpar);
    my = BASICCELL_MY(simpar);

	CellType *cells = VORO_BASICCELL(simpar); 

	for(i=0;i<mx*my;i++) {
        cells[i].link = NULL;
        cells[i].nmem = 0;
    }
	for(i=0;i<VORO_NP(simpar);i++){
        int ix,iy;
        ix = ((bp[i].x-xmin)/cellsize);
        iy = ((bp[i].y-ymin)/cellsize);
        int index = ix+mx*iy;
        if(index >=mx*my){
            DEBUGPRINT("error detected: %d : %d %d : %g %g %g %g\n",
                    index, mx,my, bp[i].x, bp[i].y, xmin,ymin);
        }
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)(bp+i);
        cells[index].nmem ++;
        bp[i].next = tmp;
    }
	// now for the particle padding from other ranks
    paddingTreeAllParticle(simpar, cellsize);
    bp = VORORK4_TBPP(simpar);
    for(i=0;i<VORO_NPAD(simpar);i++){
        int ix,iy;
        ix = ((bp[i].x-xmin)/cellsize);
        iy = ((bp[i].y-ymin)/cellsize);
        int index = ix+mx*iy;
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)(bp+i);
        cells[index].nmem ++;
        bp[i].next = tmp;
    }
}



void updateDenW2Pressure2D(
		SimParameters *simpar, 
		postype xmin, postype ymin, postype xmax, postype ymax,
		 postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int isave = -1;


	det2d_dpqRK4(simpar, paddingAllTreeParticles);
//	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed det2d_dpqRK4...()\n", MYID(simpar));
	int i;
	for(i=0;i<nbp;i++){
		bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
	}
	/*
	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	*/
	postype cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar);
    int mx, my;
    BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
    BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
//	DEBUGPRINT("P%d has mxy = %d %d cellsize= %g\n", MYID(simpar), mx,my,cellsize);


	CellType *cells = VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*mx*my);
	mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

    int iy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(iy=1;iy<my-1;iy++){
    	int mp=1000;
    	Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp);
    	postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh);
			int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i); 
                center.csound = p[i].csound;
                center.w2 = p[i].w2;


                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
				treevorork4particletype *ibp = p[i].bp;
                ibp->volume = Area2DPolygon(vorocorner, mp);
                ibp->den = ibp->mass/ibp->volume;
				/*
                ibp->csound = sqrt(Gamma*ibp->pressure/ibp->den);
                ibp->pressure = ibp->ie/ibp->volume*(Gamma-1);
				*/
            }
            if(np>0) my_free(p);
			if(nneigh>0) {
				my_free(neighbors); my_free(neighwork);
			}

        }
        my_free(vorocorner);
    }
	for(i=0;i<nbp;i++){
		bp[i].pressure = bp[i].ie/bp[i].volume*(Gamma-1);
		bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
	}
	// finalize mkLinkedList2D by my_free memory spaces (cell & boundary ghost particles)
	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
}

double getAccVoro2D_rt(SimParameters *simpar, postype xmin, postype ymin, 
		postype xmax, postype ymax,
		postype OrderOfAccuracy, postype Courant, postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int isave = -1;
    postype Dtime = 1.e10;


	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	postype cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar);
    int mx, my;
    BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
    BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	CellType *cells = VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my);
	mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);


	postype dtold = GAS_dtold(simpar);



	float alphavis = GAS_AlphaVis(simpar);
	float betavis = GAS_BetaVis(simpar);
	float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
	float eps2vis = GAS_EPSVIS(simpar)*GAS_EPSVIS(simpar);


//	DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);

    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=1;iy<my-1;iy++){
    	int mp=1000;
    	Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp);
    	postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh);
			int i;
            for(i=0;i<np;i++){

                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i); 
                center.csound = p[i].csound;
                center.w2 = p[i].w2;

				treevorork4particletype *ibp = p[i].bp;
				ibp->dt = 1.e10; // initialization of dt

                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
//                ibp->volume = Area2DPolygon(vorocorner, mp);
				// to use registers
				postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
				ibp_vx = ibp->vx;
				ibp_vy = ibp->vy;
				ibp_csound = ibp->csound;
				ibp_pressure = ibp->pressure;
				ibp_den = ibp->den;


				Voro2D_Corner *tmp,*tmp2;
				tmp = vorocorner;
				int j;
				double die, dte,dke,fx,fy;
				die = dte = dke = fx = fy = 0;
				do {
					if(tmp->upperrelated >=0)
					{
						// to use registers
						treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
						postype jbp_vx = jbp->vx;
						postype jbp_vy = jbp->vy;
						postype jbp_csound = jbp->csound;
                        postype jbp_pressure = jbp->pressure;
//						postype jbp_pressure = getNeighborPressure(jbp,ibp); // to consider the wallflag
						postype jbp_den = jbp->den;

						Voro2D_point line; 
	                    tmp2 = tmp->upperlink; 
						line.x = tmp2->x - tmp->x; 
						line.y = tmp2->y - tmp->y;
	                    Voro2D_point dS = voro2D_norm(&line); 
						postype facearea = Vec2DLength(tmp,tmp2); 
						dS.x = facearea*dS.x;
						dS.y = facearea*dS.y;

						/*
						postype pi = getPressure2D(
								ibp_pressure,jbp_pressure, tmp, neighwork, OrderOfAccuracy);
								*/
						postype pi = VoroRK4_Pressure2D(
								ibp_pressure,jbp_pressure, tmp, neighwork, simpar, RT);

						Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

						postype dramp = sqrt(Vec2DDotP(&dr, &dr));
						Voro2D_point er;
						er.x = dr.x/dramp;
						er.y = dr.y/dramp;

						Voro2D_point ui,ua,ub;
						/*
						ui.x = Half*(jbp_vx - ibp_vx); 
						ui.y = Half*(jbp_vy - ibp_vy);
						*/
						ui = get2dUpqradRk4(ibp,jbp,dtold);

						{
							postype rvel = Vec2DDotP(&er, &ui);
							if(rvel<0){
								postype mu = rvel/(dramp/etavis + eps2vis*etavis/dramp);
								postype meanden = 0.5*(ibp_den + jbp_den);
								postype meanCsound = 0.5*(ibp_csound + jbp_csound);
								pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
								if(center.indx==130560) DEBUGPRINT("p(%ld,%ld) has pij= %g rvel= %g : w(ibp)= %g %g ::: w(jbp) = %g %g :: vi= %g %g :: vj= %g %g :: jindx= %ld\n", 
									PINDX(ibp)%256,PINDX(ibp)/256, 
									pi, rvel, ibp->w2, ibp->w2old, jbp->w2, jbp->w2old,
									ibp->vx, ibp->vy, jbp->vx,jbp->vy, PINDX(jbp));
							}
						}

						// for the internal energy 
						die += -pi * Vec2DDotP(&ui,&dS);

						// for the total energy 
						ua.x = Half*(jbp_vx + ibp_vx); 
						ua.y = Half*(jbp_vy + ibp_vy);
						dte += -pi*Vec2DDotP(&ua, &dS);

						// for the kinetic energy 
						ub.x = ibp_vx;
						ub.y = ibp_vy;
						dke += -pi*Vec2DDotP(&ub, &dS);

						// for the force 
						fx += -pi * dS.x;
						fy += -pi * dS.y;

						if(center.indx == 130560)
						{
							DEBUGPRINT("p(%ld,%ld): xy= %g %g has ax/y= %g %g pij= %g dS= %g %g pi= %g : p%lu rj= %g %g pj= %g :: jindx= (%ld,%ld)\n", 
									PINDX(ibp)%256, PINDX(ibp)/256,
									ibp->x,ibp->y,
									fx/ibp->mass,fy/ibp->mass,
									pi,dS.x,dS.y, ibp_pressure, PINDX(jbp),jbp->x, jbp->y, jbp_pressure, PINDX(jbp)%256, PINDX(jbp)/256);
						}

						Voro2D_point dv;
						dv.x = (jbp_vx - ibp_vx); 
						dv.y = (jbp_vy - ibp_vy); 

						postype VdotR = Vec2DDotP(&dv,&er);
						postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
						postype dt = 2*Courant*dramp/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d %ld x/y= %g %g : %g %g : %g %g : %g : dv.xy= %g %g\n",
                                    MYID(simpar), i, PINDX(p+i), p[i].x, p[i].y, dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y);
                            exit(9);
                        }
                        ibp->dt = MIN(ibp->dt,dt);
                        if(dt < Dtime){
                            Dtime = dt;
                            // isave = id;
                        }
                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
               	ibp->die = die;
//                ibp->dte = dte; 
//				ibp->dke = dke;
                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;

				if(center.indx == 130560) 
				{
					DEBUGPRINT("+p%ld has ax/y= %g %g xy= %g %g : vxy= %g %g ::: rho= %g\n", PINDX(ibp),ibp->ax,ibp->ay, ibp->x, ibp->y,
							ibp->vx, ibp->vy, ibp->den);
							
				}

                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g with p= %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin,
							ibp->pressure);
					exit(0);
                }

            }
            my_free(p);my_free(neighbors); my_free(neighwork);

        }
        my_free(vorocorner);
    }
	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
	{
		postype TDtime;
        MPI_Reduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MAX, 0, MPI_COMM(simpar));
        if(MYID(simpar) == 0) Dtime = TDtime;
        MPI_Bcast(&Dtime, 1, MPI_POSTYPE,  0, MPI_COMM(simpar));
	}
	return Dtime;
}
double getAccVoro2D(SimParameters *simpar, postype xmin, postype ymin, 
		postype xmax, postype ymax,
		postype OrderOfAccuracy, postype Courant, postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int isave = -1;
    postype Dtime = 1.e10;


	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	postype cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar);
    int mx, my;
    BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
    BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	CellType *cells = VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my);
	mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);


	postype dtold = GAS_dtold(simpar);



	float alphavis = GAS_AlphaVis(simpar);
	float betavis = GAS_BetaVis(simpar);
	float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
	float eps2vis = GAS_EPSVIS(simpar)*GAS_EPSVIS(simpar);


//	DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);

    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=1;iy<my-1;iy++){
    	int mp=1000;
    	Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp);
    	postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh);
			int i;
            for(i=0;i<np;i++){

                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i); 
                center.csound = p[i].csound;
                center.w2 = p[i].w2;

				treevorork4particletype *ibp = p[i].bp;
				ibp->dt = 1.e10; // initialization of dt

                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
//                ibp->volume = Area2DPolygon(vorocorner, mp);
				// to use registers
				postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
				ibp_vx = ibp->vx;
				ibp_vy = ibp->vy;
				ibp_csound = ibp->csound;
				ibp_pressure = ibp->pressure;
				ibp_den = ibp->den;


				Voro2D_Corner *tmp,*tmp2;
				tmp = vorocorner;
				int j;
				double die, dte,dke,fx,fy;
				die = dte = dke = fx = fy = 0;
				do {
					if(tmp->upperrelated >=0)
					{
						// to use registers
						treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
						postype jbp_vx = jbp->vx;
						postype jbp_vy = jbp->vy;
						postype jbp_csound = jbp->csound;
						postype jbp_pressure = getNeighborPressure(jbp,ibp); // to consider the wallflag
						postype jbp_den = jbp->den;

						Voro2D_point line; 
	                    tmp2 = tmp->upperlink; 
						line.x = tmp2->x - tmp->x; 
						line.y = tmp2->y - tmp->y;
	                    Voro2D_point dS = voro2D_norm(&line); 
						postype facearea = Vec2DLength(tmp,tmp2); 
						dS.x = facearea*dS.x;
						dS.y = facearea*dS.y;

						/*
						postype pi = getPressure2D(
								ibp_pressure,jbp_pressure, tmp, neighwork, OrderOfAccuracy);
								*/
						postype pi = VoroRK4_Pressure2D(
								ibp_pressure,jbp_pressure, tmp, neighwork, simpar, RT);

						Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

						postype dramp = sqrt(Vec2DDotP(&dr, &dr));
						Voro2D_point er;
						er.x = dr.x/dramp;
						er.y = dr.y/dramp;

						Voro2D_point ui,ua,ub;
						/*
						ui.x = Half*(jbp_vx - ibp_vx); 
						ui.y = Half*(jbp_vy - ibp_vy);
						*/
						ui = get2dUpqradRk4(ibp,jbp,dtold);

						{
							postype rvel = Vec2DDotP(&er, &ui);
							if(rvel<0){
								postype mu = rvel/(dramp/etavis + eps2vis*etavis/dramp);
								postype meanden = 0.5*(ibp_den + jbp_den);
								postype meanCsound = 0.5*(ibp_csound + jbp_csound);
								pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
								if(center.indx==130816) DEBUGPRINT("p(%ld,%ld) has pij= %g rvel= %g : w(ibp)= %g %g ::: w(jbp) = %g %g :: vi= %g %g :: vj= %g %g :: jindx= %ld\n", 
									PINDX(ibp)/256,PINDX(ibp)%256, 
									pi, rvel, ibp->w2, ibp->w2old, jbp->w2, jbp->w2old,
									ibp->vx, ibp->vy, jbp->vx,jbp->vy, PINDX(jbp));
							}
						}

						// for the internal energy 
						die += -pi * Vec2DDotP(&ui,&dS);

						// for the total energy 
						ua.x = Half*(jbp_vx + ibp_vx); 
						ua.y = Half*(jbp_vy + ibp_vy);
						dte += -pi*Vec2DDotP(&ua, &dS);

						// for the kinetic energy 
						ub.x = ibp_vx;
						ub.y = ibp_vy;
						dke += -pi*Vec2DDotP(&ub, &dS);

						// for the force 
						fx += -pi * dS.x;
						fy += -pi * dS.y;

						if(center.indx == 130816)
						{
							DEBUGPRINT("p(%ld,%ld): xy= %g %g has ax/y= %g %g pij= %g dS= %g %g pi= %g : p%lu rj= %g %g pj= %g :: jindx= (%ld,%ld)\n", 
									PINDX(ibp)/256, PINDX(ibp)%256,
									ibp->x,ibp->y,
									fx/ibp->mass,fy/ibp->mass,
									pi,dS.x,dS.y, ibp_pressure, PINDX(jbp),jbp->x, jbp->y, jbp_pressure, PINDX(jbp)/256, PINDX(jbp)%256);
						}

						Voro2D_point dv;
						dv.x = (jbp_vx - ibp_vx); 
						dv.y = (jbp_vy - ibp_vy); 

						postype VdotR = Vec2DDotP(&dv,&er);
						postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
						postype dt = 2*Courant*dramp/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d %ld x/y= %g %g : %g %g : %g %g : %g : dv.xy= %g %g\n",
                                    MYID(simpar), i, PINDX(p+i), p[i].x, p[i].y, dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y);
                            exit(9);
                        }
                        ibp->dt = MIN(ibp->dt,dt);
                        if(dt < Dtime){
                            Dtime = dt;
                            // isave = id;
                        }
                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
               	ibp->die = die;
//                ibp->dte = dte; 
//				ibp->dke = dke;
                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;

				if(center.indx == 130816) 
				{
					DEBUGPRINT("+p%ld has ax/y= %g %g xy= %g %g : vxy= %g %g ::: rho= %g\n", PINDX(ibp),ibp->ax,ibp->ay, ibp->x, ibp->y,
							ibp->vx, ibp->vy, ibp->den);
							
				}

                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g with p= %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin,
							ibp->pressure);
					exit(0);
                }

            }
            my_free(p);my_free(neighbors); my_free(neighwork);

        }
        my_free(vorocorner);
    }
	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
	{
		postype TDtime;
        MPI_Reduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MAX, 0, MPI_COMM(simpar));
        if(MYID(simpar) == 0) Dtime = TDtime;
        MPI_Bcast(&Dtime, 1, MPI_POSTYPE,  0, MPI_COMM(simpar));
	}
	return Dtime;
}
void exam2d_centroidShift(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;

	postype xmin,ymin,zmin,xmax,ymax,zmax, cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
	// xmin,ymin,xmax,ymax are the boundaries of the local domain.
	xmin = Xmin_HydroExam(simpar)-cellsize;
    ymin = Ymin_HydroExam(simpar)-cellsize;
    xmax = Xmax_HydroExam(simpar)+cellsize;
    ymax = Ymax_HydroExam(simpar)+cellsize;

	int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);
	if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingAllTreeParticles);
	// prepare the linked-list cells for mkLinkedList2D()  and tree findings below
    CellType *cells = VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my);
    // building the linked list with "paddingAllTreeParticles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
    mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);
	// repositioning the memory space 
	bp = VORORK4_TBP(simpar);
    int i;
    int iy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(iy=1;iy<my-1;iy++){ 
		int mp=2000; 
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp); 
		postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a; 
		int ix; 
		for(ix=1;ix<mx-1;ix++){ 
			int np; 
			treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np); 
			int nneigh; 
			Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh); 
			Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh); 
			int i; 
			for(i=0;i<np;i++){ 
				Voro2D_point center; 
				center.x = p[i].x; 
				center.y = p[i].y; 
				center.indx = PINDX(p+i); 
				center.w2 = p[i].w2; 
				treevorork4particletype *ibp = p[i].bp; 
				int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
                // in this case ax, ay are used for the centroid vector
                Voro2D_Corner centroid = Centroid2DPolygon(vorocorner); 
				ibp->ax = centroid.x; 
				ibp->ay = centroid.y; 
			} 
			my_free(p);my_free(neighbors); my_free(neighwork); 
		} 
		my_free(vorocorner); 
	} 
	my_free(VORORK4_TBPP(simpar)); 
	my_free(VORO_BASICCELL(simpar));
    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;
	float fshift = GAS_FCENTROID(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	for(i=0;i<VORO_NP(simpar);i++){ 
		bp[i].x += fshift*bp[i].ax; // move it to the centroid  with a factor of shift
		bp[i].y += fshift*bp[i].ay; 
		bp[i].x = fmod(bp[i].x+Lx,Lx); 
		bp[i].y = fmod(bp[i].y+Ly,Ly); 
	}
    // migrate particles between mpi ranks 
	migrateTreeVorork4Particles(simpar);
}
void exam2dUpdateVol(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
//  determine the minimum dpq for each particle by updating the w2ceil & w2
	if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingAllTreeParticles); 

	postype xmin,ymin,zmin,xmax,ymax,zmax, cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
	// xmin,ymin,xmax,ymax are the boundaries of the local domain.
	xmin = Xmin_HydroExam(simpar)-cellsize;
    ymin = Ymin_HydroExam(simpar)-cellsize;
    xmax = Xmax_HydroExam(simpar)+cellsize;
    ymax = Ymax_HydroExam(simpar)+cellsize;

	int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);
	// prepare the linked-list cells for mkLinkedList2D()  and tree findings below
	CellType *cells = VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my);
	// building the linked list with "paddingAllTreeParticles()" which pads the domain
	// with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
	mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);

	int i;
	int iy;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=1000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.csound = p[i].csound;
                center.w2 = p[i].w2;
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
				treevorork4particletype *ibp = p[i].bp;
				ibp->volume = Area2DPolygon(vorocorner,mp);
				ibp->den = ibp->mass / ibp->volume;
			}
			my_free(neighwork);
			my_free(neighbors); 
			my_free(p);
		}
		my_free(vorocorner);
	
	}
	my_free(VORO_BASICCELL(simpar));

	VORO_NPAD(simpar) = 0;
	my_free(VORORK4_TBPP(simpar));
}
double exam2d_vph_rk4(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	postype xmin,ymin,xmax,ymax,OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);
//	DEBUGPRINT("P%d has x/y minmax= %g %g : %g %g ::: %g\n", MYID(simpar), xmin,ymin,xmax,ymax, cellsize);


    int i;
    postype Dtime,dt;
    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) bp[i].rk4.w2backup = bp[i].w2;
	// Runge-Kutta 4-th order time evolution of r and vr
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D);
    Dtime = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){ 
		bp[i].rk4.k1x = bp[i].vx*Dtime;
        bp[i].rk4.k1y = bp[i].vy*Dtime;
        bp[i].rk4.k1vx = bp[i].ax*Dtime;
        bp[i].rk4.k1vy = bp[i].ay*Dtime;
        bp[i].rk4.k1ie = bp[i].die*Dtime;
    }
    // RK4 second
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += bp[i].rk4.k1x*0.5;
        bp[i].y += bp[i].rk4.k1y*0.5;
        bp[i].vx += bp[i].rk4.k1vx*0.5;
        bp[i].vy += bp[i].rk4.k1vy*0.5;
        bp[i].ie += bp[i].rk4.k1ie*0.5;
        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
	migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D);
    dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
        paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k2x = bp[i].vx*Dtime;
        bp[i].rk4.k2y = bp[i].vy*Dtime;
        bp[i].rk4.k2vx = bp[i].ax*Dtime;
        bp[i].rk4.k2vy = bp[i].ay*Dtime;
        bp[i].rk4.k2ie = bp[i].die*Dtime;
    }
    // RK4 third
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += (bp[i].rk4.k2x-bp[i].rk4.k1x)*0.5;
        bp[i].y += (bp[i].rk4.k2y-bp[i].rk4.k1y)*0.5;
        bp[i].vx += (bp[i].rk4.k2vx-bp[i].rk4.k1vx)*0.5;
        bp[i].vy += (bp[i].rk4.k2vy-bp[i].rk4.k1vy)*0.5;
        bp[i].ie += (bp[i].rk4.k2ie-bp[i].rk4.k1ie)*0.5;

        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D);
	dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k3x = bp[i].vx*Dtime;
        bp[i].rk4.k3y = bp[i].vy*Dtime;
        bp[i].rk4.k3vx = bp[i].ax*Dtime;
        bp[i].rk4.k3vy = bp[i].ay*Dtime;
        bp[i].rk4.k3ie = bp[i].die*Dtime;
    }
    // RK4 fourth
    for(i=0;i<VORO_NP(simpar);i++){

        bp[i].x += bp[i].rk4.k3x - bp[i].rk4.k2x*0.5;
        bp[i].y += bp[i].rk4.k3y - bp[i].rk4.k2y*0.5;
        bp[i].vx += bp[i].rk4.k3vx - bp[i].rk4.k2vx*0.5;
        bp[i].vy += bp[i].rk4.k3vy - bp[i].rk4.k2vy*0.5;
        bp[i].ie += bp[i].rk4.k3ie - bp[i].rk4.k2ie*0.5;

        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D);
    dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k4x = bp[i].vx*Dtime;
        bp[i].rk4.k4y = bp[i].vy*Dtime;
        bp[i].rk4.k4vx = bp[i].ax*Dtime;
        bp[i].rk4.k4vy = bp[i].ay*Dtime;
        bp[i].rk4.k4ie = bp[i].die*Dtime;
    }
    // migrate particles between mpi ranks
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x -= bp[i].rk4.k3x;
        bp[i].y -= bp[i].rk4.k3y;
        bp[i].vx -= bp[i].rk4.k3vx;
        bp[i].vy -= bp[i].rk4.k3vy;
        bp[i].ie -= bp[i].rk4.k3ie;
        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
	// migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

    // to finalize ke & ie updates 
	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
    // Update position and velocity using RK4
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].x  += (bp[i].rk4.k1x +2*bp[i].rk4.k2x +2*bp[i].rk4.k3x +bp[i].rk4.k4x )/6.;
        bp[i].y  += (bp[i].rk4.k1y +2*bp[i].rk4.k2y +2*bp[i].rk4.k3y +bp[i].rk4.k4y )/6.;
        bp[i].vx += (bp[i].rk4.k1vx+2*bp[i].rk4.k2vx+2*bp[i].rk4.k3vx+bp[i].rk4.k4vx)/6.;
        bp[i].vy += (bp[i].rk4.k1vy+2*bp[i].rk4.k2vy+2*bp[i].rk4.k3vy+bp[i].rk4.k4vy)/6.;
        bp[i].ie += (bp[i].rk4.k1ie+2*bp[i].rk4.k2ie+2*bp[i].rk4.k3ie+bp[i].rk4.k4ie)/6.;

        bp[i].x = fmod(bp[i].x+Lx,Lx);
        bp[i].y = fmod(bp[i].y+Ly,Ly);
    }

    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);

	//centroid shift
	if(GAS_FCENTROID(simpar)>0) 
		exam2d_centroidShift(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	// update volume & density with one-step evolved positions 
    exam2dUpdateVol(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	bp = VORORK4_TBP(simpar);
	postype dmean = SIMBOX(simpar).x.max/NX(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].den = bp[i].mass/bp[i].volume;
        bp[i].pressure = bp[i].ie/bp[i].volume * (Gamma-1);
        bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
//        bp[i].w2 = measureW2(simpar,dmean, bp[i].pressure, bp[i].den);
        bp[i].w2 = getw2forHydroParticle(simpar,(bp+i));
        bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
    }
    return Dtime;
}

void dumpRk4Data2D(SimParameters *simpar, int nstep, postype t, postype dt){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int myid, nid;
	myid = MYID(simpar);
	nid = NID(simpar);
	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);


	char outfile[190];
	sprintf(outfile,"glout.%.6d.dat",nstep);
	int i;
	int snp=0;
	FILE *wp;
	int tnp = VORO_TNP(simpar);
	for(i=0;i<nid;i++){
		if(myid ==i){
			if(myid==0){
				wp = fopen(outfile,"w"); 
				fwrite(&tnp, sizeof(int),1,wp);
			}
			else wp = fopen(outfile,"a");
			fwrite(bp, sizeof(treevorork4particletype), np, wp);
			fclose(wp);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if(myid==0){
		wp = fopen(outfile,"a");
		fwrite(&t, sizeof(postype), 1, wp);
		fwrite(&nx, sizeof(int), 1, wp);
		fwrite(&ny, sizeof(int), 1, wp);
		fwrite(&(SIMBOX(simpar).x.max), sizeof(postype), 1, wp);
		fwrite(&(SIMBOX(simpar).y.max), sizeof(postype), 1, wp);
		fwrite(&GAS_AlphaVis(simpar), sizeof(float), 1, wp);
		fwrite(&GAS_BetaVis(simpar), sizeof(float), 1, wp);
		fwrite(&GAS_ETAVIS(simpar), sizeof(float), 1, wp);
		fwrite(&GAS_EPSVIS(simpar), sizeof(float), 1, wp);
		fclose(wp);
	}
	/*
	MPI_Barrier(MPI_COMM_WORLD);
	DEBUGPRINT("P%d is now exiting the dump\n", myid);
	*/
}


void readRk4Data2D(SimParameters *simpar, postype *t, int nstep){
	char infile[190];
	sprintf(infile,"glout.%.6d.dat", nstep);
	int myid = MYID(simpar);
	if(myid ==0){
		printf("P0: is now reading starting file: %s\n", infile);
		FILE *fp = fopen(infile,"r");
		int np;
		fread(&np, sizeof(int), 1, fp);
		VORO_TNP(simpar) = VORO_NP(simpar) = np;
		VORORK4_TBP(simpar) = (treevorork4particletype*)my_malloc(sizeof(treevorork4particletype)*np);
		fread(VORORK4_TBP(simpar),sizeof(treevorork4particletype), np,fp);
		fread(t,sizeof(postype), 1,fp);
		fclose(fp);
	}
	else {
		VORO_NP(simpar) = 0;
		VORORK4_TBP(simpar) = (treevorork4particletype*)my_malloc(sizeof(treevorork4particletype)*100);;

	}
	DEBUGPRINT("P%d: after reading iput data\n", myid);
	migrateTreeVorork4Particles(simpar);
	MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Barrier(MPI_COMM(simpar));
	BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
}

double exam2d_vph(SimParameters *simpar, 
		double (*measureW2Exam2D)(SimParameters *, double, double, double),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		)
{
    postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
    treevorork4particletype *bp = VORORK4_TBP(simpar);
    int nbp = VORO_NP(simpar);
    int isave = -1;
    postype Dtime = 1.e10;
    postype Courant = GAS_COURANT(simpar);
    postype cellsize;


    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;
    postype xmin,ymin,zmin,xmax,ymax,zmax, pwidth;
    pwidth = BASICCELL_CELLWIDTH(simpar) = cellsize = HydroGridSize(simpar);
    xmin = Xmin_HydroExam(simpar)-pwidth;
    ymin = Ymin_HydroExam(simpar)-pwidth;
    xmax = Xmax_HydroExam(simpar)+pwidth;
    ymax = Ymax_HydroExam(simpar)+pwidth;

    postype dtold = GAS_dtold(simpar);

//  DEBUGPRINT("P%d has dtold = %g\n", MYID(simpar), dtold);

    int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);

    // prepare the linked-list cells for mkLinkedList2D()  and tree findings below
    CellType *cells = VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my);
    // building the linked list with "paddingTreeVorork4Particles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingTreeVorork4Particles);

    postype Gamma = GAS_GAMMA(simpar);

    float alphavis = GAS_AlphaVis(simpar);
    float betavis = GAS_BetaVis(simpar);
    float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
    float eps2vis = GAS_EPSVIS(simpar)*GAS_EPSVIS(simpar);


//  DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);


    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = findCellRk4BP2D(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = searchCellRk4Neighbors2D(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.csound = p[i].csound;
                center.w2 = p[i].w2;

                treevorork4particletype *ibp = p[i].bp;
                ibp->dt = 1.e10; // initialization of dt
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
                ibp->volume = Area2DPolygon(vorocorner, mp);
                /* These are to use register */
                postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
                ibp_vx = ibp->vx;
                ibp_vy = ibp->vy;
                ibp_csound = ibp->csound;
                ibp_pressure = ibp->pressure;
                ibp_den = ibp->den;


                Voro2D_Corner *tmp,*tmp2;
                tmp = vorocorner;
                int j;
                double die, dte,dke,fx,fy;
                die = dte = dke = fx = fy = 0;
                do {
                    if(tmp->upperrelated >=0){
                        /* These are to use register */
                        treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
                        postype jbp_vx = jbp->vx;
                        postype jbp_vy = jbp->vy;
                        postype jbp_csound = jbp->csound;
//                        postype jbp_pressure = jbp->pressure;
						postype jbp_pressure = getNeighborPressure(jbp,ibp);

                        postype jbp_den = jbp->den;

                        Voro2D_point line;
                        tmp2 = tmp->upperlink;
                        line.x = tmp2->x - tmp->x;
                        line.y = tmp2->y - tmp->y;
                        Voro2D_point dS = voro2D_norm(&line);
                        postype facearea = Vec2DLength(tmp,tmp2);
                        dS.x = facearea*dS.x;
                        dS.y = facearea*dS.y;

                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,RT);

                        Voro2D_point ui,ua,ub;
                        /*
                        ui.x = Half*(jbp_vx - ibp_vx); 
                        ui.y = Half*(jbp_vy - ibp_vy);
                        */
                        ui = get2dUpqradRk4(ibp,jbp,dtold);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &ui);
                        postype mu, meanden, meanCsound;
                        if(rvel<0){
                            mu = rvel/(dramp/etavis + eps2vis*etavis/dramp);
                            meanden = 0.5*(ibp_den + jbp_den);
                            meanCsound = 0.5*(ibp_csound + jbp_csound);
                            pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        /* for the internal energy */
                        die += -pi * Vec2DDotP(&ui,&dS);

                        /* for the total energy */
                        ua.x = Half*(jbp_vx + ibp_vx);
                        ua.y = Half*(jbp_vy + ibp_vy);
                        dte += -pi*Vec2DDotP(&ua, &dS);

                        /* for the kinetic energy */
                        ub.x = ibp_vx;
                        ub.y = ibp_vy;
                        dke += -pi*Vec2DDotP(&ub, &dS);

                        /* for the force */
                        fx += -pi * dS.x;
                        fy += -pi * dS.y;

                        Voro2D_point dv;
                        dv.x = (jbp_vx - ibp_vx);
                        dv.y = (jbp_vy - ibp_vy);

                        postype VdotR = Vec2DDotP(&dv,&er);
                        postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
                        postype dt = 2*Courant*dramp/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d %ld : %g %g : %g %g : %g : dv.xy= %g %g\n",
                                    MYID(simpar), i, PINDX(p+i), dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y);
                            exit(9);
                        }
                        ibp->dt = MIN(ibp->dt,dt);
                        if(dt < Dtime){
                            Dtime = dt;
                            // isave = id;
                        }
                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
                ibp->dte = dte;
                ibp->dke = dke;
                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;
                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin);
                }

            }
            my_free(p);my_free(neighbors); my_free(neighwork);

        }
        my_free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));

//  DEBUGPRINT("P%d is exiting here\n", MYID(simpar)); exit(9);

    postype TDtime;
    MPI_Reduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MAX, 0, MPI_COMM(simpar));
    if(MYID(simpar) == 0) Dtime = TDtime;
    MPI_Bcast(&Dtime, 1, MPI_POSTYPE,  0, MPI_COMM(simpar));

    int i;
    postype dmax=0;
//  postype eps2 = KH_EPS(simpar)*KH_EPS(simpar);
#ifdef _OPENMP
#pragma omp parallel for reduction(max: dmax)
#endif
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].vx += bp[i].ax *Dtime;
        bp[i].vy += bp[i].ay *Dtime;
        bp[i].x += bp[i].vx *Dtime;
        bp[i].y += bp[i].vy *Dtime;
        dmax = MAX(dmax, fabs(bp[i].ay));
        bp[i].x = fmod(bp[i].x+Lx,Lx);
        bp[i].y = fmod(bp[i].y+Ly,Ly);
        bp[i].w2old = bp[i].w2;
    }
    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);
    if(0){
        postype xmin,ymin,xmax,ymax;
        xmin = ymin = 1.e20;
        xmax = ymax =-1.e20;
        for(i=0;i<VORO_NP(simpar);i++){
            xmin = MIN(xmin, bp[i].x);
            ymin = MIN(ymin, bp[i].y);
            xmax = MAX(xmax, bp[i].x);
            ymax = MAX(ymax, bp[i].y);
        }
    }

// This is for centroid shift.
	if(GAS_FCENTROID(simpar) > 0 ){
	    if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingTreeVorork4Particles);
    // Update the linked list 
		mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingTreeVorork4Particles);

    // repositioning the memory space 
	    bp = VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	    for(iy=1;iy<my-1;iy++){ 
			int mp=2000; 
			Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp); 
			postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a; 
			int ix; 
			for(ix=1;ix<mx-1;ix++){ 
				int np; 
				treevorork4particletype *p = findCellRk4BP2D(simpar,ix,iy,&np); 
				int nneigh; 
				Voro2D_point *neighbors = searchCellRk4Neighbors2D(simpar,ix,iy,&nneigh); 
				Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh); 
				int i; 
				for(i=0;i<np;i++){ 
					Voro2D_point center; 
					center.x = p[i].x; 
					center.y = p[i].y; 
					center.indx = PINDX(p+i); 
					center.w2 = p[i].w2; 
					
					treevorork4particletype *ibp = p[i].bp; 
					int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
					// in this case ax, ay are used for the centroid vector
					Voro2D_Corner centroid = Centroid2DPolygon(vorocorner);
					ibp->ax = centroid.x;
					ibp->ay = centroid.y; 
				} 
				my_free(p);my_free(neighbors); my_free(neighwork); 
			} 
			my_free(vorocorner); 
		} 
		my_free(VORORK4_TBPP(simpar));
		float fshift = GAS_FCENTROID(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	    for(i=0;i<VORO_NP(simpar);i++){ 
			bp[i].x += fshift*bp[i].ax; // move it to the centroid  with a factor of shift
			bp[i].y += fshift*bp[i].ay; 
        	bp[i].x = fmod(bp[i].x+Lx,Lx);
	        bp[i].y = fmod(bp[i].y+Ly,Ly);
		} 
		// migrate particles between mpi ranks 
		migrateTreeVorork4Particles(simpar);
	}
//
//
	if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingTreeVorork4Particles);
    // Update the linked list
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingTreeVorork4Particles);


    // repositioning the memory space
    bp = VORORK4_TBP(simpar);



#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = findCellRk4BP2D(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = searchCellRk4Neighbors2D(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.w2 = p[i].w2;

                treevorork4particletype *ibp = p[i].bp;
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
//              ibp->w2 = ibp->w2ceil = center.w2ceil;
                // These are to use register 
                postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
                ibp_vx = ibp->vx;
                ibp_vy = ibp->vy;
                ibp_csound = ibp->csound;
                ibp_pressure = ibp->pressure;
                ibp_den = ibp->den;

                ibp->volume = Area2DPolygon(vorocorner, mp);




                Voro2D_Corner *tmp,*tmp2;
                tmp = vorocorner;
                int j;
//                double die, dte,dke,fx,fy;
//               die = dte = dke = fx = fy = 0;
				double die, cx,cy;
				die = cx = cy = 0;
                do {
                    if(tmp->upperrelated >=0){
                        // These are to use register 
                        treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
                        postype jbp_vx = jbp->vx;
                        postype jbp_vy = jbp->vy;
                        postype jbp_csound = jbp->csound;
//                        postype jbp_pressure = jbp->pressure;
						postype jbp_pressure = getNeighborPressure(jbp,ibp);
                        postype jbp_den = jbp->den;

                        Voro2D_point line;
                        tmp2 = tmp->upperlink;
                        line.x = tmp2->x - tmp->x;
                        line.y = tmp2->y - tmp->y;
                        Voro2D_point dS = voro2D_norm(&line);
                        postype facearea = Vec2DLength(tmp,tmp2);
                        dS.x = facearea*dS.x;
                        dS.y = facearea*dS.y;
//                      postype pi = Half*(jbp_pressure + ibp_pressure);
                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,RT);

                        Voro2D_point ui,ua,ub;
//                        ui.x = Half*(jbp_vx - ibp_vx);
 //                       ui.y = Half*(jbp_vy - ibp_vy);
                        ui = get2dUpqradRk4(ibp,jbp,dtold);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &ui);
                        if(rvel<0){
                             postype mu = rvel/(dramp/etavis + eps2vis*etavis/dramp);
                             postype meanden = 0.5*(ibp_den + jbp_den);
                             postype meanCsound = 0.5*(ibp_csound + jbp_csound);
                             pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy 
                        die += -pi * Vec2DDotP(&ui,&dS);

                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
//                ibp->ax = fx/ibp->mass;
//                ibp->ay = fy/ibp->mass;
            }
            my_free(p);my_free(neighbors); my_free(neighwork);
        }
        my_free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));
    my_free(VORO_BASICCELL(simpar));

    // update the w2ceil
    postype dmean = SIMBOX(simpar).x.max/NX(simpar);

//  DEBUGPRINT("P%d is now updating the internal energy\n", MYID(simpar));
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].ie += bp[i].die * Dtime;
        bp[i].den  = bp[i].mass/bp[i].volume;
        bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1);
        bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
        if(GAS_Kappa(simpar) <0) {
            bp[i].w2 = -GAS_Kappa(simpar);
        }
        else if (GAS_Kappa(simpar) >0){
          	bp[i].w2 = getw2forHydroParticle(simpar,(bp+i));
            bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
        }
        if(bp[i].volume ==0)
                DEBUGPRINT("P%d has wrong volume at i= %d with den= %g at xy= %g %g\n", MYID(simpar),
                        i,bp[i].den, bp[i].x, bp[i].y);
        if(isnan(bp[i].csound)){
                DEBUGPRINT("P%d has nan at the sound speed = %d with"
                        "id= %ld: pressure= %g den= %g Delta IE=%g dIE= %g Dt= %g\n",
                        MYID(simpar), i,PINDX(bp+i), bp[i].pressure, bp[i].den, bp[i].ie,bp[i].die,Dtime);
            exit(9);
        }

    }
    MPI_Barrier(MPI_COMM_WORLD);
    return Dtime;
}
//
//
//
double exam2d_vph_int_rt(SimParameters *simpar, 
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		) {
    postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
    treevorork4particletype *bp = VORORK4_TBP(simpar);
    int nbp = VORO_NP(simpar);
    int isave = -1;
    postype Dtime = 1.e10;
    postype Courant = GAS_COURANT(simpar);
    postype cellsize;


    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;
    postype xmin,ymin,zmin,xmax,ymax,zmax, pwidth;
    pwidth = BASICCELL_CELLWIDTH(simpar) = cellsize = HydroGridSize(simpar);
    xmin = Xmin_HydroExam(simpar)-pwidth;
    ymin = Ymin_HydroExam(simpar)-pwidth;
    xmax = Xmax_HydroExam(simpar)+pwidth;
    ymax = Ymax_HydroExam(simpar)+pwidth;

    postype dtold = GAS_dtold(simpar);

//  DEBUGPRINT("P%d has dtold = %g\n", MYID(simpar), dtold);

    int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);

    // prepare the linked-list cells for mkLinkedList2D()  and tree findings below
    CellType *cells = VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my);
    // building the linked list with "paddingAllTreeParticles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

    postype Gamma = GAS_GAMMA(simpar);

    float alphavis = GAS_AlphaVis(simpar);
    float betavis = GAS_BetaVis(simpar);
    float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
    float eps2vis = GAS_EPSVIS(simpar)*GAS_EPSVIS(simpar);


//  DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);


    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.csound = p[i].csound;
                center.w2 = p[i].w2;

                treevorork4particletype *ibp = p[i].bp;
                ibp->dt = 1.e10; // initialization of dt
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
                ibp->volume = Area2DPolygon(vorocorner, mp);
                // These are to use register 
                postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
                ibp_vx = ibp->vx;
                ibp_vy = ibp->vy;
                ibp_csound = ibp->csound;
                ibp_pressure = ibp->pressure;
                ibp_den = ibp->den;


                Voro2D_Corner *tmp,*tmp2;
                tmp = vorocorner;
                int j;
                double die, dte,dke,fx,fy;
                die = dte = dke = fx = fy = 0;
                do {
                    if(tmp->upperrelated >=0){
                        // These are to use register 
                        treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
                        postype jbp_vx = jbp->vx;
                        postype jbp_vy = jbp->vy;
                        postype jbp_csound = jbp->csound;
                        postype jbp_pressure = jbp->pressure; // for RT, the wall particle doesnot mirror the pressure
//						postype jbp_pressure = getNeighborPressure(jbp,ibp);
                        postype jbp_den = jbp->den;

                        Voro2D_point line;
                        tmp2 = tmp->upperlink;
                        line.x = tmp2->x - tmp->x;
                        line.y = tmp2->y - tmp->y;
                        Voro2D_point dS = voro2D_norm(&line);
                        postype facearea = Vec2DLength(tmp,tmp2);
                        dS.x = facearea*dS.x;
                        dS.y = facearea*dS.y;

                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,RT);

                        Voro2D_point ui,ua,ub;
                        //
//                        ui.x = Half*(jbp_vx - ibp_vx); 
 //                       ui.y = Half*(jbp_vy - ibp_vy);
                        ui = get2dUpqradRk4(ibp,jbp,dtold);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &ui);
                        postype mu, meanden, meanCsound;
                        if(rvel<0){
                            mu = rvel/(dramp/etavis + eps2vis*etavis/dramp);
                            meanden = 0.5*(ibp_den + jbp_den);
                            meanCsound = 0.5*(ibp_csound + jbp_csound);
                            pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy 
                        die += -pi * Vec2DDotP(&ui,&dS);

                        // for the total energy
                        ua.x = Half*(jbp_vx + ibp_vx);
                        ua.y = Half*(jbp_vy + ibp_vy);
                        dte += -pi*Vec2DDotP(&ua, &dS);

                        // for the kinetic energy 
                        ub.x = ibp_vx;
                        ub.y = ibp_vy;
                        dke += -pi*Vec2DDotP(&ub, &dS);

                        // for the force 
                        fx += -pi * dS.x;
                        fy += -pi * dS.y;

                        Voro2D_point dv;
                        dv.x = (jbp_vx - ibp_vx);
                        dv.y = (jbp_vy - ibp_vy);

                        postype VdotR = Vec2DDotP(&dv,&er);
                        postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
                        postype dt = 2*Courant*dramp/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d %ld : %g %g : %g %g : %g : dv.xy= %g %g\n",
                                    MYID(simpar), i, PINDX(p+i), dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y);
                            exit(9);
                        }
                        ibp->dt = MIN(ibp->dt,dt);
                        if(dt < Dtime){
                            Dtime = dt;
                            // isave = id;
                        }
                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
                ibp->dte = dte;
                ibp->dke = dke;
                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;
                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin);
                }

            }
            my_free(p);my_free(neighbors); my_free(neighwork);

        }
        my_free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));

//  DEBUGPRINT("P%d is exiting here\n", MYID(simpar)); exit(9);

    postype TDtime;
    MPI_Reduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MAX, 0, MPI_COMM(simpar));
    if(MYID(simpar) == 0) Dtime = TDtime;
    MPI_Bcast(&Dtime, 1, MPI_POSTYPE,  0, MPI_COMM(simpar));

    int i;
    postype dmax=0;
//  postype eps2 = KH_EPS(simpar)*KH_EPS(simpar);
#ifdef _OPENMP
#pragma omp parallel for reduction(max: dmax)
#endif
    for(i=0;i<VORO_NP(simpar);i++){
//		if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
		if(targetBP(bp+i, Lx, Ly))
		{
			bp[i].ax += GAS_ACCX(simpar);
			bp[i].ay += GAS_ACCY(simpar);
	        bp[i].vx += bp[i].ax *Dtime; 
			bp[i].vy += bp[i].ay *Dtime; 
			bp[i].x += bp[i].vx *Dtime; 
			bp[i].y += bp[i].vy *Dtime; 
			dmax = MAX(dmax, fabs(bp[i].ay)); 
			bp[i].x = fmod(bp[i].x+Lx,Lx); 
			bp[i].y = fmod(bp[i].y+Ly,Ly); 
			bp[i].w2old = bp[i].w2;
		}
		else {
			bp[i].vx = bp[i].vy = bp[i].ax = bp[i].ay = 0;
		}
    }
    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);

// This is for centroid shift.
	if(GAS_FCENTROID(simpar) > 0 ){
	    if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingAllTreeParticles);
    // Update the linked list 
		mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

    // repositioning the memory space 
	    bp = VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	    for(iy=1;iy<my-1;iy++){ 
			int mp=2000; 
			Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp); 
			postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a; 
			int ix; 
			for(ix=1;ix<mx-1;ix++){ 
				int np; 
				treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np); 
				int nneigh; 
				Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh); 
				Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh); 
				int i; 
				for(i=0;i<np;i++){ 
					Voro2D_point center; 
					center.x = p[i].x; 
					center.y = p[i].y; 
					center.indx = PINDX(p+i); 
					center.w2 = p[i].w2; 
					
					treevorork4particletype *ibp = p[i].bp; 
					int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
					// in this case ax, ay are used for the centroid vector
					Voro2D_Corner centroid = Centroid2DPolygon(vorocorner);
					ibp->ax = centroid.x;
					ibp->ay = centroid.y; 
				} 
				my_free(p);my_free(neighbors); my_free(neighwork); 
			} 
			my_free(vorocorner); 
		} 
		my_free(VORORK4_TBPP(simpar));
		float fshift = GAS_FCENTROID(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	    for(i=0;i<VORO_NP(simpar);i++){ 
//			if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
			if(targetBP(bp+i, Lx, Ly))
			{
				bp[i].x += fshift*bp[i].ax; // move it to the centroid  with a factor of shift
				bp[i].y += fshift*bp[i].ay; 
   		     	bp[i].x = fmod(bp[i].x+Lx,Lx);
		        bp[i].y = fmod(bp[i].y+Ly,Ly);
			}
		} 
		// migrate particles between mpi ranks 
		migrateTreeVorork4Particles(simpar);
	}
//
//
	if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingAllTreeParticles);
    // Update the linked list
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);


    // repositioning the memory space
    bp = VORORK4_TBP(simpar);



#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.w2 = p[i].w2;

                treevorork4particletype *ibp = p[i].bp;
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
//              ibp->w2 = ibp->w2ceil = center.w2ceil;
                // These are to use register 
                postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
                ibp_vx = ibp->vx;
                ibp_vy = ibp->vy;
                ibp_csound = ibp->csound;
                ibp_pressure = ibp->pressure;
                ibp_den = ibp->den;

                ibp->volume = Area2DPolygon(vorocorner, mp);




                Voro2D_Corner *tmp,*tmp2;
                tmp = vorocorner;
                int j;
//                double die, dte,dke,fx,fy;
//               die = dte = dke = fx = fy = 0;
				double die, cx,cy;
				die = cx = cy = 0;
                do {
                    if(tmp->upperrelated >=0){
                        // These are to use register 
                        treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
                        postype jbp_vx = jbp->vx;
                        postype jbp_vy = jbp->vy;
                        postype jbp_csound = jbp->csound;
                        postype jbp_pressure = jbp->pressure; // for RT, the wall particle doesnot mirror the pressure
//						postype jbp_pressure = getNeighborPressure(jbp,ibp);
                        postype jbp_den = jbp->den;

                        Voro2D_point line;
                        tmp2 = tmp->upperlink;
                        line.x = tmp2->x - tmp->x;
                        line.y = tmp2->y - tmp->y;
                        Voro2D_point dS = voro2D_norm(&line);
                        postype facearea = Vec2DLength(tmp,tmp2);
                        dS.x = facearea*dS.x;
                        dS.y = facearea*dS.y;
//                      postype pi = Half*(jbp_pressure + ibp_pressure);
                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,RT);

                        Voro2D_point ui,ua,ub;
//                        ui.x = Half*(jbp_vx - ibp_vx);
 //                       ui.y = Half*(jbp_vy - ibp_vy);
                        ui = get2dUpqradRk4(ibp,jbp,dtold);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &ui);
                        if(rvel<0){
                             postype mu = rvel/(dramp/etavis + eps2vis*etavis/dramp);
                             postype meanden = 0.5*(ibp_den + jbp_den);
                             postype meanCsound = 0.5*(ibp_csound + jbp_csound);
                             pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy 
                        die += -pi * Vec2DDotP(&ui,&dS);

                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
//                ibp->ax = fx/ibp->mass;
//                ibp->ay = fy/ibp->mass;
            }
            my_free(p);my_free(neighbors); my_free(neighwork);
        }
        my_free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));
    my_free(VORO_BASICCELL(simpar));

    // update the w2ceil
    postype dmean = SIMBOX(simpar).x.max/NX(simpar);
//  DEBUGPRINT("P%d is now updating the internal energy\n", MYID(simpar));
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
//		if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
		if(targetBP(bp+i, Lx, Ly))
		{
	        bp[i].ie += bp[i].die * Dtime; 
			bp[i].den  = bp[i].mass/bp[i].volume; 
			bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1); 
			bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den); 
			if(GAS_Kappa(simpar) <0) { 
				bp[i].w2 = -GAS_Kappa(simpar); 
			} 
			else if (GAS_Kappa(simpar) >0){ 
//				bp[i].w2 = measureW2(simpar,dmean, bp[i].pressure, bp[i].den); 
          		bp[i].w2 = getw2forHydroParticle(simpar,(bp+i));
				bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil); 
			} 
			if(bp[i].volume ==0) 
				DEBUGPRINT("P%d has wrong volume at i= %d with den= %g at xy= %g %g\n", MYID(simpar), 
						i,bp[i].den, bp[i].x, bp[i].y); 
			if(isnan(bp[i].csound)){ 
				DEBUGPRINT("P%d has nan at the sound speed = %d with" 
						"id= %ld: pressure= %g den= %g Delta IE=%g dIE= %g Dt= %g\n", 
						MYID(simpar), i,PINDX(bp+i), bp[i].pressure, bp[i].den, bp[i].ie,bp[i].die,Dtime); 
				exit(9);
			}
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return Dtime;
}
double exam2d_vph_int(SimParameters *simpar, 
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		) {
    postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
    treevorork4particletype *bp = VORORK4_TBP(simpar);
    int nbp = VORO_NP(simpar);
    int isave = -1;
    postype Dtime = 1.e10;
    postype Courant = GAS_COURANT(simpar);
    postype cellsize;


    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;
    postype xmin,ymin,zmin,xmax,ymax,zmax, pwidth;
    pwidth = BASICCELL_CELLWIDTH(simpar) = cellsize = HydroGridSize(simpar);
    xmin = Xmin_HydroExam(simpar)-pwidth;
    ymin = Ymin_HydroExam(simpar)-pwidth;
    xmax = Xmax_HydroExam(simpar)+pwidth;
    ymax = Ymax_HydroExam(simpar)+pwidth;

    postype dtold = GAS_dtold(simpar);

//  DEBUGPRINT("P%d has dtold = %g\n", MYID(simpar), dtold);

    int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);

    // prepare the linked-list cells for mkLinkedList2D()  and tree findings below
    CellType *cells = VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my);
    // building the linked list with "paddingAllTreeParticles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

    postype Gamma = GAS_GAMMA(simpar);

    float alphavis = GAS_AlphaVis(simpar);
    float betavis = GAS_BetaVis(simpar);
    float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
    float eps2vis = GAS_EPSVIS(simpar)*GAS_EPSVIS(simpar);


//  DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);


    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.csound = p[i].csound;
                center.w2 = p[i].w2;

                treevorork4particletype *ibp = p[i].bp;
                ibp->dt = 1.e10; // initialization of dt
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
                ibp->volume = Area2DPolygon(vorocorner, mp);
                // These are to use register 
                postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
                ibp_vx = ibp->vx;
                ibp_vy = ibp->vy;
                ibp_csound = ibp->csound;
                ibp_pressure = ibp->pressure;
                ibp_den = ibp->den;


                Voro2D_Corner *tmp,*tmp2;
                tmp = vorocorner;
                int j;
                double die, dte,dke,fx,fy;
                die = dte = dke = fx = fy = 0;
                do {
                    if(tmp->upperrelated >=0){
                        // These are to use register 
                        treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
                        postype jbp_vx = jbp->vx;
                        postype jbp_vy = jbp->vy;
                        postype jbp_csound = jbp->csound;
//                        postype jbp_pressure = jbp->pressure;
						postype jbp_pressure = getNeighborPressure(jbp,ibp);
                        postype jbp_den = jbp->den;

                        Voro2D_point line;
                        tmp2 = tmp->upperlink;
                        line.x = tmp2->x - tmp->x;
                        line.y = tmp2->y - tmp->y;
                        Voro2D_point dS = voro2D_norm(&line);
                        postype facearea = Vec2DLength(tmp,tmp2);
                        dS.x = facearea*dS.x;
                        dS.y = facearea*dS.y;

                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,RT);

                        Voro2D_point ui,ua,ub;
                        //
//                        ui.x = Half*(jbp_vx - ibp_vx); 
 //                       ui.y = Half*(jbp_vy - ibp_vy);
                        ui = get2dUpqradRk4(ibp,jbp,dtold);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &ui);
                        postype mu, meanden, meanCsound;
                        if(rvel<0){
                            mu = rvel/(dramp/etavis + eps2vis*etavis/dramp);
                            meanden = 0.5*(ibp_den + jbp_den);
                            meanCsound = 0.5*(ibp_csound + jbp_csound);
                            pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy 
                        die += -pi * Vec2DDotP(&ui,&dS);

                        // for the total energy
                        ua.x = Half*(jbp_vx + ibp_vx);
                        ua.y = Half*(jbp_vy + ibp_vy);
                        dte += -pi*Vec2DDotP(&ua, &dS);

                        // for the kinetic energy 
                        ub.x = ibp_vx;
                        ub.y = ibp_vy;
                        dke += -pi*Vec2DDotP(&ub, &dS);

                        // for the force 
                        fx += -pi * dS.x;
                        fy += -pi * dS.y;

                        Voro2D_point dv;
                        dv.x = (jbp_vx - ibp_vx);
                        dv.y = (jbp_vy - ibp_vy);

                        postype VdotR = Vec2DDotP(&dv,&er);
                        postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
                        postype dt = 2*Courant*dramp/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d %ld : %g %g : %g %g : %g : dv.xy= %g %g\n",
                                    MYID(simpar), i, PINDX(p+i), dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y);
                            exit(9);
                        }
                        ibp->dt = MIN(ibp->dt,dt);
                        if(dt < Dtime){
                            Dtime = dt;
                            // isave = id;
                        }
                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
                ibp->dte = dte;
                ibp->dke = dke;
                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;
                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin);
                }

            }
            my_free(p);my_free(neighbors); my_free(neighwork);

        }
        my_free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));

//  DEBUGPRINT("P%d is exiting here\n", MYID(simpar)); exit(9);

    postype TDtime;
    MPI_Reduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MAX, 0, MPI_COMM(simpar));
    if(MYID(simpar) == 0) Dtime = TDtime;
    MPI_Bcast(&Dtime, 1, MPI_POSTYPE,  0, MPI_COMM(simpar));

    int i;
    postype dmax=0;
//  postype eps2 = KH_EPS(simpar)*KH_EPS(simpar);
#ifdef _OPENMP
#pragma omp parallel for reduction(max: dmax)
#endif
    for(i=0;i<VORO_NP(simpar);i++){
//		if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
		if(targetBP(bp+i, Lx, Ly))
		{
			bp[i].ax += GAS_ACCX(simpar);
			bp[i].ay += GAS_ACCY(simpar);
	        bp[i].vx += bp[i].ax *Dtime; 
			bp[i].vy += bp[i].ay *Dtime; 
			bp[i].x += bp[i].vx *Dtime; 
			bp[i].y += bp[i].vy *Dtime; 
			dmax = MAX(dmax, fabs(bp[i].ay)); 
			bp[i].x = fmod(bp[i].x+Lx,Lx); 
			bp[i].y = fmod(bp[i].y+Ly,Ly); 
			bp[i].w2old = bp[i].w2;
		}
		else {
			bp[i].vx = bp[i].vy = bp[i].ax = bp[i].ay = 0;
		}
    }
    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);

// This is for centroid shift.
	if(GAS_FCENTROID(simpar) > 0 ){
	    if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingAllTreeParticles);
    // Update the linked list 
		mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

    // repositioning the memory space 
	    bp = VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	    for(iy=1;iy<my-1;iy++){ 
			int mp=2000; 
			Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp); 
			postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a; 
			int ix; 
			for(ix=1;ix<mx-1;ix++){ 
				int np; 
				treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np); 
				int nneigh; 
				Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh); 
				Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh); 
				int i; 
				for(i=0;i<np;i++){ 
					Voro2D_point center; 
					center.x = p[i].x; 
					center.y = p[i].y; 
					center.indx = PINDX(p+i); 
					center.w2 = p[i].w2; 
					
					treevorork4particletype *ibp = p[i].bp; 
					int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
					// in this case ax, ay are used for the centroid vector
					Voro2D_Corner centroid = Centroid2DPolygon(vorocorner);
					ibp->ax = centroid.x;
					ibp->ay = centroid.y; 
				} 
				my_free(p);my_free(neighbors); my_free(neighwork); 
			} 
			my_free(vorocorner); 
		} 
		my_free(VORORK4_TBPP(simpar));
		float fshift = GAS_FCENTROID(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
	    for(i=0;i<VORO_NP(simpar);i++){ 
//			if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
			if(targetBP(bp+i, Lx, Ly))
			{
				bp[i].x += fshift*bp[i].ax; // move it to the centroid  with a factor of shift
				bp[i].y += fshift*bp[i].ay; 
   		     	bp[i].x = fmod(bp[i].x+Lx,Lx);
		        bp[i].y = fmod(bp[i].y+Ly,Ly);
			}
		} 
		// migrate particles between mpi ranks 
		migrateTreeVorork4Particles(simpar);
	}
//
//
	if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingAllTreeParticles);
    // Update the linked list
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);


    // repositioning the memory space
    bp = VORORK4_TBP(simpar);



#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)my_malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)my_malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.w2 = p[i].w2;

                treevorork4particletype *ibp = p[i].bp;
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
//              ibp->w2 = ibp->w2ceil = center.w2ceil;
                // These are to use register 
                postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
                ibp_vx = ibp->vx;
                ibp_vy = ibp->vy;
                ibp_csound = ibp->csound;
                ibp_pressure = ibp->pressure;
                ibp_den = ibp->den;

                ibp->volume = Area2DPolygon(vorocorner, mp);




                Voro2D_Corner *tmp,*tmp2;
                tmp = vorocorner;
                int j;
//                double die, dte,dke,fx,fy;
//               die = dte = dke = fx = fy = 0;
				double die, cx,cy;
				die = cx = cy = 0;
                do {
                    if(tmp->upperrelated >=0){
                        // These are to use register 
                        treevorork4particletype *jbp =  (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
                        postype jbp_vx = jbp->vx;
                        postype jbp_vy = jbp->vy;
                        postype jbp_csound = jbp->csound;
//                        postype jbp_pressure = jbp->pressure;
						postype jbp_pressure = getNeighborPressure(jbp,ibp);
                        postype jbp_den = jbp->den;

                        Voro2D_point line;
                        tmp2 = tmp->upperlink;
                        line.x = tmp2->x - tmp->x;
                        line.y = tmp2->y - tmp->y;
                        Voro2D_point dS = voro2D_norm(&line);
                        postype facearea = Vec2DLength(tmp,tmp2);
                        dS.x = facearea*dS.x;
                        dS.y = facearea*dS.y;
//                      postype pi = Half*(jbp_pressure + ibp_pressure);
                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,RT);

                        Voro2D_point ui,ua,ub;
//                        ui.x = Half*(jbp_vx - ibp_vx);
 //                       ui.y = Half*(jbp_vy - ibp_vy);
                        ui = get2dUpqradRk4(ibp,jbp,dtold);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &ui);
                        if(rvel<0){
                             postype mu = rvel/(dramp/etavis + eps2vis*etavis/dramp);
                             postype meanden = 0.5*(ibp_den + jbp_den);
                             postype meanCsound = 0.5*(ibp_csound + jbp_csound);
                             pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy 
                        die += -pi * Vec2DDotP(&ui,&dS);

                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
//                ibp->ax = fx/ibp->mass;
//                ibp->ay = fy/ibp->mass;
            }
            my_free(p);my_free(neighbors); my_free(neighwork);
        }
        my_free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));
    my_free(VORO_BASICCELL(simpar));

    // update the w2ceil
    postype dmean = SIMBOX(simpar).x.max/NX(simpar);
//  DEBUGPRINT("P%d is now updating the internal energy\n", MYID(simpar));
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
//		if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
		if(targetBP(bp+i, Lx, Ly))
		{
	        bp[i].ie += bp[i].die * Dtime; 
			bp[i].den  = bp[i].mass/bp[i].volume; 
			bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1); 
			bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den); 
			if(GAS_Kappa(simpar) <0) { 
				bp[i].w2 = -GAS_Kappa(simpar); 
			} 
			else if (GAS_Kappa(simpar) >0){ 
//				bp[i].w2 = measureW2(simpar,dmean, bp[i].pressure, bp[i].den); 
          		bp[i].w2 = getw2forHydroParticle(simpar,(bp+i));
				bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil); 
			} 
			if(bp[i].volume ==0) 
				DEBUGPRINT("P%d has wrong volume at i= %d with den= %g at xy= %g %g\n", MYID(simpar), 
						i,bp[i].den, bp[i].x, bp[i].y); 
			if(isnan(bp[i].csound)){ 
				DEBUGPRINT("P%d has nan at the sound speed = %d with" 
						"id= %ld: pressure= %g den= %g Delta IE=%g dIE= %g Dt= %g\n", 
						MYID(simpar), i,PINDX(bp+i), bp[i].pressure, bp[i].den, bp[i].ie,bp[i].die,Dtime); 
				exit(9);
			}
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return Dtime;
}

double exam2d_vph_rk4_int(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	postype xmin,ymin,xmax,ymax,OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);
	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has x/y minmax= %g %g : %g %g ::: %g\n", MYID(simpar), xmin,ymin,xmax,ymax, cellsize);


    int i;
    postype Dtime,dt;
    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) bp[i].rk4.w2backup = bp[i].w2;
	// Runge-Kutta 4-th order time evolution of r and vr
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D);
	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed updateDenW2P...()\n", MYID(simpar));
    Dtime = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){ 
		bp[i].rk4.k1x = bp[i].vx*Dtime;
        bp[i].rk4.k1y = bp[i].vy*Dtime;
        bp[i].rk4.k1vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k1vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k1ie = bp[i].die*Dtime;
    }
    // RK4 second
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += bp[i].rk4.k1x*0.5;
        bp[i].y += bp[i].rk4.k1y*0.5;
        bp[i].vx += bp[i].rk4.k1vx*0.5;
        bp[i].vy += bp[i].rk4.k1vy*0.5;
        bp[i].ie += bp[i].rk4.k1ie*0.5;
        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
	migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D);
	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed updateDenW2P...()\n", MYID(simpar));
    dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
        paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k2x = bp[i].vx*Dtime;
        bp[i].rk4.k2y = bp[i].vy*Dtime;
        bp[i].rk4.k2vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k2vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k2ie = bp[i].die*Dtime;
    }
    // RK4 third
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += (bp[i].rk4.k2x-bp[i].rk4.k1x)*0.5;
        bp[i].y += (bp[i].rk4.k2y-bp[i].rk4.k1y)*0.5;
        bp[i].vx += (bp[i].rk4.k2vx-bp[i].rk4.k1vx)*0.5;
        bp[i].vy += (bp[i].rk4.k2vy-bp[i].rk4.k1vy)*0.5;
        bp[i].ie += (bp[i].rk4.k2ie-bp[i].rk4.k1ie)*0.5;

        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D);
	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed updateDenW2P...()\n", MYID(simpar));
	dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k3x = bp[i].vx*Dtime;
        bp[i].rk4.k3y = bp[i].vy*Dtime;
        bp[i].rk4.k3vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k3vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k3ie = bp[i].die*Dtime;
    }
    // RK4 fourth
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += bp[i].rk4.k3x - bp[i].rk4.k2x*0.5; 
		bp[i].y += bp[i].rk4.k3y - bp[i].rk4.k2y*0.5; 
		bp[i].vx += bp[i].rk4.k3vx - bp[i].rk4.k2vx*0.5; 
		bp[i].vy += bp[i].rk4.k3vy - bp[i].rk4.k2vy*0.5; 
		bp[i].ie += bp[i].rk4.k3ie - bp[i].rk4.k2ie*0.5;

        bp[i].x = fmod(bp[i].x+Lx, Lx); 
		bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D);
	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed updateDenW2P...()\n", MYID(simpar));
    dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k4x = bp[i].vx*Dtime;
        bp[i].rk4.k4y = bp[i].vy*Dtime;
        bp[i].rk4.k4vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k4vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k4ie = bp[i].die*Dtime;
    }
    // migrate particles between mpi ranks
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x -= bp[i].rk4.k3x;
        bp[i].y -= bp[i].rk4.k3y;
        bp[i].vx -= bp[i].rk4.k3vx;
        bp[i].vy -= bp[i].rk4.k3vy;
        bp[i].ie -= bp[i].rk4.k3ie;
        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
	// migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

    // to finalize ke & ie updates 
	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
    // Update position and velocity using RK4
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
//		if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
		if(targetBP(bp+i, Lx, Ly))
		{
			bp[i].x  += (bp[i].rk4.k1x +2*bp[i].rk4.k2x +2*bp[i].rk4.k3x +bp[i].rk4.k4x )/6.; 
			bp[i].y  += (bp[i].rk4.k1y +2*bp[i].rk4.k2y +2*bp[i].rk4.k3y +bp[i].rk4.k4y )/6.; 
			bp[i].vx += (bp[i].rk4.k1vx+2*bp[i].rk4.k2vx+2*bp[i].rk4.k3vx+bp[i].rk4.k4vx)/6.; 
			bp[i].vy += (bp[i].rk4.k1vy+2*bp[i].rk4.k2vy+2*bp[i].rk4.k3vy+bp[i].rk4.k4vy)/6.; 
			bp[i].ie += (bp[i].rk4.k1ie+2*bp[i].rk4.k2ie+2*bp[i].rk4.k3ie+bp[i].rk4.k4ie)/6.;

	        bp[i].x = fmod(bp[i].x+Lx,Lx); 
			bp[i].y = fmod(bp[i].y+Ly,Ly);
		}
    }

    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);

	//centroid shift
	if(GAS_FCENTROID(simpar)>0) 
		exam2d_centroidShift(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	// update volume & density with one-step evolved positions 
    exam2dUpdateVol(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);
	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed exam2dUpdateVol...()\n", MYID(simpar));

	bp = VORORK4_TBP(simpar);
	postype dmean = SIMBOX(simpar).x.max/NX(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
//		if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
		if(targetBP(bp+i, Lx, Ly))
		{
			bp[i].den = bp[i].mass/bp[i].volume; 
			bp[i].pressure = bp[i].ie/bp[i].volume * (Gamma-1); 
			bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den); 
			if(GAS_Kappa(simpar) <0) { 
				bp[i].w2 = -GAS_Kappa(simpar); 
			} 
			else if (GAS_Kappa(simpar) >0){ 
//				bp[i].w2 = measureW2(simpar,dmean, bp[i].pressure, bp[i].den); 
          		bp[i].w2 = getw2forHydroParticle(simpar,(bp+i));
				bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil); 
			} 
			if(bp[i].volume ==0) 
				DEBUGPRINT("P%d has wrong volume at i= %d with den= %g at xy= %g %g\n", MYID(simpar), 
						i,bp[i].den, bp[i].x, bp[i].y); 
			if(isnan(bp[i].csound)){ 
				DEBUGPRINT("P%d has nan at the sound speed = %d with" 
						"id= %zu: pressure= %g den= %g Delta IE=%g dIE= %g Dt= %g\n", 
						MYID(simpar), i,PINDX(bp+i), bp[i].pressure, bp[i].den, bp[i].ie,bp[i].die,Dtime); 
				exit(9);
			}
		}
		else {
			bp[i].vx = bp[i].vy = 0;
			bp[i].w2= bp[i].rk4.w2backup;
		}
    }
    return Dtime;
}

double exam2d_vph_rk4_int_rt(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
//			void pddingAllTreeParticles(SimParameters *, postype))
			void (*)(SimParameters *, postype))
		){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	postype xmin,ymin,xmax,ymax,OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);
//	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has x/y minmax= %g %g : %g %g ::: %g\n", MYID(simpar), xmin,ymin,xmax,ymax, cellsize);


    int i;
    postype Dtime,dt;
    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) bp[i].rk4.w2backup = bp[i].w2;
	// Runge-Kutta 4-th order time evolution of r and vr
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D);
//	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed updateDenW2P...()\n", MYID(simpar));
    Dtime = getAccVoro2D_rt(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){ 
		bp[i].rk4.k1x = bp[i].vx*Dtime;
        bp[i].rk4.k1y = bp[i].vy*Dtime;
        bp[i].rk4.k1vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k1vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k1ie = bp[i].die*Dtime;
    }
    // RK4 second
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += bp[i].rk4.k1x*0.5;
        bp[i].y += bp[i].rk4.k1y*0.5;
        bp[i].vx += bp[i].rk4.k1vx*0.5;
        bp[i].vy += bp[i].rk4.k1vy*0.5;
        bp[i].ie += bp[i].rk4.k1ie*0.5;
        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
	migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D);
    dt = getAccVoro2D_rt(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
        paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k2x = bp[i].vx*Dtime;
        bp[i].rk4.k2y = bp[i].vy*Dtime;
        bp[i].rk4.k2vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k2vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k2ie = bp[i].die*Dtime;
    }
    // RK4 third
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += (bp[i].rk4.k2x-bp[i].rk4.k1x)*0.5;
        bp[i].y += (bp[i].rk4.k2y-bp[i].rk4.k1y)*0.5;
        bp[i].vx += (bp[i].rk4.k2vx-bp[i].rk4.k1vx)*0.5;
        bp[i].vy += (bp[i].rk4.k2vy-bp[i].rk4.k1vy)*0.5;
        bp[i].ie += (bp[i].rk4.k2ie-bp[i].rk4.k1ie)*0.5;

        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D);
	dt = getAccVoro2D_rt(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k3x = bp[i].vx*Dtime;
        bp[i].rk4.k3y = bp[i].vy*Dtime;
        bp[i].rk4.k3vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k3vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k3ie = bp[i].die*Dtime;
    }
    // RK4 fourth
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x += bp[i].rk4.k3x - bp[i].rk4.k2x*0.5; 
		bp[i].y += bp[i].rk4.k3y - bp[i].rk4.k2y*0.5; 
		bp[i].vx += bp[i].rk4.k3vx - bp[i].rk4.k2vx*0.5; 
		bp[i].vy += bp[i].rk4.k3vy - bp[i].rk4.k2vy*0.5; 
		bp[i].ie += bp[i].rk4.k3ie - bp[i].rk4.k2ie*0.5;

        bp[i].x = fmod(bp[i].x+Lx, Lx); 
		bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D);
    dt = getAccVoro2D_rt(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP,mkLinkedList2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].rk4.k4x = bp[i].vx*Dtime;
        bp[i].rk4.k4y = bp[i].vy*Dtime;
        bp[i].rk4.k4vx = (bp[i].ax + GAS_ACCX(simpar))*Dtime;
        bp[i].rk4.k4vy = (bp[i].ay + GAS_ACCY(simpar))*Dtime;
        bp[i].rk4.k4ie = bp[i].die*Dtime;
    }
    // migrate particles between mpi ranks
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].x -= bp[i].rk4.k3x;
        bp[i].y -= bp[i].rk4.k3y;
        bp[i].vx -= bp[i].rk4.k3vx;
        bp[i].vy -= bp[i].rk4.k3vy;
        bp[i].ie -= bp[i].rk4.k3ie;
        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
	// migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

    // to finalize ke & ie updates 
	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
    // Update position and velocity using RK4
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly)) {
			bp[i].x  += (bp[i].rk4.k1x +2*bp[i].rk4.k2x +2*bp[i].rk4.k3x +bp[i].rk4.k4x )/6.; 
			bp[i].y  += (bp[i].rk4.k1y +2*bp[i].rk4.k2y +2*bp[i].rk4.k3y +bp[i].rk4.k4y )/6.; 
			bp[i].vx += (bp[i].rk4.k1vx+2*bp[i].rk4.k2vx+2*bp[i].rk4.k3vx+bp[i].rk4.k4vx)/6.; 
			bp[i].vy += (bp[i].rk4.k1vy+2*bp[i].rk4.k2vy+2*bp[i].rk4.k3vy+bp[i].rk4.k4vy)/6.; 
			bp[i].ie += (bp[i].rk4.k1ie+2*bp[i].rk4.k2ie+2*bp[i].rk4.k3ie+bp[i].rk4.k4ie)/6.;
			if(bp[i].ie<0) bp[i].ie = 1.e-5;
	        bp[i].x = fmod(bp[i].x+Lx,Lx); 
			bp[i].y = fmod(bp[i].y+Ly,Ly);
		}
    }

    // migrate particles between mpi ranks 
    migrateTreeVorork4Particles(simpar);

	//centroid shift
	if(GAS_FCENTROID(simpar)>0) 
		exam2d_centroidShift(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	// update volume & density with one-step evolved positions 
    exam2dUpdateVol(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	bp = VORORK4_TBP(simpar);
	postype dmean = SIMBOX(simpar).x.max/NX(simpar);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
			bp[i].den = bp[i].mass/bp[i].volume; 
			bp[i].pressure = bp[i].ie/bp[i].volume * (Gamma-1); 
			bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den); 
			if(GAS_Kappa(simpar) <0) { 
				bp[i].w2 = -GAS_Kappa(simpar); 
			} 
			else if (GAS_Kappa(simpar) >0){ 
//  				bp[i].w2 = measureW2(simpar,dmean, bp[i].pressure, bp[i].den); 
          		bp[i].w2 = getw2forHydroParticle(simpar,(bp+i));
				bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil); 
			} 
			if(bp[i].volume ==0) 
				DEBUGPRINT("P%d has wrong volume at i= %d with den= %g at xy= %g %g\n", MYID(simpar), 
						i,bp[i].den, bp[i].x, bp[i].y); 
			if(isnan(bp[i].csound)){ 
				DEBUGPRINT("P%d has nan at the sound speed = %d with" 
						"id= %ld: pressure= %g den= %g Delta IE=%g dIE= %g Dt= %g\n", 
						MYID(simpar), i,PINDX(bp+i), bp[i].pressure, bp[i].den, bp[i].ie,bp[i].die,Dtime); 
				exit(9);
			}
		}
		else {
			bp[i].vx = bp[i].vy = 0;
			bp[i].w2= bp[i].rk4.w2backup;
		}
    }
    return Dtime;
}
