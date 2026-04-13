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

#ifdef USE_CUDA
#include "exam_gpu.h"
/* Forward declaration — full signature matches exam_gpu_extract.c */
double getAccVoro2DBlend_GPU(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    postype OrderOfAccuracy, postype Courant, postype Gamma,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
    treevorork4particletype *(*find2DCellBP)(SimParameters *, int, int, int *),
    void (*mkLinkedList2D)(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)));
double getAccVoro2DBlend_GPU_validate(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    postype OrderOfAccuracy, postype Courant, postype Gamma,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
    treevorork4particletype *(*find2DCellBP)(SimParameters *, int, int, int *),
    void (*mkLinkedList2D)(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)));
/* GPU nearest-neighbor (replaces det2d_dpqRK4 for Voronoi av_modes) */
void det2d_dpqRK4_GPU(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    void mkLinkedList2D(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)));
/* LagMFM GPU wrappers (exam_gpu_extract.c) */
void updateDenW2Pressure2D_LagMFM_GPU(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    postype Gamma,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    void mkLinkedList2D(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)),
    postype Dtime);
double getAccVoro2D_LagMFM_GPU(
    SimParameters *simpar,
    postype xmin, postype ymin, postype xmax, postype ymax,
    postype OrderOfAccuracy, postype Courant, postype Gamma,
    void (*paddingAllTreeParticles)(SimParameters *, postype),
    void mkLinkedList2D(SimParameters *, postype, postype, postype, postype, postype,
        void (*)(SimParameters *, postype)));
#endif


inline postype getw2forHydroParticle(SimParameters *simpar, treevorork4particletype *bp,
		postype Dtime){
//	postype w1old = sqrt(bp->w2old);
	postype gamma = GAS_GAMMA(simpar);
	int ndim = NDIM(simpar);
	if(1){
		postype w2new= GAS_Kappa(simpar)*GAS_Kappa(simpar)*
			GAS_dMean(simpar)*GAS_dMean(simpar)
			*pow(bp->pressure*GAS_invw2Scale(simpar), GAS_w2Power(simpar));
		return w2new;
	}
	else {
		postype w1 = sqrt(bp->w2);
		postype exponent = 1./(ndim*gamma);
		postype dw = GAS_w2Power(simpar)*(pow(bp->pressure/bp->avgNeighboringPressure, exponent)-1)*w1;

		postype shift = bp->csound*Dtime;
		if(PINDX(bp)==16384) 
			DEBUGPRINT("P%d has w1= %g, exponent= %g, dw= %g, shift= %g, Cs= %g Dtime= %g\n", 
					MYID(simpar), w1, exponent, dw, shift, bp->csound, Dtime); 
		if(dw > 0) dw = MIN(dw, shift);
		else dw = MAX(dw, -shift);
		w1 = w1+dw;
		postype w2 = w1*w1;
		return w2;
	}
}

/* Apply w2 control mechanisms: floor, relaxation, rate limiter, ceiling.
 * In 2D, w2ceil is a hard geometric constraint (cell degenerates if exceeded),
 * so ceiling MUST be applied last — rate limiter must not push w2 above w2ceil.
 * (This differs from 1D where ceiling is soft and rate limiter comes after.)
 * Call AFTER getw2forHydroParticle sets bp->w2hydro. */
static inline void applyW2Controls(SimParameters *simpar, treevorork4particletype *bp, postype Dtime){
	postype w2 = bp->w2hydro;

	/* 1. Floor: w >= w2_floor_frac * dMean */
	postype w2ff = GAS_W2FLOORFRAC(simpar);
	if(w2ff > 0){
		postype wf = w2ff * GAS_dMean(simpar);
		postype w2floor = wf * wf;
		if(w2 < w2floor) w2 = w2floor;
	}

	/* 2. Relaxation: exponential smoothing toward target */
	postype w2rt = GAS_W2RELAXTAU(simpar);
	if(w2rt > 0 && Dtime > 0){
		postype ci = fmax(bp->csound, 1e-10);
		postype tau = w2rt * GAS_dMean(simpar) / ci;
		postype alpha_r = fmin(1.0, Dtime / tau);
		postype w2_old = (bp->w2old > 0) ? bp->w2old : w2;
		w2 = w2_old + alpha_r * (w2 - w2_old);
	}

	/* 3. Rate limiter: |Δw2/w2_old| ≤ w2_rate_max */
	postype w2rm = GAS_W2RATEMAX(simpar);
	if(w2rm > 0 && bp->w2old > 0){
		postype w2_hi = bp->w2old * (1.0 + w2rm);
		postype w2_lo = bp->w2old * (1.0 - w2rm);
		if(w2_lo < 0) w2_lo = 0;
		if(w2 > w2_hi) w2 = w2_hi;
		if(w2 < w2_lo) w2 = w2_lo;
	}

	bp->w2hydro = w2;

	/* 4. Ceiling: w2 ≤ w2ceil (MUST be last — hard geometric constraint in 2D) */
	bp->w2 = MIN(w2, bp->w2ceil);
}

// *_w2Masure2D are obsolete and instead use getw2forHydroParticle.
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

double cyl_w2Measure2D(SimParameters *simpar, double dmean, double pressure,double rho){
#define p0_w2 71.42857
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
    PosType distmax = -1.e20;
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
    PosType distmax = -1.e20;
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

/* 5-argument nearOpen callback for find_GNear k-NN search */
int near2dOpen_kNN(void *p, TStruct *tree, int npneigh, PosType maxdist, int K){
	if(npneigh >= K){
		TPtlStruct *pt = (TPtlStruct*)p;
		PosType dx = pt->x - tree->monox, dy = pt->y - tree->monoy;
		if(sqrt(dx*dx + dy*dy) - tree->nodesize > maxdist) return NO;
	}
	return YES;
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
	TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode*2);

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
    TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode*2);
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
	/* Stride-aware: p_size may differ from sizeof(treevorork4particletype)
	   when stress particles are active (av_mode >= 1). */
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
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
	TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode*2);

	/* Refresh after padding (may realloc) */
	bp_raw = (char*)VORORK4_TBP(simpar);

	postype xmin,ymin,xmax,ymax;
	xmin = ymin = 1.e20;
	xmax = ymax = -1.e20;
	// this is for the local-domain particles.
	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		ptl[i].x = bpi->x;
		ptl[i].y = bpi->y;
		ptl[i].mass = bpi->mass; /* in this case mass */
		ptl[i].sibling = ptl + i+1;
		ptl[i].type = TYPE_PTL;
		ptl[i].indx = i;
		xmin = MIN(xmin, bpi->x);
		ymin = MIN(ymin, bpi->y);
		xmax = MAX(xmax, bpi->x);
		ymax = MAX(ymax, bpi->y);
	}

	// this is for the boundary ghost particles.
	char *pp_raw = (char*)VORORK4_TBPP(simpar);
	for(i=0;i<mp;i++){
		treevorork4particletype *ppi = (treevorork4particletype*)(pp_raw + i*p_size);
		ptl[np+i].x = ppi->x;
		ptl[np+i].y = ppi->y;
		ptl[np+i].mass = ppi->mass; /* in this case mass */
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
	bp_raw = (char*)VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		bpi->w2ceil = find_GNearest(ptl+i,  tree, nearest2dOpen, ex2d_dist);
		if(GAS_Kappa(simpar) >= 0)
			bpi->w2 = MIN(bpi->w2, bpi->w2ceil);
	}

	my_free(ptl);
	my_free(tree);
	my_free(VORORK4_TBPP(simpar));
	if(0){
		postype mindpq,maxdpq;
		mindpq = 1.e20;
		maxdpq = -1.e20;
		for(i=0;i<np;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			mindpq = MIN(mindpq, bpi->w2ceil);
			maxdpq = MAX(maxdpq, bpi->w2ceil);
		}
		DEBUGPRINT("P%d passed the finding of dpq: min/max= %g %g\n",MYID(simpar), mindpq,maxdpq);
	}
}

/* Build k-d tree for k-NN search. Tree, ptl, and padding are kept alive
   for the caller (updateDen_kNN + getAcc_kNN) to reuse.
   Caller is responsible for freeing tree_out, ptl_out, and VORORK4_TBPP. */
void buildTree2D(SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		int (*skipBP)(treevorork4particletype *),
		TStruct **tree_out, TPtlStruct **ptl_out, int *nptl_out){
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;

	int mp;
	{
		postype cellsize = HydroGridSize(simpar);
		paddingAllTreeParticles(simpar, cellsize);
		mp = VORO_NPAD(simpar);
	}

	DEBUGPRINT("P%d buildTree2D_kNN: after padding np=%d mp=%d\n",
		MYID(simpar), np, mp);

	/* Refresh after padding (may realloc) */
	bp_raw = (char*)VORORK4_TBP(simpar);

	TPtlStruct *ptl = (TPtlStruct*)my_malloc(sizeof(TPtlStruct)*(np+mp));
	int nnode = np+mp;
	TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode*2);

	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		ptl[i].x = bpi->x;
		ptl[i].y = bpi->y;
		ptl[i].mass = bpi->mass;
		ptl[i].sibling = ptl + i+1;
		ptl[i].type = TYPE_PTL;
		ptl[i].indx = i;
		ptl[i].bp = (linkedlisttype*)bpi;
	}

	/* Boundary ghost particles */
	char *pp_raw = (char*)VORORK4_TBPP(simpar);
	for(i=0;i<mp;i++){
		treevorork4particletype *ppi = (treevorork4particletype*)(pp_raw + i*p_size);
		ptl[np+i].x = ppi->x;
		ptl[np+i].y = ppi->y;
		ptl[np+i].mass = ppi->mass;
		ptl[np+i].sibling = ptl + np+i+1;
		ptl[np+i].type = TYPE_PTL;
		ptl[np+i].indx = np+i;
		ptl[np+i].bp = (linkedlisttype*)ppi; /* back-reference to padding particle */
	}
	ptl[np+mp-1].sibling = NULL;
	tree->daughter = ptl;
	tree->Nparticle = np+mp;

	Make_GNN_Tree(tree, nnode, ex2d_idivision, ex2d_findCentroid,
			ex2d_findCellSize, RECURSIVE);

	/* Compute w2ceil for local particles */
	bp_raw = (char*)VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		if(skipBP && skipBP(bpi)) continue;
		bpi->w2ceil = find_GNearest(ptl+i, tree, nearest2dOpen, ex2d_dist);
		if(GAS_Kappa(simpar) >= 0)
			bpi->w2 = MIN(bpi->w2, bpi->w2ceil);
	}

	*tree_out = tree;
	*ptl_out = ptl;
	*nptl_out = np + mp;
}

/* k-NN search around a single query point.
   Returns Voro2D_point neighbor array (malloc'd) and count.
   maxr is set to the distance of the K-th neighbor. */
int searchKNN2D(TPtlStruct *query, int K, TStruct *tree,
		Voro2D_point *neighbors, PosType *maxr){
	void *bpneighbor[MAX_NUM_NEAR];
	int nfound;

	nfound = find_GNear(query, K, tree, maxr, bpneighbor,
			near2dOpen_kNN, ex2d_dist, insGnear);

	/* Debug: if maxr=0, do a manual tree traversal with counters */
	{
		static int dbg_maxr0 = 0;
		if(*maxr < 1e-10 && nfound >= K && dbg_maxr0 < 1){
			dbg_maxr0++;
			fprintf(stderr, "searchKNN2D maxr=0: query=(%g,%g) K=%d nfound=%d maxr=%g\n",
				query->x, query->y, K, nfound, *maxr);
			/* Manual traversal mimicking find_GNear */
			int n_tree=0, n_ptl=0, n_open=0, n_skip=0, n_insert=0;
			int n_dist0=0, n_same_ptr=0;
			void *dbg_ptr = (void*)tree;
			PosType dbg_maxdist2 = 1.E23, dbg_maxdist = 1.E23;
			int dbg_npneigh = 0;
			void *first_ptl_ptr = NULL;
			while(dbg_ptr != NULL){
				if(((TYPE*)dbg_ptr)->type == TYPE_TREE){
					n_tree++;
					/* Always open (same as near2dOpen_kNN with nodesize=3e10) */
					int open_res = near2dOpen_kNN(query, (TStruct*)dbg_ptr,
						dbg_npneigh, dbg_maxdist, K);
					if(open_res == YES){
						n_open++;
						dbg_ptr = (void*)(((TStruct*)dbg_ptr)->daughter);
					} else {
						n_skip++;
						dbg_ptr = (void*)(((TStruct*)dbg_ptr)->sibling);
					}
				} else {
					n_ptl++;
					PosType d2 = ex2d_dist(query, dbg_ptr);
					if(d2 < 1e-20) n_dist0++;
					if(dbg_ptr == (void*)query) n_same_ptr++;
					if(first_ptl_ptr == NULL) first_ptl_ptr = dbg_ptr;
					/* Simulate insertion logic */
					if(dbg_npneigh < K || d2 < dbg_maxdist2){
						n_insert++;
						dbg_npneigh++;
						if(dbg_npneigh > K) dbg_npneigh = K;
						/* Track worst case maxdist2 */
						if(d2 > 0 && d2 < dbg_maxdist2) dbg_maxdist2 = d2;
						if(d2 == 0 && dbg_npneigh == 1) dbg_maxdist2 = 0;
					}
					dbg_maxdist = sqrt(dbg_maxdist2);
					dbg_ptr = (void*)(((TPtlStruct*)dbg_ptr)->sibling);
				}
				if(n_ptl + n_tree > 50000) {
					fprintf(stderr, "  ABORT: too many visits\n");
					break;
				}
			}
			fprintf(stderr, "  Manual traversal: tree_nodes=%d ptl_nodes=%d opens=%d skips=%d\n",
				n_tree, n_ptl, n_open, n_skip);
			fprintf(stderr, "  inserts=%d dist0=%d same_ptr_as_query=%d first_ptl=%p query=%p\n",
				n_insert, n_dist0, n_same_ptr, first_ptl_ptr, (void*)query);
			/* Print first few bpneighbor from actual find_GNear */
			int di;
			for(di=0;di<3 && di<nfound;di++){
				TPtlStruct *nb2 = (TPtlStruct*)bpneighbor[di];
				fprintf(stderr, "  bpneigh[%d]: pos=(%g,%g) bp=%p\n",
					di, nb2->x, nb2->y, (void*)nb2);
			}
		}
	}

	int i;
	for(i=0;i<nfound;i++){
		TPtlStruct *nb = (TPtlStruct*)bpneighbor[i];
		treevorork4particletype *tt = (treevorork4particletype*)(nb->bp);
		neighbors[i].x = tt->x;
		neighbors[i].y = tt->y;
		neighbors[i].indx = PINDX(tt);
		neighbors[i].csound = tt->csound;
		neighbors[i].pressure = tt->pressure;
		neighbors[i].w2 = tt->w2;
		neighbors[i].bp = (void*)tt;
		neighbors[i].u4if.indx = MAX_INDEX;
	}
	return nfound;
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
    TStruct *tree = (TStruct *)my_malloc(sizeof(TStruct)*nnode*2);
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
    neigh = (Voro2D_point*)malloc(sizeof(Voro2D_point)*np);
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
    neigh = realloc(neigh,sizeof(treevorork4particletype)*np);
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
    res = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*np);
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
    res = realloc(res,sizeof(treevorork4particletype)*np);
    return res;
}

/*
void new_mkLinkedList2D_rt(
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
        if(index >=mx*my || index <0){
            DEBUGPRINT("error detected: %d : %d %d : %g %g %g %g\n",
                    index, mx,my, bp[i].x, bp[i].y, xmin,ymin);
            continue;
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
*/

void mkLinkedList2D_rt(
		SimParameters *simpar,
		postype cellsize,
		postype xmin,
		postype ymin,
		postype xmax,
		postype ymax,
		void (*paddingTreeAllParticle)(SimParameters *, postype)
		){
	/* p_size: actual byte size per particle (handles treevorostressrk4particletype
	   when av_mode >= 1, where n_size is overridden in RunRT) */
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
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
	/* Build linked list for local particles (stride-aware) */
	for(i=0;i<np;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
        int ix,iy;
        ix = ((bpi->x-xmin)/cellsize);
        iy = ((bpi->y-ymin)/cellsize);
        int index = ix+mx*iy;
        if(index >=mx*my || index <0){
            DEBUGPRINT("error detected: %d : %d %d : %g %g %g %g\n",
                    index, mx,my, bpi->x, bpi->y, xmin,ymin);
            continue;
        }
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)bpi;
        cells[index].nmem ++;
        bpi->next = tmp;
    }
	// now for the particle padding from other ranks
    paddingTreeAllParticle(simpar, cellsize);

	postype getRT_Den(postype , postype , postype , postype ,postype );
	int ny = NY(simpar);

	/* Top wall mirror (y = Ly): compact padding, create mirrors
	   Stride-aware version of the original algorithm. */
	if(Ymax_HydroExam(simpar) >= (Ly - EPS)){
		char *bpp_raw = (char*)VORORK4_TBPP(simpar);
		int mpad = 0;
		for(i=0;i<VORO_NPAD(simpar);i++){
			treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
			if(bppi->y < Ly) {
				if(mpad != i) memcpy(bpp_raw + mpad*p_size, bpp_raw + i*p_size, p_size);
				mpad ++;
			}
		}
		int nmirror = 0;
		for(i=0;i<np;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			if((Ly-bpi->y) < cellsize) nmirror++;
		}
		for(i=0;i<mpad;i++){
			treevorork4particletype *bppi = (treevorork4particletype*)((char*)VORORK4_TBPP(simpar) + i*p_size);
			if((Ly-bppi->y) < cellsize) nmirror++;
		}
		if(VORO_NPAD(simpar) != (mpad+nmirror)) {
			VORORK4_TBPP(simpar) = (treevorork4particletype *)
					realloc(VORORK4_TBPP(simpar), p_size*(mpad+nmirror));
		}
		bpp_raw = (char*)VORORK4_TBPP(simpar);
		char *cursor = bpp_raw + mpad*p_size;
		for(i=0;i<np;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			if((Ly-bpi->y) < cellsize) {
				memcpy(cursor, bpi, p_size);
				treevorork4particletype *gp = (treevorork4particletype*)cursor;
				postype dy = 2*(Ly-bpi->y);
				gp->u4if.indx = MAX_INDEX;
				gp->y += dy;
				gp->vy = -bpi->vy;
				gp->den = getRT_Den(rho1, rho2, deltay, ycen, gp->y);
				gp->pressure += GAS_ACCY(simpar)* getRT_Den(rho1,rho2,deltay,ycen,Ly) * dy;
				gp->ie = gp->pressure*gp->volume/(Gamma-1);
				CLEAR_FLAG(gp);
				SET_FLAG(gp,BoundaryGhostflag);
				SET_FLAG(gp,VOROflag);
				cursor += p_size;
			}
		}
		for(i=0;i<mpad;i++){
			treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
			if((Ly-bppi->y) < cellsize) {
				memcpy(cursor, bppi, p_size);
				treevorork4particletype *gp = (treevorork4particletype*)cursor;
				postype dy = 2*(Ly-bppi->y);
				gp->u4if.indx = MAX_INDEX;
				gp->y += dy;
				gp->vy = -bppi->vy;
				gp->den = getRT_Den(rho1, rho2, deltay, ycen, gp->y);
				gp->pressure += GAS_ACCY(simpar)* getRT_Den(rho1,rho2,deltay,ycen,Ly) * dy;
				gp->ie = gp->pressure*gp->volume/(Gamma-1);
				CLEAR_FLAG(gp);
				SET_FLAG(gp,BoundaryGhostflag);
				SET_FLAG(gp,VOROflag);
				cursor += p_size;
			}
		}
		VORO_NPAD(simpar) = mpad + nmirror;
	}
	/* Bottom wall mirror (y = 0) */
	if(Ymin_HydroExam(simpar)< EPS)
	{
		char *bpp_raw = (char*)VORORK4_TBPP(simpar);
		int mpad = 0;
		for(i=0;i<VORO_NPAD(simpar);i++){
			treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
			if(bppi->y >= 0.) {
				if(mpad != i) memcpy(bpp_raw + mpad*p_size, bpp_raw + i*p_size, p_size);
				mpad ++;
			}
		}
		int nmirror = 0;
		for(i=0;i<np;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			if(bpi->y <= cellsize) nmirror++;
		}
		for(i=0;i<mpad;i++){
			treevorork4particletype *bppi = (treevorork4particletype*)((char*)VORORK4_TBPP(simpar) + i*p_size);
			if(bppi->y <= cellsize) nmirror++;
		}
		if(VORO_NPAD(simpar) != (mpad+nmirror)) {
			VORORK4_TBPP(simpar) = (treevorork4particletype *) realloc(VORORK4_TBPP(simpar),
						p_size*(mpad+nmirror));
		}
		bpp_raw = (char*)VORORK4_TBPP(simpar);
		char *cursor = bpp_raw + mpad*p_size;
		for(i=0;i<np;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			if(bpi->y <= cellsize) {
				memcpy(cursor, bpi, p_size);
				treevorork4particletype *gp = (treevorork4particletype*)cursor;
				gp->u4if.indx = MAX_INDEX;
				postype dy = -bpi->y;
				gp->y = dy;
				gp->vy = -bpi->vy;
				gp->den = getRT_Den(rho1, rho2, deltay, ycen, gp->y);
				gp->pressure -= GAS_ACCY(simpar)* getRT_Den(rho1,rho2,deltay,ycen,0.) * 2*(-dy);
				gp->ie = gp->pressure*gp->volume/(Gamma-1);
				CLEAR_FLAG(gp);
				SET_FLAG(gp,BoundaryGhostflag);
				SET_FLAG(gp,VOROflag);
				cursor += p_size;
			}
		}
		for(i=0;i<mpad;i++){
			treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
			if(bppi->y <= cellsize) {
				memcpy(cursor, bppi, p_size);
				treevorork4particletype *gp = (treevorork4particletype*)cursor;
				gp->u4if.indx = MAX_INDEX;
				postype dy = -bppi->y;
				gp->y = dy;
				gp->vy = -bppi->vy;
				gp->den = getRT_Den(rho1, rho2, deltay, ycen, gp->y);
				gp->pressure -= GAS_ACCY(simpar)* getRT_Den(rho1,rho2,deltay,ycen,0.) * 2*(-dy);
				gp->ie = gp->pressure*gp->volume/(Gamma-1);
				CLEAR_FLAG(gp);
				SET_FLAG(gp,BoundaryGhostflag);
				SET_FLAG(gp,VOROflag);
				cursor += p_size;
			}
		}
		VORO_NPAD(simpar) = mpad + nmirror;
	}
	/* Insert padding + mirror particles into linked list (stride-aware) */
	{
		char *bpp_raw = (char*)VORORK4_TBPP(simpar);
		for(i=0;i<VORO_NPAD(simpar);i++){
			treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
			int ix,iy;
			ix = ((bppi->x-xmin)/cellsize);
			iy = ((bppi->y-ymin)/cellsize);
			int index = ix+mx*iy;
			struct linkedlisttype *tmp = cells[index].link;
			cells[index].link = (struct linkedlisttype*)bppi;
			cells[index].nmem ++;
			bppi->next = tmp;
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
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	int np = VORO_NP(simpar);
    int i,j,k;
    int mx,my;

    mx = BASICCELL_MX(simpar);
    my = BASICCELL_MY(simpar);

	CellType *cells = VORO_BASICCELL(simpar);

	for(i=0;i<mx*my;i++) {
        cells[i].link = NULL;
        cells[i].nmem = 0;
    }
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
        treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
        int ix,iy;
        ix = ((bpi->x-xmin)/cellsize);
        iy = ((bpi->y-ymin)/cellsize);
        int index = ix+mx*iy;
        if(index >=mx*my || index < 0){
            DEBUGPRINT("error detected: %d : %d %d : %g %g %g %g\n",
                    index, mx,my, bpi->x, bpi->y, xmin,ymin);
            continue;
        }
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)bpi;
        cells[index].nmem ++;
        bpi->next = tmp;
    }
	// now for the particle padding from other ranks
    paddingTreeAllParticle(simpar, cellsize);
    char *bpp_raw = (char*)VORORK4_TBPP(simpar);
    for(i=0;i<VORO_NPAD(simpar);i++){
        treevorork4particletype *bppi = (treevorork4particletype*)(bpp_raw + i*p_size);
        int ix,iy;
        ix = ((bppi->x-xmin)/cellsize);
        iy = ((bppi->y-ymin)/cellsize);
        int index = ix+mx*iy;
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)bppi;
        cells[index].nmem ++;
        bppi->next = tmp;
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
			void (*)(SimParameters *, postype)),
		postype Dtime
		){
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int isave = -1;


	det2d_dpqRK4(simpar, paddingAllTreeParticles);
	int i;
	/*
	for(i=0;i<nbp;i++){
		bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
	}
	*/
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)bp;
	for(i=0;i<nbp;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		if(GAS_Kappa(simpar) <0) { 
			bpi->w2 = -GAS_Kappa(simpar); 
		} 
		else if (GAS_Kappa(simpar) >0){ 
			bpi->w2hydro = getw2forHydroParticle(simpar, bpi, Dtime);
			applyW2Controls(simpar, bpi, Dtime);
		}
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


	CellType *cells = (VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*mx*my));
	mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

    int iy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(iy=0;iy<my;iy++){
    	int mp=1000;
    	Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
    	postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
        for(ix=0;ix<mx;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
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
//                ibp->volume = Area2DPolygon(vorocorner, mp);
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);
                ibp->den = ibp->mass/ibp->volume;
				/*
                ibp->csound = sqrt(Gamma*ibp->pressure/ibp->den);
                ibp->pressure = ibp->ie/ibp->volume*(Gamma-1);
				*/
            }
            if(np>0) free(p);
			if(nneigh>0) {
				free(neighbors); free(neighwork);
			}

        }
        free(vorocorner);
    }
	for(i=0;i<nbp;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		bpi->pressure = bpi->ie/bpi->volume*(Gamma-1);
		bpi->csound = sqrt(Gamma*bpi->pressure/bpi->den);
	}
	// finalize mkLinkedList2D by my_free memory spaces (cell & boundary ghost particles)
	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
}

/* ================================================================
 *  updateDenW2Pressure2D_kNN:
 *  k-NN based density/pressure computation. No cell grid needed.
 *  Uses pre-built tree from buildTree2D.
 * ================================================================ */
void updateDenW2Pressure2D_kNN(SimParameters *simpar,
		postype xmin, postype ymin, postype xmax, postype ymax, postype Gamma,
		TStruct *tree, TPtlStruct *ptl, int nptl, postype Dtime,
		int (*skipBP)(treevorork4particletype *)){
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	int nbp = VORO_NP(simpar);
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int i;

	/* Debug: check ie at function entry (disabled for performance) */

	/* Update w2 from kappa setting */
	for(i=0;i<nbp;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		if(GAS_Kappa(simpar) <0){
			bpi->w2 = -GAS_Kappa(simpar);
		}
		else if(GAS_Kappa(simpar) >0){
			bpi->w2hydro = getw2forHydroParticle(simpar, bpi, Dtime);
			applyW2Controls(simpar, bpi, Dtime);
		}
	}

	/* Debug: ie after w2 loop (disabled for performance) */

	int K_START = 32;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<nbp;i++){
		int skip_this = 0;
		if(skipBP){
			treevorork4particletype *bpi_chk = (treevorork4particletype*)(bp_raw + i*p_size);
			if(skipBP(bpi_chk)) skip_this = 1;
		}
		if(!skip_this){

		int mp = 1000;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		Voro2D_point neighbors[MAX_NUM_NEAR];
		Voro2D_point neighwork[MAX_NUM_NEAR];
		memset(neighwork, 0, sizeof(neighwork));
		PosType maxr;
		int K = K_START;
		int nfound, ip;

		/* Adaptive kNN with retry (no goto/continue — Intel OpenMP safe) */
		int knn_ok = 0;
		while(!knn_ok){
			nfound = searchKNN2D(ptl+i, K, tree, neighbors, &maxr);
			if(nfound < 3 || maxr < 1e-20) break;

			Voro2D_point center;
			treevorork4particletype *ibp_tmp = (treevorork4particletype*)(ptl[i].bp);
			center.x = ibp_tmp->x;
			center.y = ibp_tmp->y;
			center.indx = PINDX(ibp_tmp);
			center.csound = ibp_tmp->csound;
			center.w2 = ibp_tmp->w2;

			ip = Voro2D_FindVC(&center, neighbors, neighwork, nfound, vorocorner, mp, boxsize);

			/* Completeness check */
			Voro2D_Corner *vtmp = vorocorner;
			postype maxvert2 = 0;
			do {
				postype dx = vtmp->x - center.x;
				postype dy = vtmp->y - center.y;
				postype d2 = dx*dx + dy*dy;
				if(d2 > maxvert2) maxvert2 = d2;
				vtmp = vtmp->upperlink;
			} while(vtmp != vorocorner);
			if(sqrt(maxvert2) > 0.95*maxr && K < MAX_NUM_NEAR){
				K = MIN(K*2, MAX_NUM_NEAR);
			} else {
				knn_ok = 1;
			}
		}

		if(knn_ok){
		treevorork4particletype *ibp = (treevorork4particletype*)(ptl[i].bp);

		get2dAreaAvgNeighorPressure(ibp, vorocorner, neighwork, bp);
		ibp->den = ibp->mass / ibp->volume;

		/* Debug: check particles with large volume */
		if(ibp->volume > 0.01){
			/* Count corners and big-box faces */
			int ncorners = 0, nbigbox = 0;
			Voro2D_Corner *dtmp = vorocorner;
			do {
				ncorners++;
				if(dtmp->upperrelated < 0) nbigbox++;
				dtmp = dtmp->upperlink;
			} while(dtmp != vorocorner);
/*			DEBUGPRINT("P%d kNN big cell: PINDX=%ld pos=(%g,%g) vol=%g nfound=%d K=%d maxr=%g ncorners=%d nbigbox=%d ip=%d\n",
				MYID(simpar), PINDX(ibp), ibp->x, ibp->y, ibp->volume,
				nfound, K, maxr, ncorners, nbigbox, ip); */
		}
		} /* if(knn_ok) */

		free(vorocorner);
		} /* if(!skip_this) */
	}

	/* Debug: check ie after kNN Voronoi loop */
	{
		int bad = 0;
		for(i=0;i<nbp;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			if(bpi->ie < 0){
				if(bad < 3)
					DEBUGPRINT("P%d kNN POST-VORO bad ie: i=%d PINDX=%ld ie=%g vol=%g pos=(%g,%g)\n",
						MYID(simpar), i, PINDX(bpi), bpi->ie, bpi->volume, bpi->x, bpi->y);
				bad++;
			}
		}
		if(bad > 0) DEBUGPRINT("P%d kNN POST-VORO: %d particles with negative ie\n", MYID(simpar), bad);
	}

	/* Global pressure/csound update (stride-aware for av_mode >= 1) */
	{
		int neg_ie = 0, neg_vol = 0, has_wall = 0;
		for(i=0;i<nbp;i++){
			treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
			if(bpi->ie < 0 && neg_ie < 3){
				DEBUGPRINT("P%d kNN pre-update: i=%d PINDX=%ld ie=%g vol=%g mass=%g pos=(%g,%g)\n",
						MYID(simpar), i, PINDX(bpi), bpi->ie, bpi->volume, bpi->mass,
						bpi->x, bpi->y);
				neg_ie++;
			}
			if(bpi->volume <= 0 && neg_vol < 3){
				DEBUGPRINT("P%d kNN bad vol: i=%d PINDX=%ld vol=%g\n",
						MYID(simpar), i, PINDX(bpi), bpi->volume);
				neg_vol++;
			}
		}
	}
	for(i=0;i<nbp;i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		bpi->pressure = bpi->ie/bpi->volume*(Gamma-1);
		if(bpi->pressure <= 0){
			static int pfloor_warn = 0;
			if(pfloor_warn < 10)
				DEBUGPRINT("P%d P-floor: i=%d PINDX=%ld ie=%g vol=%g P=%g pos=(%g,%g)\n",
					MYID(simpar), i, PINDX(bpi), bpi->ie, bpi->volume, bpi->pressure,
					bpi->x, bpi->y);
			pfloor_warn++;
			postype p_old = bpi->csound * bpi->csound * bpi->den / Gamma;
			postype p_floor = (p_old > 0) ? 0.01 * p_old : 1e-6;
			if(p_floor < 1e-6) p_floor = 1e-6;
			bpi->pressure = p_floor;
			bpi->ie = p_floor * bpi->volume / (Gamma-1);
		}
		bpi->csound = sqrt(Gamma*bpi->pressure/bpi->den);
		if(isnan(bpi->csound)){
			DEBUGPRINT("P%d kNN_den NaN at i=%d PINDX=%ld: ie=%g vol=%g den=%g P=%g mass=%g pos=(%g,%g)\n",
					MYID(simpar), i, PINDX(bpi), bpi->ie, bpi->volume, bpi->den, bpi->pressure,
					bpi->mass, bpi->x, bpi->y);
		}
	}
}

/* ================================================================
 *  updateDenW2Pressure2DBlend:
 *  Same as updateDenW2Pressure2D but for treevorostressrk4particletype.
 *  Also computes velocity gradient (∇⊗v) via Gauss divergence theorem
 *  and NS stress tensor τ per cell when av_mode >= 1.
 * ================================================================ */
void updateDenW2Pressure2DBlend(
		SimParameters *simpar,
		postype xmin, postype ymin, postype xmax, postype ymax,
		postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		postype Dtime
		){
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int av_mode = GAS_AVMODE(simpar);
	postype nu_phys_val = GAS_VISCOSITY(simpar);
	int use_stress = (av_mode >= 1 || nu_phys_val > 0);

	postype cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar);
#ifdef USE_CUDA
	det2d_dpqRK4_GPU(simpar, xmin, ymin, xmax, ymax,
		paddingAllTreeParticles, mkLinkedList2D);
#else
	det2d_dpqRK4(simpar, paddingAllTreeParticles);
#endif
	int i;
	for(i=0;i<nbp;i++){
		if(GAS_Kappa(simpar) <0) {
			bp[i].w2 = -GAS_Kappa(simpar);
		}
		else if (GAS_Kappa(simpar) >0){
			bp[i].w2hydro = getw2forHydroParticle(simpar,(treevorork4particletype*)(bp+i),Dtime);
			applyW2Controls(simpar, (treevorork4particletype*)(bp+i), Dtime);
		}
	}
	int mx, my;
	BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
	BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	CellType *cells = (VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*mx*my));
	mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax, paddingAllTreeParticles);

#ifdef USE_CUDA
	size_t _fused_p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *_fused_real_base = (char*)VORORK4_TBP(simpar);
	char *_fused_pad_base  = (char*)VORORK4_TBPP(simpar);
	int _fused_npad = VORO_NPAD(simpar);
	int _fused_nthreads = 1;
#ifdef _OPENMP
	_fused_nthreads = omp_get_max_threads();
#endif

	/* === GPU Tessellation Path === */
	{
		GPUContext *gctx = gpu_get_context();
		int n_total = nbp + _fused_npad;

		/* Lazy init */
		if (!gctx->initialized) {
			int max_p = (int)(n_total * 1.2) + 1024;
			int max_f = max_p * 20;
			gpu_init(gctx, max_p, max_f, MYID(simpar));
		}

		/* Realloc if needed */
		if (n_total > gctx->max_particles) {
			gpu_free(gctx);
			int max_p = (int)(n_total * 1.2) + 1024;
			int max_f = max_p * 20;
			gpu_init(gctx, max_p, max_f, MYID(simpar));
		}

		double _tess_t0 = MPI_Wtime();

		/* Build CellCSR from linked list */
		gpu_build_cell_csr(gctx, (void*)simpar, mx, my, _fused_p_size);

		/* Fill SoA (particles + padding) */
		{
			size_t p_size = _fused_p_size;
			char *real_base = _fused_real_base;
			char *pad_base  = _fused_pad_base;
			ParticleSoA *h = &gctx->h_parts;
			h->n_total = n_total;
			int ii;
			for (ii = 0; ii < n_total; ii++) {
				char *raw;
				if (ii < nbp)
					raw = real_base + ii * p_size;
				else
					raw = pad_base + (ii - nbp) * p_size;
				treevorork4particletype *bpi = (treevorork4particletype *)raw;
				h->x[ii]    = bpi->x;    h->y[ii]    = bpi->y;
				h->vx[ii]   = bpi->vx;   h->vy[ii]   = bpi->vy;
				h->mass[ii] = bpi->mass;
				h->indx[ii] = (long long)PINDX(bpi);
				h->den[ii]      = bpi->den;
				h->pressure[ii] = bpi->pressure;
				h->csound[ii]   = bpi->csound;
				h->volume[ii]   = bpi->volume;
				h->ie[ii]       = bpi->ie;
				h->w2[ii]       = bpi->w2;
				h->w2old[ii]    = bpi->w2old;
				if (use_stress) {
					treevorostressrk4particletype *sbpi =
						(treevorostressrk4particletype *)raw;
					h->alpha_cd[ii] = sbpi->stress.alpha_cd;
				}
			}
		}

		/* Upload particles + cells to GPU */
		gctx->d_parts.n_total = n_total;
		gpu_upload_particles(gctx, n_total, use_stress);
		gpu_upload_cells(gctx, mx * my, gctx->h_cells.n_entries);

		postype OoA = VoroAccuracyOrder(simpar);

		/* Run GPU tessellation: tessellate + gradients + P/cs/stress + face CSR */
		gpu_tessellate_and_build_faces(gctx, nbp, mx, my,
			cellsize, xmin, ymin, boxsize, OoA, av_mode, Gamma,
			nu_phys_val, GAS_PRANDTL(simpar), GAS_CDAMAX(simpar));

		/* Download results back to host */
		gpu_download_tess_results(gctx, nbp, use_stress);

		/* Write back to AoS particle structs */
		gpu_writeback_tess_results((void*)simpar, gctx, nbp, use_stress);

		double _tess_t1 = MPI_Wtime();
		if (MYID(simpar) == 0) {
			fprintf(stderr, "[GPU tess] nbp=%d npad=%d faces=%d time=%.3f ms\n",
				nbp, _fused_npad, gctx->d_faces.n_faces_total,
				(_tess_t1 - _tess_t0) * 1000.0);
		}
	}

#ifdef GPU_TESS_COMPARE
	/* === GPU FACE COMPARE: download GPU FaceCSR, let CPU run, compare === */
	{
		static int _face_cmp_call = 0;
		_face_cmp_call++;
		GPUContext *_fcctx = gpu_get_context();
		int _fc_nfaces = _fcctx->d_faces.n_faces_total;

		/* Allocate host buffers for GPU face CSR */
		int *_gf_offset = NULL;
		double *_gf_c1x=NULL, *_gf_c1y=NULL, *_gf_c2x=NULL, *_gf_c2y=NULL;
		int *_gf_neigh=NULL, *_gf_kp=NULL, *_gf_km=NULL, *_gf_ghost=NULL;

		if (_face_cmp_call <= 2 && MYID(simpar) == 0) {
			_gf_offset = (int*)malloc((nbp+1)*sizeof(int));
			_gf_c1x = (double*)malloc(_fc_nfaces*sizeof(double));
			_gf_c1y = (double*)malloc(_fc_nfaces*sizeof(double));
			_gf_c2x = (double*)malloc(_fc_nfaces*sizeof(double));
			_gf_c2y = (double*)malloc(_fc_nfaces*sizeof(double));
			_gf_neigh = (int*)malloc(_fc_nfaces*sizeof(int));
			_gf_kp = (int*)malloc(_fc_nfaces*sizeof(int));
			_gf_km = (int*)malloc(_fc_nfaces*sizeof(int));
			_gf_ghost = (int*)malloc(_fc_nfaces*sizeof(int));
			gpu_download_face_csr(_fcctx, nbp, _fc_nfaces,
				_gf_offset, _gf_c1x, _gf_c1y, _gf_c2x, _gf_c2y,
				_gf_neigh, _gf_kp, _gf_km, _gf_ghost);
			fprintf(stderr, "[FACE_CMP] call=%d gpu_nfaces=%d\n", _face_cmp_call, _fc_nfaces);
		}
		/* Reset tess_faces_on_device so CPU fused path and force use CPU faces */
		_fcctx->tess_faces_on_device = 0;
#else
	/* GPU did tessellation + P/cs/stress + face CSR. Skip CPU path. */
	goto _gpu_tess_skip_cpu;
#endif

	/* Fused face extraction init (for CPU fallback path) */
	gpu_fused_begin(nbp, _fused_npad, _fused_p_size, _fused_nthreads,
	                _fused_real_base, _fused_pad_base);
#endif

	int iy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(iy=0;iy<my;iy++){
#ifdef USE_CUDA
		int _fused_tid = 0;
#ifdef _OPENMP
		_fused_tid = omp_get_thread_num();
#endif
#endif
		int mp=1000;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		int ix;
		for(ix=0;ix<mx;ix++){
			int np;
			treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
			int nneigh;
			Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
			Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			int i;
			for(i=0;i<np;i++){
				Voro2D_point center;
				center.x = p[i].x;
				center.y = p[i].y;
				center.indx = PINDX(p+i);
				center.csound = p[i].csound;
				center.w2 = p[i].w2;

				int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
				treevorork4particletype *ibp_rk4 = p[i].bp;
				treevorostressrk4particletype *ibp = (treevorostressrk4particletype*)ibp_rk4;
				get2dAreaAvgNeighorPressure(ibp_rk4, vorocorner, neighwork, (treevorork4particletype*)bp);
				ibp_rk4->den = ibp_rk4->mass/ibp_rk4->volume;

#ifdef USE_CUDA
				/* Fallback face extraction when use_stress==false (no gradient traversal).
				   When use_stress==true, extraction is piggybacked on the gradient loop below. */
				if(!use_stress){
					int _fi = gpu_ptr_to_soa_index(ibp_rk4, _fused_real_base,
					            _fused_pad_base, nbp, _fused_npad, _fused_p_size);
					if(_fi >= 0 && _fi < nbp){
						Voro2D_Corner *_et = vorocorner;
						do {
							if(_et->upperrelated >= 0){
								Voro2D_Corner *_et2 = _et->upperlink;
								void *_jptr = neighwork[_et->upperrelated].bp;
								int _fj = gpu_ptr_to_soa_index(_jptr, _fused_real_base,
								            _fused_pad_base, nbp, _fused_npad, _fused_p_size);
								int _ghost = (PINDX((treevorork4particletype*)_jptr) == MAX_INDEX);
								int _fkp = -1, _fkm = -1;
								if(_et2->upperrelated >= 0)
									_fkp = gpu_ptr_to_soa_index(neighwork[_et2->upperrelated].bp,
									         _fused_real_base, _fused_pad_base, nbp, _fused_npad, _fused_p_size);
								if(_et->lowerrelated >= 0)
									_fkm = gpu_ptr_to_soa_index(neighwork[_et->lowerrelated].bp,
									         _fused_real_base, _fused_pad_base, nbp, _fused_npad, _fused_p_size);
								gpu_fused_add_face(_fused_tid, _fi,
									_et->x, _et->y, _et2->x, _et2->y,
									_fj, _fkp, _fkm, _ghost);
							}
							_et = _et->upperlink;
						} while(_et != vorocorner);
					}
				}
#endif

				/* Compute cell-averaged gradients via Green-Gauss:
				 *   ⟨∇Q⟩_i = (1/V_i) Σ_faces Q_face·dS
				 * Face-averaged values use M(n,m) stencil.
				 * Computes: ∇⊗v (velocity gradient), ∇P (MUSCL), ∇ρ (Arepo MUSCL).
				 *
				 * For av_mode==5, a second pass computes per-cell Arepo-style
				 * (Springel 2010, Barth-Jespersen) TVD slope limiter alpha_i
				 * for each of {rho, vx, vy, P} over ALL faces of cell i, and
				 * rescales the stored gradients in-place so that the MUSCL
				 * reconstruction at any face cannot create a new extremum
				 * relative to the cell-averaged values of i's neighbors.
				 * Other av_modes use unlimited gradients (for NS stress). */
				if(use_stress){
					postype sum_gUxx=0, sum_gUxy=0, sum_gUyx=0, sum_gUyy=0;
					postype sum_dPdx=0, sum_dPdy=0;
					postype sum_dRhodx=0, sum_dRhody=0;
					postype OoA = VoroAccuracyOrder(simpar);
					Voro2D_Corner *tmp = vorocorner;
					do {
						if(tmp->upperrelated >= 0){
							treevorostressrk4particletype *jbp = (treevorostressrk4particletype*)(neighwork[tmp->upperrelated].bp);
							int jbp_grad_ghost = (PINDX((treevorork4particletype*)jbp) == MAX_INDEX);
							Voro2D_Corner *tmp2 = tmp->upperlink;
							/* Face area vector dS (outward normal * area) */
							postype dSx = tmp2->y - tmp->y;
							postype dSy = -(tmp2->x - tmp->x);

							/* Laguerre-aware face interpolation weight */
							postype dx_ij = jbp->x - ibp->x;
							postype dy_ij = jbp->y - ibp->y;
							postype d2_ij = dx_ij*dx_ij + dy_ij*dy_ij;
							postype wfrac = (d2_ij > 0) ?
								0.5 + 0.5*(ibp_rk4->w2 - ((treevorork4particletype*)jbp)->w2)/d2_ij : 0.5;

							postype vx_face, vy_face, P_face, Rho_face;
							if(jbp_grad_ghost){
								/* Wall ghost: simple weighted average */
								vx_face = ibp->vx + wfrac*(jbp->vx - ibp->vx);
								vy_face = ibp->vy + wfrac*(jbp->vy - ibp->vy);
								P_face  = ibp->pressure + wfrac*(jbp->pressure - ibp->pressure);
								Rho_face= ibp->den + wfrac*(jbp->den - ibp->den);
							} else if(OoA > 0 && (tmp->upperlink)->upperrelated >= 0 && tmp->lowerrelated >= 0){
								/* M(n,m) 4-point stencil: face-averaged value */
								treevorostressrk4particletype *kp = (treevorostressrk4particletype*)(neighwork[(tmp->upperlink)->upperrelated].bp);
								treevorostressrk4particletype *km = (treevorostressrk4particletype*)(neighwork[tmp->lowerrelated].bp);
								postype OoA3 = OoA/3.;
								vx_face = ibp->vx + wfrac*(jbp->vx - ibp->vx)
								        + OoA3*(kp->vx + km->vx - ibp->vx - jbp->vx);
								vy_face = ibp->vy + wfrac*(jbp->vy - ibp->vy)
								        + OoA3*(kp->vy + km->vy - ibp->vy - jbp->vy);
								P_face  = ibp->pressure + wfrac*(jbp->pressure - ibp->pressure)
								        + OoA3*(kp->pressure + km->pressure - ibp->pressure - jbp->pressure);
								Rho_face= ibp->den + wfrac*(jbp->den - ibp->den)
								        + OoA3*(kp->den + km->den - ibp->den - jbp->den);
							} else {
								vx_face = ibp->vx + wfrac*(jbp->vx - ibp->vx);
								vy_face = ibp->vy + wfrac*(jbp->vy - ibp->vy);
								P_face  = ibp->pressure + wfrac*(jbp->pressure - ibp->pressure);
								Rho_face= ibp->den + wfrac*(jbp->den - ibp->den);
							}

							/* Accumulate Green-Gauss: ∮ Q_face · dS */
							sum_gUxx += vx_face * dSx;
							sum_gUxy += vx_face * dSy;
							sum_gUyx += vy_face * dSx;
							sum_gUyy += vy_face * dSy;
							sum_dPdx += P_face * dSx;
							sum_dPdy += P_face * dSy;
							sum_dRhodx += Rho_face * dSx;
							sum_dRhody += Rho_face * dSy;

#ifdef USE_CUDA
							/* Fused CSR face extraction (piggyback on gradient traversal) */
							{
								int _fi = gpu_ptr_to_soa_index(ibp_rk4, _fused_real_base,
								            _fused_pad_base, nbp, _fused_npad, _fused_p_size);
								if(_fi >= 0 && _fi < nbp){
									int _fj = gpu_ptr_to_soa_index((void*)jbp, _fused_real_base,
									            _fused_pad_base, nbp, _fused_npad, _fused_p_size);
									int _fkp = -1, _fkm = -1;
									if(tmp2->upperrelated >= 0)
										_fkp = gpu_ptr_to_soa_index(neighwork[tmp2->upperrelated].bp,
										         _fused_real_base, _fused_pad_base, nbp, _fused_npad, _fused_p_size);
									if(tmp->lowerrelated >= 0)
										_fkm = gpu_ptr_to_soa_index(neighwork[tmp->lowerrelated].bp,
										         _fused_real_base, _fused_pad_base, nbp, _fused_npad, _fused_p_size);
									gpu_fused_add_face(_fused_tid, _fi,
										tmp->x, tmp->y, tmp2->x, tmp2->y,
										_fj, _fkp, _fkm, jbp_grad_ghost);
								}
							}
#endif
						}
						tmp = tmp->upperlink;
					} while(tmp != vorocorner);

					postype invVol = 1.0/ibp_rk4->volume;
					ibp->stress.gUxx = sum_gUxx * invVol;
					ibp->stress.gUxy = sum_gUxy * invVol;
					ibp->stress.gUyx = sum_gUyx * invVol;
					ibp->stress.gUyy = sum_gUyy * invVol;
					ibp->stress.divv = ibp->stress.gUxx + ibp->stress.gUyy;
					ibp->stress.dPdx = sum_dPdx * invVol;
					ibp->stress.dPdy = sum_dPdy * invVol;
					ibp->stress.dRhodx = sum_dRhodx * invVol;
					ibp->stress.dRhody = sum_dRhody * invVol;

					/* ====== Springel 2010 TVD slope limiter (av_mode==5) ======
					 * Pass A: find min/max of {rho, vx, vy, P} over all real
					 *         neighbors of cell i (wall ghosts excluded).
					 * Pass B: for each face centroid x_f, compute trial delta =
					 *         grad Q . (x_f - x_i) and tighten alpha_Q so the
					 *         reconstructed value Q_i + alpha_Q*delta stays in
					 *         [Q_min, Q_max].  Finally rescale the stored
					 *         gradients in-place. */
					if(av_mode == 5){
						postype rho_i = ibp->den, vx_i = ibp->vx, vy_i = ibp->vy, P_i = ibp->pressure;
						postype rho_max=rho_i, rho_min=rho_i;
						postype vx_max=vx_i, vx_min=vx_i;
						postype vy_max=vy_i, vy_min=vy_i;
						postype P_max=P_i, P_min=P_i;
						Voro2D_Corner *tmpA = vorocorner;
						do {
							if(tmpA->upperrelated >= 0){
								treevorostressrk4particletype *jL = (treevorostressrk4particletype*)(neighwork[tmpA->upperrelated].bp);
								int ghostL = (PINDX((treevorork4particletype*)jL) == MAX_INDEX);
								if(!ghostL){
									if(jL->den      > rho_max) rho_max = jL->den;
									if(jL->den      < rho_min) rho_min = jL->den;
									if(jL->vx       > vx_max)  vx_max  = jL->vx;
									if(jL->vx       < vx_min)  vx_min  = jL->vx;
									if(jL->vy       > vy_max)  vy_max  = jL->vy;
									if(jL->vy       < vy_min)  vy_min  = jL->vy;
									if(jL->pressure > P_max)   P_max   = jL->pressure;
									if(jL->pressure < P_min)   P_min   = jL->pressure;
								}
							}
							tmpA = tmpA->upperlink;
						} while(tmpA != vorocorner);

						postype alpha_rho=1.0, alpha_vx=1.0, alpha_vy=1.0, alpha_P=1.0;
						Voro2D_Corner *tmpB = vorocorner;
						do {
							if(tmpB->upperrelated >= 0){
								Voro2D_Corner *tmpB2 = tmpB->upperlink;
								/* tmpB->x,tmpB->y are already RELATIVE to ibp
								   (Voro2D_FindVC uses Vec2DSub), so the face
								   midpoint IS the displacement from cell center. */
								postype dxf = 0.5*(tmpB->x + tmpB2->x);
								postype dyf = 0.5*(tmpB->y + tmpB2->y);

								postype drhoF = ibp->stress.dRhodx*dxf + ibp->stress.dRhody*dyf;
								postype dvxF  = ibp->stress.gUxx*dxf + ibp->stress.gUxy*dyf;
								postype dvyF  = ibp->stress.gUyx*dxf + ibp->stress.gUyy*dyf;
								postype dPF   = ibp->stress.dPdx*dxf + ibp->stress.dPdy*dyf;

								#define TVD_LIMIT(alpha, delta, Vi, Vmax, Vmin) do { \
									if(delta > 1e-30){ \
										postype r = (Vmax - Vi)/delta; \
										if(r < 0) r = 0; \
										if(r < alpha) alpha = r; \
									} else if(delta < -1e-30){ \
										postype r = (Vmin - Vi)/delta; \
										if(r < 0) r = 0; \
										if(r < alpha) alpha = r; \
									} \
								} while(0)

								TVD_LIMIT(alpha_rho, drhoF, rho_i, rho_max, rho_min);
								TVD_LIMIT(alpha_vx,  dvxF,  vx_i,  vx_max,  vx_min);
								TVD_LIMIT(alpha_vy,  dvyF,  vy_i,  vy_max,  vy_min);
								TVD_LIMIT(alpha_P,   dPF,   P_i,   P_max,   P_min);
								#undef TVD_LIMIT
							}
							tmpB = tmpB->upperlink;
						} while(tmpB != vorocorner);

						/* Apply limiter in-place (safe: av_mode==5 does not use
						   stress.gUxx/gUxy/... for NS stress computation). */
						ibp->stress.dRhodx *= alpha_rho;
						ibp->stress.dRhody *= alpha_rho;
						ibp->stress.gUxx   *= alpha_vx;
						ibp->stress.gUxy   *= alpha_vx;
						ibp->stress.gUyx   *= alpha_vy;
						ibp->stress.gUyy   *= alpha_vy;
						ibp->stress.dPdx   *= alpha_P;
						ibp->stress.dPdy   *= alpha_P;
					}
				}
			}
			if(np>0) free(p);
			if(nneigh>0) {
				free(neighbors); free(neighwork);
			}
		}
		free(vorocorner);
	}

	/* Update pressure, csound, and NS stress tensor */
	for(i=0;i<nbp;i++){
		bp[i].pressure = bp[i].ie/bp[i].volume*(Gamma-1);
		if(bp[i].pressure <= 0){
			bp[i].pressure = 1e-6;
			bp[i].ie = bp[i].pressure * bp[i].volume / (Gamma-1);
		}
		bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);

		/* Deviatoric NS stress convention: τ_code = -ν ρ S_dev
		   with S_dev = ∇v + ∇v^T - (2/d)·(∇·v) I.
		   We use d = 3 (3D slab / quasi-3D) — the (2/3)·divv trace
		   subtraction treats the 2D plane as a slab of a 3D fluid with
		   ∂_z = 0.  For truly 2D Navier-Stokes (d = 2) the coefficient
		   would be (2/2)·divv = divv.  Impact on nearly-incompressible
		   KH (divv ≈ 0) is negligible; for compressible RT/Sod-2D it
		   matters at the ~ν·divv level.  Do NOT change without
		   updating the documented convention in Docs/cfd_comp.tex.

		   Sign: τ_code is defined with the opposite sign of the
		   physical τ_NS so that the force loop accumulates +τ_NS·dS
		   (viscosity damps shear) via tau_dot_dS = -(τ_code·dS). */
		if(av_mode == 0 && nu_phys_val > 0){
			/* av_mode=0 + nu_phys: NS stress with fixed nu = nu_phys (no CD10) */
			postype divv = bp[i].stress.divv;
			bp[i].stress.tauxx = -nu_phys_val * bp[i].den * (2.0*bp[i].stress.gUxx - (2.0/3.0)*divv);
			bp[i].stress.tauxy = -nu_phys_val * bp[i].den * (bp[i].stress.gUxy + bp[i].stress.gUyx);
			bp[i].stress.tauyy = -nu_phys_val * bp[i].den * (2.0*bp[i].stress.gUyy - (2.0/3.0)*divv);
		} else if(av_mode == 1){
			/* NS stress: τ = -ν ρ (∇v + ∇v^T - 2/3 (∇·v) I)
			   ν = max(nu_phys, α_CD·h·c_s) */
			postype h = sqrt(bp[i].volume);
			postype nu_cd = bp[i].stress.alpha_cd * h * bp[i].csound;
			postype nu_phys = GAS_VISCOSITY(simpar);
			postype nu = (nu_phys > 0) ? fmax(nu_phys, nu_cd) : nu_cd;
			postype divv = bp[i].stress.divv;
			bp[i].stress.tauxx = -nu * bp[i].den * (2.0*bp[i].stress.gUxx - (2.0/3.0)*divv);
			bp[i].stress.tauxy = -nu * bp[i].den * (bp[i].stress.gUxy + bp[i].stress.gUyx);
			bp[i].stress.tauyy = -nu * bp[i].den * (2.0*bp[i].stress.gUyy - (2.0/3.0)*divv);
		} else if(av_mode >= 2){
			/* Design intent — NOT a bug.  av_mode>=2 uses CD10 (α_CD)
			   as the dissipation mechanism in three places:
			     (a) blend factor f_pq = α_max/cd_amax (getAccVoro2DBlend)
			     (b) HLLC viscous pressure Π_CD10 = ½α vsig ρ̄ (-Δv_n)
			     (c) weight (1-f_pq) for M(n,m) path
			   NS stress would double-count α_CD dissipation, so τ is
			   intentionally zeroed here as defensive cleanup.  The
			   force loop gate at getAccVoro2DBlend already excludes
			   av_mode>=2 from applying tau_dot_dS — this zeroing is
			   redundant with that gate but guards against stale reads
			   from diagnostics. */
			bp[i].stress.tauxx = 0;
			bp[i].stress.tauxy = 0;
			bp[i].stress.tauyy = 0;
		}
	}

#ifdef GPU_TESS_COMPARE
	/* === Compare GPU face CSR vs CPU fused face CSR === */
	if (_face_cmp_call <= 2 && MYID(simpar) == 0 && _gf_offset) {
		/* CPU fused faces are now in gpu_get_context()->h_faces */
		FaceCSR *cf = &_fcctx->h_faces;
		int cpu_nfaces = cf->n_faces_total;
		fprintf(stderr, "[FACE_CMP] cpu_nfaces=%d gpu_nfaces=%d\n", cpu_nfaces, _fc_nfaces);

		/* Compare per-particle face counts */
		int count_mismatch = 0, first_mm = -1;
		for(int ii=0; ii<nbp; ii++){
			int gpu_fc = _gf_offset[ii+1] - _gf_offset[ii];
			int cpu_fc = cf->face_offset[ii+1] - cf->face_offset[ii];
			if(gpu_fc != cpu_fc){
				if(first_mm < 0) first_mm = ii;
				count_mismatch++;
			}
		}
		fprintf(stderr, "[FACE_CMP] face_count_mismatch=%d (of %d particles) first=%d\n",
			count_mismatch, nbp, first_mm);

		/* For each particle, compare faces by matching neighbor_idx */
		int n_neigh_err=0, n_kp_err=0, n_km_err=0, n_ghost_err=0;
		int n_c1_err=0, n_c2_err=0;
		int first_neigh_err=-1;
		for(int ii=0; ii<nbp && ii<nbp; ii++){
			int gbase = _gf_offset[ii], gend = _gf_offset[ii+1];
			int cbase = cf->face_offset[ii], cend = cf->face_offset[ii+1];
			if(gend - gbase != cend - cbase) continue; /* skip count mismatches */
			int nf = gend - gbase;
			/* Match GPU faces to CPU faces by neighbor_idx */
			for(int gf=0; gf<nf; gf++){
				int gn = _gf_neigh[gbase+gf];
				/* Find matching CPU face */
				int found = 0;
				for(int cf2=0; cf2<nf; cf2++){
					if(cf->neighbor_idx[cbase+cf2] == gn){
						found = 1;
						int gi = gbase+gf, ci = cbase+cf2;
						/* Compare kp/km */
						if(_gf_kp[gi] != cf->kp_idx[ci]) n_kp_err++;
						if(_gf_km[gi] != cf->km_idx[ci]) n_km_err++;
						if(_gf_ghost[gi] != cf->is_ghost[ci]) n_ghost_err++;
						/* Compare corners (relative, should be close) */
						double dc1 = fabs(_gf_c1x[gi]-cf->c1x[ci])+fabs(_gf_c1y[gi]-cf->c1y[ci]);
						double dc2 = fabs(_gf_c2x[gi]-cf->c2x[ci])+fabs(_gf_c2y[gi]-cf->c2y[ci]);
						if(dc1 > 1e-10) n_c1_err++;
						if(dc2 > 1e-10) n_c2_err++;
						/* Print first few errors */
						if(n_kp_err+n_km_err<=5 && (_gf_kp[gi]!=cf->kp_idx[ci] || _gf_km[gi]!=cf->km_idx[ci])){
							fprintf(stderr, "[FACE_CMP] i=%d neigh=%d: gpu_kp=%d cpu_kp=%d gpu_km=%d cpu_km=%d\n",
								ii, gn, _gf_kp[gi], cf->kp_idx[ci], _gf_km[gi], cf->km_idx[ci]);
						}
						if(n_c1_err+n_c2_err<=3 && (dc1>1e-10 || dc2>1e-10)){
							fprintf(stderr, "[FACE_CMP] i=%d neigh=%d: corner err dc1=%g dc2=%g\n"
								"  gpu c1=(%g,%g) c2=(%g,%g)\n  cpu c1=(%g,%g) c2=(%g,%g)\n",
								ii, gn, dc1, dc2,
								_gf_c1x[gi],_gf_c1y[gi],_gf_c2x[gi],_gf_c2y[gi],
								cf->c1x[ci],cf->c1y[ci],cf->c2x[ci],cf->c2y[ci]);
						}
						break;
					}
				}
				if(!found){
					n_neigh_err++;
					if(first_neigh_err<0) first_neigh_err = ii;
				}
			}
		}
		fprintf(stderr, "[FACE_CMP] neigh_unmatched=%d kp_err=%d km_err=%d ghost_err=%d c1_err=%d c2_err=%d\n",
			n_neigh_err, n_kp_err, n_km_err, n_ghost_err, n_c1_err, n_c2_err);
		if(first_neigh_err>=0)
			fprintf(stderr, "[FACE_CMP] first_neigh_unmatched at particle %d\n", first_neigh_err);

		/* Print first particle's faces for eyeball check */
		if(nbp > 0){
			int gbase = _gf_offset[0], gend = _gf_offset[1];
			int cbase = cf->face_offset[0], cend = cf->face_offset[1];
			fprintf(stderr, "[FACE_CMP] particle 0: gpu_nf=%d cpu_nf=%d\n", gend-gbase, cend-cbase);
			for(int f=gbase;f<gend && f<gbase+4;f++)
				fprintf(stderr, "[FACE_CMP]   gpu f=%d: neigh=%d kp=%d km=%d ghost=%d c1=(%g,%g) c2=(%g,%g)\n",
					f-gbase, _gf_neigh[f], _gf_kp[f], _gf_km[f], _gf_ghost[f],
					_gf_c1x[f], _gf_c1y[f], _gf_c2x[f], _gf_c2y[f]);
			for(int f=cbase;f<cend && f<cbase+4;f++)
				fprintf(stderr, "[FACE_CMP]   cpu f=%d: neigh=%d kp=%d km=%d ghost=%d c1=(%g,%g) c2=(%g,%g)\n",
					f-cbase, cf->neighbor_idx[f], cf->kp_idx[f], cf->km_idx[f], cf->is_ghost[f],
					cf->c1x[f], cf->c1y[f], cf->c2x[f], cf->c2y[f]);
		}

		free(_gf_offset); free(_gf_c1x); free(_gf_c1y); free(_gf_c2x); free(_gf_c2y);
		free(_gf_neigh); free(_gf_kp); free(_gf_km); free(_gf_ghost);
	}
	} /* close GPU_TESS_COMPARE block */
#endif

#ifdef USE_CUDA
_gpu_tess_skip_cpu:
#endif
	my_free(VORO_BASICCELL(simpar));
#ifdef USE_CUDA
	if (!gpu_get_context()->tess_faces_on_device) {
		gpu_fused_end();
	}
	/* Keep VORORK4_TBPP alive — GPU wrapper will free it after SoA fill */
#else
	my_free(VORORK4_TBPP(simpar));
#endif
}

/* ================================================================
 *  LagMFM helpers (av_mode=4): Wendland C2 kernel (2D) + 2x2 inverse
 *  Used by meshless path: kernel density, gradients, effective faces.
 *  Kernel support: r < 2h (q = r/h in [0,2])
 * ================================================================ */
static inline postype mfm_W_wendland2d(postype r, postype h)
{
	/* 2D Wendland C2 kernel, compact support 2h, normalization 7/(4 pi h^2).
	 * W(q) = (7/(4 pi h^2)) * (1 - q/2)^4 * (1 + 2 q),  q = r/h,  q in [0,2]
	 */
	if(h <= 0) return 0;
	postype q = r / h;
	if(q >= 2.0) return 0;
	postype one_minus_half_q = 1.0 - 0.5*q;
	postype t2 = one_minus_half_q*one_minus_half_q;
	postype t4 = t2*t2;
	return (7.0/(4.0*M_PI*h*h)) * t4 * (1.0 + 2.0*q);
}

static inline postype mfm_dWdr_wendland2d(postype r, postype h)
{
	/* dW/dr for 2D Wendland C2.
	 * W(r) = C * (1 - r/(2h))^4 * (1 + 2 r/h), C = 7/(4 pi h^2)
	 * Let u = r/(2h), v = 1 - u. W = C v^4 (1 + 4 u)  [since 2 r/h = 4 u]
	 * dW/du = C [ -4 v^3 (1 + 4 u) + v^4 * 4 ]
	 *       = 4 C v^3 [ v - (1 + 4 u) ]
	 *       = 4 C v^3 [ -5 u ]
	 *       = -20 C u v^3
	 * dW/dr = dW/du * du/dr = -20 C u v^3 / (2 h) = -10 C u v^3 / h
	 */
	if(h <= 0) return 0;
	postype q = r / h;
	if(q >= 2.0) return 0;
	postype u = 0.5 * q;          /* r/(2h) */
	postype v = 1.0 - u;
	postype v3 = v*v*v;
	postype C = 7.0/(4.0*M_PI*h*h);
	return -10.0 * C * u * v3 / h;
}

static inline int mfm_invert2x2(postype Exx, postype Exy, postype Eyx, postype Eyy,
		postype *Ixx, postype *Ixy, postype *Iyx, postype *Iyy)
{
	/* D-B5 fix: Invert a 2x2 matrix [[Exx,Exy],[Eyx,Eyy]] -> store into Ixx..Iyy.
	 * Returns 1 on success, 0 if near-singular.
	 * On failure, writes a ZERO matrix (not identity).  Callers then get
	 *   ∇Q = (E_inv) · S = 0,
	 * which degrades the particle to first-order (central state) — a safe
	 * and self-degrading fallback.  Identity would have given badly scaled
	 * gradients (~ O(1) instead of O(1/h²)) and corrupted the force loop.
	 * We also use a relative tolerance ‖E‖² · eps² to avoid triggering on
	 * legitimately small but non-singular E matrices in low-density regions.
	 */
	postype det = Exx*Eyy - Exy*Eyx;
	postype adet = det >= 0 ? det : -det;
	postype nrm2 = Exx*Exx + Exy*Exy + Eyx*Eyx + Eyy*Eyy;
	postype tol  = 1.0e-14 * nrm2 + 1.0e-30;
	if(adet < tol){
		*Ixx = 0.0; *Ixy = 0.0;
		*Iyx = 0.0; *Iyy = 0.0;
		return 0;
	}
	postype invdet = 1.0/det;
	*Ixx =  Eyy*invdet;
	*Ixy = -Exy*invdet;
	*Iyx = -Eyx*invdet;
	*Iyy =  Exx*invdet;
	return 1;
}

/* LagMFM tuning knobs.
 *
 *  LAGMFM_ETA       : h_i = eta * sqrt(V_i).  With eta=1.8 we get ~30
 *                     neighbors per particle for uniform layouts.
 *  LAGMFM_H_ITER    : max Newton iterations to converge h_i so that
 *                     V_i = 1/n_i(h_i) is self-consistent (D-B1 fix).
 *  LAGMFM_H_TOL     : relative convergence tolerance on h_i.
 *  LAGMFM_WALL_NOSLIP: 0 = free-slip wall (tangential v_face = v_i_tan),
 *                     1 = no-slip wall (v_face = 0).   D-B6 option.
 */
#define LAGMFM_ETA           ((postype)1.8)
#define LAGMFM_H_ITER        3
#define LAGMFM_H_TOL         ((postype)0.005)
#define LAGMFM_WALL_NOSLIP   0

/* ================================================================
 *  updateDenW2Pressure2D_LagMFM (av_mode=4):
 *    Kernel-based density and matrix-weighted gradients (Hopkins 2015).
 *    Replaces Voronoi tessellation with a direct neighbor walk on the
 *    existing cell linked list.  Produces:
 *      - ibp->den, ibp->volume, ibp->csound, ibp->pressure
 *      - stress.gUxx..gUyy, divv, dPdx, dPdy
 *      - stress.E_inv_xx..E_inv_yy, stress.h_mfm
 *    NS stress tau is filled the same way as av_mode=1.
 * ================================================================ */
#ifdef USE_CUDA
/* w2 refresh helper — called from integrator AFTER GPU density kernel,
   which already computed w2ceil via nearest-neighbor in the neighbor loop.
   This handles w2 = f(P) from density results + applyW2Controls. */
static void lagmfm_w2_post_density(SimParameters *simpar, postype Dtime)
{
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int i;
	for(i=0;i<nbp;i++){
		if(GAS_Kappa(simpar) < 0){
			bp[i].w2 = -GAS_Kappa(simpar);
		} else if(GAS_Kappa(simpar) > 0){
			bp[i].w2hydro = getw2forHydroParticle(simpar,(treevorork4particletype*)(bp+i),Dtime);
			applyW2Controls(simpar, (treevorork4particletype*)(bp+i), Dtime);
		} else {
			/* Kappa=0: w2ceil from GPU density kernel, clip w2 */
			bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
		}
	}
}
#endif
void updateDenW2Pressure2D_LagMFM(
		SimParameters *simpar,
		postype xmin, postype ymin, postype xmax, postype ymax,
		postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		postype Dtime
		){
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype nu_phys_val = GAS_VISCOSITY(simpar);

	/* Refresh w2 for Laguerre weight (reused for dpq + w2ceil tracking even
	 * though LagMFM uses h_mfm, not w2, as kernel bandwidth). */
	det2d_dpqRK4(simpar, paddingAllTreeParticles);
	int i;
	for(i=0;i<nbp;i++){
		if(GAS_Kappa(simpar) < 0){
			bp[i].w2 = -GAS_Kappa(simpar);
		} else if(GAS_Kappa(simpar) > 0){
			bp[i].w2hydro = getw2forHydroParticle(simpar,(treevorork4particletype*)(bp+i),Dtime);
			applyW2Controls(simpar, (treevorork4particletype*)(bp+i), Dtime);
		}
	}

	/* Build cell linked list */
	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	int mx, my;
	BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
	BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	CellType *cells = (VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*mx*my));
	mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax, paddingAllTreeParticles);

	postype invcs = 1.0/cellsize;
	postype eta   = LAGMFM_ETA;

	/* -----------------------------------------------------------------
	 * Pass 1: kernel density + matrix-weighted gradients
	 * -----------------------------------------------------------------
	 *   n_i        = Σ_j  W(r_ij, h_i)
	 *   V_i        = 1 / n_i
	 *   ρ_i        = m_i · n_i
	 *   E_i^{αβ}   = Σ_j W_ij (x_j-x_i)^α (x_j-x_i)^β
	 *   ψ_j(x_i)   = W_ij · Σ_β E_inv^{αβ} (x_j-x_i)^β   (used in force)
	 *   ∇Q|_i^α    = Σ_β E_inv^{αβ} · Σ_j W_ij (Q_j-Q_i) (x_j-x_i)^β
	 *
	 * h_i is set from the previous volume:   h_i = eta * sqrt(V_i_prev).
	 * On first call V_i comes from kh_findVol (Voronoi) so the
	 * bootstrap is well behaved.
	 * ----------------------------------------------------------------- */
	int iy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(iy=0; iy<my; iy++){
		int ix;
		for(ix=0; ix<mx; ix++){
			struct linkedlisttype *li = cells[ix + mx*iy].link;
			while(li){
				treevorostressrk4particletype *ibp =
						(treevorostressrk4particletype*)li;
				int iwall = (PINDX((treevorork4particletype*)ibp) == MAX_INDEX);
				if(iwall || IS_FLAG_ON(li,BoundaryGhostflag)){
					li = li->next;
					continue;
				}

				/* D-B1 fix: adaptive h_i via Newton iteration.
				 * Start with h_i = eta * sqrt(V_prev), then iterate:
				 *   1) compute n_i = Σ_j W(r_ij, h_i) and moments E, S
				 *   2) V_new = 1/n_i, h_new = eta * sqrt(V_new)
				 *   3) stop if |h_new - h_i| < tol * h_i (converged)
				 * Only the LAST iteration's E, S values are used downstream,
				 * so on non-convergence we accept the final iteration's state.
				 * nbr (cell walk radius) is recomputed each iteration since
				 * a grown h_i needs a wider stencil. */
				postype Vprev = ibp->volume > 0 ? ibp->volume : cellsize*cellsize;
				postype h_i = eta * sqrt(Vprev);

				postype n_i = 0;
				postype Exx=0, Exy=0, Eyy=0;
				postype Sxx_vx=0, Sxy_vx=0;
				postype Sxx_vy=0, Sxy_vy=0;
				postype Sxx_P =0, Sxy_P =0;

				int h_iter_done = 0;
				int h_it;
				for(h_it = 0; h_it < LAGMFM_H_ITER; h_it++){
					/* Reset accumulators each iteration. */
					n_i = 0;
					Exx = Exy = Eyy = 0;
					Sxx_vx = Sxy_vx = 0;
					Sxx_vy = Sxy_vy = 0;
					Sxx_P  = Sxy_P  = 0;

					postype two_h = 2.0 * h_i;
					postype two_h2 = two_h * two_h;
					int nbr = (int)ceil(two_h * invcs);
					if(nbr < 1) nbr = 1;
					if(nbr > 10) nbr = 10;

					int jy;
					for(jy = iy-nbr; jy <= iy+nbr; jy++){
						if(jy < 0 || jy >= my) continue;
						int jx;
						for(jx = ix-nbr; jx <= ix+nbr; jx++){
							if(jx < 0 || jx >= mx) continue;
							struct linkedlisttype *lj = cells[jx + mx*jy].link;
							while(lj){
								treevorostressrk4particletype *jbp =
										(treevorostressrk4particletype*)lj;
								postype dx = jbp->x - ibp->x;
								postype dy = jbp->y - ibp->y;
								postype r2 = dx*dx + dy*dy;
								if(r2 < two_h2){
									postype r = sqrt(r2);
									postype Wij = mfm_W_wendland2d(r, h_i);
									n_i += Wij;
									/* E matrix is symmetric in 2D */
									Exx += Wij * dx * dx;
									Exy += Wij * dx * dy;
									Eyy += Wij * dy * dy;
									/* First moments of (Q_j - Q_i) */
									postype dvx = jbp->vx - ibp->vx;
									postype dvy = jbp->vy - ibp->vy;
									postype dP  = jbp->pressure - ibp->pressure;
									Sxx_vx += Wij * dvx * dx;
									Sxy_vx += Wij * dvx * dy;
									Sxx_vy += Wij * dvy * dx;
									Sxy_vy += Wij * dvy * dy;
									Sxx_P  += Wij * dP  * dx;
									Sxy_P  += Wij * dP  * dy;
								}
								lj = lj->next;
							}
						}
					}

					if(n_i <= 0) break;  /* isolated */

					postype V_new = 1.0 / n_i;
					postype h_new = eta * sqrt(V_new);
					postype dh    = h_new - h_i;
					postype adh   = dh >= 0 ? dh : -dh;
					if(adh < LAGMFM_H_TOL * h_i){
						h_i = h_new;
						h_iter_done = 1;
						break;
					}
					h_i = h_new;
				}
				(void)h_iter_done;  /* reserved for optional diagnostics */

				/* Finalize density / volume */
				if(n_i <= 0){
					/* Isolated particle — retain previous state */
					li = li->next;
					continue;
				}
				postype Vi = 1.0 / n_i;
				ibp->volume = Vi;
				ibp->den    = ibp->mass * n_i;

				/* Invert E to get the matrix-weighted gradient operator */
				postype Einv_xx, Einv_xy, Einv_yx, Einv_yy;
				mfm_invert2x2(Exx, Exy, Exy, Eyy,
						&Einv_xx, &Einv_xy, &Einv_yx, &Einv_yy);

				ibp->stress.E_inv_xx = Einv_xx;
				ibp->stress.E_inv_xy = Einv_xy;
				ibp->stress.E_inv_yx = Einv_yx;
				ibp->stress.E_inv_yy = Einv_yy;
				ibp->stress.h_mfm    = h_i;

				/* Gradients:  ∇Q|_i^α = Σ_β (E_inv)^{αβ} Σ_j W_ij dQ_ij dx_ij^β
				 * (no 1/n_i factor — it cancels between E and S because we
				 *  dropped the common 1/n_i prefactor from both.) */
				postype gUxx = Einv_xx*Sxx_vx + Einv_xy*Sxy_vx;
				postype gUxy = Einv_yx*Sxx_vx + Einv_yy*Sxy_vx;
				postype gUyx = Einv_xx*Sxx_vy + Einv_xy*Sxy_vy;
				postype gUyy = Einv_yx*Sxx_vy + Einv_yy*Sxy_vy;
				postype dPdx = Einv_xx*Sxx_P  + Einv_xy*Sxy_P;
				postype dPdy = Einv_yx*Sxx_P  + Einv_yy*Sxy_P;

				ibp->stress.gUxx = gUxx;
				ibp->stress.gUxy = gUxy;
				ibp->stress.gUyx = gUyx;
				ibp->stress.gUyy = gUyy;
				ibp->stress.divv = gUxx + gUyy;
				ibp->stress.dPdx = dPdx;
				ibp->stress.dPdy = dPdy;

				li = li->next;
			}
		}
	}

	/* -----------------------------------------------------------------
	 * Pass 2: pressure, sound speed, NS stress tensor (τ = -ν ρ S_dev)
	 *
	 * A-M1 analog: deviatoric NS stress with (2/d)·divv trace subtraction.
	 *   We use d = 3 (3D slab / quasi-3D convention).  For truly 2D NS
	 *   the coefficient would be (2/2)·divv.  Impact on KH (divv ≈ 0) is
	 *   negligible, but matters at the ~ν·divv level for compressible
	 *   RT/Sod-2D.  Convention matches av_mode=1 (see updateDenW2Pressure2DBlend).
	 * Sign: τ_code = -ν ρ S_dev (negative of physical τ_NS).  The force
	 *   loop computes τ_NS·A via tau_A = -(τ_code·A).
	 *
	 * A-M2 analog: stride-aware access via TVORORK4_DDINFO.  bp is typed
	 *   treevorostressrk4particletype (av_mode >= 1 allocates this), so a
	 *   plain bp[i] works, but we honor n_size for future-proofing.
	 * ----------------------------------------------------------------- */
	{
		size_t p_size2 = TVORORK4_DDINFO(simpar)[0].n_size;
		char *bp_raw2  = (char*)VORORK4_TBP(simpar);
		for(i=0;i<nbp;i++){
			treevorostressrk4particletype *bpi =
				(treevorostressrk4particletype*)(bp_raw2 + i*p_size2);
			bpi->pressure = bpi->ie/bpi->volume*(Gamma-1);
			if(bpi->pressure <= 0){
				bpi->pressure = 1e-6;
				bpi->ie = bpi->pressure * bpi->volume / (Gamma-1);
			}
			bpi->csound = sqrt(Gamma*bpi->pressure/bpi->den);

			/* NS stress (3D-slab trace convention; matches av_mode=1). */
			postype h = sqrt(bpi->volume);
			postype nu_cd = bpi->stress.alpha_cd * h * bpi->csound;
			postype nu = (nu_phys_val > 0) ? fmax(nu_phys_val, nu_cd) : nu_cd;
			postype divv = bpi->stress.divv;
			bpi->stress.tauxx = -nu * bpi->den * (2.0*bpi->stress.gUxx - (2.0/3.0)*divv);
			bpi->stress.tauxy = -nu * bpi->den * (bpi->stress.gUxy + bpi->stress.gUyx);
			bpi->stress.tauyy = -nu * bpi->den * (2.0*bpi->stress.gUyy - (2.0/3.0)*divv);
		}
	}

	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
}

/* ================================================================
 *  hllc_face_2d: HLLC approximate Riemann solver
 *  Input: left/right states in face-normal direction
 *  Output: star-region pressure and normal velocity
 * ================================================================ */
static inline void hllc_face_2d(
		postype rhoL, postype pL, postype vnL, postype cL,
		postype rhoR, postype pR, postype vnR, postype cR,
		postype Gamma,
		postype *pstar, postype *vnstar)
{
	postype ZL = rhoL * cL, ZR = rhoR * cR;
	postype GP1 = Gamma + 1.0;

	/* Acoustic estimate (PVRS) */
	postype p_pvrs = (ZR*pL + ZL*pR + ZL*ZR*(vnL - vnR)) / (ZL + ZR);
	if(p_pvrs < 0) p_pvrs = 0;

	/* Shock corrections */
	postype qL = 1.0, qR = 1.0;
	if(p_pvrs > pL)
		qL = sqrt(1.0 + GP1/(2.0*Gamma) * (p_pvrs/pL - 1.0));
	if(p_pvrs > pR)
		qR = sqrt(1.0 + GP1/(2.0*Gamma) * (p_pvrs/pR - 1.0));

	postype WL = rhoL * cL * qL;
	postype WR = rhoR * cR * qR;
	postype Ws = WL + WR;

	*pstar  = (WR*pL + WL*pR + WL*WR*(vnL - vnR)) / Ws;
	*vnstar = (WL*vnL + WR*vnR + pL - pR) / Ws;
}

/* ================================================================
 *  hllc_face_2d_rest_frame: HLLC in the face rest frame (Arepo-style)
 *
 *  Boost the normal velocities into the face frame w_n (face moves
 *  with velocity w in the lab frame), solve HLLC there, then unboost.
 *  Pressure is Galilean-invariant, so p_star is returned unchanged.
 *  The normal star velocity is returned in the lab frame.
 *
 *  Key property (Springel 2010 §32): for Lagrangian particles where
 *  w_ij = (v_i + v_j)/2 is close to both v_i and v_j, the face-frame
 *  relative velocities are small, so HLLC dissipation is small and
 *  contact discontinuities do not suffer from the "fake compression"
 *  artifact seen with lab-frame Godunov on moving meshes.
 * ================================================================ */
static inline void hllc_face_2d_rest_frame(
		postype rhoL, postype pL, postype vnL_lab, postype cL,
		postype rhoR, postype pR, postype vnR_lab, postype cR,
		postype wn,
		postype Gamma,
		postype *pstar, postype *vnstar_lab)
{
	postype vnL = vnL_lab - wn;
	postype vnR = vnR_lab - wn;
	postype pst, vnst;
	hllc_face_2d(rhoL, pL, vnL, cL,
	             rhoR, pR, vnR, cR,
	             Gamma, &pst, &vnst);
	*pstar       = pst;        /* pressure is frame-invariant */
	*vnstar_lab  = vnst + wn;  /* boost back to lab frame */
}

/* ================================================================
 *  update_alpha_cd_2d: Cullen-Dehnen (2010) viscosity switch
 *  Called once per full RK4 timestep after final state combination.
 * ================================================================ */
static void update_alpha_cd_2d(SimParameters *simpar, postype dt){
	/* Stride-aware access: honor TVORORK4_DDINFO n_size in case future
	   extensions append fields beyond Stress. */
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw = (char*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype cd_amax = GAS_CDAMAX(simpar);
	postype cd_ell  = GAS_CDELL(simpar);
	postype cd_amin = GAS_CDAMIN(simpar);

	/* Revision 4 defensive guard: skip if CD10 is disabled or dt is
	 * invalid (happens on first integrator call before dt estimation
	 * or if parameters are corrupt). Prevents NaN from divisions. */
	if(cd_amax <= 0 || cd_ell <= 0 || dt <= 0) return;

	int i;
	for(i=0;i<nbp;i++){
		treevorostressrk4particletype *bpi =
			(treevorostressrk4particletype*)(bp_raw + i*p_size);
		postype h = sqrt(bpi->volume);
		postype A = fmax(0.0, -(bpi->stress.divv - bpi->stress.divv_old) / dt);
		postype vs = bpi->stress.vsig_max;
		postype alpha_loc = cd_amax * h*h * A / (vs*vs + h*h * A + 1e-30);
		postype tau = h / (2.0 * cd_ell * vs + 1e-30);
		if(alpha_loc > bpi->stress.alpha_cd)
			bpi->stress.alpha_cd = alpha_loc;
		else
			bpi->stress.alpha_cd = alpha_loc + (bpi->stress.alpha_cd - alpha_loc) * exp(-dt/tau);
		if(bpi->stress.alpha_cd < cd_amin) bpi->stress.alpha_cd = cd_amin;
		if(bpi->stress.alpha_cd > cd_amax) bpi->stress.alpha_cd = cd_amax;
		bpi->stress.divv_old = bpi->stress.divv;
		bpi->stress.vsig_max = 0;
	}
}

/* ================================================================
 *  getAccVoro2DBlend: Two-tier blended force computation.
 *  Per-face blend between M(n,m)+NS-stress and HLLC+CD10.
 * ================================================================ */
double getAccVoro2DBlend(SimParameters *simpar, postype xmin, postype ymin,
		postype xmax, postype ymax,
		postype OrderOfAccuracy, postype Courant, postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype))
		){
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype Dtime = 1.e10;
	int nx = NX(simpar);
	int av_mode = GAS_AVMODE(simpar);
	int use_muscl = GAS_USEMUSCL(simpar);
	postype cd_amax = GAS_CDAMAX(simpar);
	postype blend_theta = GAS_BLENDTHETA(simpar);
	postype nu_phys = GAS_VISCOSITY(simpar);
	postype prandtl = GAS_PRANDTL(simpar);

	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	postype cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar);
	int mx, my;
	BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
	BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
	mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);

	postype dtold = GAS_dtold(simpar);

	/* Monaghan AV parameters (used when av_mode==0) */
	float alphavis = GAS_AlphaVis(simpar);
	float betavis = GAS_BetaVis(simpar);
	float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
	float epsvis = GAS_EPSVIS(simpar);

	int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
	for(iy=1;iy<my-1;iy++){
		int mp=1000;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		int ix;
		for(ix=1;ix<mx-1;ix++){
			int np;
			treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
			int nneigh;
			Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
			Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			int i;
			for(i=0;i<np;i++){
				Voro2D_point center;
				center.x = p[i].x;
				center.y = p[i].y;
				center.indx = PINDX(p+i);
				center.csound = p[i].csound;
				center.w2 = p[i].w2;

				treevorork4particletype *ibp_rk4 = p[i].bp;
				treevorostressrk4particletype *ibp = (treevorostressrk4particletype*)ibp_rk4;
				ibp_rk4->dt = 1.e10;

				int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
				postype ibp_vx = ibp->vx;
				postype ibp_vy = ibp->vy;
				postype ibp_csound = ibp->csound;
				postype ibp_pressure = ibp->pressure;
				postype ibp_den = ibp->den;

				Voro2D_Corner *tmp,*tmp2;
				tmp = vorocorner;
				double die, dte, dke, fx, fy;
				die = dte = dke = fx = fy = 0;
				do {
					if(tmp->upperrelated >= 0){
						treevorostressrk4particletype *jbp = (treevorostressrk4particletype*)(neighwork[tmp->upperrelated].bp);
						int jbp_is_ghost = (PINDX((treevorork4particletype*)jbp) == MAX_INDEX);
						postype jbp_vx = jbp->vx;
						postype jbp_vy = jbp->vy;
						postype jbp_csound = jbp->csound;
						postype jbp_pressure = jbp->pressure;
						postype jbp_den = jbp->den;

						/* Face geometry */
						tmp2 = tmp->upperlink;
						Voro2D_point line;
						line.x = tmp2->x - tmp->x;
						line.y = tmp2->y - tmp->y;
						Voro2D_point dS = voro2D_norm(&line);
						postype facearea = Vec2DLength(tmp,tmp2);
						dS.x = facearea*dS.x;
						dS.y = facearea*dS.y;

						/* Inter-particle direction */
						Voro2D_point dr = EunhaVec2DSub((treevorork4particletype*)jbp, ibp_rk4);
						postype dramp = sqrt(Vec2DDotP(&dr, &dr));
						Voro2D_point er;
						er.x = dr.x/dramp;
						er.y = dr.y/dramp;

						postype pi_total;
						postype tau_dot_dS_x = 0, tau_dot_dS_y = 0;

						if(av_mode == 0 && nu_phys <= 0){
							/* Original Monaghan AV path (no NS stress) */
							postype pi = VoroRK4_Pressure2D(
									ibp_pressure, jbp_pressure, tmp, neighwork, simpar, RT);

							Voro2D_point uij;
							uij.x = jbp_vx - ibp_vx;
							uij.y = jbp_vy - ibp_vy;
							postype rvel = Vec2DDotP(&er, &uij);
							if(rvel < 0){
								postype wcomp = sqrt(ibp_rk4->w2)+sqrt(((treevorork4particletype*)jbp)->w2);
								postype scaleFactor = (wcomp==0 ? etavis: wcomp);
								postype drampScale = dramp/scaleFactor;
								postype mu = rvel/(drampScale + epsvis/drampScale);
								postype meanden = 0.5*(ibp_den + jbp_den);
								postype meanCsound = 0.5*(ibp_csound + jbp_csound);
								pi = pi + (-alphavis*meanCsound*mu + betavis*mu*mu)*meanden;
							}
							pi_total = pi;

						} else if(jbp_is_ghost) {
							/* Ghost face in blend path: use M(n,m) pressure
							   + Monaghan AV for wall damping (same as av_mode=0).
							   Skip NS stress and HLLC — ghost stress fields are
							   not properly mirrored and HLLC sees artificial collision. */
							postype pi = VoroStressRK4_Pressure2D(
									ibp_pressure, jbp_pressure, tmp, neighwork, simpar, RT);
							Voro2D_point uij;
							uij.x = jbp_vx - ibp_vx;
							uij.y = jbp_vy - ibp_vy;
							postype rvel = Vec2DDotP(&er, &uij);
							if(rvel < 0){
								postype wcomp = sqrt(ibp_rk4->w2)+sqrt(((treevorork4particletype*)jbp)->w2);
								postype scaleFactor = (wcomp==0 ? etavis: wcomp);
								postype drampScale = dramp/scaleFactor;
								postype mu = rvel/(drampScale + epsvis/drampScale);
								postype meanden = 0.5*(ibp_den + jbp_den);
								postype meanCsound = 0.5*(ibp_csound + jbp_csound);
								pi = pi + (-alphavis*meanCsound*mu + betavis*mu*mu)*meanden;
							}
							pi_total = pi;

						} else if(av_mode == 3) {
							/* ========== DEPRECATED: Pure HLLC (lab frame) ==========
							 * WARNING (Revision 4, expert review):
							 *   This mode is mathematically guaranteed to fail on
							 *   Lagrangian Voronoi/Laguerre meshes in flows with
							 *   contact discontinuities (e.g. KH, RT). The contact
							 *   wave is a linearly degenerate characteristic, so a
							 *   Riemann solver alone carries zero entropy flux. In
							 *   an Eulerian code the upwind advection term provides
							 *   implicit numerical dissipation, but Lagrangian cells
							 *   move with the fluid and have no advective term.
							 *   Consequence: grid-scale noise at contacts grows
							 *   unbounded → dt collapse (confirmed by KH Step 4/5).
							 *   Proof: Toro 2009 §10.2 (HLLC contact wave), LeVeque
							 *   2002 §9.5 (Godunov on linearly degenerate fields),
							 *   Springel 2010 §32 (Arepo face rest-frame motivation).
							 * USE av_mode=5 INSTEAD (face-rest-frame HLLC + mesh
							 * regularization, which is the proven Arepo recipe).
							 * No M(n,m), no NS stress, no CD10 blending.
							 * HLLC Riemann solver is the sole dissipation mechanism. */
							postype ds_mag_inv = 1.0/(facearea + 1e-30);
							postype nx_hat = dS.x * ds_mag_inv;
							postype ny_hat = dS.y * ds_mag_inv;
							postype vnL = ibp_vx*nx_hat + ibp_vy*ny_hat;
							postype vnR = jbp_vx*nx_hat + jbp_vy*ny_hat;
							postype pL = ibp_pressure;
							postype pR = jbp_pressure;

							/* MUSCL reconstruction at face midpoint.
							   tmp->x,tmp->y are stored RELATIVE to generator ibp
							   (Voro2D_FindVC uses Vec2DSub). */
							postype xf_rel_i = 0.5*(tmp->x + tmp2->x);
							postype yf_rel_i = 0.5*(tmp->y + tmp2->y);
							postype dx_ij_rel = jbp->x - ibp->x;
							postype dy_ij_rel = jbp->y - ibp->y;
							postype dx_iF = xf_rel_i;
							postype dy_iF = yf_rel_i;
							postype dx_jF = xf_rel_i - dx_ij_rel;
							postype dy_jF = yf_rel_i - dy_ij_rel;

							/* Pressure reconstruction */
							postype dp_i = ibp->stress.dPdx*dx_iF + ibp->stress.dPdy*dy_iF;
							postype dp_j = jbp->stress.dPdx*dx_jF + jbp->stress.dPdy*dy_jF;
							pL += dp_i;
							pR += dp_j;

							/* Velocity reconstruction */
							postype dvnL = (ibp->stress.gUxx*nx_hat + ibp->stress.gUxy*ny_hat)*dx_iF
							             + (ibp->stress.gUyx*nx_hat + ibp->stress.gUyy*ny_hat)*dy_iF;
							postype dvnR = (jbp->stress.gUxx*nx_hat + jbp->stress.gUxy*ny_hat)*dx_jF
							             + (jbp->stress.gUyx*nx_hat + jbp->stress.gUyy*ny_hat)*dy_jF;
							vnL += dvnL;
							vnR += dvnR;

							/* Slope limiter: clamp to prevent new extrema */
							postype pmin = fmin(ibp_pressure, jbp_pressure);
							postype pmax = fmax(ibp_pressure, jbp_pressure);
							if(pL < pmin) pL = pmin; if(pL > pmax) pL = pmax;
							if(pR < pmin) pR = pmin; if(pR > pmax) pR = pmax;
							postype vnmin = fmin(vnL - dvnL, vnR - dvnR);
							postype vnmax = fmax(vnL - dvnL, vnR - dvnR);
							if(vnL < vnmin) vnL = vnmin; if(vnL > vnmax) vnL = vnmax;
							if(vnR < vnmin) vnR = vnmin; if(vnR > vnmax) vnR = vnmax;

							/* Positivity: floor on reconstructed pressure */
							if(pL < 1e-10) pL = 1e-10;
							if(pR < 1e-10) pR = 1e-10;

							postype pst, vnst;
							hllc_face_2d(ibp_den, pL, vnL, ibp_csound,
							             jbp_den, pR, vnR, jbp_csound,
							             Gamma, &pst, &vnst);
							pi_total = pst;

							/* CD10 viscous pressure (compression only) */
							if(cd_amax > 0){
								postype dvn_hllc = vnR - vnL;
								if(dvn_hllc < 0){
									postype alpha_face = 0.5*(ibp->stress.alpha_cd + jbp->stress.alpha_cd);
									postype rho_mean = 0.5*(ibp_den + jbp_den);
									postype vsig_cd = ibp_csound + jbp_csound - fmin(0.0, dvn_hllc);
									pi_total += 0.5 * alpha_face * vsig_cd * rho_mean * (-dvn_hllc);
								}
							}

							/* Optional Monaghan AV for contact noise control */
							if(alphavis > 0){
								Voro2D_point uij;
								uij.x = jbp_vx - ibp_vx;
								uij.y = jbp_vy - ibp_vy;
								postype rvel = Vec2DDotP(&er, &uij);
								if(rvel < 0){
									postype wcomp = sqrt(ibp_rk4->w2)+sqrt(((treevorork4particletype*)jbp)->w2);
									postype scaleFactor = (wcomp==0 ? etavis: wcomp);
									postype drampScale = dramp/scaleFactor;
									postype mu = rvel/(drampScale + epsvis/drampScale);
									postype meanden = 0.5*(ibp_den + jbp_den);
									postype meanCsound = 0.5*(ibp_csound + jbp_csound);
									pi_total += (-alphavis*meanCsound*mu + betavis*mu*mu)*meanden;
								}
							}
							/* tau_dot_dS_x/y remain 0 — no NS stress */

						} else if(av_mode == 5) {
							/* ========== AREPO-LAGUERRE: face rest frame HLLC + MUSCL ========== */
							/* Pure HLLC on every (real-real) face, solved in the
							   face rest frame.  MUSCL reconstruction of P and
							   v_n using cell-averaged gradients from
							   updateDenW2Pressure2DBlend, with a pair-wise
							   clamp limiter (prevents new extrema between i and
							   j).  Pressure floor for positivity.  No M(n,m)
							   blending, no NS stress, no explicit AV — the face
							   rest frame boost is what makes contact
							   discontinuities survive without the "fake
							   compression" artefact of lab-frame Godunov.

							   E-C1 fix: face velocity is taken from
							   get2dUpqradRk4, which is the same Laguerre-aware
							   (wfrac + AA-corrected) value used by the energy
							   path below.  This keeps momentum and internal
							   energy updates consistent with each other on
							   Laguerre meshes (w2_i ≠ w2_j). */
							postype ds_mag_inv = 1.0/(facearea + 1e-30);
							postype nx_hat = dS.x * ds_mag_inv;
							postype ny_hat = dS.y * ds_mag_inv;

							/* Face velocity in lab frame: v_face = v_i + uradix_ui
							   with uradix_ui = (v_face - v_i) from get2dUpqradRk4. */
							Voro2D_point uradix_ui_boost =
								get2dUpqradRk4(ibp_rk4,
									(treevorork4particletype*)jbp, dtold);
							postype wx = ibp_vx + uradix_ui_boost.x;
							postype wy = ibp_vy + uradix_ui_boost.y;
							postype wn = wx*nx_hat + wy*ny_hat;

							postype vnL_lab = ibp_vx*nx_hat + ibp_vy*ny_hat;
							postype vnR_lab = jbp_vx*nx_hat + jbp_vy*ny_hat;

							postype pL = ibp_pressure;
							postype pR = jbp_pressure;
							postype rhoL = ibp_den;
							postype rhoR = jbp_den;

							if(use_muscl){
								/* Springel-Arepo MUSCL: reconstruct {rho, vn, P}
								   at the face centroid using Green-Gauss
								   gradients that were already scaled by the
								   per-cell Barth-Jespersen/Springel TVD limiter
								   in updateDenW2Pressure2DBlend.  NOTE: Voronoi
								   corners tmp->x,tmp->y are stored in LOCAL
								   coordinates relative to the generator ibp
								   (see Voro2D_FindVC which calls Vec2DSub).
								   So xf is already (face - ibp), and the
								   face-relative-to-j is xf - (jbp-ibp). */
								postype xf_rel_i = 0.5*(tmp->x + tmp2->x);
								postype yf_rel_i = 0.5*(tmp->y + tmp2->y);
								postype dx_ij_rel = jbp->x - ibp->x;
								postype dy_ij_rel = jbp->y - ibp->y;
								postype dx_iF = xf_rel_i;
								postype dy_iF = yf_rel_i;
								postype dx_jF = xf_rel_i - dx_ij_rel;
								postype dy_jF = yf_rel_i - dy_ij_rel;

								postype drho_i = ibp->stress.dRhodx*dx_iF + ibp->stress.dRhody*dy_iF;
								postype drho_j = jbp->stress.dRhodx*dx_jF + jbp->stress.dRhody*dy_jF;
								postype dp_i = ibp->stress.dPdx*dx_iF + ibp->stress.dPdy*dy_iF;
								postype dp_j = jbp->stress.dPdx*dx_jF + jbp->stress.dPdy*dy_jF;
								postype dvxL = ibp->stress.gUxx*dx_iF + ibp->stress.gUxy*dy_iF;
								postype dvyL = ibp->stress.gUyx*dx_iF + ibp->stress.gUyy*dy_iF;
								postype dvnL0 = nx_hat*dvxL + ny_hat*dvyL;
								postype dvxR = jbp->stress.gUxx*dx_jF + jbp->stress.gUxy*dy_jF;
								postype dvyR = jbp->stress.gUyx*dx_jF + jbp->stress.gUyy*dy_jF;
								postype dvnR0 = nx_hat*dvxR + ny_hat*dvyR;

								rhoL += drho_i;
								rhoR += drho_j;
								pL += dp_i;
								pR += dp_j;
								vnL_lab += dvnL0;
								vnR_lab += dvnR0;
							}

							/* Positivity floor */
							if(pL < 1e-10) pL = 1e-10;
							if(pR < 1e-10) pR = 1e-10;
							if(rhoL < 1e-10) rhoL = 1e-10;
							if(rhoR < 1e-10) rhoR = 1e-10;

							/* Face rest frame HLLC */
							postype pst, vnst_lab;
							hllc_face_2d_rest_frame(rhoL, pL, vnL_lab, ibp_csound,
							                        rhoR, pR, vnR_lab, jbp_csound,
							                        wn, Gamma, &pst, &vnst_lab);
							pi_total = pst;

							/* Optional Monaghan AV for grid-scale noise control
							   (uses GAS_AlphaVis; disabled when alphavis==0).
							   Applied only on compression to damp cell jitter. */
							if(alphavis > 0){
								Voro2D_point uij;
								uij.x = jbp_vx - ibp_vx;
								uij.y = jbp_vy - ibp_vy;
								postype rvel = Vec2DDotP(&er, &uij);
								if(rvel < 0){
									postype wcomp = sqrt(ibp_rk4->w2)+sqrt(((treevorork4particletype*)jbp)->w2);
									postype scaleFactor = (wcomp==0 ? etavis: wcomp);
									postype drampScale = dramp/scaleFactor;
									postype mu = rvel/(drampScale + epsvis/drampScale);
									postype meanden = 0.5*(ibp_den + jbp_den);
									postype meanCsound = 0.5*(ibp_csound + jbp_csound);
									pi_total += (-alphavis*meanCsound*mu + betavis*mu*mu)*meanden;
								}
							}
							/* tau_dot_dS_x/y remain 0 — no NS stress */

						} else {
							/* Two-tier blended path (real-real faces only) */
							postype OoA = OrderOfAccuracy;

							/* === Tier 1: M(n,m) pressure + NS stress === */
							postype p_mnm = VoroStressRK4_Pressure2D(
									ibp_pressure, jbp_pressure, tmp, neighwork, simpar, RT);

							/* Viscous traction vector: τ·dS (full tensor, not just normal projection) */
							tau_dot_dS_x = 0; tau_dot_dS_y = 0;
							if(av_mode == 1 || (av_mode == 0 && nu_phys > 0)){
								/* M(n,m) face-average of τ components */
								postype txx_face, txy_face, tyy_face;
								if(OoA > 0 && (tmp->upperlink)->upperrelated >= 0 && tmp->lowerrelated >= 0){
									txx_face = VoroStressRK4_Scalar2D(
										ibp->stress.tauxx, jbp->stress.tauxx, tmp, neighwork, simpar, RT, stress.tauxx);
									txy_face = VoroStressRK4_Scalar2D(
										ibp->stress.tauxy, jbp->stress.tauxy, tmp, neighwork, simpar, RT, stress.tauxy);
									tyy_face = VoroStressRK4_Scalar2D(
										ibp->stress.tauyy, jbp->stress.tauyy, tmp, neighwork, simpar, RT, stress.tauyy);
								} else {
									txx_face = 0.5*(ibp->stress.tauxx + jbp->stress.tauxx);
									txy_face = 0.5*(ibp->stress.tauxy + jbp->stress.tauxy);
									tyy_face = 0.5*(ibp->stress.tauyy + jbp->stress.tauyy);
								}
								/* Viscous traction: τ_code = -τ_NS, so
								   τ_NS·dS = -(τ_code·dS).  Negate here so that
								   subsequent "+tau_dot_dS" gives correct
								   +τ_NS·dS direction (viscosity damps shear). */
								tau_dot_dS_x = -(txx_face * dS.x + txy_face * dS.y);
								tau_dot_dS_y = -(txy_face * dS.x + tyy_face * dS.y);
							}

							if(av_mode == 0 || av_mode == 1){
								/* M(n,m) pressure + NS stress + optional Monaghan AV */
								pi_total = p_mnm;
								if(alphavis > 0){
									Voro2D_point uij;
									uij.x = jbp_vx - ibp_vx;
									uij.y = jbp_vy - ibp_vy;
									postype rvel = Vec2DDotP(&er, &uij);
									if(rvel < 0){
										postype wcomp = sqrt(ibp_rk4->w2)+sqrt(((treevorork4particletype*)jbp)->w2);
										postype scaleFactor = (wcomp==0 ? etavis: wcomp);
										postype drampScale = dramp/scaleFactor;
										postype mu = rvel/(drampScale + epsvis/drampScale);
										postype meanden = 0.5*(ibp_den + jbp_den);
										postype meanCsound = 0.5*(ibp_csound + jbp_csound);
										pi_total += (-alphavis*meanCsound*mu + betavis*mu*mu)*meanden;
									}
								}
								/* tau_dot_dS already computed above */
							} else {
								/* === av_mode == 2: Two-tier blend ===
								 * Revision 4 (expert review):
								 *   FIX1: use hllc_face_2d_rest_frame (Springel 2010 §32).
								 *         Lab-frame HLLC at Lagrangian contacts gives ~zero
								 *         dissipation and produces surface-tension artifacts.
								 *   FIX2: gate CD10 viscous pressure when f_pq > 0.5
								 *         (HLLC alone provides dissipation; avoid double
                                 *         dissipation positive feedback).
								 *   FIX3: contact-aware MUSCL — skip velocity reconstruction
								 *         when density jump > 10% but pressure continuous
								 *         (< 2%). Shear gradients at smooth contacts otherwise
								 *         create spurious relative v_n for HLLC to dissipate.
								 */

								/* Blending factor: symmetric f_pq = f_qp */
								postype alpha_max_pq = fmax(ibp->stress.alpha_cd, jbp->stress.alpha_cd);
								postype theta_pq = 0;
								if(ibp->stress.divv < -1e-10 || jbp->stress.divv < -1e-10)
									theta_pq = blend_theta;
								postype f_pq = fmin(1.0, alpha_max_pq/(cd_amax + 1e-30) + theta_pq);

								postype p_hllc = p_mnm;  /* fallback */
								postype Pi_cd10 = 0;

								if(f_pq > 1e-6){
									/* HLLC path */
									postype vnL, vnR, pL, pR, rhoL, rhoR, cL, cR;

									/* Project velocities onto face normal */
									postype ds_mag_inv = 1.0/(facearea + 1e-30);
									postype nx_hat = dS.x * ds_mag_inv;
									postype ny_hat = dS.y * ds_mag_inv;
									vnL = ibp_vx*nx_hat + ibp_vy*ny_hat;
									vnR = jbp_vx*nx_hat + jbp_vy*ny_hat;

									/* FIX3: detect pure contact (ρ jump, P continuous).
									 *       At such faces, velocity MUSCL must be
									 *       suppressed to avoid spurious dvn from
									 *       transverse shear gradient. */
									postype rho_mean_det = 0.5*(ibp_den + jbp_den);
									postype p_mean_det   = 0.5*(ibp_pressure + jbp_pressure);
									postype drho_rel = fabs(ibp_den - jbp_den) /
									                   (rho_mean_det + 1e-30);
									postype dp_rel   = fabs(ibp_pressure - jbp_pressure) /
									                   (p_mean_det + 1e-30);
									#define CONTACT_DRHO_THRESH 0.1
									#define CONTACT_DP_THRESH   0.02
									int is_pure_contact = (drho_rel > CONTACT_DRHO_THRESH) && (dp_rel < CONTACT_DP_THRESH);

									if(use_muscl){
										/* MUSCL reconstruction using cell-averaged gradients.
										 * ∇P, ∇v computed via Green-Gauss in updateDenW2Pressure2DBlend.
										 * Reconstruct at face midpoint: Q_L = Q̄_i + ⟨∇Q⟩_i · (x_f - x_i).
										 * tmp->x,tmp->y are stored RELATIVE to generator ibp
										 * (Voro2D_FindVC uses Vec2DSub). */
										postype xf_rel_i = 0.5*(tmp->x + tmp2->x);
										postype yf_rel_i = 0.5*(tmp->y + tmp2->y);
										postype dx_ij_rel = jbp->x - ibp->x;
										postype dy_ij_rel = jbp->y - ibp->y;

										/* Displacement: cell center → face midpoint */
										postype dx_iF = xf_rel_i;
										postype dy_iF = yf_rel_i;
										postype dx_jF = xf_rel_i - dx_ij_rel;
										postype dy_jF = yf_rel_i - dy_ij_rel;

										/* Pressure reconstruction */
										postype dp_i = ibp->stress.dPdx*dx_iF + ibp->stress.dPdy*dy_iF;
										postype dp_j = jbp->stress.dPdx*dx_jF + jbp->stress.dPdy*dy_jF;
										pL = ibp_pressure + dp_i;
										pR = jbp_pressure + dp_j;

										/* Velocity reconstruction (project gradient onto face normal).
										 * FIX3: skip at pure contacts — shear gradient
										 * ∇v is nonzero even with v·n̂ = 0, which fabricates
										 * relative v_n that HLLC would then dissipate. */
										postype dvnL = 0, dvnR = 0;
										if(!is_pure_contact){
											dvnL = (ibp->stress.gUxx*nx_hat + ibp->stress.gUxy*ny_hat)*dx_iF
											     + (ibp->stress.gUyx*nx_hat + ibp->stress.gUyy*ny_hat)*dy_iF;
											dvnR = (jbp->stress.gUxx*nx_hat + jbp->stress.gUxy*ny_hat)*dx_jF
											     + (jbp->stress.gUyx*nx_hat + jbp->stress.gUyy*ny_hat)*dy_jF;
											vnL += dvnL;
											vnR += dvnR;
										}

										/* Simple limiter: clamp to prevent new extrema */
										postype pmin = fmin(ibp_pressure, jbp_pressure);
										postype pmax = fmax(ibp_pressure, jbp_pressure);
										if(pL < pmin) pL = pmin; if(pL > pmax) pL = pmax;
										if(pR < pmin) pR = pmin; if(pR > pmax) pR = pmax;
										postype vnmin = fmin(vnL - dvnL, vnR - dvnR);
										postype vnmax = fmax(vnL - dvnL, vnR - dvnR);
										if(vnL < vnmin) vnL = vnmin; if(vnL > vnmax) vnL = vnmax;
										if(vnR < vnmin) vnR = vnmin; if(vnR > vnmax) vnR = vnmax;
									}

									pL = use_muscl ? pL : ibp_pressure;
									pR = use_muscl ? pR : jbp_pressure;
									rhoL = ibp_den;
									rhoR = jbp_den;
									cL = ibp_csound;
									cR = jbp_csound;

									/* FIX1: Face rest frame HLLC (Arepo-style).
									 *       w_n = ((v_i + v_j)/2) · n̂ boosts the Riemann
									 *       problem into the face's rest frame, where
									 *       contact discontinuities show proper relative
									 *       velocity. Unboost pstar (invariant) and
									 *       vnstar (lab = rest + w_n). */
									postype wx = 0.5*(ibp_vx + jbp_vx);
									postype wy = 0.5*(ibp_vy + jbp_vy);
									postype wn = wx*nx_hat + wy*ny_hat;

									postype pst, vnst_lab;
									hllc_face_2d_rest_frame(rhoL, pL, vnL, cL,
									                        rhoR, pR, vnR, cR,
									                        wn, Gamma,
									                        &pst, &vnst_lab);
									p_hllc = pst;

									/* FIX2: CD10 viscous pressure is active ONLY when
									 * HLLC weight is low (f_pq <= 0.5). When HLLC
									 * dominates the face it provides its own
									 * dissipation; stacking CD10 on top causes
									 * double heating → dt collapse positive feedback. */
									if(f_pq <= 0.5){
										postype dvn = vnR - vnL;
										if(dvn < 0){ /* approaching (in face rest frame same sign) */
											postype alpha_face = 0.5*(ibp->stress.alpha_cd + jbp->stress.alpha_cd);
											postype rho_mean = 0.5*(rhoL + rhoR);
											postype vsig = cL + cR - fmin(0.0, dvn);
											Pi_cd10 = 0.5 * alpha_face * vsig * rho_mean * (-dvn);
										}
									}
								}

								/* Blend: pressure is blended, NS traction scaled by (1-f) */
								pi_total = (1.0-f_pq)*p_mnm + f_pq*(p_hllc + Pi_cd10);
								tau_dot_dS_x *= (1.0-f_pq);
								tau_dot_dS_y *= (1.0-f_pq);
							}
						}

						/* === Accumulate forces and energy rates === */
						/* Internal energy: pressure work + viscous heating + heat conduction */
						/* get2dUpqradRk4 returns v_face - v_i (Laguerre wfrac + AA correction) */
						Voro2D_point uradix_ui = get2dUpqradRk4(ibp_rk4, (treevorork4particletype*)jbp, dtold);
						/* Pressure work: -p * (v_face - v_i) · dS */
						die += -pi_total * Vec2DDotP(&uradix_ui, &dS);
						/* Viscous heating: (τ·dS) · (v_face - v_i) */
						die += tau_dot_dS_x * uradix_ui.x + tau_dot_dS_y * uradix_ui.y;

						/* Heat conduction: Q = χ ρ_face (Tj-Ti)/d_ij * facearea */
						if(nu_phys > 0 && prandtl > 0 && !jbp_is_ghost){
							postype chi = nu_phys / prandtl;
							postype Ti = ibp_pressure / ibp_den;
							postype Tj = jbp_pressure / jbp_den;
							postype rho_face = 0.5*(ibp_den + jbp_den);
							die += chi * rho_face * (Tj - Ti) / dramp * facearea;
						}

						/* Total energy: ua = v_face = v_i + (v_face - v_i)
						   ensures dte = die + dke consistency */
						Voro2D_point ua;
						ua.x = ibp_vx + uradix_ui.x;
						ua.y = ibp_vy + uradix_ui.y;
						dte += -pi_total * Vec2DDotP(&ua, &dS)
						     + tau_dot_dS_x * ua.x + tau_dot_dS_y * ua.y;

						/* Kinetic energy */
						Voro2D_point ub;
						ub.x = ibp_vx;
						ub.y = ibp_vy;
						dke += -pi_total * Vec2DDotP(&ub, &dS)
						     + tau_dot_dS_x * ub.x + tau_dot_dS_y * ub.y;

						/* Force: -p·dS + τ·dS */
						fx += -pi_total * dS.x + tau_dot_dS_x;
						fy += -pi_total * dS.y + tau_dot_dS_y;

						/* Signal velocity for CFL and CD10 */
						Voro2D_point dv;
						dv.x = jbp_vx - ibp_vx;
						dv.y = jbp_vy - ibp_vy;
						postype VdotR = Vec2DDotP(&dv, &er);
						postype vsig = jbp_csound + ibp_csound - MIN(0, VdotR);
						/* CFL floor: prevent dt→0 from particle clustering
						   or geometric proximity of ghost mirrors */
						postype heff = 0.25*sqrt(ibp_rk4->volume);
						postype dramp_cfl = fmax(dramp, heff);
						postype dt = 2*Courant*dramp_cfl/vsig;

						if(isnan(dt)){
							DEBUGPRINT("P%d blend: nan dt %d p%ld xy= %g %g : j=%ld jxy= %g %g dramp= %g vsig= %g cs_i= %g cs_j= %g pi_tot= %g\n",
									MYID(simpar), i, (long)PINDX(p+i), p[i].x, p[i].y,
									(long)PINDX((treevorork4particletype*)jbp), jbp->x, jbp->y,
									dramp, vsig, ibp_csound, jbp_csound, pi_total);
							dt = 1e-10;
						}
						if(dt < 1e-6){
							fprintf(stderr,"[HYDRO_CFL] P%d i=%d dt=%g dramp=%g dramp_cfl=%g vsig=%g csound_i=%g csound_j=%g vol=%g x=%g y=%g\n",
								MYID(simpar), i, dt, dramp, dramp_cfl, vsig, ibp_csound, jbp_csound, ibp_rk4->volume, ibp_rk4->x, ibp_rk4->y);
						}
						ibp_rk4->dt = MIN(ibp_rk4->dt, dt);
						if(dt < Dtime) Dtime = dt;

						/* Track vsig_max for CD10 (home cell only).
						   The face i-j also appears in j's face list where j
						   is home; by face-symmetry of vsig, j's vsig_max
						   will be set there. Writing only ibp->stress.vsig_max
						   removes the OMP critical-section hot spot. */
						if(av_mode >= 1){
							/* For ghost faces, use sound-speed-only signal velocity
							   (reversed velocity inflates the full vsig artificially) */
							postype vsig_cd = jbp_is_ghost ?
								(ibp_csound + jbp_csound) : vsig;
							if(vsig_cd > ibp->stress.vsig_max) ibp->stress.vsig_max = vsig_cd;
						}
					}
					tmp = tmp->upperlink;
				} while(tmp != vorocorner);

				/* Viscous CFL: dt_visc = 0.5 h^2 / max(nu_phys, nu_cd, chi)
				   Evaluated once per particle (depends only on i's state). */
				if(av_mode >= 1 || nu_phys > 0){
					postype h_i = sqrt(ibp_rk4->volume);
					postype nu_cd_i = ibp->stress.alpha_cd * h_i * ibp_csound;
					postype chi = (nu_phys > 0 && prandtl > 0) ? nu_phys / prandtl : 0;
					postype nu_eff = fmax(fmax(nu_phys, nu_cd_i), chi);
					if(nu_eff > 0){
						postype dt_visc = 0.5 * h_i * h_i / nu_eff;
						if(dt_visc < 1e-6){
							fprintf(stderr,"[VISC_CFL] P%d i=%d dt_visc=%g h_i=%g nu_eff=%g alpha_cd=%g csound=%g den=%g P=%g ie=%g vol=%g x=%g y=%g\n",
								MYID(simpar), i, dt_visc, h_i, nu_eff, ibp->stress.alpha_cd, ibp_csound, ibp_den, ibp_pressure, ibp_rk4->ie, ibp_rk4->volume, ibp_rk4->x, ibp_rk4->y);
						}
						ibp_rk4->dt = MIN(ibp_rk4->dt, dt_visc);
						if(dt_visc < Dtime) Dtime = dt_visc;
					}
				}

				ibp_rk4->die = die;
				ibp_rk4->ax = fx/ibp_rk4->mass;
				ibp_rk4->ay = fy/ibp_rk4->mass;

				if(isnan(fx)){
					DEBUGPRINT("P%d blend: nan fx %d : %g %g %ld p= %g\n",
							MYID(simpar), i, ibp->x, ibp->y, PINDX(ibp), ibp_pressure);
					exit(0);
				}
			}
			free(p); free(neighbors); free(neighwork);
		}
		free(vorocorner);
	}
	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
	{
		postype TDtime;
		MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
		Dtime = TDtime;
	}
	return Dtime;
}

/* ================================================================
 *  getAccVoro2D_LagMFM (av_mode=4):
 *    Meshless force computation using MFM effective faces
 *    (Hopkins 2015 GIZMO).  Lagrangian limit: mass per particle is
 *    fixed, face flux carries pressure work only.
 *
 *    Effective face:
 *      A_ij^α = V_i² W_ij · (Ẽ_inv_i^{αβ} dx^β)
 *             + V_j² W_ji · (Ẽ_inv_j^{αβ} dx^β)
 *    (where Ẽ_inv is the inverse of the unnormalized E = Σ W dx⊗dx
 *     already stored in stress.E_inv_xx..yy by updateDenW2Pressure2D_LagMFM.)
 *
 *    Riemann:   pstar, vnstar from hllc_face_2d (first-order).
 *    Force on i: fx += -pstar * A_ij_x + τ·A_ij_x (viscous)
 *    Energy on i: die += -pstar * (vnstar - vn_i) * |A_ij|
 * ================================================================ */
double getAccVoro2D_LagMFM(SimParameters *simpar, postype xmin, postype ymin,
		postype xmax, postype ymax,
		postype OrderOfAccuracy, postype Courant, postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype))
		){
	treevorostressrk4particletype *bp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype Dtime = 1.e10;
	postype nu_phys = GAS_VISCOSITY(simpar);
	float alphavis = GAS_AlphaVis(simpar);
	float betavis  = GAS_BetaVis(simpar);
	float etavis   = GAS_ETAVIS(simpar) * SIMBOX(simpar).x.max / NX(simpar);
	float epsvis   = GAS_EPSVIS(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	int mx, my;
	BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
	BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	CellType *cells = (VORO_BASICCELL(simpar) = (CellType*)my_malloc(sizeof(CellType)*mx*my));
	mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax, paddingAllTreeParticles);

	postype invcs = 1.0/cellsize;

	int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
	for(iy=0; iy<my; iy++){
		int ix;
		for(ix=0; ix<mx; ix++){
			struct linkedlisttype *li = cells[ix + mx*iy].link;
			while(li){
				treevorostressrk4particletype *ibp = (treevorostressrk4particletype*)li;
				int iwall = (PINDX((treevorork4particletype*)ibp) == MAX_INDEX);
				if(iwall || IS_FLAG_ON(li,BoundaryGhostflag)){
					li = li->next;
					continue;
				}

				postype ibp_vx  = ibp->vx;
				postype ibp_vy  = ibp->vy;
				postype ibp_cs  = ibp->csound;
				postype ibp_P   = ibp->pressure;
				postype ibp_rho = ibp->den;
				postype Vi      = ibp->volume;
				postype h_i     = ibp->stress.h_mfm;
				if(h_i <= 0) h_i = LAGMFM_ETA * sqrt(Vi > 0 ? Vi : cellsize*cellsize);
				postype Ei_xx = ibp->stress.E_inv_xx;
				postype Ei_xy = ibp->stress.E_inv_xy;
				postype Ei_yx = ibp->stress.E_inv_yx;
				postype Ei_yy = ibp->stress.E_inv_yy;

				postype fx=0, fy=0, die=0;
				postype dt_cell = 1.e10;

				/* Neighbor search: up to 2 max(h_i, h_neighbor_max).
				 * We don't know h_j a priori, so use 2*h_i and validate via
				 * per-pair 2*max(h_i,h_j) cutoff inside the inner loop. */
				postype two_h_i = 2.0 * h_i;
				int nbr = (int)ceil(two_h_i * invcs) + 1;
				if(nbr < 1) nbr = 1;
				if(nbr > 10) nbr = 10;

				int jy;
				for(jy = iy-nbr; jy <= iy+nbr; jy++){
					if(jy < 0 || jy >= my) continue;
					int jx;
					for(jx = ix-nbr; jx <= ix+nbr; jx++){
						if(jx < 0 || jx >= mx) continue;
						struct linkedlisttype *lj = cells[jx + mx*jy].link;
						while(lj){
							treevorostressrk4particletype *jbp =
									(treevorostressrk4particletype*)lj;
							if((treevorork4particletype*)jbp == (treevorork4particletype*)ibp){
								lj = lj->next;
								continue;
							}
							int jwall = (PINDX((treevorork4particletype*)jbp) == MAX_INDEX);

							postype dx = jbp->x - ibp->x;
							postype dy = jbp->y - ibp->y;
							postype r2 = dx*dx + dy*dy;

							postype h_j = jbp->stress.h_mfm;
							if(h_j <= 0) h_j = h_i;
							postype h_max = h_i > h_j ? h_i : h_j;
							postype two_h_max = 2.0*h_max;
							if(r2 > two_h_max*two_h_max){
								lj = lj->next;
								continue;
							}
							postype r = sqrt(r2);

							/* Effective face area vector (antisymmetric in i↔j).
							 * Wall-ghost: cannot use j's E_inv (unmirrored), so
							 * use only the i-side contribution. */
							postype W_ij = mfm_W_wendland2d(r, h_i);
							postype W_ji = mfm_W_wendland2d(r, h_j);
							postype Vi2 = Vi*Vi;
							postype Vj  = jbp->volume;
							postype Vj2 = Vj*Vj;
							postype Ej_xx = jbp->stress.E_inv_xx;
							postype Ej_xy = jbp->stress.E_inv_xy;
							postype Ej_yx = jbp->stress.E_inv_yx;
							postype Ej_yy = jbp->stress.E_inv_yy;

							postype Ax_i = Vi2 * W_ij * (Ei_xx*dx + Ei_xy*dy);
							postype Ay_i = Vi2 * W_ij * (Ei_yx*dx + Ei_yy*dy);
							postype Ax_j, Ay_j;
							if(jwall){
								Ax_j = 0;
								Ay_j = 0;
							} else {
								Ax_j = Vj2 * W_ji * (Ej_xx*dx + Ej_xy*dy);
								Ay_j = Vj2 * W_ji * (Ej_yx*dx + Ej_yy*dy);
							}
							postype Ax = Ax_i + Ax_j;
							postype Ay = Ay_i + Ay_j;
							postype Amag2 = Ax*Ax + Ay*Ay;
							if(Amag2 < 1e-30){
								lj = lj->next;
								continue;
							}
							postype Amag = sqrt(Amag2);
							postype nx_hat = Ax / Amag;
							postype ny_hat = Ay / Amag;

							postype jbp_vx = jbp->vx;
							postype jbp_vy = jbp->vy;
							postype jbp_cs = jbp->csound;
							postype jbp_P  = jbp->pressure;
							postype jbp_rho = jbp->den;

							/* ============================================
							 * D-B2 MUSCL reconstruction.
							 *
							 * Face midpoint approximation:
							 *   x_face ≈ 0.5*(x_i + x_j)
							 * so dx_iF = 0.5*dx, dx_jF = -0.5*dx.
							 *
							 * Reconstruct P, v at x_face from each side
							 * using the matrix-weighted gradients stored in
							 * stress (computed in updateDenW2Pressure2D_LagMFM):
							 *   P_L = P_i + ∇P_i · (x_face - x_i)
							 *   P_R = P_j + ∇P_j · (x_face - x_j)
							 * and analogously for vx, vy.
							 *
							 * Density is NOT reconstructed (ρ_L = ρ_i, ρ_R = ρ_j)
							 * to avoid needing ∇ρ storage — same simplification
							 * as av_mode=5 (documented as E-C2 future work).
							 *
							 * Pair-wise clamp limiter (Barth-Jespersen style
							 * restricted to the face pair): the reconstructed
							 * value is clamped into [min(Q_i,Q_j), max(Q_i,Q_j)]
							 * so that no new extrema appear at the face.
							 * ============================================ */
							postype hdx = 0.5*dx;
							postype hdy = 0.5*dy;

							postype pL_rec = ibp_P + ibp->stress.dPdx*hdx + ibp->stress.dPdy*hdy;
							postype pR_rec = jbp_P - jbp->stress.dPdx*hdx - jbp->stress.dPdy*hdy;
							postype vxL_rec = ibp_vx + ibp->stress.gUxx*hdx + ibp->stress.gUxy*hdy;
							postype vxR_rec = jbp_vx - jbp->stress.gUxx*hdx - jbp->stress.gUxy*hdy;
							postype vyL_rec = ibp_vy + ibp->stress.gUyx*hdx + ibp->stress.gUyy*hdy;
							postype vyR_rec = jbp_vy - jbp->stress.gUyx*hdx - jbp->stress.gUyy*hdy;

							/* Pair-wise clamp (no new extrema between i and j) */
							#define LAGMFM_CLAMP(Qrec,Qi,Qj) do { \
								postype qmn = (Qi) < (Qj) ? (Qi) : (Qj); \
								postype qmx = (Qi) > (Qj) ? (Qi) : (Qj); \
								if(Qrec < qmn) Qrec = qmn; \
								if(Qrec > qmx) Qrec = qmx; \
							} while(0)
							LAGMFM_CLAMP(pL_rec,  ibp_P,  jbp_P);
							LAGMFM_CLAMP(pR_rec,  ibp_P,  jbp_P);
							LAGMFM_CLAMP(vxL_rec, ibp_vx, jbp_vx);
							LAGMFM_CLAMP(vxR_rec, ibp_vx, jbp_vx);
							LAGMFM_CLAMP(vyL_rec, ibp_vy, jbp_vy);
							LAGMFM_CLAMP(vyR_rec, ibp_vy, jbp_vy);
							#undef LAGMFM_CLAMP

							/* Wall ghost: j-side gradients are not mirrored,
							 * so skip reconstruction on wall faces (use central
							 * state for the non-wall side only). */
							if(jwall){
								pL_rec  = ibp_P;
								vxL_rec = ibp_vx;
								vyL_rec = ibp_vy;
								pR_rec  = jbp_P;
								vxR_rec = jbp_vx;
								vyR_rec = jbp_vy;
							}

							/* Project reconstructed velocities to the normal */
							postype vnL_lab = vxL_rec*nx_hat + vyL_rec*ny_hat;
							postype vnR_lab = vxR_rec*nx_hat + vyR_rec*ny_hat;

							postype pL = pL_rec;
							postype pR = pR_rec;
							postype rhoL = ibp_rho;
							postype rhoR = jbp_rho;
							postype cL = ibp_cs;
							postype cR = jbp_cs;

							/* ============================================
							 * D-B3 Face rest frame HLLC.
							 *
							 * Lagrangian LagMFM: the effective face moves with
							 * the average particle velocity.  The natural face
							 * velocity is
							 *   w = 0.5*(v_i + v_j)         (lab frame)
							 * so that in the face rest frame both sides have
							 * small relative velocities, yielding small HLLC
							 * dissipation at contact discontinuities.
							 * ============================================ */
							postype wx = 0.5*(ibp_vx + jbp_vx);
							postype wy = 0.5*(ibp_vy + jbp_vy);
							postype wn = wx*nx_hat + wy*ny_hat;

							postype pstar, vnstar_lab;
							if(jwall){
								/* Wall face: bypass Riemann, use central
								 * pressure and v_n = 0 (wall enforces no
								 * normal flux). */
								pstar       = 0.5*(pL + pR);
								vnstar_lab  = 0;
							} else {
								hllc_face_2d_rest_frame(rhoL, pL, vnL_lab, cL,
								                        rhoR, pR, vnR_lab, cR,
								                        wn, Gamma,
								                        &pstar, &vnstar_lab);
							}

							/* Optional Monaghan AV (small, for grid noise) */
							if(alphavis > 0 && !jwall){
								postype uij_x = jbp_vx - ibp_vx;
								postype uij_y = jbp_vy - ibp_vy;
								postype er_x  = dx/r;
								postype er_y  = dy/r;
								postype rvel  = er_x*uij_x + er_y*uij_y;
								if(rvel < 0){
									postype wcomp = sqrt(ibp->w2) +
									                sqrt(((treevorork4particletype*)jbp)->w2);
									postype scaleFactor = (wcomp==0 ? etavis : wcomp);
									postype drampScale  = r/scaleFactor;
									postype mu = rvel/(drampScale + epsvis/drampScale);
									postype meanden = 0.5*(ibp_rho + jbp_rho);
									postype meanCs  = 0.5*(ibp_cs  + jbp_cs);
									pstar += (-alphavis*meanCs*mu + betavis*mu*mu)*meanden;
								}
							}

							/* NS viscous traction: τ_face · A_ij (full tensor).
							 * τ in stress.tauxx..tauyy is -νρS (our sign). */
							postype tau_A_x = 0, tau_A_y = 0;
							if(!jwall){
								postype txx_f = 0.5*(ibp->stress.tauxx + jbp->stress.tauxx);
								postype txy_f = 0.5*(ibp->stress.tauxy + jbp->stress.tauxy);
								postype tyy_f = 0.5*(ibp->stress.tauyy + jbp->stress.tauyy);
								/* τ_NS = -τ_code */
								tau_A_x = -(txx_f*Ax + txy_f*Ay);
								tau_A_y = -(txy_f*Ax + tyy_f*Ay);
							}

							/* Force: -pstar * A + τ_NS·A */
							fx += -pstar * Ax + tau_A_x;
							fy += -pstar * Ay + tau_A_y;

							/* ============================================
							 * D-B4 Full face velocity + viscous heating.
							 *
							 * Newton 3rd law & energy conservation check:
							 *
							 *   Per face,  A_ij is antisymmetric (A_ji = -A_ij)
							 *   and τ_face, pstar are symmetric in (i,j).
							 *   Hence  F_ij + F_ji = 0  (checked).
							 *
							 *   Total viscous heating per face:
							 *     die_i + die_j = (τ_NS·A_ij)·(v_face - v_i)
							 *                   + (τ_NS·A_ji)·(v_face - v_j)
							 *                   = (τ_NS·A_ij)·(v_j - v_i).
							 *   Total kinetic-energy change per face:
							 *     dKE_i + dKE_j = F_ij·v_i + F_ji·v_j
							 *                   = (τ_NS·A_ij)·(v_i - v_j).
							 *   Sum: 0  ⇒ energy conserved.
							 *   Since τ_NS is symmetric negative semidefinite
							 *   on the shear it dissipates v_j-v_i in the
							 *   correct direction — viscous heating ≥ 0.
							 *
							 * Full face velocity (normal + tangential):
							 *   v_face = vnstar · n̂ + v_t · t̂
							 * where the tangential component is approximated
							 * by the average (t̂ = (-ny_hat, nx_hat)):
							 *   v_t = 0.5*((v_i + v_j)·t̂).
							 * This captures tangential HLLC passive advection
							 * and yields correct viscous work along shear.
							 * ============================================ */
							postype tx_hat = -ny_hat;
							postype ty_hat =  nx_hat;
							postype vt_avg = 0.5*((ibp_vx + jbp_vx)*tx_hat
							                     +(ibp_vy + jbp_vy)*ty_hat);
							postype vfx = vnstar_lab * nx_hat + vt_avg * tx_hat;
							postype vfy = vnstar_lab * ny_hat + vt_avg * ty_hat;

							/* D-B6 Wall slip option:
							 *   Free-slip (default):  v_face keeps the
							 *     particle's tangential velocity
							 *     (approximated by the avg above, which for
							 *      a wall pair is 0.5*v_i_tan in the limit
							 *      v_j = mirrored v_i_tan).
							 *   No-slip (LAGMFM_WALL_NOSLIP=1): zero the
							 *     face velocity entirely at walls. */
							if(jwall){
#if LAGMFM_WALL_NOSLIP
								vfx = 0;
								vfy = 0;
#else
								/* Free-slip: tangential = v_i, normal = 0 */
								postype vti = ibp_vx*tx_hat + ibp_vy*ty_hat;
								vfx = vti*tx_hat;
								vfy = vti*ty_hat;
#endif
							}

							die += -pstar*((vfx - ibp_vx)*Ax + (vfy - ibp_vy)*Ay);
							/* Viscous heating: (τ_NS·A) · (v_face - v_i) */
							die += tau_A_x*(vfx - ibp_vx) + tau_A_y*(vfy - ibp_vy);

							/* Signal velocity for CFL */
							postype dvx = jbp_vx - ibp_vx;
							postype dvy = jbp_vy - ibp_vy;
							postype er_x = dx/r;
							postype er_y = dy/r;
							postype VdotR = dvx*er_x + dvy*er_y;
							postype vsig = cL + cR - MIN(0, VdotR);
							postype heff = 0.25*sqrt(Vi);
							postype dramp_cfl = fmax(r, heff);
							postype dt_face = 2.0*Courant*dramp_cfl/vsig;
							if(dt_face < dt_cell) dt_cell = dt_face;

							/* Track vsig_max for CD10 */
							if(!jwall){
								if(vsig > ibp->stress.vsig_max)
									ibp->stress.vsig_max = vsig;
							}

							lj = lj->next;
						}
					}
				}

				/* Viscous CFL */
				{
					postype h_eff = sqrt(Vi);
					postype nu_cd = ibp->stress.alpha_cd * h_eff * ibp_cs;
					postype nu_eff = fmax(nu_phys, nu_cd);
					if(nu_eff > 0){
						postype dt_visc = 0.5*h_eff*h_eff/nu_eff;
						if(dt_visc < dt_cell) dt_cell = dt_visc;
					}
				}

				ibp->die = die / ibp->mass;  /* store as rate per unit mass? */
				/* Actually die accumulated above is m_i * de/dt.  For consistency
				 * with the blend code convention (ibp_rk4->die = dE/dt), store
				 * die as the total rate, not per-mass.  We must match the
				 * integrator which does: ie += die*Dtime. */
				ibp->die = die;
				ibp->ax  = fx / ibp->mass;
				ibp->ay  = fy / ibp->mass;
				ibp->dt  = dt_cell;
				if(dt_cell < Dtime) Dtime = dt_cell;

				li = li->next;
			}
		}
	}

	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
	{
		postype TDtime;
		MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
		Dtime = TDtime;
	}
	return Dtime;
}

double getAccVoro2D_rt(SimParameters *simpar, postype xmin, postype ymin,
		postype xmax, postype ymax,
		postype OrderOfAccuracy, postype Courant, postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype))
		){
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int isave = -1;
    postype Dtime = 1.e10;
	int nx=NX(simpar);


	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	postype cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar);
    int mx, my;
    BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
    BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
	mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);


	postype dtold = GAS_dtold(simpar);



	float alphavis = GAS_AlphaVis(simpar);
	float betavis = GAS_BetaVis(simpar);
	float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
	float epsvis = GAS_EPSVIS(simpar);


//	DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);

    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=0;iy<my;iy++){
    	int mp=1000;
    	Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
    	postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
        for(ix=0;ix<mx;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
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
					/*
					if(tmp->upperrelated <0) {
						tmp = vorocorner;
						do { 
							treevorork4particletype *jbp = (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
							DEBUGPRINT("P%d has p%d :neighbor %d : %g %g : %g %g\n", 
									MYID(simpar), center.indx, tmp->upperrelated, ibp->x,ibp->y,jbp->x, jbp->y);
						} while(tmp);
						exit(99);
					}
					*/
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
								ibp_pressure,jbp_pressure, tmp, neighwork, simpar, LAG);

						Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

						postype dramp = sqrt(Vec2DDotP(&dr, &dr));
						Voro2D_point er;
						er.x = dr.x/dramp;
						er.y = dr.y/dramp;

						Voro2D_point uij,ua,ub;


						uij.x = (jbp_vx - ibp_vx); 
						uij.y = (jbp_vy - ibp_vy);

						// artificial bulk viscosity
						{
							postype rvel = Vec2DDotP(&er, &uij);
							if(rvel<0){
								postype wcomp = 0.5*(sqrt(ibp->w2hydro)+sqrt(jbp->w2hydro));
								postype scaleFactor = (wcomp ==0 ? etavis: wcomp);
								postype drampScale = dramp/wcomp;
								postype mu = rvel/(drampScale + epsvis/drampScale);
								postype meanden = 0.5*(ibp_den + jbp_den);
								postype meanCsound = 0.5*(ibp_csound + jbp_csound);

								/*
								if(center.indx == 58443){
									DEBUGPRINT("p(%ld,%ld): has viscos: wcomp= %g drampScale= %g mu= %g meanden= %g meanCsound= %g :: pi= %g  dpi= %g\n", 
											PINDX(ibp)%256, PINDX(ibp)/256, wcomp, drampScale, mu, meanden, meanCsound, 
											pi, (-alphavis * meanCsound*mu + betavis*mu*mu)*meanden);
								}
								*/

								pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
							}
						}

						// for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);

						die += -pi * Vec2DDotP(&uradix_ui,&dS);

						// for the total energy
						ua.x = ibp_vx + uradix_ui.x;
						ua.y = ibp_vy + uradix_ui.y;
						dte += -pi*Vec2DDotP(&ua, &dS);

						// for the kinetic energy 
						ub.x = ibp_vx;
						ub.y = ibp_vy;
						dke += -pi*Vec2DDotP(&ub, &dS);

						// for the force 
						fx += -pi * dS.x;
						fy += -pi * dS.y;

						/*
						if(center.indx == 58443)
						{
							DEBUGPRINT("p(%ld,%ld): xy= %g %g has ax/y= %g %g pij= %g dS= %g %g pi= %g : p%lu rj= %g %g pj= %g :: jindx= (%ld,%ld)\n", 
									PINDX(ibp)%256, PINDX(ibp)/256,
									ibp->x,ibp->y,
									fx/ibp->mass,fy/ibp->mass,
									pi,dS.x,dS.y, ibp_pressure, PINDX(jbp),jbp->x, jbp->y, jbp_pressure, PINDX(jbp)%256, PINDX(jbp)/256);
						}
						*/

						Voro2D_point dv;
						dv.x = (jbp_vx - ibp_vx); 
						dv.y = (jbp_vy - ibp_vy); 

						postype VdotR = Vec2DDotP(&dv,&er);
						postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
						postype dramp_cfl = fmax(dramp, 0.25*sqrt(ibp->volume));
						postype dt = 2*Courant*dramp_cfl/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d p(%ld,%ld) x/y= %g %g : %g %g : %g %g : %g : dv.xy= %g %g for p(%ld,%ld)\n",
                                    MYID(simpar), i, PINDX(p+i)%nx, PINDX(p+i)/nx,
									p[i].x, p[i].y, dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y,
									PINDX(jbp)%nx, PINDX(jbp)/nx);
//							DEBUGGING(simpar);
							dt = 1e-10;
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

                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;

				/*
				if(center.indx == 58443) 
				{
					DEBUGPRINT("+p%ld has ax/y= %g %g xy= %g %g : vxy= %g %g ::: rho= %g\n", PINDX(ibp),ibp->ax,ibp->ay, ibp->x, ibp->y,
							ibp->vx, ibp->vy, ibp->den);
					exit(0);
				}
				*/

                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g with p= %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin,
							ibp->pressure);
					exit(0);
                }

            }
            free(p);free(neighbors); free(neighwork);

        }
        free(vorocorner);
    }
	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
	{
		postype TDtime;
        MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
        Dtime = TDtime;
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
	int nx = NX(simpar);


	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	postype cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar);
    int mx, my;
    BASICCELL_MX(simpar) = mx = ceil((xmax-xmin)/cellsize);
    BASICCELL_MY(simpar) = my = ceil((ymax-ymin)/cellsize);
	CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
	mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);


	postype dtold = GAS_dtold(simpar);



	float alphavis = GAS_AlphaVis(simpar);
	float betavis = GAS_BetaVis(simpar);
	float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
	float epsvis = GAS_EPSVIS(simpar);


//	DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);

    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=1;iy<my-1;iy++){
    	int mp=1000;
    	Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
    	postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighboringBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
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

						Voro2D_point uij,ua,ub;
						uij.x = (jbp_vx - ibp_vx); 
						uij.y = (jbp_vy - ibp_vy);

						{
							postype rvel = Vec2DDotP(&er, &uij);
							if(rvel<0){
								postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
								postype scaleFactor = (wcomp ==0 ? etavis: wcomp);
								postype drampScale = dramp/scaleFactor;
								postype mu = rvel/(drampScale + epsvis/drampScale);
								postype meanden = 0.5*(ibp_den + jbp_den);
								postype meanCsound = 0.5*(ibp_csound + jbp_csound);
								pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
							}
						}

						// for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

						// for the total energy
						ua.x = ibp_vx + uradix_ui.x;
						ua.y = ibp_vy + uradix_ui.y;
						dte += -pi*Vec2DDotP(&ua, &dS);

						// for the kinetic energy
						ub.x = ibp_vx;
						ub.y = ibp_vy;
						dke += -pi*Vec2DDotP(&ub, &dS);

						// for the force
						fx += -pi * dS.x;
						fy += -pi * dS.y;

						/*
						if(center.indx == 130816)
						{
							DEBUGPRINT("p(%ld,%ld): xy= %g %g has ax/y= %g %g pij= %g dS= %g %g pi= %g : p%lu rj= %g %g pj= %g :: jindx= (%ld,%ld)\n", 
									PINDX(ibp)/256, PINDX(ibp)%256,
									ibp->x,ibp->y,
									fx/ibp->mass,fy/ibp->mass,
									pi,dS.x,dS.y, ibp_pressure, PINDX(jbp),jbp->x, jbp->y, jbp_pressure, PINDX(jbp)/256, PINDX(jbp)%256);
						}
						*/

						Voro2D_point dv;
						dv.x = (jbp_vx - ibp_vx); 
						dv.y = (jbp_vy - ibp_vy); 

						postype VdotR = Vec2DDotP(&dv,&er);
						postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
						postype dramp_cfl = fmax(dramp, 0.25*sqrt(ibp->volume));
						postype dt = 2*Courant*dramp_cfl/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d %ld x/y= %g %g : %g %g : %g %g : %g : dv.xy= %g %g\n",
                                    MYID(simpar), i, PINDX(p+i), p[i].x, p[i].y, dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y);
                            dt = 1e-10;
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
                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;

				/*
				if(center.indx == 130816) 
				{
					DEBUGPRINT("+p%ld has ax/y= %g %g xy= %g %g : vxy= %g %g ::: rho= %g\n", PINDX(ibp),ibp->ax,ibp->ay, ibp->x, ibp->y,
							ibp->vx, ibp->vy, ibp->den);
							
				}
				*/

                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g with p= %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin,
							ibp->pressure);
					exit(0);
                }

            }
            free(p);free(neighbors); free(neighwork);

        }
        free(vorocorner);
    }
	my_free(VORO_BASICCELL(simpar));
	my_free(VORORK4_TBPP(simpar));
	{
		postype TDtime;
        MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN,  MPI_COMM(simpar));
        Dtime = TDtime;
	}
	return Dtime;
}

/* ================================================================
 *  getAccVoro2D_kNN:
 *  k-NN based force computation. Reuses pre-built tree.
 *  No cell grid needed.
 * ================================================================ */
double getAccVoro2D_kNN(SimParameters *simpar,
		postype xmin, postype ymin, postype xmax, postype ymax,
		postype OrderOfAccuracy, postype Courant, postype Gamma,
		TStruct *tree, TPtlStruct *ptl, int nptl,
		int (*skipBP)(treevorork4particletype *)){
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype Dtime = 1.e10;

	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	postype dtold = GAS_dtold(simpar);

	float alphavis = GAS_AlphaVis(simpar);
	float betavis = GAS_BetaVis(simpar);
	float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
	float epsvis = GAS_EPSVIS(simpar);

	int K_START = 32;
	size_t p_size_a = TVORORK4_DDINFO(simpar)[0].n_size;
	char *bp_raw_a = (char*)VORORK4_TBP(simpar);

#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
	for(int ii=0;ii<nbp;ii++){
		int skip_this = 0;
		if(skipBP){
			treevorork4particletype *bpi_chk = (treevorork4particletype*)(bp_raw_a + ii*p_size_a);
			if(skipBP(bpi_chk)) skip_this = 1;
		}
		if(!skip_this){

		int mp = 1000;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		Voro2D_point neighbors[MAX_NUM_NEAR];
		Voro2D_point neighwork[MAX_NUM_NEAR];
		memset(neighwork, 0, sizeof(neighwork));
		PosType maxr;
		int K = K_START;
		int nfound, ip;

		/* Adaptive kNN with retry (no goto/continue — Intel OpenMP safe) */
		int knn_ok = 0;
		while(!knn_ok){
			nfound = searchKNN2D(ptl+ii, K, tree, neighbors, &maxr);
			if(nfound < 3 || maxr < 1e-20) break;

			treevorork4particletype *ibp_tmp = (treevorork4particletype*)(ptl[ii].bp);
			Voro2D_point center;
			center.x = ibp_tmp->x;
			center.y = ibp_tmp->y;
			center.indx = PINDX(ibp_tmp);
			center.csound = ibp_tmp->csound;
			center.w2 = ibp_tmp->w2;

			ibp_tmp->dt = 1.e10;
			ip = Voro2D_FindVC(&center, neighbors, neighwork, nfound, vorocorner, mp, boxsize);

			/* Completeness check */
			Voro2D_Corner *vtmp = vorocorner;
			postype maxvert2 = 0;
			do {
				postype dx = vtmp->x - center.x;
				postype dy = vtmp->y - center.y;
				postype d2 = dx*dx + dy*dy;
				if(d2 > maxvert2) maxvert2 = d2;
				vtmp = vtmp->upperlink;
			} while(vtmp != vorocorner);
			if(sqrt(maxvert2) > 0.95*maxr && K < MAX_NUM_NEAR){
				K = MIN(K*2, MAX_NUM_NEAR);
			} else {
				knn_ok = 1;
			}
		}

		if(knn_ok){
		treevorork4particletype *ibp = (treevorork4particletype*)(ptl[ii].bp);

		postype ibp_vx = ibp->vx;
		postype ibp_vy = ibp->vy;
		postype ibp_csound = ibp->csound;
		postype ibp_pressure = ibp->pressure;
		postype ibp_den = ibp->den;

		Voro2D_Corner *tmp, *tmp2;
		tmp = vorocorner;
		double die, dte, dke, fx, fy;
		die = dte = dke = fx = fy = 0;
		do {
			if(tmp->upperrelated >= 0){
				treevorork4particletype *jbp = (treevorork4particletype*)(neighwork[tmp->upperrelated].bp);
				postype jbp_vx = jbp->vx;
				postype jbp_vy = jbp->vy;
				postype jbp_csound = jbp->csound;
				postype jbp_pressure = jbp->pressure;
				postype jbp_den = jbp->den;

				Voro2D_point line;
				tmp2 = tmp->upperlink;
				line.x = tmp2->x - tmp->x;
				line.y = tmp2->y - tmp->y;
				Voro2D_point dS = voro2D_norm(&line);
				postype facearea = Vec2DLength(tmp, tmp2);
				dS.x = facearea*dS.x;
				dS.y = facearea*dS.y;

				postype pi = VoroRK4_Pressure2D(
						ibp_pressure, jbp_pressure, tmp, neighwork, simpar, RT);

				Voro2D_point dr = EunhaVec2DSub(jbp, ibp);
				postype dramp = sqrt(Vec2DDotP(&dr, &dr));
				/* Skip degenerate face: particles at same position */
				if(dramp < 1e-20){
					tmp = tmp->upperlink;
					continue;
				}
				Voro2D_point er;
				er.x = dr.x/dramp;
				er.y = dr.y/dramp;

				Voro2D_point uij;
				uij.x = (jbp_vx - ibp_vx);
				uij.y = (jbp_vy - ibp_vy);

				/* Skip AV for boundary ghost neighbors — the velocity
				   difference from wall reflection is not a physical shock */
				if(PINDX(jbp) != MAX_INDEX){
					postype rvel = Vec2DDotP(&er, &uij);
					if(rvel < 0){
						postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
						postype scaleFactor = (wcomp == 0 ? etavis : wcomp);
						postype drampScale = dramp/scaleFactor;
						postype mu = rvel/(drampScale + epsvis/drampScale);
						postype meanden = 0.5*(ibp_den + jbp_den);
						postype meanCsound = 0.5*(ibp_csound + jbp_csound);
						pi = pi + (-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
					}
				}

				/* internal energy */
				Voro2D_point uradix_ui = get2dUpqradRk4(ibp, jbp, dtold);
				die += -pi * Vec2DDotP(&uradix_ui, &dS);

				/* total energy */
				Voro2D_point ua;
				ua.x = ibp_vx + uradix_ui.x;
				ua.y = ibp_vy + uradix_ui.y;
				dte += -pi*Vec2DDotP(&ua, &dS);

				/* kinetic energy */
				Voro2D_point ub;
				ub.x = ibp_vx;
				ub.y = ibp_vy;
				dke += -pi*Vec2DDotP(&ub, &dS);

				/* force */
				fx += -pi * dS.x;
				fy += -pi * dS.y;

				/* CFL timestep */
				Voro2D_point dv;
				dv.x = (jbp_vx - ibp_vx);
				dv.y = (jbp_vy - ibp_vy);
				postype VdotR = Vec2DDotP(&dv, &er);
				postype vsig = (jbp_csound + ibp_csound - MIN(0, VdotR));
				/* CFL floor: prevent dt→0 from particle clustering
				   or geometric proximity of ghost mirrors */
				postype dramp_cfl = fmax(dramp, 0.25*sqrt(ibp->volume));
				postype dt = 2*Courant*dramp_cfl/vsig;
				if(isnan(dt)){
					DEBUGPRINT("P%d kNN has error dt %d %ld : %g %g : %g %g : %g\n",
							MYID(simpar), ii, PINDX(ibp), dramp, vsig, jbp_csound, ibp_csound, VdotR);
					dt = 1e-10;  /* Fallback: dont crash, use tiny dt */
				}
				ibp->dt = MIN(ibp->dt, dt);
				if(dt < Dtime) Dtime = dt;
			}
			tmp = tmp->upperlink;
		} while(tmp != vorocorner);
		ibp->die = die;
		ibp->ax = fx/ibp->mass;
		ibp->ay = fy/ibp->mass;

		if(isnan(fx)){
			static int nan_warn = 0;
			if(nan_warn < 10)
				DEBUGPRINT("P%d kNN has nan %d : %g %g %ld\n",
						MYID(simpar), ii, ibp->x, ibp->y, PINDX(ibp));
			nan_warn++;
			ibp->ax = 0; ibp->ay = 0; ibp->die = 0;
		}
		} /* if(knn_ok) */
		free(vorocorner);
		} /* if(!skip_this) */
	}
	{
		postype TDtime;
		MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
		Dtime = TDtime;
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
	/* Stride-aware: p_size may differ from sizeof(treevorork4particletype)
	   when stress particles are active (av_mode >= 1). */
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
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
    CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
    // building the linked list with "paddingAllTreeParticles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
    mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);
	// repositioning the memory space
	char *bp_raw = (char*)VORORK4_TBP(simpar);
    int i;
    int iy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(iy=1;iy<my-1;iy++){
		int mp=2000;
		Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
		postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
		for(ix=1;ix<mx-1;ix++){
			int np;
			treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
			int nneigh;
			Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
			Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
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
			free(p);free(neighbors); free(neighwork);
		}
		free(vorocorner);
	}
	my_free(VORORK4_TBPP(simpar));
	my_free(VORO_BASICCELL(simpar));
    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;
	float fshift = GAS_FCENTROID(simpar);
	bp_raw = (char*)VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)(bp_raw + i*p_size);
		bpi->x += fshift*bpi->ax; // move it to the centroid  with a factor of shift
		bpi->y += fshift*bpi->ay;
		bpi->x = fmod(bpi->x+Lx,Lx);
		bpi->y = fmod(bpi->y+Ly,Ly);
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

	postype xmin,ymin,zmin,xmax,ymax,zmax, cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
    ymin = Ymin_HydroExam(simpar)-cellsize;
    xmax = Xmax_HydroExam(simpar)+cellsize;
    ymax = Ymax_HydroExam(simpar)+cellsize;
//  determine the minimum dpq for each particle by updating the w2ceil & w2
#ifdef USE_CUDA
	if(GAS_Kappa(simpar)>0) det2d_dpqRK4_GPU(simpar, xmin, ymin, xmax, ymax,
		paddingAllTreeParticles, mkLinkedList2D);
#else
	if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingAllTreeParticles);
#endif

	int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);
	// prepare the linked-list cells for mkLinkedList2D()  and tree findings below
	CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
	// building the linked list with "paddingAllTreeParticles()" which pads the domain
	// with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
	mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);

	int i;
	int iy;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(iy=0;iy<my;iy++){
        int mp=1000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=0;ix<mx;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
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
//				ibp->volume = Area2DPolygon(vorocorner,mp);
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);
				ibp->den = ibp->mass / ibp->volume;
			}
			free(neighwork);
			free(neighbors); 
			free(p);
		}
		free(vorocorner);
	
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
			mkLinkedList2D, 1);
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
			mkLinkedList2D, Dtime);
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
			mkLinkedList2D, Dtime);
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
			mkLinkedList2D, Dtime);
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

	// split getw2forHydroParticle and pressure update because the former uses the latter.
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(GAS_Kappa(simpar) <0) { 
			bp[i].w2 = -GAS_Kappa(simpar); 
		} 
		else if (GAS_Kappa(simpar) >0){ 
			bp[i].w2hydro = getw2forHydroParticle(simpar,(bp+i),Dtime);
			applyW2Controls(simpar, bp+i, Dtime);
		}
    }
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].den = bp[i].mass/bp[i].volume;
        bp[i].pressure = bp[i].ie/bp[i].volume * (Gamma-1);
        bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
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


void readRk4Data2D(SimParameters *simpar, postype *t, postype *dt, int nstep){
	float dtold;
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
		fread(dt,sizeof(postype), 1,fp);
		fread(&GAS_dtold(simpar),sizeof(float), 1,fp);
		fclose(fp);
	}
	else {
		VORO_NP(simpar) = 0;
		VORORK4_TBP(simpar) = (treevorork4particletype*)my_malloc(sizeof(treevorork4particletype)*100);;

	}
	DEBUGPRINT("P%d: after reading iput data\n", myid);
	migrateTreeVorork4Particles(simpar);
	MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(dt, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Bcast(&GAS_dtold(simpar), 1, MPI_FLOAT,0, MPI_COMM(simpar));
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
    CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
    // building the linked list with "paddingTreeVorork4Particles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingTreeVorork4Particles);

    postype Gamma = GAS_GAMMA(simpar);

    float alphavis = GAS_AlphaVis(simpar);
    float betavis = GAS_BetaVis(simpar);
    float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
    float epsvis = GAS_EPSVIS(simpar);


//  DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);


    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = findCellRk4BP2D(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = searchCellRk4Neighbors2D(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
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
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);
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

                        Voro2D_point uij,ua,ub;
                        uij.x = (jbp_vx - ibp_vx); 
                        uij.y = (jbp_vy - ibp_vy);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &uij);
                        postype mu, meanden, meanCsound;
                        if(rvel<0){
							postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
							postype scaleFactor = (wcomp ==0 ? etavis: wcomp);
							postype drampScale = dramp/scaleFactor;
                            mu = rvel/(drampScale + epsvis/drampScale);
                            meanden = 0.5*(ibp_den + jbp_den);
                            meanCsound = 0.5*(ibp_csound + jbp_csound);
                            pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        /* for the internal energy */
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

                        /* for the total energy */
                        ua.x = ibp_vx + uradix_ui.x;
                        ua.y = ibp_vy + uradix_ui.y;
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
                        postype dramp_cfl = fmax(dramp, 0.25*sqrt(ibp->volume));
                        postype dt = 2*Courant*dramp_cfl/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d %ld : %g %g : %g %g : %g : dv.xy= %g %g\n",
                                    MYID(simpar), i, PINDX(p+i), dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y);
                            dt = 1e-10;
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
				/*
                ibp->dte = dte;
                ibp->dke = dke;
				*/
                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;
                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin);
                }

            }
            free(p);free(neighbors); free(neighwork);

        }
        free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));

//  DEBUGPRINT("P%d is exiting here\n", MYID(simpar)); exit(9);

    postype TDtime;
    MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN, MPI_COMM(simpar));
    Dtime = TDtime;

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
			Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp); 
			postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a; 
			int ix; 
			for(ix=1;ix<mx-1;ix++){ 
				int np; 
				treevorork4particletype *p = findCellRk4BP2D(simpar,ix,iy,&np); 
				int nneigh; 
				Voro2D_point *neighbors = searchCellRk4Neighbors2D(simpar,ix,iy,&nneigh); 
				Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh); 
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
				free(p);free(neighbors); free(neighwork); 
			} 
			free(vorocorner); 
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
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = findCellRk4BP2D(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = searchCellRk4Neighbors2D(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
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

//                ibp->volume = Area2DPolygon(vorocorner, mp);
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);




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

                        Voro2D_point uij,ua,ub;
                        uij.x = (jbp_vx - ibp_vx);
                        uij.y = (jbp_vy - ibp_vy);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &uij);
                        if(rvel<0){
							 postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
							 postype scaleFactor = (wcomp ==0 ? etavis: wcomp);
							 postype drampScale = dramp/scaleFactor;
                             postype mu = rvel/(drampScale + epsvis/drampScale);
                             postype meanden = 0.5*(ibp_den + jbp_den);
                             postype meanCsound = 0.5*(ibp_csound + jbp_csound);
                             pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
//                ibp->ax = fx/ibp->mass;
//                ibp->ay = fy/ibp->mass;
            }
            free(p);free(neighbors); free(neighwork);
        }
        free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));
    my_free(VORO_BASICCELL(simpar));

    // update the w2ceil
    postype dmean = SIMBOX(simpar).x.max/NX(simpar);

//  DEBUGPRINT("P%d is now updating the internal energy\n", MYID(simpar));
	// split getw2forHydroParticle and pressure update because the former uses the latter.
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
        if(GAS_Kappa(simpar) <0) {
            bp[i].w2 = -GAS_Kappa(simpar);
        }
        else if (GAS_Kappa(simpar) >0){
          	bp[i].w2hydro = getw2forHydroParticle(simpar,(bp+i),Dtime);
            applyW2Controls(simpar, bp+i, Dtime);
        }
    }
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].ie += bp[i].die * Dtime;
        bp[i].den  = bp[i].mass/bp[i].volume;
        bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1);
        bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
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
    CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
    // building the linked list with "paddingAllTreeParticles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

    postype Gamma = GAS_GAMMA(simpar);

    float alphavis = GAS_AlphaVis(simpar);
    float betavis = GAS_BetaVis(simpar);
    float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
    float epsvis = GAS_EPSVIS(simpar);


//  DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);


    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
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
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);
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

                        Voro2D_point uij,ua,ub;
                        uij.x = (jbp_vx - ibp_vx); 
                        uij.y = (jbp_vy - ibp_vy);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &uij);
                        postype mu, meanden, meanCsound;
                        if(rvel<0){ 
							postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
							postype scaleFactor = (wcomp ==0 ? etavis: wcomp);
							postype drampScale = dramp/scaleFactor;
                            mu = rvel/(drampScale + epsvis/drampScale);
                            meanden = 0.5*(ibp_den + jbp_den);
                            meanCsound = 0.5*(ibp_csound + jbp_csound);
                            pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

                        // for the total energy
                        ua.x = ibp_vx + uradix_ui.x;
                        ua.y = ibp_vy + uradix_ui.y;
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
                        postype dramp_cfl = fmax(dramp, 0.25*sqrt(ibp->volume));
                        postype dt = 2*Courant*dramp_cfl/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d %ld : %g %g : %g %g : %g : dv.xy= %g %g\n",
                                    MYID(simpar), i, PINDX(p+i), dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y);
                            dt = 1e-10;
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
				/*
                ibp->dte = dte;
                ibp->dke = dke;
				*/
                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;
                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin);
                }

            }
            free(p);free(neighbors); free(neighwork);

        }
        free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));

//  DEBUGPRINT("P%d is exiting here\n", MYID(simpar)); exit(9);

    postype TDtime;
    MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN,  MPI_COMM(simpar));
    Dtime = TDtime;

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
			Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp); 
			postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a; 
			int ix; 
			for(ix=1;ix<mx-1;ix++){ 
				int np; 
				treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np); 
				int nneigh; 
				Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh); 
				Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh); 
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
				free(p);free(neighbors); free(neighwork); 
			} 
			free(vorocorner); 
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
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
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

//                ibp->volume = Area2DPolygon(vorocorner, mp);
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);




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

                        Voro2D_point uij,ua,ub;
                        uij.x = (jbp_vx - ibp_vx);
                        uij.y = (jbp_vy - ibp_vy);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &uij);
                        if(rvel<0){
							 postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
							 postype scaleFactor = (wcomp ==0 ? etavis: wcomp);
							 postype drampScale = dramp/scaleFactor;
                             postype mu = rvel/(drampScale + epsvis/drampScale);
                             postype meanden = 0.5*(ibp_den + jbp_den);
                             postype meanCsound = 0.5*(ibp_csound + jbp_csound);
                             pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
//                ibp->ax = fx/ibp->mass;
//                ibp->ay = fy/ibp->mass;
            }
            free(p);free(neighbors); free(neighwork);
        }
        free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));
    my_free(VORO_BASICCELL(simpar));

    // update the w2ceil
    postype dmean = SIMBOX(simpar).x.max/NX(simpar);
//  DEBUGPRINT("P%d is now updating the internal energy\n", MYID(simpar));
	// split getw2forHydroParticle and pressure update because the former uses the latter.
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
			if(GAS_Kappa(simpar) <0) { 
				bp[i].w2 = -GAS_Kappa(simpar); 
			} 
			else if (GAS_Kappa(simpar) >0){ 
          		bp[i].w2hydro = getw2forHydroParticle(simpar,(bp+i),Dtime);
				applyW2Controls(simpar, bp+i, Dtime);
			} 
		}
	}
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
	        bp[i].ie += bp[i].die * Dtime; 
			bp[i].den  = bp[i].mass/bp[i].volume; 
			bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1); 
			bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den); 
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
    CellType *cells = (VORO_BASICCELL(simpar)= (CellType*)my_malloc(sizeof(CellType)*mx*my));
    // building the linked list with "paddingAllTreeParticles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
    mkLinkedList2D(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

    postype Gamma = GAS_GAMMA(simpar);

    float alphavis = GAS_AlphaVis(simpar);
    float betavis = GAS_BetaVis(simpar);
    float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
    float epsvis = GAS_EPSVIS(simpar);


//  DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);


    int iy;
#ifdef _OPENMP
#pragma omp parallel for reduction(min:Dtime)
#endif
    for(iy=1;iy<my-1;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
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
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);
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

                        Voro2D_point uij,ua,ub;
                        //
                        uij.x = (jbp_vx - ibp_vx); 
                        uij.y = (jbp_vy - ibp_vy);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &uij);
                        postype mu, meanden, meanCsound;
                        if(rvel<0){ 
							postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
							postype scaleFactor = (wcomp ==0 ? etavis: wcomp);
							postype drampScale = dramp/scaleFactor;
                            mu = rvel/(drampScale + epsvis/drampScale);
                            meanden = 0.5*(ibp_den + jbp_den);
                            meanCsound = 0.5*(ibp_csound + jbp_csound);
                            pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

                        // for the total energy
                        ua.x = ibp_vx + uradix_ui.x;
                        ua.y = ibp_vy + uradix_ui.y;
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
                        postype dramp_cfl = fmax(dramp, 0.25*sqrt(ibp->volume));
                        postype dt = 2*Courant*dramp_cfl/vsig;
                        if(isnan(dt)){
                            DEBUGPRINT("P%d has error dt %d %ld : %g %g : %g %g : %g : dv.xy= %g %g\n",
                                    MYID(simpar), i, PINDX(p+i), dramp, vsig,jbp_csound, ibp_csound, VdotR, dv.x, dv.y);
                            dt = 1e-10;
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
				/*
                ibp->dte = dte;
                ibp->dke = dke;
				*/
                ibp->ax = fx/ibp->mass;
                ibp->ay = fy/ibp->mass;
                if(isnan(fx)){
                    DEBUGPRINT("P%d has nan %d : %d %d : %g %g %ld in xymin= %g %g\n",
                            MYID(simpar), i, ix,iy, ibp->x,ibp->y, PINDX(ibp),xmin, ymin);
                }

            }
            free(p);free(neighbors); free(neighwork);

        }
        free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));

//  DEBUGPRINT("P%d is exiting here\n", MYID(simpar)); exit(9);

    postype TDtime;
    MPI_Allreduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MIN,  MPI_COMM(simpar));
    Dtime = TDtime;

    int i;
    postype dmax=0;
//  postype eps2 = KH_EPS(simpar)*KH_EPS(simpar);
#ifdef _OPENMP
#pragma omp parallel for reduction(max: dmax)
#endif
    for(i=0;i<VORO_NP(simpar);i++){
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
			Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp); 
			postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a; 
			int ix; 
			for(ix=1;ix<mx-1;ix++){ 
				int np; 
				treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np); 
				int nneigh; 
				Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh); 
				Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh); 
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
				free(p);free(neighbors); free(neighwork); 
			} 
			free(vorocorner); 
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
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=1;ix<mx-1;ix++){
            int np;
            treevorork4particletype *p = find2DCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = find2DNeighborBP(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
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

//                ibp->volume = Area2DPolygon(vorocorner, mp);
                get2dAreaAvgNeighorPressure(ibp,vorocorner, neighwork,bp);




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

                        Voro2D_point uij,ua,ub;
                        uij.x = (jbp_vx - ibp_vx);
                        uij.y = (jbp_vy - ibp_vy);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;


                        postype rvel = Vec2DDotP(&er, &uij);
                        if(rvel<0){ 
							postype wcomp = sqrt(ibp->w2)+sqrt(jbp->w2);
							postype scaleFactor = (wcomp ==0 ? etavis: wcomp);
							postype drampScale = dramp/scaleFactor; 
							postype mu = rvel/(drampScale + epsvis/drampScale); 
							postype meanden = 0.5*(ibp_den + jbp_den); 
							postype meanCsound = 0.5*(ibp_csound + jbp_csound); 
							pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                        }

                        // for the internal energy
                        Voro2D_point uradix_ui = get2dUpqradRk4(ibp,jbp,dtold);
						die += -pi * Vec2DDotP(&uradix_ui,&dS);

                    }
                    tmp = tmp->upperlink;
                } while( tmp != vorocorner);
                ibp->die = die;
//                ibp->ax = fx/ibp->mass;
//                ibp->ay = fy/ibp->mass;
            }
            free(p);free(neighbors); free(neighwork);
        }
        free(vorocorner);
    }
    my_free(VORORK4_TBPP(simpar));
    my_free(VORO_BASICCELL(simpar));

    // update the w2ceil
    postype dmean = SIMBOX(simpar).x.max/NX(simpar);
//  DEBUGPRINT("P%d is now updating the internal energy\n", MYID(simpar));
	// split getw2forHydroParticle and pressure update because the former uses the latter.
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
//		if(bp[i].y >=0.1*Ly && bp[i].y < 0.9*Ly)
		if(targetBP(bp+i, Lx, Ly))
		{
			if(GAS_Kappa(simpar) <0) { 
				bp[i].w2 = -GAS_Kappa(simpar); 
			} 
			else if (GAS_Kappa(simpar) >0){ 
          		bp[i].w2hydro = getw2forHydroParticle(simpar,(bp+i),Dtime);
				applyW2Controls(simpar, bp+i, Dtime);
			} 
		}
	}
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
//	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has x/y minmax= %g %g : %g %g ::: %g\n", MYID(simpar), xmin,ymin,xmax,ymax, cellsize);

	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);


    int i;
    postype Dtime,dt;
    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) {
		bp[i].w2old = bp[i].w2;
		bp[i].rk4.w2backup = bp[i].w2;
	}
	// Runge-Kutta 4-th order time evolution of r and vr
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, 1);
//	MPI_Barrier(MPI_COMM(simpar)); if(MYID(simpar)==0) DEBUGPRINT("P%d has passed updateDenW2P...()\n", MYID(simpar));
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
			mkLinkedList2D, Dtime);
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
			mkLinkedList2D, Dtime);
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
			mkLinkedList2D, Dtime);
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
	// split getw2forHydroParticle and pressure update because the former uses the latter.
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
			if(GAS_Kappa(simpar) <0) { 
				bp[i].w2 = -GAS_Kappa(simpar); 
			} 
			else if (GAS_Kappa(simpar) >0){ 
          		bp[i].w2hydro = getw2forHydroParticle(simpar,(bp+i),Dtime);
				applyW2Controls(simpar, bp+i, Dtime);
			} 
		}
	}
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
			bp[i].den = bp[i].mass/bp[i].volume; 
			bp[i].pressure = bp[i].ie/bp[i].volume * (Gamma-1); 
			bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den); 
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

/* ================================================================
 *  periodic_postStage_kNN / periodic_postStage_blend:
 *  Default postStage callbacks for periodic boundary models (KH etc.).
 *  Wrap x,y and migrate particles between MPI ranks.
 * ================================================================ */
void periodic_postStage_kNN(SimParameters *simpar){
	postype Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	postype Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	size_t p_size = TVORORK4_DDINFO(simpar)[0].n_size;
	char *raw = (char*)VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;
	for(i=0;i<np;i++){
		treevorork4particletype *bp = (treevorork4particletype*)(raw + i*p_size);
		bp->x = fmod(bp->x + Lx, Lx);
		bp->y = fmod(bp->y + Ly, Ly);
	}
	migrateTreeVorork4Particles(simpar);
}

void periodic_postStage_blend(SimParameters *simpar){
	postype Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	postype Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i;
	for(i=0;i<np;i++){
		sbp[i].x = fmod(sbp[i].x + Lx, Lx);
		sbp[i].y = fmod(sbp[i].y + Ly, Ly);
	}
	migrateTreeVorork4Particles(simpar);
}

/* ================================================================
 *  exam2d_vph_rk4_int_blend:
 *  RK4 integrator with NS stress + HLLC + CD10 two-tier blending.
 *  Uses treevorostressrk4particletype for stress fields.
 *  Same RK4 structure as exam2d_vph_rk4_int, but calls
 *  updateDenW2Pressure2DBlend and getAccVoro2DBlend.
 * ================================================================ */
double exam2d_vph_rk4_int_blend(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorork4particletype *(*find2DCellBP)(SimParameters *, int , int , int *),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		void (*postStage)(SimParameters *)
		){
	/* All array indexing uses treevorostressrk4particletype* (sbp) for correct
	   stride when particles are stress-sized.  Common fields (x, y, vx, vy,
	   rk4, VoroQ, etc.) are at the same offsets in both struct types. */
	treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	postype xmin,ymin,xmax,ymax,OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);
	int av_mode = GAS_AVMODE(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);

	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);

	int i;
	postype Dtime,dt;
	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) {
		sbp[i].w2old = sbp[i].w2;
		sbp[i].rk4.w2backup = sbp[i].w2;
	}

	/* Reset vsig_max at start */
	if(av_mode >= 1){
		for(i=0;i<VORO_NP(simpar);i++)
			sbp[i].stress.vsig_max = 0;
	}

	/* Revision 4 defensive init: first-call guarantee that CD10 state
	 * fields (divv_old, alpha_cd) are not garbage, regardless of
	 * allocator behaviour.  my_malloc may or may not zero-init; this
	 * removes the dependency. */
	{
		static int first_call = 1;
		if(first_call && av_mode >= 1){
			postype cd_amin_init = GAS_CDAMIN(simpar);
			for(i=0;i<VORO_NP(simpar);i++){
				sbp[i].stress.divv_old = 0;
				sbp[i].stress.alpha_cd = cd_amin_init;
			}
			first_call = 0;
		}
	}

	double _t0, _t_update=0, _t_force=0, _t_post=0, _t_fin=0;

	/* === K1 evaluation === */
	_t0 = MPI_Wtime();
	updateDenW2Pressure2DBlend(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, 1);
	_t_update += MPI_Wtime() - _t0;
	_t0 = MPI_Wtime();
#if defined(USE_CUDA) && defined(GPU_VALIDATE)
	Dtime = getAccVoro2DBlend_GPU_validate(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
#elif defined(USE_CUDA)
	Dtime = getAccVoro2DBlend_GPU(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
#else
	Dtime = getAccVoro2DBlend(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
#endif
	_t_force += MPI_Wtime() - _t0;

	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].rk4.k1x  = sbp[i].vx*Dtime;
		sbp[i].rk4.k1y  = sbp[i].vy*Dtime;
		sbp[i].rk4.k1vx = (sbp[i].ax + GAS_ACCX(simpar))*Dtime;
		sbp[i].rk4.k1vy = (sbp[i].ay + GAS_ACCY(simpar))*Dtime;
		sbp[i].rk4.k1ie = sbp[i].die*Dtime;
	}

	/* === K2 evaluation === */
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].x  += sbp[i].rk4.k1x*0.5;
		sbp[i].y  += sbp[i].rk4.k1y*0.5;
		sbp[i].vx += sbp[i].rk4.k1vx*0.5;
		sbp[i].vy += sbp[i].rk4.k1vy*0.5;
		sbp[i].ie += sbp[i].rk4.k1ie*0.5;
	}
	_t0 = MPI_Wtime();
	postStage(simpar);
	_t_post += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;
	if(av_mode >= 1)
		for(i=0;i<VORO_NP(simpar);i++) sbp[i].stress.vsig_max = 0;
	_t0 = MPI_Wtime();
	updateDenW2Pressure2DBlend(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
	_t_update += MPI_Wtime() - _t0;
	_t0 = MPI_Wtime();
#if defined(USE_CUDA) && defined(GPU_VALIDATE)
	dt = getAccVoro2DBlend_GPU_validate(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
#elif defined(USE_CUDA)
	dt = getAccVoro2DBlend_GPU(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
#else
	dt = getAccVoro2DBlend(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
#endif
	_t_force += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].rk4.k2x  = sbp[i].vx*Dtime;
		sbp[i].rk4.k2y  = sbp[i].vy*Dtime;
		sbp[i].rk4.k2vx = (sbp[i].ax + GAS_ACCX(simpar))*Dtime;
		sbp[i].rk4.k2vy = (sbp[i].ay + GAS_ACCY(simpar))*Dtime;
		sbp[i].rk4.k2ie = sbp[i].die*Dtime;
	}

	/* === K3 evaluation === */
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].x  += (sbp[i].rk4.k2x -sbp[i].rk4.k1x )*0.5;
		sbp[i].y  += (sbp[i].rk4.k2y -sbp[i].rk4.k1y )*0.5;
		sbp[i].vx += (sbp[i].rk4.k2vx-sbp[i].rk4.k1vx)*0.5;
		sbp[i].vy += (sbp[i].rk4.k2vy-sbp[i].rk4.k1vy)*0.5;
		sbp[i].ie += (sbp[i].rk4.k2ie-sbp[i].rk4.k1ie)*0.5;
	}
	_t0 = MPI_Wtime();
	postStage(simpar);
	_t_post += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;
	if(av_mode >= 1)
		for(i=0;i<VORO_NP(simpar);i++) sbp[i].stress.vsig_max = 0;
	_t0 = MPI_Wtime();
	updateDenW2Pressure2DBlend(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
	_t_update += MPI_Wtime() - _t0;
	_t0 = MPI_Wtime();
#if defined(USE_CUDA) && defined(GPU_VALIDATE)
	dt = getAccVoro2DBlend_GPU_validate(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
#elif defined(USE_CUDA)
	dt = getAccVoro2DBlend_GPU(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
#else
	dt = getAccVoro2DBlend(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
#endif
	_t_force += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].rk4.k3x  = sbp[i].vx*Dtime;
		sbp[i].rk4.k3y  = sbp[i].vy*Dtime;
		sbp[i].rk4.k3vx = (sbp[i].ax + GAS_ACCX(simpar))*Dtime;
		sbp[i].rk4.k3vy = (sbp[i].ay + GAS_ACCY(simpar))*Dtime;
		sbp[i].rk4.k3ie = sbp[i].die*Dtime;
	}

	/* === K4 evaluation === */
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].x  += sbp[i].rk4.k3x  - sbp[i].rk4.k2x*0.5;
		sbp[i].y  += sbp[i].rk4.k3y  - sbp[i].rk4.k2y*0.5;
		sbp[i].vx += sbp[i].rk4.k3vx - sbp[i].rk4.k2vx*0.5;
		sbp[i].vy += sbp[i].rk4.k3vy - sbp[i].rk4.k2vy*0.5;
		sbp[i].ie += sbp[i].rk4.k3ie - sbp[i].rk4.k2ie*0.5;
	}
	_t0 = MPI_Wtime();
	postStage(simpar);
	_t_post += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;
	if(av_mode >= 1)
		for(i=0;i<VORO_NP(simpar);i++) sbp[i].stress.vsig_max = 0;
	_t0 = MPI_Wtime();
	updateDenW2Pressure2DBlend(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
	_t_update += MPI_Wtime() - _t0;
	_t0 = MPI_Wtime();
#if defined(USE_CUDA) && defined(GPU_VALIDATE)
	dt = getAccVoro2DBlend_GPU_validate(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
#elif defined(USE_CUDA)
	dt = getAccVoro2DBlend_GPU(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
#else
	dt = getAccVoro2DBlend(simpar, xmin, ymin, xmax, ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
			find2DNeighborBP, find2DCellBP, mkLinkedList2D);
#endif
	_t_force += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].rk4.k4x  = sbp[i].vx*Dtime;
		sbp[i].rk4.k4y  = sbp[i].vy*Dtime;
		sbp[i].rk4.k4vx = (sbp[i].ax + GAS_ACCX(simpar))*Dtime;
		sbp[i].rk4.k4vy = (sbp[i].ay + GAS_ACCY(simpar))*Dtime;
		sbp[i].rk4.k4ie = sbp[i].die*Dtime;
	}
	_t0 = MPI_Wtime();

	/* === Undo K4 shift and prepare for final combination === */
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].x  -= sbp[i].rk4.k3x;
		sbp[i].y  -= sbp[i].rk4.k3y;
		sbp[i].vx -= sbp[i].rk4.k3vx;
		sbp[i].vy -= sbp[i].rk4.k3vy;
		sbp[i].ie -= sbp[i].rk4.k3ie;
	}
	_t0 = MPI_Wtime();
	postStage(simpar);
	_t_post += MPI_Wtime() - _t0;
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;

	/* === Final RK4 combination === */
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
			sbp[i].x  += (sbp[i].rk4.k1x +2*sbp[i].rk4.k2x +2*sbp[i].rk4.k3x +sbp[i].rk4.k4x )/6.;
			sbp[i].y  += (sbp[i].rk4.k1y +2*sbp[i].rk4.k2y +2*sbp[i].rk4.k3y +sbp[i].rk4.k4y )/6.;
			sbp[i].vx += (sbp[i].rk4.k1vx+2*sbp[i].rk4.k2vx+2*sbp[i].rk4.k3vx+sbp[i].rk4.k4vx)/6.;
			sbp[i].vy += (sbp[i].rk4.k1vy+2*sbp[i].rk4.k2vy+2*sbp[i].rk4.k3vy+sbp[i].rk4.k4vy)/6.;
			sbp[i].ie += (sbp[i].rk4.k1ie+2*sbp[i].rk4.k2ie+2*sbp[i].rk4.k3ie+sbp[i].rk4.k4ie)/6.;
		}
	}
	postStage(simpar);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	migrateTreeVorork4Particles(simpar);

	/* Arepo-Laguerre (av_mode=5): optional centroidal steering for
	   mesh regularization.  Nudges each generator toward its Voronoi
	   cell centroid by fraction GAS_FCENTROID, preventing cell
	   degeneracy in strong shear flows (KH).  Safe no-op when
	   GAS_FCENTROID == 0. */
	if(av_mode == 5 && GAS_FCENTROID(simpar) > 0){
		exam2d_centroidShift(simpar, paddingAllTreeParticles,
				find2DNeighborBP, find2DCellBP, mkLinkedList2D);
	}

	/* Update volume & density */
	exam2dUpdateVol(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,mkLinkedList2D);

	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	/* Update w2 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
			if(GAS_Kappa(simpar) <0){
				sbp[i].w2 = -GAS_Kappa(simpar);
			}
			else if (GAS_Kappa(simpar) >0){
				sbp[i].w2hydro = getw2forHydroParticle(simpar,(treevorork4particletype*)(sbp+i),Dtime);
				applyW2Controls(simpar, (treevorork4particletype*)(sbp+i), Dtime);
			}
		}
	}

	/* Update pressure, csound, NaN check */
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
			sbp[i].den = sbp[i].mass/sbp[i].volume;
			sbp[i].pressure = sbp[i].ie/sbp[i].volume * (Gamma-1);
			if(sbp[i].pressure <= 0){
				sbp[i].pressure = 1e-6;
				sbp[i].ie = sbp[i].pressure * sbp[i].volume / (Gamma-1);
			}
			sbp[i].csound = sqrt(Gamma*sbp[i].pressure/sbp[i].den);
		}
		else {
			sbp[i].vx = sbp[i].vy = 0;
			sbp[i].w2 = sbp[i].rk4.w2backup;
		}
	}

	/* CD10 alpha update after full RK4 step */
	if(av_mode >= 1){
		update_alpha_cd_2d(simpar, Dtime);
	}

	_t_fin = MPI_Wtime() - _t0;
	{
		static int _step_timer = 0;
		_step_timer++;
		double _t_total = _t_update + _t_force + _t_post + _t_fin;
		if(MYID(simpar) == 0)
			printf("[CPU step %d] update=%.1fs force=%.1fs post=%.1fs final=%.1fs | total=%.1fs\n",
				_step_timer, _t_update, _t_force, _t_post, _t_fin, _t_total);
		fflush(stdout);
	}

	return Dtime;
}

/* ================================================================
 *  exam2d_vph_rk4_int_lagmfm (av_mode=4):
 *    RK4 integrator calling the kernel-based MFM density/force path.
 *    Mirrors exam2d_vph_rk4_int_blend but uses
 *      updateDenW2Pressure2D_LagMFM   (kernel density + gradients)
 *      getAccVoro2D_LagMFM             (MFM effective-face force)
 * ================================================================ */
double exam2d_vph_rk4_int_lagmfm(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		void mkLinkedList2D(SimParameters *, postype, postype , postype , postype , postype,
			void (*)(SimParameters *, postype)),
		void (*postStage)(SimParameters *)
		){
	treevorostressrk4particletype *sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	postype xmin,ymin,xmax,ymax,OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);
	int av_mode = GAS_AVMODE(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);

	int i;
	postype Dtime, dt;
	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) {
		sbp[i].w2old = sbp[i].w2;
		sbp[i].rk4.w2backup = sbp[i].w2;
	}
	for(i=0;i<VORO_NP(simpar);i++)
		sbp[i].stress.vsig_max = 0;

	/* === K1 === */
#ifdef USE_CUDA
	updateDenW2Pressure2D_LagMFM_GPU(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, mkLinkedList2D, 1);
	lagmfm_w2_post_density(simpar, 1);
	Dtime = getAccVoro2D_LagMFM_GPU(simpar, xmin,ymin,xmax,ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, mkLinkedList2D);
#else
	updateDenW2Pressure2D_LagMFM(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, mkLinkedList2D, 1);
	Dtime = getAccVoro2D_LagMFM(simpar, xmin,ymin,xmax,ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, mkLinkedList2D);
#endif
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].rk4.k1x  = sbp[i].vx*Dtime;
		sbp[i].rk4.k1y  = sbp[i].vy*Dtime;
		sbp[i].rk4.k1vx = (sbp[i].ax + GAS_ACCX(simpar))*Dtime;
		sbp[i].rk4.k1vy = (sbp[i].ay + GAS_ACCY(simpar))*Dtime;
		sbp[i].rk4.k1ie = sbp[i].die*Dtime;
	}

	/* === K2 === */
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].x  += sbp[i].rk4.k1x*0.5;
		sbp[i].y  += sbp[i].rk4.k1y*0.5;
		sbp[i].vx += sbp[i].rk4.k1vx*0.5;
		sbp[i].vy += sbp[i].rk4.k1vy*0.5;
		sbp[i].ie += sbp[i].rk4.k1ie*0.5;
	}
	postStage(simpar);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;
	for(i=0;i<VORO_NP(simpar);i++) sbp[i].stress.vsig_max = 0;
#ifdef USE_CUDA
	updateDenW2Pressure2D_LagMFM_GPU(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, mkLinkedList2D, Dtime);
	lagmfm_w2_post_density(simpar, Dtime);
	dt = getAccVoro2D_LagMFM_GPU(simpar, xmin,ymin,xmax,ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, mkLinkedList2D);
#else
	updateDenW2Pressure2D_LagMFM(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, mkLinkedList2D, Dtime);
	dt = getAccVoro2D_LagMFM(simpar, xmin,ymin,xmax,ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, mkLinkedList2D);
#endif
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].rk4.k2x  = sbp[i].vx*Dtime;
		sbp[i].rk4.k2y  = sbp[i].vy*Dtime;
		sbp[i].rk4.k2vx = (sbp[i].ax + GAS_ACCX(simpar))*Dtime;
		sbp[i].rk4.k2vy = (sbp[i].ay + GAS_ACCY(simpar))*Dtime;
		sbp[i].rk4.k2ie = sbp[i].die*Dtime;
	}

	/* === K3 === */
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].x  += (sbp[i].rk4.k2x -sbp[i].rk4.k1x )*0.5;
		sbp[i].y  += (sbp[i].rk4.k2y -sbp[i].rk4.k1y )*0.5;
		sbp[i].vx += (sbp[i].rk4.k2vx-sbp[i].rk4.k1vx)*0.5;
		sbp[i].vy += (sbp[i].rk4.k2vy-sbp[i].rk4.k1vy)*0.5;
		sbp[i].ie += (sbp[i].rk4.k2ie-sbp[i].rk4.k1ie)*0.5;
	}
	postStage(simpar);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;
	for(i=0;i<VORO_NP(simpar);i++) sbp[i].stress.vsig_max = 0;
#ifdef USE_CUDA
	updateDenW2Pressure2D_LagMFM_GPU(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, mkLinkedList2D, Dtime);
	lagmfm_w2_post_density(simpar, Dtime);
	dt = getAccVoro2D_LagMFM_GPU(simpar, xmin,ymin,xmax,ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, mkLinkedList2D);
#else
	updateDenW2Pressure2D_LagMFM(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, mkLinkedList2D, Dtime);
	dt = getAccVoro2D_LagMFM(simpar, xmin,ymin,xmax,ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, mkLinkedList2D);
#endif
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].rk4.k3x  = sbp[i].vx*Dtime;
		sbp[i].rk4.k3y  = sbp[i].vy*Dtime;
		sbp[i].rk4.k3vx = (sbp[i].ax + GAS_ACCX(simpar))*Dtime;
		sbp[i].rk4.k3vy = (sbp[i].ay + GAS_ACCY(simpar))*Dtime;
		sbp[i].rk4.k3ie = sbp[i].die*Dtime;
	}

	/* === K4 === */
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].x  += sbp[i].rk4.k3x  - sbp[i].rk4.k2x*0.5;
		sbp[i].y  += sbp[i].rk4.k3y  - sbp[i].rk4.k2y*0.5;
		sbp[i].vx += sbp[i].rk4.k3vx - sbp[i].rk4.k2vx*0.5;
		sbp[i].vy += sbp[i].rk4.k3vy - sbp[i].rk4.k2vy*0.5;
		sbp[i].ie += sbp[i].rk4.k3ie - sbp[i].rk4.k2ie*0.5;
	}
	postStage(simpar);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;
	for(i=0;i<VORO_NP(simpar);i++) sbp[i].stress.vsig_max = 0;
#ifdef USE_CUDA
	updateDenW2Pressure2D_LagMFM_GPU(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, mkLinkedList2D, Dtime);
	lagmfm_w2_post_density(simpar, Dtime);
	dt = getAccVoro2D_LagMFM_GPU(simpar, xmin,ymin,xmax,ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, mkLinkedList2D);
#else
	updateDenW2Pressure2D_LagMFM(simpar, xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles, mkLinkedList2D, Dtime);
	dt = getAccVoro2D_LagMFM(simpar, xmin,ymin,xmax,ymax,
			OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles, mkLinkedList2D);
#endif
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].rk4.k4x  = sbp[i].vx*Dtime;
		sbp[i].rk4.k4y  = sbp[i].vy*Dtime;
		sbp[i].rk4.k4vx = (sbp[i].ax + GAS_ACCX(simpar))*Dtime;
		sbp[i].rk4.k4vy = (sbp[i].ay + GAS_ACCY(simpar))*Dtime;
		sbp[i].rk4.k4ie = sbp[i].die*Dtime;
	}

	/* === Undo K4 shift, then final combination === */
	for(i=0;i<VORO_NP(simpar);i++){
		sbp[i].x  -= sbp[i].rk4.k3x;
		sbp[i].y  -= sbp[i].rk4.k3y;
		sbp[i].vx -= sbp[i].rk4.k3vx;
		sbp[i].vy -= sbp[i].rk4.k3vy;
		sbp[i].ie -= sbp[i].rk4.k3ie;
	}
	postStage(simpar);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++) sbp[i].w2 = sbp[i].rk4.w2backup;

	for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
			sbp[i].x  += (sbp[i].rk4.k1x +2*sbp[i].rk4.k2x +2*sbp[i].rk4.k3x +sbp[i].rk4.k4x )/6.;
			sbp[i].y  += (sbp[i].rk4.k1y +2*sbp[i].rk4.k2y +2*sbp[i].rk4.k3y +sbp[i].rk4.k4y )/6.;
			sbp[i].vx += (sbp[i].rk4.k1vx+2*sbp[i].rk4.k2vx+2*sbp[i].rk4.k3vx+sbp[i].rk4.k4vx)/6.;
			sbp[i].vy += (sbp[i].rk4.k1vy+2*sbp[i].rk4.k2vy+2*sbp[i].rk4.k3vy+sbp[i].rk4.k4vy)/6.;
			sbp[i].ie += (sbp[i].rk4.k1ie+2*sbp[i].rk4.k2ie+2*sbp[i].rk4.k3ie+sbp[i].rk4.k4ie)/6.;
		}
	}
	postStage(simpar);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	migrateTreeVorork4Particles(simpar);
	sbp = (treevorostressrk4particletype*)VORORK4_TBP(simpar);

	/* Refresh w2 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
			if(GAS_Kappa(simpar) < 0){
				sbp[i].w2 = -GAS_Kappa(simpar);
			}
			else if(GAS_Kappa(simpar) > 0){
				sbp[i].w2hydro = getw2forHydroParticle(simpar,(treevorork4particletype*)(sbp+i),Dtime);
				applyW2Controls(simpar, (treevorork4particletype*)(sbp+i), Dtime);
			}
		}
	}

	/* Pressure/csound refresh (volume already updated by last density pass) */
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP((treevorork4particletype*)(sbp+i), Lx, Ly)){
			if(sbp[i].volume <= 0) sbp[i].volume = cellsize*cellsize;
			sbp[i].den = sbp[i].mass/sbp[i].volume;
			sbp[i].pressure = sbp[i].ie/sbp[i].volume * (Gamma-1);
			if(sbp[i].pressure <= 0){
				sbp[i].pressure = 1e-6;
				sbp[i].ie = sbp[i].pressure * sbp[i].volume / (Gamma-1);
			}
			sbp[i].csound = sqrt(Gamma*sbp[i].pressure/sbp[i].den);
		}
		else {
			sbp[i].vx = sbp[i].vy = 0;
			sbp[i].w2 = sbp[i].rk4.w2backup;
		}
	}

	/* CD10 alpha update */
	update_alpha_cd_2d(simpar, Dtime);

	return Dtime;
}

/* ================================================================
 *  exam2d_vph_rk4_int_kNN:
 *  RK4 integrator using k-NN tree search instead of cell-based search.
 *  No cell callbacks needed — uses buildTree2D + kNN functions.
 * ================================================================ */
double exam2d_vph_rk4_int_kNN(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		int (*targetBP)(treevorork4particletype*, postype, postype),
		int (*skipBP)(treevorork4particletype *),
		void (*postStage)(SimParameters *)
		){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	size_t p_sz = TVORORK4_DDINFO(simpar)[0].n_size;
	postype xmin, ymin, xmax, ymax, OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);

	int i;
	postype Dtime, dt;
	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;

	/* Debug: check ie at very start of RK4 */
	DEBUGPRINT("P%d RK4 start: np=%d bp=%p bp[0] ie=%g vol=%g pos=(%g,%g) PINDX=%ld\n",
		MYID(simpar), VORO_NP(simpar), (void*)bp, bp->ie, bp->volume, bp->x, bp->y, PINDX(bp));

	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->w2old = bpi->w2;
		bpi->rk4.w2backup = bpi->w2;
	}

	/* ---- RK4 Stage 1 ---- */
	{
		TStruct *tree; TPtlStruct *ptl; int nptl;
		buildTree2D(simpar, paddingAllTreeParticles, skipBP, &tree, &ptl, &nptl);
		bp = VORORK4_TBP(simpar);
		updateDenW2Pressure2D_kNN(simpar, xmin, ymin, xmax, ymax, Gamma,
				tree, ptl, nptl, 1, skipBP);
		Dtime = getAccVoro2D_kNN(simpar, xmin, ymin, xmax, ymax,
				OrderOfAccuracy, Courant, Gamma, tree, ptl, nptl, skipBP);
		my_free(ptl); my_free(tree); my_free(VORORK4_TBPP(simpar));
	}
	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->rk4.k1x = bpi->vx*Dtime;
		bpi->rk4.k1y = bpi->vy*Dtime;
		bpi->rk4.k1vx = (bpi->ax + GAS_ACCX(simpar))*Dtime;
		bpi->rk4.k1vy = (bpi->ay + GAS_ACCY(simpar))*Dtime;
		bpi->rk4.k1ie = bpi->die*Dtime;
	}
	/* RK4 second half-step */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->x += bpi->rk4.k1x*0.5;
		bpi->y += bpi->rk4.k1y*0.5;
		bpi->vx += bpi->rk4.k1vx*0.5;
		bpi->vy += bpi->rk4.k1vy*0.5;
		bpi->ie += bpi->rk4.k1ie*0.5;
	}
	postStage(simpar);
	bp = VORORK4_TBP(simpar);
	/* ---- RK4 Stage 2 ---- */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->w2 = bpi->rk4.w2backup;
	}
	{
		TStruct *tree; TPtlStruct *ptl; int nptl;
		buildTree2D(simpar, paddingAllTreeParticles, skipBP, &tree, &ptl, &nptl);
		updateDenW2Pressure2D_kNN(simpar, xmin, ymin, xmax, ymax, Gamma,
				tree, ptl, nptl, Dtime, skipBP);
		dt = getAccVoro2D_kNN(simpar, xmin, ymin, xmax, ymax,
				OrderOfAccuracy, Courant, Gamma, tree, ptl, nptl, skipBP);
		my_free(ptl); my_free(tree); my_free(VORORK4_TBPP(simpar));
	}

	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->rk4.k2x = bpi->vx*Dtime;
		bpi->rk4.k2y = bpi->vy*Dtime;
		bpi->rk4.k2vx = (bpi->ax + GAS_ACCX(simpar))*Dtime;
		bpi->rk4.k2vy = (bpi->ay + GAS_ACCY(simpar))*Dtime;
		bpi->rk4.k2ie = bpi->die*Dtime;
	}
	/* RK4 third half-step */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->x += (bpi->rk4.k2x-bpi->rk4.k1x)*0.5;
		bpi->y += (bpi->rk4.k2y-bpi->rk4.k1y)*0.5;
		bpi->vx += (bpi->rk4.k2vx-bpi->rk4.k1vx)*0.5;
		bpi->vy += (bpi->rk4.k2vy-bpi->rk4.k1vy)*0.5;
		bpi->ie += (bpi->rk4.k2ie-bpi->rk4.k1ie)*0.5;
	}
	postStage(simpar);
	bp = VORORK4_TBP(simpar);

	/* ---- RK4 Stage 3 ---- */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->w2 = bpi->rk4.w2backup;
	}
	{
		TStruct *tree; TPtlStruct *ptl; int nptl;
		buildTree2D(simpar, paddingAllTreeParticles, skipBP, &tree, &ptl, &nptl);
		updateDenW2Pressure2D_kNN(simpar, xmin, ymin, xmax, ymax, Gamma,
				tree, ptl, nptl, Dtime, skipBP);
		dt = getAccVoro2D_kNN(simpar, xmin, ymin, xmax, ymax,
				OrderOfAccuracy, Courant, Gamma, tree, ptl, nptl, skipBP);
		my_free(ptl); my_free(tree); my_free(VORORK4_TBPP(simpar));
	}

	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->rk4.k3x = bpi->vx*Dtime;
		bpi->rk4.k3y = bpi->vy*Dtime;
		bpi->rk4.k3vx = (bpi->ax + GAS_ACCX(simpar))*Dtime;
		bpi->rk4.k3vy = (bpi->ay + GAS_ACCY(simpar))*Dtime;
		bpi->rk4.k3ie = bpi->die*Dtime;
	}
	/* RK4 fourth full-step */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->x += bpi->rk4.k3x - bpi->rk4.k2x*0.5;
		bpi->y += bpi->rk4.k3y - bpi->rk4.k2y*0.5;
		bpi->vx += bpi->rk4.k3vx - bpi->rk4.k2vx*0.5;
		bpi->vy += bpi->rk4.k3vy - bpi->rk4.k2vy*0.5;
		bpi->ie += bpi->rk4.k3ie - bpi->rk4.k2ie*0.5;
	}
	postStage(simpar);
	bp = VORORK4_TBP(simpar);

	/* ---- RK4 Stage 4 ---- */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->w2 = bpi->rk4.w2backup;
	}
	{
		TStruct *tree; TPtlStruct *ptl; int nptl;
		buildTree2D(simpar, paddingAllTreeParticles, skipBP, &tree, &ptl, &nptl);
		updateDenW2Pressure2D_kNN(simpar, xmin, ymin, xmax, ymax, Gamma,
				tree, ptl, nptl, Dtime, skipBP);
		dt = getAccVoro2D_kNN(simpar, xmin, ymin, xmax, ymax,
				OrderOfAccuracy, Courant, Gamma, tree, ptl, nptl, skipBP);
		my_free(ptl); my_free(tree); my_free(VORORK4_TBPP(simpar));
	}

	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->rk4.k4x = bpi->vx*Dtime;
		bpi->rk4.k4y = bpi->vy*Dtime;
		bpi->rk4.k4vx = (bpi->ax + GAS_ACCX(simpar))*Dtime;
		bpi->rk4.k4vy = (bpi->ay + GAS_ACCY(simpar))*Dtime;
		bpi->rk4.k4ie = bpi->die*Dtime;
	}
	/* Rollback k3 step */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->x -= bpi->rk4.k3x;
		bpi->y -= bpi->rk4.k3y;
		bpi->vx -= bpi->rk4.k3vx;
		bpi->vy -= bpi->rk4.k3vy;
		bpi->ie -= bpi->rk4.k3ie;
	}
	postStage(simpar);
	bp = VORORK4_TBP(simpar);

	/* Final RK4 combination */
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		bpi->w2 = bpi->rk4.w2backup;
	}
	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		if(targetBP(bpi, Lx, Ly)){
			bpi->x  += (bpi->rk4.k1x +2*bpi->rk4.k2x +2*bpi->rk4.k3x +bpi->rk4.k4x )/6.;
			bpi->y  += (bpi->rk4.k1y +2*bpi->rk4.k2y +2*bpi->rk4.k3y +bpi->rk4.k4y )/6.;
			bpi->vx += (bpi->rk4.k1vx+2*bpi->rk4.k2vx+2*bpi->rk4.k3vx+bpi->rk4.k4vx)/6.;
			bpi->vy += (bpi->rk4.k1vy+2*bpi->rk4.k2vy+2*bpi->rk4.k3vy+bpi->rk4.k4vy)/6.;
			bpi->ie += (bpi->rk4.k1ie+2*bpi->rk4.k2ie+2*bpi->rk4.k3ie+bpi->rk4.k4ie)/6.;
		}
	}
	postStage(simpar);
	bp = VORORK4_TBP(simpar);
	migrateTreeVorork4Particles(simpar);

	/* Update volume & density with final positions — uses cell-based method */
	exam2dUpdateVol(simpar, paddingAllTreeParticles,
			searchCellRk4Neighbors2D, findCellRk4BP2D, mkLinkedList2D_oExam);

	bp = VORORK4_TBP(simpar);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		if(targetBP(bpi, Lx, Ly)){
			if(GAS_Kappa(simpar) < 0){
				bpi->w2 = -GAS_Kappa(simpar);
			}
			else if(GAS_Kappa(simpar) > 0){
				bpi->w2hydro = getw2forHydroParticle(simpar, bpi, Dtime);
				applyW2Controls(simpar, bpi, Dtime);
			}
		}
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NP(simpar);i++){
		treevorork4particletype *bpi = (treevorork4particletype*)((char*)bp + i*p_sz);
		if(targetBP(bpi, Lx, Ly)){
			bpi->den = bpi->mass/bpi->volume;
			bpi->pressure = bpi->ie/bpi->volume * (Gamma-1);
			if(bpi->pressure <= 0){
				bpi->pressure = 1e-6;
				bpi->ie = bpi->pressure * bpi->volume / (Gamma-1);
			}
			bpi->csound = sqrt(Gamma*bpi->pressure/bpi->den);
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
//
	int nx,ny;
	nx = NX(simpar);
	ny = NY(simpar);


    int i;
    postype Dtime,dt;
    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) {
		bp[i].w2old = bp[i].w2;
		bp[i].rk4.w2backup = bp[i].w2;
	}
	// Runge-Kutta 4-th order time evolution of r and vr
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D,1);
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
		if(bp[i].ie <0) DEBUGPRINT("P%d has p(%ld,%ld) : ie= %g, die= %g :: x,y= %g %g vx,y = %g %g\n", MYID(simpar), PINDX(bp+i)%nx, PINDX(bp+i)/nx, bp[i].ie, bp[i].die, bp[i].x, bp[i].y, bp[i].vx, bp[i].vy);
        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }

	migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
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
		if(bp[i].ie <0) DEBUGPRINT("P%d has p(%ld,%ld) : ie= %g, die= %g %g :: x,y= %g %g vx,y = %g %g\n", MYID(simpar), PINDX(bp+i)%nx, PINDX(bp+i)/nx, bp[i].ie, bp[i].rk4.k2ie, bp[i].rk4.k1ie, bp[i].x, bp[i].y, bp[i].vx, bp[i].vy);

        bp[i].x = fmod(bp[i].x+Lx, Lx);
        bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);
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
		if(bp[i].ie <0) DEBUGPRINT("P%d has p(%ld,%ld) : ie= %g, die= %g %g :: x,y= %g %g vx,y = %g %g\n", MYID(simpar), PINDX(bp+i)%nx, PINDX(bp+i)/nx, bp[i].ie, bp[i].rk4.k3ie, bp[i].rk4.k2ie, bp[i].x, bp[i].y, bp[i].vx, bp[i].vy);

        bp[i].x = fmod(bp[i].x+Lx, Lx);
		bp[i].y = fmod(bp[i].y+Ly, Ly);
    }
    // migrate particles between mpi ranks
    migrateTreeVorork4Particles(simpar);
	bp = VORORK4_TBP(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP,
			mkLinkedList2D, Dtime);



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
	//
	//

    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly)) {
			bp[i].x  += (bp[i].rk4.k1x +2*bp[i].rk4.k2x +2*bp[i].rk4.k3x +bp[i].rk4.k4x )/6.; 
			bp[i].y  += (bp[i].rk4.k1y +2*bp[i].rk4.k2y +2*bp[i].rk4.k3y +bp[i].rk4.k4y )/6.; 
			bp[i].vx += (bp[i].rk4.k1vx+2*bp[i].rk4.k2vx+2*bp[i].rk4.k3vx+bp[i].rk4.k4vx)/6.; 
			bp[i].vy += (bp[i].rk4.k1vy+2*bp[i].rk4.k2vy+2*bp[i].rk4.k3vy+bp[i].rk4.k4vy)/6.; 
			bp[i].ie += (bp[i].rk4.k1ie+2*bp[i].rk4.k2ie+2*bp[i].rk4.k3ie+bp[i].rk4.k4ie)/6.;
		   if(bp[i].ie <0) DEBUGPRINT("P%d has p(%ld,%ld) : ie= %g, die= %g %g %g %g :: x,y= %g %g vx,y = %g %g\n", MYID(simpar), PINDX(bp+i)%nx, PINDX(bp+i)/nx, bp[i].ie, 
				   bp[i].rk4.k1ie, bp[i].rk4.k2ie, 
				   bp[i].rk4.k3ie, bp[i].rk4.k4ie, 
				   bp[i].x, bp[i].y, bp[i].vx, bp[i].vy);
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
	// split getw2forHydroParticle and pressure update because the former uses the latter.
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
			if(GAS_Kappa(simpar) <0) { 
				bp[i].w2 = -GAS_Kappa(simpar); 
			} 
			else if (GAS_Kappa(simpar) >0){ 
          		bp[i].w2hydro = getw2forHydroParticle(simpar,(bp+i),Dtime);
				applyW2Controls(simpar, bp+i, Dtime);
			} 
		}
	}
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(targetBP(bp+i, Lx, Ly))
		{
			bp[i].den = bp[i].mass/bp[i].volume; 
			bp[i].pressure = bp[i].ie/bp[i].volume * (Gamma-1); 
			bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den); 
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
