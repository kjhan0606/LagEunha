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
#include "nnost.h"
#include "gnnost.h"
#include "exam.h"
#include "exam2d.h"




void det2d_dpqStress(
		SimParameters *simpar, 
		void (*paddingAllTreeParticles)(SimParameters *, postype) ){
	treevorostressparticletype *bp = VORO_TBP(simpar);
	int np = VORO_NP(simpar);
	size_t i;

	{
		postype cellsize = HydroGridSize(simpar);
    	paddingAllTreeParticles(simpar, cellsize);
	}
	int mp = VORO_NPAD(simpar);

	TPtlStruct *ptl = (TPtlStruct*)malloc(sizeof(TPtlStruct)*(np+mp));
	int nnode = np+mp;
	TStruct *tree = (TStruct *)malloc(sizeof(TStruct)*nnode);

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

	free(ptl);
	free(tree);
	free(VORO_TBPP(simpar));
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

void det3d_dpqStress(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype)
		){
    treevorostressparticletype *bp = VORO_TBP(simpar);
    int np = VORO_NP(simpar);

	{
		postype cellsize = HydroGridSize(simpar);
    	paddingAllTreeParticles(simpar, cellsize);

	}
	int mp = VORO_NPAD(simpar);

    TPtlStruct *ptl = (TPtlStruct*)malloc(sizeof(TPtlStruct)*(np+mp));
    int nnode = np+mp;
    TStruct *tree = (TStruct *)malloc(sizeof(TStruct)*nnode);
    size_t i;

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
    free(ptl);
    free(tree);
	free(VORO_TBPP(simpar));
}
void det2d_dpqStressRK4(
		SimParameters *simpar, 
		void (*paddingAllTreeParticles)(SimParameters *, postype) ){
	treevorostressrk4particletype *bp = VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	size_t i;

	{
		postype cellsize = HydroGridSize(simpar);
    	paddingAllTreeParticles(simpar, cellsize);
	}
	int mp = VORO_NPAD(simpar);

	TPtlStruct *ptl = (TPtlStruct*)malloc(sizeof(TPtlStruct)*(np+mp));
	int nnode = np+mp;
	TStruct *tree = (TStruct *)malloc(sizeof(TStruct)*nnode);

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

	free(ptl);
	free(tree);
	free(VORORK4_TBPP(simpar));
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

void det3d_dpqStressRK4(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype)
		){
    treevorostressrk4particletype *bp = VORORK4_TBP(simpar);
    int np = VORO_NP(simpar);

	{
		postype cellsize = HydroGridSize(simpar);
    	paddingAllTreeParticles(simpar, cellsize);

	}
	int mp = VORO_NPAD(simpar);

    TPtlStruct *ptl = (TPtlStruct*)malloc(sizeof(TPtlStruct)*(np+mp));
    int nnode = np+mp;
    TStruct *tree = (TStruct *)malloc(sizeof(TStruct)*nnode);
    size_t i;

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
    free(ptl);
    free(tree);
	free(VORORK4_TBPP(simpar));
}

Voro2D_point *searchCellStressRk4Neighbors2D(
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
    neigh = (Voro2D_point*)malloc(sizeof(Voro2D_point)*np);
    np = 0;
    for(j=iy-1;j<=iy+1;j++){
		if( j<0 || j >= my) continue;
        for(i=ix-1;i<=ix+1;i++){
			if(i<0 || i >= mx) continue;
            struct linkedlisttype *tmp = cells[i+mx*j].link;
            while(tmp){
                treevorostressrk4particletype *tt = (treevorostressrk4particletype*)tmp;
                neigh[np].x = tt->x;
                neigh[np].y = tt->y;
                neigh[np].indx = PINDX(tt);
                neigh[np].csound = tt->csound;
                neigh[np].pressure = tt->pressure;
                neigh[np].w2 = tt->w2;
                neigh[np].bp = (void*)tt;
                np++;
                tmp = tmp->next;
            }

        }
    }
    *nneigh = np;
    return neigh;
}

treevorostressrk4particletype *findCellStressRk4BP2D(SimParameters *simpar, int ix, int iy, int *mp){
    int i,j,k;
    int np;
    treevorostressrk4particletype *res;
    int mx,my;
    mx = BASICCELL_MX(simpar);
    my = BASICCELL_MY(simpar);
    CellType *cells = VORO_BASICCELL(simpar);

    np = cells[ix+mx*iy].nmem;
    res = (treevorostressrk4particletype*)malloc(sizeof(treevorostressrk4particletype)*np);
    np = 0;
    size_t ipixel = ix+mx*iy;
    struct linkedlisttype *tmp = cells[ipixel].link;
    while(tmp){
        if(!IS_FLAG_ON(tmp,BoundaryGhostflag)){
            treevorostressrk4particletype *tt = (treevorostressrk4particletype*)tmp;
//            if(tt->u4if.Flag[ENDIAN_OFFSET]> 0)
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


void mkLinkedList2DStress(
		SimParameters *simpar, 
		postype cellsize, 
		postype xmin, 
		postype ymin,
		postype xmax, 
		postype ymax, 
		void (*paddingTreeAllParticle)(SimParameters *, postype)
		){ 
	treevorostressrk4particletype *bp = VORORK4_TBP(simpar); 
	int np = VORO_NP(simpar);
    int i,j,k;
    int mx,my;

    mx = BASICCELL_MX(simpar);
    my = BASICCELL_MY(simpar);


	
//	DEBUGPRINT("P%d takes the values cellsize= %g xymin= %g %g xymax= %g %g\n", MYID(simpar), cellsize, xmin,ymin,xmax,ymax);


	CellType *cells = VORO_BASICCELL(simpar); 

	for(i=0;i<mx*my;i++) {
        cells[i].link = NULL;
        cells[i].nmem = 0;
    }
	for(i=0;i<VORO_NP(simpar);i++){
        int ix,iy;
        ix = ((bp[i].x-xmin)/cellsize);
        iy = ((bp[i].y-ymin)/cellsize);
        size_t index = ix+mx*iy;
        if(index >=mx*my){
            DEBUGPRINT("error detected: %ld : %d %d : %g %g %g %g\n",
                    index, mx,my, bp[i].x, bp[i].y, xmin,ymin);
        }
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)(bp+i);
        cells[index].nmem ++;
        bp[i].next = tmp;
    }
    paddingTreeAllParticle(simpar, cellsize);
    bp = VORORK4_TBPP(simpar);
    postype Xmin,Ymin,Xmax,Ymax;
    Xmin = Ymin = 1.e10;
    Xmax = Ymax = -1.e10;
    for(i=0;i<VORO_NPAD(simpar);i++){
        int ix,iy;
        ix = ((bp[i].x-xmin)/cellsize);
        iy = ((bp[i].y-ymin)/cellsize);
        Xmin = MIN(Xmin,bp[i].x);
        Ymin = MIN(Ymin,bp[i].y);
        Xmax = MAX(Xmax,bp[i].x);
        Ymax = MAX(Ymax,bp[i].y);
        size_t index = ix+mx*iy;
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)(bp+i);
        cells[index].nmem ++;
        bp[i].next = tmp;
    }
	/*
	DEBUGPRINT("P%d has padding particles : %g %g : %g %g ::: Lx/y = %g %g\n", MYID(simpar),
			Xmin,Xmax,Ymin,Ymax, 
			SIMBOX(simpar).x.max, SIMBOX(simpar).y.max);
			*/
}



void updateDenW2Pressure2DStress(
		SimParameters *simpar, 
		postype xmin, postype ymin, postype xmax, postype ymax,
		 postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorostressrk4particletype *(*find2DCellBP)(SimParameters *, int , int , int *)
		){
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	treevorostressrk4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int isave = -1;


	det2d_dpqStressRK4(simpar, paddingAllTreeParticles);
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


	CellType *cells = VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*mx*my);
	mkLinkedList2DStress(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingAllTreeParticles);

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
            treevorostressrk4particletype *p = find2DCellBP(simpar,ix,iy,&np);
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
				treevorostressrk4particletype *ibp = p[i].bp;
                ibp->volume = Area2DPolygon(vorocorner, mp);
                ibp->den = ibp->mass/ibp->volume;
				/*
                ibp->csound = sqrt(Gamma*ibp->pressure/ibp->den);
                ibp->pressure = ibp->ie/ibp->volume*(Gamma-1);
				*/
            }
            free(p);free(neighbors); free(neighwork);

        }
        free(vorocorner);
    }
	for(i=0;i<nbp;i++){
		bp[i].pressure = bp[i].ie/bp[i].volume*(Gamma-1);
		bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
	}
	// finalize mkLinkedList2DStress by freeing memory spaces (cell & boundary ghost particles)
	free(VORO_BASICCELL(simpar));
	free(VORORK4_TBPP(simpar));
}

double getAccVoro2DStress(SimParameters *simpar, postype xmin, postype ymin, 
		postype xmax, postype ymax,
		postype OrderOfAccuracy, postype Courant, postype Gamma,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighboringBP)(SimParameters *, int, int, int *),
		treevorostressrk4particletype *(*find2DCellBP)(SimParameters *, int , int , int *)
		){
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	treevorostressrk4particletype *bp = VORORK4_TBP(simpar);
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
	CellType *cells = VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*mx*my);
	mkLinkedList2DStress(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);


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
    for(iy=0;iy<my;iy++){
    	int mp=1000;
    	Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
    	postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
        for(ix=0;ix<mx;ix++){
            int np;
            treevorostressrk4particletype *p = find2DCellBP(simpar,ix,iy,&np);
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

				treevorostressrk4particletype *ibp = p[i].bp;
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
						treevorostressrk4particletype *jbp =  (treevorostressrk4particletype*)(neighwork[tmp->upperrelated].bp);
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
								ibp_pressure,jbp_pressure, tmp, neighwork, simpar, KH);

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

						/*
						if(center.indx == 43648){
							DEBUGPRINT("p43648 has fx/y= %g %g pi= %g dS= %g %g\n", fx,fy,pi,dS.x,dS.y);
						}
						*/

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
	free(VORO_BASICCELL(simpar));
	free(VORORK4_TBPP(simpar));
	{
		postype TDtime;
        MPI_Reduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MAX, 0, MPI_COMM(simpar));
        if(MYID(simpar) == 0) Dtime = TDtime;
        MPI_Bcast(&Dtime, 1, MPI_POSTYPE,  0, MPI_COMM(simpar));
	}
	return Dtime;
}
void exam2d_findVolStress(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorostressrk4particletype *(*find2DCellBP)(SimParameters *, int , int , int *)
		){
	treevorostressrk4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
//  determine the minimum dpq for each particle by updating the w2ceil & w2
	det2d_dpqStressRK4(simpar,paddingAllTreeParticles); 
	/*
	DEBUGPRINT("P%d passed the initial minimum dpq %g %g %g\n", MYID(simpar), 
			bp[0].w2ceil, bp[nbp/2].w2ceil, bp[nbp-1].w2ceil);
			*/


	postype xmin,ymin,zmin,xmax,ymax,zmax, cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
	// xmin,ymin,xmax,ymax are the boundaries of the local domain.
	xmin = Xmin_HydroExam(simpar)-cellsize;
    ymin = Ymin_HydroExam(simpar)-cellsize;
    xmax = Xmax_HydroExam(simpar)+cellsize;
    ymax = Ymax_HydroExam(simpar)+cellsize;
	DEBUGPRINT("P%d has cell info. %g : %g %g %g %g\n", MYID(simpar), cellsize, xmin,ymin,xmax,ymax);

	int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);
	// prepare the linked-list cells for mkLinkedList2DStress()  and tree findings below
	CellType *cells = VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*mx*my);
	// building the linked list with "paddingAllTreeParticles()" which pads the domain
	// with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
	mkLinkedList2DStress(simpar, cellsize,xmin,ymin,xmax,ymax, paddingAllTreeParticles);

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
            treevorostressrk4particletype *p = find2DCellBP(simpar,ix,iy,&np);
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
				treevorostressrk4particletype *ibp = p[i].bp;
				ibp->volume = Area2DPolygon(vorocorner,mp);
				ibp->den = ibp->mass / ibp->volume;
			}
			free(neighwork);
			free(neighbors); 
			free(p);
		}
		free(vorocorner);
	
	}
	free(VORO_BASICCELL(simpar));

	VORO_NPAD(simpar) = 0;
	free(VORORK4_TBPP(simpar));
}
double exam2d_vph_StressRk4(
		SimParameters *simpar,
		void (*paddingAllTreeParticles)(SimParameters *, postype),
		double (*measureW2)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*find2DNeighborBP)(SimParameters *, int, int, int *),
		treevorostressrk4particletype *(*find2DCellBP)(SimParameters *, int , int , int *)
		){
	treevorostressrk4particletype *bp = VORORK4_TBP(simpar);
	postype xmin,ymin,xmax,ymax,OrderOfAccuracy;
	postype Courant = GAS_COURANT(simpar);
	postype Gamma = GAS_GAMMA(simpar);

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = Xmin_HydroExam(simpar)-cellsize;
	ymin = Ymin_HydroExam(simpar)-cellsize;
	xmax = Xmax_HydroExam(simpar)+cellsize;
	ymax = Ymax_HydroExam(simpar)+cellsize;

	OrderOfAccuracy = VoroAccuracyOrder(simpar);
	DEBUGPRINT("P%d has x/y minmax= %g %g : %g %g ::: %g\n", MYID(simpar), xmin,ymin,xmax,ymax, cellsize);


    int i;
    postype Dtime,dt;
    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) bp[i].rk4.w2backup = bp[i].w2;
	// Runge-Kutta 4-th order time evolution of r and vr
	updateDenW2Pressure2DStress(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP);
    Dtime = getAccVoro2DStress(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma, paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP);

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
	updateDenW2Pressure2DStress(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP);
    dt = getAccVoro2DStress(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
        paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP);
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
	updateDenW2Pressure2DStress(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP);
	dt = getAccVoro2DStress(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP);
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
	updateDenW2Pressure2DStress(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingAllTreeParticles,find2DNeighborBP, find2DCellBP);
    dt = getAccVoro2DStress(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
		paddingAllTreeParticles,
        find2DNeighborBP, find2DCellBP);
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
	bp = VORORK4_TBP(simpar);
	// update density with one-step evolved positions 
    exam2d_findVolStress(simpar,paddingAllTreeParticles,find2DNeighborBP, find2DCellBP);

	postype dmean = SIMBOX(simpar).x.max/NX(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].den = bp[i].mass/bp[i].volume;
        bp[i].pressure = bp[i].ie/bp[i].volume * (Gamma-1);
        bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
        bp[i].w2 = measureW2(simpar,dmean, bp[i].pressure, bp[i].den);
        bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
    }
    return Dtime;
}

void dumpStressRk4Data2D(SimParameters *simpar, int nstep, postype t, postype dt){
	treevorostressrk4particletype *bp = VORORK4_TBP(simpar);
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
			fwrite(bp, sizeof(treevorostressrk4particletype), np, wp);
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


void readRk4StressData2D(SimParameters *simpar, postype *t, int nstep){
	char infile[190];
	sprintf(infile,"glout.%.6d.dat", nstep);
	int myid = MYID(simpar);
	if(myid ==0){
		printf("P0: is now reading starting file: %s\n", infile);
		FILE *fp = fopen(infile,"r");
		int np;
		fread(&np, sizeof(int), 1, fp);
		VORO_TNP(simpar) = VORO_NP(simpar) = np;
		VORORK4_TBP(simpar) = (treevorostressrk4particletype*)malloc(sizeof(treevorostressrk4particletype)*np);
		fread(VORORK4_TBP(simpar),sizeof(treevorostressrk4particletype), np,fp);
		fread(t,sizeof(postype), 1,fp);
		fclose(fp);
	}
	else {
		VORO_NP(simpar) = 0;
		VORORK4_TBP(simpar) = (treevorostressrk4particletype*)malloc(sizeof(treevorostressrk4particletype)*100);;

	}
	DEBUGPRINT("P%d: after reading iput data\n", myid);
    migrateTreeVorork4Particles(simpar);
	MPI_Bcast(t, 1, MPI_POSTYPE,0, MPI_COMM(simpar));
	MPI_Barrier(MPI_COMM(simpar));
	BASICCELL_CELLWIDTH(simpar) = HydroGridSize(simpar);
}

double exam2d_vph_Stress(SimParameters *simpar, 
		double (*measureW2Exam2D)(SimParameters *, double, double, double))
{
    postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
    treevorostressrk4particletype *bp = VORORK4_TBP(simpar);
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

    // prepare the linked-list cells for mkLinkedList2DStress()  and tree findings below
    CellType *cells = VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*mx*my);
    // building the linked list with "paddingTreeVorork4Particles()" which pads the domain
    // with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
    mkLinkedList2DStress(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingTreeVorork4Particles);

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
    for(iy=0;iy<my;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=0;ix<mx;ix++){
            int np;
            treevorostressrk4particletype *p = findCellStressRk4BP2D(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = searchCellStressRk4Neighbors2D(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.csound = p[i].csound;
                center.w2 = p[i].w2;

                treevorostressrk4particletype *ibp = p[i].bp;
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
                        treevorostressrk4particletype *jbp =  (treevorostressrk4particletype*)(neighwork[tmp->upperrelated].bp);
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

                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,KH);

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
                        /*
                        if(center.indx == 15872){
                            DEBUGPRINT("p2 has p%ld fx/y= %g %g pi= %g (mu= %g rvel= %g etavis= %g dramp= %g eps2vis= %g) "
                                    "dS= %g %g\n", 
                                    PINDX(jbp),fx,fy,pi,mu,rvel, etavis, dramp, eps2vis,
                                    dS.x,dS.y);
                        }
                        */


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
            free(p);free(neighbors); free(neighwork);

        }
        free(vorocorner);
    }
    free(VORORK4_TBPP(simpar));

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
    {
        postype xmin,ymin,xmax,ymax;
        xmin = ymin = 1.e20;
        xmax = ymax =-1.e20;
        for(i=0;i<VORO_NP(simpar);i++){
            xmin = MIN(xmin, bp[i].x);
            ymin = MIN(ymin, bp[i].y);
            xmax = MAX(xmax, bp[i].x);
            ymax = MAX(ymax, bp[i].y);
        }
//      DEBUGPRINT("P%d has xmin/max= %g %g :: ymin/max= %g %g : Lx/y= %g %g\n", MYID(simpar), xmin,xmax,ymin,ymax, Lx,Ly);
    }

    if(GAS_Kappa(simpar)>0) det2d_dpqStressRK4(simpar,paddingTreeVorork4Particles);
    // Update the linked list 
    mkLinkedList2DStress(simpar, cellsize, xmin,ymin,xmax,ymax,  paddingTreeVorork4Particles);


    // repositioning the memory space 
    bp = VORORK4_TBP(simpar);

#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(iy=0;iy<my;iy++){
        int mp=2000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=0;ix<mx;ix++){
            int np;
            treevorostressrk4particletype *p = findCellStressRk4BP2D(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = searchCellStressRk4Neighbors2D(simpar,ix,iy,&nneigh);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = PINDX(p+i);
                center.w2 = p[i].w2;

                treevorostressrk4particletype *ibp = p[i].bp;
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
                double die, dte,dke,fx,fy;
                die = dte = dke = fx = fy = 0;
                do {
                    if(tmp->upperrelated >=0){
                        // These are to use register 
                        treevorostressrk4particletype *jbp =  (treevorostressrk4particletype*)(neighwork[tmp->upperrelated].bp);
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
//                      postype pi = Half*(jbp_pressure + ibp_pressure);
                        postype pi = VoroRK4_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,KH);

                        Voro2D_point ui,ua,ub;
//                        ui.x = Half*(jbp_vx - ibp_vx);
 //                       ui.y = Half*(jbp_vy - ibp_vy);
                        ui = get2dUpqradStressRk4(ibp,jbp,dtold);

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

						/*
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
						*/

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
    free(VORORK4_TBPP(simpar));
    free(VORO_BASICCELL(simpar));

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
//      bp[i].w2old = bp[i].w2;
        if(GAS_Kappa(simpar) <0) {
            bp[i].w2 = -GAS_Kappa(simpar);
        }
        else if (GAS_Kappa(simpar) >0){
            bp[i].w2 = measureW2Exam2D(simpar,dmean, bp[i].pressure, bp[i].den);
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
//  DEBUGPRINT("+P%d has np = %ld\n", MYID(simpar), VORO_NP(simpar));
    return Dtime;
}
