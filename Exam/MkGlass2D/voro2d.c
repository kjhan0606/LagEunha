#include<math.h>
#include "eunha.h"
#include "voro.h"
#include "cosmology.h"
#include "nnost.h"
#include "gnnost.h"
#include "gl2d.h"
#include "exam.h"
#include "exam2d.h"
#include<omp.h>

static postype p0 = 1;
static postype rho0 = 1;
postype gl2d_w2Measure(SimParameters *simpar, double dmean, double pressure,double rho){
//	double res = 0.55*dmean*pow(pressure/p0*rho0/rho,0.5);
	double res = GL2D_Kappa(simpar)*dmean*pow(pressure/p0*rho0/rho,0.5);
    double res2 = res*res;
    return res2;
}

postype gl2d_Voro_Pressure2D( postype pi, postype pj, Voro2D_Corner *tmp, 
		Voro2D_point *neighwork, SimParameters *simpar){
	postype res;
	res = ( 0.5-GL2D_OA(simpar)/3.)*(pi+pj);
	postype pkm, pkp;
	pkm = ((treevorork4particletype*) (neighwork[(tmp->upperlink)->upperrelated].bp))->pressure;
	pkp = ((treevorork4particletype*) (neighwork[(tmp           )->lowerrelated].bp))->pressure;
	res = res + GL2D_OA(simpar)/3.* (pkm+pkp);
	return res;
}

void gl2d_findVol(SimParameters *simpar){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
//  determine the minimum dpq for each particle by updating the w2ceil & w2
	det2d_dpq(simpar,paddingTreeVoroParticles); 
	/*
	DEBUGPRINT("P%d passed the initial minimum dpq %g %g %g\n", MYID(simpar), 
			bp[0].w2ceil, bp[nbp/2].w2ceil, bp[nbp-1].w2ceil);
			*/


	postype xmin,ymin,zmin,xmax,ymax,zmax, cellsize;
	cellsize = BASICCELL_CELLWIDTH(simpar) = GL2D_GridSize(simpar);
	// xmin,ymin,xmax,ymax are the boundaries of the local domain.
	xmin = GL2D_XMIN(simpar)-cellsize;
    ymin = GL2D_YMIN(simpar)-cellsize;
    xmax = GL2D_XMAX(simpar)+cellsize;
    ymax = GL2D_YMAX(simpar)+cellsize;
	DEBUGPRINT("P%d has cell info. %g : %g %g %g %g\n", MYID(simpar), cellsize, xmin,ymin,xmax,ymax);

	int mx = BASICCELL_MX(simpar) = ceil((xmax-xmin)/cellsize);
    int my = BASICCELL_MY(simpar) = ceil((ymax-ymin)/cellsize);
	// prepare the linked-list cells for mkLinkedList2D()  and tree findings below
	CellType *cells = VORO_BASICCELL(simpar)= (CellType*)malloc(sizeof(CellType)*mx*my);
	// building the linked list with "paddingTreeVoroParticles()" which pads the domain
	// with the tree voro particles defined in ../Exam/mpirks.exam2d.c.
	mkLinkedList2D(simpar, cellsize,xmin,ymin,xmax,ymax, paddingTreeVoroParticles);




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
            treevorork4particletype *p = gl2d_Voro2D_FindCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = gl2d_Voro2D_FindNeighbor(simpar,ix,iy,&nneigh);
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


double gl2d_vph2D_rk4(SimParameters *simpar){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	postype xmin,ymin,xmax,ymax,OrderOfAccuracy;

	postype cellsize = BASICCELL_CELLWIDTH(simpar);
	xmin = GL2D_XMIN(simpar)-cellsize;
	ymin = GL2D_YMIN(simpar)-cellsize;
	xmax = GL2D_XMAX(simpar)+cellsize;
	ymax = GL2D_YMAX(simpar)+cellsize;
	OrderOfAccuracy = GL2D_OA(simpar);


    int i;
    postype Dtime,dt;
    postype Lx = SIMBOX(simpar).x.max;
    postype Ly = SIMBOX(simpar).y.max;

	for(i=0;i<VORO_NP(simpar);i++) bp[i].rk4.w2backup = bp[i].w2;
	// Runge-Kutta 4-th order time evolution of r and vr
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingTreeVoroParticles, searchCellNeighbors2D, findCellBP2D);
    Dtime = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma, paddingTreeVoroParticles,
		searchCellNeighbors2D, findCellBP2D);

	bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
        if(bp[i].u4if.Flag[ENDIAN_OFFSET] >0){
            bp[i].rk4.k1x = bp[i].vx*Dtime;
            bp[i].rk4.k1y = bp[i].vy*Dtime;
            bp[i].rk4.k1vx = bp[i].ax*Dtime;
            bp[i].rk4.k1vy = bp[i].ay*Dtime;
            bp[i].x += bp[i].rk4.k1x*0.5;
            bp[i].y += bp[i].rk4.k1y*0.5;
			bp[i].x = fmod(bp[i].x+Lx, Lx);
			bp[i].y = fmod(bp[i].y+Ly, Ly);
        }
    }
	migrateTreeVoroParticles(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingTreeVoroParticles, searchCellNeighbors2D,findCellBP2D);
    dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
        paddingTreeVoroParticles,
		searchCellNeighbors2D, findCellBP2D);
	bp = VORORK4_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
        if(bp[i].u4if.Flag[ENDIAN_OFFSET] >0){
            bp[i].rk4.k2x = (bp[i].vx + bp[i].rk4.k1vx*0.5)*Dtime;
            bp[i].rk4.k2y = (bp[i].vy + bp[i].rk4.k1vy*0.5)*Dtime;
            bp[i].rk4.k2vx = bp[i].ax*Dtime;
            bp[i].rk4.k2vy = bp[i].ay*Dtime;
            bp[i].x += (bp[i].rk4.k2x-bp[i].rk4.k1x)*0.5;
            bp[i].y += (bp[i].rk4.k2y-bp[i].rk4.k1y)*0.5;
			bp[i].x = fmod(bp[i].x+Lx, Lx);
			bp[i].y = fmod(bp[i].y+Ly, Ly);
        }
    }
    // migrate particles between mpi ranks 
    migrateTreeVoroParticles(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingTreeVoroParticles, searchCellNeighbors2D,  findCellBP2D);
	dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
        paddingTreeVoroParticles,
		searchCellNeighbors2D, findCellBP2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
        if(bp[i].u4if.Flag[ENDIAN_OFFSET] >0){
            bp[i].rk4.k3x = (bp[i].vx + bp[i].rk4.k2vx*0.5)*Dtime;
            bp[i].rk4.k3y = (bp[i].vy + bp[i].rk4.k2vy*0.5)*Dtime;
            bp[i].rk4.k3vx = bp[i].ax*Dtime;
            bp[i].rk4.k3vy = bp[i].ay*Dtime;
            bp[i].x += bp[i].rk4.k3x-bp[i].rk4.k2x*0.5;
            bp[i].y += bp[i].rk4.k3y-bp[i].rk4.k2y*0.5;
			bp[i].x = fmod(bp[i].x+Lx, Lx);
			bp[i].y = fmod(bp[i].y+Ly, Ly);
        }
    }
    // migrate particles between mpi ranks
    migrateTreeVoroParticles(simpar);

	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingTreeVoroParticles, searchCellNeighbors2D,  findCellBP2D);
    dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
        paddingTreeVoroParticles,
		searchCellNeighbors2D, findCellBP2D);
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
        if(bp[i].u4if.Flag[ENDIAN_OFFSET] >0){
            bp[i].rk4.k4x = (bp[i].vx + bp[i].rk4.k3vx)*Dtime;
            bp[i].rk4.k4y = (bp[i].vy + bp[i].rk4.k3vy)*Dtime;
            bp[i].rk4.k4vx = bp[i].ax*Dtime;
            bp[i].rk4.k4vy = bp[i].ay*Dtime;
            bp[i].x -= (bp[i].rk4.k3x);
            bp[i].y -= (bp[i].rk4.k3y);
			bp[i].x = fmod(bp[i].x+Lx, Lx);
			bp[i].y = fmod(bp[i].y+Ly, Ly);
        }
    }
	// migrate particles between mpi ranks 
    migrateTreeVoroParticles(simpar);

    // to finalize ke & ie updates 
	for(i=0;i<VORO_NP(simpar);i++) bp[i].w2= bp[i].rk4.w2backup;
	updateDenW2Pressure2D(simpar,xmin,ymin,xmax,ymax,
			Gamma, paddingTreeVoroParticles, searchCellNeighbors2D,  findCellBP2D);
    dt = getAccVoro2D(simpar, xmin, ymin, xmax, ymax,
        OrderOfAccuracy, Courant, Gamma,
        paddingTreeVoroParticles,
		searchCellNeighbors2D, findCellBP2D);

    // Update position and velocity using RK4
    bp = VORORK4_TBP(simpar);
    for(i=0;i<VORO_NP(simpar);i++){
        if(bp[i].u4if.Flag[ENDIAN_OFFSET] >0)
		{
            bp[i].x  += (bp[i].rk4.k1x +2*bp[i].rk4.k2x +2*bp[i].rk4.k3x +bp[i].rk4.k4x )/6.;
            bp[i].y  += (bp[i].rk4.k1y +2*bp[i].rk4.k2y +2*bp[i].rk4.k3y +bp[i].rk4.k4y )/6.;
            bp[i].x = fmod(bp[i].x+Lx,Lx);
            bp[i].y = fmod(bp[i].y+Ly,Ly);

            bp[i].vx += (bp[i].rk4.k1vx+2*bp[i].rk4.k2vx+2*bp[i].rk4.k3vx+bp[i].rk4.k4vx)/6.;
            bp[i].vy += (bp[i].rk4.k1vy+2*bp[i].rk4.k2vy+2*bp[i].rk4.k3vy+bp[i].rk4.k4vy)/6.;


        }
    }

    /* migrate particles between mpi ranks */
    migrateTreeVoroParticles(simpar);
    bp = VORORK4_TBP(simpar);

	// This is comment out to skip the evolution of internal energy.
    for(i=0;i<VORO_NP(simpar);i++){
        if(bp[i].u4if.Flag[ENDIAN_OFFSET] >0){
            bp[i].ke  = Half*bp[i].mass*(bp[i].vx*bp[i].vx + bp[i].vy*bp[i].vy);
            bp[i].ie += bp[i].die * Dtime;
            bp[i].te = bp[i].ie + bp[i].ke;
            bp[i].ie = MAX(0., bp[i].ie);
            if(bp[i].volume ==0)
                DEBUGPRINT("P%d has wrong volume at i= %d with den= %g at xy= %g %g\n", MYID(simpar),
                        i,bp[i].den, bp[i].x, bp[i].y);
        }
    }
	/* update density with one-step evolved positions */
	gl2d_findVol(simpar);

	postype dmean = SIMBOX(simpar).x.max/NX(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		bp[i].den = bp[i].mass/bp[i].volume;
		bp[i].pressure = bp[i].ie/bp[i].volume * (Gamma-1);
		bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
		bp[i].w2 = gl2d_w2Measure(simpar,dmean, bp[i].pressure, bp[i].den);
		bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
		bp[i].vx = bp[i].vy = 0; // inserted to drag the motion.
	}
    return Dtime;
}



double gl2d_vph2D(SimParameters *simpar){
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int isave = -1;
    postype Dtime = 1.e10;



	det2d_dpq(simpar,paddingTreeVoroParticles);

	gl2d_MkLinkedList(simpar);
	/*
	if(RT_SIMBOX(simpar).x.min !=0. || RT_SIMBOX(simpar).y.min !=0.){
		DEBUGPRINT("Error in the simulation minimum %g %g\n",
				RT_SIMBOX(simpar).x.min, RT_SIMBOX(simpar).y.min);
		MPI_Finalize();
		exit(0);
	}
	*/
	postype Lx = SIMBOX(simpar).x.max;
	postype Ly = SIMBOX(simpar).y.max;
	postype xmin,ymin,zmin,xmax,ymax,zmax, pwidth;
	pwidth = BASICCELL_CELLWIDTH(simpar);
	xmin = GL2D_XMIN(simpar)-pwidth;
    ymin = GL2D_YMIN(simpar)-pwidth;
    xmax = GL2D_XMAX(simpar)+pwidth;
    ymax = GL2D_YMAX(simpar)+pwidth;
    int mx = BASICCELL_MX(simpar);
    int my = BASICCELL_MY(simpar);

	float alphavis = GAS_AlphaVis(simpar);
	float betavis = GAS_BetaVis(simpar);
	float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
	float eps2vis = GAS_EPSVIS(simpar)*GAS_EPSVIS(simpar);


	DEBUGPRINT("P%d has viscosity parameters %g %g\n", MYID(simpar),alphavis, betavis);


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
            treevorork4particletype *p = gl2d_Voro2D_FindCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = gl2d_Voro2D_FindNeighbor(simpar,ix,iy,&nneigh);
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
				ibp->w2 = ibp->w2ceil = center.w2ceil;
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

						postype pi = Voro_Pressure2D(
								ibp_pressure,jbp_pressure, tmp, neighwork, simpar,GL2D);

						Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

						postype dramp = sqrt(Vec2DDotP(&dr, &dr));
						Voro2D_point er;
						er.x = dr.x/dramp;
						er.y = dr.y/dramp;

						Voro2D_point ui,ua,ub;
						ui.x = Half*(jbp_vx - ibp_vx); 
						ui.y = Half*(jbp_vy - ibp_vy);

						if(0){
							/* for the bulk viscosity */
							postype rvel = Vec2DDotP(&er, &ui);
							if(rvel<0){
								postype chi = 2*rvel/(dramp/2);
								pi = pi - Zeta*chi;
							}
						}
						else {
							postype rvel = Vec2DDotP(&er, &ui);
							postype mu = rvel/(etavis/dramp + eps2vis*dramp/etavis);
							if(rvel<0){
								postype meanden = 0.5*(ibp_den + jbp_den);
								postype meanCsound = 0.5*(ibp_csound + jbp_csound);
								pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
							}
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
            free(p);free(neighbors); free(neighwork);

        }
    	free(vorocorner);
    }
	/*
	if(Dtime > 1.E9){
		Dtime = bp[0].dt;
	}
	*/
	if(0){
		int i;
		for(i=0;i<VORO_NP(simpar);i++){
			if(PINDX(VORORK4_TBP(simpar)+i) == 16799){
				DEBUGPRINT("P%d has pindx= %ld %g %g\n", MYID(simpar), PINDX(VORORK4_TBP(simpar)+i),
						(VORORK4_TBP(simpar)+i)->x,
						(VORORK4_TBP(simpar)+i)->y);
			}
		}
	}

	free(VORORK4_TBPP(simpar));

	postype TDtime;
	MPI_Reduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MAX, 0, MPI_COMM(simpar));
	if(MYID(simpar) == 0) Dtime = TDtime;
	MPI_Bcast(&Dtime, 1, MPI_POSTYPE,  0, MPI_COMM(simpar));
	
	int i;
	postype dmax=0;
	postype eps2 = GL2D_EPS(simpar)*GL2D_EPS(simpar);
#ifdef _OPENMP
#pragma omp parallel for reduction(max: dmax)
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(bp[i].u4if.Flag[ENDIAN_OFFSET] >0)
		{
			bp[i].x += (bp[i].vx + 0.5*bp[i].ax*Dtime)*Dtime; 
			bp[i].y += (bp[i].vy + 0.5*bp[i].ay*Dtime)*Dtime; 
			dmax = MAX(dmax, fabs(bp[i].ay)); 
			bp[i].x = fmod(bp[i].x+Lx,Lx); 
			bp[i].y = fmod(bp[i].y+Ly,Ly); 
			bp[i].vx += bp[i].ax *Dtime; 
			bp[i].vy += bp[i].ay *Dtime;
			/*
			bp[i].vx = 0;
			bp[i].vy = 0;
			*/
		}
    }
	DEBUGPRINT("P%d has maximum displacement in y %g\n", MYID(simpar), dmax*Dtime*Dtime*0.5);
	/* migrate particles between mpi ranks */
	migrateTreeVoroParticles(simpar);
	/* Update the linked list */
	gl2d_MkLinkedList(simpar);


	/* repositioning the memory space */
	bp = VORORK4_TBP(simpar);

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
            treevorork4particletype *p = gl2d_Voro2D_FindCellBP(simpar,ix,iy,&np);
            int nneigh;
            Voro2D_point *neighbors = gl2d_Voro2D_FindNeighbor(simpar,ix,iy,&nneigh);
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
				ibp->w2 = ibp->w2ceil = center.w2ceil;
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
                        postype pi = Voro_Pressure2D(
                                ibp_pressure,jbp_pressure, tmp, neighwork, simpar,GL2D);

                        Voro2D_point dr = EunhaVec2DSub(jbp,ibp);

                        /*
                        if(dr.x > 0.5*Lx) dr.x = dr.x-Lx;
                        else if(dr.x < -0.5*Lx) dr.x = Lx+dr.x;
                        if(dr.y > 0.5*Ly) dr.y = dr.y-Ly;
                        else if(dr.y < -0.5*Ly) dr.y = Ly+dr.y;
                        */

                        postype dramp = sqrt(Vec2DDotP(&dr, &dr));
                        Voro2D_point er;
                        er.x = dr.x/dramp;
                        er.y = dr.y/dramp;

                        Voro2D_point ui,ua,ub;
                        ui.x = Half*(jbp_vx - ibp_vx);
                        ui.y = Half*(jbp_vy - ibp_vy);
                        /*
                        if(PINDX(ibp)==0)
                            DEBUGPRINT("-P%d has indx=0 with pi= %g ij pressure= %g %g jid= %ld\n",
                                MYID(simpar), pi, ibp_pressure, jbp_pressure, PINDX(jbp));
                                */

                        if(0){
							/* for the bulk viscosity */
                            postype rvel = Vec2DDotP(&er, &ui);
                            if(rvel<0){
                                postype chi = 2*rvel/(dramp/2);
                                pi = pi - Zeta*chi;
                            }
                        }
                        else {
                            postype rvel = Vec2DDotP(&er, &ui);
                            postype mu = rvel/(etavis/dramp + eps2vis*dramp/etavis);
                            if(rvel<0){
                                postype meanden = 0.5*(ibp_den + jbp_den);
                                postype meanCsound = 0.5*(ibp_csound + jbp_csound);
                                pi = pi +(-alphavis * meanCsound*mu + betavis*mu*mu)*meanden;
                            }
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
            }
            free(p);free(neighbors); free(neighwork);

        }
    	free(vorocorner);
    }
	DEBUGPRINT("P%d is now updating the internal energy\n", MYID(simpar));
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		if(bp[i].u4if.Flag[ENDIAN_OFFSET] >0){
			if(0){ 
				bp[i].ke  = Half*bp[i].mass*(bp[i].vx*bp[i].vx + bp[i].vy*bp[i].vy); 
				bp[i].ie += bp[i].die * Dtime;
				bp[i].te = bp[i].ie + bp[i].ke;
				bp[i].ie = MAX(0., bp[i].ie);
			} 
			if(1){
				bp[i].ke  = Half*bp[i].mass*(bp[i].vx*bp[i].vx + bp[i].vy*bp[i].vy); 
				bp[i].ie += bp[i].die * Dtime; 
				bp[i].te = bp[i].ie + bp[i].ke; 
				bp[i].w2 = 0.1*(bp[i].pressure * Lx/NX(simpar));
				bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
			}
			if(bp[i].volume ==0) 
				DEBUGPRINT("P%d has wrong volume at i= %d with den= %g at xy= %g %g\n", MYID(simpar), 
						i,bp[i].den, bp[i].x, bp[i].y);
	        bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1); 
			bp[i].den  = bp[i].mass/bp[i].volume;
			bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
			if(isnan(bp[i].csound)){
				DEBUGPRINT("P%d has nan at the sound speed = %d with id= %ld: pressure= %g den= %g Delta IE*Dt= %g\n",
						MYID(simpar), i,PINDX(bp+i), bp[i].pressure, bp[i].den, bp[i].ie-bp[i].die*Dtime);
				exit(9);
			}
		}

    }
	migrateTreeVoroParticles(simpar);

	free(VORORK4_TBPP(simpar));
	MPI_Barrier(MPI_COMM_WORLD);
	DEBUGPRINT("+P%d has np = %ld\n", MYID(simpar), VORO_NP(simpar));

    return Dtime;
}

