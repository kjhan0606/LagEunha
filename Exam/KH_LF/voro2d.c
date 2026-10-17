#include "eunha.h"
#include "voro.h"
#include "kh.h"
#include "nnost.h"
#include "gnnost.h"
#include "exam.h"
#include "exam2d.h"
#include "cosmology.h"
#include<omp.h>

double kh_w2Measure(SimParameters *, double , double ,double );

void initfindmass(SimParameters *simpar){
	treevoroparticletype *bp = VORO_TBP(simpar);
	int nbp = VORO_NP(simpar);
	postype xmin,ymin,zmin,xmax,ymax,zmax, pwidth;
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	pwidth = BASICCELL_CELLWIDTH(simpar);
	xmin = KH_XMIN(simpar)-pwidth;
    ymin = KH_YMIN(simpar)-pwidth;
    xmax = KH_XMAX(simpar)+pwidth;
    ymax = KH_YMAX(simpar)+pwidth;

    int mx = BASICCELL_MX(simpar);
    int my = BASICCELL_MY(simpar);


	/*
	DEBUGPRINT("P%d has %g : %g %g : %g %g : %d %d\n", 
			MYID(simpar), pwidth, xmin,ymin,xmax,ymax,
			mx,my);
			*/

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
            treevoroparticletype *p = Voro2D_FindCellBP(simpar,ix,iy,&np,bp);
            int nneigh;
            Voro2D_point *neighbors = Voro2D_FindNeighbor(simpar,ix,iy,&nneigh, bp);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.indx = p[i].u4if.indx;
                center.w2 = p[i].w2;
				/*
				DEBUGPRINT("P%d has w2  %g %g\n", MYID(simpar), center.w2, neighbors[0].w2);
				exit(99);
				if(p[i].u4if.indx == 0){
					printf("p0 %g %g %ld\n", p[i].x, p[i].y, PINDX(p+i));
					int k;
					for(k=0;k<nneigh;k++){
						printf("pj %g %g %d\n", neighbors[k].x, neighbors[k].y, neighbors[k].indx);
					}
				}
				*/

                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
				treevoroparticletype *ibp = p[i].bp;
				ibp->volume = Area2DPolygon(vorocorner,mp);
				ibp->mass = ibp->den * ibp->volume;
			}
			free(neighwork);
			free(neighbors); 
			free(p);
		}
		free(vorocorner);
	
	}
}

double vph2D(SimParameters *simpar){
	postype boxsize = BOXSIZE(simpar)/NX(simpar)*5;
	postype Courant = GAS_COURANT(simpar);
	treevoroparticletype *bp = VORO_TBP(simpar);
	int nbp = VORO_NP(simpar);
	int isave = -1;
    postype Dtime = 1.e10;




	// MkLinkedList will reallocate VORO_CELL and make space for VORO_TBPP(simpar)
	MkLinkedList(simpar);

	postype Lx = KH_SIMBOX(simpar).x.max;
	postype Ly = KH_SIMBOX(simpar).y.max;
	postype xmin,ymin,zmin,xmax,ymax,zmax, pwidth;
	pwidth = BASICCELL_CELLWIDTH(simpar);
	xmin = KH_XMIN(simpar)-pwidth;
    ymin = KH_YMIN(simpar)-pwidth;
    xmax = KH_XMAX(simpar)+pwidth;
    ymax = KH_YMAX(simpar)+pwidth;
    int mx = BASICCELL_MX(simpar);
    int my = BASICCELL_MY(simpar);
	postype Gamma = GAS_GAMMA(simpar);

	float alphavis = GAS_AlphaVis(simpar);
    float betavis = GAS_BetaVis(simpar);
    float etavis = GAS_ETAVIS(simpar) * Lx/NX(simpar);
    float eps2vis = GAS_EPSVIS(simpar)*GAS_EPSVIS(simpar);


	postype dtold = GAS_dtold(simpar);


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
            treevoroparticletype *p = Voro2D_FindCellBP(simpar,ix,iy,&np,bp);
            int nneigh;
            Voro2D_point *neighbors = Voro2D_FindNeighbor(simpar,ix,iy,&nneigh, bp);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
//                center.indx = p[i].u4if.indx;
                center.indx = PINDX(p+i);
                center.csound = p[i].csound;
                center.w2 = p[i].w2;

				treevoroparticletype *ibp = p[i].bp;
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
						treevoroparticletype *jbp =  (treevoroparticletype*)(neighwork[tmp->upperrelated].bp);
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
								ibp_pressure,jbp_pressure, tmp, neighwork, simpar, KH
								);
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
						ui = get2dUpqrad(ibp,jbp,dtold);

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
						ibp->dt = MIN(ibp->dt,dt);
						if(dt < Dtime){
							Dtime = dt;
							// isave = id;
						}
						if(isnan(fx)){
							DEBUGPRINT("P%d has nan pi= %g dS.x/y= %g %g ui.x/y= %g %g pi= %g pj= %g dtold= %g\n",
									MYID(simpar), pi,dS.x, dS.y, ui.x, ui.y, ibp_pressure, jbp_pressure, dtold);
							exit(0);
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

//	DEBUGPRINT("P%d is exiting here\n", MYID(simpar)); exit(9);

//	exit(90);
	/*
	if(Dtime > 1.E9){
		Dtime = bp[0].dt;
	}
	*/
	free(VORO_TBPP(simpar));
	postype TDtime;
	MPI_Reduce(&Dtime, &TDtime, 1, MPI_POSTYPE, MPI_MAX, 0, MPI_COMM(simpar));
	if(MYID(simpar) == 0) Dtime = TDtime;
	MPI_Bcast(&Dtime, 1, MPI_POSTYPE,  0, MPI_COMM(simpar));
	
	int i;


#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
        bp[i].vx += bp[i].ax *Dtime;
        bp[i].vy += bp[i].ay *Dtime;
//		bp[i].ie += bp[i].die * Dtime;

        bp[i].x += bp[i].vx *Dtime;
        bp[i].y += bp[i].vy *Dtime;
		bp[i].x = fmod(bp[i].x+Lx,Lx);
		bp[i].y = fmod(bp[i].y+Ly,Ly);

		bp[i].w2old = bp[i].w2;

    }
	/* migrate particles between mpi ranks */
	KH_TreeAllParticleMigrate(simpar);

	// update w2ceil and w2
	if(GAS_Kappa(simpar)>0) det2d_dpq(simpar,paddingTreeVoroParticles);



	// Update the linked list 
	MkLinkedList(simpar);

	bp = VORO_TBP(simpar);

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
            treevoroparticletype *p = Voro2D_FindCellBP(simpar,ix,iy,&np,bp);
            int nneigh;
            Voro2D_point *neighbors = Voro2D_FindNeighbor(simpar,ix,iy,&nneigh, bp);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
//                center.indx = p[i].u4if.indx;
                center.indx = PINDX(p+i);
                center.w2 = p[i].w2;


				treevoroparticletype *ibp = p[i].bp;
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);
				// These are to use register 
				postype ibp_vx,ibp_vy,ibp_csound,ibp_pressure,ibp_den;
				ibp_vx = ibp->vx;
				ibp_vy = ibp->vy;
				ibp_csound = ibp->csound;
				ibp_pressure = ibp->pressure;
				ibp_den = ibp->den;


                ibp->volume = Area2DPolygon(vorocorner, mp);


//				/*
				Voro2D_Corner *tmp,*tmp2;
				tmp = vorocorner;
				int j;
				double die, dte,dke,fx,fy;
				die = dte = dke = fx = fy = 0;
				do {
					if(tmp->upperrelated >=0){
//						int jid = neighwork[tmp->upperrelated].indx;
						// These are to use register 
						treevoroparticletype *jbp =  (treevoroparticletype*)(neighwork[tmp->upperrelated].bp);
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
								ibp_pressure,jbp_pressure, tmp, neighwork, simpar, KH
								);

						Voro2D_point ui,ua,ub;
//						ui.x = Half*(jbp_vx - ibp_vx); 
//						ui.y = Half*(jbp_vy - ibp_vy);
						ui = get2dUpqrad(ibp,jbp,dtold);

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


					}
					tmp = tmp->upperlink;
				} while( tmp != vorocorner);
				ibp->die = die;

				ibp->ax = fx/ibp->mass;
				ibp->ay = fy/ibp->mass;
//				*/

            }
            free(p);free(neighbors); free(neighwork);

        }
    	free(vorocorner);
    }
//	DEBUGPRINT("P%d is before passing free(voro_tbpp): %p\n", MYID(simpar), VORO_TBPP(simpar));
	free(VORO_TBPP(simpar));

	postype dmean = KH_SIMBOX(simpar).x.max/NX(simpar);

#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<VORO_NP(simpar);i++){
		bp[i].ie += bp[i].die * Dtime;
        bp[i].den  = bp[i].mass/bp[i].volume;
   	    bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1);
		bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
//		bp[i].w2old = bp[i].w2;
		if(GAS_Kappa(simpar) <0) {
			bp[i].w2 = - GAS_Kappa(simpar);
		}
		else if(GAS_Kappa(simpar)>0) {
			bp[i].w2 = kh_w2Measure(simpar,dmean, bp[i].pressure, bp[i].den);
        	bp[i].w2 = MIN(bp[i].w2, bp[i].w2ceil);
		}

    }
	postype minpressure=1.e23;
	for(i=0;i<VORO_NP(simpar);i++){
		minpressure = MIN(minpressure, bp[i].pressure);
	}
	postype tminpressure;
	MPI_Reduce(&minpressure, &tminpressure, 1, MPI_POSTYPE, MPI_MIN, 0, MPI_COMM_WORLD);
	if(MYID(simpar)==0) DEBUGPRINT("The minimum pressure is %g\n", tminpressure);


//	DEBUGPRINT("+P%d has np = %ld\n", MYID(simpar), VORO_NP(simpar));

    return Dtime;
}

