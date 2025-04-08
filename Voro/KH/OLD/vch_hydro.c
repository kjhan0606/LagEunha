#include "voro.h"
#include "kh.h"
#include<omp.h>

void initfindmass(Voro2D_GasParticle *bp, int nbp){
	int iy;
#pragma omp parallel for 
    for(iy=0;iy<Ny;iy++){
        int mp=1000;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        ptype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=0;ix<Nx;ix++){
            int np;
            Voro2D_GasParticle *p = Voro2D_FindCellBP(ix,iy,&np,bp);
            int nneigh;
            Voro2D_point *neighbors = Voro2D_FindNeighbor(ix,iy,&nneigh, bp);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.id = p[i].id;

                int id = p[i].id;
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);


                bp[id].volume = Area2DPolygon(vorocorner, mp);
				bp[id].mass = bp[id].den*bp[id].volume;
			}
			free(p);free(neighbors); free(neighwork);
		}
		free(vorocorner);
	
	}
}

double vph2D(Voro2D_GasParticle *bp, int nbp){
	int isave = -1;
    ptype Dtime = 1.e10;
	if(0){
		int i;
		for(i=0;i<nbp;i++) bp[i].dt = 1.e10;
	}
	MkLinkedList(bp,nbp);
	printf("passed mklinkedlist\n");

    int iy;
#pragma omp parallel for reduction(min:Dtime)
    for(iy=0;iy<Ny;iy++){
    	int mp=1000;
    	Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
    	ptype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
        for(ix=0;ix<Nx;ix++){
            int np;
            Voro2D_GasParticle *p = Voro2D_FindCellBP(ix,iy,&np,bp);
            int nneigh;
            Voro2D_point *neighbors = Voro2D_FindNeighbor(ix,iy,&nneigh, bp);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.id = p[i].id;
                center.csound = p[i].csound;

				bp[p[i].id].dt = 1.e10; // initialization of dt 


                int id = p[i].id;
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);


                bp[id].volume = Area2DPolygon(vorocorner, mp);


				Voro2D_Corner *tmp,*tmp2;
				tmp = vorocorner;
				int j;
				double die, dte,dke,fx,fy;
				die = dte = dke = fx = fy = 0;
				do {
					if(tmp->upperrelated >=0){
						int jid = neighwork[tmp->upperrelated].id;
						Voro2D_point line; 
	                    tmp2 = tmp->upperlink; 
						line.x = tmp2->x - tmp->x; 
						line.y = tmp2->y - tmp->y;
	                    Voro2D_point dS = voro2D_norm(&line); 
						ptype facearea = Vec2DLength(tmp,tmp2); 
						dS.x = facearea*dS.x;
						dS.y = facearea*dS.y;
						ptype pi = Half*(bp[jid].pressure + bp[id].pressure); 

						Voro2D_point dr = Vec2DSub(bp+jid,bp+id);

						if(dr.x > 0.5*Lx) dr.x = dr.x-Lx;
						else if(dr.x < -0.5*Lx) dr.x = Lx+dr.x;
						if(dr.y > 0.5*Ly) dr.y = dr.y-Ly;
						else if(dr.y < -0.5*Ly) dr.y = Ly+dr.y;

						ptype dramp = sqrt(Vec2DDotP(&dr, &dr));
						Voro2D_point er;
						er.x = dr.x/dramp;
						er.y = dr.y/dramp;

						Voro2D_point ui,ua,ub;
						ui.x = Half*(bp[jid].vx - bp[id].vx); 
						ui.y = Half*(bp[jid].vy - bp[id].vy);

						if(0){
							/* for the bulk viscosity */
							ptype rvel = Vec2DDotP(&er, &ui);
							if(rvel<0){
								ptype chi = 2*rvel/(dramp/2);
								pi = pi - Zeta*chi;
							}
						}
						else {
							ptype rvel = Vec2DDotP(&er, &ui);
							if(rvel<0){
								ptype meanden = 0.5*(bp[id].den + bp[jid].den);
								ptype meanCsound = 0.5*(bp[id].csound + bp[jid].csound);
								pi = pi +(-alphavis * meanCsound*rvel + betavis*rvel*rvel)*meanden;

							}
						}

						/* for the internal energy */
						die += -pi * Vec2DDotP(&ui,&dS);

						/* for the total energy */
						ua.x = Half*(bp[jid].vx + bp[id].vx); 
						ua.y = Half*(bp[jid].vy + bp[id].vy);
						dte += -pi*Vec2DDotP(&ua, &dS);

						/* for the kinetic energy */
						ub.x = bp[id].vx;
						ub.y = bp[id].vy;
						dke += -pi*Vec2DDotP(&ub, &dS);

						/* for the force */
						fx += -pi * dS.x;
						fy += -pi * dS.y;


						Voro2D_point dv;
						dv.x = (bp[jid].vx - bp[id].vx); 
						dv.y = (bp[jid].vy - bp[id].vy); 

						ptype VdotR = Vec2DDotP(&dv,&er);
						ptype vsig = (bp[jid].csound + bp[id].csound - MIN(0, VdotR));
						ptype dt = 2*Courant*dramp/vsig;
						bp[id].dt = MIN(bp[id].dt,dt);
						if(dt < Dtime){
							Dtime = dt;
							isave = id;
						}
					}
					tmp = tmp->upperlink;
				} while( tmp != vorocorner);
				bp[id].die = die;
				bp[id].dte = dte;
				bp[id].dke = dke;
				bp[id].ax = fx/bp[id].mass;
				bp[id].ay = fy/bp[id].mass;

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
	
	int i;
#pragma omp parallel for 
    for(i=0;i<nbp;i++){
        bp[i].x += (bp[i].vx + 0.5*bp[i].ax*Dtime)*Dtime;
        bp[i].y += (bp[i].vy + 0.5*bp[i].ay*Dtime)*Dtime;
		bp[i].x = fmod(bp[i].x+Lx,Lx);
		bp[i].y = fmod(bp[i].y+Ly,Ly);

        bp[i].vx += bp[i].ax *Dtime;
        bp[i].vy += bp[i].ay *Dtime;

		/*
        bp[i].ke  = Half*bp[i].mass*(bp[i].vx*bp[i].vx + bp[i].vy*bp[i].vy);

		bp[i].ie += bp[i].die * Dtime;
        bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1);

		bp[i].te = bp[i].ie + bp[i].ke;

        bp[i].den  = bp[i].mass/bp[i].volume;
		bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
        bp[i].poverrhogam = bp[i].pressure / pow(bp[i].den, Gamma);
		*/

    }
	/* Update the linked list */
	MkLinkedList(bp,nbp);

#pragma omp parallel for 
    for(iy=0;iy<Ny;iy++){
    	int mp=1000;
    	Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
    	ptype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
		int ix;
        for(ix=0;ix<Nx;ix++){
            int np;
            Voro2D_GasParticle *p = Voro2D_FindCellBP(ix,iy,&np,bp);
            int nneigh;
            Voro2D_point *neighbors = Voro2D_FindNeighbor(ix,iy,&nneigh, bp);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
			int i;
            for(i=0;i<np;i++){
                Voro2D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.id = p[i].id;

                int id = p[i].id;
                int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize);


                bp[id].volume = Area2DPolygon(vorocorner, mp);


				Voro2D_Corner *tmp,*tmp2;
				tmp = vorocorner;
				int j;
				double die, dte,dke,fx,fy;
				die = dte = dke = fx = fy = 0;
				do {
					if(tmp->upperrelated >=0){
						int jid = neighwork[tmp->upperrelated].id;
						Voro2D_point line; 
	                    tmp2 = tmp->upperlink; 
						line.x = tmp2->x - tmp->x; 
						line.y = tmp2->y - tmp->y;
	                    Voro2D_point dS = voro2D_norm(&line); 
						ptype facearea = Vec2DLength(tmp,tmp2); 
						dS.x = facearea*dS.x;
						dS.y = facearea*dS.y;
						ptype pi = Half*(bp[jid].pressure + bp[id].pressure); 

						Voro2D_point ui,ua,ub;
						ui.x = Half*(bp[jid].vx - bp[id].vx); 
						ui.y = Half*(bp[jid].vy - bp[id].vy);

						Voro2D_point dr = Vec2DSub(bp+jid,bp+id);
						if(dr.x > 0.5*Lx) dr.x = dr.x-Lx;
						else if(dr.x < -0.5*Lx) dr.x = Lx+dr.x;
						if(dr.y > 0.5*Ly) dr.y = dr.y-Ly;
						else if(dr.y < -0.5*Ly) dr.y = Ly+dr.y;

						ptype dramp = sqrt(Vec2DDotP(&dr, &dr));
						Voro2D_point er;
						er.x = dr.x/dramp;
						er.y = dr.y/dramp;

						if(0){
							/* for the bulk viscosity */
							ptype rvel = Vec2DDotP(&er, &ui);
							if(rvel<0){
								ptype chi = 2*rvel/(dramp/2);
								pi = pi - Zeta*chi;
							}
						}
						else {
							ptype rvel = Vec2DDotP(&er, &ui);
							if(rvel<0){
								ptype meanden = 0.5*(bp[id].den + bp[jid].den);
								ptype meanCsound = 0.5*(bp[id].csound + bp[jid].csound);
								pi = pi +(-alphavis * meanCsound*rvel + betavis*rvel*rvel)*meanden;

							}
						}


						/* for the internal energy */
						die += -pi * Vec2DDotP(&ui,&dS);


						/* for the total energy */
						ua.x = Half*(bp[jid].vx + bp[id].vx); 
						ua.y = Half*(bp[jid].vy + bp[id].vy);
						dte += -pi*Vec2DDotP(&ua, &dS);

						/* for the kinetic energy */
						ub.x = bp[id].vx;
						ub.y = bp[id].vy;
						dke += -pi*Vec2DDotP(&ub, &dS);

						/* for the force */
						fx += -pi * dS.x;
						fy += -pi * dS.y;


						/*
						Voro2D_point dv;
						dv.x = (bp[jid].vx - bp[id].vx); 
						dv.y = (bp[jid].vy - bp[id].vy); 

						ptype VdotR = Vec2DDotP(&dv,&dr);
						ptype vsig = (bp[jid].csound + bp[id].csound - MIN(0, VdotR));
						ptype dt = 2*Courant*dramp/vsig;
						bp[id].dt = MIN(bp[id].dt,dt);
						if(dt < Dtime){
							Dtime = dt;
							isave = id;
						}
						*/
					}
					tmp = tmp->upperlink;
				} while( tmp != vorocorner);
				bp[id].die = die;
				bp[id].dte = dte;
				bp[id].dke = dke;

				bp[id].ax = fx/bp[id].mass;
				bp[id].ay = fy/bp[id].mass;

            }
            free(p);free(neighbors); free(neighwork);

        }
    	free(vorocorner);
    }
#pragma omp parallel for 
    for(i=0;i<nbp;i++){

		if(1){
	        bp[i].ke  = Half*bp[i].mass*(bp[i].vx*bp[i].vx + bp[i].vy*bp[i].vy);
			bp[i].ie += bp[i].die * Dtime;
			bp[i].te = bp[i].ie + bp[i].ke;
		}
		else {
			ptype dte = bp[i].dte*Dtime;
			ptype dke = Half*(bp[i].vx*bp[i].vx + bp[i].vy*bp[i].vy)*bp[i].mass - bp[i].ke;
			ptype die = dte - dke;
			bp[i].ie += die;
			bp[i].ke  = Half*(bp[i].vx*bp[i].vx + bp[i].vy*bp[i].vy)*bp[i].mass;
		}

        bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1);
        bp[i].den  = bp[i].mass/bp[i].volume;
		bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
		/*
        bp[i].poverrhogam = bp[i].pressure / pow(bp[i].den, Gamma);
        bp[i].x += bp[i].vx *Dtime + 0.5*bp[i].ax*Dtime*Dtime;
        bp[i].y += bp[i].vy *Dtime + 0.5*bp[i].ay*Dtime*Dtime;
		bp[i].x = FMOD(bp[i].x+Lx,Lx);
		bp[i].y = FMOD(bp[i].y+Ly,Ly);
        bp[i].vx += bp[i].ax *Dtime;
        bp[i].vy += bp[i].ay *Dtime;
		*/

    }
	if(0){
		int id = Nxp*Nxp/4-1; 
		printf("bp id= %d: xy= %g %g den= %g\n", id, bp[id].x, bp[id].y, bp[id].den);
		id = id + Nxp-1;
		printf("bp id= %d: xy= %g %g den= %g\n", id, bp[id].x, bp[id].y, bp[id].den);
		id = Nxp*Nxp/4;
		printf("bp id= %d: xy= %g %g den= %g\n", id, bp[id].x, bp[id].y, bp[id].den);
		id = id + Nxp-1;
		printf("bp id= %d: xy= %g %g den= %g\n", id, bp[id].x, bp[id].y, bp[id].den);
		id = Nxp*Nxp/4 + 1;
		printf("bp id= %d: xy= %g %g den= %g\n", id, bp[id].x, bp[id].y, bp[id].den);
		id = id + Nxp-1;
		printf("bp id= %d: xy= %g %g den= %g\n", id, bp[id].x, bp[id].y, bp[id].den);
	}


    return Dtime;
}

