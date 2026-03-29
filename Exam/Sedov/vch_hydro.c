#include "voro.h"
#include "sedov.h"
#ifdef _OPENMP
#include<omp.h>
#endif


double vph3D(Voro3D_GasParticle **Bp, int *mbp, ptype avgVolume){
	int isave = -1;
    ptype Dtime = 1.e10;

	Voro3D_GasParticle *bp = *Bp;
	int nbp = *mbp;


	MkLinkedList(bp,nbp);

	int nthreads=1;
	int mp=1024; 
	Voro3D_Vertex *tvorovertex;
	Voro3D_point *tneighbors,*tneighwork;
	Voro3D_GasParticle *tp;
	int *tintwork;
	int npmax = get_npmax();

#ifdef _OPENMP
#pragma omp parallel 
#endif
	{
#ifdef _OPENMP
#pragma omp master
#endif
		{
			nthreads = omp_get_num_threads();
		}
	}

	tvorovertex = (Voro3D_Vertex*)malloc(sizeof(Voro3D_Vertex)*mp*nthreads); 
	tneighbors = (Voro3D_point*)malloc(sizeof(Voro3D_point)*npmax*27*nthreads); 
	tneighwork = (Voro3D_point*)malloc(sizeof(Voro3D_point)*npmax*27*nthreads); 
	tp = (Voro3D_GasParticle*)malloc(sizeof(Voro3D_GasParticle)*npmax*nthreads); 
	tintwork = (int*)malloc(sizeof(int)*mp*nthreads); 

#ifdef _OPENMP
#pragma omp parallel 
#endif
	{
		int ithread = omp_get_thread_num();
		Voro3D_Vertex *vorovertex = tvorovertex+ithread*mp;
		Voro3D_point *neighbors = tneighbors+ithread*27*npmax;
		Voro3D_point *neighwork = tneighwork+ithread*27*npmax;
		Voro3D_GasParticle *p = tp+ithread*npmax;
		int *intwork = tintwork+ithread*mp;
		int ixyz;
#ifdef _OPENMP
#pragma omp for reduction(min:Dtime) schedule(dynamic)
#endif
		for(ixyz=0;ixyz<Nx*Ny*Nz;ixyz++){
			int iz = ixyz/(Nx*Ny);
			int iy = ixyz%(Nx*Ny)/Nx;
			int ix = ixyz%(Nx*Ny)%Nx;
			ptype dlx,dly,dl,dvx,dvy,dv,ax,ay,a; 

			int np; 
			Voro3D_FindCellBP(p,ix,iy,iz,&np,bp); 
			int nneigh;
            Voro3D_FindNeighbor(neighbors,ix,iy,iz,&nneigh, bp); 
			int i; 
			for(i=0;i<np;i++){ 
				Voro3D_point center; 
				center.x = p[i].x; 
				center.y = p[i].y; 
				center.z = p[i].z; 
				center.id = p[i].id;
				center.csound = p[i].csound;

                int id = p[i].id; 
				int ishrink = 0;
				int ip = Voro3D_FindVC(&center,neighbors, neighwork,nneigh, vorovertex,mp,boxsize, 
						 ishrink, intwork);


				bp[id].volume = Voro3D_Volume_Polyhedron(vorovertex,ip);


				Voro3D_Vertex *tmp,*tmp2;
				tmp = vorovertex;
				int iv,j;
				double die, dte,dke,fx,fy,fz;
				for(iv=0;iv<ip;iv++){
					vorovertex[iv].considered[0] = vorovertex[iv].considered[1] = vorovertex[iv].considered[2] = No;
				}
				die = dte = dke = fx = fy = fz = 0;
				for(iv=0;iv<ip;iv++){
					Voro3D_Vertex *start = vorovertex+iv;
					for(j=0;j<3;j++){
						if(start->considered[j]==No){
							int jid = gIDMakingFace3D(neighwork,start, j); 
							Voro3D_Vertex *link = start->link[j];
							Voro3D_point dS = Voro3D_norm_polygon(start, link);
							ptype pij = Half*(bp[jid].pressure + bp[id].pressure);
							Voro3D_point uij,ua,ub; 
							uij.x = Half*(bp[jid].vx - bp[id].vx); 
							uij.y = Half*(bp[jid].vy - bp[id].vy); 
							uij.z = Half*(bp[jid].vz - bp[id].vz);
							Voro3D_point dr = Vec3DSub(bp+jid,bp+id);
                            if(dr.x > 0.5*Lx) dr.x = dr.x-Lx; 
							else if(dr.x < -0.5*Lx) dr.x = Lx+dr.x; 
							if(dr.y > 0.5*Ly) dr.y = dr.y-Ly; 
							else if(dr.y < -0.5*Ly) dr.y = Ly+dr.y; 
							if(dr.z > 0.5*Lz) dr.z = dr.z-Lz; 
							else if(dr.z < -0.5*Lz) dr.z = Lz+dr.z; 

                            ptype dramp = sqrt(Vec3DDotP(&dr, &dr)); 
							Voro3D_point er; 
							er.x = dr.x/dramp; 
							er.y = dr.y/dramp; 
							er.z = dr.z/dramp;
							ptype rvel;
							if(0){
								/* for the bulk viscosity */
								rvel = Vec3DDotP(&er, &uij);
								if(rvel<0){
									ptype chi = 2*rvel/(dramp/2);
									pij = pij - Zeta*chi;
								}
							}
							else {
								rvel = Vec3DDotP(&er, &uij);
								if(rvel<0){
									ptype meanden = 0.5*(bp[id].den + bp[jid].den);
									ptype meanCsound = 0.5*(bp[id].csound + bp[jid].csound);
									pij = pij +(-alphavis * meanCsound*rvel + betavis*rvel*rvel)*meanden;
	
								}
							}

							/* for the internal energy */
							die += -pij * Vec3DDotP(&uij,&dS);

							/* for the total energy */
							ua.x = Half*(bp[jid].vx + bp[id].vx); 
							ua.y = Half*(bp[jid].vy + bp[id].vy);
							ua.z = Half*(bp[jid].vz + bp[id].vz);
							dte += -pij*Vec3DDotP(&ua, &dS);
	
							/* for the kinetic energy */
							ub.x = bp[id].vx;
							ub.y = bp[id].vy;
							ub.z = bp[id].vz;
							dke += -pij*Vec3DDotP(&ub, &dS);
	
							/* for the force */
							fx += -pij * dS.x;
							fy += -pij * dS.y;
							fz += -pij * dS.z;
	
	
							Voro3D_point dv;
							dv.x = (bp[jid].vx - bp[id].vx); 
							dv.y = (bp[jid].vy - bp[id].vy); 
							dv.z = (bp[jid].vz - bp[id].vz); 
	
							ptype VdotR = Vec3DDotP(&dv,&dr);
							ptype vsig = (bp[jid].csound + bp[id].csound - MIN(0, VdotR));
							ptype accel = sqrt(fx*fx + fy*fy + fz*fz)/bp[id].mass;
							ptype dt1 = 2*Courant*dramp/vsig;
							ptype dt2 = sqrt(0.05*2*dramp/accel);
							ptype dt3 = (0.1*dramp/Vec3DDotP(&uij,&uij));
							if(rvel >0) {
								dt1 = dt2 = 1.e20;
							}
							ptype dt = MIN(dt1,dt2);
							dt = MIN(dt,dt3);
							bp[id].dt = MIN(bp[id].dt,dt);
							if(dt < Dtime){
								Dtime = dt;
								isave = id;
							}
						}
					}
				}
				bp[id].die = die;
				bp[id].dte = dte;
				bp[id].dke = dke;
				bp[id].ax = fx/bp[id].mass;
				bp[id].ay = fy/bp[id].mass; 
				bp[id].az = fz/bp[id].mass; 
			}
		} 
	}
	/*
	free(tvorovertex);
 	free(tneighbors);
 	free(tneighwork);
 	free(tp);
	*/

	int i;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<nbp;i++){
		ptype dx,dy,dz;
        bp[i].x += (dx=(bp[i].vx + 0.5*bp[i].ax*Dtime)*Dtime);
        bp[i].y += (dy=(bp[i].vy + 0.5*bp[i].ay*Dtime)*Dtime);
        bp[i].z += (dz=(bp[i].vz + 0.5*bp[i].az*Dtime)*Dtime);
		bp[i].x = fmod(bp[i].x+Lx,Lx);
		bp[i].y = fmod(bp[i].y+Ly,Ly);
		bp[i].z = fmod(bp[i].z+Lz,Lz);

        bp[i].vx += bp[i].ax *Dtime;
        bp[i].vy += bp[i].ay *Dtime;
        bp[i].vz += bp[i].az *Dtime;

		if(bp[i].id == Nxp/2-1+Nxp*(Nyp/2+Nyp*(Nzp/2))){
			printf("neighbor x/y/z= %g %g %g : dx/dy/dz= %g %g %g: %g %g %g\n", bp[i].x, bp[i].y, bp[i].z, dx,dy,dz, bp[i].ax, bp[i].ay, bp[i].az);
		}
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
	{
		int icenter = Nxp/2 + Nxp*(Nyp/2 + Nyp*(Nzp/2));
		bp[icenter].x = Lx/2;
		bp[icenter].y = Ly/2;
		bp[icenter].z = Lz/2;
		bp[icenter].vx = bp[icenter].vy = bp[icenter].vz = 0;
	}
	/* Update the linked list */
	MkLinkedList(bp,nbp);

	npmax = get_npmax();


	tvorovertex = (Voro3D_Vertex*)realloc(tvorovertex,sizeof(Voro3D_Vertex)*mp*nthreads); 
	tneighbors = (Voro3D_point*)realloc(tneighbors,sizeof(Voro3D_point)*npmax*27*nthreads); 
	tneighwork = (Voro3D_point*)realloc(tneighwork,sizeof(Voro3D_point)*npmax*27*nthreads); 
	tp = (Voro3D_GasParticle*)realloc(tp,sizeof(Voro3D_GasParticle)*npmax*nthreads); 
	tintwork = (int*)realloc(tintwork,sizeof(int)*mp*nthreads); 

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		int ithread = omp_get_thread_num();
		Voro3D_Vertex *vorovertex = tvorovertex+ithread*mp;
		Voro3D_point *neighbors = tneighbors+ithread*27*npmax;
		Voro3D_point *neighwork = tneighwork+ithread*27*npmax;
		Voro3D_GasParticle *p = tp+ithread*npmax;
		int *intwork = tintwork+ithread*mp;
		int ixyz;
#ifdef _OPENMP
#pragma omp for  schedule(dynamic)
#endif
		for(ixyz=0;ixyz<Nx*Ny*Nz;ixyz++){
			int iz = ixyz/(Nx*Ny);
			int iy = ixyz%(Nx*Ny)/Nx;
			int ix = ixyz%(Nx*Ny)%Nx;
			ptype dlx,dly,dl,dvx,dvy,dv,ax,ay,a; 
			int np;
            Voro3D_FindCellBP(p,ix,iy,iz,&np,bp); 
			int nneigh; 
			Voro3D_FindNeighbor(neighbors,ix,iy,iz,&nneigh, bp); 
			int i; 
			for(i=0;i<np;i++){ 
				Voro3D_point center;
                center.x = p[i].x; 
				center.y = p[i].y; 
				center.z = p[i].z; 
				center.id = p[i].id;
				center.csound = p[i].csound;

                int id = p[i].id; 
				int ishrink = 0;
				int ip = Voro3D_FindVC(&center,neighbors, neighwork,nneigh, vorovertex,mp,boxsize, 
						ishrink, intwork);

				bp[id].volume = Voro3D_Volume_Polyhedron(vorovertex,ip);

				Voro3D_Vertex *tmp,*tmp2;
				tmp = vorovertex;
				int iv,j;
				double die, dte,dke,fx,fy,fz;
				for(iv=0;iv<ip;iv++){
					vorovertex[iv].considered[0] = vorovertex[iv].considered[1] = vorovertex[iv].considered[2] = No;
				}
				die = dte = dke = fx = fy = fz = 0;
				bp[id].minsize2 = 1.e10;
				for(iv=0;iv<ip;iv++){
					Voro3D_Vertex *start = vorovertex+iv;
					ptype dist2 = vorovertex[iv].x*vorovertex[iv].x + vorovertex[iv].y*vorovertex[iv].y +vorovertex[iv].z*vorovertex[iv].z;
					bp[id].minsize2 = MIN(bp[id].minsize2, dist2);
					for(j=0;j<3;j++){
						if(start->considered[j]==No){
							int jid = gIDMakingFace3D(neighwork,start, j); 
							Voro3D_Vertex *link = start->link[j];
							Voro3D_point dS = Voro3D_norm_polygon(start, link);
							ptype pij = Half*(bp[jid].pressure + bp[id].pressure);
							Voro3D_point uij,ua,ub; 
							uij.x = Half*(bp[jid].vx - bp[id].vx); 
							uij.y = Half*(bp[jid].vy - bp[id].vy); 
							uij.z = Half*(bp[jid].vz - bp[id].vz);
							Voro3D_point dr = Vec3DSub(bp+jid,bp+id);
                            if(dr.x > 0.5*Lx) dr.x = dr.x-Lx; 
							else if(dr.x < -0.5*Lx) dr.x = Lx+dr.x; 
							if(dr.y > 0.5*Ly) dr.y = dr.y-Ly; 
							else if(dr.y < -0.5*Ly) dr.y = Ly+dr.y; 
							if(dr.z > 0.5*Lz) dr.z = dr.z-Lz; 
							else if(dr.z < -0.5*Lz) dr.z = Lz+dr.z; 

                            ptype dramp = sqrt(Vec3DDotP(&dr, &dr)); 
							Voro3D_point er; 
							er.x = dr.x/dramp; 
							er.y = dr.y/dramp; 
							er.z = dr.z/dramp;
							ptype rvel;
							if(0){
								/* for the bulk viscosity */
								rvel = Vec3DDotP(&er, &uij);
								if(rvel<0){
									ptype chi = 2*rvel/(dramp/2);
									pij = pij - Zeta*chi;
								}
							}
							else {
								rvel = Vec3DDotP(&er, &uij);
								if(rvel<0){
									ptype meanden = 0.5*(bp[id].den + bp[jid].den);
									ptype meanCsound = 0.5*(bp[id].csound + bp[jid].csound);
									pij = pij +(-alphavis * meanCsound*rvel + betavis*rvel*rvel)*meanden;
	
								}
							}
							/* for the internal energy */
							die += -pij * Vec3DDotP(&uij,&dS);


							/* for the total energy */
							ua.x = Half*(bp[jid].vx + bp[id].vx); 
							ua.y = Half*(bp[jid].vy + bp[id].vy);
							ua.z = Half*(bp[jid].vz + bp[id].vz);
							dte += -pij*Vec3DDotP(&ua, &dS);
	
							/* for the kinetic energy */
							ub.x = bp[id].vx;
							ub.y = bp[id].vy;
							ub.z = bp[id].vz;
							dke += -pij*Vec3DDotP(&ub, &dS);
	
							/* for the force */
							fx += -pij * dS.x;
							fy += -pij * dS.y;
							fz += -pij * dS.z;
						}
					}
				}


				bp[id].die = die;
				bp[id].dte = dte;
				bp[id].dke = dke;

				bp[id].ax = fx/bp[id].mass;
				bp[id].ay = fy/bp[id].mass; 
				bp[id].az = fz/bp[id].mass; 
			} 
		} 
	} 
	free(tneighbors);
	free(tvorovertex);
	free(tneighwork);
	free(tp);
	free(tintwork);
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(i=0;i<nbp;i++){

		if(1){
	        bp[i].ke  = Half*bp[i].mass*(bp[i].vx*bp[i].vx + bp[i].vy*bp[i].vy +bp[i].vz*bp[i].vz);
			bp[i].ie += bp[i].die * Dtime;
			bp[i].te = bp[i].ie + bp[i].ke;
			bp[i].ie = MAX(0., bp[i].ie);
		}
		else {
			ptype dte = bp[i].dte*Dtime;
			ptype dke = Half*(bp[i].vx*bp[i].vx + bp[i].vy*bp[i].vy + bp[i].vz*bp[i].vz)*bp[i].mass - bp[i].ke;
			ptype die = dte - dke;
			bp[i].ie += die;
			bp[i].ke  = Half*(bp[i].vx*bp[i].vx + bp[i].vy*bp[i].vy)*bp[i].mass;
		}

        bp[i].pressure = bp[i].ie / bp[i].volume * (Gamma-1);
        bp[i].den  = bp[i].mass/bp[i].volume;
		bp[i].csound = sqrt(Gamma*bp[i].pressure/bp[i].den);
		/*
        bp[i].x += bp[i].vx *Dtime + 0.5*bp[i].ax*Dtime*Dtime;
        bp[i].y += bp[i].vy *Dtime + 0.5*bp[i].ay*Dtime*Dtime;
		bp[i].x = FMOD(bp[i].x+Lx,Lx);
		bp[i].y = FMOD(bp[i].y+Ly,Ly);
        bp[i].vx += bp[i].ax *Dtime;
        bp[i].vy += bp[i].ay *Dtime;
		*/

    }
	if(0){
		int icount = 0;
		int ii;
		for(ii=0;ii<nbp;ii++){
			if(bp[ii].volume > 27* avgVolume){
				icount ++;
				int j;
				*Bp = bp = (Voro3D_GasParticle*)realloc(bp, sizeof(Voro3D_GasParticle)*(*mbp +26));
				ptype dx,dy,dz;
				dx = dy = dz = sqrt(bp[ii].minsize2);
				int i = *mbp;

				for(j=0;j<27;j++){
					int ix,iy,iz;
					iz = j/9 -1;
					iy = (j%9)/3-1;
					ix = (j%9)%3-1;

					ptype r = sqrt(iz*iz + iy*iy + ix*ix);
					ptype ratio = 1.L/r;

					if(j == 13) continue;
					bp[i].x = bp[ii].x + dx*ix*ratio;
					bp[i].y = bp[ii].y + dy*iy*ratio;
					bp[i].z = bp[ii].z + dz*iz*ratio;
					bp[i].mass = bp[ii].mass/52.;
					bp[i].ie = bp[ii].ie/52.;
					bp[i].ke = bp[ii].ke/52.;
					bp[i].te = bp[ii].te/52.;
					bp[i].volume = bp[ii].volume/52.;
					bp[i].pressure = bp[ii].pressure;
					bp[i].csound = bp[ii].csound;
					bp[i].den = bp[ii].den;
					bp[i].vx = bp[ii].vx;
					bp[i].vy = bp[ii].vy;
					bp[i].vz = bp[ii].vz;
					bp[i].id = i;
					i ++;
				}
				bp[ii].mass /= 2;
				bp[ii].ke /= 2;
				bp[ii].te /= 2;
				bp[ii].ie /= 2;
				bp[ii].volume /= 2;
				*mbp += 26;
			}

		}
		if(icount >0) printf("Total %d particles are split. Total np= %d\n", icount, *mbp);

	}
	else if(0){
		int icount = 0;
		int ii;
		for(ii=0;ii<nbp;ii++){
			if(bp[ii].volume > 7* avgVolume){
				icount ++;
				int j;
				*Bp = bp = (Voro3D_GasParticle*)realloc(bp, sizeof(Voro3D_GasParticle)*(*mbp +6));
				ptype dx,dy,dz;
				dx = dy = dz = sqrt(bp[ii].minsize2);
				int i = *mbp;

				for(j=0;j<6;j++){
					int ix,iy,iz;
					if(j==0) {
						iz = 0; ix = -1; iy = 0;
					}
					else if(j==1){
						iz = 0; ix =  1; iy = 0;
					}
					else if(j==2){
						iz = 0; ix =  0; iy = -1;
					}
					else if(j==3){
						iz = 0; ix =  0; iy =  1;
					}
					else if(j==4){
						iz = -1; ix =  0; iy =  0;
					}
					else if(j==5){
						iz = 1; ix =  0; iy =  0;
					}
					if(j == 13) continue;
					bp[i] = bp[ii];
					bp[i].x = bp[ii].x + dx*ix;
					bp[i].y = bp[ii].y + dy*iy;
					bp[i].z = bp[ii].z + dz*iz;
					bp[i].mass = bp[ii].mass/7.;
					bp[i].ie = bp[ii].ie/7.;
					bp[i].ke = bp[ii].ke/7.;
					bp[i].te = bp[ii].te/7.;
					bp[i].volume = bp[ii].volume/7.;
					bp[i].id = i;
					i ++;
				}
				bp[ii].mass /= 7;
				bp[ii].ke /= 7;
				bp[ii].te /= 7;
				bp[ii].ie /= 7;
				bp[ii].volume /= 7;
				*mbp += 6;
			}

		}
		if(icount >0) printf("Total %d particles are split. Total np= %d\n", icount, *mbp);

	}

    return Dtime;
}

