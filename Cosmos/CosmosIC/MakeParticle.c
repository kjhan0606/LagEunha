#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<omp.h>
#include<mpi.h>
#include "eunha.h"
#include "cosmosic.h"

#define MAX(a,b) ( (a)>(b)?(a):(b))

void Set_Particle_2LPT_Vel(SimParameters *simpar,  GridInfo *fx, GridInfo *fy, GridInfo *fz, enum mtype itype){
	ptrdiff_t ix,iy,iz,jx,jy,jz, i,j,k;
	/*
	if(itype == 1){
		ix = rint(SIM_LXMIN(simpar,dm));
		iy = rint(SIM_LYMIN(simpar,dm));
		iz = rint(SIM_LZMIN(simpar,dm));
		jx = rint(SIM_LXMAX(simpar,dm));
		jy = rint(SIM_LYMAX(simpar,dm));
		jz = rint(SIM_LZMAX(simpar,dm));
	}
	else if(itype ==2){
		ix = rint(SIM_LXMIN(simpar,sph));
		iy = rint(SIM_LYMIN(simpar,sph));
		iz = rint(SIM_LZMIN(simpar,sph));
		jx = rint(SIM_LXMAX(simpar,sph));
		jy = rint(SIM_LYMAX(simpar,sph));
		jz = rint(SIM_LZMAX(simpar,sph));
	}
	*/
	ix = fx->ix; iy = fx->iy; iz = fx->iz;
	jx = fx->jx; jy = fx->jy; jz = fx->jz;
	DenType *ffx, *ffy, *ffz;
	/*
	PosType boxsize = COSBOXSIZE(simpar);
	*/
	PosType boxsize = NX(simpar);

	ffx = (DenType*)( fx+1);
	ffy = (DenType*)( fy+1);
	ffz = (DenType*)( fz+1);
	size_t mx = (jx-ix + 1) ;
	size_t my = (jy-iy + 1) ;

	size_t np = (jx-ix + 1)*(jy-iy + 1)*(jz-iz + 1);

	VelOnly  *vel = (VelOnly*)malloc(sizeof(VelOnly)*np);

	np = 0;

	double velmax = 0;

	for(k=iz;k<=jz;k++){
		size_t zi = k-iz;
		for(j=iy;j<=jy;j++){
			size_t yi = j-iy;
			for(i=ix;i<=jx;i++){
				size_t xi = i-ix;
				vel[np].vx = ffx[xi + mx*(yi + my * zi)];
				vel[np].vy = ffy[xi + mx*(yi + my * zi)];
				vel[np].vz = ffz[xi + mx*(yi + my * zi)];
				velmax = MAX(velmax, vel[np].vx);
				velmax = MAX(velmax, vel[np].vy);
				velmax = MAX(velmax, vel[np].vz);
				np++;
			}
		}
	}
#ifdef DEBUG
	DEBUGPRINT("P%d has velmax = %g : %d %d : %d %d : %d %d ::: %ld\n", MYID(simpar), velmax, 
			fx->ix, fx->jx, fx->iy, fx->jy, fx->iz,fx->jz, np);
#endif
	/*
	free(fx);free(fy);free(fz);
	*/
	float pfact = (Evol_PFACT(simpar) = AI(simpar) * NX(simpar));
	if(itype== totm || itype ==cdm) {
		DM_NP(simpar) = np;
		DM_BP(simpar) = (dmparticletype*)vel;
	}
	else if(itype ==sph) {
		SPH_NP(simpar) = np;
		SPH_BP(simpar) = (sphparticletype*)vel;
	}
}
void Set_Particle_2LPT_Pos(SimParameters *simpar,  GridInfo *fx, GridInfo *fy, GridInfo *fz, enum mtype itype){
	ptrdiff_t ix,iy,iz,jx,jy,jz, i,j,k;
	/*
	if(itype == cdm){
		ix = rint(SIM_LXMIN(simpar,dm));
		iy = rint(SIM_LYMIN(simpar,dm));
		iz = rint(SIM_LZMIN(simpar,dm));
		jx = rint(SIM_LXMAX(simpar,dm));
		jy = rint(SIM_LYMAX(simpar,dm));
		jz = rint(SIM_LZMAX(simpar,dm));
	}
	else if(itype ==sph){
		ix = rint(SIM_LXMIN(simpar,sph));
		iy = rint(SIM_LYMIN(simpar,sph));
		iz = rint(SIM_LZMIN(simpar,sph));
		jx = rint(SIM_LXMAX(simpar,sph));
		jy = rint(SIM_LYMAX(simpar,sph));
		jz = rint(SIM_LZMAX(simpar,sph));
	}
	*/
	ix = fx->ix; iy = fx->iy; iz = fx->iz;
	jx = fx->jx; jy = fx->jy; jz = fx->jz;
	DenType *ffx, *ffy, *ffz;
	/*
	PosType boxsize = COSBOXSIZE(simpar);
	PosType boxsize = NX(simpar);
	*/

	ffx = (DenType*)( fx+1);
	ffy = (DenType*)( fy+1);
	ffz = (DenType*)( fz+1);
	size_t mx = (jx-ix + 1) ;
	size_t my = (jy-iy + 1) ;


#ifdef DEBUG
	DEBUGPRINT("P%d is here %p %p %p ::: %d\n",MYID(simpar), ffx,ffy,ffz, itype);
#endif
	if(itype==totm || itype ==cdm) TwoLPTPosAssign(simpar, DM, dm);
	else if(itype ==sph) TwoLPTPosAssign(simpar, SPH, sph);
#ifdef DEBUG
	if(0){
		dmparticletype *bp = DM_BP(simpar);
		DEBUGPRINT("P%d: now %g %g %g :: %g %g %g\n",
				MYID(simpar), 
				XofP(simpar, bp), YofP(simpar, bp), ZofP(simpar, bp),
				XofP(simpar, bp+DM_NP(simpar)-1), YofP(simpar, bp+DM_NP(simpar)-1), ZofP(simpar, bp+DM_NP(simpar)-1));
		for(i=0;i<DM_NP(simpar);i++){
			if(PINDX(bp+i) == 1){
				DEBUGPRINT("p%d found %g %g %g : %g %g %g :: %ld\n",MYID(simpar),
						XofP(simpar, bp+i),YofP(simpar,bp+i), ZofP(simpar,bp+i),
						(DM_BP(simpar)+i)->vx, (DM_BP(simpar)+i)->vy, (DM_BP(simpar)+i)->vz, 
						PINDX(DM_BP(simpar)+i));
	
			}
		}
		MPI_Finalize();
		exit(99);
	}
#endif
}

void Set_Particle_Zeld(SimParameters *simpar,  GridInfo *fx, GridInfo *fy, GridInfo *fz, enum mtype itype){
	ptrdiff_t ix,iy,iz,jx,jy,jz, i,j,k;
	/*
	if(itype == cdm){
		ix = rint(SIM_LXMIN(simpar,dm));
		iy = rint(SIM_LYMIN(simpar,dm));
		iz = rint(SIM_LZMIN(simpar,dm));
		jx = rint(SIM_LXMAX(simpar,dm));
		jy = rint(SIM_LYMAX(simpar,dm));
		jz = rint(SIM_LZMAX(simpar,dm));
	}
	else if(itype ==sph){
		ix = rint(SIM_LXMIN(simpar,sph));
		iy = rint(SIM_LYMIN(simpar,sph));
		iz = rint(SIM_LZMIN(simpar,sph));
		jx = rint(SIM_LXMAX(simpar,sph));
		jy = rint(SIM_LYMAX(simpar,sph));
		jz = rint(SIM_LZMAX(simpar,sph));
	}
	*/
	ix = fx->ix; iy = fx->iy; iz = fx->iz;
	jx = fx->jx; jy = fx->jy; jz = fx->jz;
	float vfact = COS_VAMP1(simpar, itype);
	DenType *ffx, *ffy, *ffz;
	PosType boxsize = COSBOXSIZE(simpar);

	ffx = (DenType*)( fx+1);
	ffy = (DenType*)( fy+1);
	ffz = (DenType*)( fz+1);
	size_t mx = (jx-ix+1) ;
	size_t my = (jy-iy+1) ;

	size_t np = (jx-ix+1)*(jy-iy+1)*(jz-iz+1);

	VelOnly  *vel = (VelOnly*)malloc(sizeof(VelOnly)*np);

	np = 0;

	for(k=iz;k<=jz;k++){
		size_t zi = k-iz;
		for(j=iy;j<=jy;j++){
			size_t yi = j-iy;
			for(i=ix;i<=jx;i++){
				size_t xi = i-ix;
				vel[np].vx = vfact*ffx[xi + mx*(yi + my * zi)];
				vel[np].vy = vfact*ffy[xi + mx*(yi + my * zi)];
				vel[np].vz = vfact*ffz[xi + mx*(yi + my * zi)];
				np++;
			}
		}
	}
	free(fx);free(fy);free(fz);
	float pfact = (Evol_PFACT(simpar) = AI(simpar) * NX(simpar));
	if(itype ==totm || itype ==cdm) {
		DM_NP(simpar) = np;
		ZeldPosAssign(dmparticletype, simpar, vel);
		DM_BP(simpar) = (dmparticletype*)vel;
	}
	else if(itype ==sph) {
		SPH_NP(simpar) = np;
		ZeldPosAssign(sphparticletype, simpar, vel);
		SPH_BP(simpar) = (sphparticletype*)vel;
	}
}
