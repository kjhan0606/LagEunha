#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "cosmosic.h"
#include "parallelIO.h"
#include "flrw.h"
#define MAX(a,b) ( (a)>(b)?(a):(b))

GridInfo *makeFftDomainVolume(SimParameters *simpar, size_t densize){
    GridInfo *a = (GridInfo*)malloc(sizeof(GridInfo)+sizeof(DenType)*densize);
    a->nx = NX(simpar); a->ny = NY(simpar); a->nz = NZ(simpar);
    a->ix = a->iy = 0;
    a->jx = NX(simpar)-1;
    a->jy = NY(simpar)-1;
    a->iz = LOCAL_Z_START(simpar);
    a->jz = a->iz + LOCAL_NZ(simpar)-1;
    a->npix = (ptrdiff_t)(a->jx-a->ix+1) * (ptrdiff_t)(a->jy-a->iy+1) * (ptrdiff_t)(a->jz-a->iz+1);
    return a;
}
GridInfo *makeRmsDomainVolume(SimParameters *simpar, enum mtype itype, int nbuff){
    GridInfo *a = (GridInfo*)malloc(sizeof(GridInfo));
    a->nx = NX(simpar); a->ny = NY(simpar); a->nz = NZ(simpar);

    if(itype <=cdm) {
        a->ix = rint(SIM_LXMIN(simpar,dm)) - nbuff;
        a->iy = rint(SIM_LYMIN(simpar,dm)) - nbuff;
        a->iz = rint(SIM_LZMIN(simpar,dm)) - nbuff;
        a->jx = rint(SIM_LXMAX(simpar,dm)) - 1 + nbuff;
        a->jy = rint(SIM_LYMAX(simpar,dm)) - 1 + nbuff;
        a->jz = rint(SIM_LZMAX(simpar,dm)) - 1 + nbuff;
    }
    else{
        a->ix = rint(SIM_LXMIN(simpar,sph)) - nbuff;
        a->iy = rint(SIM_LYMIN(simpar,sph)) - nbuff;
        a->iz = rint(SIM_LZMIN(simpar,sph)) - nbuff;
        a->jx = rint(SIM_LXMAX(simpar,sph)) - 1 + nbuff;
        a->jy = rint(SIM_LYMAX(simpar,sph)) - 1 + nbuff;
        a->jz = rint(SIM_LZMAX(simpar,sph)) - 1 + nbuff;
    }

    a->npix = (ptrdiff_t)(a->jx-a->ix+1) *(ptrdiff_t)(a->jy-a->iy+1) *(ptrdiff_t)(a->jz-a->iz+1);
    a = realloc(a, sizeof(GridInfo)+sizeof(DenType)*a->npix);
    return a;
}


void readGraficFile(SimParameters *simpar){
	size_t i,j,k;
	enum mtype itype;
	expectDMandSphNumParticle(simpar);
/* initialization of some parameters */
	AI(simpar) = 1;
	ZINIT(simpar) = AMAX(simpar)-1;
	ANOW(simpar) = 1;
	OMEI(simpar) = Omega_matter(simpar, 1.);

	for(itype=totm;itype<star;itype++){
		if(flagHydro(simpar)=='Y') continue;
		if(itype ==cdm && DM_TNP(simpar) ==0) continue;
        else if(itype ==sph && SPH_TNP(simpar) ==0) continue;

        size_t densize = 2*LOCAL_GRID_SIZE(simpar); // factor 2 is needed due to complex number
        densize = MAX(10, densize);
		int nx = NX(simpar);
		int ny = NY(simpar);
		int local_nz = LOCAL_NZ(simpar);
		int local_z_start = LOCAL_Z_START(simpar);
		GridInfo *ffx, *ffy, *ffz;
		GridInfo *fx, *fy, *fz;

// This is for the velocity.
		ffx = makeFftDomainVolume(simpar,densize);
		pGraficRead(ffx,simpar,itype, 3); // read grafic vx
		fx = makeRmsDomainVolume(simpar, itype,0);
		gmigrate(ffx,fx,DM_DDINFO(simpar), NDDINFO(simpar));

		ffy = makeFftDomainVolume(simpar,densize);
		pGraficRead(ffy,simpar,itype, 4); // read grafic vy
		fy = makeRmsDomainVolume(simpar, itype,0);
		if(0){
			double maxval = 0.;
			DenType *den = (DenType*)(fy+1);
			for(k=0;k<local_nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++){
				maxval = MAX(maxval, den[i+nx*(j+ny*k)]);
			}
			DEBUGPRINT("+P%d has fftw maxval of composite = %g ::: %g -- %g\n",
					MYID(simpar), maxval, SIM_LXMIN(simpar,dm), SIM_LXMAX(simpar,dm));
			DEBUGPRINT("P%d has %d %d %d -- %d %d %d\n",MYID(simpar), ffx->ix,ffx->iy,ffx->iz,
					ffx->jx,ffx->jy, ffx->jz);
			MPI_Finalize();
			exit(99);
		}
		gmigrate(ffy,fy,DM_DDINFO(simpar), NDDINFO(simpar));

        ffz = makeFftDomainVolume(simpar,densize);
		pGraficRead(ffz,simpar,itype, 5); // read grafic vz
		fz = makeRmsDomainVolume(simpar, itype,0);
		gmigrate(ffz,fz,DM_DDINFO(simpar), NDDINFO(simpar));
	
		Set_Particle_2LPT_Vel(simpar, fx,fy,fz,itype);

		free(fx);free(fy);free(fz);

// This is for the position.
		ffx = makeFftDomainVolume(simpar,densize);
		pGraficRead(ffx,simpar,itype, 0); // read grafic dx
		fx = makeRmsDomainVolume(simpar, itype,0);
		gmigrate(ffx,fx,DM_DDINFO(simpar), NDDINFO(simpar));


		ffy = makeFftDomainVolume(simpar,densize);
		pGraficRead(ffy,simpar,itype, 1); // read grafic dy
		fy = makeRmsDomainVolume(simpar, itype,0);
		gmigrate(ffy,fy,DM_DDINFO(simpar), NDDINFO(simpar));

		ffz = makeFftDomainVolume(simpar,densize);
		pGraficRead(ffz,simpar,itype, 2); // read grafic dz
		fz = makeRmsDomainVolume(simpar, itype,0);
		gmigrate(ffz,fz,DM_DDINFO(simpar), NDDINFO(simpar));

		Set_Particle_2LPT_Pos(simpar, fx,fy,fz,itype);

		free(fx);free(fy);free(fz);
		MPI_Finalize();exit(0);

	}

}
