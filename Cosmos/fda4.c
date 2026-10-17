#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
//#include "mpirms.h"
#include "cosmology.h"
#include "pmmain.h"
#include "timerutil.h"
#include "fda4.h"
#include "CosmosEvolFactor.h"

/*
float cputime0[2];
float cputime1[2];
*/


#define MIN(a,b) ( (a)<(b)?(a):(b))
#define MAX(a,b) ( (a)>(b)?(a):(b))

#ifdef XYZDBL
typedef double precisiontype;
#else
typedef float precisiontype;
#endif


void fda4(SimParameters *simpar, GridInfo *dengrid, EvolFact evolfact){
	ptrdiff_t i,j,k;
	ptrdiff_t local_ny_after_transpose = LOCAL_NY_AFTER_TRANSPOSE(simpar);
	ptrdiff_t local_y_start_after_transpose = LOCAL_Y_START_AFTER_TRANSPOSE(simpar);
	ptrdiff_t nx,ny,nz,mx, my;
	double rnx,rny,rnz;
	rnx = nx =NX(simpar);
	rny = ny =NY(simpar);
	rnz = nz =NZ(simpar);
	int myid = MYID(simpar);
	DenType *den = (DenType*)( dengrid + 1);
	PosType xs = dengrid->ix;
    PosType ys = dengrid->iy;
	PosType zs = dengrid->iz;
	float ffact1, ffact2;
	ffact1 = (4./3.)/2.;
	ffact2 = (-1./3.)/4.;
	float fact1, fact2;
	if(PMSTATUS(simpar) == PUSH){ 
		fact1 = evolfact.fact1_push; 
		fact2 = evolfact.fact2_push; 
		float a1, a2; 
		a1 = ANOW(simpar); a2 = a1 + ASTEP(simpar)/2; 
		fact2 = fact2/GetFact1(simpar,a1,a2); 
	} 
	else if (PMSTATUS(simpar) == PULL){ 
		fact1 = evolfact.fact1_pull; 
		fact2 = evolfact.fact2_pull; 
	} 
	else { 
		fact1 = evolfact.fact1; 
		fact2 = evolfact.fact2; 
	}



	mx = (dengrid->jx - dengrid->ix +1);
	my = (dengrid->jy - dengrid->iy +1);

#ifdef DEBUG
	{
		double xmin,ymin,zmin,xmax,ymax,zmax;
		xmin = ymin = zmin = 1.e20;
		xmax = ymax = zmax = -1.e20;
		dmparticletype *bp = DM_BP(simpar);
		for(i=0;i<DM_NP(simpar);i++){
			xmin = MIN(xmin, XofP(simpar, bp+i));
			ymin = MIN(ymin, YofP(simpar, bp+i));
			zmin = MIN(zmin, ZofP(simpar, bp+i));
			xmax = MAX(xmax, XofP(simpar, bp+i));
			ymax = MAX(ymax, YofP(simpar, bp+i));
			zmax = MAX(zmax, ZofP(simpar, bp+i));
		}
		DEBUGPRINT("P%d has %g %g %g %g %g %g: %d %d %d %d %d %d ::: %g %g\n",MYID(simpar),
				xmin,xmax,ymin,ymax,zmin,zmax, dengrid->ix,dengrid->jx,dengrid->iy,dengrid->jy,
				dengrid->iz,dengrid->jz, fact1, fact2);
	}
#endif

#ifdef USE_GPU
    {
        void call_cuda_fda4(int,int,int ,int , int , float *,
                pmparticletype *, long long , float , float );
        call_cuda_fda4(myid,nx,mx,ny, local_nz, DEN,dmp, simpar.dm.np, start_z, fact2);
    }

#else

#ifndef FASTFDA4
    if(DM_NP(simpar)>0) {
#       ifdef _OPENMP
#       pragma omp parallel for private(i)
#       endif
        SLOW_FDA4(simpar,DM,den,precisiontype, xs,ys,zs);
    }
    if(SPH_NP(simpar)>0) {
#       ifdef _OPENMP
#       pragma omp parallel for private(i)
#       endif
        SLOW_FDA4(simpar,SPH,den,precisiontype, xs,ys,zs);
    }
    if(STAR_NP(simpar)>0) { 
#       ifdef _OPENMP
#       pragma omp parallel for private(i)
#       endif
        SLOW_FDA4(simpar,STAR,den,precisiontype, xs,ys,zs);
    }
    if(AGN_NP(simpar)>0) { 
#       ifdef _OPENMP
#       pragma omp parallel for private(i)
#       endif
        SLOW_FDA4(simpar,AGN,den,precisiontype, xs,ys,zs);
    }
#else
    if(DM_NP(simpar)>0) {
#       ifdef _OPENMP
#       pragma omp parallel for private(i)
#       endif
        FAST_FDA4(simpar,DM,den,precisiontype, xs,ys,zs);
    }
    if(SPH_NP(simpar)>0) {
#       ifdef _OPENMP
#       pragma omp parallel for private(i)
#       endif
        FAST_FDA4(simpar,SPH,den,precisiontype, xs,ys,zs);
    }
    if(STAR_NP(simpar)>0) { 
#       ifdef _OPENMP
#       pragma omp parallel for private(i)
#       endif
        FAST_FDA4(simpar,STAR,den,precisiontype, xs,ys,zs);
    }
    if(AGN_NP(simpar)>0) { 
#       ifdef _OPENMP
#       pragma omp parallel for private(i)
#       endif
        FAST_FDA4(simpar,AGN,den,precisiontype, xs,ys,zs);
    }
#endif


#endif /* End of the USE_GPU */
}
    
