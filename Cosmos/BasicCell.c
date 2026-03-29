#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
#include "BasicCell.h"

void BuildLinkedList(SimParameters *simpar, ptrdiff_t mx, ptrdiff_t my, ptrdiff_t mz, 
		PosType xmin, PosType ymin, PosType zmin,PosType zmax){ 
	TreeLinkedCell *BasicCell = BASICCELL(simpar); 
	ptrdiff_t i,j,k,ix,iy,iz; 
	for(i=0;i<mx*my*mz;i++) { 
		BasicCell[i].link = NULL; 
		BasicCell[i].nmem = 0; 
	} 
#ifdef _OPENMP
#pragma omp parallel
	{
		int ithread = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		int izs,izf;
		int nz = (zmax-zmin)/CellWidth;
		int threadwidth = (nz+nthreads-1)/nthreads;
		izs = threadwidth * ithread;
		izf = threadwidth * (ithread + 1);
		if(ithread == nthreads-1)izf = nz;

		OMP_LinkParticles(simpar, DM,dm, BasicCell,mx,my,xmin,ymin,zmin,izs,izf); 
		OMP_LinkParticles(simpar, SPH,sph, BasicCell,mx,my,xmin,ymin,zmin,izs,izf); 
		OMP_LinkParticles(simpar, STAR,star, BasicCell,mx,my,xmin,ymin,zmin,izs,izf); 
		OMP_LinkParticles(simpar, AGN,agn, BasicCell,mx,my,xmin,ymin,zmin,izs,izf);
	}
#else
	LinkParticles(simpar, DM,dm, BasicCell,mx,my,xmin,ymin,zmin); 
	LinkParticles(simpar, SPH,sph, BasicCell,mx,my,xmin,ymin,zmin); 
	LinkParticles(simpar, STAR,star, BasicCell,mx,my,xmin,ymin,zmin); 
	LinkParticles(simpar, AGN,agn, BasicCell,mx,my,xmin,ymin,zmin);
#endif
}

