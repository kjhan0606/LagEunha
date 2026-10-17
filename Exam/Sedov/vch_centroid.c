#include "voro.h"
#include "sedov.h"
#ifdef _OPENMP
#include<omp.h>
#endif


double vph3Dcentroid(Voro3D_GasParticle **Bp, int *mbp, ptype avgVolume){
	int isave = -1;
    ptype Dtime = 1.e10;

	Voro3D_GasParticle *bp = *Bp;
	int nbp = *mbp;


	MkLinkedList(bp,nbp);

	int nthreads=1;
	int mp=1024; 
	Voro3D_Vertex *tvorovertex;
	Voro3D_point *tneighbors,*tneighwork;
	Voro3D_point *nextcom = (Voro3D_point*)malloc(sizeof(Voro3D_point)*nbp);;
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

                int id = p[i].id; 
				int ishrink = 0;
				int ip = Voro3D_FindVC(&center,neighbors, neighwork,nneigh, vorovertex,mp,boxsize, 
						 ishrink, intwork);
				bp[id].volume = Voro3D_Volume_Polyhedron(vorovertex,ip);
        		bp[id].den  = bp[id].mass/bp[id].volume;

				nextcom[id].x = nextcom[id].y = nextcom[id].z = 0;
				int j;
				for(j=0;j<ip;j++){
					nextcom[id].x += vorovertex[j].x;
					nextcom[id].y += vorovertex[j].y;
					nextcom[id].z += vorovertex[j].z;
				}
				nextcom[id].x /= ip;
				nextcom[id].y /= ip;
				nextcom[id].z /= ip;

			} 
		} 
	} 
	free(tneighbors);
	free(tvorovertex);
	free(tneighwork);
	free(tp);
	free(tintwork);
	int i;
	for(i=0;i<nbp;i++){
		bp[i].x += 0.2*nextcom[i].x;
		bp[i].y += 0.2*nextcom[i].y;
		bp[i].z += 0.2*nextcom[i].z;
		bp[i].x = fmod(bp[i].x+Lx,Lx);
		bp[i].y = fmod(bp[i].y+Ly,Ly);
		bp[i].z = fmod(bp[i].z+Lz,Lz);
	}
	free(nextcom);
	int icenter = Nxp/2 + Nxp*(Nyp/2 + Nyp*(Nzp/2));
	bp[icenter].x = Lx*0.5L;
	bp[icenter].y = Ly*0.5L;
	bp[icenter].z = Lz*0.5L;

	Dtime = 1.e-6;
    return Dtime;
}

