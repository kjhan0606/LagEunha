#include "voro.h"
#include "sedov.h"
#ifdef _OPENMP
#include <omp.h>
#endif


int main(int argc,char **argv){
	int np;
	int nlog = 4;
	ptype time,dt;
	FILE *fp = fopen(argv[1],"r");
	FILE *wp = fopen(argv[2],"w");
	fread(&np,sizeof(int),1,fp);
	Voro3D_GasParticle *bp = (Voro3D_GasParticle*)malloc(sizeof(Voro3D_GasParticle)*(np));
    fread(bp,sizeof(Voro3D_GasParticle), np,fp);
    fread(&time,sizeof(ptype), 1,fp);
    fread(&dt,sizeof(ptype), 1,fp);
    fclose(fp);

	MkLinkedList(bp,np);
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
#pragma omp for schedule(dynamic)
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
				bp[id].den = findSDen3D(vorovertex, ip, neighwork, intwork, mp, bp);
			}
		}
	}
	Voro3D_point center;
    center.x = Lx/2 + 0.5*Lx/Nxp;
    center.y = Ly/2 + 0.5*Ly/Nyp;
    center.z = Lz/2 + 0.5*Lz/Nzp;

	int i;
    for(i=0;i<np;i++){
        ptype tmpx, tmpy, tmpz;
        tmpx = (bp[i].x-center.x);
        tmpy = (bp[i].y-center.y);
        tmpz = (bp[i].z-center.z);
        ptype dist;
        dist = sqrt(tmpx*tmpx+tmpy*tmpy+tmpz*tmpz);
        if(dist < 1.) fprintf(wp,"%g %g\n", dist, bp[i].den);
    }
    fclose(wp);

}

