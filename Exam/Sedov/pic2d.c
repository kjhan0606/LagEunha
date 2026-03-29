#include "voro.h"
#include "pic2d.h"




CellType *cells=NULL;

ptype cellsize = Ly/Ny;

int colorizeit(float *, int, int, char *);

int makemap(Voro2D_GasParticle *bp, int np, int icount){
	int nx = Nximg;
	int ny = Nyimg;
	float *map = (float*)malloc(sizeof(float)*nx*ny);
	int i,j;
	int ii,jj;
//	for(i=0;i<nx*ny;i++) map[i] = 0;
	ptype imgcellsize = Lx/nx;

	for(j=0;j<ny;j++){
		ptype yp = (j+0.5)*imgcellsize;
		for(i=0;i<nx;i++){
			ptype xp = (i+0.5)*imgcellsize;
			ptype nearden = 0;
			ptype idist = 1.e20;
			int jy,ix;
			ix = xp/cellsize;
			jy = yp/cellsize;
			for(jj=jy-1;jj<=jy+1;jj++){
				int jp = (jj+Ny)%Ny;
				for(ii=ix-1;ii<=ix+1;ii++){
					int ip = (ii+Nx)%Nx;
					size_t ipixel = ip + Nx*jp;
					Voro2D_GasParticle *tmp = cells[ipixel].bp;
					while(tmp){
						float distx = fabs(tmp->x - xp);
						float disty = fabs(tmp->y - yp);
						if(distx>0.5*Lx) distx = Lx-distx;
						if(disty>0.5*Ly) disty = Ly-disty;
						float dist2 = distx*distx + disty*disty;
						if(dist2 < idist) {
							idist = dist2;
							nearden = tmp->den;
						}
						tmp = tmp->next;
					}

				}
			}
			map[i+nx*j] = nearden;

		}
	}


	char outfile[189]; 
	sprintf(outfile,"khmap.%.6d.ppm", icount); 
	colorizeit(map, nx,ny,outfile);

	free(map);
}
int maketscmap(Voro2D_GasParticle *bp, int np, int icount){
	int nx = Nximg;
	int ny = Nyimg;
	float cellsize = Ly/ny;
	float *map = (float*)malloc(sizeof(float)*nx*ny);
	int i;
	for(i=0;i<nx*ny;i++) map[i] = 0;
	ptype mscale = (Nxp/Lx)*(Nxp/Lx);

	for(i=0;i<np;i++){
		ptype pmass = mscale*bp[i].mass;
		ptype xp = bp[i].x/cellsize;
		ptype yp = bp[i].y/cellsize;
		int nearx = rint(xp);
		int neary = rint(yp);
		int ic = (nearx + nx) %nx;
		int jc = (neary + ny) %ny;
		ptype xmin = xp -nearx;
		ptype ymin = yp -neary;
		int icc = (nx+nearx + (int)(copysign(1., xmin))) %nx;
		int jcc = (ny+neary + (int)(copysign(1., ymin))) %ny;
		int iccc = (nx+nearx - (int)(copysign(1., xmin))) %nx;
		int jccc = (ny+neary - (int)(copysign(1., ymin))) %ny;
		float xd1 = fabs(xmin);
		float yd1 = fabs(ymin);
		float wx1 = (0.75-xd1*xd1)*pmass;
		float wy1 = (0.75-yd1*yd1);
		float wx3 = 0.5*pmass * (0.25+xd1*(xd1-1.));
		float wx2 = wx3 + pmass * xd1;
		float wy3 = 0.5 * (0.25+yd1*(yd1-1.));
		float wy2 = wy3 + yd1;

		map[ic+nx*jc] += wx1*wy1;
		map[icc+nx*jc] += wx2*wy1;
		map[iccc+nx*jc] += wx3*wy1;
		map[ic+nx*jcc] += wx1*wy2;
		map[icc+nx*jcc] += wx2*wy2;
		map[iccc+nx*jcc] += wx3*wy2;
		map[ic+nx*jccc] += wx1*wy3;
		map[icc+nx*jccc] += wx2*wy3;
		map[iccc+nx*jccc] += wx3*wy3;

	}

	char outfile[189];
	sprintf(outfile,"khmap.%.6d.ppm", icount);
	colorizeit(map, nx,ny,outfile);
	free(map);

}

void MkLinkedList(Voro3D_GasParticle *bp, int np){
	int i,j,k;
	if(cells == NULL) cells = (CellType*)malloc(sizeof(CellType)*Nx*Ny);
	cellsize = Lx/Nx;
	for(i=0;i<Nx*Ny;i++) {
		cells[i].bp = NULL;
		cells[i].np = 0;
	}
	for(i=0;i<np;i++){
		int ix,iy;
		ix = (bp[i].x/cellsize);
		iy = (bp[i].y/cellsize);
		size_t index = ix+Nx*iy;
		Voro2D_GasParticle *tmp = cells[index].bp;
		cells[index].bp = bp+i;
		cells[index].np ++;
		bp[i].next = tmp;
	}
}



Voro2D_point *Voro2D_FindNeighbor(int ix, int iy, int *nneigh, Voro2D_GasParticle *bp){
	int i,j,k;
	int np,mp;
	Voro2D_point *neigh;

	np = 0;
	for(j=iy-2;j<=iy+2;j++){
		int jj = (j+Ny)%Ny;
		for(i=ix-2;i<=ix+2;i++){
			int ii = (i+Nx)%Nx;
			np += cells[ii+Nx*jj].np;
		}
	}
	neigh = (Voro2D_point*)malloc(sizeof(Voro2D_point)*np);
	np = 0;
	for(j=iy-2;j<=iy+2;j++){
		int jflag=0;
		int jj = (j+Ny)%Ny;
		if(j<0) jflag = -1;
		else if(j>=Ny) jflag = 1;
		for(i=ix-2;i<=ix+2;i++){
			int iflag=0;
			int ii = (i+Nx)%Nx;
			if(i<0) iflag = -1;
			else if(i>=Nx) iflag = 1;
			Voro2D_GasParticle *tmp = cells[ii+Nx*jj].bp;
			while(tmp){
				neigh[np].x = tmp->x + iflag *Lx;
				neigh[np].y = tmp->y + jflag *Ly;
				neigh[np].id = tmp-bp;
				np++;
				tmp = tmp->next;
			}

		}
	}
	*nneigh = np;
	return neigh;
}

Voro2D_GasParticle *Voro2D_FindCellBP(int ix, int iy, int *mp, Voro2D_GasParticle *bp){
	int i,j,k;
	int np;
	Voro2D_GasParticle *res;

	int iix,iiy;
	iix = (ix+Nx)%Nx;
	iiy = (iy+Ny)%Ny;

	np = cells[iix+Nx*iiy].np;
	res = (Voro2D_GasParticle*)malloc(sizeof(Voro2D_GasParticle)*np);
	np = 0;
	int iflag,jflag;
	iflag = jflag = 0;
	if(iix < 0) iflag = -1;
	else if(iix >= Nx) iflag = 1;
	if(iiy < 0) jflag = -1;
	else if(iiy >= Ny) jflag = 1;
	size_t ipixel = iix+Nx*iiy;
	Voro2D_GasParticle *tmp = cells[ipixel].bp;
	while(tmp){
		res[np] = *tmp;
		res[np].id = tmp-bp;
		res[np].x += iflag *Lx;
		res[np].y += jflag *Ly;
		np++;
		tmp = tmp->next;
	}
	*mp = np;
	return res;
}
