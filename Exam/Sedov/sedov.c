#include "voro.h"
#include "sedov.h"




CellType *cells=NULL;

ptype cellsize = Ly/Ny;

int colorizeit(float *, int, int, char *);

int makemap(Voro3D_GasParticle *bp, int np, int icount){
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
					Voro3D_GasParticle *tmp = cells[ipixel].bp;
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
int maketscmap(Voro3D_GasParticle *bp, int np, int icount){
	int nx = Nximg;
	int ny = Nyimg;
	float cellsize = Ly/ny;
	float *map = (float*)malloc(sizeof(float)*nx*ny);
	int i;
	for(i=0;i<nx*ny;i++) map[i] = 0;
	ptype mscale = (Nxp/Lx)*(Nxp/Lx);
	int Np = Nxp*Nyp;

	for(i=0;i<np;i++){
		if(bp[i].id/(Np) != (Nzp/2) ) continue;
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

void refinecenter(Voro3D_GasParticle *bp, int np){
	int i,j,k;
	int ix,iy,iz;
	int oindx = Nxp/2 +Nxp*(Nyp/2 + Nyp*(Nzp/2));
	if(1){
		ptype x,y,z, meansep = Lx/Nxp;
		ptype r,ratio;
		for(k=-1;k<2;k++){
			for(j=-1;j<2;j++){
				for(i=-1;i<2;i++){
					if(i==0 && j==0 && k==0) continue;
					int indx = oindx +k*Nxp*Nyp + j*Nxp + i;
					r = sqrt(k*k + j*j + i*i);
					ratio = 1.L/r;
					bp[indx].x = Lx*0.5L + i*meansep * ratio;
					bp[indx].y = Ly*0.5L + j*meansep * ratio;
					bp[indx].z = Lz*0.5L + k*meansep * ratio;
				}
			}
		}
	}
	if(1){
		ptype x,y,z, meansep = Lx/Nxp;
		ptype r,ratio;
		for(k=-2;k<3;k++){
			for(j=-2;j<3;j++){
				for(i=-2;i<3;i++){
					if(i==0 && j==0 && k==0) continue;
					int indx = oindx +k*Nxp*Nyp + j*Nxp + i;
					r = sqrt(k*k + j*j + i*i);
					if(r > sqrt(3*1.L) && r <= sqrt(3*4.L)) {
						ratio = 2.L/r;
						bp[indx].x = Lx*0.5L + i*meansep * ratio;
						bp[indx].y = Ly*0.5L + j*meansep * ratio;
						bp[indx].z = Lz*0.5L + k*meansep * ratio;
					}
				}
			}
		}
	}
}

Voro3D_GasParticle *mkinitial(int *mp){
	int i,j,k,np;
	ptype pixsize = Lx/Nxp;
	Voro3D_GasParticle *res;
	res = (Voro3D_GasParticle*)malloc(sizeof(Voro3D_GasParticle)*Nxp*Nyp*Nzp);
	ptype meanvol = Lx*Ly*Lz/Nxp/Nyp/Nzp;
	for(k=0;k<Nzp;k++){
		ptype z = (ptype)(k+0.0)*Lz/(ptype)Nzp;
		for(j=0;j<Nyp;j++){ 
			ptype y = (ptype)(j+0.0)*Ly/(ptype)Nyp;
			for(i=0;i<Nxp;i++){
				ptype x = (ptype)(i+0.0)*Lx/(ptype)Nxp; 
				size_t indx = i+Nxp*(j+Nyp*k); 
				res[indx].id  = indx; 
				res[indx].x  = fmod(x+Lx,Lx);
				res[indx].y  = fmod(y+Ly,Ly);
				res[indx].z  = fmod(z+Lz,Lz);
				res[indx].vx  = res[indx].vy = res[indx].vz = 0; 
				res[indx].den  = 1; 
				res[indx].mass  = res[indx].den*meanvol;
				res[indx].pressure  = bckPressure;
				res[indx].ie  = res[indx].pressure *meanvol/(Gamma-1); 
				res[indx].ke  = 0;
				res[indx].te  = res[indx].ie + res[i].ke;
				res[indx].csound  = sqrt(Gamma*res[indx].pressure/res[indx].den);
			}
		}
	}
	int icenter = (Nxp/2)+Nxp*( Nyp/2 + Nyp*(Nzp/2));
	res[icenter].x = (ptype)(Nxp/2+0.5)*Lx/Nxp;
	res[icenter].y = (ptype)(Nyp/2+0.5)*Ly/Nyp;
	res[icenter].z = (ptype)(Nzp/2+0.5)*Lz/Nzp;
	res[icenter].ie = eBlast;
	res[icenter].te = res[icenter].ie;
	res[icenter].pressure = res[icenter].ie*(Gamma-1)/meanvol;
	res[icenter].csound  = sqrt(Gamma*res[icenter].pressure/res[icenter].den);
	printf("icenter = %d\n", icenter);
	int nbp = Nxp*Nyp*Nzp;
	if(0) refinecenter(res,nbp);
	*mp = nbp;
	return res;
}


void MkLinkedList(Voro3D_GasParticle *bp, int np){
	int i,j,k;
	if(cells == NULL) cells = (CellType*)malloc(sizeof(CellType)*Nx*Ny*Nz);


	cellsize = Lx/Nx;

	for(i=0;i<Nx*Ny*Nz;i++) {
		cells[i].bp = NULL;
		cells[i].np = 0;
	}

	for(i=0;i<np;i++){
		int ix,iy,iz;
		ix = (bp[i].x/cellsize);
		iy = (bp[i].y/cellsize);
		iz = (bp[i].z/cellsize);
		size_t index = ix+Nx*(iy+Ny*iz);
		Voro3D_GasParticle *tmp = cells[index].bp;
		cells[index].bp = bp+i;
		cells[index].np ++;
		bp[i].next = tmp;
	}
}

int get_npmax(){
	int k;
	int np = 0;
	for(k=0;k<Nz*Ny*Nx;k++) {
		np = MAX(np, cells[k].np);
	}
	return np;
}


void Voro3D_FindNeighbor(Voro3D_point *neigh, int ix, int iy, int iz, int *nneigh, Voro3D_GasParticle *bp){
	int i,j,k;
	int np,mp;

	np = 0;
	for(k=iz-1;k<=iz+1;k++){ 
		int kk = (k+Nz)%Nz; 
		int kflag =0;
		if(k<0) kflag = -1;
		else if(k>=Nz) kflag = 1;
		for(j=iy-1;j<=iy+1;j++){ 
			int jflag=0; 
			int jj = (j+Ny)%Ny; 
			if(j<0) jflag = -1; 
			else if(j>=Ny) jflag = 1; 
			for(i=ix-1;i<=ix+1;i++){ 
				int iflag=0; 
				int ii = (i+Nx)%Nx; 
				if(i<0) iflag = -1; 
				else if(i>=Nx) iflag = 1; 
				Voro3D_GasParticle *tmp = cells[ii+Nx*(jj+Ny*kk)].bp; 
				while(tmp){ 
					neigh[np].x = tmp->x + iflag *Lx; 
					neigh[np].y = tmp->y + jflag *Ly; 
					neigh[np].z = tmp->z + kflag *Lz; 
					neigh[np].id = tmp-bp; 
					neigh[np].csound = tmp->csound; 
					neigh[np].pressure = tmp->pressure; 
					np++; 
					tmp = tmp->next; 
				}
			}

		}
	}
	*nneigh = np;
}

void Voro3D_FindCellBP(Voro3D_GasParticle *p, int ix, int iy, int iz, int *mp, Voro3D_GasParticle *bp){
	int i,j,k;
	int np;

	int iix,iiy,iiz;
	iix = (ix+Nx)%Nx;
	iiy = (iy+Ny)%Ny;
	iiz = (iz+Nz)%Nz;
	np = 0;
	int iflag,jflag,kflag;
	iflag = jflag = kflag =  0;
	if(iix < 0) iflag = -1;
	else if(iix >= Nx) iflag = 1;
	if(iiy < 0) jflag = -1;
	else if(iiy >= Ny) jflag = 1;
	if(iiz < 0) kflag = -1;
	else if(iiz >= Nz) kflag = 1;
	size_t ipixel = iix+Nx*(iiy + Ny*iiz);
	Voro3D_GasParticle *tmp = cells[ipixel].bp;
	while(tmp){
		p[np] = *tmp;
		p[np].id = tmp-bp;
		p[np].x += iflag *Lx;
		p[np].y += jflag *Ly;
		p[np].z += kflag *Lz;
		np++;
		tmp = tmp->next;
	}
	*mp = np;
}
