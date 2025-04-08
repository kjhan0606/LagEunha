#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "voro.h"
#include "rt.h"
#include "color.h"




static CellType *cells=NULL;

//postype cellsize = Ly/Ny;

//void colorizeit(float *, int, int, char *);

int rt_makemap(SimParameters *simpar, int icount){
	postype cellsize = RT_GridSize(simpar);
	postype Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	postype Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	int nximg = NX(simpar);
	int nyimg = NY(simpar);
	float *map = (float*)malloc(sizeof(float)*nximg*nyimg);
	float *img = (float*)malloc(sizeof(float)*nximg*nyimg);
	float *ndist = (float*)malloc(sizeof(float)*nximg*nyimg);
	int i,j;
	int ii,jj;
	for(i=0;i<nximg*nyimg;i++) map[i] = 0;
	for(i=0;i<nximg*nyimg;i++) ndist[i] = 1.e20;
	postype pixsize = Lx/nximg;
	postype xmin,ymin,xmax,ymax;
	xmin = RT_XMIN(simpar);
	ymin = RT_YMIN(simpar);
	xmax = RT_XMAX(simpar);
	ymax = RT_YMAX(simpar);
	int mx,my;
	mx = ceil((xmax-xmin)/cellsize);
	my = ceil((ymax-ymin)/cellsize);
	cells = (CellType*)malloc(sizeof(CellType)*mx*my);
	for(i=0;i<mx*my;i++){
		cells[i].link = NULL;
		cells[i].nmem = 0;
	}
	treevoroparticletype *bp = VORO_TBP(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		int ix,iy;
		ix = ((bp[i].x-xmin)/cellsize);
		iy = ((bp[i].y-ymin)/cellsize);
		int index = ix+mx*iy;
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)(bp+i);
		cells[index].nmem ++;
		bp[i].next = tmp;
	}

	for(j=0;j<nyimg;j++){
		postype yp = (j+0.5)*pixsize; // img pixel position
		if(yp < ymin || yp >= ymax) continue;
		for(i=0;i<nximg;i++){
			postype xp = (i+0.5)*pixsize; // img pixel position
			if(xp < xmin || xp >= xmax) continue;
			postype nearden = 0;
			postype idist = 1.e20;
			int jy,ix;
			ix = (xp-xmin)/cellsize;
			jy = (yp-ymin)/cellsize;
			for(jj=jy-1;jj<=jy+1;jj++){
				if(jj<0 || jj>=my) continue;
				for(ii=ix-1;ii<=ix+1;ii++){
					if(ii<0 || ii>=mx) continue;
					size_t ipixel = ii + mx*jj;
					struct linkedlisttype *tmp = cells[ipixel].link;
					while(tmp){
						treevoroparticletype *tt = (treevoroparticletype*)tmp;
						if(tt->x >= xmin && tt->x < xmax && tt->y >=ymin && tt->y < ymax){
							postype distx = fabs(tt->x - xp);
							postype disty = fabs(tt->y - yp);
							postype dist2 = distx*distx + disty*disty;
							if(dist2 < idist) {
								idist = dist2;
								nearden = tt->den;
							}
						}
						tmp = tmp->next;
					}

				}
			}
			map[i+nximg*j] = nearden;
			ndist[i+nximg*j] = idist;

		}
	}
	int id;
	if(MYID(simpar)==0){
		float *idist = (float*)malloc(sizeof(float)*nximg*nyimg);
		MPI_Status status;

		for(id=1;id<NID(simpar);id++){
			MPI_Recv(idist, nximg*nyimg,MPI_FLOAT,id,id, MPI_COMM_WORLD, &status);
			MPI_Recv(img, nximg*nyimg,MPI_FLOAT,id,id, MPI_COMM_WORLD, &status);
			for(j=0;j<nximg*nyimg;j++){
				if(idist[j] < ndist[j]){
					map[j] = img[j];
					ndist[j] = idist[j];
				}
			}
		}
		free(idist);
	}
	else {
		MPI_Send(ndist, nximg*nyimg,MPI_FLOAT,0,MYID(simpar), MPI_COMM_WORLD);
		MPI_Send(map, nximg*nyimg,MPI_FLOAT,0,MYID(simpar), MPI_COMM_WORLD);
	}

	if(MYID(simpar)==0){
		float *img2 = (float*)malloc(sizeof(float)*nximg*nyimg);;
		for(j=0;j<nyimg;j++){
			for(i=0;i<nximg;i++){
				img2[i+nximg*(nyimg-j-1)] = map[i+nximg*j];
			}
		}
		char outfile[189]; 
		sprintf(outfile,"rtmap.%.6d.ppm", icount); 
		colorizeit(5,img2, nximg,nyimg,"rt.sao", outfile);
		free(img2);
	}
	free(cells);
 	free(map); free(img);
	free(ndist);
	return 0;
}
/*
int maketscmap(SimParameters *simpar, treevoroparticletype *bp, int np, int icount){
	int nx = NX(simpar);
	int ny = NY(simpar);
	float cellsize = RT_GridSize(simpar);
	float Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	float Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	int i;
	float *map = (float*)malloc(sizeof(float)*nx*ny);
	for(i=0;i<nx*ny;i++) map[i] = 0;
	postype mscale = (nx/Lx)*(ny/Ly);

	for(i=0;i<np;i++){
		postype pmass = mscale*bp[i].mass;
		postype xp = bp[i].x/cellsize;
		postype yp = bp[i].y/cellsize;
		int nearx = rint(xp);
		int neary = rint(yp);
		int ic = (nearx + nx) %nx;
		int jc = (neary + ny) %ny;
		postype xmin = xp -nearx;
		postype ymin = yp -neary;
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
//	colorizeit(map, nx,ny,outfile);
	free(map);

	return 0;
}
*/

postype getDen(postype rho1, postype rho2, postype deltay, postype ycen,postype y){
	postype rho = rho1+(rho2-rho1)/(1+exp(-(y-ycen)/deltay));
	return rho;
}
/* for hydrostatic equilibrium, the pressure at y is obtained by integral of g(y)*rho(y) */
postype  getRTpressure(SimParameters *simpar, postype rho1, postype rho2, postype ycen, postype y){
	postype deltay = RT_Deltay(simpar);
	postype pressure = RT_Phalf(simpar) + RT_ACC(simpar)*getDen(rho1,rho2, deltay, ycen, y)*(y-ycen);
	return pressure;
}
treevoroparticletype *rt_mkinitial(SimParameters *simpar, int *mp){
	int i,j,k;
	treevoroparticletype *res;
	postype rho1 = RT_DEN1(simpar), rho2 = RT_DEN2(simpar);
	postype deltay = RT_Deltay(simpar);
	postype dvy0 = RT_Vperturb(simpar);
	postype xmin,ymin,xmax,ymax;
	xmin = RT_XMIN(simpar);
	ymin = RT_YMIN(simpar);
	xmax = RT_XMAX(simpar);
	ymax = RT_YMAX(simpar);
	int nx = NX(simpar);
	int ny = NY(simpar);
	float Lx = SIMBOX(simpar).x.max;
	float Ly = SIMBOX(simpar).y.max;
	postype ycen = 0.5*Ly;


	DEBUGPRINT("P%d has initial set nx/y= %d %d Lx/y= %g %g rmin= %g %g rmax= %g %g\n",
			MYID(simpar), nx,ny,Lx,Ly,xmin,ymin,xmax,ymax);
	DEBUGPRINT("P%d has initial set den1/2 = %g %g ycen= %g\n", 
			MYID(simpar), rho1, rho2, ycen);

	res = (treevoroparticletype*)malloc(sizeof(treevoroparticletype)*nx*ny);
	postype meanvol = Lx*Ly/nx/ny;
	size_t np = 0;
	for(j=0;j<ny;j++){
		postype rho;
		postype y = (postype)(j+0.5)*Ly/(postype)ny;
		char iregion;
		rho = getDen(rho1, rho2, deltay, ycen, y); 
		iregion = 0;
		postype vx = 0;
		for(i=0;i<nx;i++){
			postype vy;
			postype x = (postype)(i+0.5)*Lx/(postype)nx;
			size_t indx = i+nx*j;
			if(y>=0.3*Ly && y < 0.7*Ly){
				vy = dvy0*(1+cos(8*M_PI*(x+0.5*Lx)))*(1+cos(5*M_PI*(y-ycen)));
			}
			else{
				vy = 0;
			}
			if(x>=xmin && x < xmax && y>=ymin && y < ymax){
				res[np].u4if.indx  = indx;
				res[np].u4if.Flag[ENDIAN_OFFSET]  = iregion;
				res[np].x  = x;
				res[np].y  = y;
				res[np].vx  = vx;
				res[np].vy  = vy;
				res[np].mass  = rho*meanvol;
				res[np].den  = rho;
				res[np].pressure  = getRTpressure(simpar, rho1, rho2, ycen, y);
				res[np].ie  = res[indx].pressure *meanvol/(Gamma-1);
				res[np].ke  = 0.5*(res[indx].vx*res[indx].vx + res[indx].vy*res[indx].vy);
				res[np].csound  = 1;
				np++;
			}
		}
	}
	res = (treevoroparticletype*)realloc(res, sizeof(treevoroparticletype)*np);
	int nbp = np;
	*mp = nbp;
	VORO_TBP(simpar) = res;
	VORO_BP(simpar) = (voroparticletype*)res;
	VORO_NP(simpar) = nbp;
	for(i=0;i<np;i++) UNSET_P_FLAG(simpar, VORO, i, BoundaryGhostflag);
	DEBUGPRINT("P%d has np = %ld in box %g %g : %g %g\n", MYID(simpar), np, xmin,ymin,xmax,ymax);



	rt_MkLinkedList(simpar);
	void rt_initfindmass(SimParameters *);
	rt_initfindmass(simpar);

	for(i=0;i<nbp;i++){
		res[i].ie = res[i].pressure*res[i].volume/(Gamma-1);
		res[i].ke = Half*res[i].mass*( res[i].vx*res[i].vx + res[i].vy*res[i].vy);
		res[i].te = res[i].ie + res[i].ke;
		res[i].csound = sqrt(Gamma*res[i].pressure/res[i].den);
	}
	free(VORO_TBPP(simpar));
	VORO_NPAD(simpar) = 0;


	return res;
}


void rt_MkLinkedList(SimParameters *simpar){
	treevoroparticletype *bp = VORO_TBP(simpar);
	int np = VORO_NP(simpar);
	int i,j,k;
	postype xmin,ymin,xmax,ymax;
	int mx,my;
	postype pwidth = RT_GridSize(simpar);

	BASICCELL_CELLWIDTH(simpar) = pwidth;
	BASICCELL_INVCELLWIDTH(simpar) = 1./pwidth;


	xmin = RT_XMIN(simpar)-pwidth;
	ymin = RT_YMIN(simpar)-pwidth;
	xmax = RT_XMAX(simpar)+pwidth;
	ymax = RT_YMAX(simpar)+pwidth;

	mx = ceil((xmax-xmin)/pwidth);
	my = ceil((ymax-ymin)/pwidth);

	BASICCELL_MX(simpar) = mx;
	BASICCELL_MY(simpar) = my;

	cells = VORO_BASICCELL(simpar)= (CellType*)realloc(VORO_BASICCELL(simpar),sizeof(CellType)*mx*my);


	postype cellsize = RT_GridSize(simpar);

	for(i=0;i<mx*my;i++) {
		cells[i].link = NULL;
		cells[i].nmem = 0;
	}

	DEBUGPRINT("P%d has np= %ld: %d %d: %g : %g %g\n", MYID(simpar), VORO_NP(simpar),
			mx,my, cellsize, xmin,ymin);

	for(i=0;i<VORO_NP(simpar);i++){
		int ix,iy;
		ix = ((bp[i].x-xmin)/cellsize);
		iy = ((bp[i].y-ymin)/cellsize);
		size_t index = ix+mx*iy;
		/*
		if(index >=mx*my){
			printf("error detected: %ld : %d %d : %g %g %g %g\n", 
					index, mx,my, bp[i].x, bp[i].y, xmin,ymin);
		}
		*/
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)(bp+i);
		cells[index].nmem ++;
		bp[i].next = tmp;
	}
#ifdef DEBUG
	DEBUGPRINT("P%d has %g %g : %g %g\n", MYID(simpar), RT_XMIN(simpar), RT_YMIN(simpar),
			RT_XMAX(simpar), RT_YMAX(simpar));
#endif
	RT_TreeAllParticlePadding(simpar, pwidth);
	bp = VORO_TBPP(simpar);
	postype Xmin,Ymin,Xmax,Ymax;
	Xmin = Ymin = 1.e10;
	Xmax = Ymax = -1.e10;
	for(i=0;i<VORO_NPAD(simpar);i++){
		int ix,iy;
		ix = ((bp[i].x-xmin)/cellsize);
		iy = ((bp[i].y-ymin)/cellsize);
		Xmin = MIN(Xmin,bp[i].x);
		Ymin = MIN(Ymin,bp[i].y);
		Xmax = MAX(Xmax,bp[i].x);
		Ymax = MAX(Ymax,bp[i].y);
		size_t index = ix+mx*iy;
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)(bp+i);
		cells[index].nmem ++;
		bp[i].next = tmp;
	}
#ifdef DEBUG
	DEBUGPRINT("P%d now after padding: %g %g : %g %g\n", MYID(simpar),Xmin,Ymin,Xmax,Ymax);
#endif
}

Voro2D_point *rt_Voro2D_FindNeighbor(SimParameters *simpar,
		int ix, int iy, int *nneigh, treevoroparticletype *bp){
	int i,j,k;
	int np,mp;
	Voro2D_point *neigh;
	int mx,my;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	cells = VORO_BASICCELL(simpar);
	float Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	float Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;

	np = 0;
	for(j=iy-1;j<=iy+1;j++){
		if(j <0 || j >=my) continue;
		for(i=ix-1;i<=ix+1;i++){
			if(i<0 || i>=mx) continue;
			np += cells[i+mx*j].nmem;
		}
	}
	neigh = (Voro2D_point*)malloc(sizeof(Voro2D_point)*np);
	np = 0;
	for(j=iy-1;j<=iy+1;j++){
		if(j<0 || j >=my)continue;
		for(i=ix-1;i<=ix+1;i++){
			if(i<0 || i>=mx) continue;
			struct linkedlisttype *tmp = cells[i+mx*j].link;
			while(tmp){
				treevoroparticletype *tt = (treevoroparticletype*)tmp;
				neigh[np].x = tt->x;
				neigh[np].y = tt->y;
				neigh[np].indx = PINDX(tt);
				neigh[np].csound = tt->csound;
				neigh[np].pressure = tt->pressure;
				neigh[np].bp = (void*)tt;
				np++;
				tmp = tmp->next;
			}

		}
	}
	*nneigh = np;
	return neigh;
}

treevoroparticletype *rt_Voro2D_FindCellBP(SimParameters *simpar,
		int ix, int iy, int *mp, treevoroparticletype *bp){
	int i,j,k;
	int np;
	treevoroparticletype *res;
	int mx,my;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	cells = VORO_BASICCELL(simpar);
//	postype Ly = SIMBOX(simpar).y.max;


	np = cells[ix+mx*iy].nmem;
	res = (treevoroparticletype*)malloc(sizeof(treevoroparticletype)*np);
	np = 0;
	/*
	int iflag,jflag;
	iflag = jflag = 0;
	if(iix < 0) iflag = -1;
	else if(iix >= mx) iflag = 1;
	if(iiy < 0) jflag = -1;
	else if(iiy >= my) jflag = 1;
	*/
	size_t ipixel = ix+mx*iy;
	struct linkedlisttype *tmp = cells[ipixel].link;
	while(tmp){
		if(!IS_FLAG_ON(tmp,BoundaryGhostflag)){
			treevoroparticletype *tt = (treevoroparticletype*)tmp;
//			if(tt->y > 0.3*Ly && tt->y < 0.7*Ly)
			{
				res[np] = *tt;
				res[np].bp = tt;
				/*
				res[np].u4if.indx = tt-bp;
				res[np].x += iflag *Lx;
				res[np].y += jflag *Ly;
				*/
				np++;
			}
		}
		tmp = tmp->next;
	}
	*mp = np;
	return res;
}
