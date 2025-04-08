#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "voro.h"
#include "kh.h"
#include "color.h"




static CellType *cells=NULL;

//postype cellsize = Ly/Ny;

//void colorizeit(float *, int, int, char *);

int makemap(SimParameters *simpar, int icount){
	postype cellsize = KH_GridSize(simpar);
	postype Lx = KH_SIMBOX(simpar).x.max - KH_SIMBOX(simpar).x.min;
	postype Ly = KH_SIMBOX(simpar).y.max - KH_SIMBOX(simpar).y.min;
	int nximg = NX(simpar);
	int nyimg = NY(simpar);
	float *map = (float*)malloc(sizeof(float)*nximg*nyimg);
	float *img = (float*)malloc(sizeof(float)*nximg*nyimg);
	postype meanvol = Lx*Ly/nximg/nyimg;
	int i,j;
	int ii,jj;
	for(i=0;i<nximg*nyimg;i++) map[i] = 0;
	postype pixsize = Lx/nximg;
	postype xmin,ymin,xmax,ymax;
	xmin = KH_XMIN(simpar);
	ymin = KH_YMIN(simpar);
	xmax = KH_XMAX(simpar);
	ymax = KH_YMAX(simpar);
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
								nearden = tt->mass/meanvol;
							}
						}
						tmp = tmp->next;
					}

				}
			}
			map[i+nximg*j] = nearden;

		}
	}
	MPI_Reduce(map, img, nximg*nyimg, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);


	if(MYID(simpar)==0){
		char outfile[189]; 
		sprintf(outfile,"khmap.%.6d.ppm", icount); 
		colorizeit(5, img, nximg,nyimg,"kh.sao", outfile);
	}
	free(cells);
 	free(map); free(img);
	return 0;
}
int maketscmap(SimParameters *simpar, int icount){ 
	treevoroparticletype *bp = VORO_TBP(simpar);
	int np = VORO_NP(simpar);
	int nx = NX(simpar);
	int ny = NY(simpar);
	int nximg,nyimg;
	nximg = nx;
	nyimg = ny;
	float Lx = KH_SIMBOX(simpar).x.max - KH_SIMBOX(simpar).x.min;
	float Ly = KH_SIMBOX(simpar).y.max - KH_SIMBOX(simpar).y.min;
	int i;
	float *map = (float*)malloc(sizeof(float)*nximg*nyimg);
	float *img = (float*)malloc(sizeof(float)*nximg*nyimg);
	float rscale = (float)nximg/(float)nx;
	float cellsize = Lx/nximg;
	postype mscale = (nximg/Lx)*(nyimg/Ly);
	for(i=0;i<nximg*nyimg;i++) map[i] = 0;


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
	MPI_Reduce(map, img, nximg*nyimg, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	if(MYID(simpar)==0){
		char outfile[189]; 
		sprintf(outfile,"khmap.%.6d.ppm", icount); 
		colorizeit(5,img, nximg,nyimg,"kh.sao", outfile);
	}

	free(map); free(img);

	return 0;
}

treevoroparticletype *mkinitial(SimParameters *simpar, int *mp){
	int i,j,k;
	treevoroparticletype *res;
	postype rho1 = KH_Rho1(simpar), rho2 = KH_Rho2(simpar);
	postype deltarho = 0.5*(rho1-rho2);
	postype deltay = KH_Deltay(simpar);
	postype dvy0 = KH_Vperturb(simpar);
	postype U1 = KH_Vel1(simpar);
	postype U2 = KH_Vel2(simpar);
	postype Um = 0.5*(U1-U2);
	postype xmin,ymin,xmax,ymax;
	xmin = KH_XMIN(simpar);
	ymin = KH_YMIN(simpar);
	xmax = KH_XMAX(simpar);
	ymax = KH_YMAX(simpar);
	int nx = NX(simpar);
	int ny = NY(simpar);
	float Lx = KH_SIMBOX(simpar).x.max - KH_SIMBOX(simpar).x.min;
	float Ly = KH_SIMBOX(simpar).y.max - KH_SIMBOX(simpar).y.min;

	res = (treevoroparticletype*)malloc(sizeof(treevoroparticletype)*nx*ny);
	postype meanvol = Lx*Ly/nx/ny;
	size_t np = 0;
	for(j=0;j<ny;j++){
		postype rho,vx;
		postype y = (postype)(j+0.5)*Ly/(postype)ny;
		char iregion;
		if(y>=0 && y < Ly*0.25) {
			rho = rho1-deltarho*exp((y-Ly*0.25)/deltay);
			vx = U1 - Um*exp((y-Ly*0.25)/deltay);
			iregion = 0;
		}
		else if(y>=Ly*0.25 && y < Ly*0.5) {
			rho = rho2+deltarho*exp((Ly*0.25-y)/deltay);
			vx = U2 + Um*exp((Ly*0.25-y)/deltay);
			iregion = 1;
		}
		else if(y>=Ly*0.5 && y < Ly*0.75) {
			rho = rho2+deltarho*exp((y-Ly*0.75)/deltay);
			vx =  U2 + Um*exp((y-Ly*0.75)/deltay);
			iregion = 1;
		}
		else if(y>=Ly*0.75 && y < Ly*1) {
			rho = rho1-deltarho*exp((Ly*0.75-y)/deltay);
			vx = U1 - Um*exp((Ly*0.75-y)/deltay);
			iregion = 2;
		}
		for(i=0;i<nx;i++){
			postype x = (postype)(i+0.5)*Lx/(postype)nx;
			postype vy = dvy0*sin(4.*M_PI*x );
			size_t indx = i+nx*j;
			if(x>=xmin && x < xmax && y>=ymin && y < ymax){
				res[np].u4if.indx  = indx;
				res[np].u4if.Flag[ENDIAN_OFFSET]  = iregion;
				res[np].x  = x;
				res[np].y  = y;
				res[np].vx  = vx;
				res[np].vy  = vy;
				res[np].mass  = rho*meanvol;
				res[np].den  = rho;
				res[np].pressure  = KH_Pressure(simpar);
				res[np].ie  = res[indx].pressure *meanvol/(Gamma-1);
				res[np].ke  = 0.5*(res[indx].vx*res[indx].vx + res[indx].vy*res[indx].vy);
				res[np].csound  = sqrt(Gamma*res[indx].pressure/res[indx].den);
				res[np].w2  = (Lx/nx*0.1)*(Lx/nx*0.1);
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



	MkLinkedList(simpar);
	DEBUGPRINT("P%d passed MkLinkedList\n", MYID(simpar));
	void initfindmass(SimParameters *);
	initfindmass(simpar);
	DEBUGPRINT("P%d passed initfindmass\n", MYID(simpar));

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


void MkLinkedList(SimParameters *simpar){
	treevoroparticletype *bp = VORO_TBP(simpar);
	int np = VORO_NP(simpar);
	int i,j,k;
	postype xmin,ymin,xmax,ymax;
	int mx,my;
	postype pwidth = KH_GridSize(simpar);

	BASICCELL_CELLWIDTH(simpar) = pwidth;
	BASICCELL_INVCELLWIDTH(simpar) = 1./pwidth;


	xmin = KH_XMIN(simpar)-pwidth;
	ymin = KH_YMIN(simpar)-pwidth;
	xmax = KH_XMAX(simpar)+pwidth;
	ymax = KH_YMAX(simpar)+pwidth;

	mx = ceil((xmax-xmin)/pwidth);
	my = ceil((ymax-ymin)/pwidth);

	BASICCELL_MX(simpar) = mx;
	BASICCELL_MY(simpar) = my;

	cells = VORO_BASICCELL(simpar)= (CellType*)realloc(VORO_BASICCELL(simpar),sizeof(CellType)*mx*my);


	postype cellsize = KH_GridSize(simpar);

	for(i=0;i<mx*my;i++) {
		cells[i].link = NULL;
		cells[i].nmem = 0;
	}

	/*
	DEBUGPRINT("P%d has np= %ld: %d %d: %g : %g %g\n", MYID(simpar), VORO_NP(simpar),
			mx,my, cellsize, xmin,ymin);
			*/

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
//	DEBUGPRINT("P%d has %g %g : %g %g\n", MYID(simpar), KH_XMIN(simpar), KH_YMIN(simpar),
//			KH_XMAX(simpar), KH_YMAX(simpar));
#endif
	KH_TreeAllParticlePadding(simpar, pwidth);
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
	/*
#ifdef DEBUG
	DEBUGPRINT("P%d now after padding: %g %g : %g %g\n", MYID(simpar),Xmin,Ymin,Xmax,Ymax);
#endif
*/
}

Voro2D_point *Voro2D_FindNeighbor(SimParameters *simpar,
		int ix, int iy, int *nneigh, treevoroparticletype *bp){
	int i,j,k;
	int np,mp;
	Voro2D_point *neigh;
	int mx,my;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	cells = VORO_BASICCELL(simpar);
	float Lx = KH_SIMBOX(simpar).x.max - KH_SIMBOX(simpar).x.min;
	float Ly = KH_SIMBOX(simpar).y.max - KH_SIMBOX(simpar).y.min;

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
				neigh[np].w2 = tt->w2;
				neigh[np].bp = (void*)tt;
				np++;
				tmp = tmp->next;
			}

		}
	}
	*nneigh = np;
	return neigh;
}

treevoroparticletype *Voro2D_FindCellBP(SimParameters *simpar,
		int ix, int iy, int *mp, treevoroparticletype *bp){
	int i,j,k;
	int np;
	treevoroparticletype *res;
	int mx,my;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	cells = VORO_BASICCELL(simpar);
	/*
	postype Lx = KH_SIMBOX(simpar).x.max - KH_SIMBOX(simpar).x.min;
	postype Ly = KH_SIMBOX(simpar).y.max - KH_SIMBOX(simpar).y.min;
	*/


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
			res[np] = *tt;
			res[np].bp = tt;
			/*
			res[np].u4if.indx = tt-bp;
			res[np].x += iflag *Lx;
			res[np].y += jflag *Ly;
			*/
			np++;
		}
		tmp = tmp->next;
	}
	*mp = np;
	return res;
}
