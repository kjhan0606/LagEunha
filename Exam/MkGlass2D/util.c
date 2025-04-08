#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "voro.h"
//#include "ost.h"
#include "nnost.h"
#include "gnnost.h"
#include "exam.h"
#include "exam2d.h"
#include "gl2d.h"
#include "color.h"




static CellType *cells=NULL;

//postype cellsize = Ly/Ny;

//void colorizeit(float *, int, int, char *);

int gl2d_makemap(SimParameters *simpar, int icount){
	postype cellsize = GL2D_GridSize(simpar);
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
	xmin = GL2D_XMIN(simpar);
	ymin = GL2D_YMIN(simpar);
	xmax = GL2D_XMAX(simpar);
	ymax = GL2D_YMAX(simpar);
	int mx,my;
	mx = ceil((xmax-xmin)/cellsize);
	my = ceil((ymax-ymin)/cellsize);
	cells = (CellType*)malloc(sizeof(CellType)*mx*my);
	for(i=0;i<mx*my;i++){
		cells[i].link = NULL;
		cells[i].nmem = 0;
	}
	treevorork4particletype *bp = VORORK4_TBP(simpar);
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
						treevorork4particletype *tt = (treevorork4particletype*)tmp;
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
		sprintf(outfile,"glmap.%.6d.ppm", icount); 
		colorizeit(7,img2, nximg,nyimg,"mk2d.sao", outfile, 0., 2.);
		free(img2);
	}
	free(cells);
 	free(map); free(img);
	free(ndist);
	return 0;
}
int gl2d_maketscmap(SimParameters *simpar, int icount){
	int nx = NX(simpar);
	int ny = NY(simpar);
	int np = VORO_NP(simpar);
	int nximg, nyimg;
	nximg = nx;
	nyimg = ny;
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	float Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	float Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	float cellsize = Lx/nximg;
	int i;
	float *map = (float*)malloc(sizeof(float)*nximg*nyimg);
	float *tmap = (float*)malloc(sizeof(float)*nximg*nyimg);
	for(i=0;i<nximg*nyimg;i++) map[i] = 0;
	for(i=0;i<nximg*nyimg;i++) tmap[i] = 0;
	postype mscale = (nximg/Lx)*(nyimg/Ly);

	for(i=0;i<np;i++){
		postype pmass = mscale*bp[i].mass;
		postype xp = bp[i].x/cellsize;
		postype yp = bp[i].y/cellsize;
		int nearx = rint(xp);
		int neary = rint(yp);
		int ic = (nearx + nximg) %nximg;
		int jc = (neary + nyimg) %nyimg;
		postype xmin = xp -nearx;
		postype ymin = yp -neary;
		int icc = (nximg+nearx + (int)(copysign(1., xmin))) %nximg;
		int jcc = (nyimg+neary + (int)(copysign(1., ymin))) %nyimg;
		int iccc = (nximg+nearx - (int)(copysign(1., xmin))) %nximg;
		int jccc = (nyimg+neary - (int)(copysign(1., ymin))) %nyimg;
		float xd1 = fabs(xmin);
		float yd1 = fabs(ymin);
		float wx1 = (0.75-xd1*xd1)*pmass;
		float wy1 = (0.75-yd1*yd1);
		float wx3 = 0.5*pmass * (0.25+xd1*(xd1-1.));
		float wx2 = wx3 + pmass * xd1;
		float wy3 = 0.5 * (0.25+yd1*(yd1-1.));
		float wy2 = wy3 + yd1;

		map[ic+nximg*jc] += wx1*wy1;
		map[icc+nximg*jc] += wx2*wy1;
		map[iccc+nximg*jc] += wx3*wy1;
		map[ic+nximg*jcc] += wx1*wy2;
		map[icc+nximg*jcc] += wx2*wy2;
		map[iccc+nximg*jcc] += wx3*wy2;
		map[ic+nximg*jccc] += wx1*wy3;
		map[icc+nximg*jccc] += wx2*wy3;
		map[iccc+nximg*jccc] += wx3*wy3;

	}
	MPI_Reduce(map, tmap, nximg*nyimg, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	if(MYID(simpar)==0){
		char outfile[189];
		sprintf(outfile,"glmap.%.6d.ppm", icount);
		colorizeit(7,tmap, nximg,nyimg,"mk2d.sao", outfile, 0., 2.);
	}
	free(map);
	free(tmap);

	return 0;
}

treevorork4particletype *gl2d_mkinitial(SimParameters *simpar, int *mp){
	int i,j,k;
	treevorork4particletype *res;
	postype xmin,ymin,xmax,ymax;
	xmin = GL2D_XMIN(simpar);
	ymin = GL2D_YMIN(simpar);
	xmax = GL2D_XMAX(simpar);
	ymax = GL2D_YMAX(simpar);
	int nx = NX(simpar);
	int ny = NY(simpar);
	float Lx = SIMBOX(simpar).x.max;
	float Ly = SIMBOX(simpar).y.max;
	postype ycen = 0.5*Ly;

	postype pwidth = BASICCELL_CELLWIDTH(simpar) = GL2D_GridSize(simpar);


	DEBUGPRINT("P%d has initial set nx/y= %d %d Lx/y= %g %g rmin= %g %g rmax= %g %g\n",
			MYID(simpar), nx,ny,Lx,Ly,xmin,ymin,xmax,ymax);

	res = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*nx*ny);
	postype meand = Lx/nx;
	postype meanvol = meand*meand;
	postype cx,cy;
	cx = 0.5*Lx;
	cy = 0.5*Ly;
	size_t np = 0;
	float ran3(long *idum);
	long iseed = -(MYID(simpar)+4);
	postype rho = 1; /* default density */
	for(j=0;j<ny;j++){
		char iregion=1;
		postype y,x;
		y = (postype)(j)*Ly/(postype)ny;
		for(i=0;i<nx;i++){
			postype yy,xx;
			if(j%2 ==0) x = (postype)(i)*Lx/(postype)nx;
			else x = (postype)(i+0.5L)*Lx/(postype)nx;
			if(x >=xmin && x < xmax && y>=ymin && y<ymax){
				xx = x + 0.5*ran3(&iseed)*Lx/(postype)nx;
				yy = y + 0.5*ran3(&iseed)*Ly/(postype)ny;
				size_t indx = i+nx*j;
				postype vx,vy;
				vx = vy = 0;
				res[np].u4if.indx  = indx;
				res[np].u4if.Flag[ENDIAN_OFFSET]  = iregion;
				res[np].x  = xx;
				res[np].y  = yy;
				res[np].vx  = vx;
				res[np].vy  = vy;
				res[np].mass  = rho*meanvol;
				res[np].den  = rho;
				res[np].pressure  = pow(res[np].den, Gamma);
				res[np].w2  = gl2d_w2Measure(simpar,meand, res[np].pressure, res[np].den);
				res[np].ie  = res[indx].pressure *meanvol/(Gamma-1);
				res[np].ke  = 0.5*(res[indx].vx*res[indx].vx + res[indx].vy*res[indx].vy);
				res[np].csound  = sqrt(Gamma*res[np].pressure/res[np].den);
				res[np].volume = meanvol;
				np++;
			}
		}
	}
	res = (treevorork4particletype*)realloc(res, sizeof(treevorork4particletype)*np);


	int nbp = np;
	*mp = nbp;
	// voro-rk4 type tree-structure base particle.
	VORORK4_TBP(simpar) = res;
	// voro-rk4 type base particle.
	VORORK4_BP(simpar) = (vorork4particletype*)res;
	VORO_NP(simpar) = nbp;
	// declare the particles as not being boundary ghost particle
	for(i=0;i<np;i++) UNSET_P_FLAG(simpar, VORO, i, BoundaryGhostflag);
	DEBUGPRINT("P%d has np = %ld in box %g %g : %g %g\n", MYID(simpar), np, xmin,ymin,xmax,ymax);

    migrateTreeVoroParticles(simpar);

	// measure the cell volume for each particle
	void gl2d_findVol(SimParameters *);
	gl2d_findVol(simpar);
	DEBUGPRINT("P%d passed initfindmass\n", MYID(simpar));



	res = VORORK4_TBP(simpar);
	nbp = VORO_NP(simpar);
	xmin = ymin = 1.e20;
	xmax = ymax = -1.e20;
	for(i=0;i<nbp;i++){
		xmin = MIN(xmin, res[i].x);
		xmax = MAX(xmax, res[i].x);
		ymin = MIN(ymin, res[i].y);
		ymax = MAX(ymax, res[i].y);
	}
	DEBUGPRINT("P%d has x/ymin %g %g  and max %g %g\n", MYID(simpar), xmin,ymin,xmax,ymax);

	for(i=0;i<nbp;i++){
		res[i].den  = res[i].mass/res[i].volume;
		res[i].pressure  = pow(res[i].den, Gamma);
		res[i].ie = res[i].pressure*res[i].volume/(Gamma-1);
		res[i].ke = Half*res[i].mass*( res[i].vx*res[i].vx + res[i].vy*res[i].vy);
		res[i].te = res[i].ie + res[i].ke;
		res[i].csound = sqrt(Gamma*res[i].pressure/res[i].den);
		res[i].w2  = gl2d_w2Measure(simpar,meand, res[i].pressure, res[i].den);
	}
	VORO_NPAD(simpar) = 0;

	DEBUGPRINT("P%d has w2 values P= %g w2= %g w2ceil = %g den= %g meanV= %g vol= %g Cs= %g\n", 
			MYID(simpar), res[0].pressure, res[0].w2, res[0].w2ceil, res[0].den,
			meanvol, res[0].volume, res[0].csound);

	return res;
}


void gl2d_MkLinkedList(SimParameters *simpar){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i,j,k;
	postype xmin,ymin,xmax,ymax;
	int mx,my;
	postype pwidth = GL2D_GridSize(simpar);

	BASICCELL_CELLWIDTH(simpar) = pwidth;
	BASICCELL_INVCELLWIDTH(simpar) = 1./pwidth;


	xmin = GL2D_XMIN(simpar)-pwidth;
	ymin = GL2D_YMIN(simpar)-pwidth;
	xmax = GL2D_XMAX(simpar)+pwidth;
	ymax = GL2D_YMAX(simpar)+pwidth;

	mx = ceil((xmax-xmin)/pwidth);
	my = ceil((ymax-ymin)/pwidth);

	BASICCELL_MX(simpar) = mx;
	BASICCELL_MY(simpar) = my;

	cells = VORO_BASICCELL(simpar)= (CellType*)realloc(VORO_BASICCELL(simpar),sizeof(CellType)*mx*my);


	postype cellsize = GL2D_GridSize(simpar);

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
		if(index >=mx*my){
			printf("error detected: %ld : %d %d : %g %g < %g %g : %g %g >\n", 
					index, mx,my, bp[i].x, bp[i].y, xmin,ymin, xmax,ymax);
			exit(0);
		}
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)(bp+i);
		cells[index].nmem ++;
		bp[i].next = tmp;
	}
	DEBUGPRINT("P%d has %g %g : %g %g\n", MYID(simpar), GL2D_XMIN(simpar), GL2D_YMIN(simpar),
			GL2D_XMAX(simpar), GL2D_YMAX(simpar));
	paddingTreeVoroParticles(simpar, pwidth);
	bp = VORORK4_TBPP(simpar);
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

Voro2D_point *gl2d_Voro2D_FindNeighbor(SimParameters *simpar,
		int ix, int iy, int *nneigh){
	int i,j,k;
	int np,mp;
	Voro2D_point *neigh;
	int mx,my;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	cells = VORO_BASICCELL(simpar);

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
				treevorork4particletype *tt = (treevorork4particletype*)tmp;
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

treevorork4particletype *gl2d_Voro2D_FindCellBP(SimParameters *simpar,
		int ix, int iy, int *mp){
	int i,j,k;
	int np;
	treevorork4particletype *res;
	int mx,my;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	cells = VORO_BASICCELL(simpar);
//	postype Ly = SIMBOX(simpar).y.max;


	np = cells[ix+mx*iy].nmem;
	res = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*np);
	np = 0;
	size_t ipixel = ix+mx*iy;
	struct linkedlisttype *tmp = cells[ipixel].link;
	while(tmp){
		// exclude the boundary ghost particles.
		if(!IS_FLAG_ON(tmp,BoundaryGhostflag)){
			treevorork4particletype *tt = (treevorork4particletype*)tmp;
			if(tt->u4if.Flag[ENDIAN_OFFSET]> 0)
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
