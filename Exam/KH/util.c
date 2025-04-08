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
#include "kh.h"
#include "color.h"




static CellType *cells=NULL;

//postype cellsize = Ly/Ny;

//void colorizeit(float *, int, int, char *);

int kh_makemap(SimParameters *simpar, int icount){
	postype cellsize = KH_GridSize(simpar);
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
int kh_maketscmap(SimParameters *simpar, int icount){
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

treevorork4particletype *kh_mkinitial(SimParameters *simpar, int *mp){
	int i,j,k;
	treevorork4particletype *res;
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
	postype dmean = Lx/nx;

    res = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*nx*ny);
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
	void kh_findVol(SimParameters *);
	kh_findVol(simpar);
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
		res[i].w2  = kh_w2Measure(simpar,dmean, res[i].pressure, res[i].den);
	}
	VORO_NPAD(simpar) = 0;

	DEBUGPRINT("P%d has w2 values P= %g w2= %g w2ceil = %g den= %g meanV= %g vol= %g Cs= %g\n", 
			MYID(simpar), res[0].pressure, res[0].w2, res[0].w2ceil, res[0].den,
			meanvol, res[0].volume, res[0].csound);

	return res;
}


void kh_MkLinkedList(SimParameters *simpar){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
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
	DEBUGPRINT("P%d has %g %g : %g %g\n", MYID(simpar), KH_XMIN(simpar), KH_YMIN(simpar),
			KH_XMAX(simpar), KH_YMAX(simpar));
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

Voro2D_point *kh_Voro2D_FindNeighbor(SimParameters *simpar,
		int ix, int iy, int *nneigh){
	int i,j,k;
	int np,mp;
	Voro2D_point *neigh;
	int mx,my;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	cells = VORO_BASICCELL(simpar);
	/*
	float Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	float Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	*/

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

treevorork4particletype *kh_Voro2D_FindCellBP(SimParameters *simpar,
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
