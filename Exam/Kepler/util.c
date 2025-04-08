#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "voro.h"
#include "exam2d.h"
#include "kp.h"
#include "color.h"




static CellType *cells=NULL;

//postype cellsize = Ly/Ny;

//void colorizeit(float *, int, int, char *);

int kp_makemap(SimParameters *simpar, int icount){
	postype cellsize = KP_GridSize(simpar);
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
	xmin = KP_XMIN(simpar);
	ymin = KP_YMIN(simpar);
	xmax = KP_XMAX(simpar);
	ymax = KP_YMAX(simpar);
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
		sprintf(outfile,"kpmap.%.6d.ppm", icount); 
		colorizeit(7,img2, nximg,nyimg,"kp.sao", outfile, 0., 2.);
		free(img2);
	}
	free(cells);
 	free(map); free(img);
	free(ndist);
	return 0;
}
int kp_maketscmap(SimParameters *simpar, int icount){
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
		sprintf(outfile,"kpmap.%.6d.ppm", icount);
		colorizeit(7,tmap, nximg,nyimg,"kp.sao", outfile, 0., 2.);
	}
	free(map);
	free(tmap);

	return 0;
}
treevorork4particletype *kp_mkinitial(SimParameters *simpar, int *mp){
	int i,j,k;
	treevorork4particletype *res;
	postype xmin,ymin,xmax,ymax;
	xmin = KP_XMIN(simpar);
	ymin = KP_YMIN(simpar);
	xmax = KP_XMAX(simpar);
	ymax = KP_YMAX(simpar);
	int nx = NX(simpar);
	int ny = NY(simpar);
	float Lx = SIMBOX(simpar).x.max;
	float Ly = SIMBOX(simpar).y.max;
	postype ycen = 0.5*Ly;


	DEBUGPRINT("P%d has initial set nx/y= %d %d Lx/y= %g %g rmin= %g %g rmax= %g %g\n",
			MYID(simpar), nx,ny,Lx,Ly,xmin,ymin,xmax,ymax);

	int mx,my;
	postype lx,ly;
	postype ratio;

	if(MYID(simpar) ==0){
		postype t;
		char infile[190];
		sprintf(infile,"kpout.000000.dat");
        FILE *fp = fopen(infile,"r");
		if(fp == NULL){
			fprintf(stderr,"Error in reading the initial data %s\n", infile);
			exit(9);
		}
        int np;
        fread(&np, sizeof(int), 1, fp);
        VORO_TNP(simpar) = VORO_NP(simpar) = np;

		treevoroparticletype *bp  = (treevoroparticletype*)malloc(sizeof(treevoroparticletype)*np);
        fread(bp,sizeof(treevoroparticletype), np,fp);

        fread(&t,sizeof(postype), 1,fp);
        fread(&mx,sizeof(int), 1,fp);
        fread(&my,sizeof(int), 1,fp);
        fread(&lx,sizeof(postype), 1,fp);
        fread(&ly,sizeof(postype), 1,fp);

		ratio = Lx/lx;

        VORORK4_TBP(simpar) = (treevorork4particletype*)
			malloc(sizeof(treevorork4particletype)*np);
		treevorork4particletype *pp = VORORK4_TBP(simpar);
		for(i=0;i<np;i++){
			pp[i].x = ratio * bp[i].x;
			pp[i].y = ratio * bp[i].y;
			pp[i].u4if.indx  = bp[i].u4if.indx;
		}
		DEBUGPRINT("P%d read the input initial data: ratio= %g %g %g\n", MYID(simpar), ratio, bp[0].x, bp[0].y);
		free(bp);
        fclose(fp);
    }
    else {
        VORO_NP(simpar) = 0;
        VORORK4_TBP(simpar) = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*100);;

    }
    migrateTreeVoroParticles(simpar);
    MPI_Barrier(MPI_COMM(simpar));

	res = VORORK4_TBP(simpar);
	postype meanvol = Lx*Ly/nx/ny;
	postype cx,cy;
	cx = 0.5*Lx;
	cy = 0.5*Ly;
	size_t np = VORO_NP(simpar);
	postype eps2 = KP_EPS(simpar)*KP_EPS(simpar);
	for(i=0;i<VORO_NP(simpar);i++){
		int iregion;
		postype xx = res[i].x;
		postype yy = res[i].y;
		postype x = xx-cx;
		postype y = yy-cy;
		postype R = sqrt(x*x + y*y);
		postype rho;
		if(R<0.5L) rho = 0.01L+ pow(R/0.5L,3.L);
		else if(R<2.L) rho = 0.01L+1.L;
		else rho = 0.01L+pow(1.L+(R-2.L)/0.1L, -3.L);
		if(R<2.5L) iregion = 1; /* evolve only for R < Lx/2 */
		else iregion = 0;
		postype vamp = R*pow(R*R+eps2,-0.75L);
		postype vx,vy;
		if(R==0){
			vx = vy = 0;
		}
		else {
			vx = -y/R*vamp;
			vy =  x/R*vamp;
		}
		res[i].u4if.Flag[ENDIAN_OFFSET] = iregion;
		res[i].x = xx;
		res[i].y = yy;
		res[i].vx = vx;
		res[i].vy = vy;
		res[i].mass = rho*meanvol;
		res[i].volume = meanvol;
		res[i].den = rho;
		res[i].pressure = 1.E-6;
		res[i].ie = res[i].pressure *meanvol/(Gamma-1);
		res[i].ke = 0.5*(vx*vx+vy*vy);
		res[i].csound = sqrt(Gamma*res[i].pressure/rho);
	}

	int nbp = np = VORO_NP(simpar);
	*mp = nbp;
	for(i=0;i<np;i++) UNSET_P_FLAG(simpar, VORORK4, i, BoundaryGhostflag);
	DEBUGPRINT("P%d has np = %ld in box %g %g : %g %g\n", MYID(simpar), np, xmin,ymin,xmax,ymax);

	xmin = ymin = 1.e20;
	xmax = ymax = -1.e20;
	for(i=0;i<nbp;i++){
		xmin = MIN(xmin, res[i].x);
		xmax = MAX(xmax, res[i].x);
		ymin = MIN(ymin, res[i].y);
		ymax = MAX(ymax, res[i].y);
	}
	DEBUGPRINT("P%d has x/ymin %g %g  and max %g %g\n", MYID(simpar), xmin,ymin,xmax,ymax);



	/*
	kp_MkLinkedList(simpar);
	void kp_initfindmass(SimParameters *);
	kp_initfindmass(simpar);
	free(VORO_TBPP(simpar));
	*/

	for(i=0;i<nbp;i++){
		res[i].ie = res[i].pressure*res[i].volume/(Gamma-1);
		res[i].ke = Half*res[i].mass*( res[i].vx*res[i].vx + res[i].vy*res[i].vy);
		res[i].te = res[i].ie + res[i].ke;
		res[i].csound = sqrt(Gamma*res[i].pressure/res[i].den);
	}

	VORO_NPAD(simpar) = 0;

	if(0){
		postype csmin, csmax;
		int imax=0;
		csmax = 0;
		csmin = 1.e20;
		for(i=0;i<VORO_NP(simpar);i++){
			if(res[i].csound > csmax) imax = i;
			csmin = MIN(csmin,res[i].csound);
			csmax = MAX(csmax,res[i].csound);
		}
		DEBUGPRINT("P%d has min max sound speed = %g %g, %g %g %g with np= %ld\n",
				MYID(simpar), csmin,csmax, res[imax].pressure, res[imax].den,
				res[imax].volume, VORO_NP(simpar)); exit(9);
	}


	return res;
}

treevorork4particletype *old_kp_mkinitial(SimParameters *simpar, int *mp){
	int i,j,k;
	treevorork4particletype *res;
	postype xmin,ymin,xmax,ymax;
	xmin = KP_XMIN(simpar);
	ymin = KP_YMIN(simpar);
	xmax = KP_XMAX(simpar);
	ymax = KP_YMAX(simpar);
	int nx = NX(simpar);
	int ny = NY(simpar);
	float Lx = SIMBOX(simpar).x.max;
	float Ly = SIMBOX(simpar).y.max;
	postype ycen = 0.5*Ly;


	DEBUGPRINT("P%d has initial set nx/y= %d %d Lx/y= %g %g rmin= %g %g rmax= %g %g\n",
			MYID(simpar), nx,ny,Lx,Ly,xmin,ymin,xmax,ymax);

	res = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*nx*ny);
	postype meanvol = Lx*Ly/nx/ny;
	postype cx,cy;
	cx = 0.5*Lx;
	cy = 0.5*Ly;
	size_t np = 0;
	postype eps2 = KP_EPS(simpar)*KP_EPS(simpar);
	for(j=0;j<ny;j++){
		postype rho;
		postype yy = (postype)(j)*Ly/(postype)ny;
		postype y = yy-cy;
		char iregion=0;
		for(i=0;i<nx;i++){
			postype xx = (postype)(i)*Lx/(postype)nx;
			if(xx >=xmin && xx < xmax && yy>=ymin && yy<ymax){
				postype x = xx-cx;
				size_t indx = i+nx*j;
				postype R = sqrt(x*x + y*y);
				if(R< 0.5){
					rho = 0.01 + pow(R/0.5,3.L);
				}
				else if (R < 2){
					rho = 0.01 + 1;
				}
				else {
					rho = 0.01 + pow(1 + (R-2)/0.1,-3.L);
				}
				if(R<Lx*0.5) iregion = 1; /* evolution only for R<Lx/2*/
				else iregion = 0;
				postype vamp = R*pow(R*R+eps2,-0.75L);
				postype vx,vy;
				if(R==0) {
					vx = vy = 0;
				}
				else {
					vx = -y/R*vamp;
					vy =  x/R*vamp;
				}
				res[np].u4if.indx  = indx;
				res[np].u4if.Flag[ENDIAN_OFFSET]  = iregion;
				res[np].x  = xx;
				res[np].y  = yy;
				res[np].vx  = vx;
				res[np].vy  = vy;
				res[np].mass  = rho*meanvol;
				res[np].den  = rho;
				res[np].pressure  = 1.E-6;
				res[np].ie  = res[indx].pressure *meanvol/(Gamma-1);
				res[np].ke  = 0.5*(res[indx].vx*res[indx].vx + res[indx].vy*res[indx].vy);
				res[np].csound  = sqrt(Gamma*res[np].pressure/res[np].den);
				np++;
			}
		}
	}
	res = (treevorork4particletype*)realloc(res, sizeof(treevorork4particletype)*np);
	int nbp = np;
	*mp = nbp;
	VORORK4_TBP(simpar) = res;
	VORORK4_BP(simpar) = (vorork4particletype*)res;
	VORO_NP(simpar) = nbp;
	for(i=0;i<np;i++) UNSET_P_FLAG(simpar, VORO, i, BoundaryGhostflag);
	DEBUGPRINT("P%d has np = %ld in box %g %g : %g %g\n", MYID(simpar), np, xmin,ymin,xmax,ymax);

	xmin = ymin = 1.e20;
	xmax = ymax = -1.e20;
	for(i=0;i<nbp;i++){
		xmin = MIN(xmin, res[i].x);
		xmax = MAX(xmax, res[i].x);
		ymin = MIN(ymin, res[i].y);
		ymax = MAX(ymax, res[i].y);
	}
	DEBUGPRINT("P%d has x/ymin %g %g  and max %g %g\n", MYID(simpar), xmin,ymin,xmax,ymax);



	/*
	kp_MkLinkedList(simpar);
	void kp_initfindmass(SimParameters *);
	kp_initfindmass(simpar);
	free(VORO_TBPP(simpar));
	*/

	for(i=0;i<nbp;i++){
		res[i].ie = res[i].pressure*res[i].volume/(Gamma-1);
		res[i].ke = Half*res[i].mass*( res[i].vx*res[i].vx + res[i].vy*res[i].vy);
		res[i].te = res[i].ie + res[i].ke;
		res[i].csound = sqrt(Gamma*res[i].pressure/res[i].den);
	}
	VORO_NPAD(simpar) = 0;


	return res;
}


void kp_MkLinkedList(SimParameters *simpar){
	treevorork4particletype *bp = VORORK4_TBP(simpar);
	int np = VORO_NP(simpar);
	int i,j,k;
	postype xmin,ymin,xmax,ymax;
	int mx,my;
	postype pwidth = KP_GridSize(simpar);

	BASICCELL_CELLWIDTH(simpar) = pwidth;
	BASICCELL_INVCELLWIDTH(simpar) = 1./pwidth;


	xmin = KP_XMIN(simpar)-pwidth;
	ymin = KP_YMIN(simpar)-pwidth;
	xmax = KP_XMAX(simpar)+pwidth;
	ymax = KP_YMAX(simpar)+pwidth;

	mx = ceil((xmax-xmin)/pwidth);
	my = ceil((ymax-ymin)/pwidth);

	BASICCELL_MX(simpar) = mx;
	BASICCELL_MY(simpar) = my;

	cells = VORORK4_BASICCELL(simpar)= 
		(CellType*)realloc(VORORK4_BASICCELL(simpar),sizeof(CellType)*mx*my);


	postype cellsize = KP_GridSize(simpar);

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
			DEBUGPRINT("error detected: %ld : %d %d : %g %g %g %g\n", 
					index, mx,my, bp[i].x, bp[i].y, xmin,ymin);
		}
		struct linkedlisttype *tmp = cells[index].link;
		cells[index].link = (struct linkedlisttype*)(bp+i);
		cells[index].nmem ++;
		bp[i].next = tmp;
	}
	DEBUGPRINT("P%d has volume %g %g : %g %g & test values %ld %ld\n", 
			MYID(simpar), KP_XMIN(simpar), KP_YMIN(simpar),
			KP_XMAX(simpar), KP_YMAX(simpar),
			cells[1+ mx*1].nmem, cells[2+mx*1].nmem);
	paddingVoroParticles(simpar, pwidth);
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
	DEBUGPRINT("P%d now after padding: %g %g : %g %g : cell= %ld %ld\n", 
			MYID(simpar),Xmin,Ymin,Xmax,Ymax,
			cells[0].nmem, cells[1].nmem);
}

Voro2D_point *kp_Voro2D_FindNeighbor(SimParameters *simpar,
		int ix, int iy, int *nneigh, treevorork4particletype *bp){
	int i,j,k;
	int np,mp;
	Voro2D_point *neigh;
	int mx,my;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	cells = VORORK4_BASICCELL(simpar);
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
				treevorork4particletype *tt = (treevorork4particletype*)tmp;
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

treevorork4particletype *kp_Voro2D_FindCellBP(SimParameters *simpar,
		int ix, int iy, int *mp, treevorork4particletype *bp){
	int i,j,k;
	int np;
	treevorork4particletype *res;
	int mx,my;
	mx = BASICCELL_MX(simpar);
	my = BASICCELL_MY(simpar);
	cells = VORORK4_BASICCELL(simpar);
//	postype Ly = SIMBOX(simpar).y.max;


	np = cells[ix+mx*iy].nmem;
	res = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*np);
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
