#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "voro.h"
#include "rt.h"
#include "nnost.h"
#include "gnnost.h"
#include "exam.h"
#include "exam2d.h"
#include "color.h"





//postype cellsize = Ly/Ny;

//void colorizeit(float *, int, int, char *);

int rt_makemap(SimParameters *simpar, int icount){
	postype cellsize = RT_GridSize(simpar);
	postype Lx = SIMBOX(simpar).x.max - SIMBOX(simpar).x.min;
	postype Ly = SIMBOX(simpar).y.max - SIMBOX(simpar).y.min;
	int nximg = NX(simpar);
	int nyimg = NY(simpar);
	float *map = (float*)my_malloc(sizeof(float)*nximg*nyimg);
	float *img = (float*)my_malloc(sizeof(float)*nximg*nyimg);
	int i,j;
	int ii,jj;
	for(i=0;i<nximg*nyimg;i++) map[i] = 0;
	postype pixsize = Lx/nximg;
	postype xmin,ymin,xmax,ymax;
	xmin = RT_XMIN(simpar);
	ymin = RT_YMIN(simpar);
	xmax = RT_XMAX(simpar);
	ymax = RT_YMAX(simpar);
	int mx,my;
	mx = ceil((xmax-xmin)/cellsize);
	my = ceil((ymax-ymin)/cellsize);
	CellType *cells;
	cells = (CellType*)my_malloc(sizeof(CellType)*mx*my);
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
			map[i+nximg*(nyimg-j-1)] = nearden;

		}
	}

	MPI_Reduce(map, img, nximg*nyimg, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
    if(MYID(simpar)==0){
        char outfile[189];
		char outsao[180];
        sprintf(outfile,"rtmap.%.6d.ppm", icount);
        strcpy(outsao,"rt.sao");
        colorizeit(5, img, nximg,nyimg,outsao, outfile);
    }
	MPI_Barrier(MPI_COMM_WORLD);
    if(cells != NULL) my_free(cells);
	if(img != NULL) my_free(img);
    if(map != NULL) my_free(map); 

	return 0;
}
postype getRT_Den(postype rho1, postype rho2, postype deltay, postype ycen,postype y){
	postype rho = rho1+(rho2-rho1)/(1+exp(-(y-ycen)/deltay));
	return rho;
}


static postype rt_rho1, rt_rho2, rt_ycen, rt_deltay;
float qsimp_getrho(float y){
    return rt_rho1 + (rt_rho2-rt_rho1) / (1.0 + exp(-(y - rt_ycen) / rt_deltay));
}
#include "nrutil.h"
#include "nr.h"
postype nr_rt_pressure(SimParameters *simpar,
		postype rho1, postype rho2, postype ycen, postype y) {
	rt_rho1 = rho1;
	rt_rho2 = rho2;
	rt_ycen = ycen;
	rt_deltay = RT_Deltay(simpar);
	postype integral;
	if(y == ycen) return RT_Phalf(simpar);
	else if (y > ycen){
		integral = qsimp(qsimp_getrho, ycen, y);
	}
	else {
		integral = -qsimp(qsimp_getrho, y, ycen);
	}
	return RT_Phalf(simpar)+ GAS_ACCY(simpar)* (( y - ycen) + integral);
}

// Density function: rho(y) = 1 + 1 / (1 + exp(-(y-0.5)/Delta))
postype rt_getrho(postype rho1, postype rho2, postype deltay, postype ycen, postype y) {
    return rho1 + (rho2-rho1) / (1.0 + exp(-(y - ycen) / deltay));
}

// Numerical integration using the trapezoidal rule
postype integrate_rt_rho(postype rho1, postype rho2, postype deltay, postype ycen, 
		postype y0, postype y1, int n) {
    postype h = (y1 - y0) / n;
    postype sum = 0.5 * (rt_getrho(rho1, rho2, deltay, ycen, y0) + rt_getrho(rho1, rho2, deltay,ycen, y1));
    for (int i = 1; i < n; ++i) {
        postype y = y0 + i * h;
        sum += rt_getrho(rho1, rho2, deltay, ycen,y);
    }
    return h * sum;
}

// Pressure function p(y)
postype rt_pressure(SimParameters *simpar,
		postype rho1, postype rho2, postype ycen, postype y) {
	int nsteps = 20000;
	postype Pcen = RT_Phalf(simpar);
	postype g = GAS_ACCY(simpar);
	postype deltay = RT_Deltay(simpar);
    if (y == ycen) return Pcen;
    postype integral;
    if (y > ycen) {
        // Integrate from 0.5 to y
        integral = integrate_rt_rho(rho1, rho2, deltay, ycen, ycen, y, nsteps);
    } else {
        // Integrate from y to 0.5 and take the negative
        integral = -integrate_rt_rho(rho1, rho2, deltay, ycen,y, ycen,  nsteps);
    }
    return Pcen + g * ((y - ycen) + integral);
}

/* for hydrostatic equilibrium, the pressure at y is obtained by integral of g(y)*rho(y) */
postype  getRTpressure(SimParameters *simpar, postype rho1, postype rho2, postype ycen, postype y){
	postype deltay = RT_Deltay(simpar);
//	postype pressure = RT_Phalf(simpar) + GAS_ACCY(simpar)*getRT_Den(rho1,rho2, deltay, ycen, y)*(y-ycen);
	postype Const = RT_Phalf(simpar)-GAS_ACCY(simpar)*(rho2/2+(rho2-rho1)*deltay*log(2.));
	postype pressure = GAS_ACCY(simpar)*(rho2*y+(rho2-rho1)*deltay*log(1+exp(-(y-ycen)/deltay))) + Const;
	return pressure;
}
treevorork4particletype *rt_mkinitial(SimParameters *simpar, int *mp){
	int i,j,k;
	treevorork4particletype *res;
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
	postype dmean = Lx/nx;
	postype ycen = 0.5*Ly;

	GAS_dMean(simpar) = dmean;

	postype Gamma = GAS_GAMMA(simpar);


	DEBUGPRINT("P%d has initial set nx/y= %d %d Lx/y= %g %g rmin= %g %g rmax= %g %g\n",
			MYID(simpar), nx,ny,Lx,Ly,xmin,ymin,xmax,ymax);
	DEBUGPRINT("P%d has initial set den1/2 = %g %g ycen= %g\n", 
			MYID(simpar), rho1, rho2, ycen);

	res = (treevorork4particletype*)my_malloc(sizeof(treevorork4particletype)*nx*ny);
	postype meanvol = Lx*Ly/(postype)nx/(postype)ny;
	postype Pressure[ny+1];
	postype dy = Ly/ny;
	for(j=ny/2;j<=ny;j++){
		postype y = (postype)(j+0.5)*Ly/ny;
		postype rho = getRT_Den(rho1, rho2, deltay, ycen, y);
		postype dr = dy;
		if(j==ny/2) {
			Pressure[j] = RT_Phalf(simpar);
		}
		else {
			Pressure[j] = Pressure[j-1] + GAS_ACCY(simpar)* rho * dr;
		}
	}
//	GAS_invw2Scale(simpar) = 1.L/(RT_Phalf(simpar)/getRT_Den(rho1,rho2,deltay,ycen,ycen));
	GAS_invw2Scale(simpar) = 1.L/(RT_Phalf(simpar)*meanvol);
	for(j=ny/2-1;j>=0;j--){
		postype y = (postype)(j+0.5)*Ly/ny;
        postype rho = getRT_Den(rho1, rho2, deltay, ycen, y);
        postype dr = -dy;
        Pressure[j] = Pressure[j+1] + GAS_ACCY(simpar)* rho * dr;

	}
	size_t np = 0;
	for(j=0;j<ny;j++){
		postype rho;
		postype y = (postype)(j+0.5)*Ly/(postype)ny;
		char iregion;
		rho = getRT_Den(rho1, rho2, deltay, ycen, y); 
		if(j==0 || j == ny-1) iregion = -1;
		else iregion = 0;
		postype vx = 0;
		for(i=0;i<nx;i++){
			postype vy;
			postype x = (postype)(i+0.5)*Lx/(postype)nx;

			if(j%2 ==0) x = (postype)(i+0.25)*Lx/(postype)nx;
			else x = (postype)(i+0.75L)*Lx/(postype)nx;
			/*
			x = (postype)(i+0.5L)*Lx/(postype)nx;
			*/

			size_t indx = i+nx*j;
			if(y>=0.3*Ly && y < 0.7*Ly){
				vy = dvy0*(1+cos(8*M_PI*(x+0.5*Lx)))*(1+cos(5*M_PI*(y-ycen)));
			}
			else{
				vy = 0;
			}
			if(x>=xmin && x < xmax && y>=ymin && y < ymax){
				UNSET_FLAG(res+np, Wallflag);
				res[np].u4if.indx  = indx;
//				if(j==0 || j==1 || j == ny-1 || j == ny-2) SET_FLAG(res+np, Wallflag);				
//				if(y <=0.1 || y >=0.9) SET_FLAG(res+np, Wallflag);
//				else UNSET_FLAG(res+np, Wallflag);				
				res[np].x  = x;
				res[np].y  = y;
				res[np].vx  = vx;
				res[np].vy  = vy;

				res[np].mass  = rho*meanvol;
				res[np].den  = rho;
				res[np].pressure  = 0.5*(Pressure[j] +Pressure[j+1]);
				res[np].ie  = res[np].pressure *meanvol/(Gamma-1);
				/*
				res[np].ke  = 0.5*(res[np].vx*res[np].vx + res[np].vy*res[np].vy);
				*/
				res[np].csound  = 1;
				res[np].w2  = (Lx/nx*GAS_Kappa(simpar))*(Lx/nx*GAS_Kappa(simpar));
                res[np].w2old  = (Lx/nx*GAS_Kappa(simpar))*(Lx/nx*GAS_Kappa(simpar));
                res[np].die  = 0;
				if(GAS_Kappa(simpar) <0) {
		            res[np].w2 = - GAS_Kappa(simpar);
       		     res[np].w2old = - GAS_Kappa(simpar);
		        }
		        else {
			     res[np].avgNeighboringPressure = res[np].pressure;
       		     res[np].w2 = getw2forHydroParticle(simpar, (res+np),1);
           		 res[np].w2old = res[np].w2;
		            res[np].w2ceil = res[np].w2;
       		 }
				np++;
			}
		}
	}
	DEBUGPRINT("P%d has pixels with mean volume = %g\n", MYID(simpar), meanvol);
	DEBUGPRINT("P%d has p%ld pressure = %g\n", MYID(simpar),PINDX(res),res[0].pressure);
	res = (treevorork4particletype*)realloc(res, sizeof(treevorork4particletype)*np);
	int nbp = np;
	*mp = nbp;
	VORORK4_TBP(simpar) = res;
	VORORK4_BP(simpar) = (vorork4particletype*)res;
	VORO_NP(simpar) = nbp;
	for(i=0;i<np;i++) UNSET_P_FLAG(simpar, VORORK4, i, BoundaryGhostflag);
	DEBUGPRINT("P%d has np = %ld in box %g %g : %g %g\n", MYID(simpar), np, xmin,ymin,xmax,ymax);

	migrateTreeVorork4Particles(simpar);

	if(0){
		treevorork4particletype *bp = VORORK4_TBP(simpar);
		for(i=0;i<NID(simpar);i++){
			if(MYID(simpar) == i){
				for(j=0;j<np;j++) printf("P%d has p%d %g %g\n", MYID(simpar), j,bp[j].x,bp[j].y);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		exit(9);
	}

	/*
	rt_MkLinkedList(simpar);
	void rt_initfindmass(SimParameters *);
	rt_initfindmass(simpar);
	*/
	// measure the cell volume for each particle
	/*
    void kh_findVol(SimParameters *);
    kh_findVol(simpar);
	*/
	exam2dUpdateVol(simpar,paddingTreeVorork4Particles,searchCellRk4Neighbors2D,findCellRk4BP2D,
			mkLinkedList2D_rt);
    DEBUGPRINT("P%d has volume vol(0) %g\n", MYID(simpar), res[0].volume);

    res = VORORK4_TBP(simpar);
    nbp = VORO_NP(simpar);

	for(i=0;i<nbp;i++){
		res[i].mass = res[i].den * res[i].volume;
		res[i].ie = res[i].pressure*res[i].volume/(Gamma-1);
		/*
		res[i].ke = Half*res[i].mass*( res[i].vx*res[i].vx + res[i].vy*res[i].vy);
		res[i].te = res[i].ie + res[i].ke;
		*/
		res[i].csound = sqrt(Gamma*res[i].pressure/res[i].den);
	}
//	my_free(VORORK4_TBPP(simpar));
	VORO_NPAD(simpar) = 0;

	// update w2ceil and w2
    if(GAS_Kappa(simpar)>0) det2d_dpqRK4(simpar,paddingTreeVorork4Particles);
//    DEBUGPRINT("P%d passed det2d_dpqRK4\n", MYID(simpar));


	return res;
}


