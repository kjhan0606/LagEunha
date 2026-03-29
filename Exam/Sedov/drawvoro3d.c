#include "voro.h"
#include "sedov.h"
//#include "pic2d.h"

/*
CellType *mycell=NULL;

ptype cellsize = Lx/Nx;


void linkedlist(Voro2D_GasParticle *bp, int np){
    int i,j,k;
    if(mycell == NULL) mycell = (CellType*)malloc(sizeof(CellType)*Nx*Ny);


    cellsize = Lx/Nx;

    for(i=0;i<Nx*Ny;i++) {
        mycell[i].bp = NULL;
        mycell[i].np = 0;
    }

    for(i=0;i<np;i++){
        int ix,iy;
        ix = bp[i].x/cellsize;
        iy = bp[i].y/cellsize;
        size_t index = ix+Nx*iy;
        Voro2D_GasParticle *tmp = mycell[index].bp;
        mycell[index].bp = bp+i;
        mycell[index].np ++;
        bp[i].next = tmp;
    }
}
*/




char mycolor[256][20];

int setcolor(char *, int);



int sminitialize(FILE *wp, float xmin, float xmax, float ymin, float ymax){
	float bigxtick,smxtick,bigytick,smytick;
	bigxtick = (xmax-xmin)/5;
	smxtick = (xmax-xmin)/20;
	bigytick = (ymax-ymin)/5;
	smytick = (ymax-ymin)/20;
	int i;
    fprintf(wp,"erase\n");
	for(i=0;i<256;i++){
		sprintf(mycolor[i],"color%d",i);
    	fprintf(wp,"ctype = CTYPE() concat %d\n", setcolor("sedov.sao",i));
	    fprintf(wp,"ctype = CTYPE(STRING) concat '%s'\n", mycolor[i]);
	}
    fprintf(wp,"ctype default\n");
    fprintf(wp,"ltype 0\n");
    fprintf(wp,"location 5000 27000 5000 27000\n");
    fprintf(wp,"limits %g %g %g %g\n", xmin,xmax,ymin,ymax);
    fprintf(wp,"ticksize %g %g %g %g\n", smxtick, bigxtick, smytick, bigytick);
//    fprintf(wp,"box \n");
}

int smfinalize(FILE *wp, ptype time){
	int i,j,k;
//    fprintf(wp,"ctype default \nbox \n");
	for(i=0;i<256;i++){
		fprintf(wp,"del_ctype %s\n", mycolor[i]);
	}
	fclose(wp);
}
int getcolor(int );
int setpolygoncolor(int nlog, float val){
	int res;
	val = val/4.;
	if(nlog>0){
		nlog = exp(nlog);
		float ccdmax = log(1.+nlog);
		val = log((1.0+nlog*val));
		float val2 = (int)(val/ccdmax*255);
		if(val2>255) val2 = 255;
		else if(val2<0) val2 = 0;
		res = (int)val2;
	}
	else {
		float val2 =  (int)(val*255);
		if(val2>255) val2 = 255;
		else if(val2<0) val2 = 0;
		res = val2;
	}
	return res;
}
CellType *cells=NULL;

int LinkedCell(Voro3D_GasParticle *bp, int np){
	int maxnp = 0;
	int i,j,k;
	if(cells == NULL) cells = (CellType*)malloc(sizeof(CellType)*Nx*Ny);
	ptype cellsize = Lx/Nx;
	for(i=0;i<Nx*Ny;i++) {
		cells[i].bp = NULL;
		cells[i].np = 0;
	}
	for(i=0;i<np;i++){
		int ix,iy;
		ix = (bp[i].x/cellsize);
		iy = (bp[i].y/cellsize);
		size_t index = ix+Nx*iy;
		Voro3D_GasParticle *tmp = cells[index].bp;
		cells[index].bp = bp+i;
		cells[index].np ++;
		bp[i].next = tmp;
	}
	for(i=0;i<Nx*Ny;i++){
		maxnp = MAX(maxnp, cells[i].np);
	}
	return maxnp;
}
void Voro2D_FindCellBP(Voro3D_GasParticle *p, int ix, int iy, int *mp, Voro3D_GasParticle *bp){
	int i,j,k;
	int np;

	int iix,iiy;
	iix = (ix+Nx)%Nx;
	iiy = (iy+Ny)%Ny;

	np = cells[iix+Nx*iiy].np;
	np = 0;
	int iflag,jflag;
	iflag = jflag = 0;
	if(iix < 0) iflag = -1;
	else if(iix >= Nx) iflag = 1;
	if(iiy < 0) jflag = -1;
	else if(iiy >= Ny) jflag = 1;
	size_t ipixel = iix+Nx*iiy;
	Voro3D_GasParticle *tmp = cells[ipixel].bp;
	while(tmp){
		p[np] = *tmp;
		p[np].x += iflag *Lx;
		p[np].y += jflag *Ly;
		np++;
		tmp = tmp->next;
	}
	*mp = np;
	return;
}
void Voro2D_FindNeighbor(Voro3D_point *neigh, int ix, int iy, int *nneigh, Voro3D_GasParticle *bp){
	int i,j,k;
	int np,mp;

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
			Voro3D_GasParticle *tmp = cells[ii+Nx*jj].bp;
			while(tmp){
				neigh[np].x = tmp->x + iflag *Lx;
				neigh[np].y = tmp->y + jflag *Ly;
				neigh[np].z = tmp->z;
				neigh[np].id = tmp->id;
				neigh[np].csound = tmp->csound;
				np++;
				tmp = tmp->next;
			}

		}
	}
	*nneigh = np;
	return;
}
void mkcolorbar(FILE *wp, float valmin, float valmax,int nlog){
    float lgvalmin = log10(valmin);
    float lgvalmax = log10(valmax);
    float valstep = ((lgvalmax)-(lgvalmin))/128.;
    fprintf(wp,"location 28000 29500 7000 25000\n");
    fprintf(wp,"limits 0 1 %g %g\n", lgvalmin, lgvalmax);
    fprintf(wp,"ticksize 10 100 -1 10\n");
    int i;
    fprintf(wp,"set x = 0, 1, 1 \n");
    fprintf(wp,"set x = x concat reverse(x) \n");
    for(i=0;i<128;i++){
        float val = valstep*(i)+lgvalmin;
        float rval = pow(10., val);
        int icolor = setpolygoncolor(nlog, rval);
        float val2 = val + valstep;
        fprintf(wp,"set y = %g concat %g concat %g concat %g\n", val, val, val2, val2);
        fprintf(wp,"ctype %s\n", mycolor[icolor]);
        fprintf(wp,"shade 0 x y\n");
        fprintf(wp,"delete y\n");
    }
    fprintf(wp,"ctype default \n box 0 0 0 2\n");
}


int main(int argc,char **argv){
	int i,j,k;
	int np;
	int nlog = 4;
	ptype time,dt;
	FILE *fp = fopen(argv[1],"r");
	FILE *wp = fopen(argv[2],"w");
	ptype xmin,ymin,xmax,ymax;
	xmin = atof(argv[3]);
	xmax = atof(argv[4]);
	ymin = atof(argv[5]);
	ymax = atof(argv[6]);
	ptype zcut = atof(argv[7]);
	sminitialize(wp, xmin,xmax,ymin,ymax);
	fread(&np,sizeof(int),1,fp);
	Voro3D_GasParticle *bp3d = (Voro3D_GasParticle*)malloc(sizeof(Voro3D_GasParticle)*(np));
    fread(bp3d,sizeof(Voro3D_GasParticle), np,fp);
    fread(&time,sizeof(ptype), 1,fp);
    fread(&dt,sizeof(ptype), 1,fp);
    fclose(fp);

	Voro3D_GasParticle *bp = (Voro3D_GasParticle*)malloc(sizeof(Voro3D_GasParticle)*(np));
	j = 0;
	for(i=0;i<np;i++){
//		if( fabs(bp3d[i].id/(Nxp*Nyp) - (Nzp/2)) < 5 || bp3d[i].id >= Nxp*Nyp*Nzp)
		if( fabs(bp3d[i].z - zcut) < 0.2 || bp3d[i].id >= Nxp*Nyp*Nzp)
		{
			bp[j] = bp3d[i];
			j++;
		}
	}
	free(bp3d);
	np = j;

	int maxnp = LinkedCell(bp,np);
	ptype meanvol = Lx*Ly/Nxp/Nyp;

	int ix,iy,iz,iys,iyf,ixs,ixf,izs,izf;
	ixs = floor((xmin)*Nx/Lx);
	ixf = ceil((xmax)*Nx/Lx);
	iys = floor((ymin)*Ny/Ly);
	iyf = ceil((ymax)*Ny/Ly);

	xmin -= 0.01;
	ymin -= 0.01;
	xmax += 0.01;
	ymax += 0.01;

//	ptype dfact = Lx*Ly*Lz/Nxp/Nyp/Nzp/1270.; 
	ptype dfact = 1.;
	printf("Density factor is %g\n", dfact);

    Voro3D_Vertex *vorovertex = (Voro3D_Vertex*)malloc(sizeof(Voro3D_Vertex)*4024);
    Voro3D_GasParticle *p = (Voro3D_GasParticle*)malloc(sizeof(Voro3D_GasParticle)*maxnp);
    Voro3D_point *neighbors = (Voro3D_point*)malloc(sizeof(Voro3D_point)*maxnp*25);
    Voro3D_point *neighwork = (Voro3D_point*)malloc(sizeof(Voro3D_point)*maxnp*25);
    int *intwork = (int*)malloc(sizeof(int)*4024);
	for(iy=iys;iy<iyf;iy++){
        int mp=Nxp*Nyp;
        ptype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=ixs;ix<ixf;ix++){
            int mp;
            Voro2D_FindCellBP(p,ix,iy,&mp,bp);
            int nneigh;
            Voro2D_FindNeighbor(neighbors,ix,iy,&nneigh,bp);
            int i;
            for(i=0;i<mp;i++){
                Voro3D_point center;
                center.x = p[i].x;
                center.y = p[i].y;
                center.z = p[i].z;
                center.id = p[i].id;
				if(center.x<xmin || center.x >xmax || center.y<ymin || center.y >ymax) continue;


                int id = p[i].id;
				int ishrink = 0;
                int ip = Voro3D_FindVC(&center,neighbors, neighwork,nneigh, 
						vorovertex,mp,boxsize, ishrink,intwork);

				Voro3D_point neighmirror;
				neighmirror.x = 0;
				neighmirror.y = 0;
				neighmirror.z = 2*(zcut-center.z);

				if(fabs(neighmirror.z)<1.e-5) neighmirror.z= 1.e-5;
				int jp=ip;
				j = 0;
				while(j<ip){
					ptype offset;
					if(vorovertex[j].status==Active && 
							Voro3D_CutOrStay(&neighmirror,vorovertex+j,&offset) == Outside){
						Voro3D_Vertex2 start; 
						start = FindNearestBoundaryVertex(&neighmirror, vorovertex+j); 
						int icount = InactivateOutsideVertex(&neighmirror, vorovertex,ip); 
						jp = CreateNewVertices(&neighmirror,0,&start,vorovertex, ip);
						break;
					}
					else {
						j++;
					}
				}
				if(jp>ip) {
					int icolor = setpolygoncolor(nlog, p[i].den*dfact);
					fprintf(wp,"set x = %g\n", p[i].x+vorovertex[ip].x);
					for(j=ip+1;j<jp;j++){
						fprintf(wp,"set x = x concat %g\n", p[i].x+vorovertex[j].x);
					}
					fprintf(wp,"set x = x concat %g\n", p[i].x+vorovertex[ip].x);
					fprintf(wp,"set y = %g\n", p[i].y+vorovertex[ip].y);
					for(j=ip+1;j<jp;j++){
						fprintf(wp,"set y = y concat %g\n", p[i].y+vorovertex[j].y);
					}
					fprintf(wp,"set y = y concat %g\n", p[i].y+vorovertex[ip].y);
	
					fprintf(wp,"ctype %s\n", mycolor[icolor]);
					fprintf(wp,"shade 0 x y\n");
					fprintf(wp,"ctype default  connect x y\n");
					fprintf(wp,"delete x\n");
					fprintf(wp,"delete y\n");
					fprintf(wp,"ptype 0 0\n");
					fprintf(wp,"ctype white relocate %g %g dot\n", p[i].x, p[i].y);
//					else if(bp[p[i].id].iregion==1) fprintf(wp,"ctype yellow  relocate %g %g dot\n", p[i].x, p[i].y);
//					fprintf(wp,"expand 0.1 lweight 1 draw_arrow %g %g %g %g expand 1\n", p[i].x, p[i].y, p[i].x+p[i].vx*dt, p[i].y+p[i].vy*dt);
					fprintf(wp,"delete history\n");
				}

			}
		}
	}
	/*
	free(vorovertex);
	free(p);free(neighbors);free(neighwork);
	*/
	fprintf(wp,"ctype default box \n");
	ptype alpha,beta;
	alpha = alphavis;
	beta = betavis;
	fprintf(wp,"ctype white limits 0 1 0 1 relocate 0.8 0.1 label \\alpha=%g\n",alpha);
	fprintf(wp,"relocate 0.8 0.05 label \\beta=%g\n", beta);
	fprintf(wp,"relocate 0.6 0.05 label t=%6.4f\n", time);

	mkcolorbar(wp,1.e-2, 4., nlog);
	smfinalize(wp, time);
	return 0;
}

