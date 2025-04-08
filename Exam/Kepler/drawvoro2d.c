#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "voro.h"
#include "kp.h"


char mycolor[256][20];

int setcolor(char *, int);

CellType *cells=NULL;

float cellsize;
postype Lx,Ly;

int Nx,Ny;

float boxsize = 0.1;


void MkLinkedListSerial(treevorork4particletype *bp, int np){
    int i,j,k;
    if(cells == NULL) cells = (CellType*)malloc(sizeof(CellType)*Nx*Ny);


    cellsize = Lx/Nx;

    for(i=0;i<Nx*Ny;i++) {
        cells[i].link = NULL;
        cells[i].nmem = 0;
    }

    for(i=0;i<np;i++){
        int ix,iy;
        ix = (bp[i].x/cellsize);
        iy = (bp[i].y/cellsize);
        size_t index = ix+Nx*iy;
        struct linkedlisttype *tmp = cells[index].link;
        cells[index].link = (struct linkedlisttype*)(bp+i);
        cells[index].nmem ++;
        bp[i].next = tmp;
    }
}

Voro2D_point *Voro2D_FindNeighborSerial(int ix, int iy, int *nneigh, treevorork4particletype *bp){
    int i,j,k;
    int np,mp;
    Voro2D_point *neigh;

    np = 0;
    for(j=iy-2;j<=iy+2;j++){
        int jj = (j+Ny)%Ny;
        for(i=ix-2;i<=ix+2;i++){
            int ii = (i+Nx)%Nx;
            np += cells[ii+Nx*jj].nmem;
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
            struct linkedlisttype *tmp = cells[ii+Nx*jj].link;
            while(tmp){
				treevorork4particletype *tt = (treevorork4particletype*)tmp;
                neigh[np].x = tt->x + iflag *Lx;
                neigh[np].y = tt->y + jflag *Ly;
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
treevorork4particletype *Voro2D_FindCellBPSerial(int ix, int iy, int *mp, treevorork4particletype *bp){
    int i,j,k;
    int np;
    treevorork4particletype *res;

    int iix,iiy;
    iix = (ix+Nx)%Nx;
    iiy = (iy+Ny)%Ny;

    np = cells[iix+Nx*iiy].nmem;
    res = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*np);
    np = 0;
    int iflag,jflag;
    iflag = jflag = 0;
    if(iix < 0) iflag = -1;
    else if(iix >= Nx) iflag = 1;
    if(iiy < 0) jflag = -1;
    else if(iiy >= Ny) jflag = 1;
    size_t ipixel = iix+Nx*iiy;
    struct linkedlisttype *tmp = cells[ipixel].link;
    while(tmp){
		treevorork4particletype *tt = (treevorork4particletype*)tmp;
        res[np] = *tt;
        res[np].bp = tt;
        np++;
        tmp = tmp->next;
    }
    *mp = np;
    return res;
}




static int xlocmin,ylocmin, xlocmax,ylocmax;

int sminitialize(FILE *wp, float xmin, float xmax, float ymin, float ymax){
	float bigxtick,smxtick,bigytick,smytick;
	bigxtick = rint(10*((xmax-xmin)/5))/10.;
	smxtick = rint(50*((xmax-xmin)/20))/50;
	bigytick = rint(10*((ymax-ymin)/5))/10.;
	smytick = rint(50*((ymax-ymin)/20))/50.;
	int i;
    fprintf(wp,"erase\n");
	for(i=0;i<256;i++){
		sprintf(mycolor[i],"color%d",i);
    	fprintf(wp,"ctype = CTYPE() concat %d\n", setcolor("kp.sao",i));
	    fprintf(wp,"ctype = CTYPE(STRING) concat '%s'\n", mycolor[i]);
	}
	float xsize = (xmax-xmin);
	float ysize = (ymax-ymin);
	float maxsize = xsize; 
	if(ysize > xsize) maxsize = ysize;
	float scale =(27000.-5000.) /maxsize;
	if(xsize > ysize){
		xlocmin = 5000;
		xlocmax = 27000;
		ylocmin = 5000;
		ylocmax = scale*(ymax-ymin) + 5000.;
	}
	else {
		ylocmin = 5000;
		ylocmax = 27000;
		xlocmin = 5000;
		xlocmax = scale*(xmax-xmin) + 5000.;
	}


    fprintf(wp,"ctype default\n");
    fprintf(wp,"ltype 0\n");
    fprintf(wp,"location %d %d %d %d\n", xlocmin, xlocmax, ylocmin, ylocmax);
//    fprintf(wp,"location 5000 27000 5000 27000\n");
    fprintf(wp,"limits %g %g %g %g\n", xmin,xmax,ymin,ymax);
    fprintf(wp,"ticksize %g %g %g %g\n", smxtick, bigxtick, smytick, bigytick);
//    fprintf(wp,"box \n");
	return 0;
}

int smfinalize(FILE *wp, postype time){
	int i,j,k;
//    fprintf(wp,"ctype default \nbox \n");
	for(i=0;i<256;i++){
		fprintf(wp,"del_ctype %s\n", mycolor[i]);
	}
	return 0;
}

int getcolor(int );
int setpolygoncolor(int nlog, float val){
	int res;
	nlog = exp(nlog);
//	val = val-1;
    val = val/2.;
	if(nlog>0){
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
void mkcolorbar(FILE *wp, float valmin, float valmax,int nlog){
    float lgvalmin = (valmin);
    float lgvalmax = (valmax);
    float valstep = ((lgvalmax)-(lgvalmin))/128.;
//    fprintf(wp,"location 27050 28550 5000 27000\n");
    fprintf(wp,"location %d %d 5000 27000\n",xlocmax+50, xlocmax+1550);
    fprintf(wp,"limits 0 1 %g %g\n", lgvalmin, lgvalmax);
    fprintf(wp,"ticksize 10 100 0.1 0.5\n");
    int i;
    fprintf(wp,"set x = 0, 1, 1 \n");
    fprintf(wp,"set x = x concat reverse(x) \n");
    for(i=0;i<128;i++){
        float val = valstep*(i)+lgvalmin;
//        float rval = pow(10., val);
        float rval = val;
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
	int nlog = -1;
	postype time,dt;
	postype xmin,ymin,xmax,ymax;
	FILE *fp, *wp;
	if(argc ==2) {
		char infile[190];
		sprintf(infile,"kpout.%.6d.dat", atoi(argv[1]));
		fp = fopen(infile,"r");
		wp = fopen("a.sm","w");
		xmin = 2.;
		xmax = 4;
		ymin = 2;
		ymax = 4;
	}
	else if(argc == 7){
		fp = fopen(argv[1],"r");
		wp = fopen(argv[2],"w");
		xmin = atof(argv[3]);
		xmax = atof(argv[4]);
		ymin = atof(argv[5]);
		ymax = atof(argv[6]);
	}
	sminitialize(wp, xmin,xmax,ymin,ymax);
	fread(&np,sizeof(int),1,fp);
	treevorork4particletype *bp = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*(np));
    fread(bp,sizeof(treevorork4particletype), np,fp);
    fread(&time,sizeof(postype), 1,fp);
    fread(&Nx,sizeof(int), 1,fp);
    fread(&Ny,sizeof(int), 1,fp);
    if(fread(&Lx,sizeof(postype), 1,fp) ==0){
		Lx = 0.5;
		Ly = 1;
	}
	else {
    	fread(&Ly,sizeof(postype), 1,fp);
	}
	float alpha,beta;
	fread(&alpha, sizeof(float),1,fp);
	fread(&beta, sizeof(float),1,fp);
    fclose(fp);
	printf("Lx/y= %g %g and Nx/y = %d %d\n", Lx,Ly, Nx, Ny);

	int Nxp,Nyp;
	Nxp = Nx;
	Nyp = Ny;

	cellsize = Ly/Ny*4;

	MkLinkedListSerial(bp,np);



	postype meanvol = Lx*Ly/Nxp/Nyp;

	int ix,iy,iys,iyf,ixs,ixf;
	ixs = floor((xmin)*Nx/Lx);
	ixf = ceil((xmax)*Nx/Lx);
	iys = floor((ymin)*Ny/Ly);
	iyf = ceil((ymax)*Ny/Ly);

	xmin -= 0.001;
	ymin -= 0.001;
	xmax += 0.001;
	ymax += 0.001;
	for(iy=iys;iy<iyf;iy++){
        int mp=Nxp*Nyp;
        Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp);
        postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=ixs;ix<ixf;ix++){
            int np;
            treevorork4particletype *p = Voro2D_FindCellBPSerial(ix,iy,&np,bp);
            int nneigh;
            Voro2D_point *neighbors = Voro2D_FindNeighborSerial(ix,iy,&nneigh, bp);
            Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh);
            int i;
            for(i=0;i<np;i++){
				double dx = p[i].x-3.L;
				double dy = p[i].y-3.L;
				double rdist = sqrt(dx*dx+dy*dy);
				int icolor;
                Voro2D_point center; 
				center.x = p[i].x; 
				center.y = p[i].y; 
				center.indx = PINDX(p+i);
				if(center.x<xmin || center.x >xmax || center.y<ymin || center.y >ymax) continue; 
				int id = PINDX(p+i);
				int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize); 
				treevorork4particletype *ibp = p[i].bp;
				ibp->volume = Area2DPolygon(vorocorner, mp); 
				bp[id].den = bp[id].mass/bp[id].volume; 
//				ibp->den = ibp->mass/meanvol;
				/*
				if(rdist >2.|| rdist <0.5){
					icolor = setpolygoncolor(nlog, 0.);
				}
				else {
				*/
					icolor = setpolygoncolor(nlog, ibp->den);
					/*
				}
				*/
				fprintf(wp,"set x = %g\n", p[i].x+vorocorner->x);
				Voro2D_Corner *tmp = vorocorner->upperlink;
				do{
					fprintf(wp,"set x = x concat %g\n", p[i].x+tmp->x);
					tmp = tmp->upperlink;
				} while(tmp != vorocorner);
				fprintf(wp,"set x = x concat %g\n", p[i].x+tmp->x);
 				fprintf(wp,"set y = %g\n", p[i].y+vorocorner->y);
				tmp = vorocorner->upperlink;
				do{
					fprintf(wp,"set y = y concat %g\n", p[i].y+tmp->y);
					tmp = tmp->upperlink;
				} while(tmp != vorocorner);
				fprintf(wp,"set y = y concat %g\n", p[i].y+tmp->y);
				fprintf(wp,"ctype %s\n", mycolor[icolor]);
				fprintf(wp,"shade 0 x y\n");
				fprintf(wp,"ctype default  connect x y\n");
				fprintf(wp,"delete x\n");
				fprintf(wp,"delete y\n");
				fprintf(wp,"ptype 0 0\n");
//				if(bp[p[i].id].iregion==0 || bp[p[i].id].iregion==2) fprintf(wp,"ctype white relocate %g %g dot\n", p[i].x, p[i].y);
//				else if(bp[p[i].id].iregion==1) fprintf(wp,"ctype yellow  relocate %g %g dot\n", p[i].x, p[i].y);
//
//				fprintf(wp,"expand 0.1 lweight 1 draw_arrow %g %g %g %g expand 1\n", p[i].x, p[i].y, p[i].x+p[i].vx*dt, p[i].y+p[i].vy*dt);
				fprintf(wp,"delete history\n");

			}
			free(p);free(neighbors);free(neighwork);
		}
		free(vorocorner);
	}
	fprintf(wp,"ctype default box \n");
	{
		fprintf(wp,"ctype white limits 0 1 0 1 relocate 0.8 0.1 label \\alpha=%g\n",alpha);
		fprintf(wp,"relocate 0.8 0.05 label \\beta=%g\n", beta);
		fprintf(wp,"relocate 0.6 0.05 label t=%5.3f\n", time);
	}
	mkcolorbar(wp,0., 2., nlog);
	smfinalize(wp, time);
	fclose(wp);
	return 0;
}

