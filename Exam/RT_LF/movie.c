#include "voro.h"
#include "kh.h"

/*
CellType *mycell=NULL;

postype cellsize = Lx/Nx;


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
    fprintf(wp,"device postencap bb.eps\n");
    fprintf(wp,"erase\n");
	for(i=0;i<256;i++){
		sprintf(mycolor[i],"color%d",i);
    	fprintf(wp,"ctype = CTYPE() concat %d\n", setcolor("kh.sao",i));
	    fprintf(wp,"ctype = CTYPE(STRING) concat '%s'\n", mycolor[i]);
	}
    fprintf(wp,"ctype default\n");
    fprintf(wp,"ltype 0\n");
    fprintf(wp,"location 5000 30000 5000 30000\n");
    fprintf(wp,"limits %g %g %g %g\n", xmin,xmax,ymin,ymax);
    fprintf(wp,"ticksize %g %g %g %g\n", smxtick, bigxtick, smytick, bigytick);
//    fprintf(wp,"box \n");
}

int smfinalize(FILE *wp, postype time){
	int i,j,k;
    fprintf(wp,"ctype default \nbox \n");
	for(i=0;i<256;i++){
		fprintf(wp,"del_ctype %s\n", mycolor[i]);
	}
	postype alpha,beta;
	alpha = alphavis;
	beta = betavis;
	fprintf(wp,"ctype white limits 0 1 0 1 relocate 0.8 0.1 label \\alpha=%g\n",alpha);
	fprintf(wp,"relocate 0.8 0.05 label \\beta=%g\n", beta);
	fprintf(wp,"relocate 0.6 0.05 label t=%5.3f\n", time);
    fprintf(wp,"hardcopy \n");
	fclose(wp);
}
int getcolor(int );
int setpolygoncolor(int nlog, float val){
	int res;
	nlog = exp(nlog);
	val = val - 1;
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

int main(int argc,char **argv){
	int i,j,k;
	int np= Nxp*Nyp;
	int nlog = -1;
	postype time,dt;
	FILE *fpfile = fopen(argv[1],"r");
	int istep = 0;

	char infile[190];
	int iflag = fscanf(fpfile,"%s", infile);
	Voro2D_GasParticle *bp = (Voro2D_GasParticle*)malloc(sizeof(Voro2D_GasParticle)*(np)); 
	while(iflag>0) {
		FILE *fp = fopen(infile,"r");
		FILE *wp = fopen(argv[2],"w");
		postype xmin,ymin,xmax,ymax;
		/*
		xmin = atof(argv[3]);
		xmax = atof(argv[4]);
		ymin = atof(argv[5]);
		ymax = atof(argv[6]);
		*/
		sminitialize(wp, xmin,xmax,ymin,ymax);
		fread(&np,sizeof(int),1,fp);
		fread(bp,sizeof(Voro2D_GasParticle), np,fp); 
		fread(&time,sizeof(postype), 1,fp); 
//		fread(&dt,sizeof(ptype), 1,fp); 
		fclose(fp);
	
		MkLinkedList(bp,np);
		postype meanvol = Lx*Ly/Nxp/Nyp;

		xmin = 0;
		ymin = 0;
		xmax = 1;
		ymax = 1;

		int ix,iy,iys,iyf,ixs,ixf;
		ixs = floor((xmin)*Nx/Lx);
		ixf = ceil((xmax)*Nx/Lx);
		iys = floor((ymin)*Ny/Ly);
		iyf = ceil((ymax)*Ny/Ly);
	

		for(iy=iys;iy<iyf;iy++){ 
			int mp=Nxp*Nyp; 
			Voro2D_Corner *vorocorner = (Voro2D_Corner*)malloc(sizeof(Voro2D_Corner)*mp); 
			postype dlx,dly,dl,dvx,dvy,dv,ax,ay,a; 
			int ix; 
			for(ix=ixs;ix<ixf;ix++){
	            int np; 
				Voro2D_GasParticle *p = Voro2D_FindCellBP(ix,iy,&np,bp); 
				int nneigh; 
				Voro2D_point *neighbors = Voro2D_FindNeighbor(ix,iy,&nneigh, bp); 
				Voro2D_point *neighwork = (Voro2D_point*)malloc(sizeof(Voro2D_point)*nneigh); 
				int i; 
				for(i=0;i<np;i++){ 
					Voro2D_point center; 
					center.x = p[i].x; 
					center.y = p[i].y; 
					center.id = p[i].id; 
					if(center.x<xmin || center.x >xmax || center.y<ymin || center.y >ymax) continue;

	                int id = p[i].id; 
					int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize); 
					bp[id].volume = Area2DPolygon(vorocorner, mp); 
					bp[id].den = bp[id].mass/meanvol; 
					int icolor = setpolygoncolor(nlog, bp[id].den); 
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
//				fprintf(wp,"lweight 1 ctype green relocate %g %g expand 0.3 label %5.3f\n", p[i].x,p[i].y,bp[id].pressure);
//				fprintf(wp,"lweight 1 ctype yellow relocate %g %g expand 0.3 label %5.3f\n", p[i].x,p[i].y-0.0005,bp[id].volume*512*512);
					if(bp[p[i].id].u4if.Flag[ENDIAN_OFFSET]==0 || bp[p[i].id].u4if.Flag[ENDIAN_OFFSET]==2) fprintf(wp,"ctype white relocate %g %g dot\n", p[i].x, p[i].y);
					else if(bp[p[i].id].u4if.Flag[ENDIAN_OFFSET]==1) fprintf(wp,"ctype yellow  relocate %g %g dot\n", p[i].x, p[i].y);
//				fprintf(wp,"expand 0.1 lweight 1 draw_arrow %g %g %g %g expand 1\n", p[i].x, p[i].y, p[i].x+p[i].vx*dt, p[i].y+p[i].vy*dt);
					fprintf(wp,"delete history\n");

				}
				free(p);free(neighbors);free(neighwork);
			}
			free(vorocorner);
		}
		smfinalize(wp,time);
		char outfile[190];
		sprintf(outfile,"sm < %s ; convert -background white -resize 1863x2636 -units pixelsperinch -density 224.993 -density 600 -flatten bb.eps movie.%.6d.gif", argv[2],istep);
		system(outfile);
		istep ++;
		iflag = fscanf(fpfile,"%s",infile);
	} 
	return 0;
}

