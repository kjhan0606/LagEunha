#include "voro.h"
#include "kh.h"

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
	bigxtick = rint(10*((xmax-xmin)/5))/10.;
	smxtick = rint(50*((xmax-xmin)/20))/50;
	bigytick = rint(10*((ymax-ymin)/5))/10.;
	smytick = rint(50*((ymax-ymin)/20))/50.;
	int i;
    fprintf(wp,"erase\n");
	for(i=0;i<256;i++){
		sprintf(mycolor[i],"color%d",i);
    	fprintf(wp,"ctype = CTYPE() concat %d\n", setcolor("kh.sao",i));
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
}

int getcolor(int );
int setpolygoncolor(int nlog, float val){
	int res;
	nlog = exp(nlog);
	val = val-1;
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
//    fprintf(wp,"location 28000 29500 7000 25000\n");
    fprintf(wp,"location 27050 28550 5000 27000\n");
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
	ptype time,dt;
	ptype xmin,ymin,xmax,ymax;
	FILE *fp, *wp;
	if(argc ==2) {
		char infile[190];
		sprintf(infile,"khout.%.6d.dat", atoi(argv[1]));
		fp = fopen(infile,"r");
		wp = fopen("a.sm","w");
		xmin = 0.3;
		xmax = 0.8;
		ymin = 0.499;
		ymax = 0.99;
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
	Voro2D_GasParticle *bp = (Voro2D_GasParticle*)malloc(sizeof(Voro2D_GasParticle)*(np));
    fread(bp,sizeof(Voro2D_GasParticle), np,fp);
    fread(&time,sizeof(ptype), 1,fp);
    fread(&dt,sizeof(ptype), 1,fp);
    fclose(fp);

	MkLinkedList(bp,np);
	ptype meanvol = Lx*Ly/Nxp/Nyp;

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
        ptype dlx,dly,dl,dvx,dvy,dv,ax,ay,a;
        int ix;
        for(ix=ixs;ix<ixf;ix++){
            int np;
            Voro2D_GasParticle *p = Voro2D_FindCellBP(ix,iy,&np,bp);
            int nneigh;
            Voro2D_point *neighbors = Voro2D_FindNeighbor(ix,iy,&nneigh, bp);
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
				center.id = p[i].id; 
				if(center.x<xmin || center.x >xmax || center.y<ymin || center.y >ymax) continue; 
				int id = p[i].id; 
				int ip = Voro2D_FindVC(&center,neighbors, neighwork,nneigh, vorocorner,mp,boxsize); 
				bp[id].volume = Area2DPolygon(vorocorner, mp); 
//				bp[id].den = bp[id].mass/bp[id].volume; 
				bp[id].den = bp[id].mass/meanvol;
				/*
				if(rdist >2.|| rdist <0.5){
					icolor = setpolygoncolor(nlog, 0.);
				}
				else {
				*/
					icolor = setpolygoncolor(nlog, bp[id].den);
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
		ptype alpha,beta;
		alpha = alphavis;
		beta = betavis;
		fprintf(wp,"ctype white limits 0 1 0 1 relocate 0.8 0.1 label \\alpha=%g\n",alpha);
		fprintf(wp,"relocate 0.8 0.05 label \\beta=%g\n", beta);
		fprintf(wp,"relocate 0.6 0.05 label t=%5.3f\n", time);
	}
	mkcolorbar(wp,1., 2., nlog);
	smfinalize(wp, time);
	fclose(wp);
	return 0;
}

