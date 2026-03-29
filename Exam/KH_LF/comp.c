#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
#include "voro.h"
#include "kh.h"

static int xlocmin,ylocmin, xlocmax,ylocmax;



char mycolor[256][20];

int setcolor(char *, int);

CellType *cells=NULL;

float cellsize;
postype Lx,Ly;

int Nx,Ny;

float boxsize = 0.1;

int vorosort(const void *a, const void *b){
	treevoroparticletype *aa, *bb;
	aa = (treevoroparticletype*)a;
	bb = (treevoroparticletype*)b;
	if(aa->u4if.indx < bb->u4if.indx) return -1;
	else if(aa->u4if.indx > bb->u4if.indx) return 1;
	else return 0;
}
int vorork4sort(const void *a, const void *b){
	treevorork4particletype *aa, *bb;
	aa = (treevorork4particletype*)a;
	bb = (treevorork4particletype*)b;
	if(aa->u4if.indx < bb->u4if.indx) return -1;
	else if(aa->u4if.indx > bb->u4if.indx) return 1;
	else return 0;
}


int main(int argc,char **argv){
	int i,j,k;
	int np;
	int nlog = -1;
	postype time,dt;
	postype xmin,ymin,xmax,ymax;
	float dtold;
	FILE *fp, *fp2;

	fp = fopen(argv[1], "r");
	fread(&np,sizeof(int),1,fp);
	treevoroparticletype *ap = (treevoroparticletype*)malloc(sizeof(treevoroparticletype)*(np));
    fread(ap,sizeof(treevoroparticletype), np,fp);
	fclose(fp);
	fp = fopen(argv[2], "r");
	fread(&np,sizeof(int),1,fp);
	treevorork4particletype *bp = (treevorork4particletype*)malloc(sizeof(treevorork4particletype)*(np));
    fread(bp,sizeof(treevorork4particletype), np,fp);
	fclose(fp);
	qsort(ap, np, sizeof(treevoroparticletype), vorosort);
	qsort(bp, np, sizeof(treevorork4particletype), vorork4sort);
	for(i=0;i<np;i++){
		postype ax = fabs((ap[i].vx-bp[i].vx)/bp[i].vx);
		postype ay = fabs((ap[i].vy-bp[i].vy)/bp[i].vy);
		if(ax > 0.004 || ay >0.004){
			printf("p%ld has difference: %g %g : %g %g\n", PINDX(ap+i), ap[i].vx, bp[i].vx, ap[i].vy, bp[i].vy);
		}
	}
}
