#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>

#include "color.h"

//#include "voro.h"
//#include "kh.h"


//#define MAX(a,b) ( (a) > (b) ? (a): (b))


typedef struct Img{
	int reverse;
	int nx, ny;
	float nlog, ccdmin0, ccdmax0;
	float *ccd;
	unsigned char *r,*g,*b, *work;
	char saofile[190];
} Img;

void f2pgm(Img *);
void pgm2col(Img *); 
void pickcol(unsigned char, unsigned char *, unsigned char *, unsigned char*);
void readpgm(unsigned char *, int , int , unsigned char *, unsigned char *, unsigned char *);
void readsao(int , FILE *);

void initializeImg(Img *img, int nx, int ny, char *saofile){

	img->nx = nx;
	img->ny = ny;
	img->ccd = (float*)malloc(sizeof(float)*nx*ny);
	img->r = (unsigned char*)malloc(sizeof(unsigned char)*nx*ny);
	img->g = (unsigned char*)malloc(sizeof(unsigned char)*nx*ny);
	img->b = (unsigned char*)malloc(sizeof(unsigned char)*nx*ny);
	img->work = (unsigned char*)malloc(sizeof(unsigned char)*nx*ny);
	img->ccdmin0 = img->ccdmax0 = 9999.;
	img->nlog = 13;
	img->reverse = 0;
	sprintf(&(img->saofile[0]),"%s", saofile);
}
void ImgHeader(int *nx, int *ny, char *infile){
	int chip;
	FILE *fp= fopen(infile,"r");
	fread(nx, sizeof(int), 1, fp);
	fread(ny, sizeof(int), 1, fp);
	fclose(fp);
}
void HR5pixel(Img *img, char *infile){
	int i,j,chip;
	int nx, ny;
	FILE *fp = fopen(infile, "r");
	/*
	fread(&chip, sizeof(int), 1, fp);
	fread(&nx, sizeof(int), 1, fp);
	fread(&chip, sizeof(int), 1, fp);
	fread(&chip, sizeof(int), 1, fp);
	fread(&ny, sizeof(int), 1, fp);
	fread(&chip, sizeof(int), 1, fp);
	fread(&chip, sizeof(int), 1, fp);
	*/
	float *org = (float*)malloc(sizeof(float)*img->nx*img->ny);
	fread(&nx, sizeof(int), 1, fp);
	fread(&ny, sizeof(int), 1, fp);
//	fread(img->ccd, sizeof(float), img->nx*img->ny, fp);
	fread(org, sizeof(float), img->nx*img->ny, fp);
	for(j=0;j<img->ny;j++){
		int jj = img->ny-j-1;
		for(i=0;i<img->nx;i++){
			img->ccd[i+img->nx*jj] = org[i+img->nx*j];
		}
	}
	free(org);
	/*
	fread(&chip, sizeof(int), 1, fp);
	*/
	fclose(fp);
}

void colorizeit(
//		float *map, int nx, int ny, char *saofile, char *outfile
		int iarg, ...
		){
	va_list ap;
	va_start(ap,iarg);
	float *map = va_arg(ap, float*);
	int nx = va_arg(ap, int);
	int ny = va_arg(ap, int);
	char *saofile = va_arg(ap, char *);
	char *outfile = va_arg(ap, char *);
	Img gas, star, dm, temp, composit, smetal;

	unsigned int *r,*g, *b;
	int i,j,k;
	char infile[190];
	FILE *fp;
	int nsize;


//	ImgHeader(&nx,&ny, outfile);
	nsize = nx*ny;

	{
		initializeImg(&gas, nx, ny, saofile);
		gas.ccdmax0 = 2;
		gas.ccdmin0 = 1.;
		gas.nlog = -1;
		gas.ccd = map;
		gas.nx = nx;
		gas.ny = ny;

		{/* for the input range */
			if(iarg ==7){
				gas.ccdmin0 = va_arg(ap, double );
				gas.ccdmax0 = va_arg(ap, double );
			}
			else {
				gas.ccdmin0 = 1;
				gas.ccdmax0 = 2;
			}
			va_end(ap);
		}
//		HR5pixel(&gas, argv[2]);
		f2pgm(&gas);
		pgm2col(&gas);
	
		r = (unsigned int*)malloc(sizeof(unsigned int)*nsize);
		b = (unsigned int*)malloc(sizeof(unsigned int)*nsize);
		g = (unsigned int*)malloc(sizeof(unsigned int)*nsize);
		for(i=0;i<nsize;i++){
			{
				if(gas.r[i]> 255) gas.r[i] = 255;
				if(gas.g[i]> 255) gas.g[i] = 255;
				if(gas.b[i]> 255) gas.b[i] = 255;
			}
			r[i] = gas.r[i];
			g[i] = gas.g[i];
			b[i] = gas.b[i];
		}
	}

	FILE *out = fopen(outfile,"w");
	fprintf(out,"P6\n%d %d\n255\n",nx,ny);
	for(i=0;i<nsize;i++){
		fwrite(&(r[i]),1,1,out);
		fwrite(&(g[i]),1,1,out);
		fwrite(&(b[i]),1,1,out);
	}
	fclose(out);
}


void f2pgm(Img *img) {
	int nx, ny;
	float nlog;
	float ccdmin0 , ccdmax0, *ccd;
	unsigned char *lccd;
	int i, j, ix, iy;
	int tmpccd;
	float ccdmin, ccdmax;

	nx = img->nx;
	ny = img->ny;
	nlog = img->nlog;
	ccdmin0 = img->ccdmin0;
	ccdmax0 = img->ccdmax0;
	ccd = img->ccd;
	lccd = img->work;

	nlog = exp(nlog);
	nlog = -1;
	int iab=0;
	for(i=0;i<nx*ny;i++){
		if(isnan(ccd[i]) || isinf(ccd[i])) iab ++;
	}
	if(iab>0) fprintf(stderr,"Error in the pixel value : %d\n",iab);

	ccdmin = ccd[0];
	ccdmax = ccd[0];
	for(i=0; i<nx*ny; i++) {
		if( ccdmax < ccd[i] ) ccdmax = ccd[i];
		if( ccdmin > ccd[i] ) ccdmin = ccd[i];
	}
	fprintf(stderr,"ccdmax %g , ccdmin %g nx/ny %d/%d\n",ccdmax,ccdmin,nx,ny);
	if( ccdmax0 == 9999. ) {
		ccdmax0 = ccdmax;
	}
	else {
		ccdmax = ccdmax0;
	}
	if( ccdmin0 == 9999. ) {
		ccdmin0 = ccdmin;
	}
	else {
		ccdmin = ccdmin0;
	}

/* renormalize to a number between 0 and 1 */
	for(i=0; i<nx*ny; i++)
		ccd[i] = (ccd[i] - ccdmin)/(ccdmax - ccdmin);

/* Do log transformation */
	if( nlog > 0 ) {
		ccdmax = log(1.0+nlog);
		for(i=0; i<nx*ny; i++) {
			ccd[i] = log((1.0 + nlog*ccd[i]));
			tmpccd = (int) (ccd[i]/ccdmax*255);
			if( tmpccd > 255 ) tmpccd = 255;
			else if( tmpccd < 0 ) tmpccd = 0;
			lccd[i] = (unsigned char) tmpccd;
		}
	}
	else {
		for(i=0; i<nx*ny; i++) {
			tmpccd = (int) (ccd[i]*255);
			if( tmpccd > 255 ) tmpccd = 255;
			else if( tmpccd < 0 ) tmpccd = 0;
			lccd[i] = (unsigned char) tmpccd;
		}
	}
//	return 1;
}

float rlevel[500], rbright[500];
float glevel[500], gbright[500];
float blevel[500], bbright[500];
float rslope[500], rx[500];
float gslope[500], gx[500];
float bslope[500], bx[500];

unsigned char rtab[256], gtab[256], btab[256];

int saoflag = 0;
int getcolor(int gray){
	int res = rtab[gray] + 256*(gtab[gray]+256*btab[gray]);
	return res;
}
int setcolor(char *saofile, int gray){
	int reverse = 0;
	FILE *colfile = fopen(saofile,"r");

	if(saoflag ==0) readsao(reverse,colfile);
	int res = rtab[gray] + 256*(gtab[gray]+256*btab[gray]);
	fclose(colfile);
	return res;
}

void pgm2col(Img *img){ 
	char *saofile;
	int reverse;
	unsigned char *ccd, *r,*g,*b;
	int nx, ny;

	nx = img->nx;
	ny = img->ny;
	saofile = (img->saofile);
	reverse = img->reverse;
	ccd = (img->work);
	r = (img->r);
	g = (img->g);
	b = (img->b);
	int i, j;
	FILE *pgmfile, *colfile;

	colfile = fopen(saofile,"r");

	readsao(reverse,colfile);
	readpgm(ccd, nx,ny, r,g,b);
}

void readpgm(unsigned char *ccd, int nx, int ny, unsigned char *rr, unsigned char *gg, unsigned char *bb) {
	int i, size;
	char line[80];

	size = nx*ny;
	unsigned char *r, *g, *b;
	r = rr;
	g = gg;
	b = bb;
	for(i=0; i<size; i++) {
		unsigned char gray;
		
		gray = ccd[i]; *r = rtab[gray]; *g = gtab[gray]; *b = btab[gray];
		r++;
		g++; 
		b++;
	}
}

int nr, ng, nb;

void readsao(int reverse, FILE *saofile) {
	int i;
	int ir, ig, ib, len;
	int n0;
	unsigned char r, g, b;
	char line[80], *p, n[15], *p0;

#ifdef SIMPLE
    i = 0;
	while( fscanf(saofile,"%f %f %f\n",&rbright[i],&gbright[i],&bbright[i]) != EOF) {
		i++;
	}
	n0 = i;
	for(i=0; i<n0; i++) {
		rbright[i] *= 255; gbright[i] *= 255; bbright[i] *= 255;
		rlevel[i] = (255.*i)/n0; glevel[i] = (255.*i)/n0; blevel[i] = (255.*i)/n0;
	}
	fclose(saofile);
	nr = ng = nb = n0;
		
#else
	fscanf(saofile,"%d\n",&nr);
	for(i=0; i<nr; i++) {
		fscanf(saofile,"%f %f\n",rlevel+i,rbright+i);
		rlevel[i] *= 255; rbright[i] *= 255;
	}
	fscanf(saofile,"%d\n",&ng);
	for(i=0; i<ng; i++) {
		fscanf(saofile,"%f %f\n",glevel+i,gbright+i);
		glevel[i] *= 255; gbright[i] *= 255;
	}
	fscanf(saofile,"%d\n",&nb);
	for(i=0; i<nb; i++) {
		fscanf(saofile,"%f %f\n",blevel+i,bbright+i);
		blevel[i] *= 255; bbright[i] *= 255;
	}
#endif
	for(i=0; i<nr-1; i++) {
		rslope[i] = (rbright[i+1] - rbright[i])/(rlevel[i+1] - rlevel[i]);
		rx[i]     = rbright[i] - rslope[i]*rlevel[i];
	}
	for(i=0; i<ng-1; i++) {
		gslope[i] = (gbright[i+1] - gbright[i])/(glevel[i+1] - glevel[i]);
		gx[i]     = gbright[i] - gslope[i]*glevel[i];
	}
	for(i=0; i<nb-1; i++) {
		bslope[i] = (bbright[i+1] - bbright[i])/(blevel[i+1] - blevel[i]);
		bx[i]     = bbright[i] - bslope[i]*blevel[i];
	}
	if( reverse ) {
		for(i=0; i<256; i++) {
			pickcol((unsigned char) i,&r,&g,&b);
			rtab[i] = 255 - r; gtab[i] = 255 - g; btab[i] = 255 - b;
		}
	}
	else {
		for(i=0; i<256; i++) {
			pickcol((unsigned char) i,&r,&g,&b);
			rtab[i] = r; gtab[i] = g; btab[i] = b;
		}
	}
}

void pickcol(gray, r, g, b)
unsigned char gray, *r, *g, *b;
{
	int ir, ig, ib;
	int rt, gt, bt;

	for(ir=0; ir < nr; ir++) { if ( rlevel[ir] > gray ) break; }
	for(ig=0; ig < ng; ig++) { if ( glevel[ig] > gray ) break; }
	for(ib=0; ib < nb; ib++) { if ( blevel[ib] > gray ) break; }
	ir--; ig--; ib--;
	if( ir >= nr-1 ) ir = nr-2;
	if( ig >= ng-1 ) ig = ng-2;
	if( ib >= nb-1 ) ib = nb-2;
	rt = rx[ir] + rslope[ir]*gray;
	gt = gx[ig] + gslope[ig]*gray;
	bt = bx[ib] + bslope[ib]*gray;
	if( rt > 255 ) rt = 255;
	if( gt > 255 ) gt = 255;
	if( bt > 255 ) bt = 255;
	if( rt < 0 ) rt = 0;
	if( gt < 0 ) gt = 0;
	if( bt < 0 ) bt = 0;
	*r = (unsigned char) rt;
	*g = (unsigned char) gt;
	*b = (unsigned char) bt;
}

