#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdarg.h>
#include <math.h>

#include "color.h"


typedef struct Img{
	int reverse;
	int nx, ny;
	float nlog, ccdmin0, ccdmax0;
	float *ccd;
	unsigned char *r,*g,*b, *work;
	char *saofile;
} Img;

void f2pgm(Img *);
void pgm2col(Img *); 
void pickcol(unsigned char, unsigned char *, unsigned char *, unsigned char*,
		int,int,int, float *, float *, float *, float *, float *, float *,
		float *, float *, float*);
void readpgm(unsigned char *, int , int , unsigned char *, unsigned char *, unsigned char *);
void readsao(int , FILE *);

void finalizeImg(Img *img){
	free(img->r);
	free(img->g);
	free(img->b);
	free(img->work);
}

void initializeImg(Img *img, float *map, int nx, int ny, char *saofile){
	img->nx = nx;
	img->ny = ny;
	img->r = (unsigned char*)malloc(sizeof(unsigned char)*nx*ny);
	img->g = (unsigned char*)malloc(sizeof(unsigned char)*nx*ny);
	img->b = (unsigned char*)malloc(sizeof(unsigned char)*nx*ny);
	img->work = (unsigned char*)malloc(sizeof(unsigned char)*nx*ny);
	img->ccd = map;
	img->ccdmin0 = img->ccdmax0 = 9999.;
	img->nlog = 13;
	img->reverse = 0;
	img->saofile = saofile;
//	snprintf(img->saofile, sizeof(img->saofile), "%s", saofile);
}

void test_colorizeit( float *map, int nx, int ny, char *saofile, char *outfile){
	/*
	Img gas;
	int i,j,k;
	int nsize;

	printf("saofile/outfile %s : %s\n", saofile, outfile);
	nsize = nx*ny;

	{
		initializeImg(&gas, map, nx, ny, saofile);
		gas.ccdmax0 = 2;
		gas.ccdmin0 = 1.;
		gas.nlog = -1;
		gas.nx = nx;
		gas.ny = ny;

		gas.ccdmin0 = 1;
		gas.ccdmax0 = 2;
		f2pgm(&gas);
		pgm2col(&gas);
		printf("ending pgm2col\n");
	
//		r = (unsigned int*)malloc(sizeof(unsigned int)*nsize);
//		g = (unsigned int*)malloc(sizeof(unsigned int)*nsize);
//		b = (unsigned int*)malloc(sizeof(unsigned int)*nsize);

		for(i=0;i<nsize;i++){
			if(gas.r[i]> 255) gas.r[i] = 255;
			if(gas.g[i]> 255) gas.g[i] = 255;
			if(gas.b[i]> 255) gas.b[i] = 255;

//			r[i] = gas.r[i];
//			g[i] = gas.g[i];
//			b[i] = gas.b[i];
		}
		printf("ending cap\n");
		FILE *out = fopen(outfile,"w");
		fprintf(out,"P6\n%d %d\n255\n",nx,ny);
		for(i=0;i<nsize;i++){
			fwrite(&(gas.r[i]),1,1,out);
			fwrite(&(gas.g[i]),1,1,out);
			fwrite(&(gas.b[i]),1,1,out);
		}
		fclose(out);
		printf("ending dump\n");
		finalizeImg(&gas);
		printf("finalizing the colorizing\n");
	}
	*/

}
void colorizeit(
		int iarg, ...
		){
	va_list ap;
	va_start(ap,iarg);
	float *map = va_arg(ap, float*);
	int nx = va_arg(ap, int);
	int ny = va_arg(ap, int);
	char saofile[190]="\0";
	char outfile[190]="\0"; 
	char *s = va_arg(ap, char *);
	strcat(saofile,s);
	s = va_arg(ap, char *);
	strcat(outfile,s);


	Img gas, star, dm, temp, composit, smetal;

	unsigned int *r,*g, *b;
	int i,j,k;
	char infile[190];
	FILE *fp;
	int nsize;


//	ImgHeader(&nx,&ny, outfile);
	nsize = nx*ny;

	{
		initializeImg(&gas, map, nx, ny, saofile);
		gas.ccdmax0 = 2;
		gas.ccdmin0 = 1.;
		gas.nlog = -1;
		gas.nx = nx;
		gas.ny = ny;

		{
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
		f2pgm(&gas);
		pgm2col(&gas);
	
//		r = (unsigned int*)malloc(sizeof(unsigned int)*nsize);
//		g = (unsigned int*)malloc(sizeof(unsigned int)*nsize);
//		b = (unsigned int*)malloc(sizeof(unsigned int)*nsize);

		for(i=0;i<nsize;i++){
			if(gas.r[i]> 255) gas.r[i] = 255;
			if(gas.g[i]> 255) gas.g[i] = 255;
			if(gas.b[i]> 255) gas.b[i] = 255;

//			r[i] = gas.r[i];
//			g[i] = gas.g[i];
//			b[i] = gas.b[i];
		}
		FILE *out = fopen(outfile,"w");
		fprintf(out,"P6\n%d %d\n255\n",nx,ny);
		for(i=0;i<nsize;i++){
			fwrite(&(gas.r[i]),1,1,out);
			fwrite(&(gas.g[i]),1,1,out);
			fwrite(&(gas.b[i]),1,1,out);
		}
		fclose(out);
		finalizeImg(&gas);
	}

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
	fprintf(stdout,"ccdmax %g , ccdmin %g nx/ny %d/%d\n",ccdmax,ccdmin,nx,ny);
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


static unsigned char rtab[256], gtab[256], btab[256];

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
	fclose(colfile);
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


void readsao(int reverse, FILE *saofile) {
	int nr, ng, nb;
	int i;
	int ir, ig, ib, len;
	int n0;
	unsigned char r, g, b;
	char line[80], *p, n[15], *p0;
	float *all = (float*)malloc(sizeof(float)*500*12);
	float *rlevel, *rbright;
	float *glevel, *gbright;
	float *blevel, *bbright;
	float *rslope, *gslope, *bslope;
	float *rx, *gx, *bx;
	rlevel = all;
	glevel = rlevel + 500;
	blevel = glevel + 500;
	rbright = blevel + 500;
	gbright = rbright + 500;
	bbright = gbright + 500;
	rslope = bbright + 500;
	gslope = rslope + 500;
	bslope = gslope + 500;
	rx = bslope + 500;
	gx = rx + 500;
	bx = gx + 500;

	/*
	float rlevel[500], rbright[500];
	float glevel[500], gbright[500];
	float blevel[500], bbright[500];
	float rslope[500], rx[500];
	float gslope[500], gx[500];
	float bslope[500], bx[500];
	*/

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
			pickcol((unsigned char) i,&r,&g,&b, nr,ng,nb,
					rlevel, glevel, blevel, rx,gx,bx, rslope, gslope, bslope);
			rtab[i] = 255 - r; gtab[i] = 255 - g; btab[i] = 255 - b;
		}
	}
	else {
		for(i=0; i<256; i++) {
			pickcol((unsigned char) i,&r,&g,&b, nr,ng,nb,
					rlevel, glevel, blevel, rx,gx,bx, rslope, gslope, bslope);
			rtab[i] = r; gtab[i] = g; btab[i] = b;
		}
	}
	free(all);
}

void pickcol(unsigned char gray, unsigned char *r, unsigned char *g, unsigned char *b, 
		int nr, int ng, int nb, 
		float *rlevel, float *glevel, float *blevel, float *rx,float *gx,float *bx, 
		float *rslope, float *gslope, float *bslope)
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

