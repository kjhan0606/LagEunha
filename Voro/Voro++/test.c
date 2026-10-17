#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>


int myvoro(int, int *, double *, double *, double *, int *, double *, double *);
int view_myvoro(int, int *, double *, double *, double *, int *, double *, double *);

int main(int argc, char **argv){
	int i,j,k;
	double  *x,*y,*z;
	int np,mp, *id;

	int *iface;
	double *area, *norm;


	np = 1000;
	x = (double*)malloc(sizeof(double)*np);
	y = (double*)malloc(sizeof(double)*np);
	z = (double*)malloc(sizeof(double)*np);
	area = (double*)malloc(sizeof(double)*np);
	norm = (double*)malloc(sizeof(double)*np*3);
	iface = (int*)malloc(sizeof(int)*np);
	id = (int*)malloc(sizeof(int)*np);

	if(1){
		np = 50;
		for(i=0;i<np;i++){
			x[i] = (((i+1)*942321)%2013)/2013.L * 2 - 1;
			y[i] = (((i+1)*122321)%9014)/9014.L * 2 - 1;
			z[i] = (((i+1)*382321)%7019)/7019.L * 2 - 1;
			id[i] = i + 1;
		}
	}
	else if(0){
		np = 6;
		x[0] = y[0] = 0; z[0] = 0.5; id[0] =  1;
		x[1] = y[1] = 0; z[1] = -0.5; id[1] = 2;
		x[2] = z[2] = 0; y[2] = 0.5; id[2] = 3;
		x[3] = z[3] = 0; y[3] = -0.5; id[3] = 4;
		z[4] = y[4] = 0; x[4] = 0.5; id[4] = 5;
		z[5] = y[5] = 0; x[5] = -0.5; id[5] = 6;
	}

//	mp = view_myvoro(np, id, x,y,z, iface, area, norm);
	mp = myvoro(np, id, x,y,z, iface, area, norm);
	printf("%d\n", mp);
}
