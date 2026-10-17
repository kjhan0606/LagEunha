#include "voro.h"



int main(int argc, char **argv){
	int i,j,k;
	Voro2D_point  *r;
	int np,mp, *id;

	int *iface;
	double *area, *norm;


	np = 1000;
	r = (Voro2D_point*)malloc(sizeof(Voro2D_point)*np);

	if(1){
		np = 10;
		for(i=0;i<np;i++){
			r[i].x = (((i+1)*942321)%2013)/2013.L * 2 - 1;
			r[i].y = (((i+1)*122321)%9014)/9014.L * 2 - 1;
			r[i].id = i + 1;
		}
	}
	else if(0){
		np = 6;
	}

//	mp = view_myvoro(np, id, x,y,z, iface, area, norm);
//	mp = myvoro(np, id, x,y,z, iface, area, norm);
	Voro2D_node voronodes[np*3];
	mp = np* 3;
	ptype boxsize = 1;
	int ip = Voro2D_FindPolyhedra(r,r,np, voronodes,mp,boxsize);
	printf("%d\n", ip);
}
