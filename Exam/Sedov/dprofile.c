#include "voro.h"
#include "sedov.h"


int main(int argc,char **argv){
	int i,j,k;
	int np;
	int nlog = 4;
	ptype time,dt;
	FILE *fp = fopen(argv[1],"r");
	FILE *wp = fopen(argv[2],"w");
	ptype xmin,ymin,xmax,ymax;
	fread(&np,sizeof(int),1,fp);
	Voro3D_GasParticle *bp = (Voro3D_GasParticle*)malloc(sizeof(Voro3D_GasParticle)*(np));
    fread(bp,sizeof(Voro3D_GasParticle), np,fp);
    fread(&time,sizeof(ptype), 1,fp);
    fread(&dt,sizeof(ptype), 1,fp);
    fclose(fp);

	ptype dfact = Lx*Ly*Lz/Nxp/Nyp/Nzp/1270.; 
	dfact = 1;
	printf("Density factor is %g\n", dfact); 
	Voro3D_point center;
	center.x = Lx/2 + 0.5*Lx/Nxp;
	center.y = Ly/2 + 0.5*Ly/Nyp;
	center.z = Lz/2 + 0.5*Lz/Nzp;

	for(i=0;i<np;i++){
		ptype tmpx, tmpy, tmpz;
		tmpx = (bp[i].x-center.x);
		tmpy = (bp[i].y-center.y);
		tmpz = (bp[i].z-center.z);
		ptype dist;
		dist = sqrt(tmpx*tmpx+tmpy*tmpy+tmpz*tmpz);
		if(dist < 1.) fprintf(wp,"%g %g\n", dist, bp[i].den*dfact);
	}
	fclose(wp);
	return 0;
}

