#include "voro.h"
#include "kh.h"

#define Nneigh 50
Voro2D_point neigh[Nneigh];

int makemap(Voro2D_GasParticle *, int, int);

void outfile(Voro2D_GasParticle *bp, int np, int nstep, ptype t, ptype dt){
	char outfile[190];
	sprintf(outfile,"khout.%.6d.dat",nstep);
	FILE *wp = fopen(outfile,"w");
	fwrite(&np, sizeof(int), 1, wp);
	fwrite(bp, sizeof(Voro2D_GasParticle), np, wp);
	fwrite(&t, sizeof(ptype), 1, wp);
	fwrite(&dt, sizeof(ptype), 1, wp);
	fclose(wp);
}


Voro2D_GasParticle *readdata(int *np, ptype *t, ptype *dt,int nstep){
	char infile[190];
	sprintf(infile,"khout.%.6d.dat", nstep);
	FILE *fp = fopen(infile,"r");
	fread(np, sizeof(int), 1, fp);
	Voro2D_GasParticle *res = (Voro2D_GasParticle*)malloc(sizeof(Voro2D_GasParticle)*(*np));
	fread(res,sizeof(Voro2D_GasParticle), *np,fp);
	fread(t,sizeof(ptype), 1,fp);
	fread(dt,sizeof(ptype), 1,fp);
	fclose(fp);
	return res;
}


int main(int argc, char **argv){
	ptype t,dt;
	int np;
	Voro2D_GasParticle *bp;
	int icount = 0;;
	int iflag,jflag;



	if(argc ==1) {
		printf("making initial conditions\n");
		t = 0;
		bp = mkinitial(&np);
		printf("made initial conditions\n");
	}
	else if(argc ==2){
		ptype dt;
		int nstep = atoi(argv[1]);
		bp = readdata(&np,&t,&dt, nstep);
		icount = nstep+1;
	}

	do {
//		MkLinkedList(bp,np);
		ptype dt = vph2D(bp,np);
		t += dt;
		printf("Time is %g with icount = %d\n",t, icount);
		fflush(stdout);


		iflag = t * 10.;
		jflag = (t-dt) * 10.;


//		if(icount%50 ==0) 
		if(iflag != jflag)		
		{
			outfile(bp, np, icount, t, dt);
			makemap(bp, np, icount);
		}


		icount ++;

	} while(t<10.);
}
