#include "voro.h"
#include "sedov.h"

#define Nneigh 50
Voro3D_point neigh[Nneigh];

int makemap(Voro3D_GasParticle *, int, int);

void outfile(Voro3D_GasParticle *bp, int np, int nstep, ptype t, ptype dt){
	char outfile[190];
	sprintf(outfile,"sedovout.%.6d.dat",nstep);
	FILE *wp = fopen(outfile,"w");
	fwrite(&np, sizeof(int), 1, wp);
	fwrite(bp, sizeof(Voro3D_GasParticle), np, wp);
	fwrite(&t, sizeof(ptype), 1, wp);
	fwrite(&dt, sizeof(ptype), 1, wp);
	fclose(wp);
}


Voro3D_GasParticle *readdata(int *np, ptype *t, ptype *dt,int nstep){
	char infile[190];
	sprintf(infile,"sedovout.%.6d.dat", nstep);
	FILE *fp = fopen(infile,"r");
	fread(np, sizeof(int), 1, fp);
	Voro3D_GasParticle *res = (Voro3D_GasParticle*)malloc(sizeof(Voro3D_GasParticle)*(*np));
	fread(res,sizeof(Voro3D_GasParticle), *np,fp);
	fread(t,sizeof(ptype), 1,fp);
	fread(dt,sizeof(ptype), 1,fp);
	fclose(fp);
	return res;
}


int main(int argc, char **argv){
	ptype t,dt;
	int np;
	Voro3D_GasParticle *bp;
	int icount = 0;;
	int iflag,jflag;



	if(argc ==1) {
		t = 0;
		bp = mkinitial(&np);
	}
	else if(argc ==2){
		ptype dt;
		int nstep = atoi(argv[1]);
		bp = readdata(&np,&t,&dt, nstep);
		icount = nstep+1;
	}

	ptype avgvolume = (Lx*Ly*Lz)/np;

	do {
		ptype dt;
//		MkLinkedList(bp,np);
		if(icount%10 ==0){
			double vph3Dcentroid(Voro3D_GasParticle **, int *, ptype );
			dt = vph3Dcentroid(&bp,&np, avgvolume);
		}
		else {
			dt = vph3D(&bp,&np, avgvolume);
		}
		bp[Nxp/2+Nxp*(Nyp/2+Nyp*(Nzp/2))].x = Lx*0.5L;
		bp[Nxp/2+Nxp*(Nyp/2+Nyp*(Nzp/2))].y = Ly*0.5L;
		bp[Nxp/2+Nxp*(Nyp/2+Nyp*(Nzp/2))].z = Lz*0.5L;
		t += dt;
		printf("Time is %g with icount = %d\n",t, icount);
		fflush(stdout);


		iflag = t * 100000.;
		jflag = (t-dt) * 100000.;


		if(icount%50 ==0) 
//		if(iflag != jflag)		
		{
			outfile(bp, np, icount, t, dt);
//			makemap(bp, np, icount);
		}


		icount ++;

	} while(t<1.e9);
}
