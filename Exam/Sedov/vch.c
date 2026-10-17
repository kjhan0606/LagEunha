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
		if(nstep ==0){
			ptype meanvol = Lx*Ly*Lz/Nxp/Nyp/Nzp;
			int i;
			for(i=0;i<np;i++){
				bp[i].vx = bp[i].vy = bp[i].vz = 0.;
				bp[i].ax = bp[i].ay = bp[i].az = 0.;
				bp[i].pressure  = bckPressure; 
				bp[i].den  = 1;
				bp[i].ie  = bp[i].pressure *bp[i].volume/(Gamma-1); 
				bp[i].mass  = bp[i].den * bp[i].volume;
			}
			int iblast = 0;
			ptype sumfact = 0;
			for(i=0;i<np;i++){
				ptype tmpx = bp[i].x - Lx/2.;
				ptype tmpy = bp[i].y - Ly/2.;
				ptype tmpz = bp[i].z - Lz/2.;
				ptype dist = sqrt(tmpx*tmpx+tmpy*tmpy+tmpz*tmpz);
				if(dist < 0.02) {
					bp[i].ie = exp(-dist*dist/(0.0070*0.0070));
					sumfact += exp(-dist*dist/(0.0070*0.0070));
					iblast ++;
				}
			}
			for(i=0;i<np;i++){
				ptype tmpx = bp[i].x - Lx/2.;
				ptype tmpy = bp[i].y - Ly/2.;
				ptype tmpz = bp[i].z - Lz/2.;
				ptype dist = sqrt(tmpx*tmpx+tmpy*tmpy+tmpz*tmpz);
				if(dist < 0.02) {
					bp[i].ie = bp[i].ie/sumfact * eBlast;
					bp[i].pressure = bp[i].ie/meanvol*(Gamma-1);
				}
			}
			printf("Total %d particles are blown \n", iblast);

//			int icenter = Nxp/2 + Nxp*(Nyp/2 + Nyp*(Nzp/2));
//			bp[icenter].ie = eBlast;
			t = 0;
			dt = 0;
		}
		icount = nstep+1;
	}

	ptype avgvolume = (Lx*Ly*Lz)/np;

	do {
//		MkLinkedList(bp,np);
		ptype dt = vph3D(&bp,&np, avgvolume);
		t += dt;
		printf("Time is %g with icount = %d\n",t, icount);
		fflush(stdout);


		iflag = t * 100000.;
		jflag = (t-dt) * 100000.;


//		if(icount%100 ==0) 
		if(iflag != jflag)		
		{
			outfile(bp, np, icount, t, dt);
//			makemap(bp, np, icount);
		}


		icount ++;

	} while(t<10.);
}
