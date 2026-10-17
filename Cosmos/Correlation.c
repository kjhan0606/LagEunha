#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
//#include "mpirms.h"
#include "flags.h"



void CorrMeasure(SimParameters *simpar, GridInfo *grid){
	MPI_Comm comm = MPI_COMM(simpar);
	DenType *den = (DenType*)(grid+1);
	ptrdiff_t nx = NX(simpar);
	ptrdiff_t ny = NY(simpar);
	ptrdiff_t nz = NZ(simpar);
	ptrdiff_t mx = 2*(nx/2+1);
	ptrdiff_t local_ny_after_transpose = LOCAL_NY_AFTER_TRANSPOSE(simpar);
	ptrdiff_t local_y_start_after_transpose = LOCAL_Y_START_AFTER_TRANSPOSE(simpar);
	ptrdiff_t i,j,k, ngh;
	double rng = nx;
	ngh = nx/2;
	if(nx != ny || nx != nz){
		DEBUGPRINT0("###############################\n");
		DEBUGPRINT("Error nx/ny/nz = %ld %ld %ld\n",nx,ny,nz);
		DEBUGPRINT0("###############################\n");
		MPI_Finalize();exit(0);
	}
	double rngc = nx*ny*nz;
	double rngm6 = 1.L/(rngc*rngc);
	double pk[nx+1], npk[nx+1], tpk[nx+1], tnpk[nx+1], cor[nx+1];
	double tscdenconvx, tscdenconvy, tscdenconvz;

	for(i=0;i<nx+1;i++) pk[i] = npk[i] = tpk[i] = tnpk[i] = cor[i] = 0;

	for(j=0;j<local_ny_after_transpose;j++){
		ptrdiff_t wj = j + local_y_start_after_transpose;
		if(wj > ngh) wj -= ny;
		if(wj ==0) tscdenconvy = 1;
		else  {
			tscdenconvy = sin(M_PI*wj/ny)/(M_PI*wj/ny);
		}
		tscdenconvy = tscdenconvy*tscdenconvy*tscdenconvy;
		for(k=0;k<nz;k++){
			ptrdiff_t wk = k;
			if(wk > ngh) wk -= nz;
			if(wk ==0) tscdenconvz = 1;
			else {
				tscdenconvz = sin(M_PI*wk/nz)/(M_PI*wk/nz);
			}
			tscdenconvz = tscdenconvz*tscdenconvz*tscdenconvz;
			for(i=0;i<ngh;i++){
				double wave = sqrt(i*i + wj*wj + wk*wk);
				if(i==0) tscdenconvx = 1;
				else {
					tscdenconvx = sin(M_PI*i/(double)nx)/(M_PI*i/(double)nx);
				}
				tscdenconvx = tscdenconvx*tscdenconvx*tscdenconvx;
				double delksq = (den[2*i+mx*(k+nz*j)]*den[2*i+mx*(k+nz*j)]+
					den[2*i+1+mx*(k+nz*j)]*den[2*i+1+mx*(k+nz*j)]) * rngm6;
				double tscdenconv = tscdenconvx*tscdenconvy*tscdenconvz;
				delksq = delksq/(tscdenconv*tscdenconv);
				ptrdiff_t kbin1 = wave;
				double wt1 = wave-kbin1;
				ptrdiff_t kbin2 = kbin1 + (wt1>0? 1:-1); 
				wt1 = fabs(wt1);
				if(kbin1 <= nx) {
					pk[kbin1] += (1-wt1)*delksq;
					npk[kbin1] += (1-wt1);
				}
				if(kbin2 <= nx){
					pk[kbin2] += wt1*delksq;
					npk[kbin2] += wt1;
				}
			}
		}
	}
	MPI_Reduce(pk,tpk, nx+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM(simpar) );
	MPI_Reduce(npk,tnpk, nx+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM(simpar) );
	int now=0;
	if(MYID(simpar)==0){
		double rth = 8.L/BOXSIZE(simpar)*nx;
		double std = 0;
		double fact = 2*M_PI*rth/rng;
		double corn1,corn2;
		for(i=1;i<nx+1;i++){
			if(tnpk[i]>0) tpk[i] /= tnpk[i];
			else tpk[i] = 0;
			double thk = fact*i;
			std += tpk[i]*i*i *pow(sin(thk)-thk*cos(thk),2.L)/pow(thk,6.L);
		}
		std = sqrt(4*M_PI*std*9.L);
		for(i=1;i<nx+1;i++){
			double seper = 2*M_PI*i/rng;
			for(k=1;k<nz+1;k++){
				cor[i] += sin(k*seper)*k*tpk[k];
			}
			cor[i] = (rng+rng)*cor[i]/i;
		}
		if(cor[3] >0 && cor[2] >0) {
			corn2 = (log(cor[3]) - log(cor[2])) / (log(4.L)- log(3.L));
			corn1 = corn2;
		}

		FILE *wp;
		if(STEPCOUNT(simpar)==1 && ANOW(simpar) == 1) wp= fopen("pmevol.map", "w");
		else wp= fopen("pmevol.map", "a");
		fprintf(wp,"%g %g %g\n",ANOW(simpar), std, corn1);
		for(i=0;i<ngh;i++){
			fprintf(wp,"%5d %5d %15.7g %15.7g\n", (int)(i+1), (int)i, tpk[i], cor[i]);
		}
		fclose(wp);
	}
}
