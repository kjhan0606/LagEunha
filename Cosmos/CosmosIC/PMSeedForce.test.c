#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<omp.h>
#include<mpi.h>
#include "eunha.h"
#include "flrw.h"
#include "mpiaux.h"


typedef struct Matterpower{
	char selectmatter[10];
}Matterpower;

Matterpower  matterpower;


void  qromb_( float (*)(float *), float *, float * , float *);
float bpowint_(float *);
float bpower_(float *);
void rinitpower_(float *, float *, float *, float *, float *, float *, float *, 
		float *, float *, float *, float *, int *, char *, int *);
void initpower_(float *, float *, float *, float *, float *, float *, float *, 
		float *, float *, float *, float *);
float gasdev(long *);


/*
float bpower(float *);
*/

float ran3(long*);

void MakeTunnel2Fortran(SimParameters *simpar, int itype){
	float rknyq, rkinv, powh,  omehs;
	float omep, omeplam, omepb, hubble, boxsize;
	int powreadflag;
	rknyq = 2*M_PI*NX(simpar)/2;
	rkinv = 1./rknyq;
	powh = 0.5*NPOW(simpar);
	omehs = 1./(OMEP(simpar)*HUBBLE(simpar)*BOXSIZE(simpar));
	omep = OMEP(simpar);
	omeplam = OMEPLAM(simpar);
	omepb = OMEPB(simpar);
	powreadflag = POWREADFLAG(simpar);
	hubble = HUBBLE(simpar);
	boxsize = BOXSIZE(simpar);
	void pkinfowrap4fortran_(float *, float *,  float *, float *, float *, float *, float *, float *, float *, int *, int *);
	pkinfowrap4fortran_(&hubble, &boxsize, &rknyq, &rkinv, &powh, &omehs, &omep, &omepb, &omeplam, &powreadflag,
			&itype);
}

void PMSeedForce(SimParameters *simpar, DenType *den, DenType *fx, DenType *fy, DenType *fz,
		DenType *den2, DenType *denpad, int itype){
	float damp1,damp2, vamp1, vamp2;
	ptrdiff_t i,j,k;
	ptrdiff_t nx = NX(simpar);
	ptrdiff_t mx = 2*(nx/2+1);
	ptrdiff_t ny = NY(simpar);
	ptrdiff_t nz = NZ(simpar);
	ptrdiff_t local_nz = LOCAL_NZ(simpar);
	ptrdiff_t local_z_start = LOCAL_Z_START(simpar);
	float bini, bnow;

	int nspace = NSPACE(simpar);
	float boxsize = COSBOXSIZE(simpar);
	float amax = AMAX(simpar);
	float npow = NPOW(simpar);
	float hubble = HUBBLE(simpar);
	long iseed = ISEED(simpar);
	float omepb = OMEPB(simpar);
	float omep = OMEP(simpar);
	float omeplam = OMEPLAM(simpar);
	char powfilename[190], asciipkfile[190];
	int seedflag;
	int nid = NID(simpar);
	int myid = MYID(simpar);
	float bias = BIAS8(simpar);
	float wlam0 = WLAM0(simpar);
	float wlam1 = WLAM1(simpar);
	MPI_Comm comm = MPI_COMM(simpar);
	MPI_Status status;
	float xnow, xini;
	int mlocalzstart, mlocalnz;



	sprintf(powfilename, "%s", POWFILENAME(simpar));
	sprintf(asciipkfile, "%s", INPAPKFILENAME(simpar));


	if(iseed >=0){
		seedflag = 1;
		iseed = -1;
	}
	else seedflag = 0;

	double rincube = 1.L/(double)nx/(double)ny/(double)nz;
	double rngc = (double)nx * (double)ny * (double)nz;
	float rnx, rny, rnz, xyz, xyy;
	rnx = nx;
	rny = ny;
	rnz = nz;
	float twopi = 2*M_PI;
	float rk;
	int knyq = nx/2;
	double rknyq2 = knyq*knyq;
	double rknyq  = twopi*knyq;
	double rkinv = 1./rknyq;
	float ai, zinit, anow;
	float powh = 0.5*npow;
	AI(simpar) = ai = 1;
	ZINIT(simpar) = zinit = amax-1;
	ANOW(simpar) = anow = ai;
	float omei,zp1;
	zp1 = 1 + zinit;

	OMEI(simpar) = omei = Omega_matter(simpar, 1.);
	float omehs = 1./(omep*hubble*boxsize);
	float pmas = 1./(nspace*nspace*nspace);

	/*
	float pmass = pmass/rngc;
	*/

	float wavemax = 100000;
	float pamp,w1,w2;
#ifdef EH98
	pamp = 0;
	w1 = 0.001*boxsize;
	for(i=0;i<5;i++){
		float pamp1;
		w2 = w1*10.;
		qromb_(bpowint_, &w1, &w2, &pamp1);
		pamp += pamp1;
		w1 = w2;
	}
#else
	MakeTunnel2Fortran(simpar, 0);
	sprintf(matterpower.selectmatter,"total");
#ifdef DEBUG
	if(MYID(simpar)==0) DEBUGPRINT("P%d: now I am here: matter: %s :: %d : %s\n",
			MYID(simpar), matterpower.selectmatter, POWREADFLAG(simpar), powfilename);
#endif
	if(POWREADFLAG(simpar) != 0 ){
		float asciirk[9000], asciipk[9000];
		int nascii;
		StartParallelRW(myid, (nid/10), comm);
		rinitpower_(&omep, &omepb, &omeplam, &wlam0, &hubble, &npow, &zinit, &bias, 
				&pamp, &boxsize, &rnx,  &POWREADFLAG(simpar),powfilename, &myid);
		CloseParallelRW(myid, nid, (nid/10), comm);
		if(POWREADFLAG(simpar) ==2){
			if(myid==0) {
				FILE *fp = fopen(INPAPKFILENAME(simpar),"r");
				int mmp;
				i = 0;
				while( (mmp = fscanf(fp,"%g %g\n", asciirk+i, asciipk+i))>0){
					if(i >8192) {
						DEBUGPRINT0("\n\nExit from reading external ascii powerspectrum\n\n");
						exit(9823);
					}
					i++;
				}
				fclose(fp);
				nascii = i;
			}
			MPI_Bcast(&nascii, 1, MPI_INT, 0, comm);
			MPI_Bcast(asciirk, nascii, MPI_FLOAT, 0, comm);
			MPI_Bcast(asciipk, nascii, MPI_FLOAT, 0, comm);
			pamp = 0;
			w1 = 0.001*boxsize;
			for(i=0;i<5;i++){
				float pamp1;
				w2 = w1*10.;
				qromb_(bpowint_, &w1, &w2, &pamp1);
				pamp += pamp1;
				w1 = w2;
			}
		}
	}
	else {
		initpower_(&omep, &omepb, &omeplam, &wlam0, &hubble, &npow, &zinit, &bias, &pamp, &boxsize, &rnx);
	}

#ifdef DEBUG
	DEBUGPRINT("P%d here %g %g \n",MYID(simpar), WLAM0(simpar), WLAM1(simpar)); MPI_Barrier(MPI_COMM(simpar));
#endif

	if(itype ==0) sprintf(matterpower.selectmatter,"total");
	if(itype ==1) sprintf(matterpower.selectmatter,"cdm");
	if(itype ==2) sprintf(matterpower.selectmatter,"baryon");
	else sprintf(matterpower.selectmatter,"total");

#ifdef DEBUG
	DEBUGPRINT("P%d here %g %g \n",MYID(simpar), WLAM0(simpar), WLAM1(simpar)); MPI_Barrier(MPI_COMM(simpar));
#endif
	MakeTunnel2Fortran(simpar, itype);
	if(POWREADFLAG(simpar) != 0){
		StartParallelRW(myid, (nid/10), comm);
		rinitpower_(&omep, &omepb, &omeplam, &wlam0, &hubble, &npow, &zinit, &bias,
				&pamp, &boxsize, &rnx, &POWREADFLAG(simpar), powfilename, &myid);
		CloseParallelRW(myid, nid, (nid/10), comm);
	}
	else {
		initpower_(&omep, &omepb, &omeplam, &wlam0, &hubble, &npow, &zinit, &bias, &pamp, &boxsize, &rnx);
	}
#endif


	if(myid ==0) printf(" <1>= %g \n", pamp);
	float omepk = 1. - omep - omeplam;

	if(WLAM1(simpar) !=0 ){
		bini = growthall(simpar,1.);
		bnow = growthall(simpar, amax);
		COS_PAMP(simpar,itype)  = pamp = (1./bias)/sqrt(pamp)*(bini/bnow);
		COS_DAMP1(simpar,itype) = damp1 = 1;
		COS_VAMP1(simpar,itype) = vamp1 = (growthall(simpar,1.+0.005)-
				growthall(simpar, 1.-0.005)) /(0.01)*1./bini;
	}
	else if(wlam0 != -1) {
		bini = growthgen(simpar, 1.);
		bnow = growthgen(simpar, amax);
		COS_PAMP(simpar,itype)  = pamp = (1./bias)/sqrt(pamp)*(bini/bnow);
		COS_DAMP1(simpar,itype) = damp1 = 1;
		COS_VAMP1(simpar,itype) = vamp1 = (growthgen(simpar, 1.+0.005)-
				growthgen(simpar,1.-0.005)) /(0.01)*1./bini;
	}
	/* II. SCDM Model */
	else if(omep ==1 && fabs(omepk) < 1.E-3){
		xnow = (1-omep)/omep;
		xini = xnow / (1+zinit);
		bini = 1.+3./xini+3.*sqrt(1.+xini)/pow(xini,1.5)*log(sqrt(1.+xini)-sqrt(xini));
		bnow = 1.+3./xnow+3.*sqrt(1.+xnow)/pow(xnow,1.5)*log(sqrt(1.+xnow)-sqrt(xnow));
		COS_PAMP(simpar,itype)  = pamp = (1./bias)/sqrt(pamp)*(bini/bnow);
		COS_DAMP1(simpar,itype) = damp1 = 1;
		COS_VAMP1(simpar,itype) = vamp1 = pow(omei, 0.6);
	}
	/* III. Zero curvature with cosmological constant model */
	else if(omep < 1 && fabs(omepk) < 1.E-3) {
		bini = growth(simpar, 1.);
		bnow = growth(simpar, amax);
		COS_PAMP(simpar,itype)  = pamp = (1./bias)/sqrt(pamp)*bini/bnow;
		COS_DAMP1(simpar,itype) = damp1 = 1;
		COS_VAMP1(simpar,itype) = vamp1 = (growth(simpar, 1.+0.005)-growth(simpar, 1.-0.005)) /(0.01)*1./bini;
	}
	/* IV. Non-zero curvature with cosmological constant model */
	else if(omeplam >0 && fabs(omepk) >1.E-3) {
		bini = growth2(simpar, 1.);
		bnow = growth2(simpar, amax);
		COS_PAMP(simpar,itype)  = pamp = (1./bias)/sqrt(pamp)*bini/bnow;
		COS_DAMP1(simpar,itype) = damp1 = 1.;
		COS_VAMP1(simpar,itype) = vamp1 = (growth2(simpar, 1.+0.005)-growth2(simpar, 1.-0.005))
			/(0.01)*1./bini;
	}
	if(myid ==0) printf("Omep= %g bini= %g bnow= %g pamp= %g\n",omep,bini, bnow, vamp1);

	float powerampoftotalmatter = pamp;

	float efold = 0.585*sqrt(2.);
	float sfac = -pow(M_PI*efold/boxsize,2);
	/* transform to the velocity potential */
	pamp = pamp * rngc * sqrt(0.5);
	/* P(k) = (pamp*bpower1(rk))**2 * boxsize**3; here rk is k and rk in the
	 * code is 2pi*u */

	COS_DAMP2(simpar,itype) = damp2 = -3./7.;
	COS_VAMP2(simpar,itype) = vamp2 = damp2*2*pow(omei, 4./7.);
#ifdef DEBUG
	DEBUGPRINT("P%d has amplitudes damp1/2= %g %g vamp1/2= %g %g for type %d\n", MYID(simpar), COS_DAMP1(simpar,itype),COS_DAMP2(simpar,itype),
			COS_VAMP1(simpar,itype), COS_VAMP2(simpar,itype),itype);
#endif

	if(MYID(simpar) == -1){
		int kkk = 1;
		while(kkk) {
			kkk = 2;
		}
	}


	long iiseed = iseed*(myid+1);

	if(seedflag !=0){
		for(k=0;k<local_nz;k++){
			ptrdiff_t kz = (k+ local_z_start);
			if(kz > nz/2) kz = kz - nz;
			for(j=0;j<ny;j++){
				ptrdiff_t ky = j;
				if(ky > ny/2) ky = ky - ny;
				for(i=0;i<=nx/2;i++){
					double xksq = i*i + ky*ky + kz*kz;
					if(xksq ==0 || xksq > rknyq2)
					{
						den[2*i+mx*(j+ny*k)] = 0;
						den[2*i+1+mx*(j+ny*k)] = 0;
					}
					else {
						rk = twopi*sqrt(xksq);
						den[2*i+mx*(j+ny*k)] = gasdev(&iiseed);
						den[2*i+1+mx*(j+ny*k)] = gasdev(&iiseed);
					}
				}
			} 
		}
	}
	else {
		for(k=0;k<local_nz;k++){ 
			ptrdiff_t kz = (k+ local_z_start); 
			if(kz > nz/2) kz = kz - nz; 
			for(j=0;j<ny;j++){
				ptrdiff_t ky = j;
				if(ky>ny/2) ky = ky - ny;
				for(i=0;i<=nx/2;i++){
					double xksq = i*i + ky*ky + kz*kz;
					if(xksq ==0 || xksq > rknyq2) {
						den[2*i+mx*(j+ny*k)] = 0;
						den[2*i+1+mx*(j+ny*k)] = 0;
					}
					else {
						rk = twopi*sqrt(xksq);
						float velpot = pamp*bpower_(&rk);
						/*
						if(isnan(velpot)) DEBUGPRINT("P%d has strange value %g::: %g %g\n",MYID(simpar), velpot, xksq, rk);
						*/
						den[2*i+mx*(j+ny*k)] = velpot*gasdev(&iiseed);
						den[2*i+1+mx*(j+ny*k)] = velpot*gasdev(&iiseed);
					}
				}
			}
		}
	}
	if(myid ==0) printf("last #= %g\n", ran3(&iiseed));


/*
cccc
c***  Make an array, den(*,*,local_nz+1:local_nz+2) which is needed to
ccc   follow the complex conjugate relation for a real array.
ccc   mpi_gatherv is gathering values with displacement of local_nz.
ccc   0 is the root node number.
ccc   Gathering all the den(*,*,local_nz+1:local_nz+2) to 
ccc     den(*,*,local_nz+1:local_nz+2) of 0'th processor.
ccc   rcounts is the size of sending buffer and also displs is displacement of
ccc   receive buffer saved to 0'th node.

c---  Satisfy complex conjugation relation for a real array
ccc   Only on rank 0 processor.
 * */
	for(k=0;k<local_nz;k++) for(j=0;j<ny;j++) for(i=0;i<2;i++){
		denpad[i+2*(j+ny*(k+local_z_start))] = den[i+mx*(j+ny*k)];
	}
	if(myid==0){
		for(i=1;i<nid;i++){
			ptrdiff_t mlocalzstart, mlocalnz;
			MPI_Recv(&mlocalzstart, 1, MPI_LONG, i, 0, comm, &status);
			MPI_Recv(&mlocalnz, 1, MPI_LONG, i, 1, comm, &status);
			MPI_Recv(denpad+2*ny*mlocalzstart, 2*ny*mlocalnz, MPI_DEN_T, i, 2, comm, &status);
		}
	}
	else {
		MPI_Send(&local_z_start, 1, MPI_LONG, 0 , 0, comm);
		MPI_Send(&local_nz, 1, MPI_LONG, 0 , 1, comm);
		MPI_Send(denpad+2*ny*local_z_start, 2*ny*local_nz, MPI_DEN_T, 0 , 2, comm);
	}
	if(myid==0){
		for(k=1; k<nz/2;k++){
			denpad[  2*ny*(nz-k)] =  denpad[  2*ny*k];
			denpad[1+2*ny*(nz-k)] = -denpad[1+2*ny*k];
			for(j=1;j<ny/2;j++){
				denpad[  2*(ny-j+ny*(nz-k))] =  denpad[  2*(j+ny*k)];
				denpad[1+2*(ny-j+ny*(nz-k))] = -denpad[1+2*(j+ny*k)];
			}
		}
		for(j=1;j<ny/2;j++){
			denpad[  2*(ny-j)] =  denpad[  2*j];
			denpad[1+2*(ny-j)] = -denpad[1+2*j];
		}
		for(k=1; k<nz/2;k++){
			for(j=ny/2+1;j<ny;j++){
				denpad[  2*(ny-j+ny*(nz-k))] =  denpad[  2*(j+ny*k)];
				denpad[1+2*(ny-j+ny*(nz-k))] = -denpad[1+2*(j+ny*k)];
			}
		}
	}
	int sendsize = 2*ny*nz;
	MPI_Bcast(denpad, sendsize, MPI_DEN_T, 0, comm);
	for(k=0;k<local_nz;k++) for(j=0;j<ny;j++) for(i=0;i<2;i++){
		den[i+mx*(j+ny*k)] = denpad[i+2*(j+ny*(local_z_start+k))];
	}
#ifdef DEBUG
	if(1){
		fftwf_mpi_execute_dft_c2r(FFTW_NB_PLAN(simpar),(fftwf_complex*)den, den);
		for(k=0;k<local_nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++){
			den[i+mx*(j+ny*k)] *= rincube;
		}
#define MIN(a,b) ( (a)<(b) ? (a): (b) )
#define MAX(a,b) ( (a)>(b) ? (a): (b) )
		double valmin,valmax;
		valmin = 1.e20L; valmax= -1e20L;
		for(i=0;i<nx*ny*local_nz;i++){
			valmin = MIN(valmin,den[i]);
			valmax = MAX(valmax,den[i]);
		}
		DEBUGPRINT("P%d has assigned value %lg %lg\n", MYID(simpar), valmin,valmax);
		FILE *wp;
		if(MYID(simpar)==0){
			wp = fopen("IDen.dat","w");
			int nnx,nny,nnz;
			nnx = nx;
			nny = ny;
			nnz = nz;
			fwrite(&nnx,sizeof(int),1,wp);
			fwrite(&nny,sizeof(int),1,wp);
			fwrite(&nnz,sizeof(int),1,wp);
			fclose(wp);
		}
		for(i=0;i<NID(simpar);i++){
			if(MYID(simpar)==i){
				wp = fopen("IDen.dat","a");
				for(k=0;k<local_nz;k++) for(j=0;j<ny;j++){
					fwrite(den+mx*(j+ny*k),sizeof(float), nx,wp);
				}
				fclose(wp);
			}
			MPI_Barrier(MPI_COMM(simpar));
		}
		MPI_Finalize();
		exit(99);
	}
#endif

	MPI_Barrier(MPI_COMM(simpar));
	if(MYID(simpar)==0) DEBUGPRINT0("just after broadcasting the conjugate information\n");
	float growthfactor;
	COS_GROWTH(simpar)  = growthfactor = bini/bnow;

	if(PORDER(simpar) ==2) {
		void TwoLPT(SimParameters *, DenType *, DenType *, DenType *, DenType *);
		TwoLPT(simpar, den, fx,fy,den2);
	}
	else if(PORDER(simpar) ==1) {
		void Zeldovich(SimParameters *, DenType *, DenType *, DenType * , DenType *);
		Zeldovich(simpar, den, fx,fy,fz);
	}
}
