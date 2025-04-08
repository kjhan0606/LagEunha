#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "../OST/ost.h"

#define NSPLINE (1<<12)

typedef double FORCEARRAYTYPE;
typedef struct DIFFFORCE{ 
	FORCEARRAYTYPE diff[3]; 
	FORCEARRAYTYPE slope[3];
} DIFFFORCE;

/* These are the GPU settings.                           */
/* ##########     WARNING       ########################
 *  * ##########     This should be calculated manually. ##
 *   * ##################################################### */
#define RANGE 6
#define MIN(a,b) ((a)>(b) ? (b): (a))
#define RMAX 6.
#define RMAXSQ 36.
#define LOG10RMAX (0.7781512) /* = alog10(RMAX) */
#define LOG10RMIN (-3.0) /* log10(RMIN) */
#define RMINSQ (1.E-6) /* == RMIN*RMIN */
#define LOG10RMAXmLOG10RMIN (3.7781512) /* == log10(RMAX)-log10(RMIN) */
#define invLOG10RMAXmLOG10RMIN (0.264680) /* == 1./(log10(RMAX)-log10(RMIN)) */



/*
#ifdef DEFINE_TREE_COMMON
#define EXTERN
#else
#define EXTERN extern
#endif
*/
/*
 * EXTERN  DIFFFORCE forcecorrect[NSPLINE];
 * #define (forcecorrectdiff(i,j)) forcecorrect[i].diff[j]
 * #define (forcecorrectslope(i,j)) forcecorrect[i].slope[j]
 * */
/*
EXTERN  FORCEARRAYTYPE diff[NSPLINE][3],slope[NSPLINE][3];

EXTERN  TREEPRECISION ran2nran,invran2nran;
#undef EXTERN
*/
FORCEARRAYTYPE diff[NSPLINE][3],slope[NSPLINE][3];
#define forcecorrectdiff(a,b) diff[a][b]
#define forcecorrectslope(a,b) slope[a][b]





double ran2nran,invran2nran;

void inline gettreeforce(ax,ay,az, npar, xpar,ypar,zpar,mpar,nnode,xnode,ynode,znode,xxnode,yynode,zznode,
		xynode,xznode,yznode,trnode,mnode)
float *ax,*ay,*az;
int npar, nnode; 
float *xpar,*ypar,*zpar,*mpar; 
float *xnode,*ynode,*znode; 
float *xxnode,*yynode, *zznode, *xynode,*xznode,*yznode,*trnode,*mnode; 
{
	int i,j,k;
	float dist;
	float aax,aay,aaz;
	aax = aay = aaz = 0;
#pragma simd
#pragma vector aligned
	for(i=0;i<npar;i++){
		float dist,fplmf;
		int ntmp;
		dist = sqrt(xpar[i]*xpar[i] + ypar[i]*ypar[i] + zpar[i]*zpar[i]);
		ntmp = dist*invran2nran;
		fplmf = forcecorrectdiff(ntmp,0) + forcecorrectslope(ntmp,0) * (dist-ntmp*ran2nran);
		fplmf *= mpar[i];
		aax += xpar[i]* fplmf;
		aay += ypar[i]* fplmf;
		aaz += zpar[i]* fplmf;
	}
#pragma simd
#pragma vector aligned
	for(i=0;i<nnode;i++){
		float dist;
		int ntmp;
		dist = sqrt(xpar[i]*xpar[i] + ypar[i]*ypar[i] + zpar[i]*zpar[i]);
		ntmp = dist*invran2nran;
		float qxx1,qxx2,qxx3,qxy, tmpdist,fplmf1,fplmf2,fplmf3, tmptmp, twofplmf2;
		qxx1 = xxnode[i]*xnode[i] + xynode[i]*ynode[i] + xznode[i]*znode[i];
		qxx2 = xynode[i]*xnode[i] + yynode[i]*ynode[i] + xznode[i]*znode[i];
		qxx3 = xznode[i]*xnode[i] + yznode[i]*ynode[i] + xznode[i]*znode[i];
		qxy = xxnode[i] * xnode[i]*xnode[i] + yynode[i]*ynode[i]*ynode[i] + zznode[i]*znode[i]*znode[i] +
			2*(xynode[i]*xnode[i]*ynode[i]+xznode[i]*xnode[i]*znode[i]+yznode[i]*ynode[i]*znode[i]);
		tmpdist = (dist-ntmp*ran2nran);
		ntmp = MIN(ntmp, NSPLINE-1);
		fplmf1 = forcecorrectdiff(ntmp,0) + forcecorrectslope(ntmp,0)*tmpdist;
		fplmf2 = forcecorrectdiff(ntmp,1) + forcecorrectslope(ntmp,1)*tmpdist;
		fplmf3 = forcecorrectdiff(ntmp,2) + forcecorrectslope(ntmp,2)*tmpdist;
		tmptmp = mnode[i]*fplmf1 + trnode[i]*fplmf2 + qxy*fplmf3;
		twofplmf2 = 2*fplmf2;
		aax += xnode[i]*tmptmp+qxx1*twofplmf2;
		aay += ynode[i]*tmptmp+qxx2*twofplmf2;
		aaz += znode[i]*tmptmp+qxx3*twofplmf2;
	}
	*ax = aax;
	*ay = aay;
	*az = aaz;
}



double LL;
double AA[8] = {
	1.L/256.L/256.L, 
	0.78224E+00L, 
	0.37971E-06L, 
	-0.60338E+00L, 
	-0.34419E-07L, 
	-0.39741E+01L, 
	-0.10607E+01L, 
	-0.38145E+00L
};


double sNewton(double x, double epsilon){
	double e,g;
	g = -1.L/sqrt(x*x +epsilon*epsilon);
	return g;
}

double pmfit(double r, double epsilon){
	double fr, tmp1, e;
	tmp1 = cosh(AA[1]*r);
	fr = r/pow(r*r+EPSILON*EPSILON,1.5L) 
		- ( 
				(1L/(r*r)*tanh(AA[1]*r)-AA[1]/(tmp1*tmp1*r) 
				 -2L*AA[2]/AA[0]*r*exp(AA[3]*r*r) 
				 -2L*AA[2]*AA[3]/AA[0]*r*r*r*exp(AA[3]*r*r) 
				 +AA[4]/AA[0]*(1L+AA[5]*r*r+AA[6]*r*r*r*r)*exp(AA[7]*r*r))
			);
	fr = fr/LL;

	return fr;
}


void i_force_spline(int nx, float Epsilon){
	double r,e,dr,rp2,rp1,rm1,rm2;
	int i,j,k;
	double fp, fpp, fr;
	double epsilon = Epsilon;
	LL = nx*nx;
	ran2nran = (double)RANGE/(double)NSPLINE;
	invran2nran = 1./ran2nran;
	for(i=0;i<3;i++) forcecorrectslope(0,i) = forcecorrectdiff(0,i) = 0.;
	dr = ran2nran*0.2L;
#ifdef _OPENMP
#pragma omp parallel for private(i,r,rm1,rm2,rp1,rp2,fp,fpp,fr)
#endif
	for(i=1;i<NSPLINE;i++){
		r = (double)i *ran2nran;
		rm1 = r - dr; rm2 = r - 2*dr; rp1 = r + dr; rp2 = r + 2*dr; 
		fp = (pmfit(rp1,epsilon) - pmfit(rm1,epsilon))/(dr+dr); 
		fpp = (pmfit(rp2,epsilon) - 2*pmfit(r,epsilon) + pmfit(rm2,epsilon))/(4*dr*dr); 
		fr = pmfit(r,epsilon); 
		forcecorrectdiff(i,0) = -(double)(fr/r); 
		forcecorrectdiff(i,1) = -(double)(0.5L*(fp*r-fr)/r/r/r); 
		forcecorrectdiff(i,2) = -(double)(0.5L*(fpp*r*r - 3*fp*r + 3*fr)/r/r/r/r/r);
	}
#ifdef _OPENMP
#pragma omp parallel for private(i,r,rm1,rm2,rp1,rp2,fp,fpp,fr)
#endif
	for(i=0;i<NSPLINE;i++){ 
		forcecorrectslope(i,0) = (forcecorrectdiff(i+1,0)-forcecorrectdiff(i,0))/ran2nran; 
		forcecorrectslope(i,1) = (forcecorrectdiff(i+1,1)-forcecorrectdiff(i,1))/ran2nran; 
		forcecorrectslope(i,2) = (forcecorrectdiff(i+1,2)-forcecorrectdiff(i,2))/ran2nran; 
	} 
	for(i=0;i<3;i++){ 
		forcecorrectdiff(NSPLINE-1,i)  =0 ; 
		forcecorrectslope(NSPLINE-1,i) =0 ; 
	}
	return;
}
void i_potent_spline(float Epsilon){ 
	double x; 
	double xstep; 
	double epsilon = Epsilon;
	int i,j,k; 
	double e,dx; 
	double xp2,xp1,xm1,xm2; 
	double  gp, gpp; 
	ran2nran=(double)RANGE/(double)NSPLINE; 
	invran2nran=1./ran2nran; 
	for(i=0;i<3;i++){ 
		forcecorrectslope(0,i) = forcecorrectdiff(0,i) = 0.; 
	} 
	dx = ran2nran;
#ifdef _OPENMP
#pragma omp parallel for private(i,x,xm1,xm2,xp1,xp2,gp,gpp)
#endif
	for(i=1;i<NSPLINE;i++){ 
		x = (double)(i) * ran2nran; 
		xm1 = x - dx; xm2 = x - 2*dx; xp1 = x + dx; xp2 = x + 2*dx; 
		gp = (sNewton(xp1,epsilon) - sNewton(xm1,epsilon))/(2*dx); 
		gpp = (sNewton(xp2,epsilon) - 2*sNewton(x,epsilon) + sNewton(xm2,epsilon))/(4.*dx*dx); 
		forcecorrectdiff(i,0) = sNewton(x,epsilon); 
		forcecorrectdiff(i,1) = 0.5L*gp/x; 
		forcecorrectdiff(i,2) = 0.5L*(gpp/x/x - gp/x/x/x);
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif 
	for(i=1;i<NSPLINE;i++){ 
		forcecorrectslope(i,0) = (forcecorrectdiff(i+1,0)-forcecorrectdiff(i,0))/dx; 
		forcecorrectslope(i,1) = (forcecorrectdiff(i+1,1)-forcecorrectdiff(i,1))/dx; 
		forcecorrectslope(i,2) = (forcecorrectdiff(i+1,2)-forcecorrectdiff(i,2))/dx; 
	} 
	return;
}


