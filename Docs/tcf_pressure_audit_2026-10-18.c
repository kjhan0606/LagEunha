#include <math.h>
static double tcf_pstar(double rhoL, double pL, double uL,
		double rhoR, double pR, double uR, double Gamma)
{
	double cL = sqrt(Gamma*pL/rhoL), cR = sqrt(Gamma*pR/rhoR);
	if(uR - uL >= 2*(cL + cR)/(Gamma - 1)) return 0;
	double cmax = cL > cR ? cL : cR;
	double SL = (uL < uR ? uL : uR) - cmax, SR = (uL > uR ? uL : uR) + cmax;
	double aL = rhoL*(SL - uL), aR = rhoR*(SR - uR);
	double P = (pL*aR - pR*aL + aL*aR*(uR - uL))/(aR - aL);
	return P > 0 ? P : 0;
}
void batch(int N,const double*a,double*out){for(int i=0;i<N;i++){const double*p=a+11*i;out[i]=tcf_pstar(p[0],p[1],p[2],p[4],p[5],p[6],p[9]);}}
