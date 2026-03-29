/*
 * http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html
 */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>

#include "init.h"




int main(int argc, char **argv){
	int i,j,k;
	double Cs1,Cs2;
	double u2,u3,u4,u5;
	double P1,P2,P3, P4, P5;
	double rho2, rho3,rho1,rho5,rho4;
	double beta;
	double sodGamma;

	double uR, uL, x0;


	if(0){ /* standard setting */
		Lden = 1;
		LP = 1;
		uL = 0;
		Rden = 0.125;
		RP = 0.1;
		uR = 0;
		pgamma = 1.4;
		x0 = 0.5;
		xmin = 0;
		xmax = 1;
	}
	else { /* GIZMO setting */
		/*
		Lden = 1;
		LP = 1;
		Rden = 0.25;
		RP = 0.1795;
		pgamma = 1.4;
		xmin = -10.;
		xmax = 10.;
		*/


		uR = 0;
		uL = 0;
		x0 = 0.;
	}

	double time;

	time = atof(argv[1]);


	sodGamma = (pgamma-1.)/(pgamma+1L);
	beta = (pgamma-1)/(2*pgamma);
	uL = sqrt(pgamma*LP/Lden);
	uR = sqrt(pgamma*RP/Rden);

	u2 = uL;


	rho1 = Lden;
	rho5 = Rden;
	P1 = LP;
	P5 = RP;

	P3 = P5+0.01L;
	do{
		u4 = (P3-P5)*sqrt((1-sodGamma)/(Rden*(P3+sodGamma*P5)));
		u3 = (pow(P1,beta)-pow(P3,beta))*sqrt( (1-sodGamma*sodGamma)*pow(P1,1/pgamma)/(sodGamma*sodGamma*Lden));
		if(u3 > u4){
			P3 = P3 +0.00001L;
		}
		else {
			break;
		}
	}while(fabs((u3-u4)/u4)>1.e-5);

	double m2 = (pgamma-1)/(pgamma+1);
	double vpost = 2*sqrt(pgamma)/(pgamma-1)*(1-pow(P3,(pgamma-1)/2/pgamma));
	double Ppost = P3;
	double Pshock = P3;
	double rhopostoverRden = (Ppost/RP+m2)/(1+m2*Ppost/RP);
	double vshock = vpost*rhopostoverRden/(rhopostoverRden-1);
	double rhopost = rhopostoverRden*Rden;
	double rhomiddle = Lden*pow(Ppost/LP,1/pgamma);
	double cleft = uL;
	double x1 = x0 - cleft*time;
	double x3 = x0 + vpost*time;
	double x4 = x0 + vshock*time;

	double x2 = x0-time/m2*cleft*( pow(rhomiddle/Lden,(pgamma-1)/2) + m2-1);

	/*
	printf("%g %g %g %g %g\n",x1,x2,x0,x3,x4);
	printf("%g %g %g %g \n",rho1,rhomiddle,rhopost,rho5);
	*/

	for(i=0;i<512;i++){
		double x = xmin + (xmax-xmin)*(i+0.5)/512.l;
		double rho,P,vel;
		if(x < x1){
			rho = Lden;
			P = LP;
			vel = 0;
		}
		else if(x <x2){
			double Csound = m2*(x0-x)/time + (1-m2)*cleft;
			rho = Lden*pow(Csound/cleft,2/(pgamma-1));
			vel = (1-m2)*(cleft-(x0-x)/time);
			P = LP *pow(rho/Lden,pgamma);
		}
		else if(x < x3){
			rho = rhomiddle;
			vel = vpost;
			P = Ppost;
		}
		else if(x<x4){
			rho = rhopost;
			vel = vpost;
			P = Pshock;
		}
		else{
			rho = Rden;
			P = RP;
			vel = 0;
		}
		double poverrhogamma = P/pow(rho,pgamma);
		printf("%g %g %g %g %g\n",x,rho,P,vel, poverrhogamma);

	}
}

