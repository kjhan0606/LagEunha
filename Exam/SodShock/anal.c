/*
 * http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html
 */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>




int main(int argc, char **argv){
	int i,j,k;
	double Cs1,Cs2;
	double u2,u3,u4,u5,uL,uR;
	double PL, P1,P2,P3, P4, P5,PR;
	double rhoL,rhoR,rho2, rho3,rho1,rho5,rho4;
	double Gamma,gamma,beta;


	if(1){ /* standard setting */
		rhoL = 1;
		PL = 1;
		uL = 0;
		rhoR = 0.125;
		PR = 0.1;
		uR = 0;
		gamma = 1.4;
	}
	else { /* GIZMO setting */
		rhoL = 1;
		PL = 1;
		uL = 0;
		rhoR = 0.1795;
		PR = 0.25;
		uR = 0;
		gamma = 1.4;
	}

	double time;

	time = atof(argv[1]);


	Gamma = (gamma-1.)/(gamma+1);
	beta = (gamma-1)/(2*gamma);
	uL = sqrt(gamma*PL/rhoL);
	uR = sqrt(gamma*PR/rhoR);

	u2 = uL;


	rho1 = rhoL;
	rho5 = rhoR;
	P1 = PL;
	P5 = PR;

	P3 = P5+0.01;
	do{
		u4 = (P3-P5)*sqrt((1-Gamma)/(rhoR*(P3+Gamma*P5)));
		u3 = (pow(P1,beta)-pow(P3,beta))*sqrt( (1-Gamma*Gamma)*pow(P1,1/gamma)/(Gamma*Gamma*rhoL));
		if(u3 > u4){
			P3 = P3 +0.0001;
		}
		else {
			break;
		}
	}while(fabs(u3-u4)>1.e-5);

	double m2 = (gamma-1)/(gamma+1);
	double vpost = 2*sqrt(gamma)/(gamma-1)*(1-pow(P3,(gamma-1)/2/gamma));
	double Ppost = P3;
	double Pshock = P3;
	double rhopostoverrhoR = (Ppost/PR+m2)/(1+m2*Ppost/PR);
	double vshock = vpost*rhopostoverrhoR/(rhopostoverrhoR-1);
	double rhopost = rhopostoverrhoR*rhoR;
	double rhomiddle = rhoL*pow(Ppost/PL,1/gamma);
	double cleft = uL;
	double x0 = 0.5;
	double x1 = x0 - cleft*time;
	double x3 = x0 + vpost*time;
	double x4 = x0 + vshock*time;

	double x2 = x0-time/m2*cleft*( pow(rhomiddle/rhoL,(gamma-1)/2) + m2-1);

	/*
	printf("%g %g %g %g %g\n",x1,x2,x0,x3,x4);
	printf("%g %g %g %g \n",rho1,rhomiddle,rhopost,rho5);
	*/

	for(i=0;i<512;i++){
		double x = (i+0.5)/512.l;
		double rho,P,vel;
		if(x < x1){
			rho = rhoL;
			P = PL;
			vel = 0;
		}
		else if(x <x2){
			double Csound = m2*(x0-x)/time + (1-m2)*cleft;
			rho = rhoL*pow(Csound/cleft,2./(gamma-1));
			vel = (1-m2)*(cleft-(x0-x)/time);
			P = PL *pow(rho/rhoL,gamma);
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
			rho = rhoR;
			P = PR;
			vel = 0;
		}
		double poverrhogamma = P/pow(rho,gamma);
		printf("%g %g %g %g %g\n",x,rho,P,vel, poverrhogamma);

	}
}

