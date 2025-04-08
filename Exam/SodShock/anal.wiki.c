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


	rhoL = 1;
	PL = 1;
	uL = 0;

	rhoR = 0.125;
	PR = 0.1;
	uR = 0;

	gamma = 1.4;

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
	rho4 = rho5*(P4+Gamma*P5)/(P5+Gamma*P4);


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

	P4 = P3;

	u5  = u3 - (P3-P5)/sqrt(rho5/2*( (gamma+1)*P3 + (gamma-1)*P5 ));


	printf("%g %g %g %g\n", u2*time, u3*time,u4*time,u5*time);

	
	rho3 = rho1*pow(P3/P1,1/gamma);
	double x0 = 0.5;

	u2 = 2/(gamma+1)*(uL + (x - x0)/time);
	rho2 = rho1*pow(1-(gamma-1)/2*u2/uL,2/(gamma-1));
	P2 = P1*pow(1-(gamma-1)/2*u2/uL,2*gamma/(gamma-1));
	/*



	double x1 = x0 - cleft*time;
	double x2 = x0 - vel*time;
	double x3 = x0 + vpost*time;
	double x4 = x0 + vshock*time;

	for(i=0;i<128;i++){
		double x = (i+0.5)/128.l;




		double csound = gamma*gamma*(x0-x)/time + (1-gamma*gamma)*cleft;
		double vel = (1-gamma*gamma)*(cleft - (x0-x)/time);
		double rho = rhoL*pow(csound/cleft,2./(gamma-1));
		double P = PL*pow(rho/rhoL,gamma);



		printf("%g %g %g %g %g\n",x,rho,vel,P,csound);
	}
	*/
}

