#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>


typedef double postype;

int main(int argc, char **argv){


	postype eps2 = 0.1*0.1;
	postype x,y, r;
	postype vx,vy;
	postype time,Dtime;
	int i,j,k;
	time = 0;
	x = 1.5;
	y = 0;
	r = sqrt(x*x + y*y);
	postype vamp = r*pow(r*r+eps2,-0.75L);
	vx = -y/r*vamp; 
	vy =  x/r*vamp;
	do{
		r = sqrt(x*x + y*y);
		vamp = sqrt(vx*vx + vy*vy);
		postype ax = -x*pow(r*r+eps2,-1.5L);
		postype ay = -y*pow(r*r+eps2,-1.5L); 
		Dtime = 2*M_PI*r/vamp/100.;
		x += (vx + 0.5*ax*Dtime)*Dtime; 
		y += (vy + 0.5*ay*Dtime)*Dtime; 
		vx += ax *Dtime; 
		vy += ay *Dtime;
		time += Dtime;
		printf("%g %g\n", x,y);
	} while(time <120);

	return 0;

}
