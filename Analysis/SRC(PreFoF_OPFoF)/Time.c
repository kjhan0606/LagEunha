#include<stdio.h>
#include<sys/time.h>
#include<sys/times.h>
#include<unistd.h>

struct timeval tv;
struct timezone tz;

double WTIME(){
	(void)gettimeofday(&tv,&tz);
	return tv.tv_sec + tv.tv_usec*1.E-6;
}
float WALLCLOCK(){
	static int iflag = 0;
	static double time0 = 0.;
	double timenow,etime;
	if(iflag == 0){
		time0 = WTIME();
		iflag = 1;
		return 0.;
	}
	else {
		timenow = WTIME();
		etime = timenow - time0;
		time0 = timenow;
		return etime;
	}
}
