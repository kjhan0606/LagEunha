#include <stdio.h>
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>


float cputime0[1024];
float cputime1[1024];

struct timeval tv;
float gettime()
{
	static int startflag=1;
	static double tsecs0, tsecs1;
	
	if( startflag ) {
		(void) gettimeofday(&tv,0);
		tsecs0 = tv.tv_sec + tv.tv_usec*1.0e-6;
		startflag = 0;
	}
	(void) gettimeofday(&tv,0);
	tsecs1 = tv.tv_sec + tv.tv_usec*1.0e-6;

	return ((float) (tsecs1-tsecs0));

}
void starttimer(int A){ cputime0[A] = gettime(); }
void stoptimer(int A){ cputime1[A] = gettime(); }
float timelapse(int A){ return (cputime1[A] - cputime0[A]); }


