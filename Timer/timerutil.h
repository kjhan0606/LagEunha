#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>

float gettime();

#define TIMER_START(A) do{ cputime0[A] = gettime();} while(0)
#define TIMER_STOP(A)  do{ cputime1[A] = gettime();} while(0)
#define ELAPSED_TIME(A) (cputime1[A] - cputime0[A])

void starttimer(int );
void stoptimer(int );
float timelapse(int );
