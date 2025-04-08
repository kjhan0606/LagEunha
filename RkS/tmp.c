#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>

typedef struct hello {
    float (*someFunction)();
} hello;

float foo(long i) {
    return sqrt(i*1.0000001L+0.4);
}
float foo1(long i) {
    return -sqrt(i+1);
}

/*
hello Hello(int (*foo)()) {
    struct hello aHello;
    aHello.someFunction = foo;
    return aHello;
}
*/

#define np 100000000L

int main()
{
	long i;
	double aa[np];

    struct hello aHello;
	for(i=0;i<np;i++){
		aHello.someFunction=foo;
		aa[i] = aHello.someFunction(i);
		aHello.someFunction=foo1;
		aa[i] += aHello.someFunction(i);
	}
	double bb;
	for(i=0;i<np;i++){
		bb += aa[i];
	}
	printf("total num = %g\n",bb);

    return 0;
} 
