#include<stdio.h>
#include<stdlib.h>

void printBits(char num){ 
	int size = sizeof(char);
	char maxPow = 1<<(size*8-1);
	printf("MAX POW : %d\n",maxPow);
	int i=0,j;
	for(;i<size;++i){ 
		for(;i<size*8;++i){ 
			printf("%u ",num&maxPow ? 1 : 0);
			num = num<<1; 
		} 
	}
	printf("\n");
}

int main(){
	char flag[8];
	int i;

	for(i=0;i<8;i++) flag[i] = 255;
	for(i=0;i<8;i++) printBits(flag[i]);
}
