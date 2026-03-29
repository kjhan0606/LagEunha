#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>


typedef union indxflag{
	size_t indx;
	unsigned char Flag[8];
} indxflag;



typedef struct bptype{
	int i,j,k;
	indxflag u4if;
	float x, y, z;
} bptype;

#define ENDIA_OFFSET 7


#define FLAG_RESET (0x00)
#define DMflag (1<<0)
#define SPHflag (1<<1)
#define STARflag (1<<2)
#define AGNflag (1<<3)
#define FoFflag (1<<4)
#define BoundaryGhostflag (1<<5)
#define NULL2 (1<<6)
#define NULL3 (1<<7)


#define CLEAR_FLAG(p)    ((p)->u4if.Flag[ENDIA_OFFSET] =   FLAG_RESET )
#define SET_FLAG(p, flag)  ((p)->u4if.Flag[ENDIA_OFFSET] |= flag)
#define UNSET_FLAG(p, flag)  ((p)->u4if.Flag[ENDIA_OFFSET] &= (~flag))
#define IS_FLAG(p, flag)  ((p)->u4if.Flag[ENDIA_OFFSET] &  flag)
#define TOGGLE_FLAG(p, flag)  ((p)->u4if.Flag[ENDIA_OFFSET] ^= flag)


#define PINDX(p) (((p)->u4if.indx <<8) >>8)

#define CHANGEINDX(p, a) do{\
	    unsigned char _b = (p)->u4if.Flag[ENDIA_OFFSET];\
	    (p)->u4if.indx = (a);\
	    (p)->u4if.Flag[ENDIA_OFFSET] = _b;\
}while(0)


int main(int argc, char **argv){
	bptype aa[2], *a;
	a = aa;


	(a+1)->u4if.indx = 128L;
	(a+1)->x = 0.48232;
	(a+1)->y = 1928.48232;
	(a+1)->z = -128370.48232;

	CHANGEINDX(a+1, 132383L);
	CLEAR_FLAG(a+1);
	SET_FLAG(a+1, DMflag);
	SET_FLAG(a+1, BoundaryGhostflag);
	SET_FLAG(a+1, NULL3);
	printf("%ld %d %d %d\n",PINDX(a+1), IS_FLAG(a+1, DMflag), IS_FLAG(a+1, BoundaryGhostflag), IS_FLAG(a+1, NULL3) );
	CHANGEINDX(a+1, 10);
	printf("%ld %d %d %d\n",PINDX(a+1), IS_FLAG(a+1, DMflag), IS_FLAG(a+1, BoundaryGhostflag), IS_FLAG(a+1, NULL3) );
}
