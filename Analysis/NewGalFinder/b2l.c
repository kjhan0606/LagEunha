#include<stdio.h>
#ifdef ENDIANCHANGE
size_t Fwrite(void *a,size_t b,size_t c, FILE *fp){
	char *A;
	char t1,t2,t3,t4;
	size_t i,nmem;
	A = (char *)a;
	for(i=0;i<b*c;i+=4){
		t1 = A[i];
		t2 = A[i+1];
		t3 = A[i+2];
		t4 = A[i+3];
		A[i] =t4;
		A[i+1] =t3;
		A[i+2] =t2;
		A[i+3] =t1;
	}
	nmem =  fwrite(a,b,c,fp);
	A = (char *)a;
	for(i=0;i<b*c;i+=4){
		t1 = A[i];
		t2 = A[i+1];
		t3 = A[i+2];
		t4 = A[i+3];
		A[i] =t4;
		A[i+1] =t3;
		A[i+2] =t2;
		A[i+3] =t1;
	}
	return nmem;
}
size_t  Fread(void *a,size_t b,size_t c, FILE *fp){
	char *A;
	char t1,t2,t3,t4;
	size_t i,nmem;
	nmem = fread(a,b,c,fp);
	A = (char *)a;
	for(i=0;i<b*nmem;i+=4){
		t1 = A[i];
		t2 = A[i+1];
		t3 = A[i+2];
		t4 = A[i+3];
		A[i] =t4;
		A[i+1] =t3;
		A[i+2] =t2;
		A[i+3] =t1;
	}
	return nmem;
}
#else
size_t Fread(void *a,size_t b,size_t c, FILE *fp){
	return fread(a,b,c,fp);
}
size_t Fwrite(void *a,size_t b,size_t c, FILE *fp){
	return fwrite(a,b,c,fp);
}
#endif
