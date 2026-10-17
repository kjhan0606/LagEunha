#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
typedef struct Halo{
        int np;
        float x,y,z;
        float vx,vy,vz;
}Halo;
Halo *halo;
#define NH 500
float unitmass;
#define pow3(a) ((a)*(a)*(a))

int main(int argc, char *argv[]){
    int i,j,k;
    int mp;

    FILE *fp;
    char infile[180];
        float size,hubble,npower,omep,omepl,bias,amax,astep,anow;
    int nx,nspace;


    halo=(Halo*) malloc(sizeof(Halo)*NH);
    fp=fopen(argv[1],"r");
        fread(&size,sizeof(float),1,fp);
        fread(&hubble,sizeof(float),1,fp);
        fread(&npower,sizeof(float),1,fp);
        fread(&omep,sizeof(float),1,fp);
        fread(&omepl,sizeof(float),1,fp);
        fread(&bias,sizeof(float),1,fp);
        fread(&nx,sizeof(int),1,fp);
        fread(&nspace,sizeof(int),1,fp);
        fread(&amax,sizeof(int),1,fp);
        fread(&astep,sizeof(int),1,fp);
        fread(&anow,sizeof(int),1,fp);
    unitmass = 2.7755E11*omep/pow3(nx/nspace/size);

    while((mp=fread(halo,sizeof(Halo),NH,fp))){
        for(i=0;i<mp;i++)
                printf("%g %g %g %g %g %g %g\n",
                halo[i].np*unitmass,halo[i].x,halo[i].y,halo[i].z,halo[i].vx,halo[i].vy,halo[i].vz);
    }

}
