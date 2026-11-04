#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../gfs_viscosity_implicit.h"
int main(int argc,char **argv){
    if(argc!=2)return 2;FILE *in=fopen(argv[1],"rb");uint64_t h[4];double dt;
    if(!in||fread(h,sizeof(h),1,in)!=1||h[0]!=UINT64_C(0x564953434641494c)||fread(&dt,sizeof(dt),1,in)!=1)return 2;
    int n=h[1],ne=h[2],dim=h[3];if(n<1||n>10000000||ne<0||ne>100000000||dim<1||dim>3)return 2;
    double *mass=calloc(n,sizeof(double)),*x=calloc(dim*n,sizeof(double)),*v=calloc(dim*n,sizeof(double)),*U=calloc(n,sizeof(double));
    GfsViscosityEdge *edges=calloc(ne,sizeof(*edges));
    if(!mass||!x||!v||!U||!edges||fread(mass,sizeof(double),n,in)!=(size_t)n||fread(x,sizeof(double),dim*n,in)!=(size_t)(dim*n)||fread(v,sizeof(double),dim*n,in)!=(size_t)(dim*n)||fread(U,sizeof(double),n,in)!=(size_t)n||fread(edges,sizeof(*edges),ne,in)!=(size_t)ne||fgetc(in)!=EOF)return 2;
    fclose(in);GfsViscosityWorkspace *w=gfs_viscosity_workspace_create();GfsViscosityStats s={0};
    int rc=gfs_viscosity_be(w,dim,n,ne,edges,mass,x,v,U,dt,&s);
    printf("{\"return_code\":%d,\"cells\":%d,\"edges\":%d,\"active_iterations\":%d,\"CG_iterations\":%d,\"energy_defect\":%.17g,\"equation_residual\":%.17g,\"heat\":%.17g}\n",rc,n,ne,s.active_iterations,s.iterations,s.energy_defect,s.relative_equation_residual,s.heat);
    gfs_viscosity_workspace_free(w);free(mass);free(x);free(v);free(U);free(edges);return rc;
}
